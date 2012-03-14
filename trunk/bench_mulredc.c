#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h> /* for LONG_MAX */
#include <assert.h>
#include <time.h>
#include <string.h>

#if TIME_WITH_SYS_TIME
# include <sys/time.h>
# include <time.h>
#else
# if HAVE_SYS_TIME_H
#  include <sys/time.h>
# else
#  include <time.h>
# endif
#endif

#define LOOPCOUNT 10000000UL
#define MAXSIZE 20

int tune_mul[MAXSIZE+1], tune_sqr[MAXSIZE+1];

#include <gmp.h>
#ifdef USE_ASM_REDC
#include "mulredc.h"
#endif
#include "mpmod.h"

#ifdef HAVE___GMPN_REDC_1
#ifndef __gmpn_redc_1
void __gmpn_redc_1 (mp_ptr, mp_ptr, mp_srcptr, mp_size_t, mp_limb_t);
#endif
#endif

#ifdef HAVE___GMPN_REDC_2
#ifndef __gmpn_redc_2
void __gmpn_redc_2 (mp_ptr, mp_ptr, mp_srcptr, mp_size_t, mp_srcptr);
#endif
#endif

#ifdef HAVE___GMPN_REDC_N
#ifndef __gmpn_redc_N
void __gmpn_redc_n (mp_ptr, mp_ptr, mp_srcptr, mp_size_t, mp_srcptr);
#endif
#endif

/* cputime () gives the elapsed time in milliseconds */

#if defined (_WIN32)
/* First case - GetProcessTimes () is the only known way of getting process
 * time (as opposed to calendar time) under mingw32 */

#include <windows.h>

long
cputime ()
{
  FILETIME lpCreationTime, lpExitTime, lpKernelTime, lpUserTime;
  ULARGE_INTEGER n;

  HANDLE hProcess = GetCurrentProcess();
  
  GetProcessTimes (hProcess, &lpCreationTime, &lpExitTime, &lpKernelTime,
      &lpUserTime);

  /* copy FILETIME to a ULARGE_INTEGER as recommended by MSDN docs */
  n.u.LowPart = lpUserTime.dwLowDateTime;
  n.u.HighPart = lpUserTime.dwHighDateTime;

  /* lpUserTime is in units of 100 ns. Return time in milliseconds */
  return (long) (n.QuadPart / 10000);
}

#elif defined (HAVE_GETRUSAGE)
/* Next case: getrusage () has higher resolution than clock () and so is
   preferred. */

#ifdef HAVE_SYS_TYPES_H
# include <sys/types.h>
#endif
#ifdef HAVE_SYS_RESOURCE_H
# include <sys/resource.h>
#endif

long
cputime ()
{
  struct rusage rus;

  getrusage (RUSAGE_SELF, &rus);
  /* This overflows a 32 bit signed int after 2147483s = 24.85 days */
  return rus.ru_utime.tv_sec * 1000L + rus.ru_utime.tv_usec / 1000L;
}

#else
/* Resort to clock (), which on some systems may return calendar time. */

long
cputime ()
{
  /* Return time in milliseconds */
  return (long) (clock () * (1000. / (double) CLOCKS_PER_SEC));
}

#endif /* defining cputime () */

void mp_print(mp_limb_t *x, int N) {
  int i;
  for (i = 0; i < N-1; ++i)
    gmp_printf("%Nd + W*(", x + i, 1);
  gmp_printf("%Nd", x + (N-1), 1);
  for (i = 0; i < N-1; ++i)
    printf(")");
  printf("\n");
}

static void
ecm_redc_1_svoboda (mp_ptr rp, mp_ptr tmp, mp_srcptr np, mp_size_t nn,
                    mp_limb_t invm, mp_srcptr sp)
{
  mp_size_t j;
  mp_limb_t t0, cy;

  /* instead of adding {np, nn} * (invm * tmp[0] mod B), we add
     {sp, nn} * tmp[0], where {np, nn} * invm = B * {sp, nn} - 1 */
  for (j = 0; j < nn - 1; j++, tmp++)
    rp[j + 1] = mpn_addmul_1 (tmp + 1, sp, nn, tmp[0]);
  /* for the last step, we reduce with {np, nn} */
  t0 = mpn_addmul_1 (tmp, np, nn, tmp[0] * invm);
  tmp ++;

  rp[0] = tmp[0];
  cy = mpn_add_n (rp + 1, rp + 1, tmp + 1, nn - 1);
  rp[nn-1] += t0;
  cy += rp[nn-1] < t0;
  if (cy != 0)
    mpn_sub_n (rp, rp, np, nn); /* a borrow should always occur here */
}

void bench(mp_size_t N)
{
  mp_limb_t *x, *y, *z, *m, *invm, *tmp, *svoboda1;
  unsigned long i;
  unsigned long iter;
  long tmul, tsqr, tredc_1, t_mulredc_1, tsvoboda1 = 0, t_sqrredc_1;
  long tmul_best = LONG_MAX, tsqr_best = LONG_MAX, tredc_best = LONG_MAX;
  mpz_t M, B;
#ifdef USE_ASM_REDC
  long t2;
#endif
#ifdef HAVE_NATIVE_MULREDC1_N
  long t3 = 0;
#endif
#ifdef HAVE___GMPN_REDC_2
  long tredc_2, t_mulredc_2, t_sqrredc_2;
#endif
#ifdef HAVE___GMPN_REDC_N
  long tredc_n, t_mulredc_n, t_sqrredc_n;
#endif
  
  x = (mp_limb_t *) malloc(N*sizeof(mp_limb_t));
  y = (mp_limb_t *) malloc(N*sizeof(mp_limb_t));
  z = (mp_limb_t *) malloc((2*N)*sizeof(mp_limb_t));
  m = (mp_limb_t *) malloc(N*sizeof(mp_limb_t));
  tmp = (mp_limb_t *) malloc((2*N+2)*sizeof(mp_limb_t));
  invm = (mp_limb_t *) malloc(N*sizeof(mp_limb_t));
  svoboda1 = (mp_limb_t *) malloc(N*sizeof(mp_limb_t));
 
  mpn_random(m, N);
  m[0] |= 1UL;
  if (m[N-1] == 0) 
    m[N-1] = 1UL;

  mpz_init (M);
  mpz_init (B);
  mpz_set_ui (M, m[1]);
  mpz_mul_2exp (M, M, GMP_NUMB_BITS);
  mpz_add_ui (M, M, m[0]);
  mpz_set_ui (B, 1);
  mpz_mul_2exp (B, B, 2 * GMP_NUMB_BITS);
  mpz_invert (M, M, B);
  mpz_sub (M, B, M);

  for (i = 0; i < (unsigned) N; i++)
    invm[i] = mpz_getlimbn(M, i);

  tmp[N] = mpn_mul_1 (tmp, m, N, invm[0]); /* {tmp,N+1} should be = -1 mod B */
  mpn_add_1 (tmp, tmp, N + 1, 1); /* now = 0 mod B */
  mpn_copyi (svoboda1, tmp + 1, N);

  mpz_clear (M);
  mpz_clear (B);

  mpn_random(x, N);
  mpn_random(y, N);

  /* we set 'iter' to get about 100ms for each test */
  tmul = cputime();
  i = 0;
  iter = 1;
  do
    {
      iter = 2 * iter;
      for (; i < iter; i++)
        mpn_mul_n (tmp, x, y, N);
    }
  while (cputime() - tmul < 100);
  iter = (long) (((double) iter * 100.0) / (double) (cputime() - tmul));

  tmul = cputime();
  for (i = 0; i < iter; ++i)
    mpn_mul_n(tmp, x, y, N);
  tmul = cputime()-tmul;

  tsqr = cputime();
  for (i = 0; i < iter; ++i)
    mpn_sqr (tmp, x, N);
  tsqr = cputime()-tsqr;

#ifdef HAVE___GMPN_REDC_1
  mpn_mul_n(tmp, x, y, N);
  tredc_1 = cputime();
  for (i = 0; i < iter; ++i)
    __gmpn_redc_1 (z, tmp, m, N, invm[0]);
  tredc_1 = cputime()-tredc_1;
  if (tredc_1 < tredc_best)
    tredc_best = tredc_1;
#endif

  if (N > 1) /* Svoboda only works for N > 1 */
    {
      mpn_mul_n(tmp, x, y, N);
      tsvoboda1 = cputime();
      for (i = 0; i < iter; ++i)
        ecm_redc_1_svoboda (z, tmp, m, N, invm[0], svoboda1);
      tsvoboda1 = cputime()-tsvoboda1;
      if (tsvoboda1 < tredc_best)
        tredc_best = tsvoboda1;
    }

#ifdef HAVE___GMPN_REDC_2
  mpn_mul_n(tmp, x, y, N);
  tredc_2 = cputime();
  for (i = 0; i < iter; ++i)
    __gmpn_redc_2 (z, tmp, m, N, invm);
  tredc_2 = cputime()-tredc_2;
  if (tredc_2 < tredc_best)
    tredc_best = tredc_2;
#endif

#ifdef HAVE___GMPN_REDC_N
  mpn_mul_n(tmp, x, y, N);
  tredc_n = cputime();
  for (i = 0; i < iter; ++i)
    __gmpn_redc_n (z, tmp, m, N, invm);
  tredc_n = cputime()-tredc_n;
  if (tredc_n < tredc_best)
    tredc_best = tredc_n;
#endif

#ifdef USE_ASM_REDC
  /* Mixed mul and redc */
  t2 = cputime();
  switch (N) {
   case 1:
    for (i=0; i < iter; ++i) {
      mulredc1(z, x[0], y[0], m[0], invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 2:
    for (i=0; i < iter; ++i) {
      mulredc2(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 3:
    for (i=0; i < iter; ++i) {
      mulredc3(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 4:
    for (i=0; i < iter; ++i) {
      mulredc4(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 5:
    for (i=0; i < iter; ++i) {
      mulredc5(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 6:
    for (i=0; i < iter; ++i) {
      mulredc6(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 7:
    for (i=0; i < iter; ++i) {
      mulredc7(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 8:
    for (i=0; i < iter; ++i) {
      mulredc8(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 9:
    for (i=0; i < iter; ++i) {
      mulredc9(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 10:
    for (i=0; i < iter; ++i) {
      mulredc10(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 11:
    for (i=0; i < iter; ++i) {
      mulredc11(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 12:
    for (i=0; i < iter; ++i) {
      mulredc12(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 13:
    for (i=0; i < iter; ++i) {
      mulredc13(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 14:
    for (i=0; i < iter; ++i) {
      mulredc14(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 15:
    for (i=0; i < iter; ++i) {
      mulredc15(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 16:
    for (i=0; i < iter; ++i) {
      mulredc16(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 17:
    for (i=0; i < iter; ++i) {
      mulredc17(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 18:
    for (i=0; i < iter; ++i) {
      mulredc18(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 19:
    for (i=0; i < iter; ++i) {
      mulredc19(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 20:
    for (i=0; i < iter; ++i) {
      mulredc20(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   default:
    for (i=0; i < iter; ++i) {
      mulredc20(z, x, y, m,  invm[0]);
      x[0] += tmp[0];
    }
  }
  t2 = cputime()-t2;
  if (t2 < tmul_best)
    {
      tmul_best = t2;
      tune_mul[N] = MPMOD_MULREDC;
    }
  if (t2 < tsqr_best)
    {
      tsqr_best = t2;
      tune_sqr[N] = MPMOD_MULREDC;
    }
#endif
  
  /* Mul followed by mpn_redc_1 */
#ifdef HAVE___GMPN_REDC_1
  t_mulredc_1 = cputime();
  for (i = 0; i < iter; ++i) {
    mpn_mul_n(tmp, x, y, N);
    __gmpn_redc_1 (z, tmp, m, N, invm[0]);
    x[0] += tmp[0];
  }
  t_mulredc_1 = cputime()-t_mulredc_1;
  if (t_mulredc_1 < tmul_best)
    {
      tune_mul[N] = MPMOD_MUL_REDC1;
      tmul_best = t_mulredc_1;
    }
#endif
  
  /* Mul followed by mpn_redc_2 */
#ifdef HAVE___GMPN_REDC_2
  t_mulredc_2 = cputime();
  for (i = 0; i < iter; ++i) {
    mpn_mul_n(tmp, x, y, N);
    __gmpn_redc_2 (z, tmp, m, N, invm);
    x[0] += tmp[0];
  }
  t_mulredc_2 = cputime()-t_mulredc_2;
  if (t_mulredc_2 < tmul_best)
    {
      tune_mul[N] = MPMOD_MUL_REDC2;
      tmul_best = t_mulredc_2;
    }
#endif
  
  /* Mul followed by mpn_redc_n */
#ifdef HAVE___GMPN_REDC_N
  t_mulredc_n = cputime();
  for (i = 0; i < iter; ++i)
    {
      mpn_mul_n (tmp, x, y, N);
      __gmpn_redc_n (z, tmp, m, N, invm);
    }
  t_mulredc_n = cputime()-t_mulredc_n;
  if (t_mulredc_n < tmul_best)
    {
      tune_mul[N] = MPMOD_MUL_REDCN;
      tmul_best = t_mulredc_n;
    }
#endif
  
  /* Sqr followed by mpn_redc_1 */
#ifdef HAVE___GMPN_REDC_1
  t_sqrredc_1 = cputime();
  for (i = 0; i < iter; ++i) {
    mpn_sqr(tmp, x, N);
    __gmpn_redc_1 (z, tmp, m, N, invm[0]);
    x[0] += tmp[0];
  }
  t_sqrredc_1 = cputime()-t_sqrredc_1;
  if (t_sqrredc_1 < tsqr_best)
    {
      tune_sqr[N] = MPMOD_MUL_REDC1;
      tsqr_best = t_sqrredc_1;
    }
#endif
  
  /* Sqr followed by mpn_redc_2 */
#ifdef HAVE___GMPN_REDC_2
  t_sqrredc_2 = cputime();
  for (i = 0; i < iter; ++i) {
    mpn_sqr(tmp, x, N);
    __gmpn_redc_2 (z, tmp, m, N, invm);
    x[0] += tmp[0];
  }
  t_sqrredc_2 = cputime()-t_sqrredc_2;
  if (t_sqrredc_2 < tsqr_best)
    {
      tune_sqr[N] = MPMOD_MUL_REDC2;
      tsqr_best = t_sqrredc_2;
    }
#endif
  
  /* Sqr followed by mpn_redc_n */
#ifdef HAVE___GMPN_REDC_N
  t_sqrredc_n = cputime();
  for (i = 0; i < iter; ++i)
    {
      mpn_sqr (tmp, x, N);
      __gmpn_redc_n (z, tmp, m, N, invm);
    }
  t_sqrredc_n = cputime()-t_sqrredc_n;
  if (t_sqrredc_n < tsqr_best)
    {
      tune_sqr[N] = MPMOD_MUL_REDCN;
      tsqr_best = t_sqrredc_n;
    }
#endif
  
#ifdef HAVE_NATIVE_MULREDC1_N
  /* mulredc1 */
  t3 = cputime();
  switch (N) {
   case 1:
    for (i=0; i<LOOPCOUNT; ++i) {
      mulredc1(z, x[0], y[0], m[0], invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 2:
    for (i=0; i<LOOPCOUNT; ++i) {
      mulredc1_2(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 3:
    for (i=0; i<LOOPCOUNT; ++i) {
      mulredc1_3(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 4:
    for (i=0; i<LOOPCOUNT; ++i) {
      mulredc1_4(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 5:
    for (i=0; i<LOOPCOUNT; ++i) {
      mulredc1_5(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 6:
    for (i=0; i<LOOPCOUNT; ++i) {
      mulredc1_6(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 7:
    for (i=0; i<LOOPCOUNT; ++i) {
      mulredc1_7(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 8:
    for (i=0; i<LOOPCOUNT; ++i) {
      mulredc1_8(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 9:
    for (i=0; i<LOOPCOUNT; ++i) {
      mulredc1_9(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 10:
    for (i=0; i<LOOPCOUNT; ++i) {
      mulredc1_10(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 11:
    for (i=0; i<LOOPCOUNT; ++i) {
      mulredc1_11(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 12:
    for (i=0; i<LOOPCOUNT; ++i) {
      mulredc1_12(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 13:
    for (i=0; i<LOOPCOUNT; ++i) {
      mulredc1_13(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 14:
    for (i=0; i<LOOPCOUNT; ++i) {
      mulredc1_14(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 15:
    for (i=0; i<LOOPCOUNT; ++i) {
      mulredc1_15(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 16:
    for (i=0; i<LOOPCOUNT; ++i) {
      mulredc1_16(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 17:
    for (i=0; i<LOOPCOUNT; ++i) {
      mulredc1_17(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 18:
    for (i=0; i<LOOPCOUNT; ++i) {
      mulredc1_18(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 19:
    for (i=0; i<LOOPCOUNT; ++i) {
      mulredc1_19(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 20:
    for (i=0; i<LOOPCOUNT; ++i) {
      mulredc1_20(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   default: ;
  }
  t3 = cputime() - t3;
#endif /* ifdef HAVE_NATIVE_MULREDC1_N */

  fprintf (stderr, "******************\n");
  fprintf (stderr, "Time in microseconds per call, size=%lu (iter=%lu):\n",
           N, iter);

  /* basic operations */
  fprintf (stderr, "mpn_mul_n  = %.3f\n",
           (double) tmul * 1e3 / (double) iter);
  fprintf (stderr, "mpn_sqr    = %.3f\n",
           (double) tsqr * 1e3 / (double) iter);
#ifdef HAVE___GMPN_REDC_1
  fprintf (stderr, "mpn_redc_1 = %.3f",
           (double) tredc_1 * 1e3 / (double) iter);
  if (tredc_1 == tredc_best)
    fprintf (stderr, " *");
  fprintf (stderr, "\n");
#endif
  if (N > 1)
    {
      fprintf (stderr, "svoboda1   = %.3f",
               (double) tsvoboda1 * 1e3 / (double) iter);
      if (tsvoboda1 == tredc_best)
        fprintf (stderr, " *");
      fprintf (stderr, "\n");
    }
#ifdef HAVE___GMPN_REDC_2
  fprintf (stderr, "mpn_redc_2 = %.3f",
           (double) tredc_2 * 1e3 / (double) iter);
  if (tredc_2 == tredc_best)
    fprintf (stderr, " *");
  fprintf (stderr, "\n");
#endif
#ifdef HAVE___GMPN_REDC_N
  fprintf (stderr, "mpn_redc_n = %.3f",
           (double) tredc_n * 1e3 / (double) iter);
  if (tredc_n == tredc_best)
    fprintf (stderr, " *");
  fprintf (stderr, "\n");
#endif

  fprintf (stderr, "\n");

  /* modular multiplication */
#ifdef USE_ASM_REDC
  fprintf (stderr, "mulredc    = %.3f",
           (double) t2 * 1e3 / (double) iter);
  if (tmul_best == t2)
    fprintf (stderr, " *");
  fprintf (stderr, "\n");
#endif
#ifdef HAVE___GMPN_REDC_1
  fprintf (stderr, "mul+redc_1 = %.3f",
           (double) t_mulredc_1 * 1e3 / (double) iter);
  if (tmul_best == t_mulredc_1)
    fprintf (stderr, " *");
  fprintf (stderr, "\n");
#endif
#ifdef HAVE___GMPN_REDC_2
  fprintf (stderr, "mul+redc_2 = %.3f",
           (double) t_mulredc_2 * 1e3 / (double) iter);
  if (tmul_best == t_mulredc_2)
    fprintf (stderr, " *");
  fprintf (stderr, "\n");
#endif
#ifdef HAVE___GMPN_REDC_N
  fprintf (stderr, "mul+redc_n = %.3f",
           (double) t_mulredc_n * 1e3 / (double) iter);
  if (tmul_best == t_mulredc_n)
    fprintf (stderr, " *");
  fprintf (stderr, "\n");
#endif

  fprintf (stderr, "\n");

  /* modular squaring */
#ifdef USE_ASM_REDC
  fprintf (stderr, "mulredc    = %.3f",
           (double) t2 * 1e3 / (double) iter);
  if (tsqr_best == t2)
    fprintf (stderr, " *");
  fprintf (stderr, "\n");
#endif
#ifdef HAVE___GMPN_REDC_1
  fprintf (stderr, "sqr+redc_1 = %.3f",
           (double) t_sqrredc_1 * 1e3 / (double) iter);
  if (tsqr_best == t_sqrredc_1)
    fprintf (stderr, " *");
  fprintf (stderr, "\n");
#endif
#ifdef HAVE___GMPN_REDC_2
  fprintf (stderr, "sqr+redc_2 = %.3f",
           (double) t_sqrredc_2 * 1e3 / (double) iter);
  if (tsqr_best == t_sqrredc_2)
    fprintf (stderr, " *");
  fprintf (stderr, "\n");
#endif
#ifdef HAVE___GMPN_REDC_N
  fprintf (stderr, "sqr+redc_n = %.3f",
           (double) t_sqrredc_n * 1e3 / (double) iter);
  if (tsqr_best == t_sqrredc_n)
    fprintf (stderr, " *");
  fprintf (stderr, "\n");
#endif

#ifdef HAVE_NATIVE_MULREDC1_N
  /* multiplication of n limbs by one limb */
  fprintf (stderr, "mulredc1   = %.3f\n",
           (double) t3 * 1e3 / (double) LOOPCOUNT);
#endif
  fflush (stderr);
  
  free (tmp);
  free (x);
  free (y);
  free (z);
  free (m);
  free (invm);
  free (svoboda1);
}
  

int main(int argc, char** argv)
{
  int i;
  int minsize = 1, maxsize = MAXSIZE;

  if (argc > 1)
    minsize = atoi (argv[1]);
  if (argc > 2)
    maxsize = atoi (argv[2]);
  
  for (i = minsize; i <= maxsize; ++i)
    bench(i);

  printf ("/* 0:mulredc 1:mul+redc_1 2:mul+redc_2 3:mul+redc_n */\n");
  printf ("#define TUNE_MULREDC_TABLE {0");
  for (i = 1; i <= maxsize; i++)
    printf (",%d", tune_mul[i]);
  printf ("}\n");
  printf ("/* 0:mulredc 1:sqr+redc_1 2:sqr+redc_2 3:sqr+redc_n */\n");
  printf ("#define TUNE_SQRREDC_TABLE {0");
  for (i = 1; i <= maxsize; i++)
    printf (",%d", tune_sqr[i]);
  printf ("}\n");
  fflush (stdout);

  return 0;
}
