#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#ifndef _MSC_VER
# include <sys/time.h>
# include <sys/times.h>
# include <sys/param.h>
# include <sys/resource.h>
# include <unistd.h>
#else
#include "getrusage.h"
#endif
# ifndef HZ
#  ifdef CLK_TCK
#   define HZ CLK_TCK
#  else
#   define HZ 100
#  endif
# endif

#define LOOPCOUNT 10000000UL
#define MAXSIZE 20

int tune_mul[MAXSIZE+1], tune_sqr[MAXSIZE+1];

#include <gmp.h>
#include "mulredc.h"
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

double CPUTime()
{
  double ret;
#if HAVE_GETRUSAGE
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  ret = (double) usage.ru_utime.tv_sec + 
        (double) usage.ru_utime.tv_usec / 1000000.;
#elif defined(USE_CLOCK)
  clock_t t;

  t = clock();
  ret = (double)t / (double)CLOCKS_PER_SEC;
#else
  struct tms t; 

  times(&t);
  ret = t.tms_utime * 1. / HZ;
#endif
  return ret;
}

void mp_print(mp_limb_t *x, int N) {
  int i;
  for (i = 0; i < N-1; ++i)
    printf("%lu + W*(", x[i]);
  printf("%lu", x[N-1]);
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
  mp_limb_t *x, *y, *z, *m, *invm, cy, *tmp, *svoboda1;
  unsigned long i;
  const unsigned long iter = LOOPCOUNT/N;
  double t2, t3 = 0., tmul, tsqr, tredc_1, t_mulredc_1;
  double t_sqrredc_1, tmul_best, tsqr_best, tsvoboda1;
  mpz_t M, B;
#ifdef HAVE___GMPN_REDC_2
  double tredc_2, t_mulredc_2, t_sqrredc_2;
#endif
#ifdef HAVE___GMPN_REDC_N
  double tredc_n, t_mulredc_n, t_sqrredc_n;
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

  tmul = CPUTime();
  for (i = 0; i < iter; ++i)
    mpn_mul_n(tmp, x, y, N);
  tmul = CPUTime()-tmul;

  tsqr = CPUTime();
  for (i = 0; i < iter; ++i)
    mpn_sqr (tmp, x, N);
  tsqr = CPUTime()-tsqr;

#ifdef HAVE___GMPN_REDC_1
  mpn_mul_n(tmp, x, y, N);
  tredc_1 = CPUTime();
  for (i = 0; i < iter; ++i)
    __gmpn_redc_1 (z, tmp, m, N, invm[0]);
  tredc_1 = CPUTime()-tredc_1;
#endif

  if (N > 1) /* Svoboda only works for N > 1 */
    {
      mpn_mul_n(tmp, x, y, N);
      tsvoboda1 = CPUTime();
      for (i = 0; i < iter; ++i)
        ecm_redc_1_svoboda (z, tmp, m, N, invm[0], svoboda1);
      tsvoboda1 = CPUTime()-tsvoboda1;
    }

#ifdef HAVE___GMPN_REDC_2
  mpn_mul_n(tmp, x, y, N);
  tredc_2 = CPUTime();
  for (i = 0; i < iter; ++i)
    __gmpn_redc_2 (z, tmp, m, N, invm);
  tredc_2 = CPUTime()-tredc_2;
#endif

#ifdef HAVE___GMPN_REDC_N
  mpn_mul_n(tmp, x, y, N);
  tredc_n = CPUTime();
  for (i = 0; i < iter; ++i)
    __gmpn_redc_n (z, tmp, m, N, invm);
  tredc_n = CPUTime()-tredc_n;
#endif

  /* Mixed mul and redc */
  t2 = CPUTime();
  switch (N) {
   case 1:
    for (i=0; i < iter; ++i) {
      cy += mulredc1(z, x[0], y[0], m[0], invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 2:
    for (i=0; i < iter; ++i) {
      cy += mulredc2(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 3:
    for (i=0; i < iter; ++i) {
      cy += mulredc3(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 4:
    for (i=0; i < iter; ++i) {
      cy += mulredc4(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 5:
    for (i=0; i < iter; ++i) {
      cy += mulredc5(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 6:
    for (i=0; i < iter; ++i) {
      cy += mulredc6(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 7:
    for (i=0; i < iter; ++i) {
      cy += mulredc7(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 8:
    for (i=0; i < iter; ++i) {
      cy += mulredc8(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 9:
    for (i=0; i < iter; ++i) {
      cy += mulredc9(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 10:
    for (i=0; i < iter; ++i) {
      cy += mulredc10(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 11:
    for (i=0; i < iter; ++i) {
      cy += mulredc11(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 12:
    for (i=0; i < iter; ++i) {
      cy += mulredc12(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 13:
    for (i=0; i < iter; ++i) {
      cy += mulredc13(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 14:
    for (i=0; i < iter; ++i) {
      cy += mulredc14(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 15:
    for (i=0; i < iter; ++i) {
      cy += mulredc15(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 16:
    for (i=0; i < iter; ++i) {
      cy += mulredc16(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 17:
    for (i=0; i < iter; ++i) {
      cy += mulredc17(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 18:
    for (i=0; i < iter; ++i) {
      cy += mulredc18(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 19:
    for (i=0; i < iter; ++i) {
      cy += mulredc19(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 20:
    for (i=0; i < iter; ++i) {
      cy += mulredc20(z, x, y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   default:
    for (i=0; i < iter; ++i) {
      cy += mulredc20(z, x, y, m,  invm[0]);
      x[0] += tmp[0];
    }
  }
  t2 = CPUTime()-t2;
  tmul_best = t2;
  tune_mul[N] = MPMOD_MULREDC;
  tsqr_best = t2;
  tune_sqr[N] = MPMOD_MULREDC;
  
  /* Mul followed by mpn_redc_1 */
#ifdef HAVE___GMPN_REDC_1
  t_mulredc_1 = CPUTime();
  for (i = 0; i < iter; ++i) {
    mpn_mul_n(tmp, x, y, N);
    __gmpn_redc_1 (z, tmp, m, N, invm[0]);
    x[0] += tmp[0];
  }
  t_mulredc_1 = CPUTime()-t_mulredc_1;
  if (t_mulredc_1 < tmul_best)
    {
      tune_mul[N] = MPMOD_MUL_REDC1;
      tmul_best = t_mulredc_1;
    }
#endif
  
  /* Mul followed by mpn_redc_2 */
#ifdef HAVE___GMPN_REDC_2
  t_mulredc_2 = CPUTime();
  for (i = 0; i < iter; ++i) {
    mpn_mul_n(tmp, x, y, N);
    __gmpn_redc_2 (z, tmp, m, N, invm);
    x[0] += tmp[0];
  }
  t_mulredc_2 = CPUTime()-t_mulredc_2;
  if (t_mulredc_2 < tmul_best)
    {
      tune_mul[N] = MPMOD_MUL_REDC2;
      tmul_best = t_mulredc_2;
    }
#endif
  
  /* Mul followed by mpn_redc_n */
#ifdef HAVE___GMPN_REDC_N
  t_mulredc_n = CPUTime();
  for (i = 0; i < iter; ++i)
    {
      mpn_mul_n (tmp, x, y, N);
      __gmpn_redc_n (z, tmp, m, N, invm);
    }
  t_mulredc_n = CPUTime()-t_mulredc_n;
  if (t_mulredc_n < tmul_best)
    {
      tune_mul[N] = MPMOD_MUL_REDCN;
      tmul_best = t_mulredc_n;
    }
#endif
  
  /* Sqr followed by mpn_redc_1 */
#if defined(HAVE___GMPN_REDC_1)
  t_sqrredc_1 = CPUTime();
  for (i = 0; i < iter; ++i) {
    mpn_sqr(tmp, x, N);
    __gmpn_redc_1 (z, tmp, m, N, invm[0]);
    x[0] += tmp[0];
  }
  t_sqrredc_1 = CPUTime()-t_sqrredc_1;
  if (t_sqrredc_1 < tsqr_best)
    {
      tune_sqr[N] = MPMOD_MUL_REDC1;
      tsqr_best = t_sqrredc_1;
    }
#endif
  
  /* Sqr followed by mpn_redc_2 */
#if defined(HAVE___GMPN_REDC_2)
  t_sqrredc_2 = CPUTime();
  for (i = 0; i < iter; ++i) {
    mpn_sqr(tmp, x, N);
    __gmpn_redc_2 (z, tmp, m, N, invm);
    x[0] += tmp[0];
  }
  t_sqrredc_2 = CPUTime()-t_sqrredc_2;
  if (t_sqrredc_2 < tsqr_best)
    {
      tune_sqr[N] = MPMOD_MUL_REDC2;
      tsqr_best = t_sqrredc_2;
    }
#endif
  
  /* Sqr followed by mpn_redc_n */
#if defined(HAVE___GMPN_REDC_N)
  t_sqrredc_n = CPUTime();
  for (i = 0; i < iter; ++i)
    {
      mpn_sqr (tmp, x, N);
      __gmpn_redc_n (z, tmp, m, N, invm);
    }
  t_sqrredc_n = CPUTime()-t_sqrredc_n;
  if (t_sqrredc_n < tsqr_best)
    {
      tune_sqr[N] = MPMOD_MUL_REDCN;
      tsqr_best = t_sqrredc_n;
    }
#endif
  
#if defined(HAVE_NATIVE_MULREDC1_N)
  /* mulredc1 */
  t3 = CPUTime();
  switch (N) {
   case 1:
    for (i=0; i<LOOPCOUNT; ++i) {
      cy += mulredc1(z, x[0], y[0], m[0], invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 2:
    for (i=0; i<LOOPCOUNT; ++i) {
      cy += mulredc1_2(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 3:
    for (i=0; i<LOOPCOUNT; ++i) {
      cy += mulredc1_3(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 4:
    for (i=0; i<LOOPCOUNT; ++i) {
      cy += mulredc1_4(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 5:
    for (i=0; i<LOOPCOUNT; ++i) {
      cy += mulredc1_5(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 6:
    for (i=0; i<LOOPCOUNT; ++i) {
      cy += mulredc1_6(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 7:
    for (i=0; i<LOOPCOUNT; ++i) {
      cy += mulredc1_7(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 8:
    for (i=0; i<LOOPCOUNT; ++i) {
      cy += mulredc1_8(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 9:
    for (i=0; i<LOOPCOUNT; ++i) {
      cy += mulredc1_9(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 10:
    for (i=0; i<LOOPCOUNT; ++i) {
      cy += mulredc1_10(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 11:
    for (i=0; i<LOOPCOUNT; ++i) {
      cy += mulredc1_11(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 12:
    for (i=0; i<LOOPCOUNT; ++i) {
      cy += mulredc1_12(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 13:
    for (i=0; i<LOOPCOUNT; ++i) {
      cy += mulredc1_13(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 14:
    for (i=0; i<LOOPCOUNT; ++i) {
      cy += mulredc1_14(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 15:
    for (i=0; i<LOOPCOUNT; ++i) {
      cy += mulredc1_15(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 16:
    for (i=0; i<LOOPCOUNT; ++i) {
      cy += mulredc1_16(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 17:
    for (i=0; i<LOOPCOUNT; ++i) {
      cy += mulredc1_17(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 18:
    for (i=0; i<LOOPCOUNT; ++i) {
      cy += mulredc1_18(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 19:
    for (i=0; i<LOOPCOUNT; ++i) {
      cy += mulredc1_19(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   case 20:
    for (i=0; i<LOOPCOUNT; ++i) {
      cy += mulredc1_20(z, x[0], y, m, invm[0]);
      x[0] += tmp[0];
    }
    break;
   default: ;
  }
  t3 = CPUTime() - t3;
  t3 *= 1000000.;
#endif /* if defined(HAVE_NATIVE_MULREDC1_N) */

  tmul *= 1000000.;
  tsqr *= 1000000.;
  tredc_1 *= 1000000.;
  tsvoboda1 *= 1000000.;
  t2 *= 1000000.;
  printf("******************\nTime in microseconds per call, size=%lu\n", N);

  /* basic operations */
  printf("mpn_mul_n  = %f\n", tmul/iter);
  printf("mpn_sqr    = %f\n", tsqr/iter);
#ifdef HAVE___GMPN_REDC_1
  printf("mpn_redc_1 = %f\n", tredc_1/iter);
#endif
  if (N > 1)
    printf("svoboda1   = %f\n", tsvoboda1/iter);
#ifdef HAVE___GMPN_REDC_2
  tredc_2 *= 1000000.;
  printf("mpn_redc_2 = %f\n", tredc_2/iter);
#endif
#ifdef HAVE___GMPN_REDC_N
  tredc_n *= 1000000.;
  printf("mpn_redc_n = %f\n", tredc_n/iter);
#endif

  /* modular multiplication and squaring */
  printf("mulredc    = %f\n", t2/iter);
#if defined(HAVE___GMPN_REDC_1)
  t_mulredc_1 *= 1000000.;
  printf("mul+redc_1 = %f\n", t_mulredc_1/iter);
#endif
#if defined(HAVE___GMPN_REDC_2)
  t_mulredc_2 *= 1000000.;
  printf("mul+redc_2 = %f\n", t_mulredc_2/iter);
#endif
#if defined(HAVE___GMPN_REDC_N)
  t_mulredc_n *= 1000000.;
  printf("mul+redc_n = %f\n", t_mulredc_n/iter);
#endif
#if defined(HAVE___GMPN_REDC_1)
  t_sqrredc_1 *= 1000000.;
  printf("sqr+redc_1 = %f\n", t_sqrredc_1/iter);
#endif
#if defined(HAVE___GMPN_REDC_2)
  t_sqrredc_2 *= 1000000.;
  printf("sqr+redc_2 = %f\n", t_sqrredc_2/iter);
#endif
#if defined(HAVE___GMPN_REDC_N)
  t_sqrredc_n *= 1000000.;
  printf("sqr+redc_n = %f\n", t_sqrredc_n/iter);
#endif

  /* multiplication of n limbs by one limb */
  printf ("mulredc1   = %f\n", t3/LOOPCOUNT);
  
  free(tmp);
  free(x); free(y); free(z); free(m);
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

  return 0;
}
