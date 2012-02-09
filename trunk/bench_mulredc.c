#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
# include <sys/time.h>
# include <sys/times.h>
# include <sys/param.h>
# include <sys/resource.h>
# include <string.h>
# include <unistd.h>
# ifndef HZ
#  ifdef CLK_TCK
#   define HZ CLK_TCK
#  else
#   define HZ 100
#  endif
# endif

#define LOOPCOUNT 10000000UL

#include <gmp.h>
#include "mulredc.h"

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

void bench(mp_size_t N)
{
  mp_limb_t *x, *y, *z, *m, invm[2], cy, *tmp;
  unsigned long i;
  const unsigned long iter = LOOPCOUNT/N;
  double t2, t3 = 0., tmul, tsqr, tredc_1, tredc_2;
  mpz_t M, B;
#if defined(HAVE_ASM_REDC3)
  double t1;
  mp_limb_t cy2;
#endif
  
  x = (mp_limb_t *) malloc(N*sizeof(mp_limb_t));
  y = (mp_limb_t *) malloc(N*sizeof(mp_limb_t));
  z = (mp_limb_t *) malloc((N+1)*sizeof(mp_limb_t));
  m = (mp_limb_t *) malloc(N*sizeof(mp_limb_t));
  tmp = (mp_limb_t *) malloc((2*N+2)*sizeof(mp_limb_t));
 
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

  invm[0] = mpz_getlimbn (M, 0);
  invm[1] = mpz_getlimbn (M, 1);

  mpz_clear (M);
  mpz_clear (B);

  mpn_random(x, N);
  mpn_random(y, N);

  tmul = CPUTime();
  for (i = 0; i < iter; ++i)
    mpn_mul_n(tmp, x, y, N);
  tmul = CPUTime()-tmul;

#ifdef HAVE_MPN_SQR
  tsqr = CPUTime();
  for (i = 0; i < iter; ++i)
    mpn_sqr(tmp, x, N);
  tsqr = CPUTime()-tsqr;
#endif

#ifdef HAVE___GMPN_REDC_1
  tredc_1 = CPUTime();
  mpn_mul_n(tmp, x, y, N);
  for (i = 0; i < iter; ++i)
    __gmpn_redc_1 (z, tmp, m, N, invm[0]);
  tredc_1 = CPUTime()-tredc_1;
#endif

#ifdef HAVE___GMPN_REDC_2
  tredc_2 = CPUTime();
  mpn_mul_n(tmp, x, y, N);
  for (i = 0; i < iter; ++i)
    __gmpn_redc_2 (z, tmp, m, N, invm);
  tredc_2 = CPUTime()-tredc_2;
#endif

  /* Mul followed by ecm_redc3 */
#ifdef HAVE_ASM_REDC3
  t1 = CPUTime();
  for (i = 0; i < iter; ++i) {
    mpn_mul_n(tmp, x, y, N);
    ecm_redc3(tmp, m, N, invm[0]);
    x[0] += tmp[0];
  }
  t1 = CPUTime()-t1;
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
  tredc_2 *= 1000000.;
  t2 *= 1000000.;
  printf("******************\nTime in microseconds per call, size=%lu\n", N);

  printf("mpn_mul_n  = %f\n", tmul/iter);
#ifdef HAVE_MPN_SQR
  printf("mpn_sqr    = %f\n", tsqr/iter);
#endif
#ifdef HAVE___GMPN_REDC_1
  printf("mpn_redc_1 = %f\n", tredc_1/iter);
#endif
#ifdef HAVE___GMPN_REDC_2
  printf("mpn_redc_2 = %f\n", tredc_2/iter);
#endif
#if defined(HAVE_ASM_REDC3)
  t1 *= 1000000.;
  printf("mul+redc   = %f\n", t1/iter);
#endif
  printf("mulredc    = %f\n", t2/iter);
  printf("mulredc1   = %f\n", t3/LOOPCOUNT);
  
  free(tmp);
  free(x); free(y); free(z); free(m);
}
  

int main(int argc, char** argv)
{
  int i;
  int minsize = 1, maxsize = 20;

  if (argc > 1)
    minsize = atoi (argv[1]);
  if (argc > 2)
    maxsize = atoi (argv[2]);
  
    for (i = minsize; i <= maxsize; ++i) {
      bench(i);
    }

  return 0;
}


#if 0

W := 2^64;

x0:= 12580274668139321508;
x1:= 9205793975152560417;
x2:= 7857372727033793057;
x := x0 + W*(x1 + W*x2);

y0:= 13688385828267279103;
y1:= 10575011835742767258;
y2:= 8802048318027595690;
y := y0 + W*(y1 + W*y2);
  
m0:= 2981542467342508025;
m1:= 5964669706257742025;
m2:= 18446744073678090270;
m := m0 + W*(m1 + W*m2);
  
invm := 9419286575570128311;



#endif
