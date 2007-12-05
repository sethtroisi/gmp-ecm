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

#define LOOPCOUNT 1000000000

#include <gmp.h>

extern mp_limb_t mulredc(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y,
             const mp_limb_t *m, mp_size_t n, mp_limb_t inv_m);

extern mp_limb_t mulredc1(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y,
             const mp_limb_t *m, mp_limb_t inv_m);
	  
extern mp_limb_t mulredc3(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y,
             const mp_limb_t *m, mp_limb_t inv_m);

extern mp_limb_t mulredc7(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y,
             const mp_limb_t *m, mp_limb_t inv_m);

extern void ecm_redc3(mp_limb_t * z, const mp_limb_t * x, size_t n, mp_limb_t m);

double CPUTime()
{
  double ret;
#ifndef USE_CLOCK
  struct tms t; 

  times(&t);
  ret = t.tms_utime * 1. / HZ;
#else
  clock_t t;

  t = clock();
  ret = (double)t / (double)CLOCKS_PER_SEC;
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
  mp_limb_t *x, *y, *z, *m, invm, cy, cy2, *tmp;
  int i, j;
  double t1, t2;
  
  x = (mp_limb_t *) malloc(N*sizeof(mp_limb_t));
  y = (mp_limb_t *) malloc(N*sizeof(mp_limb_t));
  z = (mp_limb_t *) malloc((N+1)*sizeof(mp_limb_t));
  m = (mp_limb_t *) malloc(N*sizeof(mp_limb_t));
  tmp = (mp_limb_t *) malloc((2*N+2)*sizeof(mp_limb_t));
 
  mpn_random(m, N);
  m[0] |= 1UL;
  if (m[N-1] == 0) 
    m[N-1] = 1UL;

  invm = 1UL;
  for (i = 0; i < 10; ++i)
    invm = (2*invm-m[0]*invm*invm);
  invm = -invm;

  mpn_random(x, N);
  mpn_random(y, N);

  // Mul followed by ecm_redc3
  t1 = CPUTime();
  for (i=0; i<LOOPCOUNT/N; ++i) {
    mpn_mul_n(tmp, x, y, N);
    ecm_redc3(tmp, m, N, invm);
    cy2 = mpn_add_n (tmp, tmp + N, tmp, N);
    x[0] += tmp[0];
  }
  t1 = CPUTime()-t1;
  
  // Mixed mul and redc
  t2 = CPUTime();
    // Mixed mul and redc
  switch (N) {
   case 1:
    for (i=0; i<LOOPCOUNT/N; ++i) {
      cy += mulredc1(z, x[0], y[0], m[0], invm);
      x[0] += tmp[0];
    }
    break;
   case 2:
    for (i=0; i<LOOPCOUNT/N; ++i) {
      cy += mulredc2(z, x, y, m, invm);
      x[0] += tmp[0];
    }
    break;
   case 3:
    for (i=0; i<LOOPCOUNT/N; ++i) {
      cy += mulredc3(z, x, y, m, invm);
      x[0] += tmp[0];
    }
    break;
   case 4:
    for (i=0; i<LOOPCOUNT/N; ++i) {
      cy += mulredc4(z, x, y, m, invm);
      x[0] += tmp[0];
    }
    break;
   case 5:
    for (i=0; i<LOOPCOUNT/N; ++i) {
      cy += mulredc5(z, x, y, m, invm);
      x[0] += tmp[0];
    }
    break;
   case 6:
    for (i=0; i<LOOPCOUNT/N; ++i) {
      cy += mulredc6(z, x, y, m, invm);
      x[0] += tmp[0];
    }
    break;
   case 7:
    for (i=0; i<LOOPCOUNT/N; ++i) {
      cy += mulredc7(z, x, y, m, invm);
      x[0] += tmp[0];
    }
    break;
   case 8:
    for (i=0; i<LOOPCOUNT/N; ++i) {
      cy += mulredc8(z, x, y, m, invm);
      x[0] += tmp[0];
    }
    break;
   case 9:
    for (i=0; i<LOOPCOUNT/N; ++i) {
      cy += mulredc9(z, x, y, m, invm);
      x[0] += tmp[0];
    }
    break;
   case 10:
    for (i=0; i<LOOPCOUNT/N; ++i) {
      cy += mulredc10(z, x, y, m, invm);
      x[0] += tmp[0];
    }
    break;
   case 11:
    for (i=0; i<LOOPCOUNT/N; ++i) {
      cy += mulredc11(z, x, y, m, invm);
      x[0] += tmp[0];
    }
    break;
   case 12:
    for (i=0; i<LOOPCOUNT/N; ++i) {
      cy += mulredc12(z, x, y, m, invm);
      x[0] += tmp[0];
    }
    break;
   case 13:
    for (i=0; i<LOOPCOUNT/N; ++i) {
      cy += mulredc13(z, x, y, m, invm);
      x[0] += tmp[0];
    }
    break;
   case 14:
    for (i=0; i<LOOPCOUNT/N; ++i) {
      cy += mulredc14(z, x, y, m, invm);
      x[0] += tmp[0];
    }
    break;
   case 15:
    for (i=0; i<LOOPCOUNT/N; ++i) {
      cy += mulredc15(z, x, y, m, invm);
      x[0] += tmp[0];
    }
    break;
   case 16:
    for (i=0; i<LOOPCOUNT/N; ++i) {
      cy += mulredc16(z, x, y, m, invm);
      x[0] += tmp[0];
    }
    break;
   case 17:
    for (i=0; i<LOOPCOUNT/N; ++i) {
      cy += mulredc17(z, x, y, m, invm);
      x[0] += tmp[0];
    }
    break;
   case 18:
    for (i=0; i<LOOPCOUNT/N; ++i) {
      cy += mulredc18(z, x, y, m, invm);
      x[0] += tmp[0];
    }
    break;
   case 19:
    for (i=0; i<LOOPCOUNT/N; ++i) {
      cy += mulredc19(z, x, y, m, invm);
      x[0] += tmp[0];
    }
    break;
   case 20:
    for (i=0; i<LOOPCOUNT/N; ++i) {
      cy += mulredc20(z, x, y, m, invm);
      x[0] += tmp[0];
    }
    break;
   default:
    for (i=0; i<LOOPCOUNT/N; ++i) {
      cy += mulredc20(z, x, y, m,  invm);
      x[0] += tmp[0];
    }
  }
  t2 = CPUTime()-t2;
  
  printf("******************\nN=%d\n",N);
  printf("mul+redc = %f\nmulredc  = %f\n", t1, t2);
  printf("Ratio %f, Normalized time=%f\n", t1/t2, t1/N);

  
  free(tmp);
  free(x); free(y); free(z); free(m);
}
  

int main(int argc, char** argv)
{
  mp_limb_t x[3], y[3], z[3], m[3], invm, carry;
  int nn, i;
  int minsize = 1, maxsize = 20;

  nn = 7;
  
  if (argc > 1)
    minsize = atoi (argv[1]);
  if (argc > 2)
    maxsize = atoi (argv[2]);
  
//  for (;;) {
    for (i = minsize; i <= maxsize; ++i) {
      bench(i);
    }
#if 0
    bench(nn); nn = 2*nn+3;
    bench(nn); nn = 2*nn+3;
    bench(nn); nn = 2*nn+3;
    bench(nn); nn = 2*nn+3;
    bench(nn); nn = 2*nn+3;
    bench(nn); nn = 2*nn+3;
    bench(nn); nn = 2*nn+3;
    bench(nn); nn = 2*nn+3;
    bench(nn); nn = 2*nn+3;
    bench(nn); nn = 2*nn+3;
#endif
//  }
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
