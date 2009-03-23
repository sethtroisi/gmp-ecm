#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
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


#include <gmp.h>

double CPUTime()
{
  struct tms t;
  double ret;

  times(&t);
  ret=t.tms_utime*1./HZ;
  return(ret);
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
  for (i=0; i<10000000/N; ++i) {
    mpn_mul_n(tmp, x, y, N);
    ecm_redc3(tmp, m, N, invm);
    cy2 = mpn_add_n (tmp, tmp + N, tmp, N);
    x[0] += tmp[0];
  }
  t1 = CPUTime()-t1;
  
  // Mixed mul and redc
  t2 = CPUTime();
  for (i=0; i<10000000/N; ++i) {
    // Mixed mul and redc
    switch (N) {
     case 1:
      cy = mulredc1(z, x[0], y[0], m[0], invm);
      break;
     case 2:
      cy = mulredc2(z, x, y, m, invm);
      break;
     case 3:
      cy = mulredc3(z, x, y, m, invm);
      break;
     case 4:
      cy = mulredc4(z, x, y, m, invm);
      break;
     case 5:
      cy = mulredc5(z, x, y, m, invm);
      break;
     case 6:
      cy = mulredc6(z, x, y, m, invm);
      break;
     case 7:
      cy = mulredc7(z, x, y, m, invm);
      break;
     case 8:
      cy = mulredc8(z, x, y, m, invm);
      break;
     case 9:
      cy = mulredc9(z, x, y, m, invm);
      break;
     case 10:
      cy = mulredc10(z, x, y, m, invm);
      break;
     case 11:
      cy = mulredc11(z, x, y, m, invm);
      break;
     case 12:
      cy = mulredc12(z, x, y, m, invm);
      break;
     case 13:
      cy = mulredc13(z, x, y, m, invm);
      break;
     case 14:
      cy = mulredc14(z, x, y, m, invm);
      break;
     case 15:
      cy = mulredc15(z, x, y, m, invm);
      break;
     case 16:
      cy = mulredc16(z, x, y, m, invm);
      break;
     case 17:
      cy = mulredc17(z, x, y, m, invm);
      break;
     case 18:
      cy = mulredc18(z, x, y, m, invm);
      break;
     case 19:
      cy = mulredc19(z, x, y, m, invm);
      break;
     case 20:
      cy = mulredc20(z, x, y, m, invm);
      break;
     default:
      cy = mulredc20(z, x, y, m, invm);
    }

    x[0] += tmp[0];
  }
  t2 = CPUTime()-t2;
  
  printf("******************\nN=%d\n",N);
  printf("mul+redc = %f\nmulredc  = %f\n", t1, t2);
  printf("Normalized time=%f\n", t1/N);

  
  free(tmp);
  free(x); free(y); free(z); free(m);
}
  

int main(int argc, char** argv)
{
  mp_limb_t x[3], y[3], z[3], m[3], invm, carry;
  int nn, i;

  nn = 7;
  
//  for (;;) {
    for (i = 1; i <= 20; ++i) {
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
