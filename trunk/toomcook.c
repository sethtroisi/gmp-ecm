/* This version works for all sizes, but cannot handle overlapping arrays */
/* Returns the number of multiplies performed */

#include <gmp.h>
#include "ecm.h"

#define A0 A[i]
#define A1 A[l+i]
#define A2 A[2*l+i]
#define B0 B[i]
#define B1 B[l+i]
#define B2 B[2*l+i]
#define C0 C[i]
#define C1 C[l+i]
#define C2 C[2*l+i]
#define C3 C[3*l+i]
#define C4 C[4*l+i]
#define C5 C[4*l+k+i]
#define t0 t[i]
#define t1 t[l+i]
#define t2 t[2*l+i]
#define t3 t[3*l+i]


void mpz_divby3_1op(mpz_t R);

void mpz_divby3_1op(mpz_t R) {
  if (R->_mp_size != 0) {
    mp_size_t abssize = abs(R->_mp_size);
    mpn_divexact_by3(R->_mp_d, R->_mp_d, abssize);
    if (R->_mp_d[abssize-1] == 0)
      R->_mp_size += (R->_mp_size < 0) ? 1 : -1;
  }
}

int toomcook3(mpz_t *C, mpz_t *A, mpz_t *B, int len, mpz_t *t) {
int i, l, k, r;
  
  if (len == 0) 
    return 0;
  
  if (len == 1) {
    mpz_mul(C[0], A[0], B[0]);
    return 1;
  }
  
  if (len == 2) {
    mpz_add(t[0], A[0], A[1]); /* t0 = A_0 + A_1 */
    mpz_add(C[1], B[0], B[1]); /* C1 = B_0 + B_1 */
    mpz_mul(C[1], C[1], t[0]); /* C1 = A_0*B_0 + A_0*B_1 + A_1*B_0 + A_1*B_1 */
    mpz_mul(C[0], A[0], B[0]); /* C0 = A_0 * B_0 */
    mpz_mul(C[2], A[1], B[1]); /* C2 = A_1 * B_1 */
    mpz_sub(C[1], C[1], C[0]); /* C1 = A_0*B_1 + A_1*B_0 + A_1*B_1 */
    mpz_sub(C[1], C[1], C[2]); /* C1 = A_0*B_1 + A_1*B_0 */
    return 3;
  }

#define LEN_4_SHORTCUT
#ifdef LEN_4_SHORTCUT
  /* A 2,2,0 split (12 muls) is less efficient than Karatsuba (9 muls) 
     for len==4 */
  if (len == 4) {
    karatsuba(C, A, B, len, t);
    return 9;
  }
#endif
  
  l = (len+2)/3;
  k = len - 2*l;
  
  for (i=0; i<k; i++) {
    mpz_add(t0, A0, A2);
    mpz_add(t1, B0, B2);
    mpz_sub(t2, t0, A1); /* t2 = A0 - A1 + A2 = A(-1) */
    mpz_sub(t3, t1, B1); /* t3 = B0 - B1 + B2 = B(-1) */
    mpz_add(t0, t0, A1); /* t0 = A0 + A1 + A2 = A(1) */
    mpz_add(t1, t1, B1); /* t1 = B0 + B1 + B2 = B(1) */
  }
  for (; i<l; i++) { /* A3 and B3 are smaller than the rest */
    mpz_add(t0, A0, A1);
    mpz_sub(t2, A0, A1);
    mpz_add(t1, B0, B1);
    mpz_sub(t3, B0, B1);
  }
  
  r = toomcook3(C+2*l, t, t+l, l, t+4*l);
  /* C2 = C(1), len(C1) = 2*l-1 */

  r += toomcook3(t, t+2*l, t+3*l, l, t+4*l);
  /* t0 = C(-1), len(t0) = 2*l-1 */
  
  for (i=0; i<k; i++) {
    mpz_mul_2exp(C0, A2, 1); /* C0 = A(2), C1 = B(2) */
    mpz_add(C0, C0, A1);
    mpz_mul_2exp(C0, C0, 1);
    mpz_add(C0, C0, A0);
    mpz_mul_2exp(C1, B2, 1);
    mpz_add(C1, C1, B1);
    mpz_mul_2exp(C1, C1, 1);
    mpz_add(C1, C1, B0);
  }
  for (; i<l; i++) {
    mpz_mul_2exp(C0, A1, 1);
    mpz_add(C0, C0, A0);
    mpz_mul_2exp(C1, B1, 1);
    mpz_add(C1, C1, B0);
  }  
  
  r += toomcook3(t+2*l, C, C+l, l, t+4*l);
  /* t2 = C(2), len(t2) = 2*l-1 */
  
  r += toomcook3(C, A, B, l, t+4*l);
  /* C0 = C(0), len(C0) = 2*l-1 */
  
  r += toomcook3(C+4*l, A+2*l, B+2*l, k, t+4*l);
  /* C4 = C(inf), len(C4) = 2*k-1 */
  
  /* C0: C_0  C2: C(1)  C4: C_4  t0: C(-1)  t2: C(2) */
  
  for (i=0; i<2*k-1; i++) {
    mpz_sub(t[4*l], C2, t0);
    mpz_add(C2, C2, t0);        /* C2 = C(1)+C(-1) = 2*(C_0 + C_2 + C_4) */
    mpz_set(t0, t[4*l]);        /* t0 = C(1)-C(-1) = 2*(C_1 + C_3) */
    mpz_fdiv_q_2exp(C2, C2, 1); /* C2 = C_0 + C_2 + C_4 */
    mpz_sub(C2, C2, C0);        /* C2 = C_2 + C_4 */
    mpz_sub(C2, C2, C4);        /* C2 = C_2 */
    
    mpz_sub(t2, t2, C0);        /* t2 = 2*C_1 + 4*C_2 + 8*C_3 + 16*C_4 */
    mpz_sub(t2, t2, t0);        /* t2 = 4*C_2 + 6*C_3 + 16*C_4 */
    mpz_fdiv_q_2exp(t2, t2, 1); /* t2 = 2*C_2 + 3*C_3 + 8*C_4 */
    mpz_mul_2exp(t[4*l], C2, 1);
    mpz_sub(t2, t2, t[4*l]);    /* t2 = 3*C_3 + 8*C_4 */
    mpz_mul_2exp(t[4*l], C4, 3);
    mpz_sub(t2, t2, t[4*l]);    /* t2 = 3*C_3 */
    mpz_divby3_1op(t2);         /* t2 = C_3 */
    mpz_fdiv_q_2exp(t0, t0, 1); /* t0 = C_1 + C_3 */
    mpz_sub(t0, t0, t2);        /* t0 = C_1 */
  }
  for (; i<2*l-1; i++) {
    mpz_sub(t[4*l], C2, t0);
    mpz_add(C2, C2, t0);
    mpz_set(t0, t[4*l]);
    mpz_fdiv_q_2exp(C2, C2, 1);
    mpz_sub(C2, C2, C0);
    
    mpz_sub(t2, t2, C0);
    mpz_sub(t2, t2, t0);        /* t2 = 4*C_2 + 6*C_3 + 16*C_4 */
    mpz_fdiv_q_2exp(t2, t2, 1);
    mpz_mul_2exp(t[4*l], C2, 1);
    mpz_sub(t2, t2, t[4*l]);     /* t2 = 3*C_3 + 8*C_4 */
    mpz_divby3_1op(t2);         /* t2 = C_3 */
    mpz_fdiv_q_2exp(t0, t0, 1); /* t0 = C_1 + C_3 */
    mpz_sub(t0, t0, t2);        /* t0 = C_1 */
  }
  
  for (i=0; i<l-1; i++)
    mpz_add(C1, C1, t0);
  mpz_set(C1, t0);
  for (i=l; i<2*l-1; i++)
    mpz_add(C1, C1, t0);

  for (i=0; i<l-1; i++)
    mpz_add(C3, C3, t2);
  mpz_set(C3, t2);
  for (i=l; i<l+k-1; i++)
    mpz_add(C3, C3, t2);
  
  return r;
}
