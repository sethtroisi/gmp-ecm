#include "../ecm-impl.h"
#include "../ecm.h"

#include "hecm.h"

/* wrapper for GMP-ECM stage2 for a curve in Weierstrass form

   y^2 = x^2 + A * x + B

   where B is implicitly defined by y^2 - (x^2 + A * x) mod n.
*/
int
ecmfactor2 (mpz_t f, mpz_t n, mpz_t A, mpz_t x, mpz_t y, mpz_t B2)
{
  ecm_params q;
  int res;

  ecm_init (q);

  q->sigma_is_A = -1; // indicates that we give a curve in Weierstrass form 
  mpz_set (q->sigma, A);
  mpz_set (q->x, x);
  mpz_set (q->go, y);
  mpz_set (q->B2, B2);

  res = ecm_factor (f, n, 0.0, q);

  ecm_clear (q);

  return res;
}



int hecm2 (mpz_t f, mpmod_t n, curve* T1, curve* T2, mpz_t B2) {
  int test;
  mpres_t g;
  mpz_t x,y,A;
  mpres_init (g,n);
  mpz_init(x);
  mpz_init(y);
  mpz_init(A);



  // stage 2 for the first elliptic curve
  mpres_get_z (A,T1->A,n);
  mpres_get_z (x,T1->x,n);
  mpres_get_z (y,T1->y,n);
  test = ecmfactor2 (f, n->orig_modulus, A, x, y, B2);
  mpres_set_z (g,f,n);
  
  if (test != ECM_NO_FACTOR_FOUND) {
    if (  mpres_is_zero(g,n) == 1 ) { // f=0. 
      mpres_clear (g,n);
      mpz_clear(x);
      mpz_clear(y);
      mpz_clear(A);
      return HECM_FOUND_ZERO_CURVE_1;
    }
    else { // f != 0,n   i.e. f is a real divisor of n
      mpres_clear (g,n);
      mpz_clear(x);
      mpz_clear(y);
      mpz_clear(A);
      return HECM_FACTOR_FOUND_STEP2;
    }
  }

  mpres_get_z (A,T2->A,n);
  mpres_get_z (x,T2->x,n);
  mpres_get_z (y,T2->y,n);
  // stage 2 for the second elliptic curve
  test = ecmfactor2 (f, n->orig_modulus, A, x, y, B2);
  mpres_set_z (g,f,n);
  
  if (test != ECM_NO_FACTOR_FOUND) {
    if (  mpres_is_zero(g,n) == 1 ) { // f=0. 
      mpres_clear (g,n);
      mpz_clear(x);
      mpz_clear(y);
      mpz_clear(A);
      return HECM_FOUND_ZERO_CURVE_2;
    }
    else { // f != 0,n   i.e. f is a real divisor of n
      mpres_clear (g,n);
      mpz_clear(x);
      mpz_clear(y);
      mpz_clear(A);
      return HECM_FACTOR_FOUND_STEP2;
    }
  }


    mpres_clear (g,n);
    mpz_clear(x);
    mpz_clear(y);
    mpz_clear(A);
  return HECM_NO_FACTOR_FOUND;
}
