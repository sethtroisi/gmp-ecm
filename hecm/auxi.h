#ifndef _auxi_H
#define _auxi_H

#include "../ecm-impl.h"

void prodTreeCalculkPlus (mpz_t k,double B1);
void prodTreeCalculkInterPlus (mpz_t Pr,double B1, int n, int* pstop);

void calculk (mpz_t k,double B1);
void prodTreeCalculk (mpz_t k,double B1);
void prodTreeCalculkInter (mpz_t Pr,double B1, int n, int* pstop);





/*
  We have to work in finite algebra over Z/nZ
  The algebra is k[Y]/pol(Y) where degree pol=DEGREE_ALGEBRA
  A special case is when Pol=Y^2+a
*/

#define DEGREE_ALGEBRA 2

                        

typedef mpres_t mpalgres_t[DEGREE_ALGEBRA];

// polynomial for creating the algebra K[x]/P(x)
struct mpalgpol
{
  int kind_of_algebra; // if COMPOSITE_ALGEBRA then degree=4!
  mpres_t coeff[DEGREE_ALGEBRA];
  mpres_t t1;
  mpalgres_t aTemp1;
  mpalgres_t aTemp2;
};
typedef struct mpalgpol mpalgpol_t[1];

void mpalgpol_init (mpalgpol_t, mpmod_t);
void mpalgpol_clear (mpalgpol_t, mpmod_t);


void mpalgres_init (mpalgres_t, mpmod_t);
void mpalgres_clear (mpalgres_t, mpmod_t);

void mpalgres_gcd (mpalgres_t, mpalgres_t, mpmod_t);


void mpalgres_set (mpalgres_t, mpalgres_t ,mpmod_t );
void mpalgres_set_zero (mpalgres_t ,mpmod_t );
void mpalgres_set_ui (mpalgres_t ,unsigned long ,mpmod_t );
void mpalgres_set_mpres (mpalgres_t ,mpres_t ,mpmod_t );

int mpalgres_is_zero (mpalgres_t ,mpalgpol_t ,mpmod_t);
int mpalgres_degree (mpalgres_t ,mpalgpol_t, mpmod_t);

void mpalgres_neg (mpalgres_t, mpalgres_t, mpalgpol_t, mpmod_t);

int mpalgres_invert (mpalgres_t, mpalgres_t, mpalgpol_t ,mpmod_t, mpz_t);

void mpalgres_mul (mpalgres_t, mpalgres_t, mpalgres_t, mpalgpol_t, mpmod_t);
void mpalgres_add (mpalgres_t, mpalgres_t, mpalgres_t, mpalgpol_t, mpmod_t);
void mpalgres_sub (mpalgres_t, mpalgres_t, mpalgres_t, mpalgpol_t, mpmod_t);


void mpalgres_mul_ui (mpalgres_t, mpalgres_t, unsigned long,mpalgpol_t,mpmod_t);
void mpalgres_mul_mpres (mpalgres_t, mpalgres_t, mpres_t, mpalgpol_t, mpmod_t);

void mpalgres_add_mpres (mpalgres_t, mpalgres_t, mpres_t, mpalgpol_t, mpmod_t);
void mpalgres_sub_mpres (mpalgres_t, mpalgres_t, mpres_t, mpalgpol_t, mpmod_t);
void mpalgres_add_ui (mpalgres_t, mpalgres_t, unsigned long,mpalgpol_t,mpmod_t);
void mpalgres_sub_ui (mpalgres_t, mpalgres_t, unsigned long,mpalgpol_t,mpmod_t);








#endif
