#ifndef _JACOBI_H
#define _JACOBI_H

#include "../ecm-impl.h"
#include "generation.h"
#include "auxi.h"

#include <math.h>


#define MULT_JACOBI_FAIL 0 // The multiplication on the Jacobi curve failled
#define MULT_JACOBI 1     // The multiplication on the Jacobi curve is a success



struct coorJacobi_s {
  mpres_t U;
  mpres_t V;
  mpres_t W;
  mpres_t Y;
};
typedef struct coorJacobi_s coorJacobi[1];

void coorJacobi_init (coorJacobi P, mpmod_t n );
void coorJacobi_clear (coorJacobi P, mpmod_t n );
void coorJacobi_set (coorJacobi P,coorJacobi Q, mpmod_t n );


void mulJacobiEntiers (mpz_t a,mpz_t b,int k,mpz_t x,mpz_t y,mpz_t z);


int mulJacobi2 (mpz_t f,mpmod_t n,int k,mpres_t x,mpres_t y,mpres_t Y,mpres_t ep,mpres_t dep);


void doubleJacobi2DebFin(mpz_t f,mpmod_t n,mpres_t x,mpres_t y,mpres_t Y,mpres_t ep,mpres_t dep);


void doubleJacobi2Deb(mpmod_t n,coorJacobi P2,mpres_t Y,mpres_t ep,mpres_t dep);

void addJacobi2(mpmod_t n,coorJacobi P1,coorJacobi P2,coorJacobi P3,mpres_t ep,mpres_t dep);

void addJacobi2fin(mpmod_t n,coorJacobi P1,coorJacobi P2,mpres_t X3,mpres_t Y3,mpres_t Z3,mpres_t ep,mpres_t dep);

#endif
