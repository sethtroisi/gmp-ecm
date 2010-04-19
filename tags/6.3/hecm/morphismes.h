#ifndef _MORPHISMES_H
#define _MORPHISMES_H

#include "../ecm-impl.h"
#include "auxi.h"

#define MORPHISM_ERROR -1 
#define MORPHISM_FAIL 0   // The computation of the morphisms failled
#define MORPHISM_STEP1 1  // We found a factor at the end of stage 1
#define MORPHISM 2        // The computation of the morphisms is a success
#define MORPHISM_FOUND_ZERO_CURVE_1 3 // We found the zero of the 1st curve
#define MORPHISM_FOUND_ZERO_CURVE_2 4 // We found the zero of the 2nd curve
#define MORPHISM_FOUND_ZERO_CURVE_1_AND_2 5 // We found the zero of both curves
#define MORPHISM_FOUND_ZERO_CURVE 6 // We found the zero of one of the curve



struct DivMumfordU_s {
  int degree; // 0,1 ou 2
  mpres_t u0;
  mpres_t u1;
};
typedef struct DivMumfordU_s DivMumfordU[1];

void DivMumfordU_init (DivMumfordU divU, mpmod_t n );
void DivMumfordU_clear (DivMumfordU divU, mpmod_t n );


struct DivMumfordV_s {
  mpres_t V0;
  mpres_t V1;
  mpres_t v1v0;
};
typedef struct DivMumfordV_s DivMumfordV[1];

void DivMumfordV_init (DivMumfordV divV, mpmod_t n );
void DivMumfordV_clear (DivMumfordV divV, mpmod_t n );



#include "generation.h"
#include "ariKS.h"



int DivMumford (mpz_t f,mpmod_t n,DivMumfordU divU,DivMumfordV divV,ksPoint P,thetaCst th,curveHyperEll cHEll);

int DivMumfordDegree2 (mpz_t f,mpmod_t n,DivMumfordU divU,DivMumfordV divV,ksPoint P,thetaCst th,curveHyperEll cHEll, mpres_t T13p, mpres_t T14p, mpres_t T16p);

int degree2_case_y1_equal_zero (mpz_t f,mpmod_t n,curve *T1, curve *T2,DivMumfordU divU, DivMumfordV divV, curveHyperEll cHEll);


int jac_to_EllW (mpz_t f,mpmod_t n,curve *T1, curve *T2,DivMumfordU divU, DivMumfordV divV,curveHyperEll cHEll);

int jac_to_EllW_Degree1 (mpz_t f,mpmod_t n,curve *T1, curve *T2,DivMumfordU divU, curveHyperEll cHEll);

int jac_to_EllW_Degree2 (mpz_t f,mpmod_t n,curve *T1, curve *T2,DivMumfordU divU, DivMumfordV divV,curveHyperEll cHEll);

void HEll_to_EllW(mpmod_t n,mpalgpol_t pol,mpalgres_t aX2,mpalgres_t aY2,mpalgres_t aZ2,mpalgres_t aX,mpalgres_t aY,mpalgres_t aZ,curveHyperEll cHEll );

int coeff_EllW_EllWtwist(mpz_t f,mpmod_t n,mpalgpol_t pol,curve *T,mpres_t a6,mpres_t A2,mpalgres_t R, mpalgres_t a4twist,mpalgres_t B,mpalgres_t invB,curveHyperEll cHEll );

int coeff_EllW(mpz_t f,mpmod_t n,curve *T,mpres_t A2, mpres_t a6,curveHyperEll cHEll );


void EllW_to_EllWshort(mpmod_t n,mpalgpol_t pol,mpalgres_t R,mpres_t A2,mpalgres_t x,mpalgres_t y,mpalgres_t z,mpalgres_t aX,mpalgres_t aY,mpalgres_t aZ);


int HEll_EllW_degree1 (mpz_t f,mpmod_t n,mpres_t x,mpres_t z,curve *T,mpres_t A2, mpres_t a6,curveHyperEll cHEll );

int double_short_weierstrass (mpz_t f,mpmod_t n,curve *T);

int addW_and_scale (mpz_t f,mpmod_t n,mpalgpol_t pol,curve *T,mpres_t a6,mpalgres_t X1 ,mpalgres_t Y1 ,mpalgres_t Z1,mpalgres_t X2 ,mpalgres_t Y2 ,mpalgres_t Z2,mpalgres_t a4twist ,mpalgres_t R);

#endif
