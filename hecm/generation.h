#ifndef _GENERATION__H
#define _GENERATION__H

#include "../ecm-impl.h"

#include "Jacobi.h"
#include "auxi.h"


#define GENERATION_THESE_SMALL_PARAM_TOO_BIG -1
         // the parameters that were used, give constants that are too big
#define GENERATION_NOT_CURVE 0 
         // There is no curve with these (fixed) parameters
#define GENERATION_A_CURVE 1
         // A curve is generated with these fixed parameters
#define GENERATION_CORRECT_CURVE 2
         // The generated curve can be used
#define GENERATION_FACTOR_FOUND 3
         // We found a factor
#define GENERATION_FAIL 4 
         // The generation have failed (maybe the generated curve can't be used)
#define GENERATION_PARAM_TOO_BIG 5 
         // There is no more curve with small enough parameters


#define POW_MAX 4 // nJacobi < 2^5
#define N_JACOBI_MIN 0
#define N_JACOBI_MAX_p1 (1<<POW_MAX) // max+1


struct paraGenCurve_s {
  mpres_t s;
  mpz_t a,b; // s=a/b mod n;
  mpres_t x;
  mpres_t y;
  int nJacobi;
};
typedef struct paraGenCurve_s paraGenCurve[1];

void paraGenCurve_init (paraGenCurve para, mpmod_t n );
void paraGenCurve_clear (paraGenCurve para, mpmod_t n );


struct ksThetaCst_s {
  mpres_t be;
  mpres_t ga;
  mpres_t t5p;
  mpres_t t6p;
  mpres_t t7p;
  mpres_t t10p;
  mpres_t Rac;
  mpres_t p; // "x0p*y0p" = B^2/(A*D)
};
typedef struct ksThetaCst_s thetaCst[1];

void thetaCst_init (thetaCst th, mpmod_t n );
void thetaCst_clear (thetaCst th, mpmod_t n );


struct curveHyperElliptic_s {
  mpres_t la;
  mpres_t mu;
  mpres_t nu;
  mpres_t q;
};
typedef struct curveHyperElliptic_s curveHyperEll[1];

void curveHyperEll_init (curveHyperEll cHEll, mpmod_t n );
void curveHyperEll_clear (curveHyperEll cHEll, mpmod_t n );



#include "ariKS.h"
#include "hecm.h"



int generateCurve (mpz_t f,mpmod_t n,paraGenCurve para,thetaCst th,curveHyperEll cHEll,ksPoint P,ksCstPourMul cMul,optionsHECM options); 


int generateNormalCurve (mpz_t f,mpmod_t n,paraGenCurve para,thetaCst th,curveHyperEll cHEll,ksPoint P,ksCstPourMul cMul,optionsHECM options); 

int generateNormalCurveSpecified (mpz_t f,mpmod_t n,paraGenCurve para,thetaCst th,curveHyperEll cHEll,ksPoint P,ksCstPourMul cMul);

int generateOneNormalCurve(mpz_t f,mpmod_t n,paraGenCurve para,thetaCst th,curveHyperEll cHEll,ksPoint P,ksCstPourMul cMul); 



#define INVERSE_NOT_COOR 0
#define INVERSE_COOR 1

int inverseCoorPointKS (mpz_t f,mpmod_t n,ksPoint P,ksSmallConstPourMul cMul);



int generateOneCurveSmallParam (mpz_t f,mpmod_t n,paraGenCurve para,thetaCst th,curveHyperEll cHEll,ksPoint P,ksSmallConstPourMul cMul,mpz_t a, mpz_t b, mpz_t x, mpz_t y, mpz_t z);

int generateCurveSmallParam (mpz_t f,mpmod_t n,paraGenCurve para,thetaCst th,curveHyperEll cHEll,ksPoint P,ksSmallConstPourMul cMul,mpz_t a, mpz_t b,int nJacobi,optionsHECM options);

int generateCurveSmallParamSpecified (mpz_t f,mpmod_t n,paraGenCurve para,thetaCst th,curveHyperEll cHEll,ksPoint P,ksSmallConstPourMul cMul);



#define NEXT_PARAM_ERROR 0
#define NEXT_PARAM_CAN_BE_USED 1
#define NEXT_SMALL_PARAM_TOO_BIG 3

#define N_JACOBI_SMALL_PARAM_MAX 20

int nextSmallParam (mpz_t a, mpz_t b,const long Ha,const long Hb);

int nextParam (mpz_t f,mpmod_t n,paraGenCurve para,optionsHECM options);


#endif
