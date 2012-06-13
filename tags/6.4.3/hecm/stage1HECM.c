#include "../ecm-impl.h"
#include "hecm.h"
#include "generation.h"
#include "ariKS.h"
#include "morphismes.h"
#include "auxi.h"

/*
  stage 1 with the normal parametrization
  input: n        the number to factored
         f        variable for potential factor
	 k        =lcm(2,..,B1)
	 para     parameters of the curve
	 T1,T2    for the underlying elliptic curve
	 options  options given by the user
  We begin by generating the curve and the point 
   (if needed and authorized by changing the parameters)
  We compute k*P
  We send the points on the elliptic curves
  We check the existence of factors of n
*/
int hecm1Normal (mpz_t f,mpmod_t n,
		 mpz_t k,
		 paraGenCurve para,
		 curve *T1,curve *T2,
		 optionsHECM options) {


  int test;

  curveHyperEll cHEll;
  ksPoint P;
  DivMumfordU divU;
  DivMumfordV divV;
  thetaCst th;
  ksCstPourMul cMul;

  mpres_t g;
  mpz_t s;


  curveHyperEll_init (cHEll,n);
  ksPoint_init (P,n);
  thetaCst_init (th,n);
  ksCstPourMul_init (cMul,n);
  mpres_init (g,n);


  // generate the curve
  if (options->curveSpecified == TRUE) { // we want fixed s or nJacobi 
    if (para->nJacobi == 0) {
      para->nJacobi = (rand()%(N_JACOBI_MAX_p1-N_JACOBI_MIN))+N_JACOBI_MIN;
    }
    
    test = generateNormalCurveSpecified(f,n,para,th,cHEll,P,cMul);
    
  }
  else {
    test = generateNormalCurve(f,n,para,th,cHEll,P,cMul,options);
  }

  if (test == GENERATION_FACTOR_FOUND ) {
    // f contient bien un diviseur de n ie f != 0 mod n
    curveHyperEll_clear (cHEll,n);
    ksPoint_clear (P,n);
    ksCstPourMul_clear (cMul,n);
    thetaCst_clear (th,n);
    mpres_clear (g,n);
    return HECM_FACTOR_FOUND_GENERATION;
  }
  else if (test == GENERATION_FAIL) {
    curveHyperEll_clear (cHEll,n);
    ksPoint_clear (P,n);
    ksCstPourMul_clear (cMul,n);
    thetaCst_clear (th,n);
    mpres_clear (g,n);
    return HECM_GENERATION_FAIL;
  }

  // check if s is in the good interval
  mpz_init (s);
  mpres_get_z (s,para->s,n);
  if ( mpz_cmp (s,options->heightMax) > 0 ) {
    curveHyperEll_clear (cHEll,n);
    ksPoint_clear (P,n);
    ksCstPourMul_clear (cMul,n);
    thetaCst_clear (th,n);
    mpz_clear(s);
    mpres_clear (g,n);
    return HECM_PARAM_TOO_BIG;
  }
  mpz_clear (s);



  // compute k*P
  mulKS(P,cMul,th->be,th->ga,k,n);




  ksCstPourMul_clear (cMul,n);


  DivMumfordU_init (divU,n);
  DivMumfordV_init (divV,n);


  // send the point on the hyperelliptic curve
  test = DivMumford (f,n,divU,divV,P,th,cHEll);
  mpres_set_z (g,f,n);

  ksPoint_clear (P,n);
  thetaCst_clear (th,n);

  if (test == MORPHISM_FAIL ) {
    if (  mpres_is_zero(g,n) == 1 ) { // f=0. Let try again
      mpres_clear (g,n);
      curveHyperEll_clear (cHEll,n);
      DivMumfordU_clear (divU,n);
      DivMumfordV_clear (divV,n);
      return HECM_NO_FACTOR_FOUND;  
    }
    else { // f != 0,n   i.e. f is a real divisor of n
      mpres_clear (g,n);
      curveHyperEll_clear (cHEll,n);
      DivMumfordU_clear (divU,n);
      DivMumfordV_clear (divV,n);
      return HECM_FACTOR_FOUND_MORPHISM;
    } 
  }


  // send the point on the elliptic curves
  //  test = jacCEllJacobi (f,n,T1,T2,divU,divV,cHEll);
  test = jac_to_EllW (f,n,T1,T2,divU,divV,cHEll);

  mpres_set_z (g,f,n);

  curveHyperEll_clear (cHEll,n);
  DivMumfordU_clear (divU,n);
  DivMumfordV_clear (divV,n);





  if (test == MORPHISM_FAIL ) {
    if (  mpres_is_zero(g,n) == 1 ) { // f=0. Let's try again
      mpres_clear (g,n);
      return HECM_NO_FACTOR_FOUND;  
    }
    else { // f != 0,n   i.e. f is a real divisor of n
      mpres_clear (g,n);
      return HECM_FACTOR_FOUND_MORPHISM;
    }
  }
  else if (test == MORPHISM_ERROR) {
    mpres_clear (g,n);
    return HECM_ERROR;
  }
  else if (test == MORPHISM_FOUND_ZERO_CURVE_1) {
    mpres_clear (g,n);
    return HECM_FOUND_ZERO_CURVE_1;
  }
  else if (test == MORPHISM_FOUND_ZERO_CURVE_2) {
    mpres_clear (g,n);
    return HECM_FOUND_ZERO_CURVE_2;
  }
  else if (test == MORPHISM_FOUND_ZERO_CURVE_1_AND_2) {
    mpres_clear (g,n);
    return HECM_FOUND_ZERO_CURVE_1_AND_2;
  }
  else if (test == MORPHISM_STEP1) {
    mpres_clear (g,n);
    return HECM_FACTOR_FOUND_STEP1;
  }
  else {
    mpres_clear (g,n);
    return HECM_NO_FACTOR_FOUND;
  } 

}





/*
  stage 1 with the small parameters
  input: n        the number to factored
         f        variable for potential factor
	 k        =lcm(2,..,B1)
	 para     parameters of the curve
	 T1,T2    for the underlying elliptic curve
	 options  options given by the user
  We begin by generating the curve and the point 
   (if needed and authorized by changing the parameters)
  We compute k*P
  We send the points on the elliptic curves
  We check the existence of factors of n
*/
int hecm1LowParam (mpz_t f,mpmod_t n,
		   mpz_t k,
		   paraGenCurve para ,
		   curve *T1,curve *T2,
		   optionsHECM options) {


  int test;

  curveHyperEll cHEll;
  ksPoint P;
  DivMumfordU divU;
  DivMumfordV divV;
  thetaCst th;
  ksSmallConstPourMul cMul;

  mpres_t g;


  curveHyperEll_init (cHEll,n);
  ksPoint_init (P,n);
  thetaCst_init (th,n);
  mpres_init (g,n);

  
  
  if (options->curveSpecified == TRUE) {
    test = generateCurveSmallParamSpecified (f,n,para,th,cHEll,P,cMul);
  }
  else {
    test = generateCurveSmallParam (f,n,para,th,cHEll,P,cMul,para->a,para->b,para->nJacobi,options);
  }


  if (test == GENERATION_FACTOR_FOUND ) {
    // f contient bien un diviseur de n
    curveHyperEll_clear (cHEll,n);
    ksPoint_clear (P,n);
    thetaCst_clear (th,n);
    mpres_clear (g,n);
    return HECM_FACTOR_FOUND_GENERATION;
  }
  else if (test == GENERATION_FAIL) {
    curveHyperEll_clear (cHEll,n);
    ksPoint_clear (P,n);
    thetaCst_clear (th,n);
    mpres_clear (g,n);
    return HECM_GENERATION_FAIL;
  }
  else if (test == GENERATION_PARAM_TOO_BIG) {
    curveHyperEll_clear (cHEll,n);
    ksPoint_clear (P,n);
    thetaCst_clear (th,n);
    mpres_clear (g,n);
    return HECM_PARAM_TOO_BIG;
  }


  if ( ( mpz_cmp (para->a,options->heightMax) > 0 ) || ( mpz_cmp (para->b,options->heightMax) > 0 ) ) {
    curveHyperEll_clear (cHEll,n);
    ksPoint_clear (P,n);
    thetaCst_clear (th,n);
    mpres_clear (g,n);
    return HECM_PARAM_TOO_BIG;
  }




  // compute k*P
  mulKSsmallParam (P,cMul,th->be,th->ga,k,n);




  DivMumfordU_init (divU,n);
  DivMumfordV_init (divV,n);



  // send the point on the hyperelliptic curve
  test = DivMumford (f,n,divU,divV,P,th,cHEll);
  mpres_set_z (g,f,n);

  ksPoint_clear (P,n);
  thetaCst_clear (th,n);

  if (test == MORPHISM_FAIL ) {
    if (  mpres_is_zero(g,n) == 1 ) { // f=0. Let's try again
      mpres_clear (g,n);
      curveHyperEll_clear (cHEll,n);
      DivMumfordU_clear (divU,n);
      DivMumfordV_clear (divV,n);
      return HECM_NO_FACTOR_FOUND;  
    }
    else { // f != 0,n   i.e. f is a real divisor of n
      mpres_clear (g,n);
      curveHyperEll_clear (cHEll,n);
      DivMumfordU_clear (divU,n);
      DivMumfordV_clear (divV,n);
      return HECM_FACTOR_FOUND_MORPHISM;
    } 
  }



  // send the point on the elliptic curves
  //  test = jacCEllJacobi (f,n,T1,T2,divU,divV,cHEll);
  test = jac_to_EllW (f,n,T1,T2,divU,divV,cHEll);

  mpres_set_z (g,f,n);

  curveHyperEll_clear (cHEll,n);
  DivMumfordU_clear (divU,n);
  DivMumfordV_clear (divV,n);



  if (test == MORPHISM_FAIL ) {
    if (  mpres_is_zero(g,n) == 1 ) { // f=0. Let's try again
      mpres_clear (g,n);
      return HECM_NO_FACTOR_FOUND;  
    }
    else { // f != 0,n  i.e. f is a real divisor of n
      mpres_clear (g,n);
      return HECM_FACTOR_FOUND_MORPHISM;
    }
  }
  else if (test == MORPHISM_ERROR) {
    mpres_clear (g,n);
    return HECM_ERROR;
  }
  else if (test == MORPHISM_FOUND_ZERO_CURVE_1) {
    mpres_clear (g,n);
    return HECM_FOUND_ZERO_CURVE_1;
  }
  else if (test == MORPHISM_FOUND_ZERO_CURVE_2) {
    mpres_clear (g,n);
    return HECM_FOUND_ZERO_CURVE_2;
  }
  else if (test == MORPHISM_FOUND_ZERO_CURVE_1_AND_2) {
    mpres_clear (g,n);
    return HECM_FOUND_ZERO_CURVE_1_AND_2;
  }
  else if (test == MORPHISM_STEP1) {
    mpres_clear (g,n);
    return HECM_FACTOR_FOUND_STEP1;
  }
  else {
    mpres_clear (g,n);
    return HECM_NO_FACTOR_FOUND;
  } 

}
