#ifndef _HECM_H
#define _HECM_H

#include "../ecm-impl.h"

#define HECM_FOUND_N -2
#define HECM_ERROR -1 
#define HECM_NO_FACTOR_FOUND 0 
#define HECM_FACTOR_FOUND_STEP1 1
#define HECM_FACTOR_FOUND_STEP2 2
#define HECM_FACTOR_FOUND_GENERATION 3
#define HECM_FACTOR_FOUND_MORPHISM 4
#define HECM_FOUND_ZERO_CURVE_1 5
#define HECM_FOUND_ZERO_CURVE_2 6
#define HECM_FOUND_ZERO_CURVE_1_AND_2 7
#define HECM_GENERATION_FAIL 8
#define HECM_PARAM_TOO_BIG 9

#define TRUE 1
#define FALSE 0


struct optionsHECM_s {
  int smallParam;             // use small parameters (default = TRUE)
  int curveSpecified;         // use a specified curve
  int initialCurveSpecified;  // the first curve is specified
  int verbose;                // set the verbose mode
  unsigned int nbtests;       // number of tests to do
  mpz_t heightMin;            // minimal height for the parameters
  mpz_t heightMax;            // maximal height for the parameters
};
typedef struct optionsHECM_s optionsHECM[1];

void optionsHECM_init (optionsHECM options);
void optionsHECM_clear (optionsHECM options);

#include "generation.h"

/* stage 1 */
int hecm1Normal (mpz_t f ,mpmod_t n ,mpz_t k ,paraGenCurve para ,curve *T1,curve *T2,optionsHECM options);
int hecm1LowParam (mpz_t f ,mpmod_t n ,mpz_t k ,paraGenCurve para ,curve *T1,curve *T2,optionsHECM options);

/* stage 2 */
int ecmfactor2 (mpz_t f, mpz_t n, mpz_t A, mpz_t x, mpz_t y, mpz_t B2);
int hecm2 (mpz_t f, mpmod_t n, curve* T1, curve* T2, mpz_t B2);



#endif
