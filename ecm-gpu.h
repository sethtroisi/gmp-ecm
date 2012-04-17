#ifndef _ECM_GPU_H
#define _ECM_GPU_H 1

#ifndef _DO_NOT_INCLUDE_ECM_IMPL_H
#include "ecm-impl.h"
#endif

#ifdef WITH_GPU

#define VERSION_GPUECM "0.1"

#ifndef NB_DIGITS
  #define NB_DIGITS 32 //by default
#endif

#ifndef DIGITS
  #define DIGITS 0
#endif

#if (DIGITS==0)
  #define SIZE_DIGIT 32
  typedef unsigned int digit_t;
  typedef int carry_t;
#endif

#define VOL volatile 
//#define VOL

#define ECM_GPU_MAX_BITS SIZE_DIGIT*NB_DIGITS
typedef digit_t VOL biguint_t[NB_DIGITS];
typedef digit_t VOL dbiguint_t[NB_DIGITS+1];
typedef carry_t VOL dbigint_t[NB_DIGITS+1];

/* Uncomment the next line to print the number of remaining iterations. */
//#define PRINT_REMAINING_ITER
#endif

#ifndef _DO_NOT_INCLUDE_ECM_IMPL_H

/* cudawrapper.c */
#define gpu_ecm __ECM(gpu_ecm)
int gpu_ecm (mpz_t, int, mpz_t, mpz_t, mpz_t, double *, double, mpz_t, 
             mpz_t, double, unsigned long, const int, int, int, int, int, int,
             FILE*, FILE*, char*, char *, double, double, gmp_randstate_t, 
             int (*)(void), mpz_t, double *, int, int*, unsigned int*);
#define gpu_ecm_stage1 __ECM(gpu_ecm_stage1)
int gpu_ecm_stage1 (mpz_t, mpz_t, unsigned int, unsigned int, float*);

#endif

#endif /* _ECM_GPU_H */
