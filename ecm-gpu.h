#ifndef _ECM_GPU_H
#define _ECM_GPU_H 1

#ifndef _DO_NOT_INCLUDE_ECM_IMPL_H
#include "ecm-impl.h"
#endif

#ifdef WITH_GPU

#ifndef ECM_GPU_NB_DIGITS
  #define ECM_GPU_NB_DIGITS 32 //by default
#endif

#ifndef ECM_GPU_DIGITS
  #define ECM_GPU_DIGITS 0
#endif

#if (ECM_GPU_DIGITS==0)
  #define ECM_GPU_SIZE_DIGIT 32
  typedef unsigned int digit_t;
  typedef int carry_t;
#endif

#define VOL volatile 
//#define VOL

#ifndef ECM_GPU_CURVES_BY_BLOCK
#define ECM_GPU_CURVES_BY_BLOCK 32
#endif

#define ECM_GPU_MAX_BITS ECM_GPU_SIZE_DIGIT * ECM_GPU_NB_DIGITS
typedef digit_t VOL biguint_t[ECM_GPU_NB_DIGITS];
typedef carry_t VOL bigint_t[ECM_GPU_NB_DIGITS];

/* Uncomment the next line to print the number of remaining iterations. */
//#define PRINT_REMAINING_ITER
#endif

#ifndef _DO_NOT_INCLUDE_ECM_IMPL_H

/* cudawrapper.c */
#define gpu_ecm __ECM(gpu_ecm)
#ifdef WITH_GPU
int gpu_ecm (mpz_t, mpz_t, int, mpz_t, mpz_t, mpz_t, double *, double, mpz_t,
             mpz_t, unsigned long, const int, int, int, int, int, int,
             FILE*, FILE*, char*, char *, double, int (*)(void), mpz_t, 
             double *, int, int*, unsigned int*);
#else
int gpu_ecm ();
#endif
#define gpu_ecm_stage1 __ECM(gpu_ecm_stage1)
int gpu_ecm_stage1 (mpz_t *, int *, mpz_t, mpz_t, unsigned int, unsigned int, 
                    float *, int);

#endif

#endif /* _ECM_GPU_H */
