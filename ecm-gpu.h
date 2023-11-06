#ifndef _ECM_GPU_H
#define _ECM_GPU_H 1

#ifndef _DO_NOT_INCLUDE_ECM_IMPL_H
#include "ecm-impl.h"
#endif

#ifdef WITH_GPU

// Absolute limit of CGBN support
#define ECM_GPU_CGBN_MAX_BITS 32*1024

#ifndef ECM_GPU_CURVES_BY_BLOCK
  #define ECM_GPU_CURVES_BY_BLOCK 32
#endif

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

#endif

#endif /* _ECM_GPU_H */
