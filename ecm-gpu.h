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

#ifdef WITH_GPU

/* cudawrapper.c */
#define gpu_ecm __ECM(gpu_ecm)
#define gpu_pm1 __ECM(gpu_pm1)

int gpu_ecm (mpz_t, const ecm_params, ecm_params, mpz_t, double);
int gpu_pm1 (mpz_t, const ecm_params, ecm_params, mpz_t, double);

#endif /* WITH_GPU */

#endif /* _DO_NOT_INCLUDE_ECM_IMPL_H */

#endif /* _ECM_GPU_H */
