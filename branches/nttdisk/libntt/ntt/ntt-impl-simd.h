#ifndef _NTT_IMPL_SIMD_H
#define _NTT_IMPL_SIMD_H

#include "ntt-impl.h"

/* SIMD includes */
#if defined(HAVE_SSE2)
#include "ntt-impl-sse2.h"
#elif defined(HAVE_SSE41)
#include "ntt-impl-sse41.h"
#endif

extern const nttconfig_t V(ntt2simd_config);
extern const nttconfig_t V(ntt3simd_config);
extern const nttconfig_t V(ntt4simd_config);
extern const nttconfig_t V(ntt5simd_config);
extern const nttconfig_t V(ntt7simd_config);
extern const nttconfig_t V(ntt8simd_config);
extern const nttconfig_t V(ntt9simd_config);
extern const nttconfig_t V(ntt15simd_config);
extern const nttconfig_t V(ntt16simd_config);
extern const nttconfig_t V(ntt35simd_config);
extern const nttconfig_t V(ntt40simd_config);

#endif /* _NTT_IMPL_SIMD_H */
