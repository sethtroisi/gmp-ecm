#ifndef _NTT_IMPL_SIMD_H
#define _NTT_IMPL_SIMD_H

#include "ntt-impl.h"

/* SIMD includes */
#if defined(HAVE_SSE2)
#include "ntt-impl-sse2.h"
#elif defined(HAVE_SSE42)
#include "ntt-impl-sse42.h"
#elif defined(HAVE_AVX)
#include "ntt-impl-avx.h"
#elif defined(HAVE_AVX2)
#include "ntt-impl-avx2.h"
#elif defined(HAVE_FMA)
#include "ntt-impl-fma.h"
#endif

#ifdef __cplusplus
extern "C" {
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

/* code generator */

#define DECLARE_CORE_ROUTINES_SIMD(N)				\
								\
static void							\
ntt##N##_run_simd(spv_t in, spv_size_t stride, 			\
    		spv_size_t dist, spv_size_t num_transforms, 	\
		sp_t p, spv_t ntt_const)			\
{								\
  spv_size_t i;							\
  spv_size_t num_simd = SP_SIMD_VSIZE * 			\
    			(num_transforms / SP_SIMD_VSIZE);	\
								\
  for (i = 0; i < num_simd; i += SP_SIMD_VSIZE)			\
    ntt##N##_run_core_simd(in + i * dist, stride, dist, 	\
			in + i * dist, stride, dist, p, 	\
			ntt_const, SP_SIMD_VSIZE);		\
								\
  if (i < num_transforms)					\
    ntt##N##_run_core_simd(in + i * dist, stride, dist, 	\
			in + i * dist, stride, dist, p, 	\
			ntt_const, num_transforms - i);		\
}								\
								\
								\
static void							\
ntt##N##_run_simd_interleaved(					\
    	spv_t in, spv_size_t stride, spv_size_t dist,		\
    	spv_size_t num_transforms, spv_t p, 			\
	spv_size_t vsize, spv_t ntt_const)			\
{								\
  spv_size_t i, j;						\
  spv_size_t alloc = 2 * NC;					\
								\
  for (i = 0; i < num_transforms; i++)				\
    for (j = 0; j < vsize; j += SP_SIMD_VSIZE)			\
      ntt##N##_run_core_simd_interleaved(			\
	  	in + i * dist + j, stride, 			\
                in + i * dist + j, stride, 			\
		sp_simd_load(p + j), 				\
	  	(sp_simd_t *)(ntt_const + j * alloc));		\
}								\
								\
static void							\
ntt##N##_twiddle_run_simd(spv_t in, spv_size_t stride, 		\
    			spv_size_t dist, spv_size_t num_transforms, \
			sp_t p, spv_t ntt_const, spv_t w)	\
{								\
  spv_size_t i, j;						\
  spv_size_t num_simd = SP_SIMD_VSIZE * 			\
    		(num_transforms / SP_SIMD_VSIZE);		\
								\
  for (i = j = 0; i < num_simd; i += SP_SIMD_VSIZE,		\
		  	j += 2*((N)-1)*SP_SIMD_VSIZE)		\
    ntt##N##_twiddle_run_core_simd(				\
		in + i * dist, stride, dist,			\
		in + i * dist, stride, dist,			\
		(sp_simd_t *)(w + j), p, ntt_const,		\
		SP_SIMD_VSIZE);					\
								\
  if (i < num_transforms)					\
    ntt##N##_twiddle_run_core_simd(				\
		in + i * dist, stride, dist,			\
		in + i * dist, stride, dist,			\
		(sp_simd_t *)(w + j), p, ntt_const,		\
		num_transforms - i);				\
}								\
								\
static void							\
ntt##N##_twiddle_run_simd_interleaved(				\
    	spv_t in, spv_size_t stride, spv_size_t dist,		\
    	spv_size_t num_transforms, spv_t p, 			\
	spv_size_t vsize, spv_t ntt_const, spv_t w)		\
{								\
  spv_size_t i, j;						\
  spv_size_t alloc = 2 * NC;					\
  spv_size_t twiddle_alloc = 2 * ((N)-1);			\
  spv_size_t max_twiddle_alloc = num_transforms * twiddle_alloc; \
								\
  for (j = 0; j < vsize; j += SP_SIMD_VSIZE)			\
    for (i = 0; i < num_transforms; i++)			\
      {								\
	sp_simd_t * w0 = (sp_simd_t *)(w + 			\
	    		j * max_twiddle_alloc +			\
	    		SP_SIMD_VSIZE * i * twiddle_alloc);	\
								\
	ntt##N##_twiddle_run_core_simd_interleaved(		\
	  		in + i * dist + j, stride,		\
			in + i * dist + j, stride,		\
			w0,					\
	  		sp_simd_load(p + j), 			\
			(sp_simd_t *)(ntt_const + j * alloc));	\
      }								\
}								\
								\
								\
static void							\
ntt##N##_pfa_run_simd(						\
    	spv_t in, spv_size_t stride, spv_size_t dist,		\
    	spv_size_t num_transforms, sp_t p, 			\
	spv_size_t cofactor, spv_t ntt_const)			\
{								\
  spv_size_t i, j;						\
  spv_size_t inc = cofactor * stride;				\
  spv_size_t inc2 = (N) * stride;				\
  spv_size_t ntt_size = (N) * cofactor * stride;		\
  spv_size_t num_simd = SP_SIMD_VSIZE * 			\
    		(cofactor / SP_SIMD_VSIZE);			\
								\
  for (i = 0; i < num_transforms; i++)				\
    {								\
      spv_size_t incstart = 0;					\
								\
      for (j = 0; j < num_simd; j += SP_SIMD_VSIZE, 		\
  			incstart += SP_SIMD_VSIZE * inc2)	\
	ntt##N##_pfa_run_core_simd(in + i * dist, 		\
	    			incstart, inc, inc2, ntt_size,	\
				p, ntt_const, SP_SIMD_VSIZE);	\
								\
      if (j < cofactor)						\
	ntt##N##_pfa_run_core_simd(in + i * dist, 		\
	    			incstart, inc, inc2, ntt_size, 	\
				p, ntt_const, cofactor - j);	\
    }								\
}								\
								\
static void							\
ntt##N##_pfa_run_simd_interleaved(				\
	spv_t in, spv_size_t stride, spv_size_t dist,		\
    	spv_size_t num_transforms, spv_t p, 			\
	spv_size_t cofactor, spv_size_t vsize,			\
	spv_t ntt_const)					\
{								\
  spv_size_t i, j, k;						\
  spv_size_t alloc = 2 * NC;					\
  spv_size_t inc = cofactor * stride;				\
  spv_size_t inc2 = (N) * stride;				\
  spv_size_t ntt_size = (N) * cofactor * stride;		\
								\
  for (k = 0; k < vsize; k += SP_SIMD_VSIZE)			\
    {								\
      spv_size_t incstart = 0;					\
								\
      for (j = 0; j < cofactor; j++, incstart += inc2)		\
	for (i = 0; i < num_transforms; i++)			\
	  ntt##N##_pfa_run_core_simd_interleaved(		\
	      		in + i * dist + k, incstart, 		\
	      		inc, ntt_size, sp_simd_load(p + k), 	\
			(sp_simd_t *)(ntt_const + k * alloc));	\
    }								\
}								\
								\
const nttconfig_t V(ntt##N##simd_config) = 			\
{								\
  N,								\
  NC,								\
  ntt##N##_fixed_const,						\
  X(ntt##N##_init),						\
								\
  ntt##N##_run_simd,						\
  ntt##N##_pfa_run_simd,					\
  ntt##N##_twiddle_run_simd,					\
								\
  ntt##N##_run_simd_interleaved,				\
  ntt##N##_pfa_run_simd_interleaved,				\
  ntt##N##_twiddle_run_simd_interleaved				\
};


#ifdef __cplusplus
}
#endif

#endif /* _NTT_IMPL_SIMD_H */
