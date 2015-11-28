#ifndef _NTT_IMPL_SSE2_H
#define _NTT_IMPL_SSE2_H

#include "ntt-impl-scalar.h"
#include <immintrin.h>

#ifdef __cplusplus
extern "C" {
#endif

#define HAVE_SIMD
#define SP_SIMD_VSIZE (128 / SP_TYPE_BITS)

/* mangled routine names */
#define V(name) MANGLE_SSE2(X(name))
#define SP_SIMD_NAME_SUFFIX_STR SP_NAME_SUFFIX_STR "sse2"

#define pload(addr)  _mm_load_si128((__m128i const *)(addr))
#define ploadu(addr)  _mm_loadu_si128((__m128i const *)(addr))
#define pload_lo32(addr)  (__m128i)_mm_load_ss((float *)(addr))
#define pload_lo64(addr)  (__m128i)_mm_load_sd((double const *)(addr))
#define pload_hi64(x, addr)  (__m128i)_mm_loadh_pd((__m128d)x, (double const *)(addr))
#define pstore(x, addr) _mm_store_si128((__m128i *)(addr), x)
#define pstoreu(x, addr) _mm_storeu_si128((__m128i *)(addr), x)
#define pstore_lo32(x, addr)  _mm_store_ss((float *)(addr), (__m128)x)
#define pstore_lo64(x, addr)  _mm_store_sd((double *)(addr), (__m128d)x)
#define pstore_hi64(x, addr)  _mm_storeh_pd((double *)(addr), (__m128d)x)
#define pand _mm_and_si128
#define pxor _mm_xor_si128
#define psetzero() _mm_setzero_si128()
#define paddd _mm_add_epi32
#define paddq _mm_add_epi64
#define psubd _mm_sub_epi32
#define psubq _mm_sub_epi64
#define pmuludq _mm_mul_epu32
#define pslld _mm_slli_epi32
#define psllq _mm_slli_epi64
#define psrld _mm_srli_epi32
#define psrlq _mm_srli_epi64
#define pshufd _mm_shuffle_epi32
#define pcmpgtd _mm_cmpgt_epi32
#define punpcklo32 _mm_unpacklo_epi32
#define punpcklo64 _mm_unpacklo_epi64
#define pcvt_i32 _mm_cvtsi32_si128

#if defined(_WIN64) || defined(__x86_64__)
#define pcvt_i64 _mm_cvtsi64_si128
#define pstore_i64(out, in) out = _mm_cvtsi128_si64(in)
#else
#define pcvt_i64(x) _mm_loadl_epi64((__m128i const *)&(x))
#define pstore_i64(out, in) _mm_storel_epi64((__m128i *)&(out), in)
#endif


typedef __m128i sp_simd_t;


#define sp_simd_load(x) pload(x)

static INLINE sp_simd_t sp_simd_gather(spv_t x, spv_size_t dist,
    					spv_size_t vsize)
{
#if SP_TYPE_BITS == 32

  switch (vsize)
    {
      case 4:
	if (dist == 1)
	  return ploadu(x);
	else
	  {
	    sp_simd_t t0 = pload_lo32(x + 0);
	    sp_simd_t t1 = pload_lo32(x + dist);
	    sp_simd_t t2 = pload_lo32(x + 2 * dist);
	    sp_simd_t t3 = pload_lo32(x + 3 * dist);
	    sp_simd_t r0 = punpcklo32(t0, t1);
	    sp_simd_t r1 = punpcklo32(t2, t3);
	    return punpcklo64(r0, r1);
	  }

      case 3:
	{
	  sp_simd_t t0 = pload_lo32(x + 0);
	  sp_simd_t t1 = pload_lo32(x + dist);
	  sp_simd_t t2 = pload_lo32(x + 2 * dist);
	  sp_simd_t r0 = punpcklo32(t0, t1);
	  return punpcklo64(r0, t2);
	}

      case 2:
	{
	  sp_simd_t t0 = pload_lo32(x + 0);
	  sp_simd_t t1 = pload_lo32(x + dist);
	  return punpcklo32(t0, t1);
	}

      default:
	return pload_lo32(x + 0);
    }

#else

  switch (vsize)
    {
      case 2:
	{
	  sp_simd_t t = pload_lo64(x + 0);
	  return pload_hi64(t, x + dist);
	}

      default:
	return pload_lo64(x + 0);
    }

#endif
}

static INLINE sp_simd_t sp_simd_pfa_gather(spv_t x, spv_size_t start_off, 
					spv_size_t inc, spv_size_t n,
					spv_size_t vsize)
{
#if SP_TYPE_BITS == 32
  switch (vsize)
    {
      case 4:
	{
	  spv_size_t j0 = start_off;
	  spv_size_t j1 = sp_array_inc(j0, inc, n);
	  spv_size_t j2 = sp_array_inc(j0, 2 * inc, n);
	  spv_size_t j3 = sp_array_inc(j0, 3 * inc, n);
	  sp_simd_t t0 = pload_lo32(x + j0);
	  sp_simd_t t1 = pload_lo32(x + j1);
	  sp_simd_t t2 = pload_lo32(x + j2);
	  sp_simd_t t3 = pload_lo32(x + j3);
	  sp_simd_t r0 = punpcklo32(t0, t1);
	  sp_simd_t r1 = punpcklo32(t2, t3);
	  return punpcklo64(r0, r1);
	}

      case 3:
	{
	  spv_size_t j0 = start_off;
	  spv_size_t j1 = sp_array_inc(j0, inc, n);
	  spv_size_t j2 = sp_array_inc(j0, 2 * inc, n);
	  sp_simd_t t0 = pload_lo32(x + j0);
	  sp_simd_t t1 = pload_lo32(x + j1);
	  sp_simd_t t2 = pload_lo32(x + j2);
	  sp_simd_t r0 = punpcklo32(t0, t1);
	  return punpcklo64(r0, t2);
	}

      case 2:
	{
	  spv_size_t j0 = start_off;
	  spv_size_t j1 = sp_array_inc(j0, inc, n);
	  sp_simd_t t0 = pload_lo32(x + j0);
	  sp_simd_t t1 = pload_lo32(x + j1);
	  return punpcklo32(t0, t1);
	}

      default:
	{
	  spv_size_t j0 = start_off;
	  return pload_lo32(x + j0);
	}
    }

#else

  switch (vsize)
    {
      case 2:
	{
	  spv_size_t j0 = start_off;
	  spv_size_t j1 = sp_array_inc(j0, inc, n);
	  sp_simd_t t = pload_lo64(x + j0);
	  return pload_hi64(t, x + j1);
	}

      default:
	{
	  spv_size_t j0 = start_off;
	  return pload_lo64(x + j0);
	}
    }

#endif
}

#define sp_simd_store(t, x) pstore(t, x)

static INLINE void sp_simd_scatter(sp_simd_t t, spv_t x, spv_size_t dist,
    					spv_size_t vsize)
{
#if SP_TYPE_BITS == 32

  switch (vsize)
    {
      case 4:
	if (dist == 1)
	  pstoreu(t, x);
	else
	  {
	    sp_simd_t t1 = _mm_srli_si128(t, 4);
	    sp_simd_t t2 = _mm_srli_si128(t, 8);
	    sp_simd_t t3 = _mm_srli_si128(t, 12);
	    pstore_lo32(t, x + 0);
	    pstore_lo32(t1, x + dist);
	    pstore_lo32(t2, x + 2 * dist);
	    pstore_lo32(t3, x + 3 * dist);
	  }
	break;

      case 3:
	{
	  sp_simd_t t1 = _mm_srli_si128(t, 4);
	  sp_simd_t t2 = _mm_srli_si128(t, 8);
	  pstore_lo32(t, x + 0);
	  pstore_lo32(t1, x + dist);
	  pstore_lo32(t2, x + 2 * dist);
      	}
	break;

      case 2:
	{
	  sp_simd_t t1 = _mm_srli_si128(t, 4);
	  pstore_lo32(t, x + 0);
	  pstore_lo32(t1, x + dist);
      	}
	break;

      default:
	pstore_lo32(t, x + 0);
	break;
    }

#else

  switch (vsize)
    {
      case 2:
	pstore_lo64(t, x + 0);
	pstore_hi64(t, x + dist);
	break;

      default:
	pstore_lo64(t, x + 0);
	break;
    }

#endif
}

static INLINE void sp_simd_pfa_scatter(sp_simd_t t, spv_t x, 
    				spv_size_t start_off, 
				spv_size_t inc, spv_size_t n,
				spv_size_t vsize)
{
#if SP_TYPE_BITS == 32

  switch (vsize)
    {
      case 4:
	{
	  spv_size_t j0 = start_off;
	  spv_size_t j1 = sp_array_inc(j0, inc, n);
	  spv_size_t j2 = sp_array_inc(j0, 2 * inc, n);
	  spv_size_t j3 = sp_array_inc(j0, 3 * inc, n);
	  sp_simd_t t1 = _mm_srli_si128(t, 4);
	  sp_simd_t t2 = _mm_srli_si128(t, 8);
      	  sp_simd_t t3 = _mm_srli_si128(t, 12);
	  pstore_lo32(t, x + j0);
	  pstore_lo32(t1, x + j1);
	  pstore_lo32(t2, x + j2);
	  pstore_lo32(t3, x + j3);
	}
	break;

      case 3:
	{
	  spv_size_t j0 = start_off;
	  spv_size_t j1 = sp_array_inc(j0, inc, n);
	  spv_size_t j2 = sp_array_inc(j0, 2 * inc, n);
	  sp_simd_t t1 = _mm_srli_si128(t, 4);
	  sp_simd_t t2 = _mm_srli_si128(t, 8);
	  pstore_lo32(t, x + j0);
	  pstore_lo32(t1, x + j1);
	  pstore_lo32(t2, x + j2);
	}
	break;

      case 2:
	{
	  spv_size_t j0 = start_off;
	  spv_size_t j1 = sp_array_inc(j0, inc, n);
	  sp_simd_t t1 = _mm_srli_si128(t, 4);
	  pstore_lo32(t, x + j0);
	  pstore_lo32(t1, x + j1);
	}
	break;

      default:
	{
	  spv_size_t j0 = start_off;
	  pstore_lo32(t, x + j0);
	}
	break;
    }

#else

  switch (vsize)
    {
      case 2:
	{
	  spv_size_t j0 = start_off;
	  spv_size_t j1 = sp_array_inc(j0, inc, n);
	  pstore_lo64(t, x + j0);
	  pstore_hi64(t, x + j1);
	}

      default:
	{
	  spv_size_t j0 = start_off;
	  pstore_lo64(t, x + j0);
	}
    }

#endif
}

static INLINE sp_simd_t sp_ntt_add_simd_core(
    			sp_simd_t a, sp_simd_t b, sp_simd_t vp)
{
  sp_simd_t t0, t1;

#if SP_TYPE_BITS == 32
  t0 = paddd(a, b);
  t0 = psubd(t0, vp);
  t1 = pcmpgtd(psetzero(), t0);
  t1 = pand(t1, vp);
  return paddd(t0, t1);

#else
  t0 = paddq(a, b);
  t0 = psubq(t0, vp);
  t1 = pcmpgtd(psetzero(), t0);
  t1 = pshufd(t1, 0xf5);
  t1 = pand(t1, vp);
  return paddq(t0, t1);

#endif
}

static INLINE sp_simd_t sp_ntt_add_partial_simd_core(
    			sp_simd_t a, sp_simd_t b, sp_simd_t vp)
{
#ifdef HAVE_PARTIAL_MOD

  #if SP_TYPE_BITS == 32
  return paddd(a, b);
  #else
  return paddq(a, b);
  #endif

#else

  return sp_ntt_add_simd_core(a, b, vp);

#endif
}

static INLINE sp_simd_t sp_ntt_add_simd(sp_simd_t a, sp_simd_t b, sp_t p)
{
  #ifdef HAVE_PARTIAL_MOD
  p = 2 * p;
  #endif

  return sp_ntt_add_simd_core(a, b, 
#if SP_TYPE_BITS == 32
			pshufd(pcvt_i32(p), 0x00));
#else
			pshufd(pcvt_i64(p), 0x44));
#endif
}

static INLINE sp_simd_t sp_ntt_add_simd0(sp_simd_t a, sp_simd_t b, sp_simd_t p)
{
  #ifdef HAVE_PARTIAL_MOD
    #if SP_TYPE_BITS == 32
    p = paddd(p, p);
    #else
    p = paddq(p, p);
    #endif
  #endif

  return sp_ntt_add_simd_core(a, b, p);
}

static INLINE sp_simd_t sp_ntt_add_partial_simd(sp_simd_t a, sp_simd_t b, sp_t p)
{
  #ifdef HAVE_PARTIAL_MOD
  p = 2 * p;
  #endif

  return sp_ntt_add_partial_simd_core(a, b, 
#if SP_TYPE_BITS == 32
			pshufd(pcvt_i32(p), 0x00));
#else
			pshufd(pcvt_i64(p), 0x44));
#endif
}

static INLINE sp_simd_t sp_ntt_add_partial_simd0(sp_simd_t a, sp_simd_t b, sp_simd_t p)
{
  #ifdef HAVE_PARTIAL_MOD
    #if SP_TYPE_BITS == 32
    p = paddd(p, p);
    #else
    p = paddq(p, p);
    #endif
  #endif

  return sp_ntt_add_partial_simd_core(a, b, p);
}

static INLINE sp_simd_t sp_ntt_sub_simd_core(
    			sp_simd_t a, sp_simd_t b, sp_simd_t vp)
{
  sp_simd_t t0, t1;

#if SP_TYPE_BITS == 32
  t0 = psubd(a, b);
  t1 = pcmpgtd(psetzero(), t0);
  t1 = pand(t1, vp);
  return paddd(t0, t1);

#else
  t0 = psubq(a, b);
  t1 = pcmpgtd(psetzero(), t0);
  t1 = pshufd(t1, 0xf5);
  t1 = pand(t1, vp);
  return paddq(t0, t1);

#endif
}

static INLINE sp_simd_t sp_ntt_sub_partial_simd_core(
    			sp_simd_t a, sp_simd_t b, sp_simd_t vp)
{
#ifdef HAVE_PARTIAL_MOD

  #if SP_TYPE_BITS == 32
  return paddd(psubd(a, b), vp);
  #else
  return paddq(psubq(a, b), vp);
  #endif

#else
  return sp_ntt_sub_simd_core(a, b, vp);

#endif
}

static INLINE sp_simd_t sp_ntt_sub_simd(sp_simd_t a, sp_simd_t b, sp_t p)
{
  #ifdef HAVE_PARTIAL_MOD
  p = 2 * p;
  #endif

  return sp_ntt_sub_simd_core(a, b, 
#if SP_TYPE_BITS == 32
			pshufd(pcvt_i32(p), 0x00));
#else
			pshufd(pcvt_i64(p), 0x44));
#endif
}

static INLINE sp_simd_t sp_ntt_sub_simd0(sp_simd_t a, sp_simd_t b, sp_simd_t p)
{
  #ifdef HAVE_PARTIAL_MOD
    #if SP_TYPE_BITS == 32
    p = paddd(p, p);
    #else
    p = paddq(p, p);
    #endif
  #endif

  return sp_ntt_sub_simd_core(a, b, p);
}

static INLINE sp_simd_t sp_ntt_sub_partial_simd(sp_simd_t a, sp_simd_t b, sp_t p)
{
  #ifdef HAVE_PARTIAL_MOD
  p = 2 * p;
  #endif

  return sp_ntt_sub_partial_simd_core(a, b, 
#if SP_TYPE_BITS == 32
			pshufd(pcvt_i32(p), 0x00));
#else
			pshufd(pcvt_i64(p), 0x44));
#endif
}

static INLINE sp_simd_t sp_ntt_sub_partial_simd0(sp_simd_t a, sp_simd_t b, sp_simd_t p)
{
  #ifdef HAVE_PARTIAL_MOD
    #if SP_TYPE_BITS == 32
    p = paddd(p, p);
    #else
    p = paddq(p, p);
    #endif
  #endif

  return sp_ntt_sub_partial_simd_core(a, b, p);
}


ATTRIBUTE_ALWAYS_INLINE
static INLINE sp_simd_t sp_ntt_mul_simd_core(
		sp_simd_t a, 
		sp_simd_t vw, sp_simd_t vw2,
		sp_simd_t vwi, sp_simd_t vwi2, 
		sp_simd_t vp, sp_simd_t vp2)
{
#if SP_TYPE_BITS == 32

  sp_simd_t t0, t1, t2, t3;

  t0 = pmuludq(a, vwi);
  t1 = pshufd(a, 0x31);
  t2 = pmuludq(t1, vwi2);

  t3 = pmuludq(a, vw);
  t1 = pmuludq(t1, vw2);

  t0 = psrlq(t0, 32);
  t2 = psrlq(t2, 32);
  t0 = pmuludq(t0, vp);
  t2 = pmuludq(t2, vp2);

  t3 = psubq(t3, t0);
  t1 = psubq(t1, t2);

  t3 = pshufd(t3, 0x08);
  t1 = pshufd(t1, 0x08);
  t3 = punpcklo32(t3, t1);

  #ifdef HAVE_PARTIAL_MOD
  return t3;
  #else
  return sp_ntt_sub_simd_core(t3, vp, vp);
  #endif

#elif GMP_LIMB_BITS == 32   /* 64-bit sp_t on a 32-bit machine */

  sp_simd_t t0, t1, t2, t3, t4, t5;

  sp_simd_t vmask = pshufd(pcvt_i32(0xffffffff), 0x44);
  sp_simd_t a2 = pshufd(a, 0xf5);

  t0 = pmuludq(a, vwi);
  t1 = pmuludq(a, vwi2);
  t2 = pmuludq(a2, vwi);
  t3 = pmuludq(a2, vwi2);

  t0 = psrlq(t0, 32);
  t4 = pand(t1, vmask);
  t1 = psrlq(t1, 32);
  t5 = pand(t2, vmask);
  t2 = psrlq(t2, 32);

  t5 = paddq(t5, t0);
  t5 = paddq(t5, t4);
  t5 = psrlq(t5, 32);
  t3 = paddq(t3, t1);
  t3 = paddq(t3, t2);
  t3 = paddq(t3, t5);  /* t3 = hi64(a * winv) */

  t0 = pmuludq(a, vw);
  t1 = pmuludq(a, vw2);
  t2 = pmuludq(a2, vw);

  t2 = paddq(t2, t1);
  t2 = psllq(t2, 32);
  t2 = paddq(t2, t0); /* t2 = lo64(a * w) */

  t0 = pshufd(t3, 0xf5);
  t1 = pmuludq(t3, vp);
  t4 = pmuludq(t3, vp2);
  t5 = pmuludq(t0, vp);

  t5 = paddq(t5, t4);
  t5 = psllq(t5, 32);
  t1 = paddq(t1, t5); /* t1 = lo64(t3 * p) */

  t2 = psubq(t2, t1);
  #ifdef HAVE_PARTIAL_MOD
  return t2;
  #else
  return sp_ntt_sub_simd_core(t2, vp, vp);
  #endif

#else
	/* not supported, don't call */
#endif
}

ATTRIBUTE_ALWAYS_INLINE
static INLINE sp_simd_t sp_ntt_mul_simd(
				sp_simd_t a, sp_t w, sp_t w_inv, sp_t p)
{
#if SP_TYPE_BITS == 32

  return sp_ntt_mul_simd_core(a, 
                        pshufd(pcvt_i32(w), 0x00),
                        pshufd(pcvt_i32(w), 0x00),
                        pshufd(pcvt_i32(w_inv), 0x00),
                        pshufd(pcvt_i32(w_inv), 0x00),
                        pshufd(pcvt_i32(p), 0x00),
                        pshufd(pcvt_i32(p), 0x00));

#elif GMP_LIMB_BITS == 32

  return sp_ntt_mul_simd_core(a,
                        pshufd(pcvt_i64(w), 0x44),
                        pshufd(pcvt_i64(w), 0x55),
                        pshufd(pcvt_i64(w_inv), 0x44),
                        pshufd(pcvt_i64(w_inv), 0x55),
                        pshufd(pcvt_i64(p), 0x44),
                        pshufd(pcvt_i64(p), 0x55));

#else

  /* use CPU 64-bit multiplier */

  sp_simd_t t0, t1;
  sp_t a0, a1;

  pstore_i64(a0, a);
  pstore_i64(a1, pshufd(a, 0x0e));

  a0 = sp_ntt_mul(a0, w, w_inv, p);
  a1 = sp_ntt_mul(a1, w, w_inv, p);

  t0 = pcvt_i64(a0);
  t1 = pcvt_i64(a1);
  return punpcklo64(t0, t1);

#endif
}

ATTRIBUTE_ALWAYS_INLINE
static INLINE sp_simd_t sp_ntt_mul_simd0(
				sp_simd_t a, sp_simd_t *c, sp_simd_t p)
{
#if SP_TYPE_BITS == 32

  return sp_ntt_mul_simd_core(a, 
                        pload(c),
                        pshufd(pload(c), 0x31),
                        pload(c+1),
                        pshufd(pload(c+1), 0x31),
                        p,
                        pshufd(p, 0x31));

#elif GMP_LIMB_BITS == 32

  return sp_ntt_mul_simd_core(a,
                        pload(c),
                        pshufd(pload(c), 0xf5),
                        pload(c+1),
                        pshufd(pload(c+1), 0xf5),
                        p,
                        pshufd(p, 0xf5));

#else

  /* use CPU 64-bit multiplier */

  sp_simd_t t0, t1;
  sp_t a0, a1;
  sp_t *pscalar = (sp_t *)&p;
  sp_t *cscalar = (sp_t *)c;

  pstore_i64(a0, a);
  pstore_i64(a1, pshufd(a, 0x0e));

  a0 = sp_ntt_mul(a0, cscalar[0], cscalar[2], pscalar[0]);
  a1 = sp_ntt_mul(a1, cscalar[1], cscalar[3], pscalar[1]);

  t0 = pcvt_i64(a0);
  t1 = pcvt_i64(a1);
  return punpcklo64(t0, t1);

#endif
}


ATTRIBUTE_ALWAYS_INLINE
static INLINE sp_simd_t sp_ntt_twiddle_mul_simd(sp_simd_t a, 
					sp_simd_t *w, sp_t p)
{
#if SP_TYPE_BITS == 32

  return sp_ntt_mul_simd_core(a, 
                        pload(w),
                        pshufd(pload(w), 0x31),
                        pload(w+1),
                        pshufd(pload(w+1), 0x31),
                        pshufd(pcvt_i32(p), 0x00),
                        pshufd(pcvt_i32(p), 0x00));

#elif GMP_LIMB_BITS == 32

  return sp_ntt_mul_simd_core(a, 
                        pload(w),
                        pshufd(pload(w), 0xf5),
                        pload(w+1),
                        pshufd(pload(w+1), 0xf5),
                        pshufd(pcvt_i64(p), 0x44),
                        pshufd(pcvt_i64(p), 0x55));

#else

  /* use CPU 64-bit multiplier */

  sp_simd_t t0, t1;
  sp_t a0, a1;
  sp_t *wscalar = (sp_t *)w;

  pstore_i64(a0, a);
  pstore_i64(a1, pshufd(a, 0x0e));

  a0 = sp_ntt_mul(a0, wscalar[0], wscalar[2], p);
  a1 = sp_ntt_mul(a1, wscalar[1], wscalar[3], p);

  t0 = pcvt_i64(a0);
  t1 = pcvt_i64(a1);
  return punpcklo64(t0, t1);

#endif
}

ATTRIBUTE_ALWAYS_INLINE
static INLINE sp_simd_t sp_ntt_twiddle_mul_simd0(
				sp_simd_t a, sp_simd_t *c, sp_simd_t p)
{
  return sp_ntt_mul_simd0(a, c, p);
}

#ifdef __cplusplus
}
#endif

#endif /* _NTT_IMPL_SSE2_H */
