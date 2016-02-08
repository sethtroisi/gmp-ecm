#ifndef _NTT_IMPL_AVX2_H
#define _NTT_IMPL_AVX2_H

#include "ntt-impl-scalar.h"
#include <immintrin.h>

#ifdef __cplusplus
extern "C" {
#endif

#define HAVE_SIMD
#define SP_SIMD_VSIZE (256 / SP_TYPE_BITS)

/* mangled routine names */
#define V(name) MANGLE_AVX2(X(name))
#define SP_SIMD_NAME_SUFFIX_STR SP_NAME_SUFFIX_STR "avx2"

#define pload(addr)  _mm256_load_si256((__m256i const *)(addr))
#define ploadu(addr)  _mm256_lddqu_si256((__m256i const *)(addr))
#define pstore(x, addr) _mm256_store_si256((__m256i *)(addr), x)
#define pstoreu(x, addr) _mm256_storeu_si256((__m256i *)(addr), x)

#define pand _mm256_and_si256
#define pxor _mm256_xor_si256
#define psetzero() _mm256_setzero_si256()
#define paddd _mm256_add_epi32
#define paddq _mm256_add_epi64
#define psubd _mm256_sub_epi32
#define psubq _mm256_sub_epi64
#define pmind _mm256_min_epu32
#define pmuludq _mm256_mul_epu32
#define pmullod _mm256_mullo_epi32
#define pslld _mm256_slli_epi32
#define psllq _mm256_slli_epi64
#define psrld _mm256_srli_epi32
#define psrlq _mm256_srli_epi64
#define pshufd _mm256_shuffle_epi32
#define pcmpgtd _mm256_cmpgt_epi32
#define pcmpgtq _mm256_cmpgt_epi64
#define punpcklo32 _mm256_unpacklo_epi32
#define punpcklo64 _mm256_unpacklo_epi64
#define pbroadcast32 _mm256_set1_epi32
#define pbroadcast64 _mm256_set1_epi64x
#define pcvt_i32(x) _mm256_set_epi32(0,0,0,0,0,0,0,(int)(x))

/* using FPU blend operation on 64-bit integer bitfields 
   works correctly and seems to run faster */
#define pcmovq(a,b,c) (__m256i)_mm256_blendv_pd((__m256d)a, (__m256d)b, (__m256d)c)


typedef __m256i sp_simd_t;

#define sp_simd_load(x) pload(x)

static INLINE sp_simd_t sp_simd_gather(spv_t x, spv_size_t dist,
    					spv_size_t vsize)
{
  return ploadu(x); /* XXX wrong */
}

static INLINE sp_simd_t sp_simd_pfa_gather(spv_t x, spv_size_t start_off, 
					spv_size_t inc, spv_size_t n,
					spv_size_t vsize)
{
  return ploadu(x); /* XXX wrong */
}

#define sp_simd_store(t, x) pstore(t, x)

static INLINE void sp_simd_scatter(sp_simd_t t, spv_t x, spv_size_t dist,
    					spv_size_t vsize)
{
  pstoreu(t, x); /* XXX wrong */
}

static INLINE void sp_simd_pfa_scatter(sp_simd_t t, spv_t x, 
    				spv_size_t start_off, 
				spv_size_t inc, spv_size_t n,
				spv_size_t vsize)
{
  pstoreu(t, x); /* XXX wrong */
}

static INLINE sp_simd_t sp_ntt_add_simd_core(
    			sp_simd_t a, sp_simd_t b, sp_simd_t vp)
{
  sp_simd_t t0, t1;

#if SP_TYPE_BITS == 32
  t0 = paddd(a, b);
  t1 = psubd(t0, vp);
  return pmind(t0, t1);

#else
  t0 = paddq(a, b);
  t1 = psubq(t0, vp);
  return pcmovq(t1, t0, t1);

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
			pbroadcast32(p));
#else
			pbroadcast64(p));
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
			pbroadcast32(p));
#else
			pbroadcast64(p));
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
  t1 = paddd(t0, vp);
  return pmind(t0, t1);

#else
  t0 = psubq(a, b);
  t1 = paddq(t0, vp);
  return pcmovq(t0, t1, t0);

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
			pbroadcast32(p));
#else
			pbroadcast64(p));
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
			pbroadcast32(p));
#else
			pbroadcast64(p));
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
			sp_simd_t a, sp_simd_t vw, 
			sp_simd_t vwi, sp_simd_t vwi2, 
			sp_simd_t vp)
{
#if SP_TYPE_BITS == 32

  sp_simd_t t0, t1, t2;

  t0 = pmuludq(a, vwi);
  t1 = pshufd(a, 0x31);
  t2 = pmuludq(t1, vwi2);
  t0 = pshufd(t0, 0x0d);
  t2 = pshufd(t2, 0x0d);
  t0 = punpcklo32(t0, t2);

  t0 = psubd(pmullod(a, vw), pmullod(t0, vp));

  #ifdef HAVE_PARTIAL_MOD
  return t0;
  #else
  return sp_ntt_sub_simd_core(t0, vp, vp);
  #endif

#else

  sp_simd_t t0, t1, t2, t3, t4, t5;

  sp_simd_t vmask = _mm256_set_epi32(0,-1,0,-1,0,-1,0,-1);
  sp_simd_t vp2 = pshufd(vp, 0xf5);
  sp_simd_t vw2 = pshufd(vw, 0xf5);
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

#endif
}

ATTRIBUTE_ALWAYS_INLINE
static INLINE sp_simd_t sp_ntt_mul_simd(
				sp_simd_t a, sp_t w, sp_t w_inv, sp_t p)
{
#if SP_TYPE_BITS == 32

  return sp_ntt_mul_simd_core(a,
		      pbroadcast32(w),
		      pbroadcast32(w_inv),
		      pbroadcast32(w_inv),
		      pbroadcast32(p));

#else

  return sp_ntt_mul_simd_core(a,
		      pbroadcast64(w),
		      pbroadcast64(w_inv),
		      pshufd(pbroadcast64(w_inv), 0x11),
		      pbroadcast64(p));
#endif
}

ATTRIBUTE_ALWAYS_INLINE
static INLINE sp_simd_t sp_ntt_mul_simd0(
				sp_simd_t a, sp_simd_t *c, sp_simd_t p)
{
#if SP_TYPE_BITS == 32

  return sp_ntt_mul_simd_core(a,
		      pload(c + 0),
		      pload(c + 1),
		      pshufd(pload(c + 1), 0x31),
		      p);

#else

  return sp_ntt_mul_simd_core(a,
		      pload(c + 0),
		      pload(c + 1),
		      pshufd(pload(c + 1), 0xf5),
		      p);
#endif
}

ATTRIBUTE_ALWAYS_INLINE
static INLINE sp_simd_t sp_ntt_twiddle_mul_simd(sp_simd_t a, 
					sp_simd_t *w, sp_t p)
{
#if SP_TYPE_BITS == 32

  return sp_ntt_mul_simd_core(a, 
      			pload(w), 
			pload(w + 1),
			pshufd(pload(w + 1), 0x31),
			pbroadcast32(p));

#else

  return sp_ntt_mul_simd_core(a, 
      			pload(w), 
			pload(w + 1),
			pshufd(pload(w + 1), 0xf5),
			pbroadcast64(p));
#endif
}

ATTRIBUTE_ALWAYS_INLINE
static INLINE sp_simd_t sp_ntt_twiddle_mul_simd0(
				sp_simd_t a, sp_simd_t *w, sp_simd_t p)
{
  return sp_ntt_mul_simd0(a, w, p);
}

#ifdef __cplusplus
}
#endif

#endif /* _NTT_IMPL_AVX2_H */
