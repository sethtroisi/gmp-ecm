#ifndef _NTT_IMPL_AVX_H
#define _NTT_IMPL_AVX_H

#include "ntt-impl-scalar.h"
#include <immintrin.h>

#ifdef __cplusplus
extern "C" {
#endif

#define HAVE_SIMD
#define SP_SIMD_VSIZE (256 / SP_TYPE_BITS)

/* mangled routine names */
#define V(name) MANGLE_AVX(X(name))
#define SP_SIMD_NAME_SUFFIX_STR SP_NAME_SUFFIX_STR "avx"

#define pload(addr)  _mm256_load_pd((double const *)(addr))
#define ploadu(addr)  _mm256_loadu_pd((double const *)(addr))
#define pload_lo64(addr)  _mm_load_sd((double const *)(addr))
#define pload_hi64(x, addr)  _mm_loadh_pd((__m128d)x, (double const *)(addr))
#define pstore(x, addr) _mm256_store_pd((double *)(addr), x)
#define pstoreu(x, addr) _mm256_storeu_pd((double *)(addr), x)
#define padd _mm256_add_pd
#define psub _mm256_sub_pd
#define pmul _mm256_mul_pd
#define psetzero _mm256_setzero_pd
#define pfloor _mm256_floor_pd
#define pbroadcast _mm256_set1_pd
#define pcmov _mm256_blendv_pd

typedef __m256d sp_simd_t;

#define sp_simd_load(x) pload(x)

static INLINE sp_simd_t sp_simd_gather(spv_t x, spv_size_t dist,
    					spv_size_t vsize)
{
  switch (vsize)
    {
      case 4:
	if (dist == 1)
	  return ploadu(x);
	else
	  {
	    __m128d t0, t1;
	    sp_simd_t res;

	    t0 = _mm_load_sd(x);
	    t0 = _mm_loadh_pd(t0, x + dist);
	    t1 = _mm_load_sd(x + 2 * dist);
	    t1 = _mm_loadh_pd(t1, x + 3 * dist);
	    res = _mm256_castpd128_pd256(t0);
	    return _mm256_insertf128_pd(res, t1, 1);
	  }

      case 3:
	{
      	  __m128d t0, t1;
	  sp_simd_t res;

	  t0 = _mm_load_sd(x);
	  t0 = _mm_loadh_pd(t0, x + dist);
	  t1 = _mm_load_sd(x + 2 * dist);
	  res = _mm256_castpd128_pd256(t0);
	  return _mm256_insertf128_pd(res, t1, 1);
	}

      case 2:
	{
	  __m128d t0;

	  t0 = _mm_load_sd(x);
	  t0 = _mm_loadh_pd(t0, x + dist);
	  return _mm256_insertf128_pd(psetzero(), t0, 0);
	}

      default:
	{
      	  return _mm256_insertf128_pd(psetzero(), _mm_load_sd(x), 0);
	}
    }
}

static INLINE sp_simd_t sp_simd_pfa_gather(spv_t x, spv_size_t start_off, 
					spv_size_t inc, spv_size_t n,
					spv_size_t vsize)
{
  switch (vsize)
    {
      case 4:
	{
      	  __m128d t0, t1;
	  sp_simd_t res;
	  spv_size_t j0 = start_off;
	  spv_size_t j1 = sp_array_inc(j0, inc, n);
	  spv_size_t j2 = sp_array_inc(j0, 2 * inc, n);
	  spv_size_t j3 = sp_array_inc(j0, 3 * inc, n);

	  t0 = _mm_load_sd(x + j0);
	  t0 = _mm_loadh_pd(t0, x + j1);
	  t1 = _mm_load_sd(x + j2);
	  t1 = _mm_loadh_pd(t1, x + j3);
	  res = _mm256_castpd128_pd256(t0);
	  return _mm256_insertf128_pd(res, t1, 1);
	}

      case 3:
	{
      	  __m128d t0, t1;
	  sp_simd_t res;
	  spv_size_t j0 = start_off;
	  spv_size_t j1 = sp_array_inc(j0, inc, n);
	  spv_size_t j2 = sp_array_inc(j0, 2 * inc, n);

	  t0 = _mm_load_sd(x + j0);
	  t0 = _mm_loadh_pd(t0, x + j1);
	  t1 = _mm_load_sd(x + j2);
	  res = _mm256_castpd128_pd256(t0);
	  return _mm256_insertf128_pd(res, t1, 1);
	}

      case 2:
	{
	  __m128d t0;
	  spv_size_t j0 = start_off;
	  spv_size_t j1 = sp_array_inc(j0, inc, n);

	  t0 = _mm_load_sd(x + j0);
	  t0 = _mm_loadh_pd(t0, x + j1);
	  return _mm256_insertf128_pd(psetzero(), t0, 0);
	}

      default:
	{
	  spv_size_t j0 = start_off;
      	  return _mm256_insertf128_pd(psetzero(), _mm_load_sd(x + j0), 0);
	}
    }
}

#define sp_simd_store(t, x) pstore(t, x);

static INLINE void sp_simd_scatter(sp_simd_t t, spv_t x, spv_size_t dist,
    					spv_size_t vsize)
{
  switch (vsize)
    {
      case 4:
	if (dist == 1)
	  pstoreu(t, x);
	else
	  {
	    __m128d hi = _mm256_extractf128_pd(t, 1);
	    __m128d lo = _mm256_castpd256_pd128(t);
	    _mm_store_sd(x, lo);
	    _mm_storeh_pd(x + dist, lo);
	    _mm_store_sd(x + 2 * dist, hi);
	    _mm_storeh_pd(x + 3 * dist, hi);
	  }
	break;

      case 3:
	{
	  __m128d hi = _mm256_extractf128_pd(t, 1);
	  __m128d lo = _mm256_castpd256_pd128(t);
      	  _mm_store_sd(x, lo);
	  _mm_storeh_pd(x + dist, lo);
	  _mm_store_sd(x + 2 * dist, hi);
      	}
	break;

      case 2:
	{
	  __m128d lo = _mm256_castpd256_pd128(t);
      	  _mm_store_sd(x, lo);
	  _mm_storeh_pd(x + dist, lo);
      	}
	break;

      default:
	{
	  __m128d lo = _mm256_castpd256_pd128(t);
      	  _mm_store_sd(x, lo);
	}
	break;
    }
}

static INLINE void sp_simd_pfa_scatter(sp_simd_t t, spv_t x, 
    				spv_size_t start_off, 
				spv_size_t inc, spv_size_t n,
				spv_size_t vsize)
{
  switch (vsize)
    {
      case 4:
	{
      	  __m128d hi = _mm256_extractf128_pd(t, 1);
	  __m128d lo = _mm256_castpd256_pd128(t);
	  spv_size_t j0 = start_off;
	  spv_size_t j1 = sp_array_inc(j0, inc, n);
	  spv_size_t j2 = sp_array_inc(j0, 2 * inc, n);
	  spv_size_t j3 = sp_array_inc(j0, 3 * inc, n);

	  _mm_store_sd(x + j0, lo);
	  _mm_storeh_pd(x + j1, lo);
	  _mm_store_sd(x + j2, hi);
	  _mm_storeh_pd(x + j3, hi);
	}
	break;

      case 3:
	{
	  __m128d hi = _mm256_extractf128_pd(t, 1);
	  __m128d lo = _mm256_castpd256_pd128(t);
	  spv_size_t j0 = start_off;
	  spv_size_t j1 = sp_array_inc(j0, inc, n);
	  spv_size_t j2 = sp_array_inc(j0, 2 * inc, n);

      	  _mm_store_sd(x + j0, lo);
	  _mm_storeh_pd(x + j1, lo);
	  _mm_store_sd(x + j2, hi);
      	}
	break;

      case 2:
	{
	  __m128d lo = _mm256_castpd256_pd128(t);
	  spv_size_t j0 = start_off;
	  spv_size_t j1 = sp_array_inc(j0, inc, n);

      	  _mm_store_sd(x + j0, lo);
	  _mm_storeh_pd(x + j1, lo);
      	}
	break;

      default:
	{
	  __m128d lo = _mm256_castpd256_pd128(t);
	  spv_size_t j0 = start_off;

      	  _mm_store_sd(x + j0, lo);
	}
	break;
    }
}

static INLINE sp_simd_t sp_ntt_add_simd_core(
    			sp_simd_t a, sp_simd_t b, sp_simd_t vp)
{
  sp_simd_t t0 = padd(a, b);
  sp_simd_t t1 = psub(t0, vp);
  return pcmov(t1, t0, t1);
}

static INLINE sp_simd_t sp_ntt_add_simd(sp_simd_t a, sp_simd_t b, sp_t p)
{
  return sp_ntt_add_simd_core(a, b, pbroadcast(p));
}

static INLINE sp_simd_t sp_ntt_add_simd0(sp_simd_t a, sp_simd_t b, sp_simd_t p)
{
  return sp_ntt_add_simd_core(a, b, p);
}

static INLINE sp_simd_t sp_ntt_add_partial_simd(sp_simd_t a, sp_simd_t b, sp_t p)
{
  return sp_ntt_add_simd_core(a, b, pbroadcast(p));
}

static INLINE sp_simd_t sp_ntt_add_partial_simd0(sp_simd_t a, sp_simd_t b, sp_simd_t p)
{
  return sp_ntt_add_simd_core(a, b, p);
}

static INLINE sp_simd_t sp_ntt_sub_simd_core(
    			sp_simd_t a, sp_simd_t b, sp_simd_t vp)
{
  sp_simd_t t0 = psub(a, b);
  sp_simd_t t1 = padd(t0, vp);
  return pcmov(t0, t1, t0);
}

static INLINE sp_simd_t sp_ntt_sub_simd(sp_simd_t a, sp_simd_t b, sp_t p)
{
  return sp_ntt_sub_simd_core(a, b, pbroadcast(p));
}

static INLINE sp_simd_t sp_ntt_sub_simd0(sp_simd_t a, sp_simd_t b, sp_simd_t p)
{
  return sp_ntt_sub_simd_core(a, b, p);
}

static INLINE sp_simd_t sp_ntt_sub_partial_simd(sp_simd_t a, sp_simd_t b, sp_t p)
{
  return sp_ntt_sub_simd_core(a, b, pbroadcast(p));
}

static INLINE sp_simd_t sp_ntt_sub_partial_simd0(sp_simd_t a, sp_simd_t b, sp_simd_t p)
{
  return sp_ntt_sub_simd_core(a, b, p);
}


ATTRIBUTE_ALWAYS_INLINE
static INLINE sp_simd_t sp_ntt_mul_simd_core_recip(
				sp_simd_t a, sp_simd_t w, sp_simd_t whi, 
				sp_simd_t p, sp_simd_t vrecip)
{
  sp_simd_t t0, t1, t2;
  sp_simd_t ahi, alo;
  sp_simd_t awhi, awlo;
  sp_simd_t q, qhi, qlo;
  sp_simd_t pqhi, pqlo;
  sp_simd_t phi, plo;
  sp_simd_t splitter = pbroadcast(134217729.0);
  sp_simd_t wlo = psub(w, whi);

  t0 = pmul(splitter, a);
  ahi = psub(t0, psub(t0, a));
  alo = psub(a, ahi);

  awhi = pmul(a, w);
  t0 = psub(pmul(ahi, whi), awhi);
  t0 = padd(t0, pmul(ahi, wlo));
  t0 = padd(t0, pmul(alo, whi));
  awlo = padd(t0, pmul(alo, wlo));

  q = pfloor(pmul(awhi, vrecip));
  t0 = pmul(splitter, q);
  qhi = psub(t0, psub(t0, q));
  qlo = psub(q, qhi);

  t0 = pmul(splitter, p);
  phi = psub(t0, psub(t0, p));
  plo = psub(p, phi);

  pqhi = pmul(p, q);
  t0 = psub(pmul(phi, qhi), pqhi);
  t0 = padd(t0, pmul(phi, qlo));
  t0 = padd(t0, pmul(plo, qhi));
  pqlo = padd(t0, pmul(plo, qlo));

  t0 = psub(awhi, pqhi);
  t1 = psub(t0, awhi);
  t1 = padd(psub(awhi, psub(t0, t1)), padd(pqhi, t1));
  t1 = psub(padd(t1, awlo), pqlo);
  t0 = padd(t0, t1);

  t1 = psub(t0, p);
  t0 = pcmov(t1, t0, t1);
  t1 = padd(t0, p);
  t0 = pcmov(t0, t1, t0);
  return t0;
}

ATTRIBUTE_ALWAYS_INLINE
static INLINE sp_simd_t sp_ntt_mul_simd(
				sp_simd_t a, sp_t w, sp_t whi, sp_t p)
{
  return sp_ntt_mul_simd_core_recip(a, pbroadcast(w),
      			pbroadcast(whi), pbroadcast(p),
			pbroadcast(1.0 / p));
}

ATTRIBUTE_ALWAYS_INLINE
static INLINE sp_simd_t sp_ntt_mul_simd0(sp_simd_t a, 
					sp_simd_t *w, sp_simd_t p)
{
  sp_simd_t vrecip = _mm256_div_pd(pbroadcast(1.0), p);

  return sp_ntt_mul_simd_core_recip(a, pload(w), pload(w+1), p, vrecip);
}

ATTRIBUTE_ALWAYS_INLINE
static INLINE sp_simd_t sp_ntt_twiddle_mul_simd(sp_simd_t a, 
					sp_simd_t *w, sp_t p)
{
  return sp_ntt_mul_simd_core_recip(a, pload(w), pload(w+1),
      			pbroadcast(p), pbroadcast(1.0 / p));
}

ATTRIBUTE_ALWAYS_INLINE
static INLINE sp_simd_t sp_ntt_twiddle_mul_simd0(sp_simd_t a, 
					sp_simd_t *w, sp_simd_t p)
{
  return sp_ntt_mul_simd0(a, w, p);
}

#ifdef __cplusplus
}
#endif

#endif /* _NTT_IMPL_AVX_H */
