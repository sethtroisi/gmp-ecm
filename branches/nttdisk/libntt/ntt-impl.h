#ifndef _NTT_IMPL_H
#define _NTT_IMPL_H

#include "sp.h"

typedef const uint8_t * (*get_fixed_ntt_const_t)(void);

typedef void (*nttdata_init_t)(spv_t out, 
				sp_t p, sp_t d,
				sp_t primroot, sp_t order);

typedef void (*ntt_run_t)(spv_t x, spv_size_t num_transforms, 
			  sp_t p, spv_t ntt_const);

typedef void (*ntt_pfa_run_t)(spv_t x, spv_size_t cofactor, 
			  sp_t p, spv_t ntt_const);

typedef void (*ntt_twiddle_run_t)(spv_t x, spv_t w,
			  spv_size_t stride,
			  spv_size_t num_transforms, 
			  sp_t p, spv_t ntt_const);


/* a copy of sp_add, but operating on array offsets. These
   are always size_t types, and may have a size different 
   from an sp_t */

static inline spv_size_t sp_array_inc(spv_size_t a, spv_size_t b, spv_size_t m) 
{
#if (defined(__GNUC__) || defined(__ICL)) && \
    (defined(__x86_64__) || defined(__i386__))

  spv_size_t t = a - m, tr = a + b;

  __asm__ (
    "add %2, %1    # sp_array_inc: t += b\n\t"
    "cmovc %1, %0  # sp_array_inc: if (cy) tr = t \n\t"
    : "+r" (tr), "+&r" (t)
    : "g" (b)
    : "cc"
  );

  return tr;

#elif defined(_MSC_VER) && !defined(_WIN64)

  __asm
    {
        mov     eax, a
        add     eax, b
        mov     edx, eax
        sub     edx, m
        cmovnc  eax, edx
    }

#else

  spv_size_t t = a + b;
  if (t >= m)
    t -= m;
  return t;

#endif
}

/* if the modulus is 2 bits or more smaller than the machine
   word size, the core NTT routines use a redundant representation 
   of the transform elements. Modular multiplies do not do their
   final modular reduction, and inputs to the multiplies are not 
   reduced mod p either. The transform results are in [0, 2p) */

#if SP_NUMB_BITS < SP_TYPE_BITS - 1
#define HAVE_PARTIAL_MOD
#endif

static inline sp_t sp_ntt_add(sp_t a, sp_t b, sp_t p)
{
#ifdef HAVE_PARTIAL_MOD
  return sp_add(a, b, 2 * p);
#else
  return sp_add(a, b, p);
#endif
}

static inline sp_t sp_ntt_add_partial(sp_t a, sp_t b, sp_t p)
{
#ifdef HAVE_PARTIAL_MOD
  return a + b;
#else
  return sp_ntt_add(a, b, p);
#endif
}

static inline sp_t sp_ntt_sub(sp_t a, sp_t b, sp_t p)
{
#ifdef HAVE_PARTIAL_MOD
  return sp_sub(a, b, 2 * p);
#else
  return sp_sub(a, b, p);
#endif
}

static inline sp_t sp_ntt_sub_partial(sp_t a, sp_t b, sp_t p)
{
#ifdef HAVE_PARTIAL_MOD
  return a - b + 2 * p;
#else
  return sp_ntt_sub(a, b, p);
#endif
}

/*--------------------------------------------------------------------*/
/* perform modular multiplication when the multiplier
   and its generalized inverse are both known and precomputed */

sp_t sp_ntt_reciprocal (sp_t w, sp_t p);

/* low half multiply */

static inline sp_t sp_mul_lo(sp_t a, sp_t b)
{
#if SP_TYPE_BITS <= GMP_LIMB_BITS

  return a * b;

#else  /* build product from smaller multiplies */

  mp_limb_t a1 = (mp_limb_t)((a) >> GMP_LIMB_BITS);
  mp_limb_t a0 = (mp_limb_t)(a);
  mp_limb_t b1 = (mp_limb_t)((b) >> GMP_LIMB_BITS);
  mp_limb_t b0 = (mp_limb_t)(b);
  mp_limb_t a0b0_hi, a0b0_lo;

  umul_ppmm(a0b0_hi, a0b0_lo, a0, b0);

  return (sp_t)(a0b0_hi + a0 * b1 + 
      		a1 * b0) << GMP_LIMB_BITS | a0b0_lo;
#endif
}

/* high half multiply */

static inline sp_t sp_mul_hi(sp_t a, sp_t b)
{
#if SP_TYPE_BITS == GMP_LIMB_BITS  /* use GMP functions */

  mp_limb_t hi;
  ATTRIBUTE_UNUSED mp_limb_t lo;
  umul_ppmm(hi, lo, a, b);

  return hi;

#elif SP_TYPE_BITS < GMP_LIMB_BITS  /* ordinary multiply */

  mp_limb_t prod = (mp_limb_t)(a) * (mp_limb_t)(b);
  return (sp_t)(prod >> SP_TYPE_BITS);

#else  /* build product from smaller multiplies */

  mp_limb_t a1 = (mp_limb_t)((a) >> GMP_LIMB_BITS);
  mp_limb_t a0 = (mp_limb_t)(a);
  mp_limb_t b1 = (mp_limb_t)((b) >> GMP_LIMB_BITS);
  mp_limb_t b0 = (mp_limb_t)(b);
  mp_limb_t a0b0_hi, a0b0_lo;
  mp_limb_t a0b1_hi, a0b1_lo;
  mp_limb_t a1b0_hi, a1b0_lo;
  mp_limb_t a1b1_hi, a1b1_lo;
  mp_limb_t cy;

  umul_ppmm(a0b0_hi, a0b0_lo, a0, b0);
  umul_ppmm(a0b1_hi, a0b1_lo, a0, b1);
  umul_ppmm(a1b0_hi, a1b0_lo, a1, b0);
  umul_ppmm(a1b1_hi, a1b1_lo, a1, b1);

  add_ssaaaa(cy, a0b0_hi,
      	     0, a0b0_hi,
	     0, a0b1_lo);
  add_ssaaaa(cy, a0b0_hi,
      	     cy, a0b0_hi,
	     0, a1b0_lo);
  add_ssaaaa(a1b1_hi, a1b1_lo,
      	     a1b1_hi, a1b1_lo,
	     0, cy);
  add_ssaaaa(a1b1_hi, a1b1_lo,
      	     a1b1_hi, a1b1_lo,
	     0, a0b1_hi);
  add_ssaaaa(a1b1_hi, a1b1_lo,
      	     a1b1_hi, a1b1_lo,
	     0, a1b0_hi);

  return (sp_t)a1b1_hi << GMP_LIMB_BITS | a1b1_lo;

#endif
}

static inline sp_t sp_ntt_mul(sp_t x, sp_t w, sp_t w_inv, sp_t p)
{
  sp_t q = sp_mul_hi(x, w_inv);

#ifdef HAVE_PARTIAL_MOD
  return x * w - q * p;
#else
  sp_t r = x * w - q * p;
  return sp_sub(r, p, p);
#endif
}

/*------------------- definitions for SIMD transforms ----------------*/

#ifdef HAVE_SSE2

typedef __m128i sp_simd_t;

#define SP_SIMD_VSIZE (128 / SP_TYPE_BITS)

static inline sp_simd_t sp_simd_load(spv_t x)
{
#if SP_TYPE_BITS == 32

  return ploadu(x);

#else

  sp_simd_t t = pload_lo64(x);
  return pload_hi64(t, x + 1);

#endif
}

static inline sp_simd_t sp_simd_gather(spv_t x, spv_size_t stride)
{
#if SP_TYPE_BITS == 32

  sp_simd_t t0 = pload_lo32(x + 0);
  sp_simd_t t1 = pload_lo32(x + stride);
  sp_simd_t t2 = pload_lo32(x + 2 * stride);
  sp_simd_t t3 = pload_lo32(x + 3 * stride);
  sp_simd_t r0 = punpcklo32(t0, t1);
  sp_simd_t r1 = punpcklo32(t2, t3);
  return punpcklo64(r0, r1);

#else

  sp_simd_t t = pload_lo64(x + 0);
  return pload_hi64(t, x + stride);

#endif
}

static inline sp_simd_t sp_simd_pfa_gather(spv_t x, spv_size_t start_off, 
					spv_size_t inc, spv_size_t n)
{
#if SP_TYPE_BITS == 32

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

#else

  spv_size_t j0 = start_off;
  spv_size_t j1 = sp_array_inc(j0, inc, n);
  sp_simd_t t = pload_lo64(x + j0);
  return pload_hi64(t, x + j1);

#endif
}

static inline void sp_simd_store(sp_simd_t t, spv_t x)
{
#if SP_TYPE_BITS == 32

  pstoreu(t, x);

#else

  pstore_lo64(t, x);
  pstore_hi64(t, x + 1);

#endif
}

static inline void sp_simd_scatter(sp_simd_t t, spv_t x, spv_size_t stride)
{
#if SP_TYPE_BITS == 32

  pstore_lo32(t, x + 0);
  t = _mm_srli_si128(t, 4);
  pstore_lo32(t, x + stride);
  t = _mm_srli_si128(t, 4);
  pstore_lo32(t, x + 2 * stride);
  t = _mm_srli_si128(t, 4);
  pstore_lo32(t, x + 3 * stride);

#else

  pstore_lo64(t, x + 0);
  pstore_hi64(t, x + stride);

#endif
}

static inline void sp_simd_pfa_scatter(sp_simd_t t, spv_t x, 
    				spv_size_t start_off, 
				spv_size_t inc, spv_size_t n)
{
#if SP_TYPE_BITS == 32

  spv_size_t j0 = start_off;
  spv_size_t j1 = sp_array_inc(j0, inc, n);
  spv_size_t j2 = sp_array_inc(j0, 2 * inc, n);
  spv_size_t j3 = sp_array_inc(j0, 3 * inc, n);
  pstore_lo32(t, x + j0);
  t = _mm_srli_si128(t, 4);
  pstore_lo32(t, x + j1);
  t = _mm_srli_si128(t, 4);
  pstore_lo32(t, x + j2);
  t = _mm_srli_si128(t, 4);
  pstore_lo32(t, x + j3);

#else

  spv_size_t j0 = start_off;
  spv_size_t j1 = sp_array_inc(j0, inc, n);
  pstore_lo64(t, x + j0);
  pstore_hi64(t, x + j1);

#endif
}

static inline sp_simd_t sp_ntt_add_simd(sp_simd_t a, sp_simd_t b, sp_t p)
{
  sp_simd_t vp, t0, t1;

  #ifdef HAVE_PARTIAL_MOD
  p = 2 * p;
  #endif

#if SP_TYPE_BITS == 32
  vp = pshufd(pcvt_i32(p), 0x00);
  t0 = paddd(a, b);
  t0 = psubd(t0, vp);
  t1 = pcmpgtd(psetzero(), t0);
  t1 = pand(t1, vp);
  return paddd(t0, t1);

#else
  vp = pshufd(pcvt_i64(p), 0x44);
  t0 = paddq(a, b);
  t0 = psubq(t0, vp);
  t1 = pcmpgtd(psetzero(), t0);
  t1 = pshufd(t1, 0xf5);
  t1 = pand(t1, vp);
  return paddq(t0, t1);

#endif
}

static inline sp_simd_t sp_ntt_add_partial_simd(sp_simd_t a, sp_simd_t b, sp_t p)
{
#ifdef HAVE_PARTIAL_MOD

  #if SP_TYPE_BITS == 32
  return paddd(a, b);
  #else
  return paddq(a, b);
  #endif

#else

  return sp_ntt_add_simd(a, b, p);

#endif
}

static inline sp_simd_t sp_ntt_sub_simd(sp_simd_t a, sp_simd_t b, sp_t p)
{
  sp_simd_t vp, t0, t1;

  #ifdef HAVE_PARTIAL_MOD
  p = 2 * p;
  #endif

#if SP_TYPE_BITS == 32
  vp = pshufd(pcvt_i32(p), 0x00);
  t0 = psubd(a, b);
  t1 = pcmpgtd(psetzero(), t0);
  t1 = pand(t1, vp);
  return paddd(t0, t1);

#else
  vp = pshufd(pcvt_i64(p), 0x44);
  t0 = psubq(a, b);
  t1 = pcmpgtd(psetzero(), t0);
  t1 = pshufd(t1, 0xf5);
  t1 = pand(t1, vp);
  return paddq(t0, t1);

#endif
}

static inline sp_simd_t sp_ntt_sub_partial_simd(sp_simd_t a, sp_simd_t b, sp_t p)
{
#ifdef HAVE_PARTIAL_MOD

  #if SP_TYPE_BITS == 32

  sp_simd_t vp = pshufd(pcvt_i32(2 * p), 0x00);
  return paddd(psubd(a, b), vp);

  #else

  sp_simd_t vp = pshufd(pcvt_i64(2 * p), 0x44);
  return paddq(psubq(a, b), vp);

  #endif

#else

  return sp_ntt_sub_simd(a, b, p);

#endif
}

static inline sp_simd_t sp_mul_simd(sp_simd_t a, sp_t w,
				sp_t p, sp_t d)
{
#if SP_NUMB_BITS < 32

  sp_simd_t t0, t1, t2, t3, vp, vp2, vd, vw;

  vp = pshufd(pcvt_i32(p), 0x00);
  vp2 = pshufd(pcvt_i32(p), 0x44);
  vw = pshufd(pcvt_i32(w), 0x00);
  vd = pshufd(pcvt_i32(d), 0x00);

  t0 = a;
  t1 = vw;
  t2 = pshufd(t0, 0x31);
  t3 = pshufd(t1, 0x31);
  t0 = pmuludq(t0, t1);
  t2 = pmuludq(t2, t3);
  t1 = psrlq(t0, 2 * SP_NUMB_BITS - SP_TYPE_BITS);
  t3 = psrlq(t2, 2 * SP_NUMB_BITS - SP_TYPE_BITS);
  t1 = pmuludq(t1, vd);
  t3 = pmuludq(t3, vd);

#if SP_NUMB_BITS < 31
  t1 = psrlq(t1, 33);
  t3 = psrlq(t3, 33);
  t1 = pmuludq(t1, vp);
  t3 = pmuludq(t3, vp);
  t0 = psubq(t0, t1);
  t2 = psubq(t2, t3);
#else
  t1 = pshufd(t1, 0xf5);
  t3 = pshufd(t3, 0xf5);
  t1 = pmuludq(t1, vp);
  t3 = pmuludq(t3, vp);
  t0 = psubq(t0, t1);
  t2 = psubq(t2, t3);

  t0 = psubq(t0, vp2);
  t2 = psubq(t2, vp2);
  t1 = pshufd(t0, 0xf5);
  t3 = pshufd(t2, 0xf5);
  t1 = pand(t1, vp2);
  t3 = pand(t3, vp2);
  t0 = paddq(t0, t1);
  t2 = paddq(t2, t3);
#endif

  t0 = pshufd(t0, 0x08);
  t1 = pshufd(t2, 0x08);
  t0 = punpcklo32(t0, t1);

  t0 = psubd(t0, vp);
  t1 = pcmpgtd(psetzero(), t0);
  t1 = pand(t1, vp);
  return paddd(t0, t1);

#elif GMP_LIMB_BITS == 32 /* 64-bit sp_t on a 32-bit machine */

  sp_simd_t t0, t1, t2, t3, t4, t5, vp, vd, vw;

  vp = pshufd(pcvt_i64(p), 0x44);
  vw = pshufd(pcvt_i64(w), 0x44);
  vd = pshufd(pcvt_i64(d), 0x44);

  t4 = a;
  t5 = vw;

  t3 = pshufd(t4, 0x31);
  t3 = pmuludq(t3, t5);
  t2 = pshufd(t5, 0x31);
  t2 = pmuludq(t2, t4);
  t1 = pshufd(t4, 0x31);
  t4 = pmuludq(t4, t5);
  t5 = pshufd(t5, 0x31);
  t5 = pmuludq(t5, t1);
  t3 = paddq(t2, t3);

  t0 = t4;
  t4 = psrlq(t4, 32);
  t1 = t3;
  t3 = psllq(t3, 32);
  t3 = paddq(t3, t0);
  t4 = paddq(t4, t1);
  t4 = psrlq(t4, 32);
  t4 = paddq(t4, t5);

  t0 = t3;
  t4 = psllq(t4, 2*(SP_TYPE_BITS - SP_NUMB_BITS));
  t0 = psrlq(t0, 2*SP_NUMB_BITS - SP_TYPE_BITS);
  t4 = paddq(t4, t0);

  t5 = pshufd(t4, 0x31);
  t5 = pmuludq(t5, vd);
  t2 = pshufd(vd, 0x31);
  t2 = pmuludq(t2, t4);
  t1 = pshufd(t4, 0x31);
  t4 = pmuludq(t4, vd);
  t0 = pshufd(vd, 0x31);
  t1 = pmuludq(t1, t0);
  t5 = paddq(t5, t2);

  t4 = psrlq(t4, 32);
  t5 = paddq(t5, t4);
  t5 = psrlq(t5, 32);
  t1 = paddq(t1, t5);
  t1 = psrlq(t1, 1);

  t5 = pshufd(t1, 0x31);
  t5 = pmuludq(t5, vp);
  t2 = pshufd(vp, 0x31);
  t2 = pmuludq(t2, t1);
  t1 = pmuludq(t1, vp);
  t5 = paddq(t5, t2);
  t5 = psllq(t5, 32);
  t1 = paddq(t1, t5);

  t3 = psubq(t3, t1);

  t0 = psubq(t3, vp);
  t1 = pcmpgtd(psetzero(), t0);
  t1 = pshufd(t1, 0xf5);
  t1 = pand(t1, vp);
  return paddq(t0, t1);

#else

  /* there's no way the SSE2 unit can keep up with a
  64-bit multiplier in the ALU */

  sp_simd_t t0, t1;
  sp_t a0, a1;

  pstore_i64(a0, a);
  pstore_i64(a1, pshufd(a, 0x0e));

  a0 = sp_mul(a0, w, p, d);
  a1 = sp_mul(a1, w, p, d);

  t0 = pcvt_i64(a0);
  t1 = pcvt_i64(a1);
  return punpcklo64(t0, t1);

#endif
}
 

static inline sp_simd_t sp_ntt_mul_simd_core(
				sp_simd_t a, sp_simd_t vw, 
				sp_simd_t vwi, sp_t p)
{
#if SP_TYPE_BITS == 32

  sp_simd_t t0, t1, t2, t3;

  sp_simd_t vp = pshufd(pcvt_i32(p), 0x00);

  t0 = pmuludq(a, vwi);
  t1 = pshufd(a, 0x31);
  t2 = pmuludq(t1, vwi);

  t3 = pmuludq(a, vw);
  t1 = pmuludq(t1, vw);

  t0 = psrlq(t0, 32);
  t2 = psrlq(t2, 32);
  t0 = pmuludq(t0, vp);
  t2 = pmuludq(t2, vp);

  t3 = psubq(t3, t0);
  t1 = psubq(t1, t2);

  t3 = pshufd(t3, 0x08);
  t1 = pshufd(t1, 0x08);
  t3 = punpcklo32(t3, t1);

  #ifdef HAVE_PARTIAL_MOD
  return t3;
  #else
  return sp_ntt_sub_simd(t3, vp, p);
  #endif

#elif GMP_LIMB_BITS == 32   /* 64-bit sp_t on a 32-bit machine */

  sp_simd_t vmask;
  sp_simd_t t0, t1, t2, t3, t4, t5, t6;

  sp_simd_t vp = pshufd(pcvt_i64(p), 0x44);

  t0 = pmuludq(a, vwi);
  t4 = pshufd(a, 0xf5);
  t1 = pmuludq(t4, vwi);
  t5 = pshufd(vwi, 0xf5);
  t2 = pmuludq(t5, a);
  t3 = pmuludq(t4, t5);

  t4 = psrlq(t1, 32);
  t5 = psrlq(t2, 32);
  t3 = paddq(t3, t4);
  t3 = paddq(t3, t5);

  t4 = psrlq(t0, 32);
  t1 = pand(t1, vmask);
  t2 = pand(t2, vmask);
  t4 = paddq(t4, t1);
  t4 = paddq(t4, t2);
  t4 = psrlq(t4, 32);
  t3 = paddq(t3, t4);  /* t3 = hi64(a * winv) */

  t0 = pmuludq(a, vw);
  t4 = pshufd(a, 0xf5);
  t1 = pmuludq(t4, vw);
  t5 = pshufd(vw, 0xf5);
  t2 = pmuludq(t5, a);

  t1 = psllq(t1, 32);
  t2 = psllq(t2, 32);
  t6 = paddq(t0, t1);
  t6 = paddq(t6, t2); /* t6 = lo64(a * w) */

  t0 = pmuludq(t3, vp);
  t4 = pshufd(t3, 0xf5);
  t1 = pmuludq(t4, vp);
  t5 = pshufd(vp, 0xf5);
  t2 = pmuludq(t5, t3);

  t1 = psllq(t1, 32);
  t2 = psllq(t2, 32);
  t0 = paddq(t0, t1);
  t0 = paddq(t0, t2); /* t0 = lo64(t3 * p) */

  t6 = psubq(t6, t0);
  #ifdef HAVE_PARTIAL_MOD
  return t6;
  #else
  return sp_ntt_sub_simd(t6, vp, p);
  #endif

#endif
}

static inline sp_simd_t sp_ntt_mul_simd(sp_simd_t a, sp_t w, 
					sp_t w_inv, sp_t p)
{
#if SP_TYPE_BITS == 32
   return sp_ntt_mul_simd_core(a,
  			pshufd(pcvt_i32(w), 0x00),
			pshufd(pcvt_i32(w_inv), 0x00),
			p);

#elif GMP_LIMB_BITS == 32   /* 64-bit sp_t on a 32-bit machine */
   return sp_ntt_mul_simd_core(a,
			pshufd(pcvt_i64(w), 0x44),
			pshufd(pcvt_i64(w_inv), 0x44);
			p);

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

/* twiddle multiplies get a separate SIMD version; their data
   is assumed to reside in (aligned) memory, with twiddle factors
   separated from their inverses, then packed into sp_simd_t vectors
   and concatenated */

static inline sp_simd_t sp_ntt_twiddle_mul_simd(sp_simd_t a, 
					sp_simd_t *w, sp_t p)
{
#if SP_TYPE_BITS == 32 || GMP_LIMB_BITS == 32
   return sp_ntt_mul_simd_core(a, pload(w), pload(w + 1), p);

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

#endif  /*-----------------------------------------------------------*/

typedef struct
{
  uint32_t size;
  uint32_t num_ntt_const;
  get_fixed_ntt_const_t get_fixed_ntt_const;
  nttdata_init_t nttdata_init;
  ntt_run_t ntt_run;
  ntt_pfa_run_t ntt_pfa_run;
  ntt_twiddle_run_t ntt_twiddle_run;
} nttconfig_t;

/* function pointers and precomputed data needed by
   all versions of a codelet */

typedef struct
{
  const nttconfig_t *config;
  spv_t ntt_const;
} codelet_data_t;


/* an NTT is built up of one or more passes through
   the input data */

typedef enum
{
  PASS_TYPE_DIRECT,
  PASS_TYPE_PFA,
  PASS_TYPE_TWIDDLE,
  PASS_TYPE_TWIDDLE_PRE
} pass_type_t;

#define MAX_PFA_CODELETS 6
#define MAX_PASSES 10

typedef struct
{
  pass_type_t pass_type;

  union
  {
    struct
    {
      codelet_data_t *codelet;
    } direct;

    struct
    {
      uint32_t num_codelets;
      codelet_data_t *codelets[MAX_PFA_CODELETS];
    } pfa;

    struct
    {
      spv_size_t stride;
      spv_size_t row_size;
      codelet_data_t *codelet;
      spv_t twiddle;
      spv_t twiddle_inv;
    } twiddle_pre;

  } d;

} nttpass_t;

/* central repository for all NTT data that shares a
   modulus and primitive root */
typedef struct
{
  uint32_t num_codelets;
  codelet_data_t *codelets;
  spv_t codelet_const;

  nttpass_t *passes;
} nttdata_t;

/* guides for constructing transforms */
typedef struct
{
  uint32_t codelet_size;
  pass_type_t pass_type;
} nttplan_t;

void
nttdata_init_generic(const nttconfig_t *c,
            spv_t out, sp_t p, sp_t d, 
            sp_t primroot, sp_t order);

extern const nttconfig_t ntt3_config;
extern const nttconfig_t ntt4_config;
extern const nttconfig_t ntt5_config;
extern const nttconfig_t ntt7_config;
extern const nttconfig_t ntt8_config;
extern const nttconfig_t ntt9_config;
extern const nttconfig_t ntt15_config;
extern const nttconfig_t ntt16_config;
extern const nttconfig_t ntt35_config;
extern const nttconfig_t ntt40_config;

/* external interface */

void * ntt_init(sp_t size, sp_t primroot, sp_t p, sp_t d);
void ntt_free(void *data);

uint32_t planner_init(spm_t spm, sp_t size, spm_t existing);
void planner_free(nttpass_t *passes, uint32_t num_passes);

#endif /* _NTT_IMPL_H */
