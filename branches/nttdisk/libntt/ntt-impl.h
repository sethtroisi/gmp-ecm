#ifndef _NTT_IMPL_H
#define _NTT_IMPL_H

#include "sp.h"

typedef uint32_t (*get_num_ntt_const_t)(void);

typedef void (*nttdata_init_t)(spv_t out, 
				sp_t p, sp_t d,
				sp_t primroot, sp_t order);

typedef void (*ntt_run_t)(spv_t x, spv_size_t stride,
			  sp_t p, sp_t d, spv_t ntt_const);

typedef void (*ntt_pfa_run_t)(spv_t x, spv_size_t stride,
			  spv_size_t cofactor, 
			  sp_t p, sp_t d, spv_t ntt_const);


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

/*------------------- definitions for SIMD transforms ----------------*/

#ifdef HAVE_SSE2

typedef __m128i sp_simd_t;

#define SP_SIMD_VSIZE (128 / SP_TYPE_BITS)

static inline sp_simd_t sp_simd_gather(spv_t x, spv_size_t start_off, 
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

static inline void sp_simd_scatter(sp_simd_t t, spv_t x, 
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

static inline sp_simd_t sp_simd_add(sp_simd_t a, sp_simd_t b, sp_t p)
{
#if SP_TYPE_BITS == 32

  sp_simd_t vp = pshufd(pcvt_i32(p), 0x00);
  sp_simd_t t0, t1;

  t0 = paddd(a, b);
  t0 = psubd(t0, vp);
  t1 = pcmpgtd(psetzero(), t0);
  t1 = pand(t1, vp);
  return paddd(t0, t1);

#else

  sp_simd_t vp = pshufd(pcvt_i64(p), 0x44);
  sp_simd_t t0, t1;

  t0 = paddq(a, b);
  t0 = psubq(t0, vp);
  t1 = pcmpgtd(psetzero(), t0);
  t1 = pshufd(t1, 0xf5);
  t1 = pand(t1, vp);
  return paddq(t0, t1);

#endif
}

static inline sp_simd_t sp_simd_sub(sp_simd_t a, sp_simd_t b, sp_t p)
{
#if SP_TYPE_BITS == 32

  sp_simd_t vp = pshufd(pcvt_i32(p), 0x00);
  sp_simd_t t0, t1;

  t0 = psubd(a, b);
  t1 = pcmpgtd(psetzero(), t0);
  t1 = pand(t1, vp);
  return paddd(t0, t1);

#else

  sp_simd_t vp = pshufd(pcvt_i64(p), 0x44);
  sp_simd_t t0, t1;

  t0 = psubq(a, b);
  t1 = pcmpgtd(psetzero(), t0);
  t1 = pshufd(t1, 0xf5);
  t1 = pand(t1, vp);
  return paddq(t0, t1);

#endif
}


static inline sp_simd_t sp_simd_mul(sp_simd_t a, sp_t b, sp_t p, sp_t d)
{
#if SP_TYPE_BITS == 32

  __m128i t0, t1, t2, t3, vp, vp2, vd;

  vp = pshufd(pcvt_i32(p), 0x00);
  vp2 = pshufd(pcvt_i32(p), 0x44);
  vd = pshufd(pcvt_i32(d), 0x00);
  t1 = pshufd(pcvt_i32(b), 0x44);

  t2 = pshufd(a, 0x31);
  t0 = pmuludq(a, t1);
  t2 = pmuludq(t2, t1);
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

#elif GMP_LIMB_BITS == 32   /* 64-bit sp_t on a 32-bit machine */

  __m128i t0, t1, t2, t3, t4, t5, vp, vd;

  vp = pshufd(pcvt_i64(p), 0x44);
  vd = pshufd(pcvt_i64(d), 0x44);
  t4 = a;
  t5 = pshufd(pcvt_i64(b), 0x44);

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
  t3 = psubq(t3, vp);
  t1 = pcmpgtd(psetzero(), t3);
  t1 = pshufd(t1, 0xf5);
  t1 = pand(t1, vp);
  return paddq(t3, t1);

#else

  /* there's no way the SSE2 unit can keep up with a
     64-bit multiplier in the ALU */

  sp_simd_t t0, t1;
  sp_t a0, a1;

  pstore_i64(a0, a);
  pstore_i64(a1, pshufd(a, 0x0e));

  a0 = sp_mul(a0, b, p, d);
  a1 = sp_mul(a1, b, p, d);

  t0 = pcvt_i64(a0);
  t1 = pcvt_i64(a1);
  return punpcklo64(t0, t1);

#endif
}


#endif  /*-----------------------------------------------------------*/

typedef struct
{
  uint32_t size;
  get_num_ntt_const_t get_num_ntt_const;
  nttdata_init_t nttdata_init;
  ntt_run_t ntt_run;
  ntt_pfa_run_t ntt_pfa_run;
} nttconfig_t;

extern const nttconfig_t ntt2_config;
extern const nttconfig_t ntt3_config;
extern const nttconfig_t ntt4_config;
extern const nttconfig_t ntt5_config;
extern const nttconfig_t ntt6_config;
extern const nttconfig_t ntt7_config;
extern const nttconfig_t ntt8_config;
extern const nttconfig_t ntt9_config;
extern const nttconfig_t ntt15_config;

#endif /* _NTT_IMPL_H */
