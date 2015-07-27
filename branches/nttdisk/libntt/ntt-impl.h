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

#if SP_NUMB_BITS == SP_TYPE_BITS - 2
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
  PASS_TYPE_TWIDDLE
} pass_type_t;

#define MAX_PASSES 10

typedef struct
{
  pass_type_t pass_type;
  codelet_data_t *codelet;

  union
  {
    struct
    {
      spv_size_t num_transforms;
    } direct;

    struct
    {
      spv_size_t cofactor;
    } pfa;

    struct
    {
      spv_size_t num_transforms;
      spv_size_t stride;
      spv_t w;
    } twiddle;

  } d;

} nttpass_t;

/* central repository for all NTT data that shares a
   modulus and primitive root */
typedef struct
{
  uint32_t num_codelets;
  codelet_data_t *codelets;
  spv_t codelet_const;

  uint32_t num_passes;
  nttpass_t passes[MAX_PASSES];
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

uint32_t 
ntt_build_passes(nttdata_t *data,
    		nttplan_t *plans, uint32_t num_plans,
		sp_t size, sp_t p, sp_t primroot, sp_t order, sp_t d);

/* external interface */

void * ntt_init(sp_t size, sp_t primroot, sp_t p, sp_t d);
void ntt_free(void *data);
void ntt_reset(void *data);
void ntt_run(spv_t x, sp_t p, void *data);

uint32_t planner_init(spm_t spm, sp_t size, spm_t existing);
void planner_free(nttpass_t *passes, uint32_t num_passes);

/* SIMD includes */
#ifdef HAVE_SSE2
#include "ntt-impl-sse2.h"
#endif

#endif /* _NTT_IMPL_H */
