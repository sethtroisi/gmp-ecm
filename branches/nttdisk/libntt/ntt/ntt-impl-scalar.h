#ifndef _NTT_IMPL_SCALAR_H
#define _NTT_IMPL_SCALAR_H

#include "ntt-impl.h"

#ifdef __cplusplus
extern "C" {
#endif

static INLINE sp_t sp_ntt_add(sp_t a, sp_t b, sp_t p)
{
#ifdef HAVE_PARTIAL_MOD
  return sp_add(a, b, 2 * p);
#else
  return sp_add(a, b, p);
#endif
}

static INLINE sp_t sp_ntt_add_partial(sp_t a, sp_t b, sp_t p)
{
#ifdef HAVE_PARTIAL_MOD
  return a + b;
#else
  return sp_ntt_add(a, b, p);
#endif
}

static INLINE sp_t sp_ntt_sub(sp_t a, sp_t b, sp_t p)
{
#ifdef HAVE_PARTIAL_MOD
  return sp_sub(a, b, 2 * p);
#else
  return sp_sub(a, b, p);
#endif
}

static INLINE sp_t sp_ntt_sub_partial(sp_t a, sp_t b, sp_t p)
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

#if SP_NUMB_BITS == 50 /* floating point operands */

static INLINE sp_t sp_ntt_mul(sp_t x, sp_t w, sp_t whi, sp_t p)
{
  sp_t recip = 1.0 / p;  /* common subexpression */

  sp_t hi, lo;
  {
    sp_t ahi, alo, bhi, blo;
    sp_split(ahi, alo, x);
    bhi = whi;
    blo = w - whi;   /* save the overhead of splitting w */
    hi = x * w;
    lo = ((ahi * bhi - hi) + ahi * blo + alo * bhi) + alo * blo; \
  }
  return sp_udiv_rem(hi, lo, p, recip);
}

#else /* integer operands */

/* low half multiply */

static INLINE sp_t sp_mul_lo(sp_t a, sp_t b)
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

static INLINE sp_t sp_mul_hi(sp_t a, sp_t b)
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

static INLINE sp_t sp_ntt_mul(sp_t x, sp_t w, sp_t w_inv, sp_t p)
{
  sp_t q = sp_mul_hi(x, w_inv);
  sp_t r = x * w - q * p;

#ifdef HAVE_PARTIAL_MOD
  return r;
#else
  return sp_sub(r, p, p);
#endif
}

#endif

extern const nttconfig_t X(ntt2_config);
extern const nttconfig_t X(ntt3_config);
extern const nttconfig_t X(ntt4_config);
extern const nttconfig_t X(ntt5_config);
extern const nttconfig_t X(ntt7_config);
extern const nttconfig_t X(ntt8_config);
extern const nttconfig_t X(ntt9_config);
extern const nttconfig_t X(ntt15_config);
extern const nttconfig_t X(ntt16_config);
extern const nttconfig_t X(ntt35_config);
extern const nttconfig_t X(ntt40_config);

#ifdef __cplusplus
}
#endif

#endif /* _NTT_IMPL_SCALAR_H */
