/* sp.h - header file for the sp library

  Copyright 2005, 2008, 2010 Dave Newman, Jason Papadopoulos and
                             Paul Zimmermann.
  Copyright 1991, 1993, 1994, 1995, 1996, 1997, 1999, 2000, 2001, 2002, 2003,
  2004, 2005, 2010 Free Software Foundation, Inc. (for parts from gmp-impl.h).

  This file is part of the SP library.
  
  The SP Library is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or (at your
  option) any later version.

  The SP Library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
  License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with the SP Library; see the file COPYING.LIB.  If not, write to
  the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
  MA 02110-1301, USA.
*/

#ifndef _SP_H
#define _SP_H

#include "config.h"
#include <stdlib.h>

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h> /* needed for size_t */
#endif

/* name some types whose bit widths are unambiguous */

#if defined(HAVE_STDINT_H)
#include <stdint.h>
#elif defined(HAVE_INTTYPES_H)
#include <inttypes.h>
#elif defined(_MSC_VER)
typedef unsigned __int32 uint32_t;
typedef unsigned __int64 uint64_t;
#endif

/* include (hopefully standard) SSE2 interface, along
   with mnemonics for the small part of the ISA that 
   we need */

#ifdef HAVE_SSE2
#include <emmintrin.h>

#define pload  _mm_load_si128
#define pstore _mm_store_si128
#define pand _mm_and_si128
#define psetzero() _mm_setzero_si128()
#define pcvt_i32 _mm_cvtsi32_si128
#define pcvt_i64 _mm_cvtsi64_si128
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

#endif

#ifndef TUNE
#include "ecm-params.h"
#else
extern size_t NTT_GFP_TWIDDLE_DIF_BREAKOVER;
extern size_t NTT_GFP_TWIDDLE_DIT_BREAKOVER;
extern size_t MUL_NTT_THRESHOLD;
extern size_t PREREVERTDIVISION_NTT_THRESHOLD;
extern size_t POLYINVERT_NTT_THRESHOLD;
extern size_t POLYEVALT_NTT_THRESHOLD;
extern size_t MPZSPV_NORMALISE_STRIDE;
#endif

#include <gmp.h>

#ifdef WANT_ASSERT
#include <assert.h>
#define ASSERT(expr)   assert (expr)
#else
#define ASSERT(expr)   do {} while (0)
#endif

/* Basic defs for the data types used in Number Theoretic Transforms

   The choice of basic data type has to balance computational
   efficiency with the ability to find enough small prime moduli
   to be useful. We use the GMP concept of 'nails' when choosing
   the number of bits in a prime modulus, since a modular reduction
   using B-bit arithmetic runs much faster when the modulus has
   B-2 bits or less.
   
   64-bit machines with a 64-bit OS should use 62-bit primes and 
   64-bit operations, except that Microsoft's compiler defaults 
   to 32 bits for an unsigned long even on 64-bit versions of 
   Windows, and the GMP limb size is 64 bits just to make things 
   more confusing. 
   
   For 32-bit machines the most efficient arithmetic derives from 
   using 30-bit moduli, but for large transform sizes we will run 
   out of moduli too soon. Using 31-bit moduli delays the inevitable 
   for long enough to be useful for larger problems, at a cost of 
   ~10% in slower modular reductions, but even then we will still
   run out of primes for transform lengths around 2^20, i.e. not 
   very much. So to deal with the largest problems we will need to 
   support 62-bit moduli on a 32-bit machine.

   A completely different possibility for PowerPC machines is to
   use 50-bit moduli and use the FPU for modular reductions. This 
   works on that architecture because the PowerPC standard includes
   fused multiply-add (FMA) instructions that perform a double-precision
   multiply and add with a single rounding operation, and the multiply
   result is required to have 106 bits of precision before its bottom
   half is discarded. The upshot is that we can compute a 101-bit
   signed product of two 51-bit signed inputs in just two FPU 
   instructions, which modern-day PowerPC processors can execute with
   very high throughput (likely much higher than using integer 
   multiplies).

   The defines below choose 30, 31, or 32-63 bit prime moduli, a choice 
   that is decoupled from a 32-bit or 64-bit machine word size. We 
   still include a bunch of GMP plumbing to take advantage of custom
   macros for arithmetic at the native machine word size */

/* types that longlong.h needs */
typedef mp_limb_t UWtype;
typedef unsigned int UHWtype;
#if (defined(_PA_RISC1_1) && defined(__GNUC__))
/* this seems to be needed, otherwise umul_ppmm() does not work properly */
typedef mp_limb_t USItype __attribute__ ((mode (SI)));
typedef mp_limb_t UDItype __attribute__ ((mode (DI)));
#else
typedef mp_limb_t USItype;
typedef mp_limb_t UDItype;
#endif

#ifndef W_TYPE_SIZE
#define W_TYPE_SIZE GMP_LIMB_BITS
#endif

#define LONGLONG_STANDALONE
#include "longlong.h"

/*********
 * TYPES *
 *********/

/* SP */

/* the type for both a small prime modulus, and a residue
 * modulo a small prime. The default type is the largest
 * supported by GMP arithmetic; redefine SP_NUMB_BITS
 * to override
 *
 * For a residue x modulo a sp p, we require 0 <= x < p */

#ifndef SP_NUMB_BITS
   #if GMP_LIMB_BITS == 64
      #define SP_NUMB_BITS 62
   #elif GMP_LIMB_BITS == 32
      #define SP_NUMB_BITS 31
   #endif
#endif

#if SP_NUMB_BITS < 30 || SP_NUMB_BITS > 63
#error "invalid choice of small prime size"
#endif

#if SP_NUMB_BITS < 32
typedef uint32_t sp_t;
#define SP_TYPE_BITS 32
#else
typedef uint64_t sp_t;
#define SP_TYPE_BITS 64
#endif

#define SP_MIN ((sp_t)1 << (SP_NUMB_BITS - 1))
#define SP_MAX ((sp_t)(-1) >> (SP_TYPE_BITS - SP_NUMB_BITS))

/* vector of residues modulo a common small prime */
typedef sp_t * spv_t;

/* length of a spv, and the type used for array offsets */
typedef unsigned long spv_size_t;

typedef struct
{
  spv_t ntt_roots;
  spv_size_t twiddle_size;
  spv_t twiddle;
} __sp_nttdata;

typedef __sp_nttdata sp_nttdata_t[1];

#define MAX_NTT_BLOCK_SIZE 128

/* Which steps to perform in convolution product funtions:
   forward transform, pair-wise multiplication, inverse transform */
#define NTT_MUL_STEP_FFT1 1
#define NTT_MUL_STEP_FFT2 2
#define NTT_MUL_STEP_MUL 4
#define NTT_MUL_STEP_IFFT 8

/* SPM */

/* small prime modulus - this contains some precomputed constants to
 * calculate modulo a sp */
typedef struct
{
  sp_t sp;		/* value of the sp */
  sp_t mul_c;		/* constant used for reduction mod sp */
  sp_t prim_root;
  sp_t inv_prim_root;
  sp_nttdata_t nttdata;
  sp_nttdata_t inttdata;
  spv_t scratch;
#if SP_TYPE_BITS > GMP_LIMB_BITS
  mpz_t mp_sp;
#endif
} __spm_struct;

typedef __spm_struct * spm_t;

/* MPZSPM */

typedef mpz_t * mpzv_t;

typedef struct
  {
    /* number of small primes needed to represent each coeff */
    unsigned int sp_num;
    spv_size_t max_ntt_size;
    
    mpz_t modulus;
    
    /* spm data */
    spm_t *spm;
    
    /* precomputed crt constants, see mpzspm.c */
    mpzv_t crt1, crt2;
    sp_t *crt3, **crt4, *crt5;

    /* product tree to speed up conversion from mpz to sp */
    mpzv_t *T;            /* product tree */
    unsigned int d;       /* ceil(log(sp_num)/log(2)) */
  } __mpzspm_struct;

typedef __mpzspm_struct * mpzspm_t;

/* MPZSPV */

/* sp representation of a mpz polynomial */

typedef spv_t * mpzspv_t;

#define MAX(x,y) (((x)<(y))?(y):(x))
#define MIN(x,y) (((x)<(y))?(x):(y))

#define SIZ(x) ((x)->_mp_size)
#define PTR(x) ((x)->_mp_d)

/* expanding macros and then turning them 
   into strings requires two levels of macro-izing */

#define _(x) #x
#define STRING(x) _(x)

/*************
 * FUNCTIONS *
 *************/

/* general */

static inline unsigned int
ceil_log_2 (spv_size_t x)
{
  unsigned int a = 0;
  
  x--;
  while (x)
    {
      a++;
      x >>= 1;
    }
  return a;
}

/* Conversion functions sp_t <-> mpz_t. Using mpz_*_ui() functions is not
   portable as those take unsigned long's, but on some systems 
   (e.g. 64 bit Windows with Visual C), unsigned long has 32 bits */

static inline void 
mpz_set_sp (mpz_t m, const sp_t n)
{
#if SP_TYPE_BITS == 32

  mpz_set_ui(m, (unsigned long)(uint32_t)n);

#else /* 64-bit sp_t */

  #if GMP_LIMB_BITS == 64  /* 64-bit GMP limb, 64-bit sp_t */
  if (sizeof(sp_t) > sizeof(unsigned long))
    {
       mpz_set_ui (m, (unsigned long)(uint32_t)(n >> 32));
       mpz_mul_2exp (m, m, 32);
       mpz_add_ui (m, m, (unsigned long)(uint32_t)n);
    }
  else
    mpz_set_ui(m, (unsigned long)n);

  #else                    /* 32-bit GMP limb, 64-bit sp_t */
  mpz_set_ui (m, (unsigned long)(uint32_t)(n >> 32));
  mpz_mul_2exp (m, m, 32);
  mpz_add_ui (m, m, (unsigned long)(uint32_t)n);
  #endif

#endif
}

static inline sp_t 
mpz_get_sp (const mpz_t n)
{
#if SP_TYPE_BITS == 32

  return mpz_get_ui(n);

#else /* 64-bit sp_t */

  sp_t m = mpz_getlimbn(n, 0);

  #if GMP_LIMB_BITS == 32  /* 32-bit GMP limb, 64-bit sp_t */
  if (mpz_size(n) >= 2)
    m |= (sp_t)mpz_getlimbn(n, 1) << 32;
  #endif

  return m;
#endif
}


void * sp_aligned_malloc (size_t len);
void sp_aligned_free (void *newptr);
sp_t sp_reciprocal(sp_t p);

/* sp */

/* Routines for arithmetic on residues modulo a small prime
 *
 * All functions return values in the range 0 <= x < p.
 *
 * The variable name of the modulus is 'p' if the input must be prime,
 *                                     'm' if we also allow composites. */


static inline sp_t sp_sub(sp_t a, sp_t b, sp_t m) 
{
#if (defined(__GNUC__) || defined(__ICL)) && \
    (defined(__x86_64__) || defined(__i386__)) && \
    (SP_TYPE_BITS <= GMP_LIMB_BITS)
  sp_t t = 0, tr = a;

  __asm__ (
    "sub %2, %0   # sp_sub: tr -= b\n\t"
    "cmovc %3, %1 # sp_sub: if (a < b) t = m\n\t"
    : "+&r" (tr), "+r" (t)
    : "g" (b), "g" (m)
    : "cc"
  );

  return tr + t;
#elif defined(_MSC_VER) && !defined(_WIN64) && (SP_TYPE_BITS == 32)

  __asm
    {
        mov     eax, a
        xor     edx, edx
        sub     eax, b
        cmovb   edx, m
        add     eax, edx
    }
#else
  if (a >= b)
    return a - b;
  else
    return a - b + m;
#endif
}

static inline sp_t sp_add(sp_t a, sp_t b, sp_t m) 
{
#if (defined(__GNUC__) || defined(__ICL)) && \
    (defined(__x86_64__) || defined(__i386__)) && \
    (SP_TYPE_BITS <= GMP_LIMB_BITS)
  sp_t t = a - m, tr = a + b;

  __asm__ (
    "add %2, %1    # sp_add: t += b\n\t"
    "cmovc %1, %0  # sp_add: if (cy) tr = t \n\t"
    : "+r" (tr), "+&r" (t)
    : "g" (b)
    : "cc"
  );

  return tr;
#elif defined(_MSC_VER) && !defined(_WIN64) && (SP_TYPE_BITS == 32)
  __asm
    {
        mov     eax, a
        add     eax, b
        mov     edx, eax
        sub     edx, m
        cmovnc  eax, edx
    }
#else
  sp_t t = a + b;
  if (t >= m)
    t -= m;
  return t;
#endif
}

/* widening multiply */

#if SP_TYPE_BITS == GMP_LIMB_BITS  /* use GMP functions */

  #define sp_wide_mul(hi, lo, a, b) umul_ppmm(hi, lo, a, b)
  #define sp_wide_sqr(hi, lo, a) sp_wide_mul(hi, lo, a, a)

#elif SP_TYPE_BITS < GMP_LIMB_BITS  /* ordinary multiply */

  #define sp_wide_mul(hi, lo, a, b)		    \
  	{					    \
	  mp_limb_t prod = (mp_limb_t)(a) *	    \
	  		   (mp_limb_t)(b);	    \
	  hi = (sp_t)(prod >> (GMP_LIMB_BITS / 2)); \
	  lo = (sp_t)prod;			    \
	}

  #define sp_wide_sqr(hi, lo, a) sp_wide_mul(hi, lo, a, a)

#else  /* worst case: build product from smaller multiplies */

  #define sp_wide_mul(hi, lo, a, b)		    		\
  	{					    		\
	  mp_limb_t a1 = (mp_limb_t)((a) >> GMP_LIMB_BITS);	\
	  mp_limb_t a0 = (mp_limb_t)(a);			\
	  mp_limb_t b1 = (mp_limb_t)((b) >> GMP_LIMB_BITS);	\
	  mp_limb_t b0 = (mp_limb_t)(b);			\
	  mp_limb_t a0b0_hi, a0b0_lo;				\
	  mp_limb_t a0b1_hi, a0b1_lo;				\
	  mp_limb_t a1b0_hi, a1b0_lo;				\
	  mp_limb_t a1b1_hi, a1b1_lo;				\
								\
	  umul_ppmm(a0b0_hi, a0b0_lo, a0, b0);			\
	  umul_ppmm(a0b1_hi, a0b1_lo, a0, b1);			\
	  umul_ppmm(a1b0_hi, a1b0_lo, a1, b0);			\
	  umul_ppmm(a1b1_hi, a1b1_lo, a1, b1);			\
								\
	  /* proceed a word at a time; the high word		\
	     of each product except a0b0 have at least		\
	     two bits clear, so carries propagate one		\
	     word at most */					\
								\
	  add_ssaaaa(a0b1_hi, a0b1_lo, 				\
	             a0b1_hi, a0b1_lo,				\
		     (mp_limb_t)0, a0b0_hi);			\
	  add_ssaaaa(a0b1_hi, a0b1_lo, 				\
	             a0b1_hi, a0b1_lo,				\
		     a1b0_hi, a1b0_lo);				\
	  add_ssaaaa(a1b1_hi, a1b1_lo, 				\
	             a1b1_hi, a1b1_lo,				\
		     (mp_limb_t)0, a0b1_hi);			\
								\
	  hi = (sp_t)a1b1_hi << GMP_LIMB_BITS | a1b1_lo;	\
	  lo = (sp_t)a0b1_lo << GMP_LIMB_BITS | a0b0_lo;	\
	}

  #define sp_wide_sqr(hi, lo, a)		    		\
  	{					    		\
	  mp_limb_t a1 = (mp_limb_t)((a) >> GMP_LIMB_BITS);	\
	  mp_limb_t a0 = (mp_limb_t)(a);			\
	  mp_limb_t a0a0_hi, a0a0_lo;				\
	  mp_limb_t a0a1_hi, a0a1_lo;				\
	  mp_limb_t a1a1_hi, a1a1_lo;				\
								\
	  umul_ppmm(a0a0_hi, a0a0_lo, a0, a0);			\
	  umul_ppmm(a0a1_hi, a0a1_lo, a0, a1);			\
	  umul_ppmm(a1a1_hi, a1a1_lo, a1, a1);			\
								\
	  add_ssaaaa(a0a1_hi, a0a1_lo, 				\
	             a0a1_hi, a0a1_lo,				\
		     a0a1_hi, a0a1_lo);				\
	  add_ssaaaa(a0a1_hi, a0a1_lo, 				\
	             a0a1_hi, a0a1_lo,				\
		     (mp_limb_t)0, a0a0_hi);			\
	  add_ssaaaa(a1a1_hi, a1a1_lo, 				\
	             a1a1_hi, a1a1_lo,				\
		     (mp_limb_t)0, a0a1_hi);			\
								\
	  hi = (sp_t)a1a1_hi << GMP_LIMB_BITS | a1a1_lo;	\
	  lo = (sp_t)a0a1_lo << GMP_LIMB_BITS | a0a0_lo;	\
	}

#endif

/* modular reduction */

#if SP_NUMB_BITS <= SP_TYPE_BITS - 2

    /* having a small modulus allows the reciprocal
     * to be one bit larger, which guarantees that the
     * initial remainder fits in a word and also that at
     * most one correction is necessary */

static inline sp_t sp_udiv_rem(sp_t nh, sp_t nl, sp_t d, sp_t di)
{
  sp_t r, q1, q2, tmp;

  q1 = nh << (2*(SP_TYPE_BITS - SP_NUMB_BITS)) |
	    nl >> (2*SP_NUMB_BITS - SP_TYPE_BITS);
  sp_wide_mul (q2, tmp, q1, di);
  r = nl - d * (q2 >> 1);
  return sp_sub(r, d, d);
}

#else    /* big modulus; no shortcuts allowed */

static inline sp_t sp_udiv_rem(sp_t nh, sp_t nl, sp_t d, sp_t di)
{
  sp_t q1, q2, tmp, dqh, dql;
  q1 = nh << (2*(SP_TYPE_BITS - SP_NUMB_BITS)) |
	    nl >> (2*SP_NUMB_BITS - SP_TYPE_BITS);
  sp_wide_mul (q2, tmp, q1, di);
  sp_wide_mul (dqh, dql, q2, d);

  tmp = nl;
  nl = tmp - dql;
  nh = nh - dqh - (nl > tmp);
  if (nh)
	  nl -= d;
  nl = sp_sub(nl, d, d);
  return sp_sub(nl, d, d);
}

#endif

/* x*y mod m */
static inline sp_t
sp_mul (sp_t x, sp_t y, sp_t m, sp_t d)
{
  sp_t u, v;
  sp_wide_mul (u, v, x, y);
  return sp_udiv_rem (u, v, m, d);
}

/* x*y mod m */
static inline sp_t
sp_sqr (sp_t x, sp_t m, sp_t d)
{
  sp_t u, v;
  sp_wide_sqr (u, v, x);
  return sp_udiv_rem (u, v, m, d);
}

#define sp_neg(x,m) ((x) == (sp_t) 0 ? (sp_t) 0 : (m) - (x))

/* Returns x^a % m, uses a right-to-left powering ladder */

static inline sp_t
sp_pow (sp_t x, sp_t a, sp_t m, sp_t d)
{
  sp_t partial = 1;

  while (1)
    {
      if (a & 1)
	partial = sp_mul (x, partial, m, d);

      a >>= 1;

      if (!a)
	return partial;

      x = sp_sqr (x, m, d);
    }
}

/* 1/x mod p */
#define sp_inv(x,p,d) sp_pow (x, (p) - 2, p, d)

int sp_spp (sp_t, sp_t, sp_t);
int sp_prime (sp_t);

/* spm */

spm_t spm_init (spv_size_t, sp_t);
void spm_clear (spm_t);

/* spv */

/* ASSIGNMENT */

void spv_set (spv_t, spv_t, spv_size_t);
void spv_rev (spv_t, spv_t, spv_size_t);
void spv_set_sp (spv_t, sp_t, spv_size_t);
void spv_set_zero (spv_t, spv_size_t);

/* ARITHMETIC */

/* add */
void spv_add (spv_t, spv_t, spv_t, spv_size_t, sp_t);
void spv_add_sp (spv_t, spv_t, sp_t, spv_size_t, sp_t);

/* subtract */
void spv_sub (spv_t, spv_t, spv_t, spv_size_t, sp_t);
void spv_sub_sp (spv_t, spv_t, sp_t, spv_size_t, sp_t);
void spv_neg (spv_t, spv_t, spv_size_t, sp_t);

/* pointwise multiplication */
void spv_pwmul (spv_t, spv_t, spv_t, spv_size_t, sp_t, sp_t);
void spv_pwmul_rev (spv_t, spv_t, spv_t, spv_size_t, sp_t, sp_t);
void spv_mul_sp (spv_t, spv_t, sp_t, spv_size_t, sp_t, sp_t);

void spv_random (spv_t, spv_size_t, sp_t);
int spv_cmp (spv_t, spv_t, spv_size_t);

/* ntt_gfp */

void spv_ntt_gfp_dif (spv_t, spv_size_t, spm_t);
void spv_ntt_gfp_dit (spv_t, spv_size_t, spm_t);

/* mpzspm */

spv_size_t mpzspm_max_len (mpz_t);
mpzspm_t mpzspm_init (spv_size_t, mpz_t);
void mpzspm_clear (mpzspm_t);

/* mpzspv */

/* we use the remainder tree for products of 2^I0_THRESHOLD moduli or more,
   and the naive method for fewer moduli. We must have I0_THRESHOLD >= 1. */
#define I0_THRESHOLD 7

mpzspv_t mpzspv_init (spv_size_t, mpzspm_t);
void mpzspv_clear (mpzspv_t, mpzspm_t);
int mpzspv_verify (mpzspv_t, spv_size_t, spv_size_t, mpzspm_t);
void mpzspv_set (mpzspv_t, spv_size_t, mpzspv_t, spv_size_t, spv_size_t,
    mpzspm_t);
void mpzspv_revcopy (mpzspv_t, spv_size_t, mpzspv_t, spv_size_t, spv_size_t,
    mpzspm_t);
void mpzspv_set_sp (mpzspv_t, spv_size_t, sp_t, spv_size_t, mpzspm_t);
void mpzspv_from_mpzv (mpzspv_t, const spv_size_t, const mpzv_t, 
		       const spv_size_t, mpzspm_t);
void mpzspv_reverse (mpzspv_t, spv_size_t, spv_size_t, mpzspm_t);
void mpzspv_neg (mpzspv_t, spv_size_t, mpzspv_t, spv_size_t, spv_size_t,
    mpzspm_t);
void mpzspv_add (mpzspv_t, spv_size_t, mpzspv_t, spv_size_t, mpzspv_t,
    spv_size_t, spv_size_t, mpzspm_t);
void mpzspv_to_mpzv (mpzspv_t, spv_size_t, mpzv_t, spv_size_t, mpzspm_t);
void mpzspv_normalise (mpzspv_t, spv_size_t, spv_size_t, mpzspm_t);
void mpzspv_pwmul (mpzspv_t, spv_size_t, mpzspv_t, spv_size_t, mpzspv_t, 
    spv_size_t, spv_size_t, mpzspm_t);
void mpzspv_to_ntt (mpzspv_t, spv_size_t, spv_size_t, spv_size_t, int,
    mpzspm_t);
void mpzspv_from_ntt (mpzspv_t, spv_size_t, spv_size_t, spv_size_t, mpzspm_t);
void mpzspv_mul_ntt (mpzspv_t, spv_size_t, mpzspv_t, spv_size_t, spv_size_t, 
    mpzspv_t, spv_size_t, spv_size_t, spv_size_t, int, spv_size_t, mpzspm_t, 
    int);
void mpzspv_random (mpzspv_t, spv_size_t, spv_size_t, mpzspm_t);
void mpzspv_to_dct1 (mpzspv_t, mpzspv_t, spv_size_t, spv_size_t, mpzspv_t, 
    mpzspm_t);
void mpzspv_mul_by_dct (mpzspv_t, const mpzspv_t, spv_size_t, const mpzspm_t, 
    int);
void mpzspv_sqr_reciprocal (mpzspv_t, spv_size_t, const mpzspm_t);

#endif /* _SP_H */
