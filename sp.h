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

/**************
 * GMP_IMPL.H *
 **************/

#ifdef WANT_ASSERT
#include <assert.h>
#define ASSERT(expr)   assert (expr)
#else
#define ASSERT(expr)   do {} while (0)
#endif

/* the following was inspired by longlong.h and gmp-impl.h;
 * note that a small prime must be the size of a GMP limb */
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

#ifndef ULONG_MAX
#define ULONG_MAX __GMP_ULONG_MAX
#endif

#define LONGLONG_STANDALONE
#include "longlong.h"

/* we use the remainder tree for products of 2^I0_THRESHOLD moduli or more,
   and the naive method for fewer moduli. We must have I0_THRESHOLD >= 1. */
#define I0_THRESHOLD 7

/*********
 * TYPES *
 *********/

/* SP */

/* the type for both a small prime, and a residue modulo a small prime.
 * Small primes must be 1 bit smaller than the word size for 32-bit
 * systems (otherwise there may not be enough suitable primes), but 
 * may be 2+ bits smaller when the word size exceeds 32 bits (and this
 * simplifies modular reductions)
 *
 * For a residue x modulo a sp p, we require 0 <= x < p */
typedef UWtype sp_t;

#if W_TYPE_SIZE <= 32
#define SP_NUMB_BITS (W_TYPE_SIZE - 1)
#else
#define SP_NUMB_BITS (W_TYPE_SIZE - 2)
#endif

#define SP_MIN ((sp_t)1 << (SP_NUMB_BITS - 1))
#define SP_MAX ((sp_t)(-1) >> (W_TYPE_SIZE - SP_NUMB_BITS))

/* vector of residues modulo a common small prime */
typedef sp_t * spv_t;

/* length of a spv */
typedef unsigned long spv_size_t;

typedef struct
{
  spv_t ntt_roots;
  spv_size_t twiddle_size;
  spv_t twiddle;
} __sp_nttdata;

typedef __sp_nttdata sp_nttdata_t[1];

#define MAX_NTT_BLOCK_SIZE 128

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
   (e.g. 64 bit Windows with Visual C), unsigned long has 32 bits while
   sp_t should use 64 */

static inline void 
mpz_set_sp (mpz_t m, const sp_t n)
{
  /* Is sizeof() a safe way of determining whether the conversion 
     is lossless? */
  if (sizeof (sp_t) <= sizeof (unsigned long))
    {
      mpz_set_ui (m, (unsigned long) n);
    }
  else if (sizeof (sp_t) == 8 && sizeof (unsigned long) == 4)
    {
      /* We want to right-shift by 32 bits on a 64 bit system here.
	 Putting a shift amount of 32 as a constant causes a compiler 
	 warning on 32 bit systems. So we put sizeof (sp_t) * 4
	 which always evaluates to 32 in this branch of the code, and 
	 does not cause a compiler warning if sp_t is only 4 bytes wide. */
      mpz_set_ui (m, (unsigned long) (n >> (sizeof (sp_t) * 4)));
      mpz_mul_2exp (m, m, 32UL);
      mpz_add_ui (m, m, (unsigned long int) (n & 4294967295UL));
    }
  else
    {
      abort ();
    }
}

static inline sp_t 
mpz_get_sp (const mpz_t n)
{
  if (sizeof (sp_t) == sizeof (unsigned long))
    {
      return (sp_t) mpz_get_ui (n);
    }
  else if (sizeof (sp_t) == sizeof (mp_limb_t))
    {
      /* mpz_get_ui() returns the least significant bits of the absolute 
         value of its argument that fit in an unsigned long.
         In the current GMP implementation with sign/magnitude 
	 representation, mpz_getlimbn() also returns the least sigificant
         bits of the absolute value. To allow for a future change to 
	 2's-complement representation in GMP, we should explicitly
         use mpz_abs() to a temp var here. */
      return (sp_t) mpz_getlimbn (n, 0);
    }
  else
    {
      abort ();
    }
}


void * sp_aligned_malloc (size_t len);
void sp_aligned_free (void *newptr);

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
    (defined(__x86_64__) || defined(__i386__))
  sp_t t = 0, tr = a;

  __asm__ (
    "sub %2, %0   # sp_sub: tr -= b\n\t"
    "cmovc %3, %1 # sp_sub: if (a < b) t = m\n\t"
    : "+&r" (tr), "+r" (t)
    : "g" (b), "g" (m)
    : "cc"
  );

  return tr + t;
#elif defined(_MSC_VER) && !defined(_WIN64)
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
    (defined(__x86_64__) || defined(__i386__))
  sp_t t = a - m, tr = a + b;

  __asm__ (
    "add %2, %1    # sp_add: t += b\n\t"
    "cmovc %1, %0  # sp_add: if (cy) tr = t \n\t"
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
#elif SP_NUMB_BITS <= W_TYPE_SIZE - 1
  sp_t t = a + b;
  if (t >= m)
    t -= m;
  return t;
#else
  return sp_sub(a, m - b, m);
#endif
}

/* functions used for modular reduction */

#if SP_NUMB_BITS <= W_TYPE_SIZE - 2

	/* having a small modulus allows the reciprocal
	 * to be one bit larger, which guarantees that the
	 * initial remainder fits in a word and also that at
	 * most one correction is necessary */

#define sp_reciprocal(invxl,xl)              \
  do {                                       \
    mp_limb_t dummy;                         \
    udiv_qrnnd (invxl, dummy,                \
		(sp_t) 1 << (2 * SP_NUMB_BITS + 1 -	\
		W_TYPE_SIZE), 0, xl);        \
  } while (0)

static inline sp_t sp_udiv_rem(sp_t nh, sp_t nl, sp_t d, sp_t di)
{
  sp_t r;
  mp_limb_t q1, q2, tmp;
  q1 = nh << (2*(W_TYPE_SIZE - SP_NUMB_BITS)) |
	    nl >> (2*SP_NUMB_BITS - W_TYPE_SIZE);
  umul_ppmm (q2, tmp, q1, di);
  r = nl - d * (q2 >> 1);
  return sp_sub(r, d, d);
}

#else    /* big modulus; no shortcuts allowed */

#define sp_reciprocal(invxl,xl)              \
  do {                                       \
    mp_limb_t dummy;                         \
    udiv_qrnnd (invxl, dummy,                \
		(sp_t) 1 << (2 * SP_NUMB_BITS -	\
		W_TYPE_SIZE), 0, xl);        \
  } while (0)

static inline sp_t sp_udiv_rem(sp_t nh, sp_t nl, sp_t d, sp_t di)
{
  mp_limb_t q1, q2, tmp, dqh, dql;
  q1 = nh << (2*(W_TYPE_SIZE - SP_NUMB_BITS)) |
	    nl >> (2*SP_NUMB_BITS - W_TYPE_SIZE);
  umul_ppmm (q2, tmp, q1, di);
  umul_ppmm (dqh, dql, q2, d);

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
  umul_ppmm (u, v, x, y);
  return sp_udiv_rem (u, v, m, d);
}

/* x*y mod m */
static inline sp_t
sp_sqr (sp_t x, sp_t m, sp_t d)
{
  sp_t u, v;
  umul_ppmm (u, v, x, x);
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

/* x / 2 mod m */
#define sp_div_2(x,m) (((x) & 1) ? (m) - (((m) - (x)) >> 1) : ((x) >> 1))
  
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
void mpzspv_random (mpzspv_t, spv_size_t, spv_size_t, mpzspm_t);

#endif /* _SP_H */
