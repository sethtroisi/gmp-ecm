/* sp.h - header file for the sp library

  Copyright 2005 Dave Newman.

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
  the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
  MA 02111-1307, USA.
*/

#ifndef __HAVE_SP_H
#define __HAVE_SP_H

#include "config.h"

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h> /* needed for size_t */
#endif

#ifndef TUNE
#include "ecm-params.h"
#else
extern size_t SPV_NTT_GFP_DIF_RECURSIVE_THRESHOLD;
extern size_t SPV_NTT_GFP_DIT_RECURSIVE_THRESHOLD;
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

#if WANT_ASSERT
#include <assert.h>
#define ASSERT(expr)   assert (expr)
#else
#define ASSERT(expr)   do {} while (0)
#endif

typedef unsigned long UWtype;
typedef unsigned int UHWtype;
typedef unsigned long USItype;
typedef unsigned long UDItype;

#ifndef BITS_PER_MP_LIMB
#define BITS_PER_MP_LIMB __GMP_BITS_PER_MP_LIMB
#endif

#ifndef W_TYPE_SIZE
#define W_TYPE_SIZE BITS_PER_MP_LIMB
#endif

#ifndef ULONG_MAX
#define ULONG_MAX __GMP_ULONG_MAX
#endif

#include "longlong.h"

/* Use a library function for invert_limb, if available. */
#if ! defined (invert_limb) && HAVE_NATIVE_mpn_invert_limb
#define mpn_invert_limb  __MPN(invert_limb)
mp_limb_t mpn_invert_limb _PROTO ((mp_limb_t)) ATTRIBUTE_CONST;
#define invert_limb(invxl,xl)  (invxl = mpn_invert_limb (xl))
#endif

#ifndef invert_limb
#define invert_limb(invxl,xl)                   \
  do {                                          \
    mp_limb_t dummy;                            \
    ASSERT ((xl) != 0);                         \
    if (xl << 1 == 0)                           \
      invxl = ~(mp_limb_t) 0;                   \
    else                                        \
      udiv_qrnnd (invxl, dummy, -xl, 0, xl);    \
  } while (0)
#endif

/* Divide the two-limb number in (NH,,NL) by D, with DI being the largest
   limb not larger than (2**(2*BITS_PER_MP_LIMB))/D - (2**BITS_PER_MP_LIMB).
   If this would yield overflow, DI should be the largest possible number
   (i.e., only ones).  For correct operation, the most significant bit of D
   has to be set.  Put the quotient in Q and the remainder in R.  */
#define udiv_qrnnd_preinv(q, r, nh, nl, d, di)                            \
  do {                                                                    \
    mp_limb_t _q, _ql, _r;                                                \
    mp_limb_t _xh, _xl;                                                   \
    ASSERT ((d) != 0);                                                    \
    umul_ppmm (_q, _ql, (nh), (di));                                      \
    _q += (nh);                 /* DI is 2**BITS_PER_MP_LIMB too small */ \
    umul_ppmm (_xh, _xl, _q, (d));                                        \
    sub_ddmmss (_xh, _r, (nh), (nl), _xh, _xl);                           \
    if (_xh != 0)                                                         \
      {                                                                   \
	sub_ddmmss (_xh, _r, _xh, _r, 0, (d));                            \
	_q += 1;                                                          \
	if (_xh != 0)                                                     \
	  {                                                               \
	    _r -= (d);                                                    \
	    _q += 1;                                                      \
	  }                                                               \
      }                                                                   \
    if (_r >= (d))                                                        \
      {                                                                   \
	_r -= (d);                                                        \
	_q += 1;                                                          \
      }                                                                   \
    (r) = _r;                                                             \
    (q) = _q;                                                             \
  } while (0)
/* Exactly like udiv_qrnnd_preinv, but branch-free.  It is not clear which
   version to use.  */
#define udiv_qrnnd_preinv2norm(q, r, nh, nl, d, di) \
  do {									\
    mp_limb_t _n2, _n10, _n1, _nadj, _q1;				\
    mp_limb_t _xh, _xl;							\
    _n2 = (nh);								\
    _n10 = (nl);							\
    _n1 = ((mp_limb_signed_t) _n10 >> (BITS_PER_MP_LIMB - 1));		\
    _nadj = _n10 + (_n1 & (d));						\
    umul_ppmm (_xh, _xl, di, _n2 - _n1);				\
    add_ssaaaa (_xh, _xl, _xh, _xl, 0, _nadj);				\
    _q1 = ~(_n2 + _xh);							\
    umul_ppmm (_xh, _xl, _q1, d);					\
    add_ssaaaa (_xh, _xl, _xh, _xl, nh, nl);				\
    _xh -= (d);								\
    (r) = _xl + ((d) & _xh);						\
    (q) = _xh - _q1;							\
  } while (0)

#define MAX(x,y) (((x)<(y))?(y):(x))
#define MIN(x,y) (((x)<(y))?(x):(y))

#define SIZ(x) ((x)->_mp_size)
#define PTR(x) ((x)->_mp_d)


/*********
 * TYPES *
 *********/

/* SP */

/* the type for both a small prime, and a residue modulo a small prime
 *  - for a sp, we require the top bit to be set
 *  - for a residue x modulo a sp p, we require 0 <= x < p */
typedef UWtype sp_t;

#define SP_NUMB_BITS (8 * sizeof (sp_t))
#define SP_MAX ULONG_MAX
#define SP_MIN (ULONG_MAX >> 1)

/* SPM */

/* small prime modulus - this contains some precomputed constants to
 * calculate modulo a sp */
typedef struct
{
  sp_t sp;		/* value of the sp */
  sp_t mul_c;		/* constant used for reduction mod sp */
  sp_t prim_root;       /* primitive root */
  sp_t inv_prim_root;	/* inverse of prim_root */
} __spm_struct;

typedef __spm_struct *spm_t;


/* SPV */

/* vector of residues modulo a common small prime */
typedef sp_t * spv_t;

/* length of a spv */
typedef unsigned long spv_size_t;


/* MPZSPM */

typedef mpz_t * mpzv_t;

typedef struct
  {
    /* number of small primes needed to represent each coeff */
    unsigned int sp_num;
    spv_size_t max_ntt_size;
    
    mpz_t modulus;
    
    /* spm data */
    __spm_struct *spm;
    
    /* precomputed crt constants, see mpzspm.c */
    mpzv_t crt1, crt2;
    sp_t *crt3, **crt4, *crt5;
  } __mpzspm_struct;

typedef __mpzspm_struct * mpzspm_t;

/* MPZSPV */

/* sp representation of a mpz polynomial */

typedef spv_t * mpzspv_t;


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


/* sp */

/* Routines for arithmetic on residues modulo a small prime
 *
 * All functions return values in the range 0 <= x < p.
 *
 * The variable name of the modulus is 'p' if the input must be prime,
 *                                     'm' if we also allow composites. */

#define add_sssaaaaaa(s2, s1, s0, a2, a1, a0, b2, b1, b0)		\
  __asm__ ("addl %8,%2\n\tadcl %6,%1\n\tadcl %4,%0"			\
	: "=r" ((USItype)(s2)), "=&r" ((USItype)(s1)), "=&r" ((USItype)(s0)) \
	: "0" ((USItype)(a2)), "g" ((USItype)(b2)),			\
	  "1" ((USItype)(a1)), "g" ((USItype)(b1)),			\
          "2" ((USItype)(a0)), "g" ((USItype)(b0)))


/* x*y mod m */
static inline sp_t
sp_mul (sp_t x, sp_t y, sp_t m, sp_t d)
{
  sp_t z, u, v, w;
  umul_ppmm (u, v, x, y);
  udiv_qrnnd_preinv2norm (w, z, u, v, m, d);
  
  return z;
}

/* x*y mod m */
static inline sp_t
sp_sqr (sp_t x, sp_t m, sp_t d)
{
  sp_t z, u, v, w;
  umul_ppmm (u, v, x, x);
  udiv_qrnnd_preinv2norm (w, z, u, v, m, d);
  
  return z;
}

#if 0
static inline sp_t
sp_montmul (sp_t x, sp_t y, sp_t p, sp_t d)
{
  sp_t a, b, u, v, m;
  umul_ppmm (u, v, x, y);
  m = v * d;
  umul_ppmm (a, b, m, p);
  add_ssaaaa (u, v, u, v, a, b);
  return (u < p) ? u : u - p;
}

static inline sp_t
sp_montsqr (sp_t x, sp_t p, sp_t d)
{
  sp_t a, b, u, v, m;
  umul_ppmm (u, v, x, x);
  m = v * d;
  umul_ppmm (a, b, m, p);
  add_ssaaaa (u, v, u, v, a, b);
  return u - (u < p) ? 0 : p;
}
#endif

#define sp_add(x,y,m) (((x)<(m)-(y)) ? ((x)+(y)) : ((x)+(y)-(m)))
#define sp_sub(x,y,m) (((x)>=(y)) ? ((x)-(y)) : ((x)-(y)+(m)))

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
  
/* x * 2 ^ a mod m */
#if 1
static inline sp_t
sp_mul_2exp (sp_t x, sp_t a, sp_t m, sp_t d)
{
  while (a--)
    x = sp_add (x, x, m);
  
  return x;
}
#else
#define sp_mul_2exp(x,a,m,d) sp_mul(x,1<<(a),m,d)
#endif

int sp_spp (sp_t, sp_t, sp_t);
int sp_prime (sp_t);

/* spm */

spm_t spm_init (sp_t);
void spm_clear (spm_t);

/* spv */

/* ASSIGNMENT */

void spv_set (spv_t, spv_t, spv_size_t);
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
void spv_mul_sp (spv_t, spv_t, sp_t, spv_size_t, sp_t, sp_t);
/* void spv_addmul_sp (spv_t, spv_t, sp_t, spv_size_t, sp_t, sp_t); */
/* void spv_submul_sp (spv_t, spv_t, sp_t, spv_size_t, sp_t, sp_t); */

/* polynomial multiplication */
/* void spv_mul_basecase (spv_t, spv_t, spv_t, spv_size_t, spv_size_t,
    sp_t, sp_t); */
/* void spv_mul_karatsuba (spv_t, spv_t, spv_t, spv_t, spv_size_t,
    sp_t, sp_t); */
/* void spv_mul_toomcook3 (spv_t, spv_t, spv_t, spv_t, spv_size_t, spm_t); */
/* void spv_mul_toomcook4 (spv_t, spv_t, spv_t, spv_t, spv_size_t, spm_t); */
/* void spv_mul (spv_t, spv_t, spv_size_t, spv_t, spv_size_t, spv_size_t,
    spv_size_t, int, spm_t); */
/* void spv_sqr (spv_t, spv_t, spv_size_t, int, spm_t); */
void spv_random (spv_t, spv_size_t, sp_t);
int spv_cmp (spv_t, spv_t, spv_size_t);

/* ntt_gfp */

void spv_ntt_scramble (spv_t, spv_size_t);
void spv_ntt_gfp_dif (spv_t, spv_size_t, sp_t, sp_t, sp_t);
void spv_ntt_gfp_dit (spv_t, spv_size_t, sp_t, sp_t, sp_t);
/* void spv_mul_ntt_gfp (spv_t, spv_t, spv_t, spv_size_t, spm_t); */
/* void spv_sqr_ntt_gfp (spv_t, spv_t, spv_size_t, spm_t); */

/* mpzspm */

mpzspm_t mpzspm_init (spv_size_t, mpz_t);
void mpzspm_clear (mpzspm_t);

/* mpzspv */

mpzspv_t mpzspv_init (spv_size_t, mpzspm_t);
void mpzspv_clear (mpzspv_t, mpzspm_t);
int mpzspv_verify (mpzspv_t, spv_size_t, spv_size_t, mpzspm_t);
void mpzspv_set (mpzspv_t, spv_size_t, mpzspv_t, spv_size_t, spv_size_t,
    mpzspm_t);
void mpzspv_set_sp (mpzspv_t, spv_size_t, sp_t, spv_size_t, mpzspm_t);
void mpzspv_from_mpzv (mpzspv_t, spv_size_t, mpzv_t, spv_size_t, mpzspm_t);
void mpzspv_reverse (mpzspv_t, spv_size_t, spv_size_t, mpzspm_t);
void mpzspv_neg (mpzspv_t, spv_size_t, mpzspv_t, spv_size_t, spv_size_t,
    mpzspm_t);
void mpzspv_to_mpzv (mpzspv_t, spv_size_t, mpzv_t, spv_size_t, mpzspm_t);
void mpzspv_normalise (mpzspv_t, spv_size_t, spv_size_t, mpzspm_t);
void mpzspv_pwmul (mpzspv_t, spv_size_t, mpzspv_t, spv_size_t, mpzspv_t, 
    spv_size_t, spv_size_t, mpzspm_t);
void mpzspv_to_ntt (mpzspv_t, spv_size_t, spv_size_t, spv_size_t, int,
    mpzspm_t);
void mpzspv_from_ntt (mpzspv_t, spv_size_t, spv_size_t, spv_size_t, mpzspm_t);
void mpzspv_random (mpzspv_t, spv_size_t, spv_size_t, mpzspm_t);

#endif /* __HAVE_SP_H */
