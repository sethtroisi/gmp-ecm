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

#include <gmp.h>

#define HAVE_NTT

/**************
 * GMP_IMPL.H *
 **************/

#ifdef __GNUC__
#define INLINE inline
#else
#define INLINE
#endif

#if WANT_ASSERT
#include <assert.h>
#define ASSERT(expr)   assert (expr)
#else
#define ASSERT(expr)   do {} while (0)
#endif

typedef unsigned long UWtype;
typedef unsigned int UHWtype;
typedef unsigned long USItype;
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


/* SPM */

/* small prime modulus - this contains some precomputed constants to
 * calculate modulo a sp */
typedef struct
{
  sp_t sp;		/* value of the sp */
  sp_t mul_c;		/* constant used for reduction mod sp */
  sp_t prim_root;    /* primitive root */
} __spm_struct;

typedef __spm_struct *spm_t;


/* SPV */

/* vector of residues modulo a common small prime */
typedef sp_t * spv_t;

/* length of a spv */
typedef unsigned long spv_size_t;


/* SPP */

/* small prime polynomial */
typedef struct
  {
    /* allocated length of the polynomial */
    spv_size_t alloc_len;
    
    /* current length of the polynomial
     * this is an upper bound on the number of nonzero coeffs and is
     * the value used for all computations */
    spv_size_t len;

    /* location of the coefficients (constant term first) */
    spv_t spv;
    
    /* modulus */
    spm_t spm;
  } __spp_struct;

typedef __spp_struct spp_t[1];

/* MPZSPM */

typedef mpz_t * mpzp_t;

typedef struct
  {
    /* number of small primes needed to represent each coeff */
    unsigned int sp_num;
    
    mpz_t modulus;
    
    /* spm data */
    __spm_struct *spm;
    
    /* precomputed crt constants, see mpzspm.c */
    mpzp_t crt1, crt2;
    sp_t *crt3, **crt4, *crt5;
  } __mpzspm_struct;

typedef __mpzspm_struct * mpzspm_t;

/* MPZSPP */

/* sp representation of a mpz polynomial */
typedef struct
  {
    /* allocated length of each polynomial */
    spv_size_t alloc_len;   
    
    /* current length of each polynomial
     * this is an upper bound on the number of nonzero coeffs and is
     * the value used for all computations */
    spv_size_t len;

    /* pointers to the location of the coefficients
     * of each polynomial (constant term first) */
    spv_t *spv;
    
    /* flags to tell whether or not the coeffs are normalised / transformed */
    int normalised;
    int ntt;
    
    /* is there an implied monic monomial? if so then store its position,
     * otherwise monic_pos = 0 */
    spv_size_t monic_pos;
    
    /* prime and crt info */
    mpzspm_t mpzspm;
  } __mpzspp_struct;

typedef __mpzspp_struct * mpzspp_t;


/*************
 * CONSTANTS *
 *************/

/* poly length at which to start using scramble + ntt_dif + scramble instead
 * of ntt_dit */
#ifndef DIT_DIF_THRESHOLD
#define DIT_DIF_THRESHOLD 32768
#endif

/* size of the CPU's L1 cache */
#ifndef CACHE_SIZE
#define CACHE_SIZE (64 * 1024)
#endif

/* poly length at which to start using ntts for mul */
#ifndef MUL_NTT_THRESHOLD
#define MUL_NTT_THRESHOLD 1024
#endif

/* poly length at which to start using ntts for PrerevertDivision */
#ifndef PREREVERT_DIVISION_NTT_THRESHOLD
#define PREREVERT_DIVISION_NTT_THRESHOLD 1024
#endif

/* poly length at which to start using ntts for PolyInvert */
#ifndef POLYINVERT_NTT_THRESHOLD
#define POLYINVERT_NTT_THRESHOLD 1024
#endif

/*************
 * FUNCTIONS *
 *************/

/* general */

static INLINE unsigned int
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
static INLINE sp_t
sp_mul (sp_t x, sp_t y, sp_t m, sp_t d)
{
  sp_t z, u, v, w;
  umul_ppmm (u, v, x, y);
  udiv_qrnnd_preinv2norm (w, z, u, v, m, d);
  
  return z;
}

/* x*y mod m */
static INLINE sp_t
sp_sqr (sp_t x, sp_t m, sp_t d)
{
  sp_t z, u, v, w;
  umul_ppmm (u, v, x, x);
  udiv_qrnnd_preinv2norm (w, z, u, v, m, d);
  
  return z;
}

#if 0
static INLINE sp_t
sp_montmul (sp_t x, sp_t y, sp_t p, sp_t d)
{
  sp_t a, b, u, v, m;
  umul_ppmm (u, v, x, y);
  m = v * d;
  umul_ppmm (a, b, m, p);
  add_ssaaaa (u, v, u, v, a, b);
  return (u < p) ? u : u - p;
}

static INLINE sp_t
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

static INLINE sp_t
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
static INLINE sp_t
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
void spv_addmul_sp (spv_t, spv_t, sp_t, spv_size_t, sp_t, sp_t);
void spv_submul_sp (spv_t, spv_t, sp_t, spv_size_t, sp_t, sp_t);

/* polynomial multiplication */
void spv_mul_basecase (spv_t, spv_t, spv_t, spv_size_t, spv_size_t,
    sp_t, sp_t);
void spv_mul_karatsuba (spv_t, spv_t, spv_t, spv_t, spv_size_t, sp_t, sp_t);
// void spv_mul_toomcook3 (spv_t, spv_t, spv_t, spv_t, spv_size_t, spm_t);
// void spv_mul_toomcook4 (spv_t, spv_t, spv_t, spv_t, spv_size_t, spm_t);
void spv_mul (spv_t, spv_t, spv_size_t, spv_t, spv_size_t, spv_size_t,
    spv_size_t, int, spm_t);
void spv_sqr (spv_t, spv_t, spv_size_t, int, spm_t);
void spv_random (spv_t, spv_size_t, sp_t);
int spv_cmp (spv_t, spv_t, spv_size_t);

#if 0
/* spp */

/* INITIALISATION */

void spp_init (spp_t x, spm_t spm);
void spp_init2 (spp_t x, spv_size_t len, spm_t spm);
void spp_realloc (spp_t x, spv_size_t len);
void spp_clear (spp_t x);

/* ASSIGNMENT */

void spp_set (spp_t r, spp_t x);
void spp_swap (spp_t x, spp_t y);

void spp_random (spp_t x, spv_size_t len);

/* ARITHMETIC */

/* add */
void spp_add (spp_t r, spp_t x, spp_t y);

/* subtract */
void spp_sub (spp_t r, spp_t x, spp_t y);

/* polynomial product */
void spp_mul (spp_t r, spp_t x, spp_t y, int monic);
void spp_sqr (spp_t r, spp_t x, int monic);

void spp_mul_sp (spp_t r, spp_t x, sp_t c);

/* COMPARISONS */
int spp_cmp (spp_t x, spp_t y);
#endif

/* ntt_gfp */

void spv_ntt_scramble (spv_t, spv_size_t);
void spv_ntt_gfp_dif (spv_t, spv_size_t, sp_t, sp_t, sp_t);
void spv_ntt_gfp_dit (spv_t, spv_size_t, sp_t, sp_t, sp_t);
void spv_mul_ntt_gfp (spv_t, spv_t, spv_t, spv_size_t, spm_t);
void spv_sqr_ntt_gfp (spv_t, spv_t, spv_size_t, spm_t);

/* mpzspm */

mpzspm_t mpzspm_init (mpz_t, spv_size_t);
void mpzspm_clear (mpzspm_t);

/* mpzspp */

#define MPZSPP_REALLOC(x,n) do{if((x)->alloc_len<n)mpzspp_realloc(x,n);}while(0)

mpzspp_t mpzspp_init (mpzspm_t);
void mpzspp_clear (mpzspp_t);
void mpzspp_realloc (mpzspp_t, spv_size_t);
void mpzspp_set_mpzp (mpzspp_t, mpzp_t, spv_size_t, spv_size_t);
void mpzspp_get_mpzp (mpzspp_t, mpzp_t, spv_size_t, spv_size_t);
void mpzspp_mul (mpzspp_t, mpzspp_t, mpzspp_t, int);
void mpzspp_mul_partial (mpzspp_t, mpzspp_t, mpzspp_t, spv_size_t, spv_size_t,
    int);
void mpzspp_normalise (mpzspp_t, spv_size_t, spv_size_t);
void mpzspp_pwmul (mpzspp_t, mpzspp_t, mpzspp_t);
void mpzspp_to_ntt (mpzspp_t, spv_size_t, int);
void mpzspp_from_ntt (mpzspp_t);

#endif /* __HAVE_SP_H */
