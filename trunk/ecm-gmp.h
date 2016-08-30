/* Part of file gmp-impl.h from GNU MP.

Copyright 1991, 1993, 1994, 1995, 1996, 1997, 1999, 2000, 2001, 2002 Free
Software Foundation, Inc.

This file contains modified code from the GNU MP Library.

The GNU MP Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The GNU MP Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MP Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
MA 02110-1301, USA. */

#ifndef _ECM_GMP_H
#define _ECM_GMP_H 1

#include "config.h"
#include "basicdefs.h"
#include <gmp.h>

#define ABSIZ(x) ABS (SIZ (x))
#define ALLOC(x) ((x)->_mp_alloc)
#define PTR(x) ((x)->_mp_d)
#define SIZ(x) ((x)->_mp_size)
#define TMP_DECL
#define TMP_ALLOC(x) alloca(x)
#define TMP_MARK
#define TMP_FREE
#define TMP_ALLOC_TYPE(n,type) ((type *) TMP_ALLOC ((n) * sizeof (type)))
#define TMP_ALLOC_LIMBS(n)     TMP_ALLOC_TYPE(n,mp_limb_t)

#ifndef MPZ_REALLOC
#define MPZ_REALLOC(z,n) ((n) > ALLOC(z) ? _mpz_realloc(z,n) : PTR(z))
#endif

#ifndef MPN_COPY
#include <string.h> /* for memcpy */
#define MPN_COPY(d,s,n) memcpy((d),(s),(n)*sizeof(mp_limb_t))
#endif

#ifndef MPN_NORMALIZE
#define MPN_NORMALIZE(DST, NLIMBS) \
  do {									\
    while (NLIMBS > 0)							\
      {									\
	if ((DST)[(NLIMBS) - 1] != 0)					\
	  break;							\
	NLIMBS--;							\
      }									\
  } while (0)
#endif

#ifndef MPN_ZERO
#define MPN_ZERO(dst, n)			\
  do {						\
    if ((n) != 0)				\
      {						\
	mp_ptr __dst = (dst);			\
	mp_size_t __n = (n);			\
	do					\
	  *__dst++ = 0;				\
	while (--__n);				\
      }						\
  } while (0)
#endif

/* Return non-zero if xp,xsize and yp,ysize overlap.
   If xp+xsize<=yp there's no overlap, or if yp+ysize<=xp there's no
   overlap.  If both these are false, there's an overlap. */
#define MPN_OVERLAP_P(xp, xsize, yp, ysize) \
  ((xp) + (xsize) > (yp) && (yp) + (ysize) > (xp))

/* Return non-zero if xp,xsize and yp,ysize are either identical or not
   overlapping.  Return zero if they're partially overlapping. */
#define MPN_SAME_OR_SEPARATE_P(xp, yp, size)    \
  MPN_SAME_OR_SEPARATE2_P(xp, size, yp, size)
#define MPN_SAME_OR_SEPARATE2_P(xp, xsize, yp, ysize)           \
  ((xp) == (yp) || ! MPN_OVERLAP_P (xp, xsize, yp, ysize))

#ifndef mpn_com_n
#define mpn_com_n(d,s,n)                                \
  do {                                                  \
    mp_ptr     __d = (d);                               \
    mp_srcptr  __s = (s);                               \
    mp_size_t  __n = (n);                               \
    ASSERT (__n >= 1);                                  \
    ASSERT (MPN_SAME_OR_SEPARATE_P (__d, __s, __n));    \
    do                                                  \
      *__d++ = (~ *__s++) & GMP_NUMB_MASK;              \
    while (--__n);                                      \
  } while (0)
#endif

#ifdef HAVE___GMPN_ADD_NC
#ifndef __gmpn_add_nc
__GMP_DECLSPEC mp_limb_t __gmpn_add_nc (mp_ptr, mp_srcptr, mp_srcptr,
    mp_size_t, mp_limb_t);
#endif
#endif

#define ECM_VERSION_NUM(a,b,c) (((a) << 16L) | ((b) << 8) | (c))

#if !defined( __MPIR_RELEASE ) && ECM_VERSION_NUM(__GNU_MP_VERSION,__GNU_MP_VERSION_MINOR,__GNU_MP_VERSION_PATCHLEVEL) >= ECM_VERSION_NUM(5,1,0)
#define MPN_REDC12_RETURNS_CARRY 1
#endif

/* GMP currently does not define prototypes for these, but MPIR does */
#if defined(HAVE___GMPN_REDC_1) && !defined( __MPIR_RELEASE )
#ifdef MPN_REDC12_RETURNS_CARRY
  mp_limb_t __gmpn_redc_1 (mp_ptr, mp_ptr, mp_srcptr, mp_size_t, mp_limb_t);
#else
  void __gmpn_redc_1 (mp_ptr, mp_ptr, mp_srcptr, mp_size_t, mp_limb_t);
#endif
#endif

#if defined(HAVE___GMPN_REDC_2) && !defined( __MPIR_RELEASE )
#ifdef MPN_REDC12_RETURNS_CARRY
  mp_limb_t __gmpn_redc_2 (mp_ptr, mp_ptr, mp_srcptr, mp_size_t, mp_srcptr);
#else
  void __gmpn_redc_2 (mp_ptr, mp_ptr, mp_srcptr, mp_size_t, mp_srcptr);
#endif
#endif

#if defined(HAVE___GMPN_REDC_N)
  void __gmpn_redc_n (mp_ptr, mp_ptr, mp_srcptr, mp_size_t, mp_srcptr);
#endif

#if defined(HAVE___GMPN_MULLO_N)
  void __gmpn_mullo_n (mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);
#endif

#if defined( __MPIR_RELEASE ) && __MPIR_RELEASE >= 20600

#if __MPIR_RELEASE == 20600
#error MPIR 2.6 does not support GMP-ECM, please use an alternative version
#endif
/* WARNING - the following two definitions map the internal interface
   of the new FFT in MPIR 2.6 (and later) to the GMP FFT interface -
   they work in this context but the parameters for mpn_fft_next_size
   and fft_adjust_limbs have different semantics, which means that
   these definitions may fail if used in other circumstances 
*/
#  define mpn_fft_best_k(n, k)             (0) 
#  define mpn_fft_next_size(n, k)          fft_adjust_limbs(n)

#else /* we are using the real GMP */

#define mpn_mul_fft __gmpn_mul_fft
mp_limb_t __gmpn_mul_fft (mp_ptr, mp_size_t, mp_srcptr, mp_size_t, mp_srcptr, 
                          mp_size_t, int);

#define mpn_mulmod_bnm1 __gmpn_mulmod_bnm1
void mpn_mulmod_bnm1 (mp_ptr, mp_size_t, mp_srcptr, mp_size_t, mp_srcptr, 
                      mp_size_t, mp_ptr);

#define mpn_mul_fft_full __gmpn_mul_fft_full
void __gmpn_mul_fft_full (mp_ptr, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t);

#define mpn_fft_next_size __gmpn_fft_next_size
mp_size_t __gmpn_fft_next_size (mp_size_t, int);

#define mpn_mulmod_bnm1_next_size __gmpn_mulmod_bnm1_next_size
mp_size_t mpn_mulmod_bnm1_next_size (mp_size_t);

#define mpn_fft_best_k __gmpn_fft_best_k
int __gmpn_fft_best_k (mp_size_t, int);

#endif

static inline void 
mpz_set_uint64 (mpz_t m, const uint64_t n)
{
#if GMP_LIMB_BITS == 64  /* 64-bit GMP limb */
  if (sizeof(mp_limb_t) > sizeof(unsigned long))
    {
       mpz_set_ui (m, (unsigned long)(uint32_t)(n >> 32));
       mpz_mul_2exp (m, m, 32);
       mpz_add_ui (m, m, (unsigned long)(uint32_t)n);
    }
  else
    mpz_set_ui(m, (unsigned long)n);

#else                    /* 32-bit GMP limb */
  mpz_set_ui (m, (unsigned long)(uint32_t)(n >> 32));
  mpz_mul_2exp (m, m, 32);
  mpz_add_ui (m, m, (unsigned long)(uint32_t)n);
#endif
}

static inline void 
mpz_set_int64 (mpz_t m, const int64_t n)
{
  if (n < 0)
    {
      mpz_set_uint64(m, -n);
      mpz_neg(m, m);
    }
  else
    mpz_set_uint64(m, n);
}


#endif /* _ECM_GMP_H */
