/* Part of file gmp-impl.h from GNU MP.

Copyright 1991, 1993, 1994, 1995, 1996, 1997, 1999, 2000, 2001, 2002 Free
Software Foundation, Inc.

This file contains modified code from the GNU MP Library.

The GNU MP Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

The GNU MP Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MP Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#ifndef alloca
#ifdef __GNUC__
# define alloca __builtin_alloca
#elif defined (__DECC)
# define alloca(x) __ALLOCA(x)
#elif defined (_MSC_VER)
# include <malloc.h>
# define alloca _alloca
#elif HAVE_ALLOCA_H || defined (sun)
# include <alloca.h>
#elif defined (_AIX) || defined (_IBMR2)
#pragma alloca
#else
  char *alloca ();
#endif
#endif

#define ABSIZ(x) ABS (SIZ (x))
#define ALLOC(x) ((x)->_mp_alloc)
#define PTR(x) ((x)->_mp_d)
#define SIZ(x) ((x)->_mp_size)
#define TMP_DECL(m)
#define TMP_ALLOC(x) alloca(x)
#define TMP_MARK(m)
#define TMP_FREE(m)
#define TMP_ALLOC_TYPE(n,type) ((type *) TMP_ALLOC ((n) * sizeof (type)))
#define TMP_ALLOC_LIMBS(n)     TMP_ALLOC_TYPE(n,mp_limb_t)

#ifndef MPZ_REALLOC
#define MPZ_REALLOC(z,n) ((n) > ALLOC(z) ? _mpz_realloc(z,n) : PTR(z))
#endif

#ifndef MPN_COPY
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
    ASSERT ((n) >= 0);				\
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

#ifdef HAVE_NATIVE_mpn_add_nc
#ifndef __gmpn_add_nc
mp_limb_t __gmpn_add_nc (mp_ptr, mp_srcptr, mp_srcptr, mp_size_t, mp_limb_t);
#endif
#endif

/* fft stuff */
#ifndef mpn_fft_best_k
#define mpn_fft_best_k __MPN(fft_best_k)
int     mpn_fft_best_k (mp_size_t, int);
#endif

#ifndef   mpn_fft_next_size
#define   mpn_fft_next_size __MPN(fft_next_size)
mp_size_t mpn_fft_next_size (mp_size_t, int);
#endif

#ifndef mpn_mul_fft
#define mpn_mul_fft  __MPN(mul_fft)
void    mpn_mul_fft (mp_ptr, mp_size_t, mp_srcptr, mp_size_t, mp_srcptr,
                     mp_size_t, int);
#endif
