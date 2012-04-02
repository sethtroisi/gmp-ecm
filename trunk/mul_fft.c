/* This code was extracted from GMP 5.0.4, with the copyright notice below.

   Schoenhage's fast multiplication modulo 2^N+1.

   Contributed by Paul Zimmermann.

   THE FUNCTIONS IN THIS FILE ARE INTERNAL WITH MUTABLE INTERFACES.  IT IS ONLY
   SAFE TO REACH THEM THROUGH DOCUMENTED INTERFACES.  IN FACT, IT IS ALMOST
   GUARANTEED THAT THEY WILL CHANGE OR DISAPPEAR IN A FUTURE GNU MP RELEASE.

Copyright 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008,
2009, 2010 Free Software Foundation, Inc.

This file is part of the GNU MP Library.

The GNU MP Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The GNU MP Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MP Library.  If not, see http://www.gnu.org/licenses/.  */

#define HAVE_NATIVE_mpn_add_n_sub_n 0

#include "ecm-impl.h"
#include "ecm-gmp.h"
#define TRACE(x)
#define ASSERT_MPN_ZERO_P(x,y)

#ifndef HAVE___GMPN_MUL_FFT_FULL
/* multiply {n, nl} by {m, ml}, and put the result in {op, nl+ml} */
void
mpn_mul_fft_full (mp_ptr op,
		  mp_srcptr n, mp_size_t nl,
		  mp_srcptr m, mp_size_t ml)
{
  mp_ptr pad_op;
  mp_size_t pl, pl2, pl3, l;
  int k2, k3;
  int sqr = (n == m && nl == ml);
  int cc, c2, oldcc;

  pl = nl + ml; /* total number of limbs of the result */

  /* perform a fft mod 2^(2N)+1 and one mod 2^(3N)+1.
     We must have pl3 = 3/2 * pl2, with pl2 a multiple of 2^k2, and
     pl3 a multiple of 2^k3. Since k3 >= k2, both are multiples of 2^k2,
     and pl2 must be an even multiple of 2^k2. Thus (pl2,pl3) =
     (2*j*2^k2,3*j*2^k2), which works for 3*j <= pl/2^k2 <= 5*j.
     We need that consecutive intervals overlap, i.e. 5*j >= 3*(j+1),
     which requires j>=2. Thus this scheme requires pl >= 6 * 2^FFT_FIRST_K. */

  /*  ASSERT_ALWAYS(pl >= 6 * (1 << FFT_FIRST_K)); */

  pl2 = (2 * pl - 1) / 5; /* ceil (2pl/5) - 1 */
  do
    {
      pl2++;
      k2 = mpn_fft_best_k (pl2, sqr); /* best fft size for pl2 limbs */
      pl2 = mpn_fft_next_size (pl2, k2);
      pl3 = 3 * pl2 / 2; /* since k>=FFT_FIRST_K=4, pl2 is a multiple of 2^4,
			    thus pl2 / 2 is exact */
      k3 = mpn_fft_best_k (pl3, sqr);
    }
  while (mpn_fft_next_size (pl3, k3) != pl3);

  TRACE (printf ("mpn_mul_fft_full nl=%ld ml=%ld -> pl2=%ld pl3=%ld k=%d\n",
		 nl, ml, pl2, pl3, k2));

  ASSERT_ALWAYS(pl3 <= pl);
  cc = mpn_mul_fft (op, pl3, n, nl, m, ml, k3);     /* mu */
  ASSERT(cc == 0);
  pad_op = malloc (pl2 * sizeof(mp_limb_t));
  if (pad_op == 0)
    {
      outputf (OUTPUT_ERROR, "Error: not enough memory in mpn_mul_fft_full\n");
      return;
    }
  cc = mpn_mul_fft (pad_op, pl2, n, nl, m, ml, k2); /* lambda */
  cc = -cc + mpn_sub_n (pad_op, pad_op, op, pl2);    /* lambda - low(mu) */
  /* 0 <= cc <= 1 */
  ASSERT(0 <= cc && cc <= 1);
  l = pl3 - pl2; /* l = pl2 / 2 since pl3 = 3/2 * pl2 */
  c2 = mpn_add_n (pad_op, pad_op, op + pl2, l);
  cc = mpn_add_1 (pad_op + l, pad_op + l, l, (mp_limb_t) c2) - cc;
  ASSERT(-1 <= cc && cc <= 1);
  if (cc < 0)
    cc = mpn_add_1 (pad_op, pad_op, pl2, (mp_limb_t) -cc);
  ASSERT(0 <= cc && cc <= 1);
  /* now lambda-mu = {pad_op, pl2} - cc mod 2^(pl2*GMP_NUMB_BITS)+1 */
  oldcc = cc;
#if HAVE_NATIVE_mpn_add_n_sub_n
  c2 = mpn_add_n_sub_n (pad_op + l, pad_op, pad_op, pad_op + l, l);
  /* c2 & 1 is the borrow, c2 & 2 is the carry */
  cc += c2 >> 1; /* carry out from high <- low + high */
  c2 = c2 & 1; /* borrow out from low <- low - high */
#else
  {
    mp_ptr tmp;
    TMP_DECL;

    TMP_MARK;
    tmp = TMP_ALLOC_LIMBS (l);
    MPN_COPY (tmp, pad_op, l);
    c2 = mpn_sub_n (pad_op,      pad_op, pad_op + l, l);
    cc += mpn_add_n (pad_op + l, tmp,    pad_op + l, l);
    TMP_FREE;
  }
#endif
  c2 += oldcc;
  /* first normalize {pad_op, pl2} before dividing by 2: c2 is the borrow
     at pad_op + l, cc is the carry at pad_op + pl2 */
  /* 0 <= cc <= 2 */
  cc -= mpn_sub_1 (pad_op + l, pad_op + l, l, (mp_limb_t) c2);
  /* -1 <= cc <= 2 */
  if (cc > 0)
    cc = -mpn_sub_1 (pad_op, pad_op, pl2, (mp_limb_t) cc);
  /* now -1 <= cc <= 0 */
  if (cc < 0)
    cc = mpn_add_1 (pad_op, pad_op, pl2, (mp_limb_t) -cc);
  /* now {pad_op, pl2} is normalized, with 0 <= cc <= 1 */
  if (pad_op[0] & 1) /* if odd, add 2^(pl2*GMP_NUMB_BITS)+1 */
    cc += 1 + mpn_add_1 (pad_op, pad_op, pl2, (mp_limb_t) 1);
  /* now 0 <= cc <= 2, but cc=2 cannot occur since it would give a carry
     out below */
  mpn_rshift (pad_op, pad_op, pl2, 1); /* divide by two */
  if (cc) /* then cc=1 */
    pad_op [pl2 - 1] |= (mp_limb_t) 1 << (GMP_NUMB_BITS - 1);
  /* now {pad_op,pl2}-cc = (lambda-mu)/(1-2^(l*GMP_NUMB_BITS))
     mod 2^(pl2*GMP_NUMB_BITS) + 1 */
  c2 = mpn_add_n (op, op, pad_op, pl2); /* no need to add cc (is 0) */
  /* since pl2+pl3 >= pl, necessary the extra limbs (including cc) are zero */
  MPN_COPY (op + pl3, pad_op, pl - pl3);
  ASSERT_MPN_ZERO_P (pad_op + pl - pl3, pl2 + pl3 - pl);
  free (pad_op);
  /* since the final result has at most pl limbs, no carry out below */
  mpn_add_1 (op + pl2, op + pl2, pl - pl2, (mp_limb_t) c2);
}
#endif
