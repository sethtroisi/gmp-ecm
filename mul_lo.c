/* Low-half short product (quadratic and Mulders' algorithms).

  Copyright 2003, 2005 Paul Zimmermann and Alexander Kruppa.

  This file is part of the ECM Library.

  The ECM Library is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or (at your
  option) any later version.

  The ECM Library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
  License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with the ECM Library; see the file COPYING.LIB.  If not, write to
  the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
  MA 02111-1307, USA.
*/

#include <stdio.h>
#include "gmp.h"
#include "ecm-impl.h"

/* puts in {rp, n} the low part of {np, n} times {mp, n}, i.e. equivalent to:

   mp_ptr tp;
   TMP_DECL(marker);
   TMP_MARK(marker);
   tp = TMP_ALLOC_LIMBS (2 * n);
   mpn_mul_n (tp, np, mp, n);
   MPN_COPY (rp, tp, n);
   TMP_FREE(marker);
 */
static INLINE void
ecm_mul_lo_basecase (mp_ptr rp, mp_srcptr np, mp_srcptr mp, mp_size_t n)
{
  mpn_mul_1 (rp, np, n, mp[0]);
  for (; --n;)
    mpn_addmul_1 (++rp, np, n, (++mp)[0]);
}

#define MPN_MUL_LO_THRESHOLD 32
static mp_size_t threshold[MPN_MUL_LO_THRESHOLD] =
/* {0,0,1,1,1,1,1,1,1,0,1,1,1,1,1,11,11,11,13,13,13,13,15,15,17,17,15,16,17,17,17,1}; */
{0,0,0,0,0,1,1,0,0,1,1,1,10,9,1,1,12,1,1,1,1,15,16,14,18,19,17,18,19,20,21,22};


void
ecm_mul_lo_n (mp_ptr rp, mp_srcptr np, mp_srcptr mp, mp_size_t n)
{
  mp_size_t k;

  if (n < MPN_MUL_LO_THRESHOLD)
    {
      switch (k = threshold[n])
        {
        case 0:
          {
            mpn_mul_n (rp, np, mp, n);
            return;
          }
        case 1:
          {
            ecm_mul_lo_basecase (rp, np, mp, n);
            return;
          }
          /* else go through */
        }
    }
  else
    k = (mp_size_t) (0.75 * (double) n);

  mpn_mul_n (rp, np, mp, k);
  rp += k;
  n -= k;
  ecm_mul_lo_n (rp + n, np + k, mp, n);
  mpn_add_n (rp, rp, rp + n, n);
  ecm_mul_lo_n (rp + n, np, mp + k, n);
  mpn_add_n (rp, rp, rp + n, n);
}
