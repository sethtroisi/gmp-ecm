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
  the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
  MA 02110-1301, USA.
*/

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
void
ecm_mul_lo_basecase (mp_ptr rp, mp_srcptr np, mp_srcptr mp, mp_size_t n)
{
  mpn_mul_1 (rp, np, n, mp[0]);
  for (; --n;)
    mpn_addmul_1 (++rp, np, n, (++mp)[0]);
}

#ifdef MPN_MUL_LO_THRESHOLD_TABLE
size_t mpn_mul_lo_threshold[MPN_MUL_LO_THRESHOLD] = MPN_MUL_LO_THRESHOLD_TABLE;
#else
size_t mpn_mul_lo_threshold[MPN_MUL_LO_THRESHOLD];
#endif


void
ecm_mul_lo_n (mp_ptr rp, mp_srcptr np, mp_srcptr mp, mp_size_t n)
{
  mp_size_t k;

  if (n < MPN_MUL_LO_THRESHOLD)
    {
      switch (k = mpn_mul_lo_threshold[n])
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
