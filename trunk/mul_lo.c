/* Low-half short product.

  Copyright 2003 Paul Zimmermann and Alexander Kruppa.

  This program is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 2 of the License, or (at your
  option) any later version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
  more details.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.
*/

#include "gmp.h"
#ifdef WANT_GMP_IMPL
#include "gmp-impl.h"
#endif
#include "ecm.h"

void mpn_mul_lo_basecase (mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);

/* puts in {rp, n} the low part of {np, n} times {mp, n}
   i.e. equivalent to:
   tp = TMP_ALLOC_LIMBS (2 * n);
   mpn_mul_n (tp, np, mp, n);
   MPN_COPY (rp, tp, n);
   Remark: it is assumed at least 2*n limbs are allocated starting from rp.
 */
INLINE void
mpn_mul_lo_basecase (mp_ptr rp, mp_srcptr np, mp_srcptr mp, mp_size_t n)
{
  mpn_mul_1 (rp, np, n, mp[0]);
  for (; --n;)
    mpn_addmul_1 (++rp, np, n, (++mp)[0]);
}

#define MPN_MUL_LO_THRESHOLD (2 * MUL_KARATSUBA_THRESHOLD)

INLINE void
mpn_mul_lo_n (mp_ptr rp, mp_srcptr np, mp_srcptr mp, mp_size_t n)
{
  if (n < MPN_MUL_LO_THRESHOLD)
    mpn_mul_lo_basecase (rp, np, mp, n);
  else
    {
      mp_size_t k = (mp_size_t ) (0.75 * (double) n);

      mpn_mul_n (rp, np, mp, k);
      rp += k;
      n -= k;
      mpn_mul_lo_n (rp + n, np + k, mp, n);
      mpn_add_n (rp, rp, rp + n, n);
      mpn_mul_lo_n (rp + n, np, mp + k, n);
      mpn_add_n (rp, rp, rp + n, n);
    }
}
