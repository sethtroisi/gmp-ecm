/* Polynomial multiplication using GMP's integer multiplication code

  Original code: Copyright 2004 Dave Newman
  david (dot) newman [at] jesus {dot} ox <dot> ac $dot$ uk
  Modified by Paul Zimmermann.

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

#include <stdlib.h>
#include <assert.h>

#include "gmp.h"

#ifdef WANT_GMP_IMPL
#include "gmp-impl.h"
#else
#include "ecm-gmp.h" /* for MPZ_REALLOC and MPN_COPY */
#endif /* WANT_GMP_IMPL */

#include "ecm.h"

/* Notes:
    - this code aligns the coeffs at limb boundaries - if instead we aligned
      at byte boundaries then we could save up to 3*l bytes in T0 and T1,
      but tests have shown this doesn't give any significant speed increase,
      even for large degree polynomials.
    - this code requires that all coefficients A[] and B[] are nonnegative.
*/    
int
kronecker_schonhage (listz_t R, listz_t A, listz_t B, unsigned int l,
                     listz_t T)
{
  unsigned long i;
  mp_size_t s, t = 0, size_t0, size_tmp;
  mp_ptr t0_ptr, t1_ptr, t2_ptr, r_ptr;

  s = mpz_sizeinbase (A[0], 2);
  if ((double) l * (double) s < 1e6)
    return toomcook4 (R, A, B, l, T);

  for (i = 0; i < l; i++)
    {
      if ((s = mpz_sizeinbase (A[i], 2)) > t)
        t = s;
      if ((s = mpz_sizeinbase (B[i], 2)) > t)
        t = s;
    }
  
  /* max number of bits in a coeff of T[0] * T[1] will be
     2 * t + ceil(log_2(l)) */
  s = t * 2;
  for (i = l - 1; i; s++, i >>= 1); /* ceil(log_2(l)) = 1+floor(log_2(l-1)) */
  
  /* work out the corresponding number of limbs */
  s = 1 + (s - 1) / GMP_NUMB_BITS;

  size_t0 = s * l;

  /* allocate one double-buffer to save malloc/MPN_ZERO/free calls */
  t0_ptr = (mp_ptr) malloc (2 * size_t0 * sizeof (mp_limb_t));
  t1_ptr = t0_ptr + size_t0;
    
  MPN_ZERO (t0_ptr, size_t0 + size_t0);

  for (i = 0; i < l; i++)
    {
      MPN_COPY (t0_ptr + i * s, PTR(A[i]), SIZ(A[i]));
      MPN_COPY (t1_ptr + i * s, PTR(B[i]), SIZ(B[i]));
    }

  t2_ptr = (mp_ptr) malloc (2 * size_t0 * sizeof (mp_limb_t));

  mpn_mul_n (t2_ptr, t0_ptr, t1_ptr, size_t0);
  
  for (i = 0; i < 2 * l - 1; i++)
    {
      size_tmp = s;
      MPN_NORMALIZE(t2_ptr + i * s, size_tmp);
      r_ptr = MPZ_REALLOC (R[i], size_tmp);
      MPN_COPY (r_ptr, t2_ptr + i * s, size_tmp);
      SIZ(R[i]) = size_tmp;
    }

  free (t0_ptr);
  free (t2_ptr);
  
  /* we don't have a measure of how many multiplies we've done */
  return 0;
}

