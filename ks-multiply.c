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
  if ((double) l * (double) s < KS_MUL_THRESHOLD)
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

  /* Note: s * (l - 1) + ceil(t/GMP_NUMB_BITS) should be faster,
     but no significant speedup was observed */
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

/* Given a[0..m] and c[0..l], puts in b[0..n] the coefficients
   of degree m to n+m of rev(a)*c, i.e.
   b[0] = a[0]*c[0] + ... + a[i]*c[i] with i = min(m, l)
   ...
   b[k] = a[0]*c[k] + ... + a[i]*c[i+k] with i = min(m, l-k)
   ...
   b[n] = a[0]*c[n] + ... + a[i]*c[i+n] with i = min(m, l-n) [=l-n].

   Assumes n <= l.
*/
unsigned int
TMulKS (listz_t b, unsigned int n,
        listz_t a, unsigned int m, listz_t c, unsigned int l, mpz_t modulus)
{
  unsigned long i, s = 0, t, k;
  mp_ptr ap, bp, cp;
  mp_size_t an, bn, cn;
#ifdef DEBUG
  int st = cputime ();
  printf ("n=%u m=%u l=%u bits=%u n*bits=%u: ", n, m, l,
          mpz_sizeinbase (modulus, 2), n * mpz_sizeinbase (modulus, 2));
#endif

  ASSERT (n <= l); /* otherwise the upper coefficients of b are 0 */
  if (l > n + m)
    l = n + m; /* otherwise, c has too many coeffs */

  /* compute max bits of a[] and c[] */
  for (i = 0; i <= m; i++)
    {
      if (mpz_sgn (a[i]) < 0)
        mpz_mod (a[i], a[i], modulus);
      if ((t = mpz_sizeinbase (a[i], 2)) > s)
        s = t;
    }
  for (i = 0; i <= l; i++)
    {
      if (mpz_sgn (c[i]) < 0)
        mpz_mod (c[i], c[i], modulus);
      if ((t = mpz_sizeinbase (c[i], 2)) > s)
        s = t;
    }

#define FFT_WRAP
#ifdef FFT_WRAP
  s ++; /* need one extra bit to determine sign of low(b) - high(b) */
#endif

  /* max coeff has 2*s+ceil(log2(max(m+1,l+1))) bits,
   i.e. 2*s + 1 + floor(log2(max(m,l))) */
  for (s = 2 * s, i = (m > l) ? m : l; i; s++, i >>= 1);

  /* corresponding number of limbs */
  s = 1 + (s - 1) / GMP_NUMB_BITS;

  an = (m + 1) * s;
  cn = (l + 1) * s;
  bn = an + cn;

  /* a[0..m] needs (m+1) * s limbs */
  ap = (mp_ptr) malloc (an * sizeof (mp_limb_t));
  cp = (mp_ptr) malloc (cn * sizeof (mp_limb_t));

  MPN_ZERO (ap, an);
  MPN_ZERO (cp, cn);

  /* a is reverted */
  for (i = 0; i <= m; i++)
    MPN_COPY (ap + (m - i) * s, PTR(a[i]), SIZ(a[i]));
  for (i = 0; i <= l; i++)
    MPN_COPY (cp + i * s, PTR(c[i]), SIZ(c[i]));

#ifdef FFT_WRAP
  /* the product rev(a) * c has m+l+1 coefficients.
     We throw away the first m and the last l-n <= m.
     If we compute mod (m+n+1) * s limbs, we are ok */
  k = mpn_fft_best_k ((m + n + 1) * s, 0);
  bn = mpn_fft_next_size ((m + n + 1) * s, k);
  bp = (mp_ptr) malloc ((bn + 1) * sizeof (mp_limb_t));
  mpn_mul_fft (bp, bn, ap, an, cp, cn, k);
  if (bp[m * s - 1] >> (GMP_NUMB_BITS - 1)) /* lo(b)-hi(b) is negative */
    mpn_add_1 (bp + m * s, bp + m * s, (n + 1) * s, (mp_limb_t) 1);
#else
  bp = (mp_ptr) malloc (bn * sizeof (mp_limb_t));
  if (an >= cn)
    mpn_mul (bp, ap, an, cp, cn);
  else
    mpn_mul (bp, cp, cn, ap, an);
#endif

  /* recover coefficients of degree m to n+m of product in b[0..n] */
  bp += m * s;
  for (i = 0; i <= n; i++)
    {
      t = s;
      MPN_NORMALIZE(bp, t);
      if (t==0) abort();
      _mpz_realloc (b[i], t);
      MPN_COPY (PTR(b[i]), bp, t);
      SIZ(b[i]) = t;
      bp += s;
    }
  bp -= (m + n + 1) * s;

  free (ap);
  free (cp);
  free (bp);

#ifdef DEBUG
  printf ("%dms\n", cputime () - st);
#endif
  
  return 0;
}
