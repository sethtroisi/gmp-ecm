/* Polynomial multiplication using GMP's integer multiplication code

  Original code: Copyright 2004 Dave Newman <david.newman@jesus.ox.ac.uk>
  Modified by Paul Zimmermann, 2004, 2005.

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

/* Puts in R[0..2l-2] the product of A[0..l-1] and B[0..l-1].
   T must have as much space as for toomcook4 (it is only used when that
   function is called).
   Notes:
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
  t0_ptr = (mp_ptr) xmalloc (2 * size_t0 * sizeof (mp_limb_t));
  t1_ptr = t0_ptr + size_t0;
    
  MPN_ZERO (t0_ptr, size_t0 + size_t0);

  for (i = 0; i < l; i++)
    {
      MPN_COPY (t0_ptr + i * s, PTR(A[i]), SIZ(A[i]));
      MPN_COPY (t1_ptr + i * s, PTR(B[i]), SIZ(B[i]));
    }

  t2_ptr = (mp_ptr) xmalloc (2 * size_t0 * sizeof (mp_limb_t));

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

   If rev=0, consider a instead of rev(a).

   Assumes n <= l.
*/
unsigned int
TMulKS (listz_t b, unsigned int n, listz_t a, unsigned int m,
        listz_t c, unsigned int l, mpz_t modulus, int rev)
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
  ap = (mp_ptr) xmalloc (an * sizeof (mp_limb_t));
  cp = (mp_ptr) xmalloc (cn * sizeof (mp_limb_t));

  MPN_ZERO (ap, an);
  MPN_ZERO (cp, cn);

  /* a is reverted */
  for (i = 0; i <= m; i++)
    MPN_COPY (ap + ((rev) ? (m - i) : i) * s, PTR(a[i]), SIZ(a[i]));
  for (i = 0; i <= l; i++)
    MPN_COPY (cp + i * s, PTR(c[i]), SIZ(c[i]));

#ifdef FFT_WRAP
  /* the product rev(a) * c has m+l+1 coefficients.
     We throw away the first m and the last l-n <= m.
     If we compute mod (m+n+1) * s limbs, we are ok */
  k = mpn_fft_best_k ((m + n + 1) * s, 0);
  bn = mpn_fft_next_size ((m + n + 1) * s, k);
  bp = (mp_ptr) xmalloc ((bn + 1) * sizeof (mp_limb_t));
  mpn_mul_fft (bp, bn, ap, an, cp, cn, k);
  if (m && bp[m * s - 1] >> (GMP_NUMB_BITS - 1)) /* lo(b)-hi(b) is negative */
    mpn_add_1 (bp + m * s, bp + m * s, (n + 1) * s, (mp_limb_t) 1);
#else
  bp = (mp_ptr) xmalloc (bn * sizeof (mp_limb_t));
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

#ifdef DEBUG
void
mpn_print (mp_ptr np, mp_size_t nn)
{
  mp_size_t i;
  for (i = 0; i < nn; i++)
    printf ("+%lu*B^%u", np[i], i);
  printf ("\n");
}
#endif

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

/* multiply in R[] A[0]+A[1]*x+...+A[k-1]*x^(k-1)
                by B[0]+B[1]*x+...+B[l-1]*x^(l-1) modulo n,
   wrapping around coefficients of the product up from degree m >= m0.
   Return m.
   Assumes k >= l.
*/
unsigned int
ks_wrapmul (listz_t R, unsigned int m0,
            listz_t A, unsigned int k,
            listz_t B, unsigned int l)
{
  unsigned long i, fft_k, m;
  mp_size_t s, t = 0, size_t0, size_t1, size_tmp;
  mp_ptr t0_ptr, t1_ptr, t2_ptr, r_ptr, tp;
  int negative;

#ifdef DEBUG
  if (k < l) abort();
#endif

  for (i = 0; i < k; i++)
    if ((s = mpz_sizeinbase (A[i], 2)) > t)
      t = s;
  for (i = 0; i < l; i++)
    if ((s = mpz_sizeinbase (B[i], 2)) > t)
      t = s;
  
  s = t * 2 + 1; /* one extra sign bit */
  for (i = k - 1; i; s++, i >>= 1);
  
  s = 1 + (s - 1) / GMP_NUMB_BITS;

  size_t0 = s * k;
  size_t1 = s * l;

  /* allocate one double-buffer to save malloc/MPN_ZERO/free calls */
  t0_ptr = (mp_ptr) xmalloc (size_t0 * sizeof (mp_limb_t));
  t1_ptr = (mp_ptr) xmalloc (size_t1 * sizeof (mp_limb_t));
    
  MPN_ZERO (t0_ptr, size_t0);
  MPN_ZERO (t1_ptr, size_t1);

  for (i = 0; i < k; i++)
    MPN_COPY (t0_ptr + i * s, PTR(A[i]), SIZ(A[i]));
  for (i = 0; i < l; i++)
    MPN_COPY (t1_ptr + i * s, PTR(B[i]), SIZ(B[i]));

  fft_k = mpn_fft_best_k (m0 * s, 0);
  i = mpn_fft_next_size (m0 * s, fft_k);
  /* the following loop ensures we don't cut in the middle of a
     coefficient */
  while (i % s)
    i = mpn_fft_next_size (i + 1, fft_k);
  m = i / s;

  t2_ptr = (mp_ptr) xmalloc ((i + 1) * sizeof (mp_limb_t));

  mpn_mul_fft (t2_ptr, i, t0_ptr, size_t0, t1_ptr, size_t1, fft_k);
  
  for (i = 0, tp = t2_ptr, negative = 0; i < m; i++)
    {
      size_tmp = s;
      if (negative) /* previous was negative, add 1 */
	mpn_add_1 (tp, tp, s, (mp_limb_t) 1);
      /* no need to check return value of mpn_add_1: if 1, then {tp, s}
         is now identically 0, and should remain so */
      MPN_NORMALIZE(tp, size_tmp);
      if ((size_tmp == s) && (tp[s - 1] >> (mp_bits_per_limb - 1)))
	{
	  negative = 1;
	  mpn_com_n (tp, tp, s);
	  mpn_add_1 (tp, tp, s, (mp_limb_t) 1);
	}
      else
	negative = 0;
      r_ptr = MPZ_REALLOC (R[i], size_tmp);
      MPN_COPY (r_ptr, tp, size_tmp);
      SIZ(R[i]) = (negative) ? -size_tmp : size_tmp;
      tp += s;
    }

  free (t0_ptr);
  free (t1_ptr);
  free (t2_ptr);
  
  return m;
}
