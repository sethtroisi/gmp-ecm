/* Polynomial multiplication using GMP's integer multiplication code

  Copyright 2004, 2005 Dave Newman, Paul Zimmermann.

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
#include <stdlib.h>
#include <assert.h>

#include "gmp.h"
#include "ecm-gmp.h" /* for MPZ_REALLOC and MPN_COPY */
#include "ecm-impl.h"

#ifdef HAVE_FFT
#define FFT_WRAP
#endif

/* Puts in R[0..2l-2] the product of A[0..l-1] and B[0..l-1].
   T must have as much space as for toomcook4 (it is only used when that
   function is called).
   Notes:
    - this code aligns the coeffs at limb boundaries - if instead we aligned
      at byte boundaries then we could save up to 3*l bytes in T0 and T1,
      but tests have shown this doesn't give any significant speed increase,
      even for large degree polynomials.
    - this code requires that all coefficients A[] and B[] are nonnegative.
    Return non-zero if an error occurred.
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
    {
      toomcook4 (R, A, B, l, T);
      return 0;
    }

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
  if (t0_ptr == NULL)
    return 1;
  t1_ptr = t0_ptr + size_t0;
    
  MPN_ZERO (t0_ptr, size_t0 + size_t0);

  for (i = 0; i < l; i++)
    {
      if (SIZ(A[i]))
        MPN_COPY (t0_ptr + i * s, PTR(A[i]), SIZ(A[i]));
      if (SIZ(B[i]))
        MPN_COPY (t1_ptr + i * s, PTR(B[i]), SIZ(B[i]));
    }

  t2_ptr = (mp_ptr) malloc (2 * size_t0 * sizeof (mp_limb_t));
  if (t2_ptr == NULL)
    {
      free (t0_ptr);
      return 1;
    }

  mpn_mul_n (t2_ptr, t0_ptr, t1_ptr, size_t0);
  
  for (i = 0; i < 2 * l - 1; i++)
    {
      size_tmp = s;
      MPN_NORMALIZE(t2_ptr + i * s, size_tmp);
      r_ptr = MPZ_REALLOC (R[i], size_tmp);
      if (size_tmp)
        MPN_COPY (r_ptr, t2_ptr + i * s, size_tmp);
      SIZ(R[i]) = size_tmp;
    }

  free (t0_ptr);
  free (t2_ptr);

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

   Return non-zero if an error occurred.
*/
int
TMulKS (listz_t b, unsigned int n, listz_t a, unsigned int m,
        listz_t c, unsigned int l, mpz_t modulus, int rev)
{
  unsigned long i, s = 0, t, k;
  mp_ptr ap, bp, cp;
  mp_size_t an, bn, cn;
  int ret = 0; /* default return value */
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
  if (ap == NULL)
    {
      ret = 1;
      goto TMulKS_end;
    }
  cp = (mp_ptr) malloc (cn * sizeof (mp_limb_t));
  if (cp == NULL)
    {
      ret = 1;
      goto TMulKS_free_ap;
    }

  MPN_ZERO (ap, an);
  MPN_ZERO (cp, cn);

  /* a is reverted */
  for (i = 0; i <= m; i++)
    if (SIZ(a[i]))
      MPN_COPY (ap + ((rev) ? (m - i) : i) * s, PTR(a[i]), SIZ(a[i]));
  for (i = 0; i <= l; i++)
    if (SIZ(c[i]))
      MPN_COPY (cp + i * s, PTR(c[i]), SIZ(c[i]));

#ifdef FFT_WRAP
  /* the product rev(a) * c has m+l+1 coefficients.
     We throw away the first m and the last l-n <= m.
     If we compute mod (m+n+1) * s limbs, we are ok */
  k = mpn_fft_best_k ((m + n + 1) * s, 0);
  bn = mpn_fft_next_size ((m + n + 1) * s, k);
  bp = (mp_ptr) malloc ((bn + 1) * sizeof (mp_limb_t));
  if (bp == NULL)
    {
      ret = 1;
      goto TMulKS_free_cp;
    }
  mpn_mul_fft (bp, bn, ap, an, cp, cn, k);
  if (m && bp[m * s - 1] >> (GMP_NUMB_BITS - 1)) /* lo(b)-hi(b) is negative */
    mpn_add_1 (bp + m * s, bp + m * s, (n + 1) * s, (mp_limb_t) 1);
#else
  bp = (mp_ptr) malloc (bn * sizeof (mp_limb_t));
  if (bp == NULL)
    {
      ret = 1;
      goto TMulKS_free_cp;
    }
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
      _mpz_realloc (b[i], t);
      if (t)
        MPN_COPY (PTR(b[i]), bp, t);
      SIZ(b[i]) = t;
      bp += s;
    }
  bp -= (m + n + 1) * s;

  free (bp);
 TMulKS_free_cp:
  free (cp);
 TMulKS_free_ap:
  free (ap);

#ifdef DEBUG
  printf ("%dms\n", cputime () - st);
#endif
  
 TMulKS_end:
  return ret;
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

unsigned int
ks_wrapmul_m (unsigned int m0, unsigned int k, mpz_t n)
{
  mp_size_t t, s;
  unsigned long i, fft_k, m;

  t = mpz_sizeinbase (n, 2);
  s = t * 2 + 1;
  for (i = k - 1; i; s++, i >>= 1);
  s = 1 + (s - 1) / GMP_NUMB_BITS;
  fft_k = mpn_fft_best_k (m0 * s, 0);
  i = mpn_fft_next_size (m0 * s, fft_k);
  while (i % s)
    i = mpn_fft_next_size (i + 1, fft_k);
  m = i / s;
  return m;
}

/* multiply in R[] A[0]+A[1]*x+...+A[k-1]*x^(k-1)
                by B[0]+B[1]*x+...+B[l-1]*x^(l-1) modulo n,
   wrapping around coefficients of the product up from degree m >= m0.
   Assumes k >= l.
   R is assumed to have 2*m0-3+list_mul_mem(m0-1) allocated cells.
   Return m (or 0 if an error occurred).
*/
unsigned int
ks_wrapmul (listz_t R, unsigned int m0,
            listz_t A, unsigned int k,
            listz_t B, unsigned int l,
	    mpz_t n)
{
  unsigned long i, fft_k, m, t;
  mp_size_t s, size_t0, size_t1, size_tmp;
  mp_ptr t0_ptr, t1_ptr, t2_ptr, r_ptr, tp;
  int negative;

#ifdef DEBUG
  if (k < l) abort();
#endif

  t = mpz_sizeinbase (n, 2);
  for (i = 0; i < k; i++)
    if (mpz_sgn (A[i]) < 0 || mpz_sizeinbase (A[i], 2) > t)
      mpz_mod (A[i], A[i], n);
  for (i = 0; i < l; i++)
    if (mpz_sgn (B[i]) < 0 || mpz_sizeinbase (B[i], 2) > t)
      mpz_mod (B[i], B[i], n);
  
  s = t * 2 + 1; /* one extra sign bit */
  for (i = k - 1; i; s++, i >>= 1);
  
  s = 1 + (s - 1) / GMP_NUMB_BITS;

  size_t0 = s * k;
  size_t1 = s * l;

  /* allocate one double-buffer to save malloc/MPN_ZERO/free calls */
  t0_ptr = (mp_ptr) malloc (size_t0 * sizeof (mp_limb_t));
  if (t0_ptr == NULL)
    return 0;
  t1_ptr = (mp_ptr) malloc (size_t1 * sizeof (mp_limb_t));
  if (t1_ptr == NULL)
    {
      free (t0_ptr);
      return 0;
    }
    
  MPN_ZERO (t0_ptr, size_t0);
  MPN_ZERO (t1_ptr, size_t1);

  for (i = 0; i < k; i++)
    if (SIZ(A[i]))
      MPN_COPY (t0_ptr + i * s, PTR(A[i]), SIZ(A[i]));
  for (i = 0; i < l; i++)
    if (SIZ(B[i]))
      MPN_COPY (t1_ptr + i * s, PTR(B[i]), SIZ(B[i]));

  fft_k = mpn_fft_best_k (m0 * s, 0);
  i = mpn_fft_next_size (m0 * s, fft_k);
  /* the following loop ensures we don't cut in the middle of a
     coefficient */
  while (i % s)
    i = mpn_fft_next_size (i + 1, fft_k);
  m = i / s;
#ifdef DEBUG  
  if (m > 2 * m0 - 3 + list_mul_mem (m0 - 1))
    {
      fprintf (stderr, "Error: m > 2*m0-3+list_mul_mem(m0-1) in ks_wrapmul: m=%u m0=%u\n",
	       m, m0);
      exit (1);
    }
#endif

  t2_ptr = (mp_ptr) malloc ((i + 1) * sizeof (mp_limb_t));
  if (t2_ptr == NULL)
    {
      free (t0_ptr);
      free (t1_ptr);
      return 0;
    }

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
      if (size_tmp)
        MPN_COPY (r_ptr, tp, size_tmp);
      SIZ(R[i]) = (negative) ? -size_tmp : size_tmp;
      tp += s;
    }

  free (t0_ptr);
  free (t1_ptr);
  free (t2_ptr);
  
  return m;
}
