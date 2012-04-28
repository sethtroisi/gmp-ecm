/* Polynomial multiplication using GMP's integer multiplication code

Copyright 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2012 Dave Newman,
Paul Zimmermann, Alexander Kruppa.

This file is part of the ECM Library.

The ECM Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The ECM Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the ECM Library; see the file COPYING.LIB.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#include <stdlib.h>
#include "ecm-gmp.h" /* for MPZ_REALLOC and MPN_COPY */
#include "ecm-impl.h"

#define FFT_WRAP /* always defined since mpn_mul_fft is included */

/* Copy at r+i*s the content of A[i*stride] for 0 <= i < l
   Assume all A[i*stride] are non-negative, and their size is <= s.
 */
static void
pack (mp_ptr r, mpz_t *A, mp_size_t l, mp_size_t stride, mp_size_t s)
{
  mp_size_t i, j, m;

  for (i = 0, j = 0; i < l; i++, j += stride, r += s)
    {
      m = SIZ(A[j]);
      ASSERT((0 <= m) && (m <= s));
      if (m)
        MPN_COPY (r, PTR(A[j]), m);
      if (m < s)
        MPN_ZERO (r + m, s - m);
    }
}

/* put in R[i*stride] for 0 <= i < l the content of {t+i*s, s} */
void
unpack (mpz_t *R, mp_size_t stride, mp_ptr t, mp_size_t l, mp_size_t s)
{
  mp_size_t i, j, size_tmp;
  mp_ptr r_ptr;

  for (i = 0, j = 0; i < l; i++, t += s, j += stride)
    {
      size_tmp = s;
      MPN_NORMALIZE(t, size_tmp); /* compute the actual size */
      r_ptr = MPZ_REALLOC (R[j], size_tmp);
      if (size_tmp)
        MPN_COPY (r_ptr, t, size_tmp);
      SIZ(R[j]) = size_tmp;
    }
}

void
list_mul_n_basecase (listz_t R, listz_t A, listz_t B, unsigned int n)
{
  unsigned int i, j;

  if (n == 1)
    {
      mpz_mul (R[0], A[0], B[0]);
      return;
    }

  if (n == 2) /* Karatsuba */
    {
      mpz_add (R[0], A[0], A[1]);
      mpz_add (R[2], B[0], B[1]);
      mpz_mul (R[1], R[0], R[2]);
      mpz_mul (R[0], A[0], B[0]);
      mpz_mul (R[2], A[1], B[1]);
      mpz_sub (R[1], R[1], R[0]);
      mpz_sub (R[1], R[1], R[2]);
      return;
    }

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      {
        if (i == 0 || j == n - 1)
          mpz_mul (R[i+j], A[i], B[j]);
        else
          mpz_addmul (R[i+j], A[i], B[j]);
      }
}

static uint64_t
gcd_uint64 (uint64_t a, uint64_t b)
{
  uint64_t t;

  ASSERT (b != 0);

  if (a >= b)
    a %= b;

  while (a > 0)
    {
      /* Here 0 < a < b */
      t = b % a; /* 0 <= t < a */
      b = a;
      a = t; /* 0 <= a < b */
    }

  return b;
}

/* Given f[0]...f[t] that contain respectively f(0), ..., f(t),
   put in f[0]...f[t] the coefficients of f.
   Assumes t <= MAX_T.
*/
static void
list_mul_tc_interpolate (mpz_t *f, int t)
{
#define MAX_T 15
  uint64_t M[MAX_T+1][MAX_T+1], g, h;
  int i, j, k;

  ASSERT_ALWAYS (t <= MAX_T); /* Ensures that all M[i][j] fit in uint64_t,
                                 and similarly for all intermediate
                                 computations on M[i][j]. This avoids the
                                 use of mpz_t to store the M[i][j]. */

  /* initialize M[i][j] = i^j */
  for (i = 0; i <= t; i++)
    for (j = 0; j <= t; j++)
      M[i][j] = (j == 0) ? 1 : i * M[i][j-1];

  /* forward Gauss: zero the under-diagonal coefficients while going down */
  for (i = 1; i <= t; i++)
    for (j = 0; j < i; j++)
      if (M[i][j] != 0)
        {
          g = gcd_uint64 (M[i][j], M[j][j]);
          h = M[i][j] / g;
          g = M[j][j] / g;
          /* f[i] <- g*f[i] - h*f[j] */
          mpz_mul_ui (f[i], f[i], g);
          mpz_submul_ui (f[i], f[j], h);
          for (k = j; k <= t; k++)
            M[i][k] = g * M[i][k] - h * M[j][k];
        }

  /* now zero upper-diagonal coefficients while going up */
  for (i = t; i >= 0; i--)
    {
      for (j = i + 1; j <= t; j++)
        /* f[i] = f[i] - M[i][j] * f[j] */
        mpz_submul_ui (f[i], f[j], M[i][j]);
      ASSERT (mpz_divisible_ui_p (f[i], M[i][i]));
      mpz_divexact_ui (f[i], f[i], M[i][i]);
    }
}

/* Generic Toom-Cook implementation: stores in f[0..r+s-2] the coefficients
   of g*h, where g has r coefficients and h has s coefficients, and their
   coefficients are in g[0..r-1] and h[0..s-1].
   Assume f differs from g and h, and f[0..r+s-2] are already allocated.
   Assume r + s - 2 <= MAX_T (MAX_T = 15 ensures all the matrix coefficients
   in the inversion fit into uint64_t, and those used in mpz_mul_ui calls
   fit into uint32_t).
   Assume r, s >= 1.
*/
void
list_mul_tc (listz_t f, listz_t g, unsigned int r, listz_t h, unsigned int s)
{
  int t, i, j;
  mpz_t tmp;

  ASSERT_ALWAYS(r >= 1 && s >= 1);

  r -= 1;    /* degree of g */
  s -= 1;    /* degree of h */
  t = r + s; /* product has t+1 coefficients */

  ASSERT_ALWAYS (t <= MAX_T);

  mpz_init (tmp);

  /* first store g(i)*h(i) in f[i] for 0 <= i <= t */
  for (i = 0; i <= t; i++)
    {
      /* f[i] <- g(i) */
      mpz_set (f[i], g[r]);
      for (j = r - 1; j >= 0; j--)
        {
          mpz_mul_ui (f[i], f[i], i);
          mpz_add (f[i], f[i], g[j]);
        }
      /* tmp <- h(i) */
      mpz_set (tmp, h[s]);
      for (j = s - 1; j >= 0; j--)
        {
          mpz_mul_ui (tmp, tmp, i);
          mpz_add (tmp, tmp, h[j]);
        }
      /* f[i] <- g(i)*h(i) */
      mpz_mul (f[i], f[i], tmp);
    }

  list_mul_tc_interpolate (f, t);

  mpz_clear (tmp);
}

/* Classical one-point Kronecker-Schoenhage substitution.
   Notes:
    - this code aligns the coeffs at limb boundaries - if instead we aligned
      at byte boundaries then we could save up to 3*n bytes,
      but tests have shown this doesn't give any significant speed increase,
      even for large degree polynomials.
     - this code requires that all coefficients A[] and B[] are nonnegative. */
void
list_mul_n_KS1 (listz_t R, listz_t A, listz_t B, unsigned int l)
{
  unsigned long i;
  mp_size_t s, t = 0, size_t0;
  mp_ptr t0_ptr, t1_ptr, t2_ptr;

  /* compute the largest bit-size t of the A[i] and B[i] */
  for (i = 0; i < l; i++)
    {
      if ((s = mpz_sizeinbase (A[i], 2)) > t)
        t = s;
      if ((s = mpz_sizeinbase (B[i], 2)) > t)
        t = s;
    }
  /* For n > 0, s = sizeinbase (n, 2)     ==> n < 2^s. 
     For n = 0, s = sizeinbase (n, 2) = 1 ==> n < 2^s.
     Hence all A[i], B[i] < 2^t */
  
  /* Each coeff of A(x)*B(x) < l * 2^(2*t), so max number of bits in a 
     coeff of the product will be 2 * t + ceil(log_2(l)) */
  s = 2 * t;
  for (i = l; i > 1; s++, i = (i + 1) >> 1);
  
  /* work out the corresponding number of limbs */
  s = 1 + (s - 1) / GMP_NUMB_BITS;

  size_t0 = s * l;

  /* allocate a single buffer to save malloc/MPN_ZERO/free calls */
  t0_ptr = (mp_ptr) malloc (4 * size_t0 * sizeof (mp_limb_t));
  if (t0_ptr == NULL)
    {
      outputf (OUTPUT_ERROR, "Out of memory in list_mult_n()\n");
      exit (1);
    }
  t1_ptr = t0_ptr + size_t0;
  t2_ptr = t1_ptr + size_t0;
    
  pack (t0_ptr, A, l, 1, s);
  pack (t1_ptr, B, l, 1, s);

  mpn_mul_n (t2_ptr, t0_ptr, t1_ptr, size_t0);

  unpack (R, 1, t2_ptr, 2 * l - 1, s);

  free (t0_ptr);
}

/* Two-point Kronecker substitition.
   Reference: Algorithm 2 from "Faster polynomial multiplication via multipoint
   Kronecker substitution", David Harvey, Journal of Symbolic Computation,
   number 44 (2009), pages 1502-1510.
   Assume n >= 2.
   Notes:
    - this code aligns the coeffs at limb boundaries - if instead we aligned
      at byte boundaries then we could save up to 3*n bytes,
      but tests have shown this doesn't give any significant speed increase,
      even for large degree polynomials.
     - this code requires that all coefficients A[] and B[] are nonnegative.
*/
void
list_mul_n_KS2 (listz_t R, listz_t A, listz_t B, unsigned int n)
{
  unsigned long i;
  mp_size_t s, s2, t = 0, l, h, ns2;
  mp_ptr tmp, A0, A1, B0, B1, C0, C1;
  int sA, sB;

  ASSERT_ALWAYS (n >= 2);

  /* compute the largest bit-size t of the A[i] and B[i] */
  for (i = 0; i < n; i++)
    {
      if ((s = mpz_sizeinbase (A[i], 2)) > t)
        t = s;
      if ((s = mpz_sizeinbase (B[i], 2)) > t)
        t = s;
    }
  /* For n > 0, s = sizeinbase (n, 2)     ==> n < 2^s. 
     For n = 0, s = sizeinbase (n, 2) = 1 ==> n < 2^s.
     Hence all A[i], B[i] < 2^t */
  
  /* Each coeff of A(x)*B(x) < n * 2^(2*t), so max number of bits in a 
     coeff of the product will be 2 * t + ceil(log_2(n)) */
  s = 2 * t;
  for (i = n; i > 1; s++, i = (i + 1) >> 1);
  
  /* work out the corresponding number of limbs */
  s = 1 + (s - 1) / GMP_NUMB_BITS;

  /* ensure s is even */
  s = s + (s & 1);
  s2 = s >> 1;
  ns2 = n * s2;

  l = n / 2;
  h = n - l;

  /* allocate a single buffer to save malloc/MPN_ZERO/free calls */
  tmp = (mp_ptr) malloc (8 * ns2 * sizeof (mp_limb_t));
  if (tmp == NULL)
    {
      outputf (OUTPUT_ERROR, "Out of memory in list_mult_n()\n");
      exit (1);
    }

  A0 = tmp;
  A1 = A0 + ns2;
  B0 = A1 + ns2;
  B1 = B0 + ns2;
  C0 = B1 + ns2;
  C1 = C0 + 2 * ns2;

  pack (A0, A, h, 2, s); /* A0 = Aeven(S) where S = 2^(s*GMP_NUMB_BITS) */
  /* A0 has in fact only n * s2 significant limbs:
     if n=2h, h*s = n*s2
     if n=2h-1, the last chunk from A0 has at most s2 limbs */
  MPN_ZERO(B0, s2);
  pack (B0 + s2, A + 1, l, 2, s);
  /* for the same reason as above, we have at most l*s-s2 significant limbs
     at B0+s2, thus at most l*s <= n*s2 at B0 */
  if ((sA = mpn_cmp (A0, B0, ns2)) >= 0)
    mpn_sub_n (A1, A0, B0, ns2);
  else
    mpn_sub_n (A1, B0, A0, ns2);
  mpn_add_n (A0, A0, B0, ns2);
  /* now A0 is X+ with the notations of Algorithm, A1 is sA*X- */

  pack (B0, B, h, 2, s);
  MPN_ZERO(C0, s2);
  pack (C0 + s2, B + 1, l, 2, s);
  if ((sB = mpn_cmp (B0, C0, ns2)) >= 0)
    mpn_sub_n (B1, B0, C0, ns2);
  else
    mpn_sub_n (B1, C0, B0, ns2);
  mpn_add_n (B0, B0, C0, ns2);
  /* B0 is Y+, B1 is sB*Y- with the notations of Algorithm 2 */

  mpn_mul_n (C0, A0, B0, ns2); /* C0 is Z+ = X+ * Y+ */
  mpn_mul_n (C1, A1, B1, ns2); /* C1 is sA * sB * Z- */

  if (sA * sB >= 0)
    {
      mpn_add_n (A0, C0, C1, 2 * ns2);
      mpn_sub_n (B0, C0, C1, 2 * ns2);
    }
  else
    {
      mpn_sub_n (A0, C0, C1, 2 * ns2);
      mpn_add_n (B0, C0, C1, 2 * ns2);
    }
  mpn_rshift (A0, A0, 4 * ns2, 1); /* global division by 2 */

  /* If A[] and B[] have n coefficients, the product has 2n-1 coefficients.
     The even part has n coefficients and the odd part n-1 coefficients */
  unpack (R, 2, A0, n, s);
  unpack (R + 1, 2, B0 + s2, n - 1, s);

  free (tmp);
}

/* Puts in R[0..2n-2] the product of A[0..n-1] and B[0..n-1], seen as
   polynomials.
*/
void
list_mult_n (listz_t R, listz_t A, listz_t B, unsigned int n)
{
  int T[TUNE_LIST_MUL_N_MAX_SIZE] = LIST_MUL_TABLE, best;

  /* See tune_list_mul_n() in tune.c:
     0 : list_mul_n_basecase
     1 : list_mul_tc
     2 : list_mul_n_KS1
     3 : list_mul_n_KS2 */
  best = (n < TUNE_LIST_MUL_N_MAX_SIZE) ? T[n] : 3;
    
  if (best == 0)
    list_mul_n_basecase (R, A, B, n);
  else if (best == 1)
    list_mul_tc (R, A, n, B, n);
  else if (best == 2)
    list_mul_n_KS1 (R, A, B, n);
  else
    list_mul_n_KS2 (R, A, B, n);
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

#undef TEST_OLD_S

int
TMulKS (listz_t b, unsigned int n, listz_t a, unsigned int m,
        listz_t c, unsigned int l, mpz_t modulus, int rev)
{
  unsigned long i, s = 0, t, k;
  mp_ptr ap, bp, cp;
  mp_size_t an, bn, cn;
  int ret = 0; /* default return value */
#ifdef TEST_OLD_S
  unsigned long s_old = 0, k_old;
  mp_size_t bn_old;
#endif
#ifdef DEBUG
  long st = cputime ();
  fprintf (ECM_STDOUT, "n=%u m=%u l=%u bits=%u n*bits=%u: ", n, m, l,
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

#ifdef TEST_OLD_S
  /* We used max(m,l) before. We compute the corresponding s for 
     comparison. */
  for (s_old = 2 * s, i = (m > l) ? m : l; i; s_old++, i >>= 1);
#endif

  /* max coeff has 2*s+ceil(log2(min(m+1,l+1))) bits,
   i.e. 2*s + 1 + floor(log2(min(m,l))) */
  for (s = 2 * s, i = (m < l) ? m : l; i; s++, i >>= 1);

  /* corresponding number of limbs */
  s = 1 + (s - 1) / GMP_NUMB_BITS;
#ifdef TEST_OLD_S
  s_old = 1 + (s_old - 1) / GMP_NUMB_BITS;
#endif

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
#ifdef TEST_OLD_S
  k_old = mpn_fft_best_k ((m + n + 1) * s_old, 0);
  if (k != k_old)
    outputf (OUTPUT_ERROR, 
             "Got different FFT transform length, k = %lu, k_old : %lu\n",
             k, k_old);    
  bn_old = mpn_fft_next_size ((m + n + 1) * s_old, k_old);
  if (bn != bn_old)
    outputf (OUTPUT_ERROR, "Got different FFT size, bn = %d, bn_old : %d\n",
             (int) bn, (int) bn_old);
#endif
  
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
      MPZ_REALLOC (b[i], (mp_size_t) t);
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
  fprintf (ECM_STDOUT, "%ums\n", elltime (st, cputime ()));
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
    fprintf (ECM_STDOUT, "+%lu*B^%u", np[i], i);
  fprintf (ECM_STDOUT, "\n");
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

  ASSERT(k >= l);

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
  ASSERT(m <= 2 * m0 - 3 + list_mul_mem (m0 - 1));

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
      if ((size_tmp == s) && (tp[s - 1] >> (GMP_NUMB_BITS - 1)))
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
