/* Arithmetic on lists of residues modulo n.

Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2012
Paul Zimmermann and Alexander Kruppa.

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
#include "ecm-impl.h"

#ifdef DEBUG
#define ASSERTD(x) assert(x)
#else
#define ASSERTD(x)
#endif

extern unsigned int Fermat;

/* returns a bound on the auxiliary memory needed by list_mult_n */
int
list_mul_mem (unsigned int len)
{
  return 2 * len;
}

/* creates a list of n integers, return NULL if error */
listz_t
init_list (unsigned int n)
{
  listz_t p;
  unsigned int i;

  p = (mpz_t*) malloc (n * sizeof (mpz_t));
  if (p == NULL)
    return NULL;
  for (i = 0; i < n; i++)
    mpz_init (p[i]);
  return p;
}

/* creates a list of n integers, return NULL if error. Allocates each
   mpz_t to the size of N bits */
listz_t
init_list2 (unsigned int n, unsigned int N)
{
  listz_t p;
  unsigned int i;

  p = (mpz_t*) malloc (n * sizeof (mpz_t));
  if (p == NULL)
    return NULL;
  for (i = 0; i < n; i++)
    mpz_init2 (p[i], N);
  return p;
}

/* clears a list of n integers */
void
clear_list (listz_t p, unsigned int n)
{
  unsigned int i;

  if (p == NULL)
    return;
  for (i = 0; i < n; i++)
    mpz_clear (p[i]);
  free (p);
}

#ifdef DEBUG
/* prints a list of n coefficients as a polynomial */
void
print_list (listz_t p, unsigned int n)
{
  unsigned int i;

  for (i = 0; i < n; i++)
    {
      if (i > 0 && mpz_cmp_ui (p[i], 0) >= 0)
        fprintf (ECM_STDOUT, "+");
      mpz_out_str (ECM_STDOUT, 10, p[i]);
      fprintf (ECM_STDOUT, "*x^%u", i);
    }
  fprintf (ECM_STDOUT, "\n");
}

static int
list_check (listz_t a, unsigned int l, mpz_t n)
{
  unsigned int i;

  for (i = 0; i < l; i++)
    if (mpz_cmp_ui (a[i], 0) < 0 || mpz_cmp (n, a[i]) <= 0)
      {
        fprintf (ECM_STDOUT, "l=%u i=%u\n", l, i);
        mpz_out_str (ECM_STDOUT, 10, a[i]);
	fprintf (ECM_STDOUT, "\n");
        return 0;
      }
  return 1;
}
#endif /* DEBUG */

/* Read all entries in list from stream. 
   Return 0 on success, ECM_ERROR on error */
int
list_inp_raw (listz_t a, FILE *f, unsigned int n)
{
  unsigned int i;
  
  for (i = 0; i < n; i++)
    if (mpz_inp_raw (a[i], f) == 0)
      return ECM_ERROR;
  
  return 0;
}

/* Write all entries in list to stream. 
   Return 0 on success, ECM_ERROR on error */
int
list_out_raw (FILE *f, listz_t a, unsigned int n)
{
  unsigned int i;
  
  for (i = 0; i < n; i++)
    if (mpz_out_raw (f, a[i]) == 0)
      return ECM_ERROR;
  
  return 0;
}

/* p <- q */
void
list_set (listz_t p, listz_t q, unsigned int n)
{
  unsigned int i;

  for (i = 0; i < n; i++)
    mpz_set (p[i], q[i]);
}

/* p[0] <-> p[n-1], p[1] <-> p[n-2], ... */
void
list_revert (listz_t p, unsigned int n)
{
  unsigned int i;

  for (i = 0; i < n - 1 - i; i++)
    mpz_swap (p[i], p[n - 1 - i]);
}

void
list_swap (listz_t p, listz_t q, unsigned int n)
{
  unsigned int i;

  for (i = 0; i < n; i++)
    mpz_swap (p[i], q[i]);
}

/* p <- -q, keeps residues normalized */
void
list_neg (listz_t p, listz_t q, unsigned int l, mpz_t n)
{
  unsigned int i;

  for (i = 0; i < l; i++)
    {
      if (mpz_sgn (q[i]))
        mpz_sub (p[i], n, q[i]);
      else
        mpz_set_ui (p[i], 0);
    }
}

/* p <- q modulo mod */
void
list_mod (listz_t p, listz_t q, unsigned int n, mpz_t mod)
{
  unsigned int i;

  for (i = 0; i < n; i++)
    mpz_mod (p[i], q[i], mod);
}

/* p <- q + r */
void
list_add (listz_t p, listz_t q, listz_t r, unsigned int l)
{
  unsigned int i;

  for (i = 0; i < l; i++)
    mpz_add (p[i], q[i], r[i]);
}

/* p <- q - r */
void
list_sub (listz_t p, listz_t q, listz_t r, unsigned int l)
{
  unsigned int i;

  for (i = 0; i < l; i++)
    mpz_sub (p[i], q[i], r[i]);
}

/* Multiply up the integers in l, modulo n. Each entry becomes the
   product (mod n) of itself and all previous entries */
   
void 
list_mulup (listz_t l, unsigned int k, mpz_t n, mpz_t t)
{
  unsigned int i;
  
  for (i = 1; i < k; i++)
    {
      mpz_mul (t, l[i - 1], l[i]);
      mpz_mod (l[i], t, n);
    }
}

/* p <- 0 */
void
list_zero (listz_t p, unsigned int n)
{
  unsigned int i;

  for (i = 0; i < n; i++)
    mpz_set_ui (p[i], 0);
}

/* puts in a[K-1]..a[2K-2] the K high terms of the product 
   of b[0..K-1] and c[0..K-1].
   Assumes K >= 1, and a[0..2K-2] exist.
   Needs space for list_mul_mem(K) in t.
*/
void
list_mul_high (listz_t a, listz_t b, listz_t c, unsigned int K)
{
  list_mult_n (a, b, c, K);
}

/* multiplies b[0]+...+b[k-1]*x^(k-1)+x^k by c[0]+...+c[l-1]*x^(l-1)+x^l
   and puts the results in a[0]+...+a[k+l-1]*x^(k+l-1)
   [the leading monomial x^(k+l) is implicit].
   If monic is 0, don't consider x^k in b (and x^l in c).
   Assumes k = l or k = l+1.
   The auxiliary array t contains at least list_mul_mem(l) entries.
   a and t should not overlap.
*/
void
list_mul (listz_t a, listz_t b, unsigned int k,
          listz_t c, unsigned int l, int monic, listz_t t)
{
  unsigned int i, po2;

  ASSERT(k == l || k == l + 1);
  
  for (po2 = l; (po2 & 1) == 0; po2 >>= 1);
  po2 = (po2 == 1);

#ifdef DEBUG
  if (Fermat && !(po2 && l == k))
    fprintf (ECM_STDOUT, "list_mul: Fermat number, but poly lengths %d and %d\n", k, l);
#endif

  if (po2 && Fermat)
    {
      if (monic && l == k)
        {
          F_mul (a, b, c, l, MONIC, Fermat, t);
          monic = 0;
        }
      else
        F_mul (a, b, c, l, DEFAULT, Fermat, t);
    }
  else
    list_mult_n (a, b, c, l); /* set a[0]...a[2l-2] */

  if (k > l) /* multiply b[l]*x^l by c[0]+...+c[l-1]*x^(l-1) */
    {
      for (i = 0; i < l - 1; i++)
        mpz_addmul (a[l+i], b[l], c[i]);
      mpz_mul (a[2*l-1], b[l], c[l-1]);
    }

  /* deal with x^k and x^l */
  if (monic)
    {
      mpz_set_ui (a[k + l - 1], 0);
      
      /* Single pass over a[] */

      /* a += b * x^l + c * x^k, so a[i] += b[i-l]; a[i] += c[i-k] 
         if 0 <= i-l < k  or  0 <= i-k < l, respectively */
      if (k > l) /* case k = l+1 */
        mpz_add (a[l], a[l], b[0]);
      for (i = k; i < k + l; i++)
        {
          mpz_add (a[i], a[i], b[i-l]); /* i-l < k */
          mpz_add (a[i], a[i], c[i-k]); /* i-k < l */
        }
    }
}

/*
  Multiplies b[0..k-1] by c[0..k-1], stores the result in a[0..2k-2],
  and stores the reduced product in a2[0..2k-2].
  (Here, there is no implicit monic leading monomial.)
  Requires at least list_mul_mem(k) cells in t.
 */
void
list_mulmod (listz_t a2, listz_t a, listz_t b, listz_t c, unsigned int k,
              listz_t t, mpz_t n)
{
  int i;

  /* keep the semicolon on a separate line to silence a warning with clang */
  for (i = k; (i & 1) == 0; i >>= 1)
    ;
  
  ASSERTD(list_check(b,k,n));
  ASSERTD(list_check(c,k,n));
  if (i == 1 && Fermat)
    F_mul (a, b, c, k, DEFAULT, Fermat, t);
  else
    list_mult_n (a, b, c, k); /* set a[0]...a[2l-2] */

  list_mod (a2, a, 2 * k - 1, n);
}

/* puts in G[0]..G[k-1] the coefficients from (x+a[0])...(x+a[k-1])
   Warning: doesn't fill the coefficient 1 of G[k], which is implicit.
   Needs k + list_mul_mem(k/2) cells in T.
   G == a is allowed. T must not overlap with anything else.
*/
void
PolyFromRoots (listz_t G, listz_t a, unsigned int k, listz_t T, mpz_t n)
{
  unsigned int l, m;

  ASSERT (T != G && T != a);
  ASSERT (k >= 1);

  if (k == 1)
    {
      /* we consider x + a[0], which mean we consider negated roots */
      mpz_mod (G[0], a[0], n);
      return;
    }

  m = k / 2; /* m >= 1 */
  l = k - m; /* l >= 1 */
  
  PolyFromRoots (G, a, l, T, n);
  PolyFromRoots (G + l, a + l, m, T, n);
  list_mul (T, G, l, G + l, m, 1, T + k);
  list_mod (G, T, k, n);
}

/* puts in G[0]..G[k-1] the coefficients from (x+a[0])...(x+a[k-1])
   Warning: doesn't fill the coefficient 1 of G[k], which is implicit.
   Needs k + list_mul_mem(k/2) cells in T.
   The product tree is stored in:
   G[0..k-1]       (degree k)
   Tree[0][0..k-1] (degree k/2)
   Tree[1][0..k-1] (degree k/4), ...,
   Tree[lgk-1][0..k-1] (degree 1)
   (then we should have initially Tree[lgk-1] = a).

   The parameter dolvl signals that only level 'dolvl' of
   the tree should be computed (dolvl < 0 means all levels).

   Either Tree <> NULL and TreeFile == NULL, and we write the tree to memory,
   or Tree == NULL and TreeFile <> NULL, and we write the tree to disk.
*/
int
PolyFromRoots_Tree (listz_t G, listz_t a, unsigned int k, listz_t T, 
               int dolvl, mpz_t n, listz_t *Tree, FILE *TreeFile, 
               unsigned int sh)
{
  unsigned int l, m;
  listz_t H1, *NextTree;

  ASSERT (k >= 1);

  if (k == 1)
    {
      /* we consider x + a[0], which mean we consider negated roots */
      mpz_mod (G[0], a[0], n);
      return 0;
    }

  if (Tree == NULL) /* -treefile case */
    {
      H1 = G;
      NextTree = NULL;
    }
  else
    {
      H1 = Tree[0] + sh;
      NextTree = Tree + 1;
    }

  m = k / 2;
  l = k - m;
  
  if (dolvl != 0) /* either dolvl < 0 and we need to compute all levels,
                     or dolvl > 0 and we need first to compute lower levels */
    {
      PolyFromRoots_Tree (H1, a, l, T, dolvl - 1, n, NextTree, TreeFile, sh);
      PolyFromRoots_Tree (H1 + l, a + l, m, T, dolvl - 1, n, NextTree, 
                          TreeFile, sh + l);
    }
  if (dolvl <= 0)
    {
      /* Write this level to disk, if requested */
      if (TreeFile != NULL)
        {
          if (list_out_raw (TreeFile, H1, l) == ECM_ERROR ||
              list_out_raw (TreeFile, H1 + l, m) == ECM_ERROR)
            {
              outputf (OUTPUT_ERROR, "Error writing product tree of F\n");
              return ECM_ERROR;
            }
        }
      list_mul (T, H1, l, H1 + l, m, 1, T + k);
      list_mod (G, T, k, n);
    }
  
  return 0; 
}

/* puts in q[0..K-1] the quotient of x^(2K-2) by B
   where B = b[0]+b[1]*x+...+b[K-1]*x^(K-1) with b[K-1]=1.
*/
void
PolyInvert (listz_t q, listz_t b, unsigned int K, listz_t t, mpz_t n)
{
  if (K == 1)
    {
      mpz_set_ui (q[0], 1);
      return;
    }
  else
    {
      int k, l, po2, use_middle_product = 0;

      use_middle_product = 1;

      k = K / 2;
      l = K - k;

      for (po2 = K; (po2 & 1) == 0; po2 >>= 1);
      po2 = (po2 == 1 && Fermat != 0);

      /* first determine l most-significant coeffs of Q */
      PolyInvert (q + k, b + k, l, t, n); /* Q1 = {q+k, l} */

      /* now Q1 * B = x^(2K-2) + O(x^(2K-2-l)) = x^(2K-2) + O(x^(K+k-2)).
         We need the coefficients of degree K-1 to K+k-2 of Q1*B */

      ASSERTD(list_check(q+k,l,n) && list_check(b,l,n));
      if (po2 == 0 && use_middle_product)
        {
          TMulKS (t, k - 1, q + k, l - 1, b, K - 1, n, 0);
          list_neg (t, t, k, n);
        }
      else if (po2)
        {
          list_revert (q + k, l);
          /* This expects the leading monomials explicitly in q[2k-1] and b[k+l-1] */
          F_mul_trans (t, q + k, b, K / 2, K, Fermat, t + k);
          list_revert (q + k, l);
          list_neg (t, t, k, n);
        }
      else
        {
          list_mult_n (t, q + k, b, l); /* t[0..2l-1] = Q1 * B0 */
          list_neg (t, t + l - 1, k, n);
      
          if (k > 1)
            {
              list_mul (t + k, q + k, l - 1, b + l, k - 1, 1,
			t + k + K - 2); /* Q1 * B1 */
              list_sub (t + 1, t + 1, t + k, k - 1);
            }
        }
      list_mod (t, t, k, n); /* high(1-B*Q1) */

      ASSERTD(list_check(t,k,n) && list_check(q+l,k,n));
      if (po2)
        F_mul (t + k, t, q + l, k, DEFAULT, Fermat, t + 3 * k);
      else
        list_mult_n (t + k, t, q + l, k);
      list_mod (q, t + 2 * k - 1, k, n);
    }
}

/*
  Returns in a[0]+a[1]*x+...+a[K-1]*x^(K-1)
  the remainder of the division of
  A = a[0]+a[1]*x+...+a[2K-2]*x^(2K-2)
  by B = b[0]+b[1]*x+...+b[K-1]*x^(K-1)+b[K]*x^K with b[K]=1 *explicit*.
  (We have A = Q*B + R with deg(Q)=K-2 and deg(R)=K-1.)
  Assumes invb[0]+invb[1]*x+...+invb[K-2]*x^(K-2) equals Quo(x^(2K-2), B).
  Assumes K >= 2.
  Requires 2K-1 + list_mul_mem(K) cells in t.

  Notations: R = r[0..K-1], A = a[0..2K-2], low(A) = a[0..K-1],
  high(A) = a[K..2K-2], Q = t[0..K-2]
  Return non-zero iff an error occurred.
*/
int
PrerevertDivision (listz_t a, listz_t b, listz_t invb,
                   unsigned int K, listz_t t, mpz_t n)
{
  int po2, wrap;
  listz_t t2 = NULL;

  wrap = ks_wrapmul_m (K + 1, K + 1, n) <= 2 * K - 1 + list_mul_mem (K);

  /* Q <- high(high(A) * INVB) with a short product */
  for (po2 = K; (po2 & 1) == 0; po2 >>= 1);
  po2 = (po2 == 1);
  if (Fermat && po2)
    {
      mpz_set_ui (a[2 * K - 1], 0);
      if (K <= 4 * Fermat)
        {
          F_mul (t, a + K, invb, K, DEFAULT, Fermat, t + 2 * K);
          /* Put Q in T, as we still need high(A) later on */
          list_mod (t, t + K - 2, K, n);
        }
      else
        {
          F_mul (t, a + K, invb, K, DEFAULT, Fermat, t + 2 * K);
          list_mod (a + K, t + K - 2, K, n);
        }
    }
  else /* non-Fermat case */
    {
      list_mul_high (t, a + K, invb, K - 1);
      /* the high part of A * INVB is now in {t+K-2, K-1} */
      if (wrap)
	{
	  t2 = init_list2 (K - 1, mpz_sizeinbase (n, 2));
          ASSERT_ALWAYS(t2 != NULL);
	  list_mod (t2, t + K - 2, K - 1, n);
	}
      else /* we can store in high(A) which is no longer needed */
	list_mod (a + K, t + K - 2, K - 1, n);
    }

  /* the quotient Q = trunc(A / B) has degree K-2, i.e. K-1 terms */

  /* T <- low(Q * B) with a short product */
  mpz_set_ui (a[2 * K - 1], 0);
  if (Fermat && po2)
    {
      if (K <= 4 * Fermat)
        {
          /* Multiply without zero padding, result is (mod x^K - 1) */
          F_mul (t + K, t, b, K, NOPAD, Fermat, t + 2 * K);
          /* Take the leading monomial x^K of B into account */
          list_add (t, t + K, t, K);
          /* Subtract high(A) */
          list_sub(t, t, a + K, K);
        }
      else
        F_mul (t, a + K, b, K, DEFAULT, Fermat, t + 2 * K);
    }
  else /* non-Fermat case */
    {
      if (wrap)
        /* Q = {t2, K-1}, B = {b, K+1}
           We know that Q*B vanishes with the coefficients of degree
           K to 2K-2 of {A, 2K-1} */
        {
          unsigned int m;
          m = ks_wrapmul (t, K + 1, b, K + 1, t2, K - 1, n);
          clear_list (t2, K - 1);
          /* coefficients of degree m..2K-2 wrap around,
             i.e. were added to 0..2K-2-m */
          if (m < 2 * K - 1) /* otherwise product is exact */
            list_sub (t, t, a + m, 2 * K - 1 - m);
        }
      else
        list_mult_n (t, a + K, b, K);
    }

  /* now {t, K} contains the low K terms from Q*B */
  list_sub (a, a, t, K);
  list_mod (a, a, K, n);

  return 0;
}
