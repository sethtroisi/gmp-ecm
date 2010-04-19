/* Arithmetic on lists of residues modulo n.

  Copyright 2001, 2002, 2003, 2004, 2005 Paul Zimmermann and Alexander Kruppa.

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

#include <stdlib.h>
#include "ecm-impl.h"

#ifdef DEBUG
#define ASSERTD(x) assert(x)
#else
#define ASSERTD(x)
#endif

#if (MULT == KS)
  #define LIST_MULT_N kronecker_schonhage
  #define WRAP /* use wrap-around multiplication for low short product */
#elif (MULT == TOOM4)
  #define LIST_MULT_N toomcook4
#elif (MULT == TOOM3)
  #define LIST_MULT_N toomcook3
#elif (MULT == KARA)
  #define LIST_MULT_N karatsuba
#else
  #error "MULT is neither KS, TOOM4, nor TOOM3, nor KARA"
#endif

extern unsigned int Fermat;

/* returns a bound on the auxiliary memory needed by LIST_MULT_N */
int
list_mul_mem (unsigned int len)
{
  unsigned int mem;

  mem = 2 * len;
#if defined(TOOMCOOK3) || defined(TOOMCOOK4)
  while (len > 3)
    {
      mem += 2;
      len = (len + 2) / 3; /* ceil(len/3) */
    }
  mem += 4;
#endif
  return mem;
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

/* p[i] <- q[i] * r mod m */
void
list_mul_z (listz_t p, listz_t q, mpz_t r, unsigned int n, mpz_t m)
{
  unsigned int i;

  for (i = 0; i < n; i++)
    {
      mpz_mul (p[i], q[i], r);
      mpz_mod (p[i], p[i], m);
    }
}

/* p <- gcd(n, l[0]*l[1]*...*l[k-1],
   returns non-zero iff p is non trivial.
   Clobbers l[0] */
int
list_gcd (mpz_t p, listz_t l, unsigned int k, mpz_t n)
{
  unsigned int i;
  
  for (i = 1; i < k; i++)
    {
      mpz_mul (l[0], l[0], l[i]);
      mpz_mod (l[0], l[0], n);
    }
  mpz_gcd (p, l[0], n);

  return mpz_cmp_ui (p, 1);
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

#ifndef KS_MULTIPLY
/* puts in a[0]..a[K-1] the K low terms of the product 
   of b[0..K-1] and c[0..K-1].
   Assumes K >= 1, and a[0..2K-2] exist.
   Needs space for list_mul_mem(K) in t.
*/
static void
list_mul_low (listz_t a, listz_t b, listz_t c, unsigned int K, listz_t t,
	      mpz_t n)
{
  unsigned int p, q;

  ASSERT(K > 0);
  switch (K)
    {
    case 1:
      mpz_mul (a[0], b[0], c[0]);
      return;
    case 2:
      mpz_mul (a[0], b[0], c[0]);
      mpz_mul (a[1], b[0], c[1]);
      mpz_addmul (a[1], b[1], c[0]);
      return;
    case 3:
      karatsuba (a, b, c, 2, t);
      mpz_addmul (a[2], b[2], c[0]);
      mpz_addmul (a[2], b[0], c[2]);
      return;
    default:
      /* MULT is 2 for Karatsuba, 3 for Toom3, 4 for Toom4 */
      for (p = 1; MULT * p <= K; p *= MULT); /* p = greatest power of MULT <=K */
      p = (K / p) * p;
      ASSERTD(list_check(b,p,n) && list_check(c,p,n));
      LIST_MULT_N (a, b, c, p, t);
      if ((q = K - p))
        {
          list_mul_low (t, b + p, c, q, t + 2 * q - 1, n);
          list_add (a + p, a + p, t, q);
          list_mul_low (t, c + p, b, q, t + 2 * q - 1, n);
          list_add (a + p, a + p, t, q);
        }
    }
}
#endif

/* puts in a[K-1]..a[2K-2] the K high terms of the product 
   of b[0..K-1] and c[0..K-1].
   Assumes K >= 1, and a[0..2K-2] exist.
   Needs space for list_mul_mem(K) in t.
*/
void
list_mul_high (listz_t a, listz_t b, listz_t c, unsigned int K, listz_t t)
{
#ifdef KS_MULTIPLY /* ks is faster */
  LIST_MULT_N (a, b, c, K, t);
#else
  unsigned int p, q;

  ASSERT(K > 0);
  switch (K)
    {
    case 1:
      mpz_mul (a[0], b[0], c[0]);
      return;
      
    case 2:
      mpz_mul (a[2], b[1], c[1]);
      mpz_mul (a[1], b[1], c[0]);
      mpz_addmul (a[1], b[0], c[1]);
      return;

    case 3:
      karatsuba (a + 2, b + 1, c + 1, 2, t);
      mpz_addmul (a[2], b[0], c[2]);
      mpz_addmul (a[2], b[2], c[0]);
      return;

    default:
      /* MULT is 2 for Karatsuba, 3 for Toom3, 4 for Toom4 */
      for (p = 1; MULT * p <= K; p *= MULT);
      p = (K / p) * p;
      q = K - p;
      LIST_MULT_N (a + 2 * q, b + q, c + q, p, t);
      if (q)
        {
          list_mul_high (t, b + p, c, q, t + 2 * q - 1);
          list_add (a + K - 1, a + K - 1, t + q - 1, q);
          list_mul_high (t, c + p, b, q, t + 2 * q - 1);
          list_add (a + K - 1, a + K - 1, t + q - 1, q);
        }
    }
#endif
}

/* Puts in a[0..2K-2] the product of b[0..K-1] and c[0..K-1].
   The auxiliary memory M(K) necessary in T satisfies:
   M(1)=0, M(K) = max(3*l-1,2*l-2+M(l)) <= 2*K-1 where l = ceil(K/2).
   Assumes K >= 1.
*/
void
karatsuba (listz_t a, listz_t b, listz_t c, unsigned int K, listz_t t)
{
  if (K == 1)
    {
      mpz_mul (a[0], b[0], c[0]);
    }
  else if (K == 2) /* basic Karatsuba scheme */
    {
      mpz_add (t[0], b[0], b[1]); /* t0 = b_0 + b_1 */
      mpz_add (a[1], c[0], c[1]); /* a1 = c_0 + c_1 */
      mpz_mul (a[1], a[1], t[0]); /* a1 = b_0*c_0 + b_0*c_1 + b_1*c_0 + b_1*c_1 */
      mpz_mul (a[0], b[0], c[0]); /* a0 = b_0 * c_0 */
      mpz_mul (a[2], b[1], c[1]); /* a2 = b_1 * c_1 */
      mpz_sub (a[1], a[1], a[0]); /* a1 = b_0*c_1 + b_1*c_0 + b_1*c_1 */
      mpz_sub (a[1], a[1], a[2]); /* a1 = b_0*c_1 + b_1*c_0 */
    }
  else if (K == 3)
    {
      /* implement Weimerskirch/Paar trick in 6 muls and 13 adds
         http://www.crypto.ruhr-uni-bochum.de/Publikationen/texte/kaweb.pdf */
      /* diagonal terms */
      mpz_mul (a[0], b[0], c[0]);
      mpz_mul (a[2], b[1], c[1]);
      mpz_mul (a[4], b[2], c[2]);
      /* (0,1) rectangular term */
      mpz_add (t[0], b[0], b[1]);
      mpz_add (t[1], c[0], c[1]);
      mpz_mul (a[1], t[0], t[1]);
      mpz_sub (a[1], a[1], a[0]);
      mpz_sub (a[1], a[1], a[2]);
      /* (1,2) rectangular term */
      mpz_add (t[0], b[1], b[2]);
      mpz_add (t[1], c[1], c[2]);
      mpz_mul (a[3], t[0], t[1]);
      mpz_sub (a[3], a[3], a[2]);
      mpz_sub (a[3], a[3], a[4]);
      /* (0,2) rectangular term */
      mpz_add (t[0], b[0], b[2]);
      mpz_add (t[1], c[0], c[2]);
      mpz_mul (t[2], t[0], t[1]);
      mpz_sub (t[2], t[2], a[0]);
      mpz_sub (t[2], t[2], a[4]);
      mpz_add (a[2], a[2], t[2]);
    }
  else
    { 
      unsigned int i, k, l;
      listz_t z;

      k = K / 2;
      l = K - k;

      z = t + 2 * l - 1;

      /* improved code with 7*k-3 additions, 
         contributed by Philip McLaughlin <mpbjr@qwest.net> */
      for (i = 0; i < k; i++)
        {
          mpz_sub (z[i], b[i], b[l+i]);
          mpz_sub (a[i], c[i], c[l+i]);
        }

      if (l > k) /* case K odd */
        {
          mpz_set (z[k], b[k]);
          mpz_set (a[k], c[k]);
        }

      /* as b[0..l-1] + b[l..K-1] is stored in t[2l-1..3l-2], we need
         here at least 3l-1 entries in t */

      karatsuba (t, z, a, l, a + l); /* fills t[0..2l-2] */
       
      /* trick: save t[2l-2] in a[2l-1] to enable M(K) <= 2*K-1 */
      z = t + 2 * l - 2;
      mpz_set (a[2*l-1], t[2*l-2]);

      karatsuba (a, b, c, l, z); /* fill a[0..2l-2] */
      karatsuba (a + 2 * l, b + l, c + l, k, z); /* fills a[2l..2K-2] */

      mpz_set (t[2*l-2], a[2*l-1]); /* restore t[2*l-2] */
      mpz_set_ui (a[2*l-1], 0);

      /*
	      l          l-1     1    l          2k-1-l
        _________________________________________________
	|    a0    |     a1    |0|    a2    |     a3    |
        -------------------------------------------------
              l          l-1
        ________________________
	|    t0    |     t1    |
        ------------------------

	We want to replace [a1, a2] by [a1 + a0 + a2 - t0, a2 + a1 + a3 - t1]
	i.e. [a12 + a0 - t0, a12 + a3 - t1] where a12 = a1 + a2.
       */

      list_add (a + 2 * l, a + 2 * l, a + l, l-1); /* a[2l..3l-1] <- a1+a2 */
      if (k > 1)
        {
          list_add (a + l, a + 2 * l, a, l); /* a[l..2l-1] <- a0 + a1 + a2 */
          list_add (a + 2 * l, a + 2 * l, a + 3 * l, 2 * k - 1 - l);
        }
      else /* k=1, i.e. K=2 or K=3, and a2 has only one entry */
        {
          mpz_add (a[l], a[2*l], a[0]);
          if (K == 3)
            mpz_set (a[l+1], a[1]);
        }

      list_sub (a + l, a + l, t, 2 * l - 1);
    }
}

/* multiplies b[0]+...+b[k-1]*x^(k-1)+x^k by c[0]+...+c[l-1]*x^(l-1)+x^l
   and puts the results in a[0]+...+a[k+l-1]*x^(k+l-1)
   [the leading monomial x^(k+l) is implicit].
   If monic_b (resp. monic_c) is 0, don't consider x^k in b (resp. x^l in c).
   Assumes k = l or k = l+1.
   The auxiliary array t contains at least list_mul_mem(l) entries.
   a and t should not overlap.
*/
void
list_mul (listz_t a, listz_t b, unsigned int k, int monic_b,
          listz_t c, unsigned int l, int monic_c, listz_t t)
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
      if (monic_b && monic_c && l == k)
        {
          F_mul (a, b, c, l, MONIC, Fermat, t);
          monic_b = monic_c = 0;
        }
      else
        F_mul (a, b, c, l, DEFAULT, Fermat, t);
    }
  else
    LIST_MULT_N (a, b, c, l, t); /* set a[0]...a[2l-2] */

  if (k > l) /* multiply b[l]*x^l by c[0]+...+c[l-1]*x^(l-1) */
    {
      for (i = 0; i < l - 1; i++)
        mpz_addmul (a[l+i], b[l], c[i]);
      mpz_mul (a[2*l-1], b[l], c[l-1]);
    }

  /* deal with x^k and x^l */
  if (monic_b || monic_c)
    {
      mpz_set_ui (a[k + l - 1], 0);
      
      if (monic_b && monic_c) /* Single pass over a[] */
        {
          /* a += b * x^l + c * x^k, so a[i] += b[i-l]; a[i] += c[i-k] 
             if 0 <= i-l < k  or  0 <= i-k < l, respectively */
          if (k > l)
            mpz_add (a[l], a[l], b[0]);
          for (i = k; i < k + l; i++)
            {
              mpz_add (a[i], a[i], b[i-l]); /* i-l < k */
              mpz_add (a[i], a[i], c[i-k]); /* i-k < l */
            }
        }
      else if (monic_c) /* add b * x^l */
        list_add (a + l, a + l, b, k);

      else /* only monic_b, add x^k * c */
        list_add (a + k, a + k, c, l);
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

  for (i = k; (i & 1) == 0; i >>= 1);
  
  ASSERTD(list_check(b,k,n));
  ASSERTD(list_check(c,k,n));
  if (i == 1 && Fermat)
    F_mul (a, b, c, k, DEFAULT, Fermat, t);
  else
    LIST_MULT_N (a, b, c, k, t); /* set a[0]...a[2l-2] */

  list_mod (a2, a, 2 * k - 1, n);
}

/* puts in G[0]..G[k-1] the coefficients from (x-a[0])...(x-a[k-1])
   Warning: doesn't fill the coefficient 1 of G[k], which is implicit.
   Needs k + list_mul_mem(k/2) cells in T.
   G == a is permissible. T must not overlap with anything else
*/
void
PolyFromRoots (listz_t G, listz_t a, unsigned int k, listz_t T, mpz_t n)
{
  unsigned int l, m;

  ASSERT (T != G && T != a);

  if (k <= 1)
    {
      if (k == 1)
	{
#if NEGATED_ROOTS == 1
	  mpz_mod (G[0], a[0], n);
#else
	  mpz_neg (T[0], a[0]);
	  mpz_mod (G[0], T[0], n);
#endif
	}
      return;
    }

#if 0
  /* (x-a[0]) * (x-a[1]) = x^2 - (a[0]+a[1]) * x + a[0]*a[1]
     however we construct (x+a[0]) * (x+a[1]) instead, i.e. the
     polynomial with the opposite roots. This has no consequence if
     we do it for all polynomials: if F(x) and G(x) have a common root,
     then so do F(-x) and G(-x). This saves one negation.
     
     NOT ANY MORE! We added the mpz_neg()'s  so building the poly from
     its roots works as one would expect in the FastPM1Stage2 algorithm!
  */
  if (k == 2)
    {
      mpz_mul (T[0], a[0], a[1]);
      mpz_add (T[1], a[1], a[0]); /* mpz_add may allocate extra limb */
#if NEGATED_ROOTS == 0
      mpz_neg (T[1], T[1]);
#endif
      mpz_mod (G[1], T[1], n);
      mpz_mod (G[0], T[0], n);
      return;
    }
#endif
  
  m = k / 2; /* m == 1 */
  l = k - m; /* l == 2 */
  
  PolyFromRoots (G, a, l, T, n);
  PolyFromRoots (G + l, a + l, m, T, n);
#if 0
  {
    unsigned int i;
    outputf (OUTPUT_RESVERBOSE, "N=%Zd;", n);
    outputf (OUTPUT_RESVERBOSE, "(x^%u ", l);
    for (i = 0; i < l; i++)
      outputf (OUTPUT_RESVERBOSE, "+ (Mod(%Zd, N) * x^%u) ", G[i], i);
    outputf (OUTPUT_RESVERBOSE, ") * ( x^%u", m);
    for (i = 0; i < m; i++)
      outputf (OUTPUT_RESVERBOSE, "+ (Mod(%Zd, N) * x^%u) ", G[i + l], i);
    outputf (OUTPUT_RESVERBOSE, ")  ");
  }
#endif
  list_mul (T, G, l, 1, G + l, m, 1, T + k);
  list_mod (G, T, k, n);
#if 0
  {
    unsigned int i;
    outputf (OUTPUT_RESVERBOSE, " == ( x^%u ", k);
    for (i = 0; i < k; i++)
      outputf (OUTPUT_RESVERBOSE, "+ (Mod(%Zd, N) * x^%u) ", G[i], i);
    outputf (OUTPUT_RESVERBOSE, ")\n");
  } 
#endif
}

/* puts in G[0]..G[k-1] the coefficients from (x-a[0])...(x-a[k-1])
   Warning: doesn't fill the coefficient 1 of G[k], which is implicit.
   Needs k + list_mul_mem(k/2) cells in T.
   The product tree is stored in:
   G[0..k-1]       (degree k)
   Tree[0][0..k-1] (degree k/2)
   Tree[1][0..k-1] (degree k/4), ...,
   Tree[lgk-1][0..k-1] (degree 1)
   (then we should have initially Tree[lgk-1] = a).

   depth is the depth (0 at root), dolvl signals that only this level of
   the tree should be computed (-1: all levels)
*/
int
PolyFromRoots_Tree (listz_t G, listz_t a, unsigned int k, listz_t T, 
               int dolvl, mpz_t n, listz_t *Tree, FILE *TreeFile, 
               unsigned int sh)
{
  unsigned int l, m;
  listz_t H1, *NextTree;

  if (k <= 1)
    {
      if (k == 1)
        {
#if NEGATED_ROOTS == 1
	  mpz_mod (G[0], a[0], n);
#else
	  mpz_neg (a[0], a[0]);
	  mpz_mod (G[0], a[0], n);
          mpz_neg (a[0], a[0]);
#endif
        }
      return 0;
    }

  if (Tree == NULL)
    {
      H1 = G;
      NextTree = NULL;
    }
  else
    {
      H1 = Tree[0] + sh;
      NextTree = Tree + 1;
    }

#if 0
  /* (x-a[0]) * (x-a[1]) = x^2 - (a[0]+a[1]) * x + a[0]*a[1]
     however we construct (x+a[0]) * (x+a[1]) instead, i.e. the
     polynomial with the opposite roots. This has no consequence if
     we do it for all polynomials: if F(x) and G(x) have a common root,
     then so do F(-x) and G(-x). This saves one negation.
       
       NOT ANY MORE! We added the mpz_neg() below so building the poly from
       its roots works as one would expect in the FastPM1Stage2 algorithm!
  */
  if (k == 2)
    {
      mpz_set (H1[0], a[0]);
      mpz_set (H1[1], a[1]);
      mpz_mul (T[0], a[0], a[1]);
      mpz_add (G[1], a[1], a[0]);
#if NEGATED_ROOTS == 0
      mpz_neg (G[1], G[1]);
#endif
      mpz_mod (G[1], G[1], n);
      mpz_mod (G[0], T[0], n);
      return 0;
    }
#endif

  m = k / 2;
  l = k - m;
  
  if (dolvl < 0 || dolvl > 0)
    {
      PolyFromRoots_Tree (H1, a, l, T, dolvl - 1, n, NextTree, TreeFile, sh);
      PolyFromRoots_Tree (H1 + l, a + l, m, T, dolvl - 1, n, NextTree, 
                          TreeFile, sh + l);
    }
  if (dolvl < 0 || dolvl == 0)
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
      list_mul (T, H1, l, 1, H1 + l, m, 1, T + k);
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

#ifdef KS_MULTIPLY
      use_middle_product = 1;
#endif

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
          LIST_MULT_N (t, q + k, b, l, t + 2 * l - 1); /* t[0..2l-1] = Q1 * B0 */
          list_neg (t, t + l - 1, k, n);
      
          if (k > 1)
            {
              list_mul (t + k, q + k, l - 1, 1, b + l, k - 1, 1,
			t + k + K - 2); /* Q1 * B1 */
              list_sub (t + 1, t + 1, t + k, k - 1);
            }
        }
      list_mod (t, t, k, n); /* high(1-B*Q1) */

      ASSERTD(list_check(t,k,n) && list_check(q+l,k,n));
      if (po2)
        F_mul (t + k, t, q + l, k, DEFAULT, Fermat, t + 3 * k);
      else
        LIST_MULT_N (t + k, t, q + l, k, t + 3 * k - 1);
      list_mod (q, t + 2 * k - 1, k, n);
    }
}

/*
  divides a[0]+a[1]*x+...+a[2K-1]*x^(2K-1)
  By b[0]+b[1]*x+...+b[K-1]*x^(K-1)+x^K
  i.e. a polynomial of 2K coefficients divided by a monic polynomial
  with K+1 coefficients (b[K]=1 is implicit).
  Puts the quotient in q[0]+q[1]*x+...+q[K-1]*x^(K-1)
  and the remainder in a[0]+a[1]*x+...+a[K-1]*x^(K-1)
  Needs space for list_mul_mem(K) coefficients in t.
  If top is non-zero, a[0]..a[K-1] are reduced mod n.
*/
void
RecursiveDivision (listz_t q, listz_t a, listz_t b, unsigned int K,
                   listz_t t, mpz_t n, int top)
{
  if (K == 1) /* a0+a1*x = a1*(b0+x) + a0-a1*b0 */
    {
      mpz_mod (a[1], a[1], n);
      mpz_mul (q[0], a[1], b[0]);
      mpz_mod (q[0], q[0], n);
      mpz_sub (a[0], a[0], q[0]);
      if (top)
        mpz_mod (a[0], a[0], n);
      mpz_set (q[0], a[1]);
    }
  else
    {
      unsigned int k, l, i, po2;

      k = K / 2;
      l = K - k;
      for (po2 = K; (po2 && 1) == 0; po2 >>= 1);
      po2 = (po2 == 1);

      /* first perform a (2l) / l division */
      RecursiveDivision (q + k, a + 2 * k, b + k, l, t, n, 0);
      /* subtract q[k..k+l-1] * b[0..k-1] */
      ASSERTD(list_check(q+l,k,n) && list_check(b,k,n));
      if (po2 && Fermat)
        F_mul (t, q + l, b, k, DEFAULT, Fermat, t + K); /* sets t[0..2*k-2]*/
      else
        LIST_MULT_N (t, q + l, b, k, t + K - 1); /* sets t[0..2*k-2] */
      list_sub (a + l, a + l, t, 2 * k - 1);
      if (k < l) /* don't forget to subtract q[k] * b[0..k-1] */
        {
	  for (i=0; i<k; i++)
	    {
	      mpz_mul (t[0], q[k], b[i]); /* TODO: need to reduce t[0]? */
	      mpz_sub (a[k+i], a[k+i], t[0]);
	    }
        }
      /* remainder is in a[0..K+k-1] */

      /* then perform a (2k) / k division */
      RecursiveDivision (q, a + l, b + l, k, t, n, 0);
      /* subtract q[0..k-1] * b[0..l-1] */
      ASSERTD(list_check(q,k,n) && list_check(b,k,n));
      if (po2 && Fermat)
        F_mul (t, q, b, k, DEFAULT, Fermat, t + K);
      else
        LIST_MULT_N (t, q, b, k, t + K - 1);
      list_sub (a, a, t, 2 * k - 1);
      if (k < l) /* don't forget to subtract q[0..k-1] * b[k] */
        {
          for (i=0; i<k; i++)
            {
              mpz_mul (t[0], q[i], b[k]); /* TODO: need to reduce t[0]? */
              mpz_sub (a[k+i], a[k+i], t[0]);
            }
        }

      /* normalizes the remainder wrt n */
      if (top)
        list_mod (a, a, K, n);
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
#ifdef WRAP
  wrap = ks_wrapmul_m (K + 1, K + 1, n) <= 2 * K - 1 + list_mul_mem (K);
#else
  wrap = 0;
#endif

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
      list_mul_high (t, a + K, invb, K - 1, t + 2 * K - 3);
      /* the high part of A * INVB is now in {t+K-2, K-1} */
      if (wrap)
	{
	  MEMORY_TAG;
	  t2 = init_list2 (K - 1, mpz_sizeinbase (n, 2));
	  MEMORY_UNTAG;
	  if (t2 == NULL)
	    {
	      fprintf (ECM_STDERR, "Error, not enough memory\n");
	      return ECM_ERROR;
	    }
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
#ifdef KS_MULTIPLY /* ks is faster */
      if (wrap)
        /* Q = {t2, K-1}, B = {b, K+1}
           We know that Q*B vanishes with the coefficients of degree
           K to 2K-2 of {A, 2K-1} */
        {
          unsigned int m;
          m = ks_wrapmul (t, K + 1, b, K + 1, t2, K - 1, n);
          clear_list (t2, K - 1);
          /* coefficients of degree m..2K-2 wrap around,
             i.e. were subtracted to 0..2K-2-m */
          if (m < 2 * K - 1) /* otherwise product is exact */
            list_add (t, t, a + m, 2 * K - 1 - m);
        }
      else
        LIST_MULT_N (t, a + K, b, K, t + 2 * K - 1);
#else
      list_mul_low (t, a + K, b, K, t + 2 * K - 1, n);
#endif
    }

  /* now {t, K} contains the low K terms from Q*B */
  list_sub (a, a, t, K);
  list_mod (a, a, K, n);

  return 0;
}

/* Puts in inv[0..l-1] the inverses of a[0..l-1] (mod n), using 3*(l-1) 
   multiplies and one gcdext.
   Returns 1 if a factor was found (stored in t), 0 otherwise.
*/
int
list_invert (listz_t inv, listz_t a, unsigned long l, mpz_t t, mpmod_t modulus)
{
  unsigned long i;
  
  if (l == 0)
    return 0;
  
  mpz_set (inv[0], a[0]);
  
  for (i = 1; i < l; i++)
    {
      mpz_mul (t, inv[i-1], a[i]);
      mpz_mod (inv[i], t, modulus->orig_modulus); /* inv[i] = a[0]*...*a[i] */
    }
  
  mpz_gcdext (t, inv[l-1], NULL, inv[l-1], modulus->orig_modulus);
  
  if (mpz_cmp_ui (t, 1) != 0)
    return 1;
  
  for (i = l-1; i > 0; i--)
    {
      mpz_mul (t, inv[i], inv[i-1]); /* t = (a[0]*...*a[i])^(-1) * (a[0]*...*a[i-1]) = a[i]^(-1) */
      mpz_mul (inv[i-1], inv[i], a[i]); /* inv[i-1] = (a[0]*...*a[i])^(-1) * a[i] = (a[0]*...*a[i-1])^(-1) */
      mpz_mod (inv[i-1], inv[i-1], modulus->orig_modulus);
      mpz_mod (inv[i], t, modulus->orig_modulus);
    }
  
  return 0;
}
