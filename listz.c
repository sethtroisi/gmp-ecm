/* Arithmetic on lists of integers.

  Copyright 2001, 2002, 2003 Paul Zimmermann and Alexander Kruppa.

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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "ecm.h"

#if (MULT == TOOM4)
#define LIST_MULT_N toomcook4
#elif (MULT == TOOM3)
#define LIST_MULT_N toomcook3
#elif (MULT == KARA)
#define LIST_MULT_N karatsuba
#else
#error "MULT is neither TOOM4, nor TOOM3, nor KARA"
#endif

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

/* creates a list of n integers */
listz_t
init_list (unsigned int n)
{
  listz_t p;
  unsigned int i;

  p = (mpz_t*) malloc (n * sizeof(mpz_t));
  if (p == NULL)
    {
      fprintf (stderr, "Error: not enough memory\n");
      exit (EXIT_FAILURE);
    }
  for (i = 0; i < n; i++)
    mpz_init (p[i]);
  return p;
}

/* clears a list of n integers */
void
clear_list (listz_t p, unsigned int n)
{
  unsigned int i;

  for (i = 0; i < n; i++)
    mpz_clear (p[i]);
  free (p);
}

#define POLYFORM

#ifdef DEBUG
/* prints a list of n coefficients */
void
print_list (listz_t p, unsigned int n)
{
  unsigned int i;

  for (i = 0; i < n; i++)
    {
#ifdef POLYFORM
      if (i > 0 && mpz_cmp_ui (p[i], 0) >= 0)
        printf ("+");
#endif
      mpz_out_str (stdout, 10, p[i]);
#ifdef POLYFORM
      printf ("*x^%u", i);
#else
      printf (" ");
#endif
    }
  putchar ('\n');
}

#define SIZ(x) ((x)->_mp_size)
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define ABSIZ(x) ABS (SIZ (x))

unsigned int
print_max_size (listz_t p, unsigned int n)
{
  unsigned int i, s = 0;

  for (i = 0; i < n; i++)
    if (ABSIZ(p[i]) > s)
      s = ABSIZ(p[i]);
  return s;
}
#endif

/* p <- q */
void
list_set (listz_t p, listz_t q, unsigned int n)
{
  unsigned int i;

  for (i = 0; i < n; i++)
    mpz_set (p[i], q[i]);
}

#if 0
/* p[0] <-> p[n], p[1] <-> p[n-1], ... */
void
list_revert (listz_t p, unsigned int n)
{
  unsigned int i;
  for (i = 0; i < n - i; i++)
    mpz_swap (p[i], p[n - i]);
}
#endif

/* p <- -q */
void
list_neg (listz_t p, listz_t q, unsigned int n)
{
  unsigned int i;

  for (i = 0; i < n; i++)
    mpz_neg (p[i], q[i]);
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
list_add (listz_t p, listz_t q, listz_t r, unsigned int n)
{
  unsigned int i;

  for (i = 0; i < n; i++)
    mpz_add (p[i], q[i], r[i]);
}

/* p <- q - r */
void
list_sub (listz_t p, listz_t q, listz_t r, unsigned int n)
{
  unsigned int i;

  for (i = 0; i < n; i++)
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
   returns non-zero iff p is non trivial */
int
list_gcd (mpz_t p, listz_t l, unsigned int k, mpz_t n)
{
  unsigned int i;
  
  for (i=1; i<k; i++)
    {
      mpz_mul (l[0], l[0], l[i]);
      mpz_mod (l[0], l[0], n);
    }
  mpz_gcd (p, l[0], n);

  return mpz_cmp_ui (p, 1);
}

/* p <- 0 */
void
list_zero (listz_t p, unsigned int n)
{
  unsigned int i;

  for (i = 0; i < n; i++)
    mpz_set_ui (p[i], 0);
}

/* returns non-zero iff p is 0 */
int
list_zerop (listz_t p, unsigned int n)
{
  unsigned int i;
  int iszero = 1;

  for (i = 0; iszero && (i < n); i++)
    iszero = iszero && (mpz_cmp_ui (p[i], 0) == 0);

  return iszero;
}

/* puts in a[0]..a[K-1] the K low terms of the product 
   of b[0..K-1] and c[0..K-1].
   Returns the number of scalar multiplies.
   Assumes K >= 1, and a[0..2K-2] exist.
   Needs space for list_mul_mem(K) in t.
*/
static int
list_mul_low (listz_t a, listz_t b, listz_t c, unsigned int K, listz_t t)
{
  unsigned int p, q, muls;

  switch (K)
    {
    case 1:
      mpz_mul (a[0], b[0], c[0]);
      return 1;
    case 2:
      mpz_mul (a[0], b[0], c[0]);
      mpz_mul (a[1], b[0], c[1]);
      mpz_addmul (a[1], b[1], c[0]);
      return 3;
    case 3:
      karatsuba (a, b, c, 2, t);
      mpz_addmul (a[2], b[2], c[0]);
      mpz_addmul (a[2], b[0], c[2]);
      return 5;
    default:
      /* MULT is 2 for Karatsuba, 3 for Toom3, 4 for Toom4 */
      for (p = 1; MULT * p <= K; p *= MULT);
      p = (K / p) * p;
      muls = LIST_MULT_N (a, b, c, p, t);
      if ((q = K - p))
        {
          muls += list_mul_low (t, b + p, c, q, t + 2 * q - 1);
          list_add (a + p, a + p, t, q);
          muls += list_mul_low (t, c + p, b, q, t + 2 * q - 1);
          list_add (a + p, a + p, t, q);
        }
      return muls;
    }
}

/* puts in a[K-1]..a[2K-2] the K high terms of the product 
   of b[0..K-1] and c[0..K-1].
   Returns the number of scalar multiplies.
   Assumes K >= 1, and a[0..2K-2] exist.
   Needs space for list_mul_mem(K) in t.
*/
int
list_mul_high (listz_t a, listz_t b, listz_t c, unsigned int K, listz_t t)
{
  unsigned int p, q, muls;

  switch (K)
    {
    case 1:
      mpz_mul (a[0], b[0], c[0]);
      return 1;
      
    case 2:
      mpz_mul (a[2], b[1], c[1]);
      mpz_mul (a[1], b[1], c[0]);
      mpz_addmul (a[1], b[0], c[1]);
      return 3;

    case 3:
      karatsuba (a + 2, b + 1, c + 1, 2, t);
      mpz_addmul (a[2], b[0], c[2]);
      mpz_addmul (a[2], b[2], c[0]);
      return 5;

    default:
      /* MULT is 2 for Karatsuba, 3 for Toom3, 4 for Toom4 */
      for (p = 1; MULT * p <= K; p *= MULT);
      p = (K / p) * p;
      q = K - p;
      muls = LIST_MULT_N (a + 2 * q, b + q, c + q, p, t);
      if (q)
        {
          muls += list_mul_high (t, b + p, c, q, t + 2 * q - 1);
          list_add (a + K - 1, a + K - 1, t + q - 1, q);
          muls += list_mul_high (t, c + p, b, q, t + 2 * q - 1);
          list_add (a + K - 1, a + K - 1, t + q - 1, q);
        }
      return muls;
    }
}

/* Puts in a[0..2K-2] the product of b[0..K-1] and c[0..K-1].
   The auxiliary memory M(K) necessary in T satisfies:
   M(1)=0, M(K) = max(3*l-1,2*l-2+M(l)) <= 2*K-1 where l = ceil(K/2).
   Returns the number of scalar multiplies.
   Assumes K >= 1.
*/
int
karatsuba (listz_t a, listz_t b, listz_t c, unsigned int K, listz_t t)
{
  if (K == 1)
    {
      mpz_mul (a[0], b[0], c[0]);
      return 1;
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
      return 3;
    }
  else
    { 
      unsigned int i, k, l, muls;
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

      muls = karatsuba (t, z, a, l, a + l); /* fills t[0..2l-2] */
       
      /* trick: save t[2l-2] in a[2l-1] to enable M(K) <= 2*K-1 */
      z = t + 2 * l - 2;
      mpz_set (a[2*l-1], t[2*l-2]);

      muls += karatsuba (a, b, c, l, z); /* fill a[0..2l-2] */
      muls += karatsuba (a + 2 * l, b + l, c + l, k, z); /* fills a[2l..2K-2] */

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

      return muls;
    }
}

/* multiplies b[0]+...+b[k-1]*x^(k-1)+x^k by c[0]+...+c[l-1]*x^(l-1)+x^l
   and puts the results in a[0]+...+a[k+l-1]*x^(k+l-1)
   [the leading monomial x^(k+l) is implicit].
   If monic_b (resp. monic_c) is 0, don't consider x^k in b (resp. x^l in c).
   Assumes k = l or k = l+1.
   The auxiliary array t contains at least list_mul_mem(l) entries.
   a and t should not overlap.
   Return the number of scalar multiplies.
*/
int
list_mul (listz_t a, listz_t b, unsigned int k, int monic_b,
          listz_t c, unsigned int l, int monic_c,
          listz_t t)
{
  unsigned int i, muls;

  assert (k >= l);
  muls = LIST_MULT_N (a, b, c, l, t); /* set a[0]...a[2l-2] */

  if (k > l) /* multiply b[l]*x^l by c[0]+...+c[l-1]*x^(l-1) */
    {
      for (i = 0; i < l - 1; i++)
        mpz_addmul (a[l+i], b[l], c[i]);
      mpz_mul (a[2*l-1], b[l], c[l-1]);
      muls += l;
    }

  /* deal with x^k and x^l */
  if (monic_b || monic_c)
    mpz_set_ui (a[k + l - 1], 0);

  if (monic_c) /* add b * x^l */
    list_add (a + l, a + l, b, k);

  if (monic_b) /* add x^k * c */
    list_add (a + k, a + k, c, l);
  
  return muls;
}

/*
  Multiplies b[0..k-1] by c[0..k-1], and stores the result in a[0..2k-2].
  (Here, there is no implicit monic leading monomial.)
  Requires at least list_mul_mem(k) cells in t.
  Return the number of scalar multiplies.
 */
int
list_mulmod (listz_t a, listz_t b, listz_t c, unsigned int k, listz_t t,
             mpz_t n)
{
  int muls;

  muls = LIST_MULT_N (a, b, c, k, t);
  list_mod (a, a, 2*k - 1, n);

  return muls;
}

/*
  Multiplies b[0..k-1] by c[0..k-1], stores the result in a[0..2k-2],
  and stores the reduced product in a2[0..2k-2].
  (Here, there is no implicit monic leading monomial.)
  Requires at least list_mul_mem(k) cells in t.
  Return the number of scalar multiplies.
 */
int
list_mulmod2 (listz_t a2, listz_t a, listz_t b, listz_t c, unsigned int k,
              listz_t t, mpz_t n)
{
  int muls;

  muls = LIST_MULT_N (a, b, c, k, t);
  list_mod (a2, a, 2 * k - 1, n);

  return muls;
}

/* puts in G[0]..G[k-1] the coefficients from (x-a[0])...(x-a[k-1])
   Warning: doesn't fill the coefficient 1 of G[k], which is implicit.
   Needs k + list_mul_mem(k/2) cells in T.
   Return the number of multiplies.
   If Tree <> NULL, the product tree is stored in:
   G[0..k-1]       (degree k)
   Tree[0][0..k-1] (degree k/2)
   Tree[1][0..k-1] (degree k/4), ...,
   Tree[lgk-1][0..k-1] (degree 1)
   (then we should have initially Tree[lgk-1] = a).
*/
int
PolyFromRoots (listz_t G, listz_t a, unsigned int k, listz_t T, int verbose,
               mpz_t n, char F, listz_t *Tree, unsigned int sh)
{
  unsigned int l, m, st;
  unsigned long muls;
  listz_t H1, *NextTree;

   if (k <= 1)
     {
       /* if Tree=NULL, then G=a and nothing to do */
       if (Tree != NULL)
         mpz_set (G[0], a[0]);
       return 0;
     }

   st = cputime ();

   m = k / 2;
   l = k - m;

   if ((Tree == NULL) || verbose) /* top call */
     {
       H1 = (Tree == NULL) ? G : Tree[0]; /* target for rec. calls */
       NextTree = Tree;
     }
   else
     {
       H1 = Tree[1] + sh;
       NextTree = Tree + 1;
     }

    /* (x-a[0]) * (x-a[1]) = x^2 - (a[0]+a[1]) * x + a[0]*a[1]
       however we construct (x+a[0]) * (x+a[1]) instead, i.e. the
       polynomial with the opposite roots. This has no consequence if
       we do it for all polynomials: if F(x) and G(x) have a common root,
       then so do F(-x) and G(-x). This saves one negation.
     */
   if (k == 2)
     {
       if (Tree != NULL)
         {
           mpz_set (H1[0], a[0]);
           mpz_set (H1[1], a[1]);
         }
       mpz_mul (T[0], a[0], a[1]);
       mpz_add (G[1], a[1], a[0]);
       mpz_mod (G[1], G[1], n);
       mpz_mod (G[0], T[0], n);
       return 1;
     }

   muls = PolyFromRoots (H1, a, l, T, 0, n, F, NextTree, sh);
   muls += PolyFromRoots (H1 + l, a + l, m, T, 0, n, F, NextTree, sh + l);
   muls += list_mul (T, H1, l, 1, H1 + l, m, 1, T + k);
   list_mod (G, T, k, n);
   
   if (verbose >= 2)
     printf ("Building %c from its roots took %ums and %lu muls\n", F,
	      cputime() - st, muls);

   return muls;
}

/* puts in q[0..K-1] the quotient of x^(2K-1) by B
   where B = b[0]+b[1]*x+...+b[K-1]*x^(K-1) with b[K-1]=1.
   Return the number of scalar multiplies performed.
*/
int
PolyInvert (listz_t q, listz_t b, unsigned int K, listz_t t, mpz_t n)
{
  if (K == 1)
    {
      mpz_set_ui (q[0], 1);
      return 0;
    }
  else
    {
      int k, l, muls;

      k = K / 2;
      l = K - k;

      muls = PolyInvert (q + k, b + k, l, t, n); /* Q1 = {q+k, l} */

      muls += LIST_MULT_N (t, q + k, b, l, t + 2 * l - 1); /* Q1 * B0 */
      list_neg (t, t + l - 1, k);
      
      if (k > 1)
        {
          muls += list_mul (t + k, q + k, l - 1, 1, b + l, k - 1, 1,
                            t + k + K - 2); /* Q1 * B1 */
          list_sub (t + 1, t + 1, t + k, k - 1);
        }
      list_mod (t, t, k, n); /* high(1-B*Q1) */
      muls += LIST_MULT_N (t + k, t, q + l, k, t + 3 * k - 1);
      list_mod (q, t + 2 * k - 1, k, n);

      return muls;
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
  Return the number of scalar multiplies performed.
  If top is non-zero, a[0]..a[K-1] are reduced mod n.
*/
int
RecursiveDivision (listz_t q, listz_t a, listz_t b, unsigned int K,
                   listz_t t, mpz_t n, int top)
{
  if (K == 1) /* a0+a1*x = a1*(b0+x) + a0-a1*b0 */
    {
      mpz_mod (a[1], a[1], n);
      mpz_mulmod (q[0], a[1], b[0], n);
      mpz_sub (a[0], a[0], q[0]);
      if (top)
        mpz_mod (a[0], a[0], n);
      mpz_set (q[0], a[1]);
      return 1;
    }
  else
    {
      unsigned int k, l, i, muls;

      k = K / 2;
      l = K - k;

      /* first perform a (2l) / l division */
      muls = RecursiveDivision (q + k, a + 2 * k, b + k, l, t, n, 0);
      /* subtract q[k..k+l-1] * b[0..k-1] */
      muls += LIST_MULT_N (t, q + l, b, k, t + K - 1); /* sets t[0..2*k-2] */
      list_sub (a + l, a + l, t, 2 * k - 1);
      if (k < l) /* don't forget to subtract q[k] * b[0..k-1] */
        {
	  for (i=0; i<k; i++)
	    {
	      mpz_mul (t[0], q[k], b[i]); /* TODO: need to reduce t[0]? */
	      mpz_sub (a[k+i], a[k+i], t[0]);
	    }
          muls += k;
        }
      /* remainder is in a[0..K+k-1] */

      /* then perform a (2k) / k division */
      muls += RecursiveDivision (q, a + l, b + l, k, t, n, 0);
      /* subtract q[0..k-1] * b[0..l-1] */
      muls += LIST_MULT_N (t, q, b, k, t + K - 1);
      list_sub (a, a, t, 2 * k - 1);
      if (k < l) /* don't forget to subtract q[0..k-1] * b[k] */
        {
          for (i=0; i<k; i++)
            {
              mpz_mul (t[0], q[i], b[k]); /* TODO: need to reduce t[0]? */
              mpz_sub (a[k+i], a[k+i], t[0]);
            }
          muls += k;
        }

      /* normalizes the remainder wrt n */
      if (top)
        list_mod (a, a, K, n);

      return muls;
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
*/
int
PrerevertDivision (listz_t a, listz_t b, listz_t invb,
                   unsigned int K, listz_t t, mpz_t n)
{
  int muls;

  /* Q <- high(A * INVB) with a short product */
  muls = list_mul_high (t, a + K, invb, K - 1, t + 2 * K - 3);
  list_mod (a + K, t + K - 2, K - 1, n);

  /* the quotient Q has degree K-2, i.e. K-1 terms */

  /* T <- low(Q * B) with a short product */
  mpz_set_ui (a[2 * K - 1], 0);
  muls += list_mul_low (t, a + K, b, K, t + 2 * K - 1);

  /* now {t, K} contains the low K terms from Q*B */
  list_sub (a, a, t, K);
  list_mod (a, a, K, n);

  return muls;
}

/* a <- a mod b, with quotient q stored in a[d]*x + a[d-1],
   a has degree d, b has degree (d-1). */
int
list_mod1 (mpz_t p, listz_t a, listz_t b, unsigned int d, mpz_t n, listz_t t)
{
  unsigned int i;

  mpz_gcdext (p, t[0], NULL, b[d-1], n);
  if (mpz_cmp_ui (p, 1) != 0)
    return 1;
  /* t[0] = 1/b[d-1] mod n */

  /* subtract a[d]*t[0]*x*b from a */
  mpz_mulmod (a[d], a[d], t[0], n);
  for (i=0; i<d-1; i++)
    {
      mpz_mulmod (t[1], b[i], a[d], n);
      mpz_sub (a[i+1], a[i+1], t[1]);
    }
  /* now deg(a) = d-1 */
  
  /* subtract a[d-1]*t[0]*b from a */
  mpz_mulmod (a[d-1], a[d-1], t[0], n);
  for (i=0; i<d-1; i++)
    {
      mpz_mulmod (t[1], b[i], a[d-1], n);
      mpz_sub (a[i], a[i], t[1]);
    }
  /* now deg(a) = d-2, q = a[d]*x + a[d-1] */

  return 0;
}

/* a[0..d] <- a[0..d] - (q[1]*x + q[0]) * b[0..d-1] */
void
poly_submul2 (listz_t a, listz_t q, listz_t b, unsigned int d, mpz_t n, mpz_t t)
{
  unsigned int i;

  mpz_set_ui (a[d+1], 0);
  mpz_set_ui (a[d], 0);
  for (i=0; i<d; i++)
    {
      mpz_mulmod (t, b[i], q[0], n);
      mpz_sub (a[i], a[i], t);
      mpz_mulmod (t, b[i], q[1], n);
      mpz_sub (a[i+1], a[i+1], t);
    }
}

/* Puts in inv[0..l-1] the inverses of a[0..l-1] (mod n), using 3*(l-1) 
   multiplies and one gcdext.
   Returns 1 if a factor was found (stored in t), 0 otherwise.
*/
int
list_invert (listz_t inv, listz_t a, unsigned int l, mpz_t t, mpz_t n)
{
  unsigned int i;
  
  if (l == 0)
    return 0;
  
  mpz_set (inv[0], a[0]);
  
  for (i = 1; i < l; i++)
    {
      mpz_mul (t, inv[i-1], a[i]);
      mpz_mod (inv[i], t, n); /* inv[i] = a[0]*...*a[i] */
    }
  
  mpz_gcdext (t, inv[l-1], NULL, inv[l-1], n);
  
  if (mpz_cmp_ui (t, 1) != 0)
    return 1;
  
  for (i = l-1; i > 0; i--)
    {
      mpz_mul (t, inv[i], inv[i-1]); /* t = (a[0]*...*a[i])^(-1) * (a[0]*...*a[i-1]) = a[i]^(-1) */
      mpz_mul (inv[i-1], inv[i], a[i]); /* inv[i-1] = (a[0]*...*a[i])^(-1) * a[i] = (a[0]*...*a[i-1])^(-1) */
      mpz_mod (inv[i-1], inv[i-1], n);
      mpz_mod (inv[i], t, n);
    }
  
  return 0;
}
