/* Polynomial arithmetic.

  Copyright (C) 2001 Paul Zimmermann,
  LORIA/INRIA Lorraine, zimmerma@loria.fr
  See http://www.loria.fr/~zimmerma/records/ecmnet.html

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

/* creates a (monic) polynomial of degree n */
mpz_t*
init_poly (unsigned int n)
{
  mpz_t *p;
  int i;

  p = (mpz_t*) malloc (n * sizeof(mpz_t));
  for (i = 0; i < n; i++)
    mpz_init (p[i]);
  return p;
}

/* clears a polynomial of degree n */
void
clear_poly (mpz_t *p, unsigned int n)
{
  unsigned int i;

  for (i = 0; i < n; i++)
    mpz_clear (p[i]);
  free (p);
}

/* prints a list of n coefficients */
void
print_list (mpz_t *p, unsigned int n)
{
  unsigned int i;

  for (i = 0; i < n; i++)
    {
      mpz_out_str (stdout, 10, p[i]);
      printf (" ");
    }
  putchar ('\n');
}

/* prints a polynomial of degree n */
void
print_poly (mpz_t *p, unsigned int n, int monic)
{
  unsigned int i;

  for (i = 0; i < n; i++)
    {
      if ((i > 0) && (mpz_cmp_ui(p[i], 0) >= 0))
	printf ("+");
      mpz_out_str (stdout, 10, p[i]);
      printf ("*x^%u", i);
    }
  if (monic)
    printf ("+x^%u", n);
  putchar ('\n');
}

void
polyset (mpz_t *p, mpz_t *q, unsigned int n)
{
  unsigned int i;

  for (i = 0; i < n; i++)
    mpz_set (p[i], q[i]);
}

void
polymod (mpz_t *p, mpz_t *q, unsigned int n, mpz_t mod)
{
  unsigned int i;

  for (i = 0; i < n; i++)
    mpz_mod (p[i], q[i], mod);
}

void
polyadd (mpz_t *p, mpz_t *q, mpz_t *r, unsigned int n)
{
  unsigned int i;

  for (i = 0; i < n; i++)
    mpz_add (p[i], q[i], r[i]);
}

void
polysub (mpz_t *p, mpz_t *q, mpz_t *r, unsigned int n)
{
  unsigned int i;

  for (i = 0; i < n; i++)
    mpz_sub (p[i], q[i], r[i]);
}

void
poly_zero (mpz_t *p, unsigned int n)
{
  unsigned int i;

  for (i = 0; i < n; i++)
    mpz_set_ui (p[i], 0);
}

int
poly_iszero (mpz_t *p, unsigned int n)
{
  unsigned int i;
  int iszero = 1;

  for (i = 0; iszero && (i < n); i++)
    iszero = iszero && (mpz_cmp_ui (p[i], 0) == 0);

  return iszero;
}

/* Puts in a[0..2K-2] the product of b[0..K-1] and b[0..K-1].
   The auxiliary memory M(K) necessary in T satisfies:
   M(1)=0, M(K) = max(3*l-1,2*l-2+M(l)) <= 2*K-1 where l = ceil(K/2)
*/
void
karatsuba (mpz_t *a, mpz_t *b, mpz_t *c, unsigned int K, mpz_t *t)
{

   if (K == 1)
     mpz_mul (a[0], b[0], c[0]);
   else
     { 
       unsigned int i, k, l;
       mpz_t *z;
       
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

       polyadd (a + 2 * l, a + 2 * l, a + l, l-1); /* a[2l..3l-1] <- a1 + a2 */
       if (k > 1)
	 {
	   polyadd (a + l, a + 2 * l, a, l); /* a[l..2l-1] <- a0 + a1 + a2 */
	   polyadd (a + 2 * l, a + 2 * l, a + 3 * l, 2 * k - 1 - l);
	 }
       else /* k=1, i.e. K=2 or K=3, and a2 has only one entry */
	 {
	   mpz_add (a[l], a[2*l], a[0]);
	   if (K==3) mpz_set (a[l+1], a[1]);
	 }

       polysub (a + l, a + l, t, 2 * l - 1);
     }
}

/* multiplies b[0]+...+b[k-1]*x^(k-1)+x^k by c[0]+...+c[l-1]*x^(l-1)+x^l
   and puts the results in a[0]+...+a[k+l-1]*x^(k+l-1)
   [the leading monomial x^(k+l) is implicit].
   Assumes k >= l.
   The auxiliary array t contains at least k+l entries.
   If monic is zero, do not consider the monomials x^k from b and x^l from c.
   a and t should not overlap.
*/
void
polymul (mpz_t *a, mpz_t *b, unsigned int k, mpz_t *c, unsigned int l,
	 mpz_t *t, int monic)
{
  unsigned int i;

  assert (k >= l);
  karatsuba (a, b, c, l, t); /* set a[0]...a[2l-2] */

  for (i = l; i + l <= k; i += l)
    {
      karatsuba (t, b + i, c, l, t + 2 * l - 1);
      /* a[0..i+l-2] are already set */
      polyadd (a + i, a + i, t, l - 1);
      polyset (a + i + l - 1, t + l - 1, l); 
    }

  /* last block may be incomplete */
  if (i < k)
    {
      unsigned int m = k - i; /* m < l */

      polymul (t, c, l, b + i, m, t + l + m - 1, 0); /* set t[0..m+l-2] */
      /* a[0..i+l-2] are already set */
      polyadd (a + i, a + i, t, l - 1);
      polyset (a + i + l - 1, t + l - 1, m);
    }

  if (monic != 0)
    {
      mpz_set_ui (a[k+l-1], 0);
      /* add b * x^l */
      polyadd (a + l, a + l, b, k);
      /* add x^k * c */
      polyadd (a + k, a + k, c, l);
    }
}

/*
  Multiplies b[0..k-1] by c[0..k-1], and stores the result in a[0..2k-2].
  (Here, there is no implicit monic leading monomial.)
 */
void
polymulmod (mpz_t *a, mpz_t *b, mpz_t *c, unsigned int k, mpz_t *t, mpz_t n)
{
  karatsuba (a, b, c, k, t);
  polymod (a, a, 2*k - 1, n);
}

/* puts in G the coefficients from (x-G[0])...(x-G[k-1])
   using 2*k cells in T */
void
buildG (mpz_t *G, unsigned int k, mpz_t *T, int verbose, mpz_t n, char F)
{
   unsigned int l, m, st;

   if (k <= 1)
     return;

    /* (x-G[0]) * (x-G[1]) = x^2 - (G[0]+G[1]) * x + G[0]*G[1]
       however we construct (x+G[0]) * (x+G[1]) instead, i.e. the
       polynomial with the opposite roots. This has no consequence if
       we do it for all polynomials: if F(x) and G(x) have a common root,
       then so do F(-x) and G(-x). This saves one negation.
     */
   if (k == 2)
     {
       mpz_mul (T[0], G[0], G[1]);
       mpz_add (G[1], G[1], G[0]);
       mpz_mod (G[1], G[1], n);
       mpz_mod (G[0], T[0], n);
       return;
     }

   st = cputime ();

   m = k / 2;
   l = k - m;

   buildG (G, l, T, 0, n, F);
   buildG (G + l, m, T, 0, n, F);
   polymul (T, G, l, G + l, m, T + k, 1);
   polymod (G, T, k, n);
   
   if (verbose)
     fprintf (stderr, "Building %c from its roots took %dms\n", F,
	      cputime() - st);
}

#define mpz_mulmod(a,b,c,n) \
        mpz_mul (a, b, c); \
        mpz_mod (a, a, n);

/*
  divides a[0]+a[1]*x+...+a[2K-1]*x^(2K-1)
  by b[0]+b[1]*x+...+b[L-1]*x^(L-1)+x^L
  puts the quotient in q[0]+q[1]*x+...+q[K-1]*x^(K-1)
  and the remainder in a[0]+a[1]*x+...+a[K-1]*x^(K-1)
  Needs space for 2K coefficients in t.
*/
void
RecursiveDivision (mpz_t *q, mpz_t *a, mpz_t *b, unsigned int K, mpz_t *t,
		   mpz_t n)
{
  if (K == 1) /* a0+a1*x = a1*(b0+x) + a0-a1*b0 */
    {
      mpz_mulmod (q[0], a[1], b[0], n);
      mpz_sub (a[0], a[0], q[0]);
      mpz_set (q[0], a[1]);
    }
  else
    {
      unsigned int k, l, i;

      k = K / 2;
      l = K - k;

      /* first perform a (2l) / l division */
      RecursiveDivision (q + k, a + 2 * k, b + k, l, t, n);
      /* subtract q[k..k+l-1] * b[0..k-1] */
      karatsuba (t, q + l, b, k, t + K - 1); /* sets t[0..2*k-2] */
      polysub (a + l, a + l, t, 2 * k - 1);
      if (k < l) /* don't forget to subtract q[k] * b[0..k-1] */
	  for (i=0; i<k; i++)
	    {
	      mpz_mul (t[0], q[k], b[i]);
	      mpz_sub (a[k+i], a[k+i], t[0]);
	    }
      /* remainder is in a[0..K+k-1] */

      /* then perform a (2k) / k division */
      RecursiveDivision (q, a + l, b + l, k, t, n);
      /* subtract q[0..k-1] * b[0..l-1] */
      karatsuba (t, q, b, k, t + K - 1);
      polysub (a, a, t, 2 * k - 1);
      if (k < l) /* don't forget to subtract q[0..k-1] * b[k] */
	for (i=0; i<k; i++)
	  {
	    mpz_mul (t[0], q[i], b[k]);
	    mpz_sub (a[k+i], a[k+i], t[0]);
	  }
    }
  /* normalizes the remainder wrt n */
  polymod (a, a, K, n);
}

/*
  divides a[0]+a[1]*x+...+a[2K-1]*x^(3K-1)
  by b[0]+b[1]*x+...+b[2K-1]*x^(2K-1)+x^(2K)
  puts the quotient in q[0]+q[1]*x+...+q[K-1]*x^(K-1)
  and the remainder in a[0]+a[1]*x+...+a[2K-1]*x^(2K-1)
  K must be a power of two.
  Needs space for 4K coefficients in t.
*/
void
Div3by2 (mpz_t *q, mpz_t *a, mpz_t *b, unsigned int K, 
			 mpz_t *t, mpz_t n)
{
  if (K==1) { 
    /* a0+a1*x+a2*x^2 = a2*(b0+b1*x+x^2) + (a1-a2*b1)*x + (a0-a2*b0) */
    mpz_mulmod (q[0], a[2], b[1], n);
    mpz_sub(a[1], a[1], q[0]);
    mpz_mulmod (q[0], a[2], b[0], n);
    mpz_sub(a[0], a[0], q[0]);
    mpz_set(q[0], a[2]);
  }
  else { 
    unsigned int i=0;

    RecursiveDivision (q, a+K, b+K, K, t, n);
    karatsuba (t, q, b, K, t + 2 * K); /* needs 2K memory in t+2K */
    for (i=0;i<=2*K-2;i++) {
      mpz_sub (a[i], a[i], t[i]);
      mpz_mod (a[i], a[i], n);
    }
    /* we could also have a special version of karatsuba that directly
       subtract the result to a */
  }
}

/* a <- a mod b, with quotient q stored in a[d]*x + a[d-1] */
int polymod1 (mpz_t p, mpz_t *a, mpz_t *b, unsigned int d, mpz_t n, mpz_t *t)
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
poly_submul2 (mpz_t *a, mpz_t *q, mpz_t *b, unsigned int d, mpz_t n, mpz_t t)
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

/* Input: a[0..d] is a polynomial of degree d,
          b[0..d-1] is a polynomial of degree d-1.
   Output: R = [[R11,R12],[R21,R22]] a 2x2 matrix,
           such that (a', b') = R * (a, b) 
           with deg(b') = floor(da/2) < deg(a').
           If k = deg(b) - deg(b'), then deg(R11)=k-2, deg(R12)=deg(R21)=k-1,
           and deg(R22)=k.
   Need space for 2 coefficients in t[].
   If flag is zero, does not compute R11, R12, R21, R22.
*/
int
hgcd_naive (mpz_t p, mpz_t *R11, mpz_t *R12, mpz_t *R21, mpz_t *R22,
            mpz_t *a, mpz_t *b, unsigned int d, mpz_t n, mpz_t *t, int flag)
{
  unsigned int m, k;

  m = d / 2;
  k = (d - 1) - m; /* degree excess */
  if (k == 0)
    {
      if (flag)
        mpz_set_ui (R22[0], 1);
    }
  else
    {
      unsigned int i, k;

      /* a <- a mod b */
      if (polymod1 (p, a, b, d, n, t))
        return 1;

      /* now deg(a) = d-2, q = a[d]*x + a[d-1] */

      if (hgcd_naive (p, R12, R11, R22, R21, b, a, d - 1, n, t, flag))
        return 1;

      if (flag)
        {

          /* deg(R12) = k-3, deg(R11)=deg(R22)=k-2, deg(R21)=k-1 */
      
          /* R <- R * Matrix([[0,1],[1,-q]])
             i.e. R11 <- R11, R12 <- R12 - R11 * q,
                  R21 <- R21, R22 <- R22 - R21 * q */

          poly_submul2 (R12, a + d - 1, R11, k - 1, n, t[0]);

          poly_submul2 (R22, a + d - 1, R21, k, n, t[0]);

        }
    }

  return 0;
}

#define HGCD_THRESHOLD 1

/* Input: a[0..d] is a polynomial of degree d,
          b[0..d-1] is a polynomial of degree d-1.
   Output: R = [[R11,R12],[R21,R22]] a 2x2 matrix,
           such that (a', b') = R * (a, b) 
           with deg(b') = floor(da/2) < deg(a').
           If k = deg(b) - deg(b'), then deg(R11)=k-2, deg(R12)=deg(R21)=k-1,
           and deg(R22)=k.
*/
int
hgcd (mpz_t p, mpz_t *R11, mpz_t *R12, mpz_t *R21, mpz_t *R22,
      mpz_t *a, mpz_t *b, unsigned int d, mpz_t n, mpz_t *t, int flag)
{
  if (d <= HGCD_THRESHOLD)
    return hgcd_naive (p, R11, R12, R21, R22, a, b, d, n, t, flag);
  else
    {
      unsigned int m, k, m2, d2, k2, k3;
      mpz_t *t2, *t4, *S11, *S12, *S21, *S22;

      m = d / 2; /* target degree */
      k = (d - 1) - m; /* degree excess */
      assert (k <= m);
      t2 = t + 2 * m;
      t4 = t + 4 * m;
      hgcd (p, R11, R12, R21, R22, a + m, b + m, d - m, n, t, 1);
      d2 = m + (d - m) / 2; /* new degree of b */
      k2 = (d - m - 1) - (d - m) / 2;
      /* R11 has degree k2-2, R12 and R21 degree k2-1, and R22 degree k2 */

      /* replace [ahi, alo] by [ahi, R11*alo + R12*blo] */
      polymul (t,  b, m, R12, k,     t2, 0); /* sets t[0..m+k-2] */
      polymul (t2, a, m, R11, k - 1, t2, 0);
      polyadd (t, t, t2, m + k - 2);
      polymul (t2, a, m, R21, k,     t4, 0);
      polyset (a, t, m);
      polyadd (a + m, a + m, t + m, k - 1);
      polymul (t,  b, m, R22, k + 1, t4, 0); /* sets t[0..m+k-1] */
      polyadd (b, t, t2, m);
      polyadd (t + m, t + m, t2 + m, k);
      polyadd (b + m, b + m, t + m, k);

      if (poly_iszero (b, d2 + 1))
        return 0;

      if (polymod1 (p, a, b, d2 + 1, n, t))
        return 1;

      /* now a has degree d2-1, with quotient q = a[d2+1]*x + a[d2] */
      /* replace R by matrix([[0,1],[1,-q]]) * R, i.e. 
         R11 <- R21,         R12 <- R22,
         R21 <- R11 - q R21, R22 <- R12 - q R22 */

      poly_submul2 (R11, a + d2, R21, k2, n, t[0]);
      poly_submul2 (R12, a + d2, R22, k2 + 1, n, t[0]);
      /* swap R11 and R21, R12 and R22 */
      
      m2 = m / 2;
      S11 = t;
      S12 = t + k3;
      hgcd (p, S11, S12, S21, S22, b + m2, a + m2, d2 - m2, n, t, flag);
    }

  return 0;
}

/* Input: a[0..d]   is a polynomial of degree d.
          b[0..d-1] is a polynomial of degree d-1.
   Output: g = gcd(a, b) is in a.
   Return value: degree(g), or -1 if a factor was found, in which case
                 this factor is stored in p.
*/
int
polygcd (mpz_t p, mpz_t *a, mpz_t *b, unsigned int d, mpz_t n, mpz_t *t)
{
  int k;

  if (poly_iszero(b, d)) /* gcd is a */
      return d;

  k = d / 2; /* target degree for hgcd */

  if (hgcd_naive (p, t, t, t, t, a, b, d, n, t, 0))
    return -1;

  /* now deg(a) = k+1, deg(b) = k */
  
  if (polymod1 (p, a, b, k+1, n, t))
    return -1;

  /* now deg(b) = k, deg(a) = k-1 */

  k = polygcd (p, b, a, k, n, t);
  
  if (k < 0)
    return k;

  polyset (a, b, k + 1);

  return k;
}
