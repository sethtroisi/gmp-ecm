/* Arithmetic on lists of integers.

  Copyright (C) 2001-2002 Paul Zimmermann,
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

//#define TOOMCOOK3
#define TOOMCOOK4

#if defined(TOOMCOOK4)
#define LIST_MULT_N toomcook4
#elif defined(TOOMCOOK3)
#define LIST_MULT_N toomcook3
#else
#define LIST_MULT_N karatsuba
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

/* prints a list of n coefficients */
void
print_list (listz_t p, unsigned int n)
{
  unsigned int i;

  for (i = 0; i < n; i++)
    {
      mpz_out_str (stdout, 10, p[i]);
      printf (" ");
    }
  putchar ('\n');
}

/* p <- q */
void
list_set (listz_t p, listz_t q, unsigned int n)
{
  unsigned int i;

  for (i = 0; i < n; i++)
    mpz_set (p[i], q[i]);
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

/* Puts in a[0..2K-2] the product of b[0..K-1] and b[0..K-1].
   The auxiliary memory M(K) necessary in T satisfies:
   M(1)=0, M(K) = max(3*l-1,2*l-2+M(l)) <= 2*K-1 where l = ceil(K/2)
*/
void
karatsuba (listz_t a, listz_t b, listz_t c, unsigned int K, listz_t t)
{

   if (K == 1)
     mpz_mul (a[0], b[0], c[0]);
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
	   if (K==3) mpz_set (a[l+1], a[1]);
	 }

       list_sub (a + l, a + l, t, 2 * l - 1);
     }
}

/* multiplies b[0]+...+b[k-1]*x^(k-1)+x^k by c[0]+...+c[l-1]*x^(l-1)+x^l
   and puts the results in a[0]+...+a[k+l-1]*x^(k+l-1)
   [the leading monomial x^(k+l) is implicit].
   Assumes k = l or k = l+1.
   The auxiliary array t contains at least list_mul_mem(l) entries.
   a and t should not overlap.
*/
void
list_mul (listz_t a, listz_t b, unsigned int k, listz_t c, unsigned int l,
          listz_t t)
{
  unsigned int i;

  assert (k >= l);
  LIST_MULT_N (a, b, c, l, t); /* set a[0]...a[2l-2] */

  if (k > l) /* multiply b[l]*x^l by c[0]+...+c[l-1]*x^(l-1) */
    {
      for (i=0; i<l-1; i++)
        {
          mpz_mul (a[2*l-1], b[l], c[i]);
          mpz_add (a[l+i], a[l+i], a[2*l-1]);
        }
      mpz_mul (a[2*l-1], b[l], c[l-1]);
    }

  /* deal with x^k and x^l */
  mpz_set_ui (a[k+l-1], 0);
  /* add b * x^l */
  list_add (a + l, a + l, b, k);
  /* add x^k * c */
  list_add (a + k, a + k, c, l);
}

/*
  Multiplies b[0..k-1] by c[0..k-1], and stores the result in a[0..2k-2].
  (Here, there is no implicit monic leading monomial.)
  Requires at least list_mul_mem(k) cells in t.
 */
void
list_mulmod (listz_t a, listz_t b, listz_t c, unsigned int k, listz_t t,
             mpz_t n)
{
  LIST_MULT_N (a, b, c, k, t);
  list_mod (a, a, 2*k - 1, n);
}

/* puts in G the coefficients from (x-G[0])...(x-G[k-1])
   Needs at least k+list_mul_mem(k/2) >= 2k cells in T */
void
buildG (listz_t G, unsigned int k, listz_t T, int verbose, mpz_t n, char F)
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
   list_mul (T, G, l, G + l, m, T + k);
   list_mod (G, T, k, n);
   
   if (verbose >= 2)
     printf ("Building %c from its roots took %dms\n", F,
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
  Needs space for list_mul_mem(K) coefficients in t.
*/
void
RecursiveDivision (listz_t q, listz_t a, listz_t b, unsigned int K, listz_t t,
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
      LIST_MULT_N (t, q + l, b, k, t + K - 1); /* sets t[0..2*k-2] */
      list_sub (a + l, a + l, t, 2 * k - 1);
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
      LIST_MULT_N (t, q, b, k, t + K - 1);
      list_sub (a, a, t, 2 * k - 1);
      if (k < l) /* don't forget to subtract q[0..k-1] * b[k] */
	for (i=0; i<k; i++)
	  {
	    mpz_mul (t[0], q[i], b[k]);
	    mpz_sub (a[k+i], a[k+i], t[0]);
	  }
    }
  /* normalizes the remainder wrt n */
  list_mod (a, a, K, n);
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
Div3by2 (listz_t q, listz_t a, listz_t b, unsigned int K, 
			 listz_t t, mpz_t n)
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
    LIST_MULT_N (t, q, b, K, t + 2 * K); /* needs 2K memory in t+2K */
    for (i=0;i<=2*K-2;i++) {
      mpz_sub (a[i], a[i], t[i]);
      mpz_mod (a[i], a[i], n);
    }
    /* we could also have a special version of karatsuba/toomcook that directly
       subtract the result to a */
  }
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
