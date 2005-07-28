/* Elliptic Curve Method: toplevel and stage 1 routines.

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
  the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
  MA 02111-1307, USA.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h> /* for ULONG_MAX */
#include <gmp.h>
#include "ecm.h"
#include "ecm-impl.h"

/******************************************************************************
*                                                                             *
*                            Elliptic Curve Method                            *
*                                                                             *
******************************************************************************/

#define mpz_mulmod5(r,s1,s2,m,t) { mpz_mul(t,s1,s2); mpz_mod(r, t, m); }

/* Computes curve parameter A and a starting point (x:1) from a given 
   sigma value.
   If a factor of n was found during the process, returns 1 (and factor 
   in A), 0 otherwise.
*/
static int
get_curve_from_sigma (mpz_t f, mpres_t A, mpres_t x, mpz_t sigma, mpmod_t n)
{
  mpres_t t, u, v, b, z;
  
  mpres_init (t, n);
  mpres_init (u, n);
  mpres_init (v, n);
  mpres_init (b, n);
  mpres_init (z, n);

  mpres_set_z  (u, sigma, n);
  mpres_mul_ui (v, u, 4, n);   /* v = (4*sigma) mod n */
  mpres_mul (t, u, u, n);
  mpres_sub_ui (u, t, 5, n);       /* u = (sigma^2-5) mod n */
  mpres_mul (t, u, u, n);
  mpres_mul (x, t, u, n);          /* x = (u^3) mod n */
  mpres_mul (t, v, v, n);
  mpres_mul (z, t, v, n);          /* z = (v^3) mod n */
  mpres_mul (t, x, v, n);
  mpres_mul_ui (b, t, 4, n);       /* b = (4*x*v) mod n */
  mpres_mul_ui (t, u, 3, n);
  mpres_sub (u, v, u, n);          /* u' = v-u */
  mpres_add (v, t, v, n);          /* v' = (3*u+v) mod n */
  mpres_mul (t, u, u, n);
  mpres_mul (u, t, u, n);          /* u'' = ((v-u)^3) mod n */
  mpres_mul (A, u, v, n);          /* a = (u'' * v') mod n = 
                                      ((v-u)^3 * (3*u+v)) mod n */
  
  /* Normalize b and z to 1 */
  mpres_mul (v, b, z, n);
  if (!mpres_invert (u, v, n)) /* u = (b*z)^(-1) (mod n) */
    {
      mpres_gcd (f, v, n);
      mpres_clear (t, n);
      mpres_clear (u, n);
      mpres_clear (v, n);
      mpres_clear (b, n);
      mpres_clear (z, n);
      return 1;
    }
  
  mpres_mul (v, u, b, n);   /* v = z^(-1) (mod n) */
  mpres_mul (x, x, v, n);   /* x = x * z^(-1) */
  
  mpres_mul (v, u, z, n);   /* v = b^(-1) (mod n) */
  mpres_mul (t, A, v, n);
  mpres_sub_ui (A, t, 2, n);
  
  mpres_clear (t, n);
  mpres_clear (u, n);
  mpres_clear (v, n);
  mpres_clear (b, n);
  mpres_clear (z, n);

  return 0;
}

/* switch from Montgomery's form g*y^2 = x^3 + a*x^2 + x
   to Weierstrass' form          Y^2 = X^3 + A*X + B
   by change of variables x -> g*X-a/3, y -> g*Y.
   We have A = (3-a^2)/(3g^2), X = (3x+a)/(3g), Y = y/g.
*/
static int 
montgomery_to_weierstrass (mpz_t f, mpres_t x, mpres_t y, mpres_t A, mpmod_t n)
{
  mpres_t g;
  
  mpres_init (g, n);
  mpres_add (g, x, A, n);
  mpres_mul (g, g, x, n);
  mpres_add_ui (g, g, 1, n);
  mpres_mul (g, g, x, n);    /* g = x^3+a*x^2+x (y=1) */
  mpres_mul_ui (y, g, 3, n);
  mpres_mul (y, y, g, n);    /* y = 3g^2 */
  if (!mpres_invert (y, y, n)) /* y = 1/(3g^2) temporarily */
    {
      mpres_gcd (f, y, n);
      mpres_clear (g, n);
      return 1;
    }
  
  /* update x */
  mpres_mul_ui (x, x, 3, n); /* 3x */
  mpres_add (x, x, A, n);    /* 3x+a */
  mpres_mul (x, x, g, n);    /* (3x+a)*g */
  mpres_mul (x, x, y, n);    /* (3x+a)/(3g) */

  /* update A */
  mpres_mul (A, A, A, n);    /* a^2 */
  mpres_sub_ui (A, A, 3, n);
  mpres_neg (A, A, n);       /* 3-a^2 */
  mpres_mul (A, A, y, n);    /* (3-a^2)/(3g^2) */

  /* update y */
  mpres_mul_ui (g, g, 3, n); /* 3g */
  mpres_mul (y, y, g, n);    /* (3g)/(3g^2) = 1/g */
  
  mpres_clear (g, n);
  return 0;
}

/* adds Q=(x2:z2) and R=(x1:z1) and puts the result in (x3:z3),
     using 6 muls (4 muls and 2 squares), and 6 add/sub.
   One assumes that Q-R=P or R-Q=P where P=(x:z).
     - n : number to factor
     - u, v, w : auxiliary variables
   Modifies: x3, z3, u, v, w.
   (x3,z3) may be identical to (x2,z2) and to (x,z)
*/
void
add3 (mpres_t x3, mpres_t z3, mpres_t x2, mpres_t z2, mpres_t x1, mpres_t z1, 
      mpres_t x, mpres_t z, mpmod_t n, mpres_t u, mpres_t v, mpres_t w)
{
  mpres_sub (u, x2, z2, n);
  mpres_add (v, x1, z1, n);      /* u = x2-z2, v = x1+z1 */

  mpres_mul (u, u, v, n);        /* u = (x2-z2)*(x1+z1) */

  mpres_add (w, x2, z2, n);
  mpres_sub (v, x1, z1, n);      /* w = x2+z2, v = x1-z1 */

  mpres_mul (v, w, v, n);        /* v = (x2+z2)*(x1-z1) */

  mpres_add (w, u, v, n);        /* w = 2*(x1*x2-z1*z2) */
  mpres_sub (v, u, v, n);        /* v = 2*(x2*z1-x1*z2) */

  mpres_mul (w, w, w, n);        /* w = 4*(x1*x2-z1*z2)^2 */
  mpres_mul (v, v, v, n);        /* v = 4*(x2*z1-x1*z2)^2 */

  if (x == x3) /* same variable: in-place variant */
    {
      /* x3 <- w * z mod n
	 z3 <- x * v mod n */
      mpres_mul (z3, w, z, n);
      mpres_mul (x3, x, v, n);
      mpres_swap (x3, z3, n);
    }
  else
    {
      mpres_mul (x3, w, z, n);   /* x3 = 4*z*(x1*x2-z1*z2)^2 mod n */
      mpres_mul (z3, x, v, n);   /* z3 = 4*x*(x2*z1-x1*z2)^2 mod n */
    }
  /* mul += 6; */
}

/* computes 2P=(x2:z2) from P=(x1:z1), with 5 muls (3 muls and 2 squares)
   and 4 add/sub.
     - n : number to factor
     - b : (a+2)/4 mod n
     - t, u, v, w : auxiliary variables
*/
void
duplicate (mpres_t x2, mpres_t z2, mpres_t x1, mpres_t z1, mpmod_t n, 
           mpres_t b, mpres_t u, mpres_t v, mpres_t w)
{
  mpres_add (u, x1, z1, n);
  mpres_mul (u, u, u, n);   /* u = (x1+z1)^2 mod n */
  mpres_sub (v, x1, z1, n);
  mpres_mul (v, v, v, n);   /* v = (x1-z1)^2 mod n */
  mpres_mul (x2, u, v, n);  /* x2 = u*v = (x1^2 - z1^2)^2 mod n */
  mpres_sub (w, u, v, n);   /* w = u-v = 4*x1*z1 */
  mpres_mul (u, w, b, n);   /* u = w*b = ((A+2)/4*(4*x1*z1)) mod n */
  mpres_add (u, u, v, n);   /* u = (x1-z1)^2+(A+2)/4*(4*x1*z1) */
  mpres_mul (z2, w, u, n);  /* z2 = ((4*x1*z1)*((x1-z1)^2+(A+2)/4*(4*x1*z1))) mod n */
}

/* multiply P=(x:z) by e and puts the result in (x:z). */
void
ecm_mul (mpres_t x, mpres_t z, mpz_t e, mpmod_t n, mpres_t b)
{
  size_t l;
  int negated = 0;
  mpres_t x0, z0, x1, z1, u, v, w;

  /* In Montgomery coordinates, the point at infinity is (0::0) */
  if (mpz_sgn (e) == 0)
    {
      mpz_set_ui (x, 0);
      mpz_set_ui (z, 0);
      return;
    }

  /* The negative of a point (x:y:z) is (x:-y:u). Since we do not compute
     y, e*(x::z) == (-e)*(x::z). */
  if (mpz_sgn (e) < 0)
    {
      negated = 1;
      mpz_neg (e, e);
    }

  if (mpz_cmp_ui (e, 1) == 0)
    goto ecm_mul_end;

  mpres_init (x0, n);
  mpres_init (z0, n);
  mpres_init (x1, n);
  mpres_init (z1, n);
  mpres_init (u, n);
  mpres_init (v, n);
  mpres_init (w, n);

  l = mpz_sizeinbase (e, 2) - 1; /* l >= 1 */

  mpres_set (x0, x, n);
  mpres_set (z0, z, n);
  duplicate (x1, z1, x0, z0, n, b, u, v, w);

  /* invariant: (P1,P0) = ((k+1)P, kP) where k = floor(e/2^l) */

  while (l-- > 0)
    {
      if (mpz_tstbit (e, l)) /* k, k+1 -> 2k+1, 2k+2 */
        {
          add3 (x0, z0, x0, z0, x1, z1, x, z, n, u, v, w); /* 2k+1 */
          duplicate (x1, z1, x1, z1, n, b, u, v, w); /* 2k+2 */
        }
      else /* k, k+1 -> 2k, 2k+1 */
        {
          add3 (x1, z1, x1, z1, x0, z0, x, z, n, u, v, w); /* 2k+1 */
          duplicate (x0, z0, x0, z0, n, b, u, v, w); /* 2k */
        }
    }

  mpres_set (x, x0, n);
  mpres_set (z, z0, n);

  mpres_clear (x0, n);
  mpres_clear (z0, n);
  mpres_clear (x1, n);
  mpres_clear (z1, n);
  mpres_clear (u, n);
  mpres_clear (v, n);
  mpres_clear (w, n);

ecm_mul_end:

  /* Undo negation to avoid changing the caller's e value */
  if (negated)
    mpz_neg (e, e);
}

#define ADD 6.0 /* number of multiplications in an addition */
#define DUP 5.0 /* number of multiplications in a duplicate */

#define START(d,v) ((d)/1.6180339887498948482-128.0+(v))

/* returns the number of modular multiplications for computing
   V_n from V_r * V_{n-r} - V_{n-2r}.
   ADD is the cost of an addition
   DUP is the cost of a duplicate
*/
static double
lucas_cost (unsigned long n, double v)
{
  unsigned long d, e, r;
  double c; /* cost */

  d = n;
  r = (unsigned long) ((double) d / v + 0.5);
  if (r >= n)
    return (ADD * (double) n);
  d = n - r;
  e = 2 * r - n;
  c = DUP + ADD; /* initial duplicate and final addition */
  while (d != e)
    {
      if (d < e)
        {
          r = d;
          d = e;
          e = r;
        }
      if (4 * d <= 5 * e && ((d + e) % 3) == 0)
        { /* condition 1 */
          d = (2 * d - e) / 3;
          e = (e - d) / 2;
          c += 3.0 * ADD; /* 3 additions */
        }
      else if (4 * d <= 5 * e && (d - e) % 6 == 0)
        { /* condition 2 */
          d = (d - e) / 2;
          c += ADD + DUP; /* one addition, one duplicate */
        }
      else if (d <= 4 * e)
        { /* condition 3 */
          d -= e;
          c += ADD; /* one addition */
        }
      else if ((d + e) % 2 == 0)
        { /* condition 4 */
          d = (d - e) / 2;
          c += ADD + DUP; /* one addition, one duplicate */
        }
      /* now d+e is odd */
      else if (d % 2 == 0)
        { /* condition 5 */
          d /= 2;
          c += ADD + DUP; /* one addition, one duplicate */
        }
      /* now d is odd and e is even */
      else if (d % 3 == 0)
        { /* condition 6 */
          d = d / 3 - e;
          c += 3.0 * ADD + DUP; /* three additions, one duplicate */
        }
      else if ((d + e) % 3 == 0)
        { /* condition 7 */
          d = (d - 2 * e) / 3;
          c += 3.0 * ADD + DUP; /* three additions, one duplicate */
        }
      else if ((d - e) % 3 == 0)
        { /* condition 8 */
          d = (d - e) / 3;
          c += 3.0 * ADD + DUP; /* three additions, one duplicate */
        }
      else /* necessarily e is even: catches all cases */
        { /* condition 9 */
          e /= 2;
          c += ADD + DUP; /* one addition, one duplicate */
        }
    }
  
  return c;
}


/* computes kP from P=(xA:zA) and puts the result in (xA:zA). Assumes k>2. */
static void
prac (mpres_t xA, mpres_t zA, unsigned long k, mpmod_t n, mpres_t b,
      mpres_t u, mpres_t v, mpres_t w, mpres_t xB, mpres_t zB, mpres_t xC, 
      mpres_t zC, mpres_t xT, mpres_t zT, mpres_t xT2, mpres_t zT2)
{
  unsigned long d, e, r, i = 0;
  double c, cmin;
  __mpz_struct *tmp;
#define NV 10  
  static double val[NV] =
    { 1.61803398875, 1.72360679775, 1.618347119656, 1.617914406529,
      1.612429949509, 1.632839806089, 1.620181980807, 1.580178728295,
      1.617214616534, 1.38196601125 };
  
  /* chooses the best value of v */
  for (d = 0, cmin = ADD * (double) k; d < NV; d++)
    {
      c = lucas_cost (k, val[d]);
      if (c < cmin)
        {
          cmin = c;
          i = d;
        }
    }
  d = k;
  r = (unsigned long) ((double) d / val[i] + 0.5);
  
  /* first iteration always begins by Condition 3, then a swap */
  d = k - r;
  e = 2 * r - k;
  mpres_set (xB, xA, n);
  mpres_set (zB, zA, n); /* B=A */
  mpres_set (xC, xA, n);
  mpres_set (zC, zA, n); /* C=A */
  duplicate (xA, zA, xA, zA, n, b, u, v, w); /* A = 2*A */
  while (d != e)
    {
      if (d < e)
        {
          r = d;
          d = e;
          e = r;
          mpres_swap (xA, xB, n);
          mpres_swap (zA, zB, n);
        }
      /* do the first line of Table 4 whose condition qualifies */
      if (4 * d <= 5 * e && ((d + e) % 3) == 0)
        { /* condition 1 */
          d = (2 * d - e) / 3;
          e = (e - d) / 2;
          add3 (xT, zT, xA, zA, xB, zB, xC, zC, n, u, v, w); /* T = f(A,B,C) */
          add3 (xT2, zT2, xT, zT, xA, zA, xB, zB, n, u, v, w); /* T2 = f(T,A,B) */
          add3 (xB, zB, xB, zB, xT, zT, xA, zA, n, u, v, w); /* B = f(B,T,A) */
          mpres_swap (xA, xT2, n);
          mpres_swap (zA, zT2, n); /* swap A and T2 */
        }
      else if (4 * d <= 5 * e && (d - e) % 6 == 0)
        { /* condition 2 */
          d = (d - e) / 2;
          add3 (xB, zB, xA, zA, xB, zB, xC, zC, n, u, v, w); /* B = f(A,B,C) */
          duplicate (xA, zA, xA, zA, n, b, u, v, w); /* A = 2*A */
        }
      else if (d <= (4 * e))
        { /* condition 3 */
          d -= e;
          add3 (xT, zT, xB, zB, xA, zA, xC, zC, n, u, v, w); /* T = f(B,A,C) */
          /* circular permutation (B,T,C) */
          tmp = xB;
          xB = xT;
          xT = xC;
          xC = tmp;
          tmp = zB;
          zB = zT;
          zT = zC;
          zC = tmp;
        }
      else if ((d + e) % 2 == 0)
        { /* condition 4 */
          d = (d - e) / 2;
          add3 (xB, zB, xB, zB, xA, zA, xC, zC, n, u, v, w); /* B = f(B,A,C) */
          duplicate (xA, zA, xA, zA, n, b, u, v, w); /* A = 2*A */
        }
      /* now d+e is odd */
      else if (d % 2 == 0)
        { /* condition 5 */
          d /= 2;
          add3 (xC, zC, xC, zC, xA, zA, xB, zB, n, u, v, w); /* C = f(C,A,B) */
          duplicate (xA, zA, xA, zA, n, b, u, v, w); /* A = 2*A */
        }
      /* now d is odd, e is even */
      else if (d % 3 == 0)
        { /* condition 6 */
          d = d / 3 - e;
          duplicate (xT, zT, xA, zA, n, b, u, v, w); /* T = 2*A */
          add3 (xT2, zT2, xA, zA, xB, zB, xC, zC, n, u, v, w); /* T2 = f(A,B,C) */
          add3 (xA, zA, xT, zT, xA, zA, xA, zA, n, u, v, w); /* A = f(T,A,A) */
          add3 (xT, zT, xT, zT, xT2, zT2, xC, zC, n, u, v, w); /* T = f(T,T2,C) */
          /* circular permutation (C,B,T) */
          tmp = xC;
          xC = xB;
          xB = xT;
          xT = tmp;
          tmp = zC;
          zC = zB;
          zB = zT;
          zT = tmp;
        }
      else if ((d + e) % 3 == 0)
        { /* condition 7 */
          d = (d - 2 * e) / 3;
          add3 (xT, zT, xA, zA, xB, zB, xC, zC, n, u, v, w); /* T = f(A,B,C) */
          add3 (xB, zB, xT, zT, xA, zA, xB, zB, n, u, v, w); /* B = f(T,A,B) */
          duplicate (xT, zT, xA, zA, n, b, u, v, w);
          add3 (xA, zA, xA, zA, xT, zT, xA, zA, n, u, v, w); /* A = 3*A */
        }
      else if ((d - e) % 3 == 0)
        { /* condition 8 */
          d = (d - e) / 3;
          add3 (xT, zT, xA, zA, xB, zB, xC, zC, n, u, v, w); /* T = f(A,B,C) */
          add3 (xC, zC, xC, zC, xA, zA, xB, zB, n, u, v, w); /* C = f(A,C,B) */
          mpres_swap (xB, xT, n);
          mpres_swap (zB, zT, n); /* swap B and T */
          duplicate (xT, zT, xA, zA, n, b, u, v, w);
          add3 (xA, zA, xA, zA, xT, zT, xA, zA, n, u, v, w); /* A = 3*A */
        }
      else /* necessarily e is even here */
        { /* condition 9 */
          e /= 2;
          add3 (xC, zC, xC, zC, xB, zB, xA, zA, n, u, v, w); /* C = f(C,B,A) */
          duplicate (xB, zB, xB, zB, n, b, u, v, w); /* B = 2*B */
        }
    }
  
  add3 (xA, zA, xA, zA, xB, zB, xC, zC, n, u, v, w);

  ASSERT(d == 1);
}


/* Input: x is initial point
          A is curve parameter in Montgomery's form:
          g*y^2*z = x^3 + a*x^2*z + x*z^2
          n is the number to factor
	  B1 is the stage 1 bound
   Output: If a factor is found, it is returned in x.
           Otherwise, x contains the x-coordinate of the point computed
           in stage 1 (with z coordinate normalized to 1).
   Return value: non-zero iff a factor is found
*/
static int
ecm_stage1 (mpz_t f, mpres_t x, mpres_t A, mpmod_t n, double B1, double B1done,
	    mpz_t go)
{
  mpres_t b, z, u, v, w, xB, zB, xC, zC, xT, zT, xT2, zT2;
  double q, r;
  int ret = 0;

  mpres_init (b, n);
  mpres_init (z, n);
  mpres_init (u, n);
  mpres_init (v, n);
  mpres_init (w, n);
  mpres_init (xB, n);
  mpres_init (zB, n);
  mpres_init (xC, n);
  mpres_init (zC, n);
  mpres_init (xT, n);
  mpres_init (zT, n);
  mpres_init (xT2, n);
  mpres_init (zT2, n);

  mpres_set_ui (z, 1, n);

  mpres_add_ui (b, A, 2, n);
  mpres_div_2exp (b, b, 2, n); /* b == (A0+2)*B/4, where B=2^(k*GMP_NUMB_LIMB)
                                  for MODMULN or REDC, B=1 otherwise */
#ifndef FULL_REDUCTION
  mpres_semi_normalize (b, mpz_size (n->orig_modulus));
#endif

  /* preload group order */
  if (go != NULL)
    ecm_mul (x, z, go, n, b);

  /* prac() wants multiplicands > 2 */
  for (r = 2.0; r <= B1; r *= 2.0)
    if (r > B1done)
      duplicate (x, z, x, z, n, b, u, v, w);
  
  /* We'll do 3 manually, too (that's what ecm4 did..) */
  for (r = 3.0; r <= B1; r *= 3.0)
    if (r > B1done)
      {
        duplicate (xB, zB, x, z, n, b, u, v, w);
        add3 (x, z, x, z, xB, zB, x, z, n, u, v, w);
      }
  
  q = getprime (2.0); /* Puts 3.0 into q. Next call gives 5.0 */
  for (q = getprime (q); q <= B1; q = getprime (q))
    {
      for (r = q; r <= B1; r *= q)
	if (r > B1done)
	  prac (x, z, (unsigned long) q, n, b, u, v, w, xB, zB, xC, zC, xT,
		zT, xT2, zT2);
    }
  getprime (FREE_PRIME_TABLE); /* free the prime tables, and reinitialize */

  /* Normalize z to 1 */
#ifndef FULL_REDUCTION
  mpres_normalize (z); /* needed for gcd */
#endif
  if (!mpres_invert (u, z, n)) /* Factor found? */
    {
      mpres_gcd (f, z, n);
      ret = 1;
    }
  mpres_mul (x, x, u, n);

  mpres_clear (zT2, n);
  mpres_clear (xT2, n);
  mpres_clear (zT, n);
  mpres_clear (xT, n);
  mpres_clear (zC, n);
  mpres_clear (xC, n);
  mpres_clear (zB, n);
  mpres_clear (xB, n);
  mpres_clear (w, n);
  mpres_clear (v, n);
  mpres_clear (u, n);
  mpres_clear (z, n);
  mpres_clear (b, n);

  return ret;
}

/* choose "optimal" S according to step 2 range B2 */
int
choose_S (mpz_t B2len)
{
  if (mpz_cmp_d (B2len, 1e7) < 0)
    return 1;   /* x^1 */
  else if (mpz_cmp_d (B2len, 1e8) < 0)
    return 2;   /* x^2 */
  else if (mpz_cmp_d (B2len, 1e9) < 0)
    return -3;  /* Dickson(3) */
  else if (mpz_cmp_d (B2len, 1e10) < 0)
    return -6;  /* Dickson(6) */
  else if (mpz_cmp_d (B2len, 3e11) < 0)
    return -12; /* Dickson(12) */
  else
    return -30; /* Dickson(30) */
}

static void
print_expcurves (mpz_t B2min, mpz_t effB2, unsigned long dF, unsigned long k, 
                 int S, int clear)
{
  double prob;
  int i;
  char sep;

  if (!test_verbose (OUTPUT_VERBOSE))
    return;

  rhoinit (256, 10);
  outputf (OUTPUT_VERBOSE, "Expected number of curves to find a factor "
           "of n digits:\n20\t25\t30\t35\t40\t45\t50\t55\t60\t65\n");
  for (i = 20; i <= 65; i += 5)
    {
      sep = (i < 65) ? '\t' : '\n';
      prob = ecmprob (mpz_get_d (B2min), mpz_get_d (effB2),
                      pow (10., i - .5), (double) dF * dF * k, S);
      if (prob > 1. / 10000000)
        outputf (OUTPUT_VERBOSE, "%.0f%c", floor (1. / prob + .5), sep);
      else if (prob > 0.)
        outputf (OUTPUT_VERBOSE, "%.2g%c", floor (1. / prob + .5), sep);
      else
        outputf (OUTPUT_VERBOSE, "Inf%c", sep);
    }

  if (clear)
    rhoinit (1, 0); /* Free memory of rhotable */
}

static void
print_exptime (mpz_t B2min, mpz_t effB2, unsigned long dF, unsigned long k, 
               int S, double tottime, int clear)
{
  double prob, exptime;
  int i;
  char sep;
  
  if (!test_verbose (OUTPUT_VERBOSE))
    return;
  
  rhoinit (256, 10);
  outputf (OUTPUT_VERBOSE, "Expected time to find a factor of n digits:\n"
    "20\t25\t30\t35\t40\t45\t50\t55\t60\t65\n");
  for (i = 20; i <= 65; i += 5)
    {
      sep = (i < 65) ? '\t' : '\n';
      prob = ecmprob (mpz_get_d (B2min), mpz_get_d (effB2), 
                      pow (10., i - .5), (double) dF * dF * k, S);
      exptime = (prob > 0.) ? tottime / prob : HUGE_VAL;
      outputf (OUTPUT_TRACE, "Digits: %d, Total time: %.0f, probability: "
               "%g, expected time: %.0f\n", i, tottime, prob, exptime);
      if (exptime < 1000.)
        outputf (OUTPUT_VERBOSE, "%.0fms%c", exptime, sep);
      else if (exptime < 60000.) /* One minute */
        outputf (OUTPUT_VERBOSE, "%.2fs%c", exptime / 1000., sep);
      else if (exptime < 3600000.) /* One hour */
        outputf (OUTPUT_VERBOSE, "%.2fm%c", exptime / 60000., sep);
      else if (exptime < 86400000.) /* One day */
        outputf (OUTPUT_VERBOSE, "%.2fh%c", exptime / 3600000., sep);
      else if (exptime < 31536000000.) /* One year */
        outputf (OUTPUT_VERBOSE, "%.2fd%c", exptime / 86400000., sep);
      else if (exptime < 31536000000000.) /* One thousand years */
        outputf (OUTPUT_VERBOSE, "%.2fy%c", exptime / 31536000000., sep);
      else if (exptime < 31536000000000000.) /* One million years */
        outputf (OUTPUT_VERBOSE, "%.0fy%c", exptime / 31536000000., sep);
      else if (prob > 0.)
        outputf (OUTPUT_VERBOSE, "%.1gy%c", exptime / 31536000000., sep);
      else 
        outputf (OUTPUT_VERBOSE, "Inf%c", sep);
    }

  if (clear)
    rhoinit (1, 0); /* Free memory of rhotable */
}

/* Input: x is starting point or zero
          sigma is sigma value (if x is set to zero) or 
            A parameter (if x is non-zero) of curve
          n is the number to factor
          go is the initial group order to preload  
          B1, B2 are the stage 1/stage 2 bounds, respectively
          B2min the lower bound for stage 2
          B2scale is the stage 2 scale factor
          k is the number of blocks to do in stage 2
          S is the degree of the Suyama-Brent extension for stage 2
          verbose is verbosity level: 0 no output, 1 normal output,
            2 diagnostic output.
          sigma_is_a: If true, the sigma parameter contains the curve's A value
   Output: f is the factor found.
   Return value: ECM_FACTOR_FOUND_STEPn if a factor was found,
                 ECM_NO_FACTOR_FOUND if no factor was found,
		 ECM_ERROR in case of error.
*/
int
ecm (mpz_t f, mpz_t x, mpz_t sigma, mpz_t n, mpz_t go, double B1done,
     double B1, mpz_t B2min_parm, mpz_t B2_parm, double B2scale, 
     unsigned long k, const int S, int verbose, int repr, int sigma_is_A, 
     FILE *os, FILE* es, char *TreeFilename, double maxmem, 
     gmp_randstate_t rng)
{
  int youpi = ECM_NO_FACTOR_FOUND;
  int base2 = 0;  /* If n is of form 2^n[+-]1, set base to [+-]n */
  int Fermat = 0; /* If base2 > 0 is a power of 2, set Fermat to base2 */
  int po2 = 0;    /* Whether we should use power-of-2 poly degree */
  unsigned int st;
  mpmod_t modulus;
  curve P;
  mpz_t B2min, B2; /* Local B2, B2min to avoid changing caller's values */
  unsigned long dF;
  root_params_t root_params;

  set_verbose (verbose);
  ECM_STDOUT = (os == NULL) ? stdout : os;
  ECM_STDERR = (es == NULL) ? stdout : es;

  /* if n is even, return 2 */
  if (mpz_divisible_2exp_p (n, 1))
    {
      mpz_set_ui (f, 2);
      return ECM_FACTOR_FOUND_STEP1;
    }

  /* now n is odd */

  /* check that B1 is not too large */
  if (B1 > (double) ULONG_MAX)
    {
      outputf (OUTPUT_ERROR, "Error, maximal step 1 bound for ECM is %lu.\n", 
               ULONG_MAX);
      return ECM_ERROR;
    }

  st = cputime ();

  /* See what kind of number we have as that may influence optimal parameter 
     selection. Test for base 2 number */

  if (repr != ECM_MOD_NOBASE2)
    base2 = (abs (repr) >= 16) ? repr : isbase2 (n, BASE2_THRESHOLD);

  /* For for a Fermat number (base2 a positive power of 2) */
  for (Fermat = base2; Fermat > 0 && (Fermat & 1) == 0; Fermat >>= 1);
  if (Fermat == 1) 
    {
      Fermat = base2;
      po2 = 1;
    }
  else
      Fermat = 0;

  if (base2)
    {
      if (mpmod_init_BASE2 (modulus, base2, n) == ECM_ERROR)
        return ECM_ERROR;
    }
  else if (repr == ECM_MOD_MPZ)
    mpmod_init_MPZ (modulus, n);
  else if (repr == ECM_MOD_MODMULN)
    mpmod_init_MODMULN (modulus, n);
  else if (repr == ECM_MOD_REDC)
    mpmod_init_REDC (modulus, n);
  else /* automatic choice of general reduction method */
    mpmod_init (modulus, n, -1);

  mpres_init (P.x, modulus);
  mpres_init (P.y, modulus);
  mpres_init (P.A, modulus);

  mpres_set_z (P.x, x, modulus);
  mpres_set_ui (P.y, 1, modulus);
  
  mpz_init_set (B2min, B2min_parm);
  mpz_init_set (B2, B2_parm);
  
  mpz_init (root_params.i0);

  /* set second stage bound B2: when using polynomial multiplication of
     complexity n^alpha, stage 2 has complexity about B2^(alpha/2), and
     we want stage 2 to take about half of stage 1, thus we choose
     B2 = (c*B1)^(2/alpha). Experimentally, c=1/4 seems to work well.
     For Toom-Cook 3, this gives alpha=log(5)/log(3), and B2 ~ (c*B1)^1.365.
     For Toom-Cook 4, this gives alpha=log(7)/log(4), and B2 ~ (c*B1)^1.424. */

  /* We take the cost of P+1 stage 1 to be about twice that of P-1.
     Since nai"ve P+1 and ECM cost respectively 2 and 11 multiplies per
     addition and duplicate, and both are optimized with PRAC, we can
     assume the ratio remains about 11/2. */

  /* Also scale B2 by what the user said (or by the default scaling of 1.0) */

  if (ECM_IS_DEFAULT_B2(B2))
    mpz_set_d (B2, B2scale * pow (11.0 / 6.0 * B1, 1.424828748));

  /* set B2min */
  if (mpz_sgn (B2min) < 0)
    mpz_set_d (B2min, B1);

  /* Let bestD determine parameters for root generation and the 
     effective B2 */

#ifdef HAVE_NTT
  po2 = 1;
#endif

  if (bestD (&root_params, &k, &dF, B2min, B2, po2, maxmem, 
             (TreeFilename != NULL), modulus) == ECM_ERROR)
    {
      youpi = ECM_ERROR;
      goto end_of_ecm;
    }

  /* Set default degree for Brent-Suyama extension */
  /* We try to keep the time used by the Brent-Suyama extension
     at about 10% of the stage 2 time */
  /* Degree S Dickson polys and x^S are equally fast for ECM, so we go for
     the better Dickson polys whenever possible. For S == 1, 2, they behave
     identically. */

  root_params.S = S;
  if (root_params.S == ECM_DEFAULT_S)
    {
      if (Fermat > 0)
        {
          /* For Fermat numbers, default is 1 (no Brent-Suyama) */
          root_params.S = 1;
        }
      else
        {
          mpz_t t;
          mpz_init (t);
          mpz_sub (t, B2, B2min);
          root_params.S = choose_S (t);
          mpz_clear (t);
        }
    }
  
  if (sigma_is_A == 0)
    {
      /* if sigma=0, generate it at random */
      if (mpz_sgn (sigma) == 0)
        {
          mpz_urandomb (sigma, rng, 32);
          mpz_add_ui (sigma, sigma, 6);
        }

      /* sigma contains sigma value, A and x values must be computed */
      /* If x is specified by user, keep that value */
      youpi = get_curve_from_sigma (f, P.A, P.x, sigma, modulus);
      if (youpi != ECM_NO_FACTOR_FOUND)
	  goto end_of_ecm;
    }
  else
    {
      /* sigma contains the A value */
      mpres_set_z (P.A, sigma, modulus);
      /* TODO: make a valid, random starting point in case none was given */
      /* Problem: this may be as hard as factoring as we'd need to determine
         whether x^3 + a*x^2 + x is a quadratic residue or not */
      /* For now, we'll just chicken out. */
      if (mpz_sgn (x) == 0)
        {
          outputf (OUTPUT_ERROR, 
                   "Error, -A requires a starting point (-x0 x).\n");
	  youpi = ECM_ERROR;
	  goto end_of_ecm;
        }
    }

  /* If a nonzero value is given in x, then we use it as the starting point */
  if (mpz_sgn (x) != 0)
      mpres_set_z (P.x, x, modulus);

  /* Now that the parameters are decided, print info */

  outputf (OUTPUT_NORMAL, "Using B1=%1.0f, B2=", B1);
  if (mpz_cmp_d (B2min, B1) == 0)
      outputf (OUTPUT_NORMAL, "%Zd", B2);
  else
      outputf (OUTPUT_NORMAL, "%Zd-%Zd", B2min, B2);
  outputf (OUTPUT_NORMAL, ", polynomial ");
  if (root_params.S > 0)
      outputf (OUTPUT_NORMAL, "x^%u", root_params.S);
  else
      outputf (OUTPUT_NORMAL, "Dickson(%u)", -root_params.S);
  outputf (OUTPUT_NORMAL, ", sigma=%Zd\n", sigma);

#if 0
  outputf (OUTPUT_VERBOSE, "b2=%1.0f, dF=%lu, k=%lu, d=%lu, d2=%lu, i0=%Zd\n", 
           b2, dF, k, root_params.d1, root_params.d2, root_params.i0);
#else
  outputf (OUTPUT_VERBOSE, "dF=%lu, k=%lu, d=%lu, d2=%lu, i0=%Zd\n", 
           dF, k, root_params.d1, root_params.d2, root_params.i0);
#endif

  if (test_verbose (OUTPUT_RESVERBOSE))
    {
      mpz_t t;

      mpz_init (t);
      mpres_get_z (t, P.A, modulus);
      outputf (OUTPUT_RESVERBOSE, "a=%Zd\n", t);
      mpres_get_z (t, P.x, modulus);
      outputf (OUTPUT_RESVERBOSE, "starting point: x=%Zd\n", t);
      mpz_clear (t);
    }

  if (go != NULL && mpz_cmp_ui (go, 1) > 0)
    outputf (OUTPUT_VERBOSE, "initial group order: %Zd\n", go);

  print_expcurves (B2min, B2, dF, k, root_params.S, 0);

#ifdef HAVE_GWNUM
  /* Right now, we only do base 2 numbers with GWNUM */

  if (base2 != 0 && B1 >= B1done)
      youpi = gw_ecm_stage1 (f, &P, modulus, B1, &B1done, go);

  /* At this point B1 == B1done unless interrupted, or no GWNUM ecm_stage1
     is available */
#endif

  if (B1 > B1done)
    youpi = ecm_stage1 (f, P.x, P.A, modulus, B1, B1done, go);
  
  outputf (OUTPUT_NORMAL, "Step 1 took %ums\n", elltime (st, cputime ()));

  /* Store end-of-stage-1 residue in x in case we write it to a save file, 
     before P.x is converted to Weierstrass form */
  
  mpres_get_z (x, P.x, modulus);

  if (youpi != ECM_NO_FACTOR_FOUND)
    goto end_of_ecm;

  if (test_verbose (OUTPUT_RESVERBOSE)) 
    {
      mpz_t t;
      
      mpz_init (t);
      mpres_get_z (t, P.x, modulus);
      outputf (OUTPUT_RESVERBOSE, "x=%Zd\n", t);
      mpz_clear (t);
    }

#ifdef MONT_ROOTS
  /* If we want to use Montgomery form for generating the roots for stage 2, 
     convert to Weierstrass form only if S != 1 */
  if (root_params.S != 1)
    {
      mpz_t t;
      outputf (OUTPUT_DEVVERBOSE, "ecm: Converting to Weierstrass form\n");
      youpi = montgomery_to_weierstrass (f, P.x, P.y, P.A, modulus);
      if (test_verbose (OUTPUT_TRACE))
        {
          mpz_init (t);
          mpres_get_z (t, P.x, modulus);
          outputf (OUTPUT_TRACE, "P = (%Zd, ", t);
          mpres_get_z (t, P.y, modulus);
          outputf (OUTPUT_TRACE, "%Zd)\n", t);
          mpz_clear (t);
        }
    }
#else
  youpi = montgomery_to_weierstrass (f, P.x, P.y, P.A, modulus);
#endif

  if (youpi == ECM_NO_FACTOR_FOUND && mpz_cmp (B2, B2min) >= 0)
    youpi = stage2 (f, &P, modulus, dF, k, &root_params, ECM_ECM, 
                    TreeFilename);
  
  if (youpi == ECM_NO_FACTOR_FOUND)
    print_exptime (B2min, B2, dF, k, root_params.S, 
                   elltime (st, cputime ()), 1);

end_of_ecm:
  mpres_clear (P.A, modulus);
  mpres_clear (P.y, modulus);
  mpres_clear (P.x, modulus);
  mpmod_clear (modulus);
  mpz_clear (root_params.i0);
  mpz_clear (B2);
  mpz_clear (B2min);

  return youpi;
}
