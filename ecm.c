/* Elliptic Curve Method implementation.

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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h> /* for ULONG_MAX */
#include "gmp.h"
#include "ecm.h"

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
     - t, u, v, w : auxiliary variables
   Modifies: x3, z3, t, u, v, w.
   (x3,z3) may be identical to (x2,z2) and to (x,z)
*/
static void
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
static void
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
      else if (d % 2 == 0)
        { /* condition 5 */
          d /= 2;
          c += ADD + DUP; /* one addition, one duplicate */
        }
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
      else if (e % 2 == 0)
        { /* condition 9 */
          e /= 2;
          c += ADD + DUP; /* one addition, one duplicate */
        }
      else
        {
          fprintf (stderr, "lucas_cost: no condition qualifies for d=%lu e=%lu\n", d, e);
          exit (EXIT_FAILURE);
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
  static double val[] =
    { 1.61803398875, 1.72360679775, 1.618347119656, 1.617914406529,
      1.612429949509, 1.632839806089, 1.620181980807, 1.580178728295,
      1.617214616534, 1.38196601125, 0.0 };
  
  /* chooses the best value of v */
  for (d = 0, cmin = ADD * (double) k; (c = val[d]) != 0.0; d++)
    {
      c = lucas_cost (k, c);
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
      else if (d % 2 == 0)
        { /* condition 5 */
          d /= 2;
          add3 (xC, zC, xC, zC, xA, zA, xB, zB, n, u, v, w); /* C = f(C,A,B) */
          duplicate (xA, zA, xA, zA, n, b, u, v, w); /* A = 2*A */
        }
      else if (d % 3 == 0)
        { /* condition 6 */
          d = d / 3 - e;
          duplicate (xT, zT, xA, zA, n, b, u, v, w); /* T1 = 2*A */
          add3 (xT2, zT2, xA, zA, xB, zB, xC, zC, n, u, v, w); /* T2 = f(A,B,C) */
          add3 (xA, zA, xT, zT, xA, zA, xA, zA, n, u, v, w); /* A = f(T1,A,A) */
          add3 (xT, zT, xT, zT, xT2, zT2, xC, zC, n, u, v, w); /* T1 = f(T1,T2,C) */
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
          add3 (xT, zT, xA, zA, xB, zB, xC, zC, n, u, v, w); /* T1 = f(A,B,C) */
          add3 (xB, zB, xT, zT, xA, zA, xB, zB, n, u, v, w); /* B = f(T1,A,B) */
          duplicate (xT, zT, xA, zA, n, b, u, v, w);
          add3 (xA, zA, xA, zA, xT, zT, xA, zA, n, u, v, w); /* A = 3*A */
        }
      else if ((d - e) % 3 == 0)
        { /* condition 8 */
          d = (d - e) / 3;
          add3 (xT, zT, xA, zA, xB, zB, xC, zC, n, u, v, w); /* T1 = f(A,B,C) */
          add3 (xC, zC, xC, zC, xA, zA, xB, zB, n, u, v, w); /* C = f(A,C,B) */
          mpres_swap (xB, xT, n);
          mpres_swap (zB, zT, n); /* swap B and T */
          duplicate (xT, zT, xA, zA, n, b, u, v, w);
          add3 (xA, zA, xA, zA, xT, zT, xA, zA, n, u, v, w); /* A = 3*A */
        }
      else if (e % 2 == 0)
        { /* condition 9 */
          e /= 2;
          add3 (xC, zC, xC, zC, xB, zB, xA, zA, n, u, v, w); /* C = f(C,B,A) */
          duplicate (xB, zB, xB, zB, n, b, u, v, w); /* B = 2*B */
        }
      else
        {
          fprintf (stderr, "no condition qualifies for d=%lu e=%lu\n", d, e);
          exit (EXIT_FAILURE);
        }
    }
  
  add3 (xA, zA, xA, zA, xB, zB, xC, zC, n, u, v, w);
  
#ifdef DEBUG
  if (d != 1)
    {
      fprintf (stderr, "d!=1 at the end of PRAC\n");
      exit (EXIT_FAILURE);
    }
#endif
}


/* Input: x is initial point
          A is curve parameter in Montgomery's form:
          g*y^2*z = x^3 + a*x^2*z + x*z^2
          n is the number to factor
	  B1 is the stage 1 bound
	  B2 is the stage 2 bound
          verbose is verbosity level: 0 no output, 1 normal output,
            2 diagnostic output
   Output: If a factor is found, it is returned in x.
           Otherwise, x contains the x-coordinate of the point computed
           in stage 1 (with z coordinate normalized to 1).
   Return value: non-zero iff a factor is found
*/
static int
ecm_stage1 (mpz_t f, mpres_t x, mpres_t A, mpmod_t n, double B1, double B1done)
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
  mpres_div_2exp (b, b, 2, n); /* b == (A+2)/4 */

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
  /* check that B1 is not too large */
  if (B1 > (double) ULONG_MAX)
    {
      fprintf (stderr, "Error, maximal step1 bound B1 for ECM is %lu\n", ULONG_MAX);
      exit (EXIT_FAILURE);
    }
  for (q = getprime (q); q <= B1; q = getprime (q))
    {
      for (r = q; r <= B1; r *= q)
	if (r > B1done)
	  prac (x, z, (unsigned long) q, n, b, u, v, w, xB, zB, xC, zC, xT, zT, xT2, zT2);
    }
  getprime (0.0); /* free the prime tables, and reinitialize */


  /* Normalize z to 1 */
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

/* Input: x is starting point or zero
          sigma is sigma value (if x is set to zero) or 
            A parameter (if x is non-zero) of curve
          n is the number to factor
          B1, B2 are the stage 1/stage 2 bounds, respectively
          k is the number of blocks to do in stage 2
          S is the degree of the Suyama-Brent extension for stage 2
          verbose is verbosity level: 0 no output, 1 normal output,
            2 diagnostic output.
          sigma_is_a: If true, the sigma parameter contains the curve's A value
   Output: f is the factor found.
   Return value: non-zero iff a factor was found.
*/
int
ecm (mpz_t f, mpz_t x, mpz_t sigma, mpz_t n, double B1done, double B1,
     double B2min, double B2, double B2scale, unsigned int k, int S, int verbose,
     int repr, int sigma_is_A)
{
  int youpi = 0, st;
  mpmod_t modulus;
  curve P;

  st = cputime ();

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

  if (B2 == 0.0)
    B2 = pow (11.0 / 6.0 * B1, 1.424828748);

  /* Scale B2 by what the user said (or by the default scaling of 1.0) */
  B2 *= B2scale;

  /* Set default degree for Brent-Suyama extension */
  /* We try to keep the time used by the Brent-Suyama extension
     at about 10% of the stage 2 time */
  /* Degree S Dickson polys and x^S are equally fast for ECM, so we go for
     the better Dickson polys whenever possible. For S == 1, 2, they behave
     identically. */

  if (S == 0)
    {
      if (B2 - B2min < 1.e7)
        S = 1;    /* x^1 */
      else if (B2 - B2min < 1.e8)
        S = 2;   /* x^2 */
      else if (B2 - B2min < 8.e8)
        S = -3;   /* Dickson(3) */
      else if (B2 - B2min < 4.e9)
        S = -6;  /* Dickson(6) */
      else if (B2 - B2min < 3.e10)
        S = -12;  /* Dickson(12) */
      else
        S = -30;  /* Dickson(30) */
    }
  
  if (repr == 1)
    mpmod_init_MPZ (modulus, n);
  else if (repr == 2)
    mpmod_init_MODMULN (modulus, n);
  else if (repr == 3)
    mpmod_init_REDC (modulus, n);
  else if (repr > 16)
    mpmod_init_BASE2 (modulus, repr, n);
  else
    mpmod_init (modulus, n, repr, verbose);

  mpres_init (P.x, modulus);
  mpres_init (P.A, modulus);

  if (sigma_is_A == 0)
    {
      /* sigma contains sigma value, A value must be computed */
      if (verbose >= 1)
        {
          printf ("Using B1=%1.0f, B2=", B1);
	  if (B2min <= B1)
	    printf ("%1.0f", B2);
	  else
	    printf ("%1.0f-%1.0f", B2min, B2);
	  printf (", polynomial ");
	  if (S > 0)
	    printf ("x^%u", S);
	  else
	    printf ("Dickson(%u)", -S);
	  printf (", sigma=");
          mpz_out_str (stdout, 10, sigma);
          printf ("\n");
          fflush (stdout);
        }
      if ((youpi = get_curve_from_sigma (f, P.A, P.x, sigma, modulus)))
        goto end_of_ecm;
    }
  else
    {
      /* sigma contains the A value */
      mpres_set_z (P.A, sigma, modulus);
      /* TODO: make a valid, random starting point in case none was given */
      /* For now, we'll just chicken out. */
      if (mpz_sgn (x) == 0)
        {
          fprintf (stderr, "The -A option requires that a starting point is given as well by the -x0\noption\n");
          exit (EXIT_FAILURE);
        }
    }

  /* If a nonzero value is given in x, then we use it as the starting point */
  if (mpz_sgn (x) != 0)
      mpres_set_z (P.x, x, modulus);
  
  if (verbose >= 2)
    {
      printf ("a=");
      mpres_out_str (stdout, 10, P.A, modulus);
      printf ("\nstarting point: x=");
      mpres_out_str (stdout, 10, P.x, modulus);
      printf ("\n");
    }

  if (B1 > B1done)
    youpi = ecm_stage1 (f, P.x, P.A, modulus, B1, B1done);
  
  st = cputime () - st;
  
  if (verbose >= 1)
    {
      printf ("Step 1 took %dms\n", st);
      fflush (stdout);
    }

  /* Store end-of-stage-1 residue in x in case we write it to a save file, 
     before P.x is converted to Weierstrass form */
  
  mpres_get_z (x, P.x, modulus);

  if (youpi != 0) /* a factor was found */
    goto end_of_ecm;

  if (verbose >= 2) 
    {
      printf ("x=");
      mpres_out_str (stdout, 10, P.x, modulus);
      printf ("\n");
      fflush (stdout);
    }

  mpres_init (P.y, modulus);

  youpi = montgomery_to_weierstrass (f, P.x, P.y, P.A, modulus);

  if (youpi == 0)
    youpi = stage2 (f, &P, modulus, B2min, B2, k, S, verbose, EC_METHOD, st);
  
  mpres_clear (P.y, modulus);

end_of_ecm:
  mpres_clear (P.x, modulus);
  mpres_clear (P.A, modulus);
  mpmod_clear (modulus);

  return youpi;
}
