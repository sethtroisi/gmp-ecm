/* Elliptic Curve Method implementation.

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

#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "ecm.h"

/******************************************************************************
*                                                                             *
*                            Elliptic Curve Method                            *
*                                                                             *
******************************************************************************/

#define mpz_mulmod(r,s1,s2,m,t) { mpz_mul(t,s1,s2); mpz_mod(r, t, m); }


/* Computes curve parameter A and a starting point (x:1) from a given 
   sigma value.
   If a factor of n was found during the process, returns 1 (and factor 
   in A), 0 otherwise.
*/
int
get_curve_from_sigma (mpz_t A, mpz_t x, mpz_t sigma, mpz_t n)
{
  mpz_t t, u, v, b, z;
  
  mpz_init (t);
  mpz_init (u);
  mpz_init (v);
  mpz_init (b);
  mpz_init (z);

  mpz_mul_2exp (t, sigma, 2);
  mpz_mod (v, t, n);            /* v = (4*sigma) mod n */
  mpz_mul (t, sigma, sigma);
  mpz_sub_ui (t, t, 5);
  mpz_mod (u, t, n);            /* u = (sigma^2-5) mod n */
  mpz_mul (t, u, u);
  mpz_mul (t, t, u);
  mpz_mod (x, t, n);            /* x = (u^3) mod n */
  mpz_mul (t, v, v);
  mpz_mul (t, t, v);
  mpz_mod (z, t, n);            /* z = (v^3) mod n */
  mpz_mul (t, x, v);
  mpz_mul_2exp (t, t, 2);
  mpz_mod (b, t, n);            /* b = (4*x*v) mod n */
  mpz_mul_ui (t, u, 3);
  mpz_sub (u, v, u);            /* u' = v-u */
  mpz_add (t, t, v);
  mpz_mod (v, t, n);            /* v' = (3*u+v) mod n */
  mpz_mul (t, u, u);
  mpz_mul (t, t, u);
  mpz_mod (u, t, n);            /* u'' = ((v-u)^3) mod n */
  mpz_mulmod (A, u, v, n, t);   /* a = (u'' * v') mod n = 
                                   ((v-u)^3 * (3*u+v)) mod n */
  
  /* Normalize b and z to 1 */
  mpz_mulmod (v, b, z, n, t);
  mpz_gcdext (t, u, (__mpz_struct *) NULL, v, n); /* u = (b*z)^(-1) (mod n) */
  if (mpz_cmp_ui (t, 1) != 0)
    {
      mpz_set (A, t);
      mpz_clear (t);
      mpz_clear (u);
      mpz_clear (v);
      mpz_clear (b);
      mpz_clear (z);
      return 1;
    }
  
  mpz_mulmod (v, u, b, n, t);   /* v = z^(-1) (mod n) */
  mpz_mulmod (x, x, v, n, t);   /* x = x * z^(-1) */
  
  mpz_mulmod (v, u, z, n, t);   /* v = b^(-1) (mod n) */
  mpz_mul (t, A, v);
  mpz_sub_ui (t, t, 2);
  mpz_mod (A, t, n);
  
  mpz_clear (t);
  mpz_clear (u);
  mpz_clear (v);
  mpz_clear (b);
  mpz_clear (z);

  return 0;
}

/* adds Q=(x2:z2) and R=(x1:z1) and puts the result in (x3:z3),
     using 5/6 mul, 6 add/sub and 6 mod. 
   One assumes that Q-R=P or R-Q=P where P=(x:z).
     - n : number to factor
     - t, u, v, w : auxiliary variables
   Modifies: x3, z3, t, u, v, w.
   (x3,z3) may be identical to (x2,z2) and to (x,z)
*/
void
add3 (mpz_t x3, mpz_t z3, mpz_t x2, mpz_t z2, mpz_t x1, mpz_t z1, mpz_t x,
      mpz_t z, mpz_t n, mpz_t t, mpz_t u, mpz_t v, mpz_t w)
{
  mpz_sub (u, x2, z2);
  mpz_add (v, x1, z1);          /* u = x2-z2, v = x1+z1 */

  mpz_mulmod (u, u, v, n, t);   /* u = (x2-z2)*(x1+z1) */

  mpz_add (w, x2, z2);
  mpz_sub (v, x1, z1);          /* w = x2+z2, v = x1-z1 */

  mpz_mulmod (v, w, v, n, t);   /* v = (x2+z2)*(x1-z1) */

  mpz_add (w, u, v);            /* w = 2*(x1*x2-z1*z2) */
  mpz_sub (v, u, v);            /* v = 2*(x2*z1-x1*z2) */

  mpz_mulmod (w, w, w, n, t);   /* w = 4*(x1*x2-z1*z2)^2 */
  mpz_mulmod (v, v, v, n, t);   /* v = 4*(x2*z1-x1*z2)^2 */

  if (x == x3) /* mpz_cmp() ? */
    {
      mpz_set (u, x);
      mpz_mulmod (w, w, z, n, t);
      mpz_mulmod (z3, u, v, n, t);
      mpz_set (x3, w);
    }
  else
    {
      mpz_mulmod(x3, w, z, n, t); /* x3 = 4*z*(x1*x2-z1*z2)^2 mod n */
      mpz_mulmod(z3, x, v, n, t); /* z3 = 4*x*(x2*z1-x1*z2)^2 mod n */
    }
  /* mul += 6; */
}

/* computes 2P=(x2:z2) from P=(x1:z1), with 5 mul, 4 add/sub, 5 mod.
     - n : number to factor
     - b : (a+2)/4 mod n
     - t, u, v, w : auxiliary variables
*/
void
duplicate (mpz_t x2, mpz_t z2, mpz_t x1, mpz_t z1, mpz_t n, mpz_t b,
           mpz_t t, mpz_t u, mpz_t v, mpz_t w)
{
  mpz_add (u, x1, z1);
  mpz_mulmod (u, u, u, n, t);   /* u = (x1+z1)^2 mod n */
  mpz_sub (v, x1, z1);
  mpz_mulmod (v, v, v, n, t);   /* v = (x1-z1)^2 mod n */
  mpz_mulmod (x2, u, v, n, t);  /* x2 = u*v = (x1^2 - z1^2)^2 mod n */
  mpz_sub (w, u, v);            /* w = u-v = 4*x1*z1 */
  mpz_mulmod (u, w, b, n, t);   /* u = w*b = ((A+2)/4*(4*x1*z1)) mod n */
  mpz_add (u, u, v);            /* u = (x1-z1)^2+(A+2)/4*(4*x1*z1) */
  mpz_mulmod (z2, w, u, n, t);  /* z2 = ((4*x1*z1)*((x1-z1)^2+(A+2)/4*(4*x1*z1))) mod n */
}

#define ADD 6 /* number of multiplications in an addition */
#define DUP 5 /* number of multiplications in a duplicate */

#define START(d,v) ((d)/1.6180339887498948482-128.0+(v))

/* returns the number of modular multiplications for computing
   V_n from V_r * V_{n-r} - V_{n-2r}.
   ADD is the cost of an addition
   DUP is the cost of a duplicate
*/
unsigned int
lucas_cost (unsigned n, double v)
{
  unsigned int c, d, e, r;

  d = n;
  r = (unsigned int) ((double) d / v + 0.5);
  if (r >= n)
    return (ADD * n);
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
          c += 3 * ADD; /* 3 additions */
        }
      else if (4 * d <= 5 * e && (d - e) % 6 == 0)
        { /* condition 2 */
          d = (d - e) / 2;
          c += ADD+DUP; /* one addition, one duplicate */
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
          c += ADD+DUP; /* one addition, one duplicate */
        }
      else if (d % 3 == 0)
        { /* condition 6 */
          d = d / 3 - e;
          c += 3 * ADD + DUP; /* three additions, one duplicate */
        }
      else if ((d + e) % 3 == 0)
        { /* condition 7 */
          d = (d - 2 * e) / 3;
          c += 3 * ADD + DUP; /* three additions, one duplicate */
        }
      else if ((d - e) % 3 == 0)
        { /* condition 8 */
          d = (d - e) / 3;
          c += 3 * ADD + DUP; /* three additions, one duplicate */
        }
      else if (e % 2 == 0)
        { /* condition 9 */
          e /= 2;
          c += ADD + DUP; /* one addition, one duplicate */
        }
      else
        {
          printf ("lucas_cost: no condition qualifies for d=%u e=%u\n", d, e);
          exit (EXIT_FAILURE);
        }
    }
  
  return (c);
}


#define NV 10

/* computes kP from P=(xA:zA) and puts the result in (xA:zA). Assumes k>2. */
void
prac (mpz_t xA, mpz_t zA, unsigned int k, mpz_t n, mpz_t b, mpz_t t, mpz_t u,
      mpz_t v, mpz_t w, mpz_t xB, mpz_t zB, mpz_t xC, mpz_t zC, mpz_t xT,
      mpz_t zT, mpz_t xT2, mpz_t zT2)
{
  unsigned int d, e, r, i = 0;
  __mpz_struct *tmp;
  static double val[NV] =
    { 1.61803398875, 1.72360679775, 1.618347119656, 1.617914406529,
      1.612429949509, 1.632839806089, 1.620181980807, 1.580178728295,
      1.617214616534, 1.38196601125 };
  
  /* chooses the best value of v */
  for (d = 0, r = ADD * k; d < NV; d++)
    {
      e = lucas_cost (k, val[d]);
      if (e < r)
        {
          r = e;
          i = d;
        }
    }
  d = k;
  r = (int) ((double) d / val[i] + 0.5);
  
  /* first iteration always begins by Condition 3, then a swap */
  d = k - r;
  e = 2 * r - k;
  mpz_set (xB, xA);
  mpz_set (zB, zA); /* B=A */
  mpz_set (xC, xA);
  mpz_set (zC, zA); /* C=A */
  duplicate (xA, zA, xA, zA, n, b, t, u, v, w); /* A = 2*A */
  while (d != e)
    {
      if (d < e)
        {
          r = d;
          d = e;
          e = r;
          mpz_swap (xA, xB);
          mpz_swap (zA, zB);
        }
      /* do the first line of Table 4 whose condition qualifies */
      if (4 * d <= 5 * e && ((d + e) % 3) == 0)
        { /* condition 1 */
          d = (2 * d - e) / 3;
          e = (e - d) / 2;
          add3 (xT, zT, xA, zA, xB, zB, xC, zC, n, t, u, v, w); /* T = f(A,B,C) */
          add3 (xT2, zT2, xT, zT, xA, zA, xB, zB, n, t, u, v, w); /* T2 = f(T,A,B) */
          add3 (xB, zB, xB, zB, xT, zT, xA, zA, n, t, u, v, w); /* B = f(B,T,A) */
          mpz_swap (xA, xT2);
          mpz_swap (zA, zT2); /* swap A and T2 */
        }
      else if (4 * d <= 5 * e && (d - e) % 6 == 0)
        { /* condition 2 */
          d = (d - e) / 2;
          add3 (xB, zB, xA, zA, xB, zB, xC, zC, n, t, u, v, w); /* B = f(A,B,C) */
          duplicate (xA, zA, xA, zA, n, b, t, u, v, w); /* A = 2*A */
        }
      else if (d <= (4 * e))
        { /* condition 3 */
          d -= e;
          add3 (xT, zT, xB, zB, xA, zA, xC, zC, n, t, u, v, w); /* T = f(B,A,C) */
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
          add3 (xB, zB, xB, zB, xA, zA, xC, zC, n, t, u, v, w); /* B = f(B,A,C) */
          duplicate (xA, zA, xA, zA, n, b, t, u, v, w); /* A = 2*A */
        }
      else if (d % 2 == 0)
        { /* condition 5 */
          d /= 2;
          add3 (xC, zC, xC, zC, xA, zA, xB, zB, n, t, u, v, w); /* C = f(C,A,B) */
          duplicate (xA, zA, xA, zA, n, b, t, u, v, w); /* A = 2*A */
        }
      else if (d % 3 == 0)
        { /* condition 6 */
          d = d / 3 - e;
          duplicate (xT, zT, xA, zA, n, b, t, u, v, w); /* T1 = 2*A */
          add3 (xT2, zT2, xA, zA, xB, zB, xC, zC, n, t, u, v, w); /* T2 = f(A,B,C) */
          add3 (xA, zA, xT, zT, xA, zA, xA, zA, n, t, u, v, w); /* A = f(T1,A,A) */
          add3 (xT, zT, xT, zT, xT2, zT2, xC, zC, n, t, u, v, w); /* T1 = f(T1,T2,C) */
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
          add3 (xT, zT, xA, zA, xB, zB, xC, zC, n, t, u, v, w); /* T1 = f(A,B,C) */
          add3 (xB, zB, xT, zT, xA, zA, xB, zB, n, t, u, v, w); /* B = f(T1,A,B) */
          duplicate (xT, zT, xA, zA, n, b, t, u, v, w);
          add3 (xA, zA, xA, zA, xT, zT, xA, zA, n, t, u, v, w); /* A = 3*A */
        }
      else if ((d - e) % 3 == 0)
        { /* condition 8 */
          d = (d - e) / 3;
          add3 (xT, zT, xA, zA, xB, zB, xC, zC, n, t, u, v, w); /* T1 = f(A,B,C) */
          add3 (xC, zC, xC, zC, xA, zA, xB, zB, n, t, u, v, w); /* C = f(A,C,B) */
          mpz_swap (xB, xT);
          mpz_swap (zB, zT); /* swap B and T */
          duplicate (xT, zT, xA, zA, n, b, t, u, v, w);
          add3 (xA, zA, xA, zA, xT, zT, xA, zA, n, t, u, v, w); /* A = 3*A */
        }
      else if (e % 2 == 0)
        { /* condition 9 */
          e /= 2;
          add3 (xC, zC, xC, zC, xB, zB, xA, zA, n, t, u, v, w); /* C = f(C,B,A) */
          duplicate (xB, zB, xB, zB, n, b, t, u, v, w); /* B = 2*B */
        }
      else
        {
          fprintf (stderr, "no condition qualifies for d=%u e=%u\n", d, e);
          exit (EXIT_FAILURE);
        }
    }
  
  add3 (xA, zA, xA, zA, xB, zB, xC, zC, n, t, u, v, w);
  
#ifdef DEBUG
  if (d != 1)
    {
      fprintf (stderr, "d!=1 at the end of PRAC\n");
      exit (EXIT_FAILURE);
    }
#endif
}


/* Input: x is initial point
          A is curve parameter A
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
int
ecm_stage1 (mpz_t x, mpz_t A, mpz_t n, double B1, int verbose)
{
  mpz_t b, z, t, u, v, w, xB, zB, xC, zC, xT, zT, xT2, zT2;
  double q, r;
  int ret = 0;

  mpz_init2 (b, mpz_size (n) + 1);
  mpz_init2 (z, mpz_size (n) + 1);
  mpz_init2 (t, 3 * mpz_size (n));
  mpz_init2 (u, mpz_size (n) + 1);
  mpz_init2 (v, mpz_size (n) + 1);
  mpz_init2 (w, mpz_size (n) + 1);
  mpz_init2 (xB, mpz_size (n) + 1);
  mpz_init2 (zB, mpz_size (n) + 1);
  mpz_init2 (xC, mpz_size (n) + 1);
  mpz_init2 (zC, mpz_size (n) + 1);
  mpz_init2 (xT, mpz_size (n) + 1);
  mpz_init2 (zT, mpz_size (n) + 1);
  mpz_init2 (xT2, mpz_size (n) + 1);
  mpz_init2 (zT2, mpz_size (n) + 1);
  
  mpz_set_ui (z, 1);
  
  mpz_add_ui (b, A, 2);
  if (mpz_odd_p (b))
    mpz_add (b, b, n); /* Assumes n is odd */
  mpz_tdiv_q_2exp (b, b, 1);
  if (mpz_odd_p (b))
    mpz_add (b, b, n);
  mpz_tdiv_q_2exp (b, b, 1); /* Now b == (A+2)/4 */

  /* prac() wants multiplicands > 2 */
  for (r = 2.0; r <= B1; r *= 2.0)
    duplicate (x, z, x, z, n, b, t, u, v, w);
  
  /* We'll do 3 manually, too (that's what ecm4 did..) */
  for (r = 3.0; r <= B1; r *= 3.0)
    {
      duplicate (xB, zB, x, z, n, b, t, u, v, w);
      add3 (x, z, x, z, xB, zB, x, z, n, t, u, v, w);
    }
  
  q = getprime (2.0); /* Puts 3.0 into q. Next call gives 5.0 */
  for (q = getprime (q); q <= B1; q = getprime (q))
    for (r = q; r <= B1; r *= q)
      prac (x, z, (int) q, n, b, t, u, v, w, xB, zB, xC, zC, xT, zT, xT2, zT2);

  mpz_gcdext (t, u, (__mpz_struct *) NULL, z, n);
  if (mpz_cmp_ui (t, 1) > 0) /* Factor found? */
    {
      ret = 1;
      mpz_set (x, t);
    }
  else
    mpz_mulmod (x, x, u, n, t);

  mpz_clear (zT2);
  mpz_clear (xT2);
  mpz_clear (zT);
  mpz_clear (xT);
  mpz_clear (zC);
  mpz_clear (xC);
  mpz_clear (zB);
  mpz_clear (xB);
  mpz_clear (w);
  mpz_clear (v);
  mpz_clear (u);
  mpz_clear (t);
  mpz_clear (z);
  mpz_clear (b);
  
  return ret;
}

/* Input: p is starting point or zero
          sigma is sigma value (if p is set to zero) or 
            A parameter (if p is non-zero) of curve
          n is the number to factor
          B1, B2 are the stage 1/stage 2 bounds, respectively
          k is the number of blocks to do in stage 2
          S is the degree of the Suyama-Brent extension for stage 2
          verbose is verbosity level: 0 no output, 1 normal output,
            2 diagnostic output
*/
int
ecm (mpz_t p, mpz_t sigma, mpz_t n, double B1, double B2, unsigned int k,
     unsigned int S, int verbose)
{
  int youpi, st;
  mpz_t A;

  st = cputime ();
  
  mpz_init (A);
  if (mpz_sgn (p) == 0)
    {
      /* sigma contains sigma value, starting point and A value must be 
         computed */
      if (verbose >= 1)
        {
          printf ("Using sigma=");
          mpz_out_str (stdout, 10, sigma);
          printf("\n");
        }
      get_curve_from_sigma (A, p, sigma, n);
    }
  else
    {
      /* sigma contains the A value, p contains starting point */
      mpz_set (A, sigma);
    }
  
  if (verbose >= 2)
    {
      printf("A=");
      mpz_out_str(stdout, 10, A);
      printf("\nstarting point: x=");
      mpz_out_str(stdout, 10, p);
      printf("\n");
    }
  
  youpi = ecm_stage1 (p, A, n, B1, verbose);

  if (verbose >= 1)
    {
      printf ("Stage 1 took %dms\n", cputime () - st);
      fflush (stdout);
    }

  if (youpi != 0) /* a factor was found */
    return 1;

  if (verbose >= 2)
    gmp_printf ("x=%Zd\n", p);

//  return (B2 > B1) ? stage2 (p, n, B2, k, S, verbose, 1, EC_METHOD) : 0;
  return 0;
}
