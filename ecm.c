/* Elliptic Curve Method: toplevel and stage 1 routines.

Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011,
2012, 2016, 2024 Paul Zimmermann, Alexander Kruppa, Cyril Bouvier,
David Cleaver, Philip McLaughlin.

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ecm-impl.h"
#include "getprime_r.h"
#include <math.h>

#ifdef HAVE_ADDLAWS
#include "addlaws.h"
#endif

#ifdef HAVE_LIMITS_H
# include <limits.h>
#else
# define ULONG_MAX __GMP_ULONG_MAX
#endif

#ifdef TIMING_CRT
extern int mpzspv_from_mpzv_slow_time, mpzspv_to_mpzv_time,
  mpzspv_normalise_time;
#endif

/* the following factor takes into account the smaller expected smoothness
   for Montgomery's curves (batch mode) with respect to Suyama's curves */
/* For param 1 we use A=4d-2 with d a square (see main.c). In that
   case, Cyril Bouvier and Razvan Barbulescu have shown that the average
   expected torsion is that of a generic Suyama curve multiplied by the
   constant 2^(1/3)/(3*3^(1/128)) */
#define EXTRA_SMOOTHNESS_SQUARE 0.416384512396064
/* For A=4d-2 (param 3) for d a random integer, the average expected torsion 
   is that of a generic Suyama curve multiplied by the constant 
   1/(3*3^(1/128)) */
#define EXTRA_SMOOTHNESS_32BITS_D 0.330484606500389
/******************************************************************************
*                                                                             *
*                            Elliptic Curve Method                            *
*                                                                             *
******************************************************************************/

void duplicate (mpres_t, mpres_t, mpres_t, mpres_t, mpmod_t, mpres_t, 
                mpres_t, mpres_t, mpres_t) ATTRIBUTE_HOT;
void add3 (mpres_t, mpres_t, mpres_t, mpres_t, mpres_t, mpres_t, mpres_t, 
           mpres_t, mpmod_t, mpres_t, mpres_t, mpres_t) ATTRIBUTE_HOT;

#define mpz_mulmod5(r,s1,s2,m,t) { mpz_mul(t,s1,s2); mpz_mod(r, t, m); }

/* PBMcL additions  for Lucas chain codes */

/* 3-bit chain start sequences. Covers all possible Lucas chains */
#define CHAIN_START_5_8_13 0x7
#define CHAIN_START_5_8_11 0x6
#define CHAIN_START_5_8_10 0x5
#define CHAIN_START_5_7    0x4
#define CHAIN_START_5_6    0x3
#define CHAIN_START_4_7    0x2
#define CHAIN_START_4_5    0x1
#define CHAIN_START_4_6    0x0 /* precludes a completely zero code */

typedef struct
{
	/* note: "parent" is the immediately previous chain element */
	uint64_t	value;			/* integer value of this chain element */
	uint8_t	comp_offset_1;	/* larger summand (summand_1) component index counting back from parent (parent = 0) */
	uint8_t	comp_offset_2;	/* smaller summand (summand_2) component index counting back from parent */
	uint8_t	dif_offset;		/* component index of (summand_1 - summand_2) counting back from parent */
								/* note: dif_offset = 0 will indicate that this is a doubled element */
} chain_element;

/* prototypes */
static void max_continuation( chain_element *, uint8_t *, uint8_t );
static uint8_t generate_Lucas_chain( uint64_t, uint64_t, chain_element * );

/* end PBMcL additions */

/* switch from Montgomery's form g*y^2 = x^3 + a*x^2 + x
   to Weierstrass' form          Y^2 = X^3 + A*X + B
   by change of variables x -> g*X-a/3, y -> g*Y.
   We have A = (3-a^2)/(3g^2), X = (3x+a)/(3g), Y = y/g.
   If a factor is found during the modular inverse, returns 
   ECM_FACTOR_FOUND_STEP1 and the factor in f, otherwise returns
   ECM_NO_FACTOR_FOUND.

   The input value of y is the y0 value in the initial equation:
   g*y0^2 = x0^3 + a*x0^2 + x0.
*/
int 
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
      return ECM_FACTOR_FOUND_STEP1;
    }
  
  /* update x */
  mpres_mul_ui (x, x, 3, n); /* 3x */
  mpres_add (x, x, A, n);    /* 3x+a */
  mpres_mul (x, x, g, n);    /* (3x+a)*g */
  mpres_mul (x, x, y, n);    /* (3x+a)/(3g) */

  /* update A */
  mpres_sqr (A, A, n);    /* a^2 */
  mpres_sub_ui (A, A, 3, n);
  mpres_neg (A, A, n);       /* 3-a^2 */
  mpres_mul (A, A, y, n);    /* (3-a^2)/(3g^2) */

  /* update y */
  mpres_mul_ui (g, g, 3, n); /* 3g */
  mpres_mul (y, y, g, n);    /* (3g)/(3g^2) = 1/g */
  
  mpres_clear (g, n);
  return ECM_NO_FACTOR_FOUND;
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

  mpres_sqr (w, w, n);           /* w = 4*(x1*x2-z1*z2)^2 */
  mpres_sqr (v, v, n);           /* v = 4*(x2*z1-x1*z2)^2 */

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
  mpres_sqr (u, u, n);      /* u = (x1+z1)^2 mod n */
  mpres_sub (v, x1, z1, n);
  mpres_sqr (v, v, n);      /* v = (x1-z1)^2 mod n */
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
      if (ecm_tstbit (e, l)) /* k, k+1 -> 2k+1, 2k+2 */
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

/* returns the number of modular multiplications for computing
   V_n from V_r * V_{n-r} - V_{n-2r}.
   ADD is the cost of an addition
   DUP is the cost of a duplicate
*/
static double
lucas_cost (ecm_uint n, double v)
{
  ecm_uint d, e, r;
  double c; /* cost */

  d = n;
  r = (ecm_uint) ((double) d * v + 0.5);
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
      if (d - e <= e / 4 && ((d + e) % 3) == 0)
        { /* condition 1 */
          d = (2 * d - e) / 3;
          e = (e - d) / 2;
          c += 3.0 * ADD; /* 3 additions */
        }
      else if (d - e <= e / 4 && (d - e) % 6 == 0)
        { /* condition 2 */
          d = (d - e) / 2;
          c += ADD + DUP; /* one addition, one duplicate */
        }
      else if ((d + 3) / 4 <= e)
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


/* computes kP from P=(xA:zA) and puts the result in (xA:zA). Assumes k>2. 
   WARNING! The calls to add3() assume that the two input points are distinct,
   which is not neccessarily satisfied. The result can be that in rare cases
   the point at infinity (z==0) results when it shouldn't. A test case is 
   echo 33554520197234177 | ./ecm -sigma 2046841451 373 1
   which finds the prime even though it shouldn't (23^2=529 divides order).
   This is not a problem for ECM since at worst we'll find a factor we 
   shouldn't have found. For other purposes (i.e. primality proving) this 
   would have to be fixed first.
*/

static void
prac (mpres_t xA, mpres_t zA, ecm_uint k, mpmod_t n, mpres_t b,
      mpres_t u, mpres_t v, mpres_t w, mpres_t xB, mpres_t zB, mpres_t xC, 
      mpres_t zC, mpres_t xT, mpres_t zT, mpres_t xT2, mpres_t zT2)
{
  ecm_uint d, e, r, i = 0, nv;
  double c, cmin;
  __mpz_struct *tmp;
#define NV 10  
  /* 1/val[0] = the golden ratio (1+sqrt(5))/2, and 1/val[i] for i>0
     is the real number whose continued fraction expansion is all 1s
     except for a 2 in i+1-st place */
  static double val[NV] =
    { 0.61803398874989485, 0.72360679774997897, 0.58017872829546410,
      0.63283980608870629, 0.61242994950949500, 0.62018198080741576,
      0.61721461653440386, 0.61834711965622806, 0.61791440652881789,
      0.61807966846989581};

  /* for small n, it makes no sense to try 10 different Lucas chains */
  nv = mpz_size ((mpz_ptr) n);
  if (nv > NV)
    nv = NV;

  if (nv > 1)
    {
      /* chooses the best value of v */
      for (d = 0, cmin = ADD * (double) k; d < nv; d++)
        {
          c = lucas_cost (k, val[d]);
          if (c < cmin)
            {
              cmin = c;
              i = d;
            }
        }
    }

  d = k;
  r = (ecm_uint) ((double) d * val[i] + 0.5);
  
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
      if (d - e <= e / 4 && ((d + e) % 3) == 0)
        { /* condition 1 */
          d = (2 * d - e) / 3;
          e = (e - d) / 2;
          add3 (xT, zT, xA, zA, xB, zB, xC, zC, n, u, v, w); /* T = f(A,B,C) */
          add3 (xT2, zT2, xT, zT, xA, zA, xB, zB, n, u, v, w); /* T2 = f(T,A,B) */
          add3 (xB, zB, xB, zB, xT, zT, xA, zA, n, u, v, w); /* B = f(B,T,A) */
          mpres_swap (xA, xT2, n);
          mpres_swap (zA, zT2, n); /* swap A and T2 */
        }
      else if (d - e <= e / 4 && (d - e) % 6 == 0)
        { /* condition 2 */
          d = (d - e) / 2;
          add3 (xB, zB, xA, zA, xB, zB, xC, zC, n, u, v, w); /* B = f(A,B,C) */
          duplicate (xA, zA, xA, zA, n, b, u, v, w); /* A = 2*A */
        }
      else if ((d + 3) / 4 <= e)
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

/* functions for using optimal (or near-optimal) Lucas chains
   instead of PRAC. Added by P McLaughlin 2023 */


/* To verify that Lucas chains are being generated
   correctly, configure with --enable-assert */

/* generate Lucas chain for a prime >= 2.
 *
 * inputs: (1) prime p (the final chain value)
 *         (2) chain_code for p read from input file (starting with p = 11)
 *
 * Outputs: (1) Lchain - an array of Lucas chain elements leading to p
 *          (2) chain length
*/
static
uint8_t generate_Lucas_chain( uint64_t prime, uint64_t chain_code, chain_element *Lchain )
{
	uint8_t code_fragment, chain_length, i, k;
	uint64_t dif;

	/* the code file starts at p = 11, so handle 2, 3, 5, or 7 separately */
	if( prime < 11 )
	{
		if(prime == 2)
		{
			chain_length = 1;
			return chain_length;
		}
		else if( prime == 3 )
		{
			chain_length = 2;
			return chain_length;
		}
		else if( prime == 5 )
		{
			Lchain[3].value = 5;
			Lchain[3].comp_offset_1 = 0;
			Lchain[3].comp_offset_2 = 1;
			Lchain[3].dif_offset = 2;

			chain_length = 3;
			return chain_length;
		}
		else if( prime == 7 )
		{
			Lchain[3].value = 4;
			Lchain[3].comp_offset_1 = 1;
			Lchain[3].comp_offset_2 = 1;
			Lchain[3].dif_offset = 0;

			Lchain[4].value = 7;
			Lchain[4].comp_offset_1 = 0;
			Lchain[4].comp_offset_2 = 1;
			Lchain[4].dif_offset = 3;

			chain_length = 4;
			return chain_length;
		}
		else
		{
			printf("ERROR: generate_Lucas_chain entered with prime = %lu < 11 but != 2, 3, 5 or 7\n", prime);
			return 0;
		}
	}

	/* first 3 bits of code give the next two or three chain components */
	code_fragment = (uint8_t)(chain_code & 0x7);
	chain_code >>= 3;
	switch( code_fragment )
	{
	case CHAIN_START_5_8_13:
		Lchain[3].value = 5;
		Lchain[3].comp_offset_1 = 0;
		Lchain[3].comp_offset_2 = 1;
		Lchain[3].dif_offset = 2;

		Lchain[4].value = 8;
		Lchain[4].comp_offset_1 = 0;
		Lchain[4].comp_offset_2 = 1;
		Lchain[4].dif_offset = 2;

		Lchain[5].value = 13;
		Lchain[5].comp_offset_1 = 0;
		Lchain[5].comp_offset_2 = 1;
		Lchain[5].dif_offset = 2;

		chain_length = 5;
		break;

	case CHAIN_START_5_8_11:
		Lchain[3].value = 5;
		Lchain[3].comp_offset_1 = 0;
		Lchain[3].comp_offset_2 = 1;
		Lchain[3].dif_offset = 2;

		Lchain[4].value = 8;
		Lchain[4].comp_offset_1 = 0;
		Lchain[4].comp_offset_2 = 1;
		Lchain[4].dif_offset = 2;

		Lchain[5].value = 11;
		Lchain[5].comp_offset_1 = 0;
		Lchain[5].comp_offset_2 = 2;
		Lchain[5].dif_offset = 1;

		chain_length = 5;
		break;

	case CHAIN_START_5_8_10:
		Lchain[3].value = 5;
		Lchain[3].comp_offset_1 = 0;
		Lchain[3].comp_offset_2 = 1;
		Lchain[3].dif_offset = 2;

		Lchain[4].value = 8;
		Lchain[4].comp_offset_1 = 0;
		Lchain[4].comp_offset_2 = 1;
		Lchain[4].dif_offset = 2;

		Lchain[5].value = 10;
		Lchain[5].comp_offset_1 = 1;
		Lchain[5].comp_offset_2 = 1;
		Lchain[5].dif_offset = 0;

		chain_length = 5;
		break;

	case CHAIN_START_5_7:
		Lchain[3].value = 5;
		Lchain[3].comp_offset_1 = 0;
		Lchain[3].comp_offset_2 = 1;
		Lchain[3].dif_offset = 2;

		Lchain[4].value = 7;
		Lchain[4].comp_offset_1 = 0;
		Lchain[4].comp_offset_2 = 2;
		Lchain[4].dif_offset = 1;

		chain_length = 4;
		break;

	case CHAIN_START_5_6:
		Lchain[3].value = 5;
		Lchain[3].comp_offset_1 = 0;
		Lchain[3].comp_offset_2 = 1;
		Lchain[3].dif_offset = 2;

		Lchain[4].value = 6;
		Lchain[4].comp_offset_1 = 1;
		Lchain[4].comp_offset_2 = 1;
		Lchain[4].dif_offset = 0;

		chain_length = 4;
		break;

	case CHAIN_START_4_7:
		Lchain[3].value = 4;
		Lchain[3].comp_offset_1 = 1;
		Lchain[3].comp_offset_2 = 1;
		Lchain[3].dif_offset = 0;

		Lchain[4].value = 7;
		Lchain[4].comp_offset_1 = 0;
		Lchain[4].comp_offset_2 = 1;
		Lchain[4].dif_offset = 3;

		chain_length = 4;
		break;

	case CHAIN_START_4_6:
		Lchain[3].value = 4;
		Lchain[3].comp_offset_1 = 1;
		Lchain[3].comp_offset_2 = 1;
		Lchain[3].dif_offset = 0;

		Lchain[4].value = 6;
		Lchain[4].comp_offset_1 = 1;
		Lchain[4].comp_offset_2 = 1;
		Lchain[4].dif_offset = 0;

		chain_length = 4;
		break;

	case CHAIN_START_4_5:
		Lchain[3].value = 4;
		Lchain[3].comp_offset_1 = 1;
		Lchain[3].comp_offset_2 = 1;
		Lchain[3].dif_offset = 0;

		Lchain[4].value = 5;
		Lchain[4].comp_offset_1 = 0;
		Lchain[4].comp_offset_2 = 3;
		Lchain[4].dif_offset = 1;

		chain_length = 4;
		break;

	default: /* should never happen */
		printf("ERROR: bad chain code start value in generate_Lchain = %u\n", code_fragment);
		return 0;
	}

	/* rebuild chain from code fragments */
	while( chain_code != 0 )
	{
		code_fragment = (uint8_t)( chain_code & 0xF );
		chain_code >>= 4;
		switch( code_fragment )
		{
			case 0: /* step type 1 or 4 */
			{
				i = (uint8_t)( chain_code & 0x3 );
				chain_code >>= 2;
				if( i == 3 )
				{
					 /* a rare double, but does occur */
					Lchain[ chain_length+1 ].value = 2*Lchain[ chain_length-2 ].value;
					Lchain[ chain_length+1 ].comp_offset_1 = 2;
					Lchain[ chain_length+1 ].comp_offset_2 = 2;
					Lchain[ chain_length+1 ].dif_offset = 0;
					chain_length++;
				}
				else
				{
					i = 12*(i + 1);
					max_continuation( Lchain, &chain_length, i );
				}
				break;
			}
			case 1:
			{
				Lchain[ chain_length+1 ].value = Lchain[ chain_length ].value + Lchain[ chain_length-2 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 0;
				Lchain[ chain_length+1 ].comp_offset_2 = 2;
				dif = Lchain[ chain_length ].value - Lchain[ chain_length-2 ].value;
				k = 1;
				while( dif < Lchain[ chain_length-k ].value )
					k++;
				ASSERT( (dif == Lchain[ chain_length-k ].value) && (k < 15) );

				Lchain[ chain_length+1 ].dif_offset = k;
				chain_length++;
				break;
			}
			case 2:
			{
				Lchain[ chain_length+1 ].value = 2*Lchain[ chain_length-1 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 1;
				Lchain[ chain_length+1 ].comp_offset_2 = 1;
				Lchain[ chain_length+1 ].dif_offset = 0;
				chain_length++;
				break;
			}
			case 3:
			{
				i = 1;
				max_continuation( Lchain, &chain_length, i );
				Lchain[ chain_length+1 ].value = Lchain[ chain_length ].value + Lchain[ chain_length-2 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 0;
				Lchain[ chain_length+1 ].comp_offset_2 = 2;
				Lchain[ chain_length+1 ].dif_offset = 1;
				chain_length++;
				break;
			}
			case 4:
			{
				i = 1;
				max_continuation( Lchain, &chain_length, i );
				Lchain[ chain_length+1 ].value = 2*Lchain[ chain_length-1 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 1;
				Lchain[ chain_length+1 ].comp_offset_2 = 1;
				Lchain[ chain_length+1 ].dif_offset = 0;
				chain_length++;
				break;
			}
			case 5:
			{
				i = 2;
				max_continuation( Lchain, &chain_length, i );
				Lchain[ chain_length+1 ].value = Lchain[ chain_length ].value + Lchain[ chain_length-2 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 0;
				Lchain[ chain_length+1 ].comp_offset_2 = 2;
				Lchain[ chain_length+1 ].dif_offset = 1;
				chain_length++;
				break;
			}
			case 6:
			{
				i = 2;
				max_continuation( Lchain, &chain_length, i );
				Lchain[ chain_length+1 ].value = 2*Lchain[ chain_length-1 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 1;
				Lchain[ chain_length+1 ].comp_offset_2 = 1;
				Lchain[ chain_length+1 ].dif_offset = 0;
				chain_length++;
				break;
			}
			case 7:
			{
				i = 3;
				max_continuation( Lchain, &chain_length, i );
				Lchain[ chain_length+1 ].value = Lchain[ chain_length ].value + Lchain[ chain_length-2 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 0;
				Lchain[ chain_length+1 ].comp_offset_2 = 2;
				Lchain[ chain_length+1 ].dif_offset = 1;
				chain_length++;
				break;
			}
			case 8:
			{
				i = 3;
				max_continuation( Lchain, &chain_length, i );
				Lchain[ chain_length+1 ].value = 2*Lchain[ chain_length-1 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 1;
				Lchain[ chain_length+1 ].comp_offset_2 = 1;
				Lchain[ chain_length+1 ].dif_offset = 0;
				chain_length++;
				break;
			}
			case 9:
			{
				i = (uint8_t)( chain_code & 0x3 );
				chain_code >>= 2;
				i += 4;
				max_continuation( Lchain, &chain_length, i );
				Lchain[ chain_length+1 ].value = Lchain[ chain_length ].value + Lchain[ chain_length-2 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 0;
				Lchain[ chain_length+1 ].comp_offset_2 = 2;
				Lchain[ chain_length+1 ].dif_offset = 1;
				chain_length++;
				break;
			}
			case 10:
			{
				i = (uint8_t)( chain_code & 0x3 );
				chain_code >>= 2;
				i += 4;
				max_continuation( Lchain, &chain_length, i );
				Lchain[ chain_length+1 ].value = 2*Lchain[ chain_length-1 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 1;
				Lchain[ chain_length+1 ].comp_offset_2 = 1;
				Lchain[ chain_length+1 ].dif_offset = 0;
				chain_length++;
				break;
			}
			case 11:
			{
				i = (uint8_t)( chain_code & 0x3 );
				chain_code >>= 2;
				i += 8;
				max_continuation( Lchain, &chain_length, i );
				Lchain[ chain_length+1 ].value = Lchain[ chain_length ].value + Lchain[ chain_length-2 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 0;
				Lchain[ chain_length+1 ].comp_offset_2 = 2;
				Lchain[ chain_length+1 ].dif_offset = 1;
				chain_length++;
				break;
			}
			case 12:
			{
				i = (uint8_t)( chain_code & 0x3 );
				chain_code >>= 2;
				i += 8;
				max_continuation( Lchain, &chain_length, i );
				Lchain[ chain_length+1 ].value = 2*Lchain[ chain_length-1 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 1;
				Lchain[ chain_length+1 ].comp_offset_2 = 1;
				Lchain[ chain_length+1 ].dif_offset = 0;
				chain_length++;
				break;
			}
			case 13:
			{
				i = (uint8_t)( chain_code & 0x3 );
				chain_code >>= 2;

				switch( i )
				{
				case 0: /* chain encoding did not require this step until chain length = 44 */
					Lchain[ chain_length+1 ].value = Lchain[ chain_length ].value + Lchain[ chain_length-7 ].value;
					Lchain[ chain_length+1 ].comp_offset_1 = 0;
					Lchain[ chain_length+1 ].comp_offset_2 = 7;
					dif = Lchain[ chain_length ].value - Lchain[ chain_length-7 ].value;
					k = 1;
					while( dif < Lchain[ chain_length-k ].value )
						k++;
					ASSERT( (dif == Lchain[ chain_length-k ].value) && (k < 15) );

					Lchain[ chain_length+1 ].dif_offset = k;
					chain_length++;
					break;

				case 1:
					Lchain[ chain_length+1 ].value = Lchain[ chain_length ].value + Lchain[ chain_length-8 ].value;
					Lchain[ chain_length+1 ].comp_offset_1 = 0;
					Lchain[ chain_length+1 ].comp_offset_2 = 8;
					dif = Lchain[ chain_length ].value - Lchain[ chain_length-8 ].value;
					k = 1;
					while( dif < Lchain[ chain_length-k ].value )
						k++;
					ASSERT( (dif == Lchain[ chain_length-k ].value) && (k < 15) );

					Lchain[ chain_length+1 ].dif_offset = k;
					chain_length++;
					break;

				case 2: /* example chain 1 2 3 5 8 10 12 13 */
					Lchain[ chain_length+1 ].value = Lchain[ chain_length - 2 ].value + Lchain[ chain_length - 3 ].value;
					Lchain[ chain_length+1 ].comp_offset_1 = 2;
					Lchain[ chain_length+1 ].comp_offset_2 = 3;
					dif = Lchain[ chain_length - 2 ].value - Lchain[ chain_length - 3 ].value;
					k = 1;
					while( dif < Lchain[ chain_length-k ].value )
						k++;
					ASSERT( (dif == Lchain[ chain_length-k ].value) && (k < 15) );

					Lchain[ chain_length+1 ].dif_offset = k;
					chain_length++;
					break;

				case 3: /* can occur but don't have an example (yet) */
					Lchain[ chain_length+1 ].value = Lchain[ chain_length - 2 ].value + Lchain[ chain_length - 4 ].value;
					Lchain[ chain_length+1 ].comp_offset_1 = 2;
					Lchain[ chain_length+1 ].comp_offset_2 = 4;
					dif = Lchain[ chain_length - 2 ].value - Lchain[ chain_length - 4 ].value;
					k = 1;
					while( dif < Lchain[ chain_length-k ].value )
						k++;
					ASSERT( (dif == Lchain[ chain_length-k ].value) && (k < 15) );

					Lchain[ chain_length+1 ].dif_offset = k;
					chain_length++;
					break;

				default: /* should never happen */
					printf("ERROR: unimplemented code fragment 0xiD, i = %u\n", i);
					return 0;
				}
				break;
			}
			case 14:
			{
				i = (uint8_t)( chain_code & 0x3 );
				chain_code >>= 2;

				Lchain[ chain_length+1 ].value = Lchain[ chain_length ].value + Lchain[ chain_length-3-i ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 0;
				Lchain[ chain_length+1 ].comp_offset_2 = 3+i;
				dif = Lchain[ chain_length ].value - Lchain[ chain_length-3-i ].value;
				k = 1;
				while( dif < Lchain[ chain_length-k ].value )
					k++;
				ASSERT( (dif == Lchain[ chain_length-k ].value) && (k < 15) );

				Lchain[ chain_length+1 ].dif_offset = k;
				chain_length++;
				break;
			}
			case 15:
			{
				i = (uint8_t)( chain_code & 0x3 );
				chain_code >>= 2;

				Lchain[ chain_length+1 ].value = Lchain[ chain_length-1 ].value + Lchain[ chain_length-2-i ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 1;
				Lchain[ chain_length+1 ].comp_offset_2 = 2+i;
				dif = Lchain[ chain_length-1 ].value - Lchain[ chain_length-2-i ].value;
				k = 2;
				while( dif < Lchain[ chain_length-k ].value )
					k++;
				ASSERT( (dif == Lchain[ chain_length-k ].value) && (k < 15) );

				Lchain[ chain_length+1 ].dif_offset = k;
				chain_length++;
				break;
			}
		}
	}

	/* finish the chain with maximum continuations until prime is reached */
	i = 1;
	while( Lchain[chain_length].value < prime )
	{
		max_continuation( Lchain, &chain_length, i );
	}

	ASSERT( Lchain[ chain_length ].value == prime );

	return chain_length;
}

/* extend the chain with maximum elements for i steps */
static
void max_continuation( chain_element *Lchain, uint8_t *chain_length, uint8_t i )
{
	uint8_t k;
	uint64_t dif;

	Lchain[ *chain_length+1 ].value = Lchain[ *chain_length ].value + Lchain[ *chain_length-1 ].value;
	Lchain[ *chain_length+1 ].comp_offset_1 = 0;
	Lchain[ *chain_length+1 ].comp_offset_2 = 1;
	dif = Lchain[ *chain_length ].value - Lchain[ *chain_length-1 ].value;
	k = 2;
	while( dif < Lchain[ *chain_length-k ].value )
		k++;
	ASSERT( (dif == Lchain[ *chain_length-k ].value) && (k < 15) );

	Lchain[ *chain_length+1 ].dif_offset = k;
	(*chain_length)++;
	if( i > 1 )
	{
		for( k = 1; k < i; k++)
		{
			Lchain[ *chain_length+1 ].value = Lchain[ *chain_length ].value + Lchain[ *chain_length-1 ].value;
			Lchain[ *chain_length+1 ].comp_offset_1 = 0;
			Lchain[ *chain_length+1 ].comp_offset_2 = 1;
			Lchain[ *chain_length+1 ].dif_offset = 2;
			(*chain_length)++;
		}
	}
}

/* Input: x is initial point
          A is curve parameter in Montgomery's form:
          g*y^2*z = x^3 + a*x^2*z + x*z^2
          n is the number to factor
	  B1 is the stage 1 bound
   Output: If a factor is found, it is returned in f.
           Otherwise, x contains the x-coordinate of the point computed
           in stage 1 (with z coordinate normalized to 1).
	   B1done is set to B1 if stage 1 completed normally,
	   or to the largest prime processed if interrupted, but never
	   to a smaller value than B1done was upon function entry.
   Return value: ECM_FACTOR_FOUND_STEP1 if a factor, otherwise 
           ECM_NO_FACTOR_FOUND
*/
static int
ecm_stage1 (mpz_t f, mpres_t x, mpres_t A, mpmod_t n, double B1, 
            double *B1done, mpz_t go, int (*stop_asap)(void), 
            char *chkfilename)
{
  mpres_t b, z, u, v, w, xB, zB, xC, zC, xT, zT, xT2, zT2;
  uint64_t p, r, last_chkpnt_p;
  int ret = ECM_NO_FACTOR_FOUND;
  long last_chkpnt_time;
  prime_info_t prime_info;

  /* Lucas chain code file mods */
  uint64_t chain_code;
  uint8_t dum, chain_length;
  int32_t i;
  FILE *chain_code_file = NULL;
  chain_element Lchain[64];
  uint8_t using_code_file; /* logical */

  /* Elliptic curve states as we follow the Lucas chain */
  mpres_t LCS_x[16];
  mpres_t LCS_z[16];
  uint8_t base_indx, next_indx, s1_indx, s2_indx, dif_indx;
  /* end mods */

  prime_info_init (prime_info);

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
  
  /* Lucas chain code file mods */
  for(i = 0; i < 16; i++)
  {
    mpres_init ( LCS_x[i], n);
    mpres_init ( LCS_z[i], n);
  }
  using_code_file = 0;

  if( B1 > *B1done )
  {
    chain_code_file = fopen("Lchain_codes.dat", "r");
    if(chain_code_file != (FILE *)NULL )
    {
      using_code_file = 1;
/*      outputf (OUTPUT_NORMAL, "Using Lucas chain codes\n"); */

      /* the first 3 chain elements are always value = 1, 2, and 3, respectively */
	  Lchain[0].value = 1;

      Lchain[1].value = 2;
      Lchain[1].comp_offset_1 = 0;
      Lchain[1].comp_offset_2 = 0;
      Lchain[1].dif_offset = 0;

      Lchain[2].value = 3;
      Lchain[2].comp_offset_1 = 0;
      Lchain[2].comp_offset_2 = 1;
      Lchain[2].dif_offset = 1;
    }
    else
    {
      outputf (OUTPUT_NORMAL, "Lchain_codes.dat file failed to open, using prac\n");
    }
  }

  chain_length = 0; /* not necessary but stops a compiler warning */
  /* end mods */

  last_chkpnt_time = cputime ();

  mpres_set_ui (z, 1, n);

  mpres_add_ui (b, A, 2, n);
  mpres_div_2exp (b, b, 2, n); /* b == (A0+2)*B/4 */

  /* preload group order */
  if (go != NULL)
    ecm_mul (x, z, go, n, b);

  /* prac() wants multiplicands > 2 */
  for (r = 2; r <= B1; r *= 2)
    if (r > *B1done)
      duplicate (x, z, x, z, n, b, u, v, w);
  
  /* We'll do 3 manually, too (that's what ecm4 did..) */
  for (r = 3; r <= B1; r *= 3)
    if (r > *B1done)
      {
        duplicate (xB, zB, x, z, n, b, u, v, w);
        add3 (x, z, x, z, xB, zB, x, z, n, u, v, w);
      }
  
  if( using_code_file )
  {
	  /* copy x & z to LCS_x[0] & lCS_z[0] */
	  mpres_set_z (LCS_x[0], x, n);
	  mpres_set_z (LCS_z[0], z, n);
	  base_indx = 0;
  }

   last_chkpnt_p = 3;
  p = getprime_mt (prime_info); /* Puts 3 into p. Next call gives 5 */
  for (p = getprime_mt (prime_info); p <= B1; p = getprime_mt (prime_info))
    {
      if(using_code_file)
      {
        if(p >= 11) /* code file starts at p = 11 */
        dum = fread((int8_t *)&chain_code, sizeof(uint64_t), 1, chain_code_file);

        if(dum || (p < 11))
        {
          /* the following call doesn't need a chain code for p < 11 */
          chain_length = generate_Lucas_chain( p, chain_code, Lchain );
        }
        else /* reached end of code file; revert to use prac */
        {
          /* copy final LCS_x,z values back to x,z */
          mpres_set_z (x, LCS_x[base_indx], n);
          mpres_set_z (z, LCS_z[base_indx], n);

          fclose(chain_code_file);
          using_code_file = 0;
          outputf (OUTPUT_NORMAL, "Reached Lchain_codes.dat EOF at p = %lu, reverting to use prac\n", p);
        }
      } /* end if( using_code_file) */

       for (r = p; r <= B1; r *= p)
       {
         if (r > *B1done)
         {
           if(using_code_file)
           {
             next_indx = (base_indx + 1) & 0xF;
             duplicate (LCS_x[next_indx], LCS_z[next_indx],
                        LCS_x[base_indx], LCS_z[base_indx], n, b, u, v, w);
             base_indx = next_indx;

             for(i=2; i<=chain_length; i++)
             {
               next_indx = (base_indx + 1) & 0xF;
               s1_indx = (base_indx - Lchain[i].comp_offset_1) & 0xF;
               if(Lchain[i].dif_offset == 0)
               {
                 duplicate (LCS_x[next_indx], LCS_z[next_indx],
                 LCS_x[s1_indx], LCS_z[s1_indx], n, b, u, v, w);
               }
               else
               {
                 s2_indx = (base_indx - Lchain[i].comp_offset_2) & 0xF;
                 dif_indx = (base_indx - Lchain[i].dif_offset) & 0xF;
                 add3 (LCS_x[next_indx], LCS_z[next_indx], LCS_x[s1_indx],
                       LCS_z[s1_indx], LCS_x[s2_indx], LCS_z[s2_indx],
                       LCS_x[dif_indx], LCS_z[dif_indx], n, u, v, w);
               }
               base_indx = next_indx;
             }
           }
           else
             prac (x, z, (ecm_uint) p, n, b, u, v, w, xB, zB, xC, zC, xT, zT, xT2, zT2);
         }
 
         if(using_code_file)
         {
           if (mpres_is_zero (LCS_z[base_indx], n))
           {
             outputf (OUTPUT_VERBOSE, "Reached point at infinity, %.0f divides "
                      "group orders\n", p);
             break;
           }
         }
         else
         {
           if (mpres_is_zero (z, n))
           {
             outputf (OUTPUT_VERBOSE, "Reached point at infinity, %.0f divides "
                      "group orders\n", p);
             break;
           }
         }

         if (stop_asap != NULL && (*stop_asap) ())
         {
           outputf (OUTPUT_NORMAL, "Interrupted at prime %.0f\n", p);
           break;
         }

         if (chkfilename != NULL && p > last_chkpnt_p + 10000 &&
             elltime (last_chkpnt_time, cputime ()) > CHKPNT_PERIOD)
         {
           if(using_code_file) 
             writechkfile (chkfilename, ECM_ECM, MAX(p, *B1done), n, A, LCS_x[base_indx], NULL, LCS_z[base_indx]);
           else
             writechkfile (chkfilename, ECM_ECM, MAX(p, *B1done), n, A, x, NULL, z);
           last_chkpnt_p = p;
           last_chkpnt_time = cputime ();
         }
       }
    }
  
  /* If stage 1 finished normally, p is the smallest prime >B1 here.
     In that case, set to B1 */
  if (p > B1)
      p = B1;
  
  if (p > *B1done)
      *B1done = p;

  if(using_code_file)
  {
    /* copy final LCS_x,z values back to x,z */
    mpres_set_z (x, LCS_x[base_indx], n);
    mpres_set_z (z, LCS_z[base_indx], n);

    fclose(chain_code_file);
  }

  if (chkfilename != NULL)
    writechkfile (chkfilename, ECM_ECM, *B1done, n, A, x, NULL, z);

  prime_info_clear (prime_info);

  if (!mpres_invert (u, z, n)) /* Factor found? */
    {
      mpres_gcd (f, z, n);
      ret = ECM_FACTOR_FOUND_STEP1;
    }
  mpres_mul (x, x, u, n);

  for(i = 15; i >= 0; i--)
  {
    mpres_clear ( LCS_z[i], n);
    mpres_clear ( LCS_x[i], n);
  }

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

#define DEBUG_EC_W 0

#ifdef HAVE_ADDLAWS
/* Input: when Etype == ECM_EC_TYPE_WEIERSTRASS*:
            (x, y) is initial point
            A is curve parameter in Weierstrass's form:
            Y^2 = X^3 + A*X + B, where B = y^2-(x^3+A*x) is implicit
	  when Etype == ECM_EC_TYPE_TWISTED_HESSIAN:
	    (x, y) is initial point
	    A=a/d is curve parameter in Hessian form: a*X^3+Y^3+Z^3=d*X*Y*Z
          n is the number to factor
	  B1 is the stage 1 bound
	  batch_s = prod(p^e <= B1) if != 1
   Output: If a factor is found, it is returned in f.
           Otherwise, (x, y) contains the point computed in stage 1.
	   B1done is set to B1 if stage 1 completed normally,
	   or to the largest prime processed if interrupted, but never
	   to a smaller value than B1done was upon function entry.
   Return value: ECM_FACTOR_FOUND_STEP1 if a factor is found, otherwise 
           ECM_NO_FACTOR_FOUND
*/
static int
ecm_stage1_W (mpz_t f, ell_curve_t E, ell_point_t P, mpmod_t n, 
	      double B1, double *B1done, mpz_t batch_s, mpz_t go, 
	      int (*stop_asap)(void), char *chkfilename)
{
    mpres_t xB;
    ell_point_t Q;
    uint64_t p = 0, r, last_chkpnt_p;
    int ret = ECM_NO_FACTOR_FOUND, status;
    long last_chkpnt_time;
    prime_info_t prime_info;

    prime_info_init (prime_info);
    
    mpres_init (xB, n);

    ell_point_init(Q, E, n);
    
    last_chkpnt_time = cputime ();

#if DEBUG_EC_W >= 2
    gmp_printf("N:=%Zd;\n", n->orig_modulus);
    printf("E:="); ell_curve_print(E, n);
    printf("E:=[E[4], E[5]];\n");
    printf("P:="); ell_point_print(P, E, n); printf("; Q:=P;\n");
#endif
    /* preload group order */
    if (go != NULL){
	if (ell_point_mul (f, Q, go, P, E, n) == 0){
	    ret = ECM_FACTOR_FOUND_STEP1;
	    goto end_of_stage1_w;
        }
	ell_point_set(P, Q, E, n);
    }
#if DEBUG_EC_W >= 1
    gmp_printf("go:=%Zd;\n", go);
    printf("goP:="); ell_point_print(P, E, n); printf(";\n");
#endif
    if(mpz_cmp_ui(batch_s, 1) == 0){
        outputf (OUTPUT_VERBOSE, "Using traditional approach to Step 1\n");
	for (r = 2; r <= B1; r *= 2)
	    if (r > *B1done){
		if(ell_point_duplicate (f, Q, P, E, n) == 0){
		    ret = ECM_FACTOR_FOUND_STEP1;
		    goto end_of_stage1_w;
		}
		ell_point_set(P, Q, E, n);
#if DEBUG_EC_W >= 2
		printf("P%ld:=", (long)r); ell_point_print(P, E, n); printf(";\n");
		printf("Q:=EcmMult(2, Q, E, N);\n");
		printf("(Q[1]*P%ld[3]-Q[3]*P%ld[1]) mod N;\n",(long)r,(long)r);
#endif
	    }
	
	last_chkpnt_p = 3;
	for (p = getprime_mt (prime_info); p <= B1; p = getprime_mt (prime_info)){
	    for (r = p; r <= B1; r *= p){
		if (r > *B1done){
		    mpz_set_ui(f, (ecm_uint)p);
		    status = ell_point_mul (f, Q, f, P, E, n);
		    if(status == 0){
		    }
		    else if(E->law == ECM_LAW_HOMOGENEOUS){
			if(E->type == ECM_EC_TYPE_TWISTED_HESSIAN)
			    mpres_gcd(f, Q->x, n);
			else
			    mpres_gcd(f, Q->z, n);
			//			gmp_printf("gcd=%Zd\n", f);
			if(mpz_cmp(f, n->orig_modulus) < 0 
			   && mpz_cmp_ui(f, 1) > 0)
			    status = 0;
		    }
		    if(status == 0){
			ret = ECM_FACTOR_FOUND_STEP1;
			goto end_of_stage1_w;
		    }
#if DEBUG_EC_W >= 2
		    printf("R%ld:=", (long)r); ell_point_print(Q, E, n);
		    printf(";\nQ:=EcmMult(%ld, Q, E, N);\n", (long)p);
		    printf("(Q[1]*R%ld[3]-Q[3]*R%ld[1]) mod N;\n",(long)r,(long)r);
#endif
		    ell_point_set(P, Q, E, n);
		}
	    }
	    if (ell_point_is_zero (P, E, n)){
		outputf (OUTPUT_VERBOSE, "Reached point at infinity, "
			 "%.0f divides group orders\n", p);
		break;
	    }
	    
	    if (stop_asap != NULL && (*stop_asap) ()){
		outputf (OUTPUT_NORMAL, "Interrupted at prime %.0f\n", p);
		break;
	    }
	    
	    if (chkfilename != NULL && p > last_chkpnt_p + 10000 && 
		elltime (last_chkpnt_time, cputime ()) > CHKPNT_PERIOD){
		writechkfile (chkfilename, ECM_ECM, MAX(p, *B1done), 
			      n, E->a4, P->x, P->y, P->z);
		last_chkpnt_p = p;
		last_chkpnt_time = cputime ();
	    }
	}
    }
    else{
#if USE_ADD_SUB_CHAINS == 0 /* keeping it simple */
	if (ell_point_mul (f, Q, batch_s, P, E, n) == 0){
	    ret = ECM_FACTOR_FOUND_STEP1;
	    goto end_of_stage1_w;
        }
#else
	/* batch mode and special coding... */
	short *S = NULL;
	size_t iS;
	int w;
	add_sub_unpack(&w, &S, &iS, batch_s);
	if (ell_point_mul_add_sub_with_S(f, Q, P, E, n, w, S, iS) == 0){
	    ret = ECM_FACTOR_FOUND_STEP1;
        }
#endif
	ell_point_set(P, Q, E, n);
	p = B1;
    }
 end_of_stage1_w:
    /* If stage 1 finished normally, p is the smallest prime > B1 here.
       In that case, set to B1 */
    if (p > B1)
	p = B1;
    
    if (p > *B1done)
	*B1done = p;
    
    if (chkfilename != NULL)
	writechkfile (chkfilename, ECM_ECM, *B1done, n, E->a4, P->x, P->y,P->z);
    prime_info_clear (prime_info);
#if DEBUG_EC_W >= 2
    printf("lastP="); ell_point_print(P, E, n); printf("\n");
#endif
    if(ret != ECM_FACTOR_FOUND_STEP1){
	if(ell_point_is_zero(P, E, n) == 1){
	    /* too bad */
	    ell_point_set_to_zero(P, E, n);
	    mpz_set(f, n->orig_modulus);
	    ret = ECM_FACTOR_FOUND_STEP1;
	}
	else{
	    /* for affine case, z = 1 anyway */
	    if(E->law == ECM_LAW_HOMOGENEOUS){
		if (!mpres_invert (xB, P->z, n)){ /* Factor found? */
		    mpres_gcd (f, P->z, n);
		    gmp_printf("# factor found during normalization: %Zd\n", f);
		    ret = ECM_FACTOR_FOUND_STEP1;
		}
		else{
		    /* normalize to get (x:y:1) valid in W or H form... */
#if DEBUG_EC_W >= 2
		    mpres_get_z(f, xB, n); gmp_printf("1/z=%Zd\n", f);
#endif
		    mpres_mul (P->x, P->x, xB, n);
		    mpres_mul (P->y, P->y, xB, n);
#if DEBUG_EC_W >= 2
		    mpres_get_z(f, P->x, n); gmp_printf("x/z=%Zd\n", f);
		    mpres_get_z(f, P->y, n); gmp_printf("y/z=%Zd\n", f);
#endif
		    mpres_set_ui (P->z, 1, n);
		}
	    }
	}
    }

    mpres_clear (xB, n);
    ell_point_clear(Q, E, n);
    
    return ret;
}
#endif

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

#define DIGITS_START 35
#define DIGITS_INCR   5
#define DIGITS_END   80

void
print_expcurves (double B1, const mpz_t B2, unsigned long dF, unsigned long k, 
                 int S, int param)
{
  double prob;
  int i, j;
  char sep, outs[128], flt[16];
  double smoothness_correction;

  if (param == ECM_PARAM_SUYAMA || param == ECM_PARAM_BATCH_2)
      smoothness_correction = 1.0; 
  else if (param == ECM_PARAM_BATCH_SQUARE)
      smoothness_correction = EXTRA_SMOOTHNESS_SQUARE;
  else if (param == ECM_PARAM_BATCH_32BITS_D)
      smoothness_correction = EXTRA_SMOOTHNESS_32BITS_D;
  else /* This case should never happen */
      smoothness_correction = 0.0; 

  for (i = DIGITS_START, j = 0; i <= DIGITS_END; i += DIGITS_INCR)
    j += sprintf (outs + j, "%u%c", i, (i < DIGITS_END) ? '\t' : '\n');
  outs[j] = '\0';
  outputf (OUTPUT_VERBOSE, "Expected number of curves to find a factor "
           "of n digits (assuming one exists):\n%s", outs);
  for (i = DIGITS_START; i <= DIGITS_END; i += DIGITS_INCR)
    {
      sep = (i < DIGITS_END) ? '\t' : '\n';
      prob = ecmprob (B1, mpz_get_d (B2),
                      /* smoothness depends on the parametrization */
                      pow (10., i - .5) / smoothness_correction,
                      (double) dF * dF * k, S);
      if (prob > 1. / 10000000)
        outputf (OUTPUT_VERBOSE, "%.0f%c", floor (1. / prob + .5), sep);
      else if (prob > 0.)
        {
          /* on Windows: 2.6e+009   or   2.6e+025  (8 characters in string) */
          /* on Linux  : 2.6e+09    or   2.6e+25   (7 characters in string) */
          /* This will normalize the output so that the Windows ouptut will match the Linux output */
          if (sprintf (flt, "%.2g", floor (1. / prob + .5)) == 8)
            memmove (&flt[5], &flt[6], strlen(flt) - 5);
          outputf (OUTPUT_VERBOSE, "%s%c", flt, sep);
        }
      else
        outputf (OUTPUT_VERBOSE, "Inf%c", sep);
    }
}

void
print_exptime (double B1, const mpz_t B2, unsigned long dF, unsigned long k, 
               int S, double tottime, int param)
{
  double prob, exptime;
  int i, j;
  char sep, outs[128];
  double smoothness_correction;

  if (param == ECM_PARAM_SUYAMA || param == ECM_PARAM_BATCH_2)
      smoothness_correction = 1.0; 
  else if (param == ECM_PARAM_BATCH_SQUARE)
      smoothness_correction = EXTRA_SMOOTHNESS_SQUARE;
  else if (param == ECM_PARAM_BATCH_32BITS_D)
      smoothness_correction = EXTRA_SMOOTHNESS_32BITS_D;
  else /* This case should never happen */
      smoothness_correction = 0.0; 
  
  for (i = DIGITS_START, j = 0; i <= DIGITS_END; i += DIGITS_INCR)
    j += sprintf (outs + j, "%u%c", i, (i < DIGITS_END) ? '\t' : '\n');
  outs[j] = '\0';
  outputf (OUTPUT_VERBOSE, "Expected time to find a factor of n digits:\n%s",
           outs);
  for (i = DIGITS_START; i <= DIGITS_END; i += DIGITS_INCR)
    {
      sep = (i < DIGITS_END) ? '\t' : '\n';
      prob = ecmprob (B1, mpz_get_d (B2),
                      /* in batch mode, the extra smoothness is smaller */
                      pow (10., i - .5) / smoothness_correction,
                      (double) dF * dF * k, S);
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
}

/* y should be NULL for P+1, and P-1, it contains the y coordinate for the
   Weierstrass form for ECM (when sigma_is_A = -1). */
void
print_B1_B2_poly (int verbosity, int method, double B1, double B1done, 
		  const mpz_t B2min_param, const mpz_t B2min, const mpz_t B2,
                  int S, const mpz_t sigma, int sigma_is_A, int Etype,
		  const mpz_t y, int param, unsigned int nb_curves)
{
  ASSERT ((method == ECM_ECM) || (y == NULL));
  ASSERT ((-1 <= sigma_is_A) && (sigma_is_A <= 1));
  ASSERT (param != ECM_PARAM_DEFAULT || sigma_is_A == 1 || sigma_is_A == -1);

  if (test_verbose (verbosity))
  {
      outputf (verbosity, "Using ");
      if (ECM_IS_DEFAULT_B1_DONE(B1done))
	  outputf (verbosity, "B1=%1.0f, ", B1);
      else
	  outputf (verbosity, "B1=%1.0f-%1.0f, ", B1done, B1);
      if (mpz_sgn (B2min_param) < 0)
	  outputf (verbosity, "B2=%Zd", B2);
      else
	  outputf (verbosity, "B2=%Zd-%Zd", B2min, B2);
      
      if (S > 0)
	  outputf (verbosity, ", polynomial x^%u", S);
      else if (S < 0)
	  outputf (verbosity, ", polynomial Dickson(%u)", -S);
      
      /* don't print in resume case, since x0 is saved in resume file */
      if (method == ECM_ECM)
        {
	    if (sigma_is_A == 1)
		outputf (verbosity, ", A=%Zd", sigma);
	    else if (sigma_is_A == 0)
	      {
		if (nb_curves > 1) 
		  {
                    mpz_t end;
                    mpz_init(end);
		    mpz_add_ui (end, sigma, nb_curves-1);
		    outputf (verbosity, ", sigma=%d:%Zd-%Zd", param, sigma, end);
                    mpz_clear(end);
		    outputf (verbosity, " (%u curves)", nb_curves);
		  }
		else
		    outputf (verbosity, ", sigma=%d:%Zd", param, sigma);
	      }
	    else{
		if (Etype == ECM_EC_TYPE_WEIERSTRASS)
		  outputf (verbosity, ", Weierstrass(A=%Zd,y=%Zd)", sigma, y);
		else if (Etype == ECM_EC_TYPE_HESSIAN)
		  outputf (verbosity, ", Hessian(D=%Zd,y=%Zd)", sigma, y);
		else if (Etype == ECM_EC_TYPE_TWISTED_HESSIAN)
		    outputf (verbosity, ", twisted Hessian(y=%Zd)", y);
	    }
        }
      else if (ECM_IS_DEFAULT_B1_DONE(B1done))
        /* in case of P-1 or P+1, we store the initial point in sigma */
        outputf (verbosity, ", x0=%Zd", sigma);
      
      outputf (verbosity, "\n");
  }
}

/* Compute parameters for stage 2 */
int
set_stage_2_params (mpz_t B2, const mpz_t B2_parm,
                    mpz_t B2min, const mpz_t B2min_parm,
                    root_params_t *root_params, double B1,
                    unsigned long *k, const int S, int use_ntt,
                    unsigned long *dF, char *TreeFilename, double maxmem, 
                    mpmod_t modulus)
{
  mpz_set (B2min, B2min_parm);
  mpz_set (B2, B2_parm);
  
  mpz_init (root_params->i0);

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
    mpz_set_d (B2, pow (ECM_COST * B1, DEFAULT_B2_EXPONENT));

  /* set B2min */
  if (mpz_sgn (B2min) < 0)
    mpz_set_d (B2min, B1);

  /* Let bestD determine parameters for root generation and the 
     effective B2 */

  root_params->d2 = 0; /* Enable automatic choice of d2 */
  if (bestD (root_params, k, dF, B2min, B2, use_ntt, use_ntt, maxmem,
             (TreeFilename != NULL), modulus) == ECM_ERROR)
    return ECM_ERROR;

  /* Set default degree for Brent-Suyama extension */
  /* We try to keep the time used by the Brent-Suyama extension
     at about 10% of the stage 2 time */
  /* Degree S Dickson polys and x^S are equally fast for ECM, so we go for
     the better Dickson polys whenever possible. For S == 1, 2, they behave
     identically. */

  root_params->S = S;
  if (root_params->S == ECM_DEFAULT_S)
    {
      mpz_t t;
      mpz_init (t);
      mpz_sub (t, B2, B2min);
      root_params->S = choose_S (t);
      mpz_clear (t);
    }
  return ECM_NO_FACTOR_FOUND;
}

/* Input: x is starting point or zero
          y is used for Weierstrass curves (even in Step1)
          sigma is sigma value (if x is set to zero) or 
            A parameter (if x is non-zero) of curve
          n is the number to factor
          go is the initial group order to preload  
          B1, B2 are the stage 1/stage 2 bounds, respectively
          B2min the lower bound for stage 2
          k is the number of blocks to do in stage 2
          S is the degree of the Suyama-Brent extension for stage 2
          verbose is verbosity level: 0 no output, 1 normal output,
            2 diagnostic output.
	  sigma_is_A: If true, the sigma parameter contains the curve's A value
	  Etype
	  zE is a curve that is used when a special torsion group was used; in
	    that case, (x, y) must be a point on E.
   Output: f is the factor found.
   Return value: ECM_FACTOR_FOUND_STEPn if a factor was found,
                 ECM_NO_FACTOR_FOUND if no factor was found,
		 ECM_ERROR in case of error.
   (x, y) contains the new point at the end of Stage 1.
*/
int
ecm (mpz_t f, mpz_t x, mpz_t y, int param, mpz_t sigma, mpz_t n, mpz_t go, 
     double *B1done, double B1, mpz_t B2min_parm, mpz_t B2_parm,
     unsigned long k, const int S, int verbose, int repr, int nobase2step2, 
     int use_ntt, int sigma_is_A, ell_curve_t zE,
     FILE *os, FILE* es, char *chkfilename, char
     *TreeFilename, double maxmem, double stage1time, gmp_randstate_t rng, int
     (*stop_asap)(void), mpz_t batch_s, double *batch_last_B1_used,
     ATTRIBUTE_UNUSED double gw_k, ATTRIBUTE_UNUSED unsigned long gw_b,
     ATTRIBUTE_UNUSED unsigned long gw_n, ATTRIBUTE_UNUSED signed long gw_c,
     ATTRIBUTE_UNUSED signed long gw_cl_flag)
{
  int youpi = ECM_NO_FACTOR_FOUND;
  long st;
  mpmod_t modulus;
  curve P;
  ell_curve_t E;
#ifdef HAVE_ADDLAWS
  ell_point_t PE;
#endif
  mpz_t B2min, B2; /* Local B2, B2min to avoid changing caller's values */
  unsigned long dF;
#ifdef HAVE_GWNUM  
  unsigned long kbnc_bitsize;
  unsigned long gw_knbc_threshold;
  unsigned long gw_spec_div_threshold;
  unsigned long gw_generic_threshold;
#endif  
  root_params_t root_params;

  /*  1: sigma contains A from Montgomery form By^2 = x^3 + Ax^2 + x
      0: sigma contains sigma
     -1: sigma contains A from Weierstrass form y^2 = x^3 + Ax + B,
         and y contains y */
  ASSERT((-1 <= sigma_is_A) && (sigma_is_A <= 1));
  ASSERT((GMP_NUMB_BITS == 32) || (GMP_NUMB_BITS == 64));

  set_verbose (verbose);
  ECM_STDOUT = (os == NULL) ? stdout : os;
  ECM_STDERR = (es == NULL) ? stdout : es;

#ifdef MPRESN_NO_ADJUSTMENT
  /* When no adjustment is made in mpresn_ functions, N should be smaller
     than B^n/16 */
  if (mpz_sizeinbase (n, 2) > mpz_size (n) * GMP_NUMB_BITS - 4)
    {
      outputf (OUTPUT_ERROR, "Error, N should be smaller than B^n/16\n");
      return ECM_ERROR;
    }
#endif

  /* if a batch mode is requested by the user, this implies ECM_MOD_MODMULN */
  if (repr == ECM_MOD_DEFAULT && IS_BATCH_MODE(param))
    repr = ECM_MOD_MODMULN;

  /* choose the arithmetic used before the parametrization, since for divisors
     of 2^n+/-1 the default choice param=1 might not be optimal */
  if (mpmod_init (modulus, n, repr) != 0)
    return ECM_ERROR;

  repr = modulus->repr;

  /* If the parametrization is not given, choose it. */
  if (param == ECM_PARAM_DEFAULT)
    param = get_default_param (sigma_is_A, *B1done, repr);
  /* when dealing with several input numbers, if we had already computed
     batch_s, but the new number uses the base-2 representation, then we
     are forced to use ECM_PARAM_SUYAMA, and we reset batch_s to 1 to avoid
     the error "-bsaves/-bloads makes sense in batch mode only" below */
  if (param == ECM_PARAM_SUYAMA)
    mpz_set_ui (batch_s, 1);

  /* In batch mode, 
        we force repr=MODMULN, 
        B1done should be either the default value or greater than B1 
        x should be either 0 (undetermined) or 2 */
  if (IS_BATCH_MODE(param))
    {
      if (repr != ECM_MOD_MODMULN)
        {
          outputf (OUTPUT_ERROR, "Error, with param %d, repr should be " 
                                 "ECM_MOD_MODMULN.\n", param);
          return ECM_ERROR;
        }

      if (!ECM_IS_DEFAULT_B1_DONE(*B1done) && *B1done < B1)
        {
          outputf (OUTPUT_ERROR, "Error, cannot resume with param %d, except " 
		                 "for doing only stage 2\n", param);
          return ECM_ERROR;
        }

      if (sigma_is_A >= 0 && mpz_sgn (x) != 0 && mpz_cmp_ui (x, 2) != 0)
        {
          if (ECM_IS_DEFAULT_B1_DONE(*B1done))
            {
              outputf (OUTPUT_ERROR, "Error, x0 should be equal to 2 with this "
                                     "parametrization\n");
              return ECM_ERROR;
            }
        }
    }

  /* check that if ECM_PARAM_BATCH_SQUARE is used, GMP_NUMB_BITS == 64 */
  if (param == ECM_PARAM_BATCH_SQUARE && GMP_NUMB_BITS == 32)
    {
      outputf (OUTPUT_ERROR, "Error, parametrization ECM_PARAM_BATCH_SQUARE "
                             "works only with GMP_NUMB_BITS=64\n");
      return ECM_ERROR;
    }

  /* check that B1 is not too large */
  if (B1 > (double) ECM_UINT_MAX)
    {
      outputf (OUTPUT_ERROR, "Error, maximal step 1 bound for ECM is %lu.\n", 
               ECM_UINT_MAX);
      return ECM_ERROR;
    }

  /* See what kind of number we have as that may influence optimal parameter 
     selection. Test for base 2 number. Note: this was already done by
     mpmod_init. */

  mpres_init (P.x, modulus);
  mpres_init (P.y, modulus);
  mpres_init (P.A, modulus);

  /* In resume case, set (P.x,P.y) from x. */
  int resume = mpz_cmp_ui (x, 0);
  if (resume)
  {
    /* Call get_curve_from_param0() before setting P.x otherwise it would
       reset P.x. */
    if (param == ECM_PARAM_SUYAMA)
      get_curve_from_param0 (f, P.A, P.x, sigma, modulus);
    else if (param == ECM_PARAM_BATCH_SQUARE)
      get_curve_from_param1 (P.A, P.x, sigma, modulus);
    else if (param == ECM_PARAM_BATCH_2)
      get_curve_from_param2 (f, P.A, P.x, sigma, modulus);
    else
      get_curve_from_param3 (P.A, P.x, sigma, modulus);
    mpres_set_z (P.x, x, modulus);
    mpres_set_ui (P.y, 1, modulus);
  }

#ifdef HAVE_ADDLAWS
  ell_curve_set_z (E, zE, modulus);
#else
  E->type = ECM_EC_TYPE_MONTGOMERY;
#endif

  mpz_init (B2);
  mpz_init (B2min);
  youpi = set_stage_2_params (B2, B2_parm, B2min, B2min_parm,
			      &root_params, B1, &k, S, use_ntt,
			      &dF, TreeFilename, maxmem, modulus);

  /* if the user gave B2, print that B2 on the Using B1=..., B2=... line */
  if(!ECM_IS_DEFAULT_B2(B2_parm))
    mpz_set (B2, B2_parm);

  if (youpi == ECM_ERROR)
      goto end_of_ecm;

  if (!resume && sigma_is_A == 0)
    {
      if (mpz_sgn (sigma) == 0)
        {
          youpi = get_curve_from_random_parameter (f, P.A, P.x, sigma, param,
                                                   modulus, rng);

          if (youpi == ECM_ERROR)
            {
              outputf (OUTPUT_ERROR, "Error, invalid parametrization.\n");
	            goto end_of_ecm;
            }
        }
      else /* Compute A and x0 from given sigma values */
        {
          if (param == ECM_PARAM_SUYAMA)
              youpi = get_curve_from_param0 (f, P.A, P.x, sigma, modulus);
          else if (param == ECM_PARAM_BATCH_SQUARE)
              youpi = get_curve_from_param1 (P.A, P.x, sigma, modulus);
          else if (param == ECM_PARAM_BATCH_2)
              youpi = get_curve_from_param2 (f, P.A, P.x, sigma, modulus);
          else if (param == ECM_PARAM_BATCH_32BITS_D)
              youpi = get_curve_from_param3 (P.A, P.x, sigma, modulus);
          else
            {
              outputf (OUTPUT_ERROR, "Error, invalid parametrization.\n");
              youpi = ECM_ERROR;
	            goto end_of_ecm;
            }
      
          if (youpi != ECM_NO_FACTOR_FOUND)
            {
              if (youpi == ECM_ERROR)
                  outputf (OUTPUT_ERROR, "Error, invalid value of sigma.\n");

	            goto end_of_ecm;
            }
        }
    }
  else if (sigma_is_A == 1)
    {
      /* sigma contains the A value */
      mpres_set_z (P.A, sigma, modulus);
      /* TODO: make a valid, random starting point in case none was given */
      /* Problem: this may be as hard as factoring as we'd need to determine
         whether x^3 + a*x^2 + x is a quadratic residue or not */
      /* For now, we'll just chicken out. */
      /* Except for batch mode where we know that x0=2 */
      if (mpz_sgn (x) == 0)
        {
          if (IS_BATCH_MODE(param))
            mpres_set_ui (P.x, 2, modulus);
          else
            {
              outputf (OUTPUT_ERROR, 
                          "Error, -A requires a starting point (-x0 x).\n");
              youpi = ECM_ERROR;
	            goto end_of_ecm;
            }
        }
      else
          mpres_set_z (P.x, x, modulus);
    }

  /* Print B1, B2, polynomial and sigma */
  print_B1_B2_poly (OUTPUT_NORMAL, ECM_ECM, B1, *B1done, B2min_parm, B2min, 
		    B2, root_params.S, sigma, sigma_is_A, E->type,
		    y, param, 0);

#if 0
  outputf (OUTPUT_VERBOSE, "b2=%1.0f, dF=%lu, k=%lu, d=%lu, d2=%lu, i0=%Zd\n", 
           b2, dF, k, root_params.d1, root_params.d2, root_params.i0);
#else
  outputf (OUTPUT_VERBOSE, "dF=%lu, k=%lu, d=%lu, d2=%lu, i0=%Zd\n", 
           dF, k, root_params.d1, root_params.d2, root_params.i0);
#endif

  if (sigma_is_A == -1) /* Weierstrass or Hessian form. 
			   It is known that all curves in Weierstrass form do
			   not admit a Montgomery form. However, we could
			   be interested in performing some plain Step 1
			   on some special curves.
			*/
    {
	if (param != ECM_PARAM_TORSION)
	  {
	      mpres_set_z (P.A, sigma, modulus); /* sigma contains A */
	      mpres_set_z (P.x, x,    modulus);
	      mpres_set_z (P.y, y,    modulus);
	  }
	else
	  {
#ifdef HAVE_ADDLAWS
	      if(E->type == ECM_EC_TYPE_WEIERSTRASS)
#endif
              mpres_set_z(P.A, zE->a4, modulus);
#ifdef HAVE_ADDLAWS
	      else if(E->type == ECM_EC_TYPE_MONTGOMERY)
		  mpres_set_z(P.A, zE->a2, modulus);
#endif
              mpres_set_z (P.x, x, modulus);
              mpres_set_z (P.y, y, modulus);
	  }
    }

  if (test_verbose (OUTPUT_RESVERBOSE))
    {
      mpz_t t;

      mpz_init (t);
      mpres_get_z (t, P.A, modulus);
      outputf (OUTPUT_RESVERBOSE, "A=%Zd\n", t);
      mpres_get_z (t, P.x, modulus);
      outputf (OUTPUT_RESVERBOSE, "starting point: x0=%Zd\n", t);
#ifdef HAVE_ADDLAWS
      if (E->type == ECM_EC_TYPE_WEIERSTRASS)
	{
          mpres_get_z (t, P.y, modulus);
	  outputf (OUTPUT_RESVERBOSE, " y0=%Zd\n", t);
	}
#endif
      mpz_clear (t);
    }

  if (go != NULL && mpz_cmp_ui (go, 1) > 0)
    outputf (OUTPUT_VERBOSE, "initial group order: %Zd\n", go);

  if (test_verbose (OUTPUT_VERBOSE))
    {
      if (mpz_cmp_d (B2min, B1) != 0)
        {
          outputf (OUTPUT_VERBOSE, 
            "Can't compute success probabilities for B1 <> B2min\n");
        }
      else if (param == ECM_PARAM_DEFAULT)
        {
          outputf (OUTPUT_VERBOSE, "Can't compute success probabilities " 
                                   "for this parametrization.\n");
        }
      else
        {
          rhoinit (256, 10);
          print_expcurves (B1, B2, dF, k, root_params.S, param);
        }
    }

  /* Compute s for the batch mode */
  if (IS_BATCH_MODE(param) && ECM_IS_DEFAULT_B1_DONE(*B1done) &&
      (B1 != *batch_last_B1_used || mpz_cmp_ui (batch_s, 1) <= 0))
    {
      *batch_last_B1_used = B1;

      st = cputime ();
      /* construct the batch exponent */
      compute_s (batch_s, B1, NULL);
      outputf (OUTPUT_VERBOSE, "Computing batch product (of %" PRIu64
                               " bits) of primes up to B1=%1.0f took %ldms\n",
                               mpz_sizeinbase (batch_s, 2), B1,
                               elltime (st, cputime ()));
    }

  st = cputime ();

#ifdef HAVE_GWNUM
  /* gwnum only used when param == 0 and command line did not include -force-no-gwnum */
  if ((param == ECM_PARAM_SUYAMA) && (gw_cl_flag >= 0))
  {
    /* set thresholds */
    if (gw_cl_flag > 0) /* -force-gwnum specified in command line */
    {
      /* set all thresholds to minimum */
      gw_knbc_threshold = 350;
      gw_spec_div_threshold = 350;
      gw_generic_threshold = 350;
    }
    else /* use default thresholds */
    {
      gw_knbc_threshold = GWNUM_KBNC_THRESHOLD;
      gw_spec_div_threshold = GWNUM_SPEC_DIV_THRESHOLD;
      gw_generic_threshold = GWNUM_GENERIC_THRESHOLD;
    }

    /* We will use GWNUM for numbers N which have a kbnc form (b != 0)
       when N > 2^gw_knbc_threshold */
    /* we will check if N = (2^m +/- 1)/divisors (special division). In that case
       gwnum will be used only if m >= gw_spec_div_threshold */
    /* We will use GWNUM for numbers N of any other form (b == 0) when N > 2^gw_generic_threshold */

    /* check for special division case */
    if (modulus->repr == ECM_MOD_BASE2)
    {
      if ( abs (modulus->bits) >= (unsigned long) gw_spec_div_threshold )
      {
        /* use gwnum; kbnc params may not be set */
        gw_k = 1.0;
        gw_b = 2;
        gw_n = abs (modulus->bits);
        if ( modulus->bits < 0 ) gw_c = -1;
        else gw_c = 1;
      }
      else 
      {
        /* don't use gwnum */
        gw_b = 0;
        gw_k = 0.0;
      }
    }
    else if ( gw_b != 0 )
    {
      kbnc_bitsize = (unsigned long) ( (log(gw_k) + (double)gw_n * log((double)gw_b) )/log(2.0) );
      if (kbnc_bitsize < (unsigned long) gw_knbc_threshold)
      {
        /* don't use gwnum */
        gw_b = 0;
        gw_k = 0.0;
      }
    }
    else /* gw_b == 0 */
    {
      /* gwnum generic is not faster until input number is much larger -
         roughly 600 decimal digits (TBD) */
      if (mpz_sizeinbase(modulus->orig_modulus, 2) >= gw_generic_threshold)
      {
        gw_k = 1.0; /* indicate that input is large enough for gwnum generic */
        if (verbose > OUTPUT_NORMAL)
          printf ("Did not find a gwnum poly; will use gwnum generic.\n");
      }
      else
      {
        gw_k = 0.0; /* indicate that input is too small for gwnum generic */
        if (verbose > OUTPUT_NORMAL)
          printf ("Did not find a gwnum poly; input number is too small for gwnum generic.\n");
      }
    }

    if ((gw_b != 0 || (gw_b == 0 && gw_k >= 1.0) ) && B1 >= *B1done)
      youpi = gw_ecm_stage1 (f, &P, modulus, B1, B1done, go, gw_k, gw_b, gw_n, gw_c);
  }
  /* end if ((param == ECM_PARAM_SUYAMA) && (gw_cl_flag >= 0)) */

  /* At this point B1 == *B1done unless interrupted, or no GWNUM ecm_stage1
     is available */

  if (youpi != ECM_NO_FACTOR_FOUND)
    {
      outputf (OUTPUT_NORMAL, "Step 1 took %ldms\n", elltime (st, cputime ()));
      goto end_of_ecm_rhotable;
    }
#endif /* HAVE_GWNUM */

  if (B1 > *B1done || mpz_cmp_ui (go, 1) > 0)
    {
        if (IS_BATCH_MODE(param))
        /* FIXME: go, stop_asap and chkfilename are ignored in batch mode */
	    youpi = ecm_stage1_batch (f, P.x, P.A, modulus, B1, B1done, 
				      param, batch_s);
        else{
#ifdef HAVE_ADDLAWS
	    if(E->type == ECM_EC_TYPE_MONTGOMERY)
#endif
            youpi = ecm_stage1 (f, P.x, P.A, modulus, B1, B1done, go, 
                                stop_asap, chkfilename);
#ifdef HAVE_ADDLAWS
	    else{
		ell_point_init (PE, E, modulus);
		mpres_set (PE->x, P.x, modulus);
		mpres_set (PE->y, P.y, modulus);
		youpi = ecm_stage1_W (f, E, PE, modulus, B1, B1done, batch_s,
				      go, stop_asap, chkfilename);
		mpres_set (P.x, PE->x, modulus);
		mpres_set (P.y, PE->y, modulus);
                ell_point_clear (PE, E, modulus);
	    }
#endif
	}
    }
  else if (mpz_sgn (x) != 0)
    {
	/* when x <> 0, we initialize P to (x:y) */
	mpres_set_z (P.x, x, modulus);
	mpres_set_z (P.y, y, modulus);
    }

  if (stage1time > 0.)
    {
      const long st2 = elltime (st, cputime ());
      const long s1t = (long) (stage1time * 1000.);
      outputf (OUTPUT_NORMAL, 
               "Step 1 took %ldms (%ld in this run, %ld from previous runs)\n", 
               st2 + s1t, st2, s1t);
    }
  else
    outputf (OUTPUT_NORMAL, "Step 1 took %ldms\n", elltime (st, cputime ()));

  /* Store end-of-stage-1 residue in x in case we write it to a save file, 
     before P.x is (perhaps) converted to Weierstrass form */
  
  mpres_get_z (x, P.x, modulus);
#ifdef HAVE_ADDLAWS
  if (E->type == ECM_EC_TYPE_WEIERSTRASS 
      || E->type == ECM_EC_TYPE_HESSIAN 
      || E->type == ECM_EC_TYPE_TWISTED_HESSIAN)
    mpres_get_z (y, P.y, modulus);  
#endif

  if (youpi != ECM_NO_FACTOR_FOUND)
    goto end_of_ecm_rhotable;

  if (test_verbose (OUTPUT_RESVERBOSE)) 
    {
      mpz_t t;
      
      mpz_init (t);
      mpres_get_z (t, P.x, modulus);
      outputf (OUTPUT_RESVERBOSE, "x=%Zd\n", t);
#ifdef HAVE_ADDLAWS
      if (E->type == ECM_EC_TYPE_WEIERSTRASS)
	{
	  mpres_get_z (t, P.y, modulus);
	  outputf (OUTPUT_RESVERBOSE, "y=%Zd\n", t);
	}
#endif
      mpz_clear (t);
    }

  /* In case of a signal, we'll exit after the residue is printed. If no save
     file is specified, the user may still resume from the residue */
  if (stop_asap != NULL && (*stop_asap) ())
    goto end_of_ecm_rhotable;

  /* If using 2^k +/-1 modulus and 'nobase2step2' flag is set,
     set default (-nobase2) modular method and remap P.x, P.y, and P.A */
  if (modulus->repr == ECM_MOD_BASE2 && nobase2step2)
    {
      mpz_t x_t, y_t, A_t;

      mpz_init (x_t);
      mpz_init (y_t);
      mpz_init (A_t);

      mpz_mod (x_t, P.x, modulus->orig_modulus);
      mpz_mod (y_t, P.y, modulus->orig_modulus);
      mpz_mod (A_t, P.A, modulus->orig_modulus);

      mpmod_clear (modulus);

      repr = ECM_MOD_NOBASE2;
      if (mpmod_init (modulus, n, repr) != 0) /* reset modulus for nobase2 */
        return ECM_ERROR;

      /* remap x, y, and A for new modular method */
      mpres_set_z (P.x, x_t, modulus);
      mpres_set_z (P.y, y_t, modulus);
      mpres_set_z (P.A, A_t, modulus);

      mpz_clear (x_t);
      mpz_clear (y_t);
      mpz_clear (A_t);
    }

#ifdef HAVE_ADDLAWS
  if (E->type == ECM_EC_TYPE_MONTGOMERY)
#endif
  youpi = montgomery_to_weierstrass (f, P.x, P.y, P.A, modulus);
#ifdef HAVE_ADDLAWS
  else if (E->type == ECM_EC_TYPE_HESSIAN)
    {
      youpi = hessian_to_weierstrass (f, P.x, P.y, P.A, modulus);
      if(youpi == ECM_NO_FACTOR_FOUND)
        /* due to that non-trivial kernel(?) */
        youpi = mult_by_3(f, P.x, P.y, P.A, modulus);
    }
  else if (E->type == ECM_EC_TYPE_TWISTED_HESSIAN)
    {
	mpz_t c, rm;
	mpz_init(c);
	mpz_init(rm);
	mpres_get_z(rm, E->a4, modulus);
	mpz_rootrem(c, rm, rm, 3);
	if(mpz_sgn(rm) != 0){
	    printf("ECM_EC_TYPE_TWISTED_HESSIAN: not a cube!\n");
	    exit(-1);
	}
	mpres_set_z(P.A, c, modulus);
	mpz_clear(c);
	mpz_clear(rm);
	youpi = twisted_hessian_to_weierstrass (f, P.x, P.y, P.A, E->a6, modulus);
	if(youpi == ECM_NO_FACTOR_FOUND){
	    /* due to that non-trivial kernel(?) */
	    youpi = mult_by_3(f, P.x, P.y, P.A, modulus);
	}
    }
#endif
  
  if (test_verbose (OUTPUT_RESVERBOSE) && youpi == ECM_NO_FACTOR_FOUND && 
      mpz_cmp (B2, B2min) >= 0)
    {
      mpz_t t;

      mpz_init (t);
      mpres_get_z (t, P.x, modulus);
      outputf (OUTPUT_RESVERBOSE, "After switch to Weierstrass form, "
      "P=(%Zd", t);
      mpres_get_z (t, P.y, modulus);
      outputf (OUTPUT_RESVERBOSE, ", %Zd)\n", t);
      mpres_get_z (t, P.A, modulus);
      outputf (OUTPUT_RESVERBOSE, "on curve Y^2 = X^3 + %Zd * X + b\n", t);
      mpz_clear (t);
    }

  P.disc = 0; /* FIXME: should disappear one day */
  
  if (youpi == ECM_NO_FACTOR_FOUND && mpz_cmp (B2, B2min) >= 0)
    youpi = stage2 (f, &P, modulus, dF, k, &root_params, use_ntt, 
                    TreeFilename, 0, stop_asap);
#ifdef TIMING_CRT
  printf ("mpzspv_from_mpzv_slow: %dms\n", mpzspv_from_mpzv_slow_time);
  printf ("mpzspv_to_mpzv: %dms\n", mpzspv_to_mpzv_time);
  printf ("mpzspv_normalise: %dms\n", mpzspv_normalise_time);
#endif
  
end_of_ecm_rhotable:
  if (test_verbose (OUTPUT_VERBOSE))
    {
      if (mpz_cmp_d (B2min, B1) == 0 && param != ECM_PARAM_DEFAULT)
        {
          if (youpi == ECM_NO_FACTOR_FOUND && 
              (stop_asap == NULL || !(*stop_asap)()))
            print_exptime (B1, B2, dF, k, root_params.S, 
                           (long) (stage1time * 1000.) + 
                           elltime (st, cputime ()), param);
          rhoinit (1, 0); /* Free memory of rhotable */
        }
    }

end_of_ecm:
#ifdef HAVE_ADDLAWS
  ell_curve_clear(E, modulus);
#endif
  mpres_clear (P.A, modulus);
  mpres_clear (P.y, modulus);
  mpres_clear (P.x, modulus);
  mpz_clear (root_params.i0);
  mpz_clear (B2);
  mpz_clear (B2min);
  mpmod_clear (modulus);
  return youpi;
}
