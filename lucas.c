/* Auxiliary functions to evaluate Lucas sequences.

  Copyright (C) 2002 Alexander Kruppa and Paul Zimmermann.

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

  References:

  A p+1 Method of Factoring, H. C. Williams, Mathematics of Computation,
  volume 39, number 159, pages 225-234, 1982.

  Evaluating recurrences of form X_{m+n} = f(X_m, X_n, X_{m-n}) via
  Lucas chains, Peter L. Montgomery, December 1983, revised January 1992.
*/

#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"

#define ADD 1
#define DUP 1

#define duplicate pp1_duplicate
#define add3      pp1_add3

/* P <- V_2(Q) */
inline void
duplicate (mpz_t P, mpz_t Q, mpz_t n)
{
  mpz_mul (P, Q, Q);
  mpz_sub_ui (P, P, 2);
  mpz_mod (P, P, n);
}

/* P <- V_{m+n} where Q = V_m, R = V_n, S = V_{m-n}.
   t is an auxiliary variable.
   Warning: P may equal Q, R or S.
*/
inline void
add3 (mpz_t P, mpz_t Q, mpz_t R, mpz_t S, mpz_t n, mpz_t t)
{
  mpz_mul (t, Q, R);
  mpz_sub (t, t, S);
  mpz_mod (P, t, n);
}

/* returns the number of modular multiplications for computing
   V_n from V_r * V_{n-r} - V_{n-2r}.
   ADD is the cost of an addition
   DUP is the cost of a duplicate
*/
unsigned int
lucas_cost_pp1 (unsigned n, double v)
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
          fprintf (stderr, "lucas_cost: no condition qualifies for d=%u e=%u\n", d, e);
          exit (1);
        }
    }
  
  return (c);
}


#define NV 4

// #define SWAP(x,y) { __mpz_struct *tmp = x; x = y; y = tmp; }
#define SWAP mpz_swap

/* computes V_k(P) from P=A and puts the result in P=A. Assumes k>2.
   Uses auxiliary variables t, B, C, T, T2.
*/
void
pp1_mul_prac (mpz_t A, unsigned long k, mpz_t n, mpz_t t, mpz_t B, mpz_t C,
          mpz_t T, mpz_t T2)
{
  unsigned int d, e, r, i = 0;
  static double val[NV] =
  { 1.61803398875, 1.72360679775, 1.618347119656, 1.617914406529};
  
  /* chooses the best value of v */
  for (d = 0, r = ADD * k; d < NV; d++)
    {
      e = lucas_cost_pp1 (k, val[d]);
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
  mpz_set (B, A); /* B=A */
  mpz_set (C, A); /* C=A */
  duplicate (A, A, n); /* A = 2*A */
  while (d != e)
    {
      if (d < e)
        {
          r = d;
          d = e;
          e = r;
          mpz_swap (A, B);
        }
      /* do the first line of Table 4 whose condition qualifies */
      if (4 * d <= 5 * e && ((d + e) % 3) == 0)
        { /* condition 1 */
          d = (2 * d - e) / 3;
          e = (e - d) / 2;
          add3 (T,  A, B, C, n, t); /* T = f(A,B,C) */
          add3 (T2, T, A, B, n, t); /* T2 = f(T,A,B) */
          add3 (B,  B, T, A, n, t); /* B = f(B,T,A) */
          mpz_swap (A, T2);         /* swap A and T2 */
        }
      else if (4 * d <= 5 * e && (d - e) % 6 == 0)
        { /* condition 2 */
          d = (d - e) / 2;
          add3 (B, A, B, C, n, t); /* B = f(A,B,C) */
          duplicate (A, A, n);     /* A = 2*A */
        }
      else if (d <= (4 * e))
        { /* condition 3 */
          d -= e;
          add3 (C, B, A, C, n, t); /* C = f(B,A,C) */
          SWAP (B, C);
        }
      else if ((d + e) % 2 == 0)
        { /* condition 4 */
          d = (d - e) / 2;
          add3 (B, B, A, C, n, t); /* B = f(B,A,C) */
          duplicate (A, A, n);     /* A = 2*A */
        }
      else if (d % 2 == 0)
        { /* condition 5 */
          d /= 2;
          add3 (C, C, A, B, n, t); /* C = f(C,A,B) */
          duplicate (A, A, n);     /* A = 2*A */
        }
      else if (d % 3 == 0)
        { /* condition 6 */
          d = d / 3 - e;
          duplicate (T, A, n);       /* T = 2*A */
          add3 (T2, A, B, C, n, t);  /* T2 = f(A,B,C) */
          add3 (A,  T, A, A, n, t);  /* A = f(T,A,A) */
          add3 (C,  T, T2, C, n, t); /* C = f(T,T2,C) */
          SWAP (B, C);
        }
      else if ((d + e) % 3 == 0)
        { /* condition 7 */
          d = (d - 2 * e) / 3;
          add3 (T, A, B, C, n, t); /* T1 = f(A,B,C) */
          add3 (B, T, A, B, n, t); /* B = f(T1,A,B) */
          duplicate (T, A, n);
          add3 (A, A, T, A, n, t); /* A = 3*A */
        }
      else if ((d - e) % 3 == 0)
        { /* condition 8: never happens? */
          d = (d - e) / 3;
          add3 (T, A, B, C, n, t); /* T1 = f(A,B,C) */
          add3 (C, C, A, B, n, t); /* C = f(A,C,B) */
          SWAP (B, T);         /* swap B and T */
          duplicate (T, A, n);
          add3 (A, A, T, A, n, t); /* A = 3*A */
        }
      else if (e % 2 == 0)
        { /* condition 9: never happens? */
          e /= 2;
          add3 (C, C, B, A, n, t); /* C = f(C,B,A) */
          duplicate (B, B, n);     /* B = 2*B */
        }
      else
        {
          fprintf (stderr, "no condition qualifies for d=%u e=%u\n", d, e);
          exit (1);
        }
    }
  
  add3 (A, A, B, C, n, t);
  
  if (d != 1)
    {
      fprintf (stderr, "d!=1 at the end of PRAC\n");
      exit (1);
    }
}
