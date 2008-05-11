/* Auxiliary functions to evaluate Lucas sequences.

  Copyright 2002, 2003, 2005 Paul Zimmermann and Alexander Kruppa.

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

  References:

  A p+1 Method of Factoring, H. C. Williams, Mathematics of Computation,
  volume 39, number 159, pages 225-234, 1982.

  Evaluating recurrences of form X_{m+n} = f(X_m, X_n, X_{m-n}) via
  Lucas chains, Peter L. Montgomery, December 1983, revised January 1992.
*/

#include "ecm-impl.h"

/* the following constant were tested experimentally: in theory, we should
   have ADD = 2*DUP in the basecase range, and ADD = 3/2*DUP in the FFT range.
   Only the ratio ADD/DUP counts. */
#define ADD 3 /* cost of add3:      one modular multiply */
#define DUP 2 /* cost of duplicate: one square */

/* P <- V_2(Q) */
static void
pp1_duplicate (mpres_t P, mpres_t Q, mpmod_t n)
{
  mpres_mul (P, Q, Q, n);
  mpres_sub_ui (P, P, 2, n);
}

/* P <- V_{m+n} where Q = V_m, R = V_n, S = V_{m-n}.
   t is an auxiliary variable.
   Warning: P may equal Q, R or S.
*/
static void
pp1_add3 (mpres_t P, mpres_t Q, mpres_t R, mpres_t S, mpmod_t n, mpres_t t)
{
  mpres_mul (t, Q, R, n);
  mpres_sub (P, t, S, n);
}

/* returns the number of modular multiplications for computing
   V_n from V_r * V_{n-r} - V_{n-2r}.
   ADD is the cost of an addition
   DUP is the cost of a duplicate
*/
static unsigned int
lucas_cost_pp1 (unsigned long n, double v)
{
  unsigned int c;
  unsigned long d, e, r;

  d = n;
  r = (unsigned long) ((double) d / v + 0.5);
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
      if (d - e <= e / 4 && ((d + e) % 3) == 0)
        { /* condition 1 */
          d = (2 * d - e) / 3;
          e = (e - d) / 2;
          c += 3 * ADD; /* 3 additions */
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
      /* now d is odd and e even */
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
      else /* necessarily e is even */
        { /* condition 9 */
          e /= 2;
          c += ADD + DUP; /* one addition, one duplicate */
        }
    }
  
  return c;
}


#define NV 4

/* #define SWAP(x,y) { __mpz_struct *tmp = x; x = y; y = tmp; } */
#define SWAP mpres_swap

/* computes V_k(P) from P=A and puts the result in P=A. Assumes k>2.
   Uses auxiliary variables t, B, C, T, T2.
*/
void
pp1_mul_prac (mpres_t A, unsigned long k, mpmod_t n, mpres_t t, mpres_t B,
              mpres_t C, mpres_t T, mpres_t T2)
{
  unsigned long d, e, r, i = 0;
  static double val[NV] =
    {0.61803398874989485, 0.5801787282954641, 0.6179144065288179 , 0.6180796684698958};
    /* 1/GR,              5/(GR+7) (2),       1429/(GR+2311) (8),  3739/(6051-GR) (9) */

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
  r = (unsigned long) ((double) d * val[i] + 0.5);
  
  /* first iteration always begins by Condition 3, then a swap */
  d = k - r;
  e = 2 * r - k;
  mpres_set (B, A, n); /* B=A */
  mpres_set (C, A, n); /* C=A */
  pp1_duplicate (A, A, n); /* A = 2*A */
  while (d != e)
    {
      if (d < e)
        {
          r = d;
          d = e;
          e = r;
          mpres_swap (A, B, n);
        }
      /* do the first line of Table 4 whose condition qualifies */
      if (d - e <= e / 4 && ((d + e) % 3) == 0)
        { /* condition 1 */
          d = (2 * d - e) / 3;
          e = (e - d) / 2;
          pp1_add3 (T,  A, B, C, n, t); /* T = f(A,B,C) */
          pp1_add3 (T2, T, A, B, n, t); /* T2 = f(T,A,B) */
          pp1_add3 (B,  B, T, A, n, t); /* B = f(B,T,A) */
          mpres_swap (A, T2, n);    /* swap A and T2 */
        }
      else if (d - e <= e / 4 && (d - e) % 6 == 0)
        { /* condition 2 */
          d = (d - e) / 2;
          pp1_add3 (B, A, B, C, n, t); /* B = f(A,B,C) */
          pp1_duplicate (A, A, n);     /* A = 2*A */
        }
      else if ((d + 3) / 4 <= e) /* <==>  (d <= 4 * e) */
        { /* condition 3 */
          d -= e;
          pp1_add3 (C, B, A, C, n, t); /* C = f(B,A,C) */
          SWAP (B, C, n);
        }
      else if ((d + e) % 2 == 0)
        { /* condition 4 */
          d = (d - e) / 2;
          pp1_add3 (B, B, A, C, n, t); /* B = f(B,A,C) */
          pp1_duplicate (A, A, n);     /* A = 2*A */
        }
      /* d+e is now odd */
      else if (d % 2 == 0)
        { /* condition 5 */
          d /= 2;
          pp1_add3 (C, C, A, B, n, t); /* C = f(C,A,B) */
          pp1_duplicate (A, A, n);     /* A = 2*A */
        }
      /* d is odd, e even */
      else if (d % 3 == 0)
        { /* condition 6 */
          d = d / 3 - e;
          pp1_duplicate (T, A, n);       /* T = 2*A */
          pp1_add3 (T2, A, B, C, n, t);  /* T2 = f(A,B,C) */
          pp1_add3 (A,  T, A, A, n, t);  /* A = f(T,A,A) */
          pp1_add3 (C,  T, T2, C, n, t); /* C = f(T,T2,C) */
          SWAP (B, C, n);
        }
      else if ((d + e) % 3 == 0) /* d+e <= val[i]*k < k < 2^32 */
        { /* condition 7 */
          d = (d - 2 * e) / 3;
          pp1_add3 (T, A, B, C, n, t); /* T1 = f(A,B,C) */
          pp1_add3 (B, T, A, B, n, t); /* B = f(T1,A,B) */
          pp1_duplicate (T, A, n);
          pp1_add3 (A, A, T, A, n, t); /* A = 3*A */
        }
      else if ((d - e) % 3 == 0)
        { /* condition 8: never happens? */
          d = (d - e) / 3;
          pp1_add3 (T, A, B, C, n, t); /* T1 = f(A,B,C) */
          pp1_add3 (C, C, A, B, n, t); /* C = f(A,C,B) */
          SWAP (B, T, n);          /* swap B and T */
          pp1_duplicate (T, A, n);
          pp1_add3 (A, A, T, A, n, t); /* A = 3*A */
        }
      else /* necessarily e is even */
        { /* condition 9: never happens? */
          e /= 2;
          pp1_add3 (C, C, B, A, n, t); /* C = f(C,B,A) */
          pp1_duplicate (B, B, n);     /* B = 2*B */
        }
    }
  
  pp1_add3 (A, A, B, C, n, t);

  ASSERT(d == 1);
}
