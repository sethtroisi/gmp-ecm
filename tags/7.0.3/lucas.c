/* Auxiliary functions to evaluate Lucas sequences.

Copyright 2002, 2003, 2005, 2006, 2008, 2011, 2012, 2015
Paul Zimmermann, Alexander Kruppa, Dave Newman.

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

/* References:

A p+1 Method of Factoring, H. C. Williams, Mathematics of Computation,
volume 39, number 159, pages 225-234, 1982.

Evaluating recurrences of form X_{m+n} = f(X_m, X_n, X_{m-n}) via
Lucas chains, Peter L. Montgomery, December 1983, revised January 1992. */

#include "ecm-impl.h"

/* P <- V_2(Q) */
static void
pp1_duplicate (mpres_t P, mpres_t Q, mpmod_t n)
{
  mpres_sqr (P, Q, n);
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

/* computes V_k(P) from P=A and puts the result in P=A. Assumes k>2.
   Uses auxiliary variables t, B, C, T, T2.
*/
void
pp1_mul_prac (mpres_t A, ecm_uint k, mpmod_t n, mpres_t t, mpres_t B,
              mpres_t C, mpres_t T, mpres_t T2)
{
  ecm_uint d, e, r;
  static double val = 0.61803398874989485; /* 1/(golden ratio) */

  /* Note: we used to use several (4) values of "val", but:
     (1) the code to estimate the best value was buggy;
     (2) even after fixing the bug, the overhead to choose the
         best value was larger than the corresponding gain (for a c155
         and B1=10^7). */

  d = k;
  r = (ecm_uint) ((double) d * val + 0.5);
  
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
          mpres_swap (B, C, n);
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
          mpres_swap (B, C, n);
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
          mpres_swap (B, T, n);          /* swap B and T */
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
