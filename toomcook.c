/* Implementation of the Toom-Cook 3-way and 4-way algorithms for 
   polynomial convolution products. 
  
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

#include <stdio.h> /* to avoid #undef all prototypes with FILE in ecm-impl.h */
#include "gmp.h"
#include "ecm-impl.h"

#define A0 A[i]
#define A1 A[l+i]
#define A2 A[2*l+i]
#define B0 B[i]
#define B1 B[l+i]
#define B2 B[2*l+i]
#define C0 C[i]
#define C1 C[l+i]
#define C2 C[2*l+i]
#define C3 C[3*l+i]
#define C4 C[4*l+i]
#define t0 t[i]
#define t2 t[2*l+i-1]
#define T  t[4*l-2]

/* Puts in C[0..2len-2] the product of A[0..len-1] and B[0..len-1].
   This version works for all input sizes, but cannot handle input arrays
   overlapping with output.
   Assumes len >= 1.

   The auxiliary memory M(len) necessary in t satisfies:
   M(0) = 0, M(1) = 0, M(2) = 1, M(3) = 3,
   otherwise M(len) = 2*(2*l-1) + max(M(l), 1)
   with l = ceil(len/3).

   We prove M(len) <= 2*len + 2 * k with k = ceil(log[3](len)) by induction:
        4*l-2 + max(M(l), 1) <= 4*l-2 + max(2*l + 2 * (k-1), 1)
                             <= 6*l - 2 + 2 * (k-1)
                             <= 2*(len+2) - 2 + 2 * (k-1)
                             <= 2*len + 2 * k
*/

void
toomcook3 (listz_t C, listz_t A, listz_t B, unsigned int len, listz_t t)
{
  int i, l, k;

  if (len <= 2 || len == 4)
    {
      karatsuba (C, A, B, len, t);
      return;
    }

  l = (len + 2) / 3;            /* ceil(len/3) */
  k = len - 2 * l;              /* smaller part */

  for (i = 0; i < k; i++)       /* uses t[0..3*l+k-1] */
    {
      mpz_add (C0, A0, A2);
      mpz_sub (C2, C0, A1);     /* C2 = A0 - A1 + A2 = A(-1) */
      mpz_add (C0, C0, A1);     /* C0 = A0 + A1 + A2 = A(1) */
      mpz_add (C1, B0, B2);
      mpz_sub (C3, C1, B1);     /* C3 = B0 - B1 + B2 = B(-1) */
      mpz_add (C1, C1, B1);     /* C1 = B0 + B1 + B2 = B(1) */
    }
  for (; i < l; i++)            /* uses t[0..4*l-1] */
    {                           /* A2 and B2 are smaller than the rest */
      mpz_add (C0, A0, A1);
      mpz_sub (C2, A0, A1);
      mpz_add (C1, B0, B1);
      mpz_sub (C3, B0, B1);
    }

  toomcook3 (t, C + 2 * l, C + 3 * l, l, &T);
  /* t0 = C2*C3 = A(-1)*B(-1) = C(-1), len(t0) = 2*l-1 */

  for (i = 0; i < k; i++)
    {
      mpz_mul_2exp (C2, A2, 1); /* C2 = A(2), C3 = B(2) */
      mpz_add (C2, C2, A1);
      mpz_mul_2exp (C2, C2, 1);
      mpz_add (C2, C2, A0);
      mpz_mul_2exp (C3, B2, 1);
      mpz_add (C3, C3, B1);
      mpz_mul_2exp (C3, C3, 1);
      mpz_add (C3, C3, B0);
    }
  for (; i < l; i++)
    {
      mpz_mul_2exp (C2, A1, 1);
      mpz_add (C2, C2, A0);
      mpz_mul_2exp (C3, B1, 1);
      mpz_add (C3, C3, B0);
    }

  toomcook3 (t + 2 * l - 1, C + 2 * l, C + 3 * l, l, &T);
  /* t2 = C2*C3 = A(2)*B(2) = C(2), len(t2) = 2*l-1 */

  toomcook3 (C + 2 * l, C, C + l, l, &T);
  /* C2 = C0*C1 = A(1)*B(1) = C(1), len(C1) = 2*l-1 */

  toomcook3 (C, A, B, l, &T);
  /* C0 = A(0)*B(0) = C(0), len(C0) = 2*l-1 */

  toomcook3 (C + 4 * l, A + 2 * l, B + 2 * l, k, &T);
  /* C4 = A(inf)*B(inf) = C(inf), len(C4) = 2*k-1 */

  /* C0: C_0  C2: C(1)  C4: C_4  t0: C(-1)  t2: C(2) */

  /* C_3 = A_1 * B_2 + A_2 * B_1, len(C_3) = l+k-1
     We need not bother to compute C_3[2l-1] if k<l, as it will be zero */

  for (i = 0; i < 2 * l - 1; i++)       /* uses t[0..4l] */
    {
      mpz_sub (T, C2, t0);
      mpz_add (C2, C2, t0);             /* C2 = C(1)+C(-1) = 2*(C_0 + C_2 + C_4) */
      mpz_set (t0, T);                  /* t0 = C(1)-C(-1) = 2*(C_1 + C_3) */
      mpz_fdiv_q_2exp (C2, C2, 1);      /* C2 = C_0 + C_2 + C_4 */
      mpz_sub (C2, C2, C0);             /* C2 = C_2 + C_4 */
      if (i < l + k - 1)
        {
          mpz_sub (t2, t2, C0);         /* t2 = 2*C_1 + 4*C_2 + 8*C_3 + 16*C_4 */
          mpz_sub (t2, t2, t0);         /* t2 = 4*C_2 + 6*C_3 + 16*C_4 */
          mpz_fdiv_q_2exp (t2, t2, 1);  /* t2 = 2*C_2 + 3*C_3 + 8*C_4 */
          if (i < 2 * k - 1)
            {
              mpz_mul_2exp (T, C4, 3);  /* T = 8*C_4 */
              mpz_sub (t2, t2, T);      /* t2 = 2*C_2 + 3*C_3 */
              mpz_sub (C2, C2, C4);     /* C2 = C_2 */
            }
          mpz_mul_2exp (T, C2, 1);
          mpz_sub (t2, t2, T);          /* t2 = 3*C_3 */
          mpz_divby3_1op (t2);          /* t2 = C_3 */
          mpz_fdiv_q_2exp (t0, t0, 1);  /* t0 = C_1 + C_3 */
          mpz_sub (t0, t0, t2);         /* t0 = C_1 */
        }
      else                              /* C_3 == 0, t0 == 2*C_1 */
        {
          mpz_fdiv_q_2exp (t0, t0, 1);
        }
    }

  for (i = 0; i < l - 1; i++)
    mpz_add (C1, C1, t0);
  mpz_set (C1, t0);
  for (i = l; i < 2 * l - 1; i++)
    mpz_add (C1, C1, t0);

  for (i = 0; i < l - 1; i++)
    mpz_add (C3, C3, t2);
  mpz_set (C3, t2);
  for (i = l; i < l + k - 1; i++)
    mpz_add (C3, C3, t2);
}

#define A3 A[3*l+i]
#define B3 B[3*l+i]
#define C5 C[5*l+i]
#define C6 C[6*l+i]
#define t4 t[4*l+i-2]
#undef T
#define T t[6*l-3]

/* the number of multiplies satisfies

   T4(n) = 6 * T4(l) + T4(n-3*l) where l = ceil(n/4)
*/
void
toomcook4 (listz_t C, listz_t A, listz_t B, unsigned int len, listz_t t)
{
  unsigned int l, k, i;

  /* toomcook4 cannot handle len = 1, 2, 5
     since we need k := len - 3*ceil(len/4) >= 0.
     For len=3, toomcook4 would use 6 multiplies, toomcook3 uses only 5.
     For len=6, toomcook4 would use 18 multiplies, toomcook3 only 15.
     For len=9, toomcook4 would use 30 multiplies, toomcook3 only 25.
     Further values where toomcook3 is faster are 17,18,26,27,77,78,79,80,81.
   */

  if (len <= 2)
    {
      karatsuba (C, A, B, len, t);
      return;
    }

  if (len == 3 || len == 5 || len == 6 || len == 9 || len == 17 || len == 18 ||
      (25 <= len && len <= 27) || (77 <= len && len <= 81))
    {
      toomcook3 (C, A, B, len, t);
      return;
    }

  l = (len + 3) / 4;            /* l = ceil(len/4) */
  k = len - 3 * l;              /* k = smaller part. len = 3*l + k, k <= l */

  for (i = 0; i < l; i++)
    {
      /*** Evaluate A(2), A(-2), 8*A(1/2) ***/
      mpz_mul_2exp (C0, A0, 1);
      mpz_add (C0, C0, A1);
      mpz_mul_2exp (C0, C0, 1);
      mpz_add (C0, C0, A2);
      mpz_mul_2exp (C0, C0, 1);
      if (i < k)
        {
          mpz_add (C0, C0, A3);         /* C[0 .. l-1] = 8*A(1/2) */

          mpz_mul_2exp (C2, A3, 2);
          mpz_add (C2, C2, A1);
          mpz_mul_2exp (C2, C2, 1);     /* C[2l .. 3l-1] = 8*A_3 + 2*A_1 */
        }
      else
        mpz_mul_2exp (C2, A1, 1);
      mpz_mul_2exp (T, A2, 2);
      mpz_add (T, T, A0);       /* T = 4*A_2 + A0 */
      mpz_sub (C4, T, C2);      /* C[4l .. 5l-1] = A(-2) */
      mpz_add (C2, C2, T);      /* C[2l .. 3l-1] = A(2) */
#ifdef DEBUG
      gmp_printf ("8*A(1/2)[%d] = %Zd\n", i, C0);
      gmp_printf ("A(2)[%d] = %Zd\n", i, C2);
      gmp_printf ("A(-2)[%d] = %Zd\n", i, C4);
#endif

      /*** Evaluate B(2), B(-2), 8*B(1/2) ***/
      mpz_mul_2exp (C1, B0, 1);
      mpz_add (C1, C1, B1);
      mpz_mul_2exp (C1, C1, 1);
      mpz_add (C1, C1, B2);
      mpz_mul_2exp (C1, C1, 1);
      if (i < k)
        {
          mpz_add (C1, C1, B3); /* C[l .. 2l-1] = 8*B(1/2) */

          mpz_mul_2exp (C3, B3, 2);
          mpz_add (C3, C3, B1);
          mpz_mul_2exp (C3, C3, 1);     /* C[3l .. 3l+k-1] = 8*B_3 + 2*B_1 */
        }
      else
        mpz_mul_2exp (C3, B1, 1);
      mpz_mul_2exp (T, B2, 2);
      mpz_add (T, T, B0);       /* T = 4*B_2 + B0 */
      mpz_sub (C5, T, C3);      /* C[5l .. 5l+k-1] = B(-2) */
      mpz_add (C3, C3, T);      /* C[3l .. 3l+k-1] = B(2) */
#ifdef DEBUG
      gmp_printf ("8*B(1/2)[%d] = %Zd\n", i, C1);
      gmp_printf ("B(2)[%d] = %Zd\n", i, C3);
      gmp_printf ("B(-2)[%d] = %Zd\n", i, C5);
#endif
    }

  toomcook4 (t, C, C + l, l, &T);   /* t0 = 8*A(1/2) * 8*B(1/2) = 64*C(1/2) */
  toomcook4 (t + 2 * l - 1, C + 2 * l, C + 3 * l, l, &T);  /* t2 = A(2) * B(2) = C(2) */
  toomcook4 (t + 4 * l - 2, C + 4 * l, C + 5 * l, l, &T);  /* t4 = A(-2) * B(-2) = C(-2) */

  for (i = 0; i < l; i++)
    {
      mpz_add (C0, A0, A2);
      mpz_add (C1, B0, B2);
      if (i < k)
        {
          mpz_add (T, A1, A3);
          mpz_sub (C2, C0, T);  /* C2 = A(-1) */
          mpz_add (C0, C0, T);  /* C0 = A(1) */
          mpz_add (T, B1, B3);
          mpz_sub (C3, C1, T);  /* C3 = B(-1) */
          mpz_add (C1, C1, T);  /* C1 = B(1) */
        }
      else
        {
          mpz_sub (C2, C0, A1);
          mpz_add (C0, C0, A1);
          mpz_sub (C3, C1, B1);
          mpz_add (C1, C1, B1);
        }

#ifdef DEBUG
      gmp_printf ("A(1)[%d] = %Zd\n", i, C0);
      gmp_printf ("A(-1)[%d] = %Zd\n", i, C2);
      gmp_printf ("B(1)[%d] = %Zd\n", i, C1);
      gmp_printf ("B(-1)[%d] = %Zd\n", i, C3);
#endif
    }

  toomcook4 (C + 4 * l, C + 2 * l, C + 3 * l, l, &T);
  /* C4 = A(-1) * B(-1) = C(-1) */

  toomcook4 (C + 2 * l, C, C + l, l, &T);
  /* C2 = A(1) * B(1) = C(1) */

  toomcook4 (C, A, B, l, &T);
  /* C0 = A_0 * B_0 = C_0 */
  
  toomcook4 (C + 6 * l, A + 3 * l, B + 3 * l, k, &T);
  /* C6 = A_3 * B_3 = C_6 */

  for (i = 0; i < 2 * l - 1; i++)
    {
#ifdef DEBUG
      gmp_printf ("C(0)[%d] = %Zd\n", i, C0);
      gmp_printf ("C(1)[%d] = %Zd\n", i, C2);
      gmp_printf ("C(-1)[%d] = %Zd\n", i, C4);
      gmp_printf ("C(2)[%d] = %Zd\n", i, t2);
      gmp_printf ("C(-2)[%d] = %Zd\n", i, t4);
      gmp_printf ("64*C(1/2)[%d] = %Zd\n", i, t0);
      if (i < 2 * k - 1)
        gmp_printf ("C(inf)[%d] = %Zd\n", i, C6);

      gmp_printf ("C_0[%d] = %Zd\n", i, C0);
#endif

      mpz_add (t0, t0, t2);             /* t0 = 65 34 20 16 20 34 65 */

      mpz_sub (T, C2, C4);              /* T = 2*C_odd(1) = 0 2 0 2 0 2 0 */
      mpz_add (C2, C2, C4);             /* C2 = 2*C_even(1) */
      mpz_fdiv_q_2exp (C2, C2, 1);      /* C2 = C_even(1) */

      mpz_add (C4, t2, t4);             /* C4 = 2*C_even(2) */
      mpz_fdiv_q_2exp (C4, C4, 1);      /* C4 = C_even(2) */
      mpz_sub (t4, t2, t4);             /* t4 = 2*C_odd(2) */
      mpz_fdiv_q_2exp (t4, t4, 2);      /* t4 = C_odd(2)/2 = C_1 + 4*C_3 + 16*C_5 */
      mpz_fdiv_q_2exp (t2, T, 1);       /* t2 = C_odd(1) */

      mpz_sub (t0, t0, T);              /* t0 = 65 32 20 14 20 32 65 */
      mpz_mul_2exp (T, T, 4);
      mpz_sub (t0, t0, T);              /* t0 = 65 0 20 -18 20 0 65 */

      if (i < 2 * k - 1)
        {
          mpz_add (T, C0, C6);          /* T = C_0 + C_6 */
          mpz_sub (C2, C2, T);          /* C2 = C_2 + C_4 */
          mpz_sub (t0, t0, T);          /* t0 = 64 0 20 -18 20 0 64 */
          mpz_mul_2exp (T, T, 5);
        }
      else
        {
          mpz_sub (C2, C2, C0);         /* C2 = C_2 + C_4 */
          mpz_sub (t0, t0, C0);         /* t0 = 64 0 20 -18 20 0 */
          mpz_mul_2exp (T, C0, 5);
        }
      mpz_fdiv_q_2exp (t0, t0, 1);      /* t0 = 32 0 10 -9 10 0 32 */
      mpz_sub (t0, t0, T);              /* t0 = 0 0 10 -9 10 0 0 */
      mpz_sub (t0, t0, C2);             /* t0 = 0 0 9 -9 9 0 0 */
      mpz_divexact_ui (t0, t0, 9);      /* t0 = 0 0 1 -1 1 0 0 */
      mpz_sub (t0, C2, t0);             /* t0 = C_3 */
      mpz_sub (t2, t2, t0);             /* t2 = C_1 + C_5 */
      mpz_mul_2exp (T, t0, 2);          /* T = 4*C_3 */
      mpz_sub (t4, t4, T);              /* t4 = C_1 + 16*C_5 */
      mpz_sub (t4, t4, t2);             /* t4 = 15*C_5 */
      mpz_divexact_ui (t4, t4, 15);     /* t4 = C_5 */
      mpz_sub (t2, t2, t4);             /* t2 = C_1 */

      mpz_sub (C4, C4, C0);             /* C4 = 4*C_2 + 16*C_4 + 64*C_6 */
      mpz_fdiv_q_2exp (C4, C4, 2);      /* C4 = C_2 + 4*C_4 + 16*C_6 */
      if (i < 2 * k - 1)
        {
          mpz_mul_2exp (T, C6, 4);
          mpz_sub (C4, C4, T);          /* C4 = C_2 + 4*C_4 */
        }
      mpz_sub (C4, C4, C2);             /* C4 = 3*C_4 */
      mpz_divby3_1op (C4);              /* C4 = C_4 */
      mpz_sub (C2, C2, C4);             /* C2 = C_2 */
#ifdef DEBUG
      gmp_printf ("C_1[%d] = %Zd\n", i, t2);
      gmp_printf ("C_2[%d] = %Zd\n", i, C2);
      gmp_printf ("C_3[%d] = %Zd\n", i, t0);
      gmp_printf ("C_4[%d] = %Zd\n", i, C4);
      gmp_printf ("C_5[%d] = %Zd\n", i, t4);
      if (i < 2 * k - 1)
        gmp_printf ("C_6[%d] = %Zd\n", i, C6);
#endif
    }

  for (i = 0; i < l - 1; i++)
    mpz_add (C1, C1, t2);
  mpz_set (C1, t2);
  for (i = l; i < 2 * l - 1; i++)
    mpz_add (C1, C1, t2);

  for (i = 0; i < l - 1; i++)
    mpz_add (C3, C3, t0);
  mpz_set (C3, t0);
  for (i = l; i < 2 * l - 1; i++)
    mpz_add (C3, C3, t0);

  for (i = 0; i < l - 1; i++)
    mpz_add (C5, C5, t4);
  mpz_set (C5, t4);
  for (i = l; i < l + k - 1; i++)
    mpz_add (C5, C5, t4);
}
