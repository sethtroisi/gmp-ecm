/* Williams 'P+1' algorithm.

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

#include <math.h>
#include <stdio.h>
#include "gmp.h"
#include "ecm.h"

void             pp1_mul_ui (mpz_t, mpz_t, unsigned long, mpz_t, mpz_t, mpz_t);
int  count_significant_bits (mp_limb_t);
void pp1_check_factor       (mpz_t, mpz_t);
int          pp1_stage1     (mpz_t, mpz_t, mpz_t, double, double,
                             unsigned long *);

/******************************************************************************
*                                                                             *
*                                  Stage 1                                    *
*                                                                             *
******************************************************************************/

/* prime powers are accumulated up to about n^L1 */
#define L1 1

/* P1 <- V_e(P0), using P, Q as auxiliary variables,
   where V_{2k}(P0) = V_k(P0)^2 - 2
         V_{2k-1}(P0) = V_k(P0)*V_{k-1}(P0) - P0.
   (More generally V_{m+n} = V_m * V_n - V_{m-n}.)
   Warning: P1 and P0 may be equal.
   Return the number of multiplies mod n.
*/
int
pp1_mul (mpz_t P1, mpz_t P0, mpz_t e, mpz_t n, mpz_t P, mpz_t Q)
{
  unsigned long i, muls;
  mp_size_t size_e;

  if (mpz_cmp_ui (e, 0) == 0)
    {
      mpz_set_ui (P1, 2);
      return 0;
    }

  if (mpz_cmp_ui (e, 1) == 0)
    {
      mpz_set (P1, P0);
      return 0;
    }

  /* now e >= 2 */
  mpz_sub_ui (e, e, 1);
  mpz_mul (P, P0, P0);
  muls = 1;
  mpz_sub_ui (P, P, 2);
  mpz_mod (P, P, n); /* P = V_2(P0) = P0^2-2 */
  mpz_set (Q, P0);   /* Q = V_1(P0) = P0 */

  /* invariant: (P, Q) = (V_{k+1}(P0), V_k(P0)), start with k=1 */
  size_e = mpz_sizeinbase (e, 2);
  for (i = size_e - 1; i > 0;)
    {
      if (mpz_tstbit (e, --i)) /* k -> 2k+1 */
        {
          if (i) /* Q is not needed for last iteration */
            {
              mpz_mul (Q, P, Q);
              muls ++;
              mpz_sub (Q, Q, P0);
              mpz_mod (Q, Q, n);
            }
          mpz_mul (P, P, P);
          muls ++;
          mpz_sub_ui (P, P, 2);
        }
      else /* k -> 2k */
        {
          mpz_mul (P, P, Q);
          muls ++;
          mpz_sub (P, P, P0);
          if (i) /* Q is not needed for last iteration */
            {
              mpz_mul (Q, Q, Q);
              muls ++;
              mpz_sub_ui (Q, Q, 2);
              mpz_mod (Q, Q, n);
            }
        }
      if (i)
        mpz_mod (P, P, n);
    }

  mpz_mod (P1, P, n);

  /* the number of multiplies is 1 + 2*#loops - 1 (last loop) */
  return 2 * (size_e - 1);
}

int
count_significant_bits (mp_limb_t e)
{
  int i = 0;
  while (e)
    {
      e >>= 1;
      i ++;
    }
  return i;
}

/* Input:  P0 is the initial point (sigma)
           n is the number to factor
           B1 is the stage 1 bound
   Output: a is the factor found, or the value at end of stage 1
   Return value: non-zero iff a factor was found.
*/
int
pp1_stage1 (mpz_t f, mpz_t P0, mpz_t n, double B1, double B1done,
            unsigned long *muls)
{
  double B0, p, q, r;
  mpz_t g;
  mpres_t P, Q;
  mpres_t R, S, T;
  int youpi;
  unsigned int max_size;

  mpz_init (g);
  mpres_init (P, nn);
  mpres_init (Q, nn);
  mpres_init (R, nn);
  mpres_init (S, nn);
  mpres_init (T, nn);

  B0 = sqrt (B1);

  max_size = L1 * mpz_sizeinbase (n, 2);

  /* suggestion from Peter Montgomery: start with exponent n-1,
     since any prime divisor of b^m-1 which does not divide any
     algebraic factor of b^m-1 must be of the form km+1 [Williams82].
     Do this only when n is composite, otherwise all tests with prime
     n factor of a Cunningham number will succeed in stage 1. */
  if (mpz_probab_prime_p (n, 1) == 0)
    {
      mpz_sub_ui (g, n, 1);
      *muls += pp1_mul (P0, P0, g, n, P, Q);
    }
  else
    mpz_set_ui (g, 1);

  /* first loop through small primes <= sqrt(B1) */
  for (p = 2.0; p <= B0; p = getprime(p))
    {
      for (q = 1, r = p; r <= B1; r *= p)
        if (r > B1done) q *= p;
      mpz_mul_d (g, g, q, Q);
      if (mpz_sizeinbase (g, 2) >= max_size)
	{
	  *muls += pp1_mul (P0, P0, g, n, P, Q);
	  mpz_set_ui (g, 1);
	}
    }

  *muls += pp1_mul (P0, P0, g, n, P, Q);

  /* then all primes > sqrt(B1) and taken with exponent 1 */
  for (; p <= B1; p = getprime(p))
    if (p > B1done)
      *muls += pp1_mul_prac (P0, (unsigned long) p, n, P, Q, R, S, T);

  getprime (0.0); /* free the prime tables, and reinitialize */

  mpres_clear (P, nn);
  mpres_clear (Q, nn);
  mpres_clear (R, nn);
  mpres_clear (S, nn);
  mpres_clear (T, nn);

  mpz_sub_ui (g, P0, 2);
  mpz_gcd (f, g, n);
  youpi = mpz_cmp_ui (f, 1);

  mpz_clear (g);
  
  return youpi;
}

/* put in seed a valid random seed for P+1 */
void
pp1_random_seed (mpz_t seed, mpz_t n, gmp_randstate_t randstate)
{
  mpz_t q;

  /* need gcd(p^2-4, n) = 1. */
  mpz_init (q);
  do
    {
      mpz_urandomb (seed, randstate, 32);
      mpz_add_ui (seed, seed, 1);
      mpz_mul (q, seed, seed);
      mpz_sub_ui (q, q, 4);
      mpz_gcd (q, q, n);
    }
  while (mpz_cmp_ui (q, 1) != 0);
  mpz_clear (q);
}

/* checks if the factor p was found by P+1 or P-1 (when prime).
   a is the initial seed.
*/
void
pp1_check_factor (mpz_t a, mpz_t p)
{
  if (mpz_probab_prime_p (p, 25))
    {
      mpz_mul (a, a, a);
      mpz_sub_ui (a, a, 4);
      if (mpz_jacobi (a, p) == 1)
        printf ("[factor found by P-1]\n");
    }
}

/******************************************************************************
*                                                                             *
*                               Williams P+1                                  *
*                                                                             *
******************************************************************************/

/* Input: p is the initial generator (sigma)
          n is the number to factor
	  B1 is the stage 1 bound
	  B2 is the stage 2 bound
          k is the number of blocks for stage 2
          verbose is the verbose level: 0=quiet, 1=normal, 2=verbose
   Output: p is the factor found
   Return value: non-zero iff a factor is found (1 for stage 1, 2 for stage 2)
*/
int
pp1 (mpz_t f, mpz_t p, mpz_t n, double B1, double B2, double B1done, 
     unsigned int k, unsigned int S, int verbose)
{
  int youpi = 0, st;
  mpz_t a;
  curve P;
  unsigned long muls = 0;

  st = cputime ();

  mpz_init (a);
  mpz_set (a, p);

  if (verbose >= 1)
    {
      printf ("Using seed=");
      mpz_out_str (stdout, 10, a);
      printf ("\n");
      fflush (stdout);
    }

  if (B1 > B1done)
    youpi = pp1_stage1 (f, p, n, B1, B1done, &muls);

  if (verbose >= 1)
    {
      printf ("Stage 1 took %dms for %lu muls\n", cputime() - st, muls);
      fflush (stdout);
    }

  if (youpi != 0) /* a factor was found */
    goto end;

  mpz_init_set (P.x, p);
  youpi = stage2 (f, &P, n, B2, k, S, verbose, PP1_METHOD, B1);
  mpz_clear (P.x);

 end:
  if (youpi != 0)
    pp1_check_factor (a, f);
  mpz_clear (a);

  return youpi;
}
