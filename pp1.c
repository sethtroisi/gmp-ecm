/* The 'P+1' algorithm.

  Copyright 2002, 2003, 2004, 2005 Paul Zimmermann and Alexander Kruppa.

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
#include <stdlib.h>
#include <stdio.h>
#include <limits.h> /* for ULONG_MAX */
#include "gmp.h"
#include "ecm.h"
#include "ecm-impl.h"

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
*/
static void
pp1_mul (mpres_t P1, mpres_t P0, mpz_t e, mpmod_t n, mpres_t P, mpres_t Q)
{
  unsigned long i;
  mp_size_t size_e;

  if (mpz_cmp_ui (e, 0) == 0)
    {
      mpres_set_ui (P1, 2, n);
      return;
    }

  if (mpz_cmp_ui (e, 1) == 0)
    {
      mpres_set (P1, P0, n);
      return;
    }
  
  /* now e >= 2 */
  mpz_sub_ui (e, e, 1);
  mpres_mul (P, P0, P0, n);
  mpres_sub_ui (P, P, 2, n); /* P = V_2(P0) = P0^2-2 */
  mpres_set (Q, P0, n);      /* Q = V_1(P0) = P0 */

  /* invariant: (P, Q) = (V_{k+1}(P0), V_k(P0)), start with k=1 */
  size_e = mpz_sizeinbase (e, 2);
  for (i = size_e - 1; i > 0;)
    {
      if (mpz_tstbit (e, --i)) /* k -> 2k+1 */
        {
          if (i) /* Q is not needed for last iteration */
            {
              mpres_mul (Q, P, Q, n);
              mpres_sub (Q, Q, P0, n);
            }
          mpres_mul (P, P, P, n);
          mpres_sub_ui (P, P, 2, n);
        }
      else /* k -> 2k */
        {
          mpres_mul (P, P, Q, n);
          mpres_sub (P, P, P0, n);
          if (i) /* Q is not needed for last iteration */
            {
              mpres_mul (Q, Q, Q, n);
              mpres_sub_ui (Q, Q, 2, n);
            }
        }
    }

  mpres_set (P1, P, n);
}

/* Input:  P0 is the initial point (sigma)
           n is the number to factor
           B1 is the stage 1 bound
	   B1done: stage 1 was already done up to that limit
	   go: if <> 1, group order to preload
   Output: a is the factor found, or the value at end of stage 1
   Return value: non-zero iff a factor was found.
*/
static int
pp1_stage1 (mpz_t f, mpres_t P0, mpmod_t n, double B1, double B1done, mpz_t go)
{
  double B0, p, q, r;
  mpz_t g;
  mpres_t P, Q;
  mpres_t R, S, T;
  int youpi;
  unsigned int max_size, size_n;

  mpz_init (g);
  mpres_init (P, n);
  mpres_init (Q, n);
  mpres_init (R, n);
  mpres_init (S, n);
  mpres_init (T, n);

  B0 = ceil (sqrt (B1));

  size_n = mpz_sizeinbase (n->orig_modulus, 2);
  max_size = L1 * size_n;

  if (mpz_cmp_ui (go, 1) > 0)
    pp1_mul (P0, P0, go, n, P, Q);

  /* suggestion from Peter Montgomery: start with exponent n^2-1,
     as factors of Lucas and Fibonacci number are either +/-1 (mod index),
     and so is n. Therefore, index will appear as a factor
     of n^2-1 and be included in stage 1.
     Do this only when n is composite, otherwise all tests with prime
     n factor of a Cunningham number will succeed in stage 1.

     As in P-1, for small overhead, use that trick only when lg(n) <= sqrt(B1).
  */
  if ((double) size_n <= B0 &&
      mpz_probab_prime_p (n->orig_modulus, PROBAB_PRIME_TESTS) == 0)
    {
      mpz_mul (g, n->orig_modulus, n->orig_modulus);
      mpz_sub_ui (g, g, 1);
      pp1_mul (P0, P0, g, n, P, Q);
    }

  mpz_set_ui (g, 1);

  /* first loop through small primes <= sqrt(B1) */
  for (p = 2.0; p <= B0; p = getprime (p))
    {
      for (q = 1, r = p; r <= B1; r *= p)
        if (r > B1done) q *= p;
      mpz_mul_d (g, g, q, Q);
      if (mpz_sizeinbase (g, 2) >= max_size)
	{
	  pp1_mul (P0, P0, g, n, P, Q);
	  mpz_set_ui (g, 1);
	}
      /* No need to save incrementals here, or to show screen output, since this happens pretty quickly */
    }

  pp1_mul (P0, P0, g, n, P, Q);

  /* then all primes > sqrt(B1) and taken with exponent 1 */
  /* since pp1_mul_prac takes an unsigned long, we have to check
     that B1 <= MAX_ULONG */
  if (B1 > (double) ULONG_MAX)
    {
      fprintf (stderr, "Error, maximal step1 bound B1 for P+1 is %lu\n", ULONG_MAX);
      exit (EXIT_FAILURE);
    }
  for (; p <= B1; p = getprime (p))
    {
      if (p > B1done)
        pp1_mul_prac (P0, (unsigned long) p, n, P, Q, R, S, T);
    }

  getprime (FREE_PRIME_TABLE); /* free the prime tables, and reinitialize */

  mpres_clear (Q, n);
  mpres_clear (R, n);
  mpres_clear (S, n);
  mpres_clear (T, n);
  mpz_clear (g);

  mpres_sub_ui (P, P0, 2, n);
  mpres_gcd (f, P, n);
  youpi = mpz_cmp_ui (f, 1);

  mpres_clear (P, n);
  
  return youpi;
}

/* checks if the factor p was found by P+1 or P-1 (when prime).
   a is the initial seed.
*/
static void
pp1_check_factor (mpz_t a, mpz_t p)
{
  if (mpz_probab_prime_p (p, PROBAB_PRIME_TESTS))
    {
      mpz_mul (a, a, a);
      mpz_sub_ui (a, a, 4);
      if (mpz_jacobi (a, p) == 1)
        printf ("[factor found by P-1]\n");
    }
}

/******************************************************************************
*                                                                             *
*                                  Stage 2                                    *
*                                                                             *
******************************************************************************/

/* puts in F[0..dF-1] the successive values of 
   V_j(P) for Williams P+1
   For P+1, we have V_{j+6} = V_j * V_6 - V_{j-6}.
   for 0 < j = 1 mod 6 < d, j and d coprime.
   Return non-zero iff a factor was found.
*/
int
pp1_rootsF (listz_t F, unsigned int d1, unsigned int d2, unsigned int dF, 
            mpres_t *x, listz_t t, mpmod_t modulus, int verbose)
{
  unsigned int i, j, muls = 0;
  int st, st2;
  mpres_t fd[5]; /* fd[3..4] are temp vars */
  
  if (dF == 0)
    return 0;

  st2 = st = cputime ();

  if (verbose >= 3)
    printf ("pp1_rootsF: d1 = %d, d2 = %d, dF = %d\n", d1, d2, dF);

  mpres_init (fd[0], modulus);
  mpres_init (fd[1], modulus);
  mpres_init (fd[2], modulus);
  mpres_init (fd[3], modulus);
  mpres_init (fd[4], modulus);

  mpz_set_ui (*t, d2);
  pp1_mul (fd[2], *x, *t, modulus, fd[3], fd[4]);
  mpres_get_z (F[0], fd[2], modulus);
  
  mpz_set_ui (*t, 7);
  pp1_mul (fd[0], fd[2], *t, modulus, fd[3], fd[4]);

  mpz_set_ui (*t, 6); 
  pp1_mul (fd[1], fd[2], *t, modulus, fd[3], fd[4]);
  /* for P+1, fd[0] = V_{7*d2}(P), fd[1] = V_{6*d2}(P), fd[2] = V_{d2}(P) */

  if (verbose >= 2)
    printf ("Initializing table of differences for F took %dms\n", cputime () - st2);
  i = 1;
  j = 7;
  while (i < dF)
    {
      if (gcd (j, d1) == 1)
        mpres_get_z (F[i++], fd[0], modulus);

      /* V_{m+n} = V_m * V_n - V_{m-n} */
      /* fd[0] = V_m, fd[1] = V_n, fd[2] = V_{m-n} */
      mpres_swap (fd[0], fd[2], modulus);
      mpres_mul (fd[3], fd[2], fd[1], modulus);
      mpres_sub (fd[0], fd[3], fd[0], modulus);
      /* fd[0] = V_{m+n}, fd[1] = V_n, fd[2] = V_m */
      j += 6;
      muls++;
    }
  mpres_clear (fd[0], modulus);
  mpres_clear (fd[1], modulus);
  mpres_clear (fd[2], modulus);
  mpres_clear (fd[3], modulus);
  mpres_clear (fd[4], modulus);

  if (verbose >= 2)
    {
      printf ("Computing roots of F took %dms", cputime () - st);
      if (verbose > 2)
        printf (" and %d muls", muls);
      printf ("\n");
    }
  
  return 0;
}

pp1_roots_state *
pp1_rootsG_init (mpres_t *x, double s, unsigned int d, unsigned int d2, 
                 mpmod_t modulus)
{
  int st;
  mpres_t P;
  mpz_t t;
  pp1_roots_state *state;
  
  st = cputime ();
  
  state = (pp1_roots_state *) xmalloc (sizeof (pp1_roots_state));
  mpz_init (t);

  mpres_init (state->fd[0], modulus);
  mpres_init (state->fd[1], modulus);
  mpres_init (state->fd[2], modulus);
  mpres_init (state->fd[3], modulus);
  mpres_init (P, modulus);

  state->dsieve = d2;
  state->rsieve = 0; /* Assumes (d1*d2)|s */

  mpz_set_d (t, s);
  pp1_mul (state->fd[0], *x, t, modulus, state->fd[3], P);
  mpz_set_ui (t, d);
  pp1_mul (state->fd[1], *x, t, modulus, state->fd[3], P);
  mpz_set_d (t, fabs (s - (double)d));
  pp1_mul (state->fd[2], *x, t, modulus, state->fd[3], P);
  /* for P+1, fd[0] = V_s(P), fd[1] = V_d(P), fd[2] = V_{|s-d|}(P) */

  mpres_clear (P, modulus);
  mpz_clear (t);

  return state;
}

void 
pp1_rootsG_clear (pp1_roots_state *state, UNUSED mpmod_t modulus)
{
  mpres_clear (state->fd[0], modulus);
  mpres_clear (state->fd[1], modulus);
  mpres_clear (state->fd[2], modulus);
  mpres_clear (state->fd[3], modulus);
  free (state);
}

int
pp1_rootsG (listz_t G, unsigned int d, pp1_roots_state *state, 
            mpmod_t modulus, int verbose) 
{
  unsigned int i;
  int st;
  
  st = cputime ();

  for (i = 0; i < d;)
    {
      if (gcd (state->rsieve, state->dsieve) == 1)
         mpres_get_z (G[i++], state->fd[0], modulus);

      mpres_swap (state->fd[0], state->fd[2], modulus);
      mpres_mul (state->fd[3], state->fd[2], state->fd[1], modulus);
      mpres_sub (state->fd[0], state->fd[3], state->fd[0], modulus);
      state->rsieve++;
    }

  if (verbose >= 2)
    {
      printf ("Computing roots of G took %dms", cputime () - st);
      if (verbose > 2)
        printf (", %u muls", d);
      printf ("\n");
    }

  
  return 0;
}


/******************************************************************************
*                                                                             *
*                               Williams P+1                                  *
*                                                                             *
******************************************************************************/

/* Input: p is the initial generator (sigma), if 0 generate it at random.
          n is the number to factor
	  B1 is the stage 1 bound
	  B2 is the stage 2 bound
          k is the number of blocks for stage 2
          verbose is the verbose level: 0=quiet, 1=normal, 2=verbose
   Output: p is the factor found
   Return value: non-zero iff a factor is found (1 for stage 1, 2 for stage 2)
*/
int
pp1 (mpz_t f, mpz_t p, mpz_t n, mpz_t go, double B1done, double B1,
     double B2min, double B2, double B2scale, unsigned int k, unsigned int S,
     int verbose, int repr)
{
  int youpi = 0, st;
  mpres_t a;
  mpmod_t modulus;

  st = cputime ();

  if (mpz_cmp_ui (p, 0) == 0)
    {
      gmp_randstate_t state;
      gmp_randinit_default (state);
      pm1_random_seed (p, n, state);
      gmp_randclear (state);
    }

  /* Set default B2. See ecm.c for comments */
  if (IS_DEFAULT_B2(B2))
    B2 = pow (2.0 * B1 / 6.0, 1.424828748);

  /* Scale B2 by what the user said (or by the default scaling of 1.0) */
  B2 *= B2scale;

  /* Set default degree for Brent-Suyama extension */

  if (S == 0)
    S = 1;

  if (S != 1)
    {
      printf ("Warning: Brent-Suyama's extension does not work with P+1, using x^1\n");
      S = 1;
    }
  
  if (verbose >= 1)
    {
      printf ("Using ");
      if (IS_DEFAULT_B1_DONE(B1done))
        printf("B1=%1.0f", B1);
      else
        printf("B1=%1.0f-%1.0f", B1done, B1);
      if (B2min <= B1)
        printf(", B2=%1.0f", B2);
      else
        printf(", B2=%1.0f-%1.0f", B2min, B2);

      printf (", polynomial x^1");
      if (IS_DEFAULT_B1_DONE(B1done) || verbose > 1) /* don't print x0 in resume case */
	{
	  printf (", x0=");
	  mpz_out_str (stdout, 10, p);
	}
      printf ("\n");
      fflush (stdout);
    }

  if (repr == 1)
    mpmod_init_MPZ (modulus, n);
  else   if (repr == 2)
    mpmod_init_MODMULN (modulus, n);
  else if (repr == 3)
    mpmod_init_REDC (modulus, n);
  else if (repr > 16)
    mpmod_init_BASE2 (modulus, repr, n);
  else
    mpmod_init (modulus, n, repr, verbose);

  mpres_init (a, modulus);
  mpres_set_z (a, p, modulus);

  if (B1 > B1done)
    youpi = pp1_stage1 (f, a, modulus, B1, B1done, go);

  st = cputime () - st;

  if (verbose >= 1)
    {
      printf ("Step 1 took %dms\n", st);
      fflush (stdout);
    }

  if (verbose >= 2)
    {
      printf ("x=");
      mpres_out_str (stdout, 10, a, modulus);
      printf ("\n");
      fflush (stdout);
    }

  if (youpi != 0) /* a factor was found */
    goto end;

  youpi = stage2 (f, &a, modulus, B2min, B2, k, S, verbose, PP1_METHOD, st);

 end:
  if (youpi != 0 && verbose != 0)
    pp1_check_factor (p, f); /* tell user if factor was found in fact by P-1 */

  mpres_get_z (p, a, modulus);
  mpres_clear (a, modulus);
  mpmod_clear (modulus);

  return youpi;
}
