/* The 'P+1' algorithm.

  Copyright 2002, 2003 Paul Zimmermann and Alexander Kruppa.

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
#include "gmp.h"
#include "ecm.h"

void             pp1_mul_ui (mpres_t, mpres_t, unsigned long, mpmod_t, 
                             mpres_t, mpres_t);
int  count_significant_bits (mp_limb_t);
void pp1_check_factor       (mpz_t, mpz_t);
int          pp1_stage1     (mpz_t, mpres_t, mpmod_t, double, double,
                             unsigned long *, mpz_t);

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
pp1_mul (mpres_t P1, mpres_t P0, mpz_t e, mpmod_t n, mpres_t P, mpres_t Q)
{
  unsigned long i, muls;
  mp_size_t size_e;

  if (mpz_cmp_ui (e, 0) == 0)
    {
      mpres_set_ui (P1, 2, n);
      return 0;
    }

  if (mpz_cmp_ui (e, 1) == 0)
    {
      mpres_set (P1, P0, n);
      return 0;
    }
  
  /* now e >= 2 */
  mpz_sub_ui (e, e, 1);
  mpres_mul (P, P0, P0, n);
  muls = 1;
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
              muls ++;
              mpres_sub (Q, Q, P0, n);
            }
          mpres_mul (P, P, P, n);
          muls ++;
          mpres_sub_ui (P, P, 2, n);
        }
      else /* k -> 2k */
        {
          mpres_mul (P, P, Q, n);
          muls ++;
          mpres_sub (P, P, P0, n);
          if (i) /* Q is not needed for last iteration */
            {
              mpres_mul (Q, Q, Q, n);
              muls ++;
              mpres_sub_ui (Q, Q, 2, n);
            }
        }
    }

  mpres_set (P1, P, n);
  
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
pp1_stage1 (mpz_t f, mpres_t P0, mpmod_t n, double B1, double B1done,
            unsigned long *muls, mpz_t orig_n)
{
  double B0, p, q, r, percentage;
  mpz_t g;
  mpres_t P, Q;
  mpres_t R, S, T;
  int youpi;
  unsigned int max_size, size_n;
  int Counter = 0, st, st_save;

  /* Prep for stage one counter */
  fprintf (stderr, "1:000 \r");

  mpz_init (g);
  mpres_init (P, n);
  mpres_init (Q, n);
  mpres_init (R, n);
  mpres_init (S, n);
  mpres_init (T, n);

  B0 = ceil (sqrt (B1));

  size_n = mpz_sizeinbase (n->orig_modulus, 2);
  max_size = L1 * size_n;

  /* suggestion from Peter Montgomery: start with exponent n^2-1,
     as factors of Lucas and Fibonacci number are either +/-1 (mod index),
     and so is n. Therefore, index will appear as a factor
     of n^2-1 and be included in stage 1.
     Do this only when n is composite, otherwise all tests with prime
     n factor of a Cunningham number will succeed in stage 1.

     As in P-1, for small overhead, use that trick only when lg(n) <= sqrt(B1).
  */
  if ((double) size_n <= B0 && mpz_probab_prime_p (n->orig_modulus, 1) == 0)
    {
      mpz_mul (g, n->orig_modulus, n->orig_modulus);
      mpz_sub_ui (g, g, 1);
      *muls += pp1_mul (P0, P0, g, n, P, Q);
    }

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
      /* No need to save incrementals here, or to show screen output, since this happens pretty quickly */
    }

  *muls += pp1_mul (P0, P0, g, n, P, Q);

  /* update the screen mode after the cascade work is done */
  percentage = p;
  st_save = st = cputime ();
  percentage /= B1;
  percentage *= 100;
  fprintf (stderr, "1:%03d\r", (int)percentage);

  /* then all primes > sqrt(B1) and taken with exponent 1 */
  for (; p <= B1; p = getprime(p))
    {
      if (p > B1done)
        *muls += pp1_mul_prac (P0, (unsigned long) p, n, P, Q, R, S, T);
      if (++Counter == 3000)
	{
	  int st_now = cputime();
	  if (st_now - st > 30000)
	    {
	      /* Check to see if we should update our screen "percentage counter" */
	      percentage = p;
	      st = st_now;
	      percentage /= B1;
	      percentage *= 100;
	      fprintf (stderr, "1:%03d\r", (int)percentage);
  	    }
	  Counter=0;
	    /* should we save the current "ecm_wip.sav" file??? It is saved every 15 minutes */
	    /* NOTE this saving DOES NOT save the expression.  It is just a "fail-safe" measure */
  #if defined (DEBUG_AUTO_SAVE)
	  if (st_now - st_save > 2000)
  #else
	  if (st_now - st_save > 15 * 60 * 1000)
  #endif
	    {
	      st_save = st_now;
	      /* Figure out how to save here */
	      /* sigma is not needed */
	      /* A is not needed */

	      /******** THIS NEEDS TO BE MADE TO SAVE *****************/
	      /* write_temp_resumefile (PM1_METHOD, p+1, sigma, A, x, orig_n, verbose); */
	    }
	    /*  This "testing" code is here to see just how often this ++Counter loop is entered.
	    {
	      static int x;
  	      fprintf (stderr, "1:%03d  p=%.0f\r", ++x, p);
	    }
	    */
	}
    }

  getprime (0.0); /* free the prime tables, and reinitialize */

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

/* put in seed a valid random seed for P+1 */
void
pp1_random_seed (mpz_t seed, mpz_t n, gmp_randstate_t randstate)
{
  mpz_t q;

  /* need gcd(p^2-4, n) = 1. */
  mpz_init (q);
  do
    {
      mpz_urandomb (q, randstate, 32);
      mpz_add_ui (q, q, 1);
      mpz_set (seed, q);
      mpz_mul (q, q, q);
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
*                                  Stage 2                                    *
*                                                                             *
******************************************************************************/

/* puts in F[0..dF-1] the successive values of 

   V_j(P) for Williams P+1

   For P+1, we have V_{j+6} = V_j * V_6 - V_{j-6}.

   for 0 < j = 1 mod 6 < d, j and d coprime.
*/

int
pp1_rootsF (listz_t F, unsigned int d, mpres_t *x, listz_t t,
            mpmod_t modulus, int verbose, unsigned long *tot_muls)
{
  unsigned int i, j, muls = 0;
  int st, st2;
  mpres_t fd[4];
  
  st = cputime ();

  mpres_get_z (F[0], *x, modulus); /* V_1(P)=P for P+1 */
  i = 1;

  if (d > 7)
    {
      st2 = cputime ();
      mpres_init (fd[0], modulus);
      mpres_init (fd[1], modulus);
      mpres_init (fd[2], modulus);
      mpres_init (fd[3], modulus);
      mpz_set_ui (*t, 7);
      pp1_mul (fd[0], *x, *t, modulus, fd[2], fd[3]);
      mpz_set_ui (*t, 6);
      pp1_mul (fd[1], *x, *t, modulus, fd[2], fd[3]);
      mpres_set (fd[2], *x, modulus); 
      /* for P+1, fd[0] = V_7(P), fd[1] = V_6(P), fd[2] = V_{7-6}(P) */
      if (verbose >= 2)
        printf ("Initializing table of differences for F took %dms\n", cputime () - st2);
      j = 7;
      while (j < d)
        {
          if (gcd (j, d) == 1)
            mpres_get_z (F[i++], fd[0], modulus);

          mpres_swap (fd[0], fd[2], modulus);
          mpres_mul (fd[3], fd[2], fd[1], modulus);
          mpres_sub (fd[0], fd[3], fd[0], modulus);
          j += 6;
          muls++;
        }
      mpres_clear (fd[0], modulus);
      mpres_clear (fd[1], modulus);
      mpres_clear (fd[2], modulus);
      mpres_clear (fd[3], modulus);
    }

  if (verbose >= 2)
    printf ("Computing roots of F took %dms and %d muls\n", cputime () - st, 
            muls);
  
  if (tot_muls != NULL)
    *tot_muls += muls;
  
  return 0;
}

mpres_t *
pp1_rootsG_init (mpres_t *x, unsigned long s, unsigned int d, mpmod_t modulus)
{
  int st;
  mpres_t *fd, P;
  mpz_t t;
  
  st = cputime ();
  
  mpz_init (t);

  fd = (mpres_t *) malloc (4 * sizeof (mpres_t));

  mpres_init (fd[0], modulus);
  mpres_init (fd[1], modulus);
  mpres_init (fd[2], modulus);
  mpres_init (fd[3], modulus);
  mpres_init (P, modulus);

  mpz_set_ui (t, s);
  pp1_mul (fd[0], *x, t, modulus, fd[3], P);
  mpz_set_ui (t, d);
  pp1_mul (fd[1], *x, t, modulus, fd[3], P);
  mpz_set_ui (t, (s > d) ? s - d : d - s);
  pp1_mul (fd[2], *x, t, modulus, fd[3], P);
  /* for P+1, fd[0] = V_s(P), fd[1] = V_d(P), fd[2] = V_{|s-d|}(P) */

  mpres_clear (P, modulus);
  mpz_clear (t);

  return fd;
}

void 
pp1_rootsG_clear (mpres_t *fd, mpmod_t modulus)
{
  mpres_clear (fd[0], modulus);
  mpres_clear (fd[1], modulus);
  mpres_clear (fd[2], modulus);
  mpres_clear (fd[3], modulus);
  free (fd);
}

int
pp1_rootsG (listz_t G, unsigned int d, mpres_t *fd, mpmod_t modulus, 
            unsigned long *tot_muls) 
{
  unsigned int i;
  int st;
  
  st = cputime ();

  for (i = 0; i < d; i++)
    {
      mpres_get_z (G[i], fd[0], modulus);

      mpres_swap (fd[0], fd[2], modulus);
      mpres_mul (fd[3], fd[2], fd[1], modulus);
      mpres_sub (fd[0], fd[3], fd[0], modulus);
    }
  
  if (tot_muls != NULL)
    *tot_muls += d;
  
  return 0;
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
pp1 (mpz_t f, mpz_t p, mpz_t n, double B1done, double B1, double B2min, double B2,
     double B2scale, unsigned int k, unsigned int S, int verbose, int repr)
{
  int youpi = 0, st;
  mpres_t a;
  mpmod_t modulus;
  unsigned long muls = 0;

  st = cputime ();

  /* Set default B2. See ecm.c for comments */
  if (B2 == 0.0)
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
      if (B1done == 1.0)
        printf("B1=%1.0f", B1);
      else
        printf("B1=%1.0f-%1.0f", B1done, B1);
      if (B2min <= B1)
        printf(", B2=%1.0f", B2);
      else
        printf(", B2=%1.0f-%1.0f", B2min, B2);

      printf (", polynomial x^1");
      if (B1done == 1.0 || verbose > 1) /* don't print x0 in resume case */
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
    youpi = pp1_stage1 (f, a, modulus, B1, B1done, &muls, n);

  if (verbose >= 1)
    {
      printf ("Step 1 took %dms", cputime () - st);
      if (verbose >= 2)
	printf (" for %lu muls", muls);
      printf ("\n");
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

  youpi = stage2 (f, &a, modulus, B2min, B2, k, S, verbose, PP1_METHOD);

 end:
  if (youpi != 0) 
    {
      pp1_check_factor (p, f);
    }
  mpres_get_z (p, a, modulus);
  mpres_clear (a, modulus);
  mpmod_clear (modulus);

  return youpi;
}
