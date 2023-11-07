/* The 'P+1' algorithm.

Copyright 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012
Paul Zimmermann and Alexander Kruppa.

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

#include <math.h>
#include <stdlib.h>
#include <primesieve.h>
#include "ecm-impl.h"

#ifdef HAVE_LIMITS_H
# include <limits.h>
#else
# ifndef ULONG_MAX
#  define ULONG_MAX __GMP_ULONG_MAX
# endif
#endif


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
   Assume e >= 1.
*/
static void
pp1_mul (mpres_t P1, mpres_t P0, mpz_t e, mpmod_t n, mpres_t P, mpres_t Q)
{
  mp_size_t size_e, i;

  ASSERT (mpz_cmp_ui (e, 1) >= 0);

  if (mpz_cmp_ui (e, 1) == 0)
    {
      mpres_set (P1, P0, n);
      return;
    }
  
  /* now e >= 2 */
  mpz_sub_ui (e, e, 1);
  mpres_sqr (P, P0, n);
  mpres_sub_ui (P, P, 2, n); /* P = V_2(P0) = P0^2-2 */
  mpres_set (Q, P0, n);      /* Q = V_1(P0) = P0 */

  /* invariant: (P, Q) = (V_{k+1}(P0), V_k(P0)), start with k=1 */
  size_e = mpz_sizeinbase (e, 2);
  for (i = size_e - 1; i > 0;)
    {
      if (ecm_tstbit (e, --i)) /* k -> 2k+1 */
        {
          if (i) /* Q is not needed for last iteration */
            {
              mpres_mul (Q, P, Q, n);
              mpres_sub (Q, Q, P0, n);
            }
          mpres_sqr (P, P, n);
          mpres_sub_ui (P, P, 2, n);
        }
      else /* k -> 2k */
        {
          mpres_mul (P, P, Q, n);
          mpres_sub (P, P, P0, n);
          if (i) /* Q is not needed for last iteration */
            {
              mpres_sqr (Q, Q, n);
              mpres_sub_ui (Q, Q, 2, n);
            }
        }
    }

  mpres_set (P1, P, n);
  mpz_add_ui (e, e, 1); /* recover original value of e */
}

/* Input:  P0 is the initial point (sigma)
           n is the number to factor
           B1 is the stage 1 bound
	   B1done: stage 1 was already done up to that limit
	   go: if <> 1, group order to preload
   Output: a is the factor found, or the value at end of stage 1
	   B1done is set to B1 if stage 1 completed normally,
	   or to the largest prime processed if interrupted, but never
	   to a smaller value than B1done was upon function entry.
   Return value: non-zero iff a factor was found.
*/
static int
pp1_stage1 (mpz_t f, mpres_t P0, mpmod_t n, double B1, double *B1done, 
            mpz_t go, int (*stop_asap)(void), char *chkfilename)
{
  double B0, p, q, r, last_chkpnt_p;
  mpz_t g;
  mpres_t P, Q;
  mpres_t R, S, T;
  int youpi = ECM_NO_FACTOR_FOUND;
  unsigned int max_size, size_n;
  long last_chkpnt_time;
  primesieve_iterator it;

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

  last_chkpnt_p = 2.;
  last_chkpnt_time = cputime ();
  /* first loop through small primes <= sqrt(B1) */
  primesieve_init (&it);
  for (p = 2.0; p <= B0; p = (double) primesieve_next_prime (&it))
    {
      for (q = 1, r = p; r <= B1; r *= p)
        if (r > *B1done) q *= p;
      mpz_mul_d (g, g, q, Q);
      if (mpz_sizeinbase (g, 2) >= max_size)
	{
	  pp1_mul (P0, P0, g, n, P, Q);
	  mpz_set_ui (g, 1);
          if (stop_asap != NULL && (*stop_asap) ())
            {
            interrupt:
              outputf (OUTPUT_NORMAL, "Interrupted at prime %.0f\n", p);
	      if (p > *B1done)
		  *B1done = p;
              goto clear_and_exit;
            }
	}
    }

  pp1_mul (P0, P0, g, n, P, Q);

  /* All primes sqrt(B1) < p <= B1 appear with exponent 1. All primes <= B1done
     are already included with exponent at least 1, so it's safe to skip 
     ahead to B1done+1. */
 
  if (p <= *B1done) {
    primesieve_skipto (&it, *B1done, primesieve_get_max_stop());
    p = primesieve_next_prime (&it);
    ASSERT( p > *B1done );
  }

  /* then all primes > sqrt(B1) and taken with exponent 1 */
  for (; p <= B1; p = (double) primesieve_next_prime (&it))
    {
      pp1_mul_prac (P0, (ecm_uint) p, n, P, Q, R, S, T);
  
      if (stop_asap != NULL && (*stop_asap) ())
        goto interrupt;
      if (chkfilename != NULL && p > last_chkpnt_p + 10000. &&
          elltime (last_chkpnt_time, cputime ()) > CHKPNT_PERIOD)
        {
	  writechkfile (chkfilename, ECM_PP1, p, n, NULL, P0, NULL, NULL);
          last_chkpnt_p = p;
          last_chkpnt_time = cputime ();
        }
    }

  /* If stage 1 finished normally, p is the smallest prime >B1 here.
     In that case, set to B1 */
  if (p > B1)
    p = B1;
  
  if (p > *B1done)
    *B1done = p;
  
  mpres_sub_ui (P, P0, 2, n);
  mpres_gcd (f, P, n);
  youpi = mpz_cmp_ui (f, 1);

clear_and_exit:
  if (chkfilename != NULL)
    writechkfile (chkfilename, ECM_PP1, p, n, NULL, P0, NULL, NULL);
  primesieve_free_iterator (&it); /* free the prime iterator */
  mpres_clear (Q, n);
  mpres_clear (R, n);
  mpres_clear (S, n);
  mpres_clear (T, n);
  mpz_clear (g);
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
        outputf (OUTPUT_NORMAL, "[factor found by P-1]\n");
    }
}


/******************************************************************************
*                                                                             *
*                               Williams P+1                                  *
*                                                                             *
******************************************************************************/

/* Input: p is the initial generator (x0), if 0 generate it at random.
          n is the number to factor (assumed to be odd)
	  B1 is the stage 1 bound
	  B2 is the stage 2 bound
          k is the number of blocks for stage 2
          verbose is the verbosity level
   Output: f is the factor found, p is the residue at end of stage 1
   Return value: non-zero iff a factor is found (1 for stage 1, 2 for stage 2)
*/
int
pp1 (mpz_t f, mpz_t p, mpz_t n, mpz_t go, double *B1done, double B1,
     mpz_t B2min_parm, mpz_t B2_parm, unsigned long k,
     int verbose, int repr, int use_ntt, FILE *os, FILE *es,
     char *chkfilename, char *TreeFilename, double maxmem,
     gmp_randstate_t rng, int (*stop_asap)(void))
{
  int youpi = ECM_NO_FACTOR_FOUND;
  long st;
  mpres_t a;
  mpmod_t modulus;
  mpz_t B2min, B2; /* Local B2, B2min to avoid changing caller's values */
  faststage2_param_t faststage2_params;
  int twopass = 0;
  mpz_t p0;

  ASSERT (mpz_divisible_ui_p (n, 2) == 0);

  set_verbose (verbose);
  ECM_STDOUT = (os == NULL) ? stdout : os;
  ECM_STDERR = (es == NULL) ? stdout : es;

  st = cputime ();

  if (mpz_cmp_ui (p, 0) == 0)
    pp1_random_seed (p, n, rng);

  mpz_init_set (B2min, B2min_parm);
  mpz_init_set (B2, B2_parm);

  /* Set default B2. See ecm.c for comments */
  if (ECM_IS_DEFAULT_B2(B2))
    mpz_set_d (B2, pow (B1 * PP1FS2_COST, PM1FS2_DEFAULT_B2_EXPONENT));

  /* set B2min */
  if (mpz_sgn (B2min) < 0)
    mpz_set_d (B2min, B1);

  {
      long P;
      const unsigned long lmax = 1UL<<28; /* An upper bound */
      unsigned long lmax_NTT, lmax_noNTT;
      
      mpz_init (faststage2_params.m_1);
      faststage2_params.l = 0;
      faststage2_params.file_stem = TreeFilename;
      
      /* Find out what the longest transform length is we can do at all.
	 If no maxmem is given, the non-NTT can theoretically do any length. */

      lmax_NTT = 0;
      if (use_ntt)
	{
	  unsigned long t, t2 = 0;
	  /* See what transform length that the NTT can handle (due to limited 
	     primes and limited memory) */
	  t = mpzspm_max_len (n);
	  lmax_NTT = MIN (lmax, t);
	  if (maxmem != 0.)
	    {
	      t = pp1fs2_maxlen (double_to_size (maxmem), n, use_ntt, 0);
	      t = MIN (t, lmax_NTT);
	      /* Maybe the two pass variant lets us use a longer transform */
	      t2 = pp1fs2_maxlen (double_to_size (maxmem), n, use_ntt, 1);
	      t2 = MIN (t2, lmax_NTT);
	      if (t2 > t)
		{
		  t = t2;
		  twopass = 1;
		}
	      lmax_NTT = t;
	    }
	  outputf (OUTPUT_DEVVERBOSE, "NTT can handle lmax <= %lu\n", lmax_NTT);
	}

      /* See what transform length that the non-NTT code can handle */
      lmax_noNTT = lmax;
      if (maxmem != 0.)
	{
	  unsigned long t;
	  t = pp1fs2_maxlen (double_to_size (maxmem), n, 0, 0);
	  lmax_noNTT = MIN (lmax_noNTT, t);
	  outputf (OUTPUT_DEVVERBOSE, "non-NTT can handle lmax <= %lu\n", 
		   lmax_noNTT);
	}

      P = choose_P (B2min, B2, MAX(lmax_noNTT, lmax_NTT), k, 
		    &faststage2_params, B2min, B2, use_ntt, ECM_PP1);
      if (P == ECM_ERROR)
	{
          outputf (OUTPUT_ERROR, 
                   "Error: cannot choose suitable P value for your stage 2 "
                   "parameters.\nTry a shorter B2min,B2 interval.\n");
	  mpz_clear (faststage2_params.m_1);
	  return ECM_ERROR;
	}

      /* See if the selected parameters let us use NTT or not */
      if (faststage2_params.l > lmax_NTT)
	use_ntt = 0;
      
      if (maxmem != 0.)
	{
	  unsigned long MB;
	  char *s;
	  if (!use_ntt)
	    s = "out";
	  else if (twopass)
	    s = " two pass";
	  else
	    s = " one pass";

	  MB = pp1fs2_memory_use (faststage2_params.l, n, use_ntt, twopass)
	    / 1048576;
	  outputf (OUTPUT_VERBOSE, "Using lmax = %lu with%s NTT which takes "
		   "about %luMB of memory\n", faststage2_params.l, s, MB);
	}
  }

  /* Print B1, B2, polynomial and x0 */
  print_B1_B2_poly (OUTPUT_NORMAL, ECM_PP1, B1, *B1done, B2min_parm, B2min, 
		    B2, 1, p, 0, 0, NULL, 0, 0);

  /* If we do a stage 2, print its parameters */
  if (mpz_cmp (B2, B2min) >= 0)
    {
      /* can't mix 64-bit types and mpz_t on win32 for some reason */
      outputf (OUTPUT_VERBOSE, "P = %" PRIu64 ", l = %" PRIu64 ", "
            "s_1 = %" PRIu64 ", k = s_2 = %" PRIu64, 
             faststage2_params.P, faststage2_params.l,
             faststage2_params.s_1,faststage2_params.s_2);
      outputf (OUTPUT_VERBOSE, ", m_1 = %Zd\n", 
            faststage2_params.m_1);
    }

  if (test_verbose (OUTPUT_VERBOSE))
    {
      if (mpz_sgn (B2min_parm) >= 0)
        {
          outputf (OUTPUT_VERBOSE, 
            "Can't compute success probabilities for B1 <> B2min\n");
        }
      else
        {
          rhoinit (256, 10);
          /* If x0 is chosen randomly, the resulting group order will behave,
             on average, like for P-1, thus we use the same code as for P-1. */
          print_prob (B1, B2, 0, k, 1, go);
        }
    }

  mpmod_init (modulus, n, repr);
  mpres_init (a, modulus);
  mpres_set_z (a, p, modulus);

  /* since pp1_mul_prac takes an ecm_uint, we have to check
     that B1 <= ECM_UINT_MAX */
  if (B1 > (double) ECM_UINT_MAX)
    {
      outputf (OUTPUT_ERROR, "Error, maximal step1 bound for P+1 is %lu\n", 
               ECM_UINT_MAX);
      youpi = ECM_ERROR;
      goto clear_and_exit;
    }

  if (B1 > *B1done || mpz_cmp_ui (go, 1) > 0)
    youpi = pp1_stage1 (f, a, modulus, B1, B1done, go, stop_asap, 
                        chkfilename);

  outputf (OUTPUT_NORMAL, "Step 1 took %ldms\n", elltime (st, cputime ()));
  if (test_verbose (OUTPUT_RESVERBOSE))
    {
      mpz_t t;
      
      mpz_init (t);
      mpres_get_z (t, a, modulus);
      outputf (OUTPUT_RESVERBOSE, "x=%Zd\n", t);
      mpz_clear (t);
    }

  mpz_init_set (p0, p);

  /* store in p the residue at end of stage 1 */
  mpres_get_z (p, a, modulus);

  if (stop_asap != NULL && (*stop_asap) ())
    goto clear_and_exit;
      
  if (youpi == ECM_NO_FACTOR_FOUND && mpz_cmp (B2, B2min) >= 0)
    {
      if (use_ntt)
        youpi = pp1fs2_ntt (f, a, modulus, &faststage2_params, twopass);
      else 
        youpi = pp1fs2 (f, a, modulus, &faststage2_params);
    }

  if (youpi > 0 && test_verbose (OUTPUT_NORMAL))
    pp1_check_factor (p0, f); /* tell user if factor was found by P-1 */

  mpz_clear (p0);

 clear_and_exit:
  mpres_clear (a, modulus);
  mpmod_clear (modulus);
  mpz_clear (faststage2_params.m_1);
  mpz_clear (B2);
  mpz_clear (B2min);

  return youpi;
}
