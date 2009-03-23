/* The 'P+1' algorithm.

  Copyright 2002, 2003, 2004, 2005 Paul Zimmermann and Alexander Kruppa.

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

#include <math.h>
#include <stdlib.h>
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
*/
static void
pp1_mul (mpres_t P1, mpres_t P0, mpz_t e, mpmod_t n, mpres_t P, mpres_t Q)
{
  mp_size_t size_e;
  unsigned long i;
  int sign;

  sign = mpz_sgn (e);
  mpz_abs (e, e);
  if (sign == 0)
    {
      mpres_set_ui (P1, 2, n);
      goto unnegate;
    }

  if (mpz_cmp_ui (e, 1) == 0)
    {
      mpres_set (P1, P0, n);
      goto unnegate;
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
  mpz_add_ui (e, e, 1); /* recover original value of e */
unnegate:
  if (sign == -1)
    mpz_neg (e, e);
  
  return;
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
  for (p = 2.0; p <= B0; p = getprime ())
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
              outputf (OUTPUT_NORMAL, "Interrupted at prime %.0f\n", p);
	      if (p > *B1done)
		  *B1done = p;
              goto clear_and_exit;
            }
	}
    }

  pp1_mul (P0, P0, g, n, P, Q);

#if 1
  /* All primes sqrt(B1) < p <= B1 appear in exponent 1. All primes <= B1done
     are already included in exponent of at least 1, so it's save to skip 
     ahead to B1done+1 */
  
  if (*B1done > p)
    {
      getprime_seek ((*B1done) + 1.);
      p = getprime ();
    }
#endif

  /* then all primes > sqrt(B1) and taken with exponent 1 */
  for (; p <= B1; p = getprime ())
    {
      pp1_mul_prac (P0, (unsigned long) p, n, P, Q, R, S, T);
  
      if (stop_asap != NULL && (*stop_asap) ())
        {
          outputf (OUTPUT_NORMAL, "Interrupted at prime %.0f\n", p);
	  if (p > *B1done)
	      *B1done = p;
          goto clear_and_exit;
        }
      if (chkfilename != NULL && p > last_chkpnt_p + 10000. &&
          elltime (last_chkpnt_time, cputime ()) > CHKPNT_PERIOD)
        {
          writechkfile (chkfilename, ECM_PP1, p, n, NULL, P0, NULL);
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
#ifndef FULL_REDUCTION
  mpres_normalize (P); /* needed for gcd */
#endif
  mpres_gcd (f, P, n);
  youpi = mpz_cmp_ui (f, 1);

clear_and_exit:
  if (chkfilename != NULL)
    writechkfile (chkfilename, ECM_PP1, p, n, NULL, P0, NULL);
  getprime_clear (); /* free the prime tables, and reinitialize */
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
*                                  Stage 2                                    *
*                                                                             *
******************************************************************************/

/* let alpha, beta be the roots of x^2-Px+1=0
   set a, b such that alpha^e = a*alpha+b (idem for beta),
   i.e. a*x+b = rem(x^e, x^2-Px+1).
   Since (x-alpha)*(x-beta) = x^2-Px+1, we have alpha*beta = 1
   and alpha+beta = P, i.e. 1/alpha = beta = -alpha + P.
   It seems that if x^e % (x^2-Px+1) = a*x+b, then
   x^{-e+1} % (x^2-Px+1) = b*x+a. Proof?
*/

static void
pp1_mul2 (mpres_t a, mpres_t b, mpres_t P, mpz_t e, mpmod_t n)
{
  unsigned long l;
  mpres_t t;
  mpz_t abs_e;
  const int positive_e = (mpz_sgn (e) > 0);

  if (mpz_cmp_ui (e, 0UL) == 0) /* x^0 = 1 */
    {
      mpres_set_ui (a, 0, n);
      mpres_set_ui (b, 1, n);
      return;
    }
  
  mpres_init (t, n);
  mpz_init (abs_e);
  mpz_abs (abs_e, e);

  if (positive_e)
    {
      mpres_set_ui (a, 1, n);
      mpres_set_ui (b, 0, n);
    }
  else
    {
      /* Set to -x+P */
      mpres_set_ui (a, 1, n);
      mpres_neg (a, a, n);
      mpres_set (b, P, n);
    }

  l = mpz_sizeinbase (abs_e, 2) - 1; /* number of bits of e (minus 1) */
  
  while (l--)
    {
      /* square: (ax+b)^2 = (a^2P+2ab) x + (b^2-a^2) */
      mpres_mul (t, a, a, n); /* a^2 */
      mpres_mul (a, a, b, n);
      mpres_mul_ui (a, a, 2, n); /* 2ab */
      mpres_mul (b, b, b, n); /* b^2 */
      mpres_sub (b, b, t, n); /* b^2-a^2 */
      mpres_mul (t, t, P, n); /* a^2P */
      mpres_add (a, t, a, n); /* a^2P+2ab */
      
      if (mpz_tstbit (abs_e, l))
	{
	  if (positive_e)
	    {
	      /* multiply: (ax+b)*x = (aP+b) x - a */
	      mpres_mul (t, a, P, n);
	      mpres_add (t, t, b, n);
	      mpres_neg (b, a, n);
	      mpres_set (a, t, n);
	    }
	  else
	    {
	      /* multiply: (ax+b)*(-x+P) = 
                           -ax^2+(aP-b)x+b*P == 
			   -bx + (bP + a) (mod x^2-P*x+1) */
	      mpres_mul (t, b, P, n);
	      mpres_add (t, t, a, n);
	      mpres_neg (a, b, n);
	      mpres_set (b, t, n);
	    }
	}
    }

  mpz_clear (abs_e);
  mpres_clear (t, n);
}

/* 
   Performs the following:
   for (i=0;i<m;i++)
      for (j=0;j<n;j++)
        (x[j+(n+1)*i],y[j+(n+1)*i]) += (x[j+1+(n+1)*i],y[j+1+(n+1)*i])
 */
static void
addWnm (point *X, mpres_t P, mpmod_t modulus, unsigned int m,
        unsigned int n, unsigned long *tot_muls)
{
  unsigned long i, j, k;
  mpres_t t, u;

  mpres_init (t, modulus);
  mpres_init (u, modulus);
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      { /* (a*x+b)*(c*x+d) = (Pac+ad+bc) x + (bd-ac) */
	k = (n + 1) * i + j; /* n is the polynomial degree, so each table of
				differences contains n+1 entries */
	mpres_add (t, X[k].x, X[k].y, modulus); /* a+b */
	mpres_add (u, X[k+1].x, X[k+1].y, modulus); /* c+d */
	mpres_mul (t, t, u, modulus); /* ad+bc+ac+bd */
	mpres_mul (X[k].y, X[k].y, X[k+1].y, modulus); /* bd */
	mpres_sub (t, t, X[k].y, modulus); /* ad+bc+ac */
	mpres_mul (u, X[k].x, X[k+1].x, modulus); /* ac */
	mpres_sub (X[k].y, X[k].y, u, modulus); /* bd-ac */
	mpres_sub (t, t, u, modulus); /* ad+bc */
	mpres_mul (u, u, P, modulus); /* Pac */
	mpres_add (X[k].x, t, u, modulus); /* ad+bc+Pac */
	*tot_muls += 4;
      }
  mpres_clear (t, modulus);
  mpres_clear (u, modulus);
}

/* puts in F[0..dF-1] the successive values of 
   V_f(j*d2)(P) for Williams P+1,
   where f(x) = x^S for positive S, Dickson_{S,a}(x) for negative S,
   and 0 < j = 1 mod 6 < d1, j and d1 coprime.

   Return non-zero iff a factor was found (always zero in fact).
*/
int
pp1_rootsF (listz_t F, root_params_t *root_params, unsigned long dF, 
            mpres_t *x, listz_t t, mpmod_t modulus)
{
  unsigned long i, j;
  unsigned long muls = 0;
  long st, st1;
  int youpi = ECM_NO_FACTOR_FOUND;
  mpres_t fd[3];
  mpres_t u, v; /* auxiliary variables */
  listz_t coeffs;
  ecm_roots_state_t state;
  progression_params_t *params = &state.params; /* for less typing */

  if (dF == 0)
    return youpi;

  st1 = st = cputime ();

  outputf (OUTPUT_DEVVERBOSE, 
	   "pp1_rootsF: d1 = %lu, d2 = %lu, S = %d, dF = %lu\n",
	   root_params->d1, root_params->d2, root_params->S, dF);

  mpres_init (u, modulus);
  mpres_init (v, modulus);

  if (ABS(root_params->S) == 1) /* special code with d1/6 muls */
    {
      mpres_init (fd[0], modulus);
      mpres_init (fd[1], modulus);
      mpres_init (fd[2], modulus);
  
      mpz_set_ui (*t, root_params->d2);
      pp1_mul (fd[2], *x, *t, modulus, u, v);
      mpres_get_z (F[0], fd[2], modulus);
  
      mpz_set_ui (*t, 7UL);
      pp1_mul (fd[0], fd[2], *t, modulus, u, v);

      mpz_set_ui (*t, 6UL);
      pp1_mul (fd[1], fd[2], *t, modulus, u, v);

      /* fd[0] = V_{7*d2}(P), fd[1] = V_{6*d2}(P), fd[2] = V_{d2}(P) */

      outputf (OUTPUT_VERBOSE,
	       "Initializing table of differences for F took %ldms\n",
	       elltime (st1, cputime ()));

      i = 1;
      j = 7;
      while (i < dF)
	{
	  if (gcd (j, root_params->d1) == 1) /* (d2,d1) == 1 ==> (j*d2,d1) == (j,d1) */
	    mpres_get_z (F[i++], fd[0], modulus);

	  /* V_{m+n} = V_m * V_n - V_{m-n} */
	  /* fd[0] = V_m, fd[1] = V_n, fd[2] = V_{m-n} */
	  mpres_swap (fd[0], fd[2], modulus);
	  /* fd[0] = V_{m-n}, fd[1] = V_n, fd[2] = V_m */
	  mpres_mul (u, fd[2], fd[1], modulus);
	  /* u = V_n * V_m */
	  mpres_sub (fd[0], u, fd[0], modulus);
	  /* fd[0] = V_n * V_m - V_{m-n} = V_{m+n}, hence */
	  /* fd[0] = V_{m+n}, fd[1] = V_n, fd[2] = V_m */
	  j += 6;
	  muls ++;
	}

      mpres_clear (fd[0], modulus);
      mpres_clear (fd[1], modulus);
      mpres_clear (fd[2], modulus);
    }
  else /* case |S| <> 1: this code works also for S=1, but is more
	  expensive, since it can use up to 4*(d1/6) muls */
    {
      init_roots_params (params, root_params->S, root_params->d1, 
                         root_params->d2, 1.0);
      mpz_set_ui (*t, 0UL);
      coeffs = init_progression_coeffs (*t, params->dsieve, root_params->d2, 1, 
                                        6, params->S, params->dickson_a);
      
      if (coeffs == NULL)
        return ECM_ERROR;

      state.fd = (point *) malloc (params->size_fd * sizeof (point));
      if (state.fd == NULL)
        {
	  clear_list (coeffs, params->size_fd);
	  return ECM_ERROR;
        }
      for (i = 0; i < params->size_fd; i++)
	{
          mpres_init (state.fd[i].x, modulus);
          mpres_init (state.fd[i].y, modulus);
	  /* if i = k*(S+1) + S for k>=1, we can copy x and y from i - (S+1) */
	  if (i > params->S && (i % (params->S + 1) == params->S))
	    {
	      mpres_set (state.fd[i].x, state.fd[params->S].x, modulus);
	      mpres_set (state.fd[i].y, state.fd[params->S].y, modulus);
	    }
	  else
	    pp1_mul2 (state.fd[i].x, state.fd[i].y, x[0], coeffs[i], modulus);
	}
      clear_list (coeffs, params->size_fd);

      outputf (OUTPUT_VERBOSE,
	       "Initializing table of differences for F took %ldms\n",
	       elltime (st1, cputime ()));

      /* Now for the actual calculation of the roots. */
      for (i = 0; i < dF && !youpi;)
	{
	  /* Is this a rsieve value where we computed Dickson(j * d2) * X? */
	  if (gcd (params->rsieve, params->dsieve) == 1)
	    {
	      /* Did we use every progression since the last update? */
	      if (params->next == params->nr)
		{
		  /* Yes, time to update again */
		  addWnm (state.fd, x[0], modulus, params->nr, params->S, &muls);
		  params->next = 0;
		}

	      /* Is this a j value where we want Dickson(j*d2)*X as a root? */
	      if (gcd (params->rsieve, root_params->d1) == 1)
		{
		  /* we have alpha^k = x * alpha + y
		     thus alpha^k + beta^k = x * P + 2 * y.
                     FIXME: can we avoid returning to the Lucas form?
                  */
		  mpres_mul (u, state.fd[params->next * (params->S + 1)].x, x[0],
			     modulus);
		  mpres_mul_ui (v, state.fd[params->next * (params->S + 1)].y,
				2, modulus);
		  mpres_add (u, u, v, modulus);
		  mpres_get_z (F[i++], u, modulus);
		}

	      params->next ++;
	    }
	  params->rsieve += 6;
	}

      for (i = 0; i < params->size_fd; i++)
        {
          mpres_clear (state.fd[i].x, modulus);
          mpres_clear (state.fd[i].y, modulus);
        }
      free (state.fd);
    }

  mpres_clear (u, modulus);
  mpres_clear (v, modulus);

  outputf (OUTPUT_VERBOSE, "Computing roots of F took %ldms",
	   elltime (st, cputime ()));
  outputf (OUTPUT_DEVVERBOSE, " and %d muls", muls);
  outputf (OUTPUT_VERBOSE, "\n");
  
  return youpi;
}

/* return NULL if an error occurred */
pp1_roots_state_t *
pp1_rootsG_init (mpres_t *x, root_params_t *root_params, mpmod_t modulus)
{
  mpres_t P;
  pp1_roots_state_t *state;
  progression_params_t *params; /* for less typing */
  unsigned long i;

  ASSERT (gcd (root_params->d1, root_params->d2) == 1);

  state = (pp1_roots_state_t *) malloc (sizeof (pp1_roots_state_t));
  if (state == NULL)
    return NULL;
  params = &(state->params);

  /* we don't need the sign anymore after pp1_rootsG_init */
  params->S = ABS(root_params->S); 

  if (params->S == 1)
    {
      mpz_t t;
      mpz_init (t);
      mpres_init (P, modulus);
      for (i = 0; i < 4; i++)
        mpres_init (state->tmp[i], modulus);

      params->dsieve = root_params->d2; /* needed in pp1_rootsG */
      /* We want to skip values where gcd((i0 + i) * d1, d2) != 1.
	 We can test for gcd(i0 + i, d2) instead and let pp1_rootsG()
	 advance params->rsieve in steps of 1 */
      /* params->rsieve = i0 % d2 */
      params->rsieve = mpz_fdiv_ui (root_params->i0, root_params->d2);
      
      outputf (OUTPUT_DEVVERBOSE, "pp1_rootsG_init: i0 = %Zd, state: "
               "dsieve = %d, rsieve = %d, S = %d\n", root_params->i0, 
               params->dsieve, params->rsieve, params->S);
      
      mpz_set_ui (t, root_params->d1);
      pp1_mul (state->tmp[1], *x, t, modulus, state->tmp[3], P);
      pp1_mul (state->tmp[0], state->tmp[1], root_params->i0, modulus, 
               state->tmp[3], P);
      mpz_sub_ui (t, root_params->i0, 1);
      mpz_abs (t, t);
      pp1_mul (state->tmp[2], state->tmp[1], t, modulus, state->tmp[3], P);
      /* for P+1, tmp[0] = V_s(P), tmp[1] = V_d1(P), tmp[2] = V_{|s-d1|}(P) */

      mpres_clear (P, modulus);
      mpz_clear (t);
    }
  else
    {
      listz_t coeffs;
      
      params->dickson_a = (root_params->S < 0) ? -1 : 0;
      params->nr = (root_params->d2 > 1) ? root_params->d2 - 1 : 1;
      params->size_fd = params->nr * (params->S + 1);
      params->next = 0;
      params->dsieve = 1;
      params->rsieve = 1;

      state->fd = (point *) malloc (params->size_fd * sizeof (point));
      if (state->fd == NULL)
        {
          free (state);
          return NULL;
        }
      
      coeffs = init_progression_coeffs (root_params->i0, root_params->d2, 
                 root_params->d1, 1, 1, params->S, params->dickson_a);
      if (coeffs == NULL)
        {
          free (state->fd);
          free (state);
          return NULL;
        }

      for (i = 0; i < params->size_fd; i++) 
        {
          mpres_init (state->fd[i].x, modulus);
          mpres_init (state->fd[i].y, modulus);
          /* The S-th coeff of all progressions is identical */
          if (i > params->S && i % (params->S + 1) == params->S) 
            {
              /* Simply copy from the first progression */
              mpres_set (state->fd[i].x, state->fd[params->S].x, modulus); 
              mpres_set (state->fd[i].y, state->fd[params->S].y, modulus); 
            }
          else
            pp1_mul2 (state->fd[i].x, state->fd[i].y, x[0], coeffs[i], modulus);
        }

      clear_list (coeffs, params->size_fd);
    }

  return state;
}

void 
pp1_rootsG_clear (pp1_roots_state_t *state, ATTRIBUTE_UNUSED mpmod_t modulus)
{
  unsigned long i;

  if (state->params.S == 1)
    {
      for (i = 0; i < 4; i++)
        mpres_clear (state->tmp[i], modulus);
    }
  else
    {
      for (i = 0; i < state->params.size_fd; i++)
        {
          mpres_clear (state->fd[i].x, modulus);
          mpres_clear (state->fd[i].y, modulus);
        }
      free (state->fd);
    }

  free (state);
}

int
pp1_rootsG (listz_t G, unsigned long dF, pp1_roots_state_t *state, 
	    mpmod_t modulus, mpres_t *x)
{
  unsigned long i;
  unsigned long muls = 0;
  long st;
  progression_params_t *params = &(state->params); /* for less typing */

  st = cputime ();

  /* params->S is positive: we don't need the sign anymore, since the
     polynomial is defined by the table of differences */

  if (params->S == 1)
    {
      for (i = 0; i < dF;)
        {
          if (gcd (params->rsieve, params->dsieve) == 1)
            {
              outputf (OUTPUT_TRACE, "pp1_rootsG: Taking root G[%d], rsieve = %d\n", 
                       i, params->rsieve);
              mpres_get_z (G[i++], state->tmp[0], modulus);
            }
          else
            {
              outputf (OUTPUT_TRACE, "pp1_rootsG: NOT taking root, rsieve = %d, gcd = %d\n", 
                       params->rsieve, gcd (params->rsieve, params->dsieve));
            }

          mpres_swap (state->tmp[0], state->tmp[2], modulus);
          mpres_mul (state->tmp[3], state->tmp[2], state->tmp[1], modulus);
          mpres_sub (state->tmp[0], state->tmp[3], state->tmp[0], modulus);
          params->rsieve++;
        }
    }
  else
    {
      mpres_t u, v;

      mpres_init (u, modulus);
      mpres_init (v, modulus);
      for (i = 0; i < dF;)
        {
          /* Did we use every progression since the last update? */
          if (params->next == params->nr)
            {
              /* Yes, time to update again */
              addWnm (state->fd, x[0], modulus, params->nr, params->S, &muls);
              params->next = 0;
            }
      
          /* Is this a root we should skip? (Take only if gcd == 1) */
          if (gcd (params->rsieve, params->dsieve) == 1)
            {
              mpres_mul (u, state->fd[params->next * (params->S + 1)].x, x[0],
                         modulus);
              mpres_mul_ui (v, state->fd[params->next * (params->S + 1)].y,
                            2, modulus);
              mpres_add (u, u, v, modulus);
              mpres_get_z (G[i++], u, modulus);
            }
      
          params->next ++;
          params->rsieve ++;
        }
      mpres_clear (u, modulus);
      mpres_clear (v, modulus);
    }

  outputf (OUTPUT_VERBOSE, "Computing roots of G took %ldms",
	   elltime (st, cputime ()));
  outputf (OUTPUT_DEVVERBOSE, ", %lu muls", dF);
  outputf (OUTPUT_VERBOSE, "\n");
  
  return ECM_NO_FACTOR_FOUND;
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
          verbose is the verbosity level
   Output: p is the factor found
   Return value: non-zero iff a factor is found (1 for stage 1, 2 for stage 2)
*/
int
pp1 (mpz_t f, mpz_t p, mpz_t n, mpz_t go, double *B1done, double B1,
     mpz_t B2min_parm, mpz_t B2_parm, double B2scale, unsigned long k, 
     const int S, int verbose, int repr, int use_ntt, FILE *os, FILE *es, 
     char *chkfilename, char *TreeFilename, double maxmem, 
     gmp_randstate_t rng, int (*stop_asap)(void))
{
  int youpi = ECM_NO_FACTOR_FOUND;
  int po2 = 0;    /* Whether we should use power-of-2 poly degree */
  long st;
  mpres_t a;
  mpmod_t modulus;
  mpz_t B2min, B2; /* Local B2, B2min to avoid changing caller's values */
  unsigned long dF;
  root_params_t root_params;
  faststage2_param_t faststage2_params;
  const int stage2_variant = (S == 1 || S == ECM_DEFAULT_S);
  int twopass = 0;

  set_verbose (verbose);
  ECM_STDOUT = (os == NULL) ? stdout : os;
  ECM_STDERR = (es == NULL) ? stdout : es;

  /* if n is even, return 2 */
  if (mpz_divisible_2exp_p (n, 1))
    {
      mpz_set_ui (f, 2);
      return ECM_FACTOR_FOUND_STEP1;
    }

  st = cputime ();

  if (mpz_cmp_ui (p, 0) == 0)
    pm1_random_seed (p, n, rng);

  mpz_init_set (B2min, B2min_parm);
  mpz_init_set (B2, B2_parm);

  /* Set default B2. See ecm.c for comments */
  if (ECM_IS_DEFAULT_B2(B2))
    {
      if (stage2_variant == 0)
        mpz_set_d (B2, B2scale * pow (B1 * PP1_COST, DEFAULT_B2_EXPONENT));
      else
        mpz_set_d (B2, B2scale * pow (B1 * PP1FS2_COST, 
                   PM1FS2_DEFAULT_B2_EXPONENT));
    }

  /* set B2min */
  if (mpz_sgn (B2min) < 0)
    mpz_set_d (B2min, B1);

  mpmod_init (modulus, n, repr);

  if (use_ntt)
    po2 = 1;
  
  if (stage2_variant != 0)
    {
      long P;
      const unsigned long lmax = 1UL<<28; /* An upper bound */
      unsigned long lmax_NTT, lmax_noNTT;
      
      mpz_init (faststage2_params.m_1);
      
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
  else 
    {
      mpz_init (root_params.i0);
      root_params.d2 = 0; /* Enable automatic choice of d2 */
      if (bestD (&root_params, &k, &dF, B2min, B2, po2, use_ntt, maxmem,
		 (TreeFilename != NULL), modulus) == ECM_ERROR)
	{
	  youpi = ECM_ERROR;
	  goto clear_and_exit;
	}
      
      /* Set default degree for Brent-Suyama extension */
      root_params.S = S;
      if (root_params.S == ECM_DEFAULT_S)
	{
	  if (modulus->repr == ECM_MOD_BASE2 && modulus->Fermat > 0)
	    {
	      /* For Fermat numbers, default is 1 (no Brent-Suyama) */
	      root_params.S = 1;
	    }
	  else
	    {
	      mpz_t t;
	      mpz_init (t);
	      mpz_sub (t, B2, B2min);
	      root_params.S = choose_S (t);
	      mpz_clear (t);
	    }
	}
    }

  /* Print B1, B2, polynomial and x0 */
  print_B1_B2_poly (OUTPUT_NORMAL, ECM_PP1, B1, *B1done, B2min_parm, B2min, 
		    B2, (stage2_variant == 0) ? root_params.S : 1, p, 0);

  /* If we do a stage 2, print its parameters */
  if (mpz_cmp (B2, B2min) >= 0)
    {
      if (stage2_variant != 0)
        outputf (OUTPUT_VERBOSE, "P = %lu, l = %lu, s_1 = %lu, k = s_2 = %lu, "
                 "m_1 = %Zd\n", faststage2_params.P, faststage2_params.l,
                 faststage2_params.s_1,faststage2_params.s_2,
                 faststage2_params.m_1);
      else
        outputf (OUTPUT_VERBOSE, "dF=%lu, k=%lu, d=%lu, d2=%lu, i0=%Zd\n", 
                 dF, k, root_params.d1, root_params.d2, 
                 S == 1 ? faststage2_params.m_1 : root_params.i0);
    }

  mpres_init (a, modulus);
  mpres_set_z (a, p, modulus);

  /* since pp1_mul_prac takes an unsigned long, we have to check
     that B1 <= MAX_ULONG */
  if (B1 > (double) ULONG_MAX)
    {
      outputf (OUTPUT_ERROR, "Error, maximal step1 bound for P+1 is %lu\n", 
               ULONG_MAX);
      youpi = ECM_ERROR;
      goto clear_and_exit;
    }

  if (B1 > *B1done)
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

  mpres_get_z (p, a, modulus);

  if (stop_asap != NULL && (*stop_asap) ())
    goto clear_and_exit;
      
  if (youpi == ECM_NO_FACTOR_FOUND && mpz_cmp (B2, B2min) >= 0)
    {
      if (stage2_variant != 0)
        {
          if (use_ntt)
            youpi = pp1fs2_ntt (f, a, modulus, &faststage2_params, twopass);
          else 
            youpi = pp1fs2 (f, a, modulus, &faststage2_params);
        }
      else
	youpi = stage2 (f, &a, modulus, dF, k, &root_params, ECM_PP1, 
			use_ntt, TreeFilename, stop_asap);
    }

  if (youpi > 0 && test_verbose (OUTPUT_NORMAL))
    pp1_check_factor (p, f); /* tell user if factor was found by P-1 */

 clear_and_exit:
  mpres_clear (a, modulus);
  mpmod_clear (modulus);
  if (stage2_variant != 0)
    mpz_clear (faststage2_params.m_1);
  else
    mpz_clear (root_params.i0);
  mpz_clear (B2);
  mpz_clear (B2min);

  return youpi;
}
