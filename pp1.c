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
  the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
  MA 02111-1307, USA.

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
    }

  pp1_mul (P0, P0, g, n, P, Q);

  /* then all primes > sqrt(B1) and taken with exponent 1 */
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
pp1_check_factor (mpz_t a, mpz_t p, FILE *ECM_STDOUT)
{
  if (mpz_probab_prime_p (p, PROBAB_PRIME_TESTS))
    {
      mpz_mul (a, a, a);
      mpz_sub_ui (a, a, 4);
      if (mpz_jacobi (a, p) == 1)
        fprintf (ECM_STDOUT, "[factor found by P-1]\n");
    }
}

/******************************************************************************
*                                                                             *
*                                  Stage 2                                    *
*                                                                             *
******************************************************************************/

/* let alpha, beta be the roots of x^2-Px+1=0
   set a, b such that alpha^e = a*alpha+b (idem for beta),
   i.e. a*x+b = rem(x^e, x^2-Px+1) */
static void
pp1_mul2 (mpres_t a, mpres_t b, mpres_t P, mpz_t e, mpmod_t n)
{
  unsigned long l, i;
  mpres_t t;

  if (mpz_cmp_ui (e, 0) == 0) /* x^0 = 1 */
    {
      mpres_set_ui (a, 0, n);
      mpres_set_ui (b, 1, n);
      return;
    }

  /* now e >= 1 */
  mpres_set_ui (a, 1, n);
  mpres_set_ui (b, 0, n);

  l = mpz_sizeinbase (e, 2) - 1; /* number of bits of e (minus 1) */

  mpres_init (t, n);

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

      if (mpz_tstbit (e, l)) /* multiply: (ax+b)*x = (aP+b) x - a */
        {
          mpres_mul (t, a, P, n);
          mpres_add (t, t, b, n);
          mpres_neg (b, a, n);
          mpz_swap (a, t);
        }
    }

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
pp1_rootsF (listz_t F, unsigned int d1, unsigned int d2, unsigned int dF, 
            mpres_t *x, listz_t t, int S, mpmod_t modulus)
{
  unsigned int i, j;
  unsigned long muls = 0;
  int st, st2, youpi = ECM_NO_FACTOR_FOUND;
  mpres_t fd[3];
  mpres_t u, v; /* auxiliary variables */
  listz_t coeffs;
  ecm_roots_state state;

  if (dF == 0)
    return youpi;

  st2 = st = cputime ();

  outputf (OUTPUT_DEVVERBOSE, "pp1_rootsF: d1 = %d, d2 = %d, dF = %d\n",
	   d1, d2, dF);

  mpres_init (u, modulus);
  mpres_init (v, modulus);

  if (ABS(S) == 1) /* special code with d1/6 muls */
    {
      mpres_init (fd[0], modulus);
      mpres_init (fd[1], modulus);
      mpres_init (fd[2], modulus);
  
      mpz_set_ui (*t, d2);
      pp1_mul (fd[2], *x, *t, modulus, u, v);
      mpres_get_z (F[0], fd[2], modulus);
  
      mpz_set_ui (*t, 7);
      pp1_mul (fd[0], fd[2], *t, modulus, u, v);

      mpz_set_ui (*t, 6);
      pp1_mul (fd[1], fd[2], *t, modulus, u, v);

      /* fd[0] = V_{7*d2}(P), fd[1] = V_{6*d2}(P), fd[2] = V_{d2}(P) */

      outputf (OUTPUT_VERBOSE,
	       "Initializing table of differences for F took %dms\n",
	       cputime () - st2);

      i = 1;
      j = 7;
      while (i < dF)
	{
	  if (gcd (j, d1) == 1)
	    mpres_get_z (F[i++], fd[0], modulus);

	  /* V_{m+n} = V_m * V_n - V_{m-n} */
	  /* fd[0] = V_m, fd[1] = V_n, fd[2] = V_{m-n} */
	  mpres_swap (fd[0], fd[2], modulus);
	  mpres_mul (u, fd[2], fd[1], modulus);
	  mpres_sub (fd[0], u, fd[0], modulus);
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
      init_roots_state (&state, S, d1, d2, 1.0);
      coeffs = init_progression_coeffs (0.0, state.dsieve, d2, 1, 6, state.S,
                                        state.dickson_a);
      if (coeffs == NULL)
        return ECM_ERROR;

      state.fd = (point *) malloc (state.size_fd * sizeof (point));
      if (state.fd == NULL)
        {
	  clear_list (coeffs, state.size_fd);
	  return ECM_ERROR;
        }
      for (i = 0; i < state.size_fd; i++)
	{
          mpres_init (state.fd[i].x, modulus);
          mpres_init (state.fd[i].y, modulus);
	  /* if i = k*(S+1) + S for k>=1, we can copy x and y from i - (S+1) */
	  if (i > state.S && (i % (state.S + 1) == state.S))
	    {
	      mpres_set (state.fd[i].x, state.fd[state.S].x, modulus);
	      mpres_set (state.fd[i].y, state.fd[state.S].y, modulus);
	    }
	  else
	    pp1_mul2 (state.fd[i].x, state.fd[i].y, x[0], coeffs[i], modulus);
	}
      clear_list (coeffs, state.size_fd);

      outputf (OUTPUT_VERBOSE,
	       "Initializing table of differences for F took %dms\n",
	       cputime () - st2);

      /* Now for the actual calculation of the roots. */
      for (i = 0; i < dF && !youpi;)
	{
	  /* Is this a rsieve value where we computed Dickson(j * d2) * X? */
	  if (gcd (state.rsieve, state.dsieve) == 1)
	    {
	      /* Did we use every progression since the last update? */
	      if (state.next == state.nr)
		{
		  /* Yes, time to update again */
		  addWnm (state.fd, x[0], modulus, state.nr, state.S, &muls);
		  state.next = 0;
		}

	      /* Is this a j value where we want Dickson(j*d2)*X as a root? */
	      if (gcd (state.rsieve, d1) == 1)
		{
		  /* we have alpha^k = x * alpha + y
		     thus alpha^k + beta^k = x * P + 2 * y.
                     FIXME: can we avoid returning to the Lucas form?
                  */
		  mpres_mul (u, state.fd[state.next * (state.S + 1)].x, x[0],
			     modulus);
		  mpres_mul_ui (v, state.fd[state.next * (state.S + 1)].y,
				2, modulus);
		  mpres_add (u, u, v, modulus);
		  mpres_get_z (F[i++], u, modulus);
		}

	      state.next ++;
	    }
	  state.rsieve += 6;
	}

      
    clear_fdi:
      for (i = 0; i < state.size_fd; i++)
        {
          mpres_clear (state.fd[i].x, modulus);
          mpres_clear (state.fd[i].y, modulus);
        }
    }

  mpres_clear (u, modulus);
  mpres_clear (v, modulus);

  outputf (OUTPUT_VERBOSE, "Computing roots of F took %dms", cputime () - st);
  outputf (OUTPUT_DEVVERBOSE, " and %d muls", muls);
  outputf (OUTPUT_VERBOSE, "\n");
  
  return youpi;
}

/* return NULL if an error occurred */
pp1_roots_state *
pp1_rootsG_init (mpres_t *x, double s, unsigned int d1, unsigned int d2, 
                 int S, mpmod_t modulus)
{
  int st;
  mpres_t P;
  mpz_t t;
  pp1_roots_state *state;
  unsigned long i;

  ASSERT (gcd (d1, d2) == 1);

  st = cputime ();

  state = (pp1_roots_state *) malloc (sizeof (pp1_roots_state));
  if (state == NULL)
    return NULL;

  state->S = ABS(S); /* we don't need the sign anymore after pp1_rootsG_init */

  if (state->S == 1)
    {
      mpz_init (t);
      mpres_init (P, modulus);
      for (i = 0; i < 4; i++)
        mpres_init (state->tmp[i], modulus);

      state->d = d1; /* needed in pp1_rootsG */
      state->dsieve = d2; /* needed in pp1_rootsG */
      /* We want to skip values where gcd(s + i * d1, d2) != 1 */
      /* state->rsieve = s % d2 */
      state->rsieve = (unsigned int) (s - floor (s / (double) d2) * (double) d2);

      mpz_set_d (t, s);
      pp1_mul (state->tmp[0], *x, t, modulus, state->tmp[3], P);
      mpz_set_ui (t, d1);
      pp1_mul (state->tmp[1], *x, t, modulus, state->tmp[3], P);
      mpz_set_d (t, fabs (s - (double) d1));
      pp1_mul (state->tmp[2], *x, t, modulus, state->tmp[3], P);
      /* for P+1, tmp[0] = V_s(P), tmp[1] = V_d1(P), tmp[2] = V_{|s-d1|}(P) */

      mpres_clear (P, modulus);
      mpz_clear (t);
    }
  else
    {
      listz_t coeffs;
      int dickson_a = (S < 0) ? -1 : 0;

      state->nr = (d2 > 1) ? d2 - 1 : 1;
      state->size_fd = state->nr * (state->S + 1);
      state->next = 0;
      state->dsieve = 1;
      state->rsieve = 1;

      state->fd = (point *) malloc (state->size_fd * sizeof (point));
      if (state->fd == NULL)
        {
          free (state);
          return NULL;
        }

      coeffs = init_progression_coeffs (s, d2, d1, 1, 1, state->S, dickson_a);
      if (coeffs == NULL)
        {
          free (state->fd);
          free (state);
          return NULL;
        }

      for (i = 0; i < state->size_fd; i++) 
        {
          mpres_init (state->fd[i].x, modulus);
          mpres_init (state->fd[i].y, modulus);
          /* The S-th coeff of all progressions is identical */
          if (i > state->S && i % (state->S + 1) == state->S) 
            {
              /* Simply copy from the first progression */
              mpres_set (state->fd[i].x, state->fd[state->S].x, modulus); 
              mpres_set (state->fd[i].y, state->fd[state->S].y, modulus); 
            }
          else
            pp1_mul2 (state->fd[i].x, state->fd[i].y, x[0], coeffs[i], modulus);
        }

      clear_list (coeffs, state->size_fd);
    }

  return state;
}

void 
pp1_rootsG_clear (pp1_roots_state *state, ATTRIBUTE_UNUSED mpmod_t modulus)
{
  unsigned long i;

  if (state->S == 1)
    {
      for (i = 0; i < 4; i++)
        mpres_clear (state->tmp[i], modulus);
    }
  else
    {
      for (i = 0; i < state->size_fd; i++)
        {
          mpres_clear (state->fd[i].x, modulus);
          mpres_clear (state->fd[i].y, modulus);
        }
      free (state->fd);
    }

  free (state);
}

int
pp1_rootsG (listz_t G, unsigned int dF, pp1_roots_state *state, mpmod_t modulus,
            mpres_t *x)
{
  unsigned int i;
  unsigned long muls = 0;
  int st;

  st = cputime ();

  /* state->S is positive: we don't need the sign anymore, since the
     polynomial is defined by the table of differences */

  if (state->S == 1)
    {
      for (i = 0; i < dF;)
        {
          if (gcd (state->rsieve, state->dsieve) == 1)
            mpres_get_z (G[i++], state->tmp[0], modulus);

          mpres_swap (state->tmp[0], state->tmp[2], modulus);
          mpres_mul (state->tmp[3], state->tmp[2], state->tmp[1], modulus);
          mpres_sub (state->tmp[0], state->tmp[3], state->tmp[0], modulus);
          state->rsieve = (state->rsieve + state->d) % state->dsieve;
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
          if (state->next == state->nr)
            {
              /* Yes, time to update again */
              addWnm (state->fd, x[0], modulus, state->nr, state->S, &muls);
              state->next = 0;
            }
      
          /* Is this a root we should skip? (Take only if gcd == 1) */
          if (gcd (state->rsieve, state->dsieve) == 1)
            {
              mpres_mul (u, state->fd[state->next * (state->S + 1)].x, x[0],
                         modulus);
              mpres_mul_ui (v, state->fd[state->next * (state->S + 1)].y,
                            2, modulus);
              mpres_add (u, u, v, modulus);
              mpres_get_z (G[i++], u, modulus);
            }
      
          state->next ++;
          state->rsieve ++;
        }
      mpres_clear (u, modulus);
      mpres_clear (v, modulus);
    }

  outputf (OUTPUT_VERBOSE, "Computing roots of G took %dms", cputime () - st);
  outputf (OUTPUT_DEVVERBOSE, ", %u muls", dF);
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
pp1 (mpz_t f, mpz_t p, mpz_t n, mpz_t go, double B1done, double B1,
     double B2min, double B2, double B2scale, unsigned int k, int S,
     int verbose, int repr, FILE *os, FILE *es, char *TreeFilename)
{
  int youpi = 0, st;
  mpres_t a;
  mpmod_t modulus;

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
    {
      gmp_randstate_t state;
      gmp_randinit_default (state);
      pm1_random_seed (p, n, state);
      gmp_randclear (state);
    }

  /* Set default B2. See ecm.c for comments */
  if (ECM_IS_DEFAULT_B2(B2))
    B2 = pow (2.0 * B1 / 6.0, 1.424828748);

  /* Scale B2 by what the user said (or by the default scaling of 1.0) */
  B2 *= B2scale;

  /* set B2min */
  if (B2min < 0.0)
    B2min = B1;

  /* Set default degree for Brent-Suyama extension */
  if (S == ECM_DEFAULT_S)
    S = choose_S (B2 - B2min);

  if (test_verbose (OUTPUT_NORMAL))
    {
      fprintf (ECM_STDOUT, "Using ");
      if (ECM_IS_DEFAULT_B1_DONE(B1done))
        fprintf (ECM_STDOUT, "B1=%1.0f", B1);
      else
        fprintf (ECM_STDOUT, "B1=%1.0f-%1.0f", B1done, B1);
      if (B2min <= B1)
        fprintf (ECM_STDOUT, ", B2=%1.0f", B2);
      else
        fprintf (ECM_STDOUT, ", B2=%1.0f-%1.0f", B2min, B2);

      if (S > 0)
        fprintf (ECM_STDOUT, ", polynomial x^%u", S);
      else
        fprintf (ECM_STDOUT, ", polynomial Dickson(%u)", -S);
      if (ECM_IS_DEFAULT_B1_DONE(B1done) || test_verbose (OUTPUT_VERBOSE)) /* don't print x0 in resume case */
	{
	  fprintf (ECM_STDOUT, ", x0=");
	  mpz_out_str (ECM_STDOUT, 10, p);
	}
      fprintf (ECM_STDOUT, "\n");
      fflush (ECM_STDOUT);
    }

  if (repr == 1)
    mpmod_init_MPZ (modulus, n);
  else   if (repr == 2)
    mpmod_init_MODMULN (modulus, n);
  else if (repr == 3)
    mpmod_init_REDC (modulus, n);
  else if (repr > 16)
    mpmod_init_BASE2 (modulus, repr, n);
  else /* automatic choice */
    mpmod_init (modulus, n, repr);

  mpres_init (a, modulus);
  mpres_set_z (a, p, modulus);

  /* since pp1_mul_prac takes an unsigned long, we have to check
     that B1 <= MAX_ULONG */
  if (B1 > (double) ULONG_MAX)
    {
      fprintf (ECM_STDERR, "Error, maximal step1 bound for P+1 is %lu\n", ULONG_MAX);
      youpi = ECM_ERROR;
      goto clear_pp1;
    }

  if (B1 > B1done)
    youpi = pp1_stage1 (f, a, modulus, B1, B1done, go);

  st = cputime () - st;

  outputf (OUTPUT_NORMAL, "Step 1 took %dms\n", st);
  if (test_verbose (OUTPUT_VERBOSE))
    {
      fprintf (ECM_STDOUT, "x=");
      mpres_out_str (ECM_STDOUT, 10, a, modulus);
      fprintf (ECM_STDOUT, "\n");
      fflush (ECM_STDOUT);
    }

  if (youpi == ECM_NO_FACTOR_FOUND) /* no factor found, no error */
    youpi = stage2 (f, &a, modulus, B2min, B2, k, S, ECM_PP1, st, 
                    TreeFilename);

  if (youpi > 0 && test_verbose (OUTPUT_NORMAL))
    pp1_check_factor (p, f, ECM_STDOUT); /* tell user if factor was found by P-1 */

  mpres_get_z (p, a, modulus);

 clear_pp1:
  mpres_clear (a, modulus);
  mpmod_clear (modulus);

  return youpi;
}
