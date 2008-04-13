/* Pollard 'P-1' algorithm.

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
  the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
  MA 02110-1301, USA.
*/

#include <math.h>
#include <stdlib.h>
#include "ecm-impl.h"

#define CASCADE_THRES 3
#define CASCADE_MAX 50000000.0
#ifndef POWM_THRESHOLD
#define POWM_THRESHOLD 100
#endif

typedef struct {
  unsigned int size;
  mpz_t *val;
} mul_casc;

/******************************************************************************
*                                                                             *
*                                  Stage 1                                    *
*                                                                             *
******************************************************************************/

/* prime powers are accumulated up to about n^L1 */
#define L1 16

/*** Cascaded multiply ***/

/* return NULL if an error occurred */
static mul_casc *
mulcascade_init (void)
{
  mul_casc *t;

  t = (mul_casc *) malloc (sizeof (mul_casc));
  if (t == NULL)
    {
      outputf (OUTPUT_ERROR, "mulcascade_init: could not allocate memory\n");
      return NULL;
    }
  t->val = (mpz_t*) malloc (sizeof (mpz_t));
  if (t->val == NULL)
    {
      outputf (OUTPUT_ERROR, "mulcascade_init: could not allocate memory\n");
      free (t);
      return NULL;
    }
  mpz_init (t->val[0]);
  t->size = 1;
  return t;
}

static void 
mulcascade_free (mul_casc *c)
{
  unsigned int i;

  for (i = 0; i < c->size; i++)
    mpz_clear (c->val[i]);
  free (c->val);
  free (c);
}

static mul_casc * 
mulcascade_mul_d (mul_casc *c, const double n, ATTRIBUTE_UNUSED mpz_t t)
{
  unsigned int i;

  if (mpz_sgn (c->val[0]) == 0)
    {
      mpz_set_d (c->val[0], n);
      return c;
    }

  mpz_mul_d (c->val[0], c->val[0], n, t);
  if (mpz_size (c->val[0]) <= CASCADE_THRES)
    return c;
  
  for (i = 1; i < c->size; i++) 
    {
      if (mpz_sgn (c->val[i]) == 0) 
        {
          mpz_set (c->val[i], c->val[i-1]);
          mpz_set_ui (c->val[i-1], 0);
          return c;
        }
      else
	{
          mpz_mul (c->val[i], c->val[i], c->val[i-1]);
          mpz_set_ui (c->val[i-1], 0);
        }
    }
  
  /* Allocate more space for cascade */
  
  i = c->size++;
  c->val = (mpz_t*) realloc (c->val, c->size * sizeof (mpz_t));
  mpz_init (c->val[i]);
  mpz_swap (c->val[i], c->val[i-1]);

  return c;
}

static mul_casc * 
mulcascade_mul (mul_casc *c, mpz_t n)
{
  unsigned int i;

  if (mpz_sgn (c->val[0]) == 0)
    {
      mpz_set (c->val[0], n);
      return c;
    }

  mpz_mul (c->val[0], c->val[0], n);
  if (mpz_size (c->val[0]) <= CASCADE_THRES)
    return c;
  
  for (i = 1; i < c->size; i++) 
    {
      if (mpz_sgn (c->val[i]) == 0) 
        {
          mpz_set (c->val[i], c->val[i-1]);
          mpz_set_ui (c->val[i-1], 0);
          return c;
        } else {
          mpz_mul (c->val[i], c->val[i], c->val[i-1]);
          mpz_set_ui (c->val[i-1], 0);
        }
    }
  
  /* Allocate more space for cascade */
  
  i = c->size++;
  c->val = (mpz_t*) realloc (c->val, c->size * sizeof (mpz_t));
  mpz_init (c->val[i]);
  mpz_swap (c->val[i], c->val[i-1]);

  return c;
}

static void 
mulcascade_get_z (mpz_t r, mul_casc *c) 
{
  unsigned int i;
  
  if (c->size == 0)
    {
      mpz_set_ui (r, 1); /* Empty product */
      return;
    }

  mpz_set_ui (r, 1);
  
  for (i = 0; i < c->size; i++)
    if (mpz_sgn (c->val[i]) != 0)
      mpz_mul (r, r, c->val[i]);
}


/* Input:  a is the generator (sigma)
           n is the number to factor
           B1 is the stage 1 bound
	   B1done: stage 1 was already done up to that limit
	   go is the group order to preload
   Output: f is the factor found, a is the value at end of stage 1
	   B1done is set to B1 if stage 1 completed normally,
	   or to the largest prime processed if interrupted, but never
	   to a smaller value than B1done was upon function entry.
   Return value: non-zero iff a factor was found (or an error occurred).
*/

static int
pm1_stage1 (mpz_t f, mpres_t a, mpmod_t n, double B1, double *B1done, 
            mpz_t go, int (*stop_asap)(void), char *chkfilename)
{
  double p, q, r, cascade_limit, last_chkpnt_p;
  mpz_t g, d;
  int youpi = ECM_NO_FACTOR_FOUND;
  unsigned int size_n, max_size;
  unsigned int smallbase = 0;
  mul_casc *cascade;
  long last_chkpnt_time;
  const double B0 = sqrt (B1);

  mpz_init (g);
  mpz_init (d);

  size_n = mpz_sizeinbase (n->orig_modulus, 2);
  max_size = L1 * size_n;

  mpres_get_z (g, a, n);
  if (mpz_fits_uint_p (g))
    smallbase = mpz_get_ui (g);

  /* suggestion from Peter Montgomery: start with exponent n-1,
     since any prime divisor of b^m-1 which does not divide any
     algebraic factor of b^m-1 must be of the form km+1 [Williams82].
     Do this only when n is composite, otherwise all tests with prime
     n factor of a Cunningham number will succeed in stage 1.

     Since mpz_probab_prime_p and a^(n-1) mod n require about lg(n) modular
     multiplications, and P-1 perform about B1 modular multiplications,
     to ensure small overhead, use that trick only when lg(n) <= sqrt(B1).
  */
  /* For now, this p^N-1 is left in.  We might want it out at a later time */
  if ((double) size_n <= B0 &&
      mpz_probab_prime_p (n->orig_modulus, PROBAB_PRIME_TESTS) == 0)
    {
      mpz_sub_ui (g, n->orig_modulus, 1);
      mpres_pow (a, a, g, n);
    }
  else
    mpz_set_ui (g, 1);

  /* Set a limit of roughly 10000 * log_10(N) for the primes that are 
     multiplied up in the exponent, i.e. 1M for a 100 digit number, 
     but limit to CASCADE_MAX to avoid problems with stack allocation */
  
  cascade_limit = 3000.0 * (double) size_n;

  if (cascade_limit > CASCADE_MAX)
    cascade_limit = CASCADE_MAX;
  
  if (cascade_limit > B1)
    cascade_limit = B1;

  cascade = mulcascade_init ();
  if (cascade == NULL)
    {
      youpi = ECM_ERROR;
      goto clear_pm1_stage1;
    }

  /* since B0 = sqrt(B1), we can have B0 > cascade_limit only when
     B1 > cascade_limit^2. This cannot happen when cascade_limit=B1,
     thus we need B1 > min(CASCADE_MAX, 3000*sizeinbase(n,2))^2.
     For sizeinbase(n,2) <= CASCADE_MAX/3000 (less than 5017 digits 
     for CASCADE_MAX=5e7) this means B1 > 9e6*sizeinbase(n,2)^2.
     For sizeinbase(n,2) > CASCADE_MAX/3000, this means B1 > CASCADE_MAX^2,
     i.e. B1 > 25e14 for CASCADE_MAX=5e7.
*/

  /* if the user knows that P-1 has a given divisor, he can supply it */
  if (mpz_cmp_ui (go, 1) > 0)
    cascade = mulcascade_mul (cascade, go);

  last_chkpnt_time = cputime ();
  last_chkpnt_p = 2.;
  if (B0 <= cascade_limit)
    {
      /* first loop through small primes <= sqrt(B1) */
      for (p = 2.0; p <= B0; p = getprime ())
        {
          for (q = 1, r = p; r <= B1; r *= p)
            if (r > *B1done) q *= p;
          cascade = mulcascade_mul_d (cascade, q, d);
        }

      /* then all sqrt(B1) < primes < cascade_limit and taken with 
         exponent 1 */
      for ( ; p <= cascade_limit; p = getprime ())
        if (p > *B1done)
          cascade = mulcascade_mul_d (cascade, p, d);
   
      mulcascade_get_z (g, cascade);
      mulcascade_free (cascade);

      outputf (OUTPUT_DEVVERBOSE, "Exponent has %u bits\n", 
               mpz_sizeinbase (g, 2));

      if (smallbase)
        {
	  outputf (OUTPUT_VERBOSE, "Using mpres_ui_pow, base %u\n", smallbase);
          mpres_ui_pow (a, smallbase, g, n);
        }
      else
	{
	  mpres_pow (a, a, g, n);
	}
      mpz_set_ui (g, 1);
    }
  else
    {
      for (p = 2.0; p <= cascade_limit; p = getprime ())
        {
          for (q = 1.0, r = p; r <= B1; r *= p)
            if (r > *B1done) q *= p;
          cascade = mulcascade_mul_d (cascade, q, d);
        }
      
      mulcascade_get_z (g, cascade);
      mulcascade_free (cascade);

      outputf (OUTPUT_DEVVERBOSE, "Exponent has %u bits\n", 
               mpz_sizeinbase (g, 2));

      if (smallbase)
        {
	  outputf (OUTPUT_VERBOSE, "Using mpres_ui_pow, base %u\n", smallbase);
          mpres_ui_pow (a, smallbase, g, n);
        }
      else
        {
          mpres_pow (a, a, g, n);
        }
      mpz_set_ui (g, 1);
      
      for ( ; p <= B0; p = getprime ())
        {
          for (q = 1, r = p; r <= B1; r *= p)
            if (r > *B1done) q *= p;
          mpz_mul_d (g, g, q, d);
          if (mpz_sizeinbase (g, 2) >= max_size)
            {
              mpres_pow (a, a, g, n);
              mpz_set_ui (g, 1);
            if (stop_asap != NULL && (*stop_asap) ())
              {
                outputf (OUTPUT_NORMAL, "Interrupted at prime %.0f\n", p);
                if (p > *B1done)
                  *B1done = p;
                goto clear_pm1_stage1;
              }
            }
        }
    }
  
  /* All primes sqrt(B1) < p <= B1 appear in exponent 1. All primes <= B1done
     are already included in exponent of at least 1, so it's save to skip  
     ahead to B1done+1 */
     
  if (*B1done > p)
    {
      getprime_seek ((*B1done) + 1);
      p = getprime ();
    }

  /* then remaining primes > max(sqrt(B1), cascade_limit) and taken 
     with exponent 1 */
  for (; p <= B1; p = getprime ())
  {
    mpz_mul_d (g, g, p, d);
    if (mpz_sizeinbase (g, 2) >= max_size)
      {
        mpres_pow (a, a, g, n);
        mpz_set_ui (g, 1);
        if (stop_asap != NULL && (*stop_asap) ())
          {
            outputf (OUTPUT_NORMAL, "Interrupted at prime %.0f\n", p);
             if (p > *B1done)
              *B1done = p;
            goto clear_pm1_stage1;
          }
        if (chkfilename != NULL && p > last_chkpnt_p + 10000. &&
            elltime (last_chkpnt_time, cputime ()) > CHKPNT_PERIOD)
          {
            writechkfile (chkfilename, ECM_PM1, p, n, NULL, a, NULL);
            last_chkpnt_p = p;
            last_chkpnt_time = cputime ();
          }
      }
  }

  mpres_pow (a, a, g, n);
  
  /* If stage 1 finished normally, p is the smallest prime >B1 here.
     In that case, set to B1 */
  if (p > B1)
      p = B1;
  
  if (p > *B1done)
    *B1done = p;
  
  mpres_sub_ui (a, a, 1, n);
  mpres_gcd (f, a, n);
  if (mpz_cmp_ui (f, 1) > 0)
    youpi = ECM_FACTOR_FOUND_STEP1;
  mpres_add_ui (a, a, 1, n);

 clear_pm1_stage1:
  if (chkfilename != NULL)
    writechkfile (chkfilename, ECM_PM1, *B1done, n, NULL, a, NULL);
  getprime_clear (); /* free the prime tables, and reinitialize */
  mpz_clear (d);
  mpz_clear (g);

  return youpi;
}

/******************************************************************************
*                                                                             *
*                                  Stage 2                                    *
*                                                                             *
******************************************************************************/

/* For each of the nr progressions each of S+1 entries in fd[], performs
   the update fd[k] *= fd[k+1], 0 <= k < S+1. */
static void
update_fd (mpres_t *fd, unsigned int nr, unsigned int S, mpmod_t modulus,
           unsigned long *muls)
{
  unsigned int j, k;
  
  for (j = 0; j < nr * (S + 1); j += S + 1)
    for (k = 0; k < S; k++)
      mpres_mul (fd[j + k], fd[j + k], fd[j + k + 1], modulus);

  if (muls != NULL)
    *muls += (unsigned long) nr * S;
}

/* Puts in F[0..dF-1] the successive values of 

   x^(Dickson_{S, a}(j * d2))
   
     for j == 1 mod 6 , j and d1 coprime, where Dickson_{S, a}
     is the degree S Dickson polynomial with parameter a. For a == 0, 
     Dickson_{S, a} (x) = x^S.
   Uses the x+1/x trick whenever S > 6 and even, then the Dickson 
     parameter a must be 0.
   Requires (dF+1) cells in t for the x+1/x trick.
   Returns non-zero iff a factor was found (then stored in f),
   or an error occurred.
*/

int
pm1_rootsF (mpz_t f, listz_t F, root_params_t *root_params, 
            unsigned long dF, mpres_t *x, listz_t t, mpmod_t modulus)
{
  unsigned long i;
  unsigned long muls = 0, gcds = 0;
  long st, st1;
  pm1_roots_state_t state;
  progression_params_t *params = &state.params; /* for less typing */
  listz_t coeffs;
  mpz_t ts;

  if (dF == 0)
    return 0;

  st = cputime ();

  /* Relative cost of point add during init and computing roots assumed =1 */
  init_roots_params (&state.params, root_params->S, root_params->d1, 
		     root_params->d2, 1.0);

  /* The invtrick is profitable for x^S, S even and > 6. Does not work for 
     Dickson polynomials (root_params->S < 0)! */
  if (root_params->S > 6 && (root_params->S & 1) == 0)
    {
      state.invtrick = 1;
      params->S /= 2;
      params->size_fd = params->nr * (params->S + 1);
    }
  else
    state.invtrick = 0;

  outputf (OUTPUT_DEVVERBOSE, 
	   "pm1_rootsF: state: nr = %d, dsieve = %d, size_fd = %d, S = %d, "
	   "dickson_a = %d, invtrick = %d\n", params->nr, params->dsieve, 
	   params->size_fd, params->S, params->dickson_a, state.invtrick);

  /* Init finite differences tables */
  mpz_init (ts); /* ts = 0 */
  coeffs = init_progression_coeffs (ts, params->dsieve, root_params->d2, 
				    1, 6, params->S, params->dickson_a);
  mpz_clear (ts);

  if (coeffs == NULL)
    return ECM_ERROR;

  /* Allocate memory for fd[] and compute x^coeff[]*/
  state.fd = (mpres_t *) malloc (params->size_fd * sizeof (mpres_t));
  if (state.fd == NULL)
    {
      clear_list (coeffs, params->size_fd);
      return ECM_ERROR;
    }

  for (i = 0; i < params->size_fd; i++) 
    {
      outputf (OUTPUT_TRACE, "pm1_rootsF: coeffs[%d] = %Zd\n", i, coeffs[i]);
      mpres_init (state.fd[i], modulus);
      /* The highest coefficient of all progressions is identical */
      if (i > params->S + 1 && i % (params->S + 1) == params->S)
	{
	  ASSERT (mpz_cmp (coeffs[i], coeffs[params->S]) == 0);
	  mpres_set (state.fd[i], state.fd[params->S], modulus);
	}
      else
        mpres_pow (state.fd[i], *x, coeffs[i], modulus);
    }

  clear_list (coeffs, params->size_fd);
  coeffs = NULL;
  
  st1 = cputime ();
  outputf (OUTPUT_VERBOSE,
	   "Initializing table of differences for F took %ldms\n",
	   elltime (st, st1));
  st = st1;

  /* Now for the actual calculation of the roots. */
  for (i = 0; i < dF;)
    {
      /* Is this a rsieve value where we computed x^Dickson(j * d2) ? */
      if (gcd (params->rsieve, params->dsieve) == 1)
        {
          /* Did we use every progression since the last update? */
          if (params->next == params->nr)
            {
              /* Yes, time to update again */
              update_fd (state.fd, params->nr, params->S, modulus, 
			 &muls);
              params->next = 0;
            }
          
          /* Is this a j value where we want x^Dickson(j * d2) as a root? */
          if (gcd (params->rsieve, root_params->d1) == 1)
            mpres_get_z (F[i++], state.fd[params->next * (params->S + 1)], 
                         modulus);
          params->next ++;
        }
      params->rsieve += 6;
    }

  for (i = 0; i < params->size_fd; i++)
    mpres_clear (state.fd[i], modulus);
  free (state.fd);
  state.fd = NULL;

  if (state.invtrick)
    {
      if (list_invert (t, F, dF, t[dF], modulus)) 
        {
          /* Should never happen */
	  outputf (OUTPUT_ERROR, 
		   "Found factor unexpectedly while inverting F[0]*..*F[dF]\n");
          mpz_set (f, t[dF]);
          return ECM_FACTOR_FOUND_STEP2;
        }
      
      muls += 3 * (dF - 1);
      gcds ++;
      
      for (i = 0; i < dF; i++) 
        {
          mpz_add (F[i], F[i], t[i]);
          mpz_mod (F[i], F[i], modulus->orig_modulus);
        }
    }
  
  outputf (OUTPUT_VERBOSE, "Computing roots of F took %ldms",
	   elltime (st, cputime ()));
  outputf (OUTPUT_DEVVERBOSE, ", %lu muls and %lu extgcds", muls, gcds);
  outputf (OUTPUT_VERBOSE, "\n");
  
  return ECM_NO_FACTOR_FOUND;
}

/* Perform the necessary initialisation to allow computation of

   x^(Dickson_{S, a}(s+n*d))
   
     for successive n, where Dickson_{S, a} is the degree S Dickson
     polynomial with parameter a. For a == 0, Dickson_{S, a} (x) = x^S.
   Uses the x+1/x trick whenever S > 6 and even.
   Return NULL if an error occurred.
*/

pm1_roots_state_t *
pm1_rootsG_init (mpres_t *x, root_params_t *root_params, mpmod_t modulus)
{
  unsigned int i;
  listz_t coeffs;
  pm1_roots_state_t *state;
  progression_params_t *params; /* for less typing */

  state = (pm1_roots_state_t *) malloc (sizeof (pm1_roots_state_t));
  if (state == NULL)
    return NULL;
  params = &(state->params);
  
  params->dickson_a = (root_params->S < 0) ? -1 : 0;
  params->nr = (root_params->d2 > 1) ? root_params->d2 - 1 : 1;
  params->next = 0;
  state->invtrick = (root_params->S > 6 && (root_params->S & 1) == 0);
  params->S = (state->invtrick) ? abs (root_params->S) / 2 : 
                                 abs (root_params->S);
  params->size_fd = params->nr * (params->S + 1);
  params->dsieve = 1;
  params->rsieve = 1;
  
  outputf (OUTPUT_DEVVERBOSE,
	   "pm1_rootsG_init: d1 = %lu, d2 = %lu, state: dsieve = %d, "
	   "nr = %d, size_fd = %d, S = %d, invtrick = %d\n", 
	   root_params->d1, root_params->d2, params->dsieve, 
	   params->nr, params->size_fd, params->S, state->invtrick);
  
  state->fd = (mpres_t *) malloc (params->size_fd * sizeof (mpres_t));
  if (state->fd == NULL)
    {
      free (state);
      return NULL;
    }

  /* Init for Dickson_{E,a} (i0 * d + d1 * n) */
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
      outputf (OUTPUT_TRACE, "pm1_rootsG_init: coeffs[%d] = %Zd\n", 
               i, coeffs[i]);
      mpres_init (state->fd[i], modulus);
      /* The S-th coeff of all progressions is identical */
      if (i > params->S && i % (params->S + 1) == params->S) 
        {
          ASSERT (mpz_cmp (coeffs[i], coeffs[params->S]) == 0);
          /* Simply copy from the first progression */
          mpres_set (state->fd[i], state->fd[params->S], modulus); 
        }
      else
        {
          if (mpz_sgn (coeffs[i]) < 0)
            {
              mpz_neg (coeffs[i], coeffs[i]);
              mpres_pow (state->fd[i], *x, coeffs[i], modulus);
              mpres_invert (state->fd[i],  state->fd[i],  modulus);
              mpz_neg (coeffs[i], coeffs[i]);
            }
          else
            {
              mpres_pow (state->fd[i], *x, coeffs[i], modulus);
            }
        }
    }

  clear_list (coeffs, params->size_fd);
   
  return state;
}

/* Frees all the dynamic variables allocated by pm1_rootsG_init() */

void 
pm1_rootsG_clear (pm1_roots_state_t *state, ATTRIBUTE_UNUSED mpmod_t modulus)
{
  unsigned int k;
  
  for (k = 0; k < state->params.size_fd; k++)
    mpres_clear (state->fd[k], modulus);

  free (state->fd);
  state->fd = NULL;
  
  free (state);
}

/* Puts in G the successive values of 
    
    x^(Dickson_{S, a}(s+j*k))
    
    for 1 <= j <= d, where k is the 'd' value from pm1_rootsG_init()
    and s is the 's' value of pm1_rootsG_init() or where a previous
    call to pm1_rootsG has left off.
   
   Requires (d+1) cells in t for the x+1/x trick.
   Returns non-zero iff a factor was found (then stored in f).
   No error can occur.
*/

int
pm1_rootsG (mpz_t f, listz_t G, unsigned long dF, pm1_roots_state_t *state, 
            listz_t t, mpmod_t modulus)
{
  unsigned long i;
  unsigned long muls = 0, gcds = 0;
  unsigned int st;
  progression_params_t *params = &(state->params); /* for less typing */
  
  outputf (OUTPUT_TRACE,
	   "pm1_rootsG: dF = %d, state: size_fd = %d, nr = %d, S = %d\n",
	   dF, params->size_fd, params->nr, params->S);
  
  st = cputime ();
  
  for (i = 0; i < dF;)
    {
      /* Did we use every progression since the last update? */
      if (params->next == params->nr)
        {
          /* Yes, time to update again */
	  outputf (OUTPUT_TRACE, "pm1_rootsG: Updating table at rsieve = %d\n",
		   params->rsieve);
          update_fd (state->fd, params->nr, params->S, modulus, &muls);
          params->next = 0;
        }
      
      /* Is this a root we should skip? (Take only if gcd == 1) */
      if (gcd (params->rsieve, params->dsieve) == 1)
        {
	  outputf (OUTPUT_TRACE,
		   "pm1_rootsG: Taking root G[%d] at rsieve = %d\n",
		   i, params->rsieve);
          mpres_get_z (G[i++], state->fd[params->next * (params->S + 1)], 
		       modulus);
        }
      else
	outputf (OUTPUT_TRACE, "pm1_rootsG: Skipping root at rsieve = %d\n",
		 params->rsieve);
      
      params->next ++;
      params->rsieve ++;
    }
  
  if (state->invtrick)
    {
      if (list_invert (t, G, dF, t[dF], modulus)) 
        {
	  outputf (OUTPUT_VERBOSE,
		   "Found factor while inverting G[0]*..*G[d]\n");
          mpz_set (f, t[dF]);
          return ECM_FACTOR_FOUND_STEP2;
        }

      muls += 3 * (dF - 1);
      gcds ++;
      
      for (i = 0; i < dF; i++) 
        {
          mpz_add (G[i], G[i], t[i]);
          mpz_mod (G[i], G[i], modulus->orig_modulus);
        }
    }
  
  outputf (OUTPUT_VERBOSE, "Computing roots of G took %ldms",
	   elltime (st, cputime ()));
  outputf (OUTPUT_DEVVERBOSE, ", %lu muls and %lu extgcds", muls, gcds);
  outputf (OUTPUT_VERBOSE, "\n");
  
  return ECM_NO_FACTOR_FOUND;
}


/******************************************************************************
*                                                                             *
*                                Pollard P-1                                  *
*                                                                             *
******************************************************************************/

/* Input: p is the initial generator (sigma), if 0, generate it at random.
          N is the number to factor
	  B1 is the stage 1 bound
	  B2 is the stage 2 bound
	  B1done is the stage 1 limit to which supplied residue has 
	    already been computed
          k is the number of blocks for stage 2
          verbose is the verbosity level
   Output: f is the factor found, p is the residue at end of stage 1
   Return value: non-zero iff a factor is found (1 for stage 1, 2 for stage 2)
*/
int
pm1 (mpz_t f, mpz_t p, mpz_t N, mpz_t go, double *B1done, double B1,
     mpz_t B2min_parm, mpz_t B2_parm, double B2scale, unsigned long k, 
     const int S, int verbose, int repr, int use_ntt, FILE *os, FILE *es, 
     char *chkfilename, char *TreeFilename, double maxmem, 
     gmp_randstate_t rng, int (*stop_asap)(void))
{
  int youpi = ECM_NO_FACTOR_FOUND;
  int base2 = 0;
  int Nbits, smallbase;
  int po2 = 0;    /* Whether we should use power-of-2 poly degree */
  long st;
  mpmod_t modulus;
  mpres_t x;
  mpz_t B2min, B2; /* Local B2, B2min to avoid changing caller's values */
  unsigned long dF;
  root_params_t root_params;
  faststage2_param_t faststage2_params;
  /* If stage2_variant != 0, we use the new fast stage 2 */
  const int stage2_variant = (S == 1 || S == ECM_DEFAULT_S);

  set_verbose (verbose);
  ECM_STDOUT = (os == NULL) ? stdout : os;
  ECM_STDERR = (es == NULL) ? stdout : es;

  /* if n is even, return 2 */
  if (mpz_divisible_2exp_p (N, 1))
    {
      mpz_set_ui (f, 2);
      return ECM_FACTOR_FOUND_STEP1;
    }

  st = cputime ();

  if (mpz_cmp_ui (p, 0) == 0)
    pm1_random_seed (p, N, rng);
  
  mpz_init_set (B2min, B2min_parm);
  mpz_init_set (B2, B2_parm);
  
  /* Set default B2. See ecm.c for comments */
  if (ECM_IS_DEFAULT_B2(B2))
    mpz_set_d (B2, B2scale * pow (B1 * PM1_COST, DEFAULT_B2_EXPONENT));

  /* set B2min */
  if (mpz_sgn (B2min) < 0)
    mpz_set_d (B2min, B1);

  if (repr != ECM_MOD_DEFAULT && repr != ECM_MOD_NOBASE2)
    {
      if (repr == ECM_MOD_MODMULN)
        mpmod_init_MODMULN (modulus, N);
      else if (repr == ECM_MOD_REDC)
        mpmod_init_REDC (modulus, N);
      else if (abs (repr) > 16)
        {
          if (mpmod_init_BASE2 (modulus, repr, N) == ECM_ERROR)
            return ECM_ERROR;
        }
      else
        mpmod_init_MPZ (modulus, N);
    }
  else /* automatic choice */
    {
      /* Find a good arithmetic for this number */
      Nbits = mpz_sizeinbase (N, 2);
      base2 = (repr == 0) ? isbase2 (N, BASE2_THRESHOLD) : 0;
      smallbase = mpz_fits_uint_p (p);

      /* TODO: make dependent on Nbits and base2 */
      if (base2)
        {
          mpmod_init_BASE2 (modulus, base2, N);
        }

      else if (mpz_size (N) <= 2 * POWM_THRESHOLD && smallbase && B1 <= 1e6)
      /* Below POWM_THRESHOLD, mpz_powm uses MODMULN reduction, too, but 
         without special code for small bases which makes our MODMULN
         faster. Above POWM_THRESHOLD mpz_powm uses faster mod reduction,
         at about 2*POWM_THRESHOLD it catches up with our smallbase-MODMULN
         and then is faster until REDC takes over. */
        {
	  outputf (OUTPUT_VERBOSE, "Using MODMULN\n");
          mpmod_init_MODMULN (modulus, N);
        }
      else if (Nbits > 50000 ||  (Nbits > 3500 && smallbase))
        {
	  outputf (OUTPUT_VERBOSE, "Using REDC\n");
          mpmod_init_REDC (modulus, N);
        }
      else
        {
	  outputf (OUTPUT_VERBOSE, "Using mpz_powm\n");
          mpmod_init_MPZ (modulus, N);
        }
    }
  

  /* Determine parameters (polynomial degree etc.) */

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
	  unsigned long t;
	  /* See what transform length that the NTT can handle (due to limited 
	     primes and limited memory) */
	  t = mpzspm_max_len (N);
	  lmax_NTT = MIN (lmax, t);
	  if (maxmem != 0.)
	    {
	      t = pm1fs2_maxlen ((size_t) maxmem, N, use_ntt);
	      lmax_NTT = MIN (lmax_NTT, t);
	    }
	  outputf (OUTPUT_DEVVERBOSE, "NTT can handle lmax <= %lu\n", lmax_NTT);
	}
      
      /* See what transform length that the non-NTT code can handle */
      lmax_noNTT = lmax;
      if (maxmem != 0.)
	{
	  unsigned long t;
	  t = pm1fs2_maxlen ((size_t) maxmem, N, 0);
	  lmax_noNTT = MIN (lmax_noNTT, t);
	  outputf (OUTPUT_DEVVERBOSE, "non-NTT can handle lmax <= %lu\n", 
		   lmax_noNTT);
	}
      
      P = choose_P (B2min, B2, MAX (lmax_NTT, lmax_noNTT), k, 
		    &faststage2_params, B2min, B2, use_ntt);
      if (P == ECM_ERROR)
        {
          mpz_clear (faststage2_params.m_1);
          return ECM_ERROR;
        }
      
      /* See if the selected parameters let us use NTT or not */
      if (faststage2_params.l > lmax_NTT)
	use_ntt = 0;
      
      if (maxmem != 0.)
	  outputf (OUTPUT_VERBOSE, "Using lmax = %lu with%s NTT which takes "
		   "about %luMB of memory\n", faststage2_params.l, 
		   (use_ntt) ? "" : "out", 
		   pm1fs2_memory_use (faststage2_params.l, N, use_ntt)/1048576);
    }
  else
    {
      mpz_init (root_params.i0);
      root_params.d2 = 0; /* Enable automatic choice of d2 */
      
      if (use_ntt || (modulus->repr == ECM_MOD_BASE2 && modulus->Fermat > 0))
	po2 = 1;

      if (bestD (&root_params, &k, &dF, B2min, B2, po2, use_ntt, maxmem,
                 (TreeFilename != NULL), modulus) == ECM_ERROR)
	{
	  youpi = ECM_ERROR;
	  goto clear_and_exit;
	}
  
      root_params.S = S;
      /* Set default degree for Brent-Suyama extension */
      if (root_params.S == ECM_DEFAULT_S)
	{
	  if (modulus->repr == ECM_MOD_BASE2 && modulus->Fermat > 0)
	    {
	      /* For Fermat numbers, default is 2 (no Brent-Suyama) */
	      root_params.S = 2;
	    }
	  else
	    {
	      mpz_t t;
	      mpz_init (t);
	      mpz_sub (t, B2, B2min);
	      if (mpz_cmp_d (t, 3.5e5) < 0) /* B1 < 50000 */
		root_params.S = -4; /* Dickson polys give a slightly better chance of success */
	      else if (mpz_cmp_d (t, 1.1e7) < 0) /* B1 < 500000 */
		root_params.S = -6;
	      else if (mpz_cmp_d (t, 1.25e8) < 0) /* B1 < 3000000 */
		root_params.S = 12; /* but for S>6, S-th powers are faster thanks to invtrick */
	      else if (mpz_cmp_d (t, 7.e9) < 0) /* B1 < 50000000 */
		root_params.S = 24;
	      else if (mpz_cmp_d (t, 1.9e10) < 0) /* B1 < 100000000 */
		root_params.S = 48;
	      else if (mpz_cmp_d (t, 5.e11) < 0) /* B1 < 1000000000 */
		root_params.S = 60;
	      else
		root_params.S = 120;
	      mpz_clear (t);
	    }
	}
      
      /* We need Suyama's power even and at least 2 for P-1 stage 2 to work 
	 correctly */

      if (root_params.S & 1)
	root_params.S *= 2; /* FIXME: Is this what the user would expect? */
    }
  
  /* Print B1, B2, polynomial and x0 */
  print_B1_B2_poly (OUTPUT_NORMAL, ECM_PM1, B1, *B1done, B2min_parm, B2min, 
		    B2, (stage2_variant == 0) ? root_params.S : 1, p, 0);

  if (stage2_variant != 0)
    outputf (OUTPUT_VERBOSE, "P = %lu, l = %lu, s_1 = %lu, s_2 = %lu, "
	     "m_1 = %Zd\n", faststage2_params.P, faststage2_params.l,
	     faststage2_params.s_1,faststage2_params.s_2,
	     faststage2_params.m_1);
  else
    outputf (OUTPUT_VERBOSE, "dF=%lu, k=%lu, d=%lu, d2=%lu, i0=%Zd\n", 
	     dF, k, root_params.d1, root_params.d2, root_params.i0);

  mpres_init (x, modulus);
  mpres_set_z (x, p, modulus);

  st = cputime ();

  if (B1 > *B1done)
    youpi = pm1_stage1 (f, x, modulus, B1, B1done, go, stop_asap, chkfilename);

  st = elltime (st, cputime ());

  outputf (OUTPUT_NORMAL, "Step 1 took %ldms\n", st);
  if (test_verbose (OUTPUT_RESVERBOSE))
    {
      mpz_t tx;
      mpz_init (tx);
      mpres_get_z (tx, x, modulus);
      outputf (OUTPUT_RESVERBOSE, "x=%Zd\n", tx);
      mpz_clear (tx);
    }

  if (stop_asap != NULL && (*stop_asap) ())
    goto clear_and_exit;

  if (youpi == ECM_NO_FACTOR_FOUND && mpz_cmp (B2, B2min) >= 0)
    {
      if (stage2_variant != 0)
        {
          if (use_ntt)
            youpi = pm1fs2_ntt (f, x, modulus, &faststage2_params);
          else
            youpi = pm1fs2 (f, x, modulus, &faststage2_params);
        }
      else
        youpi = stage2 (f, &x, modulus, dF, k, &root_params, ECM_PM1, 
                        use_ntt, TreeFilename, stop_asap);
    }

clear_and_exit:
  mpres_get_z (p, x, modulus);
  mpres_clear (x, modulus);
  mpmod_clear (modulus);
  if (stage2_variant != 0)
    mpz_clear (faststage2_params.m_1);
  else
    mpz_clear (root_params.i0);
  mpz_clear (B2);
  mpz_clear (B2min);

  return youpi;
}
