/* Pollard 'P-1' algorithm.

  Copyright 2001, 2002, 2003 Paul Zimmermann and Alexander Kruppa.

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
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "ecm.h"

/* #define DEBUG */

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

/* put in 'a' a valid random seed for P-1, i.e. gcd(a, n)=1 and a <> {-1,1} */
void
pm1_random_seed (mpres_t a, mpz_t n, gmp_randstate_t randstate)
{
  mpz_t q;

  mpz_init (q);
  do
    {
      mpz_urandomb (a, randstate, 32);
      mpz_gcd (q, a, n);
    }
  while (mpz_cmp_ui (q, 1) != 0 || mpz_cmp_ui (a, 1) == 0 ||
         mpz_cmp_si (a, -1) == 0);
  mpz_clear (q);
}

/*** Cascaded multiply ***/

static mul_casc *
mulcascade_init (void)
{
  mul_casc *t;
  t = (mul_casc *) malloc (sizeof (mul_casc));
  if (t == NULL)
    {
      fprintf (stderr, "mulcascade_init: could not allocate memory\n");
      exit (EXIT_FAILURE);
    }
  t->val = (mpz_t*) malloc (sizeof (mpz_t));
  if (t->val == NULL)
    {
      fprintf (stderr, "mulcascade_init: could not allocate memory\n");
      exit (EXIT_FAILURE);
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
mulcascade_mul_d (mul_casc *c, const double n, mpz_t t)
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
   Return value: non-zero iff a factor was found.
*/

static int
pm1_stage1 (mpz_t f, mpres_t a, mpmod_t n, double B1, double B1done,
	    int verbose, mpz_t go)
{
  double B0, p, q, r, cascade_limit;
  mpz_t g, d;
  int youpi;
  unsigned int size_n, max_size;
  unsigned int smallbase = 0;
  mul_casc *cascade;

  mpz_init (g);
  mpz_init (d);

  B0 = sqrt (B1);

  size_n = mpz_sizeinbase (n->orig_modulus, 2);
  max_size = L1 * size_n;

  mpres_get_z (g, a, n);
  if (mpz_fits_uint_p (g))
    {
      smallbase = mpz_get_ui (g);
    }

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

  /* since B0 = sqrt(B1), we can have B0 > cascade_limit only when
     B1 > cascade_limit^2. This cannot happen when cascade_limit=B1,
     thus we need B1 > min(CASCADE_MAX, 3000*sizeinbase(n,2))^2.
     For sizeinbase(n,2) <= CASCADE_MAX/3000 (less than 5017 digits 
     for CASCADE_MAX=5e7) this means B1 > 9e6*sizeinbase(n,2)^2.
     For sizeinbase(n,2) > CASCADE_MAX/3000, this means B1 > CASCADE_MAX^2,
     i.e. B1 > 25e14 for CASCADE_MAX=5e7.
*/

  /* if the user "knows" that P-1 has a given, he/she can "enter" it */
  if (mpz_cmp_ui (go, 1) > 0)
    cascade = mulcascade_mul (cascade, go);

  if (B0 <= cascade_limit)
    {
      /* first loop through small primes <= sqrt(B1) */
      for (p = 2.0; p <= B0; p = getprime (p))
        {
          for (q = 1, r = p; r <= B1; r *= p)
            if (r > B1done) q *= p;
          cascade = mulcascade_mul_d (cascade, q, d);
        }

      /* then all sqrt(B1) < primes < cascade_limit and taken with 
         exponent 1 */
      for ( ; p <= cascade_limit; p = getprime (p))
        if (p > B1done)
          cascade = mulcascade_mul_d (cascade, p, d);
   
      mulcascade_get_z (g, cascade);
      mulcascade_free (cascade);
#ifdef DEBUG
      printf ("Exponent has %u bits\n", mpz_sizeinbase (g, 2));
#endif
      if (smallbase)
        {
	  if (verbose > 1)
	    printf ("Using mpres_ui_pow, base %u\n", smallbase);
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
      for (p = 2.0; p <= cascade_limit; p = getprime (p))
        {
          for (q = 1.0, r = p; r <= B1; r *= p)
            if (r > B1done) q *= p;
          cascade = mulcascade_mul_d (cascade, q, d);
        }
      
      mulcascade_get_z (g, cascade);
      mulcascade_free (cascade);
#ifdef DEBUG
      printf("Exponent has %u bits\n", mpz_sizeinbase (g, 2));
#endif
      if (smallbase)
        {
	  if (verbose > 1)
	    printf ("Using mpres_ui_pow, base %u\n", smallbase);
          mpres_ui_pow (a, smallbase, g, n);
        }
      else
        {
          mpres_pow (a, a, g, n);
        }
      mpz_set_ui (g, 1);
      
      for ( ; p <= B0; p = getprime (p))
        {
          for (q = 1, r = p; r <= B1; r *= p)
            if (r > B1done) q *= p;
          mpz_mul_d (g, g, q, d);
          if (mpz_sizeinbase (g, 2) >= max_size)
            {
              mpres_pow (a, a, g, n);
              mpz_set_ui (g, 1);
            }
        }
    }
  
  /* then remaining primes > max(sqrt(B1), cascade_limit) and taken 
     with exponent 1 */
  for (; p <= B1; p = getprime (p))
  {
    if (p > B1done)
      {
        mpz_mul_d (g, g, p, d);
        if (mpz_sizeinbase (g, 2) >= max_size)
	  {
	    mpres_pow (a, a, g, n);
	    mpz_set_ui (g, 1);
	  }
      }
  }

  getprime (FREE_PRIME_TABLE); /* free the prime tables, and reinitialize */

  mpz_clear (d);
  mpres_pow (a, a, g, n);
  mpz_clear (g);
  
  mpres_sub_ui (a, a, 1, n);
  mpres_gcd (f, a, n);
  youpi = mpz_cmp_ui (f, 1);
  mpres_add_ui (a, a, 1, n);

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
  
  *muls += nr * S;
}

/* Puts in F[0..dF-1] the successive values of 

   x^(Dickson_{S, a}(j * d2))
   
     for j == 1 mod 6 , j and d1 coprime, where Dickson_{S, a}
     is the degree S Dickson polynomial with parameter a. For a == 0, 
     Dickson_{S, a} (x) = x^S.
   Returns non-zero iff a factor was found (then stored in f).
   Uses the x+1/x trick whenever S > 6 and even, then the Dickson 
     parameter a must be 0.
   Requires (dF+1) cells in t for the x+1/x trick.
*/

int
pm1_rootsF (mpz_t f, listz_t F, unsigned int d1, unsigned int d2, 
            unsigned int dF, mpres_t *x, listz_t t, int S, mpmod_t modulus, 
            int verbose)
{
  unsigned int i;
  unsigned long muls = 0, gcds = 0;
  int st;
  listz_t coeffs;
  pm1_roots_state state;
  
  if (dF == 0)
    return 0;

  if (verbose >= 2)
    st = cputime ();

  /* Relative cost of point add during init and computing roots assumed =1 */
  /* The typecast from hell: the relevant fields of ecm_roots_state and
     pm1_roots_state match in position so init_roots_state() can init a
     pm1_roots_state as well. UGLY. OOP would really help here */
  init_roots_state ((ecm_roots_state *) &state, S, d1, d2, 1.0);

  /* The invtrick is profitable for x^S, S even and > 6 */
  if (S > 6 && (S & 1) == 0)
    {
      state.invtrick = 1;
      state.S /= 2;
      state.size_fd = state.nr * (state.S + 1);
    } else
    state.invtrick = 0;

  if (verbose >= 3)
    printf ("pm1_rootsF: state: nr = %d, dsieve = %d, size_fd = %d, S = %d, "
            "dickson_a = %d, invtrick = %d\n", state.nr, state.dsieve, 
            state.size_fd, state.S, state.dickson_a, state.invtrick);
  
  /* Init finite differences tables */
  
  coeffs = init_progression_coeffs (0., state.dsieve, d2, 1, 6, state.S, 
                                    state.dickson_a);
  
  /* Allocate memory for fd[] and compute x^coeff[]*/
  
  state.fd = (mpres_t *) xmalloc (state.size_fd * sizeof (mpres_t));
  for (i = 0; i < state.size_fd; i++) 
    {
      if (verbose >= 4)
        gmp_printf ("pm1_rootsF: coeffs[%d] = %Zd\n", i, coeffs[i]);
      mpres_init (state.fd[i], modulus);
      /* The highest coefficient of all progressions is identical */
      if (i > state.S + 1 && i % (state.S + 1) == state.S)
        mpres_set (state.fd[i], state.fd[state.S], modulus);
      else
        mpres_pow (state.fd[i], *x, coeffs[i], modulus);
    }

  clear_list (coeffs, state.size_fd);
  coeffs = NULL;
  
  if (verbose >= 2)
    {
      int st1;
      
      st1 = cputime ();
      printf ("Initializing tables of differences for F took %dms\n", 
              st1 - st);
      st = st1;
    }

  /* Now for the actual calculation of the roots. */
  
  for (i = 0; i < dF;)
    {
      /* Is this a rsieve value where we computed x^Dickson(j * d2) ? */
      if (gcd (state.rsieve, state.dsieve) == 1)
        {
          /* Did we use every progression since the last update? */
          if (state.next == state.nr)
            {
              /* Yes, time to update again */
              update_fd (state.fd, state.nr, state.S, modulus, &muls);
              
              state.next = 0;
            }
          
          /* Is this a j value where we want x^Dickson(j * d2) as a root? */
          if (gcd (state.rsieve, d1) == 1)
            mpres_get_z (F[i++], state.fd[state.next * (state.S + 1)], modulus);
            
          state.next++;
        }
      state.rsieve += 6;
    }

  for (i = 0; i < state.size_fd; i++)
    mpres_clear (state.fd[i], modulus);
  free (state.fd);
  state.fd = NULL;

  if (state.invtrick)
    {
      if (list_invert (t, F, dF, t[dF], modulus)) 
        {
          if (verbose >= 2)
            printf ("Found factor while inverting F[0]*..*F[d]\n");
          mpz_set (f, t[dF]);
          return 1;
        }
      
      muls += 3 * (dF - 1);
      gcds++;
      
      for (i = 0; i < dF; i++) 
        {
          mpz_add (F[i], F[i], t[i]);
          mpz_mod (F[i], F[i], modulus->orig_modulus);
        }
    }
  
  if (verbose >= 2)
    {
      printf ("Computing roots of F took %dms", cputime () - st);
      if (verbose > 2)
        printf (", %lu muls and %lu extgcds", muls, gcds);
      printf ("\n");
    }
  
  return 0;
}

/* Perform the neccessary initialisation to allow computation of

   x^(Dickson_{S, a}(s+n*d))
   
     for successive n, where Dickson_{S, a} is the degree S Dickson
     polynomial with parameter a. For a == 0, Dickson_{S, a} (x) = x^S.
   Uses the x+1/x trick whenever S > 6 and even, then the Dickson
     parameter a must be 0.
*/

pm1_roots_state *
pm1_rootsG_init (mpres_t *x, double s, unsigned int d1, unsigned int d2, 
                 int S, int verbose, mpmod_t modulus)
{
  unsigned int i;
  int dickson_a;
  listz_t coeffs;
  pm1_roots_state *state;

  state = (pm1_roots_state *) xmalloc (sizeof (pm1_roots_state));
  
  dickson_a = (S < 0) ? -1 : 0;
  
  state->nr = (d2 > 1) ? d2 - 1 : 1;
  state->next = 0;
  state->invtrick = (S > 6 && (S & 1) == 0);
  state->S = (state->invtrick) ? abs (S) / 2 : abs (S);
  state->size_fd = state->nr * (state->S + 1);
  state->dsieve = 1;
  state->rsieve = 1;
  
  if (verbose >= 3)
    printf ("pm1_rootsG_init: d1 = %d, d2 = %d, state: dsieve = %d, nr = %d, "
            "size_fd = %d, S = %d, invtrick = %d\n", d1, d2, state->dsieve, 
            state->nr, state->size_fd, state->S, state->invtrick);
  
  state->fd = (mpres_t *) xmalloc (state->size_fd * sizeof (mpres_t));

  /* Init for Dickson_{E,a} (s + d1 * n) */
  coeffs = init_progression_coeffs (s, d2, d1, 1, 1, state->S, dickson_a);

  for (i = 0; i < state->size_fd; i++) 
    {
      /* gmp_printf ("pm1_rootsG_init: coeffs[%d] = %Zd\n", i, coeffs[i]); */
      mpres_init (state->fd[i], modulus);
      /* The S-th coeff of all progressions is identical */
      if (i > state->S + 1 && i % (state->S + 1) == state->S) 
        {
#ifdef DEBUG
          if (mpz_cmp(coeffs[i], coeffs[state->S]) != 0)
            {
              fprintf (stderr, "pm1_rootsG_init: coeffs[%d] != coeffs[%d]\n", 
                       i, state->S);
              exit (EXIT_FAILURE);
            }
#endif
          /* Simply copy from the first progression */
          mpres_set (state->fd[i], state->fd[state->S], modulus); 
        }
      else
        mpres_pow (state->fd[i], *x, coeffs[i], modulus);
    }

  clear_list (coeffs, state->size_fd);
   
  return state;
}

/* Frees all the dynamic variables allocated by pm1_rootsG_init() */

void 
pm1_rootsG_clear (pm1_roots_state *state, UNUSED mpmod_t modulus)
{
  unsigned int k;
  
  for (k = 0; k < state->size_fd; k++)
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
   
   Returns non-zero iff a factor was found (then stored in f).
   Requires (d+1) cells in t for the x+1/x trick.
*/

int
pm1_rootsG (mpz_t f, listz_t G, unsigned int dF, pm1_roots_state *state, 
            listz_t t, mpmod_t modulus, int verbose)
{
  unsigned int i;
  unsigned long muls = 0, gcds = 0;
  int st;
  
  if (verbose >= 4)
    printf ("pm1_rootsG: dF = %d, state: size_fd = %d, nr = %d, S = %d\n",
            dF, state->size_fd, state->nr, state->S);
  
  st = cputime ();
  
  for (i = 0; i < dF;)
    {
      /* Did we use every progression since the last update? */
      if (state->next == state->nr)
        {
          /* Yes, time to update again */
          if (verbose >= 4)
            printf ("pm1_rootsG: Updating table at rsieve = %d\n", state->rsieve);
          
          update_fd (state->fd, state->nr, state->S, modulus, &muls);
          state->next = 0;
        }
      
      /* Is this a root we should skip? (Take only if gcd == 1) */
      if (gcd(state->rsieve, state->dsieve) == 1)
        {
          if (verbose >= 4)
            printf ("pm1_rootsG: Taking root G[%d] at rsieve = %d\n", i, state->rsieve);
          mpres_get_z (G[i++], state->fd[state->next * (state->S + 1)], modulus);
        }
      else
        if (verbose >= 4)
          printf ("pm1_rootsG: Skipping root at rsieve = %d\n", state->rsieve);
        
      
      state->next++;
      state->rsieve++;
    }
  
  if (state->invtrick)
    {
      if (list_invert (t, G, dF, t[dF], modulus)) 
        {
          if (verbose >= 2)
            printf ("Found factor while inverting G[0]*..*G[d]\n");
          mpz_set (f, t[dF]);
          return 1;
        }

      muls += 3 * (dF - 1);
      gcds++;
      
      for (i = 0; i < dF; i++) 
        {
          mpz_add (G[i], G[i], t[i]);
          mpz_mod (G[i], G[i], modulus->orig_modulus);
        }
    }
  
  if (verbose >= 2)
    {
      printf ("Computing roots of G took %dms", cputime () - st);
      if (verbose > 2)
        printf (", %lu muls and %lu extgcds", muls, gcds);
      printf ("\n");
    }
  
  return 0;
}


/******************************************************************************
*                                                                             *
*                                Pollard P-1                                  *
*                                                                             *
******************************************************************************/

/* Input: p is the initial generator (sigma)
          n is the number to factor
	  B1 is the stage 1 bound
	  B2 is the stage 2 bound
	  B1done is the stage 1 limit to which supplied residue has 
	    already been computed
          k is the number of blocks for stage 2
          verbose is the verbose level: 0=quiet, 1=normal, 2=verbose
   Output: f is the factor found, p is the residue at end of stage 1
   Return value: non-zero iff a factor is found (1 for stage 1, 2 for stage 2)
*/
int
pm1 (mpz_t f, mpz_t p, mpz_t N, mpz_t go, double B1done, double B1,
     double B2min, double B2, double B2scale, unsigned int k, int S,
     int verbose, int repr)
{
  mpmod_t modulus;
  mpres_t x;
  int youpi = 0, st, base2, Nbits, smallbase;

  st = cputime ();
  
  /* Set default B2. See ecm.c for comments */
  if (IS_DEFAULT_B2(B2))
    B2 = pow (B1 / 6.0, 1.424828748);

  /* Scale B2 by what the user said (or by the default scaling of 1.0) */
  B2 *= B2scale;

  /* Set default degree for Brent-Suyama extension */
  
  if (S == 0)
    {
      if (B2 - B2min < 3.5e5) /* B1 < 50000 */
        S = -4; /* Dickson polys give a slightly better chance of success */
      else if (B2 - B2min < 1.1e7) /* B1 < 500000 */
        S = -6;
      else if (B2 - B2min < 1.25e8) /* B1 < 3000000 */
        S = 12; /* but for S>6, S-th powers are faster thanks to invtrick */
      else if (B2 - B2min < 7.e9) /* B1 < 50000000 */
        S = 24;
      else if (B2 - B2min < 1.9e10) /* B1 < 100000000 */
        S = 48;
      else if (B2 - B2min < 5.e11) /* B1 < 1000000000 */
        S = 60;
      else
        S = 120;
    }

  /* We need Suyama's power even and at least 2 for P-1 stage 2 to work 
     correctly */
  if (abs(S) < 2)
    S = 2;

  if (S & 1)
    S *= 2; /* FIXME: Is this what the user would expect? */
  
  if (verbose >= 1)
    {
      printf ("Using ");
      if (IS_DEFAULT_B1_DONE(B1done))
        printf("B1=%1.0f", B1);
      else
        printf("B1=%1.0f-%1.0f", B1done, B1);
      if (B2min <= B1)
        printf(", B2=%1.0f, ", B2);
      else
        printf(", B2=%1.0f-%1.0f, ", B2min, B2);
      if (S > 0)
        printf("polynomial x^%u", S);
      else
        printf("polynomial Dickson(%u)", -S);

      if (IS_DEFAULT_B1_DONE(B1done) || verbose > 1) 
	/* don't print in resume case, since x0 is saved in resume file */
	{
	  printf (", x0=");
	  mpz_out_str (stdout, 10, p);
	}
      printf ("\n");
      fflush (stdout);
    }

  if (repr > 0) /* repr = 0 is the default, -1 means nobase2 */
    {
      if (repr == 2)
        mpmod_init_MODMULN (modulus, N);
      else if (repr == 3)
        mpmod_init_REDC (modulus, N);
      else if (repr > 16)
        mpmod_init_BASE2 (modulus, repr, N);
      else
        mpmod_init_MPZ (modulus, N);
    }
  else
    {
      /* Find a good arithmetic for this number */
      Nbits = mpz_sizeinbase (N, 2);
      base2 = (repr == 0) ? isbase2 (N, BASE2_THRESHOLD) : 0;
      smallbase = mpz_fits_uint_p (p);

      /* TODO: make dependent on Nbits and base2 */
      if (base2)
        {
	  if (verbose > 1)
	    printf ("Using special division for factor of 2^%d%c1\n", 
		    abs (base2), (base2 > 0) ? '+' : '-');
          mpmod_init_BASE2 (modulus, base2, N);
        }

      else if (mpz_size (N) <= 2 * POWM_THRESHOLD && smallbase && B1 <= 1e6)
      /* Below POWM_THRESHOLD, mpz_powm uses MODMULN reduction, too, but 
         without special code for small bases which makes our MODMULN
         faster. Above POWM_THRESHOLD mpz_powm uses faster mod reduction,
         at about 2*POWM_THRESHOLD it catches up with our smallbase-MODMULN
         and then is faster until REDC takes over. */
        {
	  if (verbose > 1)
	    printf ("Using MODMULN\n");
          mpmod_init_MODMULN (modulus, N);
        }
      else if (Nbits > 50000 ||  (Nbits > 3500 && smallbase))
        {
	  if (verbose > 1)
	    printf ("Using REDC\n");
          mpmod_init_REDC (modulus, N);
        }
      else
        {
	  if (verbose > 1)
	    printf ("Using mpz_powm\n");
          mpmod_init_MPZ (modulus, N);
        }
    }
  
  mpres_init (x, modulus);
  mpres_set_z (x, p, modulus);

  if (B1 > B1done)
    youpi = pm1_stage1 (f, x, modulus, B1, B1done, verbose, go);

  st = cputime() - st;

  if (verbose >= 1)
    {
      printf ("Step 1 took %dms\n", st);
      fflush (stdout);
    }

  if (verbose >= 2)
    {
      printf ("x=");
      mpres_out_str (stdout, 10, x, modulus);
      printf ("\n");
      fflush (stdout);
    }

  if (youpi != 0) /* a factor was found */
    goto clear_and_exit;

  youpi = stage2 (f, &x, modulus, B2min, B2, k, S, verbose, PM1_METHOD, st);

clear_and_exit:
  mpres_get_z (p, x, modulus);
  mpres_clear (x, modulus);
  mpmod_clear (modulus);

  return youpi;
}
