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
#ifdef WANT_GMP_IMPL
#include "gmp-mparam.h"
#endif
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

int      pm1_stage1     (mpz_t, mpres_t, mpmod_t, double, double, int, mpz_t, mpz_t, mpz_t);
mul_casc *mulcascade_init (void);
void     mulcascade_free (mul_casc *);
mul_casc *mulcascade_mul_d (mul_casc *c, const double n, mpz_t t);
mul_casc *mulcascade_mul   (mul_casc *c, mpz_t n);
void     mulcascade_get_z (mpz_t, mul_casc *);

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

mul_casc *
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

void 
mulcascade_free (mul_casc *c)
{
  unsigned int i;
  for (i = 0; i < c->size; i++)
    mpz_clear (c->val[i]);
  free (c->val);
  free (c);
}

mul_casc * 
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

mul_casc * 
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

void 
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
   Output: f is the factor found, a is the value at end of stage 1
   Return value: non-zero iff a factor was found.
*/

int
pm1_stage1 (mpz_t f, mpres_t a, mpmod_t n, double B1, double B1done,
	    int verbose, mpz_t orig_n, mpz_t orig_X0, mpz_t go)
{
  double B0, p, q, r, cascade_limit;
  mpz_t g, d;
  int youpi;
  unsigned int size_n, max_size;
  unsigned int smallbase = 0;
  mul_casc *cascade;
  int Counter = 0, st_save;

  mpz_init (g);
  mpz_init (d);

  /* Prep for stage one counter */
  showscreenticks_change_stage(1);

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

  /* If the user "knows" that P-1 factors always have a certain form, then the user can "enter" that known factor */
  /* HOWEVER, this is done in a double, so the user must only load in at most 53 bits or so? */
  if (mpz_get_ui(go) > 1 || mpz_size (go) > 1)
    cascade = mulcascade_mul (cascade, go);

  if (B0 <= cascade_limit)
    {
      /* first loop through small primes <= sqrt(B1) */
      for (p = 2.0; p <= B0; p = getprime(p))
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
      for (p = 2.0; p <= cascade_limit; p = getprime(p))
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
      
      for ( ; p <= B0; p = getprime(p))
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
  
  /* update the screen mode after the cascade work is done */
  st_save = cputime ();
  showscreenticks(1,(int) (100.0 * (double) p / (double) B1));

  /* then remaining primes > max(sqrt(B1), cascade_limit) and taken 
     with exponent 1 */
  for (; p <= B1; p = getprime(p))
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
    if (++Counter == 250)
      {
        showscreenticks(1,(int) (100.0 * (double) p / (double) B1));
	Counter=0;
	  /* should we save the current "ecm_wip.sav" file??? It is saved every 15 minutes */
	  /* NOTE this saving DOES NOT save the expression.  It is just a "fail-safe" measure */
#if defined (DEBUG_AUTO_SAVE)
        if (cputime () - st_save > 2000)
#else
        if (cputime () - st_save > 15 * 60 * 1000)
#endif
	  {
	    st_save = cputime ();
	    /* sigma and A are not needed for the save.  Simply create "dummies" here to pass in */
/*	    mpz_t sigma, A, X;
	    mpz_init_set_ui (sigma, 0);
	    mpz_init_set_ui (A, 0);
	    mpz_init (X);
*/
	    /* Suck the X value out of a */
	    /*mpres_get_z (X, a, n);*/
/*
	    mpres_get_z (X, g, n);
		write_temp_resumefile (PM1_METHOD, p, sigma, A, X, orig_n, orig_X0, verbose);
	    mpz_clear (X);
	    mpz_clear (A);
	    mpz_clear (sigma);
*/		
	    /* Reset our "timeout" value, so we save again in 15 minutes */
	  }
	  /*  This "testing" code is here to see just how often this ++Counter loop is entered.
	  {
	    static int x;
  	    fprintf (stderr, "1:%02d  p=%.0f\r", ++x, p);
	  }
	  */
      }
  }

  getprime (0.0); /* free the prime tables, and reinitialize */

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

/* Puts in F[0..dF-1] the successive values of 

   x^(Dickson_{S, a}(j))
   
     for 0 < j = 1 mod 6 < d, j and d coprime, where Dickson_{S, a}
     is the degree S Dickson polynomial with parameter a. For a == 0, 
     Dickson_{S, a} (x) = x^S.
   Returns non-zero iff a factor was found (then stored in f).
   Uses the x+1/x trick whenever S > 6 and even, then the Dickson 
     parameter a must be 0.
   Requires (dF+1) cells in t for the x+1/x trick.
*/

int
pm1_rootsF (mpz_t f, listz_t F, unsigned int d, unsigned int dF, mpres_t *x,
        listz_t t, int S, mpmod_t modulus, int verbose, 
        unsigned long *tot_muls)
{
  unsigned int i, j, k, muls = 0;
  int st, st2;
  mpres_t *fd;
  listz_t coeffs;
  int invtrick = 0, dickson_a = 0;

  st = cputime ();

  mpres_get_z (F[0], *x, modulus); /* s^1 for P-1 */

  if (S > 6 && (S & 1) == 0) /* If we use S-th power, S > 6 */
    {
      invtrick = 1; /* then the invtrick is profitable */
      S /= 2;
    }
  else if (S < 0)
    {
      dickson_a = -1;
      S = -S;
    }

  if (d > 7)
    {
      st2 = cputime ();
      coeffs = init_list (S + 1);
      
      fin_diff_coeff (coeffs, 7, 6, S, dickson_a);
      
      fd = (mpres_t *) malloc ((S + 1) * sizeof (mpres_t));
      if (fd == NULL)
        {
          fprintf (stderr, "Error: not enough memory\n");
          exit (EXIT_FAILURE);
        }
      for (k = 0; k <= (unsigned) S; k++) 
        {
          mpres_init (fd[k], modulus);
          mpres_pow (fd[k], *x, coeffs[k], modulus);
        }

      clear_list (coeffs, S + 1);
      
      if (verbose >= 2)
        printf ("Initializing table of differences for F took %dms\n", cputime () - st2);

      for (j = 7, i = 1; i < dF; j += 6)
        {
          if (gcd (j, d) == 1)
            mpres_get_z (F[i++], fd[0], modulus);
          
          /* Compute value of f_{S, a}(7+n*6) for the next n */
          for (k = 0; k < (unsigned) S; k++)
            mpres_mul (fd[k], fd[k], fd[k+1], modulus);
          
          muls += S;
        }
      
      for (k = 0; k <= (unsigned) S; k++)
        mpres_clear (fd[k], modulus);
      free (fd);
    }

  if (invtrick)
    {
      if (list_invert (t, F, dF, t[dF], modulus)) 
        {
          if (verbose >= 2)
            printf ("Found factor while inverting F[0]*..*F[d]\n");
          mpz_set (f, t[dF]);
          return 1;
        }
      
      muls += 3 * dF;
      
      for (j = 0; j < dF; j++) 
        {
          mpz_add (F[j], F[j], t[j]);
          mpz_mod (F[j], F[j], modulus->orig_modulus);
        }
    }
  
  if (verbose >= 2)
    printf ("Computing roots of F took %dms and %d muls\n", cputime () - st, 
            muls);
  
  if (tot_muls != NULL)
    *tot_muls += muls;

  return 0;
}

/* Perform the neccessary initialisation to allow computation of

   x^(Dickson_{S, a}(s+n*d))
   
     for successive n, where Dickson_{S, a} is the degree S Dickson
     polynomial with parameter a. For a == 0, Dickson_{S, a} (x) = x^S.
   Uses the x+1/x trick whenever S > 6 and even, then the Dickson
     parameter a must be 0.
*/

mpres_t *
pm1_rootsG_init (mpres_t *x, double s, unsigned int d, int S,
                 mpmod_t modulus)
{
  unsigned int k;
  int dickson_a = 0;
  listz_t coeffs;
  mpres_t *fd;

  if (S > 6 && (S & 1) == 0)
    {
      S /= 2;
    }
  else if (S < 0)
    {
      dickson_a = -1;
      S = -S;
    }
  
  coeffs = init_list (S + 1);

  fin_diff_coeff (coeffs, s, d, S, dickson_a);
  
  fd = (mpres_t *) malloc((S + 1) * sizeof(mpres_t));
  if (fd == NULL)
    {
      fprintf (stderr, "Error: not enough memory\n");
      exit (EXIT_FAILURE);
    }

  for (k = 0; k <= (unsigned) S; k++) 
    {
      mpres_init (fd[k], modulus);
      mpres_pow (fd[k], *x, coeffs[k], modulus);
    }

  clear_list (coeffs, S + 1);
      
  return fd;  
}

/* Frees all the dynamic variables allocated by pm1_rootsG_init() */

void 
pm1_rootsG_clear (mpres_t *fd, int S, mpmod_t modulus)
{
  unsigned int k;
  
  if (S > 6 && (S & 1) == 0)
    {
      S /= 2;
    }
  S = abs (S);
  
  for (k = 0; k <= (unsigned) S; k++)
    mpres_clear (fd[k], modulus);
  
  free (fd);
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
pm1_rootsG (mpz_t f, listz_t G, unsigned int d, mpres_t *fd, listz_t t, 
        int S, mpmod_t modulus, int verbose, unsigned long *tot_muls)
{
  unsigned int i, j, muls = 0;
  int st;
  int invtrick = 0, dickson_a = 0;
  
  st = cputime ();

  if (S > 6 && (S & 1) == 0)
    {
      invtrick = 1;
      S /= 2;
    }
  else if (S < 0)
    {
      dickson_a = -1;
      S = -S;
    }

  for (i = 0; i < d; i++)
    {
      mpres_get_z (G[i], fd[0], modulus);
      for (j = 0; j < (unsigned) S; j++)
        mpres_mul(fd[j], fd[j], fd[j+1], modulus);
    }
  muls += d * S;
  
  if (invtrick)
    {
      if (list_invert (t, G, d, t[d], modulus)) 
        {
          if (verbose >= 2)
            printf ("Found factor while inverting G[0]*..*G[d]\n");
          mpz_set (f, t[d]);
          return 1;
        }
      muls += 3 * d;
    
      for (i = 0; i < d; i++) 
        {
          mpz_add (G[i], G[i], t[i]);
          mpz_mod (G[i], G[i], modulus->orig_modulus);
        }
    }
  
  if (tot_muls != NULL)
    *tot_muls += muls;
  
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
pm1 (mpz_t f, mpz_t p, mpz_t N, mpz_t go, double B1done, double B1, double B2min, double B2,
     double B2scale, unsigned int k, int S, int verbose, int repr, mpz_t orig_X0)
{
  mpmod_t modulus;
  mpres_t x;
  int youpi = 0, st, base2, Nbits, smallbase;

  st = cputime ();
  
  /* Set default B2. See ecm.c for comments */
  if (B2 == 0.0)
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
      if (B1done == 1.0)
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

      if (B1done == 1.0 || verbose > 1) 
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
    youpi = pm1_stage1 (f, x, modulus, B1, B1done, verbose, N, orig_X0, go);

  if (verbose >= 1)
    {
      printf ("Step 1 took %dms\n", cputime() - st);
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

  youpi = stage2 (f, &x, modulus, B2min, B2, k, S, verbose, PM1_METHOD);

clear_and_exit:
  mpres_get_z (p, x, modulus);
  mpres_clear (x, modulus);
  mpmod_clear (modulus);

  return youpi;
}
