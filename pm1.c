/* Pollard 'P-1' algorithm.

  Copyright (C) 2001 Paul Zimmermann,
  LORIA/INRIA Lorraine, zimmerma@loria.fr
  See http://www.loria.fr/~zimmerma/records/ecmnet.html

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
#include "gmp-mparam.h"
#include "ecm.h"

#define DEBUG

typedef struct {
  unsigned int size;
  mpz_t val[1];
} mul_casc;

int      pm1_stage1     (mpz_t, mpres_t, mpmod_t, double, double);
mul_casc *mulcascade_init (void);
void     mulcascade_free (mul_casc *);
mul_casc *mulcascade_mul_ui (mul_casc *, unsigned int);
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
#define CASCADE_THRES 3

mul_casc *
mulcascade_init (void)
{
  mul_casc *t;
  t = (mul_casc *) malloc (sizeof (unsigned int) + sizeof (mpz_t));
  if (t == NULL)
    {
      fprintf (stderr, "mulcascade_init: could not allocate memory\n");
      exit (EXIT_FAILURE);
    }
  t->size = 1;
  mpz_init_set_ui (t->val[0], 1);
  return t;
}

void 
mulcascade_free (mul_casc *c) {
  unsigned int i;
  for (i = 0; i < c->size; i++)
    mpz_clear( c->val[i] );
  free(c);
}

/* TODO mulcascade_mul_d */
mul_casc * 
mulcascade_mul_ui (mul_casc *c, unsigned int n) {
  unsigned int i;

  if (mpz_cmp_ui (c->val[0], 1) == 0) /* TODO: Compare to 0 for speed? */
    {
      mpz_set_ui(c->val[0], n);
      return c;
    }

  mpz_mul_ui (c->val[0], c->val[0], n);
  if (mpz_size (c->val[0]) <= CASCADE_THRES)
    return c;
  
  for (i = 1; i < c->size; i++) 
    {
      if (mpz_cmp_ui (c->val[i], 1) == 0) 
        {
          mpz_set (c->val[i], c->val[i-1]);
          mpz_set_ui (c->val[i-1], 1);
          return c;
        } else {
          mpz_mul (c->val[i], c->val[i], c->val[i-1]);
          mpz_set_ui (c->val[i-1], 1);
        }
    }
  
  /* Allocate more space for cascade */
  
  i = c->size++;
  c = realloc (c, sizeof (unsigned int) + c->size * sizeof (mpz_t));
  mpz_init (c->val[i]);
  mpz_set (c->val[i], c->val[i-1]);
  mpz_set_ui (c->val[i-1], 1);
  return c;
}

void 
mulcascade_get_z (mpz_t r, mul_casc *c) {
  unsigned int i;
  
  if (c->size == 0)
    {
      mpz_set_ui (r, 1); /* Empty product */
      return;
    }

  mpz_set(r, c->val[0]);
  
  for (i = 1; i < c->size; i++)
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
pm1_stage1 (mpz_t f, mpres_t a, mpmod_t n, double B1, double B1done)
{
  double B0, p, q, r;
  mpz_t g, d;
  int youpi;
  unsigned int max_size, cascade_limit;
  unsigned int smallbase = 0;
  mul_casc *cascade;

  mpz_init (g);
  mpz_init (d);

  B0 = sqrt (B1);

  max_size = L1 * mpz_sizeinbase (n->orig_modulus, 2);

  mpres_get_z (g, a, n);
  if (mpz_fits_uint_p (g))
    {
      smallbase = mpz_get_ui (g);
    }

#ifdef EXPONENTIATE_BY_N_MINUS_1
  /* suggestion from Peter Montgomery: start with exponent n-1,
     since any prime divisor of b^m-1 which does not divide any
     algebraic factor of b^m-1 must be of the form km+1 [Williams82].
     Do this only when n is composite, otherwise all tests with prime
     n factor of a Cunningham number will succeed in stage 1. */
  if (mpz_probab_prime_p (n->orig_modulus, 1) == 0)
    {
      mpz_sub_ui (g, n->orig_modulus, 1);
      mpres_pow (a, a, g, n);
    }
  else
#endif
    mpz_set_ui (g, 1);

  /* TODO set dynamically depending on size of N */
  cascade_limit = 1000000;
  if (cascade_limit > B1)
    cascade_limit = B1;

  cascade = mulcascade_init ();

  if (B0 <= cascade_limit)
    {
      /* first loop through small primes <= sqrt(B1) */
      for (p = 2.0; p <= B0; p = getprime(p))
        {
          for (q = 1, r = p; r <= B1; r *= p)
            if (r > B1done) q *= p;
          cascade = mulcascade_mul_ui (cascade, q);
        }

      /* then all sqrt(B1) < primes < cascade_limit and taken with 
         exponent 1 */
      for ( ; p <= cascade_limit; p = getprime(p))
        if (p > B1done)
          cascade = mulcascade_mul_ui (cascade, p);
      
      mulcascade_get_z (g, cascade);
      mulcascade_free (cascade);
#ifdef DEBUG
      printf ("Exponent has %d bits\n", mpz_sizeinbase (g, 2));
#endif
      if (smallbase)
        {
#ifdef DEBUG
          printf ("Using mpres_ui_pow, base %u\n", smallbase);
#endif
          mpres_ui_pow (a, smallbase, g, n);
        }
      else
        {
          mpres_pow (a, a, g, n);
        }
      mpz_set_ui (g, 1);
    } else {
      for (p = 2.0; p <= cascade_limit; p = getprime(p))
        {
          for (q = 1, r = p; r <= B1; r *= p)
            if (r > B1done) q *= p;
          cascade = mulcascade_mul_ui (cascade, q);
        }
      
      mulcascade_get_z (g, cascade);
      mulcascade_free (cascade);
#ifdef DEBUG
      printf("Exponent has %d bits\n", mpz_sizeinbase (g, 2));
#endif
      if (smallbase)
        {
#ifdef DEBUG
          printf ("Using mpres_ui_pow, base %u\n", smallbase);
#endif
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
  
  /* then remaining primes > max(sqrt(B1), cascade_limit) and taken 
     with exponent 1 */
  for (; p <= B1; p = getprime(p))
    if (p > B1done)
      {
        mpz_mul_d (g, g, p, d);
        if (mpz_sizeinbase (g, 2) >= max_size)
	  {
	    mpres_pow (a, a, g, n);
	    mpz_set_ui (g, 1);
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
pm1_rootsF (mpz_t f, listz_t F, unsigned int d, mpres_t *x, listz_t t,
        unsigned int S, mpmod_t modulus, int verbose)
{
  unsigned int i, j, k;
  int st, st2;
  mpres_t *fd;
  listz_t coeffs;
  int invtrick = 0, dickson_a = -1;

  st = cputime ();

  mpres_get_z (F[0], *x, modulus); /* s^1 for P-1 */
  i = 1;

  if (S > 6 && (S & 1) == 0)
    { /* We can't use both invtrick and proper Dickson polys */
      invtrick = 1;
      dickson_a = 0;
      S /= 2;
    }

  if (d > 7)
    {
      st2 = cputime ();
      coeffs = init_list (S + 1);
      
      fin_diff_coeff (coeffs, 7, 6, S, dickson_a);
      
      fd = (mpres_t *) malloc ((S + 1) * sizeof (mpres_t));
      for (k = 0; k <= S; k++) 
        {
          mpres_init (fd[k], modulus);
          mpres_pow (fd[k], *x, coeffs[k], modulus);
        }

      clear_list (coeffs, S + 1);
      
      if (verbose >= 2)
        printf ("Initializing table of differences for F took %dms\n", cputime () - st2);

      for (j = 7; j < d; j += 6)
        {
          if (gcd (j, d) == 1)
            mpres_get_z (F[i++], fd[0], modulus);
          
          /* Compute value of f_{S, a}(7+n*6) for the next n */
          for (k = 0; k < S; k++)
            mpres_mul (fd[k], fd[k], fd[k+1], modulus);
        }
      
      for (k = 0; k <= S; k++)
        mpres_clear (fd[k], modulus);
      free (fd);
    }

  if (invtrick)
    {
      if (list_invert (t, F, i, t[i], modulus->orig_modulus)) 
        {
          mpz_set (f, t[i]);
          return 1;
        }
      
      for (j = 0; j < i; j++) 
        {
          mpz_add (F[j], F[j], t[j]);
          mpz_mod (F[j], F[j], modulus->orig_modulus);
        }
    }
  
  if (verbose >= 2)
    printf ("Computing roots of F took %dms\n", cputime () - st);

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
pm1_rootsG_init (mpres_t *x, unsigned int s, unsigned int d, unsigned int S,
                 mpmod_t modulus)
{
  unsigned int k;
  int invtrick = 0, dickson_a = -1;
  listz_t coeffs;
  mpres_t *fd;
  
  if (S > 6 && (S & 1) == 0)
    {
      invtrick = 1;
      dickson_a = 0;
      S /= 2;
    }
  
  coeffs = init_list (S + 1);

  fin_diff_coeff(coeffs, s, d, S, dickson_a);
  
  fd = (mpres_t *) malloc((S + 1) * sizeof(mpres_t));
  for (k = 0; k <= S; k++) 
    {
      mpres_init (fd[k], modulus);
      mpres_pow (fd[k], *x, coeffs[k], modulus);
    }

  clear_list (coeffs, S + 1);
      
  return fd;  
}

/* Frees all the dynamic variables allocated by pm1_rootsG_init() */

void 
pm1_rootsG_clear (mpres_t *fd, unsigned int S, mpmod_t modulus)
{
  unsigned int k;
  
  if (S > 6 && (S & 1) == 0)
    {
      S /= 2;
    }
  
  for (k = 0; k <= S; k++)
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
pm1_rootsG (mpz_t f, listz_t G, unsigned int d, mpres_t *fd, 
        listz_t t, unsigned int S, mpmod_t modulus, int verbose)
{
  unsigned int i, j;
  int st;
  int invtrick = 0;
  
  st = cputime ();

  if (S > 6 && (S & 1) == 0) 
    {
      invtrick = 1;
      S /= 2;
    }
  
  for (i = 0; i < d; i++)
    {
      mpres_get_z (G[i], fd[0], modulus);
      for (j = 0; j < S; j++)
        mpres_mul(fd[j], fd[j], fd[j+1], modulus);
    }
  
  if (invtrick)
    {
      if (list_invert (t, G, i, t[i], modulus->orig_modulus)) 
        {
          mpz_set (f, t[i]);
          return 1;
        }
    
      for (j = 0; j < i; j++) 
        {
          mpz_add (G[j], G[j], t[j]);
          mpz_mod (G[j], G[j], modulus->orig_modulus);
        }
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
pm1 (mpz_t f, mpz_t p, mpz_t N, double B1, double B2, double B1done,
     unsigned int k, unsigned int S, int verbose, int repr)
{
  mpmod_t modulus;
  mpres_t x;
  int youpi = 0, st, base2, Nbits, smallbase;

  st = cputime ();

  if (verbose >= 1)
    {
      printf ("Using seed=");
      mpz_out_str (stdout, 10, p);
      printf ("\n");
      fflush (stdout);
    }

  if (repr)
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
      base2 = isbase2 (N, 2.0);
      smallbase = mpz_fits_uint_p (p);
      
      /* TODO: make dependant on Nbits and base2 */
      if (base2)
        {
          printf ("Using base-2 representation: 2^%d%c1\n", 
                  abs (base2), (base2 > 0) ? '+' : '-');
          mpmod_init_BASE2 (modulus, base2, N);
        }

      else if (mpz_size (N) <= 2 * POWM_THRESHOLD && smallbase)
      /* Below POWM_THRESHOLD, mpz_powm uses MODMULN reduction, too, but 
         without special code for small bases which makes our MODMULN
         faster. Above POWM_THRESHOLD mpz_powm uses faster mod reduction,
         at about 2*POWM_THRESHOLD it catches up with our smallbase-MODMULN
         and then is faster until REDC takes over. */
        {
          printf ("Using MODMULN\n");
          mpmod_init_MODMULN (modulus, N);
        }
      else if (Nbits > 50000 ||  (Nbits > 3500 && smallbase))
        {
          printf ("Using REDC\n");
          mpmod_init_REDC (modulus, N);
        }
      else
        {
          printf ("Using mpz_powm\n");
          mpmod_init_MPZ (modulus, N);
        }
    }
  
  mpres_init (x, modulus);
  mpres_set_z (x, p, modulus);

  if (B1 > B1done)
    youpi = pm1_stage1 (f, x, modulus, B1, B1done);

  if (verbose >= 1)
    {
      printf ("Stage 1 took %dms\n", cputime() - st);
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

  /* We need Suyama's power even and at least 2 for stage 2 to work 
     correctly */
  if (S < 2)
    S = 2;

  if (S & 1)
    S++;

  youpi = stage2 (f, &x, modulus, B2, k, S, verbose, PM1_METHOD, B1);

clear_and_exit:
  mpres_get_z (p, x, modulus);
  mpres_clear (x, modulus);
  mpmod_clear (modulus);

  return youpi;
}
