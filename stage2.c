/* Common stage 2 for ECM, P-1 and P+1 (improved standard continuation
   with subquadratic polynomial arithmetic).

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
  the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
  MA 02111-1307, USA.
*/

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <math.h> /* for floor */
#include "sp.h"
#if defined (_MSC_VER)
#define snprintf _snprintf
#endif
#include "gmp.h"
#include "ecm.h"
#include "ecm-impl.h"

//#undef HAVE_NTT

extern unsigned int Fermat;

/* r <- Dickson(n,a)(x) */
static void 
dickson (mpz_t r, mpz_t x, unsigned int n, int a)
{
  unsigned int i, b = 0;
  mpz_t t, u;

  if (n == 0)
    {
      mpz_set_ui (r, 1);
      return;
    }
  
  while (n > 2 && (n & 1) == 0)
    {
      b++;
      n >>= 1;
    }
  
  mpz_set (r, x);
  
  mpz_init (t);
  mpz_init (u);

  if (n > 1)
    {
      mpz_set (r, x);
      mpz_mul (r, r, r);
      mpz_sub_si (r, r, a);
      mpz_sub_si (r, r, a); /* r = dickson(x, 2, a) */
      
      mpz_set (t, x);    /* t = dickson(x, 1, a) */
      
      for (i = 2; i < n; i++)
        {
          mpz_mul_si (u, t, a);
          mpz_set (t, r);     /* t = dickson(x, i, a) */
          mpz_mul (r, r, x);
          mpz_sub (r, r, u);  /* r = dickson(x, i+1, a) */
        }
    }
  
  for ( ; b > 0; b--)
    {
      mpz_mul (t, r, r); /* t = dickson(x, n, a) ^ 2 */
      mpz_ui_pow_ui (u, abs (a), n);
      if (n & 1 && a < 0)
        mpz_neg (u, u);
      mpz_mul_2exp (u, u, 1); /* u = 2 * a^n */
      mpz_sub (r, t, u); /* r = dickson(x, 2*n, a) */
      n <<= 1;
    }
  
  mpz_clear (t);
  mpz_clear (u);
}


/* Init table to allow computation of

   Dickson_{E, a} (s + n*D), 

   for successive n, where Dickson_{E, a} is the Dickson polynomial 
   of degree E with parameter a. For a == 0, Dickson_{E, a} (x) = x^E .

   See Knuth, TAOCP vol.2, 4.6.4 and exercise 7 in 4.6.4, and
   "An FFT Extension of the Elliptic Curve Method of Factorization",
   Peter Montgomery, Dissertation, 1992, Chapter 5.

   Ternary return value.
*/

static void
fin_diff_coeff (listz_t coeffs, mpz_t s, mpz_t D, unsigned int E, 
                int dickson_a)
{
  unsigned int i, k;
  mpz_t t;
  
  mpz_init_set (t, s);
  
  for (i = 0; i <= E; i++)
    {
      if (dickson_a != 0)         /* fd[i] = dickson_{E,a} (s+i*D) */
        dickson (coeffs[i], t, E, dickson_a); 
      else                        /* fd[i] = (s+i*D)^E */
        mpz_pow_ui (coeffs[i], t, E);
      mpz_add (t, t, D);          /* t = s + i * D */
    }
  
  for (k = 1; k <= E; k++)
    for (i = E; i >= k; i--)
      mpz_sub (coeffs[i], coeffs[i], coeffs[i-1]);
  
  mpz_clear (t);
}


/* Init several disjoint progressions for the computation of 

   Dickson_{E,a} (s + e * (i + d * n * k)), 0 <= i < k * d, gcd(s+e*i, d) == 1,
                                            i == 1 (mod m)
   
   for successive n. m must divide d, e must divide s (d need not).
   
   This means there will be k sets of progressions, where each set contains
   eulerphi(d) progressions that generate the values coprime to d and with
   i == 1 (mod m).
   
   Return NULL if an error occurred.
*/

listz_t
init_progression_coeffs (mpz_t s, unsigned int d, unsigned int e, 
                         unsigned int k, unsigned int m, unsigned int E, 
                         int dickson_a)
{
  unsigned int i, j, size_fd;
  mpz_t t, dke;
  listz_t fd;

  ASSERT (d % m == 0);
  ASSERT (mpz_fdiv_ui (s, e) == 0);

  size_fd = k * phi(d) / phi(m) * (E + 1);
  outputf (OUTPUT_TRACE, "init_progression_coeffs: s = %Zd, d = %u, e = %u, "
           "k = %u, m = %u, E = %u, a = %d, size_fd = %u\n", 
           s, d, e, k, m, E, dickson_a, size_fd);

  fd = (listz_t) malloc (size_fd * sizeof (mpz_t));
  if (fd == NULL)
    return NULL;
  for (i = 0; i < size_fd; i++)
    mpz_init (fd[i]);
  mpz_init (t);
  mpz_set_ui (t, e * (1 % m));
  mpz_add (t, t, s);
  
  /* dke = d * k * e */
  mpz_init (dke);
  mpz_set_ui (dke, d);
  mpz_mul_ui (dke, dke, k);
  mpz_mul_ui (dke, dke, e);
  
  for (i = 1 % m, j = 0; i < k * d; i += m)
    {
      if (mpz_gcd_ui (NULL, t, d) == 1)
        {
          outputf (OUTPUT_TRACE, "init_progression_coeffs: initing a "
                   "progression for Dickson_{%d,%d}(%Zd + n * %Zd)\n", 
                   E, dickson_a, t, dke);
          fin_diff_coeff (fd + j, t, dke, E, dickson_a);
          j += E + 1;
        } else
          if (test_verbose (OUTPUT_TRACE))
            outputf (OUTPUT_TRACE, "init_progression_coeffs: NOT initing a "
                     "progression for Dickson_{%d,%d}(%Zd + n * %Zd), "
                     "gcd (%Zd, %u) == %u)\n", E, dickson_a, t, dke, t, d,
                     mpz_gcd_ui (NULL, t, d));
      mpz_add_ui (t, t, e * m); /* t = s + i * e */
    }

  mpz_clear (dke);
  mpz_clear (t);
  return fd;
}

void 
init_roots_state (ecm_roots_state *state, int S, unsigned int d1, 
                  unsigned int d2, double cost)
{
  ASSERT (gcd (d1, d2) == 1);
  /* If S < 0, use degree |S| Dickson poly, otherwise use x^S */
  state->S = abs (S);
  state->dickson_a = (S < 0) ? -1 : 0;

  /* We only calculate Dickson_{S, a}(j * d2) * s where
     gcd (j, dsieve) == 1 and j == 1 (mod 6)
     by doing nr = eulerphi(dsieve / 6) separate progressions. */
  /* Now choose a value for dsieve. */
  state->dsieve = 6;
  state->nr = 1;

  /* Prospective saving by sieving out multiples of 5:
     d1 / state->dsieve * state->nr / 5 roots, each one costs S point adds
     Prospective cost increase:
     4 times as many progressions to init (that is, 3 * state->nr more),
     each costs ~ S * S * log_2(5 * dsieve * d2) / 2 point adds
     The state->nr and one S cancel.
  */
  if (d1 % 5 == 0 &&
      d1 / state->dsieve / 5. * cost > 
      3. * state->S * log (5. * state->dsieve * d2) / 2.)
    {
      state->dsieve *= 5;
      state->nr *= 4;
    }

  if (d1 % 7 == 0 &&
      d1 / state->dsieve / 7. * cost > 
      5. * state->S * log (7. * state->dsieve * d2) / 2.)
    {
      state->dsieve *= 7;
      state->nr *= 6;
    }

  if (d1 % 11 == 0 &&
      d1 / state->dsieve / 11. * cost > 
      9. * state->S * log (11. * state->dsieve * d2) / 2.)
    {
      state->dsieve *= 11;
      state->nr *= 10;
    }

  state->size_fd = state->nr * (state->S + 1);
  state->next = 0;
  state->rsieve = 1;
}

/* Input:  X is the point at end of stage 1
           n is the number to factor
           B2min-B2 is the stage 2 range (we consider B2min is done)
           k0 is the number of blocks (if 0, use default)
           S is the exponent for Brent-Suyama's extension
           invtrick is non-zero iff one uses x+1/x instead of x.
           method: ECM_ECM, ECM_PM1 or ECM_PP1
           Cf "Speeding the Pollard and Elliptic Curve Methods
               of Factorization", Peter Montgomery, Math. of Comp., 1987,
               page 257: using x^(i^e)+1/x^(i^e) instead of x^(i^(2e))
               reduces the cost of Brent-Suyama's extension from 2*e
               to e+3 multiplications per value of i.
   Output: f is the factor found
   Return value: 2 (step number) iff a factor was found,
                 or ECM_ERROR if an error occurred.
*/
int
stage2 (mpz_t f, void *X, mpmod_t modulus, mpz_t B2min, mpz_t B2,
        unsigned int k0, int S, int method, int stage1time, 
        char *TreeFilename)
{
  double b2;
  unsigned int k;
  unsigned int i, d, d2, dF, sizeT;
  mpz_t n, i0, s, effB2; /* s = i0 * d */
  listz_t F, G, H, T;
  int youpi = 0;
  unsigned int st, st0;
  void *rootsG_state = NULL;
  listz_t *Tree = NULL; /* stores the product tree for F */
  unsigned int lgk; /* ceil(log(k)/log(2)) */
  listz_t invF = NULL;
  double mem;

  /* check alloc. size of f */
  mpres_realloc (f, modulus);

  if (mpz_cmp(B2, B2min) < 0)
    return 0;

  st0 = cputime ();

  /* since we consider only residues = 1 mod 6 in intervals of length d
     (d multiple of 6), each interval [i*d,(i+1)*d] covers partially itself
     and [(i-1)*d,i*d]. Thus to cover [B2min, B2] with all intervals 
     [i*d,(i+1)*d] for i0 <= i < i1 , we should  have i0*d <= B2min and 
     B2 <= (i1-1)*d */
  d = d2 = dF = 0;
  Fermat = 0;
  k = k0;
  mpz_init_set (effB2, B2);
  mpz_init (i0);
  mpz_init (s);
  if (modulus->repr == 1 && modulus->bits > 0)
    {
      for (i = modulus->bits; (i & 1) == 0; i >>= 1);
      if (1 && i == 1)
        {
          Fermat = modulus->bits;
          outputf (OUTPUT_DEVVERBOSE, "Choosing power of 2 poly length "
                   "for 2^%d+1 (%d blocks)\n", Fermat, k0);
          if (bestD (B2min, effB2, 1, &d, &d2, &k, &dF, i0) == ECM_ERROR)
            {
              youpi = ECM_ERROR;
              goto clear_s_i0;
            }
        }
    }
  if (d == 0)
    {
#if 1 || defined HAVE_NTT
      /* choose power-of-two dF */
      if (bestD (B2min, effB2, 1, &d, &d2, &k, &dF, i0) == ECM_ERROR)
#else
      if (bestD (B2min, effB2, 0, &d, &d2, &k, &dF, i0) == ECM_ERROR)
#endif
        {
          youpi = ECM_ERROR;
          goto clear_s_i0;
        }
    }
  
  mpz_mul_ui (s, i0, d); /* s = i0 * d */
  b2 = (double) dF * (double) d * (double) d2 / (double) phi (d2);

  outputf (OUTPUT_VERBOSE, "B2'=%Zd k=%u b2=%1.0f d=%u d2=%u dF=%u, "
           "i0=%Zd\n", effB2, k, b2, d, d2, dF, i0);

  lgk = ceil_log2 (dF);
  mem = 9.0 + (double) lgk;
#if (MULT == KS)
  mem += 24.0; /* estimated memory for kronecker_schonhage */
  mem += 1.0;  /* for the wrap-case in PrerevertDivision   */
#endif
  mem *= (double) dF;
  mem *= (double) mpz_size (modulus->orig_modulus);
  mem *= (double) mp_bits_per_limb / 8.0;
  if (mem < 1e4)
    outputf (OUTPUT_VERBOSE, "Estimated memory usage: %1.0f\n", mem);
  else if (mem < 1e7)
    outputf (OUTPUT_VERBOSE, "Estimated memory usage: %1.0fK\n", mem / 1e3);
  else if (mem < 1e10)
    outputf (OUTPUT_VERBOSE, "Estimated memory usage: %1.0fM\n", mem / 1e6);
  else
    outputf (OUTPUT_VERBOSE, "Estimated memory usage: %1.0fG\n", mem / 1e9);

  if (method == ECM_ECM && test_verbose (OUTPUT_VERBOSE))
    {
      double nrcurves;
      rhoinit (256, 10);
      outputf (OUTPUT_VERBOSE, "Expected number of curves to find a factor "
               "of n digits:\n20\t25\t30\t35\t40\t45\t50\t55\t60\t65\n");
      for (i = 20; i <= 65; i += 5)
        {
          nrcurves = 1. / ecmprob (mpz_get_d (B2min), mpz_get_d (effB2), 
                                 pow (10., i - .5), (double) dF * dF * k, S); 
          if (nrcurves < 10000000)
            outputf (OUTPUT_VERBOSE, "%.0f%c", 
                     floor (nrcurves + .5), i < 65 ? '\t' : '\n');
          else
            outputf (OUTPUT_VERBOSE, "%.2g%c", 
                    floor (nrcurves + .5), i < 65 ? '\t' : '\n');
        }
    }
    
  F = init_list (dF + 1);
  if (F == NULL)
    {
      youpi = ECM_ERROR;
      goto clear_s_i0;
    }

  sizeT = 3 * dF + list_mul_mem (dF);
  if (dF > 3)
    sizeT += dF;
  T = init_list (sizeT);
  if (T == NULL)
    {
      youpi = ECM_ERROR;
      goto clear_F;
    }
  H = T;

  /* needs dF+1 cells in T */
  if (method == ECM_PM1)
    youpi = pm1_rootsF (f, F, d, d2, dF, (mpres_t*) X, T, S, modulus);
  else if (method == ECM_PP1)
    youpi = pp1_rootsF (F, d, d2, dF, (mpres_t*) X, T, S, modulus);
  else 
    youpi = ecm_rootsF (f, F, d, d2, dF, (curve*) X, S, modulus);

  if (youpi != ECM_NO_FACTOR_FOUND)
    {
      if (youpi != ECM_ERROR)
	youpi = 2;
      goto clear_T;
    }

  /* ----------------------------------------------
     |   F    |  invF  |   G   |         T        |
     ----------------------------------------------
     | rootsF |  ???   |  ???  |      ???         |
     ---------------------------------------------- */

  if (TreeFilename == NULL)
    {
      Tree = (listz_t*) malloc (lgk * sizeof (listz_t));
      if (Tree == NULL)
        {
          outputf (OUTPUT_ERROR, "Error: not enough memory\n");
          youpi = ECM_ERROR;
          goto clear_T;
        }
      for (i = 0; i < lgk; i++)
        {
          Tree[i] = init_list (dF);
          if (Tree[i] == NULL)
            {
              /* clear already allocated Tree[i] */
              while (i)
              clear_list (Tree[--i], dF);
              free (Tree);
              youpi = ECM_ERROR;
              goto clear_T;
            }
        }
      /* list_set (Tree[lgk - 1], F, dF); PolyFromRoots_Tree does it */
    }
  else
    Tree = NULL;
  
#ifdef TELLEGEN_DEBUG
  outputf (OUTPUT_ALWAYS, "Roots = ");
  print_list (os, F, dF);
#endif
  mpz_init_set (n, modulus->orig_modulus);
  st = cputime ();
  if (TreeFilename != NULL)
    {
      FILE *TreeFile;
      char fullname[256];
      for (i = lgk; i > 0; i--)
        {
#ifdef HAVE_snprintf
          snprintf (fullname, 256, "%.252s.%d", TreeFilename, i - 1);
#else
          sprintf (fullname, "%.252s.%d", TreeFilename, i - 1);
#endif
          TreeFile = fopen (fullname, "wb");
          if (TreeFile == NULL)
            {
              outputf (OUTPUT_ERROR, "Error opening file for product tree of F\n");
              youpi = ECM_ERROR;
              goto free_Tree_i;
            }
          if (PolyFromRoots_Tree (F, F, dF, T, i - 1, n, NULL, TreeFile, 0)
              == ECM_ERROR)
            {
              fclose (TreeFile);
              youpi = ECM_ERROR;
              goto free_Tree_i;
            };
          if (fclose (TreeFile) != 0)
            {
              youpi = ECM_ERROR;
              goto free_Tree_i;
            }
        }
    }
  else
#ifdef HAVE_NTT
    ntt_PolyFromRoots_Tree (F, F, dF, T, n, Tree, NULL);
#else
    PolyFromRoots_Tree (F, F, dF, T, -1, n, Tree, NULL, 0);
#endif
  
  outputf (OUTPUT_VERBOSE, "Building F from its roots took %ums\n", 
           elltime (st, cputime ()));

  /* needs dF+list_mul_mem(dF/2) cells in T */

  mpz_set_ui (F[dF], 1); /* the leading monic coefficient needs to be stored
                             explicitly for PrerevertDivision */

  /* ----------------------------------------------
     |   F    |  invF  |   G   |         T        |
     ----------------------------------------------
     |  F(x)  |  ???   |  ???  |      ???         |
     ---------------------------------------------- */

  /* G*H has degree 2*dF-2, hence we must cancel dF-1 coefficients
     to get degree dF-1 */
  if (dF > 1)
    {
      /* only dF-1 coefficients of 1/F are needed to reduce G*H,
         but we need one more for TUpTree */
      invF = init_list (dF + 1);
      if (invF == NULL)
	{
	  youpi = ECM_ERROR;
	  goto free_Tree_i;
	}
      st = cputime ();
#ifdef HAVE_NTT
      ntt_PolyInvert (invF, F + 1, dF, T, n);
#else
      PolyInvert (invF, F + 1, dF, T, n);
#endif
      /* now invF[0..dF-1] = Quo(x^(2dF-1), F) */
      outputf (OUTPUT_VERBOSE, "Computing 1/F took %ums\n",
	       elltime (st, cputime ()));
      
      /* ----------------------------------------------
         |   F    |  invF  |   G   |         T        |
         ----------------------------------------------
         |  F(x)  | 1/F(x) |  ???  |      ???         |
         ---------------------------------------------- */
    }

  /* start computing G with roots at i0*d, (i0+1)*d, (i0+2)*d, ... 
     where i0*d <= B2min < (i0+1)*d */
  G = init_list (dF);
  if (G == NULL)
    {
      youpi = ECM_ERROR;
      goto clear_invF;
    }

  st = cputime ();
  if (method == ECM_PM1)
    rootsG_state = pm1_rootsG_init ((mpres_t *) X, s, d, d2, S, modulus);
  else if (method == ECM_PP1)
    rootsG_state = pp1_rootsG_init ((mpres_t *) X, s, d, d2, S, modulus);
  else /* ECM_ECM */
    rootsG_state = ecm_rootsG_init (f, (curve *) X, s, d, d2, dF, k, S, 
                                    modulus);

  /* rootsG_state=NULL if an error occurred or (ecm only) a factor was found */
  if (rootsG_state == NULL)
    {
      /* ecm: f = -1 if an error occurred */
      youpi = (method == ECM_ECM && mpz_cmp_si (f, -1)) ? 2 : ECM_ERROR;
      goto clear_G;
    }

  if (method != ECM_ECM) /* ecm_rootsG_init prints itself */
    outputf (OUTPUT_VERBOSE, "Initializing table of differences for G "
             "took %ums\n", elltime (st, cputime ()));

  for (i = 0; i < k; i++)
    {
      /* needs dF+1 cells in T+dF */
      if (method == ECM_PM1)
	youpi = pm1_rootsG (f, G, dF, (pm1_roots_state *) rootsG_state, T + dF,
			    modulus);
      else if (method == ECM_PP1)
        youpi = pp1_rootsG (G, dF, (pp1_roots_state *) rootsG_state, modulus,
                            (mpres_t *) X);
      else
	youpi = ecm_rootsG (f, G, dF, (ecm_roots_state *) rootsG_state, 
			    modulus);

      ASSERT(youpi != ECM_ERROR); /* xxx_rootsG cannot fail */
      if (youpi) /* factor found */
        {
          youpi = 2;
          goto clear_fd;
        }

  /* -----------------------------------------------
     |   F    |  invF  |   G    |         T        |
     -----------------------------------------------
     |  F(x)  | 1/F(x) | rootsG |      ???         |
     ----------------------------------------------- */

      st = cputime ();
#ifdef HAVE_NTT
      ntt_PolyFromRoots (G, G, dF, T + dF, n);
#else
      PolyFromRoots (G, G, dF, T + dF, n);
#endif

      /* needs 2*dF+list_mul_mem(dF/2) cells in T */
      outputf (OUTPUT_VERBOSE, "Building G from its roots took %ums\n", 
               elltime (st, cputime ()));

  /* -----------------------------------------------
     |   F    |  invF  |   G    |         T        |
     -----------------------------------------------
     |  F(x)  | 1/F(x) |  G(x)  |      ???         |
     ----------------------------------------------- */

      if (i == 0)
        {
          list_sub (H, G, F, dF); /* coefficients 1 of degree cancel,
                                     thus T is of degree < dF */
          list_mod (H, H, dF, n);
          /* ------------------------------------------------
             |   F    |  invF  |    G    |         T        |
             ------------------------------------------------
             |  F(x)  | 1/F(x) |  ???    |G(x)-F(x)|  ???   |
             ------------------------------------------------ */
        }
      else
	{
          /* since F and G are monic of same degree, G mod F = G - F */
          list_sub (G, G, F, dF);
          list_mod (G, G, dF, n);

          /* ------------------------------------------------
             |   F    |  invF  |    G    |         T        |
             ------------------------------------------------
             |  F(x)  | 1/F(x) |G(x)-F(x)|  H(x)  |         |
             ------------------------------------------------ */

	  st = cputime ();
	  /* previous G mod F is in H, with degree < dF, i.e. dF coefficients:
	     requires 3dF-1+list_mul_mem(dF) cells in T */
#ifdef HAVE_NTT
	  ntt_mul (T + dF, G, H, dF, T + 3 * dF, 0, n);
	  list_mod (H, T + dF, 2 * dF, n);
#else
	  list_mulmod (H, T + dF, G, H, dF, T + 3 * dF, n);
#endif

          outputf (OUTPUT_VERBOSE, "Computing G * H took %ums\n", 
                   elltime (st, cputime ()));

          /* ------------------------------------------------
             |   F    |  invF  |    G    |         T        |
             ------------------------------------------------
             |  F(x)  | 1/F(x) |G(x)-F(x)| G * H  |         |
             ------------------------------------------------ */

	  st = cputime ();
#ifdef HAVE_NTT
	  ntt_PrerevertDivision (H, F, invF + 1, dF, T + 2 * dF, n);
#else
	  if (PrerevertDivision (H, F, invF + 1, dF, T + 2 * dF, n))
	    {
	      youpi = ECM_ERROR;
	      goto clear_fd;
	    }
#endif
          outputf (OUTPUT_VERBOSE, "Reducing  G * H mod F took %ums\n", 
                   elltime (st, cputime ()));
	}
    }
  
  clear_list (F, dF + 1);
  F = NULL;
  clear_list (G, dF);
  G = NULL;
  st = cputime ();
#ifdef POLYEVALTELLEGEN
  youpi = polyeval_tellegen (T, dF, Tree, T + dF + 1, sizeT - dF - 1, invF,
			     n, TreeFilename);
  if (youpi)
    {
      outputf (OUTPUT_ERROR, "Error, not enough memory\n");
      goto clear_fd;
    }
#else
  clear_list (invF, dF + 1);
  invF = NULL;
  polyeval (T, dF, Tree, T + dF + 1, n, 0, ECM_STDERR);
#endif

  outputf (OUTPUT_VERBOSE, "Computing polyeval(F,G) took %ums\n", 
           elltime (st, cputime ()));

  youpi = list_gcd (f, T, dF, n) ? 2 : 0;
  outputf (OUTPUT_RESVERBOSE, "Product of G(f_i) = %Zd\n", T[0]);

 clear_fd:
  if (method == ECM_PM1)
    pm1_rootsG_clear ((pm1_roots_state *) rootsG_state, modulus);
  else if (method == ECM_PP1)
    pp1_rootsG_clear ((pp1_roots_state *) rootsG_state, modulus);
  else /* ECM_ECM */
    ecm_rootsG_clear ((ecm_roots_state *) rootsG_state, S, modulus);

clear_G:
  clear_list (G, dF);
 clear_invF:
  clear_list (invF, dF + 1);

free_Tree_i:
  if (Tree != NULL)
    {
      for (i = 0; i < lgk; i++)
        clear_list (Tree[i], dF);
      free (Tree);
    }

  mpz_clear (n);

clear_T:
  clear_list (T, sizeT);
clear_F:
  clear_list (F, dF + 1);

clear_s_i0:
  mpz_clear (i0);
  mpz_clear (s);

  st0 = elltime (st0, cputime ());

  outputf (OUTPUT_NORMAL, "Step 2 took %ums\n", st0);

  if (method == ECM_ECM && test_verbose (OUTPUT_VERBOSE) && youpi != ECM_ERROR)
    {
      double prob, tottime, exptime;
      rhoinit (256, 10);
      outputf (OUTPUT_VERBOSE, "Expected time to find a factor of n digits:\n"
	       "20\t25\t30\t35\t40\t45\t50\t55\t60\t65\n");
      tottime = (double) stage1time + (double) st0;
      for (i = 20; i <= 65; i += 5)
        {
          const char sep = (i < 65) ? '\t' : '\n';
          prob = ecmprob (mpz_get_d (B2min), mpz_get_d (effB2), 
                          pow (10., i - .5), (double) dF * dF * k, S);
          exptime = tottime / prob;
          /* outputf (OUTPUT_VERBOSE, "Total time: %.0f, probability: %f, expected time: %.0f\n", tottime, prob, exptime); */
          if (exptime < 1000.)
            outputf (OUTPUT_VERBOSE, "%.0fms%c", exptime, sep);
          else if (exptime < 60000.) /* One minute */
            outputf (OUTPUT_VERBOSE, "%.2fs%c", exptime / 1000., sep);
          else if (exptime < 3600000.) /* One hour */
            outputf (OUTPUT_VERBOSE, "%.2fm%c", exptime / 60000., sep);
          else if (exptime < 86400000.) /* One day */
            outputf (OUTPUT_VERBOSE, "%.2fh%c", exptime / 3600000., sep);
          else if (exptime < 31536000000.) /* One year */
            outputf (OUTPUT_VERBOSE, "%.2fd%c", exptime / 86400000., sep);
          else if (exptime < 31536000000000.) /* One thousand years */
            outputf (OUTPUT_VERBOSE, "%.2fy%c", exptime / 31536000000., sep);
          else if (exptime < 31536000000000000.) /* One million years */
            outputf (OUTPUT_VERBOSE, "%.0fy%c", exptime / 31536000000., sep);
          else if (prob > 0.)
            outputf (OUTPUT_VERBOSE, "%.1gy%c", exptime / 31536000000., sep);
          else 
            outputf (OUTPUT_VERBOSE, "Inf%c", sep);
        }
      rhoinit (1, 0); /* Free memory of rhotable */
    }
  mpz_clear (effB2);

  return youpi;
}
