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
#include "gmp.h"
#include "ecm.h"
#include "ecm-impl.h"

extern unsigned int Fermat;

/* r <- x^n */
static void
mpz_d_pow_ui (mpz_t r, double x, unsigned long int n)
{
  mpz_set_d (r, x);
  mpz_pow_ui (r, r, n);
}

/* r <- Dickson(n,a)(x) */
static void 
dickson_ui (mpz_t r, double x, unsigned int n, int a)
{
  unsigned int i, b = 0;
  mpz_t t, u, v;

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
  
  mpz_set_d (r, x);
  
  mpz_init (t);
  mpz_init (u);
  mpz_init (v);

  mpz_set_d (v, x);

  if (n > 1)
    {
      mpz_set (r, v);
      mpz_mul (r, r, r);
      mpz_sub_si (r, r, a);
      mpz_sub_si (r, r, a); /* r = dickson(x, 2, a) */
      
      mpz_set (t, v);    /* t = dickson(x, 1, a) */
      
      for (i = 2; i < n; i++)
        {
          mpz_mul_si (u, t, a);
          mpz_set (t, r);     /* t = dickson(x, i, a) */
          mpz_mul (r, r, v);
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
  mpz_clear (v);
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

static int
fin_diff_coeff (listz_t coeffs, double s, double D,
                unsigned int E, int dickson_a)
{
  unsigned int i, k;

  /* check maximal value of s + i * D does not overflow */
  if (s + (double) E * D > TWO53) /* 2^53 */
    {
      outputf (OUTPUT_ERROR, "Error, overflow in fin_diff_coeff\n");
      outputf (OUTPUT_ERROR, "Please use a smaller B1 or B2min\n");
      return ECM_ERROR;
    }
  for (i = 0; i <= E; i++)
    if (dickson_a != 0)         /* fd[i] = dickson_{E,a} (s+i*D) */
      dickson_ui (coeffs[i], s + (double) i * D, E, dickson_a); 
    else                        /* fd[i] = (s+i*D)^E */
      mpz_d_pow_ui (coeffs[i], s + (double) i * D, E);
  
  for (k = 1; k <= E; k++)
    for (i = E; i >= k; i--)
      mpz_sub (coeffs[i], coeffs[i], coeffs[i-1]);

  return ECM_NO_FACTOR_FOUND;
}


/* Init several disjoint progressions for the computation of 

   Dickson_{E,a} (s + e * (i + d * n * k)), 0 <= i < k * d, gcd(i, d) == 1,
                                            i == 1 (mod m)
   
   for successive n. m must divide d.
   
   This means there will be k sets of progressions, where each set contains
   eulerphi(d) progressions that generate the values coprime to d and 
   == 1 (mod m).
   
   Return NULL if an error occurred.
*/

listz_t
init_progression_coeffs (double s, unsigned int d, unsigned int e, 
                         unsigned int k, unsigned int m, unsigned int E, 
                         int dickson_a)
{
  unsigned int i, j, size_fd, smd;
  listz_t fd;
  const double de = (double) e;

  ASSERT (d % m == 0);

  size_fd = k * phi(d) / phi(m) * (E + 1);
  fd = (listz_t) malloc (size_fd * sizeof (mpz_t));
  if (fd == NULL)
    return NULL;
  for (i = 0; i < size_fd; i++)
    mpz_init (fd[i]);

  /* smd := s % d */
  smd = (unsigned int) (s - floor (s / (double) d) * (double) d);

  j = 0;
  for (i = 1 % m; i < k * d; i += m)
    {
      if (gcd (smd + e * (i % d), d) == 1)
        {
          if (fin_diff_coeff (fd + j, s + de * i, de * k * d, E, dickson_a)
	      == ECM_ERROR)
	    {
	      for (i = 0; i < size_fd; i++)
		mpz_clear (fd[i]);
	      free (fd);
	      return NULL;
	    }
          j += E + 1;
        }
    }

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
     The state->dsieve and one S cancel.
  */
  if (d1 % 5 == 0 &&
      d1 / state->dsieve / 5 * cost > 
      3. * state->S * log (5. * state->dsieve * d2) / 2.)
    {
      state->dsieve *= 5;
      state->nr *= 4;
    }

  if (d1 % 7 == 0 &&
      d1 / state->dsieve / 7 * cost > 
      5. * state->S * log (7. * state->dsieve * d2) / 2.)
    {
      state->dsieve *= 7;
      state->nr *= 6;
    }

  if (d1 % 11 == 0 &&
      d1 / state->dsieve / 11 * cost > 
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
stage2 (mpz_t f, void *X, mpmod_t modulus, double B2min, double B2,
        unsigned int k0, int S, int method, int stage1time, 
        char *TreeFilename)
{
  double b2, i0;
  unsigned int k;
  unsigned int i, d, d2, dF, sizeT;
  mpz_t n;
  listz_t F, G, H, T;
  int youpi = 0, st, st0;
  void *rootsG_state = NULL;
  listz_t *Tree = NULL; /* stores the product tree for F */
  unsigned int lgk; /* ceil(log(k)/log(2)) */
  listz_t invF = NULL;

  /* check alloc. size of f */
  mpres_realloc (f, modulus);

  if (B2 < B2min)
    return 0;

  st0 = cputime ();

  if (k0 == 0)
    {
      outputf (OUTPUT_ERROR, 
               "Error: number of blocks in step 2 should be positive\n");
      return ECM_ERROR;
    }

  /* since we consider only residues = 1 mod 6 in intervals of length d
     (d multiple of 6), each interval [i*d,(i+1)*d] covers partially itself
     and [(i-1)*d,i*d]. Thus to cover [B2min, B2] with all intervals 
     [i*d,(i+1)*d] for i0 <= i < i1 , we should  have i0*d <= B2min and 
     B2 <= (i1-1)*d */
  d = d2 = dF = 0;
  Fermat = 0;
  k = k0;
  if (modulus->repr == 1 && modulus->bits > 0)
    {
      for (i = modulus->bits; (i & 1) == 0; i >>= 1);
      if (1 && i == 1)
        {
          Fermat = modulus->bits;
          outputf (OUTPUT_DEVVERBOSE, "Choosing power of 2 poly length "
                   "for 2^%d+1 (%d blocks)\n", Fermat, k0);
          if (bestD_po2 (B2min, B2, &d, &d2, &k, &i0) == ECM_ERROR)
            return ECM_ERROR;
          dF = 1 << ceil_log2 (phi (d) / 2);
        }
    }
  if (d == 0)
    {
      if (bestD (B2min, B2, &d, &d2, &k, &i0) == ECM_ERROR)
        return ECM_ERROR;
      dF = phi (d) / 2;
    }
  
  /* check that i0 * d does not overflow */
  if (i0 * (double) d > TWO53) /* 2^53 */
    {
      outputf (OUTPUT_ERROR, "Error, overflow in stage 2\n"
               "Please use a smaller B1 or B2min\n");
      return ECM_ERROR;
    }

  b2 = (double) dF * (double) d * (double) d2 / (double) phi (d2);

  /* compute real B2 */
  B2 = i0 * (double) d + floor ((double) k * b2 / d / d2) * d * d2;

  outputf (OUTPUT_VERBOSE, "B2'=%1.0f k=%u b2=%1.0f d=%u d2=%u dF=%u, "
           "i0=%.0f\n", B2, k, b2, d, d2, dF, i0);

  if (method == ECM_ECM && test_verbose (OUTPUT_VERBOSE))
    {
      double nrcurves;
      rhoinit (256, 10);
      outputf (OUTPUT_VERBOSE, "Expected number of curves to find a factor "
               "of n digits:\n20\t25\t30\t35\t40\t45\t50\t55\t60\t65\n");
      for (i = 20; i <= 65; i += 5)
        {
          nrcurves = 1. / ecmprob (B2min, B2, pow (10., i - .5), 
                                 (double)dF * (double)dF * k, S); 
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
      goto clear_n;
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

  lgk = ceil_log2 (dF);
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
              youpi = ECM_ERROR;
              goto free_Tree;
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
          snprintf (fullname, 256, "%.252s.%d", TreeFilename, i - 1);
          TreeFile = fopen (fullname, "wb");
          if (TreeFile == NULL)
            {
              outputf (OUTPUT_ERROR, "Error opening file for product tree of F\n");
              youpi = ECM_ERROR;
              goto clear_T;
            }
          if (PolyFromRoots_Tree (F, F, dF, T, i - 1, n, NULL, TreeFile, 0)
              == ECM_ERROR)
            {
              fclose (TreeFile);
              youpi = ECM_ERROR;
              goto clear_T;
            };
          if (fclose (TreeFile) != 0)
            {
              youpi = ECM_ERROR;
              goto clear_T;
            }
        }
    }
  else
    PolyFromRoots_Tree (F, F, dF, T, -1, n, Tree, NULL, 0);
  
  outputf (OUTPUT_VERBOSE, "Building F from its roots took %ums\n", 
           cputime() - st);

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
      PolyInvert (invF, F + 1, dF, T, n);

      /* now invF[0..dF-1] = Quo(x^(2dF-1), F) */
      outputf (OUTPUT_VERBOSE, "Computing 1/F took %ums\n", cputime () - st);
      
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
    rootsG_state = pm1_rootsG_init ((mpres_t *) X, i0 * (double) d, d, d2, S,
				    modulus);
  else if (method == ECM_PP1)
    rootsG_state = pp1_rootsG_init ((mpres_t *) X, i0 * (double) d, d, d2, S,
				    modulus);
  else /* ECM_ECM */
    rootsG_state = ecm_rootsG_init (f, (curve *) X, i0 * (double) d, d, d2,
				    dF, k, S, modulus);

  /* rootsG_state=NULL if an error occurred or (ecm only) a factor was found */
  if (rootsG_state == NULL)
    {
      /* ecm: f = -1 if an error occurred */
      youpi = (method == ECM_ECM && mpz_cmp_si (f, -1)) ? 2 : ECM_ERROR;
      goto clear_G;
    }

  if (method != ECM_ECM) /* ecm_rootsG_init prints itself */
    outputf (OUTPUT_VERBOSE, "Initializing table of differences for G "
             "took %dms\n", cputime () - st);

  for (i = 0; i < k; i++)
    {
      /* needs dF+1 cells in T+dF */
      if (method == ECM_PM1)
	youpi = pm1_rootsG (f, G, dF, (pm1_roots_state *) rootsG_state, T + dF,
			    modulus);
      else if (method == ECM_PP1)
        youpi = pp1_rootsG (G, dF, (pp1_roots_state *) rootsG_state, modulus);
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
      PolyFromRoots (G, G, dF, T + dF, n);
      /* needs 2*dF+list_mul_mem(dF/2) cells in T */
      outputf (OUTPUT_VERBOSE, "Building G from its roots took %ums\n", 
               cputime() - st);

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
	  list_mulmod (H, T + dF, G, H, dF, T + 3 * dF, n);

          outputf (OUTPUT_VERBOSE, "Computing G * H took %ums\n", 
                   cputime () - st);

          /* ------------------------------------------------
             |   F    |  invF  |    G    |         T        |
             ------------------------------------------------
             |  F(x)  | 1/F(x) |G(x)-F(x)| G * H  |         |
             ------------------------------------------------ */

	  st = cputime ();
          if (PrerevertDivision (H, F, invF + 1, dF, T + 2 * dF, n))
	    {
	      youpi = ECM_ERROR;
	      goto clear_fd;
	    }

          outputf (OUTPUT_VERBOSE, "Reducing  G * H mod F took %ums\n", 
                   cputime () - st);
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
           cputime () - st);

  youpi = list_gcd (f, T, dF, n) ? 2 : 0;
  outputf (OUTPUT_DEVVERBOSE, "Product of G(f_i) = %Zd\n", T[0]);

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
  if (TreeFilename == NULL)
    for (i = 0; i < lgk; i++)
      clear_list (Tree[i], dF);
free_Tree:
  if (TreeFilename == NULL)
    free (Tree);

clear_T:
  clear_list (T, sizeT);
clear_F:
  clear_list (F, dF + 1);

  st0 = cputime () - st0;

  outputf (OUTPUT_NORMAL, "Step 2 took %dms\n", st0);

  if (method == ECM_ECM && test_verbose (OUTPUT_VERBOSE))
    {
      double prob, tottime, exptime;
      rhoinit (256, 10);
      outputf (OUTPUT_VERBOSE, "Expected time to find a factor of n digits:\n"
	       "20\t25\t30\t35\t40\t45\t50\t55\t60\t65\n");
      tottime = (double) stage1time + (double) st0;
      for (i = 20; i <= 65; i += 5)
        {
          const char sep = (i < 65) ? '\t' : '\n';
          prob = ecmprob (B2min, B2, pow (10., i - .5), 
                          (double)dF * (double)dF * k, S);
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
    }

 clear_n:
  mpz_clear (n);

  return youpi;
}
