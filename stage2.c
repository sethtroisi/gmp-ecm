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
  the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
  MA 02110-1301, USA.
*/

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h> /* for floor */
#include <string.h> /* for strlen */

#if HAVE_UNISTD_H
#include <unistd.h> /* for unlink */
#endif

#include "ecm-impl.h"
#include "sp.h"

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
  
  MPZ_INIT (t);
  MPZ_INIT (u);

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
  
  MPZ_INIT (t);
  mpz_set (t, s);
  
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

   Dickson_{E,a} (e * (i0 + i + d * n * k)), for 0 <= i < k * d 
                  with gcd(e * (i0 + i), d) == 1, i == 1 (mod m)
   
   for successive n. m must divide d, e must divide s (d need not).
   
   This means there will be k sets of progressions, where each set contains
   eulerphi(d) progressions that generate the values coprime to d and with
   i == 1 (mod m).
   
   Return NULL if an error occurred.
*/

listz_t
init_progression_coeffs (mpz_t i0, const unsigned long d, const unsigned long e, 
                         const unsigned int k, const unsigned int m, 
                         const unsigned int E, const int dickson_a)
{
  unsigned int i, j, size_fd;
  mpz_t t, dke;
  listz_t fd;

  ASSERT (d % m == 0);

  size_fd = k * eulerphi(d) / eulerphi(m) * (E + 1);
  fd = (listz_t) malloc (size_fd * sizeof (mpz_t));
  if (fd == NULL)
    return NULL;
  for (i = 0; i < size_fd; i++)
    MPZ_INIT (fd[i]);

  MPZ_INIT (t);
  outputf (OUTPUT_TRACE, "init_progression_coeffs: i0 = %Zd, d = %u, e = %u, "
           "k = %u, m = %u, E = %u, a = %d, size_fd = %u\n", 
           i0 ? i0 : t, d, e, k, m, E, dickson_a, size_fd);

  if (i0 != NULL)
    mpz_set (t, i0);
  mpz_mul_ui (t, t, e);
  mpz_add_ui (t, t, e * (1 % m));
  
  /* dke = d * k * e */
  MPZ_INIT (dke);
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
init_roots_state (ecm_roots_state *state, const int S, const unsigned long d1, 
                  const unsigned long d2, const double cost)
{
  ASSERT (gcd (d1, d2) == 1);
  /* If S < 0, use degree |S| Dickson poly, otherwise use x^S */
  state->S = abs (S);
  state->dickson_a = (S < 0) ? -1 : 0;

  /* We only calculate Dickson_{S, a}(j * d2) * s where
     gcd (j, dsieve) == 1 and j == 1 (mod 6)
     by doing nr = eulerphi(dsieve)/2 separate progressions. */
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

double 
memory_use (unsigned long dF, unsigned int sp_num, unsigned int Ftreelvl,
            mpmod_t modulus)
{
  double mem;
  
  /* printf ("memory_use (%lu, %d, %d, )\n", dF, sp_num, Ftreelvl); */

  mem = 9.0 * 1.5; /* Some of the T[] entries contain unreduced products. */
  mem += (double) Ftreelvl;
  mem *= (double) dF;
#if (MULT == KS)
   /* estimated memory for kronecker_schonhage /
      wrap-case in PrerevertDivision respectively */
  mem += (24.0 + 1.0) * (double) (sp_num ? MUL_NTT_THRESHOLD : dF);
#endif
  mem *= (double) (mpz_size (modulus->orig_modulus) + 3)
                  * (mp_bits_per_limb / 8.0)
         + sizeof (mpz_t);
  
  if (sp_num)
    mem += /* peak malloc in ecm_ntt.c */
         (4.0 * dF * sp_num * sizeof (sp_t))
	 
	 /* mpzspv_normalise */
	 + (MPZSPV_NORMALISE_STRIDE * ((double) sp_num * 
	 	sizeof (sp_t) + 6.0 * sizeof (sp_t) + sizeof (float)))

	 /* sp_F, sp_invF */
	 + ((1.0 + 2.0) * dF * sp_num * sizeof (sp_t));

  return mem;
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
stage2 (mpz_t f, void *X, mpmod_t modulus, unsigned long dF, unsigned long k, 
        root_params_t *root_params, int method, int use_ntt, 
        char *TreeFilename, int (*stop_asap)(void))
{
  unsigned long i, sizeT;
  mpz_t n;
  listz_t F, G, H, T;
  int youpi = ECM_NO_FACTOR_FOUND;
  long st, st0;
  void *rootsG_state = NULL;
  listz_t *Tree = NULL; /* stores the product tree for F */
  unsigned int treefiles_used = 0; /* Number of tree files currently in use */
  unsigned int lgk; /* ceil(log(k)/log(2)) */
  listz_t invF = NULL;
  double mem;
  mpzspm_t mpzspm = NULL;
  mpzspv_t sp_F = NULL, sp_invF = NULL;
  
  /* check alloc. size of f */
  mpres_realloc (f, modulus);

  st0 = cputime ();

  Fermat = 0;
  if (modulus->repr == ECM_MOD_BASE2 && modulus->Fermat > 0)
    {
      Fermat = modulus->Fermat;
      use_ntt = 0; /* don't use NTT for Fermat numbers */
    }

  if (use_ntt)
    {
      mpzspm = mpzspm_init (2 * dF, modulus->orig_modulus);
  
      if (mpzspm == NULL)
        return ECM_ERROR;

      outputf (OUTPUT_VERBOSE,
	  "Using %u small primes for NTT\n", mpzspm->sp_num);
    }

  lgk = ceil_log2 (dF);

  mem = memory_use (dF, use_ntt ? mpzspm->sp_num : 0,
      (TreeFilename == NULL) ? lgk : 0, modulus);

  if (mem < 1e4)
    outputf (OUTPUT_VERBOSE, "Estimated memory usage: %1.0f\n", mem);
  else if (mem < 1e7)
    outputf (OUTPUT_VERBOSE, "Estimated memory usage: %1.0fK\n", mem / 1024.);
  else if (mem < 1e10)
    outputf (OUTPUT_VERBOSE, "Estimated memory usage: %1.0fM\n", 
             mem / 1048576.);
  else
    outputf (OUTPUT_VERBOSE, "Estimated memory usage: %1.0fG\n", 
             mem / 1073741824.);

  MEMORY_TAG;
  F = init_list2 (dF + 1, mpz_sizeinbase (modulus->orig_modulus, 2) + 
                          3 * mp_bits_per_limb);
  MEMORY_UNTAG;
  if (F == NULL)
    {
      youpi = ECM_ERROR;
      goto clear_i0;
    }

  sizeT = 3 * dF + list_mul_mem (dF);
  if (dF > 3)
    sizeT += dF;
  MEMORY_TAG;
  T = init_list2 (sizeT, 2 * mpz_sizeinbase (modulus->orig_modulus, 2) + 
                         3 * mp_bits_per_limb);
  MEMORY_UNTAG;
  if (T == NULL)
    {
      youpi = ECM_ERROR;
      goto clear_F;
    }
  H = T;

  /* needs dF+1 cells in T */
  if (method == ECM_PM1)
    youpi = pm1_rootsF (f, F, root_params, dF, (mpres_t*) X, T, modulus);
  else if (method == ECM_PP1)
    youpi = pp1_rootsF (F, root_params, dF, (mpres_t*) X, T, modulus);
  else 
    youpi = ecm_rootsF (f, F, root_params, dF, (curve*) X, modulus);

  if (youpi != ECM_NO_FACTOR_FOUND)
    {
      if (youpi != ECM_ERROR)
	youpi = ECM_FACTOR_FOUND_STEP2;
      goto clear_T;
    }
  if (stop_asap != NULL && (*stop_asap)())
    goto clear_T;


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
          MEMORY_TAG;
          Tree[i] = init_list2 (dF, mpz_sizeinbase (modulus->orig_modulus, 2) 
                                    + mp_bits_per_limb);
          MEMORY_UNTAG;
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
      /* assume this "small" malloc will not fail in normal usage */
      char *fullname = (char *) malloc (strlen (TreeFilename) + 1 + 2 + 1);
      
      for (i = lgk; i > 0; i--)
        {
          if (stop_asap != NULL && (*stop_asap)())
            goto free_Tree_i;
          sprintf (fullname, "%s.%lu", TreeFilename, i - 1);
          
	  TreeFile = fopen (fullname, "wb");
          if (TreeFile == NULL)
            {
              outputf (OUTPUT_ERROR, 
                       "Error opening file for product tree of F\n");
              youpi = ECM_ERROR;
              goto free_Tree_i;
            }
	  treefiles_used++;
	  
	  if (use_ntt)
	    {
	      if (ntt_PolyFromRoots_Tree (F, F, dF, T, i - 1, mpzspm, NULL,
		    TreeFile) == ECM_ERROR)
                {
                  fclose (TreeFile);
                  youpi = ECM_ERROR;
                  goto free_Tree_i;
		}
            }
	  else
	    {
	      if (PolyFromRoots_Tree (F, F, dF, T, i - 1, n, NULL,
		    TreeFile, 0) == ECM_ERROR)
   	        {
                  fclose (TreeFile);
                  youpi = ECM_ERROR;
                  goto free_Tree_i;
		}
	    }

          if (fclose (TreeFile) != 0)
            {
              youpi = ECM_ERROR;
              goto free_Tree_i;
            }
        }
      free (fullname);
    }
  else
    {
      /* TODO: how to check for stop_asap() here? */
      if (use_ntt)
        ntt_PolyFromRoots_Tree (F, F, dF, T, -1, mpzspm, Tree, NULL);
      else
	PolyFromRoots_Tree (F, F, dF, T, -1, n, Tree, NULL, 0);
    }
      
  outputf (OUTPUT_VERBOSE, "Building F from its roots took %ldms\n", 
           elltime (st, cputime ()));

  if (stop_asap != NULL && (*stop_asap)())
    goto free_Tree_i;


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
      MEMORY_TAG;
      invF = init_list2 (dF + 1, mpz_sizeinbase (modulus->orig_modulus, 2) + 
                                 2 * mp_bits_per_limb);
      MEMORY_UNTAG;
      if (invF == NULL)
	{
	  youpi = ECM_ERROR;
	  goto free_Tree_i;
	}
      st = cputime ();
      
      if (use_ntt)
        {
	  sp_F = mpzspv_init (dF, mpzspm);
	  mpzspv_from_mpzv (sp_F, 0, F, dF, mpzspm);
	  mpzspv_to_ntt (sp_F, 0, dF, dF, 1, mpzspm);
	  
	  ntt_PolyInvert (invF, F + 1, dF, T, mpzspm);
	  sp_invF = mpzspv_init (2 * dF, mpzspm);
	  mpzspv_from_mpzv (sp_invF, 0, invF, dF, mpzspm);
	  mpzspv_to_ntt (sp_invF, 0, dF, 2 * dF, 0, mpzspm);
	}
      else
        PolyInvert (invF, F + 1, dF, T, n);
      
      /* now invF[0..dF-1] = Quo(x^(2dF-1), F) */
      outputf (OUTPUT_VERBOSE, "Computing 1/F took %ldms\n",
	       elltime (st, cputime ()));
      
      /* ----------------------------------------------
         |   F    |  invF  |   G   |         T        |
         ----------------------------------------------
         |  F(x)  | 1/F(x) |  ???  |      ???         |
         ---------------------------------------------- */
    }

  if (stop_asap != NULL && (*stop_asap)())
    goto clear_invF;


  /* start computing G with roots at i0*d, (i0+1)*d, (i0+2)*d, ... 
     where i0*d <= B2min < (i0+1)*d */
  MEMORY_TAG;
  G = init_list2 (dF, mpz_sizeinbase (modulus->orig_modulus, 2) + 
                      3 * mp_bits_per_limb);
  MEMORY_UNTAG;
  if (G == NULL)
    {
      youpi = ECM_ERROR;
      goto clear_invF;
    }

  st = cputime ();
  if (method == ECM_PM1)
    rootsG_state = pm1_rootsG_init ((mpres_t *) X, root_params, modulus);
  else if (method == ECM_PP1)
    rootsG_state = pp1_rootsG_init ((mpres_t *) X, root_params, modulus);
  else /* ECM_ECM */
    rootsG_state = ecm_rootsG_init (f, (curve *) X, root_params, dF, k, 
                                    modulus);

  /* rootsG_state=NULL if an error occurred or (ecm only) a factor was found */
  if (rootsG_state == NULL)
    {
      /* ecm: f = -1 if an error occurred */
      youpi = (method == ECM_ECM && mpz_cmp_si (f, -1)) ? 
              ECM_FACTOR_FOUND_STEP2 : ECM_ERROR;
      goto clear_G;
    }

  if (method != ECM_ECM) /* ecm_rootsG_init prints itself */
    outputf (OUTPUT_VERBOSE, "Initializing table of differences for G "
             "took %ldms\n", elltime (st, cputime ()));

  if (stop_asap != NULL && (*stop_asap)())
    goto clear_fd;

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
          youpi = ECM_FACTOR_FOUND_STEP2;
          goto clear_fd;
        }

    if (stop_asap != NULL && (*stop_asap)())
      goto clear_fd;

  /* -----------------------------------------------
     |   F    |  invF  |   G    |         T        |
     -----------------------------------------------
     |  F(x)  | 1/F(x) | rootsG |      ???         |
     ----------------------------------------------- */

      st = cputime ();

      if (use_ntt)
        ntt_PolyFromRoots (G, G, dF, T + dF, mpzspm);
      else
        PolyFromRoots (G, G, dF, T + dF, n);

      /* needs 2*dF+list_mul_mem(dF/2) cells in T */
      outputf (OUTPUT_VERBOSE, "Building G from its roots took %ldms\n", 
               elltime (st, cputime ()));

    if (stop_asap != NULL && (*stop_asap)())
      goto clear_fd;

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
          if (use_ntt)
	    {
	      ntt_mul (T + dF, G, H, dF, T + 3 * dF, 0, mpzspm);
	      list_mod (H, T + dF, 2 * dF, n);
	    }
	  else
	    list_mulmod (H, T + dF, G, H, dF, T + 3 * dF, n);

          outputf (OUTPUT_VERBOSE, "Computing G * H took %ldms\n", 
                   elltime (st, cputime ()));

          if (stop_asap != NULL && (*stop_asap)())
            goto clear_fd;

          /* ------------------------------------------------
             |   F    |  invF  |    G    |         T        |
             ------------------------------------------------
             |  F(x)  | 1/F(x) |G(x)-F(x)| G * H  |         |
             ------------------------------------------------ */

	  st = cputime ();

          if (use_ntt)
	    {
	      ntt_PrerevertDivision (H, F, invF + 1, sp_F, sp_invF, dF,
		  T + 2 * dF, mpzspm);
	    }
	  else
	    {
	      if (PrerevertDivision (H, F, invF + 1, dF, T + 2 * dF, n))
	        {
	          youpi = ECM_ERROR;
	          goto clear_fd;
	        }
	    }
          
	  outputf (OUTPUT_VERBOSE, "Reducing  G * H mod F took %ldms\n", 
                   elltime (st, cputime ()));

          if (stop_asap != NULL && (*stop_asap)())
            goto clear_fd;
	}
    }
  
  clear_list (F, dF + 1);
  F = NULL;
  clear_list (G, dF);
  G = NULL;
  st = cputime ();
#ifdef POLYEVALTELLEGEN
  if (use_ntt)
    youpi = ntt_polyevalT (T, dF, Tree, T + dF + 1, sp_invF,
	mpzspm, TreeFilename);
  else
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
  treefiles_used = 0; /* Polyeval deletes treefiles by itself */

  outputf (OUTPUT_VERBOSE, "Computing polyeval(F,G) took %ldms\n", 
           elltime (st, cputime ()));

  st = cputime ();
  list_mulup (T, dF, n, T[dF]);
  outputf (OUTPUT_VERBOSE, "Computing product of all F(g_i) took %ldms\n", 
           elltime (st, cputime ()));

  mpz_gcd (f, T[dF - 1], n);
  if (mpz_cmp_ui (f, 1) > 0)
    {
      youpi = ECM_FACTOR_FOUND_STEP2;
      if (method == ECM_ECM && test_verbose (OUTPUT_RESVERBOSE))
        {
          /* Find out for which i*X, (i,d)==1, a factor was found */
          /* Note that the factor we found may be composite */
          /* TBD: use binary search */
          unsigned long j, k;
          mpz_set (T[dF], f);
          for (k = 0, j = 1; k < dF; j += 6)
            {
              if (gcd (j, root_params->d1) > 1)
                continue;
              mpz_gcd (T[dF + 1], T[k], T[dF]);
              if (mpz_cmp_ui (T[dF + 1], 1) > 0)
                {
                  long i; /* We need it as a signed type here */
                  /* Find i so that $f(i d1) X = f(j d2) X$ over GF(f) */
                  i = ecm_findmatch (j, root_params, (curve *)X, modulus, f);
                  
                  mpz_add_ui (T[dF + 2], root_params->i0, abs (i));
                  outputf (OUTPUT_RESVERBOSE, 
                           "Divisor %Zd first occurs in T[%lu] = "
                           "((f(%Zd*%lu)%cf(%lu*%lu))*X)_x\n", T[dF + 1], k, 
                           T[dF + 2], root_params->d1, i < 0 ? '+' : '-',
                           j, root_params->d2);
                  
                  mpz_mul_ui (T[dF + 2], T[dF + 2], root_params->d1);
                  if (i < 0)
                    mpz_add_ui (T[dF + 2], T[dF + 2], j * root_params->d2);
                  else
                    mpz_sub_ui (T[dF + 2], T[dF + 2], j * root_params->d2);
                  mpz_abs (T[dF + 2], T[dF + 2]);
                  
                  outputf (OUTPUT_RESVERBOSE, "Maybe largest group order "
                           "factor is or divides %Zd\n", T[dF + 2]);
                  
                  /* Don't report this divisor again */
                  mpz_divexact (T[dF], T[dF], T[dF + 1]);
                }
              k++;
            }
        }
    } else {
      /* Here, mpz_cmp_ui (f, 1) == 0, i.e. no factor was found */
      outputf (OUTPUT_RESVERBOSE, "Product of G(f_i) = %Zd\n", T[0]);
    }

clear_fd:
  if (method == ECM_PM1)
    pm1_rootsG_clear ((pm1_roots_state *) rootsG_state, modulus);
  else if (method == ECM_PP1)
    pp1_rootsG_clear ((pp1_roots_state *) rootsG_state, modulus);
  else /* ECM_ECM */
    ecm_rootsG_clear ((ecm_roots_state *) rootsG_state, modulus);

clear_G:
  clear_list (G, dF);
clear_invF:
  clear_list (invF, dF + 1);

  if (use_ntt)
    {
      mpzspv_clear (sp_F, mpzspm);
      mpzspv_clear (sp_invF, mpzspm);
    }

free_Tree_i:
  if (Tree != NULL)
    {
      for (i = 0; i < lgk; i++)
        clear_list (Tree[i], dF);
      free (Tree);
    }
  if (TreeFilename != NULL && treefiles_used > 0)
    {
      /* Unlink any treefiles still in use */
      char *fullname = (char *) malloc (strlen (TreeFilename) + 1 + 2 + 1);
      for (i = 0; i < treefiles_used; i++)
        {
          sprintf (fullname, "%s.%lu", TreeFilename, i);
          outputf (OUTPUT_DEVVERBOSE, "Unlinking %s\n", fullname);
          if (unlink (fullname) != 0)
            outputf (OUTPUT_ERROR, "Could not delete %s\n", fullname);
        }
      free (fullname);
    }
  mpz_clear (n);

clear_T:
  clear_list (T, sizeT);
clear_F:
  clear_list (F, dF + 1);

clear_i0:

  if (use_ntt)
    mpzspm_clear (mpzspm);
  
  if (Fermat)
    F_clear ();
  

  if (stop_asap == NULL || !(*stop_asap)())
    {
      st0 = elltime (st0, cputime ());
      outputf (OUTPUT_NORMAL, "Step 2 took %ldms\n", st0);
    }

  return youpi;
}
