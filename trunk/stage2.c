/* Common stage 2 for ECM and 'P-1' (improved standard continuation).

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

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "gmp.h"
#include "ecm.h"
#include "cputime.h"

#define INVF /* precompute 1/F for divisions by F */

/* Init table to allow successive computation of x^((s + n*D)^E) mod N */
/* See Knuth, TAOCP vol.2, 4.6.4 and exercise 7 in 4.6.4 */
mpz_t *
fin_diff_init (mpz_t x, unsigned int s, unsigned int D, unsigned int E,
               mpz_t N) 
{
  mpz_t *fd;
  unsigned int i, k;
  
  fd = (mpz_t *) malloc ((E+1) * sizeof(mpz_t));
  for (i = 0; i <= E; i++)
    mpz_init (fd[i]);
  
  for (i = 0; i <= E; i++)
    mpz_ui_pow_ui (fd[i], s + i * D, E); /* fd[i] = (s+i*D)^E */
  
  for (k = 1; k <= E; k++)
    for (i = E; i >= k; i--)
      mpz_sub (fd[i], fd[i], fd[i-1]);
  
  for (i = 0; i <= E; i++)
    mpz_powm (fd[i], x, fd[i], N);

  return fd;
}

/* Computes x^((s + (n+1)*D)^E and stores in fd[0] */
void 
fin_diff_next (mpz_t *fd, unsigned int E, mpz_t N) 
{
  unsigned int i;
  
  for (i = 0; i < E; i++)
    {
      mpz_mul (fd[i], fd[i], fd[i+1]);
      mpz_mod (fd[i], fd[i], N);
    }
  return;
}

void 
fin_diff_clear (mpz_t *fd, unsigned int E)
{
  unsigned int i;
  
  for (i = 0; i <= E; i++)
    mpz_clear (fd[i]);
  free (fd);
  return;
}

/* puts in F[0..dF-1] the successive values of s^(j^S), for 0 < j = 1 mod 7 < d
   j and d coprime.
   Returns df = degree of F, or 0 if a factor was found.
   If invs=0, don't use the x+1/x trick.
   Requires (dF+1) cells in t.
*/
int
rootsF (listz_t F, unsigned int d, mpz_t s, mpz_t invs, listz_t t,
        unsigned int S, mpz_t n, int verbose)
{
  unsigned int i, j;
  int st, st2;
  mpz_t *fd;
  
  st = cputime ();

  mpz_set (F[0], s);
  i = 1;
  
  if (d > 7)
    {
      st2 = cputime ();
      fd = fin_diff_init (s, 7, 6, S, n);
      if (verbose >= 2)
        printf ("Initializing table of differences for F took %dms\n", cputime () - st2);
      j = 7;
      while (j < d) 
        {
          if (gcd (j, d) == 1)
            mpz_set (F[i++], fd[0]);
          fin_diff_next (fd, S, n);
          j += 6;
        }
      fin_diff_clear (fd, S);
    }

  if (mpz_cmp_ui (invs, 0) && (d/6)*S > 3*(i-1)) /* Batch inversion is cheaper */
    {
      if (list_invert (t, F, i, t[i], n)) 
        {
          mpz_set (s, t[i]);
          return 0;
        }
     
      for (j = 0; j < i; j++) 
        {
          mpz_add (F[j], F[j], t[j]);
          mpz_mod (F[j], F[j], n);
        }
     
      mpz_set (invs, t[0]); /* Save s^(-1) in invs */
    
    }
  else
    { /* fin_diff code is cheaper */

      mpz_gcdext (*t, invs, NULL, s, n);

      if (mpz_cmp_ui (*t, 1) != 0)
        {
          mpz_set (s, *t);
          return 0;
        }

      mpz_add (F[0], F[0], invs);
      mpz_mod (F[0], F[0], n);

      i = 1;
       
      if (d > 7) 
        {
          fd = fin_diff_init (invs, 7, 6, S, n);
          j = 7;
          while (j < d) 
            {
              if (gcd (j, d) == 1)
                {
                  mpz_add (F[i], F[i], fd[0]);
                  mpz_mod (F[i], F[i], n);
                  i++;
                }
              fin_diff_next (fd, S, n);
              j += 6;
            }
          fin_diff_clear (fd, S);
        }
    }
  
  if (verbose >= 2)
    printf ("Computing roots of F took %dms\n", cputime () - st);

  return i;
}


/* puts in G the successive values of t0*s^j, for 1 <= j <= d
   returns in t the value of t0*s^d, in u the value of u0*invs^d.
   If listz_t=NULL, don't use the x+1/x trick.
   Needs d+1 cells in t.
*/
void
rootsG (listz_t G, unsigned int d, listz_t fd_x, listz_t fd_invx, 
        listz_t t, unsigned int S, mpz_t n, int verbose)
{
  unsigned int i;
  int st;

  st = cputime ();

  if (fd_invx != NULL && S > 3)
    {
      mpz_set (G[0], fd_x[0]);
      mpz_set (t[0], fd_x[0]);

      for (i = 1; i < d; i++)
        {
          fin_diff_next (fd_x, S, n);
          mpz_set (G[i], fd_x[0]);
          mpz_mul (t[i], t[i-1], fd_x[0]);
          mpz_mod (t[i], t[i], n);
        }

      mpz_gcdext (t[d], t[d-1], NULL, t[d-1], n);

      for (i = d-1; i > 0; i--) 
        {
          mpz_mul (t[d], t[i], t[i-1]);
          mpz_mul (t[i-1], t[i], G[i]);
          mpz_mod (t[i-1], t[i-1], n);
          mpz_add (t[d], t[d], G[i]);
          mpz_mod (G[i], t[d], n);
        }

    }
  else /* fd_invx=NULL or S <= 2 */
    {
      for (i = 0; i < d; i++)
        {
          if (fd_invx != NULL)
            {
              mpz_add (G[i], fd_x[0], fd_invx[0]);
              mpz_mod (G[i], G[i], n);
              fin_diff_next (fd_invx, S, n);
            }
          else
            mpz_set (G[i], fd_x[0]);
          fin_diff_next (fd_x, S, n);
        }
    }

  if (verbose >= 2)
    printf ("Computing roots of G took %dms\n", cputime () - st);
}

/* Input:  x is the value at end of stage 1
           n is the number to factor
           B2 is the stage 2 bound
           k is the number of blocks
           S is the exponent for Brent-Suyama's extension
           verbose is the verbose level
           invtrick is non-zero iff one uses x+1/x instead of x.
           Cf "Speeding the Pollard and Elliptic Curve Methods
               of Factorization", Peter Montgomery, Math. of Comp., 1987,
               page 257: using x^(i^e)+1/x^(i^e) instead of x^(i^(2e))
               reduces the cost of Brent-Suyama's extension from 2*e
               to e+3 multiplications per value of i.
   Output: x is the factor found
   Return value: non-zero iff a factor was found.
*/
int
stage2 (mpz_t x, mpz_t n, double B2, unsigned int k, unsigned int S, 
        int verbose, int invtrick)
{
  double b2;
  unsigned int i, d, dF, sizeT;
  unsigned long muls;
  listz_t F, G, H, T, fd_x, fd_invx = NULL;
  polyz_t polyF, polyT;
  mpz_t invx;
  int youpi = 0, st, st0;
#ifdef INVF
  listz_t invF;
#endif

  st0 = cputime ();

  if (verbose >= 2)
    {
      printf ("starting stage 2 with x=");
      mpz_out_str (stdout, 10, x);
      putchar ('\n');
    }

  b2 = ceil(B2 / k); /* b2 = ceil(B2/k): small block size */

  d = bestD (b2);

  b2 = block_size (d);

  B2 = (double) k * b2;

  dF = phi (d) / 2;

  if (verbose >= 2)
    printf ("B2=%1.0f k=%u b2=%1.0f d=%u dF=%u\n", B2, k, b2, d, dF);

  F = init_list (dF + 1);

  mpz_init_set_ui (invx, invtrick);
  sizeT = 3 * dF - 1 + list_mul_mem (dF);
#ifdef INVF
  if (dF > 3)
    sizeT += dF - 3;
#endif
  T = init_list (sizeT);
  H = T;

  /* needs dF+1 cells in T */
  if ((i = rootsF (F, d, x, invx, T, S, n, verbose)) == 0)
    {
      youpi = 2;
      goto clear_F;
    }

  /* ----------------------------------------------
     |   F    |  invF  |   G   |         T        |
     ----------------------------------------------
     | rootsF |  ???   |  ???  |      ???         |
     ---------------------------------------------- */

  assert (i == dF);

  PolyFromRoots (F, dF, T, verbose, n, 'F'); /* needs dF+list_mul_mem(dF/2) cells in T */
  mpz_set_ui (F[dF], 1); /* the leading monic coefficient needs to be stored
                             explicitely for PrerevertDivision and polygcd */

  /* ----------------------------------------------
     |   F    |  invF  |   G   |         T        |
     ----------------------------------------------
     |  F(x)  |  ???   |  ???  |      ???         |
     ---------------------------------------------- */

#ifdef INVF
  /* G*H has degree 2*dF-2, hence we must cancel dF-1 coefficients
     to get degree dF-1 */
  if (dF > 1)
    {
      invF = init_list (dF - 1);
      st = cputime ();
#if 0
      list_zero (T, 2 * dF - 3);
      mpz_set_ui (T[2 * dF - 3], 1); /* T[0..2dF-3] = x^(2dF-3) */
      muls = RecursiveDivision (invF, T, F + 1, dF - 1, T + 2 * dF - 2, n);
#else
      muls = PolyInvert (invF, F + 2, dF - 1, T, n);
#endif
      /* now invF[0..K-2] = Quo(x^(2dF-3), F) */
      if (verbose >= 2)
        printf ("Computing 1/F took %ums and %lumuls\n", cputime() - st, muls);
      
      /* ----------------------------------------------
         |   F    |  invF  |   G   |         T        |
         ----------------------------------------------
         |  F(x)  | 1/F(x) |  ???  |      ???         |
         ---------------------------------------------- */
    }
#endif

  G = init_list (dF);
  st = cputime ();
  fd_x = fin_diff_init (x, 0, d, S, n);
  if (verbose >= 2)
    printf ("Initializing table of differences for G took %dms\n", cputime () - st);
    
  if (invtrick)
    fd_invx = fin_diff_init (invx, 0, d, S, n);

  for (i=0; i<k; i++)
    {
      /* needs dF+1 cells in T+dF */
      rootsG (G, dF, fd_x, fd_invx, T + dF, S, n, verbose);

  /* -----------------------------------------------
     |   F    |  invF  |   G    |         T        |
     -----------------------------------------------
     |  F(x)  | 1/F(x) | rootsG |      ???         |
     ----------------------------------------------- */

      PolyFromRoots (G, dF, T + dF, verbose, n, 'G'); /* needs 2*dF+list_mul_mem(dF/2) cells in T */

  /* -----------------------------------------------
     |   F    |  invF  |   G    |         T        |
     -----------------------------------------------
     |  F(x)  | 1/F(x) |  G(x)  |      ???         |
     ----------------------------------------------- */

      if (i == 0)
        {
          list_sub (H, G, F, dF); /* coefficients 1 of degree cancel,
                                     thus T is of degree < dF */
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

          /* ------------------------------------------------
             |   F    |  invF  |    G    |         T        |
             ------------------------------------------------
             |  F(x)  | 1/F(x) |G(x)-F(x)|  H(x)  |         |
             ------------------------------------------------ */

	  st = cputime ();
	  /* previous G mod F is in H, with degree < dF, i.e. dF coefficients:
	     requires 3dF-1+list_mul_mem(dF) cells in T */
	  muls = list_mulmod2 (H, T + dF, G, H, dF, T + 3 * dF - 1, n);
          if (verbose >= 2)
            printf ("Computing G * H took %ums and %lumuls\n", cputime() - st,
                    muls);

          /* ------------------------------------------------
             |   F    |  invF  |    G    |         T        |
             ------------------------------------------------
             |  F(x)  | 1/F(x) |G(x)-F(x)| G * H  |         |
             ------------------------------------------------ */

	  st = cputime ();
#ifdef INVF
          muls = PrerevertDivision (H, F, invF, dF, T + 2 * dF - 1, n);
#else
          mpz_set_ui (T[2*dF-1], 0); /* since RecursiveDivision expects a
                                        dividend of 2*dF coefficients */
	  muls = RecursiveDivision (G, H, F, dF, T + 2 * dF, n);
#endif
          if (verbose >= 2)
            printf ("Reducing G * H mod F took %ums and %lumuls\n",
                    cputime() - st, muls);
	}
    }

  st = cputime ();
  init_poly_list (polyF, dF, F);
  init_poly_list (polyT, dF - 1, T);
  if ((youpi = poly_gcd (x, polyF, polyT, n, T + dF)))
    NTL_get_factor (x);
  if (verbose >= 2)
    printf ("Computing gcd of F and G took %dms\n", cputime() - st);

  clear_list (G, dF);
  if (invtrick)
    fin_diff_clear (fd_invx, S);
  fin_diff_clear (fd_x, S);
  clear_list (T, sizeT);

#ifdef INVF
  clear_list (invF, dF - 1);
#endif

 clear_F:
  mpz_clear (invx);
  clear_list (F, dF + 1);

  if (verbose >= 1)
    printf ("Stage 2 took %dms\n", cputime() - st0);

  return youpi;
}
