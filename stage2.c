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

/* Init table to allow successive computation of x^((s + n*D)^E) mod N */
/* See Knuth, TAOCP vol.2, 4.6.4 and exercise 7 in 4.6.4 */
mpz_t *fin_diff_init(mpz_t x, unsigned int s, unsigned int D, unsigned int E,
                     mpz_t N) {
  mpz_t *fd;
  unsigned int i, k;
  
  fd = (mpz_t *) malloc((E+1) * sizeof(mpz_t));
  for (i=0; i <= E; i++)
    mpz_init(fd[i]);
  
  for (i=0; i <= E; i++)
    mpz_ui_pow_ui(fd[i], s+i*D, E);
  
  for (k=1; k<=E; k++)
    for (i=E; i>=k; i--)
      mpz_sub(fd[i], fd[i], fd[i-1]);
  
  for (i=0; i<=E; i++)
    mpz_powm(fd[i], x, fd[i], N);

  return fd;
}

/* Computes x^((s + (n+1)*D)^E and stores in fd[0] */
void 
fin_diff_next(mpz_t *fd, unsigned int E, mpz_t N) {
  unsigned int i;
  
  for (i=0; i<E; i++) {
    mpz_mul(fd[i], fd[i], fd[i+1]);
    mpz_mod(fd[i], fd[i], N);
  }
  return;
}

void 
fin_diff_clear(mpz_t *fd, unsigned int E) {
  unsigned int i;
  
  for (i=0; i<=E; i++)
    mpz_clear(fd[i]);
  free(fd);
  return;
}

/* puts in F the successive values of s^(j^S), for 0 < j < d, j=1 mod 7, 
   j and d coprime.
   Returns the degree of F, or 0 if a factor was found.
*/
int
rootsF (listz_t F, unsigned int d, mpz_t s, mpz_t invs, listz_t t, 
        unsigned int S, mpz_t n, int verbose)
  {
    unsigned int i, j;
    int st;
    mpz_t *fd;
    
    st = cputime ();

    mpz_set(F[0], s);
    i = 1;
    
    if (d > 7) 
      {
        fd = fin_diff_init(s, 7, 6, S, n);
        j = 7;
        while (j < d) 
          {
            if (gcd(j, d) == 1)
              mpz_set(F[i++], fd[0]);
            fin_diff_next(fd, S, n);
            j += 6;
          }
        fin_diff_clear(fd, S);
      }
    
    if (S > 3) /* Batch inversion is cheaper */
      {
        if (list_invert(t, F, i, t[i], n)) {
          mpz_set(s, t[i]);
          return 0;
        }

       for (j = 0; j < i; j++) {
         mpz_add(F[j], F[j], t[j]);
         mpz_mod(F[j], F[j], n);
       }

       mpz_set(invs, t[0]); /* Save s^(-1) in invs */
       
      } else { /* S <= 3: fin_diff code is cheaper */

        mpz_gcdext (*t, invs, NULL, s, n);

        if (mpz_cmp_ui (*t, 1) != 0)
          {
            mpz_set (s, *t);
            return 0;
          }

        mpz_add(F[0], F[0], invs);
        mpz_mod(F[0], F[0], n);

        i = 1;
        
        if (d > 7) 
          {
            fd = fin_diff_init(invs, 7, 6, S, n);
            j = 7;
            while (j < d) 
              {
                if (gcd(j, d) == 1)
                  {
                    mpz_add(F[i], F[i], fd[0]);
                    mpz_mod(F[i], F[i], n);
                    i++;
                  }
                fin_diff_next(fd, S, n);
                j += 6;
              }
            fin_diff_clear(fd, S);
          }
     }
    
    if (verbose >= 2)
      printf ("Computing roots of F took %dms\n", cputime () - st);

    return i;
}


/* puts in G the successive values of t0*s^j, for 1 <= j <= d
   returns in t the value of t0*s^d, in u the value of u0*invs^d
*/
void
rootsG (mpz_t *G, unsigned int d, mpz_t s, mpz_t invs, listz_t fd_x, 
        listz_t fd_invx, unsigned int S, mpz_t n, int verbose)
{
  unsigned int i;
  int st;

  st = cputime ();

  for (i = 0; i < d; i++)
    {
      mpz_add (G[i], fd_x[0], fd_invx[0]);
      mpz_mod (G[i], G[i], n);
      fin_diff_next(fd_x, S, n);
      fin_diff_next(fd_invx, S, n);
    }

  if (verbose >= 2)
    printf ("Computing roots of G took %dms\n", cputime () - st);
}

/* Input:  x is the value at end of stage 1
           n is the number to factor
           B2 is the stage 2 bound
           k is the number of blocks
   Output: x is the factor found
   Return value: non-zero iff a factor was found.
*/
int
stage2 (mpz_t x, mpz_t n, double B2, unsigned int k, unsigned int S, 
        int verbose)
{
  double b2;
  unsigned int i, d, dF, dG, sizeT;
  listz_t F, G, T, fd_x, fd_invx;
  polyz_t polyF, polyT;
  mpz_t invx;
  int youpi = 0, st, st0;

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
  dG = dF - 1;

  if (S == 0)
    S = 1;

  if (verbose >= 2)
    printf ("B2=%1.0f b2=%1.0f d=%u dF=%u dG=%u S=%u\n", B2, b2, d, dF, dG, S);

  F = init_list (dF + 1); 

  mpz_init (invx);
  sizeT = 3 * dF + list_mul_mem (dF);
  T = init_list (sizeT);

  if ((i = rootsF (F, d, x, invx, T, S, n, verbose)) == 0)
    {
      youpi = 2;
      goto clear_F;
    }

  assert (i == dF);

  buildG (F, dF, T, verbose, n, 'F'); /* needs dF+list_mul_mem(dF/2) cells in T */

  G = init_list (dG + 1);
  fd_x = fin_diff_init(x, 0, d, S, n);
  fd_invx = fin_diff_init(invx, 0, d, S, n);

  for (i=0; i<k; i++)
    {
      rootsG (G, dG, x, invx, fd_x, fd_invx, S, n, verbose);

      buildG (G, dG, T + dF, 1, n, 'G'); /* needs 2*dF+list_mul_mem(dF/2) cells in T */
      mpz_set_ui (G[dG], 1);

      if (i == 0)
	  list_set (T, G, dF);
      else
	{
	  st = cputime ();
	  /* previous G is in T, with degree < dF, i.e. dF coefficients
	     and dG = dF - 1, requires 3dF+list_mul_mem(dF) cells in T */
	  list_mulmod (T + dF, G, T, dF, T + 3 * dF, n);
          if (verbose >= 2)
            printf ("Computing G * H took %dms\n", cputime() - st);
	  st = cputime ();
          mpz_set_ui (T[3*dF-1], 0); /* since RecursiveDivision expects a
                                        dividend of 2*dF coefficients */
	  RecursiveDivision (T, T + dF, F, dF, T + 3 * dF, n);
          list_set (T, T + dF, dF);
          if (verbose >= 2)
            printf ("Reducing G * H mod F took %dms\n", cputime() - st);
	}
    }

  mpz_set_ui (F[dF], 1); /* the leading monic coefficient needs to be stored
                             explicitely for polygcd */
  st = cputime ();
  init_poly_list (polyF, dF, F);
  init_poly_list (polyT, dF - 1, T);
  if ((youpi = poly_gcd (x, polyF, polyT, n, T + dF)))
    NTL_get_factor (x);
  if (verbose >= 2)
    printf ("Computing gcd of F and G took %dms\n", cputime() - st);

  clear_list (G, dG + 1);
  fin_diff_clear(fd_invx, S);
  fin_diff_clear(fd_x, S);
  clear_list (T, sizeT);

 clear_F:
  mpz_clear (invx);
  clear_list (F, dF + 1);

  if (verbose >= 1)
    printf ("Stage 2 took %dms\n", cputime() - st0);

  return youpi;
}
