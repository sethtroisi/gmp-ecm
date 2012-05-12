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

/* puts in F the successive values of s^j, for 0 < j < d, prime j
   Return the degree of F, or 0 if a factor was found.
   Output: s is s0^6, invs = s0^(-6).
*/
int
rootsF (mpz_t *F, unsigned int d, mpz_t s, mpz_t invs, mpz_t t, mpz_t u,
	mpz_t n)
{
  unsigned int i, j;
  int st;

  st = cputime ();

  mpz_gcdext (t, invs, NULL, s, n);
  
  if (mpz_cmp_ui (t, 1) != 0)
    {
      mpz_set (s, t);
      return 0;
    }
  
  mpz_set (t, s); /* t = s0 */
  mpz_set (u, invs); /* u = 1/s0 */

  mpz_add (F[0], t, u);
  mpz_mod (F[0], F[0], n);

  mpz_powm_ui (s, s, 6, n); /* s = s0^6 */
  mpz_powm_ui (invs, invs, 6, n); /* invs = s0^(-6) */

  for (i = 1, j = 7 ; j < d; j += 6)
    {
      mpz_mul (t, t, s);
      mpz_mod (t, t, n);
      mpz_mul (u, u, invs);
      mpz_mod (u, u, n);
      if (gcd (j, d) == 1)
	{
	  mpz_add (F[i], t, u);
	  mpz_mod (F[i], F[i], n);
	  i ++;
	}
    }

  fprintf (stderr, "Computing roots of F took %dms\n", cputime () - st);

  return i;
}

/* puts in G the successive values of s^j, for 1 <= j <= d
*/
void
rootsG (mpz_t *G, unsigned int d, mpz_t s, mpz_t invs, mpz_t t, mpz_t u,
	mpz_t n)
{
  unsigned int i;
  int st;

  st = cputime ();

  for (i = 0; i < d; i++)
    {
      mpz_mul (t, t, s);
      mpz_mod (t, t, n);
      mpz_mul (u, u, invs);
      mpz_mod (u, u, n);

      /* now t = s^i, u = s^(-i) */
      mpz_add (G[i], t, u);
      mpz_mod (G[i], G[i], n);
    }

  fprintf (stderr, "Computing roots of G took %dms\n", cputime () - st);
}

/* Input:  x is the value at end of stage 1
           n is the number to factor
           B2 is the stage 2 bound
           k is the number of blocks
   Output: x is the factor found
   Return value: non-zero iff a factor was found.
*/
int
stage2 (mpz_t x, mpz_t n, double B2, unsigned int k)
{
  double b2;
  unsigned int d, dF, dG;
  mpz_t *F, *G, *T, invx, t, u;
  int i, youpi = 0, st;

  printf ("starting stage 2 with x=");
  mpz_out_str (stdout, 10, x);
  putchar ('\n');

  b2 = ceil(B2 / k); /* b2 = ceil(B2/k): small block size */

  d = bestD (b2);

  b2 = block_size (d);

  B2 = (double) k * b2;

  dF = phi (d) / 2;
  dG = dF - 1;

  printf ("B2=%1.0f b2=%1.0f d=%u dF=%u dG=%u\n", B2, b2, d, dF, dG);

  F = init_poly (dF + 1); 

  mpz_init (invx);
  mpz_init (t);
  mpz_init (u);
  if ((i = rootsF (F, d, x, invx, t, u, n)) == 0)
    {
      youpi = 1;
      goto clear_F;
    }

  assert (i == dF);
  
  T = init_poly (5 * dF - 1);
  buildG (F, dF, T, 1, n, 'F');

  mpz_powm_ui (x, x, d / 6, n);
  mpz_powm_ui (invx, invx, d / 6, n);

  /* now x = x0^d, and invx = x0^(-d) */

  G = init_poly (dF);

  mpz_set_ui (t, 1);
  mpz_set_ui (u, 1);

  for (i=0; i<k; i++)
    {
      rootsG (G, dG, x, invx, t, u, n);

      buildG (G, dG, T + dF, 1, n, 'G');
      mpz_set_ui (G[dG], 1);

      if (i == 0)
	  polyset (T, G, dF);
      else
	{
	  st = cputime ();
	  /* previous G is in T, with degree < dF, i.e. dF coefficients
	     and dG = dF - 1 */
	  polymulmod (T + dF, G, T, dF, T + 3 * dF, n);
	  fprintf (stderr, "Computing G * H took %dms\n", cputime() - st);
	  st = cputime ();
          mpz_set_ui (T[3*dF-1], 0); /* since RecursiveDivision expects a
                                        dividend of 2*dF coefficients */
	  RecursiveDivision (T, T + dF, F, dF, T + 3 * dF, n);
          polyset (T, T + dF, dF);
	  fprintf (stderr, "Reducing G * H mod F took %dms\n", cputime() - st);
	}
    }

  mpz_set_ui (F[dF], 1); /* the leading monic coefficient needs to be stored
                             explicitely for polygcd */
  st = cputime ();
  youpi = polygcd (x, F, T, dF, n, T + dF);
  fprintf (stderr, "Computing gcd of F and G took %dms\n", cputime() - st);

  clear_poly (G, dG);
  clear_poly (T, 2 * dF);

 clear_F:
  mpz_clear (invx);
  mpz_clear (t);
  mpz_clear (u);
  clear_poly (F, dF);

  return youpi;
}