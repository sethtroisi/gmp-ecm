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
#include "gmp.h"
#include "ecm.h"

/******************************************************************************
*                                                                             *
*                                  Stage 1                                    *
*                                                                             *
******************************************************************************/

/* a <- b * c where a and b are mpz, c is a double, and t an auxiliary mpz */
#if (BITS_PER_MP_LIMB >= 53)
#define mpz_mul_d(a, b, c, t) \
   mpz_mul_ui (a, b, (unsigned long int) c);
#else
#if (BITS_PER_MP_LIMB >= 32)
#define mpz_mul_d(a, b, c, t) \
   if (c < 4294967296.0) \
      mpz_mul_ui (a, b, (unsigned long int) c); \
   else { \
   mpz_set_d (t, c); \
   mpz_mul (a, b, t); }
#else
#define mpz_mul_d(a, b, c, t) \
   mpz_set_d (t, c); \
   mpz_mul (a, b, t);
#endif
#endif

/* prime powers are accumulated up to about n^L1 */
#define L1 16

/* Input:  a is the generator (sigma)
           n is the number to factor
           B1 is the stage 1 bound
   Output: a is the factor found, or the value at end of stage 1
   Return value: non-zero iff a factor was found.
*/
int
stage1 (mpz_t a, mpz_t n, double B1)
{
  double B0, p, q, r;
  mpz_t g, d;
  int youpi;
  unsigned int max_size;

  mpz_init_set_ui (g, 1);
  mpz_init (d);

  B0 = sqrt (B1);

  max_size = L1 * mpz_sizeinbase (n, 2);

  /* first loop through small primes <= sqrt(B1) */
  for (p = 2.0; p <= B0; p = getprime(p))
    {
      for (q = 1, r = p; r <= B1; q = r, r *= p);
      mpz_mul_d (g, g, q, d);
      if (mpz_sizeinbase (g, 2) >= max_size)
	{
	  mpz_powm (a, a, g, n);
	  mpz_set_ui (g, 1);
	}
    }

  /* then all primes > sqrt(B1) and taken with exponent 1 */
  for (; p <= B1; p = getprime(p))
    {
      mpz_mul_d (g, g, p, d);
      if (mpz_sizeinbase (g, 2) >= max_size)
	{
	  mpz_powm (a, a, g, n);
	  mpz_set_ui (g, 1);
	}
    }

  getprime (0.0); /* free the prime tables, and reinitialize */

  mpz_clear (d);

  mpz_powm (a, a, g, n);

  mpz_sub_ui (g, a, 1);
  mpz_gcd (g, g, n);
  if ((youpi = mpz_cmp_ui (g, 1)))
    mpz_set (a, g);

  mpz_clear (g);
  
  return youpi;
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
          k is the number of blocks for stage 2
          verbose is the verbose level: 0=quiet, 1=normal, 2=verbose
   Output: p is the factor found
   Return value: non-zero iff a factor is found (1 for stage 1, 2 for stage 2)
*/
int
pm1 (mpz_t p, mpz_t n, double B1, double B2, unsigned int k, unsigned int S, 
     int verbose)
{
  int youpi, st;

  st = cputime ();
  youpi = stage1 (p, n, B1);
  if (verbose >= 1)
    {
      printf ("Stage 1 took %dms\n", cputime() - st);
      fflush (stdout);
    }

  if (youpi != 0) /* a factor was found */
    return 1;

  return (B2 > B1) ? stage2 (p, n, B2, k, S, verbose) : 0;
}
