/* Tune program.

  Copyright 2003, 2005 Paul Zimmermann, Alexander Kruppa, Dave Newman.

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

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "ecm.h"
#include "ecm-gmp.h"
#include "ecm-impl.h"


/* we don't need any more precision */
#define GRANULARITY 1e2

/* Throughout, each function pointer points to a function
 * 
 *   double f0 (size_t limbs, unsigned int mintime);
 *
 * that runs for at least mintime ms and then returns the number of iterations
 * performed per ms. */

double
tune_mpres_mul (size_t limbs, unsigned int mintime, int repr)
{
  mpmod_t modulus;
  mpres_t x, y, z;
  mpz_t N, p, q;
  unsigned int st, k = 0;

  mpz_init (N);
  mpz_init (p);
  mpz_init (q);
  
  /* No need to generate a probable prime, just ensure N is not
     divisible by 2 or 3 */
  do
    {
      mpz_random (N, limbs);
      while (mpz_gcd_ui (NULL, N, 6) != 1)
        mpz_add_ui (N, N, 1);
    }
  while ((mp_size_t) mpz_size (N) != limbs);
  
  switch (repr)
  {
    case ECM_MOD_MPZ:
      mpmod_init_MPZ (modulus, N);
      break;
    case ECM_MOD_MODMULN:
      mpmod_init_MODMULN (modulus, N);
      break;
    case ECM_MOD_REDC:
      mpmod_init_REDC (modulus, N);
      break;
  }

  mpz_random (p, limbs);
  mpz_random (q, limbs);
  mpz_mod (p, p, N);
  mpz_mod (q, q, N);
  
  mpres_init (x, modulus);
  mpres_init (y, modulus);
  mpres_init (z, modulus);

  mpres_set_z (x, p, modulus);
  mpres_set_z (y, q, modulus);

  st = cputime ();

  do
    {
      mpres_mul (z, x, y, modulus);
      k++;
    }
  while (cputime () - st < mintime);

  st = cputime () - st;

  mpres_clear (x, modulus);
  mpres_clear (y, modulus);
  mpres_clear (z, modulus);
  mpmod_clear (modulus);
  mpz_clear (N);
  mpz_clear (p);
  mpz_clear (q);

  return (double) k / (double) st;
}

double
tune_mpres_mul_mpz (size_t n, unsigned int mintime)
{
  return tune_mpres_mul (n, mintime, ECM_MOD_MPZ);
}

double
tune_mpres_mul_modmuln (size_t n, unsigned int mintime)
{
  return tune_mpres_mul (n, mintime, ECM_MOD_MODMULN);
}

double
tune_mpres_mul_redc (size_t n, unsigned int mintime)
{
  return tune_mpres_mul (n, mintime, ECM_MOD_REDC);
}


/* Assume f0 and f1 are monotone decreasing. Return the first n in the range
 * [min_n, max_n) for which f1(n) > f0(n), or return max_n if no such n
 * exists. */
size_t
crossover (double (*f0)(size_t, unsigned int),
    double (*f1)(size_t, unsigned int), size_t min_n, size_t max_n)
{
  size_t mid_n;
  
  if (min_n == max_n)
    return min_n;

  mid_n = (max_n + min_n) / 2;
  return ((f0)(mid_n, GRANULARITY) >= (f1)(mid_n, GRANULARITY))
    ? crossover (f0, f1, mid_n + 1, max_n)
    : crossover (f0, f1, min_n, mid_n);
}

/* Return the lowest n with min_n <= n < max_n such that
 * f1(t) > f0(t) for all t in [n, n + k)
 *
 * Return max_n if no such n exists. */
size_t
crossover2 (double (*f0)(size_t, unsigned int),
    double (*f1)(size_t, unsigned int), size_t min_n, size_t max_n, size_t k)
{
  size_t n = min_n;
  size_t t;
  
  while (n < max_n)
    {
      for (t = n + k - 1; t > n; t--)
	if ((f0)(t, GRANULARITY) >= (f1)(t, GRANULARITY))
	  break;

      if (t == n)
	return n;

      n = t + 1;
    };

  return max_n;
}


int main ()
{
  size_t mpzmod_threshold, redc_threshold;
  
  mpzmod_threshold = crossover2 (tune_mpres_mul_modmuln, tune_mpres_mul_mpz,
      1, 512, 10);
  redc_threshold = crossover2 (tune_mpres_mul_mpz, tune_mpres_mul_redc,
      mpzmod_threshold, 512, 10);
  
  printf ("#define MPZMOD_THRESHOLD %u\n", mpzmod_threshold);
  printf ("#define REDC_THRESHOLD %u\n", redc_threshold);

  return 0;
}

  
