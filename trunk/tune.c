/* Tune program.

  Copyright 2003, 2005 Paul Zimmermann and Alexander Kruppa.

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
#include "ecm-gmp.h"
#include "ecm-impl.h"

#define MINTIME 1000 /* one second */

/* performs k computations of p*q mod N using representation 'repr'
   and return the total time.
*/
static unsigned int
test (mpz_t N, mpz_t p, mpz_t q, int repr, int k)
{
  mpmod_t modulus;
  mpres_t x, y, z;
  unsigned int st;

  if (repr == 1)
    mpmod_init_MPZ (modulus, N);
  else if (repr == 2)
    mpmod_init_MODMULN (modulus, N);
  else if (repr == 3)
    mpmod_init_REDC (modulus, N);

  mpres_init (x, modulus);
  mpres_init (y, modulus);
  mpres_init (z, modulus);

  mpres_set_z (x, p, modulus);
  mpres_set_z (y, q, modulus);

  st = cputime ();

  while (k--)
    mpres_mul (z, x, y, modulus);

  st = elltime (st, cputime ());

  mpres_clear (x, modulus);
  mpres_clear (y, modulus);
  mpres_clear (z, modulus);
  mpmod_clear (modulus);

  return st;
}

int
main (int argc, char *argv[])
{
  mp_size_t n, n0;
  mpz_t N, p, q;
  int k;
  unsigned int st[3];
  int mpzmod_threshold = 0;
  int redc_threshold = 0;

  printf ("MUL_KARATSUBA_THRESHOLD=%u\n", MUL_KARATSUBA_THRESHOLD);
  printf ("DIV_DC_THRESHOLD=%u\n", DIV_DC_THRESHOLD);
  printf ("MPZMOD_THRESHOLD_DEFAULT=%u\n", MPZMOD_THRESHOLD_DEFAULT);
  printf ("REDC_THRESHOLD_DEFAULT=%u\n", REDC_THRESHOLD_DEFAULT);

  n0 = (argc > 1) ? atoi (argv[1]) : 1;

  mpz_init (N);
  mpz_init (p);
  mpz_init (q);

  printf ("n\tmpzmod\tmodmuln\tredc\n");

  for (n = n0; ; n++)
    {

      printf ("%lu\t", n);

      /* no need to generate a probable prime, just ensure N is not
         divisible by 2 or 3 */
      do
        {
          mpz_random (N, n);
          while (mpz_gcd_ui (NULL, N, 6) != 1)
            mpz_add_ui (N, N, 1);
        }
      while ((mp_size_t) mpz_size (N) != n);

      mpz_random (p, n);
      mpz_mod (p, p, N);

      mpz_random (q, n);
      mpz_mod (q, q, N);

      /* first calibrate */
      for (k = 1; (st[0] = test (N, p, q, 1, k)) < MINTIME; k *= 2);

      k = (int) (((double) k * (double) MINTIME) / (double) st[0]);

      printf ("%u\t", st[0] = test (N, p, q, 1, k)); /* mpzmod */
      printf ("%u\t", st[1] = test (N, p, q, 2, k)); /* modmuln */
      printf ("%u\t", st[2] = test (N, p, q, 3, k)); /* redc */

      /* since modmuln is O(n^2), we should have asymptotically
         mpzmod faster than modmuln.
         Also, mpzmod uses plain division which is asymptotically
         about k multiplications with k > 2 (k=2 in the Karatsuba range),
         whereas redc is equivalent to two multiplications, thus
         asymptotically redc is faster than mpzmod.
      */

      if (st[0] < st[1] && st[0] <= st[2])
        printf ("\tmpzmod\n");
      else if (st[1] <= st[0] && st[1] <= st[2])
        printf ("\tmodmuln\n");
      else
        printf ("\tredc\n");

      if (st[0] >= st[1])
        mpzmod_threshold = n + 1;

      if (st[2] >= st[0])
        redc_threshold = n + 1;

      /* stop when order did not change for past 10 sizes */
      if (n > n0 + 10 && mpzmod_threshold + 10 < n && redc_threshold + 10 < n)
        break;
    }

  printf ("#define MPZMOD_THRESHOLD %u\n", mpzmod_threshold);
  printf ("#define REDC_THRESHOLD %u\n", redc_threshold);

  mpz_clear (p);
  mpz_clear (q);
  mpz_clear (N);
  
  return 0;
}
