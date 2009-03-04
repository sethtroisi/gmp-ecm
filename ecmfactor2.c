/* ecmfactor2.c - example of use of libecm.a with HECM.

  Copyright 2008 Paul Zimmermann.

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
  Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
  02111-1307, USA.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h> /* GMP header file */
#include "ecm.h" /* ecm header file */

/* wrapper for GMP-ECM stage2 for a curve in Weierstrass form

   y^2 = x^2 + A * x + B

   where B is implicitly defined by y^2 - (x^2 + A * x) mod n.
*/
int
ecmfactor2 (mpz_t f, mpz_t n, mpz_t A, mpz_t x, mpz_t y, mpz_t B2)
{
  ecm_params q;
  int res;

  gmp_printf ("Performing one curve with B2=%Zd\n", B2);

  ecm_init (q);

  q->sigma_is_A = -1; /* indicates that we give a curve in Weierstrass form */
  mpz_set (q->sigma, A);
  mpz_set (q->x, x);
  mpz_set (q->go, y);
  mpz_set (q->B2, B2);

  res = ecm_factor (f, n, 0.0, q);

  ecm_clear (q);

  return res;
}

int
main (int argc, char *argv[])
{
  mpz_t n, f, A, x, y, B2;
  int res;

  if (argc != 6)
    {
      fprintf (stderr, "Usage: ecmfactor2 <number> <A> <x> <y> <B2>\n");
      /* example: n=100000000000000000039 A=97286259809299325273
         x=53919461457074783540 y=67664780906038509384 B2=10000000 */
      exit (1);
    }

  mpz_init (n);
  /* read number on command line */
  if (mpz_set_str (n, argv[1], 10))
    {
      fprintf (stderr, "Invalid number: %s\n", argv[1]);
      exit (1);
    }

  /* read A */
  mpz_init (A);
  if (mpz_set_str (A, argv[2], 10))
    {
      fprintf (stderr, "Invalid A: %s\n", argv[2]);
      exit (1);
    }

  /* read x */
  mpz_init (x);
  if (mpz_set_str (x, argv[3], 10))
    {
      fprintf (stderr, "Invalid x: %s\n", argv[3]);
      exit (1);
    }

  /* read y */
  mpz_init (y);
  if (mpz_set_str (y, argv[4], 10))
    {
      fprintf (stderr, "Invalid y: %s\n", argv[4]);
      exit (1);
    }

  /* read stage 2 bound B2 */
  mpz_init (B2);
  if (mpz_set_str (B2, argv[5], 10))
    {
      fprintf (stderr, "Invalid B2: %s\n", argv[5]);
      exit (1);
    }

  mpz_init (f); /* for potential factor */

  res = ecmfactor2 (f, n, A, x, y, B2);

  if (res > 0)
    {
      printf ("found factor in step %u: ", res);
      mpz_out_str (stdout, 10, f);
      printf ("\n");
    }
  else if (res == ECM_NO_FACTOR_FOUND)
    printf ("found no factor\n");
  else
    printf ("error\n");

  mpz_clear (f);
  mpz_clear (n);
  mpz_clear (A);
  mpz_clear (x);
  mpz_clear (y);
  mpz_clear (B2);

  return 0;
}
