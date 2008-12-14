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

int
main (int argc, char *argv[])
{
  mpz_t n, f;
  int res;
  ecm_params q;

  if (argc != 6)
    {
      fprintf (stderr, "Usage: ecmfactor2 <number> <A> <x> <y> <B2>\n");
      /* example: n=100000000000000000039 A=2098635068150078946
         x=42553831306919029237 y=70495752344261309320 B2=12408991157
	 (with B=16991832675985765075) */
      exit (1);
    }

  ecm_init (q);
  q->sigma_is_A = -1; /* indicates that we give a curve in Weierstrass form */

  mpz_init (n);
  /* read number on command line */
  if (mpz_set_str (n, argv[1], 10))
    {
      fprintf (stderr, "Invalid number: %s\n", argv[1]);
      exit (1);
    }

  /* read A */
  if (mpz_set_str (q->sigma, argv[2], 10))
    {
      fprintf (stderr, "Invalid A: %s\n", argv[2]);
      exit (1);
    }

  /* read x */
  if (mpz_set_str (q->x, argv[3], 10))
    {
      fprintf (stderr, "Invalid x: %s\n", argv[3]);
      exit (1);
    }

  /* read y */
  if (mpz_set_str (q->go, argv[4], 10))
    {
      fprintf (stderr, "Invalid y: %s\n", argv[4]);
      exit (1);
    }

  if (mpz_set_str (q->B2, argv[5], 10))
    {
      fprintf (stderr, "Invalid B2: %s\n", argv[5]);
      exit (1);
    }

  mpz_init (f); /* for potential factor */

  gmp_printf ("Performing one curve with B2=%Zd\n", q->B2);

  res = ecm_factor (f, n, 0.0, q);

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

  ecm_clear (q);

  return 0;
}
