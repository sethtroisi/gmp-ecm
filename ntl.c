/* Interface with NTL.

  Copyright 2002, 2003 Paul Zimmermann and Alexander Kruppa.

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
#include "gmp.h"
#include "ecm.h"
#include "NTL/ZZ_pX.h"
#include "NTL/version.h"

#ifdef NTL_CLIENT
NTL_CLIENT
#endif

static int ntl_found = 0;
static mpz_t factor_found;

/* copies a NTL bigint into a GMP mpz */
void gmp_of_ntl (mpz_t a, ZZ& b)
{
  unsigned char *s;
  long n;

  n = NumBytes (b); /* number of 256-digits to represent abs(b) */
  s = (unsigned char *) malloc (n * sizeof (char));
  if (s == NULL)
    {
      fprintf (stderr, "Error: not enough memory\n");
      exit (EXIT_FAILURE);
    }
  BytesFromZZ (s, b, n);
  mpz_set_ui (a, 0);
  while (n)
    {
      mpz_mul_2exp (a, a, 8);
      mpz_add_ui (a, a, s[--n]);
    }
  free (s);
}

void
NTL_divide_handler (const ZZ_p& b)
{
  ZZ p;
  ntl_found = 2;
  p = GCD(rep(b), ZZ_p::modulus());
  gmp_of_ntl (factor_found, p);
}

void
NTL_init (void)
{
  ZZ_p::DivHandler = NTL_divide_handler;
  mpz_init (factor_found);
}

void
NTL_clear (void)
{
  mpz_clear (factor_found);
}

void
NTL_get_factor (mpz_t p)
{
  mpz_set (p, factor_found);
}

int
NTL_major_version (void)
{
  return NTL_MAJOR_VERSION;
}

int
NTL_minor_version (void)
{
  return NTL_MINOR_VERSION;
}

/* copies a GMP mpz into a NTL bigint */
void
ntl_of_gmp (ZZ& a, mpz_t b)
{
  /* simple string interface */
  char *s;

  s = mpz_get_str (NULL, 10, b);
  a = to_ZZ (s);
  free (s);
}

void
ntl_set_poly (ZZ_pX& aa, polyz_t a)
{
  unsigned int da = a -> degree;
  ZZ t;

  for (long int i = da; i >= 0; i--)
    {
      ntl_of_gmp (t, a->coeff[i]);
      SetCoeff (aa, i, to_ZZ_p(t));
    }
}

int
ntl_poly_gcd (mpz_t p, polyz_t a, polyz_t b, mpz_t n)
{
  ZZ nn;

  ntl_found = 0;
  
  ntl_of_gmp (nn, n);
  ZZ_p::init (nn);

  ZZ_pX aa, bb, gg;
  
  ntl_set_poly (aa, a);
  ntl_set_poly (bb, b);

  GCD (gg, aa, bb);

  if (deg(gg) != 0)
    { /* x^i+1/x^i = (x^d)^j+1/(x^d)^j (mod n) for some i, j */
      ntl_found = 2;
      mpz_set (factor_found, n);
    }

  return ntl_found;
}
