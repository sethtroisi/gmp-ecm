/* Polynomial arithmetic.

  Copyright (C) 2002 Paul Zimmermann,
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
#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "ecm.h"

/* initializes a polynomial of degree n */
void
init_poly (polyz_t p, int n)
{
  p->coeff = init_list (n + 1);
  p->alloc = n + 1;
}

/* initializes a polynomial of degree n from a given list */
void
init_poly_list (polyz_t p, int n, listz_t l)
{
  p->coeff = l;
  p->degree = n;
  p->alloc = n + 1;
}

/* clears a polynomial */
void
clear_poly (polyz_t p)
{
  clear_list (p->coeff, p->alloc);
}

/* returns non-zero iff p is 0 */
int
poly_zerop (polyz_t p)
{
  return list_zerop (p->coeff, p->degree + 1);
}

void
poly_set_ui (polyz_t p, unsigned long int i)
{
  assert(p->alloc >= 1);
  mpz_set_ui (p->coeff[0], i);
}

/* Input: a is a polynomial.
          b is a polynomial with deg(b) < deg(a).
   Output: g = gcd(a, b) is in a.
   Return value: deg(g), or -1 if a factor of n was found, in which case
                 this factor is stored in p.
*/
int
poly_gcd (mpz_t p, polyz_t a, polyz_t b, mpz_t n, listz_t t)
{
  int result;

#ifdef MEMORY_DEBUG
  tests_memory_reset ();
#endif
  result = ntl_poly_gcd (p, a, b, n);
#ifdef MEMORY_DEBUG
  tests_memory_start ();
#endif
  return result;
}

