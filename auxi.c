/* Auxiliary functions for GMP-ECM.

  Copyright 2002, 2003, 2004, 2005 Paul Zimmermann, Alexander Kruppa, Laurent Fousse, Jim Fougeron.

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

#include <gmp.h>
#include "ecm-ecm.h"

/******************************************************************************
*                                                                             *
*                            Auxiliary functions                              *
*                                                                             *
******************************************************************************/

/* returns the number of decimal digits of n */
unsigned int
nb_digits (const mpz_t n)
{
   mpz_t x;
   unsigned int size;
   
   mpz_init_set (x, n);
   
   for (size = 0; mpz_size (x); size++)
     mpz_tdiv_q_ui (x, x, 10);

   mpz_clear (x);

   return size;
}

