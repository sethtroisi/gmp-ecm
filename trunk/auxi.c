/* Auxiliary functions for GMP-ECM.

  Copyright 2001, 2002, 2003, 2004, 2005 Paul Zimmermann and Alexander Kruppa.

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
#include <string.h>
#include <gmp.h>
#include "ecm.h"
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
   unsigned int size;
   char *str;
 
   str = mpz_get_str (NULL, 10, n);
   size = strlen (str);
   FREE (str, size + 1); /* FIXME: use GMP's own free function
			    if/when available, alloc size is size+1 */
   return size;
}

