/* A Very Simple trial division module

  Copyright 2003, 2004, 2005 Jim Fougeron, Paul Zimmermann.

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
#include <math.h> /* log */
#include "ecm-ecm.h"

int
trial_factor (mpcandi_t *n, double maxfact, int deep)
{
  unsigned long factors = 0, exponent;
  double p;

  getprime (FREE_PRIME_TABLE);  /* free the prime tables, and reinitialize */
  /* brain dead trial factor'r but it works */
  for (p = 2.0; p <= maxfact; p = getprime (p))
    {
      for (exponent = 0; mpcandi_t_addfoundfactor_d (n, p); exponent++);
      
      if (exponent)
	{
	  printf ("********** Factor found trial div: %lu\n", (unsigned long) p);
	  printf ("Found proven prime factor of %2lu digits: %lu",
	    (unsigned long) (log (p) / log (10.0)) + 1, (unsigned long) p);
	  
	  if (exponent > 1)
	    printf ("^%lu", exponent);
	  printf ("\n");
	  
	  factors += exponent;
	  if (!deep)
	    /* We only want the first factor if not in "deep" mode */
	    break;
	}
    }
  getprime (FREE_PRIME_TABLE);  /* free the prime tables, and reinitialize */

  return factors;
}
