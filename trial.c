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
#include <string.h>
#if !defined (_MSC_VER)
#include <unistd.h>
#endif
#include <gmp.h>
#include "ecm.h"
#include "ecm-ecm.h"

int
trial_factor (mpcandi_t *n, double maxfact, int deep)
{
  int factors = 0, cnt_this_fact;
  mpz_t t;
  unsigned long Remainder;
  double p;
  char numbuf[40];

  mpz_init (t);

  getprime (FREE_PRIME_TABLE);  /* free the prime tables, and reinitialize */
  /* brain dead trial factor'r but it works */
  for (p = 2.0; p <= maxfact; p = getprime (p))
    {
      Remainder = mpz_mod_ui (t, n->n, (unsigned long) p);
      if (!Remainder)
	{
	  cnt_this_fact = 0;
	  sprintf (numbuf, "%.0f", p);
	  while (mpcandi_t_addfoundfactor_d (n, p))
	    ++cnt_this_fact;
	  printf ("********** Factor found trial div: %s\n", numbuf);
	  if (cnt_this_fact > 1)
	    printf ("Found Proven Prime   factor of %2u digits: %s^%d\n",
		    (unsigned int) strlen (numbuf), numbuf, cnt_this_fact);
	  else
	    printf ("Found Proven Prime   factor of %2u digits: %s\n",
		    (unsigned int) strlen (numbuf), numbuf); 
	  factors += cnt_this_fact;
	  if (!deep)
	    /* We only want the first factor if not in "deep" mode */
	    break;
	}
    }
  mpz_clear (t);
  getprime (FREE_PRIME_TABLE);  /* free the prime tables, and reinitialize */

  return factors;
}
