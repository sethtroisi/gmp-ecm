/* Dynamic Eratosthenes sieve.

  Copyright 2001, 2002, 2003 Paul Zimmermann and Alexander Kruppa.

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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "ecm.h"

/*
  B1          10^5    10^6     10^7     10^8        1e9        1e10
  Pi(B1)      9592    78498    664579   5761455   50847534    455052511
  prev(B1)    99991   999983   9999991  99999989  999999937   9999999967
  ecm/step1   0.1     2.02     42.7     
  getprime    0.0     0.04     0.38     3.87        43.64      536.8
*/

/* This function returns successive odd primes, starting with 3.
   To perform a loop over all primes <= B1, do the following:

      for (p = 2.0; p <= B1; p = getprime (p))
         {
            ...
         }
*/
double
getprime (double pp)
{
  static double offset = 0.0; /* offset for current primes */
  static long int current = -1; /* index of previous prime */
  static unsigned *primes = NULL; /* table of small primes up to sqrt(p) */
  static unsigned long int nprimes = 0; /* length of primes[] */
  static unsigned char *sieve = NULL; /* sieving table */
  static long int len = 0; /* length of sieving table */
  static unsigned long int *moduli = NULL; /* offset for small primes */

  if (WANT_FREE_PRIME_TABLE(pp)) /* free the tables, and reinitialize */
    {
      offset = 0.0;
      current = -1;
      free (primes);
      primes = NULL;
      nprimes = 0;
      free (sieve);
      sieve = NULL;
      len = 0;
      free (moduli);
      moduli = NULL;
      return pp;
    }

  while ((++current < len) && (sieve[current] == 0));

  if (current < len) /* most calls will end here */
    return offset + 2.0 * (double) current;

  /* otherwise we have to sieve */
  offset += 2.0 * (double) len;

  /* first enlarge sieving table if too small */
  if ((double) len * (double) len < offset)
    {
      free (sieve);
      len *= 2;
      sieve = (unsigned char *) malloc (len * sizeof(unsigned char));
      if (sieve == NULL)
        {
          fprintf (stderr, "Error: not enough memory\n");
          exit (EXIT_FAILURE);
        }
    }

  /* now enlarge small prime table if too small */
  if ((nprimes == 0) || (primes[nprimes-1] < sqrt(offset + len)))
      {
	if (nprimes == 0) /* initialization */
	  {
	    nprimes = 1;
	    primes = (unsigned *) malloc (nprimes * sizeof(unsigned long int));
            if (primes == NULL)
              {
                fprintf (stderr, "Error: not enough memory\n");
                exit (EXIT_FAILURE);
              }
	    moduli = (long unsigned int *) malloc (nprimes *
                                                   sizeof(unsigned long int));
            if (moduli == NULL)
              {
                fprintf (stderr, "Error: not enough memory\n");
                exit (EXIT_FAILURE);
              }
	    len = 1;
	    sieve = (unsigned char *) malloc(len *
                                       sizeof(unsigned char)); /* len=1 here */
            if (sieve == NULL)
              {
                fprintf (stderr, "Error: not enough memory\n");
                exit (EXIT_FAILURE);
              }
	    offset = 5.0;
	    sieve[0] = 1; /* corresponding to 5 */
	    primes[0] = 3;
	    moduli[0] = 1; /* next odd multiple of 3 is 7, i.e. next to 5 */
	    current = -1;
	    return 3.0;
	  }
	else
	  {
	    unsigned long int i, p, j, ok;

	    i = nprimes;
	    nprimes *= 2;
	    primes = (unsigned *) realloc (primes, nprimes *
                                           sizeof(unsigned long int));
	    moduli = (unsigned long int *) realloc (moduli, nprimes *
                                                    sizeof(unsigned long int));
	    for (p = primes[i-1]; i < nprimes; i++)
	      {
		/* find next (odd) prime > p */
		do
		  {
		    for (p += 2, ok = 1, j = 0; (ok != 0) && (j < i); j++)
		      ok = p % primes[j];
		  }
		while (ok == 0);
		primes[i] = p;
		/* moduli[i] is the smallest m such that offset + 2*m = k*p */
		j = (unsigned long) fmod (offset, (double) p);
		j = (j == 0) ? j : p - j; /* -offset mod p */
		if ((j % 2) != 0)
		  j += p; /* ensure j is even */
		moduli[i] = j / 2;
	      }
	  }
      }

  /* now sieve for new primes */
  {
    long int i;
    unsigned long int j, p;
    
    for (i = 0; i < len; i++)
      sieve[i] = 1;
    for (j = 0; j < nprimes; j++)
      {
	p = primes[j];
	for (i = moduli[j]; i < len; i += p)
	  sieve[i] = 0;
	moduli[j] = i - len; /* for next sieving array */
      }
  }

  current = -1;
  while ((++current < len) && (sieve[current] == 0));

  if (current < len)
    return offset + 2.0 * (double) current;
  else
    abort ();
}
