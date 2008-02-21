/* Dynamic Eratosthenes sieve.

  Copyright 2001, 2002, 2003, 2005 Paul Zimmermann and Alexander Kruppa.

  This file is part of the ECM Library.

  The ECM Library is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or (at your
  option) any later version.

  The ECM Library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
  License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with the ECM Library; see the file COPYING.LIB.  If not, write to
  the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
  MA 02110-1301, USA.
*/

#include <math.h>
#include <stdlib.h>

#ifdef OUTSIDE_LIBECM
# include "ecm-ecm.h"
#else
# include "ecm-impl.h"
#endif

/* This function returns successive odd primes, starting with 3.
   To perform a loop over all primes <= B1, do the following
   (compile this file with -DMAIN to count primes):

      for (p = 2.0; p <= B1; p = getprime ())
         {
            ...
         }

   It is slightly less efficient (1.5 to 2 times) than Dan Bernstein's
   primegen library (http://cr.yp.to/primegen.html), however it is
   fast enough for our usage here.
*/

/* We allow primes up to 2^53. This means len, the primes in primes[] etc.
   will stay well below 2^32. */

static double offset = 0.0; /* offset for current primes, must be 0 or odd */
static int current = -1; /* index of previous prime */
static unsigned int *primes = NULL; /* table of small primes up to sqrt(p) */
static unsigned int nprimes = 0; /* length of primes[] */
static unsigned char *sieve = NULL; /* sieving table */
static int len = 0; /* length of sieving table, WITHOUT sentinel */
static unsigned int *moduli = NULL; /* offset for small primes, 
				       moduli[i] = offset mod primes[i] */

/* sieve[i] == 1 if offset+2*i is a prime, otherwise sieve[i] == 0.
   sieve has len + 1 bytes allocated, the last byte is always 1 (a sentinel). 
   This allows us avoid testing for the array end in the loop that looks
   for the next prime in sieve[]. */
/* The last prime returned by getprime is offset + 2*current */
/* primes[] contains small primes to needed to sieve out composites in 
   sieve, i.e. all primes <= sqrt(offset + 2 * (len - 1)). 
   moduli[i] contains the smallest k so that offset+2*(len+k) is divisible
   by primes[i], i.e. after advancing the sieve array by len, 
   sieve[moduli[i]] is divisible by primes[i].
*/

void
getprime_clear ()
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
}

/* For p > 1, return 1 if p is prime and 0 if p is not prime.
   Requires that all primes <= sqrt(p) are in *primes */

static int
isprime_ui (unsigned int p, unsigned int *primes)
{
  int i;
  
  for (i = 0; primes[i] * primes[i] <= p; i++)
    if (p % primes[i] == 0)
      return 0;
  
  return 1;
}

double
getprime ()
{
  /* the following complex block is equivalent to:
     while ((++current < len) && (sieve[current] == 0));
     but is faster.
  */
  if (len > 0L)
  {
    unsigned char *ptr = sieve + current;
    while (*(++ptr) == 0);
    current = ptr - sieve;
  }
  else
    current = len;

  if (current < len) /* most calls will end here */
    return offset + 2.0 * (double) current;

  /* otherwise we have to advance the sieve */
  offset += 2.0 * (double) len;

  /* first enlarge sieving table if too small */
  if ((double) len * (double) len < offset && len > 0)
    {
      free (sieve);
      len *= 2;
      sieve = (unsigned char *) malloc ((len + 1) * sizeof (unsigned char));
      /* assume this "small" malloc will not fail in normal usage */
      ASSERT(sieve != NULL);
    }

  /* now enlarge small prime table if too small */
  if ((nprimes == 0) || (primes[nprimes-1] < sqrt(offset + 2*len)))
    {
      if (nprimes == 0) /* initialization */
	{
	  nprimes = 1;
	  primes = (unsigned int *) malloc (nprimes * sizeof(unsigned int));
	  /* assume this "small" malloc will not fail in normal usage */
	  ASSERT(primes != NULL);
	  moduli = (unsigned int *) malloc (nprimes * sizeof(unsigned int));
	  /* assume this "small" malloc will not fail in normal usage */
	  ASSERT(moduli != NULL);
	  len = 1;
	  sieve = (unsigned char *) malloc((len + 1) *
					   sizeof(unsigned char)); /* len=1 here */
	  /* assume this "small" malloc will not fail in normal usage */
	  ASSERT(sieve != NULL);
	  offset = 5.0;
	  sieve[0] = 1; /* corresponding to 5 */
	  sieve[1] = 1; /* place the sentinel */
	  primes[0] = 3;
	  moduli[0] = 1; /* After we advance sieve[], sieve[0] will 
			    correspond to 7 and sieve[1] to 9, which is
			    the smallest odd multiple of 3 */
	  current = -1;
	  return 3.0;
	}
      else
	{
	  /* extend the existing table of small primes */
	  unsigned int i, j;
	  
	  i = nprimes;
	  nprimes *= 2;
	  primes = (unsigned int *) realloc (primes, nprimes *
					     sizeof(unsigned int));
	  moduli = (unsigned int *) realloc (moduli, nprimes *
					     sizeof(unsigned int));
	  /* assume those "small" realloc's will not fail in normal usage */
	  ASSERT_ALWAYS(primes != NULL && moduli != NULL);
	  for (; i < nprimes; i++)
	    {
	      unsigned int p;
	      /* find next (odd) prime */
	      for (p = primes[i - 1] + 2; !isprime_ui (p, primes); p += 2);
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
    int i, p;
    unsigned int j;
        
    /* Set sieve (including sentinel at the end) to 1 */
    for (i = 0; i < len + 1; i++)
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
  while (sieve[++current] == 0);

  ASSERT(current < len); /* otherwise we found a prime gap >= sqrt(x) around x */
  return offset + 2.0 * (double) current;
}


/* Skips forward or backward in the sieve so that the next call to getprime
   returns the smallest prime >= pp */

void
getprime_seek (double pp)
{
  int i, p;
  unsigned int j;

  if (pp <= 3.)
    {
      getprime_clear ();
      return;
    }

  offset = floor (pp / 2.) * 2. + 1.; /* make sure offset is odd */

  /* Choose a large enough sieve array length */
  for (i = 2; (double) i * (double) i < offset; i *= 2);

  /* Now allocate sieving table */
  if (len > 0)
    free (sieve);
  len = i;
  sieve = (unsigned char *) malloc ((len + 1) * sizeof (unsigned char));
  /* assume this "small" malloc will not fail in normal usage */
  ASSERT_ALWAYS(sieve != NULL);

  j = 1;
  /* Find out how many small odd primes we'll need */
  for (p = 5; (double)p*(double)p <= offset + (double)(2*len); p += 2)
    {
      for (i = 3; i*i <= p && p % i != 0; i += 2);
      if (i*i <= p)
	continue;
      if ((double)p*(double)p < offset + (double)len)
	j++;
    }

  /* Allocate memory for small primes */
  if (nprimes != 0)
    {
      free (primes);
      free (moduli);
    }
  nprimes = j;
  primes = (unsigned int *) malloc (nprimes * sizeof(unsigned int));
  moduli = (unsigned int *) malloc (nprimes * sizeof(unsigned int));
  ASSERT_ALWAYS(primes != NULL && moduli != NULL);
  
  /* Fill small primes and moduli arrays */
  for (p = 3, j = 0; j < nprimes; p += 2)
    {
      for (i = 3; i*i <= p && p % i != 0; i += 2);
      if (i*i <= p)
	continue;
      primes[j] = p;
      i = (unsigned int) fmod (offset, (double)p);
      i = (i == 0) ? i : p - i; /* -offset mod p */
      if (i % 2 != 0)
	i += p; /* ensure i is even */
      moduli[j] = i / 2;
      j++;
    }
  
  /* now sieve for new primes */
  for (i = 0; i < len + 1; i++)
    sieve[i] = 1;
  
  for (j = 0; j < nprimes; j++)
    {
      p = primes[j];
      for (i = moduli[j]; i < len; i += p)
	sieve[i] = 0;
      moduli[j] = i - len; /* for next sieving array */
    }

  current = -1;
}


#ifdef MAIN
int
main (int argc, char *argv[])
{
  double p, B1, B2;
  unsigned long pi = 0;

  if (argc != 3)
    {
      fprintf (stderr, "Usage: getprime <bound> <bound>\n");
      exit (EXIT_FAILURE);
    }

  B1 = atof (argv[1]);
  B2 = atof (argv[2]);

  if (B1 > 0.)
    getprime_seek (B1);

  p = 0;
  if (B1 <= 2)
    {
      printf("2\n");
      pi++;
    }
  
  for (p = getprime (); p <= B2; p = getprime (), pi++)
    printf("%1.0f\n", p);
  //  printf ("pi(%1.0f) - pi(%1.0f - 1) = %lu\n", B1, B2, pi);

  getprime_clear ();

  return 0;
}
#endif
