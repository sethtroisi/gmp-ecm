/* Dynamic Eratosthenes sieve.

 Copyright 2001, 2002, 2003, 2005 Paul Zimmermann and Alexander Kruppa.
 (Modified wrt GMP-ECM to use 'unsigned long' instead of 'double'.)

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

#include <stdio.h>
#include <stdlib.h>
#include "getprime.h"

/* provided for in cado.h, but we want getprime.c to be standalone */
#ifndef ASSERT
#define ASSERT(x)
#endif

/* This function returns successive odd primes, starting with 3.
  To perform a loop over all primes <= B1, do the following
  (compile this file with -DMAIN to count primes):

     for (p = 2; p <= B1; p = getprime (p))
        {
           ...
        }

     getprime (0);  { free the memory used by getprime() }

  It is slightly less efficient (1.5 to 2 times) than Dan Bernstein's
  primegen library (http://cr.yp.to/primegen.html), however it is
  fast enough for our usage here.

  FIXME: An MT-safe version would be nice (and easy).
*/

unsigned long
getprime (unsigned long pp)
{
 static unsigned long offset = 0;     /* offset for current primes */
 static long current = -1;            /* index of previous prime */
 static unsigned *primes = NULL;      /* small primes up to sqrt(p) */
 static unsigned long nprimes = 0;    /* length of primes[] */
 static unsigned char *sieve = NULL;  /* sieving table */
 static long len = 0;                 /* length of sieving table */
 static unsigned long *moduli = NULL; /* offset for small primes */

 if (pp == 0) /* free the tables, and reinitialize */
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

 /* the following complex block is equivalent to:
    while ((++current < len) && (sieve[current] == 0));
    but is faster.
 */
 {
   unsigned char *ptr = sieve + current;
   unsigned char *end = sieve + len;
   while ((++ptr < end) && (*ptr == 0));
   current = ptr - sieve;
 }

 if (current < len) /* most calls will end here */
   return offset + 2 * current;

 /* otherwise we have to sieve */
 offset += 2 * len;

 /* first enlarge sieving table if too small */
 if ((unsigned long) len * len < offset)
   {
     free (sieve);
     len *= 2;
     sieve = (unsigned char *) malloc (len * sizeof (unsigned char));
     /* assume this "small" malloc will not fail in normal usage */
     ASSERT(sieve != NULL);
   }

 /* now enlarge small prime table if too small */
 if ((nprimes == 0) ||
     (primes[nprimes - 1] * primes[nprimes - 1] < offset + len))
     {
       if (nprimes == 0) /* initialization */
         {
           nprimes = 1;
           primes = (unsigned *) malloc (nprimes * sizeof(unsigned long));
           /* assume this "small" malloc will not fail in normal usage */
           ASSERT(primes != NULL);
           moduli = (long unsigned int *) malloc (nprimes *
                                                  sizeof(unsigned long));
           /* assume this "small" malloc will not fail in normal usage */
           ASSERT(moduli != NULL);
           len = 1;
           sieve = (unsigned char *) malloc(len *
                                      sizeof(unsigned char)); /* len=1 here */
           /* assume this "small" malloc will not fail in normal usage */
           ASSERT(sieve != NULL);
           offset = 5.0;
           sieve[0] = 1; /* corresponding to 5 */
           primes[0] = 3;
           moduli[0] = 1; /* next odd multiple of 3 is 7, i.e. next to 5 */
           current = -1;
           return 3.0;
         }
       else
         {
           unsigned int i, p, j, ok;

           i = nprimes;
           nprimes *= 2;
           primes = (unsigned *) realloc (primes, nprimes *
                                          sizeof(unsigned long));
           moduli = (unsigned long*) realloc (moduli, nprimes *
                                                   sizeof(unsigned long));
           /* assume those "small" realloc's will not fail in normal usage */
           ASSERT(primes != NULL && moduli != NULL);
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
               j = offset % p;
               j = (j == 0) ? j : p - j; /* -offset mod p */
               if ((j % 2) != 0)
                 j += p; /* ensure j is even */
               moduli[i] = j / 2;
             }
         }
     }

 /* now sieve for new primes */
 {
   long i;
   unsigned long j, p;

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

 ASSERT(current < len); /* otherwise we found a prime gap >= sqrt(x) around x */
 return offset + 2 * current;
}

#ifdef MAIN
int
main (int argc, char *argv[])
{
 unsigned long p, B;
 unsigned long pi = 0;

 if (argc != 2)
   {
     fprintf (stderr, "Usage: getprime <bound>\n");
     exit (EXIT_FAILURE);
   }

 B = atoi (argv[1]);

 for (pi = 0, p = 2; p <= B; p = getprime (p), pi++);
 printf ("pi(%lu)=%lu\n", B, pi);

 getprime (0); /* free the tables */

 return 0;
}
#endif
