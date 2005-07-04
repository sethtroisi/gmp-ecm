/* Random initialization for P-1 and P+1.

  Copyright 2005 Paul Zimmermann and Alexander Kruppa.

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
  the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
  MA 02111-1307, USA.
*/

#include <stdio.h>
#include <stdlib.h>
#if !defined (_MSC_VER)
#include <unistd.h> /* getpid */
#include <sys/time.h> /* gettimeofday */
#endif
#include <time.h>

#include <gmp.h>
#ifdef OUTSIDE_LIBECM
#include "ecm-ecm.h"
#else
#include "ecm-impl.h"
#endif

/* put in 'a' a valid random seed for P-1, i.e. gcd(a,n)=1 and a <> {-1,1} */
void
pm1_random_seed (mpz_t a, mpz_t n, gmp_randstate_t randstate)
{
  mpz_t q;

  mpz_init (q);
  do
    {
      mpz_urandomb (a, randstate, 32);
      mpz_gcd (q, a, n);
    }
  while (mpz_cmp_ui (q, 1) != 0 || mpz_cmp_ui (a, 1) == 0 ||
         mpz_cmp_si (a, -1) == 0);
  mpz_clear (q);
}

/* put in seed a valid random seed for P+1 */
void
pp1_random_seed (mpz_t seed, mpz_t n, gmp_randstate_t randstate)
{
  mpz_t q;

  /* need gcd(p^2-4, n) = 1. */
  mpz_init (q);
  do
    {
      mpz_urandomb (q, randstate, 32);
      mpz_add_ui (q, q, 1);
      mpz_set (seed, q);
      mpz_mul (q, q, q);
      mpz_sub_ui (q, q, 4);
      mpz_gcd (q, q, n);
    }
  while (mpz_cmp_ui (q, 1) != 0);
  mpz_clear (q);
}

/* Produces a random unsigned int value */
#if defined (_MSC_VER) || defined (__MINGW32__)
#include <windows.h>
unsigned int 
get_random_ui (void)
{
  SYSTEMTIME tv;
  GetSystemTime(&tv);
  /* This gets us 27 bits of somewhat "random" data based on the time clock.
     It would probably do the program justice if a better random mixing was done
     in the non-MinGW get_random_ui if /dev/random does not exist */
  return ((tv.wHour<<22)+(tv.wMinute<<16)+(tv.wSecond<<10)+tv.wMilliseconds) ^
         ((tv.wMilliseconds<<17)+(tv.wMinute<<11)+(tv.wHour<<6)+tv.wSecond);
}
#else
unsigned int 
get_random_ui (void)
{
  FILE *rndfd;
  struct timeval tv;
  unsigned int t;

  /* Try /dev/urandom */
  rndfd = fopen ("/dev/urandom", "r");
  if (rndfd != NULL)
    {
      if (fread (&t, sizeof(unsigned int), 1, rndfd) == 1)
        {
#ifndef OUTSIDE_LIBECM /* warning: outputf is not exported from libecm */
          outputf (OUTPUT_DEVVERBOSE, "Got seed for RNG from /dev/urandom\n");
#endif
          fclose (rndfd);
          return t;
        }
      fclose (rndfd);
    }

  if (gettimeofday (&tv, NULL) == 0)
    {
#ifndef OUTSIDE_LIBECM
      outputf (OUTPUT_DEVVERBOSE, "Got seed for RNG from gettimeofday()\n");
#endif
      return tv.tv_sec + tv.tv_usec;
    }

#ifndef OUTSIDE_LIBECM
  outputf (OUTPUT_DEVVERBOSE, "Got seed for RNG from time()+getpid()\n");
#endif

  return time (NULL) + getpid ();
}
#endif


