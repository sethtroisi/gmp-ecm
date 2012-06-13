/* Random initialization for P-1 and P+1.

Copyright 2005, 2006, 2008 Paul Zimmermann, Alexander Kruppa, Dave Newman.

This file is part of the ECM Library.

The ECM Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The ECM Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the ECM Library; see the file COPYING.LIB.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#include <stdio.h>
#include <stdlib.h>

#ifdef OUTSIDE_LIBECM
# include "ecm-ecm.h"
#else
# include "ecm-impl.h"
#endif

#ifdef HAVE_UNISTD_H
# include <unistd.h> /* getpid */
#endif

#ifdef TIME_WITH_SYS_TIME
# include <sys/time.h>
# include <time.h>
#else
# if HAVE_SYS_TIME_H
#  include <sys/time.h>
# else
#  include <time.h>
# endif
#endif

#if defined (_MSC_VER) || defined (__MINGW32__)
# include <windows.h>
# include <wincrypt.h>
#endif

#if 0  /* dirty hack until outputf gets fixed */
#ifdef outputf
# undef outputf
# define outputf(x,y) printf(y)
#endif
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

/* Produces a random unsigned long value */

#if defined (_MSC_VER) || defined (__MINGW32__)
unsigned long 
get_random_ul (void)
{
  SYSTEMTIME tv;
  HCRYPTPROV Prov;

  if (CryptAcquireContext (&Prov, NULL, NULL, PROV_RSA_FULL,
    CRYPT_VERIFYCONTEXT))
    {
      int r;
      unsigned long rnd;
    
      r = CryptGenRandom (Prov, sizeof (unsigned long), (void *) &rnd);
      CryptReleaseContext (Prov, 0);
      if (r)
        {
/* warning: outputf is not exported from libecm */
#if !defined (OUTSIDE_LIBECM) && !defined(GPUECM)
          outputf (OUTPUT_DEVVERBOSE, "Got seed for RNG from CryptGenRandom\n");
#endif
          return rnd;
        }
    }
  
/* warning: outputf is not exported from libecm */
#if !defined (OUTSIDE_LIBECM) && !defined(GPUECM)
  outputf (OUTPUT_DEVVERBOSE, "Got seed for RNG from GetSystemTime\n");
#endif

  GetSystemTime (&tv);
  /* This gets us 27 bits of somewhat "random" data based on the time clock.
     It would probably do the program justice if a better random mixing was done
     in the non-MinGW get_random_ul if /dev/random does not exist */
  return ((tv.wHour<<22)+(tv.wMinute<<16)+(tv.wSecond<<10)+tv.wMilliseconds) ^
         ((tv.wMilliseconds<<17)+(tv.wMinute<<11)+(tv.wHour<<6)+tv.wSecond);
}

#else

unsigned long 
get_random_ul (void)
{
  FILE *rndfd;
  struct timeval tv;
  unsigned long t;

  /* Try /dev/urandom */
  rndfd = fopen ("/dev/urandom", "r");
  if (rndfd != NULL)
    {
      if (fread (&t, sizeof (unsigned long), 1, rndfd) == 1)
        {
/* warning: outputf is not exported from libecm */
#if !defined (OUTSIDE_LIBECM) && !defined(GPUECM)
          outputf (OUTPUT_DEVVERBOSE, "Got seed for RNG from /dev/urandom\n");
#endif
          fclose (rndfd);
          return t;
        }
      fclose (rndfd);
    }

#ifdef HAVE_GETTIMEOFDAY
  if (gettimeofday (&tv, NULL) == 0)
    {
/* warning: outputf is not exported from libecm */
#if !defined (OUTSIDE_LIBECM) && !defined(GPUECM)
      outputf (OUTPUT_DEVVERBOSE, "Got seed for RNG from gettimeofday()\n");
#endif
      return (unsigned long) tv.tv_sec + 
             (unsigned long) tv.tv_usec * 2147483629UL;
    }
#endif

/* warning: outputf is not exported from libecm */
#if !defined (OUTSIDE_LIBECM) && !defined(GPUECM)
  outputf (OUTPUT_DEVVERBOSE, "Got seed for RNG from time()+getpid()\n");
#endif

  /* Multiply one value by a large prime to get a bit of avalance effect */
  return (unsigned long) time (NULL) + 
         (unsigned long) getpid () * 2147483629UL;
}
#endif
