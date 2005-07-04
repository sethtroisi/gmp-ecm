/* Auxiliary routines for the ecm library.

  Copyright 2001, 2002, 2003, 2004, 2005 Paul Zimmermann and Alexander Kruppa.

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
#include <time.h>
#include <stdarg.h>
#include <gmp.h>
#include "ecm-impl.h"

#define VERBOSE __ECM(verbose)
static int VERBOSE = OUTPUT_NORMAL;

unsigned int
gcd (unsigned int a, unsigned int b)
{
  unsigned int t;

  while (b != 0)
    {
      t = a % b;
      a = b;
      b = t;
    }

  return a;
}

void 
mpz_sub_si (mpz_t r, mpz_t s, int i)
{
  if (i >= 0)
    mpz_sub_ui (r, s, (unsigned int) i);
  else
    mpz_add_ui (r, s, (unsigned int) (-i));
}

/* Divide RS by 3 */
void
mpz_divby3_1op (mpz_t RS)
{
  mp_size_t abssize = mpz_size (RS);
  mp_limb_t r;
  
  if (abssize == 0)
    return;
  
  r = mpn_divexact_by3 (RS->_mp_d, RS->_mp_d, abssize);
  ASSERT(r == 0);

  if (RS->_mp_d[abssize - 1] == 0)
    RS->_mp_size -= mpz_sgn (RS);
}

/* returns ceil(log(n)/log(2)) */
unsigned int
ceil_log2 (unsigned long n)
{
  unsigned int k = 0;

  ASSERT (n > 0);

  n--;
  while (n)
    {
      k++;
      n >>= 1;
    }

  return k;
}

/* cputime() gives the CPU time in microseconds */
#if defined (ANSIONLY) || defined (USG) || defined (__SVR4) || defined (_UNICOS) || defined(__hpux)

/* multiplier to get time in milliseconds */
#define TICK (1000.0 / (double) CLOCKS_PER_SEC)

unsigned int
cputime ()
{
  return (unsigned int) clock ();
}

#else	/* ANSIONLY and others */

#if defined (__MINGW32__) || defined (_MSC_VER)

#define TICK (1000.0 / (double) CLOCKS_PER_SEC)

static unsigned int
cputime_x (void)
{
  return (unsigned int) clock ();
}

#include <windows.h>

unsigned int
cputime ()
{
  static int First = 1;
  static LARGE_INTEGER PF;
  LARGE_INTEGER i;
  double d;

  if (First == 1)
  {
    First = 0;
    if (!QueryPerformanceFrequency (&PF))
        First = -1;
  }
  if (First == -1)
    return cputime_x ();

  QueryPerformanceCounter (&i);
  d = (double)*(__int64*)&i;
  d /= *(__int64*)&PF;
  d *= (double) CLOCKS_PER_SEC;

  return (unsigned int) d;
}

#else /* __MINGW32___ or VC */

#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>

#define TICK 1.0

unsigned int
cputime ()
{
  struct rusage rus;

  getrusage (RUSAGE_SELF, &rus);
  return rus.ru_utime.tv_sec * 1000 + rus.ru_utime.tv_usec / 1000;
}
#endif  /* __MINGW32___ or VC */

#endif /* ANSIONLY and others */

/* ellapsed time (in milliseconds) between st0 and st1 (values of cputime) */
unsigned int
elltime (unsigned int st0, unsigned int st1)
{
  st1 = st1 - st0; /* correct even if wrap around, as long as the time
		      difference is less than UINT_MAX */
  return (unsigned int) ((double) st1 * TICK);
}

int 
get_verbose ()
{
  return VERBOSE;
}

/* Tests if loglevel gets printed with the current verbose setting */

int 
test_verbose (int loglevel)
{
  return (loglevel <= VERBOSE);
}

void 
set_verbose (int v)
{
  VERBOSE = v;
}

int 
inc_verbose ()
{
  VERBOSE ++;
  return VERBOSE;
}

int
outputf (int loglevel, char *format, ...)
{
  va_list ap;
  int n = 0;
  
  va_start (ap, format);

  if (loglevel != OUTPUT_ERROR && loglevel <= VERBOSE)
    {
      n = gmp_vfprintf (ECM_STDOUT, format, ap);
      fflush (ECM_STDOUT);
    }
  else if (loglevel == OUTPUT_ERROR)
    n = gmp_vfprintf (ECM_STDERR, format, ap);
  
  va_end (ap);
  
  return n;
}
