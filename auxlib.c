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

/* need stdio.h and stdarg.h for gmp.h to declare gmp_vfprintf */
#include <stdio.h>
#include <stdarg.h>
#include <gmp.h>
#include "ecm-impl.h"

#if TIME_WITH_SYS_TIME
# include <sys/time.h>
# include <time.h>
#else
# if HAVE_SYS_TIME_H
#  include <sys/time.h>
# else
#  include <time.h>
# endif
#endif

#if HAVE_LIMITS_H
# include <limits.h>
#else
# ifndef ULONG_MAX
#  define LONG_MAX (__GMP_ULONG_MAX / 2)
# endif
#endif

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

/* returns Euler's totient phi function */
unsigned long
eulerphi (unsigned long n)
{
  unsigned long phi = 1, p;

  for (p = 2; p * p <= n; p += 2)
    {
      if (n % p == 0)
	{
	  phi *= p - 1;
	  n /= p;
	  while (n % p == 0)
	    {
	      phi *= p;
	      n /= p;
	    }
	}

      if (p == 2)
	p--;
    }

  /* now n is prime or 1 */

  return (n == 1) ? phi : phi * (n - 1);
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
  ASSERT (r == 0);

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

/* cputime () gives the elapsed time in milliseconds */

#if defined (_WIN32)
/* First case - GetProcessTimes () is the only known way of getting process
 * time (as opposed to calendar time) under mingw32 */

#include <windows.h>

long
cputime ()
{
  FILETIME lpCreationTime, lpExitTime, lpKernelTime, lpUserTime;
  ULARGE_INTEGER n;

  HANDLE hProcess = GetCurrentProcess();
  
  GetProcessTimes (hProcess, &lpCreationTime, &lpExitTime, &lpKernelTime,
      &lpUserTime);

  /* copy FILETIME to a ULARGE_INTEGER as recommended by MSDN docs */
  n.u.LowPart = lpUserTime.dwLowDateTime;
  n.u.HighPart = lpUserTime.dwHighDateTime;

  /* lpUserTime is in units of 100 ns. Return time in milliseconds */
  return (long) (n.QuadPart / 10000);
}

#elif defined (HAVE_GETRUSAGE)
/* Next case: getrusage () has higher resolution than clock () and so is
   preferred. */

#ifdef HAVE_SYS_TYPES_H
# include <sys/types.h>
#endif
#ifdef HAVE_SYS_RESOURCE_H
# include <sys/resource.h>
#endif

long
cputime ()
{
  struct rusage rus;

  getrusage (RUSAGE_SELF, &rus);
  /* This overflows a 32 bit signed int after 2147483s = 24.85 days */
  return rus.ru_utime.tv_sec * 1000L + rus.ru_utime.tv_usec / 1000L;
}

#else
/* Resort to clock (), which on some systems may return calendar time. */

long
cputime ()
{
  /* Return time in milliseconds */
  return (long) (clock () * (1000. / (double) CLOCKS_PER_SEC));
}

#endif /* defining cputime () */

/* ellapsed time (in milliseconds) between st0 and st1 (values of cputime) */
long
elltime (long st0, long st1)
{
  if (st1 >= st0)
    return st1 - st0;
  else
    {
      /* A wrap around can only really happen on a system where long int is 
         32 bit and where we use clock(). So we assume that there was exactly 
         one wrap-around which "swallowed" 
         LONG_MAX * (1000. / (double) CLOCKS_PER_SEC) milliseconds. */
      
      return st1 - st0 + (long)(LONG_MAX * (1000. / (double) CLOCKS_PER_SEC));
    }
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

  MEMORY_TAG; /* For gmp_*printf's temp allocs */
  if (loglevel != OUTPUT_ERROR && loglevel <= VERBOSE)
    {
      n = gmp_vfprintf (ECM_STDOUT, format, ap);
      fflush (ECM_STDOUT);
    }
  else if (loglevel == OUTPUT_ERROR)
    n = gmp_vfprintf (ECM_STDERR, format, ap);
  MEMORY_UNTAG;
  
  
  va_end (ap);
  
  return n;
}

void
writechkfile (char *chkfilename, int method, double p, mpmod_t modulus, 
              mpres_t A, mpres_t x, mpres_t z)
{
  FILE *chkfile;
  char *methodname;
  mpz_t t;

  outputf (OUTPUT_DEVVERBOSE, "Writing checkpoint to %s at p = %.0f\n",
           chkfilename, p);

  switch (method)
    {
    case ECM_ECM : methodname = "ECM"; break;
    case ECM_PM1 : methodname = "P-1"; break;
    case ECM_PP1 : methodname = "P+1"; break;
    default: 
      outputf (OUTPUT_ERROR, "writechkfile: Invalid method\n");
      return;
    }

  chkfile = fopen (chkfilename, "w");
  if (chkfile == NULL)
    {
      outputf (OUTPUT_ERROR, "Error opening checkpoint file %s\n", 
	       chkfilename);
      return;
    }

  mpz_init (t);

  gmp_fprintf (chkfile, "METHOD=%s; B1=%.0f; N=0x%Zx;", 
	       methodname, p, modulus->orig_modulus);
  mpres_get_z (t, x, modulus);
  gmp_fprintf (chkfile, " X=0x%Zx;", t);
  if (method == ECM_ECM)
    {
      mpres_get_z (t, z, modulus);
      gmp_fprintf (chkfile, " Z=0x%Zx;", t);
      mpres_get_z (t, A, modulus);
      gmp_fprintf (chkfile, " A=0x%Zx;", t);
    }
  fprintf (chkfile, "\n");
  mpz_clear (t);
  fflush (chkfile);
  fclose (chkfile);
}
