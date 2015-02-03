/* Auxiliary arithmetic routines on unsigned long ints for the ecm library.

Copyright 2001, 2002, 2003, 2004, 2005, 2007, 2008 Paul Zimmermann and 
Alexander Kruppa.

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

#include "config.h"
#include "ecm-impl.h"

/* Returns the gcd of a and b */
unsigned long
gcd (unsigned long a, unsigned long b)
{
  unsigned long t;

  while (b != 0UL)
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
  unsigned long phi = 1UL, p;

  for (p = 2UL; p * p <= n; p += 2UL)
    {
      if (n % p == 0UL)
	{
	  phi *= p - 1UL;
	  n /= p;
	  while (n % p == 0UL)
	    {
	      phi *= p;
	      n /= p;
	    }
	}

      if (p == 2UL)
	p--;
    }

  /* now n is prime or 1 */

  return (n == 1UL) ? phi : phi * (n - 1UL);
}


/* returns ceil(log(n)/log(2)) */
unsigned int
ceil_log2 (unsigned long n)
{
  unsigned int k = 0;

  ASSERT (n > 0UL);

  n--;
  while (n)
    {
      k++;
      n >>= 1;
    }

  return k;
}

/* Returns the smallest prime factor of N. If N == 1, return 1. */
unsigned long
find_factor (const unsigned long N)
{
  unsigned long i;

  ASSERT_ALWAYS (N != 0UL);

  if (N == 1UL)
    return 1UL;

  if (N % 2UL == 0UL)
    return 2UL;

  for (i = 3UL; i*i <= N; i += 2UL)
    if (N % i == 0UL)
      return i;

  return N;
}
