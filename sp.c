/* sp.c - "small prime" functions that don't need to be inlined

  Copyright 2005, 2008 Dave Newman and Jason Papadopoulos.

  The SP Library is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or (at your
  option) any later version.

  The SP Library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
  License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with the SP Library; see the file COPYING.LIB.  If not, write to
  the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
  MA 02110-1301, USA.
*/

#include <stdio.h> /* for stderr */
#include <stdlib.h>
#include "sp.h"

/* Test if m is a base "a" strong probable prime */

int
sp_spp (sp_t a, sp_t m, sp_t d)
{
  sp_t r, s, t, e;

  if (m == a)
    return 1;
	
  /* Set e * 2^s = m-1, e odd */
  for (s = 0, e = m - 1; !(e & 1); s++, e >>= 1);

  t = sp_pow (a, e, m, d);
	
  if (t == 1)
    return 1;
	
  for (r = 0; r < s; r++)
    {
      if (t == m - 1)
        return 1;

      t = sp_sqr (t, m, d);
    }
	
  return 0;
}

/* Test if x is a prime, return 1 if it is. Note this only works on sp's, 
   i.e. we need the top bit of x set */
int
sp_prime (sp_t x)
{
  sp_t d;

  if (!(x & 1))
    return 0;
  
  if (x < SP_MIN)
    return 1;
  
  sp_reciprocal (d, x);
  
  if (SP_NUMB_BITS <= 32)
    {
      /* 32-bit primality test
       * See http://primes.utm.edu/prove/prove2_3.html */
      
      if (!sp_spp (2, x, d) || !sp_spp (7, x, d) || !sp_spp (61, x, d))
        return 0;
    }
  else
    {
      ASSERT (SP_NUMB_BITS <= 64);
      /* 64-bit primality test
       * follows from results by Jaeschke, "On strong pseudoprimes to several
       * bases" Math. Comp. 61 (1993) p916 */
      
      if (!sp_spp (2, x, d) || !sp_spp (3, x, d) || !sp_spp (5, x, d)
        || !sp_spp (7, x, d) || !sp_spp (11, x, d) || !sp_spp (13, x, d)
        || !sp_spp (17, x, d) || ! sp_spp (19, x, d) || !sp_spp (23, x, d)
        || !sp_spp (29, x, d))
          return 0;
    }
  
  return 1;
}

#define CACHE_LINE_SIZE 64

void *
sp_aligned_malloc (size_t len)
{
  void *ptr, *aligned_ptr;
  unsigned long addr;

  ptr = malloc (len + CACHE_LINE_SIZE);
  if (ptr == NULL)
    {
      fprintf (stderr, "Cannot allocate memory in sp_aligned_malloc\n");
      exit (1);
    }

  addr = (unsigned long)ptr;				
  addr = CACHE_LINE_SIZE - (addr % CACHE_LINE_SIZE);
  aligned_ptr = (void *)((char *)ptr + addr);

  *( (void **)aligned_ptr - 1 ) = ptr;
  return aligned_ptr;
}

void
sp_aligned_free (void *newptr) 
{
  void *ptr;

  if (newptr == NULL) 
    return;
  ptr = *( (void **)newptr - 1 );
  free (ptr);
}
