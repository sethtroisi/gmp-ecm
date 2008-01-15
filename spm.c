/* spm.c - "small prime modulus" functions to precompute an inverse and a
   primitive root for a small prime

  Copyright 2005 Dave Newman.

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

#include <stdlib.h>
#include "sp.h"

/* Returns the exponent of $q$ in the factorisation of $n$ */
static int
exponent (const sp_t q, sp_t n)
{
  int i;
  for (i = 0; n % q == (sp_t) 0; i++, n /= q);
  return i;
}

/* Returns i so that ord(a) = q^i. This assumes that ord(a) is indeed
   a low power of q. */
static int
ordpow (const sp_t q, sp_t a, const sp_t sp, const sp_t mul_c)
{
  int i = 0;
  for (i = 0; a != (sp_t) 1; i++, a = sp_pow (a, q, sp, mul_c));
  return i;
}

/* Compute some constants, including a primitive n'th root of unity. */
spm_t
spm_init (spv_size_t n, sp_t sp)
{
  sp_t a, b, bd, sc;
  spv_size_t q, nc, d;
  spm_t spm = (spm_t) malloc (sizeof (__spm_struct));

  ASSERT (sp % (sp_t) n == (sp_t) 1);

  spm->sp = sp;
  invert_limb (spm->mul_c, sp);

  /* find an $n$-th primitive root $a$ of unity $(mod sp)$. */

  /* Construct a $b$ whose order $(mod sp)$ is equal to $n$.
     We try different a values and test if the exponent of $q$ in $ord(a)$
     is at least as large as in $n$. If it isn't, we move to another $a$.
     If it is, we optionally exponentiate to make the exponents equal and
     test for the remaining $q$'s.
     We assume that the largest prime dividing $n$ is very small, 
     so no optimizations in factoring n are made. */
  a = 2;
  b = a;
  nc = n; /* nc is remaining cofactor of n */
  q = 2;
  sc = sp - 1;
#ifdef PARI
  printf ("/* spm_init */ n = %lu; sp = %lu; /* PARI */\n", n, sp);
  printf ("exponent(a,b) = {local(i); while(b%%a == 0,i++;b/=a); "
	  "return(i)} /* PARI */\n");
#endif
  for ( ; nc != (spv_size_t) 1; q++)
    {
      if (nc % q == (spv_size_t) 0)
        {
	  const int k = exponent (q, n); /* q^k || n */
	  int l;
#ifdef PARI
	  printf ("exponent(%lu, n) == %d /* PARI */\n", q, k);
#endif
          /* Remove all factors of $q$ from $sp-1$ */
          for (d = sp - 1; d % q == (spv_size_t) 0; d /= q);
	  bd = sp_pow (b, d, sp, spm->mul_c);
          /* Now ord(bd) = q^l, q^l || ord(a) */
	  l = ordpow (q, bd, sp, spm->mul_c);
#ifdef PARI
	  printf ("exponent(%lu, znorder(Mod(%lu, sp))) == %d /* PARI */\n", 
		  q, b, l);
#endif
          if (l < k)
            {
              /* No good, q appears in ord(a) in a lower power than in n. 
		 Try next $a$ */
              a++;
              b = a;
              nc = n;
              q = 1; /* Loop increment following "continue" will make q=2 */
	      sc = sp - 1;
              continue;
            }
          else
            {
	      /* Reduce the exponent of $q$ in $ord(b)$ until is it 
		 equal to that in $n$ */
	      for ( ; l > k; l--)
		{
#ifdef PARI
		  printf ("Exponentiating %lu by %lu\n", b, q);
#endif
		  b = sp_pow (b, q, sp, spm->mul_c);
		}
#ifdef PARI
	      printf ("New b = %lu\n", b);
#endif
            }
	  do {nc /= q;} while (nc % q == 0); /* Divide out all q from nc */
	  while (sc % q == (sp_t) 0) /* Divide out all q from sc */
	    sc /= q;
        }
    }
  
  b = sp_pow (b, sc, sp, spm->mul_c);
#ifdef PARI
  printf ("znorder(Mod(%lu, sp)) == n /* PARI */\n", b, sp, n);
#endif

  /* turn this into a primitive n'th root of unity mod p */
  spm->prim_root = b;
  spm->inv_prim_root = sp_inv (spm->prim_root, sp, spm->mul_c);

  return spm;
}

void
spm_clear (spm_t spm)
{
  free (spm);
}
