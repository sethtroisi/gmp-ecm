#include "sp.h"

sp_t X(sp_reciprocal)(sp_t p)
{
    /* integer reciprocal */

    #if SP_NUMB_BITS <= SP_TYPE_BITS - 2
    mp_limb_t shift = 2 * SP_NUMB_BITS + 1;
    #else
    mp_limb_t shift = 2 * SP_NUMB_BITS;
    #endif

#if SP_TYPE_BITS == GMP_LIMB_BITS  /* use GMP functions */
    mp_limb_t recip, dummy;

    udiv_qrnnd (recip, dummy,
		(mp_limb_t) 1 << (shift - SP_TYPE_BITS), 0, p);
    return recip;

#elif SP_TYPE_BITS < GMP_LIMB_BITS  /* ordinary division */

    return ((mp_limb_t)1 << shift) / p;

#else  /* worst case: bit-at-a-time */

    sp_t r = (sp_t)1 << (SP_NUMB_BITS - 1);
    sp_t q = 0;
    mp_limb_t i;

    for (i = 0; i < shift + 1 - SP_NUMB_BITS; i++)
      {
	q += q;
	r += r;
	if (r >= p)
	{
	  r -= p;
	  q |= 1;
	}
      }
    return q;

#endif
}

/* Test if m is a base "a" strong probable prime */

static int
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

/* Test if x is a prime, return 1 if it is */

int
X(sp_prime)(sp_t x)
{
  sp_t d;

  if (!(x & 1))
    return 0;
  
  if (x < SP_MIN)
    return 1;
  
  d = X(sp_reciprocal)(x);
  
#if SP_NUMB_BITS <= 32
  /* 32-bit primality test
   * See http://primes.utm.edu/prove/prove2_3.html */
  
  if (!sp_spp (2, x, d) || !sp_spp (7, x, d) || !sp_spp (61, x, d))
    return 0;

#else

  ASSERT (SP_NUMB_BITS <= 64);
  /* 64-bit primality test
   * follows from results by Jaeschke, "On strong pseudoprimes to several
   * bases" Math. Comp. 61 (1993) p916 */
  
  if (!sp_spp (2, x, d) || !sp_spp (3, x, d) || !sp_spp (5, x, d)
    || !sp_spp (7, x, d) || !sp_spp (11, x, d) || !sp_spp (13, x, d)
    || !sp_spp (17, x, d) || ! sp_spp (19, x, d) || !sp_spp (23, x, d)
    || !sp_spp (29, x, d))
      return 0;
#endif 

  return 1;
}
