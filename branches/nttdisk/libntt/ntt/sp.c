#include "sp.h"

sp_t X(sp_reciprocal)(sp_t p)
{
    /* integer reciprocal */

    #if SP_NUMB_BITS <= SP_TYPE_BITS - 2
    mp_limb_t shift = 2 * SP_NUMB_BITS + 1;
    #else
    mp_limb_t shift = 2 * SP_NUMB_BITS;
    #endif

#if SP_NUMB_BITS == 50 /* floating point reciprocal */

    return 1.0 / p;

#elif SP_TYPE_BITS == GMP_LIMB_BITS  /* use GMP functions */
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
  uint64_t r, s, t, e;

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

  if (x < SP_MIN)
    return 1;
  
  d = X(sp_reciprocal)(x);
  
  /* see http://mersenneforum.org/showthread.php?t=12209 */

#if SP_NUMB_BITS <= 32

  if (!sp_spp (sprp32_lookup[((uint32_t)x * 0x3AC69A35) >> 22], x, d))
    return 0;

#else

  if (!sp_spp (2, x, d) || 
      !sp_spp (325, x, d) || 
      !sp_spp (9375, x, d) ||
      !sp_spp (sprp64_lookup[(((uint64_t)x * 
                     3141592653589793239ULL) >> 42) & 0x7ff], x, d))
      return 0;
#endif 

  return 1;
}
