#include "sp.h"

int
sp_spp (sp_t a, sp_t m, sp_t d)
{
  sp_t r, s, t, e;

  if (m == a)
    return 1;
	
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

/* note this only works on sp's, i.e. we need the top bit of x set */
int sp_prime (sp_t x)
{
  sp_t d;

  if (!(x & 1))
    return 0;
  
  if (x < (ULONG_MAX >> 1))
    return 1;
  
  invert_limb (d, x);
  
  if (!sp_spp (2, x, d))
    return 0;
  
  if (!sp_spp (7, x, d))
    return 0;
  
  if (!sp_spp (61, x, d))
    return 0;
  
  return 1;
}
