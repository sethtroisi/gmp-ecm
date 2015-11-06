#include "sp.h"

void
X(spv_random) (spv_t x, spv_size_t len, sp_t m)
{
  spv_size_t i;
#if SP_NUMB_BITS == 50

  uint64_t t;
  for (i = 0; i < len; i++)
    {
      mpn_random ((mp_limb_t *)&t, sizeof(t) / sizeof(mp_limb_t));
      x[i] = t % (uint64_t)m;
    }
#else

  #if SP_TYPE_BITS == GMP_LIMB_BITS
  mpn_random ((mp_limb_t *)x, len);

  #elif SP_TYPE_BITS < GMP_LIMB_BITS
  mpn_random ((mp_limb_t *)x, len / 2);
  if (len % 2)
    {
      mp_limb_t t;
      mpn_random (&t, 1);
      x[len - 1] = (sp_t)t;
    }

  #else
  mpn_random ((mp_limb_t *)x, 2 * len);
  #endif

  for (i = 0; i < len; i++)
  #if SP_NUMB_BITS > SP_TYPE_BITS - 3
    while (x[i] >= m) 
      x[i] -= m;
  #else
    x[i] %= m;
  #endif
#endif
}
