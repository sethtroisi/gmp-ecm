#include <stdio.h>
#include "gmp.h"
#include "ecm.h"
#include "ecm-impl.h"

/* returns ECM_FACTOR_FOUND, ECM_NO_FACTOR_FOUND, or ECM_ERROR */
int
ecm_factor (mpz_t f, mpz_t n, double B1, ecm_params p)
{
  int res = ECM_NO_FACTOR_FOUND; /* return value */
  int p_is_null;

  if (p_is_null = (p == NULL))
    {
      p = malloc (sizeof(ecm_params));
    }

  if (p_is_null)
    free (p);

  return res;
}
