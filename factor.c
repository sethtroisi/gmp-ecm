#include <stdio.h>
#include "gmp.h"
#include "ecm.h"
#include "ecm-impl.h"

void
ecm_init (ecm_params q)
{
  q->method = ECM_ECM; /* default method */
  mpz_init_set_ui (q->x, 0);
  mpz_init_set_ui (q->sigma, 0);
  q->sigma_is_A = 0;
  mpz_init_set_ui (q->go, 1);
  q->B1done = ECM_DEFAULT_B1_DONE;
  mpz_init_set_si (q->B2min, -1.0); /* default: B2min will be set to B1 */
  mpz_init_set_si (q->B2, ECM_DEFAULT_B2);
  q->k = ECM_DEFAULT_K;
  q->S = ECM_DEFAULT_S; /* automatic choice of polynomial */
  q->repr = ECM_DEFAULT_REPR; /* automatic choice of representation */
  q->verbose = 0; /* no output (default in library mode) */
  q->os = stdout; /* standard output */
  q->es = stderr; /* error output */
  q->TreeFilename = NULL;
}

void
ecm_clear (ecm_params q)
{
  mpz_clear (q->x);
  mpz_clear (q->sigma);
  mpz_clear (q->go);
  mpz_clear (q->B2min);
  mpz_clear (q->B2);
}

/* returns ECM_FACTOR_FOUND, ECM_NO_FACTOR_FOUND, or ECM_ERROR */
int
ecm_factor (mpz_t f, mpz_t n, double B1, ecm_params p)
{
  int res; /* return value */
  int p_is_null;
  ecm_params q;

  if ((p_is_null = (p == NULL)))
    {
      p = q;
      ecm_init (q);
    }

  if (p->method == ECM_ECM)
    res = ecm (f, p->x, p->sigma, n, p->go, p->B1done, B1, p->B2min, p->B2, 1.0,
               p->k, p->S, p->verbose, p->repr, p->sigma_is_A, p->os, p->es,
               p->TreeFilename);
  else if (p->method == ECM_PM1)
    res = pm1 (f, p->x, n, p->go, p->B1done, B1, p->B2min, p->B2, 1.0,
               p->k, p->S, p->verbose, p->repr, p->os, p->es,
               p->TreeFilename);
  else if (p->method == ECM_PP1)
    res = pp1 (f, p->x, n, p->go, p->B1done, B1, p->B2min, p->B2, 1.0,
               p->k, p->S, p->verbose, p->repr, p->os, p->es,
               p->TreeFilename);
  else
    {
      fprintf (p->es, "Error, unknown method: %d\n", p->method);
      res = ECM_ERROR;
    }

  if (p_is_null)
    ecm_clear (q);

  return res;
}
