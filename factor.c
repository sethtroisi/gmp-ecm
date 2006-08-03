/* factor.c - public interface for libecm.

  Copyright 2005 Paul Zimmermann and Alexander Kruppa.

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

#include <stdio.h>
#include "ecm-impl.h"

void
ecm_init (ecm_params q)
{
  q->method = ECM_ECM; /* default method */
  MEMORY_TAG;
  mpz_init_set_ui (q->x, 0);
  mpz_init_set_ui (q->sigma, 0);
  q->sigma_is_A = 0;
  mpz_init_set_ui (q->go, 1);
  q->B1done = ECM_DEFAULT_B1_DONE;
  mpz_init_set_si (q->B2min, -1.0); /* default: B2min will be set to B1 */
  mpz_init_set_si (q->B2, ECM_DEFAULT_B2);
  q->k = ECM_DEFAULT_K;
  q->S = ECM_DEFAULT_S; /* automatic choice of polynomial */
  q->repr = ECM_MOD_DEFAULT; /* automatic choice of representation */
  q->verbose = 0; /* no output (default in library mode) */
  q->os = stdout; /* standard output */
  q->es = stderr; /* error output */
  q->chkfilename = NULL;
  q->TreeFilename = NULL;
  q->maxmem = 0.0;
  q->stage1time = 0.0;
  MEMORY_TAG;
  gmp_randinit_default (q->rng);
  MEMORY_TAG;
  gmp_randseed_ui (q->rng, get_random_ui ());
  MEMORY_UNTAG;
  q->use_ntt = 1;
  q->stop_asap = NULL;
}

void
ecm_clear (ecm_params q)
{
  mpz_clear (q->x);
  mpz_clear (q->sigma);
  mpz_clear (q->go);
  mpz_clear (q->B2min);
  mpz_clear (q->B2);
  gmp_randclear (q->rng);
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
    res = ecm (f, p->x, p->sigma, n, p->go, &(p->B1done), B1, p->B2min, p->B2, 
               1.0, p->k, p->S, p->verbose, p->repr, p->use_ntt, p->sigma_is_A,
	       p->os, p->es, p->chkfilename, p->TreeFilename, p->maxmem, 
	       p->stage1time, p->rng, p->stop_asap);
  else if (p->method == ECM_PM1)
    res = pm1 (f, p->x, n, p->go, &(p->B1done), B1, p->B2min, p->B2, 1.0,
               p->k, p->S, p->verbose, p->repr, p->use_ntt, p->os, p->es,
               p->chkfilename, p->TreeFilename, p->maxmem, p->rng, 
               p->stop_asap);
  else if (p->method == ECM_PP1)
    res = pp1 (f, p->x, n, p->go, &(p->B1done), B1, p->B2min, p->B2, 1.0,
               p->k, p->S, p->verbose, p->repr, p->use_ntt, p->os, p->es,
               p->chkfilename, p->TreeFilename, p->maxmem, p->rng, 
               p->stop_asap);
  else
    {
      fprintf (p->es, "Error, unknown method: %d\n", p->method);
      res = ECM_ERROR;
    }

  if (p_is_null)
    ecm_clear (q);

  return res;
}
