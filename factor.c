/* factor.c - public interface for libecm.

Copyright 2005, 2006, 2007, 2009, 2011 Paul Zimmermann, Alexander Kruppa,
David Cleaver, Cyril Bouvier.

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

#include <stdio.h>
#include <math.h>
#include "ecm-impl.h"
#include "ecm-gpu.h"


const char *
ecm_version ()
{
  static const char *version = ECM_VERSION;
  return version;
}

void
ecm_init (ecm_params q)
{
  __ell_curve_struct *ptrE = (__ell_curve_struct *) malloc(sizeof(__ell_curve_struct));

  q->method = ECM_ECM; /* default method */
  mpz_init_set_ui (q->x, 0);
  mpz_init_set_ui (q->y, 0);
  mpz_init_set_ui (q->sigma, 0);
  q->sigma_is_A = 0;
  mpz_init_set_ui (ptrE->a1, 0);
  mpz_init_set_ui (ptrE->a3, 0);
  mpz_init_set_ui (ptrE->a2, 0);
  mpz_init_set_ui (ptrE->a4, 0);
  mpz_init_set_ui (ptrE->a6, 0);
  ptrE->type = ECM_EC_TYPE_MONTGOMERY;
  ptrE->disc = 0;
  mpz_init_set_ui (ptrE->sq[0], 1);
  q->E = ptrE;
  q->param = ECM_PARAM_DEFAULT;
  mpz_init_set_ui (q->go, 1);
  q->B1done = ECM_DEFAULT_B1_DONE;
  mpz_init_set_si (q->B2min, -1.0); /* default: B2min will be set to B1 */
  mpz_init_set_si (q->B2, ECM_DEFAULT_B2);
  q->k = ECM_DEFAULT_K;
  q->S = ECM_DEFAULT_S; /* automatic choice of polynomial */
  q->repr = ECM_MOD_DEFAULT; /* automatic choice of representation */
  q->nobase2step2 = 0; /* continue special base 2 code in ecm step 2, if used */
  q->verbose = 0; /* no output (default in library mode) */
  q->os = stdout; /* standard output */
  q->es = stderr; /* error output */
  q->chkfilename = NULL;
  q->TreeFilename = NULL;
  q->maxmem = 0.0;
  q->stage1time = 0.0;
  gmp_randinit_default (q->rng);
  mpz_set_ui (q->rng->_mp_seed, 0); /* trick to tell that the random number
                                       generator has not been initialized */
  q->use_ntt = 1;
  q->stop_asap = NULL;
  q->batch_last_B1_used = 1.0;
  mpz_init_set_ui (q->batch_s, 1);
  q->gpu = 0; /* no gpu by default in library mode */
  q->gpu_device = -1; 
  q->gpu_device_init = 0; 
  q->gpu_number_of_curves = 0; 
  q->gw_k = 0.0;
  q->gw_b = 0;
  q->gw_n = 0;
  q->gw_c = 0;
}

void
ecm_clear (ecm_params q)
{
  mpz_clear (q->x);
  mpz_clear (q->y);
  mpz_clear (q->sigma);
  mpz_clear (q->go);
  mpz_clear (q->B2min);
  mpz_clear (q->B2);
  gmp_randclear (q->rng);
  mpz_clear (q->batch_s);
  mpz_clear (q->E->a1);
  mpz_clear (q->E->a3);
  mpz_clear (q->E->a2);
  mpz_clear (q->E->a4);
  mpz_clear (q->E->a6);
  mpz_clear (q->E->sq[0]);
  free (q->E);
}

/* returns ECM_FACTOR_FOUND, ECM_NO_FACTOR_FOUND, or ECM_ERROR */
int
ecm_factor (mpz_t f, mpz_t n, double B1, ecm_params p0)
{
  int res; /* return value */
  ecm_params q;
  ecm_params_ptr p;

  if (mpz_cmp_ui (n, 0) <= 0)
    {
      fprintf ((p0 == NULL) ? stderr : p0->es,
               "Error, n should be positive.\n");
      return ECM_ERROR;
    }
  else if (mpz_cmp_ui (n, 1) == 0)
    {
      mpz_set_ui (f, 1);
      return ECM_FACTOR_FOUND_STEP1;
    }
  else if (mpz_divisible_2exp_p (n, 1))
    {
      mpz_set_ui (f, 2);
      return ECM_FACTOR_FOUND_STEP1;
    }
  
  if (p0 == NULL)
    {
      p = q;
      ecm_init (q);
    }
  else
    p = p0;

  if (p->method == ECM_ECM)
    {
#ifdef WITH_GPU
      if (p->gpu == 0)
        {
#endif
            res = ecm (f, p->x, p->y, p->param, p->sigma, n, p->go,
		       &(p->B1done),
                       B1, p->B2min, p->B2, p->k, p->S, p->verbose,
                       p->repr, p->nobase2step2, p->use_ntt, 
		       p->sigma_is_A, p->E,
                       p->os, p->es, p->chkfilename, p->TreeFilename, p->maxmem,
                       p->stage1time, p->rng, p->stop_asap, p->batch_s,
                       &(p->batch_last_B1_used), p->gw_k, p->gw_b, p->gw_n,
                       p->gw_c);
#ifdef WITH_GPU
        }
      else
        {
          res = gpu_ecm (f, p->x, p->param, p->sigma, n, p->go,
                         &(p->B1done), B1, p->B2min, p->B2, p->k,
                         p->S, p->verbose, p->repr, p->nobase2step2, 
                         p->use_ntt, p->sigma_is_A, p->os, p->es,
                         p->chkfilename, p->TreeFilename, p->maxmem,
                         p->stop_asap, p->batch_s, &(p->batch_last_B1_used),
                         p->gpu_device, &(p->gpu_device_init),
                         &(p->gpu_number_of_curves));
        }
#endif
    }
  else if (p->method == ECM_PM1)
    res = pm1 (f, p->x, n, p->go, &(p->B1done), B1, p->B2min, p->B2,
               p->k, p->verbose, p->repr, p->use_ntt, p->os, p->es,
               p->chkfilename, p->TreeFilename, p->maxmem, p->rng,
               p->stop_asap);
  else if (p->method == ECM_PP1)
    res = pp1 (f, p->x, n, p->go, &(p->B1done), B1, p->B2min, p->B2,
               p->k, p->verbose, p->repr, p->use_ntt, p->os, p->es,
               p->chkfilename, p->TreeFilename, p->maxmem, p->rng,
               p->stop_asap);
  else
    {
      fprintf (p->es, "Error, unknown method: %d\n", p->method);
      res = ECM_ERROR;
    }

  if (p0 == NULL)
    ecm_clear (q);

  return res;
}
