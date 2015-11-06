#include "libntt.h"

const nttinit_t master_init[] = {
#if GMP_LIMB_BITS == 32
  (const nttinit_t)&nttinit_sp30w32,
  (const nttinit_t)&nttinit_sp31w32,
  (const nttinit_t)&nttinit_sp50w32,
  (const nttinit_t)&nttinit_sp62w32,
#else
  (const nttinit_t)&nttinit_sp30w64,
  (const nttinit_t)&nttinit_sp31w64,
  (const nttinit_t)&nttinit_sp50w64,
  (const nttinit_t)&nttinit_sp62w64,
#endif
};

void 
nttwork_clear(nttwork_t nttwork)
{
  uint32_t i;

  if (nttwork == NULL)
    return;

  for (i = 0; i < nttwork->mpzspm_num; i++)
    {
      nttwork->nttinit[i]->mpzspm_clear(nttwork->mpzspm[i]);
    }

  free (nttwork->nttinit);
  free (nttwork->mpzspm);
  mpz_clear (nttwork->modulus);
  free (nttwork);
}

nttwork_t
nttwork_init(uint32_t max_ntt_size, mpz_t modulus, uint32_t interleaved,
		uint32_t * sp_sizes, uint32_t num_sp_sizes)
{
  uint32_t i, j;
  uint32_t max_mpzspm;
  uint32_t done = 0;
  mpz_t P, S;
  nttwork_t nttwork;

  nttwork = (nttwork_t) calloc (1, sizeof (__nttwork_struct));
  if (nttwork == NULL)
    return NULL;
  
  mpz_init_set_ui (P, 1);
  mpz_init_set_ui (S, 0);

  mpz_init_set(nttwork->modulus, modulus);
  max_mpzspm = 10;
  nttwork->max_ntt_size = max_ntt_size;
  nttwork->mpzspm_num = 0;
  nttwork->nttinit = (nttinit_t *) malloc (max_mpzspm * sizeof (nttinit_t));
  nttwork->mpzspm = (void **) malloc (max_mpzspm * sizeof (void *));
  if (nttwork->mpzspm == NULL || nttwork->nttinit == NULL)
    goto cleanup;

  for (i = 0; !done && i < num_sp_sizes; i++)
    {
      for (j = 0; j < sizeof(master_init) / sizeof(master_init[0]); j++)
	{
	  nttinit_t init = master_init[j];

	  if (init->sp_bits != sp_sizes[i])
	    continue;

    	  /* add this mpzspm */

	  if (nttwork->mpzspm_num == max_mpzspm)
	    {
	      nttinit_t *tmp1;
	      void **tmp2;

	      max_mpzspm *= 2;
	      tmp1 = (nttinit_t *)realloc(nttwork->nttinit, 
				    max_mpzspm * sizeof(nttinit_t));
	      if (tmp1 == NULL)
		goto cleanup;

	      tmp2 = (void **)realloc(nttwork->mpzspm, 
				    max_mpzspm * sizeof(void *));
	      if (tmp2 == NULL)
		{
		  free (tmp1);
		  goto cleanup;
		}

	      nttwork->nttinit = tmp1;
	      nttwork->mpzspm = tmp2;
	    }

	  nttwork->nttinit[nttwork->mpzspm_num] = init;
	  nttwork->mpzspm[nttwork->mpzspm_num] = init->mpzspm_init(
	      						max_ntt_size,
							modulus, P, S, 
							interleaved, &done);
	  if (nttwork->mpzspm[nttwork->mpzspm_num] != NULL)
	    nttwork->mpzspm_num++;

	  break;
	}
    }

cleanup:
  mpz_clear (P);
  mpz_clear (S);
  if (done)
      return nttwork;

  nttwork_clear(nttwork);
  return NULL;
}

void
nttwork_ntt_init(nttwork_t nttwork, nttplangroup_t *plans_in)
{
  uint32_t i;
  uint32_t ntt_size = 1;

  for (i = 0; i < plans_in[0].num_plans; i++)
      ntt_size *= plans_in[0].plans[i].codelet_size;

  nttwork->ntt_size = ntt_size;

  for (i = 0; i < nttwork->mpzspm_num; i++)
      nttwork->nttinit[i]->mpzspm_ntt_init(
			      nttwork->mpzspm[i], 
			      nttwork->ntt_size, 
			      nttwork->max_ntt_size, 
			      plans_in + i);
}

void
nttwork_ntt_run(nttwork_t nttwork, mpz_t * x)
{
  uint32_t i;

  for (i = 0; i < nttwork->mpzspm_num; i++)
    nttwork->nttinit[i]->mpzspm_ntt_run(nttwork->mpzspm[i], 
	  				x, nttwork->ntt_size);
}

void
nttwork_ntt_reset(nttwork_t nttwork)
{
  uint32_t i;

  for (i = 0; i < nttwork->mpzspm_num; i++)
    nttwork->nttinit[i]->mpzspm_ntt_reset(nttwork->mpzspm[i]);
}

void
nttwork_random(nttwork_t nttwork)
{
  uint32_t i;

  for (i = 0; i < nttwork->mpzspm_num; i++)
    nttwork->nttinit[i]->mpzspm_random(nttwork->mpzspm[i], 
					nttwork->ntt_size);
}

double
nttwork_ntt_test(nttwork_t nttwork, uint32_t verify)
{
  uint32_t i, j;

  if (verify)
    {
      for (i = 0; i < nttwork->mpzspm_num; i++)
	  nttwork->nttinit[i]->mpzspm_test(nttwork->mpzspm[i], 
	      				nttwork->ntt_size,
	      				nttwork->max_ntt_size);
      return 0;
    }
  else
    {
      uint64_t start, stop, elapsed;

      elapsed = (uint64_t)(-1);
      for (i = 0; i < 10; i++)
	{
	  start = read_clock();

	  for (j = 0; j < 10; j++)
	    nttwork_ntt_run(nttwork, NULL);

	  stop = read_clock();
	  if (stop - start < elapsed)
    	    elapsed = stop - start;
	}

      return (double)elapsed / 10 / nttwork->ntt_size;
    }
}


