#include <stdlib.h>
#include "ntt-impl.h"

extern const nttconfig_t ntt3_config;
extern const nttconfig_t ntt4_config;
extern const nttconfig_t ntt5_config;
extern const nttconfig_t ntt7_config;
extern const nttconfig_t ntt8_config;
extern const nttconfig_t ntt9_config;
extern const nttconfig_t ntt15_config;

static const nttconfig_t * ntt_config[] = 
{
  &ntt3_config,
  &ntt4_config,
  &ntt5_config,
  &ntt7_config,
  &ntt8_config,
  &ntt9_config,
  &ntt15_config,
};

#define NUM_CODELETS (sizeof(ntt_config) / sizeof(ntt_config[0]))

/*-------------------------------------------------------------------------*/
void * ntt_init(sp_t size, sp_t primroot, sp_t p, sp_t recip)
{
  uint32_t i;
  uint32_t num_const;
  spv_t curr_const;
  nttdata_t *d;

  d = (nttdata_t *)calloc(1, sizeof(nttdata_t));
  if (d == NULL)
    return d;

  for (i = num_const = 0; i < NUM_CODELETS; i++)
    {
      const nttconfig_t *c = ntt_config[i];

      if (size % c->size != 0)
	continue;

      num_const += c->get_num_ntt_const();
      d->num_codelets++;
    }

  if (num_const == 0)
    goto error_free;

  d->codelets = (codelet_data_t *)malloc(d->num_codelets * 
		  			sizeof(codelet_data_t));
  if (d->codelets == NULL)
    goto error_free;

  d->codelet_const = curr_const = (spv_t)malloc(num_const * sizeof(sp_t));
  if (d->codelet_const == NULL)
    goto error_free;

  for (i = 0; i < NUM_CODELETS; i++)
    {
      const nttconfig_t *c = ntt_config[i];

      if (size % c->size != 0)
	continue;

      d->codelets[i].config = c;
      d->codelets[i].ntt_const = curr_const;
      c->nttdata_init(curr_const, p, recip, primroot, size);
      curr_const += c->get_num_ntt_const();
    }

  return d;

error_free:
  ntt_free(d);
  return NULL;
}

/*-------------------------------------------------------------------------*/
void ntt_free(void *data)
{
  nttdata_t *d = (nttdata_t *)data;

  if (d == NULL)
    return;

  free(d->codelets);
  free(d->codelet_const);
  free(d);
}


