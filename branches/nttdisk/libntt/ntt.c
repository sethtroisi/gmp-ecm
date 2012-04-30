#include "ntt-impl.h"

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

void * ntt_init(sp_t size, sp_t primroot, sp_t p, sp_t d)
{
}

void ntt_free(void *data)
{
  nttdata_t *d = (nttdata_t *)data;

  if (d == NULL)
    return;

  free(d->codelets);
  free(d->passes);
  free(d);
}


