#include "ntt-impl-scalar.h"

/*-------------------------------------------------------------------------*/
static const nttconfig_t * ntt_config[] = 
{
  &X(ntt2_config),
  &X(ntt3_config),
  &X(ntt4_config),
  &X(ntt5_config),
  &X(ntt7_config),
  &X(ntt8_config),
  &X(ntt9_config),
  &X(ntt15_config),
  &X(ntt16_config),
  &X(ntt35_config),
  &X(ntt40_config),
};

static const nttconfig_t ** 
ntt_config_list(void)
{
  return ntt_config;
}

const nttgroup_t X(ntt_group) =
{
    SP_NAME_SUFFIX_STR,
    sizeof(ntt_config) / sizeof(ntt_config[0]),
    1,
    ntt_config_list
};

