#include "ntt-impl-simd.h"

/*-------------------------------------------------------------------------*/
static const nttconfig_t * ntt_config_simd[] = 
{
  &V(ntt2simd_config),
  &V(ntt3simd_config),
  &V(ntt4simd_config),
  &V(ntt5simd_config),
  &V(ntt7simd_config),
  &V(ntt8simd_config),
  &V(ntt9simd_config),
  &V(ntt15simd_config),
  &V(ntt16simd_config),
  &V(ntt35simd_config),
  &V(ntt40simd_config),
};

static const nttconfig_t ** 
ntt_config_list_simd(void)
{
  return ntt_config_simd;
}

/*-------------------------------------------------------------------------*/
const nttgroup_t V(ntt_group_simd) =
{
    SP_SIMD_NAME_SUFFIX_STR,
    sizeof(ntt_config_simd) / sizeof(ntt_config_simd[0]),
    SP_SIMD_VSIZE,
    ntt_config_list_simd
};

