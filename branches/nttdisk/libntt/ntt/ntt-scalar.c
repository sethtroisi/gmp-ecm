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

/*-------------------------------------------------------------------------*/
static spv_t
alloc_twiddle(sp_t primroot, sp_t order, sp_t p, sp_t d,
    		spv_size_t rows, spv_size_t cols)
{
  spv_size_t size = rows * cols;
  spv_t res = (spv_t)sp_aligned_malloc(2 * (rows - 1) * cols * sizeof(sp_t));

  sp_t w = sp_pow(primroot, order / size, p, d);
  sp_t w_inc = 1;
  spv_size_t i, j, k;

  for (i = 0; i < cols; i++)
    {
      sp_t w0 = w_inc;

      for (j = 0; j < rows - 1; j++)
	{
	  res[2*i*(rows-1) + 2*j] = w0;
	  res[2*i*(rows-1) + 2*j+1] = X(sp_ntt_reciprocal)(w0, p);
	  w0 = sp_mul(w0, w_inc, p, d);
	}

      w_inc = sp_mul(w_inc, w, p, d);
    }

  return res;
}

/*-------------------------------------------------------------------------*/
static void
free_twiddle(spv_t t)
{
    sp_aligned_free(t);
}

/*-------------------------------------------------------------------------*/
const nttgroup_t X(ntt_group) =
{
    SP_NAME_SUFFIX_STR,
    sizeof(ntt_config) / sizeof(ntt_config[0]),
    ntt_config_list,
    alloc_twiddle,
    free_twiddle
};

