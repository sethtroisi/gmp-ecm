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
static spv_t
alloc_twiddle_simd(sp_t primroot, sp_t order, sp_t p, sp_t d,
    		spv_size_t rows, spv_size_t cols)
{
  spv_size_t alloc = 2 * (rows - 1) * SP_SIMD_VSIZE *
      			((cols + SP_SIMD_VSIZE - 1) / SP_SIMD_VSIZE);
  spv_t res = (spv_t)sp_aligned_malloc(alloc * sizeof(sp_t));

  sp_t w = sp_pow(primroot, order / (rows * cols), p, d);
  sp_t w_inc = 1;
  spv_size_t i, j, k;
  spv_size_t num_simd = SP_SIMD_VSIZE * ((cols + SP_SIMD_VSIZE - 1) / 
      					SP_SIMD_VSIZE);

  for (i = 0; i < num_simd; i += SP_SIMD_VSIZE)
    {
      for (j = 0; j < SP_SIMD_VSIZE; j++)
	{
	  sp_t w0 = w_inc;

	  if (i + j < cols)
	    {
	      for (k = 0; k < rows - 1; k++)
		{
		  res[2*i*(rows-1) + (2*k)*SP_SIMD_VSIZE + j] = w0;
		  res[2*i*(rows-1) + (2*k+1)*SP_SIMD_VSIZE + j] = X(sp_ntt_reciprocal)(w0, p);
		  w0 = sp_mul(w0, w_inc, p, d);
		}

	      w_inc = sp_mul(w_inc, w, p, d);
	    }
	  else
	    {
	      for (k = 0; k < rows - 1; k++)
		{
		  res[2*i*(rows-1) + (2*k)*SP_SIMD_VSIZE + j] = 0;
		  res[2*i*(rows-1) + (2*k+1)*SP_SIMD_VSIZE + j] = 0;
		}
	    }
	}
    }

  return res;
}

/*-------------------------------------------------------------------------*/
static void
free_twiddle_simd(spv_t t)
{
    sp_aligned_free(t);
}

/*-------------------------------------------------------------------------*/
const nttgroup_t V(ntt_group_simd) =
{
    SP_SIMD_NAME_SUFFIX_STR,
    sizeof(ntt_config_simd) / sizeof(ntt_config_simd[0]),
    ntt_config_list_simd,
    alloc_twiddle_simd,
    free_twiddle_simd
};

