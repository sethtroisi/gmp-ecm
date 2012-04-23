#include "ntt-impl.h"

static uint32_t 
ntt2_get_num_const(void)
{
  return 2;
}

static void
ntt2_init(spv_t out, sp_t p, sp_t d, 
	  sp_t primroot, sp_t order)
{
  out[0] = 1;
  out[1] = 1;
}

static void
ntt2_run(spv_t x, spv_size_t stride,
	  sp_t p, sp_t d, spv_t ntt_const)
{
  sp_t x0, x1;

  x0 = x[0 * stride];
  x1 = x[1 * stride];

  x[0 * stride] = sp_add(x0, x1, p);
  x[1 * stride] = sp_sub(x0, x1, p);
}

static void
ntt2_pfa_run(spv_t x, spv_size_t stride,
	  spv_size_t cofactor,
	  sp_t p, sp_t d, spv_t ntt_const)
{
  spv_size_t i, jstart;
  spv_size_t n = 2 * cofactor * stride;
  spv_size_t inc = cofactor * stride;
  spv_size_t inc2 = 2 * stride;

  for (i = jstart = 0; i < cofactor; i++, jstart += inc2)
    {
      spv_size_t j0, j1;

      sp_t x0, x1;
      sp_t t0, t1;

      j0 = jstart;
      j1 = sp_array_inc(j0, inc, n);

      x0 = x[j0];
      x1 = x[j1];

      t0 = sp_add(x0, x1, p);
      t1 = sp_sub(x0, x1, p);

      x[j0] = t0;
      x[j1] = t1;
    }
}

const nttconfig_t ntt2_config = 
{
  2,
  ntt2_get_num_const,
  ntt2_init,
  ntt2_run,
  ntt2_pfa_run
};

