#include "ntt/ntt-impl-scalar.h"

#define NC 2

static const uint8_t ntt2_fixed_const[NC] = {1, 1};

void
X(ntt2_init)(spv_t out, sp_t p, sp_t d, 
	  sp_t primroot, sp_t order, sp_t perm)
{
  out[0] = out[1] = 1;
}

static void 
ntt2_run_core(spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		sp_t p, spv_t ntt_const)
{
  sp_t p0, p1;
  sp_t x0, x1;

  x0 = in[0 * istride];
  x1 = in[1 * istride];

  p0 = sp_ntt_add(x0, x1, p);
  p1 = sp_ntt_sub(x0, x1, p);

  out[0 * ostride] = p0;
  out[1 * ostride] = p1;
}


static void
ntt2_twiddle_run_core(spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		spv_t w, sp_t p, spv_t ntt_const)
{
  sp_t p0, p1;
  sp_t x0, x1;

  x0 = in[0 * istride];
  x1 = in[1 * istride];

  p0 = sp_ntt_add(x0, x1, p);
  p1 = sp_ntt_sub_partial(x0, x1, p);

  p1 = sp_ntt_mul(p1, w[0], w[1], p);

  out[0 * ostride] = p0;
  out[1 * ostride] = p1;
}


static void
ntt2_pfa_run_core(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t n,
	  sp_t p, spv_t ntt_const)
{
  spv_size_t j0, j1;
  sp_t p0, p1;
  sp_t x0, x1;

  j0 = start;
  j1 = sp_array_inc(j0, inc, n);

  x0 = x[j0];
  x1 = x[j1];

  p0 = sp_ntt_add(x0, x1, p);
  p1 = sp_ntt_sub(x0, x1, p);

  x[j0] = p0;
  x[j1] = p1;
}

DECLARE_CORE_ROUTINES(2)

