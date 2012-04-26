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
ntt2_pfa_run_core(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t n,
	  sp_t p, sp_t d, spv_t ntt_const)
{
  spv_size_t j0, j1;
  sp_t x0, x1;
  sp_t t0, t1;

  j0 = start;
  j1 = sp_array_inc(j0, inc, n);

  x0 = x[j0];
  x1 = x[j1];

  t0 = sp_add(x0, x1, p);
  t1 = sp_sub(x0, x1, p);

  x[j0] = t0;
  x[j1] = t1;
}

#ifdef HAVE_SSE2
static void
ntt2_pfa_run_core_simd(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t inc2, spv_size_t n,
	  sp_t p, sp_t d, spv_t ntt_const)
{
  spv_size_t j0, j1;
  sp_simd_t x0, x1;
  sp_simd_t t0, t1;

  j0 = start;
  j1 = sp_array_inc(j0, inc, n);

  x0 = sp_simd_gather(x, j0, inc2, n);
  x1 = sp_simd_gather(x, j1, inc2, n);

  t0 = sp_simd_add(x0, x1, p);
  t1 = sp_simd_sub(x0, x1, p);

  sp_simd_scatter(t0, x, j0, inc2, n);
  sp_simd_scatter(t1, x, j1, inc2, n);
}
#endif

static void
ntt2_pfa_run(spv_t x, spv_size_t stride,
	  spv_size_t cofactor,
	  sp_t p, sp_t d, spv_t ntt_const)
{
  spv_size_t i = 0;
  spv_size_t incstart = 0;
  spv_size_t n = 2 * cofactor * stride;
  spv_size_t inc = cofactor * stride;
  spv_size_t inc2 = 2 * stride;

#ifdef HAVE_SSE2
  spv_size_t num_simd = SP_SIMD_VSIZE * (cofactor / SP_SIMD_VSIZE);

  for (i = 0; i < num_simd; i += SP_SIMD_VSIZE)
    {
      ntt2_pfa_run_core_simd(x, incstart, inc, inc2, n, p, d, ntt_const);
      incstart += SP_SIMD_VSIZE * inc2;
    }
#endif

  for (; i < cofactor; i++, incstart += inc2)
    ntt2_pfa_run_core(x, incstart, inc, n, p, d, ntt_const);

}

const nttconfig_t ntt2_config = 
{
  2,
  ntt2_get_num_const,
  ntt2_init,
  ntt2_run,
  ntt2_pfa_run
};

