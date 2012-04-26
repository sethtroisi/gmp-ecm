#include "ntt-impl.h"

static uint32_t 
ntt4_get_num_const(void)
{
  return 4;
}

void
ntt4_init(spv_t out, sp_t p, sp_t d, 
	  sp_t primroot, sp_t order)
{
  out[0] = 1;
  out[1] = 1;
  out[2] = 1;
  out[3] = sp_pow(primroot, order / 4, p, d);
}

static void
ntt4_run(spv_t x, spv_size_t stride,
	  sp_t p, sp_t d, spv_t ntt_const)
{
  sp_t x0, x1, x2, x3;
  sp_t t0, t1, t2, t3;
  sp_t p0, p1, p2, p3;

  x0 = x[0 * stride];
  x1 = x[1 * stride];
  x2 = x[2 * stride];
  x3 = x[3 * stride];

  t0 = sp_add(x0, x2, p);
  t2 = sp_sub(x0, x2, p);
  t1 = sp_add(x1, x3, p);
  t3 = sp_sub(x1, x3, p);

  t3 = sp_mul(t3, ntt_const[3], p, d);

  p0 = sp_add(t0, t1, p);
  p1 = sp_sub(t0, t1, p);
  p2 = sp_add(t2, t3, p);
  p3 = sp_sub(t2, t3, p);

  x[0 * stride] = p0;
  x[1 * stride] = p2;
  x[2 * stride] = p1;
  x[3 * stride] = p3;
}

static void
ntt4_pfa_run_core(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t n,
	  sp_t p, sp_t d, spv_t ntt_const)
{
  spv_size_t j0, j1, j2, j3;

  sp_t x0, x1, x2, x3;
  sp_t t0, t1, t2, t3;
  sp_t p0, p1, p2, p3;

  j0 = start;
  j1 = sp_array_inc(j0, inc, n);
  j2 = sp_array_inc(j0, 2 * inc, n);
  j3 = sp_array_inc(j0, 3 * inc, n);

  x0 = x[j0];
  x1 = x[j1];
  x2 = x[j2];
  x3 = x[j3];

  t0 = sp_add(x0, x2, p);
  t2 = sp_sub(x0, x2, p);
  t1 = sp_add(x1, x3, p);
  t3 = sp_sub(x1, x3, p);

  t3 = sp_mul(t3, ntt_const[3], p, d);

  p0 = sp_add(t0, t1, p);
  p1 = sp_sub(t0, t1, p);
  p2 = sp_add(t2, t3, p);
  p3 = sp_sub(t2, t3, p);

  x[j0] = p0;
  x[j1] = p2;
  x[j2] = p1;
  x[j3] = p3;
}

#ifdef HAVE_SSE2
static void
ntt4_pfa_run_core_simd(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t inc2, spv_size_t n,
	  sp_t p, sp_t d, spv_t ntt_const)
{
  spv_size_t j0, j1, j2, j3;
  sp_simd_t x0, x1, x2, x3;
  sp_simd_t t0, t1, t2, t3;
  sp_simd_t p0, p1, p2, p3;

  j0 = start;
  j1 = sp_array_inc(j0, inc, n);
  j2 = sp_array_inc(j0, 2 * inc, n);
  j3 = sp_array_inc(j0, 3 * inc, n);

  x0 = sp_simd_gather(x, j0, inc2, n);
  x1 = sp_simd_gather(x, j1, inc2, n);
  x2 = sp_simd_gather(x, j2, inc2, n);
  x3 = sp_simd_gather(x, j3, inc2, n);

  t0 = sp_simd_add(x0, x2, p);
  t2 = sp_simd_sub(x0, x2, p);
  t1 = sp_simd_add(x1, x3, p);
  t3 = sp_simd_sub(x1, x3, p);

  t3 = sp_simd_mul(t3, ntt_const[3], p, d);

  p0 = sp_simd_add(t0, t1, p);
  p1 = sp_simd_sub(t0, t1, p);
  p2 = sp_simd_add(t2, t3, p);
  p3 = sp_simd_sub(t2, t3, p);

  sp_simd_scatter(p0, x, j0, inc2, n);
  sp_simd_scatter(p2, x, j1, inc2, n);
  sp_simd_scatter(p1, x, j2, inc2, n);
  sp_simd_scatter(p3, x, j3, inc2, n);
}
#endif

static void
ntt4_pfa_run(spv_t x, spv_size_t stride,
	  spv_size_t cofactor,
	  sp_t p, sp_t d, spv_t ntt_const)
{
  spv_size_t i = 0;
  spv_size_t incstart = 0;
  spv_size_t n = 4 * cofactor * stride;
  spv_size_t inc = cofactor * stride;
  spv_size_t inc2 = 4 * stride;

#ifdef HAVE_SSE2
  spv_size_t num_simd = SP_SIMD_VSIZE * (cofactor / SP_SIMD_VSIZE);

  for (i = 0; i < num_simd; i += SP_SIMD_VSIZE)
    {
      ntt4_pfa_run_core_simd(x, incstart, inc, inc2, n, p, d, ntt_const);
      incstart += SP_SIMD_VSIZE * inc2;
    }
#endif

  for (; i < cofactor; i++, incstart += inc2)
    ntt4_pfa_run_core(x, incstart, inc, n, p, d, ntt_const);

}

const nttconfig_t ntt4_config = 
{
  4,
  ntt4_get_num_const,
  ntt4_init,
  ntt4_run,
  ntt4_pfa_run
};

