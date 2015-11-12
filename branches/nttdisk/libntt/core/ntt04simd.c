#include "ntt/ntt-impl-simd.h"

#define NC 4

static const uint8_t ntt4_fixed_const[NC] = {1, 1, 1, 0};

extern void X(ntt4_init)(spv_t out, sp_t p, sp_t d, sp_t primroot, 
                        sp_t order, sp_t perm);

static void
ntt4_run_core_simd(spv_t in, spv_size_t istride, spv_size_t idist,
		spv_t out, spv_size_t ostride, spv_size_t odist,
		sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  sp_simd_t x0, x1, x2, x3;
  sp_simd_t t0, t1, t2, t3;
  sp_simd_t p0, p1, p2, p3;

  x0 = sp_simd_gather(in + 0 * istride, idist, vsize);
  x1 = sp_simd_gather(in + 1 * istride, idist, vsize);
  x2 = sp_simd_gather(in + 2 * istride, idist, vsize);
  x3 = sp_simd_gather(in + 3 * istride, idist, vsize);

  t0 = sp_ntt_add_simd(x0, x2, p);
  t2 = sp_ntt_sub_simd(x0, x2, p);
  t1 = sp_ntt_add_simd(x1, x3, p);
  t3 = sp_ntt_sub_partial_simd(x1, x3, p);

  t3 = sp_ntt_mul_simd(t3, ntt_const[3], ntt_const[NC+3], p);

  p0 = sp_ntt_add_simd(t0, t1, p);
  p1 = sp_ntt_sub_simd(t0, t1, p);
  p2 = sp_ntt_add_simd(t2, t3, p);
  p3 = sp_ntt_sub_simd(t2, t3, p);

  sp_simd_scatter(p0, out + 0 * ostride, odist, vsize);
  sp_simd_scatter(p2, out + 1 * ostride, odist, vsize);
  sp_simd_scatter(p1, out + 2 * ostride, odist, vsize);
  sp_simd_scatter(p3, out + 3 * ostride, odist, vsize);
}

static void
ntt4_twiddle_run_core_simd(
        spv_t in, spv_size_t istride, spv_size_t idist,
		spv_t out, spv_size_t ostride, spv_size_t odist,
		sp_simd_t *w, sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  sp_simd_t x0, x1, x2, x3;
  sp_simd_t t0, t1, t2, t3;
  sp_simd_t p0, p1, p2, p3;

  x0 = sp_simd_gather(in + 0 * istride, idist, vsize);
  x1 = sp_simd_gather(in + 1 * istride, idist, vsize);
  x2 = sp_simd_gather(in + 2 * istride, idist, vsize);
  x3 = sp_simd_gather(in + 3 * istride, idist, vsize);

  t0 = sp_ntt_add_simd(x0, x2, p);
  t2 = sp_ntt_sub_simd(x0, x2, p);
  t1 = sp_ntt_add_simd(x1, x3, p);
  t3 = sp_ntt_sub_partial_simd(x1, x3, p);

  t3 = sp_ntt_mul_simd(t3, ntt_const[3], ntt_const[NC+3], p);

  p0 = sp_ntt_add_simd(t0, t1, p);
  p1 = sp_ntt_sub_partial_simd(t0, t1, p);
  p2 = sp_ntt_add_partial_simd(t2, t3, p);
  p3 = sp_ntt_sub_partial_simd(t2, t3, p);

  p2 = sp_ntt_twiddle_mul_simd(p2, w + 0, p);
  p1 = sp_ntt_twiddle_mul_simd(p1, w + 2, p);
  p3 = sp_ntt_twiddle_mul_simd(p3, w + 4, p);

  sp_simd_scatter(p0, out + 0 * ostride, odist, vsize);
  sp_simd_scatter(p2, out + 1 * ostride, odist, vsize);
  sp_simd_scatter(p1, out + 2 * ostride, odist, vsize);
  sp_simd_scatter(p3, out + 3 * ostride, odist, vsize);
}

static void
ntt4_pfa_run_core_simd(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t inc2, spv_size_t n,
	  sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  spv_size_t j0, j1, j2, j3;
  sp_simd_t x0, x1, x2, x3;
  sp_simd_t t0, t1, t2, t3;
  sp_simd_t p0, p1, p2, p3;

  j0 = start;
  j1 = sp_array_inc(j0, inc, n);
  j2 = sp_array_inc(j0, 2 * inc, n);
  j3 = sp_array_inc(j0, 3 * inc, n);

  x0 = sp_simd_pfa_gather(x, j0, inc2, n, vsize);
  x1 = sp_simd_pfa_gather(x, j1, inc2, n, vsize);
  x2 = sp_simd_pfa_gather(x, j2, inc2, n, vsize);
  x3 = sp_simd_pfa_gather(x, j3, inc2, n, vsize);

  t0 = sp_ntt_add_simd(x0, x2, p);
  t2 = sp_ntt_sub_simd(x0, x2, p);
  t1 = sp_ntt_add_simd(x1, x3, p);
  t3 = sp_ntt_sub_partial_simd(x1, x3, p);

  t3 = sp_ntt_mul_simd(t3, ntt_const[3], ntt_const[NC+3], p);

  p0 = sp_ntt_add_simd(t0, t1, p);
  p1 = sp_ntt_sub_simd(t0, t1, p);
  p2 = sp_ntt_add_simd(t2, t3, p);
  p3 = sp_ntt_sub_simd(t2, t3, p);

  sp_simd_pfa_scatter(p0, x, j0, inc2, n, vsize);
  sp_simd_pfa_scatter(p2, x, j1, inc2, n, vsize);
  sp_simd_pfa_scatter(p1, x, j2, inc2, n, vsize);
  sp_simd_pfa_scatter(p3, x, j3, inc2, n, vsize);
}

DECLARE_CORE_ROUTINES_SIMD(4)
