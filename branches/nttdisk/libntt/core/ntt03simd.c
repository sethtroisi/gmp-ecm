#include "ntt/ntt-impl-simd.h"

#define NC 3

static const uint8_t ntt3_fixed_const[NC] = {1};

extern void X(ntt3_init)(spv_t out, sp_t p, sp_t d, sp_t primroot, 
                        sp_t order, sp_t perm);


static void
ntt3_run_core_simd(spv_t in, spv_size_t istride, spv_size_t idist,
		spv_t out, spv_size_t ostride, spv_size_t odist,
		sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  sp_simd_t p0, p1, p2;
  sp_simd_t x0, x1, x2;
  sp_simd_t     t1, t2;

  x0 = sp_simd_gather(in + 0 * istride, idist, vsize);
  x1 = sp_simd_gather(in + 1 * istride, idist, vsize);
  x2 = sp_simd_gather(in + 2 * istride, idist, vsize);

  t1 = sp_ntt_add_simd(x1, x2, p);
  t2 = sp_ntt_sub_partial_simd(x1, x2, p);

  p0 = sp_ntt_add_simd(x0, t1, p);

  p1 = sp_ntt_mul_simd(t1, ntt_const[1], ntt_const[NC+1], p);
  p2 = sp_ntt_mul_simd(t2, ntt_const[2], ntt_const[NC+2], p);

  p1 = sp_ntt_add_simd(p0, p1, p);

  t1 = sp_ntt_add_simd(p1, p2, p);
  t2 = sp_ntt_sub_simd(p1, p2, p);

  sp_simd_scatter(p0, out + 0 * ostride, odist, vsize);
  sp_simd_scatter(t1, out + 1 * ostride, odist, vsize);
  sp_simd_scatter(t2, out + 2 * ostride, odist, vsize);
}


static void
ntt3_run_core_simd_interleaved(
    		spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		sp_simd_t p, sp_simd_t * c)
{
  sp_simd_t p0, p1, p2;
  sp_simd_t x0, x1, x2;
  sp_simd_t     t1, t2;

  x0 = sp_simd_load(in + 0 * istride);
  x1 = sp_simd_load(in + 1 * istride);
  x2 = sp_simd_load(in + 2 * istride);

  t1 = sp_ntt_add_simd0(x1, x2, p);
  t2 = sp_ntt_sub_partial_simd0(x1, x2, p);

  p0 = sp_ntt_add_simd0(x0, t1, p);

  p1 = sp_ntt_mul_simd0(t1, c+2, p);
  p2 = sp_ntt_mul_simd0(t2, c+4, p);

  p1 = sp_ntt_add_simd0(p0, p1, p);

  t1 = sp_ntt_add_simd0(p1, p2, p);
  t2 = sp_ntt_sub_simd0(p1, p2, p);

  sp_simd_store(p0, out + 0 * ostride);
  sp_simd_store(t1, out + 1 * ostride);
  sp_simd_store(t2, out + 2 * ostride);
}

static void
ntt3_twiddle_run_core_simd(
        spv_t in, spv_size_t istride, spv_size_t idist,
		spv_t out, spv_size_t ostride, spv_size_t odist,
		sp_simd_t *w, sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  sp_simd_t p0, p1, p2;
  sp_simd_t x0, x1, x2;
  sp_simd_t     t1, t2;

  x0 = sp_simd_gather(in + 0 * istride, idist, vsize);
  x1 = sp_simd_gather(in + 1 * istride, idist, vsize);
  x2 = sp_simd_gather(in + 2 * istride, idist, vsize);

  t1 = sp_ntt_add_simd(x1, x2, p);
  t2 = sp_ntt_sub_partial_simd(x1, x2, p);

  p0 = sp_ntt_add_simd(x0, t1, p);

  p1 = sp_ntt_mul_simd(t1, ntt_const[1], ntt_const[NC+1], p);
  p2 = sp_ntt_mul_simd(t2, ntt_const[2], ntt_const[NC+2], p);

  p1 = sp_ntt_add_simd(p0, p1, p);

  t1 = sp_ntt_add_partial_simd(p1, p2, p);
  t2 = sp_ntt_sub_partial_simd(p1, p2, p);

  t1 = sp_ntt_twiddle_mul_simd(t1, w + 0, p);
  t2 = sp_ntt_twiddle_mul_simd(t2, w + 2, p);

  sp_simd_scatter(p0, out + 0 * ostride, odist, vsize);
  sp_simd_scatter(t1, out + 1 * ostride, odist, vsize);
  sp_simd_scatter(t2, out + 2 * ostride, odist, vsize);
}

static void
ntt3_twiddle_run_core_simd_interleaved(
		spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		sp_simd_t *w, sp_simd_t p, sp_simd_t * c)
{
  sp_simd_t p0, p1, p2;
  sp_simd_t x0, x1, x2;
  sp_simd_t     t1, t2;

  x0 = sp_simd_load(in + 0 * istride);
  x1 = sp_simd_load(in + 1 * istride);
  x2 = sp_simd_load(in + 2 * istride);

  t1 = sp_ntt_add_simd0(x1, x2, p);
  t2 = sp_ntt_sub_partial_simd0(x1, x2, p);

  p0 = sp_ntt_add_simd0(x0, t1, p);

  p1 = sp_ntt_mul_simd0(t1, c+2, p);
  p2 = sp_ntt_mul_simd0(t2, c+4, p);

  p1 = sp_ntt_add_simd0(p0, p1, p);

  t1 = sp_ntt_add_partial_simd0(p1, p2, p);
  t2 = sp_ntt_sub_partial_simd0(p1, p2, p);

  t1 = sp_ntt_twiddle_mul_simd0(t1, w+0, p);
  t2 = sp_ntt_twiddle_mul_simd0(t2, w+2, p);

  sp_simd_store(p0, out + 0 * ostride);
  sp_simd_store(t1, out + 1 * ostride);
  sp_simd_store(t2, out + 2 * ostride);
}

static void
ntt3_pfa_run_core_simd(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t inc2, spv_size_t n,
	  sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  spv_size_t j0, j1, j2;
  sp_simd_t p0, p1, p2;
  sp_simd_t x0, x1, x2;
  sp_simd_t     t1, t2;

  j0 = start;
  j1 = sp_array_inc(j0, inc, n);
  j2 = sp_array_inc(j0, 2 * inc, n);

  x0 = sp_simd_pfa_gather(x, j0, inc2, n, vsize);
  x1 = sp_simd_pfa_gather(x, j1, inc2, n, vsize);
  x2 = sp_simd_pfa_gather(x, j2, inc2, n, vsize);

  t1 = sp_ntt_add_simd(x1, x2, p);
  t2 = sp_ntt_sub_partial_simd(x1, x2, p);

  p0 = sp_ntt_add_simd(x0, t1, p);

  p1 = sp_ntt_mul_simd(t1, ntt_const[1], ntt_const[NC+1], p);
  p2 = sp_ntt_mul_simd(t2, ntt_const[2], ntt_const[NC+2], p);

  p1 = sp_ntt_add_simd(p0, p1, p);

  t1 = sp_ntt_add_simd(p1, p2, p);
  t2 = sp_ntt_sub_simd(p1, p2, p);

  sp_simd_pfa_scatter(p0, x, j0, inc2, n, vsize);
  sp_simd_pfa_scatter(t1, x, j1, inc2, n, vsize);
  sp_simd_pfa_scatter(t2, x, j2, inc2, n, vsize);
}

static void
ntt3_pfa_run_core_simd_interleaved(
	  spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t n,
	  sp_simd_t p, sp_simd_t * c)
{
  spv_size_t j0, j1, j2;
  sp_simd_t p0, p1, p2;
  sp_simd_t x0, x1, x2;
  sp_simd_t     t1, t2;

  j0 = start;
  j1 = sp_array_inc(j0, inc, n);
  j2 = sp_array_inc(j0, 2 * inc, n);

  x0 = sp_simd_load(x + j0);
  x1 = sp_simd_load(x + j1);
  x2 = sp_simd_load(x + j2);

  t1 = sp_ntt_add_simd0(x1, x2, p);
  t2 = sp_ntt_sub_partial_simd0(x1, x2, p);

  p0 = sp_ntt_add_simd0(x0, t1, p);

  p1 = sp_ntt_mul_simd0(t1, c+2, p);
  p2 = sp_ntt_mul_simd0(t2, c+4, p);

  p1 = sp_ntt_add_simd0(p0, p1, p);

  t1 = sp_ntt_add_simd0(p1, p2, p);
  t2 = sp_ntt_sub_simd0(p1, p2, p);

  sp_simd_store(p0, x + j0);
  sp_simd_store(t1, x + j1);
  sp_simd_store(t2, x + j2);
}

DECLARE_CORE_ROUTINES_SIMD(3)
