#include "ntt/ntt-impl-simd.h"

#define NC 6

static const uint8_t ntt5_fixed_const[NC] = {1};

extern void X(ntt5_init)(spv_t out, sp_t p, sp_t d, sp_t primroot, 
                        sp_t order, sp_t perm);

static void
ntt5_run_core_simd(spv_t in, spv_size_t istride, spv_size_t idist,
		spv_t out, spv_size_t ostride, spv_size_t odist,
		sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  sp_simd_t p0, p1, p2, p3, p4, p5;
  sp_simd_t x0, x1, x2, x3, x4;
  sp_simd_t     t1, t2, t3, t4;

  x0 = sp_simd_gather(in + 0 * istride, idist, vsize);
  x1 = sp_simd_gather(in + 1 * istride, idist, vsize);
  x4 = sp_simd_gather(in + 2 * istride, idist, vsize);
  x2 = sp_simd_gather(in + 3 * istride, idist, vsize);
  x3 = sp_simd_gather(in + 4 * istride, idist, vsize);

  t1 = sp_ntt_add_simd(x1, x3, p);
  t3 = sp_ntt_sub_simd(x1, x3, p);
  t2 = sp_ntt_add_simd(x2, x4, p);
  t4 = sp_ntt_sub_simd(x2, x4, p);

  p1 = sp_ntt_add_simd(t1, t2, p);
  p2 = sp_ntt_sub_partial_simd(t1, t2, p);
  p3 = t3;
  p4 = t4;
  p5 = sp_ntt_add_partial_simd(t3, t4, p);

  p0 = sp_ntt_add_simd(x0, p1, p);

  p1 = sp_ntt_mul_simd(p1, ntt_const[1], ntt_const[NC+1], p);
  p2 = sp_ntt_mul_simd(p2, ntt_const[2], ntt_const[NC+2], p);
  p3 = sp_ntt_mul_simd(p3, ntt_const[3], ntt_const[NC+3], p);
  p4 = sp_ntt_mul_simd(p4, ntt_const[4], ntt_const[NC+4], p);
  p5 = sp_ntt_mul_simd(p5, ntt_const[5], ntt_const[NC+5], p);

  p1 = sp_ntt_add_simd(p0, p1, p);

  t1 = sp_ntt_add_simd(p1, p2, p);
  t2 = sp_ntt_sub_simd(p1, p2, p);
  t3 = sp_ntt_add_simd(p3, p5, p);
  t4 = sp_ntt_add_simd(p4, p5, p);

  p1 = sp_ntt_add_simd(t1, t3, p);
  p2 = sp_ntt_add_simd(t2, t4, p);
  p3 = sp_ntt_sub_simd(t1, t3, p);
  p4 = sp_ntt_sub_simd(t2, t4, p);

  sp_simd_scatter(p0, out + 0 * ostride, odist, vsize);
  sp_simd_scatter(p4, out + 1 * ostride, odist, vsize);
  sp_simd_scatter(p3, out + 2 * ostride, odist, vsize);
  sp_simd_scatter(p1, out + 3 * ostride, odist, vsize);
  sp_simd_scatter(p2, out + 4 * ostride, odist, vsize);
}

static void
ntt5_run_core_simd_interleaved(
    		spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		sp_simd_t p, sp_simd_t * c)
{
  sp_simd_t p0, p1, p2, p3, p4, p5;
  sp_simd_t x0, x1, x2, x3, x4;
  sp_simd_t     t1, t2, t3, t4;

  x0 = sp_simd_load(in + 0 * istride);
  x1 = sp_simd_load(in + 1 * istride);
  x4 = sp_simd_load(in + 2 * istride);
  x2 = sp_simd_load(in + 3 * istride);
  x3 = sp_simd_load(in + 4 * istride);

  t1 = sp_ntt_add_simd0(x1, x3, p);
  t3 = sp_ntt_sub_simd0(x1, x3, p);
  t2 = sp_ntt_add_simd0(x2, x4, p);
  t4 = sp_ntt_sub_simd0(x2, x4, p);

  p1 = sp_ntt_add_simd0(t1, t2, p);
  p2 = sp_ntt_sub_partial_simd0(t1, t2, p);
  p3 = t3;
  p4 = t4;
  p5 = sp_ntt_add_partial_simd0(t3, t4, p);

  p0 = sp_ntt_add_simd0(x0, p1, p);

  p1 = sp_ntt_mul_simd0(p1, c+2, p);
  p2 = sp_ntt_mul_simd0(p2, c+4, p);
  p3 = sp_ntt_mul_simd0(p3, c+6, p);
  p4 = sp_ntt_mul_simd0(p4, c+8, p);
  p5 = sp_ntt_mul_simd0(p5, c+10, p);

  p1 = sp_ntt_add_simd0(p0, p1, p);

  t1 = sp_ntt_add_simd0(p1, p2, p);
  t2 = sp_ntt_sub_simd0(p1, p2, p);
  t3 = sp_ntt_add_simd0(p3, p5, p);
  t4 = sp_ntt_add_simd0(p4, p5, p);

  p1 = sp_ntt_add_simd0(t1, t3, p);
  p2 = sp_ntt_add_simd0(t2, t4, p);
  p3 = sp_ntt_sub_simd0(t1, t3, p);
  p4 = sp_ntt_sub_simd0(t2, t4, p);

  sp_simd_store(p0, out + 0 * ostride);
  sp_simd_store(p4, out + 1 * ostride);
  sp_simd_store(p3, out + 2 * ostride);
  sp_simd_store(p1, out + 3 * ostride);
  sp_simd_store(p2, out + 4 * ostride);
}

static void
ntt5_twiddle_run_core_simd(
        spv_t in, spv_size_t istride, spv_size_t idist,
		spv_t out, spv_size_t ostride, spv_size_t odist,
		sp_simd_t *w, sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  sp_simd_t p0, p1, p2, p3, p4, p5;
  sp_simd_t x0, x1, x2, x3, x4;
  sp_simd_t     t1, t2, t3, t4;

  x0 = sp_simd_gather(in + 0 * istride, idist, vsize);
  x1 = sp_simd_gather(in + 1 * istride, idist, vsize);
  x4 = sp_simd_gather(in + 2 * istride, idist, vsize);
  x2 = sp_simd_gather(in + 3 * istride, idist, vsize);
  x3 = sp_simd_gather(in + 4 * istride, idist, vsize);

  t1 = sp_ntt_add_simd(x1, x3, p);
  t3 = sp_ntt_sub_simd(x1, x3, p);
  t2 = sp_ntt_add_simd(x2, x4, p);
  t4 = sp_ntt_sub_simd(x2, x4, p);

  p1 = sp_ntt_add_simd(t1, t2, p);
  p2 = sp_ntt_sub_partial_simd(t1, t2, p);
  p3 = t3;
  p4 = t4;
  p5 = sp_ntt_add_partial_simd(t3, t4, p);

  p0 = sp_ntt_add_simd(x0, p1, p);

  p1 = sp_ntt_mul_simd(p1, ntt_const[1], ntt_const[NC+1], p);
  p2 = sp_ntt_mul_simd(p2, ntt_const[2], ntt_const[NC+2], p);
  p3 = sp_ntt_mul_simd(p3, ntt_const[3], ntt_const[NC+3], p);
  p4 = sp_ntt_mul_simd(p4, ntt_const[4], ntt_const[NC+4], p);
  p5 = sp_ntt_mul_simd(p5, ntt_const[5], ntt_const[NC+5], p);

  p1 = sp_ntt_add_simd(p0, p1, p);

  t1 = sp_ntt_add_simd(p1, p2, p);
  t2 = sp_ntt_sub_simd(p1, p2, p);
  t3 = sp_ntt_add_simd(p3, p5, p);
  t4 = sp_ntt_add_simd(p4, p5, p);

  p1 = sp_ntt_add_partial_simd(t1, t3, p);
  p2 = sp_ntt_add_partial_simd(t2, t4, p);
  p3 = sp_ntt_sub_partial_simd(t1, t3, p);
  p4 = sp_ntt_sub_partial_simd(t2, t4, p);

  p4 = sp_ntt_twiddle_mul_simd(p4, w + 0, p);
  p3 = sp_ntt_twiddle_mul_simd(p3, w + 2, p);
  p1 = sp_ntt_twiddle_mul_simd(p1, w + 4, p);
  p2 = sp_ntt_twiddle_mul_simd(p2, w + 6, p);

  sp_simd_scatter(p0, out + 0 * ostride, odist, vsize);
  sp_simd_scatter(p4, out + 1 * ostride, odist, vsize);
  sp_simd_scatter(p3, out + 2 * ostride, odist, vsize);
  sp_simd_scatter(p1, out + 3 * ostride, odist, vsize);
  sp_simd_scatter(p2, out + 4 * ostride, odist, vsize);
}

static void
ntt5_twiddle_run_core_simd_interleaved(
		spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		sp_simd_t *w, sp_simd_t p, sp_simd_t *c)
{
  sp_simd_t p0, p1, p2, p3, p4, p5;
  sp_simd_t x0, x1, x2, x3, x4;
  sp_simd_t     t1, t2, t3, t4;

  x0 = sp_simd_load(in + 0 * istride);
  x1 = sp_simd_load(in + 1 * istride);
  x4 = sp_simd_load(in + 2 * istride);
  x2 = sp_simd_load(in + 3 * istride);
  x3 = sp_simd_load(in + 4 * istride);

  t1 = sp_ntt_add_simd0(x1, x3, p);
  t3 = sp_ntt_sub_simd0(x1, x3, p);
  t2 = sp_ntt_add_simd0(x2, x4, p);
  t4 = sp_ntt_sub_simd0(x2, x4, p);

  p1 = sp_ntt_add_simd0(t1, t2, p);
  p2 = sp_ntt_sub_partial_simd0(t1, t2, p);
  p3 = t3;
  p4 = t4;
  p5 = sp_ntt_add_partial_simd0(t3, t4, p);

  p0 = sp_ntt_add_simd0(x0, p1, p);

  p1 = sp_ntt_mul_simd0(p1, c+2, p);
  p2 = sp_ntt_mul_simd0(p2, c+4, p);
  p3 = sp_ntt_mul_simd0(p3, c+6, p);
  p4 = sp_ntt_mul_simd0(p4, c+8, p);
  p5 = sp_ntt_mul_simd0(p5, c+10, p);

  p1 = sp_ntt_add_simd0(p0, p1, p);

  t1 = sp_ntt_add_simd0(p1, p2, p);
  t2 = sp_ntt_sub_simd0(p1, p2, p);
  t3 = sp_ntt_add_simd0(p3, p5, p);
  t4 = sp_ntt_add_simd0(p4, p5, p);

  p1 = sp_ntt_add_partial_simd0(t1, t3, p);
  p2 = sp_ntt_add_partial_simd0(t2, t4, p);
  p3 = sp_ntt_sub_partial_simd0(t1, t3, p);
  p4 = sp_ntt_sub_partial_simd0(t2, t4, p);

  p4 = sp_ntt_twiddle_mul_simd0(p4, w+0, p);
  p3 = sp_ntt_twiddle_mul_simd0(p3, w+2, p);
  p1 = sp_ntt_twiddle_mul_simd0(p1, w+4, p);
  p2 = sp_ntt_twiddle_mul_simd0(p2, w+6, p);

  sp_simd_store(p0, out + 0 * ostride);
  sp_simd_store(p4, out + 1 * ostride);
  sp_simd_store(p3, out + 2 * ostride);
  sp_simd_store(p1, out + 3 * ostride);
  sp_simd_store(p2, out + 4 * ostride);
}

static void
ntt5_pfa_run_core_simd(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t inc2, spv_size_t n,
	  sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  spv_size_t j0, j1, j2, j3, j4;
  sp_simd_t p0, p1, p2, p3, p4, p5;
  sp_simd_t x0, x1, x2, x3, x4;
  sp_simd_t     t1, t2, t3, t4;

  j0 = start;
  j1 = sp_array_inc(j0, inc, n);
  j2 = sp_array_inc(j0, 2 * inc, n);
  j3 = sp_array_inc(j0, 3 * inc, n);
  j4 = sp_array_inc(j0, 4 * inc, n);

  x0 = sp_simd_pfa_gather(x, j0, inc2, n, vsize);
  x1 = sp_simd_pfa_gather(x, j1, inc2, n, vsize);
  x4 = sp_simd_pfa_gather(x, j2, inc2, n, vsize);
  x2 = sp_simd_pfa_gather(x, j3, inc2, n, vsize);
  x3 = sp_simd_pfa_gather(x, j4, inc2, n, vsize);

  t1 = sp_ntt_add_simd(x1, x3, p);
  t3 = sp_ntt_sub_simd(x1, x3, p);
  t2 = sp_ntt_add_simd(x2, x4, p);
  t4 = sp_ntt_sub_simd(x2, x4, p);

  p1 = sp_ntt_add_simd(t1, t2, p);
  p2 = sp_ntt_sub_partial_simd(t1, t2, p);
  p3 = t3;
  p4 = t4;
  p5 = sp_ntt_add_partial_simd(t3, t4, p);

  p0 = sp_ntt_add_simd(x0, p1, p);

  p1 = sp_ntt_mul_simd(p1, ntt_const[1], ntt_const[NC+1], p);
  p2 = sp_ntt_mul_simd(p2, ntt_const[2], ntt_const[NC+2], p);
  p3 = sp_ntt_mul_simd(p3, ntt_const[3], ntt_const[NC+3], p);
  p4 = sp_ntt_mul_simd(p4, ntt_const[4], ntt_const[NC+4], p);
  p5 = sp_ntt_mul_simd(p5, ntt_const[5], ntt_const[NC+5], p);

  p1 = sp_ntt_add_simd(p0, p1, p);

  t1 = sp_ntt_add_simd(p1, p2, p);
  t2 = sp_ntt_sub_simd(p1, p2, p);
  t3 = sp_ntt_add_simd(p3, p5, p);
  t4 = sp_ntt_add_simd(p4, p5, p);

  p1 = sp_ntt_add_simd(t1, t3, p);
  p2 = sp_ntt_add_simd(t2, t4, p);
  p3 = sp_ntt_sub_simd(t1, t3, p);
  p4 = sp_ntt_sub_simd(t2, t4, p);

  sp_simd_pfa_scatter(p0, x, j0, inc2, n, vsize);
  sp_simd_pfa_scatter(p4, x, j1, inc2, n, vsize);
  sp_simd_pfa_scatter(p3, x, j2, inc2, n, vsize);
  sp_simd_pfa_scatter(p1, x, j3, inc2, n, vsize);
  sp_simd_pfa_scatter(p2, x, j4, inc2, n, vsize);
}

static void
ntt5_pfa_run_core_simd_interleaved(
	  spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t n,
	  sp_simd_t p, sp_simd_t *c)
{
  spv_size_t j0, j1, j2, j3, j4;
  sp_simd_t p0, p1, p2, p3, p4, p5;
  sp_simd_t x0, x1, x2, x3, x4;
  sp_simd_t     t1, t2, t3, t4;

  j0 = start;
  j1 = sp_array_inc(j0, inc, n);
  j2 = sp_array_inc(j0, 2 * inc, n);
  j3 = sp_array_inc(j0, 3 * inc, n);
  j4 = sp_array_inc(j0, 4 * inc, n);

  x0 = sp_simd_load(x + j0);
  x1 = sp_simd_load(x + j1);
  x4 = sp_simd_load(x + j2);
  x2 = sp_simd_load(x + j3);
  x3 = sp_simd_load(x + j4);

  t1 = sp_ntt_add_simd0(x1, x3, p);
  t3 = sp_ntt_sub_simd0(x1, x3, p);
  t2 = sp_ntt_add_simd0(x2, x4, p);
  t4 = sp_ntt_sub_simd0(x2, x4, p);

  p1 = sp_ntt_add_simd0(t1, t2, p);
  p2 = sp_ntt_sub_partial_simd0(t1, t2, p);
  p3 = t3;
  p4 = t4;
  p5 = sp_ntt_add_partial_simd0(t3, t4, p);

  p0 = sp_ntt_add_simd0(x0, p1, p);

  p1 = sp_ntt_mul_simd0(p1, c+2, p);
  p2 = sp_ntt_mul_simd0(p2, c+4, p);
  p3 = sp_ntt_mul_simd0(p3, c+6, p);
  p4 = sp_ntt_mul_simd0(p4, c+8, p);
  p5 = sp_ntt_mul_simd0(p5, c+10, p);

  p1 = sp_ntt_add_simd0(p0, p1, p);

  t1 = sp_ntt_add_simd0(p1, p2, p);
  t2 = sp_ntt_sub_simd0(p1, p2, p);
  t3 = sp_ntt_add_simd0(p3, p5, p);
  t4 = sp_ntt_add_simd0(p4, p5, p);

  p1 = sp_ntt_add_simd0(t1, t3, p);
  p2 = sp_ntt_add_simd0(t2, t4, p);
  p3 = sp_ntt_sub_simd0(t1, t3, p);
  p4 = sp_ntt_sub_simd0(t2, t4, p);

  sp_simd_store(p0, x + j0);
  sp_simd_store(p4, x + j1);
  sp_simd_store(p3, x + j2);
  sp_simd_store(p1, x + j3);
  sp_simd_store(p2, x + j4);
}

DECLARE_CORE_ROUTINES_SIMD(5)
