#include "ntt/ntt-impl-simd.h"

#define NC 8

static const uint8_t ntt8_fixed_const[NC] = {1, 1, 1, 0, 1, 0, 0, 0};

extern void X(ntt8_init)(spv_t out, sp_t p, sp_t d, sp_t primroot, 
                        sp_t order, sp_t perm);


static void
ntt8_run_core_simd(spv_t in, spv_size_t istride, spv_size_t idist,
		spv_t out, spv_size_t ostride, spv_size_t odist,
		sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  sp_simd_t x0, x1, x2, x3, x4, x5, x6, x7;
  sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7; 
  sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7;

  x0 = sp_simd_gather(in + 0 * istride, idist, vsize);
  x1 = sp_simd_gather(in + 1 * istride, idist, vsize);
  x2 = sp_simd_gather(in + 2 * istride, idist, vsize);
  x3 = sp_simd_gather(in + 3 * istride, idist, vsize);
  x4 = sp_simd_gather(in + 4 * istride, idist, vsize);
  x5 = sp_simd_gather(in + 5 * istride, idist, vsize);
  x6 = sp_simd_gather(in + 6 * istride, idist, vsize);
  x7 = sp_simd_gather(in + 7 * istride, idist, vsize);

  t0 = sp_ntt_add_simd(x0, x4, p);
  t4 = sp_ntt_sub_simd(x0, x4, p);
  t1 = sp_ntt_add_simd(x1, x5, p);
  t5 = sp_ntt_sub_simd(x1, x5, p);
  t2 = sp_ntt_add_simd(x2, x6, p);
  t6 = sp_ntt_sub_partial_simd(x2, x6, p);
  t3 = sp_ntt_add_simd(x3, x7, p);
  t7 = sp_ntt_sub_simd(x3, x7, p);

  p0 = sp_ntt_add_simd(t0, t2, p);
  p1 = sp_ntt_sub_simd(t0, t2, p);
  p2 = sp_ntt_add_simd(t1, t3, p);
  p3 = sp_ntt_sub_partial_simd(t1, t3, p);
  p4 = t4;
  p5 = t6;
  p6 = sp_ntt_sub_partial_simd(t5, t7, p);
  p7 = sp_ntt_add_partial_simd(t5, t7, p); 

  p3 = sp_ntt_mul_simd(p3, ntt_const[3], ntt_const[NC+3], p);
  p5 = sp_ntt_mul_simd(p5, ntt_const[5], ntt_const[NC+5], p);
  p6 = sp_ntt_mul_simd(p6, ntt_const[6], ntt_const[NC+6], p);
  p7 = sp_ntt_mul_simd(p7, ntt_const[7], ntt_const[NC+7], p);

  t0 = sp_ntt_add_simd(p4, p5, p);
  t1 = sp_ntt_sub_simd(p4, p5, p);
  t2 = sp_ntt_add_simd(p6, p7, p);
  t3 = sp_ntt_sub_simd(p6, p7, p);
  t4 = sp_ntt_add_simd(t0, t2, p);
  t5 = sp_ntt_sub_simd(t0, t2, p);
  t6 = sp_ntt_add_simd(t1, t3, p);
  t7 = sp_ntt_sub_simd(t1, t3, p);

  t0 = sp_ntt_add_simd(p0, p2, p);
  t1 = sp_ntt_sub_simd(p0, p2, p);
  t2 = sp_ntt_add_simd(p1, p3, p);
  t3 = sp_ntt_sub_simd(p1, p3, p);

  sp_simd_scatter(t0, out + 0 * ostride, odist, vsize);
  sp_simd_scatter(t4, out + 1 * ostride, odist, vsize);
  sp_simd_scatter(t2, out + 2 * ostride, odist, vsize);
  sp_simd_scatter(t7, out + 3 * ostride, odist, vsize);
  sp_simd_scatter(t1, out + 4 * ostride, odist, vsize);
  sp_simd_scatter(t5, out + 5 * ostride, odist, vsize);
  sp_simd_scatter(t3, out + 6 * ostride, odist, vsize);
  sp_simd_scatter(t6, out + 7 * ostride, odist, vsize);
}

static void
ntt8_run_core_simd_interleaved(
		spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		sp_simd_t p, sp_simd_t * c)
{
  sp_simd_t x0, x1, x2, x3, x4, x5, x6, x7;
  sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7; 
  sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7;

  x0 = sp_simd_load(in + 0 * istride);
  x1 = sp_simd_load(in + 1 * istride);
  x2 = sp_simd_load(in + 2 * istride);
  x3 = sp_simd_load(in + 3 * istride);
  x4 = sp_simd_load(in + 4 * istride);
  x5 = sp_simd_load(in + 5 * istride);
  x6 = sp_simd_load(in + 6 * istride);
  x7 = sp_simd_load(in + 7 * istride);

  t0 = sp_ntt_add_simd0(x0, x4, p);
  t4 = sp_ntt_sub_simd0(x0, x4, p);
  t1 = sp_ntt_add_simd0(x1, x5, p);
  t5 = sp_ntt_sub_simd0(x1, x5, p);
  t2 = sp_ntt_add_simd0(x2, x6, p);
  t6 = sp_ntt_sub_partial_simd0(x2, x6, p);
  t3 = sp_ntt_add_simd0(x3, x7, p);
  t7 = sp_ntt_sub_simd0(x3, x7, p);

  p0 = sp_ntt_add_simd0(t0, t2, p);
  p1 = sp_ntt_sub_simd0(t0, t2, p);
  p2 = sp_ntt_add_simd0(t1, t3, p);
  p3 = sp_ntt_sub_partial_simd0(t1, t3, p);
  p4 = t4;
  p5 = t6;
  p6 = sp_ntt_sub_partial_simd0(t5, t7, p);
  p7 = sp_ntt_add_partial_simd0(t5, t7, p); 

  p3 = sp_ntt_mul_simd0(p3, c+6, p);
  p5 = sp_ntt_mul_simd0(p5, c+10, p);
  p6 = sp_ntt_mul_simd0(p6, c+12, p);
  p7 = sp_ntt_mul_simd0(p7, c+14, p);

  t0 = sp_ntt_add_simd0(p4, p5, p);
  t1 = sp_ntt_sub_simd0(p4, p5, p);
  t2 = sp_ntt_add_simd0(p6, p7, p);
  t3 = sp_ntt_sub_simd0(p6, p7, p);
  t4 = sp_ntt_add_simd0(t0, t2, p);
  t5 = sp_ntt_sub_simd0(t0, t2, p);
  t6 = sp_ntt_add_simd0(t1, t3, p);
  t7 = sp_ntt_sub_simd0(t1, t3, p);

  t0 = sp_ntt_add_simd0(p0, p2, p);
  t1 = sp_ntt_sub_simd0(p0, p2, p);
  t2 = sp_ntt_add_simd0(p1, p3, p);
  t3 = sp_ntt_sub_simd0(p1, p3, p);

  sp_simd_store(t0, out + 0 * ostride);
  sp_simd_store(t4, out + 1 * ostride);
  sp_simd_store(t2, out + 2 * ostride);
  sp_simd_store(t7, out + 3 * ostride);
  sp_simd_store(t1, out + 4 * ostride);
  sp_simd_store(t5, out + 5 * ostride);
  sp_simd_store(t3, out + 6 * ostride);
  sp_simd_store(t6, out + 7 * ostride);
}

static void
ntt8_twiddle_run_core_simd(
        spv_t in, spv_size_t istride, spv_size_t idist,
		spv_t out, spv_size_t ostride, spv_size_t odist,
		sp_simd_t *w, sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  sp_simd_t x0, x1, x2, x3, x4, x5, x6, x7;
  sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7; 
  sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7;

  x0 = sp_simd_gather(in + 0 * istride, idist, vsize);
  x1 = sp_simd_gather(in + 1 * istride, idist, vsize);
  x2 = sp_simd_gather(in + 2 * istride, idist, vsize);
  x3 = sp_simd_gather(in + 3 * istride, idist, vsize);
  x4 = sp_simd_gather(in + 4 * istride, idist, vsize);
  x5 = sp_simd_gather(in + 5 * istride, idist, vsize);
  x6 = sp_simd_gather(in + 6 * istride, idist, vsize);
  x7 = sp_simd_gather(in + 7 * istride, idist, vsize);

  t0 = sp_ntt_add_simd(x0, x4, p);
  t4 = sp_ntt_sub_simd(x0, x4, p);
  t1 = sp_ntt_add_simd(x1, x5, p);
  t5 = sp_ntt_sub_simd(x1, x5, p);
  t2 = sp_ntt_add_simd(x2, x6, p);
  t6 = sp_ntt_sub_partial_simd(x2, x6, p);
  t3 = sp_ntt_add_simd(x3, x7, p);
  t7 = sp_ntt_sub_simd(x3, x7, p);

  p0 = sp_ntt_add_simd(t0, t2, p);
  p1 = sp_ntt_sub_simd(t0, t2, p);
  p2 = sp_ntt_add_simd(t1, t3, p);
  p3 = sp_ntt_sub_partial_simd(t1, t3, p);
  p4 = t4;
  p5 = t6;
  p6 = sp_ntt_sub_partial_simd(t5, t7, p);
  p7 = sp_ntt_add_partial_simd(t5, t7, p); 

  p3 = sp_ntt_mul_simd(p3, ntt_const[3], ntt_const[NC+3], p);
  p5 = sp_ntt_mul_simd(p5, ntt_const[5], ntt_const[NC+5], p);
  p6 = sp_ntt_mul_simd(p6, ntt_const[6], ntt_const[NC+6], p);
  p7 = sp_ntt_mul_simd(p7, ntt_const[7], ntt_const[NC+7], p);

  t0 = sp_ntt_add_simd(p4, p5, p);
  t1 = sp_ntt_sub_simd(p4, p5, p);
  t2 = sp_ntt_add_simd(p6, p7, p);
  t3 = sp_ntt_sub_simd(p6, p7, p);
  t4 = sp_ntt_add_partial_simd(t0, t2, p);
  t5 = sp_ntt_sub_partial_simd(t0, t2, p);
  t6 = sp_ntt_add_partial_simd(t1, t3, p);
  t7 = sp_ntt_sub_partial_simd(t1, t3, p);

  t0 = sp_ntt_add_simd(p0, p2, p);
  t1 = sp_ntt_sub_partial_simd(p0, p2, p);
  t2 = sp_ntt_add_partial_simd(p1, p3, p);
  t3 = sp_ntt_sub_partial_simd(p1, p3, p);

  t4 = sp_ntt_twiddle_mul_simd(t4, w + 0, p);
  t2 = sp_ntt_twiddle_mul_simd(t2, w + 2, p);
  t7 = sp_ntt_twiddle_mul_simd(t7, w + 4, p);
  t1 = sp_ntt_twiddle_mul_simd(t1, w + 6, p);
  t5 = sp_ntt_twiddle_mul_simd(t5, w + 8, p);
  t3 = sp_ntt_twiddle_mul_simd(t3, w + 10, p);
  t6 = sp_ntt_twiddle_mul_simd(t6, w + 12, p);

  sp_simd_scatter(t0, out + 0 * ostride, odist, vsize);
  sp_simd_scatter(t4, out + 1 * ostride, odist, vsize);
  sp_simd_scatter(t2, out + 2 * ostride, odist, vsize);
  sp_simd_scatter(t7, out + 3 * ostride, odist, vsize);
  sp_simd_scatter(t1, out + 4 * ostride, odist, vsize);
  sp_simd_scatter(t5, out + 5 * ostride, odist, vsize);
  sp_simd_scatter(t3, out + 6 * ostride, odist, vsize);
  sp_simd_scatter(t6, out + 7 * ostride, odist, vsize);
}

static void
ntt8_twiddle_run_core_simd_interleaved(
		spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		sp_simd_t *w, sp_simd_t p, sp_simd_t * c)
{
  sp_simd_t x0, x1, x2, x3, x4, x5, x6, x7;
  sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7; 
  sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7;

  x0 = sp_simd_load(in + 0 * istride);
  x1 = sp_simd_load(in + 1 * istride);
  x2 = sp_simd_load(in + 2 * istride);
  x3 = sp_simd_load(in + 3 * istride);
  x4 = sp_simd_load(in + 4 * istride);
  x5 = sp_simd_load(in + 5 * istride);
  x6 = sp_simd_load(in + 6 * istride);
  x7 = sp_simd_load(in + 7 * istride);

  t0 = sp_ntt_add_simd0(x0, x4, p);
  t4 = sp_ntt_sub_simd0(x0, x4, p);
  t1 = sp_ntt_add_simd0(x1, x5, p);
  t5 = sp_ntt_sub_simd0(x1, x5, p);
  t2 = sp_ntt_add_simd0(x2, x6, p);
  t6 = sp_ntt_sub_partial_simd0(x2, x6, p);
  t3 = sp_ntt_add_simd0(x3, x7, p);
  t7 = sp_ntt_sub_simd0(x3, x7, p);

  p0 = sp_ntt_add_simd0(t0, t2, p);
  p1 = sp_ntt_sub_simd0(t0, t2, p);
  p2 = sp_ntt_add_simd0(t1, t3, p);
  p3 = sp_ntt_sub_partial_simd0(t1, t3, p);
  p4 = t4;
  p5 = t6;
  p6 = sp_ntt_sub_partial_simd0(t5, t7, p);
  p7 = sp_ntt_add_partial_simd0(t5, t7, p); 

  p3 = sp_ntt_mul_simd0(p3, c+6, p);
  p5 = sp_ntt_mul_simd0(p5, c+10, p);
  p6 = sp_ntt_mul_simd0(p6, c+12, p);
  p7 = sp_ntt_mul_simd0(p7, c+14, p);

  t0 = sp_ntt_add_simd0(p4, p5, p);
  t1 = sp_ntt_sub_simd0(p4, p5, p);
  t2 = sp_ntt_add_simd0(p6, p7, p);
  t3 = sp_ntt_sub_simd0(p6, p7, p);
  t4 = sp_ntt_add_partial_simd0(t0, t2, p);
  t5 = sp_ntt_sub_partial_simd0(t0, t2, p);
  t6 = sp_ntt_add_partial_simd0(t1, t3, p);
  t7 = sp_ntt_sub_partial_simd0(t1, t3, p);

  t0 = sp_ntt_add_simd0(p0, p2, p);
  t1 = sp_ntt_sub_partial_simd0(p0, p2, p);
  t2 = sp_ntt_add_partial_simd0(p1, p3, p);
  t3 = sp_ntt_sub_partial_simd0(p1, p3, p);

  t4 = sp_ntt_twiddle_mul_simd0(t4, w+0, p);
  t2 = sp_ntt_twiddle_mul_simd0(t2, w+2, p);
  t7 = sp_ntt_twiddle_mul_simd0(t7, w+4, p);
  t1 = sp_ntt_twiddle_mul_simd0(t1, w+6, p);
  t5 = sp_ntt_twiddle_mul_simd0(t5, w+8, p);
  t3 = sp_ntt_twiddle_mul_simd0(t3, w+10, p);
  t6 = sp_ntt_twiddle_mul_simd0(t6, w+12, p);

  sp_simd_store(t0, out + 0 * ostride);
  sp_simd_store(t4, out + 1 * ostride);
  sp_simd_store(t2, out + 2 * ostride);
  sp_simd_store(t7, out + 3 * ostride);
  sp_simd_store(t1, out + 4 * ostride);
  sp_simd_store(t5, out + 5 * ostride);
  sp_simd_store(t3, out + 6 * ostride);
  sp_simd_store(t6, out + 7 * ostride);
}

static void
ntt8_pfa_run_core_simd(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t inc2, spv_size_t n,
	  sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  spv_size_t j0, j1, j2, j3, j4, j5, j6, j7;
  sp_simd_t x0, x1, x2, x3, x4, x5, x6, x7;
  sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7; 
  sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7;

  j0 = start;
  j1 = sp_array_inc(j0, inc, n);
  j2 = sp_array_inc(j0, 2 * inc, n);
  j3 = sp_array_inc(j0, 3 * inc, n);
  j4 = sp_array_inc(j0, 4 * inc, n);
  j5 = sp_array_inc(j0, 5 * inc, n);
  j6 = sp_array_inc(j0, 6 * inc, n);
  j7 = sp_array_inc(j0, 7 * inc, n);

  x0 = sp_simd_pfa_gather(x, j0, inc2, n, vsize);
  x1 = sp_simd_pfa_gather(x, j1, inc2, n, vsize);
  x2 = sp_simd_pfa_gather(x, j2, inc2, n, vsize);
  x3 = sp_simd_pfa_gather(x, j3, inc2, n, vsize);
  x4 = sp_simd_pfa_gather(x, j4, inc2, n, vsize);
  x5 = sp_simd_pfa_gather(x, j5, inc2, n, vsize);
  x6 = sp_simd_pfa_gather(x, j6, inc2, n, vsize);
  x7 = sp_simd_pfa_gather(x, j7, inc2, n, vsize);

  t0 = sp_ntt_add_simd(x0, x4, p);
  t4 = sp_ntt_sub_simd(x0, x4, p);
  t1 = sp_ntt_add_simd(x1, x5, p);
  t5 = sp_ntt_sub_simd(x1, x5, p);
  t2 = sp_ntt_add_simd(x2, x6, p);
  t6 = sp_ntt_sub_partial_simd(x2, x6, p);
  t3 = sp_ntt_add_simd(x3, x7, p);
  t7 = sp_ntt_sub_simd(x3, x7, p);

  p0 = sp_ntt_add_simd(t0, t2, p);
  p1 = sp_ntt_sub_simd(t0, t2, p);
  p2 = sp_ntt_add_simd(t1, t3, p);
  p3 = sp_ntt_sub_partial_simd(t1, t3, p);
  p4 = t4;
  p5 = t6;
  p6 = sp_ntt_sub_partial_simd(t5, t7, p);
  p7 = sp_ntt_add_partial_simd(t5, t7, p); 

  p3 = sp_ntt_mul_simd(p3, ntt_const[3], ntt_const[NC+3], p);
  p5 = sp_ntt_mul_simd(p5, ntt_const[5], ntt_const[NC+5], p);
  p6 = sp_ntt_mul_simd(p6, ntt_const[6], ntt_const[NC+6], p);
  p7 = sp_ntt_mul_simd(p7, ntt_const[7], ntt_const[NC+7], p);

  t0 = sp_ntt_add_simd(p4, p5, p);
  t1 = sp_ntt_sub_simd(p4, p5, p);
  t2 = sp_ntt_add_simd(p6, p7, p);
  t3 = sp_ntt_sub_simd(p6, p7, p);
  t4 = sp_ntt_add_simd(t0, t2, p);
  t5 = sp_ntt_sub_simd(t0, t2, p);
  t6 = sp_ntt_add_simd(t1, t3, p);
  t7 = sp_ntt_sub_simd(t1, t3, p);

  t0 = sp_ntt_add_simd(p0, p2, p);
  t1 = sp_ntt_sub_simd(p0, p2, p);
  t2 = sp_ntt_add_simd(p1, p3, p);
  t3 = sp_ntt_sub_simd(p1, p3, p);

  sp_simd_pfa_scatter(t0, x, j0, inc2, n, vsize);
  sp_simd_pfa_scatter(t4, x, j1, inc2, n, vsize);
  sp_simd_pfa_scatter(t2, x, j2, inc2, n, vsize);
  sp_simd_pfa_scatter(t7, x, j3, inc2, n, vsize);
  sp_simd_pfa_scatter(t1, x, j4, inc2, n, vsize);
  sp_simd_pfa_scatter(t5, x, j5, inc2, n, vsize);
  sp_simd_pfa_scatter(t3, x, j6, inc2, n, vsize);
  sp_simd_pfa_scatter(t6, x, j7, inc2, n, vsize);
}

static void
ntt8_pfa_run_core_simd_interleaved(
	  spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t n,
	  sp_simd_t p, sp_simd_t * c)
{
  spv_size_t j0, j1, j2, j3, j4, j5, j6, j7;
  sp_simd_t x0, x1, x2, x3, x4, x5, x6, x7;
  sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7; 
  sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7;

  j0 = start;
  j1 = sp_array_inc(j0, inc, n);
  j2 = sp_array_inc(j0, 2 * inc, n);
  j3 = sp_array_inc(j0, 3 * inc, n);
  j4 = sp_array_inc(j0, 4 * inc, n);
  j5 = sp_array_inc(j0, 5 * inc, n);
  j6 = sp_array_inc(j0, 6 * inc, n);
  j7 = sp_array_inc(j0, 7 * inc, n);

  x0 = sp_simd_load(x + j0);
  x1 = sp_simd_load(x + j1);
  x2 = sp_simd_load(x + j2);
  x3 = sp_simd_load(x + j3);
  x4 = sp_simd_load(x + j4);
  x5 = sp_simd_load(x + j5);
  x6 = sp_simd_load(x + j6);
  x7 = sp_simd_load(x + j7);

  t0 = sp_ntt_add_simd0(x0, x4, p);
  t4 = sp_ntt_sub_simd0(x0, x4, p);
  t1 = sp_ntt_add_simd0(x1, x5, p);
  t5 = sp_ntt_sub_simd0(x1, x5, p);
  t2 = sp_ntt_add_simd0(x2, x6, p);
  t6 = sp_ntt_sub_partial_simd0(x2, x6, p);
  t3 = sp_ntt_add_simd0(x3, x7, p);
  t7 = sp_ntt_sub_simd0(x3, x7, p);

  p0 = sp_ntt_add_simd0(t0, t2, p);
  p1 = sp_ntt_sub_simd0(t0, t2, p);
  p2 = sp_ntt_add_simd0(t1, t3, p);
  p3 = sp_ntt_sub_partial_simd0(t1, t3, p);
  p4 = t4;
  p5 = t6;
  p6 = sp_ntt_sub_partial_simd0(t5, t7, p);
  p7 = sp_ntt_add_partial_simd0(t5, t7, p); 

  p3 = sp_ntt_mul_simd0(p3, c+6, p);
  p5 = sp_ntt_mul_simd0(p5, c+10, p);
  p6 = sp_ntt_mul_simd0(p6, c+12, p);
  p7 = sp_ntt_mul_simd0(p7, c+14, p);

  t0 = sp_ntt_add_simd0(p4, p5, p);
  t1 = sp_ntt_sub_simd0(p4, p5, p);
  t2 = sp_ntt_add_simd0(p6, p7, p);
  t3 = sp_ntt_sub_simd0(p6, p7, p);
  t4 = sp_ntt_add_simd0(t0, t2, p);
  t5 = sp_ntt_sub_simd0(t0, t2, p);
  t6 = sp_ntt_add_simd0(t1, t3, p);
  t7 = sp_ntt_sub_simd0(t1, t3, p);

  t0 = sp_ntt_add_simd0(p0, p2, p);
  t1 = sp_ntt_sub_simd0(p0, p2, p);
  t2 = sp_ntt_add_simd0(p1, p3, p);
  t3 = sp_ntt_sub_simd0(p1, p3, p);

  sp_simd_store(t0, x + j0);
  sp_simd_store(t4, x + j1);
  sp_simd_store(t2, x + j2);
  sp_simd_store(t7, x + j3);
  sp_simd_store(t1, x + j4);
  sp_simd_store(t5, x + j5);
  sp_simd_store(t3, x + j6);
  sp_simd_store(t6, x + j7);
}

DECLARE_CORE_ROUTINES_SIMD(8)
