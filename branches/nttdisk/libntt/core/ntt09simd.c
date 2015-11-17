#include "ntt/ntt-impl-simd.h"

#define NC 11

static const uint8_t ntt9_fixed_const[NC] = {1};

extern void X(ntt9_init)(spv_t out, sp_t p, sp_t d, sp_t primroot, 
                        sp_t order, sp_t perm);


static void
ntt9_run_core_simd(spv_t in, spv_size_t istride, spv_size_t idist,
		spv_t out, spv_size_t ostride, spv_size_t odist,
		sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  sp_simd_t x0, x1, x2, x3, x4, x5, x6;
  sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
  sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

  sp_simd_t x0e, x1e, t0e, t1e, p0e, p1e, p2e;

  x0 = sp_simd_gather(in + 0 * istride, idist, vsize);
  x1 = sp_simd_gather(in + 1 * istride, idist, vsize);
  x2 = sp_simd_gather(in + 2 * istride, idist, vsize);
  x0e = sp_simd_gather(in + 3 * istride, idist, vsize);
  x3 = sp_simd_gather(in + 4 * istride, idist, vsize);
  x6 = sp_simd_gather(in + 5 * istride, idist, vsize);
  x1e = sp_simd_gather(in + 6 * istride, idist, vsize);
  x5 = sp_simd_gather(in + 7 * istride, idist, vsize);
  x4 = sp_simd_gather(in + 8 * istride, idist, vsize);

  t0e = sp_ntt_add_simd(x0e, x1e, p);
  t1e = sp_ntt_sub_partial_simd(x0e, x1e, p);

  p1 = sp_ntt_add_simd(x1, x3, p);
  p1 = sp_ntt_add_simd(p1, x5, p);
  p2 = sp_ntt_add_simd(x2, x4, p);
  p2 = sp_ntt_add_simd(p2, x6, p);
  p3 = sp_ntt_sub_simd(x1, x5, p);
  p4 = sp_ntt_sub_simd(x2, x6, p);
  p5 = sp_ntt_sub_simd(x3, x5, p);
  p6 = sp_ntt_sub_simd(x4, x6, p);

  t1 = sp_ntt_add_simd(p1, p2, p);
  t2 = sp_ntt_sub_partial_simd(p1, p2, p);
  t3 = sp_ntt_sub_simd(p3, p5, p);
  t5 = sp_ntt_add_simd(t3, p6, p);
  t3 = sp_ntt_sub_simd(t3, p6, p);
  t4 = sp_ntt_add_simd(p4, p5, p);
  t6 = sp_ntt_sub_simd(p4, p5, p);

  p0e = sp_ntt_add_simd(x0, t0e, p);
  p0 = t1;
  p1 = t1;
  p2 = t2;
  p3 = t3;
  p4 = t4;
  p5 = sp_ntt_add_partial_simd(t3, t4, p);
  p6 = t5;
  p7 = t6;
  p8 = sp_ntt_add_partial_simd(t5, t6, p);

  p1 = sp_ntt_mul_simd(p1, ntt_const[1], ntt_const[NC+1], p);
  p2 = sp_ntt_mul_simd(p2, ntt_const[2], ntt_const[NC+2], p);
  t0e = sp_ntt_mul_simd(t0e, ntt_const[3], ntt_const[NC+3], p);
  t1e = sp_ntt_mul_simd(t1e, ntt_const[4], ntt_const[NC+4], p);
  p3 = sp_ntt_mul_simd(p3, ntt_const[5], ntt_const[NC+5], p);
  p4 = sp_ntt_mul_simd(p4, ntt_const[6], ntt_const[NC+6], p);
  p5 = sp_ntt_mul_simd(p5, ntt_const[7], ntt_const[NC+7], p);
  p6 = sp_ntt_mul_simd(p6, ntt_const[8], ntt_const[NC+8], p);
  p7 = sp_ntt_mul_simd(p7, ntt_const[9], ntt_const[NC+9], p);
  p8 = sp_ntt_mul_simd(p8, ntt_const[10], ntt_const[NC+10], p);

  t0e = sp_ntt_add_simd(t0e, p0e, p);
  t1 = sp_ntt_add_simd(p1, p2, p);
  t2 = sp_ntt_sub_simd(p1, p2, p);
  t3 = sp_ntt_add_simd(p3, p5, p);
  t4 = sp_ntt_add_simd(p4, p5, p);
  t5 = sp_ntt_add_simd(p6, p8, p);
  t6 = sp_ntt_add_simd(p7, p8, p);

  p1e = sp_ntt_add_simd(t0e, t1e, p);
  p2e = sp_ntt_sub_simd(t0e, t1e, p);
  p3 = sp_ntt_add_simd(t3, t5, p);
  p4 = sp_ntt_add_simd(t4, t6, p);
  p5 = sp_ntt_sub_simd(t4, t6, p);
  p5 = sp_ntt_sub_simd(p5, p3, p);
  p6 = sp_ntt_sub_simd(t5, t3, p);

  p0 = sp_ntt_add_simd(p0, p0e, p);
  t1 = sp_ntt_add_simd(t1, p0e, p);
  t2 = sp_ntt_add_simd(t2, p0e, p);
  t3 = sp_ntt_add_simd(p3, p1e, p);
  t4 = sp_ntt_add_simd(p4, p2e, p);
  t5 = sp_ntt_add_simd(p5, p1e, p);
  t6 = sp_ntt_add_simd(p6, p2e, p);
  t7 = sp_ntt_add_simd(p3, p5, p);
  t7 = sp_ntt_sub_simd(p1e, t7, p);
  t8 = sp_ntt_add_simd(p4, p6, p);
  t8 = sp_ntt_sub_simd(p2e, t8, p);

  sp_simd_scatter(p0, out + 0 * ostride, odist, vsize);
  sp_simd_scatter(t8, out + 1 * ostride, odist, vsize);
  sp_simd_scatter(t3, out + 2 * ostride, odist, vsize);
  sp_simd_scatter(t2, out + 3 * ostride, odist, vsize);
  sp_simd_scatter(t4, out + 4 * ostride, odist, vsize);
  sp_simd_scatter(t7, out + 5 * ostride, odist, vsize);
  sp_simd_scatter(t1, out + 6 * ostride, odist, vsize);
  sp_simd_scatter(t6, out + 7 * ostride, odist, vsize);
  sp_simd_scatter(t5, out + 8 * ostride, odist, vsize);
}

static void
ntt9_run_core_simd_interleaved(
    		spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		sp_simd_t p, sp_simd_t * c)
{
  sp_simd_t x0, x1, x2, x3, x4, x5, x6;
  sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
  sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

  sp_simd_t x0e, x1e, t0e, t1e, p0e, p1e, p2e;

  x0 = sp_simd_load(in + 0 * istride);
  x1 = sp_simd_load(in + 1 * istride);
  x2 = sp_simd_load(in + 2 * istride);
  x0e = sp_simd_load(in + 3 * istride);
  x3 = sp_simd_load(in + 4 * istride);
  x6 = sp_simd_load(in + 5 * istride);
  x1e = sp_simd_load(in + 6 * istride);
  x5 = sp_simd_load(in + 7 * istride);
  x4 = sp_simd_load(in + 8 * istride);

  t0e = sp_ntt_add_simd0(x0e, x1e, p);
  t1e = sp_ntt_sub_partial_simd0(x0e, x1e, p);

  p1 = sp_ntt_add_simd0(x1, x3, p);
  p1 = sp_ntt_add_simd0(p1, x5, p);
  p2 = sp_ntt_add_simd0(x2, x4, p);
  p2 = sp_ntt_add_simd0(p2, x6, p);
  p3 = sp_ntt_sub_simd0(x1, x5, p);
  p4 = sp_ntt_sub_simd0(x2, x6, p);
  p5 = sp_ntt_sub_simd0(x3, x5, p);
  p6 = sp_ntt_sub_simd0(x4, x6, p);

  t1 = sp_ntt_add_simd0(p1, p2, p);
  t2 = sp_ntt_sub_partial_simd0(p1, p2, p);
  t3 = sp_ntt_sub_simd0(p3, p5, p);
  t5 = sp_ntt_add_simd0(t3, p6, p);
  t3 = sp_ntt_sub_simd0(t3, p6, p);
  t4 = sp_ntt_add_simd0(p4, p5, p);
  t6 = sp_ntt_sub_simd0(p4, p5, p);

  p0e = sp_ntt_add_simd0(x0, t0e, p);
  p0 = t1;
  p1 = t1;
  p2 = t2;
  p3 = t3;
  p4 = t4;
  p5 = sp_ntt_add_partial_simd0(t3, t4, p);
  p6 = t5;
  p7 = t6;
  p8 = sp_ntt_add_partial_simd0(t5, t6, p);

  p1 = sp_ntt_mul_simd0(p1, c+2, p);
  p2 = sp_ntt_mul_simd0(p2, c+4, p);
  t0e = sp_ntt_mul_simd0(t0e, c+6, p);
  t1e = sp_ntt_mul_simd0(t1e, c+8, p);
  p3 = sp_ntt_mul_simd0(p3, c+10, p);
  p4 = sp_ntt_mul_simd0(p4, c+12, p);
  p5 = sp_ntt_mul_simd0(p5, c+14, p);
  p6 = sp_ntt_mul_simd0(p6, c+16, p);
  p7 = sp_ntt_mul_simd0(p7, c+18, p);
  p8 = sp_ntt_mul_simd0(p8, c+20, p);

  t0e = sp_ntt_add_simd0(t0e, p0e, p);
  t1 = sp_ntt_add_simd0(p1, p2, p);
  t2 = sp_ntt_sub_simd0(p1, p2, p);
  t3 = sp_ntt_add_simd0(p3, p5, p);
  t4 = sp_ntt_add_simd0(p4, p5, p);
  t5 = sp_ntt_add_simd0(p6, p8, p);
  t6 = sp_ntt_add_simd0(p7, p8, p);

  p1e = sp_ntt_add_simd0(t0e, t1e, p);
  p2e = sp_ntt_sub_simd0(t0e, t1e, p);
  p3 = sp_ntt_add_simd0(t3, t5, p);
  p4 = sp_ntt_add_simd0(t4, t6, p);
  p5 = sp_ntt_sub_simd0(t4, t6, p);
  p5 = sp_ntt_sub_simd0(p5, p3, p);
  p6 = sp_ntt_sub_simd0(t5, t3, p);

  p0 = sp_ntt_add_simd0(p0, p0e, p);
  t1 = sp_ntt_add_simd0(t1, p0e, p);
  t2 = sp_ntt_add_simd0(t2, p0e, p);
  t3 = sp_ntt_add_simd0(p3, p1e, p);
  t4 = sp_ntt_add_simd0(p4, p2e, p);
  t5 = sp_ntt_add_simd0(p5, p1e, p);
  t6 = sp_ntt_add_simd0(p6, p2e, p);
  t7 = sp_ntt_add_simd0(p3, p5, p);
  t7 = sp_ntt_sub_simd0(p1e, t7, p);
  t8 = sp_ntt_add_simd0(p4, p6, p);
  t8 = sp_ntt_sub_simd0(p2e, t8, p);

  sp_simd_store(p0, out + 0 * ostride);
  sp_simd_store(t8, out + 1 * ostride);
  sp_simd_store(t3, out + 2 * ostride);
  sp_simd_store(t2, out + 3 * ostride);
  sp_simd_store(t4, out + 4 * ostride);
  sp_simd_store(t7, out + 5 * ostride);
  sp_simd_store(t1, out + 6 * ostride);
  sp_simd_store(t6, out + 7 * ostride);
  sp_simd_store(t5, out + 8 * ostride);
}

static void
ntt9_twiddle_run_core_simd(
        spv_t in, spv_size_t istride, spv_size_t idist,
		spv_t out, spv_size_t ostride, spv_size_t odist,
		sp_simd_t *w, sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  sp_simd_t x0, x1, x2, x3, x4, x5, x6;
  sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
  sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

  sp_simd_t x0e, x1e, t0e, t1e, p0e, p1e, p2e;

  x0 = sp_simd_gather(in + 0 * istride, idist, vsize);
  x1 = sp_simd_gather(in + 1 * istride, idist, vsize);
  x2 = sp_simd_gather(in + 2 * istride, idist, vsize);
  x0e = sp_simd_gather(in + 3 * istride, idist, vsize);
  x3 = sp_simd_gather(in + 4 * istride, idist, vsize);
  x6 = sp_simd_gather(in + 5 * istride, idist, vsize);
  x1e = sp_simd_gather(in + 6 * istride, idist, vsize);
  x5 = sp_simd_gather(in + 7 * istride, idist, vsize);
  x4 = sp_simd_gather(in + 8 * istride, idist, vsize);

  t0e = sp_ntt_add_simd(x0e, x1e, p);
  t1e = sp_ntt_sub_partial_simd(x0e, x1e, p);

  p1 = sp_ntt_add_simd(x1, x3, p);
  p1 = sp_ntt_add_simd(p1, x5, p);
  p2 = sp_ntt_add_simd(x2, x4, p);
  p2 = sp_ntt_add_simd(p2, x6, p);
  p3 = sp_ntt_sub_simd(x1, x5, p);
  p4 = sp_ntt_sub_simd(x2, x6, p);
  p5 = sp_ntt_sub_simd(x3, x5, p);
  p6 = sp_ntt_sub_simd(x4, x6, p);

  t1 = sp_ntt_add_simd(p1, p2, p);
  t2 = sp_ntt_sub_partial_simd(p1, p2, p);
  t3 = sp_ntt_sub_simd(p3, p5, p);
  t5 = sp_ntt_add_simd(t3, p6, p);
  t3 = sp_ntt_sub_simd(t3, p6, p);
  t4 = sp_ntt_add_simd(p4, p5, p);
  t6 = sp_ntt_sub_simd(p4, p5, p);

  p0e = sp_ntt_add_simd(x0, t0e, p);
  p0 = t1;
  p1 = t1;
  p2 = t2;
  p3 = t3;
  p4 = t4;
  p5 = sp_ntt_add_partial_simd(t3, t4, p);
  p6 = t5;
  p7 = t6;
  p8 = sp_ntt_add_partial_simd(t5, t6, p);

  p1 = sp_ntt_mul_simd(p1, ntt_const[1], ntt_const[NC+1], p);
  p2 = sp_ntt_mul_simd(p2, ntt_const[2], ntt_const[NC+2], p);
  t0e = sp_ntt_mul_simd(t0e, ntt_const[3], ntt_const[NC+3], p);
  t1e = sp_ntt_mul_simd(t1e, ntt_const[4], ntt_const[NC+4], p);
  p3 = sp_ntt_mul_simd(p3, ntt_const[5], ntt_const[NC+5], p);
  p4 = sp_ntt_mul_simd(p4, ntt_const[6], ntt_const[NC+6], p);
  p5 = sp_ntt_mul_simd(p5, ntt_const[7], ntt_const[NC+7], p);
  p6 = sp_ntt_mul_simd(p6, ntt_const[8], ntt_const[NC+8], p);
  p7 = sp_ntt_mul_simd(p7, ntt_const[9], ntt_const[NC+9], p);
  p8 = sp_ntt_mul_simd(p8, ntt_const[10], ntt_const[NC+10], p);

  t0e = sp_ntt_add_simd(t0e, p0e, p);
  t1 = sp_ntt_add_simd(p1, p2, p);
  t2 = sp_ntt_sub_simd(p1, p2, p);
  t3 = sp_ntt_add_simd(p3, p5, p);
  t4 = sp_ntt_add_simd(p4, p5, p);
  t5 = sp_ntt_add_simd(p6, p8, p);
  t6 = sp_ntt_add_simd(p7, p8, p);

  p1e = sp_ntt_add_simd(t0e, t1e, p);
  p2e = sp_ntt_sub_simd(t0e, t1e, p);
  p3 = sp_ntt_add_simd(t3, t5, p);
  p4 = sp_ntt_add_simd(t4, t6, p);
  p5 = sp_ntt_sub_simd(t4, t6, p);
  p5 = sp_ntt_sub_simd(p5, p3, p);
  p6 = sp_ntt_sub_simd(t5, t3, p);

  p0 = sp_ntt_add_simd(p0, p0e, p);
  t1 = sp_ntt_add_partial_simd(t1, p0e, p);
  t2 = sp_ntt_add_partial_simd(t2, p0e, p);
  t3 = sp_ntt_add_partial_simd(p3, p1e, p);
  t4 = sp_ntt_add_partial_simd(p4, p2e, p);
  t5 = sp_ntt_add_partial_simd(p5, p1e, p);
  t6 = sp_ntt_add_partial_simd(p6, p2e, p);
  t7 = sp_ntt_add_simd(p3, p5, p);
  t7 = sp_ntt_sub_partial_simd(p1e, t7, p);
  t8 = sp_ntt_add_simd(p4, p6, p);
  t8 = sp_ntt_sub_partial_simd(p2e, t8, p);

  t8 = sp_ntt_twiddle_mul_simd(t8, w + 0, p);
  t3 = sp_ntt_twiddle_mul_simd(t3, w + 2, p);
  t2 = sp_ntt_twiddle_mul_simd(t2, w + 4, p);
  t4 = sp_ntt_twiddle_mul_simd(t4, w + 6, p);
  t7 = sp_ntt_twiddle_mul_simd(t7, w + 8, p);
  t1 = sp_ntt_twiddle_mul_simd(t1, w + 10, p);
  t6 = sp_ntt_twiddle_mul_simd(t6, w + 12, p);
  t5 = sp_ntt_twiddle_mul_simd(t5, w + 14, p);

  sp_simd_scatter(p0, out + 0 * ostride, odist, vsize);
  sp_simd_scatter(t8, out + 1 * ostride, odist, vsize);
  sp_simd_scatter(t3, out + 2 * ostride, odist, vsize);
  sp_simd_scatter(t2, out + 3 * ostride, odist, vsize);
  sp_simd_scatter(t4, out + 4 * ostride, odist, vsize);
  sp_simd_scatter(t7, out + 5 * ostride, odist, vsize);
  sp_simd_scatter(t1, out + 6 * ostride, odist, vsize);
  sp_simd_scatter(t6, out + 7 * ostride, odist, vsize);
  sp_simd_scatter(t5, out + 8 * ostride, odist, vsize);
}

static void
ntt9_twiddle_run_core_simd_interleaved(
		spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		sp_simd_t *w, sp_simd_t p, sp_simd_t * c)
{
  sp_simd_t x0, x1, x2, x3, x4, x5, x6;
  sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
  sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

  sp_simd_t x0e, x1e, t0e, t1e, p0e, p1e, p2e;

  x0 = sp_simd_load(in + 0 * istride);
  x1 = sp_simd_load(in + 1 * istride);
  x2 = sp_simd_load(in + 2 * istride);
  x0e = sp_simd_load(in + 3 * istride);
  x3 = sp_simd_load(in + 4 * istride);
  x6 = sp_simd_load(in + 5 * istride);
  x1e = sp_simd_load(in + 6 * istride);
  x5 = sp_simd_load(in + 7 * istride);
  x4 = sp_simd_load(in + 8 * istride);

  t0e = sp_ntt_add_simd0(x0e, x1e, p);
  t1e = sp_ntt_sub_partial_simd0(x0e, x1e, p);

  p1 = sp_ntt_add_simd0(x1, x3, p);
  p1 = sp_ntt_add_simd0(p1, x5, p);
  p2 = sp_ntt_add_simd0(x2, x4, p);
  p2 = sp_ntt_add_simd0(p2, x6, p);
  p3 = sp_ntt_sub_simd0(x1, x5, p);
  p4 = sp_ntt_sub_simd0(x2, x6, p);
  p5 = sp_ntt_sub_simd0(x3, x5, p);
  p6 = sp_ntt_sub_simd0(x4, x6, p);

  t1 = sp_ntt_add_simd0(p1, p2, p);
  t2 = sp_ntt_sub_partial_simd0(p1, p2, p);
  t3 = sp_ntt_sub_simd0(p3, p5, p);
  t5 = sp_ntt_add_simd0(t3, p6, p);
  t3 = sp_ntt_sub_simd0(t3, p6, p);
  t4 = sp_ntt_add_simd0(p4, p5, p);
  t6 = sp_ntt_sub_simd0(p4, p5, p);

  p0e = sp_ntt_add_simd0(x0, t0e, p);
  p0 = t1;
  p1 = t1;
  p2 = t2;
  p3 = t3;
  p4 = t4;
  p5 = sp_ntt_add_partial_simd0(t3, t4, p);
  p6 = t5;
  p7 = t6;
  p8 = sp_ntt_add_partial_simd0(t5, t6, p);

  p1 = sp_ntt_mul_simd0(p1, c+2, p);
  p2 = sp_ntt_mul_simd0(p2, c+4, p);
  t0e = sp_ntt_mul_simd0(t0e, c+6, p);
  t1e = sp_ntt_mul_simd0(t1e, c+8, p);
  p3 = sp_ntt_mul_simd0(p3, c+10, p);
  p4 = sp_ntt_mul_simd0(p4, c+12, p);
  p5 = sp_ntt_mul_simd0(p5, c+14, p);
  p6 = sp_ntt_mul_simd0(p6, c+16, p);
  p7 = sp_ntt_mul_simd0(p7, c+18, p);
  p8 = sp_ntt_mul_simd0(p8, c+20, p);

  t0e = sp_ntt_add_simd0(t0e, p0e, p);
  t1 = sp_ntt_add_simd0(p1, p2, p);
  t2 = sp_ntt_sub_simd0(p1, p2, p);
  t3 = sp_ntt_add_simd0(p3, p5, p);
  t4 = sp_ntt_add_simd0(p4, p5, p);
  t5 = sp_ntt_add_simd0(p6, p8, p);
  t6 = sp_ntt_add_simd0(p7, p8, p);

  p1e = sp_ntt_add_simd0(t0e, t1e, p);
  p2e = sp_ntt_sub_simd0(t0e, t1e, p);
  p3 = sp_ntt_add_simd0(t3, t5, p);
  p4 = sp_ntt_add_simd0(t4, t6, p);
  p5 = sp_ntt_sub_simd0(t4, t6, p);
  p5 = sp_ntt_sub_simd0(p5, p3, p);
  p6 = sp_ntt_sub_simd0(t5, t3, p);

  p0 = sp_ntt_add_simd0(p0, p0e, p);
  t1 = sp_ntt_add_partial_simd0(t1, p0e, p);
  t2 = sp_ntt_add_partial_simd0(t2, p0e, p);
  t3 = sp_ntt_add_partial_simd0(p3, p1e, p);
  t4 = sp_ntt_add_partial_simd0(p4, p2e, p);
  t5 = sp_ntt_add_partial_simd0(p5, p1e, p);
  t6 = sp_ntt_add_partial_simd0(p6, p2e, p);
  t7 = sp_ntt_add_simd0(p3, p5, p);
  t7 = sp_ntt_sub_partial_simd0(p1e, t7, p);
  t8 = sp_ntt_add_simd0(p4, p6, p);
  t8 = sp_ntt_sub_partial_simd0(p2e, t8, p);

  t8 = sp_ntt_twiddle_mul_simd0(t8, w+0, p);
  t3 = sp_ntt_twiddle_mul_simd0(t3, w+2, p);
  t2 = sp_ntt_twiddle_mul_simd0(t2, w+4, p);
  t4 = sp_ntt_twiddle_mul_simd0(t4, w+6, p);
  t7 = sp_ntt_twiddle_mul_simd0(t7, w+8, p);
  t1 = sp_ntt_twiddle_mul_simd0(t1, w+10, p);
  t6 = sp_ntt_twiddle_mul_simd0(t6, w+12, p);
  t5 = sp_ntt_twiddle_mul_simd0(t5, w+14, p);

  sp_simd_store(p0, out + 0 * ostride);
  sp_simd_store(t8, out + 1 * ostride);
  sp_simd_store(t3, out + 2 * ostride);
  sp_simd_store(t2, out + 3 * ostride);
  sp_simd_store(t4, out + 4 * ostride);
  sp_simd_store(t7, out + 5 * ostride);
  sp_simd_store(t1, out + 6 * ostride);
  sp_simd_store(t6, out + 7 * ostride);
  sp_simd_store(t5, out + 8 * ostride);
}


static void
ntt9_pfa_run_core_simd(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t inc2, spv_size_t n,
	  sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  spv_size_t j0, j1, j2, j3, j4, j5, j6, j7, j8;
  sp_simd_t x0, x1, x2, x3, x4, x5, x6;
  sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
  sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

  sp_simd_t x0e, x1e, t0e, t1e, p0e, p1e, p2e;

  j0 = start;
  j1 = sp_array_inc(j0, inc, n);
  j2 = sp_array_inc(j0, 2 * inc, n);
  j3 = sp_array_inc(j0, 3 * inc, n);
  j4 = sp_array_inc(j0, 4 * inc, n);
  j5 = sp_array_inc(j0, 5 * inc, n);
  j6 = sp_array_inc(j0, 6 * inc, n);
  j7 = sp_array_inc(j0, 7 * inc, n);
  j8 = sp_array_inc(j0, 8 * inc, n);

  x0 = sp_simd_pfa_gather(x, j0, inc2, n, vsize);
  x1 = sp_simd_pfa_gather(x, j1, inc2, n, vsize);
  x2 = sp_simd_pfa_gather(x, j2, inc2, n, vsize);
  x0e = sp_simd_pfa_gather(x, j3, inc2, n, vsize);
  x3 = sp_simd_pfa_gather(x, j4, inc2, n, vsize);
  x6 = sp_simd_pfa_gather(x, j5, inc2, n, vsize);
  x1e = sp_simd_pfa_gather(x, j6, inc2, n, vsize);
  x5 = sp_simd_pfa_gather(x, j7, inc2, n, vsize);
  x4 = sp_simd_pfa_gather(x, j8, inc2, n, vsize);

  t0e = sp_ntt_add_simd(x0e, x1e, p);
  t1e = sp_ntt_sub_partial_simd(x0e, x1e, p);

  p1 = sp_ntt_add_simd(x1, x3, p);
  p1 = sp_ntt_add_simd(p1, x5, p);
  p2 = sp_ntt_add_simd(x2, x4, p);
  p2 = sp_ntt_add_simd(p2, x6, p);
  p3 = sp_ntt_sub_simd(x1, x5, p);
  p4 = sp_ntt_sub_simd(x2, x6, p);
  p5 = sp_ntt_sub_simd(x3, x5, p);
  p6 = sp_ntt_sub_simd(x4, x6, p);

  t1 = sp_ntt_add_simd(p1, p2, p);
  t2 = sp_ntt_sub_partial_simd(p1, p2, p);
  t3 = sp_ntt_sub_simd(p3, p5, p);
  t5 = sp_ntt_add_simd(t3, p6, p);
  t3 = sp_ntt_sub_simd(t3, p6, p);
  t4 = sp_ntt_add_simd(p4, p5, p);
  t6 = sp_ntt_sub_simd(p4, p5, p);

  p0e = sp_ntt_add_simd(x0, t0e, p);
  p0 = t1;
  p1 = t1;
  p2 = t2;
  p3 = t3;
  p4 = t4;
  p5 = sp_ntt_add_partial_simd(t3, t4, p);
  p6 = t5;
  p7 = t6;
  p8 = sp_ntt_add_partial_simd(t5, t6, p);

  p1 = sp_ntt_mul_simd(p1, ntt_const[1], ntt_const[NC+1], p);
  p2 = sp_ntt_mul_simd(p2, ntt_const[2], ntt_const[NC+2], p);
  t0e = sp_ntt_mul_simd(t0e, ntt_const[3], ntt_const[NC+3], p);
  t1e = sp_ntt_mul_simd(t1e, ntt_const[4], ntt_const[NC+4], p);
  p3 = sp_ntt_mul_simd(p3, ntt_const[5], ntt_const[NC+5], p);
  p4 = sp_ntt_mul_simd(p4, ntt_const[6], ntt_const[NC+6], p);
  p5 = sp_ntt_mul_simd(p5, ntt_const[7], ntt_const[NC+7], p);
  p6 = sp_ntt_mul_simd(p6, ntt_const[8], ntt_const[NC+8], p);
  p7 = sp_ntt_mul_simd(p7, ntt_const[9], ntt_const[NC+9], p);
  p8 = sp_ntt_mul_simd(p8, ntt_const[10], ntt_const[NC+10], p);

  t0e = sp_ntt_add_simd(t0e, p0e, p);
  t1 = sp_ntt_add_simd(p1, p2, p);
  t2 = sp_ntt_sub_simd(p1, p2, p);
  t3 = sp_ntt_add_simd(p3, p5, p);
  t4 = sp_ntt_add_simd(p4, p5, p);
  t5 = sp_ntt_add_simd(p6, p8, p);
  t6 = sp_ntt_add_simd(p7, p8, p);

  p1e = sp_ntt_add_simd(t0e, t1e, p);
  p2e = sp_ntt_sub_simd(t0e, t1e, p);
  p3 = sp_ntt_add_simd(t3, t5, p);
  p4 = sp_ntt_add_simd(t4, t6, p);
  p5 = sp_ntt_sub_simd(t4, t6, p);
  p5 = sp_ntt_sub_simd(p5, p3, p);
  p6 = sp_ntt_sub_simd(t5, t3, p);

  p0 = sp_ntt_add_simd(p0, p0e, p);
  t1 = sp_ntt_add_simd(t1, p0e, p);
  t2 = sp_ntt_add_simd(t2, p0e, p);
  t3 = sp_ntt_add_simd(p3, p1e, p);
  t4 = sp_ntt_add_simd(p4, p2e, p);
  t5 = sp_ntt_add_simd(p5, p1e, p);
  t6 = sp_ntt_add_simd(p6, p2e, p);
  t7 = sp_ntt_add_simd(p3, p5, p);
  t7 = sp_ntt_sub_simd(p1e, t7, p);
  t8 = sp_ntt_add_simd(p4, p6, p);
  t8 = sp_ntt_sub_simd(p2e, t8, p);

  sp_simd_pfa_scatter(p0, x, j0, inc2, n, vsize);
  sp_simd_pfa_scatter(t8, x, j1, inc2, n, vsize);
  sp_simd_pfa_scatter(t3, x, j2, inc2, n, vsize);
  sp_simd_pfa_scatter(t2, x, j3, inc2, n, vsize);
  sp_simd_pfa_scatter(t4, x, j4, inc2, n, vsize);
  sp_simd_pfa_scatter(t7, x, j5, inc2, n, vsize);
  sp_simd_pfa_scatter(t1, x, j6, inc2, n, vsize);
  sp_simd_pfa_scatter(t6, x, j7, inc2, n, vsize);
  sp_simd_pfa_scatter(t5, x, j8, inc2, n, vsize);
}

static void
ntt9_pfa_run_core_simd_interleaved(
	  spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t n,
	  sp_simd_t p, sp_simd_t * c)
{
  spv_size_t j0, j1, j2, j3, j4, j5, j6, j7, j8;
  sp_simd_t x0, x1, x2, x3, x4, x5, x6;
  sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
  sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

  sp_simd_t x0e, x1e, t0e, t1e, p0e, p1e, p2e;

  j0 = start;
  j1 = sp_array_inc(j0, inc, n);
  j2 = sp_array_inc(j0, 2 * inc, n);
  j3 = sp_array_inc(j0, 3 * inc, n);
  j4 = sp_array_inc(j0, 4 * inc, n);
  j5 = sp_array_inc(j0, 5 * inc, n);
  j6 = sp_array_inc(j0, 6 * inc, n);
  j7 = sp_array_inc(j0, 7 * inc, n);
  j8 = sp_array_inc(j0, 8 * inc, n);

  x0 = sp_simd_load(x + j0);
  x1 = sp_simd_load(x + j1);
  x2 = sp_simd_load(x + j2);
  x0e = sp_simd_load(x + j3);
  x3 = sp_simd_load(x + j4);
  x6 = sp_simd_load(x + j5);
  x1e = sp_simd_load(x + j6);
  x5 = sp_simd_load(x + j7);
  x4 = sp_simd_load(x + j8);

  t0e = sp_ntt_add_simd0(x0e, x1e, p);
  t1e = sp_ntt_sub_partial_simd0(x0e, x1e, p);

  p1 = sp_ntt_add_simd0(x1, x3, p);
  p1 = sp_ntt_add_simd0(p1, x5, p);
  p2 = sp_ntt_add_simd0(x2, x4, p);
  p2 = sp_ntt_add_simd0(p2, x6, p);
  p3 = sp_ntt_sub_simd0(x1, x5, p);
  p4 = sp_ntt_sub_simd0(x2, x6, p);
  p5 = sp_ntt_sub_simd0(x3, x5, p);
  p6 = sp_ntt_sub_simd0(x4, x6, p);

  t1 = sp_ntt_add_simd0(p1, p2, p);
  t2 = sp_ntt_sub_partial_simd0(p1, p2, p);
  t3 = sp_ntt_sub_simd0(p3, p5, p);
  t5 = sp_ntt_add_simd0(t3, p6, p);
  t3 = sp_ntt_sub_simd0(t3, p6, p);
  t4 = sp_ntt_add_simd0(p4, p5, p);
  t6 = sp_ntt_sub_simd0(p4, p5, p);

  p0e = sp_ntt_add_simd0(x0, t0e, p);
  p0 = t1;
  p1 = t1;
  p2 = t2;
  p3 = t3;
  p4 = t4;
  p5 = sp_ntt_add_partial_simd0(t3, t4, p);
  p6 = t5;
  p7 = t6;
  p8 = sp_ntt_add_partial_simd0(t5, t6, p);

  p1 = sp_ntt_mul_simd0(p1, c+2, p);
  p2 = sp_ntt_mul_simd0(p2, c+4, p);
  t0e = sp_ntt_mul_simd0(t0e, c+6, p);
  t1e = sp_ntt_mul_simd0(t1e, c+8, p);
  p3 = sp_ntt_mul_simd0(p3, c+10, p);
  p4 = sp_ntt_mul_simd0(p4, c+12, p);
  p5 = sp_ntt_mul_simd0(p5, c+14, p);
  p6 = sp_ntt_mul_simd0(p6, c+16, p);
  p7 = sp_ntt_mul_simd0(p7, c+18, p);
  p8 = sp_ntt_mul_simd0(p8, c+20, p);

  t0e = sp_ntt_add_simd0(t0e, p0e, p);
  t1 = sp_ntt_add_simd0(p1, p2, p);
  t2 = sp_ntt_sub_simd0(p1, p2, p);
  t3 = sp_ntt_add_simd0(p3, p5, p);
  t4 = sp_ntt_add_simd0(p4, p5, p);
  t5 = sp_ntt_add_simd0(p6, p8, p);
  t6 = sp_ntt_add_simd0(p7, p8, p);

  p1e = sp_ntt_add_simd0(t0e, t1e, p);
  p2e = sp_ntt_sub_simd0(t0e, t1e, p);
  p3 = sp_ntt_add_simd0(t3, t5, p);
  p4 = sp_ntt_add_simd0(t4, t6, p);
  p5 = sp_ntt_sub_simd0(t4, t6, p);
  p5 = sp_ntt_sub_simd0(p5, p3, p);
  p6 = sp_ntt_sub_simd0(t5, t3, p);

  p0 = sp_ntt_add_simd0(p0, p0e, p);
  t1 = sp_ntt_add_simd0(t1, p0e, p);
  t2 = sp_ntt_add_simd0(t2, p0e, p);
  t3 = sp_ntt_add_simd0(p3, p1e, p);
  t4 = sp_ntt_add_simd0(p4, p2e, p);
  t5 = sp_ntt_add_simd0(p5, p1e, p);
  t6 = sp_ntt_add_simd0(p6, p2e, p);
  t7 = sp_ntt_add_simd0(p3, p5, p);
  t7 = sp_ntt_sub_simd0(p1e, t7, p);
  t8 = sp_ntt_add_simd0(p4, p6, p);
  t8 = sp_ntt_sub_simd0(p2e, t8, p);

  sp_simd_store(p0, x + j0);
  sp_simd_store(t8, x + j1);
  sp_simd_store(t3, x + j2);
  sp_simd_store(t2, x + j3);
  sp_simd_store(t4, x + j4);
  sp_simd_store(t7, x + j5);
  sp_simd_store(t1, x + j6);
  sp_simd_store(t6, x + j7);
  sp_simd_store(t5, x + j8);
}

DECLARE_CORE_ROUTINES_SIMD(9)
