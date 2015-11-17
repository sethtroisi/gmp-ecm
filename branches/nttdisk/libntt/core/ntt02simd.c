#include "ntt/ntt-impl-simd.h"

#define NC 2

static const uint8_t ntt2_fixed_const[NC] = {1, 1};

extern void X(ntt2_init)(spv_t out, sp_t p, sp_t d, sp_t primroot, 
                        sp_t order, sp_t perm);


static void
ntt2_run_core_simd(spv_t in, spv_size_t istride, spv_size_t idist,
		spv_t out, spv_size_t ostride, spv_size_t odist,
		sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  sp_simd_t p0, p1;
  sp_simd_t x0, x1;

  x0 = sp_simd_gather(in + 0 * istride, idist, vsize);
  x1 = sp_simd_gather(in + 1 * istride, idist, vsize);

  p0 = sp_ntt_add_simd(x0, x1, p);
  p1 = sp_ntt_sub_simd(x0, x1, p);

  sp_simd_scatter(p0, out + 0 * ostride, odist, vsize);
  sp_simd_scatter(p1, out + 1 * ostride, odist, vsize);
}

static void
ntt2_run_core_simd_interleaved(
        spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		sp_simd_t p, sp_simd_t * c)
{
  sp_simd_t p0, p1;
  sp_simd_t x0, x1;

  x0 = sp_simd_load(in + 0 * istride);
  x1 = sp_simd_load(in + 1 * istride);

  p0 = sp_ntt_add_simd0(x0, x1, p);
  p1 = sp_ntt_sub_simd0(x0, x1, p);

  sp_simd_store(p0, out + 0 * ostride);
  sp_simd_store(p1, out + 1 * ostride);
}

static void
ntt2_twiddle_run_core_simd(
        spv_t in, spv_size_t istride, spv_size_t idist,
		spv_t out, spv_size_t ostride, spv_size_t odist,
		sp_simd_t *w, sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  sp_simd_t p0, p1;
  sp_simd_t x0, x1;

  x0 = sp_simd_gather(in + 0 * istride, idist, vsize);
  x1 = sp_simd_gather(in + 1 * istride, idist, vsize);

  p0 = sp_ntt_add_simd(x0, x1, p);
  p1 = sp_ntt_sub_partial_simd(x0, x1, p);

  p1 = sp_ntt_twiddle_mul_simd(p1, w + 0, p);

  sp_simd_scatter(p0, out + 0 * ostride, odist, vsize);
  sp_simd_scatter(p1, out + 1 * ostride, odist, vsize);
}

static void
ntt2_twiddle_run_core_simd_interleaved(
        spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		sp_simd_t *w, sp_simd_t p, sp_simd_t * c)
{
  sp_simd_t p0, p1;
  sp_simd_t x0, x1;

  x0 = sp_simd_load(in + 0 * istride);
  x1 = sp_simd_load(in + 1 * istride);

  p0 = sp_ntt_add_simd0(x0, x1, p);
  p1 = sp_ntt_sub_partial_simd0(x0, x1, p);

  p1 = sp_ntt_twiddle_mul_simd0(p1, w+0, p);

  sp_simd_store(p0, out + 0 * ostride);
  sp_simd_store(p1, out + 1 * ostride);
}


static void
ntt2_pfa_run_core_simd(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t inc2, spv_size_t n,
	  sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  spv_size_t j0, j1;
  sp_simd_t p0, p1;
  sp_simd_t x0, x1;

  j0 = start;
  j1 = sp_array_inc(j0, inc, n);

  x0 = sp_simd_pfa_gather(x, j0, inc2, n, vsize);
  x1 = sp_simd_pfa_gather(x, j1, inc2, n, vsize);

  p0 = sp_ntt_add_simd(x0, x1, p);
  p1 = sp_ntt_sub_simd(x0, x1, p);

  sp_simd_pfa_scatter(p0, x, j0, inc2, n, vsize);
  sp_simd_pfa_scatter(p1, x, j1, inc2, n, vsize);
}


static void
ntt2_pfa_run_core_simd_interleaved(
          spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t n,
	  sp_simd_t p, sp_simd_t * c)
{
  spv_size_t j0, j1;
  sp_simd_t p0, p1;
  sp_simd_t x0, x1;

  j0 = start;
  j1 = sp_array_inc(j0, inc, n);

  x0 = sp_simd_load(x + j0);
  x1 = sp_simd_load(x + j1);

  p0 = sp_ntt_add_simd0(x0, x1, p);
  p1 = sp_ntt_sub_simd0(x0, x1, p);

  sp_simd_store(p0, x + j0);
  sp_simd_store(p1, x + j1);
}

DECLARE_CORE_ROUTINES_SIMD(2)
