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
ntt2_run_simd(spv_t in, spv_size_t istride, spv_size_t idist,
    		spv_t out, spv_size_t ostride, spv_size_t odist,
    		spv_size_t num_transforms, sp_t p, spv_t ntt_const)
{
  spv_size_t i = 0;

  spv_size_t num_simd = SP_SIMD_VSIZE * (num_transforms / SP_SIMD_VSIZE);

  for (; i < num_simd; i += SP_SIMD_VSIZE)
    ntt2_run_core_simd(in + i * idist, istride, idist, 
                        out + i * odist, ostride, odist, p, 
                        ntt_const, SP_SIMD_VSIZE);

  if (i < num_transforms)
    ntt2_run_core_simd(in + i * idist, istride, idist, 
                        out + i * odist, ostride, odist, p, 
                        ntt_const, num_transforms - i);
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
ntt2_twiddle_run_simd(spv_t in, spv_size_t istride, spv_size_t idist,
    			spv_t out, spv_size_t ostride, spv_size_t odist,
    			spv_t w, spv_size_t num_transforms, sp_t p, spv_t ntt_const)
{
  spv_size_t i = 0, j = 0;

  spv_size_t num_simd = SP_SIMD_VSIZE * (num_transforms / SP_SIMD_VSIZE);

  for (; i < num_simd; i += SP_SIMD_VSIZE,
		  	j += 2*(2-1)*SP_SIMD_VSIZE)
    ntt2_twiddle_run_core_simd(
		in + i * idist, istride, idist,
		out + i * odist, ostride, odist,
		(sp_simd_t *)(w + j), p, ntt_const,
		SP_SIMD_VSIZE);

  if (i < num_transforms)
    ntt2_twiddle_run_core_simd(
		in + i * idist, istride, idist,
		out + i * odist, ostride, odist,
		(sp_simd_t *)(w + j), p, ntt_const,
		num_transforms - i);
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
ntt2_pfa_run_simd(spv_t x, spv_size_t cofactor,
	  sp_t p, spv_t ntt_const)
{
  spv_size_t i = 0;
  spv_size_t incstart = 0;
  spv_size_t n = 2 * cofactor;
  spv_size_t inc = cofactor;
  spv_size_t inc2 = 2;

  spv_size_t num_simd = SP_SIMD_VSIZE * (cofactor / SP_SIMD_VSIZE);

  for (i = 0; i < num_simd; i += SP_SIMD_VSIZE)
    {
      ntt2_pfa_run_core_simd(x, incstart, inc, inc2, 
	  			n, p, ntt_const, SP_SIMD_VSIZE);
      incstart += SP_SIMD_VSIZE * inc2;
    }

  if (i < cofactor)
    ntt2_pfa_run_core_simd(x, incstart, inc, inc2, 
  			n, p, ntt_const, cofactor - i);

}

const nttconfig_t V(ntt2simd_config) = 
{
  2,
  NC,
  ntt2_fixed_const,
  X(ntt2_init),
  ntt2_run_simd,
  ntt2_pfa_run_simd,
  ntt2_twiddle_run_simd
};

