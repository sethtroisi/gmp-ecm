#include "ntt/ntt-impl-simd.h"

#define NC 9

static const uint8_t ntt7_fixed_const[NC] = {1};

extern void X(ntt7_init)(spv_t out, sp_t p, sp_t d, sp_t primroot, 
                        sp_t order, sp_t perm);


static void
ntt7_run_core_simd(spv_t in, spv_size_t istride, spv_size_t idist,
		spv_t out, spv_size_t ostride, spv_size_t odist,
		sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  sp_simd_t x0, x1, x2, x3, x4, x5, x6;
  sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
  sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

  x0 = sp_simd_gather(in + 0 * istride, idist, vsize);
  x1 = sp_simd_gather(in + 1 * istride, idist, vsize);
  x5 = sp_simd_gather(in + 2 * istride, idist, vsize);
  x6 = sp_simd_gather(in + 3 * istride, idist, vsize);
  x3 = sp_simd_gather(in + 4 * istride, idist, vsize);
  x2 = sp_simd_gather(in + 5 * istride, idist, vsize);
  x4 = sp_simd_gather(in + 6 * istride, idist, vsize);

  p1 = sp_ntt_add_simd(x1, x4, p);
  p2 = sp_ntt_add_simd(x2, x5, p);
  p3 = sp_ntt_add_simd(x3, x6, p);
  p4 = sp_ntt_sub_simd(x1, x4, p);
  p5 = sp_ntt_sub_simd(x2, x5, p);
  p6 = sp_ntt_sub_simd(x3, x6, p);

  t1 = sp_ntt_add_simd(p1, p2, p);
  t1 = sp_ntt_add_simd(t1, p3, p);
  t2 = sp_ntt_sub_simd(p1, p3, p);
  t3 = sp_ntt_sub_simd(p2, p3, p);
  t4 = sp_ntt_sub_simd(p4, p5, p);
  t4 = sp_ntt_add_partial_simd(t4, p6, p);
  t5 = sp_ntt_sub_simd(p4, p6, p);
  t6 = sp_ntt_add_simd(p5, p6, p);

  p0 = sp_ntt_add_simd(x0, t1, p);
  p1 = t1;
  p2 = t2;
  p3 = t3;
  p4 = sp_ntt_add_partial_simd(t2, t3, p);
  p5 = t4;
  p6 = t5;
  p7 = t6;
  p8 = sp_ntt_add_partial_simd(t5, t6, p);

  p1 = sp_ntt_mul_simd(p1, ntt_const[1], ntt_const[NC+1], p);
  p2 = sp_ntt_mul_simd(p2, ntt_const[2], ntt_const[NC+2], p);
  p3 = sp_ntt_mul_simd(p3, ntt_const[3], ntt_const[NC+3], p);
  p4 = sp_ntt_mul_simd(p4, ntt_const[4], ntt_const[NC+4], p);
  p5 = sp_ntt_mul_simd(p5, ntt_const[5], ntt_const[NC+5], p);
  p6 = sp_ntt_mul_simd(p6, ntt_const[6], ntt_const[NC+6], p);
  p7 = sp_ntt_mul_simd(p7, ntt_const[7], ntt_const[NC+7], p);
  p8 = sp_ntt_mul_simd(p8, ntt_const[8], ntt_const[NC+8], p);

  t1 = sp_ntt_add_simd(p0, p1, p);
  t2 = sp_ntt_add_simd(p2, p4, p);
  t3 = sp_ntt_add_simd(p3, p4, p);
  t4 = p5;
  t5 = sp_ntt_add_simd(p6, p8, p);
  t6 = sp_ntt_add_simd(p7, p8, p);

  p1 = sp_ntt_add_simd(t1, t2, p);
  p2 = sp_ntt_add_simd(t1, t3, p);
  p3 = sp_ntt_sub_simd(t1, t2, p);
  p3 = sp_ntt_sub_simd(p3, t3, p);
  p4 = sp_ntt_add_simd(t4, t5, p);
  p5 = sp_ntt_sub_simd(t6, t4, p);
  p6 = sp_ntt_sub_simd(t4, t5, p);
  p6 = sp_ntt_add_simd(p6, t6, p);

  t1 = sp_ntt_add_simd(p1, p4, p);
  t2 = sp_ntt_add_simd(p2, p5, p);
  t3 = sp_ntt_add_simd(p3, p6, p);
  t4 = sp_ntt_sub_simd(p1, p4, p);
  t5 = sp_ntt_sub_simd(p2, p5, p);
  t6 = sp_ntt_sub_simd(p3, p6, p);

  sp_simd_scatter(p0, out + 0 * ostride, odist, vsize);
  sp_simd_scatter(t6, out + 1 * ostride, odist, vsize);
  sp_simd_scatter(t4, out + 2 * ostride, odist, vsize);
  sp_simd_scatter(t5, out + 3 * ostride, odist, vsize);
  sp_simd_scatter(t2, out + 4 * ostride, odist, vsize);
  sp_simd_scatter(t1, out + 5 * ostride, odist, vsize);
  sp_simd_scatter(t3, out + 6 * ostride, odist, vsize);
}


static void
ntt7_run_simd(spv_t in, spv_size_t istride, spv_size_t idist,
    		spv_t out, spv_size_t ostride, spv_size_t odist,
    		spv_size_t num_transforms, sp_t p, spv_t ntt_const)
{
  spv_size_t i = 0;

  spv_size_t num_simd = SP_SIMD_VSIZE * (num_transforms / SP_SIMD_VSIZE);

  for (; i < num_simd; i += SP_SIMD_VSIZE)
    ntt7_run_core_simd(in + i * idist, istride, idist, 
                        out + i * odist, ostride, odist, p, 
			ntt_const, SP_SIMD_VSIZE);

  if (i < num_transforms)
    ntt7_run_core_simd(in + i * idist, istride, idist, 
                        out + i * odist, ostride, odist, p, 
			ntt_const, num_transforms - i);
}


static void
ntt7_twiddle_run_core_simd(
        spv_t in, spv_size_t istride, spv_size_t idist,
		spv_t out, spv_size_t ostride, spv_size_t odist,
		sp_simd_t *w, sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  sp_simd_t x0, x1, x2, x3, x4, x5, x6;
  sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
  sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

  x0 = sp_simd_gather(in + 0 * istride, idist, vsize);
  x1 = sp_simd_gather(in + 1 * istride, idist, vsize);
  x5 = sp_simd_gather(in + 2 * istride, idist, vsize);
  x6 = sp_simd_gather(in + 3 * istride, idist, vsize);
  x3 = sp_simd_gather(in + 4 * istride, idist, vsize);
  x2 = sp_simd_gather(in + 5 * istride, idist, vsize);
  x4 = sp_simd_gather(in + 6 * istride, idist, vsize);

  p1 = sp_ntt_add_simd(x1, x4, p);
  p2 = sp_ntt_add_simd(x2, x5, p);
  p3 = sp_ntt_add_simd(x3, x6, p);
  p4 = sp_ntt_sub_simd(x1, x4, p);
  p5 = sp_ntt_sub_simd(x2, x5, p);
  p6 = sp_ntt_sub_simd(x3, x6, p);

  t1 = sp_ntt_add_simd(p1, p2, p);
  t1 = sp_ntt_add_simd(t1, p3, p);
  t2 = sp_ntt_sub_simd(p1, p3, p);
  t3 = sp_ntt_sub_simd(p2, p3, p);
  t4 = sp_ntt_sub_simd(p4, p5, p);
  t4 = sp_ntt_add_partial_simd(t4, p6, p);
  t5 = sp_ntt_sub_simd(p4, p6, p);
  t6 = sp_ntt_add_simd(p5, p6, p);

  p0 = sp_ntt_add_simd(x0, t1, p);
  p1 = t1;
  p2 = t2;
  p3 = t3;
  p4 = sp_ntt_add_partial_simd(t2, t3, p);
  p5 = t4;
  p6 = t5;
  p7 = t6;
  p8 = sp_ntt_add_partial_simd(t5, t6, p);

  p1 = sp_ntt_mul_simd(p1, ntt_const[1], ntt_const[NC+1], p);
  p2 = sp_ntt_mul_simd(p2, ntt_const[2], ntt_const[NC+2], p);
  p3 = sp_ntt_mul_simd(p3, ntt_const[3], ntt_const[NC+3], p);
  p4 = sp_ntt_mul_simd(p4, ntt_const[4], ntt_const[NC+4], p);
  p5 = sp_ntt_mul_simd(p5, ntt_const[5], ntt_const[NC+5], p);
  p6 = sp_ntt_mul_simd(p6, ntt_const[6], ntt_const[NC+6], p);
  p7 = sp_ntt_mul_simd(p7, ntt_const[7], ntt_const[NC+7], p);
  p8 = sp_ntt_mul_simd(p8, ntt_const[8], ntt_const[NC+8], p);

  t1 = sp_ntt_add_simd(p0, p1, p);
  t2 = sp_ntt_add_simd(p2, p4, p);
  t3 = sp_ntt_add_simd(p3, p4, p);
  t4 = p5;
  t5 = sp_ntt_add_simd(p6, p8, p);
  t6 = sp_ntt_add_simd(p7, p8, p);

  p1 = sp_ntt_add_simd(t1, t2, p);
  p2 = sp_ntt_add_simd(t1, t3, p);
  p3 = sp_ntt_sub_simd(t1, t2, p);
  p3 = sp_ntt_sub_simd(p3, t3, p);
  p4 = sp_ntt_add_simd(t4, t5, p);
  p5 = sp_ntt_sub_simd(t6, t4, p);
  p6 = sp_ntt_sub_simd(t4, t5, p);
  p6 = sp_ntt_add_simd(p6, t6, p);

  t1 = sp_ntt_add_partial_simd(p1, p4, p);
  t2 = sp_ntt_add_partial_simd(p2, p5, p);
  t3 = sp_ntt_add_partial_simd(p3, p6, p);
  t4 = sp_ntt_sub_partial_simd(p1, p4, p);
  t5 = sp_ntt_sub_partial_simd(p2, p5, p);
  t6 = sp_ntt_sub_partial_simd(p3, p6, p);

  t6 = sp_ntt_twiddle_mul_simd(t6, w + 0, p);
  t4 = sp_ntt_twiddle_mul_simd(t4, w + 2, p);
  t5 = sp_ntt_twiddle_mul_simd(t5, w + 4, p);
  t2 = sp_ntt_twiddle_mul_simd(t2, w + 6, p);
  t1 = sp_ntt_twiddle_mul_simd(t1, w + 8, p);
  t3 = sp_ntt_twiddle_mul_simd(t3, w + 10, p);

  sp_simd_scatter(p0, out + 0 * ostride, odist, vsize);
  sp_simd_scatter(t6, out + 1 * ostride, odist, vsize);
  sp_simd_scatter(t4, out + 2 * ostride, odist, vsize);
  sp_simd_scatter(t5, out + 3 * ostride, odist, vsize);
  sp_simd_scatter(t2, out + 4 * ostride, odist, vsize);
  sp_simd_scatter(t1, out + 5 * ostride, odist, vsize);
  sp_simd_scatter(t3, out + 6 * ostride, odist, vsize);
}


static void
ntt7_twiddle_run_simd(spv_t in, spv_size_t istride, spv_size_t idist,
    			spv_t out, spv_size_t ostride, spv_size_t odist,
    			spv_t w, spv_size_t num_transforms, sp_t p, spv_t ntt_const)
{
  spv_size_t i = 0, j = 0;

  spv_size_t num_simd = SP_SIMD_VSIZE * (num_transforms / SP_SIMD_VSIZE);

  for (; i < num_simd; i += SP_SIMD_VSIZE,
		  	j += 2*(7-1)*SP_SIMD_VSIZE)
    ntt7_twiddle_run_core_simd(
		in + i * idist, istride, idist,
		out + i * odist, ostride, odist,
		(sp_simd_t *)(w + j), p, 
		ntt_const, SP_SIMD_VSIZE);

  if (i < num_transforms)
    ntt7_twiddle_run_core_simd(
		in + i * idist, istride, idist,
		out + i * odist, ostride, odist,
		(sp_simd_t *)(w + j), p, 
		ntt_const, num_transforms - i);
}


static void
ntt7_pfa_run_core_simd(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t inc2, spv_size_t n,
	  sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  spv_size_t j0, j1, j2, j3, j4, j5, j6;
  sp_simd_t x0, x1, x2, x3, x4, x5, x6;
  sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
  sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

  j0 = start;
  j1 = sp_array_inc(j0, inc, n);
  j2 = sp_array_inc(j0, 2 * inc, n);
  j3 = sp_array_inc(j0, 3 * inc, n);
  j4 = sp_array_inc(j0, 4 * inc, n);
  j5 = sp_array_inc(j0, 5 * inc, n);
  j6 = sp_array_inc(j0, 6 * inc, n);

  x0 = sp_simd_pfa_gather(x, j0, inc2, n, vsize);
  x1 = sp_simd_pfa_gather(x, j1, inc2, n, vsize);
  x5 = sp_simd_pfa_gather(x, j2, inc2, n, vsize);
  x6 = sp_simd_pfa_gather(x, j3, inc2, n, vsize);
  x3 = sp_simd_pfa_gather(x, j4, inc2, n, vsize);
  x2 = sp_simd_pfa_gather(x, j5, inc2, n, vsize);
  x4 = sp_simd_pfa_gather(x, j6, inc2, n, vsize);

  p1 = sp_ntt_add_simd(x1, x4, p);
  p2 = sp_ntt_add_simd(x2, x5, p);
  p3 = sp_ntt_add_simd(x3, x6, p);
  p4 = sp_ntt_sub_simd(x1, x4, p);
  p5 = sp_ntt_sub_simd(x2, x5, p);
  p6 = sp_ntt_sub_simd(x3, x6, p);

  t1 = sp_ntt_add_simd(p1, p2, p);
  t1 = sp_ntt_add_simd(t1, p3, p);
  t2 = sp_ntt_sub_simd(p1, p3, p);
  t3 = sp_ntt_sub_simd(p2, p3, p);
  t4 = sp_ntt_sub_simd(p4, p5, p);
  t4 = sp_ntt_add_partial_simd(t4, p6, p);
  t5 = sp_ntt_sub_simd(p4, p6, p);
  t6 = sp_ntt_add_simd(p5, p6, p);

  p0 = sp_ntt_add_simd(x0, t1, p);
  p1 = t1;
  p2 = t2;
  p3 = t3;
  p4 = sp_ntt_add_partial_simd(t2, t3, p);
  p5 = t4;
  p6 = t5;
  p7 = t6;
  p8 = sp_ntt_add_partial_simd(t5, t6, p);

  p1 = sp_ntt_mul_simd(p1, ntt_const[1], ntt_const[NC+1], p);
  p2 = sp_ntt_mul_simd(p2, ntt_const[2], ntt_const[NC+2], p);
  p3 = sp_ntt_mul_simd(p3, ntt_const[3], ntt_const[NC+3], p);
  p4 = sp_ntt_mul_simd(p4, ntt_const[4], ntt_const[NC+4], p);
  p5 = sp_ntt_mul_simd(p5, ntt_const[5], ntt_const[NC+5], p);
  p6 = sp_ntt_mul_simd(p6, ntt_const[6], ntt_const[NC+6], p);
  p7 = sp_ntt_mul_simd(p7, ntt_const[7], ntt_const[NC+7], p);
  p8 = sp_ntt_mul_simd(p8, ntt_const[8], ntt_const[NC+8], p);

  t1 = sp_ntt_add_simd(p0, p1, p);
  t2 = sp_ntt_add_simd(p2, p4, p);
  t3 = sp_ntt_add_simd(p3, p4, p);
  t4 = p5;
  t5 = sp_ntt_add_simd(p6, p8, p);
  t6 = sp_ntt_add_simd(p7, p8, p);

  p1 = sp_ntt_add_simd(t1, t2, p);
  p2 = sp_ntt_add_simd(t1, t3, p);
  p3 = sp_ntt_sub_simd(t1, t2, p);
  p3 = sp_ntt_sub_simd(p3, t3, p);
  p4 = sp_ntt_add_simd(t4, t5, p);
  p5 = sp_ntt_sub_simd(t6, t4, p);
  p6 = sp_ntt_sub_simd(t4, t5, p);
  p6 = sp_ntt_add_simd(p6, t6, p);

  t1 = sp_ntt_add_simd(p1, p4, p);
  t2 = sp_ntt_add_simd(p2, p5, p);
  t3 = sp_ntt_add_simd(p3, p6, p);
  t4 = sp_ntt_sub_simd(p1, p4, p);
  t5 = sp_ntt_sub_simd(p2, p5, p);
  t6 = sp_ntt_sub_simd(p3, p6, p);

  sp_simd_pfa_scatter(p0, x, j0, inc2, n, vsize);
  sp_simd_pfa_scatter(t6, x, j1, inc2, n, vsize);
  sp_simd_pfa_scatter(t4, x, j2, inc2, n, vsize);
  sp_simd_pfa_scatter(t5, x, j3, inc2, n, vsize);
  sp_simd_pfa_scatter(t2, x, j4, inc2, n, vsize);
  sp_simd_pfa_scatter(t1, x, j5, inc2, n, vsize);
  sp_simd_pfa_scatter(t3, x, j6, inc2, n, vsize);
}


static void
ntt7_pfa_run_simd(spv_t x, spv_size_t cofactor,
	  sp_t p, spv_t ntt_const)
{
  spv_size_t i = 0;
  spv_size_t incstart = 0;
  spv_size_t n = 7 * cofactor;
  spv_size_t inc = cofactor;
  spv_size_t inc2 = 7;

  spv_size_t num_simd = SP_SIMD_VSIZE * (cofactor / SP_SIMD_VSIZE);

  for (i = 0; i < num_simd; i += SP_SIMD_VSIZE)
    {
      ntt7_pfa_run_core_simd(x, incstart, inc, inc2, n, p, 
	  			ntt_const, SP_SIMD_VSIZE);
      incstart += SP_SIMD_VSIZE * inc2;
    }

  if (i < cofactor)
    ntt7_pfa_run_core_simd(x, incstart, inc, inc2, n, p, 
	  			ntt_const, cofactor - i);
}

const nttconfig_t V(ntt7simd_config) = 
{
  7,
  NC,
  ntt7_fixed_const,
  X(ntt7_init),
  ntt7_run_simd,
  ntt7_pfa_run_simd,
  ntt7_twiddle_run_simd,
};

