#include "ntt-impl.h"

#define NC 2

static const uint8_t fixed_const[NC] = {1};

static const uint8_t *
ntt2_get_fixed_ntt_const(void)
{
  return fixed_const;
}

void
ntt2_init(spv_t out, sp_t p, sp_t d, 
	  sp_t primroot, sp_t order, sp_t perm)
{
  out[0] = out[1] = 1;
}

static void 
ntt2_run_core(spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		sp_t p, spv_t ntt_const)
{
  sp_t p0, p1;
  sp_t x0, x1;

  x0 = in[0 * istride];
  x1 = in[1 * istride];

  p0 = sp_ntt_add(x0, x1, p);
  p1 = sp_ntt_sub(x0, x1, p);

  out[0 * ostride] = p0;
  out[1 * ostride] = p1;
}

#ifdef HAVE_SIMD
static void
ntt2_run_core_simd(spv_t in, spv_size_t istride, spv_size_t idist,
		spv_t out, spv_size_t ostride, spv_size_t odist,
		sp_t p, spv_t ntt_const)
{
  sp_simd_t p0, p1;
  sp_simd_t x0, x1;

  x0 = sp_simd_gather(in + 0 * istride, idist);
  x1 = sp_simd_gather(in + 1 * istride, idist);

  p0 = sp_ntt_add_simd(x0, x1, p);
  p1 = sp_ntt_sub_simd(x0, x1, p);

  sp_simd_scatter(p0, out + 0 * ostride, odist);
  sp_simd_scatter(p1, out + 1 * ostride, odist);
}
#endif


static void
ntt2_run(spv_t in, spv_size_t istride, spv_size_t idist,
    		spv_t out, spv_size_t ostride, spv_size_t odist,
    		spv_size_t num_transforms, sp_t p, spv_t ntt_const)
{
  spv_size_t i = 0;

#ifdef HAVE_SIMD
  spv_size_t num_simd = SP_SIMD_VSIZE * (num_transforms / SP_SIMD_VSIZE);

  for (; i < num_simd; i += SP_SIMD_VSIZE)
    ntt2_run_core_simd(in + i * idist, istride, idist, 
                        out + i * odist, ostride, odist, p, ntt_const);
#endif

  for (; i < num_transforms; i++)
    ntt2_run_core(in + i * idist, istride, 
                out + i * odist, ostride, p, ntt_const);
}


static void
ntt2_twiddle_run_core(spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		spv_t w, sp_t p, spv_t ntt_const)
{
  sp_t p0, p1;
  sp_t x0, x1;

  x0 = in[0 * istride];
  x1 = in[1 * istride];

  p0 = sp_ntt_add(x0, x1, p);
  p1 = sp_ntt_sub_partial(x0, x1, p);

  p1 = sp_ntt_mul(p1, w[0], w[1], p);

  out[0 * ostride] = p0;
  out[1 * ostride] = p1;
}

#ifdef HAVE_SIMD
static void
ntt2_twiddle_run_core_simd(
        spv_t in, spv_size_t istride, spv_size_t idist,
		spv_t out, spv_size_t ostride, spv_size_t odist,
		sp_simd_t *w, sp_t p, spv_t ntt_const)
{
  sp_simd_t p0, p1;
  sp_simd_t x0, x1;

  x0 = sp_simd_gather(in + 0 * istride, idist);
  x1 = sp_simd_gather(in + 1 * istride, idist);

  p0 = sp_ntt_add_simd(x0, x1, p);
  p1 = sp_ntt_sub_partial_simd(x0, x1, p);

  p1 = sp_ntt_twiddle_mul_simd(p1, w + 0, p);

  sp_simd_scatter(p0, out + 0 * ostride, odist);
  sp_simd_scatter(p1, out + 1 * ostride, odist);
}
#endif

static void
ntt2_twiddle_run(spv_t in, spv_size_t istride, spv_size_t idist,
    			spv_t out, spv_size_t ostride, spv_size_t odist,
    			spv_t w, spv_size_t num_transforms, sp_t p, spv_t ntt_const)
{
  spv_size_t i = 0, j = 0;

#ifdef HAVE_SIMD
  spv_size_t num_simd = SP_SIMD_VSIZE * (num_transforms / SP_SIMD_VSIZE);

  for (; i < num_simd; i += SP_SIMD_VSIZE,
		  	j += 2*(2-1)*SP_SIMD_VSIZE)
    ntt2_twiddle_run_core_simd(
		in + i * idist, istride, idist,
		out + i * odist, ostride, odist,
		(sp_simd_t *)(w + j), p, ntt_const);
#endif

  for (; i < num_transforms; i++, j += 2*(2-1))
    ntt2_twiddle_run_core(in + i * idist, istride, 
			out + i * odist, ostride,
			w + j, p, ntt_const);
}


static void
ntt2_pfa_run_core(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t n,
	  sp_t p, spv_t ntt_const)
{
  spv_size_t j0, j1;
  sp_t p0, p1;
  sp_t x0, x1;

  j0 = start;
  j1 = sp_array_inc(j0, inc, n);

  x0 = x[j0];
  x1 = x[j1];

  p0 = sp_ntt_add(x0, x1, p);
  p1 = sp_ntt_sub(x0, x1, p);

  x[j0] = p0;
  x[j1] = p1;
}

#ifdef HAVE_SIMD
static void
ntt2_pfa_run_core_simd(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t inc2, spv_size_t n,
	  sp_t p, spv_t ntt_const)
{
  spv_size_t j0, j1;
  sp_simd_t p0, p1;
  sp_simd_t x0, x1;

  j0 = start;
  j1 = sp_array_inc(j0, inc, n);

  x0 = sp_simd_pfa_gather(x, j0, inc2, n);
  x1 = sp_simd_pfa_gather(x, j1, inc2, n);

  p0 = sp_ntt_add_simd(x0, x1, p);
  p1 = sp_ntt_sub_simd(x0, x1, p);

  sp_simd_pfa_scatter(p0, x, j0, inc2, n);
  sp_simd_pfa_scatter(p1, x, j1, inc2, n);
}
#endif

static void
ntt2_pfa_run(spv_t x, spv_size_t cofactor,
	  sp_t p, spv_t ntt_const)
{
  spv_size_t i = 0;
  spv_size_t incstart = 0;
  spv_size_t n = 2 * cofactor;
  spv_size_t inc = cofactor;
  spv_size_t inc2 = 2;

#ifdef HAVE_SIMD
  spv_size_t num_simd = SP_SIMD_VSIZE * (cofactor / SP_SIMD_VSIZE);

  for (i = 0; i < num_simd; i += SP_SIMD_VSIZE)
    {
      ntt2_pfa_run_core_simd(x, incstart, inc, inc2, n, p, ntt_const);
      incstart += SP_SIMD_VSIZE * inc2;
    }
#endif

  for (; i < cofactor; i++, incstart += inc2)
    ntt2_pfa_run_core(x, incstart, inc, n, p, ntt_const);

}

const nttconfig_t ntt2_config = 
{
  2,
  NC,
  ntt2_get_fixed_ntt_const,
  ntt2_init,
  ntt2_run,
  ntt2_pfa_run,
  ntt2_twiddle_run
};

