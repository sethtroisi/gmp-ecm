#include "ntt/ntt-impl-scalar.h"

#define NC 3

static const uint8_t ntt3_fixed_const[NC] = {1};

void
X(ntt3_init)(spv_t out, sp_t p, sp_t d,
	  sp_t primroot, sp_t order, sp_t perm)
{
  X(nttdata_init_generic)(&X(ntt3_config), out, p, d, primroot, order, perm);
}

static void 
ntt3_run_core(spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		sp_t p, spv_t ntt_const)
{
  sp_t p0, p1, p2;
  sp_t x0, x1, x2;
  sp_t     t1, t2;

  x0 = in[0 * istride];
  x1 = in[1 * istride];
  x2 = in[2 * istride];

  t1 = sp_ntt_add(x1, x2, p);
  t2 = sp_ntt_sub_partial(x1, x2, p);

  p0 = sp_ntt_add(x0, t1, p);

  p1 = sp_ntt_mul(t1, ntt_const[1], ntt_const[NC+1], p);
  p2 = sp_ntt_mul(t2, ntt_const[2], ntt_const[NC+2], p);

  p1 = sp_ntt_add(p0, p1, p);

  t1 = sp_ntt_add(p1, p2, p);
  t2 = sp_ntt_sub(p1, p2, p);

  out[0 * ostride] = p0;
  out[1 * ostride] = t1;
  out[2 * ostride] = t2;
}

static void
ntt3_run(spv_t in, spv_size_t istride, spv_size_t idist,
    		spv_t out, spv_size_t ostride, spv_size_t odist,
    		spv_size_t num_transforms, sp_t p, spv_t ntt_const)
{
  spv_size_t i = 0;

  for (; i < num_transforms; i++)
    ntt3_run_core(in + i * idist, istride, 
                out + i * odist, ostride, p, ntt_const);
}


static void
ntt3_twiddle_run_core(spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		spv_t w, sp_t p, spv_t ntt_const)
{
  sp_t p0, p1, p2;
  sp_t x0, x1, x2;
  sp_t     t1, t2;

  x0 = in[0 * istride];
  x1 = in[1 * istride];
  x2 = in[2 * istride];

  t1 = sp_ntt_add(x1, x2, p);
  t2 = sp_ntt_sub_partial(x1, x2, p);

  p0 = sp_ntt_add(x0, t1, p);

  p1 = sp_ntt_mul(t1, ntt_const[1], ntt_const[NC+1], p);
  p2 = sp_ntt_mul(t2, ntt_const[2], ntt_const[NC+2], p);

  p1 = sp_ntt_add(p0, p1, p);

  t1 = sp_ntt_add_partial(p1, p2, p);
  t2 = sp_ntt_sub_partial(p1, p2, p);

  t1 = sp_ntt_mul(t1, w[0], w[1], p);
  t2 = sp_ntt_mul(t2, w[2], w[3], p);

  out[0 * ostride] = p0;
  out[1 * ostride] = t1;
  out[2 * ostride] = t2;
}


static void
ntt3_twiddle_run(spv_t in, spv_size_t istride, spv_size_t idist,
    			spv_t out, spv_size_t ostride, spv_size_t odist,
    			spv_t w, spv_size_t num_transforms, sp_t p, spv_t ntt_const)
{
  spv_size_t i = 0, j = 0;

  for (; i < num_transforms; i++, j += 2*(3-1))
    ntt3_twiddle_run_core(in + i * idist, istride, 
			out + i * odist, ostride,
			w + j, p, ntt_const);
}


static void
ntt3_pfa_run_core(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t n,
	  sp_t p, spv_t ntt_const)
{
  spv_size_t j0, j1, j2;
  sp_t p0, p1, p2;
  sp_t x0, x1, x2;
  sp_t     t1, t2;

  j0 = start;
  j1 = sp_array_inc(j0, inc, n);
  j2 = sp_array_inc(j0, 2 * inc, n);

  x0 = x[j0];
  x1 = x[j1];
  x2 = x[j2];

  t1 = sp_ntt_add(x1, x2, p);
  t2 = sp_ntt_sub_partial(x1, x2, p);

  p0 = sp_ntt_add(x0, t1, p);

  p1 = sp_ntt_mul(t1, ntt_const[1], ntt_const[NC+1], p);
  p2 = sp_ntt_mul(t2, ntt_const[2], ntt_const[NC+2], p);

  p1 = sp_ntt_add(p0, p1, p);

  t1 = sp_ntt_add(p1, p2, p);
  t2 = sp_ntt_sub(p1, p2, p);

  x[j0] = p0;
  x[j1] = t1;
  x[j2] = t2;
}

static void
ntt3_pfa_run(spv_t x, spv_size_t cofactor,
	  sp_t p, spv_t ntt_const)
{
  spv_size_t i = 0;
  spv_size_t incstart = 0;
  spv_size_t n = 3 * cofactor;
  spv_size_t inc = cofactor;
  spv_size_t inc2 = 3;

  for (; i < cofactor; i++, incstart += inc2)
    ntt3_pfa_run_core(x, incstart, inc, n, p, ntt_const);

}


const nttconfig_t X(ntt3_config) = 
{
  3,
  NC,
  ntt3_fixed_const,
  X(ntt3_init),
  ntt3_run,
  ntt3_pfa_run,
  ntt3_twiddle_run
};

