#include "ntt/ntt-impl-scalar.h"

#define NC 4

static const uint8_t ntt4_fixed_const[NC] = {1, 1, 1, 0};

void
X(ntt4_init)(spv_t out, sp_t p, sp_t d, 
	  sp_t primroot, sp_t order, sp_t perm)
{
  out[0] = 1;
  out[1] = 1;
  out[2] = 1;
  out[3] = sp_pow(primroot, order / 4, p, d);
}

static void 
ntt4_run_core(spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		sp_t p, spv_t ntt_const)
{
  sp_t x0, x1, x2, x3;
  sp_t t0, t1, t2, t3;
  sp_t p0, p1, p2, p3;

  x0 = in[0 * istride];
  x1 = in[1 * istride];
  x2 = in[2 * istride];
  x3 = in[3 * istride];

  t0 = sp_ntt_add(x0, x2, p);
  t2 = sp_ntt_sub(x0, x2, p);
  t1 = sp_ntt_add(x1, x3, p);
  t3 = sp_ntt_sub_partial(x1, x3, p);

  t3 = sp_ntt_mul(t3, ntt_const[3], ntt_const[NC+3], p);

  p0 = sp_ntt_add(t0, t1, p);
  p1 = sp_ntt_sub(t0, t1, p);
  p2 = sp_ntt_add(t2, t3, p);
  p3 = sp_ntt_sub(t2, t3, p);

  out[0 * ostride] = p0;
  out[1 * ostride] = p2;
  out[2 * ostride] = p1;
  out[3 * ostride] = p3;
}


static void
ntt4_run(spv_t in, spv_size_t istride, spv_size_t idist,
    		spv_t out, spv_size_t ostride, spv_size_t odist,
    		spv_size_t num_transforms, sp_t p, spv_t ntt_const)
{
  spv_size_t i = 0;

  for (; i < num_transforms; i++)
    ntt4_run_core(in + i * idist, istride, 
                out + i * odist, ostride, p, ntt_const);
}


static void
ntt4_twiddle_run_core(spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		spv_t w, sp_t p, spv_t ntt_const)
{
  sp_t x0, x1, x2, x3;
  sp_t t0, t1, t2, t3;
  sp_t p0, p1, p2, p3;

  x0 = in[0 * istride];
  x1 = in[1 * istride];
  x2 = in[2 * istride];
  x3 = in[3 * istride];

  t0 = sp_ntt_add(x0, x2, p);
  t2 = sp_ntt_sub(x0, x2, p);
  t1 = sp_ntt_add(x1, x3, p);
  t3 = sp_ntt_sub_partial(x1, x3, p);

  t3 = sp_ntt_mul(t3, ntt_const[3], ntt_const[NC+3], p);

  p0 = sp_ntt_add(t0, t1, p);
  p1 = sp_ntt_sub_partial(t0, t1, p);
  p2 = sp_ntt_add_partial(t2, t3, p);
  p3 = sp_ntt_sub_partial(t2, t3, p);

  p2 = sp_ntt_mul(p2, w[0], w[1], p);
  p1 = sp_ntt_mul(p1, w[2], w[3], p);
  p3 = sp_ntt_mul(p3, w[4], w[5], p);

  out[0 * ostride] = p0;
  out[1 * ostride] = p2;
  out[2 * ostride] = p1;
  out[3 * ostride] = p3;
}

static void
ntt4_twiddle_run(spv_t in, spv_size_t istride, spv_size_t idist,
    			spv_t out, spv_size_t ostride, spv_size_t odist,
    			spv_t w, spv_size_t num_transforms, sp_t p, spv_t ntt_const)
{
  spv_size_t i = 0, j = 0;

  for (; i < num_transforms; i++, j += 2*(4-1))
    ntt4_twiddle_run_core(in + i * idist, istride, 
			out + i * odist, ostride,
			w + j, p, ntt_const);
}

static void
ntt4_pfa_run_core(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t n,
	  sp_t p, spv_t ntt_const)
{
  spv_size_t j0, j1, j2, j3;

  sp_t x0, x1, x2, x3;
  sp_t t0, t1, t2, t3;
  sp_t p0, p1, p2, p3;

  j0 = start;
  j1 = sp_array_inc(j0, inc, n);
  j2 = sp_array_inc(j0, 2 * inc, n);
  j3 = sp_array_inc(j0, 3 * inc, n);

  x0 = x[j0];
  x1 = x[j1];
  x2 = x[j2];
  x3 = x[j3];

  t0 = sp_ntt_add(x0, x2, p);
  t2 = sp_ntt_sub(x0, x2, p);
  t1 = sp_ntt_add(x1, x3, p);
  t3 = sp_ntt_sub_partial(x1, x3, p);

  t3 = sp_ntt_mul(t3, ntt_const[3], ntt_const[NC+3], p);

  p0 = sp_ntt_add(t0, t1, p);
  p1 = sp_ntt_sub(t0, t1, p);
  p2 = sp_ntt_add(t2, t3, p);
  p3 = sp_ntt_sub(t2, t3, p);

  x[j0] = p0;
  x[j1] = p2;
  x[j2] = p1;
  x[j3] = p3;
}

static void
ntt4_pfa_run(spv_t x, spv_size_t cofactor,
	  sp_t p, spv_t ntt_const)
{
  spv_size_t i = 0;
  spv_size_t incstart = 0;
  spv_size_t n = 4 * cofactor;
  spv_size_t inc = cofactor;
  spv_size_t inc2 = 4;

  for (; i < cofactor; i++, incstart += inc2)
    ntt4_pfa_run_core(x, incstart, inc, n, p, ntt_const);
}


const nttconfig_t X(ntt4_config) = 
{
  4,
  NC,
  ntt4_fixed_const,
  X(ntt4_init),
  ntt4_run,
  ntt4_pfa_run,
  ntt4_twiddle_run
};

