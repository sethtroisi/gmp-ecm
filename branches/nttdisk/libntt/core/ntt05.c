#include "ntt/ntt-impl-scalar.h"

#define NC 6

static const uint8_t ntt5_fixed_const[NC] = {1};

void
X(ntt5_init)(spv_t out, sp_t p, sp_t d,
	  sp_t primroot, sp_t order, sp_t perm)
{
  X(nttdata_init_generic)(&X(ntt5_config), out, p, d, primroot, order, perm);
}

static void 
ntt5_run_core(spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		sp_t p, spv_t ntt_const)
{
  sp_t p0, p1, p2, p3, p4, p5;
  sp_t x0, x1, x2, x3, x4;
  sp_t     t1, t2, t3, t4;

  x0 = in[0 * istride];
  x1 = in[1 * istride];
  x4 = in[2 * istride];
  x2 = in[3 * istride];
  x3 = in[4 * istride];

  t1 = sp_ntt_add(x1, x3, p);
  t3 = sp_ntt_sub(x1, x3, p);
  t2 = sp_ntt_add(x2, x4, p);
  t4 = sp_ntt_sub(x2, x4, p);

  p1 = sp_ntt_add(t1, t2, p);
  p2 = sp_ntt_sub_partial(t1, t2, p);
  p3 = t3;
  p4 = t4;
  p5 = sp_ntt_add_partial(t3, t4, p);

  p0 = sp_ntt_add(x0, p1, p);

  p1 = sp_ntt_mul(p1, ntt_const[1], ntt_const[NC+1], p);
  p2 = sp_ntt_mul(p2, ntt_const[2], ntt_const[NC+2], p);
  p3 = sp_ntt_mul(p3, ntt_const[3], ntt_const[NC+3], p);
  p4 = sp_ntt_mul(p4, ntt_const[4], ntt_const[NC+4], p);
  p5 = sp_ntt_mul(p5, ntt_const[5], ntt_const[NC+5], p);

  p1 = sp_ntt_add(p0, p1, p);

  t1 = sp_ntt_add(p1, p2, p);
  t2 = sp_ntt_sub(p1, p2, p);
  t3 = sp_ntt_add(p3, p5, p);
  t4 = sp_ntt_add(p4, p5, p);

  p1 = sp_ntt_add(t1, t3, p);
  p2 = sp_ntt_add(t2, t4, p);
  p3 = sp_ntt_sub(t1, t3, p);
  p4 = sp_ntt_sub(t2, t4, p);

  out[0 * ostride] = p0;
  out[1 * ostride] = p4;
  out[2 * ostride] = p3;
  out[3 * ostride] = p1;
  out[4 * ostride] = p2;
}

static void
ntt5_run(spv_t in, spv_size_t istride, spv_size_t idist,
    		spv_t out, spv_size_t ostride, spv_size_t odist,
    		spv_size_t num_transforms, sp_t p, spv_t ntt_const)
{
  spv_size_t i = 0;

  for (; i < num_transforms; i++)
    ntt5_run_core(in + i * idist, istride, 
                out + i * odist, ostride, p, ntt_const);
}


static void
ntt5_twiddle_run_core(spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		spv_t w, sp_t p, spv_t ntt_const)
{
  sp_t p0, p1, p2, p3, p4, p5;
  sp_t x0, x1, x2, x3, x4;
  sp_t     t1, t2, t3, t4;

  x0 = in[0 * istride];
  x1 = in[1 * istride];
  x4 = in[2 * istride];
  x2 = in[3 * istride];
  x3 = in[4 * istride];

  t1 = sp_ntt_add(x1, x3, p);
  t3 = sp_ntt_sub(x1, x3, p);
  t2 = sp_ntt_add(x2, x4, p);
  t4 = sp_ntt_sub(x2, x4, p);

  p1 = sp_ntt_add(t1, t2, p);
  p2 = sp_ntt_sub_partial(t1, t2, p);
  p3 = t3;
  p4 = t4;
  p5 = sp_ntt_add_partial(t3, t4, p);

  p0 = sp_ntt_add(x0, p1, p);

  p1 = sp_ntt_mul(p1, ntt_const[1], ntt_const[NC+1], p);
  p2 = sp_ntt_mul(p2, ntt_const[2], ntt_const[NC+2], p);
  p3 = sp_ntt_mul(p3, ntt_const[3], ntt_const[NC+3], p);
  p4 = sp_ntt_mul(p4, ntt_const[4], ntt_const[NC+4], p);
  p5 = sp_ntt_mul(p5, ntt_const[5], ntt_const[NC+5], p);

  p1 = sp_ntt_add(p0, p1, p);

  t1 = sp_ntt_add(p1, p2, p);
  t2 = sp_ntt_sub(p1, p2, p);
  t3 = sp_ntt_add(p3, p5, p);
  t4 = sp_ntt_add(p4, p5, p);

  p1 = sp_ntt_add_partial(t1, t3, p);
  p2 = sp_ntt_add_partial(t2, t4, p);
  p3 = sp_ntt_sub_partial(t1, t3, p);
  p4 = sp_ntt_sub_partial(t2, t4, p);

  p4 = sp_ntt_mul(p4, w[0], w[1], p);
  p3 = sp_ntt_mul(p3, w[2], w[3], p);
  p1 = sp_ntt_mul(p1, w[4], w[5], p);
  p2 = sp_ntt_mul(p2, w[6], w[7], p);

  out[0 * ostride] = p0;
  out[1 * ostride] = p4;
  out[2 * ostride] = p3;
  out[3 * ostride] = p1;
  out[4 * ostride] = p2;
}


static void
ntt5_twiddle_run(spv_t in, spv_size_t istride, spv_size_t idist,
    			spv_t out, spv_size_t ostride, spv_size_t odist,
    			spv_t w, spv_size_t num_transforms, sp_t p, spv_t ntt_const)
{
  spv_size_t i = 0, j = 0;

  for (; i < num_transforms; i++, j += 2*(5-1))
    ntt5_twiddle_run_core(in + i * idist, istride, 
			out + i * odist, ostride,
			w + j, p, ntt_const);
}


static void
ntt5_pfa_run_core(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t n,
	  sp_t p, spv_t ntt_const)
{
  spv_size_t j0, j1, j2, j3, j4;
  sp_t p0, p1, p2, p3, p4, p5;
  sp_t x0, x1, x2, x3, x4;
  sp_t     t1, t2, t3, t4;

  j0 = start;
  j1 = sp_array_inc(j0, inc, n);
  j2 = sp_array_inc(j0, 2 * inc, n);
  j3 = sp_array_inc(j0, 3 * inc, n);
  j4 = sp_array_inc(j0, 4 * inc, n);

  x0 = x[j0];
  x1 = x[j1];
  x4 = x[j2];
  x2 = x[j3];
  x3 = x[j4];

  t1 = sp_ntt_add(x1, x3, p);
  t3 = sp_ntt_sub(x1, x3, p);
  t2 = sp_ntt_add(x2, x4, p);
  t4 = sp_ntt_sub(x2, x4, p);

  p1 = sp_ntt_add(t1, t2, p);
  p2 = sp_ntt_sub_partial(t1, t2, p);
  p3 = t3;
  p4 = t4;
  p5 = sp_ntt_add_partial(t3, t4, p);

  p0 = sp_ntt_add(x0, p1, p);

  p1 = sp_ntt_mul(p1, ntt_const[1], ntt_const[NC+1], p);
  p2 = sp_ntt_mul(p2, ntt_const[2], ntt_const[NC+2], p);
  p3 = sp_ntt_mul(p3, ntt_const[3], ntt_const[NC+3], p);
  p4 = sp_ntt_mul(p4, ntt_const[4], ntt_const[NC+4], p);
  p5 = sp_ntt_mul(p5, ntt_const[5], ntt_const[NC+5], p);

  p1 = sp_ntt_add(p0, p1, p);

  t1 = sp_ntt_add(p1, p2, p);
  t2 = sp_ntt_sub(p1, p2, p);
  t3 = sp_ntt_add(p3, p5, p);
  t4 = sp_ntt_add(p4, p5, p);

  p1 = sp_ntt_add(t1, t3, p);
  p2 = sp_ntt_add(t2, t4, p);
  p3 = sp_ntt_sub(t1, t3, p);
  p4 = sp_ntt_sub(t2, t4, p);

  x[j0] = p0;
  x[j1] = p4;
  x[j2] = p3;
  x[j3] = p1;
  x[j4] = p2;
}

static void
ntt5_pfa_run(spv_t x, spv_size_t cofactor,
	  sp_t p, spv_t ntt_const)
{
  spv_size_t i = 0;
  spv_size_t incstart = 0;
  spv_size_t n = 5 * cofactor;
  spv_size_t inc = cofactor;
  spv_size_t inc2 = 5;

  for (; i < cofactor; i++, incstart += inc2)
    ntt5_pfa_run_core(x, incstart, inc, n, p, ntt_const);
}


const nttconfig_t X(ntt5_config) = 
{
  5,
  NC,
  ntt5_fixed_const,
  X(ntt5_init),
  ntt5_run,
  ntt5_pfa_run,
  ntt5_twiddle_run
};

