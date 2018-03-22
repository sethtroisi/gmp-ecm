#include "ntt/ntt-impl-scalar.h"

#define NC 11

static const uint8_t ntt9_fixed_const[NC] = {1};

void
X(ntt9_init)(spv_t out, sp_t p, sp_t d, 
	  sp_t primroot, sp_t order, sp_t perm)
{
  X(nttdata_init_generic)(&X(ntt9_config), out, p, d, primroot, order, perm);
}

static void 
ntt9_run_core(spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		sp_t p, spv_t ntt_const)
{
  sp_t x0, x1, x2, x3, x4, x5, x6;
  sp_t     t1, t2, t3, t4, t5, t6, t7, t8;
  sp_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

  sp_t x0e, x1e, t0e, t1e, p0e, p1e, p2e;

  x0 = in[0 * istride];
  x1 = in[1 * istride];
  x2 = in[2 * istride];
  x0e = in[3 * istride];
  x3 = in[4 * istride];
  x6 = in[5 * istride];
  x1e = in[6 * istride];
  x5 = in[7 * istride];
  x4 = in[8 * istride];

  t0e = sp_ntt_add(x0e, x1e, p);
  t1e = sp_ntt_sub_partial(x0e, x1e, p);

  p1 = sp_ntt_add(x1, x3, p);
  p1 = sp_ntt_add(p1, x5, p);
  p2 = sp_ntt_add(x2, x4, p);
  p2 = sp_ntt_add(p2, x6, p);
  p3 = sp_ntt_sub(x1, x5, p);
  p4 = sp_ntt_sub(x2, x6, p);
  p5 = sp_ntt_sub(x3, x5, p);
  p6 = sp_ntt_sub(x4, x6, p);

  t1 = sp_ntt_add(p1, p2, p);
  t2 = sp_ntt_sub_partial(p1, p2, p);
  t3 = sp_ntt_sub(p3, p5, p);
  t5 = sp_ntt_add(t3, p6, p);
  t3 = sp_ntt_sub(t3, p6, p);
  t4 = sp_ntt_add(p4, p5, p);
  t6 = sp_ntt_sub(p4, p5, p);

  p0e = sp_ntt_add(x0, t0e, p);
  p0 = t1;
  p1 = t1;
  p2 = t2;
  p3 = t3;
  p4 = t4;
  p5 = sp_ntt_add_partial(t3, t4, p);
  p6 = t5;
  p7 = t6;
  p8 = sp_ntt_add_partial(t5, t6, p);

  p1 = sp_ntt_mul(p1, ntt_const[1], ntt_const[NC+1], p);
  p2 = sp_ntt_mul(p2, ntt_const[2], ntt_const[NC+2], p);
  t0e = sp_ntt_mul(t0e, ntt_const[3], ntt_const[NC+3], p);
  t1e = sp_ntt_mul(t1e, ntt_const[4], ntt_const[NC+4], p);
  p3 = sp_ntt_mul(p3, ntt_const[5], ntt_const[NC+5], p);
  p4 = sp_ntt_mul(p4, ntt_const[6], ntt_const[NC+6], p);
  p5 = sp_ntt_mul(p5, ntt_const[7], ntt_const[NC+7], p);
  p6 = sp_ntt_mul(p6, ntt_const[8], ntt_const[NC+8], p);
  p7 = sp_ntt_mul(p7, ntt_const[9], ntt_const[NC+9], p);
  p8 = sp_ntt_mul(p8, ntt_const[10], ntt_const[NC+10], p);

  t0e = sp_ntt_add(t0e, p0e, p);
  t1 = sp_ntt_add(p1, p2, p);
  t2 = sp_ntt_sub(p1, p2, p);
  t3 = sp_ntt_add(p3, p5, p);
  t4 = sp_ntt_add(p4, p5, p);
  t5 = sp_ntt_add(p6, p8, p);
  t6 = sp_ntt_add(p7, p8, p);

  p1e = sp_ntt_add(t0e, t1e, p);
  p2e = sp_ntt_sub(t0e, t1e, p);
  p3 = sp_ntt_add(t3, t5, p);
  p4 = sp_ntt_add(t4, t6, p);
  p5 = sp_ntt_sub(t4, t6, p);
  p5 = sp_ntt_sub(p5, p3, p);
  p6 = sp_ntt_sub(t5, t3, p);

  p0 = sp_ntt_add(p0, p0e, p);
  t1 = sp_ntt_add(t1, p0e, p);
  t2 = sp_ntt_add(t2, p0e, p);
  t3 = sp_ntt_add(p3, p1e, p);
  t4 = sp_ntt_add(p4, p2e, p);
  t5 = sp_ntt_add(p5, p1e, p);
  t6 = sp_ntt_add(p6, p2e, p);
  t7 = sp_ntt_add(p3, p5, p);
  t7 = sp_ntt_sub(p1e, t7, p);
  t8 = sp_ntt_add(p4, p6, p);
  t8 = sp_ntt_sub(p2e, t8, p);

  out[0 * ostride] = p0;
  out[1 * ostride] = t8;
  out[2 * ostride] = t3;
  out[3 * ostride] = t2;
  out[4 * ostride] = t4;
  out[5 * ostride] = t7;
  out[6 * ostride] = t1;
  out[7 * ostride] = t6;
  out[8 * ostride] = t5;
}

static void
ntt9_twiddle_run_core(spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		spv_t w, sp_t p, spv_t ntt_const)
{
  sp_t x0, x1, x2, x3, x4, x5, x6;
  sp_t     t1, t2, t3, t4, t5, t6, t7, t8;
  sp_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

  sp_t x0e, x1e, t0e, t1e, p0e, p1e, p2e;

  x0 = in[0 * istride];
  x1 = in[1 * istride];
  x2 = in[2 * istride];
  x0e = in[3 * istride];
  x3 = in[4 * istride];
  x6 = in[5 * istride];
  x1e = in[6 * istride];
  x5 = in[7 * istride];
  x4 = in[8 * istride];

  t0e = sp_ntt_add(x0e, x1e, p);
  t1e = sp_ntt_sub_partial(x0e, x1e, p);

  p1 = sp_ntt_add(x1, x3, p);
  p1 = sp_ntt_add(p1, x5, p);
  p2 = sp_ntt_add(x2, x4, p);
  p2 = sp_ntt_add(p2, x6, p);
  p3 = sp_ntt_sub(x1, x5, p);
  p4 = sp_ntt_sub(x2, x6, p);
  p5 = sp_ntt_sub(x3, x5, p);
  p6 = sp_ntt_sub(x4, x6, p);

  t1 = sp_ntt_add(p1, p2, p);
  t2 = sp_ntt_sub_partial(p1, p2, p);
  t3 = sp_ntt_sub(p3, p5, p);
  t5 = sp_ntt_add(t3, p6, p);
  t3 = sp_ntt_sub(t3, p6, p);
  t4 = sp_ntt_add(p4, p5, p);
  t6 = sp_ntt_sub(p4, p5, p);

  p0e = sp_ntt_add(x0, t0e, p);
  p0 = t1;
  p1 = t1;
  p2 = t2;
  p3 = t3;
  p4 = t4;
  p5 = sp_ntt_add_partial(t3, t4, p);
  p6 = t5;
  p7 = t6;
  p8 = sp_ntt_add_partial(t5, t6, p);

  p1 = sp_ntt_mul(p1, ntt_const[1], ntt_const[NC+1], p);
  p2 = sp_ntt_mul(p2, ntt_const[2], ntt_const[NC+2], p);
  t0e = sp_ntt_mul(t0e, ntt_const[3], ntt_const[NC+3], p);
  t1e = sp_ntt_mul(t1e, ntt_const[4], ntt_const[NC+4], p);
  p3 = sp_ntt_mul(p3, ntt_const[5], ntt_const[NC+5], p);
  p4 = sp_ntt_mul(p4, ntt_const[6], ntt_const[NC+6], p);
  p5 = sp_ntt_mul(p5, ntt_const[7], ntt_const[NC+7], p);
  p6 = sp_ntt_mul(p6, ntt_const[8], ntt_const[NC+8], p);
  p7 = sp_ntt_mul(p7, ntt_const[9], ntt_const[NC+9], p);
  p8 = sp_ntt_mul(p8, ntt_const[10], ntt_const[NC+10], p);

  t0e = sp_ntt_add(t0e, p0e, p);
  t1 = sp_ntt_add(p1, p2, p);
  t2 = sp_ntt_sub(p1, p2, p);
  t3 = sp_ntt_add(p3, p5, p);
  t4 = sp_ntt_add(p4, p5, p);
  t5 = sp_ntt_add(p6, p8, p);
  t6 = sp_ntt_add(p7, p8, p);

  p1e = sp_ntt_add(t0e, t1e, p);
  p2e = sp_ntt_sub(t0e, t1e, p);
  p3 = sp_ntt_add(t3, t5, p);
  p4 = sp_ntt_add(t4, t6, p);
  p5 = sp_ntt_sub(t4, t6, p);
  p5 = sp_ntt_sub(p5, p3, p);
  p6 = sp_ntt_sub(t5, t3, p);

  p0 = sp_ntt_add(p0, p0e, p);
  t1 = sp_ntt_add_partial(t1, p0e, p);
  t2 = sp_ntt_add_partial(t2, p0e, p);
  t3 = sp_ntt_add_partial(p3, p1e, p);
  t4 = sp_ntt_add_partial(p4, p2e, p);
  t5 = sp_ntt_add_partial(p5, p1e, p);
  t6 = sp_ntt_add_partial(p6, p2e, p);
  t7 = sp_ntt_add(p3, p5, p);
  t7 = sp_ntt_sub_partial(p1e, t7, p);
  t8 = sp_ntt_add(p4, p6, p);
  t8 = sp_ntt_sub_partial(p2e, t8, p);

  t8 = sp_ntt_mul(t8, w[0], w[1], p);
  t3 = sp_ntt_mul(t3, w[2], w[3], p);
  t2 = sp_ntt_mul(t2, w[4], w[5], p);
  t4 = sp_ntt_mul(t4, w[6], w[7], p);
  t7 = sp_ntt_mul(t7, w[8], w[9], p);
  t1 = sp_ntt_mul(t1, w[10], w[11], p);
  t6 = sp_ntt_mul(t6, w[12], w[13], p);
  t5 = sp_ntt_mul(t5, w[14], w[15], p);

  out[0 * ostride] = p0;
  out[1 * ostride] = t8;
  out[2 * ostride] = t3;
  out[3 * ostride] = t2;
  out[4 * ostride] = t4;
  out[5 * ostride] = t7;
  out[6 * ostride] = t1;
  out[7 * ostride] = t6;
  out[8 * ostride] = t5;
}

static void
ntt9_pfa_run_core(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t n,
	  sp_t p, spv_t ntt_const)
{
  spv_size_t j0, j1, j2, j3, j4, j5, j6, j7, j8;
  sp_t x0, x1, x2, x3, x4, x5, x6;
  sp_t     t1, t2, t3, t4, t5, t6, t7, t8;
  sp_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

  sp_t x0e, x1e, t0e, t1e, p0e, p1e, p2e;

  j0 = start;
  j1 = sp_array_inc(j0, inc, n);
  j2 = sp_array_inc(j0, 2 * inc, n);
  j3 = sp_array_inc(j0, 3 * inc, n);
  j4 = sp_array_inc(j0, 4 * inc, n);
  j5 = sp_array_inc(j0, 5 * inc, n);
  j6 = sp_array_inc(j0, 6 * inc, n);
  j7 = sp_array_inc(j0, 7 * inc, n);
  j8 = sp_array_inc(j0, 8 * inc, n);

  x0 = x[j0];
  x1 = x[j1];
  x2 = x[j2];
  x0e = x[j3];
  x3 = x[j4];
  x6 = x[j5];
  x1e = x[j6];
  x5 = x[j7];
  x4 = x[j8];

  t0e = sp_ntt_add(x0e, x1e, p);
  t1e = sp_ntt_sub_partial(x0e, x1e, p);

  p1 = sp_ntt_add(x1, x3, p);
  p1 = sp_ntt_add(p1, x5, p);
  p2 = sp_ntt_add(x2, x4, p);
  p2 = sp_ntt_add(p2, x6, p);
  p3 = sp_ntt_sub(x1, x5, p);
  p4 = sp_ntt_sub(x2, x6, p);
  p5 = sp_ntt_sub(x3, x5, p);
  p6 = sp_ntt_sub(x4, x6, p);

  t1 = sp_ntt_add(p1, p2, p);
  t2 = sp_ntt_sub_partial(p1, p2, p);
  t3 = sp_ntt_sub(p3, p5, p);
  t5 = sp_ntt_add(t3, p6, p);
  t3 = sp_ntt_sub(t3, p6, p);
  t4 = sp_ntt_add(p4, p5, p);
  t6 = sp_ntt_sub(p4, p5, p);

  p0e = sp_ntt_add(x0, t0e, p);
  p0 = t1;
  p1 = t1;
  p2 = t2;
  p3 = t3;
  p4 = t4;
  p5 = sp_ntt_add_partial(t3, t4, p);
  p6 = t5;
  p7 = t6;
  p8 = sp_ntt_add_partial(t5, t6, p);

  p1 = sp_ntt_mul(p1, ntt_const[1], ntt_const[NC+1], p);
  p2 = sp_ntt_mul(p2, ntt_const[2], ntt_const[NC+2], p);
  t0e = sp_ntt_mul(t0e, ntt_const[3], ntt_const[NC+3], p);
  t1e = sp_ntt_mul(t1e, ntt_const[4], ntt_const[NC+4], p);
  p3 = sp_ntt_mul(p3, ntt_const[5], ntt_const[NC+5], p);
  p4 = sp_ntt_mul(p4, ntt_const[6], ntt_const[NC+6], p);
  p5 = sp_ntt_mul(p5, ntt_const[7], ntt_const[NC+7], p);
  p6 = sp_ntt_mul(p6, ntt_const[8], ntt_const[NC+8], p);
  p7 = sp_ntt_mul(p7, ntt_const[9], ntt_const[NC+9], p);
  p8 = sp_ntt_mul(p8, ntt_const[10], ntt_const[NC+10], p);

  t0e = sp_ntt_add(t0e, p0e, p);
  t1 = sp_ntt_add(p1, p2, p);
  t2 = sp_ntt_sub(p1, p2, p);
  t3 = sp_ntt_add(p3, p5, p);
  t4 = sp_ntt_add(p4, p5, p);
  t5 = sp_ntt_add(p6, p8, p);
  t6 = sp_ntt_add(p7, p8, p);

  p1e = sp_ntt_add(t0e, t1e, p);
  p2e = sp_ntt_sub(t0e, t1e, p);
  p3 = sp_ntt_add(t3, t5, p);
  p4 = sp_ntt_add(t4, t6, p);
  p5 = sp_ntt_sub(t4, t6, p);
  p5 = sp_ntt_sub(p5, p3, p);
  p6 = sp_ntt_sub(t5, t3, p);

  p0 = sp_ntt_add(p0, p0e, p);
  t1 = sp_ntt_add(t1, p0e, p);
  t2 = sp_ntt_add(t2, p0e, p);
  t3 = sp_ntt_add(p3, p1e, p);
  t4 = sp_ntt_add(p4, p2e, p);
  t5 = sp_ntt_add(p5, p1e, p);
  t6 = sp_ntt_add(p6, p2e, p);
  t7 = sp_ntt_add(p3, p5, p);
  t7 = sp_ntt_sub(p1e, t7, p);
  t8 = sp_ntt_add(p4, p6, p);
  t8 = sp_ntt_sub(p2e, t8, p);

  x[j0] = p0;
  x[j1] = t8;
  x[j2] = t3;
  x[j3] = t2;
  x[j4] = t4;
  x[j5] = t7;
  x[j6] = t1;
  x[j7] = t6;
  x[j8] = t5;
}

DECLARE_CORE_ROUTINES(9)
