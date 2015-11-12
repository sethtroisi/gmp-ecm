#include "ntt/ntt-impl-scalar.h"

#define NC 9

static const uint8_t ntt7_fixed_const[NC] = {1};

void
X(ntt7_init)(spv_t out, sp_t p, sp_t d, 
	  sp_t primroot, sp_t order, sp_t perm)
{
  X(nttdata_init_generic)(&X(ntt7_config), out, p, d, primroot, order, perm);
}

static void 
ntt7_run_core(spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		sp_t p, spv_t ntt_const)
{
  sp_t x0, x1, x2, x3, x4, x5, x6;
  sp_t     t1, t2, t3, t4, t5, t6, t7, t8;
  sp_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

  x0 = in[0 * istride];
  x1 = in[1 * istride];
  x5 = in[2 * istride];
  x6 = in[3 * istride];
  x3 = in[4 * istride];
  x2 = in[5 * istride];
  x4 = in[6 * istride];

  p1 = sp_ntt_add(x1, x4, p);
  p2 = sp_ntt_add(x2, x5, p);
  p3 = sp_ntt_add(x3, x6, p);
  p4 = sp_ntt_sub(x1, x4, p);
  p5 = sp_ntt_sub(x2, x5, p);
  p6 = sp_ntt_sub(x3, x6, p);

  t1 = sp_ntt_add(p1, p2, p);
  t1 = sp_ntt_add(t1, p3, p);
  t2 = sp_ntt_sub(p1, p3, p);
  t3 = sp_ntt_sub(p2, p3, p);
  t4 = sp_ntt_sub(p4, p5, p);
  t4 = sp_ntt_add_partial(t4, p6, p);
  t5 = sp_ntt_sub(p4, p6, p);
  t6 = sp_ntt_add(p5, p6, p);

  p0 = sp_ntt_add(x0, t1, p);
  p1 = t1;
  p2 = t2;
  p3 = t3;
  p4 = sp_ntt_add_partial(t2, t3, p);
  p5 = t4;
  p6 = t5;
  p7 = t6;
  p8 = sp_ntt_add_partial(t5, t6, p);

  p1 = sp_ntt_mul(p1, ntt_const[1], ntt_const[NC+1], p);
  p2 = sp_ntt_mul(p2, ntt_const[2], ntt_const[NC+2], p);
  p3 = sp_ntt_mul(p3, ntt_const[3], ntt_const[NC+3], p);
  p4 = sp_ntt_mul(p4, ntt_const[4], ntt_const[NC+4], p);
  p5 = sp_ntt_mul(p5, ntt_const[5], ntt_const[NC+5], p);
  p6 = sp_ntt_mul(p6, ntt_const[6], ntt_const[NC+6], p);
  p7 = sp_ntt_mul(p7, ntt_const[7], ntt_const[NC+7], p);
  p8 = sp_ntt_mul(p8, ntt_const[8], ntt_const[NC+8], p);

  t1 = sp_ntt_add(p0, p1, p);
  t2 = sp_ntt_add(p2, p4, p);
  t3 = sp_ntt_add(p3, p4, p);
  t4 = p5;
  t5 = sp_ntt_add(p6, p8, p);
  t6 = sp_ntt_add(p7, p8, p);

  p1 = sp_ntt_add(t1, t2, p);
  p2 = sp_ntt_add(t1, t3, p);
  p3 = sp_ntt_sub(t1, t2, p);
  p3 = sp_ntt_sub(p3, t3, p);
  p4 = sp_ntt_add(t4, t5, p);
  p5 = sp_ntt_sub(t6, t4, p);
  p6 = sp_ntt_sub(t4, t5, p);
  p6 = sp_ntt_add(p6, t6, p);

  t1 = sp_ntt_add(p1, p4, p);
  t2 = sp_ntt_add(p2, p5, p);
  t3 = sp_ntt_add(p3, p6, p);
  t4 = sp_ntt_sub(p1, p4, p);
  t5 = sp_ntt_sub(p2, p5, p);
  t6 = sp_ntt_sub(p3, p6, p);

  out[0 * ostride] = p0;
  out[1 * ostride] = t6;
  out[2 * ostride] = t4;
  out[3 * ostride] = t5;
  out[4 * ostride] = t2;
  out[5 * ostride] = t1;
  out[6 * ostride] = t3;
}

static void
ntt7_twiddle_run_core(spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		spv_t w, sp_t p, spv_t ntt_const)
{
  sp_t x0, x1, x2, x3, x4, x5, x6;
  sp_t     t1, t2, t3, t4, t5, t6, t7, t8;
  sp_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

  x0 = in[0 * istride];
  x1 = in[1 * istride];
  x5 = in[2 * istride];
  x6 = in[3 * istride];
  x3 = in[4 * istride];
  x2 = in[5 * istride];
  x4 = in[6 * istride];

  p1 = sp_ntt_add(x1, x4, p);
  p2 = sp_ntt_add(x2, x5, p);
  p3 = sp_ntt_add(x3, x6, p);
  p4 = sp_ntt_sub(x1, x4, p);
  p5 = sp_ntt_sub(x2, x5, p);
  p6 = sp_ntt_sub(x3, x6, p);

  t1 = sp_ntt_add(p1, p2, p);
  t1 = sp_ntt_add(t1, p3, p);
  t2 = sp_ntt_sub(p1, p3, p);
  t3 = sp_ntt_sub(p2, p3, p);
  t4 = sp_ntt_sub(p4, p5, p);
  t4 = sp_ntt_add_partial(t4, p6, p);
  t5 = sp_ntt_sub(p4, p6, p);
  t6 = sp_ntt_add(p5, p6, p);

  p0 = sp_ntt_add(x0, t1, p);
  p1 = t1;
  p2 = t2;
  p3 = t3;
  p4 = sp_ntt_add_partial(t2, t3, p);
  p5 = t4;
  p6 = t5;
  p7 = t6;
  p8 = sp_ntt_add_partial(t5, t6, p);

  p1 = sp_ntt_mul(p1, ntt_const[1], ntt_const[NC+1], p);
  p2 = sp_ntt_mul(p2, ntt_const[2], ntt_const[NC+2], p);
  p3 = sp_ntt_mul(p3, ntt_const[3], ntt_const[NC+3], p);
  p4 = sp_ntt_mul(p4, ntt_const[4], ntt_const[NC+4], p);
  p5 = sp_ntt_mul(p5, ntt_const[5], ntt_const[NC+5], p);
  p6 = sp_ntt_mul(p6, ntt_const[6], ntt_const[NC+6], p);
  p7 = sp_ntt_mul(p7, ntt_const[7], ntt_const[NC+7], p);
  p8 = sp_ntt_mul(p8, ntt_const[8], ntt_const[NC+8], p);

  t1 = sp_ntt_add(p0, p1, p);
  t2 = sp_ntt_add(p2, p4, p);
  t3 = sp_ntt_add(p3, p4, p);
  t4 = p5;
  t5 = sp_ntt_add(p6, p8, p);
  t6 = sp_ntt_add(p7, p8, p);

  p1 = sp_ntt_add(t1, t2, p);
  p2 = sp_ntt_add(t1, t3, p);
  p3 = sp_ntt_sub(t1, t2, p);
  p3 = sp_ntt_sub(p3, t3, p);
  p4 = sp_ntt_add(t4, t5, p);
  p5 = sp_ntt_sub(t6, t4, p);
  p6 = sp_ntt_sub(t4, t5, p);
  p6 = sp_ntt_add(p6, t6, p);

  t1 = sp_ntt_add_partial(p1, p4, p);
  t2 = sp_ntt_add_partial(p2, p5, p);
  t3 = sp_ntt_add_partial(p3, p6, p);
  t4 = sp_ntt_sub_partial(p1, p4, p);
  t5 = sp_ntt_sub_partial(p2, p5, p);
  t6 = sp_ntt_sub_partial(p3, p6, p);

  t6 = sp_ntt_mul(t6, w[0], w[1], p);
  t4 = sp_ntt_mul(t4, w[2], w[3], p);
  t5 = sp_ntt_mul(t5, w[4], w[5], p);
  t2 = sp_ntt_mul(t2, w[6], w[7], p);
  t1 = sp_ntt_mul(t1, w[8], w[9], p);
  t3 = sp_ntt_mul(t3, w[10], w[11], p);

  out[0 * ostride] = p0;
  out[1 * ostride] = t6;
  out[2 * ostride] = t4;
  out[3 * ostride] = t5;
  out[4 * ostride] = t2;
  out[5 * ostride] = t1;
  out[6 * ostride] = t3;
}

static void
ntt7_pfa_run_core(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t n,
	  sp_t p, spv_t ntt_const)
{
  spv_size_t j0, j1, j2, j3, j4, j5, j6;
  sp_t x0, x1, x2, x3, x4, x5, x6;
  sp_t     t1, t2, t3, t4, t5, t6, t7, t8;
  sp_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

  j0 = start;
  j1 = sp_array_inc(j0, inc, n);
  j2 = sp_array_inc(j0, 2 * inc, n);
  j3 = sp_array_inc(j0, 3 * inc, n);
  j4 = sp_array_inc(j0, 4 * inc, n);
  j5 = sp_array_inc(j0, 5 * inc, n);
  j6 = sp_array_inc(j0, 6 * inc, n);

  x0 = x[j0];
  x1 = x[j1];
  x5 = x[j2];
  x6 = x[j3];
  x3 = x[j4];
  x2 = x[j5];
  x4 = x[j6];

  p1 = sp_ntt_add(x1, x4, p);
  p2 = sp_ntt_add(x2, x5, p);
  p3 = sp_ntt_add(x3, x6, p);
  p4 = sp_ntt_sub(x1, x4, p);
  p5 = sp_ntt_sub(x2, x5, p);
  p6 = sp_ntt_sub(x3, x6, p);

  t1 = sp_ntt_add(p1, p2, p);
  t1 = sp_ntt_add(t1, p3, p);
  t2 = sp_ntt_sub(p1, p3, p);
  t3 = sp_ntt_sub(p2, p3, p);
  t4 = sp_ntt_sub(p4, p5, p);
  t4 = sp_ntt_add_partial(t4, p6, p);
  t5 = sp_ntt_sub(p4, p6, p);
  t6 = sp_ntt_add(p5, p6, p);

  p0 = sp_ntt_add(x0, t1, p);
  p1 = t1;
  p2 = t2;
  p3 = t3;
  p4 = sp_ntt_add_partial(t2, t3, p);
  p5 = t4;
  p6 = t5;
  p7 = t6;
  p8 = sp_ntt_add_partial(t5, t6, p);

  p1 = sp_ntt_mul(p1, ntt_const[1], ntt_const[NC+1], p);
  p2 = sp_ntt_mul(p2, ntt_const[2], ntt_const[NC+2], p);
  p3 = sp_ntt_mul(p3, ntt_const[3], ntt_const[NC+3], p);
  p4 = sp_ntt_mul(p4, ntt_const[4], ntt_const[NC+4], p);
  p5 = sp_ntt_mul(p5, ntt_const[5], ntt_const[NC+5], p);
  p6 = sp_ntt_mul(p6, ntt_const[6], ntt_const[NC+6], p);
  p7 = sp_ntt_mul(p7, ntt_const[7], ntt_const[NC+7], p);
  p8 = sp_ntt_mul(p8, ntt_const[8], ntt_const[NC+8], p);

  t1 = sp_ntt_add(p0, p1, p);
  t2 = sp_ntt_add(p2, p4, p);
  t3 = sp_ntt_add(p3, p4, p);
  t4 = p5;
  t5 = sp_ntt_add(p6, p8, p);
  t6 = sp_ntt_add(p7, p8, p);

  p1 = sp_ntt_add(t1, t2, p);
  p2 = sp_ntt_add(t1, t3, p);
  p3 = sp_ntt_sub(t1, t2, p);
  p3 = sp_ntt_sub(p3, t3, p);
  p4 = sp_ntt_add(t4, t5, p);
  p5 = sp_ntt_sub(t6, t4, p);
  p6 = sp_ntt_sub(t4, t5, p);
  p6 = sp_ntt_add(p6, t6, p);

  t1 = sp_ntt_add(p1, p4, p);
  t2 = sp_ntt_add(p2, p5, p);
  t3 = sp_ntt_add(p3, p6, p);
  t4 = sp_ntt_sub(p1, p4, p);
  t5 = sp_ntt_sub(p2, p5, p);
  t6 = sp_ntt_sub(p3, p6, p);

  x[j0] = p0;
  x[j1] = t6;
  x[j2] = t4;
  x[j3] = t5;
  x[j4] = t2;
  x[j5] = t1;
  x[j6] = t3;
}

DECLARE_CORE_ROUTINES(7)
