#include "ntt/ntt-impl-scalar.h"

#define NC 8

static const uint8_t ntt8_fixed_const[NC] = {1, 1, 1, 0, 1, 0, 0, 0};

void
X(ntt8_init)(spv_t out, sp_t p, sp_t d, 
	  sp_t primroot, sp_t order, sp_t perm)
{
  X(nttdata_init_generic)(&X(ntt8_config), out, p, d, primroot, order, perm);
}

static void 
ntt8_run_core(spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		sp_t p, spv_t ntt_const)
{
  sp_t x0, x1, x2, x3, x4, x5, x6, x7;
  sp_t t0, t1, t2, t3, t4, t5, t6, t7; 
  sp_t p0, p1, p2, p3, p4, p5, p6, p7;

  x0 = in[0 * istride];
  x1 = in[1 * istride];
  x2 = in[2 * istride];
  x3 = in[3 * istride];
  x4 = in[4 * istride];
  x5 = in[5 * istride];
  x6 = in[6 * istride];
  x7 = in[7 * istride];

  t0 = sp_ntt_add(x0, x4, p);
  t4 = sp_ntt_sub(x0, x4, p);
  t1 = sp_ntt_add(x1, x5, p);
  t5 = sp_ntt_sub(x1, x5, p);
  t2 = sp_ntt_add(x2, x6, p);
  t6 = sp_ntt_sub_partial(x2, x6, p);
  t3 = sp_ntt_add(x3, x7, p);
  t7 = sp_ntt_sub(x3, x7, p);

  p0 = sp_ntt_add(t0, t2, p);
  p1 = sp_ntt_sub(t0, t2, p);
  p2 = sp_ntt_add(t1, t3, p);
  p3 = sp_ntt_sub_partial(t1, t3, p);
  p4 = t4;
  p5 = t6;
  p6 = sp_ntt_sub_partial(t5, t7, p);
  p7 = sp_ntt_add_partial(t5, t7, p); 

  p3 = sp_ntt_mul(p3, ntt_const[3], ntt_const[NC+3], p);
  p5 = sp_ntt_mul(p5, ntt_const[5], ntt_const[NC+5], p);
  p6 = sp_ntt_mul(p6, ntt_const[6], ntt_const[NC+6], p);
  p7 = sp_ntt_mul(p7, ntt_const[7], ntt_const[NC+7], p);

  t0 = sp_ntt_add(p4, p5, p);
  t1 = sp_ntt_sub(p4, p5, p);
  t2 = sp_ntt_add(p6, p7, p);
  t3 = sp_ntt_sub(p6, p7, p);
  t4 = sp_ntt_add(t0, t2, p);
  t5 = sp_ntt_sub(t0, t2, p);
  t6 = sp_ntt_add(t1, t3, p);
  t7 = sp_ntt_sub(t1, t3, p);

  t0 = sp_ntt_add(p0, p2, p);
  t1 = sp_ntt_sub(p0, p2, p);
  t2 = sp_ntt_add(p1, p3, p);
  t3 = sp_ntt_sub(p1, p3, p);

  out[0 * ostride] = t0;
  out[1 * ostride] = t4;
  out[2 * ostride] = t2;
  out[3 * ostride] = t7;
  out[4 * ostride] = t1;
  out[5 * ostride] = t5;
  out[6 * ostride] = t3;
  out[7 * ostride] = t6;
}

static void
ntt8_twiddle_run_core(spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		spv_t w, sp_t p, spv_t ntt_const)
{
  sp_t x0, x1, x2, x3, x4, x5, x6, x7;
  sp_t t0, t1, t2, t3, t4, t5, t6, t7; 
  sp_t p0, p1, p2, p3, p4, p5, p6, p7;

  x0 = in[0 * istride];
  x1 = in[1 * istride];
  x2 = in[2 * istride];
  x3 = in[3 * istride];
  x4 = in[4 * istride];
  x5 = in[5 * istride];
  x6 = in[6 * istride];
  x7 = in[7 * istride];

  t0 = sp_ntt_add(x0, x4, p);
  t4 = sp_ntt_sub(x0, x4, p);
  t1 = sp_ntt_add(x1, x5, p);
  t5 = sp_ntt_sub(x1, x5, p);
  t2 = sp_ntt_add(x2, x6, p);
  t6 = sp_ntt_sub_partial(x2, x6, p);
  t3 = sp_ntt_add(x3, x7, p);
  t7 = sp_ntt_sub(x3, x7, p);

  p0 = sp_ntt_add(t0, t2, p);
  p1 = sp_ntt_sub(t0, t2, p);
  p2 = sp_ntt_add(t1, t3, p);
  p3 = sp_ntt_sub_partial(t1, t3, p);
  p4 = t4;
  p5 = t6;
  p6 = sp_ntt_sub_partial(t5, t7, p);
  p7 = sp_ntt_add_partial(t5, t7, p); 

  p3 = sp_ntt_mul(p3, ntt_const[3], ntt_const[NC+3], p);
  p5 = sp_ntt_mul(p5, ntt_const[5], ntt_const[NC+5], p);
  p6 = sp_ntt_mul(p6, ntt_const[6], ntt_const[NC+6], p);
  p7 = sp_ntt_mul(p7, ntt_const[7], ntt_const[NC+7], p);

  t0 = sp_ntt_add(p4, p5, p);
  t1 = sp_ntt_sub(p4, p5, p);
  t2 = sp_ntt_add(p6, p7, p);
  t3 = sp_ntt_sub(p6, p7, p);
  t4 = sp_ntt_add_partial(t0, t2, p);
  t5 = sp_ntt_sub_partial(t0, t2, p);
  t6 = sp_ntt_add_partial(t1, t3, p);
  t7 = sp_ntt_sub_partial(t1, t3, p);

  t0 = sp_ntt_add(p0, p2, p);
  t1 = sp_ntt_sub_partial(p0, p2, p);
  t2 = sp_ntt_add_partial(p1, p3, p);
  t3 = sp_ntt_sub_partial(p1, p3, p);

  t4 = sp_ntt_mul(t4, w[0], w[1], p);
  t2 = sp_ntt_mul(t2, w[2], w[3], p);
  t7 = sp_ntt_mul(t7, w[4], w[5], p);
  t1 = sp_ntt_mul(t1, w[6], w[7], p);
  t5 = sp_ntt_mul(t5, w[8], w[9], p);
  t3 = sp_ntt_mul(t3, w[10], w[11], p);
  t6 = sp_ntt_mul(t6, w[12], w[13], p);

  out[0 * ostride] = t0;
  out[1 * ostride] = t4;
  out[2 * ostride] = t2;
  out[3 * ostride] = t7;
  out[4 * ostride] = t1;
  out[5 * ostride] = t5;
  out[6 * ostride] = t3;
  out[7 * ostride] = t6;
}

static void
ntt8_pfa_run_core(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t n,
	  sp_t p, spv_t ntt_const)
{
  spv_size_t j0, j1, j2, j3, j4, j5, j6, j7;
  sp_t x0, x1, x2, x3, x4, x5, x6, x7;
  sp_t t0, t1, t2, t3, t4, t5, t6, t7; 
  sp_t p0, p1, p2, p3, p4, p5, p6, p7;

  j0 = start;
  j1 = sp_array_inc(j0, inc, n);
  j2 = sp_array_inc(j0, 2 * inc, n);
  j3 = sp_array_inc(j0, 3 * inc, n);
  j4 = sp_array_inc(j0, 4 * inc, n);
  j5 = sp_array_inc(j0, 5 * inc, n);
  j6 = sp_array_inc(j0, 6 * inc, n);
  j7 = sp_array_inc(j0, 7 * inc, n);

  x0 = x[j0];
  x1 = x[j1];
  x2 = x[j2];
  x3 = x[j3];
  x4 = x[j4];
  x5 = x[j5];
  x6 = x[j6];
  x7 = x[j7];

  t0 = sp_ntt_add(x0, x4, p);
  t4 = sp_ntt_sub(x0, x4, p);
  t1 = sp_ntt_add(x1, x5, p);
  t5 = sp_ntt_sub(x1, x5, p);
  t2 = sp_ntt_add(x2, x6, p);
  t6 = sp_ntt_sub_partial(x2, x6, p);
  t3 = sp_ntt_add(x3, x7, p);
  t7 = sp_ntt_sub(x3, x7, p);

  p0 = sp_ntt_add(t0, t2, p);
  p1 = sp_ntt_sub(t0, t2, p);
  p2 = sp_ntt_add(t1, t3, p);
  p3 = sp_ntt_sub_partial(t1, t3, p);
  p4 = t4;
  p5 = t6;
  p6 = sp_ntt_sub_partial(t5, t7, p);
  p7 = sp_ntt_add_partial(t5, t7, p); 

  p3 = sp_ntt_mul(p3, ntt_const[3], ntt_const[NC+3], p);
  p5 = sp_ntt_mul(p5, ntt_const[5], ntt_const[NC+5], p);
  p6 = sp_ntt_mul(p6, ntt_const[6], ntt_const[NC+6], p);
  p7 = sp_ntt_mul(p7, ntt_const[7], ntt_const[NC+7], p);

  t0 = sp_ntt_add(p4, p5, p);
  t1 = sp_ntt_sub(p4, p5, p);
  t2 = sp_ntt_add(p6, p7, p);
  t3 = sp_ntt_sub(p6, p7, p);
  t4 = sp_ntt_add(t0, t2, p);
  t5 = sp_ntt_sub(t0, t2, p);
  t6 = sp_ntt_add(t1, t3, p);
  t7 = sp_ntt_sub(t1, t3, p);

  t0 = sp_ntt_add(p0, p2, p);
  t1 = sp_ntt_sub(p0, p2, p);
  t2 = sp_ntt_add(p1, p3, p);
  t3 = sp_ntt_sub(p1, p3, p);

  x[j0] = t0;
  x[j1] = t4;
  x[j2] = t2;
  x[j3] = t7;
  x[j4] = t1;
  x[j5] = t5;
  x[j6] = t3;
  x[j7] = t6;
}

DECLARE_CORE_ROUTINES(8)
