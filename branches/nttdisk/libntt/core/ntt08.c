#include "ntt-impl.h"

static uint32_t 
ntt8_get_num_const(void)
{
  return 8;
}

void
ntt8_init(spv_t out, sp_t p, sp_t d, 
	  sp_t primroot, sp_t order)
{
  sp_t w1 = sp_pow(primroot, order / 8, p, d);
  sp_t w2 = sp_sqr(w1, p, d);
  sp_t w3 = sp_mul(w2, w1, p, d);
  sp_t inv2 = sp_inv(2, p, d);

  out[0] = 1;
  out[1] = 1;
  out[2] = 1;
  out[3] = w2;
  out[4] = 1;
  out[5] = sp_mul(inv2, sp_sub(w1, w3, p), p, d);
  out[6] = w2;
  out[7] = sp_mul(inv2, sp_add(w1, w3, p), p, d);
}

static void
ntt8_run(spv_t x, spv_size_t stride,
	  sp_t p, sp_t d, spv_t ntt_const)
{
  sp_t x0, x1, x2, x3, x4, x5, x6, x7;
  sp_t t0, t1, t2, t3, t4, t5, t6, t7; 
  sp_t p0, p1, p2, p3, p4, p5, p6, p7;

  x0 = x[0 * stride];
  x1 = x[1 * stride];
  x2 = x[2 * stride];
  x3 = x[3 * stride];
  x4 = x[4 * stride];
  x5 = x[5 * stride];
  x6 = x[6 * stride];
  x7 = x[7 * stride];

  t0 = sp_add(x0, x4, p);
  t4 = sp_sub(x0, x4, p);
  t1 = sp_add(x1, x5, p);
  t5 = sp_sub(x1, x5, p);
  t2 = sp_add(x2, x6, p);
  t6 = sp_sub(x2, x6, p);
  t3 = sp_add(x3, x7, p);
  t7 = sp_sub(x3, x7, p);

  p0 = sp_add(t0, t2, p);
  p1 = sp_sub(t0, t2, p);
  p2 = sp_add(t1, t3, p);
  p3 = sp_sub(t1, t3, p);
  p4 = t4;
  p5 = sp_sub(t5, t7, p);
  p6 = t6;
  p7 = sp_add(t5, t7, p); 

  p3 = sp_mul(p3, ntt_const[3], p, d);
  p5 = sp_mul(p5, ntt_const[5], p, d);
  p6 = sp_mul(p6, ntt_const[6], p, d);
  p7 = sp_mul(p7, ntt_const[7], p, d);

  t0 = sp_add(p4, p6, p);
  t1 = sp_sub(p4, p6, p);
  t2 = sp_add(p5, p7, p);
  t3 = sp_sub(p5, p7, p);
  t4 = sp_add(t0, t2, p);
  t5 = sp_sub(t0, t2, p);
  t6 = sp_add(t1, t3, p);
  t7 = sp_sub(t1, t3, p);

  t0 = sp_add(p0, p2, p);
  t1 = sp_sub(p0, p2, p);
  t2 = sp_add(p1, p3, p);
  t3 = sp_sub(p1, p3, p);

  x[0 * stride] = t0;
  x[1 * stride] = t4;
  x[2 * stride] = t2;
  x[3 * stride] = t7;
  x[4 * stride] = t1;
  x[5 * stride] = t5;
  x[6 * stride] = t3;
  x[7 * stride] = t6;
}

const nttconfig_t ntt8_config = 
{
  8,
  ntt8_get_num_const,
  ntt8_init,
  ntt8_run,
  NULL
};

