#include "ntt-impl.h"

static uint32_t 
ntt7_get_num_const(void)
{
  return 9;
}

void
ntt7_init(spv_t out, sp_t p, sp_t d, 
	  sp_t primroot, sp_t order)
{
  sp_t w1, w2, w3, w4, w5, w6;
  sp_t h1, h2, h3, h4, h5, h6, h7, h8;
  sp_t t1, t2, t3, t4, t5, t6, t7, t8;
  sp_t inv6 = sp_inv(6, p, d);
    
  h1 = sp_pow(primroot, order / 7, p, d);
  h2 = sp_sqr(h1, p, d);
  h3 = sp_mul(h2, h1, p, d);
  h4 = sp_mul(h3, h1, p, d);
  h5 = sp_mul(h4, h1, p, d);
  h6 = sp_mul(h5, h1, p, d);

  w1 = h5;
  w2 = h4;
  w3 = h6;
  w4 = h2;
  w5 = h3;
  w6 = h1;

  t1 = sp_add(w1, w2, p);
  t1 = sp_add(t1, w3, p);
  t1 = sp_add(t1, w4, p);
  t1 = sp_add(t1, w5, p);
  t1 = sp_add(t1, w6, p);

  t2 = sp_add(w1, w1, p);
  t2 = sp_sub(t2, w2, p);
  t2 = sp_sub(t2, w3, p);
  t2 = sp_add(t2, w4, p);
  t2 = sp_add(t2, w4, p);
  t2 = sp_sub(t2, w5, p);
  t2 = sp_sub(t2, w6, p);

  t3 = sp_neg(w1, p);
  t3 = sp_add(t3, w2, p);
  t3 = sp_add(t3, w2, p);
  t3 = sp_sub(t3, w3, p);
  t3 = sp_sub(t3, w4, p);
  t3 = sp_add(t3, w5, p);
  t3 = sp_add(t3, w5, p);
  t3 = sp_sub(t3, w6, p);

  t4 = sp_sub(w1, w2, p);
  t4 = sp_add(t4, w3, p);
  t4 = sp_sub(t4, w4, p);
  t4 = sp_add(t4, w5, p);
  t4 = sp_sub(t4, w6, p);

  t5 = sp_add(w1, w1, p);
  t5 = sp_add(t5, w2, p);
  t5 = sp_sub(t5, w3, p);
  t5 = sp_sub(t5, w4, p);
  t5 = sp_sub(t5, w4, p);
  t5 = sp_sub(t5, w5, p);
  t5 = sp_add(t5, w6, p);

  t6 = sp_add(w1, w2, p);
  t6 = sp_add(t6, w2, p);
  t6 = sp_add(t6, w3, p);
  t6 = sp_sub(t6, w4, p);
  t6 = sp_sub(t6, w5, p);
  t6 = sp_sub(t6, w5, p);
  t6 = sp_sub(t6, w6, p);

  h1 = t1;
  h2 = sp_sub(t2, t3, p);
  h3 = sp_neg(sp_add(t2, t3, p), p);
  h3 = sp_sub(h3, t3, p);
  h4 = t3;
  h5 = t4;
  h6 = sp_sub(t5, t6, p);
  h7 = sp_neg(t5, p); 
  h8 = t6;

  out[0] = 1;
  out[1] = sp_sub(sp_mul(h1, inv6, p, d), 1, p);
  out[2] = sp_mul(h2, inv6, p, d);
  out[3] = sp_mul(h3, inv6, p, d);
  out[4] = sp_mul(h4, inv6, p, d);
  out[5] = sp_mul(h5, inv6, p, d);
  out[6] = sp_mul(h6, inv6, p, d);
  out[7] = sp_mul(h7, inv6, p, d);
  out[8] = sp_mul(h8, inv6, p, d);
}

static void
ntt7_run(spv_t x, spv_size_t stride,
	  sp_t p, sp_t d, spv_t ntt_const)
{
  sp_t x0, x1, x2, x3, x4, x5, x6;
  sp_t     t1, t2, t3, t4, t5, t6, t7, t8;
  sp_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

  x0 = x[0 * stride];
  x1 = x[1 * stride];
  x5 = x[2 * stride];
  x6 = x[3 * stride];
  x3 = x[4 * stride];
  x2 = x[5 * stride];
  x4 = x[6 * stride];

  p1 = sp_add(x1, x4, p);
  p2 = sp_add(x2, x5, p);
  p3 = sp_add(x3, x6, p);
  p4 = sp_sub(x1, x4, p);
  p5 = sp_sub(x2, x5, p);
  p6 = sp_sub(x3, x6, p);

  t1 = sp_add(p1, p2, p);
  t1 = sp_add(t1, p3, p);
  t2 = sp_sub(p1, p3, p);
  t3 = sp_sub(p2, p3, p);
  t4 = sp_sub(p4, p5, p);
  t4 = sp_add(t4, p6, p);
  t5 = sp_sub(p4, p6, p);
  t6 = sp_add(p5, p6, p);

  p0 = sp_add(x0, t1, p);
  p1 = t1;
  p2 = t2;
  p3 = t3;
  p4 = sp_add(t2, t3, p);
  p5 = t4;
  p6 = t5;
  p7 = t6;
  p8 = sp_add(t5, t6, p);

  p1 = sp_mul(p1, ntt_const[1], p, d);
  p2 = sp_mul(p2, ntt_const[2], p, d);
  p3 = sp_mul(p3, ntt_const[3], p, d);
  p4 = sp_mul(p4, ntt_const[4], p, d);
  p5 = sp_mul(p5, ntt_const[5], p, d);
  p6 = sp_mul(p6, ntt_const[6], p, d);
  p7 = sp_mul(p7, ntt_const[7], p, d);
  p8 = sp_mul(p8, ntt_const[8], p, d);

  t1 = sp_add(p0, p1, p);
  t2 = sp_add(p2, p4, p);
  t3 = sp_add(p3, p4, p);
  t4 = p5;
  t5 = sp_add(p6, p8, p);
  t6 = sp_add(p7, p8, p);

  p1 = sp_add(t1, t2, p);
  p2 = sp_add(t1, t3, p);
  p3 = sp_sub(t1, t2, p);
  p3 = sp_sub(p3, t3, p);
  p4 = sp_add(t4, t5, p);
  p5 = sp_sub(t6, t4, p);
  p6 = sp_sub(t4, t5, p);
  p6 = sp_add(p6, t6, p);

  t1 = sp_add(p1, p4, p);
  t2 = sp_add(p2, p5, p);
  t3 = sp_add(p3, p6, p);
  t4 = sp_sub(p1, p4, p);
  t5 = sp_sub(p2, p5, p);
  t6 = sp_sub(p3, p6, p);

  x[0 * stride] = p0;
  x[1 * stride] = t6;
  x[2 * stride] = t4;
  x[3 * stride] = t5;
  x[4 * stride] = t2;
  x[5 * stride] = t1;
  x[6 * stride] = t3;
}

const nttconfig_t ntt7_config = 
{
  7,
  ntt7_get_num_const,
  ntt7_init,
  ntt7_run,
  NULL
};

