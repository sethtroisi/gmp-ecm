#include "ntt-impl.h"

static uint32_t 
ntt5_get_num_const(void)
{
  return 6;
}

void
ntt5_init(spv_t out, sp_t p, sp_t d, 
	  sp_t primroot, sp_t order)
{
  sp_t w1, w2, w3, w4;
  sp_t h1, h2, h3, h4, h5;
  sp_t t1, t2, t3, t4;
  sp_t inv4 = sp_inv(4, p, d);
    
  h1 = sp_pow(primroot, order / 5, p, d);
  h2 = sp_mul(h1, h1, p, d);
  h3 = sp_mul(h2, h1, p, d);
  h4 = sp_mul(h3, h1, p, d);

  w1 = h3;
  w2 = h4;
  w3 = h2;
  w4 = h1;

  t1 = sp_add(w1, w2, p);
  t1 = sp_add(t1, w3, p);
  t1 = sp_add(t1, w4, p);

  t2 = sp_sub(w1, w2, p);
  t2 = sp_add(t2, w3, p);
  t2 = sp_sub(t2, w4, p);

  t3 = sp_add(w1, w1, p);
  t3 = sp_sub(t3, w3, p);
  t3 = sp_sub(t3, w3, p);

  t4 = sp_add(w2, w2, p);
  t4 = sp_sub(t4, w4, p);
  t4 = sp_sub(t4, w4, p);

  h1 = t1;
  h2 = t2;
  h3 = sp_sub(t3, t4, p);
  h4 = sp_neg(sp_add(t3, t4, p), p);
  h5 = t4;

  out[0] = 1;
  out[1] = sp_sub(sp_mul(h1, inv4, p, d), 1, p);
  out[2] = sp_mul(h2, inv4, p, d);
  out[3] = sp_mul(h3, inv4, p, d);
  out[4] = sp_mul(h4, inv4, p, d);
  out[5] = sp_mul(h5, inv4, p, d);
}

static void
ntt5_run(spv_t x, spv_size_t stride,
	  sp_t p, sp_t d, spv_t ntt_const)
{
  sp_t p0, p1, p2, p3, p4, p5;
  sp_t x0, x1, x2, x3, x4;
  sp_t     t1, t2, t3, t4;

  x0 = x[0 * stride];
  x1 = x[1 * stride];
  x4 = x[2 * stride];
  x2 = x[3 * stride];
  x3 = x[4 * stride];

  t1 = sp_add(x1, x3, p);
  t3 = sp_sub(x1, x3, p);
  t2 = sp_add(x2, x4, p);
  t4 = sp_sub(x2, x4, p);

  p1 = sp_add(t1, t2, p);
  p2 = sp_sub(t1, t2, p);
  p3 = t3;
  p4 = t4;
  p5 = sp_add(t3, t4, p);

  p0 = sp_add(x0, p1, p);

  p1 = sp_mul(p1, ntt_const[1], p, d);
  p2 = sp_mul(p2, ntt_const[2], p, d);
  p3 = sp_mul(p3, ntt_const[3], p, d);
  p4 = sp_mul(p4, ntt_const[4], p, d);
  p5 = sp_mul(p5, ntt_const[5], p, d);

  p1 = sp_add(p0, p1, p);

  t1 = sp_add(p1, p2, p);
  t2 = sp_sub(p1, p2, p);
  t3 = sp_add(p3, p5, p);
  t4 = sp_add(p4, p5, p);

  p1 = sp_add(t1, t3, p);
  p2 = sp_add(t2, t4, p);
  p3 = sp_sub(t1, t3, p);
  p4 = sp_sub(t2, t4, p);

  x[0 * stride] = p0;
  x[1 * stride] = p4;
  x[2 * stride] = p3;
  x[3 * stride] = p1;
  x[4 * stride] = p2;
}

const nttconfig_t ntt5_config = 
{
  5,
  ntt5_get_num_const,
  ntt5_init,
  ntt5_run,
  NULL
};

