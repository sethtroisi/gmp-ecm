#include "ntt-impl.h"

static uint32_t 
ntt3_get_num_const(void)
{
  return 3;
}

void
ntt3_init(spv_t out, sp_t p, sp_t d, 
	  sp_t primroot, sp_t order)
{
  sp_t w1, w2;
  sp_t h1, h2;
  sp_t inv2 = sp_inv(2, p, d);

  w1 = sp_pow(primroot, order / 3, p, d);
  w2 = sp_sqr(w1, p, d);

  h1 = sp_add(w1, w2, p);
  h2 = sp_sub(w1, w2, p);

  out[0] = 1;
  out[1] = sp_sub(sp_mul(h1, inv2, p, d), 1, p);
  out[2] = sp_mul(h2, inv2, p, d);
}

static void
ntt3_run(spv_t x, spv_size_t stride,
	  sp_t p, sp_t d, spv_t ntt_const)
{
  sp_t p0, p1, p2;
  sp_t x0, x1, x2;
  sp_t     t1, t2;

  x0 = x[0 * stride];
  x1 = x[1 * stride];
  x2 = x[2 * stride];

  t1 = sp_add(x1, x2, p);
  t2 = sp_sub(x1, x2, p);

  p0 = sp_add(x0, t1, p);

  p1 = sp_mul(t1, ntt_const[1], p, d);
  p2 = sp_mul(t2, ntt_const[2], p, d);

  p1 = sp_add(p0, p1, p);

  t1 = sp_add(p1, p2, p);
  t2 = sp_sub(p1, p2, p);

  x[0 * stride] = p0;
  x[1 * stride] = t1;
  x[2 * stride] = t2;
}

static void
ntt3_pfa_run(spv_t x, spv_size_t stride,
	  spv_size_t cofactor,
	  sp_t p, sp_t d, spv_t ntt_const)
{
  spv_size_t i, jstart;
  spv_size_t n = 3 * cofactor * stride;
  spv_size_t inc = cofactor * stride;
  spv_size_t inc3 = 3 * stride;

  for (i = jstart = 0; i < cofactor; i++, jstart += inc3)
    {
      spv_size_t j0, j1, j2;

      sp_t p0, p1, p2;
      sp_t x0, x1, x2;
      sp_t     t1, t2;

      j0 = jstart;
      j1 = sp_array_inc(j0, inc, n);
      j2 = sp_array_inc(j0, 2 * inc, n);

      x0 = x[j0];
      x1 = x[j1];
      x2 = x[j2];

      t1 = sp_add(x1, x2, p);
      t2 = sp_sub(x1, x2, p);

      p0 = sp_add(x0, t1, p);

      p1 = sp_mul(t1, ntt_const[1], p, d);
      p2 = sp_mul(t2, ntt_const[2], p, d);

      p1 = sp_add(p0, p1, p);

      t1 = sp_add(p1, p2, p);
      t2 = sp_sub(p1, p2, p);

      x[j0] = p0;
      x[j1] = t1;
      x[j2] = t2;
    }
}

const nttconfig_t ntt3_config = 
{
  3,
  ntt3_get_num_const,
  ntt3_init,
  ntt3_run,
  ntt3_pfa_run
};

