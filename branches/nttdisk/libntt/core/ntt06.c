#include "ntt-impl.h"

static uint32_t 
ntt6_get_num_const(void)
{
  return 6;
}

extern void ntt3_init(spv_t out, sp_t p, sp_t d, 
			sp_t primroot, sp_t order);

void
ntt6_init(spv_t out, sp_t p, sp_t d, 
	  sp_t primroot, sp_t order)
{
  ntt3_init(out, p, d, primroot, order);
  ntt3_init(out + 3, p, d, primroot, order);
}

static void
ntt6_run(spv_t x, spv_size_t stride,
	  sp_t p, sp_t d, spv_t ntt_const)
{
  sp_t t0, t1, t2, t3, t4, t5;

  {
    sp_t x0 = x[0 * stride];
    sp_t x1 = x[2 * stride];
    sp_t x2 = x[4 * stride];
    sp_t b0, b1, b2;

    b1 = sp_add(x1, x2, p);
    b2 = sp_sub(x1, x2, p);

    b0 = sp_add(x0, b1, p);

    b1 = sp_mul(b1, ntt_const[1], p, d);
    b2 = sp_mul(b2, ntt_const[2], p, d);

    b1 = sp_add(b0, b1, p);

    t0 = b0;
    t1 = sp_add(b1, b2, p);
    t2 = sp_sub(b1, b2, p);
  }
  {
    sp_t x5 = x[1 * stride];
    sp_t x3 = x[3 * stride];
    sp_t x4 = x[5 * stride];
    sp_t b0, b1, b2;

    b1 = sp_add(x4, x5, p);
    b2 = sp_sub(x4, x5, p);

    b0 = sp_add(x3, b1, p);

    b1 = sp_mul(b1, ntt_const[4], p, d);
    b2 = sp_mul(b2, ntt_const[5], p, d);

    b1 = sp_add(b0, b1, p);

    t3 = b0;
    t4 = sp_add(b1, b2, p);
    t5 = sp_sub(b1, b2, p);
  }
  {
    sp_t x0 = sp_add(t0, t3, p);
    sp_t x3 = sp_sub(t0, t3, p);

    x[0 * stride] = x0;
    x[3 * stride] = x3;
  }
  {
    sp_t x1 = sp_add(t1, t4, p);
    sp_t x4 = sp_sub(t1, t4, p);

    x[1 * stride] = x4;
    x[4 * stride] = x1;
  }
  {
    sp_t x2 = sp_add(t2, t5, p);
    sp_t x5 = sp_sub(t2, t5, p);

    x[2 * stride] = x2;
    x[5 * stride] = x5;
  }

}

const nttconfig_t ntt6_config = 
{
  6,
  ntt6_get_num_const,
  ntt6_init,
  ntt6_run,
  NULL
};

