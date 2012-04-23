#include "ntt-impl.h"

static uint32_t 
ntt15_get_num_const(void)
{
  return 3 * 6;
}

extern void ntt3_init(spv_t out, sp_t p, sp_t d, 
			sp_t primroot, sp_t order);

extern void ntt5_init(spv_t out, sp_t p, sp_t d, 
			sp_t primroot, sp_t order);

void
ntt15_init(spv_t out, sp_t p, sp_t d, 
	  sp_t primroot, sp_t order)
{
  uint32_t i, j;
  uint32_t num3 = 3;
  uint32_t num5 = 6;
  spv_t root3 = (spv_t)alloca(num3 * sizeof(sp_t));
  spv_t root5 = (spv_t)alloca(num5 * sizeof(sp_t));

  ntt3_init(root3, p, d, primroot, order);
  ntt5_init(root5, p, d, primroot, order);

  for (i = 0; i < num3; i++)
    for (j = 0; j < num5; j++)
      out[num5 * i + j] = sp_mul(root3[i], root5[j], p, d);
}

static void
ntt15_run(spv_t x, spv_size_t stride,
	  sp_t p, sp_t d, spv_t ntt_const)
{
  sp_t a00, a01, a02, a03, a04, 
       a05, a06, a07, a08, a09, 
       a10, a11, a12, a13, a14;

  {
    sp_t x00 = x[ 0 * stride];
    sp_t x01 = x[ 5 * stride];
    sp_t x02 = x[10 * stride];

    a01 = sp_add(x01, x02, p);
    a02 = sp_sub(x01, x02, p);

    a00 = sp_add(x00, a01, p);
  }
  {
    sp_t x03 = x[ 3 * stride];
    sp_t x04 = x[ 8 * stride];
    sp_t x05 = x[13 * stride];

    a04 = sp_add(x04, x05, p);
    a05 = sp_sub(x04, x05, p);

    a03 = sp_add(x03, a04, p);
  }
  {
    sp_t x08 = x[ 1 * stride];
    sp_t x06 = x[ 6 * stride];
    sp_t x07 = x[11 * stride];

    a07 = sp_add(x07, x08, p);
    a08 = sp_sub(x07, x08, p);

    a06 = sp_add(x06, a07, p);
  }
  {
    sp_t x11 = x[ 4 * stride];
    sp_t x09 = x[ 9 * stride];
    sp_t x10 = x[14 * stride];

    a10 = sp_add(x10, x11, p);
    a11 = sp_sub(x10, x11, p);

    a09 = sp_add(x09, a10, p);
  }
  {
    sp_t x13 = x[ 2 * stride];
    sp_t x14 = x[ 7 * stride];
    sp_t x12 = x[12 * stride];

    a13 = sp_add(x13, x14, p);
    a14 = sp_sub(x13, x14, p);

    a12 = sp_add(x12, a13, p);
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a00;
    b1 = a03;
    b4 = a06;
    b2 = a09;
    b3 = a12;

    c1 = sp_add(b1, b3, p);
    c3 = sp_sub(b1, b3, p);
    c2 = sp_add(b2, b4, p);
    c4 = sp_sub(b2, b4, p);

    b1 = sp_add(c1, c2, p);
    b2 = sp_sub(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_add(c3, c4, p);

    b0 = sp_add(b0, b1, p);

    b1 = sp_mul(b1, ntt_const[1], p, d);
    b2 = sp_mul(b2, ntt_const[2], p, d);
    b3 = sp_mul(b3, ntt_const[3], p, d);
    b4 = sp_mul(b4, ntt_const[4], p, d);
    b5 = sp_mul(b5, ntt_const[5], p, d);

    b1 = sp_add(b0, b1, p);

    c1 = sp_add(b1, b2, p);
    c2 = sp_sub(b1, b2, p);
    c3 = sp_add(b3, b5, p);
    c4 = sp_add(b4, b5, p);

    b1 = sp_add(c1, c3, p);
    b2 = sp_add(c2, c4, p);
    b3 = sp_sub(c1, c3, p);
    b4 = sp_sub(c2, c4, p);

    a00 = b0;
    a03 = b4;
    a06 = b3;
    a09 = b1;
    a12 = b2;
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a01;
    b1 = a04;
    b4 = a07;
    b2 = a10;
    b3 = a13;

    c1 = sp_add(b1, b3, p);
    c3 = sp_sub(b1, b3, p);
    c2 = sp_add(b2, b4, p);
    c4 = sp_sub(b2, b4, p);

    b1 = sp_add(c1, c2, p);
    b2 = sp_sub(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_add(c3, c4, p);

    b0 = sp_add(b0, b1, p);

    b0 = sp_mul(b0, ntt_const[6], p, d);
    b1 = sp_mul(b1, ntt_const[7], p, d);
    b2 = sp_mul(b2, ntt_const[8], p, d);
    b3 = sp_mul(b3, ntt_const[9], p, d);
    b4 = sp_mul(b4, ntt_const[10], p, d);
    b5 = sp_mul(b5, ntt_const[11], p, d);

    b1 = sp_add(b0, b1, p);

    c1 = sp_add(b1, b2, p);
    c2 = sp_sub(b1, b2, p);
    c3 = sp_add(b3, b5, p);
    c4 = sp_add(b4, b5, p);

    b1 = sp_add(c1, c3, p);
    b2 = sp_add(c2, c4, p);
    b3 = sp_sub(c1, c3, p);
    b4 = sp_sub(c2, c4, p);

    a01 = b0;
    a04 = b4;
    a07 = b3;
    a10 = b1;
    a13 = b2;
  }

  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a02;
    b1 = a05;
    b4 = a08;
    b2 = a11;
    b3 = a14;

    c1 = sp_add(b1, b3, p);
    c3 = sp_sub(b1, b3, p);
    c2 = sp_add(b2, b4, p);
    c4 = sp_sub(b2, b4, p);

    b1 = sp_add(c1, c2, p);
    b2 = sp_sub(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_add(c3, c4, p);

    b0 = sp_add(b0, b1, p);

    b0 = sp_mul(b0, ntt_const[12], p, d);
    b1 = sp_mul(b1, ntt_const[13], p, d);
    b2 = sp_mul(b2, ntt_const[14], p, d);
    b3 = sp_mul(b3, ntt_const[15], p, d);
    b4 = sp_mul(b4, ntt_const[16], p, d);
    b5 = sp_mul(b5, ntt_const[17], p, d);

    b1 = sp_add(b0, b1, p);

    c1 = sp_add(b1, b2, p);
    c2 = sp_sub(b1, b2, p);
    c3 = sp_add(b3, b5, p);
    c4 = sp_add(b4, b5, p);

    b1 = sp_add(c1, c3, p);
    b2 = sp_add(c2, c4, p);
    b3 = sp_sub(c1, c3, p);
    b4 = sp_sub(c2, c4, p);

    a02 = b0;
    a05 = b4;
    a08 = b3;
    a11 = b1;
    a14 = b2;
  }
  {
    sp_t x00, x01, x02;

    a01 = sp_add(a00, a01, p);

    x00 = a00;
    x01 = sp_add(a01, a02, p);
    x02 = sp_sub(a01, a02, p);

    x[ 0 * stride] = x00;
    x[ 5 * stride] = x02;
    x[10 * stride] = x01;
  }
  {
    sp_t x03, x04, x05;

    a04 = sp_add(a03, a04, p);

    x03 = a03;
    x04 = sp_add(a04, a05, p);
    x05 = sp_sub(a04, a05, p);

    x[ 1 * stride] = x04;
    x[ 6 * stride] = x03;
    x[11 * stride] = x05;
  }
  {
    sp_t x06, x07, x08;

    a07 = sp_add(a06, a07, p);

    x06 = a06;
    x07 = sp_add(a07, a08, p);
    x08 = sp_sub(a07, a08, p);

    x[ 2 * stride] = x08;
    x[ 7 * stride] = x07;
    x[12 * stride] = x06;
  }
  {
    sp_t x09, x10, x11;

    a10 = sp_add(a09, a10, p);

    x09 = a09;
    x10 = sp_add(a10, a11, p);
    x11 = sp_sub(a10, a11, p);

    x[ 3 * stride] = x09;
    x[ 8 * stride] = x11;
    x[13 * stride] = x10;
  }
  {
    sp_t x12, x13, x14;

    a13 = sp_add(a12, a13, p);

    x12 = a12;
    x13 = sp_add(a13, a14, p);
    x14 = sp_sub(a13, a14, p);

    x[ 4 * stride] = x13;
    x[ 9 * stride] = x12;
    x[14 * stride] = x14;
  }

}

const nttconfig_t ntt15_config = 
{
  15,
  ntt15_get_num_const,
  ntt15_init,
  ntt15_run,
  NULL
};

