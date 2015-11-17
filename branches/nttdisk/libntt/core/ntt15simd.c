#include "ntt/ntt-impl-simd.h"

#define NC 18

static const uint8_t ntt15_fixed_const[NC] = {1};

extern void X(ntt15_init)(spv_t out, sp_t p, sp_t d, sp_t primroot, 
                        sp_t order, sp_t perm);


static void
ntt15_run_core_simd(spv_t in, spv_size_t istride, spv_size_t idist,
		spv_t out, spv_size_t ostride, spv_size_t odist,
		sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  sp_simd_t a00, a01, a02, a03, a04, 
            a05, a06, a07, a08, a09, 
            a10, a11, a12, a13, a14;

  {
    sp_simd_t x00 = sp_simd_gather(in + 0 * istride, idist, vsize);
    sp_simd_t x01 = sp_simd_gather(in + 5 * istride, idist, vsize);
    sp_simd_t x02 = sp_simd_gather(in +10 * istride, idist, vsize);

    a01 = sp_ntt_add_simd(x01, x02, p);
    a02 = sp_ntt_sub_simd(x01, x02, p);

    a00 = sp_ntt_add_simd(x00, a01, p);
  }
  {
    sp_simd_t x03 = sp_simd_gather(in + 3 * istride, idist, vsize);
    sp_simd_t x04 = sp_simd_gather(in + 8 * istride, idist, vsize);
    sp_simd_t x05 = sp_simd_gather(in +13 * istride, idist, vsize);

    a04 = sp_ntt_add_simd(x04, x05, p);
    a05 = sp_ntt_sub_simd(x04, x05, p);

    a03 = sp_ntt_add_simd(x03, a04, p);
  }
  {
    sp_simd_t x08 = sp_simd_gather(in + 1 * istride, idist, vsize);
    sp_simd_t x06 = sp_simd_gather(in + 6 * istride, idist, vsize);
    sp_simd_t x07 = sp_simd_gather(in +11 * istride, idist, vsize);

    a07 = sp_ntt_add_simd(x07, x08, p);
    a08 = sp_ntt_sub_simd(x07, x08, p);

    a06 = sp_ntt_add_simd(x06, a07, p);
  }
  {
    sp_simd_t x11 = sp_simd_gather(in + 4 * istride, idist, vsize);
    sp_simd_t x09 = sp_simd_gather(in + 9 * istride, idist, vsize);
    sp_simd_t x10 = sp_simd_gather(in +14 * istride, idist, vsize);

    a10 = sp_ntt_add_simd(x10, x11, p);
    a11 = sp_ntt_sub_simd(x10, x11, p);

    a09 = sp_ntt_add_simd(x09, a10, p);
  }
  {
    sp_simd_t x13 = sp_simd_gather(in + 2 * istride, idist, vsize);
    sp_simd_t x14 = sp_simd_gather(in + 7 * istride, idist, vsize);
    sp_simd_t x12 = sp_simd_gather(in +12 * istride, idist, vsize);

    a13 = sp_ntt_add_simd(x13, x14, p);
    a14 = sp_ntt_sub_simd(x13, x14, p);

    a12 = sp_ntt_add_simd(x12, a13, p);
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a00;
    b1 = a03;
    b4 = a06;
    b2 = a09;
    b3 = a12;

    c1 = sp_ntt_add_simd(b1, b3, p);
    c3 = sp_ntt_sub_simd(b1, b3, p);
    c2 = sp_ntt_add_simd(b2, b4, p);
    c4 = sp_ntt_sub_simd(b2, b4, p);

    b1 = sp_ntt_add_simd(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd(c3, c4, p);

    b0 = sp_ntt_add_simd(b0, b1, p);

    b1 = sp_ntt_mul_simd(b1, ntt_const[1], ntt_const[NC+1], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[2], ntt_const[NC+2], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[3], ntt_const[NC+3], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[4], ntt_const[NC+4], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[5], ntt_const[NC+5], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a00 = b0;
    a03 = b4;
    a06 = b3;
    a09 = b1;
    a12 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a01;
    b1 = a04;
    b4 = a07;
    b2 = a10;
    b3 = a13;

    c1 = sp_ntt_add_simd(b1, b3, p);
    c3 = sp_ntt_sub_simd(b1, b3, p);
    c2 = sp_ntt_add_simd(b2, b4, p);
    c4 = sp_ntt_sub_simd(b2, b4, p);

    b1 = sp_ntt_add_simd(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd(c3, c4, p);

    b0 = sp_ntt_add_partial_simd(b0, b1, p);

    b0 = sp_ntt_mul_simd(b0, ntt_const[6], ntt_const[NC+6], p);
    b1 = sp_ntt_mul_simd(b1, ntt_const[7], ntt_const[NC+7], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[8], ntt_const[NC+8], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[9], ntt_const[NC+9], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[10], ntt_const[NC+10], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[11], ntt_const[NC+11], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a01 = b0;
    a04 = b4;
    a07 = b3;
    a10 = b1;
    a13 = b2;
  }

  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a02;
    b1 = a05;
    b4 = a08;
    b2 = a11;
    b3 = a14;

    c1 = sp_ntt_add_simd(b1, b3, p);
    c3 = sp_ntt_sub_simd(b1, b3, p);
    c2 = sp_ntt_add_simd(b2, b4, p);
    c4 = sp_ntt_sub_simd(b2, b4, p);

    b1 = sp_ntt_add_simd(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd(c3, c4, p);

    b0 = sp_ntt_add_partial_simd(b0, b1, p);

    b0 = sp_ntt_mul_simd(b0, ntt_const[12], ntt_const[NC+12], p);
    b1 = sp_ntt_mul_simd(b1, ntt_const[13], ntt_const[NC+13], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[14], ntt_const[NC+14], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[15], ntt_const[NC+15], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[16], ntt_const[NC+16], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[17], ntt_const[NC+17], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a02 = b0;
    a05 = b4;
    a08 = b3;
    a11 = b1;
    a14 = b2;
  }
  {
    sp_simd_t x00, x01, x02;

    a01 = sp_ntt_add_simd(a00, a01, p);

    x00 = a00;
    x01 = sp_ntt_add_simd(a01, a02, p);
    x02 = sp_ntt_sub_simd(a01, a02, p);

    sp_simd_scatter(x00, out +  0 * ostride, odist, vsize);
    sp_simd_scatter(x02, out +  5 * ostride, odist, vsize);
    sp_simd_scatter(x01, out + 10 * ostride, odist, vsize);
  }
  {
    sp_simd_t x03, x04, x05;

    a04 = sp_ntt_add_simd(a03, a04, p);

    x03 = a03;
    x04 = sp_ntt_add_simd(a04, a05, p);
    x05 = sp_ntt_sub_simd(a04, a05, p);

    sp_simd_scatter(x04, out +  1 * ostride, odist, vsize);
    sp_simd_scatter(x03, out +  6 * ostride, odist, vsize);
    sp_simd_scatter(x05, out + 11 * ostride, odist, vsize);
  }
  {
    sp_simd_t x06, x07, x08;

    a07 = sp_ntt_add_simd(a06, a07, p);

    x06 = a06;
    x07 = sp_ntt_add_simd(a07, a08, p);
    x08 = sp_ntt_sub_simd(a07, a08, p);

    sp_simd_scatter(x08, out +  2 * ostride, odist, vsize);
    sp_simd_scatter(x07, out +  7 * ostride, odist, vsize);
    sp_simd_scatter(x06, out + 12 * ostride, odist, vsize);
  }
  {
    sp_simd_t x09, x10, x11;

    a10 = sp_ntt_add_simd(a09, a10, p);

    x09 = a09;
    x10 = sp_ntt_add_simd(a10, a11, p);
    x11 = sp_ntt_sub_simd(a10, a11, p);

    sp_simd_scatter(x09, out +  3 * ostride, odist, vsize);
    sp_simd_scatter(x11, out +  8 * ostride, odist, vsize);
    sp_simd_scatter(x10, out + 13 * ostride, odist, vsize);
  }
  {
    sp_simd_t x12, x13, x14;

    a13 = sp_ntt_add_simd(a12, a13, p);

    x12 = a12;
    x13 = sp_ntt_add_simd(a13, a14, p);
    x14 = sp_ntt_sub_simd(a13, a14, p);

    sp_simd_scatter(x13, out +  4 * ostride, odist, vsize);
    sp_simd_scatter(x12, out +  9 * ostride, odist, vsize);
    sp_simd_scatter(x14, out + 14 * ostride, odist, vsize);
  }
}

static void
ntt15_run_core_simd_interleaved(
		spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		sp_simd_t p, sp_simd_t * c)
{
  sp_simd_t a00, a01, a02, a03, a04, 
            a05, a06, a07, a08, a09, 
            a10, a11, a12, a13, a14;

  {
    sp_simd_t x00 = sp_simd_load(in + 0 * istride);
    sp_simd_t x01 = sp_simd_load(in + 5 * istride);
    sp_simd_t x02 = sp_simd_load(in +10 * istride);

    a01 = sp_ntt_add_simd0(x01, x02, p);
    a02 = sp_ntt_sub_simd0(x01, x02, p);

    a00 = sp_ntt_add_simd0(x00, a01, p);
  }
  {
    sp_simd_t x03 = sp_simd_load(in + 3 * istride);
    sp_simd_t x04 = sp_simd_load(in + 8 * istride);
    sp_simd_t x05 = sp_simd_load(in +13 * istride);

    a04 = sp_ntt_add_simd0(x04, x05, p);
    a05 = sp_ntt_sub_simd0(x04, x05, p);

    a03 = sp_ntt_add_simd0(x03, a04, p);
  }
  {
    sp_simd_t x08 = sp_simd_load(in + 1 * istride);
    sp_simd_t x06 = sp_simd_load(in + 6 * istride);
    sp_simd_t x07 = sp_simd_load(in +11 * istride);

    a07 = sp_ntt_add_simd0(x07, x08, p);
    a08 = sp_ntt_sub_simd0(x07, x08, p);

    a06 = sp_ntt_add_simd0(x06, a07, p);
  }
  {
    sp_simd_t x11 = sp_simd_load(in + 4 * istride);
    sp_simd_t x09 = sp_simd_load(in + 9 * istride);
    sp_simd_t x10 = sp_simd_load(in +14 * istride);

    a10 = sp_ntt_add_simd0(x10, x11, p);
    a11 = sp_ntt_sub_simd0(x10, x11, p);

    a09 = sp_ntt_add_simd0(x09, a10, p);
  }
  {
    sp_simd_t x13 = sp_simd_load(in + 2 * istride);
    sp_simd_t x14 = sp_simd_load(in + 7 * istride);
    sp_simd_t x12 = sp_simd_load(in +12 * istride);

    a13 = sp_ntt_add_simd0(x13, x14, p);
    a14 = sp_ntt_sub_simd0(x13, x14, p);

    a12 = sp_ntt_add_simd0(x12, a13, p);
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a00;
    b1 = a03;
    b4 = a06;
    b2 = a09;
    b3 = a12;

    c1 = sp_ntt_add_simd0(b1, b3, p);
    c3 = sp_ntt_sub_simd0(b1, b3, p);
    c2 = sp_ntt_add_simd0(b2, b4, p);
    c4 = sp_ntt_sub_simd0(b2, b4, p);

    b1 = sp_ntt_add_simd0(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd0(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd0(c3, c4, p);

    b0 = sp_ntt_add_simd0(b0, b1, p);

    b1 = sp_ntt_mul_simd0(b1, c+2*1, p);
    b2 = sp_ntt_mul_simd0(b2, c+2*2, p);
    b3 = sp_ntt_mul_simd0(b3, c+2*3, p);
    b4 = sp_ntt_mul_simd0(b4, c+2*4, p);
    b5 = sp_ntt_mul_simd0(b5, c+2*5, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a00 = b0;
    a03 = b4;
    a06 = b3;
    a09 = b1;
    a12 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a01;
    b1 = a04;
    b4 = a07;
    b2 = a10;
    b3 = a13;

    c1 = sp_ntt_add_simd0(b1, b3, p);
    c3 = sp_ntt_sub_simd0(b1, b3, p);
    c2 = sp_ntt_add_simd0(b2, b4, p);
    c4 = sp_ntt_sub_simd0(b2, b4, p);

    b1 = sp_ntt_add_simd0(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd0(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd0(c3, c4, p);

    b0 = sp_ntt_add_partial_simd0(b0, b1, p);

    b0 = sp_ntt_mul_simd0(b0, c+2*6, p);
    b1 = sp_ntt_mul_simd0(b1, c+2*7, p);
    b2 = sp_ntt_mul_simd0(b2, c+2*8, p);
    b3 = sp_ntt_mul_simd0(b3, c+2*9, p);
    b4 = sp_ntt_mul_simd0(b4, c+2*10, p);
    b5 = sp_ntt_mul_simd0(b5, c+2*11, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a01 = b0;
    a04 = b4;
    a07 = b3;
    a10 = b1;
    a13 = b2;
  }

  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a02;
    b1 = a05;
    b4 = a08;
    b2 = a11;
    b3 = a14;

    c1 = sp_ntt_add_simd0(b1, b3, p);
    c3 = sp_ntt_sub_simd0(b1, b3, p);
    c2 = sp_ntt_add_simd0(b2, b4, p);
    c4 = sp_ntt_sub_simd0(b2, b4, p);

    b1 = sp_ntt_add_simd0(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd0(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd0(c3, c4, p);

    b0 = sp_ntt_add_partial_simd0(b0, b1, p);

    b0 = sp_ntt_mul_simd0(b0, c+2*12, p);
    b1 = sp_ntt_mul_simd0(b1, c+2*13, p);
    b2 = sp_ntt_mul_simd0(b2, c+2*14, p);
    b3 = sp_ntt_mul_simd0(b3, c+2*15, p);
    b4 = sp_ntt_mul_simd0(b4, c+2*16, p);
    b5 = sp_ntt_mul_simd0(b5, c+2*17, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a02 = b0;
    a05 = b4;
    a08 = b3;
    a11 = b1;
    a14 = b2;
  }
  {
    sp_simd_t x00, x01, x02;

    a01 = sp_ntt_add_simd0(a00, a01, p);

    x00 = a00;
    x01 = sp_ntt_add_simd0(a01, a02, p);
    x02 = sp_ntt_sub_simd0(a01, a02, p);

    sp_simd_store(x00, out +  0 * ostride);
    sp_simd_store(x02, out +  5 * ostride);
    sp_simd_store(x01, out + 10 * ostride);
  }
  {
    sp_simd_t x03, x04, x05;

    a04 = sp_ntt_add_simd0(a03, a04, p);

    x03 = a03;
    x04 = sp_ntt_add_simd0(a04, a05, p);
    x05 = sp_ntt_sub_simd0(a04, a05, p);

    sp_simd_store(x04, out +  1 * ostride);
    sp_simd_store(x03, out +  6 * ostride);
    sp_simd_store(x05, out + 11 * ostride);
  }
  {
    sp_simd_t x06, x07, x08;

    a07 = sp_ntt_add_simd0(a06, a07, p);

    x06 = a06;
    x07 = sp_ntt_add_simd0(a07, a08, p);
    x08 = sp_ntt_sub_simd0(a07, a08, p);

    sp_simd_store(x08, out +  2 * ostride);
    sp_simd_store(x07, out +  7 * ostride);
    sp_simd_store(x06, out + 12 * ostride);
  }
  {
    sp_simd_t x09, x10, x11;

    a10 = sp_ntt_add_simd0(a09, a10, p);

    x09 = a09;
    x10 = sp_ntt_add_simd0(a10, a11, p);
    x11 = sp_ntt_sub_simd0(a10, a11, p);

    sp_simd_store(x09, out +  3 * ostride);
    sp_simd_store(x11, out +  8 * ostride);
    sp_simd_store(x10, out + 13 * ostride);
  }
  {
    sp_simd_t x12, x13, x14;

    a13 = sp_ntt_add_simd0(a12, a13, p);

    x12 = a12;
    x13 = sp_ntt_add_simd0(a13, a14, p);
    x14 = sp_ntt_sub_simd0(a13, a14, p);

    sp_simd_store(x13, out +  4 * ostride);
    sp_simd_store(x12, out +  9 * ostride);
    sp_simd_store(x14, out + 14 * ostride);
  }
}


static void
ntt15_twiddle_run_core_simd(
        spv_t in, spv_size_t istride, spv_size_t idist,
		spv_t out, spv_size_t ostride, spv_size_t odist,
		sp_simd_t *w, sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  sp_simd_t a00, a01, a02, a03, a04, 
            a05, a06, a07, a08, a09, 
            a10, a11, a12, a13, a14;

  {
    sp_simd_t x00 = sp_simd_gather(in +  0 * istride, idist, vsize);
    sp_simd_t x01 = sp_simd_gather(in +  5 * istride, idist, vsize);
    sp_simd_t x02 = sp_simd_gather(in + 10 * istride, idist, vsize);

    a01 = sp_ntt_add_simd(x01, x02, p);
    a02 = sp_ntt_sub_simd(x01, x02, p);

    a00 = sp_ntt_add_simd(x00, a01, p);
  }
  {
    sp_simd_t x03 = sp_simd_gather(in +  3 * istride, idist, vsize);
    sp_simd_t x04 = sp_simd_gather(in +  8 * istride, idist, vsize);
    sp_simd_t x05 = sp_simd_gather(in + 13 * istride, idist, vsize);

    a04 = sp_ntt_add_simd(x04, x05, p);
    a05 = sp_ntt_sub_simd(x04, x05, p);

    a03 = sp_ntt_add_simd(x03, a04, p);
  }
  {
    sp_simd_t x08 = sp_simd_gather(in +  1 * istride, idist, vsize);
    sp_simd_t x06 = sp_simd_gather(in +  6 * istride, idist, vsize);
    sp_simd_t x07 = sp_simd_gather(in + 11 * istride, idist, vsize);

    a07 = sp_ntt_add_simd(x07, x08, p);
    a08 = sp_ntt_sub_simd(x07, x08, p);

    a06 = sp_ntt_add_simd(x06, a07, p);
  }
  {
    sp_simd_t x11 = sp_simd_gather(in +  4 * istride, idist, vsize);
    sp_simd_t x09 = sp_simd_gather(in +  9 * istride, idist, vsize);
    sp_simd_t x10 = sp_simd_gather(in + 14 * istride, idist, vsize);

    a10 = sp_ntt_add_simd(x10, x11, p);
    a11 = sp_ntt_sub_simd(x10, x11, p);

    a09 = sp_ntt_add_simd(x09, a10, p);
  }
  {
    sp_simd_t x13 = sp_simd_gather(in +  2 * istride, idist, vsize);
    sp_simd_t x14 = sp_simd_gather(in +  7 * istride, idist, vsize);
    sp_simd_t x12 = sp_simd_gather(in + 12 * istride, idist, vsize);

    a13 = sp_ntt_add_simd(x13, x14, p);
    a14 = sp_ntt_sub_simd(x13, x14, p);

    a12 = sp_ntt_add_simd(x12, a13, p);
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a00;
    b1 = a03;
    b4 = a06;
    b2 = a09;
    b3 = a12;

    c1 = sp_ntt_add_simd(b1, b3, p);
    c3 = sp_ntt_sub_simd(b1, b3, p);
    c2 = sp_ntt_add_simd(b2, b4, p);
    c4 = sp_ntt_sub_simd(b2, b4, p);

    b1 = sp_ntt_add_simd(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd(c3, c4, p);

    b0 = sp_ntt_add_simd(b0, b1, p);

    b1 = sp_ntt_mul_simd(b1, ntt_const[1], ntt_const[NC+1], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[2], ntt_const[NC+2], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[3], ntt_const[NC+3], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[4], ntt_const[NC+4], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[5], ntt_const[NC+5], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a00 = b0;
    a03 = b4;
    a06 = b3;
    a09 = b1;
    a12 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a01;
    b1 = a04;
    b4 = a07;
    b2 = a10;
    b3 = a13;

    c1 = sp_ntt_add_simd(b1, b3, p);
    c3 = sp_ntt_sub_simd(b1, b3, p);
    c2 = sp_ntt_add_simd(b2, b4, p);
    c4 = sp_ntt_sub_simd(b2, b4, p);

    b1 = sp_ntt_add_simd(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd(c3, c4, p);

    b0 = sp_ntt_add_partial_simd(b0, b1, p);

    b0 = sp_ntt_mul_simd(b0, ntt_const[6], ntt_const[NC+6], p);
    b1 = sp_ntt_mul_simd(b1, ntt_const[7], ntt_const[NC+7], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[8], ntt_const[NC+8], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[9], ntt_const[NC+9], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[10], ntt_const[NC+10], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[11], ntt_const[NC+11], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a01 = b0;
    a04 = b4;
    a07 = b3;
    a10 = b1;
    a13 = b2;
  }

  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a02;
    b1 = a05;
    b4 = a08;
    b2 = a11;
    b3 = a14;

    c1 = sp_ntt_add_simd(b1, b3, p);
    c3 = sp_ntt_sub_simd(b1, b3, p);
    c2 = sp_ntt_add_simd(b2, b4, p);
    c4 = sp_ntt_sub_simd(b2, b4, p);

    b1 = sp_ntt_add_simd(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd(c3, c4, p);

    b0 = sp_ntt_add_partial_simd(b0, b1, p);

    b0 = sp_ntt_mul_simd(b0, ntt_const[12], ntt_const[NC+12], p);
    b1 = sp_ntt_mul_simd(b1, ntt_const[13], ntt_const[NC+13], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[14], ntt_const[NC+14], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[15], ntt_const[NC+15], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[16], ntt_const[NC+16], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[17], ntt_const[NC+17], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a02 = b0;
    a05 = b4;
    a08 = b3;
    a11 = b1;
    a14 = b2;
  }
  {
    sp_simd_t x00, x01, x02;

    a01 = sp_ntt_add_simd(a00, a01, p);

    x00 = a00;
    x01 = sp_ntt_add_partial_simd(a01, a02, p);
    x02 = sp_ntt_sub_partial_simd(a01, a02, p);

    x02 = sp_ntt_twiddle_mul_simd(x02, w + 8, p);
    x01 = sp_ntt_twiddle_mul_simd(x01, w + 18, p);

    sp_simd_scatter(x00, out +  0 * ostride, odist, vsize);
    sp_simd_scatter(x02, out +  5 * ostride, odist, vsize);
    sp_simd_scatter(x01, out + 10 * ostride, odist, vsize);
  }
  {
    sp_simd_t x03, x04, x05;

    a04 = sp_ntt_add_simd(a03, a04, p);

    x03 = a03;
    x04 = sp_ntt_add_partial_simd(a04, a05, p);
    x05 = sp_ntt_sub_partial_simd(a04, a05, p);

    x04 = sp_ntt_twiddle_mul_simd(x04, w + 0, p);
    x03 = sp_ntt_twiddle_mul_simd(x03, w + 10, p);
    x05 = sp_ntt_twiddle_mul_simd(x05, w + 20, p);

    sp_simd_scatter(x04, out +  1 * ostride, odist, vsize);
    sp_simd_scatter(x03, out +  6 * ostride, odist, vsize);
    sp_simd_scatter(x05, out + 11 * ostride, odist, vsize);
  }
  {
    sp_simd_t x06, x07, x08;

    a07 = sp_ntt_add_simd(a06, a07, p);

    x06 = a06;
    x07 = sp_ntt_add_partial_simd(a07, a08, p);
    x08 = sp_ntt_sub_partial_simd(a07, a08, p);

    x08 = sp_ntt_twiddle_mul_simd(x08, w + 2, p);
    x07 = sp_ntt_twiddle_mul_simd(x07, w + 12, p);
    x06 = sp_ntt_twiddle_mul_simd(x06, w + 22, p);

    sp_simd_scatter(x08, out +  2 * ostride, odist, vsize);
    sp_simd_scatter(x07, out +  7 * ostride, odist, vsize);
    sp_simd_scatter(x06, out + 12 * ostride, odist, vsize);
  }
  {
    sp_simd_t x09, x10, x11;

    a10 = sp_ntt_add_simd(a09, a10, p);

    x09 = a09;
    x10 = sp_ntt_add_partial_simd(a10, a11, p);
    x11 = sp_ntt_sub_partial_simd(a10, a11, p);

    x09 = sp_ntt_twiddle_mul_simd(x09, w + 4, p);
    x11 = sp_ntt_twiddle_mul_simd(x11, w + 14, p);
    x10 = sp_ntt_twiddle_mul_simd(x10, w + 24, p);

    sp_simd_scatter(x09, out +  3 * ostride, odist, vsize);
    sp_simd_scatter(x11, out +  8 * ostride, odist, vsize);
    sp_simd_scatter(x10, out + 13 * ostride, odist, vsize);
  }
  {
    sp_simd_t x12, x13, x14;

    a13 = sp_ntt_add_simd(a12, a13, p);

    x12 = a12;
    x13 = sp_ntt_add_partial_simd(a13, a14, p);
    x14 = sp_ntt_sub_partial_simd(a13, a14, p);

    x13 = sp_ntt_twiddle_mul_simd(x13, w + 6, p);
    x12 = sp_ntt_twiddle_mul_simd(x12, w + 16, p);
    x14 = sp_ntt_twiddle_mul_simd(x14, w + 26, p);

    sp_simd_scatter(x13, out +  4 * ostride, odist, vsize);
    sp_simd_scatter(x12, out +  9 * ostride, odist, vsize);
    sp_simd_scatter(x14, out + 14 * ostride, odist, vsize);
  }
}

static void
ntt15_twiddle_run_core_simd_interleaved(
		spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		sp_simd_t *w, sp_simd_t p, sp_simd_t * c)
{
  sp_simd_t a00, a01, a02, a03, a04, 
            a05, a06, a07, a08, a09, 
            a10, a11, a12, a13, a14;

  {
    sp_simd_t x00 = sp_simd_load(in +  0 * istride);
    sp_simd_t x01 = sp_simd_load(in +  5 * istride);
    sp_simd_t x02 = sp_simd_load(in + 10 * istride);

    a01 = sp_ntt_add_simd0(x01, x02, p);
    a02 = sp_ntt_sub_simd0(x01, x02, p);

    a00 = sp_ntt_add_simd0(x00, a01, p);
  }
  {
    sp_simd_t x03 = sp_simd_load(in +  3 * istride);
    sp_simd_t x04 = sp_simd_load(in +  8 * istride);
    sp_simd_t x05 = sp_simd_load(in + 13 * istride);

    a04 = sp_ntt_add_simd0(x04, x05, p);
    a05 = sp_ntt_sub_simd0(x04, x05, p);

    a03 = sp_ntt_add_simd0(x03, a04, p);
  }
  {
    sp_simd_t x08 = sp_simd_load(in +  1 * istride);
    sp_simd_t x06 = sp_simd_load(in +  6 * istride);
    sp_simd_t x07 = sp_simd_load(in + 11 * istride);

    a07 = sp_ntt_add_simd0(x07, x08, p);
    a08 = sp_ntt_sub_simd0(x07, x08, p);

    a06 = sp_ntt_add_simd0(x06, a07, p);
  }
  {
    sp_simd_t x11 = sp_simd_load(in +  4 * istride);
    sp_simd_t x09 = sp_simd_load(in +  9 * istride);
    sp_simd_t x10 = sp_simd_load(in + 14 * istride);

    a10 = sp_ntt_add_simd0(x10, x11, p);
    a11 = sp_ntt_sub_simd0(x10, x11, p);

    a09 = sp_ntt_add_simd0(x09, a10, p);
  }
  {
    sp_simd_t x13 = sp_simd_load(in +  2 * istride);
    sp_simd_t x14 = sp_simd_load(in +  7 * istride);
    sp_simd_t x12 = sp_simd_load(in + 12 * istride);

    a13 = sp_ntt_add_simd0(x13, x14, p);
    a14 = sp_ntt_sub_simd0(x13, x14, p);

    a12 = sp_ntt_add_simd0(x12, a13, p);
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a00;
    b1 = a03;
    b4 = a06;
    b2 = a09;
    b3 = a12;

    c1 = sp_ntt_add_simd0(b1, b3, p);
    c3 = sp_ntt_sub_simd0(b1, b3, p);
    c2 = sp_ntt_add_simd0(b2, b4, p);
    c4 = sp_ntt_sub_simd0(b2, b4, p);

    b1 = sp_ntt_add_simd0(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd0(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd0(c3, c4, p);

    b0 = sp_ntt_add_simd0(b0, b1, p);

    b1 = sp_ntt_mul_simd0(b1, c+2*1, p);
    b2 = sp_ntt_mul_simd0(b2, c+2*2, p);
    b3 = sp_ntt_mul_simd0(b3, c+2*3, p);
    b4 = sp_ntt_mul_simd0(b4, c+2*4, p);
    b5 = sp_ntt_mul_simd0(b5, c+2*5, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a00 = b0;
    a03 = b4;
    a06 = b3;
    a09 = b1;
    a12 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a01;
    b1 = a04;
    b4 = a07;
    b2 = a10;
    b3 = a13;

    c1 = sp_ntt_add_simd0(b1, b3, p);
    c3 = sp_ntt_sub_simd0(b1, b3, p);
    c2 = sp_ntt_add_simd0(b2, b4, p);
    c4 = sp_ntt_sub_simd0(b2, b4, p);

    b1 = sp_ntt_add_simd0(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd0(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd0(c3, c4, p);

    b0 = sp_ntt_add_partial_simd0(b0, b1, p);

    b0 = sp_ntt_mul_simd0(b0, c+2*6, p);
    b1 = sp_ntt_mul_simd0(b1, c+2*7, p);
    b2 = sp_ntt_mul_simd0(b2, c+2*8, p);
    b3 = sp_ntt_mul_simd0(b3, c+2*9, p);
    b4 = sp_ntt_mul_simd0(b4, c+2*10, p);
    b5 = sp_ntt_mul_simd0(b5, c+2*11, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a01 = b0;
    a04 = b4;
    a07 = b3;
    a10 = b1;
    a13 = b2;
  }

  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a02;
    b1 = a05;
    b4 = a08;
    b2 = a11;
    b3 = a14;

    c1 = sp_ntt_add_simd0(b1, b3, p);
    c3 = sp_ntt_sub_simd0(b1, b3, p);
    c2 = sp_ntt_add_simd0(b2, b4, p);
    c4 = sp_ntt_sub_simd0(b2, b4, p);

    b1 = sp_ntt_add_simd0(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd0(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd0(c3, c4, p);

    b0 = sp_ntt_add_partial_simd0(b0, b1, p);

    b0 = sp_ntt_mul_simd0(b0, c+2*12, p);
    b1 = sp_ntt_mul_simd0(b1, c+2*13, p);
    b2 = sp_ntt_mul_simd0(b2, c+2*14, p);
    b3 = sp_ntt_mul_simd0(b3, c+2*15, p);
    b4 = sp_ntt_mul_simd0(b4, c+2*16, p);
    b5 = sp_ntt_mul_simd0(b5, c+2*17, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a02 = b0;
    a05 = b4;
    a08 = b3;
    a11 = b1;
    a14 = b2;
  }
  {
    sp_simd_t x00, x01, x02;

    a01 = sp_ntt_add_simd0(a00, a01, p);

    x00 = a00;
    x01 = sp_ntt_add_partial_simd0(a01, a02, p);
    x02 = sp_ntt_sub_partial_simd0(a01, a02, p);

    x02 = sp_ntt_twiddle_mul_simd0(x02, w+8, p);
    x01 = sp_ntt_twiddle_mul_simd0(x01, w+18, p);

    sp_simd_store(x00, out +  0 * ostride);
    sp_simd_store(x02, out +  5 * ostride);
    sp_simd_store(x01, out + 10 * ostride);
  }
  {
    sp_simd_t x03, x04, x05;

    a04 = sp_ntt_add_simd0(a03, a04, p);

    x03 = a03;
    x04 = sp_ntt_add_partial_simd0(a04, a05, p);
    x05 = sp_ntt_sub_partial_simd0(a04, a05, p);

    x04 = sp_ntt_twiddle_mul_simd0(x04, w+0, p);
    x03 = sp_ntt_twiddle_mul_simd0(x03, w+10, p);
    x05 = sp_ntt_twiddle_mul_simd0(x05, w+20, p);

    sp_simd_store(x04, out +  1 * ostride);
    sp_simd_store(x03, out +  6 * ostride);
    sp_simd_store(x05, out + 11 * ostride);
  }
  {
    sp_simd_t x06, x07, x08;

    a07 = sp_ntt_add_simd0(a06, a07, p);

    x06 = a06;
    x07 = sp_ntt_add_partial_simd0(a07, a08, p);
    x08 = sp_ntt_sub_partial_simd0(a07, a08, p);

    x08 = sp_ntt_twiddle_mul_simd0(x08, w+2, p);
    x07 = sp_ntt_twiddle_mul_simd0(x07, w+12, p);
    x06 = sp_ntt_twiddle_mul_simd0(x06, w+22, p);

    sp_simd_store(x08, out +  2 * ostride);
    sp_simd_store(x07, out +  7 * ostride);
    sp_simd_store(x06, out + 12 * ostride);
  }
  {
    sp_simd_t x09, x10, x11;

    a10 = sp_ntt_add_simd0(a09, a10, p);

    x09 = a09;
    x10 = sp_ntt_add_partial_simd0(a10, a11, p);
    x11 = sp_ntt_sub_partial_simd0(a10, a11, p);

    x09 = sp_ntt_twiddle_mul_simd0(x09, w+4, p);
    x11 = sp_ntt_twiddle_mul_simd0(x11, w+14, p);
    x10 = sp_ntt_twiddle_mul_simd0(x10, w+24, p);

    sp_simd_store(x09, out +  3 * ostride);
    sp_simd_store(x11, out +  8 * ostride);
    sp_simd_store(x10, out + 13 * ostride);
  }
  {
    sp_simd_t x12, x13, x14;

    a13 = sp_ntt_add_simd0(a12, a13, p);

    x12 = a12;
    x13 = sp_ntt_add_partial_simd0(a13, a14, p);
    x14 = sp_ntt_sub_partial_simd0(a13, a14, p);

    x13 = sp_ntt_twiddle_mul_simd0(x13, w+6, p);
    x12 = sp_ntt_twiddle_mul_simd0(x12, w+16, p);
    x14 = sp_ntt_twiddle_mul_simd0(x14, w+26, p);

    sp_simd_store(x13, out +  4 * ostride);
    sp_simd_store(x12, out +  9 * ostride);
    sp_simd_store(x14, out + 14 * ostride);
  }
}


static void
ntt15_pfa_run_core_simd(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t inc2, spv_size_t n,
	  sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  spv_size_t j00, j01, j02, j03, j04,
             j05, j06, j07, j08, j09,
	     j10, j11, j12, j13, j14;

  sp_simd_t a00, a01, a02, a03, a04, 
            a05, a06, a07, a08, a09, 
            a10, a11, a12, a13, a14;

  j00 = start;
  j01 = sp_array_inc(j00, inc, n);
  j02 = sp_array_inc(j00, 2 * inc, n);
  j03 = sp_array_inc(j00, 3 * inc, n);
  j04 = sp_array_inc(j00, 4 * inc, n);
  j05 = sp_array_inc(j00, 5 * inc, n);
  j06 = sp_array_inc(j00, 6 * inc, n);
  j07 = sp_array_inc(j00, 7 * inc, n);
  j08 = sp_array_inc(j00, 8 * inc, n);
  j09 = sp_array_inc(j00, 9 * inc, n);
  j10 = sp_array_inc(j00, 10 * inc, n);
  j11 = sp_array_inc(j00, 11 * inc, n);
  j12 = sp_array_inc(j00, 12 * inc, n);
  j13 = sp_array_inc(j00, 13 * inc, n);
  j14 = sp_array_inc(j00, 14 * inc, n);

  {
    sp_simd_t x00 = sp_simd_pfa_gather(x, j00, inc2, n, vsize);
    sp_simd_t x01 = sp_simd_pfa_gather(x, j05, inc2, n, vsize);
    sp_simd_t x02 = sp_simd_pfa_gather(x, j10, inc2, n, vsize);

    a01 = sp_ntt_add_simd(x01, x02, p);
    a02 = sp_ntt_sub_simd(x01, x02, p);

    a00 = sp_ntt_add_simd(x00, a01, p);
  }
  {
    sp_simd_t x03 = sp_simd_pfa_gather(x, j03, inc2, n, vsize);
    sp_simd_t x04 = sp_simd_pfa_gather(x, j08, inc2, n, vsize);
    sp_simd_t x05 = sp_simd_pfa_gather(x, j13, inc2, n, vsize);

    a04 = sp_ntt_add_simd(x04, x05, p);
    a05 = sp_ntt_sub_simd(x04, x05, p);

    a03 = sp_ntt_add_simd(x03, a04, p);
  }
  {
    sp_simd_t x08 = sp_simd_pfa_gather(x, j01, inc2, n, vsize);
    sp_simd_t x06 = sp_simd_pfa_gather(x, j06, inc2, n, vsize);
    sp_simd_t x07 = sp_simd_pfa_gather(x, j11, inc2, n, vsize);

    a07 = sp_ntt_add_simd(x07, x08, p);
    a08 = sp_ntt_sub_simd(x07, x08, p);

    a06 = sp_ntt_add_simd(x06, a07, p);
  }
  {
    sp_simd_t x11 = sp_simd_pfa_gather(x, j04, inc2, n, vsize);
    sp_simd_t x09 = sp_simd_pfa_gather(x, j09, inc2, n, vsize);
    sp_simd_t x10 = sp_simd_pfa_gather(x, j14, inc2, n, vsize);

    a10 = sp_ntt_add_simd(x10, x11, p);
    a11 = sp_ntt_sub_simd(x10, x11, p);

    a09 = sp_ntt_add_simd(x09, a10, p);
  }
  {
    sp_simd_t x13 = sp_simd_pfa_gather(x, j02, inc2, n, vsize);
    sp_simd_t x14 = sp_simd_pfa_gather(x, j07, inc2, n, vsize);
    sp_simd_t x12 = sp_simd_pfa_gather(x, j12, inc2, n, vsize);

    a13 = sp_ntt_add_simd(x13, x14, p);
    a14 = sp_ntt_sub_simd(x13, x14, p);

    a12 = sp_ntt_add_simd(x12, a13, p);
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a00;
    b1 = a03;
    b4 = a06;
    b2 = a09;
    b3 = a12;

    c1 = sp_ntt_add_simd(b1, b3, p);
    c3 = sp_ntt_sub_simd(b1, b3, p);
    c2 = sp_ntt_add_simd(b2, b4, p);
    c4 = sp_ntt_sub_simd(b2, b4, p);

    b1 = sp_ntt_add_simd(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd(c3, c4, p);

    b0 = sp_ntt_add_simd(b0, b1, p);

    b1 = sp_ntt_mul_simd(b1, ntt_const[1], ntt_const[NC+1], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[2], ntt_const[NC+2], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[3], ntt_const[NC+3], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[4], ntt_const[NC+4], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[5], ntt_const[NC+5], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a00 = b0;
    a03 = b4;
    a06 = b3;
    a09 = b1;
    a12 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a01;
    b1 = a04;
    b4 = a07;
    b2 = a10;
    b3 = a13;

    c1 = sp_ntt_add_simd(b1, b3, p);
    c3 = sp_ntt_sub_simd(b1, b3, p);
    c2 = sp_ntt_add_simd(b2, b4, p);
    c4 = sp_ntt_sub_simd(b2, b4, p);

    b1 = sp_ntt_add_simd(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd(c3, c4, p);

    b0 = sp_ntt_add_partial_simd(b0, b1, p);

    b0 = sp_ntt_mul_simd(b0, ntt_const[6], ntt_const[NC+6], p);
    b1 = sp_ntt_mul_simd(b1, ntt_const[7], ntt_const[NC+7], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[8], ntt_const[NC+8], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[9], ntt_const[NC+9], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[10], ntt_const[NC+10], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[11], ntt_const[NC+11], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a01 = b0;
    a04 = b4;
    a07 = b3;
    a10 = b1;
    a13 = b2;
  }

  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a02;
    b1 = a05;
    b4 = a08;
    b2 = a11;
    b3 = a14;

    c1 = sp_ntt_add_simd(b1, b3, p);
    c3 = sp_ntt_sub_simd(b1, b3, p);
    c2 = sp_ntt_add_simd(b2, b4, p);
    c4 = sp_ntt_sub_simd(b2, b4, p);

    b1 = sp_ntt_add_simd(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd(c3, c4, p);

    b0 = sp_ntt_add_partial_simd(b0, b1, p);

    b0 = sp_ntt_mul_simd(b0, ntt_const[12], ntt_const[NC+12], p);
    b1 = sp_ntt_mul_simd(b1, ntt_const[13], ntt_const[NC+13], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[14], ntt_const[NC+14], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[15], ntt_const[NC+15], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[16], ntt_const[NC+16], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[17], ntt_const[NC+17], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a02 = b0;
    a05 = b4;
    a08 = b3;
    a11 = b1;
    a14 = b2;
  }
  {
    sp_simd_t x00, x01, x02;

    a01 = sp_ntt_add_simd(a00, a01, p);

    x00 = a00;
    x01 = sp_ntt_add_simd(a01, a02, p);
    x02 = sp_ntt_sub_simd(a01, a02, p);

    sp_simd_pfa_scatter(x00, x, j00, inc2, n, vsize);
    sp_simd_pfa_scatter(x02, x, j05, inc2, n, vsize);
    sp_simd_pfa_scatter(x01, x, j10, inc2, n, vsize);
  }
  {
    sp_simd_t x03, x04, x05;

    a04 = sp_ntt_add_simd(a03, a04, p);

    x03 = a03;
    x04 = sp_ntt_add_simd(a04, a05, p);
    x05 = sp_ntt_sub_simd(a04, a05, p);

    sp_simd_pfa_scatter(x04, x, j01, inc2, n, vsize);
    sp_simd_pfa_scatter(x03, x, j06, inc2, n, vsize);
    sp_simd_pfa_scatter(x05, x, j11, inc2, n, vsize);
  }
  {
    sp_simd_t x06, x07, x08;

    a07 = sp_ntt_add_simd(a06, a07, p);

    x06 = a06;
    x07 = sp_ntt_add_simd(a07, a08, p);
    x08 = sp_ntt_sub_simd(a07, a08, p);

    sp_simd_pfa_scatter(x08, x, j02, inc2, n, vsize);
    sp_simd_pfa_scatter(x07, x, j07, inc2, n, vsize);
    sp_simd_pfa_scatter(x06, x, j12, inc2, n, vsize);
  }
  {
    sp_simd_t x09, x10, x11;

    a10 = sp_ntt_add_simd(a09, a10, p);

    x09 = a09;
    x10 = sp_ntt_add_simd(a10, a11, p);
    x11 = sp_ntt_sub_simd(a10, a11, p);

    sp_simd_pfa_scatter(x09, x, j03, inc2, n, vsize);
    sp_simd_pfa_scatter(x11, x, j08, inc2, n, vsize);
    sp_simd_pfa_scatter(x10, x, j13, inc2, n, vsize);
  }
  {
    sp_simd_t x12, x13, x14;

    a13 = sp_ntt_add_simd(a12, a13, p);

    x12 = a12;
    x13 = sp_ntt_add_simd(a13, a14, p);
    x14 = sp_ntt_sub_simd(a13, a14, p);

    sp_simd_pfa_scatter(x13, x, j04, inc2, n, vsize);
    sp_simd_pfa_scatter(x12, x, j09, inc2, n, vsize);
    sp_simd_pfa_scatter(x14, x, j14, inc2, n, vsize);
  }
}

static void
ntt15_pfa_run_core_simd_interleaved(
	  spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t n,
	  sp_simd_t p, sp_simd_t * c)
{
  spv_size_t j00, j01, j02, j03, j04,
             j05, j06, j07, j08, j09,
	     j10, j11, j12, j13, j14;

  sp_simd_t a00, a01, a02, a03, a04, 
            a05, a06, a07, a08, a09, 
            a10, a11, a12, a13, a14;

  j00 = start;
  j01 = sp_array_inc(j00, inc, n);
  j02 = sp_array_inc(j00, 2 * inc, n);
  j03 = sp_array_inc(j00, 3 * inc, n);
  j04 = sp_array_inc(j00, 4 * inc, n);
  j05 = sp_array_inc(j00, 5 * inc, n);
  j06 = sp_array_inc(j00, 6 * inc, n);
  j07 = sp_array_inc(j00, 7 * inc, n);
  j08 = sp_array_inc(j00, 8 * inc, n);
  j09 = sp_array_inc(j00, 9 * inc, n);
  j10 = sp_array_inc(j00, 10 * inc, n);
  j11 = sp_array_inc(j00, 11 * inc, n);
  j12 = sp_array_inc(j00, 12 * inc, n);
  j13 = sp_array_inc(j00, 13 * inc, n);
  j14 = sp_array_inc(j00, 14 * inc, n);

  {
    sp_simd_t x00 = sp_simd_load(x + j00);
    sp_simd_t x01 = sp_simd_load(x + j05);
    sp_simd_t x02 = sp_simd_load(x + j10);

    a01 = sp_ntt_add_simd0(x01, x02, p);
    a02 = sp_ntt_sub_simd0(x01, x02, p);

    a00 = sp_ntt_add_simd0(x00, a01, p);
  }
  {
    sp_simd_t x03 = sp_simd_load(x + j03);
    sp_simd_t x04 = sp_simd_load(x + j08);
    sp_simd_t x05 = sp_simd_load(x + j13);

    a04 = sp_ntt_add_simd0(x04, x05, p);
    a05 = sp_ntt_sub_simd0(x04, x05, p);

    a03 = sp_ntt_add_simd0(x03, a04, p);
  }
  {
    sp_simd_t x08 = sp_simd_load(x + j01);
    sp_simd_t x06 = sp_simd_load(x + j06);
    sp_simd_t x07 = sp_simd_load(x + j11);

    a07 = sp_ntt_add_simd0(x07, x08, p);
    a08 = sp_ntt_sub_simd0(x07, x08, p);

    a06 = sp_ntt_add_simd0(x06, a07, p);
  }
  {
    sp_simd_t x11 = sp_simd_load(x + j04);
    sp_simd_t x09 = sp_simd_load(x + j09);
    sp_simd_t x10 = sp_simd_load(x + j14);

    a10 = sp_ntt_add_simd0(x10, x11, p);
    a11 = sp_ntt_sub_simd0(x10, x11, p);

    a09 = sp_ntt_add_simd0(x09, a10, p);
  }
  {
    sp_simd_t x13 = sp_simd_load(x + j02);
    sp_simd_t x14 = sp_simd_load(x + j07);
    sp_simd_t x12 = sp_simd_load(x + j12);

    a13 = sp_ntt_add_simd0(x13, x14, p);
    a14 = sp_ntt_sub_simd0(x13, x14, p);

    a12 = sp_ntt_add_simd0(x12, a13, p);
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a00;
    b1 = a03;
    b4 = a06;
    b2 = a09;
    b3 = a12;

    c1 = sp_ntt_add_simd0(b1, b3, p);
    c3 = sp_ntt_sub_simd0(b1, b3, p);
    c2 = sp_ntt_add_simd0(b2, b4, p);
    c4 = sp_ntt_sub_simd0(b2, b4, p);

    b1 = sp_ntt_add_simd0(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd0(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd0(c3, c4, p);

    b0 = sp_ntt_add_simd0(b0, b1, p);

    b1 = sp_ntt_mul_simd0(b1, c+2*1, p);
    b2 = sp_ntt_mul_simd0(b2, c+2*2, p);
    b3 = sp_ntt_mul_simd0(b3, c+2*3, p);
    b4 = sp_ntt_mul_simd0(b4, c+2*4, p);
    b5 = sp_ntt_mul_simd0(b5, c+2*5, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a00 = b0;
    a03 = b4;
    a06 = b3;
    a09 = b1;
    a12 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a01;
    b1 = a04;
    b4 = a07;
    b2 = a10;
    b3 = a13;

    c1 = sp_ntt_add_simd0(b1, b3, p);
    c3 = sp_ntt_sub_simd0(b1, b3, p);
    c2 = sp_ntt_add_simd0(b2, b4, p);
    c4 = sp_ntt_sub_simd0(b2, b4, p);

    b1 = sp_ntt_add_simd0(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd0(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd0(c3, c4, p);

    b0 = sp_ntt_add_partial_simd0(b0, b1, p);

    b0 = sp_ntt_mul_simd0(b0, c+2*6, p);
    b1 = sp_ntt_mul_simd0(b1, c+2*7, p);
    b2 = sp_ntt_mul_simd0(b2, c+2*8, p);
    b3 = sp_ntt_mul_simd0(b3, c+2*9, p);
    b4 = sp_ntt_mul_simd0(b4, c+2*10, p);
    b5 = sp_ntt_mul_simd0(b5, c+2*11, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a01 = b0;
    a04 = b4;
    a07 = b3;
    a10 = b1;
    a13 = b2;
  }

  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a02;
    b1 = a05;
    b4 = a08;
    b2 = a11;
    b3 = a14;

    c1 = sp_ntt_add_simd0(b1, b3, p);
    c3 = sp_ntt_sub_simd0(b1, b3, p);
    c2 = sp_ntt_add_simd0(b2, b4, p);
    c4 = sp_ntt_sub_simd0(b2, b4, p);

    b1 = sp_ntt_add_simd0(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd0(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd0(c3, c4, p);

    b0 = sp_ntt_add_partial_simd0(b0, b1, p);

    b0 = sp_ntt_mul_simd0(b0, c+2*12, p);
    b1 = sp_ntt_mul_simd0(b1, c+2*13, p);
    b2 = sp_ntt_mul_simd0(b2, c+2*14, p);
    b3 = sp_ntt_mul_simd0(b3, c+2*15, p);
    b4 = sp_ntt_mul_simd0(b4, c+2*16, p);
    b5 = sp_ntt_mul_simd0(b5, c+2*17, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a02 = b0;
    a05 = b4;
    a08 = b3;
    a11 = b1;
    a14 = b2;
  }
  {
    sp_simd_t x00, x01, x02;

    a01 = sp_ntt_add_simd0(a00, a01, p);

    x00 = a00;
    x01 = sp_ntt_add_simd0(a01, a02, p);
    x02 = sp_ntt_sub_simd0(a01, a02, p);

    sp_simd_store(x00, x + j00);
    sp_simd_store(x02, x + j05);
    sp_simd_store(x01, x + j10);
  }
  {
    sp_simd_t x03, x04, x05;

    a04 = sp_ntt_add_simd0(a03, a04, p);

    x03 = a03;
    x04 = sp_ntt_add_simd0(a04, a05, p);
    x05 = sp_ntt_sub_simd0(a04, a05, p);

    sp_simd_store(x04, x + j01);
    sp_simd_store(x03, x + j06);
    sp_simd_store(x05, x + j11);
  }
  {
    sp_simd_t x06, x07, x08;

    a07 = sp_ntt_add_simd0(a06, a07, p);

    x06 = a06;
    x07 = sp_ntt_add_simd0(a07, a08, p);
    x08 = sp_ntt_sub_simd0(a07, a08, p);

    sp_simd_store(x08, x + j02);
    sp_simd_store(x07, x + j07);
    sp_simd_store(x06, x + j12);
  }
  {
    sp_simd_t x09, x10, x11;

    a10 = sp_ntt_add_simd0(a09, a10, p);

    x09 = a09;
    x10 = sp_ntt_add_simd0(a10, a11, p);
    x11 = sp_ntt_sub_simd0(a10, a11, p);

    sp_simd_store(x09, x + j03);
    sp_simd_store(x11, x + j08);
    sp_simd_store(x10, x + j13);
  }
  {
    sp_simd_t x12, x13, x14;

    a13 = sp_ntt_add_simd0(a12, a13, p);

    x12 = a12;
    x13 = sp_ntt_add_simd0(a13, a14, p);
    x14 = sp_ntt_sub_simd0(a13, a14, p);

    sp_simd_store(x13, x + j04);
    sp_simd_store(x12, x + j09);
    sp_simd_store(x14, x + j14);
  }
}

DECLARE_CORE_ROUTINES_SIMD(15)
