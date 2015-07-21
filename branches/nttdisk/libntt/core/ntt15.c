#include "ntt-impl.h"

#define NC 18

static const uint8_t fixed_const[NC] = {1};

static const uint8_t *
ntt15_get_fixed_ntt_const(void)
{
  return fixed_const;
}

void
ntt15_init(spv_t out, sp_t p, sp_t d, 
	  sp_t primroot, sp_t order)
{
  nttdata_init_generic(&ntt15_config, out, p, d, primroot, order);
}

static void 
ntt15_run_core(spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		sp_t p, spv_t ntt_const)
{
  sp_t a00, a01, a02, a03, a04, 
       a05, a06, a07, a08, a09, 
       a10, a11, a12, a13, a14;

  {
    sp_t x00 = in[ 0 * istride];
    sp_t x01 = in[ 5 * istride];
    sp_t x02 = in[10 * istride];

    a01 = sp_ntt_add(x01, x02, p);
    a02 = sp_ntt_sub(x01, x02, p);

    a00 = sp_ntt_add(x00, a01, p);
  }
  {
    sp_t x03 = in[ 3 * istride];
    sp_t x04 = in[ 8 * istride];
    sp_t x05 = in[13 * istride];

    a04 = sp_ntt_add(x04, x05, p);
    a05 = sp_ntt_sub(x04, x05, p);

    a03 = sp_ntt_add(x03, a04, p);
  }
  {
    sp_t x08 = in[ 1 * istride];
    sp_t x06 = in[ 6 * istride];
    sp_t x07 = in[11 * istride];

    a07 = sp_ntt_add(x07, x08, p);
    a08 = sp_ntt_sub(x07, x08, p);

    a06 = sp_ntt_add(x06, a07, p);
  }
  {
    sp_t x11 = in[ 4 * istride];
    sp_t x09 = in[ 9 * istride];
    sp_t x10 = in[14 * istride];

    a10 = sp_ntt_add(x10, x11, p);
    a11 = sp_ntt_sub(x10, x11, p);

    a09 = sp_ntt_add(x09, a10, p);
  }
  {
    sp_t x13 = in[ 2 * istride];
    sp_t x14 = in[ 7 * istride];
    sp_t x12 = in[12 * istride];

    a13 = sp_ntt_add(x13, x14, p);
    a14 = sp_ntt_sub(x13, x14, p);

    a12 = sp_ntt_add(x12, a13, p);
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a00;
    b1 = a03;
    b4 = a06;
    b2 = a09;
    b3 = a12;

    c1 = sp_ntt_add(b1, b3, p);
    c3 = sp_ntt_sub(b1, b3, p);
    c2 = sp_ntt_add(b2, b4, p);
    c4 = sp_ntt_sub(b2, b4, p);

    b1 = sp_ntt_add(c1, c2, p);
    b2 = sp_ntt_sub_partial(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial(c3, c4, p);

    b0 = sp_ntt_add(b0, b1, p);

    b1 = sp_ntt_mul(b1, ntt_const[1], ntt_const[NC+1], p);
    b2 = sp_ntt_mul(b2, ntt_const[2], ntt_const[NC+2], p);
    b3 = sp_ntt_mul(b3, ntt_const[3], ntt_const[NC+3], p);
    b4 = sp_ntt_mul(b4, ntt_const[4], ntt_const[NC+4], p);
    b5 = sp_ntt_mul(b5, ntt_const[5], ntt_const[NC+5], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

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

    c1 = sp_ntt_add(b1, b3, p);
    c3 = sp_ntt_sub(b1, b3, p);
    c2 = sp_ntt_add(b2, b4, p);
    c4 = sp_ntt_sub(b2, b4, p);

    b1 = sp_ntt_add(c1, c2, p);
    b2 = sp_ntt_sub_partial(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial(c3, c4, p);

    b0 = sp_ntt_add_partial(b0, b1, p);

    b0 = sp_ntt_mul(b0, ntt_const[6], ntt_const[NC+6], p);
    b1 = sp_ntt_mul(b1, ntt_const[7], ntt_const[NC+7], p);
    b2 = sp_ntt_mul(b2, ntt_const[8], ntt_const[NC+8], p);
    b3 = sp_ntt_mul(b3, ntt_const[9], ntt_const[NC+9], p);
    b4 = sp_ntt_mul(b4, ntt_const[10], ntt_const[NC+10], p);
    b5 = sp_ntt_mul(b5, ntt_const[11], ntt_const[NC+11], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

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

    c1 = sp_ntt_add(b1, b3, p);
    c3 = sp_ntt_sub(b1, b3, p);
    c2 = sp_ntt_add(b2, b4, p);
    c4 = sp_ntt_sub(b2, b4, p);

    b1 = sp_ntt_add(c1, c2, p);
    b2 = sp_ntt_sub_partial(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial(c3, c4, p);

    b0 = sp_ntt_add_partial(b0, b1, p);

    b0 = sp_ntt_mul(b0, ntt_const[12], ntt_const[NC+12], p);
    b1 = sp_ntt_mul(b1, ntt_const[13], ntt_const[NC+13], p);
    b2 = sp_ntt_mul(b2, ntt_const[14], ntt_const[NC+14], p);
    b3 = sp_ntt_mul(b3, ntt_const[15], ntt_const[NC+15], p);
    b4 = sp_ntt_mul(b4, ntt_const[16], ntt_const[NC+16], p);
    b5 = sp_ntt_mul(b5, ntt_const[17], ntt_const[NC+17], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

    a02 = b0;
    a05 = b4;
    a08 = b3;
    a11 = b1;
    a14 = b2;
  }
  {
    sp_t x00, x01, x02;

    a01 = sp_ntt_add(a00, a01, p);

    x00 = a00;
    x01 = sp_ntt_add(a01, a02, p);
    x02 = sp_ntt_sub(a01, a02, p);

    out[ 0 * ostride] = x00;
    out[ 5 * ostride] = x02;
    out[10 * ostride] = x01;
  }
  {
    sp_t x03, x04, x05;

    a04 = sp_ntt_add(a03, a04, p);

    x03 = a03;
    x04 = sp_ntt_add(a04, a05, p);
    x05 = sp_ntt_sub(a04, a05, p);

    out[ 1 * ostride] = x04;
    out[ 6 * ostride] = x03;
    out[11 * ostride] = x05;
  }
  {
    sp_t x06, x07, x08;

    a07 = sp_ntt_add(a06, a07, p);

    x06 = a06;
    x07 = sp_ntt_add(a07, a08, p);
    x08 = sp_ntt_sub(a07, a08, p);

    out[ 2 * ostride] = x08;
    out[ 7 * ostride] = x07;
    out[12 * ostride] = x06;
  }
  {
    sp_t x09, x10, x11;

    a10 = sp_ntt_add(a09, a10, p);

    x09 = a09;
    x10 = sp_ntt_add(a10, a11, p);
    x11 = sp_ntt_sub(a10, a11, p);

    out[ 3 * ostride] = x09;
    out[ 8 * ostride] = x11;
    out[13 * ostride] = x10;
  }
  {
    sp_t x12, x13, x14;

    a13 = sp_ntt_add(a12, a13, p);

    x12 = a12;
    x13 = sp_ntt_add(a13, a14, p);
    x14 = sp_ntt_sub(a13, a14, p);

    out[ 4 * ostride] = x13;
    out[ 9 * ostride] = x12;
    out[14 * ostride] = x14;
  }
}

#ifdef HAVE_SSE2
static void
ntt15_run_core_simd(spv_t in, spv_size_t istride, spv_size_t idist,
		spv_t out, spv_size_t ostride, spv_size_t odist,
		sp_t p, spv_t ntt_const)
{
  sp_simd_t a00, a01, a02, a03, a04, 
            a05, a06, a07, a08, a09, 
            a10, a11, a12, a13, a14;

  {
    sp_simd_t x00 = sp_simd_gather(in + 0 * istride, idist);
    sp_simd_t x01 = sp_simd_gather(in + 5 * istride, idist);
    sp_simd_t x02 = sp_simd_gather(in +10 * istride, idist);

    a01 = sp_ntt_add_simd(x01, x02, p);
    a02 = sp_ntt_sub_simd(x01, x02, p);

    a00 = sp_ntt_add_simd(x00, a01, p);
  }
  {
    sp_simd_t x03 = sp_simd_gather(in + 3 * istride, idist);
    sp_simd_t x04 = sp_simd_gather(in + 8 * istride, idist);
    sp_simd_t x05 = sp_simd_gather(in +13 * istride, idist);

    a04 = sp_ntt_add_simd(x04, x05, p);
    a05 = sp_ntt_sub_simd(x04, x05, p);

    a03 = sp_ntt_add_simd(x03, a04, p);
  }
  {
    sp_simd_t x08 = sp_simd_gather(in + 1 * istride, idist);
    sp_simd_t x06 = sp_simd_gather(in + 6 * istride, idist);
    sp_simd_t x07 = sp_simd_gather(in +11 * istride, idist);

    a07 = sp_ntt_add_simd(x07, x08, p);
    a08 = sp_ntt_sub_simd(x07, x08, p);

    a06 = sp_ntt_add_simd(x06, a07, p);
  }
  {
    sp_simd_t x11 = sp_simd_gather(in + 4 * istride, idist);
    sp_simd_t x09 = sp_simd_gather(in + 9 * istride, idist);
    sp_simd_t x10 = sp_simd_gather(in +14 * istride, idist);

    a10 = sp_ntt_add_simd(x10, x11, p);
    a11 = sp_ntt_sub_simd(x10, x11, p);

    a09 = sp_ntt_add_simd(x09, a10, p);
  }
  {
    sp_simd_t x13 = sp_simd_gather(in + 2 * istride, idist);
    sp_simd_t x14 = sp_simd_gather(in + 7 * istride, idist);
    sp_simd_t x12 = sp_simd_gather(in +12 * istride, idist);

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

    sp_simd_scatter(x00, out +  0 * ostride, odist);
    sp_simd_scatter(x02, out +  5 * ostride, odist);
    sp_simd_scatter(x01, out + 10 * ostride, odist);
  }
  {
    sp_simd_t x03, x04, x05;

    a04 = sp_ntt_add_simd(a03, a04, p);

    x03 = a03;
    x04 = sp_ntt_add_simd(a04, a05, p);
    x05 = sp_ntt_sub_simd(a04, a05, p);

    sp_simd_scatter(x04, out +  1 * ostride, odist);
    sp_simd_scatter(x03, out +  6 * ostride, odist);
    sp_simd_scatter(x05, out + 11 * ostride, odist);
  }
  {
    sp_simd_t x06, x07, x08;

    a07 = sp_ntt_add_simd(a06, a07, p);

    x06 = a06;
    x07 = sp_ntt_add_simd(a07, a08, p);
    x08 = sp_ntt_sub_simd(a07, a08, p);

    sp_simd_scatter(x08, out +  2 * ostride, odist);
    sp_simd_scatter(x07, out +  7 * ostride, odist);
    sp_simd_scatter(x06, out + 12 * ostride, odist);
  }
  {
    sp_simd_t x09, x10, x11;

    a10 = sp_ntt_add_simd(a09, a10, p);

    x09 = a09;
    x10 = sp_ntt_add_simd(a10, a11, p);
    x11 = sp_ntt_sub_simd(a10, a11, p);

    sp_simd_scatter(x09, out +  3 * ostride, odist);
    sp_simd_scatter(x11, out +  8 * ostride, odist);
    sp_simd_scatter(x10, out + 13 * ostride, odist);
  }
  {
    sp_simd_t x12, x13, x14;

    a13 = sp_ntt_add_simd(a12, a13, p);

    x12 = a12;
    x13 = sp_ntt_add_simd(a13, a14, p);
    x14 = sp_ntt_sub_simd(a13, a14, p);

    sp_simd_scatter(x13, out +  4 * ostride, odist);
    sp_simd_scatter(x12, out +  9 * ostride, odist);
    sp_simd_scatter(x14, out + 14 * ostride, odist);
  }
}
#endif

static void
ntt15_run(spv_t x, spv_size_t num_transforms,
	  sp_t p, spv_t ntt_const)
{
  spv_size_t i = 0;

#ifdef HAVE_SSE2
  spv_size_t num_simd = SP_SIMD_VSIZE * (num_transforms / SP_SIMD_VSIZE);

  for (i = 0; i < num_simd; i += 15 * SP_SIMD_VSIZE)
    ntt15_run_core_simd(x + i, 1, 15, x + i, 1, 15, p, ntt_const);
#endif

  for (; i < num_transforms; i += 15)
    ntt15_run_core(x + i, 1, x + i, 1, p, ntt_const);
}


static void
ntt15_twiddle_run_core(spv_t x, spv_t w, spv_size_t stride,
			sp_t p, spv_t ntt_const)
{
  sp_t a00, a01, a02, a03, a04, 
       a05, a06, a07, a08, a09, 
       a10, a11, a12, a13, a14;

  {
    sp_t x00 = x[ 0 * stride];
    sp_t x01 = x[ 5 * stride];
    sp_t x02 = x[10 * stride];

    a01 = sp_ntt_add(x01, x02, p);
    a02 = sp_ntt_sub(x01, x02, p);

    a00 = sp_ntt_add(x00, a01, p);
  }
  {
    sp_t x03 = x[ 3 * stride];
    sp_t x04 = x[ 8 * stride];
    sp_t x05 = x[13 * stride];

    a04 = sp_ntt_add(x04, x05, p);
    a05 = sp_ntt_sub(x04, x05, p);

    a03 = sp_ntt_add(x03, a04, p);
  }
  {
    sp_t x08 = x[ 1 * stride];
    sp_t x06 = x[ 6 * stride];
    sp_t x07 = x[11 * stride];

    a07 = sp_ntt_add(x07, x08, p);
    a08 = sp_ntt_sub(x07, x08, p);

    a06 = sp_ntt_add(x06, a07, p);
  }
  {
    sp_t x11 = x[ 4 * stride];
    sp_t x09 = x[ 9 * stride];
    sp_t x10 = x[14 * stride];

    a10 = sp_ntt_add(x10, x11, p);
    a11 = sp_ntt_sub(x10, x11, p);

    a09 = sp_ntt_add(x09, a10, p);
  }
  {
    sp_t x13 = x[ 2 * stride];
    sp_t x14 = x[ 7 * stride];
    sp_t x12 = x[12 * stride];

    a13 = sp_ntt_add(x13, x14, p);
    a14 = sp_ntt_sub(x13, x14, p);

    a12 = sp_ntt_add(x12, a13, p);
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a00;
    b1 = a03;
    b4 = a06;
    b2 = a09;
    b3 = a12;

    c1 = sp_ntt_add(b1, b3, p);
    c3 = sp_ntt_sub(b1, b3, p);
    c2 = sp_ntt_add(b2, b4, p);
    c4 = sp_ntt_sub(b2, b4, p);

    b1 = sp_ntt_add(c1, c2, p);
    b2 = sp_ntt_sub_partial(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial(c3, c4, p);

    b0 = sp_ntt_add(b0, b1, p);

    b1 = sp_ntt_mul(b1, ntt_const[1], ntt_const[NC+1], p);
    b2 = sp_ntt_mul(b2, ntt_const[2], ntt_const[NC+2], p);
    b3 = sp_ntt_mul(b3, ntt_const[3], ntt_const[NC+3], p);
    b4 = sp_ntt_mul(b4, ntt_const[4], ntt_const[NC+4], p);
    b5 = sp_ntt_mul(b5, ntt_const[5], ntt_const[NC+5], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

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

    c1 = sp_ntt_add(b1, b3, p);
    c3 = sp_ntt_sub(b1, b3, p);
    c2 = sp_ntt_add(b2, b4, p);
    c4 = sp_ntt_sub(b2, b4, p);

    b1 = sp_ntt_add(c1, c2, p);
    b2 = sp_ntt_sub_partial(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial(c3, c4, p);

    b0 = sp_ntt_add_partial(b0, b1, p);

    b0 = sp_ntt_mul(b0, ntt_const[6], ntt_const[NC+6], p);
    b1 = sp_ntt_mul(b1, ntt_const[7], ntt_const[NC+7], p);
    b2 = sp_ntt_mul(b2, ntt_const[8], ntt_const[NC+8], p);
    b3 = sp_ntt_mul(b3, ntt_const[9], ntt_const[NC+9], p);
    b4 = sp_ntt_mul(b4, ntt_const[10], ntt_const[NC+10], p);
    b5 = sp_ntt_mul(b5, ntt_const[11], ntt_const[NC+11], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

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

    c1 = sp_ntt_add(b1, b3, p);
    c3 = sp_ntt_sub(b1, b3, p);
    c2 = sp_ntt_add(b2, b4, p);
    c4 = sp_ntt_sub(b2, b4, p);

    b1 = sp_ntt_add(c1, c2, p);
    b2 = sp_ntt_sub_partial(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial(c3, c4, p);

    b0 = sp_ntt_add_partial(b0, b1, p);

    b0 = sp_ntt_mul(b0, ntt_const[12], ntt_const[NC+12], p);
    b1 = sp_ntt_mul(b1, ntt_const[13], ntt_const[NC+13], p);
    b2 = sp_ntt_mul(b2, ntt_const[14], ntt_const[NC+14], p);
    b3 = sp_ntt_mul(b3, ntt_const[15], ntt_const[NC+15], p);
    b4 = sp_ntt_mul(b4, ntt_const[16], ntt_const[NC+16], p);
    b5 = sp_ntt_mul(b5, ntt_const[17], ntt_const[NC+17], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

    a02 = b0;
    a05 = b4;
    a08 = b3;
    a11 = b1;
    a14 = b2;
  }
  {
    sp_t x00, x01, x02;

    a01 = sp_ntt_add(a00, a01, p);

    x00 = a00;
    x01 = sp_ntt_add_partial(a01, a02, p);
    x02 = sp_ntt_sub_partial(a01, a02, p);

    x02 = sp_ntt_mul(x02, w[8], w[9], p);
    x01 = sp_ntt_mul(x01, w[18], w[19], p);

    x[ 0 * stride] = x00;
    x[ 5 * stride] = x02;
    x[10 * stride] = x01;
  }
  {
    sp_t x03, x04, x05;

    a04 = sp_ntt_add(a03, a04, p);

    x03 = a03;
    x04 = sp_ntt_add_partial(a04, a05, p);
    x05 = sp_ntt_sub_partial(a04, a05, p);

    x04 = sp_ntt_mul(x04, w[0], w[1], p);
    x03 = sp_ntt_mul(x03, w[10], w[11], p);
    x05 = sp_ntt_mul(x05, w[20], w[21], p);

    x[ 1 * stride] = x04;
    x[ 6 * stride] = x03;
    x[11 * stride] = x05;
  }
  {
    sp_t x06, x07, x08;

    a07 = sp_ntt_add(a06, a07, p);

    x06 = a06;
    x07 = sp_ntt_add_partial(a07, a08, p);
    x08 = sp_ntt_sub_partial(a07, a08, p);

    x08 = sp_ntt_mul(x08, w[2], w[3], p);
    x07 = sp_ntt_mul(x07, w[12], w[13], p);
    x06 = sp_ntt_mul(x06, w[22], w[23], p);

    x[ 2 * stride] = x08;
    x[ 7 * stride] = x07;
    x[12 * stride] = x06;
  }
  {
    sp_t x09, x10, x11;

    a10 = sp_ntt_add(a09, a10, p);

    x09 = a09;
    x10 = sp_ntt_add_partial(a10, a11, p);
    x11 = sp_ntt_sub_partial(a10, a11, p);

    x09 = sp_ntt_mul(x09, w[4], w[5], p);
    x11 = sp_ntt_mul(x11, w[14], w[15], p);
    x10 = sp_ntt_mul(x10, w[24], w[25], p);

    x[ 3 * stride] = x09;
    x[ 8 * stride] = x11;
    x[13 * stride] = x10;
  }
  {
    sp_t x12, x13, x14;

    a13 = sp_ntt_add(a12, a13, p);

    x12 = a12;
    x13 = sp_ntt_add_partial(a13, a14, p);
    x14 = sp_ntt_sub_partial(a13, a14, p);

    x13 = sp_ntt_mul(x13, w[6], w[7], p);
    x12 = sp_ntt_mul(x12, w[16], w[17], p);
    x14 = sp_ntt_mul(x14, w[26], w[27], p);

    x[ 4 * stride] = x13;
    x[ 9 * stride] = x12;
    x[14 * stride] = x14;
  }
}

#ifdef HAVE_SSE2
static void
ntt15_twiddle_run_core_simd(spv_t x, sp_simd_t *w,
			spv_size_t stride,
			sp_t p, spv_t ntt_const)
{
  sp_simd_t a00, a01, a02, a03, a04, 
            a05, a06, a07, a08, a09, 
            a10, a11, a12, a13, a14;

  {
    sp_simd_t x00 = sp_simd_load(x +  0 * stride);
    sp_simd_t x01 = sp_simd_load(x +  5 * stride);
    sp_simd_t x02 = sp_simd_load(x + 10 * stride);

    a01 = sp_ntt_add_simd(x01, x02, p);
    a02 = sp_ntt_sub_simd(x01, x02, p);

    a00 = sp_ntt_add_simd(x00, a01, p);
  }
  {
    sp_simd_t x03 = sp_simd_load(x +  3 * stride);
    sp_simd_t x04 = sp_simd_load(x +  8 * stride);
    sp_simd_t x05 = sp_simd_load(x + 13 * stride);

    a04 = sp_ntt_add_simd(x04, x05, p);
    a05 = sp_ntt_sub_simd(x04, x05, p);

    a03 = sp_ntt_add_simd(x03, a04, p);
  }
  {
    sp_simd_t x08 = sp_simd_load(x +  1 * stride);
    sp_simd_t x06 = sp_simd_load(x +  6 * stride);
    sp_simd_t x07 = sp_simd_load(x + 11 * stride);

    a07 = sp_ntt_add_simd(x07, x08, p);
    a08 = sp_ntt_sub_simd(x07, x08, p);

    a06 = sp_ntt_add_simd(x06, a07, p);
  }
  {
    sp_simd_t x11 = sp_simd_load(x +  4 * stride);
    sp_simd_t x09 = sp_simd_load(x +  9 * stride);
    sp_simd_t x10 = sp_simd_load(x + 14 * stride);

    a10 = sp_ntt_add_simd(x10, x11, p);
    a11 = sp_ntt_sub_simd(x10, x11, p);

    a09 = sp_ntt_add_simd(x09, a10, p);
  }
  {
    sp_simd_t x13 = sp_simd_load(x +  2 * stride);
    sp_simd_t x14 = sp_simd_load(x +  7 * stride);
    sp_simd_t x12 = sp_simd_load(x + 12 * stride);

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

    sp_simd_store(x00, x +  0 * stride);
    sp_simd_store(x02, x +  5 * stride);
    sp_simd_store(x01, x + 10 * stride);
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

    sp_simd_store(x04, x +  1 * stride);
    sp_simd_store(x03, x +  6 * stride);
    sp_simd_store(x05, x + 11 * stride);
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

    sp_simd_store(x08, x +  2 * stride);
    sp_simd_store(x07, x +  7 * stride);
    sp_simd_store(x06, x + 12 * stride);
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

    sp_simd_store(x09, x +  3 * stride);
    sp_simd_store(x11, x +  8 * stride);
    sp_simd_store(x10, x + 13 * stride);
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

    sp_simd_store(x13, x +  4 * stride);
    sp_simd_store(x12, x +  9 * stride);
    sp_simd_store(x14, x + 14 * stride);
  }
}
#endif

static void
ntt15_twiddle_run(spv_t x, spv_t w,
	  spv_size_t stride,
	  spv_size_t num_transforms,
	  sp_t p, spv_t ntt_const)
{
  spv_size_t i = 0, j = 0;

#ifdef HAVE_SSE2
  spv_size_t num_simd = SP_SIMD_VSIZE * (num_transforms / SP_SIMD_VSIZE);

  for (i = 0; i < num_simd; i += SP_SIMD_VSIZE,
		  	j += 2*(15-1)*SP_SIMD_VSIZE)
    ntt15_twiddle_run_core_simd(x + i, (sp_simd_t *)(w + j),
				stride, p, ntt_const);
#endif

  for (; i < num_transforms; i++, j += 2*(15-1))
    ntt15_twiddle_run_core(x + i, w + j, stride, p, ntt_const);
}

static void
ntt15_pfa_run_core(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t n,
	  sp_t p, spv_t ntt_const)
{
  spv_size_t j00, j01, j02, j03, j04,
             j05, j06, j07, j08, j09,
	     j10, j11, j12, j13, j14;

  sp_t a00, a01, a02, a03, a04, 
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
    sp_t x00 = x[j00];
    sp_t x01 = x[j05];
    sp_t x02 = x[j10];

    a01 = sp_ntt_add(x01, x02, p);
    a02 = sp_ntt_sub(x01, x02, p);

    a00 = sp_ntt_add(x00, a01, p);
  }
  {
    sp_t x03 = x[j03];
    sp_t x04 = x[j08];
    sp_t x05 = x[j13];

    a04 = sp_ntt_add(x04, x05, p);
    a05 = sp_ntt_sub(x04, x05, p);

    a03 = sp_ntt_add(x03, a04, p);
  }
  {
    sp_t x08 = x[j01];
    sp_t x06 = x[j06];
    sp_t x07 = x[j11];

    a07 = sp_ntt_add(x07, x08, p);
    a08 = sp_ntt_sub(x07, x08, p);

    a06 = sp_ntt_add(x06, a07, p);
  }
  {
    sp_t x11 = x[j04];
    sp_t x09 = x[j09];
    sp_t x10 = x[j14];

    a10 = sp_ntt_add(x10, x11, p);
    a11 = sp_ntt_sub(x10, x11, p);

    a09 = sp_ntt_add(x09, a10, p);
  }
  {
    sp_t x13 = x[j02];
    sp_t x14 = x[j07];
    sp_t x12 = x[j12];

    a13 = sp_ntt_add(x13, x14, p);
    a14 = sp_ntt_sub(x13, x14, p);

    a12 = sp_ntt_add(x12, a13, p);
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a00;
    b1 = a03;
    b4 = a06;
    b2 = a09;
    b3 = a12;

    c1 = sp_ntt_add(b1, b3, p);
    c3 = sp_ntt_sub(b1, b3, p);
    c2 = sp_ntt_add(b2, b4, p);
    c4 = sp_ntt_sub(b2, b4, p);

    b1 = sp_ntt_add(c1, c2, p);
    b2 = sp_ntt_sub_partial(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial(c3, c4, p);

    b0 = sp_ntt_add(b0, b1, p);

    b1 = sp_ntt_mul(b1, ntt_const[1], ntt_const[NC+1], p);
    b2 = sp_ntt_mul(b2, ntt_const[2], ntt_const[NC+2], p);
    b3 = sp_ntt_mul(b3, ntt_const[3], ntt_const[NC+3], p);
    b4 = sp_ntt_mul(b4, ntt_const[4], ntt_const[NC+4], p);
    b5 = sp_ntt_mul(b5, ntt_const[5], ntt_const[NC+5], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

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

    c1 = sp_ntt_add(b1, b3, p);
    c3 = sp_ntt_sub(b1, b3, p);
    c2 = sp_ntt_add(b2, b4, p);
    c4 = sp_ntt_sub(b2, b4, p);

    b1 = sp_ntt_add(c1, c2, p);
    b2 = sp_ntt_sub_partial(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial(c3, c4, p);

    b0 = sp_ntt_add_partial(b0, b1, p);

    b0 = sp_ntt_mul(b0, ntt_const[6], ntt_const[NC+6], p);
    b1 = sp_ntt_mul(b1, ntt_const[7], ntt_const[NC+7], p);
    b2 = sp_ntt_mul(b2, ntt_const[8], ntt_const[NC+8], p);
    b3 = sp_ntt_mul(b3, ntt_const[9], ntt_const[NC+9], p);
    b4 = sp_ntt_mul(b4, ntt_const[10], ntt_const[NC+10], p);
    b5 = sp_ntt_mul(b5, ntt_const[11], ntt_const[NC+11], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

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

    c1 = sp_ntt_add(b1, b3, p);
    c3 = sp_ntt_sub(b1, b3, p);
    c2 = sp_ntt_add(b2, b4, p);
    c4 = sp_ntt_sub(b2, b4, p);

    b1 = sp_ntt_add(c1, c2, p);
    b2 = sp_ntt_sub_partial(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial(c3, c4, p);

    b0 = sp_ntt_add_partial(b0, b1, p);

    b0 = sp_ntt_mul(b0, ntt_const[12], ntt_const[NC+12], p);
    b1 = sp_ntt_mul(b1, ntt_const[13], ntt_const[NC+13], p);
    b2 = sp_ntt_mul(b2, ntt_const[14], ntt_const[NC+14], p);
    b3 = sp_ntt_mul(b3, ntt_const[15], ntt_const[NC+15], p);
    b4 = sp_ntt_mul(b4, ntt_const[16], ntt_const[NC+16], p);
    b5 = sp_ntt_mul(b5, ntt_const[17], ntt_const[NC+17], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

    a02 = b0;
    a05 = b4;
    a08 = b3;
    a11 = b1;
    a14 = b2;
  }
  {
    sp_t x00, x01, x02;

    a01 = sp_ntt_add(a00, a01, p);

    x00 = a00;
    x01 = sp_ntt_add(a01, a02, p);
    x02 = sp_ntt_sub(a01, a02, p);

    x[j00] = x00;
    x[j05] = x02;
    x[j10] = x01;
  }
  {
    sp_t x03, x04, x05;

    a04 = sp_ntt_add(a03, a04, p);

    x03 = a03;
    x04 = sp_ntt_add(a04, a05, p);
    x05 = sp_ntt_sub(a04, a05, p);

    x[j01] = x04;
    x[j06] = x03;
    x[j11] = x05;
  }
  {
    sp_t x06, x07, x08;

    a07 = sp_ntt_add(a06, a07, p);

    x06 = a06;
    x07 = sp_ntt_add(a07, a08, p);
    x08 = sp_ntt_sub(a07, a08, p);

    x[j02] = x08;
    x[j07] = x07;
    x[j12] = x06;
  }
  {
    sp_t x09, x10, x11;

    a10 = sp_ntt_add(a09, a10, p);

    x09 = a09;
    x10 = sp_ntt_add(a10, a11, p);
    x11 = sp_ntt_sub(a10, a11, p);

    x[j03] = x09;
    x[j08] = x11;
    x[j13] = x10;
  }
  {
    sp_t x12, x13, x14;

    a13 = sp_ntt_add(a12, a13, p);

    x12 = a12;
    x13 = sp_ntt_add(a13, a14, p);
    x14 = sp_ntt_sub(a13, a14, p);

    x[j04] = x13;
    x[j09] = x12;
    x[j14] = x14;
  }
}

#ifdef HAVE_SSE2
static void
ntt15_pfa_run_core_simd(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t inc2, spv_size_t n,
	  sp_t p, spv_t ntt_const)
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
    sp_simd_t x00 = sp_simd_pfa_gather(x, j00, inc2, n);
    sp_simd_t x01 = sp_simd_pfa_gather(x, j05, inc2, n);
    sp_simd_t x02 = sp_simd_pfa_gather(x, j10, inc2, n);

    a01 = sp_ntt_add_simd(x01, x02, p);
    a02 = sp_ntt_sub_simd(x01, x02, p);

    a00 = sp_ntt_add_simd(x00, a01, p);
  }
  {
    sp_simd_t x03 = sp_simd_pfa_gather(x, j03, inc2, n);
    sp_simd_t x04 = sp_simd_pfa_gather(x, j08, inc2, n);
    sp_simd_t x05 = sp_simd_pfa_gather(x, j13, inc2, n);

    a04 = sp_ntt_add_simd(x04, x05, p);
    a05 = sp_ntt_sub_simd(x04, x05, p);

    a03 = sp_ntt_add_simd(x03, a04, p);
  }
  {
    sp_simd_t x08 = sp_simd_pfa_gather(x, j01, inc2, n);
    sp_simd_t x06 = sp_simd_pfa_gather(x, j06, inc2, n);
    sp_simd_t x07 = sp_simd_pfa_gather(x, j11, inc2, n);

    a07 = sp_ntt_add_simd(x07, x08, p);
    a08 = sp_ntt_sub_simd(x07, x08, p);

    a06 = sp_ntt_add_simd(x06, a07, p);
  }
  {
    sp_simd_t x11 = sp_simd_pfa_gather(x, j04, inc2, n);
    sp_simd_t x09 = sp_simd_pfa_gather(x, j09, inc2, n);
    sp_simd_t x10 = sp_simd_pfa_gather(x, j14, inc2, n);

    a10 = sp_ntt_add_simd(x10, x11, p);
    a11 = sp_ntt_sub_simd(x10, x11, p);

    a09 = sp_ntt_add_simd(x09, a10, p);
  }
  {
    sp_simd_t x13 = sp_simd_pfa_gather(x, j02, inc2, n);
    sp_simd_t x14 = sp_simd_pfa_gather(x, j07, inc2, n);
    sp_simd_t x12 = sp_simd_pfa_gather(x, j12, inc2, n);

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

    sp_simd_pfa_scatter(x00, x, j00, inc2, n);
    sp_simd_pfa_scatter(x02, x, j05, inc2, n);
    sp_simd_pfa_scatter(x01, x, j10, inc2, n);
  }
  {
    sp_simd_t x03, x04, x05;

    a04 = sp_ntt_add_simd(a03, a04, p);

    x03 = a03;
    x04 = sp_ntt_add_simd(a04, a05, p);
    x05 = sp_ntt_sub_simd(a04, a05, p);

    sp_simd_pfa_scatter(x04, x, j01, inc2, n);
    sp_simd_pfa_scatter(x03, x, j06, inc2, n);
    sp_simd_pfa_scatter(x05, x, j11, inc2, n);
  }
  {
    sp_simd_t x06, x07, x08;

    a07 = sp_ntt_add_simd(a06, a07, p);

    x06 = a06;
    x07 = sp_ntt_add_simd(a07, a08, p);
    x08 = sp_ntt_sub_simd(a07, a08, p);

    sp_simd_pfa_scatter(x08, x, j02, inc2, n);
    sp_simd_pfa_scatter(x07, x, j07, inc2, n);
    sp_simd_pfa_scatter(x06, x, j12, inc2, n);
  }
  {
    sp_simd_t x09, x10, x11;

    a10 = sp_ntt_add_simd(a09, a10, p);

    x09 = a09;
    x10 = sp_ntt_add_simd(a10, a11, p);
    x11 = sp_ntt_sub_simd(a10, a11, p);

    sp_simd_pfa_scatter(x09, x, j03, inc2, n);
    sp_simd_pfa_scatter(x11, x, j08, inc2, n);
    sp_simd_pfa_scatter(x10, x, j13, inc2, n);
  }
  {
    sp_simd_t x12, x13, x14;

    a13 = sp_ntt_add_simd(a12, a13, p);

    x12 = a12;
    x13 = sp_ntt_add_simd(a13, a14, p);
    x14 = sp_ntt_sub_simd(a13, a14, p);

    sp_simd_pfa_scatter(x13, x, j04, inc2, n);
    sp_simd_pfa_scatter(x12, x, j09, inc2, n);
    sp_simd_pfa_scatter(x14, x, j14, inc2, n);
  }
}
#endif

static void
ntt15_pfa_run(spv_t x, spv_size_t cofactor,
	  sp_t p, spv_t ntt_const)
{
  spv_size_t i = 0;
  spv_size_t incstart = 0;
  spv_size_t n = 15 * cofactor;
  spv_size_t inc = cofactor;
  spv_size_t inc2 = 15;

#ifdef HAVE_SSE2
  spv_size_t num_simd = SP_SIMD_VSIZE * (cofactor / SP_SIMD_VSIZE);

  for (i = 0; i < num_simd; i += SP_SIMD_VSIZE)
    {
      ntt15_pfa_run_core_simd(x, incstart, inc, inc2, n, p, ntt_const);
      incstart += SP_SIMD_VSIZE * inc2;
    }
#endif

  for (; i < cofactor; i++, incstart += inc2)
    ntt15_pfa_run_core(x, incstart, inc, n, p, ntt_const);

}

const nttconfig_t ntt15_config = 
{
  15,
  NC, 
  ntt15_get_fixed_ntt_const,
  ntt15_init,
  ntt15_run,
  ntt15_pfa_run,
  ntt15_twiddle_run
};

