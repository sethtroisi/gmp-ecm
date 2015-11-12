#include "ntt/ntt-impl-scalar.h"

#define NC 18

static const uint8_t ntt15_fixed_const[NC] = {1};

void
X(ntt15_init)(spv_t out, sp_t p, sp_t d, 
	  sp_t primroot, sp_t order, sp_t perm)
{
  X(nttdata_init_generic)(&X(ntt15_config), out, p, d, primroot, order, perm);
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

static void
ntt15_twiddle_run_core(spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		spv_t w, sp_t p, spv_t ntt_const)
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
    x01 = sp_ntt_add_partial(a01, a02, p);
    x02 = sp_ntt_sub_partial(a01, a02, p);

    x02 = sp_ntt_mul(x02, w[8], w[9], p);
    x01 = sp_ntt_mul(x01, w[18], w[19], p);

    out[ 0 * ostride] = x00;
    out[ 5 * ostride] = x02;
    out[10 * ostride] = x01;
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

    out[ 1 * ostride] = x04;
    out[ 6 * ostride] = x03;
    out[11 * ostride] = x05;
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

    out[ 2 * ostride] = x08;
    out[ 7 * ostride] = x07;
    out[12 * ostride] = x06;
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

    out[ 3 * ostride] = x09;
    out[ 8 * ostride] = x11;
    out[13 * ostride] = x10;
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

    out[ 4 * ostride] = x13;
    out[ 9 * ostride] = x12;
    out[14 * ostride] = x14;
  }
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

DECLARE_CORE_ROUTINES(15)
