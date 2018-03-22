#include "ntt/ntt-impl-scalar.h"

#define NC 48

static const uint8_t ntt40_fixed_const[NC] = {1, 0, 0, 0, 0, 0,
					1, 0, 0, 0, 0, 0,
					1, 0, 0, 0, 0, 0,
					0, 0, 0, 0, 0, 0,
					1};
void
X(ntt40_init)(spv_t out, sp_t p, sp_t d, 
	  sp_t primroot, sp_t order, sp_t perm)
{
  X(nttdata_init_generic)(&X(ntt40_config), out, p, d, primroot, order, perm);
}

static void 
ntt40_run_core(spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		sp_t p, spv_t ntt_const)
{
  sp_t a00, a01, a02, a03, a04, a05, a06, a07,
       a08, a09, a10, a11, a12, a13, a14, a15,
       a16, a17, a18, a19, a20, a21, a22, a23,
       a24, a25, a26, a27, a28, a29, a30, a31,
       a32, a33, a34, a35, a36, a37, a38, a39;

  {
    sp_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_t x0 = in[ 0 * istride];
    sp_t x1 = in[ 5 * istride];
    sp_t x2 = in[10 * istride];
    sp_t x3 = in[15 * istride];
    sp_t x4 = in[20 * istride];
    sp_t x5 = in[25 * istride];
    sp_t x6 = in[30 * istride];
    sp_t x7 = in[35 * istride];

    t0 = sp_ntt_add(x0, x4, p);
    t4 = sp_ntt_sub(x0, x4, p);
    t1 = sp_ntt_add(x1, x5, p);
    t5 = sp_ntt_sub(x1, x5, p);
    t2 = sp_ntt_add(x2, x6, p);
    t6 = sp_ntt_sub(x2, x6, p);
    t3 = sp_ntt_add(x3, x7, p);
    t7 = sp_ntt_sub(x3, x7, p);

    a00 = sp_ntt_add(t0, t2, p);
    a01 = sp_ntt_sub(t0, t2, p);
    a02 = sp_ntt_add(t1, t3, p);
    a03 = sp_ntt_sub(t1, t3, p);
    a04 = t4;
    a05 = t6;
    a06 = sp_ntt_sub(t5, t7, p);
    a07 = sp_ntt_add(t5, t7, p); 
  }
  {
    sp_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_t x7 = in[ 3 * istride];
    sp_t x0 = in[ 8 * istride];
    sp_t x1 = in[13 * istride];
    sp_t x2 = in[18 * istride];
    sp_t x3 = in[23 * istride];
    sp_t x4 = in[28 * istride];
    sp_t x5 = in[33 * istride];
    sp_t x6 = in[38 * istride];

    t0 = sp_ntt_add(x0, x4, p);
    t4 = sp_ntt_sub(x0, x4, p);
    t1 = sp_ntt_add(x1, x5, p);
    t5 = sp_ntt_sub(x1, x5, p);
    t2 = sp_ntt_add(x2, x6, p);
    t6 = sp_ntt_sub(x2, x6, p);
    t3 = sp_ntt_add(x3, x7, p);
    t7 = sp_ntt_sub(x3, x7, p);

    a08 = sp_ntt_add(t0, t2, p);
    a09 = sp_ntt_sub(t0, t2, p);
    a10 = sp_ntt_add(t1, t3, p);
    a11 = sp_ntt_sub(t1, t3, p);
    a12 = t4;
    a13 = t6;
    a14 = sp_ntt_sub(t5, t7, p);
    a15 = sp_ntt_add(t5, t7, p); 
  }
  {
    sp_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_t x5 = in[ 1 * istride];
    sp_t x6 = in[ 6 * istride];
    sp_t x7 = in[11 * istride];
    sp_t x0 = in[16 * istride];
    sp_t x1 = in[21 * istride];
    sp_t x2 = in[26 * istride];
    sp_t x3 = in[31 * istride];
    sp_t x4 = in[36 * istride];

    t0 = sp_ntt_add(x0, x4, p);
    t4 = sp_ntt_sub(x0, x4, p);
    t1 = sp_ntt_add(x1, x5, p);
    t5 = sp_ntt_sub(x1, x5, p);
    t2 = sp_ntt_add(x2, x6, p);
    t6 = sp_ntt_sub(x2, x6, p);
    t3 = sp_ntt_add(x3, x7, p);
    t7 = sp_ntt_sub(x3, x7, p);

    a16 = sp_ntt_add(t0, t2, p);
    a17 = sp_ntt_sub(t0, t2, p);
    a18 = sp_ntt_add(t1, t3, p);
    a19 = sp_ntt_sub(t1, t3, p);
    a20 = t4;
    a21 = t6;
    a22 = sp_ntt_sub(t5, t7, p);
    a23 = sp_ntt_add(t5, t7, p); 
  }
  {
    sp_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_t x0 = in[24 * istride];
    sp_t x1 = in[29 * istride];
    sp_t x2 = in[34 * istride];
    sp_t x3 = in[39 * istride];
    sp_t x4 = in[ 4 * istride];
    sp_t x5 = in[ 9 * istride];
    sp_t x6 = in[14 * istride];
    sp_t x7 = in[19 * istride];

    t0 = sp_ntt_add(x0, x4, p);
    t4 = sp_ntt_sub(x0, x4, p);
    t1 = sp_ntt_add(x1, x5, p);
    t5 = sp_ntt_sub(x1, x5, p);
    t2 = sp_ntt_add(x2, x6, p);
    t6 = sp_ntt_sub(x2, x6, p);
    t3 = sp_ntt_add(x3, x7, p);
    t7 = sp_ntt_sub(x3, x7, p);

    a24 = sp_ntt_add(t0, t2, p);
    a25 = sp_ntt_sub(t0, t2, p);
    a26 = sp_ntt_add(t1, t3, p);
    a27 = sp_ntt_sub(t1, t3, p);
    a28 = t4;
    a29 = t6;
    a30 = sp_ntt_sub(t5, t7, p);
    a31 = sp_ntt_add(t5, t7, p); 
  }
  {
    sp_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_t x2 = in[ 2 * istride];
    sp_t x3 = in[ 7 * istride];
    sp_t x4 = in[12 * istride];
    sp_t x5 = in[17 * istride];
    sp_t x6 = in[22 * istride];
    sp_t x7 = in[27 * istride];
    sp_t x0 = in[32 * istride];
    sp_t x1 = in[37 * istride];

    t0 = sp_ntt_add(x0, x4, p);
    t4 = sp_ntt_sub(x0, x4, p);
    t1 = sp_ntt_add(x1, x5, p);
    t5 = sp_ntt_sub(x1, x5, p);
    t2 = sp_ntt_add(x2, x6, p);
    t6 = sp_ntt_sub(x2, x6, p);
    t3 = sp_ntt_add(x3, x7, p);
    t7 = sp_ntt_sub(x3, x7, p);

    a32 = sp_ntt_add(t0, t2, p);
    a33 = sp_ntt_sub(t0, t2, p);
    a34 = sp_ntt_add(t1, t3, p);
    a35 = sp_ntt_sub(t1, t3, p);
    a36 = t4;
    a37 = t6;
    a38 = sp_ntt_sub(t5, t7, p);
    a39 = sp_ntt_add(t5, t7, p); 
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a00;
    b1 = a08;
    b4 = a16;
    b2 = a24;
    b3 = a32;

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
    a08 = b4;
    a16 = b3;
    a24 = b1;
    a32 = b2;
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a01;
    b1 = a09;
    b4 = a17;
    b2 = a25;
    b3 = a33;

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
    a09 = b4;
    a17 = b3;
    a25 = b1;
    a33 = b2;
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a02;
    b1 = a10;
    b4 = a18;
    b2 = a26;
    b3 = a34;

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
    a10 = b4;
    a18 = b3;
    a26 = b1;
    a34 = b2;
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a03;
    b1 = a11;
    b4 = a19;
    b2 = a27;
    b3 = a35;

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

    b0 = sp_ntt_mul(b0, ntt_const[18], ntt_const[NC+18], p);
    b1 = sp_ntt_mul(b1, ntt_const[19], ntt_const[NC+19], p);
    b2 = sp_ntt_mul(b2, ntt_const[20], ntt_const[NC+20], p);
    b3 = sp_ntt_mul(b3, ntt_const[21], ntt_const[NC+21], p);
    b4 = sp_ntt_mul(b4, ntt_const[22], ntt_const[NC+22], p);
    b5 = sp_ntt_mul(b5, ntt_const[23], ntt_const[NC+23], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

    a03 = b0;
    a11 = b4;
    a19 = b3;
    a27 = b1;
    a35 = b2;
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a04;
    b1 = a12;
    b4 = a20;
    b2 = a28;
    b3 = a36;

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

    b1 = sp_ntt_mul(b1, ntt_const[25], ntt_const[NC+25], p);
    b2 = sp_ntt_mul(b2, ntt_const[26], ntt_const[NC+26], p);
    b3 = sp_ntt_mul(b3, ntt_const[27], ntt_const[NC+27], p);
    b4 = sp_ntt_mul(b4, ntt_const[28], ntt_const[NC+28], p);
    b5 = sp_ntt_mul(b5, ntt_const[29], ntt_const[NC+29], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

    a04 = b0;
    a12 = b4;
    a20 = b3;
    a28 = b1;
    a36 = b2;
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a05;
    b1 = a13;
    b4 = a21;
    b2 = a29;
    b3 = a37;

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

    b0 = sp_ntt_mul(b0, ntt_const[30], ntt_const[NC+30], p);
    b1 = sp_ntt_mul(b1, ntt_const[31], ntt_const[NC+31], p);
    b2 = sp_ntt_mul(b2, ntt_const[32], ntt_const[NC+32], p);
    b3 = sp_ntt_mul(b3, ntt_const[33], ntt_const[NC+33], p);
    b4 = sp_ntt_mul(b4, ntt_const[34], ntt_const[NC+34], p);
    b5 = sp_ntt_mul(b5, ntt_const[35], ntt_const[NC+35], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

    a05 = b0;
    a13 = b4;
    a21 = b3;
    a29 = b1;
    a37 = b2;
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a06;
    b1 = a14;
    b4 = a22;
    b2 = a30;
    b3 = a38;

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

    b0 = sp_ntt_mul(b0, ntt_const[36], ntt_const[NC+36], p);
    b1 = sp_ntt_mul(b1, ntt_const[37], ntt_const[NC+37], p);
    b2 = sp_ntt_mul(b2, ntt_const[38], ntt_const[NC+38], p);
    b3 = sp_ntt_mul(b3, ntt_const[39], ntt_const[NC+39], p);
    b4 = sp_ntt_mul(b4, ntt_const[40], ntt_const[NC+40], p);
    b5 = sp_ntt_mul(b5, ntt_const[41], ntt_const[NC+41], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

    a06 = b0;
    a14 = b4;
    a22 = b3;
    a30 = b1;
    a38 = b2;
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a07;
    b1 = a15;
    b4 = a23;
    b2 = a31;
    b3 = a39;

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

    b0 = sp_ntt_mul(b0, ntt_const[42], ntt_const[NC+42], p);
    b1 = sp_ntt_mul(b1, ntt_const[43], ntt_const[NC+43], p);
    b2 = sp_ntt_mul(b2, ntt_const[44], ntt_const[NC+44], p);
    b3 = sp_ntt_mul(b3, ntt_const[45], ntt_const[NC+45], p);
    b4 = sp_ntt_mul(b4, ntt_const[46], ntt_const[NC+46], p);
    b5 = sp_ntt_mul(b5, ntt_const[47], ntt_const[NC+47], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

    a07 = b0;
    a15 = b4;
    a23 = b3;
    a31 = b1;
    a39 = b2;
  }
  {
    sp_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_t p0 = a00;
    sp_t p1 = a01;
    sp_t p2 = a02;
    sp_t p3 = a03;
    sp_t p4 = a04;
    sp_t p5 = a05;
    sp_t p6 = a06;
    sp_t p7 = a07;

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

    out[ 0 * ostride] = t0;
    out[ 5 * ostride] = t5;
    out[10 * ostride] = t2;
    out[15 * ostride] = t6;
    out[20 * ostride] = t1;
    out[25 * ostride] = t4;
    out[30 * ostride] = t3;
    out[35 * ostride] = t7;
  }
  {
    sp_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_t p0 = a08;
    sp_t p1 = a09;
    sp_t p2 = a10;
    sp_t p3 = a11;
    sp_t p4 = a12;
    sp_t p5 = a13;
    sp_t p6 = a14;
    sp_t p7 = a15;

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

    out[ 1 * ostride] = t4;
    out[ 6 * ostride] = t3;
    out[11 * ostride] = t7;
    out[16 * ostride] = t0;
    out[21 * ostride] = t5;
    out[26 * ostride] = t2;
    out[31 * ostride] = t6;
    out[36 * ostride] = t1;
  }
  {
    sp_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_t p0 = a16;
    sp_t p1 = a17;
    sp_t p2 = a18;
    sp_t p3 = a19;
    sp_t p4 = a20;
    sp_t p5 = a21;
    sp_t p6 = a22;
    sp_t p7 = a23;

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

    out[ 2 * ostride] = t2;
    out[ 7 * ostride] = t6;
    out[12 * ostride] = t1;
    out[17 * ostride] = t4;
    out[22 * ostride] = t3;
    out[27 * ostride] = t7;
    out[32 * ostride] = t0;
    out[37 * ostride] = t5;
  }
  {
    sp_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_t p0 = a24;
    sp_t p1 = a25;
    sp_t p2 = a26;
    sp_t p3 = a27;
    sp_t p4 = a28;
    sp_t p5 = a29;
    sp_t p6 = a30;
    sp_t p7 = a31;

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

    out[ 3 * ostride] = t7;
    out[ 8 * ostride] = t0;
    out[13 * ostride] = t5;
    out[18 * ostride] = t2;
    out[23 * ostride] = t6;
    out[28 * ostride] = t1;
    out[33 * ostride] = t4;
    out[38 * ostride] = t3;
  }
  {
    sp_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_t p0 = a32;
    sp_t p1 = a33;
    sp_t p2 = a34;
    sp_t p3 = a35;
    sp_t p4 = a36;
    sp_t p5 = a37;
    sp_t p6 = a38;
    sp_t p7 = a39;

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

    out[ 4 * ostride] = t1;
    out[ 9 * ostride] = t4;
    out[14 * ostride] = t3;
    out[19 * ostride] = t7;
    out[24 * ostride] = t0;
    out[29 * ostride] = t5;
    out[34 * ostride] = t2;
    out[39 * ostride] = t6;
  }
}

static void
ntt40_twiddle_run_core(spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		spv_t w, sp_t p, spv_t ntt_const)
{
  sp_t a00, a01, a02, a03, a04, a05, a06, a07,
       a08, a09, a10, a11, a12, a13, a14, a15,
       a16, a17, a18, a19, a20, a21, a22, a23,
       a24, a25, a26, a27, a28, a29, a30, a31,
       a32, a33, a34, a35, a36, a37, a38, a39;

  {
    sp_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_t x0 = in[ 0 * istride];
    sp_t x1 = in[ 5 * istride];
    sp_t x2 = in[10 * istride];
    sp_t x3 = in[15 * istride];
    sp_t x4 = in[20 * istride];
    sp_t x5 = in[25 * istride];
    sp_t x6 = in[30 * istride];
    sp_t x7 = in[35 * istride];

    t0 = sp_ntt_add(x0, x4, p);
    t4 = sp_ntt_sub(x0, x4, p);
    t1 = sp_ntt_add(x1, x5, p);
    t5 = sp_ntt_sub(x1, x5, p);
    t2 = sp_ntt_add(x2, x6, p);
    t6 = sp_ntt_sub(x2, x6, p);
    t3 = sp_ntt_add(x3, x7, p);
    t7 = sp_ntt_sub(x3, x7, p);

    a00 = sp_ntt_add(t0, t2, p);
    a01 = sp_ntt_sub(t0, t2, p);
    a02 = sp_ntt_add(t1, t3, p);
    a03 = sp_ntt_sub(t1, t3, p);
    a04 = t4;
    a05 = t6;
    a06 = sp_ntt_sub(t5, t7, p);
    a07 = sp_ntt_add(t5, t7, p); 
  }
  {
    sp_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_t x7 = in[ 3 * istride];
    sp_t x0 = in[ 8 * istride];
    sp_t x1 = in[13 * istride];
    sp_t x2 = in[18 * istride];
    sp_t x3 = in[23 * istride];
    sp_t x4 = in[28 * istride];
    sp_t x5 = in[33 * istride];
    sp_t x6 = in[38 * istride];

    t0 = sp_ntt_add(x0, x4, p);
    t4 = sp_ntt_sub(x0, x4, p);
    t1 = sp_ntt_add(x1, x5, p);
    t5 = sp_ntt_sub(x1, x5, p);
    t2 = sp_ntt_add(x2, x6, p);
    t6 = sp_ntt_sub(x2, x6, p);
    t3 = sp_ntt_add(x3, x7, p);
    t7 = sp_ntt_sub(x3, x7, p);

    a08 = sp_ntt_add(t0, t2, p);
    a09 = sp_ntt_sub(t0, t2, p);
    a10 = sp_ntt_add(t1, t3, p);
    a11 = sp_ntt_sub(t1, t3, p);
    a12 = t4;
    a13 = t6;
    a14 = sp_ntt_sub(t5, t7, p);
    a15 = sp_ntt_add(t5, t7, p); 
  }
  {
    sp_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_t x5 = in[ 1 * istride];
    sp_t x6 = in[ 6 * istride];
    sp_t x7 = in[11 * istride];
    sp_t x0 = in[16 * istride];
    sp_t x1 = in[21 * istride];
    sp_t x2 = in[26 * istride];
    sp_t x3 = in[31 * istride];
    sp_t x4 = in[36 * istride];

    t0 = sp_ntt_add(x0, x4, p);
    t4 = sp_ntt_sub(x0, x4, p);
    t1 = sp_ntt_add(x1, x5, p);
    t5 = sp_ntt_sub(x1, x5, p);
    t2 = sp_ntt_add(x2, x6, p);
    t6 = sp_ntt_sub(x2, x6, p);
    t3 = sp_ntt_add(x3, x7, p);
    t7 = sp_ntt_sub(x3, x7, p);

    a16 = sp_ntt_add(t0, t2, p);
    a17 = sp_ntt_sub(t0, t2, p);
    a18 = sp_ntt_add(t1, t3, p);
    a19 = sp_ntt_sub(t1, t3, p);
    a20 = t4;
    a21 = t6;
    a22 = sp_ntt_sub(t5, t7, p);
    a23 = sp_ntt_add(t5, t7, p); 
  }
  {
    sp_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_t x0 = in[24 * istride];
    sp_t x1 = in[29 * istride];
    sp_t x2 = in[34 * istride];
    sp_t x3 = in[39 * istride];
    sp_t x4 = in[ 4 * istride];
    sp_t x5 = in[ 9 * istride];
    sp_t x6 = in[14 * istride];
    sp_t x7 = in[19 * istride];

    t0 = sp_ntt_add(x0, x4, p);
    t4 = sp_ntt_sub(x0, x4, p);
    t1 = sp_ntt_add(x1, x5, p);
    t5 = sp_ntt_sub(x1, x5, p);
    t2 = sp_ntt_add(x2, x6, p);
    t6 = sp_ntt_sub(x2, x6, p);
    t3 = sp_ntt_add(x3, x7, p);
    t7 = sp_ntt_sub(x3, x7, p);

    a24 = sp_ntt_add(t0, t2, p);
    a25 = sp_ntt_sub(t0, t2, p);
    a26 = sp_ntt_add(t1, t3, p);
    a27 = sp_ntt_sub(t1, t3, p);
    a28 = t4;
    a29 = t6;
    a30 = sp_ntt_sub(t5, t7, p);
    a31 = sp_ntt_add(t5, t7, p); 
  }
  {
    sp_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_t x2 = in[ 2 * istride];
    sp_t x3 = in[ 7 * istride];
    sp_t x4 = in[12 * istride];
    sp_t x5 = in[17 * istride];
    sp_t x6 = in[22 * istride];
    sp_t x7 = in[27 * istride];
    sp_t x0 = in[32 * istride];
    sp_t x1 = in[37 * istride];

    t0 = sp_ntt_add(x0, x4, p);
    t4 = sp_ntt_sub(x0, x4, p);
    t1 = sp_ntt_add(x1, x5, p);
    t5 = sp_ntt_sub(x1, x5, p);
    t2 = sp_ntt_add(x2, x6, p);
    t6 = sp_ntt_sub(x2, x6, p);
    t3 = sp_ntt_add(x3, x7, p);
    t7 = sp_ntt_sub(x3, x7, p);

    a32 = sp_ntt_add(t0, t2, p);
    a33 = sp_ntt_sub(t0, t2, p);
    a34 = sp_ntt_add(t1, t3, p);
    a35 = sp_ntt_sub(t1, t3, p);
    a36 = t4;
    a37 = t6;
    a38 = sp_ntt_sub(t5, t7, p);
    a39 = sp_ntt_add(t5, t7, p); 
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a00;
    b1 = a08;
    b4 = a16;
    b2 = a24;
    b3 = a32;

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
    a08 = b4;
    a16 = b3;
    a24 = b1;
    a32 = b2;
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a01;
    b1 = a09;
    b4 = a17;
    b2 = a25;
    b3 = a33;

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
    a09 = b4;
    a17 = b3;
    a25 = b1;
    a33 = b2;
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a02;
    b1 = a10;
    b4 = a18;
    b2 = a26;
    b3 = a34;

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
    a10 = b4;
    a18 = b3;
    a26 = b1;
    a34 = b2;
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a03;
    b1 = a11;
    b4 = a19;
    b2 = a27;
    b3 = a35;

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

    b0 = sp_ntt_mul(b0, ntt_const[18], ntt_const[NC+18], p);
    b1 = sp_ntt_mul(b1, ntt_const[19], ntt_const[NC+19], p);
    b2 = sp_ntt_mul(b2, ntt_const[20], ntt_const[NC+20], p);
    b3 = sp_ntt_mul(b3, ntt_const[21], ntt_const[NC+21], p);
    b4 = sp_ntt_mul(b4, ntt_const[22], ntt_const[NC+22], p);
    b5 = sp_ntt_mul(b5, ntt_const[23], ntt_const[NC+23], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

    a03 = b0;
    a11 = b4;
    a19 = b3;
    a27 = b1;
    a35 = b2;
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a04;
    b1 = a12;
    b4 = a20;
    b2 = a28;
    b3 = a36;

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

    b1 = sp_ntt_mul(b1, ntt_const[25], ntt_const[NC+25], p);
    b2 = sp_ntt_mul(b2, ntt_const[26], ntt_const[NC+26], p);
    b3 = sp_ntt_mul(b3, ntt_const[27], ntt_const[NC+27], p);
    b4 = sp_ntt_mul(b4, ntt_const[28], ntt_const[NC+28], p);
    b5 = sp_ntt_mul(b5, ntt_const[29], ntt_const[NC+29], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

    a04 = b0;
    a12 = b4;
    a20 = b3;
    a28 = b1;
    a36 = b2;
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a05;
    b1 = a13;
    b4 = a21;
    b2 = a29;
    b3 = a37;

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

    b0 = sp_ntt_mul(b0, ntt_const[30], ntt_const[NC+30], p);
    b1 = sp_ntt_mul(b1, ntt_const[31], ntt_const[NC+31], p);
    b2 = sp_ntt_mul(b2, ntt_const[32], ntt_const[NC+32], p);
    b3 = sp_ntt_mul(b3, ntt_const[33], ntt_const[NC+33], p);
    b4 = sp_ntt_mul(b4, ntt_const[34], ntt_const[NC+34], p);
    b5 = sp_ntt_mul(b5, ntt_const[35], ntt_const[NC+35], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

    a05 = b0;
    a13 = b4;
    a21 = b3;
    a29 = b1;
    a37 = b2;
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a06;
    b1 = a14;
    b4 = a22;
    b2 = a30;
    b3 = a38;

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

    b0 = sp_ntt_mul(b0, ntt_const[36], ntt_const[NC+36], p);
    b1 = sp_ntt_mul(b1, ntt_const[37], ntt_const[NC+37], p);
    b2 = sp_ntt_mul(b2, ntt_const[38], ntt_const[NC+38], p);
    b3 = sp_ntt_mul(b3, ntt_const[39], ntt_const[NC+39], p);
    b4 = sp_ntt_mul(b4, ntt_const[40], ntt_const[NC+40], p);
    b5 = sp_ntt_mul(b5, ntt_const[41], ntt_const[NC+41], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

    a06 = b0;
    a14 = b4;
    a22 = b3;
    a30 = b1;
    a38 = b2;
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a07;
    b1 = a15;
    b4 = a23;
    b2 = a31;
    b3 = a39;

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

    b0 = sp_ntt_mul(b0, ntt_const[42], ntt_const[NC+42], p);
    b1 = sp_ntt_mul(b1, ntt_const[43], ntt_const[NC+43], p);
    b2 = sp_ntt_mul(b2, ntt_const[44], ntt_const[NC+44], p);
    b3 = sp_ntt_mul(b3, ntt_const[45], ntt_const[NC+45], p);
    b4 = sp_ntt_mul(b4, ntt_const[46], ntt_const[NC+46], p);
    b5 = sp_ntt_mul(b5, ntt_const[47], ntt_const[NC+47], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

    a07 = b0;
    a15 = b4;
    a23 = b3;
    a31 = b1;
    a39 = b2;
  }
  {
    sp_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_t p0 = a00;
    sp_t p1 = a01;
    sp_t p2 = a02;
    sp_t p3 = a03;
    sp_t p4 = a04;
    sp_t p5 = a05;
    sp_t p6 = a06;
    sp_t p7 = a07;

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

    t5 = sp_ntt_mul(t5, w[8], w[9], p);
    t2 = sp_ntt_mul(t2, w[18], w[19], p);
    t6 = sp_ntt_mul(t6, w[28], w[29], p);
    t1 = sp_ntt_mul(t1, w[38], w[39], p);
    t4 = sp_ntt_mul(t4, w[48], w[49], p);
    t3 = sp_ntt_mul(t3, w[58], w[59], p);
    t7 = sp_ntt_mul(t7, w[68], w[69], p);

    out[ 0 * ostride] = t0;
    out[ 5 * ostride] = t5;
    out[10 * ostride] = t2;
    out[15 * ostride] = t6;
    out[20 * ostride] = t1;
    out[25 * ostride] = t4;
    out[30 * ostride] = t3;
    out[35 * ostride] = t7;
  }
  {
    sp_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_t p0 = a08;
    sp_t p1 = a09;
    sp_t p2 = a10;
    sp_t p3 = a11;
    sp_t p4 = a12;
    sp_t p5 = a13;
    sp_t p6 = a14;
    sp_t p7 = a15;

    t0 = sp_ntt_add(p4, p5, p);
    t1 = sp_ntt_sub(p4, p5, p);
    t2 = sp_ntt_add(p6, p7, p);
    t3 = sp_ntt_sub(p6, p7, p);
    t4 = sp_ntt_add_partial(t0, t2, p);
    t5 = sp_ntt_sub_partial(t0, t2, p);
    t6 = sp_ntt_add_partial(t1, t3, p);
    t7 = sp_ntt_sub_partial(t1, t3, p);

    t0 = sp_ntt_add_partial(p0, p2, p);
    t1 = sp_ntt_sub_partial(p0, p2, p);
    t2 = sp_ntt_add_partial(p1, p3, p);
    t3 = sp_ntt_sub_partial(p1, p3, p);

    t4 = sp_ntt_mul(t4, w[0], w[1], p);
    t3 = sp_ntt_mul(t3, w[10], w[11], p);
    t7 = sp_ntt_mul(t7, w[20], w[21], p);
    t0 = sp_ntt_mul(t0, w[30], w[31], p);
    t5 = sp_ntt_mul(t5, w[40], w[41], p);
    t2 = sp_ntt_mul(t2, w[50], w[51], p);
    t6 = sp_ntt_mul(t6, w[60], w[61], p);
    t1 = sp_ntt_mul(t1, w[70], w[71], p);

    out[ 1 * ostride] = t4;
    out[ 6 * ostride] = t3;
    out[11 * ostride] = t7;
    out[16 * ostride] = t0;
    out[21 * ostride] = t5;
    out[26 * ostride] = t2;
    out[31 * ostride] = t6;
    out[36 * ostride] = t1;
  }
  {
    sp_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_t p0 = a16;
    sp_t p1 = a17;
    sp_t p2 = a18;
    sp_t p3 = a19;
    sp_t p4 = a20;
    sp_t p5 = a21;
    sp_t p6 = a22;
    sp_t p7 = a23;

    t0 = sp_ntt_add(p4, p5, p);
    t1 = sp_ntt_sub(p4, p5, p);
    t2 = sp_ntt_add(p6, p7, p);
    t3 = sp_ntt_sub(p6, p7, p);
    t4 = sp_ntt_add_partial(t0, t2, p);
    t5 = sp_ntt_sub_partial(t0, t2, p);
    t6 = sp_ntt_add_partial(t1, t3, p);
    t7 = sp_ntt_sub_partial(t1, t3, p);

    t0 = sp_ntt_add_partial(p0, p2, p);
    t1 = sp_ntt_sub_partial(p0, p2, p);
    t2 = sp_ntt_add_partial(p1, p3, p);
    t3 = sp_ntt_sub_partial(p1, p3, p);

    t2 = sp_ntt_mul(t2, w[2], w[3], p);
    t6 = sp_ntt_mul(t6, w[12], w[13], p);
    t1 = sp_ntt_mul(t1, w[22], w[23], p);
    t4 = sp_ntt_mul(t4, w[32], w[33], p);
    t3 = sp_ntt_mul(t3, w[42], w[43], p);
    t7 = sp_ntt_mul(t7, w[52], w[53], p);
    t0 = sp_ntt_mul(t0, w[62], w[63], p);
    t5 = sp_ntt_mul(t5, w[72], w[73], p);

    out[ 2 * ostride] = t2;
    out[ 7 * ostride] = t6;
    out[12 * ostride] = t1;
    out[17 * ostride] = t4;
    out[22 * ostride] = t3;
    out[27 * ostride] = t7;
    out[32 * ostride] = t0;
    out[37 * ostride] = t5;
  }
  {
    sp_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_t p0 = a24;
    sp_t p1 = a25;
    sp_t p2 = a26;
    sp_t p3 = a27;
    sp_t p4 = a28;
    sp_t p5 = a29;
    sp_t p6 = a30;
    sp_t p7 = a31;

    t0 = sp_ntt_add(p4, p5, p);
    t1 = sp_ntt_sub(p4, p5, p);
    t2 = sp_ntt_add(p6, p7, p);
    t3 = sp_ntt_sub(p6, p7, p);
    t4 = sp_ntt_add_partial(t0, t2, p);
    t5 = sp_ntt_sub_partial(t0, t2, p);
    t6 = sp_ntt_add_partial(t1, t3, p);
    t7 = sp_ntt_sub_partial(t1, t3, p);

    t0 = sp_ntt_add_partial(p0, p2, p);
    t1 = sp_ntt_sub_partial(p0, p2, p);
    t2 = sp_ntt_add_partial(p1, p3, p);
    t3 = sp_ntt_sub_partial(p1, p3, p);

    t7 = sp_ntt_mul(t7, w[4], w[5], p);
    t0 = sp_ntt_mul(t0, w[14], w[15], p);
    t5 = sp_ntt_mul(t5, w[24], w[25], p);
    t2 = sp_ntt_mul(t2, w[34], w[35], p);
    t6 = sp_ntt_mul(t6, w[44], w[45], p);
    t1 = sp_ntt_mul(t1, w[54], w[55], p);
    t4 = sp_ntt_mul(t4, w[64], w[65], p);
    t3 = sp_ntt_mul(t3, w[74], w[75], p);

    out[ 3 * ostride] = t7;
    out[ 8 * ostride] = t0;
    out[13 * ostride] = t5;
    out[18 * ostride] = t2;
    out[23 * ostride] = t6;
    out[28 * ostride] = t1;
    out[33 * ostride] = t4;
    out[38 * ostride] = t3;
  }
  {
    sp_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_t p0 = a32;
    sp_t p1 = a33;
    sp_t p2 = a34;
    sp_t p3 = a35;
    sp_t p4 = a36;
    sp_t p5 = a37;
    sp_t p6 = a38;
    sp_t p7 = a39;

    t0 = sp_ntt_add(p4, p5, p);
    t1 = sp_ntt_sub(p4, p5, p);
    t2 = sp_ntt_add(p6, p7, p);
    t3 = sp_ntt_sub(p6, p7, p);
    t4 = sp_ntt_add_partial(t0, t2, p);
    t5 = sp_ntt_sub_partial(t0, t2, p);
    t6 = sp_ntt_add_partial(t1, t3, p);
    t7 = sp_ntt_sub_partial(t1, t3, p);

    t0 = sp_ntt_add_partial(p0, p2, p);
    t1 = sp_ntt_sub_partial(p0, p2, p);
    t2 = sp_ntt_add_partial(p1, p3, p);
    t3 = sp_ntt_sub_partial(p1, p3, p);

    t1 = sp_ntt_mul(t1, w[6], w[7], p);
    t4 = sp_ntt_mul(t4, w[16], w[17], p);
    t3 = sp_ntt_mul(t3, w[26], w[27], p);
    t7 = sp_ntt_mul(t7, w[36], w[37], p);
    t0 = sp_ntt_mul(t0, w[46], w[47], p);
    t5 = sp_ntt_mul(t5, w[56], w[57], p);
    t2 = sp_ntt_mul(t2, w[66], w[67], p);
    t6 = sp_ntt_mul(t6, w[76], w[77], p);

    out[ 4 * ostride] = t1;
    out[ 9 * ostride] = t4;
    out[14 * ostride] = t3;
    out[19 * ostride] = t7;
    out[24 * ostride] = t0;
    out[29 * ostride] = t5;
    out[34 * ostride] = t2;
    out[39 * ostride] = t6;
  }
}

static void
ntt40_pfa_run_core(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t n,
	  sp_t p, spv_t ntt_const)
{
  spv_size_t j00, j01, j02, j03, j04, j05, j06, j07,
	     j08, j09, j10, j11, j12, j13, j14, j15,
	     j16, j17, j18, j19, j20, j21, j22, j23,
	     j24, j25, j26, j27, j28, j29, j30, j31,
	     j32, j33, j34, j35, j36, j37, j38, j39;

  sp_t a00, a01, a02, a03, a04, a05, a06, a07,
       a08, a09, a10, a11, a12, a13, a14, a15,
       a16, a17, a18, a19, a20, a21, a22, a23,
       a24, a25, a26, a27, a28, a29, a30, a31,
       a32, a33, a34, a35, a36, a37, a38, a39;

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
  j15 = sp_array_inc(j00, 15 * inc, n);
  j16 = sp_array_inc(j00, 16 * inc, n);
  j17 = sp_array_inc(j00, 17 * inc, n);
  j18 = sp_array_inc(j00, 18 * inc, n);
  j19 = sp_array_inc(j00, 19 * inc, n);
  j20 = sp_array_inc(j00, 20 * inc, n);
  j21 = sp_array_inc(j00, 21 * inc, n);
  j22 = sp_array_inc(j00, 22 * inc, n);
  j23 = sp_array_inc(j00, 23 * inc, n);
  j24 = sp_array_inc(j00, 24 * inc, n);
  j25 = sp_array_inc(j00, 25 * inc, n);
  j26 = sp_array_inc(j00, 26 * inc, n);
  j27 = sp_array_inc(j00, 27 * inc, n);
  j28 = sp_array_inc(j00, 28 * inc, n);
  j29 = sp_array_inc(j00, 29 * inc, n);
  j30 = sp_array_inc(j00, 30 * inc, n);
  j31 = sp_array_inc(j00, 31 * inc, n);
  j32 = sp_array_inc(j00, 32 * inc, n);
  j33 = sp_array_inc(j00, 33 * inc, n);
  j34 = sp_array_inc(j00, 34 * inc, n);
  j35 = sp_array_inc(j00, 35 * inc, n);
  j36 = sp_array_inc(j00, 36 * inc, n);
  j37 = sp_array_inc(j00, 37 * inc, n);
  j38 = sp_array_inc(j00, 38 * inc, n);
  j39 = sp_array_inc(j00, 39 * inc, n);

  {
    sp_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_t x0 = x[j00];
    sp_t x1 = x[j05];
    sp_t x2 = x[j10];
    sp_t x3 = x[j15];
    sp_t x4 = x[j20];
    sp_t x5 = x[j25];
    sp_t x6 = x[j30];
    sp_t x7 = x[j35];

    t0 = sp_ntt_add(x0, x4, p);
    t4 = sp_ntt_sub(x0, x4, p);
    t1 = sp_ntt_add(x1, x5, p);
    t5 = sp_ntt_sub(x1, x5, p);
    t2 = sp_ntt_add(x2, x6, p);
    t6 = sp_ntt_sub(x2, x6, p);
    t3 = sp_ntt_add(x3, x7, p);
    t7 = sp_ntt_sub(x3, x7, p);

    a00 = sp_ntt_add(t0, t2, p);
    a01 = sp_ntt_sub(t0, t2, p);
    a02 = sp_ntt_add(t1, t3, p);
    a03 = sp_ntt_sub(t1, t3, p);
    a04 = t4;
    a05 = t6;
    a06 = sp_ntt_sub(t5, t7, p);
    a07 = sp_ntt_add(t5, t7, p); 
  }
  {
    sp_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_t x7 = x[j03];
    sp_t x0 = x[j08];
    sp_t x1 = x[j13];
    sp_t x2 = x[j18];
    sp_t x3 = x[j23];
    sp_t x4 = x[j28];
    sp_t x5 = x[j33];
    sp_t x6 = x[j38];

    t0 = sp_ntt_add(x0, x4, p);
    t4 = sp_ntt_sub(x0, x4, p);
    t1 = sp_ntt_add(x1, x5, p);
    t5 = sp_ntt_sub(x1, x5, p);
    t2 = sp_ntt_add(x2, x6, p);
    t6 = sp_ntt_sub(x2, x6, p);
    t3 = sp_ntt_add(x3, x7, p);
    t7 = sp_ntt_sub(x3, x7, p);

    a08 = sp_ntt_add(t0, t2, p);
    a09 = sp_ntt_sub(t0, t2, p);
    a10 = sp_ntt_add(t1, t3, p);
    a11 = sp_ntt_sub(t1, t3, p);
    a12 = t4;
    a13 = t6;
    a14 = sp_ntt_sub(t5, t7, p);
    a15 = sp_ntt_add(t5, t7, p); 
  }
  {
    sp_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_t x5 = x[j01];
    sp_t x6 = x[j06];
    sp_t x7 = x[j11];
    sp_t x0 = x[j16];
    sp_t x1 = x[j21];
    sp_t x2 = x[j26];
    sp_t x3 = x[j31];
    sp_t x4 = x[j36];

    t0 = sp_ntt_add(x0, x4, p);
    t4 = sp_ntt_sub(x0, x4, p);
    t1 = sp_ntt_add(x1, x5, p);
    t5 = sp_ntt_sub(x1, x5, p);
    t2 = sp_ntt_add(x2, x6, p);
    t6 = sp_ntt_sub(x2, x6, p);
    t3 = sp_ntt_add(x3, x7, p);
    t7 = sp_ntt_sub(x3, x7, p);

    a16 = sp_ntt_add(t0, t2, p);
    a17 = sp_ntt_sub(t0, t2, p);
    a18 = sp_ntt_add(t1, t3, p);
    a19 = sp_ntt_sub(t1, t3, p);
    a20 = t4;
    a21 = t6;
    a22 = sp_ntt_sub(t5, t7, p);
    a23 = sp_ntt_add(t5, t7, p); 
  }
  {
    sp_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_t x0 = x[j24];
    sp_t x1 = x[j29];
    sp_t x2 = x[j34];
    sp_t x3 = x[j39];
    sp_t x4 = x[j04];
    sp_t x5 = x[j09];
    sp_t x6 = x[j14];
    sp_t x7 = x[j19];

    t0 = sp_ntt_add(x0, x4, p);
    t4 = sp_ntt_sub(x0, x4, p);
    t1 = sp_ntt_add(x1, x5, p);
    t5 = sp_ntt_sub(x1, x5, p);
    t2 = sp_ntt_add(x2, x6, p);
    t6 = sp_ntt_sub(x2, x6, p);
    t3 = sp_ntt_add(x3, x7, p);
    t7 = sp_ntt_sub(x3, x7, p);

    a24 = sp_ntt_add(t0, t2, p);
    a25 = sp_ntt_sub(t0, t2, p);
    a26 = sp_ntt_add(t1, t3, p);
    a27 = sp_ntt_sub(t1, t3, p);
    a28 = t4;
    a29 = t6;
    a30 = sp_ntt_sub(t5, t7, p);
    a31 = sp_ntt_add(t5, t7, p); 
  }
  {
    sp_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_t x2 = x[j02];
    sp_t x3 = x[j07];
    sp_t x4 = x[j12];
    sp_t x5 = x[j17];
    sp_t x6 = x[j22];
    sp_t x7 = x[j27];
    sp_t x0 = x[j32];
    sp_t x1 = x[j37];

    t0 = sp_ntt_add(x0, x4, p);
    t4 = sp_ntt_sub(x0, x4, p);
    t1 = sp_ntt_add(x1, x5, p);
    t5 = sp_ntt_sub(x1, x5, p);
    t2 = sp_ntt_add(x2, x6, p);
    t6 = sp_ntt_sub(x2, x6, p);
    t3 = sp_ntt_add(x3, x7, p);
    t7 = sp_ntt_sub(x3, x7, p);

    a32 = sp_ntt_add(t0, t2, p);
    a33 = sp_ntt_sub(t0, t2, p);
    a34 = sp_ntt_add(t1, t3, p);
    a35 = sp_ntt_sub(t1, t3, p);
    a36 = t4;
    a37 = t6;
    a38 = sp_ntt_sub(t5, t7, p);
    a39 = sp_ntt_add(t5, t7, p); 
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a00;
    b1 = a08;
    b4 = a16;
    b2 = a24;
    b3 = a32;

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
    a08 = b4;
    a16 = b3;
    a24 = b1;
    a32 = b2;
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a01;
    b1 = a09;
    b4 = a17;
    b2 = a25;
    b3 = a33;

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
    a09 = b4;
    a17 = b3;
    a25 = b1;
    a33 = b2;
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a02;
    b1 = a10;
    b4 = a18;
    b2 = a26;
    b3 = a34;

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
    a10 = b4;
    a18 = b3;
    a26 = b1;
    a34 = b2;
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a03;
    b1 = a11;
    b4 = a19;
    b2 = a27;
    b3 = a35;

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

    b0 = sp_ntt_mul(b0, ntt_const[18], ntt_const[NC+18], p);
    b1 = sp_ntt_mul(b1, ntt_const[19], ntt_const[NC+19], p);
    b2 = sp_ntt_mul(b2, ntt_const[20], ntt_const[NC+20], p);
    b3 = sp_ntt_mul(b3, ntt_const[21], ntt_const[NC+21], p);
    b4 = sp_ntt_mul(b4, ntt_const[22], ntt_const[NC+22], p);
    b5 = sp_ntt_mul(b5, ntt_const[23], ntt_const[NC+23], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

    a03 = b0;
    a11 = b4;
    a19 = b3;
    a27 = b1;
    a35 = b2;
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a04;
    b1 = a12;
    b4 = a20;
    b2 = a28;
    b3 = a36;

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

    b1 = sp_ntt_mul(b1, ntt_const[25], ntt_const[NC+25], p);
    b2 = sp_ntt_mul(b2, ntt_const[26], ntt_const[NC+26], p);
    b3 = sp_ntt_mul(b3, ntt_const[27], ntt_const[NC+27], p);
    b4 = sp_ntt_mul(b4, ntt_const[28], ntt_const[NC+28], p);
    b5 = sp_ntt_mul(b5, ntt_const[29], ntt_const[NC+29], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

    a04 = b0;
    a12 = b4;
    a20 = b3;
    a28 = b1;
    a36 = b2;
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a05;
    b1 = a13;
    b4 = a21;
    b2 = a29;
    b3 = a37;

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

    b0 = sp_ntt_mul(b0, ntt_const[30], ntt_const[NC+30], p);
    b1 = sp_ntt_mul(b1, ntt_const[31], ntt_const[NC+31], p);
    b2 = sp_ntt_mul(b2, ntt_const[32], ntt_const[NC+32], p);
    b3 = sp_ntt_mul(b3, ntt_const[33], ntt_const[NC+33], p);
    b4 = sp_ntt_mul(b4, ntt_const[34], ntt_const[NC+34], p);
    b5 = sp_ntt_mul(b5, ntt_const[35], ntt_const[NC+35], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

    a05 = b0;
    a13 = b4;
    a21 = b3;
    a29 = b1;
    a37 = b2;
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a06;
    b1 = a14;
    b4 = a22;
    b2 = a30;
    b3 = a38;

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

    b0 = sp_ntt_mul(b0, ntt_const[36], ntt_const[NC+36], p);
    b1 = sp_ntt_mul(b1, ntt_const[37], ntt_const[NC+37], p);
    b2 = sp_ntt_mul(b2, ntt_const[38], ntt_const[NC+38], p);
    b3 = sp_ntt_mul(b3, ntt_const[39], ntt_const[NC+39], p);
    b4 = sp_ntt_mul(b4, ntt_const[40], ntt_const[NC+40], p);
    b5 = sp_ntt_mul(b5, ntt_const[41], ntt_const[NC+41], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

    a06 = b0;
    a14 = b4;
    a22 = b3;
    a30 = b1;
    a38 = b2;
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a07;
    b1 = a15;
    b4 = a23;
    b2 = a31;
    b3 = a39;

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

    b0 = sp_ntt_mul(b0, ntt_const[42], ntt_const[NC+42], p);
    b1 = sp_ntt_mul(b1, ntt_const[43], ntt_const[NC+43], p);
    b2 = sp_ntt_mul(b2, ntt_const[44], ntt_const[NC+44], p);
    b3 = sp_ntt_mul(b3, ntt_const[45], ntt_const[NC+45], p);
    b4 = sp_ntt_mul(b4, ntt_const[46], ntt_const[NC+46], p);
    b5 = sp_ntt_mul(b5, ntt_const[47], ntt_const[NC+47], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

    a07 = b0;
    a15 = b4;
    a23 = b3;
    a31 = b1;
    a39 = b2;
  }
  {
    sp_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_t p0 = a00;
    sp_t p1 = a01;
    sp_t p2 = a02;
    sp_t p3 = a03;
    sp_t p4 = a04;
    sp_t p5 = a05;
    sp_t p6 = a06;
    sp_t p7 = a07;

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

    x[j00] = t0;
    x[j05] = t5;
    x[j10] = t2;
    x[j15] = t6;
    x[j20] = t1;
    x[j25] = t4;
    x[j30] = t3;
    x[j35] = t7;
  }
  {
    sp_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_t p0 = a08;
    sp_t p1 = a09;
    sp_t p2 = a10;
    sp_t p3 = a11;
    sp_t p4 = a12;
    sp_t p5 = a13;
    sp_t p6 = a14;
    sp_t p7 = a15;

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

    x[j01] = t4;
    x[j06] = t3;
    x[j11] = t7;
    x[j16] = t0;
    x[j21] = t5;
    x[j26] = t2;
    x[j31] = t6;
    x[j36] = t1;
  }
  {
    sp_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_t p0 = a16;
    sp_t p1 = a17;
    sp_t p2 = a18;
    sp_t p3 = a19;
    sp_t p4 = a20;
    sp_t p5 = a21;
    sp_t p6 = a22;
    sp_t p7 = a23;

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

    x[j02] = t2;
    x[j07] = t6;
    x[j12] = t1;
    x[j17] = t4;
    x[j22] = t3;
    x[j27] = t7;
    x[j32] = t0;
    x[j37] = t5;
  }
  {
    sp_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_t p0 = a24;
    sp_t p1 = a25;
    sp_t p2 = a26;
    sp_t p3 = a27;
    sp_t p4 = a28;
    sp_t p5 = a29;
    sp_t p6 = a30;
    sp_t p7 = a31;

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

    x[j03] = t7;
    x[j08] = t0;
    x[j13] = t5;
    x[j18] = t2;
    x[j23] = t6;
    x[j28] = t1;
    x[j33] = t4;
    x[j38] = t3;
  }
  {
    sp_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_t p0 = a32;
    sp_t p1 = a33;
    sp_t p2 = a34;
    sp_t p3 = a35;
    sp_t p4 = a36;
    sp_t p5 = a37;
    sp_t p6 = a38;
    sp_t p7 = a39;

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

    x[j04] = t1;
    x[j09] = t4;
    x[j14] = t3;
    x[j19] = t7;
    x[j24] = t0;
    x[j29] = t5;
    x[j34] = t2;
    x[j39] = t6;
  }

}

DECLARE_CORE_ROUTINES(40)
