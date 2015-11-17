#include "ntt/ntt-impl-simd.h"

#define NC 48

static const uint8_t ntt40_fixed_const[NC] = {1, 0, 0, 0, 0, 0,
					1, 0, 0, 0, 0, 0,
					1, 0, 0, 0, 0, 0,
					0, 0, 0, 0, 0, 0,
					1};


extern void X(ntt40_init)(spv_t out, sp_t p, sp_t d, sp_t primroot, 
                        sp_t order, sp_t perm);


static void
ntt40_run_core_simd(spv_t in, spv_size_t istride, spv_size_t idist,
		spv_t out, spv_size_t ostride, spv_size_t odist,
		sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  sp_simd_t a00, a01, a02, a03, a04, a05, a06, a07,
     	    a08, a09, a10, a11, a12, a13, a14, a15,
	    a16, a17, a18, a19, a20, a21, a22, a23,
	    a24, a25, a26, a27, a28, a29, a30, a31,
	    a32, a33, a34, a35, a36, a37, a38, a39;

  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x0 = sp_simd_gather(in +  0 * istride, idist, vsize);
    sp_simd_t x1 = sp_simd_gather(in +  5 * istride, idist, vsize);
    sp_simd_t x2 = sp_simd_gather(in + 10 * istride, idist, vsize);
    sp_simd_t x3 = sp_simd_gather(in + 15 * istride, idist, vsize);
    sp_simd_t x4 = sp_simd_gather(in + 20 * istride, idist, vsize);
    sp_simd_t x5 = sp_simd_gather(in + 25 * istride, idist, vsize);
    sp_simd_t x6 = sp_simd_gather(in + 30 * istride, idist, vsize);
    sp_simd_t x7 = sp_simd_gather(in + 35 * istride, idist, vsize);

    t0 = sp_ntt_add_simd(x0, x4, p);
    t4 = sp_ntt_sub_simd(x0, x4, p);
    t1 = sp_ntt_add_simd(x1, x5, p);
    t5 = sp_ntt_sub_simd(x1, x5, p);
    t2 = sp_ntt_add_simd(x2, x6, p);
    t6 = sp_ntt_sub_simd(x2, x6, p);
    t3 = sp_ntt_add_simd(x3, x7, p);
    t7 = sp_ntt_sub_simd(x3, x7, p);

    a00 = sp_ntt_add_simd(t0, t2, p);
    a01 = sp_ntt_sub_simd(t0, t2, p);
    a02 = sp_ntt_add_simd(t1, t3, p);
    a03 = sp_ntt_sub_simd(t1, t3, p);
    a04 = t4;
    a05 = t6;
    a06 = sp_ntt_sub_simd(t5, t7, p);
    a07 = sp_ntt_add_simd(t5, t7, p); 
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x7 = sp_simd_gather(in +  3 * istride, idist, vsize);
    sp_simd_t x0 = sp_simd_gather(in +  8 * istride, idist, vsize);
    sp_simd_t x1 = sp_simd_gather(in + 13 * istride, idist, vsize);
    sp_simd_t x2 = sp_simd_gather(in + 18 * istride, idist, vsize);
    sp_simd_t x3 = sp_simd_gather(in + 23 * istride, idist, vsize);
    sp_simd_t x4 = sp_simd_gather(in + 28 * istride, idist, vsize);
    sp_simd_t x5 = sp_simd_gather(in + 33 * istride, idist, vsize);
    sp_simd_t x6 = sp_simd_gather(in + 38 * istride, idist, vsize);

    t0 = sp_ntt_add_simd(x0, x4, p);
    t4 = sp_ntt_sub_simd(x0, x4, p);
    t1 = sp_ntt_add_simd(x1, x5, p);
    t5 = sp_ntt_sub_simd(x1, x5, p);
    t2 = sp_ntt_add_simd(x2, x6, p);
    t6 = sp_ntt_sub_simd(x2, x6, p);
    t3 = sp_ntt_add_simd(x3, x7, p);
    t7 = sp_ntt_sub_simd(x3, x7, p);

    a08 = sp_ntt_add_simd(t0, t2, p);
    a09 = sp_ntt_sub_simd(t0, t2, p);
    a10 = sp_ntt_add_simd(t1, t3, p);
    a11 = sp_ntt_sub_simd(t1, t3, p);
    a12 = t4;
    a13 = t6;
    a14 = sp_ntt_sub_simd(t5, t7, p);
    a15 = sp_ntt_add_simd(t5, t7, p); 
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x5 = sp_simd_gather(in +  1 * istride, idist, vsize);
    sp_simd_t x6 = sp_simd_gather(in +  6 * istride, idist, vsize);
    sp_simd_t x7 = sp_simd_gather(in + 11 * istride, idist, vsize);
    sp_simd_t x0 = sp_simd_gather(in + 16 * istride, idist, vsize);
    sp_simd_t x1 = sp_simd_gather(in + 21 * istride, idist, vsize);
    sp_simd_t x2 = sp_simd_gather(in + 26 * istride, idist, vsize);
    sp_simd_t x3 = sp_simd_gather(in + 31 * istride, idist, vsize);
    sp_simd_t x4 = sp_simd_gather(in + 36 * istride, idist, vsize);

    t0 = sp_ntt_add_simd(x0, x4, p);
    t4 = sp_ntt_sub_simd(x0, x4, p);
    t1 = sp_ntt_add_simd(x1, x5, p);
    t5 = sp_ntt_sub_simd(x1, x5, p);
    t2 = sp_ntt_add_simd(x2, x6, p);
    t6 = sp_ntt_sub_simd(x2, x6, p);
    t3 = sp_ntt_add_simd(x3, x7, p);
    t7 = sp_ntt_sub_simd(x3, x7, p);

    a16 = sp_ntt_add_simd(t0, t2, p);
    a17 = sp_ntt_sub_simd(t0, t2, p);
    a18 = sp_ntt_add_simd(t1, t3, p);
    a19 = sp_ntt_sub_simd(t1, t3, p);
    a20 = t4;
    a21 = t6;
    a22 = sp_ntt_sub_simd(t5, t7, p);
    a23 = sp_ntt_add_simd(t5, t7, p); 
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x0 = sp_simd_gather(in + 24 * istride, idist, vsize);
    sp_simd_t x1 = sp_simd_gather(in + 29 * istride, idist, vsize);
    sp_simd_t x2 = sp_simd_gather(in + 34 * istride, idist, vsize);
    sp_simd_t x3 = sp_simd_gather(in + 39 * istride, idist, vsize);
    sp_simd_t x4 = sp_simd_gather(in +  4 * istride, idist, vsize);
    sp_simd_t x5 = sp_simd_gather(in +  9 * istride, idist, vsize);
    sp_simd_t x6 = sp_simd_gather(in + 14 * istride, idist, vsize);
    sp_simd_t x7 = sp_simd_gather(in + 19 * istride, idist, vsize);

    t0 = sp_ntt_add_simd(x0, x4, p);
    t4 = sp_ntt_sub_simd(x0, x4, p);
    t1 = sp_ntt_add_simd(x1, x5, p);
    t5 = sp_ntt_sub_simd(x1, x5, p);
    t2 = sp_ntt_add_simd(x2, x6, p);
    t6 = sp_ntt_sub_simd(x2, x6, p);
    t3 = sp_ntt_add_simd(x3, x7, p);
    t7 = sp_ntt_sub_simd(x3, x7, p);

    a24 = sp_ntt_add_simd(t0, t2, p);
    a25 = sp_ntt_sub_simd(t0, t2, p);
    a26 = sp_ntt_add_simd(t1, t3, p);
    a27 = sp_ntt_sub_simd(t1, t3, p);
    a28 = t4;
    a29 = t6;
    a30 = sp_ntt_sub_simd(t5, t7, p);
    a31 = sp_ntt_add_simd(t5, t7, p); 
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x2 = sp_simd_gather(in +  2 * istride, idist, vsize);
    sp_simd_t x3 = sp_simd_gather(in +  7 * istride, idist, vsize);
    sp_simd_t x4 = sp_simd_gather(in + 12 * istride, idist, vsize);
    sp_simd_t x5 = sp_simd_gather(in + 17 * istride, idist, vsize);
    sp_simd_t x6 = sp_simd_gather(in + 22 * istride, idist, vsize);
    sp_simd_t x7 = sp_simd_gather(in + 27 * istride, idist, vsize);
    sp_simd_t x0 = sp_simd_gather(in + 32 * istride, idist, vsize);
    sp_simd_t x1 = sp_simd_gather(in + 37 * istride, idist, vsize);

    t0 = sp_ntt_add_simd(x0, x4, p);
    t4 = sp_ntt_sub_simd(x0, x4, p);
    t1 = sp_ntt_add_simd(x1, x5, p);
    t5 = sp_ntt_sub_simd(x1, x5, p);
    t2 = sp_ntt_add_simd(x2, x6, p);
    t6 = sp_ntt_sub_simd(x2, x6, p);
    t3 = sp_ntt_add_simd(x3, x7, p);
    t7 = sp_ntt_sub_simd(x3, x7, p);

    a32 = sp_ntt_add_simd(t0, t2, p);
    a33 = sp_ntt_sub_simd(t0, t2, p);
    a34 = sp_ntt_add_simd(t1, t3, p);
    a35 = sp_ntt_sub_simd(t1, t3, p);
    a36 = t4;
    a37 = t6;
    a38 = sp_ntt_sub_simd(t5, t7, p);
    a39 = sp_ntt_add_simd(t5, t7, p); 
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a00;
    b1 = a08;
    b4 = a16;
    b2 = a24;
    b3 = a32;

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
    a08 = b4;
    a16 = b3;
    a24 = b1;
    a32 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a01;
    b1 = a09;
    b4 = a17;
    b2 = a25;
    b3 = a33;

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
    a09 = b4;
    a17 = b3;
    a25 = b1;
    a33 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a02;
    b1 = a10;
    b4 = a18;
    b2 = a26;
    b3 = a34;

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
    a10 = b4;
    a18 = b3;
    a26 = b1;
    a34 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a03;
    b1 = a11;
    b4 = a19;
    b2 = a27;
    b3 = a35;

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

    b0 = sp_ntt_mul_simd(b0, ntt_const[18], ntt_const[NC+18], p);
    b1 = sp_ntt_mul_simd(b1, ntt_const[19], ntt_const[NC+19], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[20], ntt_const[NC+20], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[21], ntt_const[NC+21], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[22], ntt_const[NC+22], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[23], ntt_const[NC+23], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a03 = b0;
    a11 = b4;
    a19 = b3;
    a27 = b1;
    a35 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a04;
    b1 = a12;
    b4 = a20;
    b2 = a28;
    b3 = a36;

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

    b1 = sp_ntt_mul_simd(b1, ntt_const[25], ntt_const[NC+25], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[26], ntt_const[NC+26], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[27], ntt_const[NC+27], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[28], ntt_const[NC+28], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[29], ntt_const[NC+29], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a04 = b0;
    a12 = b4;
    a20 = b3;
    a28 = b1;
    a36 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a05;
    b1 = a13;
    b4 = a21;
    b2 = a29;
    b3 = a37;

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

    b0 = sp_ntt_mul_simd(b0, ntt_const[30], ntt_const[NC+30], p);
    b1 = sp_ntt_mul_simd(b1, ntt_const[31], ntt_const[NC+31], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[32], ntt_const[NC+32], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[33], ntt_const[NC+33], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[34], ntt_const[NC+34], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[35], ntt_const[NC+35], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a05 = b0;
    a13 = b4;
    a21 = b3;
    a29 = b1;
    a37 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a06;
    b1 = a14;
    b4 = a22;
    b2 = a30;
    b3 = a38;

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

    b0 = sp_ntt_mul_simd(b0, ntt_const[36], ntt_const[NC+36], p);
    b1 = sp_ntt_mul_simd(b1, ntt_const[37], ntt_const[NC+37], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[38], ntt_const[NC+38], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[39], ntt_const[NC+39], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[40], ntt_const[NC+40], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[41], ntt_const[NC+41], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a06 = b0;
    a14 = b4;
    a22 = b3;
    a30 = b1;
    a38 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a07;
    b1 = a15;
    b4 = a23;
    b2 = a31;
    b3 = a39;

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

    b0 = sp_ntt_mul_simd(b0, ntt_const[42], ntt_const[NC+42], p);
    b1 = sp_ntt_mul_simd(b1, ntt_const[43], ntt_const[NC+43], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[44], ntt_const[NC+44], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[45], ntt_const[NC+45], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[46], ntt_const[NC+46], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[47], ntt_const[NC+47], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a07 = b0;
    a15 = b4;
    a23 = b3;
    a31 = b1;
    a39 = b2;
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t p0 = a00;
    sp_simd_t p1 = a01;
    sp_simd_t p2 = a02;
    sp_simd_t p3 = a03;
    sp_simd_t p4 = a04;
    sp_simd_t p5 = a05;
    sp_simd_t p6 = a06;
    sp_simd_t p7 = a07;

    t0 = sp_ntt_add_simd(p4, p5, p);
    t1 = sp_ntt_sub_simd(p4, p5, p);
    t2 = sp_ntt_add_simd(p6, p7, p);
    t3 = sp_ntt_sub_simd(p6, p7, p);
    t4 = sp_ntt_add_simd(t0, t2, p);
    t5 = sp_ntt_sub_simd(t0, t2, p);
    t6 = sp_ntt_add_simd(t1, t3, p);
    t7 = sp_ntt_sub_simd(t1, t3, p);

    t0 = sp_ntt_add_simd(p0, p2, p);
    t1 = sp_ntt_sub_simd(p0, p2, p);
    t2 = sp_ntt_add_simd(p1, p3, p);
    t3 = sp_ntt_sub_simd(p1, p3, p);

    sp_simd_scatter(t0, out +  0 * ostride, odist, vsize);
    sp_simd_scatter(t5, out +  5 * ostride, odist, vsize);
    sp_simd_scatter(t2, out + 10 * ostride, odist, vsize);
    sp_simd_scatter(t6, out + 15 * ostride, odist, vsize);
    sp_simd_scatter(t1, out + 20 * ostride, odist, vsize);
    sp_simd_scatter(t4, out + 25 * ostride, odist, vsize);
    sp_simd_scatter(t3, out + 30 * ostride, odist, vsize);
    sp_simd_scatter(t7, out + 35 * ostride, odist, vsize);
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t p0 = a08;
    sp_simd_t p1 = a09;
    sp_simd_t p2 = a10;
    sp_simd_t p3 = a11;
    sp_simd_t p4 = a12;
    sp_simd_t p5 = a13;
    sp_simd_t p6 = a14;
    sp_simd_t p7 = a15;

    t0 = sp_ntt_add_simd(p4, p5, p);
    t1 = sp_ntt_sub_simd(p4, p5, p);
    t2 = sp_ntt_add_simd(p6, p7, p);
    t3 = sp_ntt_sub_simd(p6, p7, p);
    t4 = sp_ntt_add_simd(t0, t2, p);
    t5 = sp_ntt_sub_simd(t0, t2, p);
    t6 = sp_ntt_add_simd(t1, t3, p);
    t7 = sp_ntt_sub_simd(t1, t3, p);

    t0 = sp_ntt_add_simd(p0, p2, p);
    t1 = sp_ntt_sub_simd(p0, p2, p);
    t2 = sp_ntt_add_simd(p1, p3, p);
    t3 = sp_ntt_sub_simd(p1, p3, p);

    sp_simd_scatter(t4, out +  1 * ostride, odist, vsize);
    sp_simd_scatter(t3, out +  6 * ostride, odist, vsize);
    sp_simd_scatter(t7, out + 11 * ostride, odist, vsize);
    sp_simd_scatter(t0, out + 16 * ostride, odist, vsize);
    sp_simd_scatter(t5, out + 21 * ostride, odist, vsize);
    sp_simd_scatter(t2, out + 26 * ostride, odist, vsize);
    sp_simd_scatter(t6, out + 31 * ostride, odist, vsize);
    sp_simd_scatter(t1, out + 36 * ostride, odist, vsize);
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t p0 = a16;
    sp_simd_t p1 = a17;
    sp_simd_t p2 = a18;
    sp_simd_t p3 = a19;
    sp_simd_t p4 = a20;
    sp_simd_t p5 = a21;
    sp_simd_t p6 = a22;
    sp_simd_t p7 = a23;

    t0 = sp_ntt_add_simd(p4, p5, p);
    t1 = sp_ntt_sub_simd(p4, p5, p);
    t2 = sp_ntt_add_simd(p6, p7, p);
    t3 = sp_ntt_sub_simd(p6, p7, p);
    t4 = sp_ntt_add_simd(t0, t2, p);
    t5 = sp_ntt_sub_simd(t0, t2, p);
    t6 = sp_ntt_add_simd(t1, t3, p);
    t7 = sp_ntt_sub_simd(t1, t3, p);

    t0 = sp_ntt_add_simd(p0, p2, p);
    t1 = sp_ntt_sub_simd(p0, p2, p);
    t2 = sp_ntt_add_simd(p1, p3, p);
    t3 = sp_ntt_sub_simd(p1, p3, p);

    sp_simd_scatter(t2, out +  2 * ostride, odist, vsize);
    sp_simd_scatter(t6, out +  7 * ostride, odist, vsize);
    sp_simd_scatter(t1, out + 12 * ostride, odist, vsize);
    sp_simd_scatter(t4, out + 17 * ostride, odist, vsize);
    sp_simd_scatter(t3, out + 22 * ostride, odist, vsize);
    sp_simd_scatter(t7, out + 27 * ostride, odist, vsize);
    sp_simd_scatter(t0, out + 32 * ostride, odist, vsize);
    sp_simd_scatter(t5, out + 37 * ostride, odist, vsize);
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t p0 = a24;
    sp_simd_t p1 = a25;
    sp_simd_t p2 = a26;
    sp_simd_t p3 = a27;
    sp_simd_t p4 = a28;
    sp_simd_t p5 = a29;
    sp_simd_t p6 = a30;
    sp_simd_t p7 = a31;

    t0 = sp_ntt_add_simd(p4, p5, p);
    t1 = sp_ntt_sub_simd(p4, p5, p);
    t2 = sp_ntt_add_simd(p6, p7, p);
    t3 = sp_ntt_sub_simd(p6, p7, p);
    t4 = sp_ntt_add_simd(t0, t2, p);
    t5 = sp_ntt_sub_simd(t0, t2, p);
    t6 = sp_ntt_add_simd(t1, t3, p);
    t7 = sp_ntt_sub_simd(t1, t3, p);

    t0 = sp_ntt_add_simd(p0, p2, p);
    t1 = sp_ntt_sub_simd(p0, p2, p);
    t2 = sp_ntt_add_simd(p1, p3, p);
    t3 = sp_ntt_sub_simd(p1, p3, p);

    sp_simd_scatter(t7, out +  3 * ostride, odist, vsize);
    sp_simd_scatter(t0, out +  8 * ostride, odist, vsize);
    sp_simd_scatter(t5, out + 13 * ostride, odist, vsize);
    sp_simd_scatter(t2, out + 18 * ostride, odist, vsize);
    sp_simd_scatter(t6, out + 23 * ostride, odist, vsize);
    sp_simd_scatter(t1, out + 28 * ostride, odist, vsize);
    sp_simd_scatter(t4, out + 33 * ostride, odist, vsize);
    sp_simd_scatter(t3, out + 38 * ostride, odist, vsize);
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t p0 = a32;
    sp_simd_t p1 = a33;
    sp_simd_t p2 = a34;
    sp_simd_t p3 = a35;
    sp_simd_t p4 = a36;
    sp_simd_t p5 = a37;
    sp_simd_t p6 = a38;
    sp_simd_t p7 = a39;

    t0 = sp_ntt_add_simd(p4, p5, p);
    t1 = sp_ntt_sub_simd(p4, p5, p);
    t2 = sp_ntt_add_simd(p6, p7, p);
    t3 = sp_ntt_sub_simd(p6, p7, p);
    t4 = sp_ntt_add_simd(t0, t2, p);
    t5 = sp_ntt_sub_simd(t0, t2, p);
    t6 = sp_ntt_add_simd(t1, t3, p);
    t7 = sp_ntt_sub_simd(t1, t3, p);

    t0 = sp_ntt_add_simd(p0, p2, p);
    t1 = sp_ntt_sub_simd(p0, p2, p);
    t2 = sp_ntt_add_simd(p1, p3, p);
    t3 = sp_ntt_sub_simd(p1, p3, p);

    sp_simd_scatter(t1, out +  4 * ostride, odist, vsize);
    sp_simd_scatter(t4, out +  9 * ostride, odist, vsize);
    sp_simd_scatter(t3, out + 14 * ostride, odist, vsize);
    sp_simd_scatter(t7, out + 19 * ostride, odist, vsize);
    sp_simd_scatter(t0, out + 24 * ostride, odist, vsize);
    sp_simd_scatter(t5, out + 29 * ostride, odist, vsize);
    sp_simd_scatter(t2, out + 34 * ostride, odist, vsize);
    sp_simd_scatter(t6, out + 39 * ostride, odist, vsize);
  }
}

static void
ntt40_run_core_simd_interleaved(
		spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		sp_simd_t p, sp_simd_t * c)
{
  sp_simd_t a00, a01, a02, a03, a04, a05, a06, a07,
     	    a08, a09, a10, a11, a12, a13, a14, a15,
	    a16, a17, a18, a19, a20, a21, a22, a23,
	    a24, a25, a26, a27, a28, a29, a30, a31,
	    a32, a33, a34, a35, a36, a37, a38, a39;

  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x0 = sp_simd_load(in +  0 * istride);
    sp_simd_t x1 = sp_simd_load(in +  5 * istride);
    sp_simd_t x2 = sp_simd_load(in + 10 * istride);
    sp_simd_t x3 = sp_simd_load(in + 15 * istride);
    sp_simd_t x4 = sp_simd_load(in + 20 * istride);
    sp_simd_t x5 = sp_simd_load(in + 25 * istride);
    sp_simd_t x6 = sp_simd_load(in + 30 * istride);
    sp_simd_t x7 = sp_simd_load(in + 35 * istride);

    t0 = sp_ntt_add_simd0(x0, x4, p);
    t4 = sp_ntt_sub_simd0(x0, x4, p);
    t1 = sp_ntt_add_simd0(x1, x5, p);
    t5 = sp_ntt_sub_simd0(x1, x5, p);
    t2 = sp_ntt_add_simd0(x2, x6, p);
    t6 = sp_ntt_sub_simd0(x2, x6, p);
    t3 = sp_ntt_add_simd0(x3, x7, p);
    t7 = sp_ntt_sub_simd0(x3, x7, p);

    a00 = sp_ntt_add_simd0(t0, t2, p);
    a01 = sp_ntt_sub_simd0(t0, t2, p);
    a02 = sp_ntt_add_simd0(t1, t3, p);
    a03 = sp_ntt_sub_simd0(t1, t3, p);
    a04 = t4;
    a05 = t6;
    a06 = sp_ntt_sub_simd0(t5, t7, p);
    a07 = sp_ntt_add_simd0(t5, t7, p); 
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x7 = sp_simd_load(in +  3 * istride);
    sp_simd_t x0 = sp_simd_load(in +  8 * istride);
    sp_simd_t x1 = sp_simd_load(in + 13 * istride);
    sp_simd_t x2 = sp_simd_load(in + 18 * istride);
    sp_simd_t x3 = sp_simd_load(in + 23 * istride);
    sp_simd_t x4 = sp_simd_load(in + 28 * istride);
    sp_simd_t x5 = sp_simd_load(in + 33 * istride);
    sp_simd_t x6 = sp_simd_load(in + 38 * istride);

    t0 = sp_ntt_add_simd0(x0, x4, p);
    t4 = sp_ntt_sub_simd0(x0, x4, p);
    t1 = sp_ntt_add_simd0(x1, x5, p);
    t5 = sp_ntt_sub_simd0(x1, x5, p);
    t2 = sp_ntt_add_simd0(x2, x6, p);
    t6 = sp_ntt_sub_simd0(x2, x6, p);
    t3 = sp_ntt_add_simd0(x3, x7, p);
    t7 = sp_ntt_sub_simd0(x3, x7, p);

    a08 = sp_ntt_add_simd0(t0, t2, p);
    a09 = sp_ntt_sub_simd0(t0, t2, p);
    a10 = sp_ntt_add_simd0(t1, t3, p);
    a11 = sp_ntt_sub_simd0(t1, t3, p);
    a12 = t4;
    a13 = t6;
    a14 = sp_ntt_sub_simd0(t5, t7, p);
    a15 = sp_ntt_add_simd0(t5, t7, p); 
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x5 = sp_simd_load(in +  1 * istride);
    sp_simd_t x6 = sp_simd_load(in +  6 * istride);
    sp_simd_t x7 = sp_simd_load(in + 11 * istride);
    sp_simd_t x0 = sp_simd_load(in + 16 * istride);
    sp_simd_t x1 = sp_simd_load(in + 21 * istride);
    sp_simd_t x2 = sp_simd_load(in + 26 * istride);
    sp_simd_t x3 = sp_simd_load(in + 31 * istride);
    sp_simd_t x4 = sp_simd_load(in + 36 * istride);

    t0 = sp_ntt_add_simd0(x0, x4, p);
    t4 = sp_ntt_sub_simd0(x0, x4, p);
    t1 = sp_ntt_add_simd0(x1, x5, p);
    t5 = sp_ntt_sub_simd0(x1, x5, p);
    t2 = sp_ntt_add_simd0(x2, x6, p);
    t6 = sp_ntt_sub_simd0(x2, x6, p);
    t3 = sp_ntt_add_simd0(x3, x7, p);
    t7 = sp_ntt_sub_simd0(x3, x7, p);

    a16 = sp_ntt_add_simd0(t0, t2, p);
    a17 = sp_ntt_sub_simd0(t0, t2, p);
    a18 = sp_ntt_add_simd0(t1, t3, p);
    a19 = sp_ntt_sub_simd0(t1, t3, p);
    a20 = t4;
    a21 = t6;
    a22 = sp_ntt_sub_simd0(t5, t7, p);
    a23 = sp_ntt_add_simd0(t5, t7, p); 
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x0 = sp_simd_load(in + 24 * istride);
    sp_simd_t x1 = sp_simd_load(in + 29 * istride);
    sp_simd_t x2 = sp_simd_load(in + 34 * istride);
    sp_simd_t x3 = sp_simd_load(in + 39 * istride);
    sp_simd_t x4 = sp_simd_load(in +  4 * istride);
    sp_simd_t x5 = sp_simd_load(in +  9 * istride);
    sp_simd_t x6 = sp_simd_load(in + 14 * istride);
    sp_simd_t x7 = sp_simd_load(in + 19 * istride);

    t0 = sp_ntt_add_simd0(x0, x4, p);
    t4 = sp_ntt_sub_simd0(x0, x4, p);
    t1 = sp_ntt_add_simd0(x1, x5, p);
    t5 = sp_ntt_sub_simd0(x1, x5, p);
    t2 = sp_ntt_add_simd0(x2, x6, p);
    t6 = sp_ntt_sub_simd0(x2, x6, p);
    t3 = sp_ntt_add_simd0(x3, x7, p);
    t7 = sp_ntt_sub_simd0(x3, x7, p);

    a24 = sp_ntt_add_simd0(t0, t2, p);
    a25 = sp_ntt_sub_simd0(t0, t2, p);
    a26 = sp_ntt_add_simd0(t1, t3, p);
    a27 = sp_ntt_sub_simd0(t1, t3, p);
    a28 = t4;
    a29 = t6;
    a30 = sp_ntt_sub_simd0(t5, t7, p);
    a31 = sp_ntt_add_simd0(t5, t7, p); 
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x2 = sp_simd_load(in +  2 * istride);
    sp_simd_t x3 = sp_simd_load(in +  7 * istride);
    sp_simd_t x4 = sp_simd_load(in + 12 * istride);
    sp_simd_t x5 = sp_simd_load(in + 17 * istride);
    sp_simd_t x6 = sp_simd_load(in + 22 * istride);
    sp_simd_t x7 = sp_simd_load(in + 27 * istride);
    sp_simd_t x0 = sp_simd_load(in + 32 * istride);
    sp_simd_t x1 = sp_simd_load(in + 37 * istride);

    t0 = sp_ntt_add_simd0(x0, x4, p);
    t4 = sp_ntt_sub_simd0(x0, x4, p);
    t1 = sp_ntt_add_simd0(x1, x5, p);
    t5 = sp_ntt_sub_simd0(x1, x5, p);
    t2 = sp_ntt_add_simd0(x2, x6, p);
    t6 = sp_ntt_sub_simd0(x2, x6, p);
    t3 = sp_ntt_add_simd0(x3, x7, p);
    t7 = sp_ntt_sub_simd0(x3, x7, p);

    a32 = sp_ntt_add_simd0(t0, t2, p);
    a33 = sp_ntt_sub_simd0(t0, t2, p);
    a34 = sp_ntt_add_simd0(t1, t3, p);
    a35 = sp_ntt_sub_simd0(t1, t3, p);
    a36 = t4;
    a37 = t6;
    a38 = sp_ntt_sub_simd0(t5, t7, p);
    a39 = sp_ntt_add_simd0(t5, t7, p); 
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a00;
    b1 = a08;
    b4 = a16;
    b2 = a24;
    b3 = a32;

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
    a08 = b4;
    a16 = b3;
    a24 = b1;
    a32 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a01;
    b1 = a09;
    b4 = a17;
    b2 = a25;
    b3 = a33;

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
    a09 = b4;
    a17 = b3;
    a25 = b1;
    a33 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a02;
    b1 = a10;
    b4 = a18;
    b2 = a26;
    b3 = a34;

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
    a10 = b4;
    a18 = b3;
    a26 = b1;
    a34 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a03;
    b1 = a11;
    b4 = a19;
    b2 = a27;
    b3 = a35;

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

    b0 = sp_ntt_mul_simd0(b0, c+2*18, p);
    b1 = sp_ntt_mul_simd0(b1, c+2*19, p);
    b2 = sp_ntt_mul_simd0(b2, c+2*20, p);
    b3 = sp_ntt_mul_simd0(b3, c+2*21, p);
    b4 = sp_ntt_mul_simd0(b4, c+2*22, p);
    b5 = sp_ntt_mul_simd0(b5, c+2*23, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a03 = b0;
    a11 = b4;
    a19 = b3;
    a27 = b1;
    a35 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a04;
    b1 = a12;
    b4 = a20;
    b2 = a28;
    b3 = a36;

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

    b1 = sp_ntt_mul_simd0(b1, c+2*25, p);
    b2 = sp_ntt_mul_simd0(b2, c+2*26, p);
    b3 = sp_ntt_mul_simd0(b3, c+2*27, p);
    b4 = sp_ntt_mul_simd0(b4, c+2*28, p);
    b5 = sp_ntt_mul_simd0(b5, c+2*29, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a04 = b0;
    a12 = b4;
    a20 = b3;
    a28 = b1;
    a36 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a05;
    b1 = a13;
    b4 = a21;
    b2 = a29;
    b3 = a37;

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

    b0 = sp_ntt_mul_simd0(b0, c+2*30, p);
    b1 = sp_ntt_mul_simd0(b1, c+2*31, p);
    b2 = sp_ntt_mul_simd0(b2, c+2*32, p);
    b3 = sp_ntt_mul_simd0(b3, c+2*33, p);
    b4 = sp_ntt_mul_simd0(b4, c+2*34, p);
    b5 = sp_ntt_mul_simd0(b5, c+2*35, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a05 = b0;
    a13 = b4;
    a21 = b3;
    a29 = b1;
    a37 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a06;
    b1 = a14;
    b4 = a22;
    b2 = a30;
    b3 = a38;

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

    b0 = sp_ntt_mul_simd0(b0, c+2*36, p);
    b1 = sp_ntt_mul_simd0(b1, c+2*37, p);
    b2 = sp_ntt_mul_simd0(b2, c+2*38, p);
    b3 = sp_ntt_mul_simd0(b3, c+2*39, p);
    b4 = sp_ntt_mul_simd0(b4, c+2*40, p);
    b5 = sp_ntt_mul_simd0(b5, c+2*41, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a06 = b0;
    a14 = b4;
    a22 = b3;
    a30 = b1;
    a38 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a07;
    b1 = a15;
    b4 = a23;
    b2 = a31;
    b3 = a39;

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

    b0 = sp_ntt_mul_simd0(b0, c+2*42, p);
    b1 = sp_ntt_mul_simd0(b1, c+2*43, p);
    b2 = sp_ntt_mul_simd0(b2, c+2*44, p);
    b3 = sp_ntt_mul_simd0(b3, c+2*45, p);
    b4 = sp_ntt_mul_simd0(b4, c+2*46, p);
    b5 = sp_ntt_mul_simd0(b5, c+2*47, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a07 = b0;
    a15 = b4;
    a23 = b3;
    a31 = b1;
    a39 = b2;
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t p0 = a00;
    sp_simd_t p1 = a01;
    sp_simd_t p2 = a02;
    sp_simd_t p3 = a03;
    sp_simd_t p4 = a04;
    sp_simd_t p5 = a05;
    sp_simd_t p6 = a06;
    sp_simd_t p7 = a07;

    t0 = sp_ntt_add_simd0(p4, p5, p);
    t1 = sp_ntt_sub_simd0(p4, p5, p);
    t2 = sp_ntt_add_simd0(p6, p7, p);
    t3 = sp_ntt_sub_simd0(p6, p7, p);
    t4 = sp_ntt_add_simd0(t0, t2, p);
    t5 = sp_ntt_sub_simd0(t0, t2, p);
    t6 = sp_ntt_add_simd0(t1, t3, p);
    t7 = sp_ntt_sub_simd0(t1, t3, p);

    t0 = sp_ntt_add_simd0(p0, p2, p);
    t1 = sp_ntt_sub_simd0(p0, p2, p);
    t2 = sp_ntt_add_simd0(p1, p3, p);
    t3 = sp_ntt_sub_simd0(p1, p3, p);

    sp_simd_store(t0, out +  0 * ostride);
    sp_simd_store(t5, out +  5 * ostride);
    sp_simd_store(t2, out + 10 * ostride);
    sp_simd_store(t6, out + 15 * ostride);
    sp_simd_store(t1, out + 20 * ostride);
    sp_simd_store(t4, out + 25 * ostride);
    sp_simd_store(t3, out + 30 * ostride);
    sp_simd_store(t7, out + 35 * ostride);
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t p0 = a08;
    sp_simd_t p1 = a09;
    sp_simd_t p2 = a10;
    sp_simd_t p3 = a11;
    sp_simd_t p4 = a12;
    sp_simd_t p5 = a13;
    sp_simd_t p6 = a14;
    sp_simd_t p7 = a15;

    t0 = sp_ntt_add_simd0(p4, p5, p);
    t1 = sp_ntt_sub_simd0(p4, p5, p);
    t2 = sp_ntt_add_simd0(p6, p7, p);
    t3 = sp_ntt_sub_simd0(p6, p7, p);
    t4 = sp_ntt_add_simd0(t0, t2, p);
    t5 = sp_ntt_sub_simd0(t0, t2, p);
    t6 = sp_ntt_add_simd0(t1, t3, p);
    t7 = sp_ntt_sub_simd0(t1, t3, p);

    t0 = sp_ntt_add_simd0(p0, p2, p);
    t1 = sp_ntt_sub_simd0(p0, p2, p);
    t2 = sp_ntt_add_simd0(p1, p3, p);
    t3 = sp_ntt_sub_simd0(p1, p3, p);

    sp_simd_store(t4, out +  1 * ostride);
    sp_simd_store(t3, out +  6 * ostride);
    sp_simd_store(t7, out + 11 * ostride);
    sp_simd_store(t0, out + 16 * ostride);
    sp_simd_store(t5, out + 21 * ostride);
    sp_simd_store(t2, out + 26 * ostride);
    sp_simd_store(t6, out + 31 * ostride);
    sp_simd_store(t1, out + 36 * ostride);
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t p0 = a16;
    sp_simd_t p1 = a17;
    sp_simd_t p2 = a18;
    sp_simd_t p3 = a19;
    sp_simd_t p4 = a20;
    sp_simd_t p5 = a21;
    sp_simd_t p6 = a22;
    sp_simd_t p7 = a23;

    t0 = sp_ntt_add_simd0(p4, p5, p);
    t1 = sp_ntt_sub_simd0(p4, p5, p);
    t2 = sp_ntt_add_simd0(p6, p7, p);
    t3 = sp_ntt_sub_simd0(p6, p7, p);
    t4 = sp_ntt_add_simd0(t0, t2, p);
    t5 = sp_ntt_sub_simd0(t0, t2, p);
    t6 = sp_ntt_add_simd0(t1, t3, p);
    t7 = sp_ntt_sub_simd0(t1, t3, p);

    t0 = sp_ntt_add_simd0(p0, p2, p);
    t1 = sp_ntt_sub_simd0(p0, p2, p);
    t2 = sp_ntt_add_simd0(p1, p3, p);
    t3 = sp_ntt_sub_simd0(p1, p3, p);

    sp_simd_store(t2, out +  2 * ostride);
    sp_simd_store(t6, out +  7 * ostride);
    sp_simd_store(t1, out + 12 * ostride);
    sp_simd_store(t4, out + 17 * ostride);
    sp_simd_store(t3, out + 22 * ostride);
    sp_simd_store(t7, out + 27 * ostride);
    sp_simd_store(t0, out + 32 * ostride);
    sp_simd_store(t5, out + 37 * ostride);
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t p0 = a24;
    sp_simd_t p1 = a25;
    sp_simd_t p2 = a26;
    sp_simd_t p3 = a27;
    sp_simd_t p4 = a28;
    sp_simd_t p5 = a29;
    sp_simd_t p6 = a30;
    sp_simd_t p7 = a31;

    t0 = sp_ntt_add_simd0(p4, p5, p);
    t1 = sp_ntt_sub_simd0(p4, p5, p);
    t2 = sp_ntt_add_simd0(p6, p7, p);
    t3 = sp_ntt_sub_simd0(p6, p7, p);
    t4 = sp_ntt_add_simd0(t0, t2, p);
    t5 = sp_ntt_sub_simd0(t0, t2, p);
    t6 = sp_ntt_add_simd0(t1, t3, p);
    t7 = sp_ntt_sub_simd0(t1, t3, p);

    t0 = sp_ntt_add_simd0(p0, p2, p);
    t1 = sp_ntt_sub_simd0(p0, p2, p);
    t2 = sp_ntt_add_simd0(p1, p3, p);
    t3 = sp_ntt_sub_simd0(p1, p3, p);

    sp_simd_store(t7, out +  3 * ostride);
    sp_simd_store(t0, out +  8 * ostride);
    sp_simd_store(t5, out + 13 * ostride);
    sp_simd_store(t2, out + 18 * ostride);
    sp_simd_store(t6, out + 23 * ostride);
    sp_simd_store(t1, out + 28 * ostride);
    sp_simd_store(t4, out + 33 * ostride);
    sp_simd_store(t3, out + 38 * ostride);
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t p0 = a32;
    sp_simd_t p1 = a33;
    sp_simd_t p2 = a34;
    sp_simd_t p3 = a35;
    sp_simd_t p4 = a36;
    sp_simd_t p5 = a37;
    sp_simd_t p6 = a38;
    sp_simd_t p7 = a39;

    t0 = sp_ntt_add_simd0(p4, p5, p);
    t1 = sp_ntt_sub_simd0(p4, p5, p);
    t2 = sp_ntt_add_simd0(p6, p7, p);
    t3 = sp_ntt_sub_simd0(p6, p7, p);
    t4 = sp_ntt_add_simd0(t0, t2, p);
    t5 = sp_ntt_sub_simd0(t0, t2, p);
    t6 = sp_ntt_add_simd0(t1, t3, p);
    t7 = sp_ntt_sub_simd0(t1, t3, p);

    t0 = sp_ntt_add_simd0(p0, p2, p);
    t1 = sp_ntt_sub_simd0(p0, p2, p);
    t2 = sp_ntt_add_simd0(p1, p3, p);
    t3 = sp_ntt_sub_simd0(p1, p3, p);

    sp_simd_store(t1, out +  4 * ostride);
    sp_simd_store(t4, out +  9 * ostride);
    sp_simd_store(t3, out + 14 * ostride);
    sp_simd_store(t7, out + 19 * ostride);
    sp_simd_store(t0, out + 24 * ostride);
    sp_simd_store(t5, out + 29 * ostride);
    sp_simd_store(t2, out + 34 * ostride);
    sp_simd_store(t6, out + 39 * ostride);
  }
}

static void
ntt40_twiddle_run_core_simd(
        spv_t in, spv_size_t istride, spv_size_t idist,
		spv_t out, spv_size_t ostride, spv_size_t odist,
		sp_simd_t *w, sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  sp_simd_t a00, a01, a02, a03, a04, a05, a06, a07,
     	    a08, a09, a10, a11, a12, a13, a14, a15,
	    a16, a17, a18, a19, a20, a21, a22, a23,
	    a24, a25, a26, a27, a28, a29, a30, a31,
	    a32, a33, a34, a35, a36, a37, a38, a39;

  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x0 = sp_simd_gather(in +  0 * istride, idist, vsize);
    sp_simd_t x1 = sp_simd_gather(in +  5 * istride, idist, vsize);
    sp_simd_t x2 = sp_simd_gather(in + 10 * istride, idist, vsize);
    sp_simd_t x3 = sp_simd_gather(in + 15 * istride, idist, vsize);
    sp_simd_t x4 = sp_simd_gather(in + 20 * istride, idist, vsize);
    sp_simd_t x5 = sp_simd_gather(in + 25 * istride, idist, vsize);
    sp_simd_t x6 = sp_simd_gather(in + 30 * istride, idist, vsize);
    sp_simd_t x7 = sp_simd_gather(in + 35 * istride, idist, vsize);

    t0 = sp_ntt_add_simd(x0, x4, p);
    t4 = sp_ntt_sub_simd(x0, x4, p);
    t1 = sp_ntt_add_simd(x1, x5, p);
    t5 = sp_ntt_sub_simd(x1, x5, p);
    t2 = sp_ntt_add_simd(x2, x6, p);
    t6 = sp_ntt_sub_simd(x2, x6, p);
    t3 = sp_ntt_add_simd(x3, x7, p);
    t7 = sp_ntt_sub_simd(x3, x7, p);

    a00 = sp_ntt_add_simd(t0, t2, p);
    a01 = sp_ntt_sub_simd(t0, t2, p);
    a02 = sp_ntt_add_simd(t1, t3, p);
    a03 = sp_ntt_sub_simd(t1, t3, p);
    a04 = t4;
    a05 = t6;
    a06 = sp_ntt_sub_simd(t5, t7, p);
    a07 = sp_ntt_add_simd(t5, t7, p); 
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x7 = sp_simd_gather(in +  3 * istride, idist, vsize);
    sp_simd_t x0 = sp_simd_gather(in +  8 * istride, idist, vsize);
    sp_simd_t x1 = sp_simd_gather(in + 13 * istride, idist, vsize);
    sp_simd_t x2 = sp_simd_gather(in + 18 * istride, idist, vsize);
    sp_simd_t x3 = sp_simd_gather(in + 23 * istride, idist, vsize);
    sp_simd_t x4 = sp_simd_gather(in + 28 * istride, idist, vsize);
    sp_simd_t x5 = sp_simd_gather(in + 33 * istride, idist, vsize);
    sp_simd_t x6 = sp_simd_gather(in + 38 * istride, idist, vsize);

    t0 = sp_ntt_add_simd(x0, x4, p);
    t4 = sp_ntt_sub_simd(x0, x4, p);
    t1 = sp_ntt_add_simd(x1, x5, p);
    t5 = sp_ntt_sub_simd(x1, x5, p);
    t2 = sp_ntt_add_simd(x2, x6, p);
    t6 = sp_ntt_sub_simd(x2, x6, p);
    t3 = sp_ntt_add_simd(x3, x7, p);
    t7 = sp_ntt_sub_simd(x3, x7, p);

    a08 = sp_ntt_add_simd(t0, t2, p);
    a09 = sp_ntt_sub_simd(t0, t2, p);
    a10 = sp_ntt_add_simd(t1, t3, p);
    a11 = sp_ntt_sub_simd(t1, t3, p);
    a12 = t4;
    a13 = t6;
    a14 = sp_ntt_sub_simd(t5, t7, p);
    a15 = sp_ntt_add_simd(t5, t7, p); 
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x5 = sp_simd_gather(in +  1 * istride, idist, vsize);
    sp_simd_t x6 = sp_simd_gather(in +  6 * istride, idist, vsize);
    sp_simd_t x7 = sp_simd_gather(in + 11 * istride, idist, vsize);
    sp_simd_t x0 = sp_simd_gather(in + 16 * istride, idist, vsize);
    sp_simd_t x1 = sp_simd_gather(in + 21 * istride, idist, vsize);
    sp_simd_t x2 = sp_simd_gather(in + 26 * istride, idist, vsize);
    sp_simd_t x3 = sp_simd_gather(in + 31 * istride, idist, vsize);
    sp_simd_t x4 = sp_simd_gather(in + 36 * istride, idist, vsize);

    t0 = sp_ntt_add_simd(x0, x4, p);
    t4 = sp_ntt_sub_simd(x0, x4, p);
    t1 = sp_ntt_add_simd(x1, x5, p);
    t5 = sp_ntt_sub_simd(x1, x5, p);
    t2 = sp_ntt_add_simd(x2, x6, p);
    t6 = sp_ntt_sub_simd(x2, x6, p);
    t3 = sp_ntt_add_simd(x3, x7, p);
    t7 = sp_ntt_sub_simd(x3, x7, p);

    a16 = sp_ntt_add_simd(t0, t2, p);
    a17 = sp_ntt_sub_simd(t0, t2, p);
    a18 = sp_ntt_add_simd(t1, t3, p);
    a19 = sp_ntt_sub_simd(t1, t3, p);
    a20 = t4;
    a21 = t6;
    a22 = sp_ntt_sub_simd(t5, t7, p);
    a23 = sp_ntt_add_simd(t5, t7, p); 
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x0 = sp_simd_gather(in + 24 * istride, idist, vsize);
    sp_simd_t x1 = sp_simd_gather(in + 29 * istride, idist, vsize);
    sp_simd_t x2 = sp_simd_gather(in + 34 * istride, idist, vsize);
    sp_simd_t x3 = sp_simd_gather(in + 39 * istride, idist, vsize);
    sp_simd_t x4 = sp_simd_gather(in +  4 * istride, idist, vsize);
    sp_simd_t x5 = sp_simd_gather(in +  9 * istride, idist, vsize);
    sp_simd_t x6 = sp_simd_gather(in + 14 * istride, idist, vsize);
    sp_simd_t x7 = sp_simd_gather(in + 19 * istride, idist, vsize);

    t0 = sp_ntt_add_simd(x0, x4, p);
    t4 = sp_ntt_sub_simd(x0, x4, p);
    t1 = sp_ntt_add_simd(x1, x5, p);
    t5 = sp_ntt_sub_simd(x1, x5, p);
    t2 = sp_ntt_add_simd(x2, x6, p);
    t6 = sp_ntt_sub_simd(x2, x6, p);
    t3 = sp_ntt_add_simd(x3, x7, p);
    t7 = sp_ntt_sub_simd(x3, x7, p);

    a24 = sp_ntt_add_simd(t0, t2, p);
    a25 = sp_ntt_sub_simd(t0, t2, p);
    a26 = sp_ntt_add_simd(t1, t3, p);
    a27 = sp_ntt_sub_simd(t1, t3, p);
    a28 = t4;
    a29 = t6;
    a30 = sp_ntt_sub_simd(t5, t7, p);
    a31 = sp_ntt_add_simd(t5, t7, p); 
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x2 = sp_simd_gather(in +  2 * istride, idist, vsize);
    sp_simd_t x3 = sp_simd_gather(in +  7 * istride, idist, vsize);
    sp_simd_t x4 = sp_simd_gather(in + 12 * istride, idist, vsize);
    sp_simd_t x5 = sp_simd_gather(in + 17 * istride, idist, vsize);
    sp_simd_t x6 = sp_simd_gather(in + 22 * istride, idist, vsize);
    sp_simd_t x7 = sp_simd_gather(in + 27 * istride, idist, vsize);
    sp_simd_t x0 = sp_simd_gather(in + 32 * istride, idist, vsize);
    sp_simd_t x1 = sp_simd_gather(in + 37 * istride, idist, vsize);

    t0 = sp_ntt_add_simd(x0, x4, p);
    t4 = sp_ntt_sub_simd(x0, x4, p);
    t1 = sp_ntt_add_simd(x1, x5, p);
    t5 = sp_ntt_sub_simd(x1, x5, p);
    t2 = sp_ntt_add_simd(x2, x6, p);
    t6 = sp_ntt_sub_simd(x2, x6, p);
    t3 = sp_ntt_add_simd(x3, x7, p);
    t7 = sp_ntt_sub_simd(x3, x7, p);

    a32 = sp_ntt_add_simd(t0, t2, p);
    a33 = sp_ntt_sub_simd(t0, t2, p);
    a34 = sp_ntt_add_simd(t1, t3, p);
    a35 = sp_ntt_sub_simd(t1, t3, p);
    a36 = t4;
    a37 = t6;
    a38 = sp_ntt_sub_simd(t5, t7, p);
    a39 = sp_ntt_add_simd(t5, t7, p); 
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a00;
    b1 = a08;
    b4 = a16;
    b2 = a24;
    b3 = a32;

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
    a08 = b4;
    a16 = b3;
    a24 = b1;
    a32 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a01;
    b1 = a09;
    b4 = a17;
    b2 = a25;
    b3 = a33;

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
    a09 = b4;
    a17 = b3;
    a25 = b1;
    a33 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a02;
    b1 = a10;
    b4 = a18;
    b2 = a26;
    b3 = a34;

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
    a10 = b4;
    a18 = b3;
    a26 = b1;
    a34 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a03;
    b1 = a11;
    b4 = a19;
    b2 = a27;
    b3 = a35;

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

    b0 = sp_ntt_mul_simd(b0, ntt_const[18], ntt_const[NC+18], p);
    b1 = sp_ntt_mul_simd(b1, ntt_const[19], ntt_const[NC+19], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[20], ntt_const[NC+20], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[21], ntt_const[NC+21], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[22], ntt_const[NC+22], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[23], ntt_const[NC+23], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a03 = b0;
    a11 = b4;
    a19 = b3;
    a27 = b1;
    a35 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a04;
    b1 = a12;
    b4 = a20;
    b2 = a28;
    b3 = a36;

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

    b1 = sp_ntt_mul_simd(b1, ntt_const[25], ntt_const[NC+25], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[26], ntt_const[NC+26], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[27], ntt_const[NC+27], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[28], ntt_const[NC+28], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[29], ntt_const[NC+29], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a04 = b0;
    a12 = b4;
    a20 = b3;
    a28 = b1;
    a36 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a05;
    b1 = a13;
    b4 = a21;
    b2 = a29;
    b3 = a37;

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

    b0 = sp_ntt_mul_simd(b0, ntt_const[30], ntt_const[NC+30], p);
    b1 = sp_ntt_mul_simd(b1, ntt_const[31], ntt_const[NC+31], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[32], ntt_const[NC+32], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[33], ntt_const[NC+33], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[34], ntt_const[NC+34], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[35], ntt_const[NC+35], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a05 = b0;
    a13 = b4;
    a21 = b3;
    a29 = b1;
    a37 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a06;
    b1 = a14;
    b4 = a22;
    b2 = a30;
    b3 = a38;

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

    b0 = sp_ntt_mul_simd(b0, ntt_const[36], ntt_const[NC+36], p);
    b1 = sp_ntt_mul_simd(b1, ntt_const[37], ntt_const[NC+37], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[38], ntt_const[NC+38], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[39], ntt_const[NC+39], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[40], ntt_const[NC+40], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[41], ntt_const[NC+41], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a06 = b0;
    a14 = b4;
    a22 = b3;
    a30 = b1;
    a38 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a07;
    b1 = a15;
    b4 = a23;
    b2 = a31;
    b3 = a39;

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

    b0 = sp_ntt_mul_simd(b0, ntt_const[42], ntt_const[NC+42], p);
    b1 = sp_ntt_mul_simd(b1, ntt_const[43], ntt_const[NC+43], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[44], ntt_const[NC+44], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[45], ntt_const[NC+45], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[46], ntt_const[NC+46], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[47], ntt_const[NC+47], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a07 = b0;
    a15 = b4;
    a23 = b3;
    a31 = b1;
    a39 = b2;
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t p0 = a00;
    sp_simd_t p1 = a01;
    sp_simd_t p2 = a02;
    sp_simd_t p3 = a03;
    sp_simd_t p4 = a04;
    sp_simd_t p5 = a05;
    sp_simd_t p6 = a06;
    sp_simd_t p7 = a07;

    t0 = sp_ntt_add_simd(p4, p5, p);
    t1 = sp_ntt_sub_simd(p4, p5, p);
    t2 = sp_ntt_add_simd(p6, p7, p);
    t3 = sp_ntt_sub_simd(p6, p7, p);
    t4 = sp_ntt_add_partial_simd(t0, t2, p);
    t5 = sp_ntt_sub_partial_simd(t0, t2, p);
    t6 = sp_ntt_add_partial_simd(t1, t3, p);
    t7 = sp_ntt_sub_partial_simd(t1, t3, p);

    t0 = sp_ntt_add_simd(p0, p2, p);
    t1 = sp_ntt_sub_partial_simd(p0, p2, p);
    t2 = sp_ntt_add_partial_simd(p1, p3, p);
    t3 = sp_ntt_sub_partial_simd(p1, p3, p);

    t5 = sp_ntt_twiddle_mul_simd(t5, w + 8, p);
    t2 = sp_ntt_twiddle_mul_simd(t2, w + 18, p);
    t6 = sp_ntt_twiddle_mul_simd(t6, w + 28, p);
    t1 = sp_ntt_twiddle_mul_simd(t1, w + 38, p);
    t4 = sp_ntt_twiddle_mul_simd(t4, w + 48, p);
    t3 = sp_ntt_twiddle_mul_simd(t3, w + 58, p);
    t7 = sp_ntt_twiddle_mul_simd(t7, w + 68, p);

    sp_simd_scatter(t0, out +  0 * ostride, odist, vsize);
    sp_simd_scatter(t5, out +  5 * ostride, odist, vsize);
    sp_simd_scatter(t2, out + 10 * ostride, odist, vsize);
    sp_simd_scatter(t6, out + 15 * ostride, odist, vsize);
    sp_simd_scatter(t1, out + 20 * ostride, odist, vsize);
    sp_simd_scatter(t4, out + 25 * ostride, odist, vsize);
    sp_simd_scatter(t3, out + 30 * ostride, odist, vsize);
    sp_simd_scatter(t7, out + 35 * ostride, odist, vsize);
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t p0 = a08;
    sp_simd_t p1 = a09;
    sp_simd_t p2 = a10;
    sp_simd_t p3 = a11;
    sp_simd_t p4 = a12;
    sp_simd_t p5 = a13;
    sp_simd_t p6 = a14;
    sp_simd_t p7 = a15;

    t0 = sp_ntt_add_simd(p4, p5, p);
    t1 = sp_ntt_sub_simd(p4, p5, p);
    t2 = sp_ntt_add_simd(p6, p7, p);
    t3 = sp_ntt_sub_simd(p6, p7, p);
    t4 = sp_ntt_add_partial_simd(t0, t2, p);
    t5 = sp_ntt_sub_partial_simd(t0, t2, p);
    t6 = sp_ntt_add_partial_simd(t1, t3, p);
    t7 = sp_ntt_sub_partial_simd(t1, t3, p);

    t0 = sp_ntt_add_partial_simd(p0, p2, p);
    t1 = sp_ntt_sub_partial_simd(p0, p2, p);
    t2 = sp_ntt_add_partial_simd(p1, p3, p);
    t3 = sp_ntt_sub_partial_simd(p1, p3, p);

    t4 = sp_ntt_twiddle_mul_simd(t4, w + 0, p);
    t3 = sp_ntt_twiddle_mul_simd(t3, w + 10, p);
    t7 = sp_ntt_twiddle_mul_simd(t7, w + 20, p);
    t0 = sp_ntt_twiddle_mul_simd(t0, w + 30, p);
    t5 = sp_ntt_twiddle_mul_simd(t5, w + 40, p);
    t2 = sp_ntt_twiddle_mul_simd(t2, w + 50, p);
    t6 = sp_ntt_twiddle_mul_simd(t6, w + 60, p);
    t1 = sp_ntt_twiddle_mul_simd(t1, w + 70, p);

    sp_simd_scatter(t4, out +  1 * ostride, odist, vsize);
    sp_simd_scatter(t3, out +  6 * ostride, odist, vsize);
    sp_simd_scatter(t7, out + 11 * ostride, odist, vsize);
    sp_simd_scatter(t0, out + 16 * ostride, odist, vsize);
    sp_simd_scatter(t5, out + 21 * ostride, odist, vsize);
    sp_simd_scatter(t2, out + 26 * ostride, odist, vsize);
    sp_simd_scatter(t6, out + 31 * ostride, odist, vsize);
    sp_simd_scatter(t1, out + 36 * ostride, odist, vsize);
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t p0 = a16;
    sp_simd_t p1 = a17;
    sp_simd_t p2 = a18;
    sp_simd_t p3 = a19;
    sp_simd_t p4 = a20;
    sp_simd_t p5 = a21;
    sp_simd_t p6 = a22;
    sp_simd_t p7 = a23;

    t0 = sp_ntt_add_simd(p4, p5, p);
    t1 = sp_ntt_sub_simd(p4, p5, p);
    t2 = sp_ntt_add_simd(p6, p7, p);
    t3 = sp_ntt_sub_simd(p6, p7, p);
    t4 = sp_ntt_add_partial_simd(t0, t2, p);
    t5 = sp_ntt_sub_partial_simd(t0, t2, p);
    t6 = sp_ntt_add_partial_simd(t1, t3, p);
    t7 = sp_ntt_sub_partial_simd(t1, t3, p);

    t0 = sp_ntt_add_partial_simd(p0, p2, p);
    t1 = sp_ntt_sub_partial_simd(p0, p2, p);
    t2 = sp_ntt_add_partial_simd(p1, p3, p);
    t3 = sp_ntt_sub_partial_simd(p1, p3, p);

    t2 = sp_ntt_twiddle_mul_simd(t2, w + 2, p);
    t6 = sp_ntt_twiddle_mul_simd(t6, w + 12, p);
    t1 = sp_ntt_twiddle_mul_simd(t1, w + 22, p);
    t4 = sp_ntt_twiddle_mul_simd(t4, w + 32, p);
    t3 = sp_ntt_twiddle_mul_simd(t3, w + 42, p);
    t7 = sp_ntt_twiddle_mul_simd(t7, w + 52, p);
    t0 = sp_ntt_twiddle_mul_simd(t0, w + 62, p);
    t5 = sp_ntt_twiddle_mul_simd(t5, w + 72, p);

    sp_simd_scatter(t2, out +  2 * ostride, odist, vsize);
    sp_simd_scatter(t6, out +  7 * ostride, odist, vsize);
    sp_simd_scatter(t1, out + 12 * ostride, odist, vsize);
    sp_simd_scatter(t4, out + 17 * ostride, odist, vsize);
    sp_simd_scatter(t3, out + 22 * ostride, odist, vsize);
    sp_simd_scatter(t7, out + 27 * ostride, odist, vsize);
    sp_simd_scatter(t0, out + 32 * ostride, odist, vsize);
    sp_simd_scatter(t5, out + 37 * ostride, odist, vsize);
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t p0 = a24;
    sp_simd_t p1 = a25;
    sp_simd_t p2 = a26;
    sp_simd_t p3 = a27;
    sp_simd_t p4 = a28;
    sp_simd_t p5 = a29;
    sp_simd_t p6 = a30;
    sp_simd_t p7 = a31;

    t0 = sp_ntt_add_simd(p4, p5, p);
    t1 = sp_ntt_sub_simd(p4, p5, p);
    t2 = sp_ntt_add_simd(p6, p7, p);
    t3 = sp_ntt_sub_simd(p6, p7, p);
    t4 = sp_ntt_add_partial_simd(t0, t2, p);
    t5 = sp_ntt_sub_partial_simd(t0, t2, p);
    t6 = sp_ntt_add_partial_simd(t1, t3, p);
    t7 = sp_ntt_sub_partial_simd(t1, t3, p);

    t0 = sp_ntt_add_partial_simd(p0, p2, p);
    t1 = sp_ntt_sub_partial_simd(p0, p2, p);
    t2 = sp_ntt_add_partial_simd(p1, p3, p);
    t3 = sp_ntt_sub_partial_simd(p1, p3, p);

    t7 = sp_ntt_twiddle_mul_simd(t7, w + 4, p);
    t0 = sp_ntt_twiddle_mul_simd(t0, w + 14, p);
    t5 = sp_ntt_twiddle_mul_simd(t5, w + 24, p);
    t2 = sp_ntt_twiddle_mul_simd(t2, w + 34, p);
    t6 = sp_ntt_twiddle_mul_simd(t6, w + 44, p);
    t1 = sp_ntt_twiddle_mul_simd(t1, w + 54, p);
    t4 = sp_ntt_twiddle_mul_simd(t4, w + 64, p);
    t3 = sp_ntt_twiddle_mul_simd(t3, w + 74, p);

    sp_simd_scatter(t7, out +  3 * ostride, odist, vsize);
    sp_simd_scatter(t0, out +  8 * ostride, odist, vsize);
    sp_simd_scatter(t5, out + 13 * ostride, odist, vsize);
    sp_simd_scatter(t2, out + 18 * ostride, odist, vsize);
    sp_simd_scatter(t6, out + 23 * ostride, odist, vsize);
    sp_simd_scatter(t1, out + 28 * ostride, odist, vsize);
    sp_simd_scatter(t4, out + 33 * ostride, odist, vsize);
    sp_simd_scatter(t3, out + 38 * ostride, odist, vsize);
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t p0 = a32;
    sp_simd_t p1 = a33;
    sp_simd_t p2 = a34;
    sp_simd_t p3 = a35;
    sp_simd_t p4 = a36;
    sp_simd_t p5 = a37;
    sp_simd_t p6 = a38;
    sp_simd_t p7 = a39;

    t0 = sp_ntt_add_simd(p4, p5, p);
    t1 = sp_ntt_sub_simd(p4, p5, p);
    t2 = sp_ntt_add_simd(p6, p7, p);
    t3 = sp_ntt_sub_simd(p6, p7, p);
    t4 = sp_ntt_add_partial_simd(t0, t2, p);
    t5 = sp_ntt_sub_partial_simd(t0, t2, p);
    t6 = sp_ntt_add_partial_simd(t1, t3, p);
    t7 = sp_ntt_sub_partial_simd(t1, t3, p);

    t0 = sp_ntt_add_partial_simd(p0, p2, p);
    t1 = sp_ntt_sub_partial_simd(p0, p2, p);
    t2 = sp_ntt_add_partial_simd(p1, p3, p);
    t3 = sp_ntt_sub_partial_simd(p1, p3, p);

    t1 = sp_ntt_twiddle_mul_simd(t1, w + 6, p);
    t4 = sp_ntt_twiddle_mul_simd(t4, w + 16, p);
    t3 = sp_ntt_twiddle_mul_simd(t3, w + 26, p);
    t7 = sp_ntt_twiddle_mul_simd(t7, w + 36, p);
    t0 = sp_ntt_twiddle_mul_simd(t0, w + 46, p);
    t5 = sp_ntt_twiddle_mul_simd(t5, w + 56, p);
    t2 = sp_ntt_twiddle_mul_simd(t2, w + 66, p);
    t6 = sp_ntt_twiddle_mul_simd(t6, w + 76, p);

    sp_simd_scatter(t1, out +  4 * ostride, odist, vsize);
    sp_simd_scatter(t4, out +  9 * ostride, odist, vsize);
    sp_simd_scatter(t3, out + 14 * ostride, odist, vsize);
    sp_simd_scatter(t7, out + 19 * ostride, odist, vsize);
    sp_simd_scatter(t0, out + 24 * ostride, odist, vsize);
    sp_simd_scatter(t5, out + 29 * ostride, odist, vsize);
    sp_simd_scatter(t2, out + 34 * ostride, odist, vsize);
    sp_simd_scatter(t6, out + 39 * ostride, odist, vsize);
  }
}

static void
ntt40_twiddle_run_core_simd_interleaved(
		spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		sp_simd_t *w, sp_simd_t p, sp_simd_t * c)
{
  sp_simd_t a00, a01, a02, a03, a04, a05, a06, a07,
     	    a08, a09, a10, a11, a12, a13, a14, a15,
	    a16, a17, a18, a19, a20, a21, a22, a23,
	    a24, a25, a26, a27, a28, a29, a30, a31,
	    a32, a33, a34, a35, a36, a37, a38, a39;

  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x0 = sp_simd_load(in +  0 * istride);
    sp_simd_t x1 = sp_simd_load(in +  5 * istride);
    sp_simd_t x2 = sp_simd_load(in + 10 * istride);
    sp_simd_t x3 = sp_simd_load(in + 15 * istride);
    sp_simd_t x4 = sp_simd_load(in + 20 * istride);
    sp_simd_t x5 = sp_simd_load(in + 25 * istride);
    sp_simd_t x6 = sp_simd_load(in + 30 * istride);
    sp_simd_t x7 = sp_simd_load(in + 35 * istride);

    t0 = sp_ntt_add_simd0(x0, x4, p);
    t4 = sp_ntt_sub_simd0(x0, x4, p);
    t1 = sp_ntt_add_simd0(x1, x5, p);
    t5 = sp_ntt_sub_simd0(x1, x5, p);
    t2 = sp_ntt_add_simd0(x2, x6, p);
    t6 = sp_ntt_sub_simd0(x2, x6, p);
    t3 = sp_ntt_add_simd0(x3, x7, p);
    t7 = sp_ntt_sub_simd0(x3, x7, p);

    a00 = sp_ntt_add_simd0(t0, t2, p);
    a01 = sp_ntt_sub_simd0(t0, t2, p);
    a02 = sp_ntt_add_simd0(t1, t3, p);
    a03 = sp_ntt_sub_simd0(t1, t3, p);
    a04 = t4;
    a05 = t6;
    a06 = sp_ntt_sub_simd0(t5, t7, p);
    a07 = sp_ntt_add_simd0(t5, t7, p); 
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x7 = sp_simd_load(in +  3 * istride);
    sp_simd_t x0 = sp_simd_load(in +  8 * istride);
    sp_simd_t x1 = sp_simd_load(in + 13 * istride);
    sp_simd_t x2 = sp_simd_load(in + 18 * istride);
    sp_simd_t x3 = sp_simd_load(in + 23 * istride);
    sp_simd_t x4 = sp_simd_load(in + 28 * istride);
    sp_simd_t x5 = sp_simd_load(in + 33 * istride);
    sp_simd_t x6 = sp_simd_load(in + 38 * istride);

    t0 = sp_ntt_add_simd0(x0, x4, p);
    t4 = sp_ntt_sub_simd0(x0, x4, p);
    t1 = sp_ntt_add_simd0(x1, x5, p);
    t5 = sp_ntt_sub_simd0(x1, x5, p);
    t2 = sp_ntt_add_simd0(x2, x6, p);
    t6 = sp_ntt_sub_simd0(x2, x6, p);
    t3 = sp_ntt_add_simd0(x3, x7, p);
    t7 = sp_ntt_sub_simd0(x3, x7, p);

    a08 = sp_ntt_add_simd0(t0, t2, p);
    a09 = sp_ntt_sub_simd0(t0, t2, p);
    a10 = sp_ntt_add_simd0(t1, t3, p);
    a11 = sp_ntt_sub_simd0(t1, t3, p);
    a12 = t4;
    a13 = t6;
    a14 = sp_ntt_sub_simd0(t5, t7, p);
    a15 = sp_ntt_add_simd0(t5, t7, p); 
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x5 = sp_simd_load(in +  1 * istride);
    sp_simd_t x6 = sp_simd_load(in +  6 * istride);
    sp_simd_t x7 = sp_simd_load(in + 11 * istride);
    sp_simd_t x0 = sp_simd_load(in + 16 * istride);
    sp_simd_t x1 = sp_simd_load(in + 21 * istride);
    sp_simd_t x2 = sp_simd_load(in + 26 * istride);
    sp_simd_t x3 = sp_simd_load(in + 31 * istride);
    sp_simd_t x4 = sp_simd_load(in + 36 * istride);

    t0 = sp_ntt_add_simd0(x0, x4, p);
    t4 = sp_ntt_sub_simd0(x0, x4, p);
    t1 = sp_ntt_add_simd0(x1, x5, p);
    t5 = sp_ntt_sub_simd0(x1, x5, p);
    t2 = sp_ntt_add_simd0(x2, x6, p);
    t6 = sp_ntt_sub_simd0(x2, x6, p);
    t3 = sp_ntt_add_simd0(x3, x7, p);
    t7 = sp_ntt_sub_simd0(x3, x7, p);

    a16 = sp_ntt_add_simd0(t0, t2, p);
    a17 = sp_ntt_sub_simd0(t0, t2, p);
    a18 = sp_ntt_add_simd0(t1, t3, p);
    a19 = sp_ntt_sub_simd0(t1, t3, p);
    a20 = t4;
    a21 = t6;
    a22 = sp_ntt_sub_simd0(t5, t7, p);
    a23 = sp_ntt_add_simd0(t5, t7, p); 
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x0 = sp_simd_load(in + 24 * istride);
    sp_simd_t x1 = sp_simd_load(in + 29 * istride);
    sp_simd_t x2 = sp_simd_load(in + 34 * istride);
    sp_simd_t x3 = sp_simd_load(in + 39 * istride);
    sp_simd_t x4 = sp_simd_load(in +  4 * istride);
    sp_simd_t x5 = sp_simd_load(in +  9 * istride);
    sp_simd_t x6 = sp_simd_load(in + 14 * istride);
    sp_simd_t x7 = sp_simd_load(in + 19 * istride);

    t0 = sp_ntt_add_simd0(x0, x4, p);
    t4 = sp_ntt_sub_simd0(x0, x4, p);
    t1 = sp_ntt_add_simd0(x1, x5, p);
    t5 = sp_ntt_sub_simd0(x1, x5, p);
    t2 = sp_ntt_add_simd0(x2, x6, p);
    t6 = sp_ntt_sub_simd0(x2, x6, p);
    t3 = sp_ntt_add_simd0(x3, x7, p);
    t7 = sp_ntt_sub_simd0(x3, x7, p);

    a24 = sp_ntt_add_simd0(t0, t2, p);
    a25 = sp_ntt_sub_simd0(t0, t2, p);
    a26 = sp_ntt_add_simd0(t1, t3, p);
    a27 = sp_ntt_sub_simd0(t1, t3, p);
    a28 = t4;
    a29 = t6;
    a30 = sp_ntt_sub_simd0(t5, t7, p);
    a31 = sp_ntt_add_simd0(t5, t7, p); 
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x2 = sp_simd_load(in +  2 * istride);
    sp_simd_t x3 = sp_simd_load(in +  7 * istride);
    sp_simd_t x4 = sp_simd_load(in + 12 * istride);
    sp_simd_t x5 = sp_simd_load(in + 17 * istride);
    sp_simd_t x6 = sp_simd_load(in + 22 * istride);
    sp_simd_t x7 = sp_simd_load(in + 27 * istride);
    sp_simd_t x0 = sp_simd_load(in + 32 * istride);
    sp_simd_t x1 = sp_simd_load(in + 37 * istride);

    t0 = sp_ntt_add_simd0(x0, x4, p);
    t4 = sp_ntt_sub_simd0(x0, x4, p);
    t1 = sp_ntt_add_simd0(x1, x5, p);
    t5 = sp_ntt_sub_simd0(x1, x5, p);
    t2 = sp_ntt_add_simd0(x2, x6, p);
    t6 = sp_ntt_sub_simd0(x2, x6, p);
    t3 = sp_ntt_add_simd0(x3, x7, p);
    t7 = sp_ntt_sub_simd0(x3, x7, p);

    a32 = sp_ntt_add_simd0(t0, t2, p);
    a33 = sp_ntt_sub_simd0(t0, t2, p);
    a34 = sp_ntt_add_simd0(t1, t3, p);
    a35 = sp_ntt_sub_simd0(t1, t3, p);
    a36 = t4;
    a37 = t6;
    a38 = sp_ntt_sub_simd0(t5, t7, p);
    a39 = sp_ntt_add_simd0(t5, t7, p); 
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a00;
    b1 = a08;
    b4 = a16;
    b2 = a24;
    b3 = a32;

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
    a08 = b4;
    a16 = b3;
    a24 = b1;
    a32 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a01;
    b1 = a09;
    b4 = a17;
    b2 = a25;
    b3 = a33;

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
    a09 = b4;
    a17 = b3;
    a25 = b1;
    a33 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a02;
    b1 = a10;
    b4 = a18;
    b2 = a26;
    b3 = a34;

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
    a10 = b4;
    a18 = b3;
    a26 = b1;
    a34 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a03;
    b1 = a11;
    b4 = a19;
    b2 = a27;
    b3 = a35;

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

    b0 = sp_ntt_mul_simd0(b0, c+2*18, p);
    b1 = sp_ntt_mul_simd0(b1, c+2*19, p);
    b2 = sp_ntt_mul_simd0(b2, c+2*20, p);
    b3 = sp_ntt_mul_simd0(b3, c+2*21, p);
    b4 = sp_ntt_mul_simd0(b4, c+2*22, p);
    b5 = sp_ntt_mul_simd0(b5, c+2*23, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a03 = b0;
    a11 = b4;
    a19 = b3;
    a27 = b1;
    a35 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a04;
    b1 = a12;
    b4 = a20;
    b2 = a28;
    b3 = a36;

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

    b1 = sp_ntt_mul_simd0(b1, c+2*25, p);
    b2 = sp_ntt_mul_simd0(b2, c+2*26, p);
    b3 = sp_ntt_mul_simd0(b3, c+2*27, p);
    b4 = sp_ntt_mul_simd0(b4, c+2*28, p);
    b5 = sp_ntt_mul_simd0(b5, c+2*29, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a04 = b0;
    a12 = b4;
    a20 = b3;
    a28 = b1;
    a36 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a05;
    b1 = a13;
    b4 = a21;
    b2 = a29;
    b3 = a37;

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

    b0 = sp_ntt_mul_simd0(b0, c+2*30, p);
    b1 = sp_ntt_mul_simd0(b1, c+2*31, p);
    b2 = sp_ntt_mul_simd0(b2, c+2*32, p);
    b3 = sp_ntt_mul_simd0(b3, c+2*33, p);
    b4 = sp_ntt_mul_simd0(b4, c+2*34, p);
    b5 = sp_ntt_mul_simd0(b5, c+2*35, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a05 = b0;
    a13 = b4;
    a21 = b3;
    a29 = b1;
    a37 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a06;
    b1 = a14;
    b4 = a22;
    b2 = a30;
    b3 = a38;

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

    b0 = sp_ntt_mul_simd0(b0, c+2*36, p);
    b1 = sp_ntt_mul_simd0(b1, c+2*37, p);
    b2 = sp_ntt_mul_simd0(b2, c+2*38, p);
    b3 = sp_ntt_mul_simd0(b3, c+2*39, p);
    b4 = sp_ntt_mul_simd0(b4, c+2*40, p);
    b5 = sp_ntt_mul_simd0(b5, c+2*41, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a06 = b0;
    a14 = b4;
    a22 = b3;
    a30 = b1;
    a38 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a07;
    b1 = a15;
    b4 = a23;
    b2 = a31;
    b3 = a39;

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

    b0 = sp_ntt_mul_simd0(b0, c+2*42, p);
    b1 = sp_ntt_mul_simd0(b1, c+2*43, p);
    b2 = sp_ntt_mul_simd0(b2, c+2*44, p);
    b3 = sp_ntt_mul_simd0(b3, c+2*45, p);
    b4 = sp_ntt_mul_simd0(b4, c+2*46, p);
    b5 = sp_ntt_mul_simd0(b5, c+2*47, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a07 = b0;
    a15 = b4;
    a23 = b3;
    a31 = b1;
    a39 = b2;
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t p0 = a00;
    sp_simd_t p1 = a01;
    sp_simd_t p2 = a02;
    sp_simd_t p3 = a03;
    sp_simd_t p4 = a04;
    sp_simd_t p5 = a05;
    sp_simd_t p6 = a06;
    sp_simd_t p7 = a07;

    t0 = sp_ntt_add_simd0(p4, p5, p);
    t1 = sp_ntt_sub_simd0(p4, p5, p);
    t2 = sp_ntt_add_simd0(p6, p7, p);
    t3 = sp_ntt_sub_simd0(p6, p7, p);
    t4 = sp_ntt_add_partial_simd0(t0, t2, p);
    t5 = sp_ntt_sub_partial_simd0(t0, t2, p);
    t6 = sp_ntt_add_partial_simd0(t1, t3, p);
    t7 = sp_ntt_sub_partial_simd0(t1, t3, p);

    t0 = sp_ntt_add_simd0(p0, p2, p);
    t1 = sp_ntt_sub_partial_simd0(p0, p2, p);
    t2 = sp_ntt_add_partial_simd0(p1, p3, p);
    t3 = sp_ntt_sub_partial_simd0(p1, p3, p);

    t5 = sp_ntt_twiddle_mul_simd0(t5, w+8, p);
    t2 = sp_ntt_twiddle_mul_simd0(t2, w+18, p);
    t6 = sp_ntt_twiddle_mul_simd0(t6, w+28, p);
    t1 = sp_ntt_twiddle_mul_simd0(t1, w+38, p);
    t4 = sp_ntt_twiddle_mul_simd0(t4, w+48, p);
    t3 = sp_ntt_twiddle_mul_simd0(t3, w+58, p);
    t7 = sp_ntt_twiddle_mul_simd0(t7, w+68, p);

    sp_simd_store(t0, out +  0 * ostride);
    sp_simd_store(t5, out +  5 * ostride);
    sp_simd_store(t2, out + 10 * ostride);
    sp_simd_store(t6, out + 15 * ostride);
    sp_simd_store(t1, out + 20 * ostride);
    sp_simd_store(t4, out + 25 * ostride);
    sp_simd_store(t3, out + 30 * ostride);
    sp_simd_store(t7, out + 35 * ostride);
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t p0 = a08;
    sp_simd_t p1 = a09;
    sp_simd_t p2 = a10;
    sp_simd_t p3 = a11;
    sp_simd_t p4 = a12;
    sp_simd_t p5 = a13;
    sp_simd_t p6 = a14;
    sp_simd_t p7 = a15;

    t0 = sp_ntt_add_simd0(p4, p5, p);
    t1 = sp_ntt_sub_simd0(p4, p5, p);
    t2 = sp_ntt_add_simd0(p6, p7, p);
    t3 = sp_ntt_sub_simd0(p6, p7, p);
    t4 = sp_ntt_add_partial_simd0(t0, t2, p);
    t5 = sp_ntt_sub_partial_simd0(t0, t2, p);
    t6 = sp_ntt_add_partial_simd0(t1, t3, p);
    t7 = sp_ntt_sub_partial_simd0(t1, t3, p);

    t0 = sp_ntt_add_partial_simd0(p0, p2, p);
    t1 = sp_ntt_sub_partial_simd0(p0, p2, p);
    t2 = sp_ntt_add_partial_simd0(p1, p3, p);
    t3 = sp_ntt_sub_partial_simd0(p1, p3, p);

    t4 = sp_ntt_twiddle_mul_simd0(t4, w+0, p);
    t3 = sp_ntt_twiddle_mul_simd0(t3, w+10, p);
    t7 = sp_ntt_twiddle_mul_simd0(t7, w+20, p);
    t0 = sp_ntt_twiddle_mul_simd0(t0, w+30, p);
    t5 = sp_ntt_twiddle_mul_simd0(t5, w+40, p);
    t2 = sp_ntt_twiddle_mul_simd0(t2, w+50, p);
    t6 = sp_ntt_twiddle_mul_simd0(t6, w+60, p);
    t1 = sp_ntt_twiddle_mul_simd0(t1, w+70, p);

    sp_simd_store(t4, out +  1 * ostride);
    sp_simd_store(t3, out +  6 * ostride);
    sp_simd_store(t7, out + 11 * ostride);
    sp_simd_store(t0, out + 16 * ostride);
    sp_simd_store(t5, out + 21 * ostride);
    sp_simd_store(t2, out + 26 * ostride);
    sp_simd_store(t6, out + 31 * ostride);
    sp_simd_store(t1, out + 36 * ostride);
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t p0 = a16;
    sp_simd_t p1 = a17;
    sp_simd_t p2 = a18;
    sp_simd_t p3 = a19;
    sp_simd_t p4 = a20;
    sp_simd_t p5 = a21;
    sp_simd_t p6 = a22;
    sp_simd_t p7 = a23;

    t0 = sp_ntt_add_simd0(p4, p5, p);
    t1 = sp_ntt_sub_simd0(p4, p5, p);
    t2 = sp_ntt_add_simd0(p6, p7, p);
    t3 = sp_ntt_sub_simd0(p6, p7, p);
    t4 = sp_ntt_add_partial_simd0(t0, t2, p);
    t5 = sp_ntt_sub_partial_simd0(t0, t2, p);
    t6 = sp_ntt_add_partial_simd0(t1, t3, p);
    t7 = sp_ntt_sub_partial_simd0(t1, t3, p);

    t0 = sp_ntt_add_partial_simd0(p0, p2, p);
    t1 = sp_ntt_sub_partial_simd0(p0, p2, p);
    t2 = sp_ntt_add_partial_simd0(p1, p3, p);
    t3 = sp_ntt_sub_partial_simd0(p1, p3, p);

    t2 = sp_ntt_twiddle_mul_simd0(t2, w+2, p);
    t6 = sp_ntt_twiddle_mul_simd0(t6, w+12, p);
    t1 = sp_ntt_twiddle_mul_simd0(t1, w+22, p);
    t4 = sp_ntt_twiddle_mul_simd0(t4, w+32, p);
    t3 = sp_ntt_twiddle_mul_simd0(t3, w+42, p);
    t7 = sp_ntt_twiddle_mul_simd0(t7, w+52, p);
    t0 = sp_ntt_twiddle_mul_simd0(t0, w+62, p);
    t5 = sp_ntt_twiddle_mul_simd0(t5, w+72, p);

    sp_simd_store(t2, out +  2 * ostride);
    sp_simd_store(t6, out +  7 * ostride);
    sp_simd_store(t1, out + 12 * ostride);
    sp_simd_store(t4, out + 17 * ostride);
    sp_simd_store(t3, out + 22 * ostride);
    sp_simd_store(t7, out + 27 * ostride);
    sp_simd_store(t0, out + 32 * ostride);
    sp_simd_store(t5, out + 37 * ostride);
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t p0 = a24;
    sp_simd_t p1 = a25;
    sp_simd_t p2 = a26;
    sp_simd_t p3 = a27;
    sp_simd_t p4 = a28;
    sp_simd_t p5 = a29;
    sp_simd_t p6 = a30;
    sp_simd_t p7 = a31;

    t0 = sp_ntt_add_simd0(p4, p5, p);
    t1 = sp_ntt_sub_simd0(p4, p5, p);
    t2 = sp_ntt_add_simd0(p6, p7, p);
    t3 = sp_ntt_sub_simd0(p6, p7, p);
    t4 = sp_ntt_add_partial_simd0(t0, t2, p);
    t5 = sp_ntt_sub_partial_simd0(t0, t2, p);
    t6 = sp_ntt_add_partial_simd0(t1, t3, p);
    t7 = sp_ntt_sub_partial_simd0(t1, t3, p);

    t0 = sp_ntt_add_partial_simd0(p0, p2, p);
    t1 = sp_ntt_sub_partial_simd0(p0, p2, p);
    t2 = sp_ntt_add_partial_simd0(p1, p3, p);
    t3 = sp_ntt_sub_partial_simd0(p1, p3, p);

    t7 = sp_ntt_twiddle_mul_simd0(t7, w+4, p);
    t0 = sp_ntt_twiddle_mul_simd0(t0, w+14, p);
    t5 = sp_ntt_twiddle_mul_simd0(t5, w+24, p);
    t2 = sp_ntt_twiddle_mul_simd0(t2, w+34, p);
    t6 = sp_ntt_twiddle_mul_simd0(t6, w+44, p);
    t1 = sp_ntt_twiddle_mul_simd0(t1, w+54, p);
    t4 = sp_ntt_twiddle_mul_simd0(t4, w+64, p);
    t3 = sp_ntt_twiddle_mul_simd0(t3, w+74, p);

    sp_simd_store(t7, out +  3 * ostride);
    sp_simd_store(t0, out +  8 * ostride);
    sp_simd_store(t5, out + 13 * ostride);
    sp_simd_store(t2, out + 18 * ostride);
    sp_simd_store(t6, out + 23 * ostride);
    sp_simd_store(t1, out + 28 * ostride);
    sp_simd_store(t4, out + 33 * ostride);
    sp_simd_store(t3, out + 38 * ostride);
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t p0 = a32;
    sp_simd_t p1 = a33;
    sp_simd_t p2 = a34;
    sp_simd_t p3 = a35;
    sp_simd_t p4 = a36;
    sp_simd_t p5 = a37;
    sp_simd_t p6 = a38;
    sp_simd_t p7 = a39;

    t0 = sp_ntt_add_simd0(p4, p5, p);
    t1 = sp_ntt_sub_simd0(p4, p5, p);
    t2 = sp_ntt_add_simd0(p6, p7, p);
    t3 = sp_ntt_sub_simd0(p6, p7, p);
    t4 = sp_ntt_add_partial_simd0(t0, t2, p);
    t5 = sp_ntt_sub_partial_simd0(t0, t2, p);
    t6 = sp_ntt_add_partial_simd0(t1, t3, p);
    t7 = sp_ntt_sub_partial_simd0(t1, t3, p);

    t0 = sp_ntt_add_partial_simd0(p0, p2, p);
    t1 = sp_ntt_sub_partial_simd0(p0, p2, p);
    t2 = sp_ntt_add_partial_simd0(p1, p3, p);
    t3 = sp_ntt_sub_partial_simd0(p1, p3, p);

    t1 = sp_ntt_twiddle_mul_simd0(t1, w+6, p);
    t4 = sp_ntt_twiddle_mul_simd0(t4, w+16, p);
    t3 = sp_ntt_twiddle_mul_simd0(t3, w+26, p);
    t7 = sp_ntt_twiddle_mul_simd0(t7, w+36, p);
    t0 = sp_ntt_twiddle_mul_simd0(t0, w+46, p);
    t5 = sp_ntt_twiddle_mul_simd0(t5, w+56, p);
    t2 = sp_ntt_twiddle_mul_simd0(t2, w+66, p);
    t6 = sp_ntt_twiddle_mul_simd0(t6, w+76, p);

    sp_simd_store(t1, out +  4 * ostride);
    sp_simd_store(t4, out +  9 * ostride);
    sp_simd_store(t3, out + 14 * ostride);
    sp_simd_store(t7, out + 19 * ostride);
    sp_simd_store(t0, out + 24 * ostride);
    sp_simd_store(t5, out + 29 * ostride);
    sp_simd_store(t2, out + 34 * ostride);
    sp_simd_store(t6, out + 39 * ostride);
  }
}

static void
ntt40_pfa_run_core_simd(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t inc2, spv_size_t n,
	  sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  spv_size_t j00, j01, j02, j03, j04, j05, j06, j07,
	     j08, j09, j10, j11, j12, j13, j14, j15,
	     j16, j17, j18, j19, j20, j21, j22, j23,
	     j24, j25, j26, j27, j28, j29, j30, j31,
	     j32, j33, j34, j35, j36, j37, j38, j39;

  sp_simd_t a00, a01, a02, a03, a04, a05, a06, a07,
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
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x0 = sp_simd_pfa_gather(x, j00, inc2, n, vsize);
    sp_simd_t x1 = sp_simd_pfa_gather(x, j05, inc2, n, vsize);
    sp_simd_t x2 = sp_simd_pfa_gather(x, j10, inc2, n, vsize);
    sp_simd_t x3 = sp_simd_pfa_gather(x, j15, inc2, n, vsize);
    sp_simd_t x4 = sp_simd_pfa_gather(x, j20, inc2, n, vsize);
    sp_simd_t x5 = sp_simd_pfa_gather(x, j25, inc2, n, vsize);
    sp_simd_t x6 = sp_simd_pfa_gather(x, j30, inc2, n, vsize);
    sp_simd_t x7 = sp_simd_pfa_gather(x, j35, inc2, n, vsize);

    t0 = sp_ntt_add_simd(x0, x4, p);
    t4 = sp_ntt_sub_simd(x0, x4, p);
    t1 = sp_ntt_add_simd(x1, x5, p);
    t5 = sp_ntt_sub_simd(x1, x5, p);
    t2 = sp_ntt_add_simd(x2, x6, p);
    t6 = sp_ntt_sub_simd(x2, x6, p);
    t3 = sp_ntt_add_simd(x3, x7, p);
    t7 = sp_ntt_sub_simd(x3, x7, p);

    a00 = sp_ntt_add_simd(t0, t2, p);
    a01 = sp_ntt_sub_simd(t0, t2, p);
    a02 = sp_ntt_add_simd(t1, t3, p);
    a03 = sp_ntt_sub_simd(t1, t3, p);
    a04 = t4;
    a05 = t6;
    a06 = sp_ntt_sub_simd(t5, t7, p);
    a07 = sp_ntt_add_simd(t5, t7, p); 
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x7 = sp_simd_pfa_gather(x, j03, inc2, n, vsize);
    sp_simd_t x0 = sp_simd_pfa_gather(x, j08, inc2, n, vsize);
    sp_simd_t x1 = sp_simd_pfa_gather(x, j13, inc2, n, vsize);
    sp_simd_t x2 = sp_simd_pfa_gather(x, j18, inc2, n, vsize);
    sp_simd_t x3 = sp_simd_pfa_gather(x, j23, inc2, n, vsize);
    sp_simd_t x4 = sp_simd_pfa_gather(x, j28, inc2, n, vsize);
    sp_simd_t x5 = sp_simd_pfa_gather(x, j33, inc2, n, vsize);
    sp_simd_t x6 = sp_simd_pfa_gather(x, j38, inc2, n, vsize);

    t0 = sp_ntt_add_simd(x0, x4, p);
    t4 = sp_ntt_sub_simd(x0, x4, p);
    t1 = sp_ntt_add_simd(x1, x5, p);
    t5 = sp_ntt_sub_simd(x1, x5, p);
    t2 = sp_ntt_add_simd(x2, x6, p);
    t6 = sp_ntt_sub_simd(x2, x6, p);
    t3 = sp_ntt_add_simd(x3, x7, p);
    t7 = sp_ntt_sub_simd(x3, x7, p);

    a08 = sp_ntt_add_simd(t0, t2, p);
    a09 = sp_ntt_sub_simd(t0, t2, p);
    a10 = sp_ntt_add_simd(t1, t3, p);
    a11 = sp_ntt_sub_simd(t1, t3, p);
    a12 = t4;
    a13 = t6;
    a14 = sp_ntt_sub_simd(t5, t7, p);
    a15 = sp_ntt_add_simd(t5, t7, p); 
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x5 = sp_simd_pfa_gather(x, j01, inc2, n, vsize);
    sp_simd_t x6 = sp_simd_pfa_gather(x, j06, inc2, n, vsize);
    sp_simd_t x7 = sp_simd_pfa_gather(x, j11, inc2, n, vsize);
    sp_simd_t x0 = sp_simd_pfa_gather(x, j16, inc2, n, vsize);
    sp_simd_t x1 = sp_simd_pfa_gather(x, j21, inc2, n, vsize);
    sp_simd_t x2 = sp_simd_pfa_gather(x, j26, inc2, n, vsize);
    sp_simd_t x3 = sp_simd_pfa_gather(x, j31, inc2, n, vsize);
    sp_simd_t x4 = sp_simd_pfa_gather(x, j36, inc2, n, vsize);

    t0 = sp_ntt_add_simd(x0, x4, p);
    t4 = sp_ntt_sub_simd(x0, x4, p);
    t1 = sp_ntt_add_simd(x1, x5, p);
    t5 = sp_ntt_sub_simd(x1, x5, p);
    t2 = sp_ntt_add_simd(x2, x6, p);
    t6 = sp_ntt_sub_simd(x2, x6, p);
    t3 = sp_ntt_add_simd(x3, x7, p);
    t7 = sp_ntt_sub_simd(x3, x7, p);

    a16 = sp_ntt_add_simd(t0, t2, p);
    a17 = sp_ntt_sub_simd(t0, t2, p);
    a18 = sp_ntt_add_simd(t1, t3, p);
    a19 = sp_ntt_sub_simd(t1, t3, p);
    a20 = t4;
    a21 = t6;
    a22 = sp_ntt_sub_simd(t5, t7, p);
    a23 = sp_ntt_add_simd(t5, t7, p); 
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x0 = sp_simd_pfa_gather(x, j24, inc2, n, vsize);
    sp_simd_t x1 = sp_simd_pfa_gather(x, j29, inc2, n, vsize);
    sp_simd_t x2 = sp_simd_pfa_gather(x, j34, inc2, n, vsize);
    sp_simd_t x3 = sp_simd_pfa_gather(x, j39, inc2, n, vsize);
    sp_simd_t x4 = sp_simd_pfa_gather(x, j04, inc2, n, vsize);
    sp_simd_t x5 = sp_simd_pfa_gather(x, j09, inc2, n, vsize);
    sp_simd_t x6 = sp_simd_pfa_gather(x, j14, inc2, n, vsize);
    sp_simd_t x7 = sp_simd_pfa_gather(x, j19, inc2, n, vsize);

    t0 = sp_ntt_add_simd(x0, x4, p);
    t4 = sp_ntt_sub_simd(x0, x4, p);
    t1 = sp_ntt_add_simd(x1, x5, p);
    t5 = sp_ntt_sub_simd(x1, x5, p);
    t2 = sp_ntt_add_simd(x2, x6, p);
    t6 = sp_ntt_sub_simd(x2, x6, p);
    t3 = sp_ntt_add_simd(x3, x7, p);
    t7 = sp_ntt_sub_simd(x3, x7, p);

    a24 = sp_ntt_add_simd(t0, t2, p);
    a25 = sp_ntt_sub_simd(t0, t2, p);
    a26 = sp_ntt_add_simd(t1, t3, p);
    a27 = sp_ntt_sub_simd(t1, t3, p);
    a28 = t4;
    a29 = t6;
    a30 = sp_ntt_sub_simd(t5, t7, p);
    a31 = sp_ntt_add_simd(t5, t7, p); 
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x2 = sp_simd_pfa_gather(x, j02, inc2, n, vsize);
    sp_simd_t x3 = sp_simd_pfa_gather(x, j07, inc2, n, vsize);
    sp_simd_t x4 = sp_simd_pfa_gather(x, j12, inc2, n, vsize);
    sp_simd_t x5 = sp_simd_pfa_gather(x, j17, inc2, n, vsize);
    sp_simd_t x6 = sp_simd_pfa_gather(x, j22, inc2, n, vsize);
    sp_simd_t x7 = sp_simd_pfa_gather(x, j27, inc2, n, vsize);
    sp_simd_t x0 = sp_simd_pfa_gather(x, j32, inc2, n, vsize);
    sp_simd_t x1 = sp_simd_pfa_gather(x, j37, inc2, n, vsize);

    t0 = sp_ntt_add_simd(x0, x4, p);
    t4 = sp_ntt_sub_simd(x0, x4, p);
    t1 = sp_ntt_add_simd(x1, x5, p);
    t5 = sp_ntt_sub_simd(x1, x5, p);
    t2 = sp_ntt_add_simd(x2, x6, p);
    t6 = sp_ntt_sub_simd(x2, x6, p);
    t3 = sp_ntt_add_simd(x3, x7, p);
    t7 = sp_ntt_sub_simd(x3, x7, p);

    a32 = sp_ntt_add_simd(t0, t2, p);
    a33 = sp_ntt_sub_simd(t0, t2, p);
    a34 = sp_ntt_add_simd(t1, t3, p);
    a35 = sp_ntt_sub_simd(t1, t3, p);
    a36 = t4;
    a37 = t6;
    a38 = sp_ntt_sub_simd(t5, t7, p);
    a39 = sp_ntt_add_simd(t5, t7, p); 
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a00;
    b1 = a08;
    b4 = a16;
    b2 = a24;
    b3 = a32;

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
    a08 = b4;
    a16 = b3;
    a24 = b1;
    a32 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a01;
    b1 = a09;
    b4 = a17;
    b2 = a25;
    b3 = a33;

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
    a09 = b4;
    a17 = b3;
    a25 = b1;
    a33 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a02;
    b1 = a10;
    b4 = a18;
    b2 = a26;
    b3 = a34;

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
    a10 = b4;
    a18 = b3;
    a26 = b1;
    a34 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a03;
    b1 = a11;
    b4 = a19;
    b2 = a27;
    b3 = a35;

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

    b0 = sp_ntt_mul_simd(b0, ntt_const[18], ntt_const[NC+18], p);
    b1 = sp_ntt_mul_simd(b1, ntt_const[19], ntt_const[NC+19], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[20], ntt_const[NC+20], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[21], ntt_const[NC+21], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[22], ntt_const[NC+22], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[23], ntt_const[NC+23], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a03 = b0;
    a11 = b4;
    a19 = b3;
    a27 = b1;
    a35 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a04;
    b1 = a12;
    b4 = a20;
    b2 = a28;
    b3 = a36;

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

    b1 = sp_ntt_mul_simd(b1, ntt_const[25], ntt_const[NC+25], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[26], ntt_const[NC+26], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[27], ntt_const[NC+27], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[28], ntt_const[NC+28], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[29], ntt_const[NC+29], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a04 = b0;
    a12 = b4;
    a20 = b3;
    a28 = b1;
    a36 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a05;
    b1 = a13;
    b4 = a21;
    b2 = a29;
    b3 = a37;

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

    b0 = sp_ntt_mul_simd(b0, ntt_const[30], ntt_const[NC+30], p);
    b1 = sp_ntt_mul_simd(b1, ntt_const[31], ntt_const[NC+31], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[32], ntt_const[NC+32], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[33], ntt_const[NC+33], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[34], ntt_const[NC+34], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[35], ntt_const[NC+35], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a05 = b0;
    a13 = b4;
    a21 = b3;
    a29 = b1;
    a37 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a06;
    b1 = a14;
    b4 = a22;
    b2 = a30;
    b3 = a38;

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

    b0 = sp_ntt_mul_simd(b0, ntt_const[36], ntt_const[NC+36], p);
    b1 = sp_ntt_mul_simd(b1, ntt_const[37], ntt_const[NC+37], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[38], ntt_const[NC+38], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[39], ntt_const[NC+39], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[40], ntt_const[NC+40], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[41], ntt_const[NC+41], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a06 = b0;
    a14 = b4;
    a22 = b3;
    a30 = b1;
    a38 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a07;
    b1 = a15;
    b4 = a23;
    b2 = a31;
    b3 = a39;

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

    b0 = sp_ntt_mul_simd(b0, ntt_const[42], ntt_const[NC+42], p);
    b1 = sp_ntt_mul_simd(b1, ntt_const[43], ntt_const[NC+43], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[44], ntt_const[NC+44], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[45], ntt_const[NC+45], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[46], ntt_const[NC+46], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[47], ntt_const[NC+47], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a07 = b0;
    a15 = b4;
    a23 = b3;
    a31 = b1;
    a39 = b2;
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t p0 = a00;
    sp_simd_t p1 = a01;
    sp_simd_t p2 = a02;
    sp_simd_t p3 = a03;
    sp_simd_t p4 = a04;
    sp_simd_t p5 = a05;
    sp_simd_t p6 = a06;
    sp_simd_t p7 = a07;

    t0 = sp_ntt_add_simd(p4, p5, p);
    t1 = sp_ntt_sub_simd(p4, p5, p);
    t2 = sp_ntt_add_simd(p6, p7, p);
    t3 = sp_ntt_sub_simd(p6, p7, p);
    t4 = sp_ntt_add_simd(t0, t2, p);
    t5 = sp_ntt_sub_simd(t0, t2, p);
    t6 = sp_ntt_add_simd(t1, t3, p);
    t7 = sp_ntt_sub_simd(t1, t3, p);

    t0 = sp_ntt_add_simd(p0, p2, p);
    t1 = sp_ntt_sub_simd(p0, p2, p);
    t2 = sp_ntt_add_simd(p1, p3, p);
    t3 = sp_ntt_sub_simd(p1, p3, p);

    sp_simd_pfa_scatter(t0, x, j00, inc2, n, vsize);
    sp_simd_pfa_scatter(t5, x, j05, inc2, n, vsize);
    sp_simd_pfa_scatter(t2, x, j10, inc2, n, vsize);
    sp_simd_pfa_scatter(t6, x, j15, inc2, n, vsize);
    sp_simd_pfa_scatter(t1, x, j20, inc2, n, vsize);
    sp_simd_pfa_scatter(t4, x, j25, inc2, n, vsize);
    sp_simd_pfa_scatter(t3, x, j30, inc2, n, vsize);
    sp_simd_pfa_scatter(t7, x, j35, inc2, n, vsize);
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t p0 = a08;
    sp_simd_t p1 = a09;
    sp_simd_t p2 = a10;
    sp_simd_t p3 = a11;
    sp_simd_t p4 = a12;
    sp_simd_t p5 = a13;
    sp_simd_t p6 = a14;
    sp_simd_t p7 = a15;

    t0 = sp_ntt_add_simd(p4, p5, p);
    t1 = sp_ntt_sub_simd(p4, p5, p);
    t2 = sp_ntt_add_simd(p6, p7, p);
    t3 = sp_ntt_sub_simd(p6, p7, p);
    t4 = sp_ntt_add_simd(t0, t2, p);
    t5 = sp_ntt_sub_simd(t0, t2, p);
    t6 = sp_ntt_add_simd(t1, t3, p);
    t7 = sp_ntt_sub_simd(t1, t3, p);

    t0 = sp_ntt_add_simd(p0, p2, p);
    t1 = sp_ntt_sub_simd(p0, p2, p);
    t2 = sp_ntt_add_simd(p1, p3, p);
    t3 = sp_ntt_sub_simd(p1, p3, p);

    sp_simd_pfa_scatter(t4, x, j01, inc2, n, vsize);
    sp_simd_pfa_scatter(t3, x, j06, inc2, n, vsize);
    sp_simd_pfa_scatter(t7, x, j11, inc2, n, vsize);
    sp_simd_pfa_scatter(t0, x, j16, inc2, n, vsize);
    sp_simd_pfa_scatter(t5, x, j21, inc2, n, vsize);
    sp_simd_pfa_scatter(t2, x, j26, inc2, n, vsize);
    sp_simd_pfa_scatter(t6, x, j31, inc2, n, vsize);
    sp_simd_pfa_scatter(t1, x, j36, inc2, n, vsize);
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t p0 = a16;
    sp_simd_t p1 = a17;
    sp_simd_t p2 = a18;
    sp_simd_t p3 = a19;
    sp_simd_t p4 = a20;
    sp_simd_t p5 = a21;
    sp_simd_t p6 = a22;
    sp_simd_t p7 = a23;

    t0 = sp_ntt_add_simd(p4, p5, p);
    t1 = sp_ntt_sub_simd(p4, p5, p);
    t2 = sp_ntt_add_simd(p6, p7, p);
    t3 = sp_ntt_sub_simd(p6, p7, p);
    t4 = sp_ntt_add_simd(t0, t2, p);
    t5 = sp_ntt_sub_simd(t0, t2, p);
    t6 = sp_ntt_add_simd(t1, t3, p);
    t7 = sp_ntt_sub_simd(t1, t3, p);

    t0 = sp_ntt_add_simd(p0, p2, p);
    t1 = sp_ntt_sub_simd(p0, p2, p);
    t2 = sp_ntt_add_simd(p1, p3, p);
    t3 = sp_ntt_sub_simd(p1, p3, p);

    sp_simd_pfa_scatter(t2, x, j02, inc2, n, vsize);
    sp_simd_pfa_scatter(t6, x, j07, inc2, n, vsize);
    sp_simd_pfa_scatter(t1, x, j12, inc2, n, vsize);
    sp_simd_pfa_scatter(t4, x, j17, inc2, n, vsize);
    sp_simd_pfa_scatter(t3, x, j22, inc2, n, vsize);
    sp_simd_pfa_scatter(t7, x, j27, inc2, n, vsize);
    sp_simd_pfa_scatter(t0, x, j32, inc2, n, vsize);
    sp_simd_pfa_scatter(t5, x, j37, inc2, n, vsize);
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t p0 = a24;
    sp_simd_t p1 = a25;
    sp_simd_t p2 = a26;
    sp_simd_t p3 = a27;
    sp_simd_t p4 = a28;
    sp_simd_t p5 = a29;
    sp_simd_t p6 = a30;
    sp_simd_t p7 = a31;

    t0 = sp_ntt_add_simd(p4, p5, p);
    t1 = sp_ntt_sub_simd(p4, p5, p);
    t2 = sp_ntt_add_simd(p6, p7, p);
    t3 = sp_ntt_sub_simd(p6, p7, p);
    t4 = sp_ntt_add_simd(t0, t2, p);
    t5 = sp_ntt_sub_simd(t0, t2, p);
    t6 = sp_ntt_add_simd(t1, t3, p);
    t7 = sp_ntt_sub_simd(t1, t3, p);

    t0 = sp_ntt_add_simd(p0, p2, p);
    t1 = sp_ntt_sub_simd(p0, p2, p);
    t2 = sp_ntt_add_simd(p1, p3, p);
    t3 = sp_ntt_sub_simd(p1, p3, p);

    sp_simd_pfa_scatter(t7, x, j03, inc2, n, vsize);
    sp_simd_pfa_scatter(t0, x, j08, inc2, n, vsize);
    sp_simd_pfa_scatter(t5, x, j13, inc2, n, vsize);
    sp_simd_pfa_scatter(t2, x, j18, inc2, n, vsize);
    sp_simd_pfa_scatter(t6, x, j23, inc2, n, vsize);
    sp_simd_pfa_scatter(t1, x, j28, inc2, n, vsize);
    sp_simd_pfa_scatter(t4, x, j33, inc2, n, vsize);
    sp_simd_pfa_scatter(t3, x, j38, inc2, n, vsize);
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t p0 = a32;
    sp_simd_t p1 = a33;
    sp_simd_t p2 = a34;
    sp_simd_t p3 = a35;
    sp_simd_t p4 = a36;
    sp_simd_t p5 = a37;
    sp_simd_t p6 = a38;
    sp_simd_t p7 = a39;

    t0 = sp_ntt_add_simd(p4, p5, p);
    t1 = sp_ntt_sub_simd(p4, p5, p);
    t2 = sp_ntt_add_simd(p6, p7, p);
    t3 = sp_ntt_sub_simd(p6, p7, p);
    t4 = sp_ntt_add_simd(t0, t2, p);
    t5 = sp_ntt_sub_simd(t0, t2, p);
    t6 = sp_ntt_add_simd(t1, t3, p);
    t7 = sp_ntt_sub_simd(t1, t3, p);

    t0 = sp_ntt_add_simd(p0, p2, p);
    t1 = sp_ntt_sub_simd(p0, p2, p);
    t2 = sp_ntt_add_simd(p1, p3, p);
    t3 = sp_ntt_sub_simd(p1, p3, p);

    sp_simd_pfa_scatter(t1, x, j04, inc2, n, vsize);
    sp_simd_pfa_scatter(t4, x, j09, inc2, n, vsize);
    sp_simd_pfa_scatter(t3, x, j14, inc2, n, vsize);
    sp_simd_pfa_scatter(t7, x, j19, inc2, n, vsize);
    sp_simd_pfa_scatter(t0, x, j24, inc2, n, vsize);
    sp_simd_pfa_scatter(t5, x, j29, inc2, n, vsize);
    sp_simd_pfa_scatter(t2, x, j34, inc2, n, vsize);
    sp_simd_pfa_scatter(t6, x, j39, inc2, n, vsize);
  }
}

static void
ntt40_pfa_run_core_simd_interleaved(
	  spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t n,
	  sp_simd_t p, sp_simd_t * c)
{
  spv_size_t j00, j01, j02, j03, j04, j05, j06, j07,
	     j08, j09, j10, j11, j12, j13, j14, j15,
	     j16, j17, j18, j19, j20, j21, j22, j23,
	     j24, j25, j26, j27, j28, j29, j30, j31,
	     j32, j33, j34, j35, j36, j37, j38, j39;

  sp_simd_t a00, a01, a02, a03, a04, a05, a06, a07,
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
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x0 = sp_simd_load(x + j00);
    sp_simd_t x1 = sp_simd_load(x + j05);
    sp_simd_t x2 = sp_simd_load(x + j10);
    sp_simd_t x3 = sp_simd_load(x + j15);
    sp_simd_t x4 = sp_simd_load(x + j20);
    sp_simd_t x5 = sp_simd_load(x + j25);
    sp_simd_t x6 = sp_simd_load(x + j30);
    sp_simd_t x7 = sp_simd_load(x + j35);

    t0 = sp_ntt_add_simd0(x0, x4, p);
    t4 = sp_ntt_sub_simd0(x0, x4, p);
    t1 = sp_ntt_add_simd0(x1, x5, p);
    t5 = sp_ntt_sub_simd0(x1, x5, p);
    t2 = sp_ntt_add_simd0(x2, x6, p);
    t6 = sp_ntt_sub_simd0(x2, x6, p);
    t3 = sp_ntt_add_simd0(x3, x7, p);
    t7 = sp_ntt_sub_simd0(x3, x7, p);

    a00 = sp_ntt_add_simd0(t0, t2, p);
    a01 = sp_ntt_sub_simd0(t0, t2, p);
    a02 = sp_ntt_add_simd0(t1, t3, p);
    a03 = sp_ntt_sub_simd0(t1, t3, p);
    a04 = t4;
    a05 = t6;
    a06 = sp_ntt_sub_simd0(t5, t7, p);
    a07 = sp_ntt_add_simd0(t5, t7, p); 
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x7 = sp_simd_load(x + j03);
    sp_simd_t x0 = sp_simd_load(x + j08);
    sp_simd_t x1 = sp_simd_load(x + j13);
    sp_simd_t x2 = sp_simd_load(x + j18);
    sp_simd_t x3 = sp_simd_load(x + j23);
    sp_simd_t x4 = sp_simd_load(x + j28);
    sp_simd_t x5 = sp_simd_load(x + j33);
    sp_simd_t x6 = sp_simd_load(x + j38);

    t0 = sp_ntt_add_simd0(x0, x4, p);
    t4 = sp_ntt_sub_simd0(x0, x4, p);
    t1 = sp_ntt_add_simd0(x1, x5, p);
    t5 = sp_ntt_sub_simd0(x1, x5, p);
    t2 = sp_ntt_add_simd0(x2, x6, p);
    t6 = sp_ntt_sub_simd0(x2, x6, p);
    t3 = sp_ntt_add_simd0(x3, x7, p);
    t7 = sp_ntt_sub_simd0(x3, x7, p);

    a08 = sp_ntt_add_simd0(t0, t2, p);
    a09 = sp_ntt_sub_simd0(t0, t2, p);
    a10 = sp_ntt_add_simd0(t1, t3, p);
    a11 = sp_ntt_sub_simd0(t1, t3, p);
    a12 = t4;
    a13 = t6;
    a14 = sp_ntt_sub_simd0(t5, t7, p);
    a15 = sp_ntt_add_simd0(t5, t7, p); 
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x5 = sp_simd_load(x + j01);
    sp_simd_t x6 = sp_simd_load(x + j06);
    sp_simd_t x7 = sp_simd_load(x + j11);
    sp_simd_t x0 = sp_simd_load(x + j16);
    sp_simd_t x1 = sp_simd_load(x + j21);
    sp_simd_t x2 = sp_simd_load(x + j26);
    sp_simd_t x3 = sp_simd_load(x + j31);
    sp_simd_t x4 = sp_simd_load(x + j36);

    t0 = sp_ntt_add_simd0(x0, x4, p);
    t4 = sp_ntt_sub_simd0(x0, x4, p);
    t1 = sp_ntt_add_simd0(x1, x5, p);
    t5 = sp_ntt_sub_simd0(x1, x5, p);
    t2 = sp_ntt_add_simd0(x2, x6, p);
    t6 = sp_ntt_sub_simd0(x2, x6, p);
    t3 = sp_ntt_add_simd0(x3, x7, p);
    t7 = sp_ntt_sub_simd0(x3, x7, p);

    a16 = sp_ntt_add_simd0(t0, t2, p);
    a17 = sp_ntt_sub_simd0(t0, t2, p);
    a18 = sp_ntt_add_simd0(t1, t3, p);
    a19 = sp_ntt_sub_simd0(t1, t3, p);
    a20 = t4;
    a21 = t6;
    a22 = sp_ntt_sub_simd0(t5, t7, p);
    a23 = sp_ntt_add_simd0(t5, t7, p); 
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x0 = sp_simd_load(x + j24);
    sp_simd_t x1 = sp_simd_load(x + j29);
    sp_simd_t x2 = sp_simd_load(x + j34);
    sp_simd_t x3 = sp_simd_load(x + j39);
    sp_simd_t x4 = sp_simd_load(x + j04);
    sp_simd_t x5 = sp_simd_load(x + j09);
    sp_simd_t x6 = sp_simd_load(x + j14);
    sp_simd_t x7 = sp_simd_load(x + j19);

    t0 = sp_ntt_add_simd0(x0, x4, p);
    t4 = sp_ntt_sub_simd0(x0, x4, p);
    t1 = sp_ntt_add_simd0(x1, x5, p);
    t5 = sp_ntt_sub_simd0(x1, x5, p);
    t2 = sp_ntt_add_simd0(x2, x6, p);
    t6 = sp_ntt_sub_simd0(x2, x6, p);
    t3 = sp_ntt_add_simd0(x3, x7, p);
    t7 = sp_ntt_sub_simd0(x3, x7, p);

    a24 = sp_ntt_add_simd0(t0, t2, p);
    a25 = sp_ntt_sub_simd0(t0, t2, p);
    a26 = sp_ntt_add_simd0(t1, t3, p);
    a27 = sp_ntt_sub_simd0(t1, t3, p);
    a28 = t4;
    a29 = t6;
    a30 = sp_ntt_sub_simd0(t5, t7, p);
    a31 = sp_ntt_add_simd0(t5, t7, p); 
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x2 = sp_simd_load(x + j02);
    sp_simd_t x3 = sp_simd_load(x + j07);
    sp_simd_t x4 = sp_simd_load(x + j12);
    sp_simd_t x5 = sp_simd_load(x + j17);
    sp_simd_t x6 = sp_simd_load(x + j22);
    sp_simd_t x7 = sp_simd_load(x + j27);
    sp_simd_t x0 = sp_simd_load(x + j32);
    sp_simd_t x1 = sp_simd_load(x + j37);

    t0 = sp_ntt_add_simd0(x0, x4, p);
    t4 = sp_ntt_sub_simd0(x0, x4, p);
    t1 = sp_ntt_add_simd0(x1, x5, p);
    t5 = sp_ntt_sub_simd0(x1, x5, p);
    t2 = sp_ntt_add_simd0(x2, x6, p);
    t6 = sp_ntt_sub_simd0(x2, x6, p);
    t3 = sp_ntt_add_simd0(x3, x7, p);
    t7 = sp_ntt_sub_simd0(x3, x7, p);

    a32 = sp_ntt_add_simd0(t0, t2, p);
    a33 = sp_ntt_sub_simd0(t0, t2, p);
    a34 = sp_ntt_add_simd0(t1, t3, p);
    a35 = sp_ntt_sub_simd0(t1, t3, p);
    a36 = t4;
    a37 = t6;
    a38 = sp_ntt_sub_simd0(t5, t7, p);
    a39 = sp_ntt_add_simd0(t5, t7, p); 
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a00;
    b1 = a08;
    b4 = a16;
    b2 = a24;
    b3 = a32;

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
    a08 = b4;
    a16 = b3;
    a24 = b1;
    a32 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a01;
    b1 = a09;
    b4 = a17;
    b2 = a25;
    b3 = a33;

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
    a09 = b4;
    a17 = b3;
    a25 = b1;
    a33 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a02;
    b1 = a10;
    b4 = a18;
    b2 = a26;
    b3 = a34;

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
    a10 = b4;
    a18 = b3;
    a26 = b1;
    a34 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a03;
    b1 = a11;
    b4 = a19;
    b2 = a27;
    b3 = a35;

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

    b0 = sp_ntt_mul_simd0(b0, c+2*18, p);
    b1 = sp_ntt_mul_simd0(b1, c+2*19, p);
    b2 = sp_ntt_mul_simd0(b2, c+2*20, p);
    b3 = sp_ntt_mul_simd0(b3, c+2*21, p);
    b4 = sp_ntt_mul_simd0(b4, c+2*22, p);
    b5 = sp_ntt_mul_simd0(b5, c+2*23, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a03 = b0;
    a11 = b4;
    a19 = b3;
    a27 = b1;
    a35 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a04;
    b1 = a12;
    b4 = a20;
    b2 = a28;
    b3 = a36;

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

    b1 = sp_ntt_mul_simd0(b1, c+2*25, p);
    b2 = sp_ntt_mul_simd0(b2, c+2*26, p);
    b3 = sp_ntt_mul_simd0(b3, c+2*27, p);
    b4 = sp_ntt_mul_simd0(b4, c+2*28, p);
    b5 = sp_ntt_mul_simd0(b5, c+2*29, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a04 = b0;
    a12 = b4;
    a20 = b3;
    a28 = b1;
    a36 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a05;
    b1 = a13;
    b4 = a21;
    b2 = a29;
    b3 = a37;

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

    b0 = sp_ntt_mul_simd0(b0, c+2*30, p);
    b1 = sp_ntt_mul_simd0(b1, c+2*31, p);
    b2 = sp_ntt_mul_simd0(b2, c+2*32, p);
    b3 = sp_ntt_mul_simd0(b3, c+2*33, p);
    b4 = sp_ntt_mul_simd0(b4, c+2*34, p);
    b5 = sp_ntt_mul_simd0(b5, c+2*35, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a05 = b0;
    a13 = b4;
    a21 = b3;
    a29 = b1;
    a37 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a06;
    b1 = a14;
    b4 = a22;
    b2 = a30;
    b3 = a38;

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

    b0 = sp_ntt_mul_simd0(b0, c+2*36, p);
    b1 = sp_ntt_mul_simd0(b1, c+2*37, p);
    b2 = sp_ntt_mul_simd0(b2, c+2*38, p);
    b3 = sp_ntt_mul_simd0(b3, c+2*39, p);
    b4 = sp_ntt_mul_simd0(b4, c+2*40, p);
    b5 = sp_ntt_mul_simd0(b5, c+2*41, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a06 = b0;
    a14 = b4;
    a22 = b3;
    a30 = b1;
    a38 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a07;
    b1 = a15;
    b4 = a23;
    b2 = a31;
    b3 = a39;

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

    b0 = sp_ntt_mul_simd0(b0, c+2*42, p);
    b1 = sp_ntt_mul_simd0(b1, c+2*43, p);
    b2 = sp_ntt_mul_simd0(b2, c+2*44, p);
    b3 = sp_ntt_mul_simd0(b3, c+2*45, p);
    b4 = sp_ntt_mul_simd0(b4, c+2*46, p);
    b5 = sp_ntt_mul_simd0(b5, c+2*47, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a07 = b0;
    a15 = b4;
    a23 = b3;
    a31 = b1;
    a39 = b2;
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t p0 = a00;
    sp_simd_t p1 = a01;
    sp_simd_t p2 = a02;
    sp_simd_t p3 = a03;
    sp_simd_t p4 = a04;
    sp_simd_t p5 = a05;
    sp_simd_t p6 = a06;
    sp_simd_t p7 = a07;

    t0 = sp_ntt_add_simd0(p4, p5, p);
    t1 = sp_ntt_sub_simd0(p4, p5, p);
    t2 = sp_ntt_add_simd0(p6, p7, p);
    t3 = sp_ntt_sub_simd0(p6, p7, p);
    t4 = sp_ntt_add_simd0(t0, t2, p);
    t5 = sp_ntt_sub_simd0(t0, t2, p);
    t6 = sp_ntt_add_simd0(t1, t3, p);
    t7 = sp_ntt_sub_simd0(t1, t3, p);

    t0 = sp_ntt_add_simd0(p0, p2, p);
    t1 = sp_ntt_sub_simd0(p0, p2, p);
    t2 = sp_ntt_add_simd0(p1, p3, p);
    t3 = sp_ntt_sub_simd0(p1, p3, p);

    sp_simd_store(t0, x + j00);
    sp_simd_store(t5, x + j05);
    sp_simd_store(t2, x + j10);
    sp_simd_store(t6, x + j15);
    sp_simd_store(t1, x + j20);
    sp_simd_store(t4, x + j25);
    sp_simd_store(t3, x + j30);
    sp_simd_store(t7, x + j35);
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t p0 = a08;
    sp_simd_t p1 = a09;
    sp_simd_t p2 = a10;
    sp_simd_t p3 = a11;
    sp_simd_t p4 = a12;
    sp_simd_t p5 = a13;
    sp_simd_t p6 = a14;
    sp_simd_t p7 = a15;

    t0 = sp_ntt_add_simd0(p4, p5, p);
    t1 = sp_ntt_sub_simd0(p4, p5, p);
    t2 = sp_ntt_add_simd0(p6, p7, p);
    t3 = sp_ntt_sub_simd0(p6, p7, p);
    t4 = sp_ntt_add_simd0(t0, t2, p);
    t5 = sp_ntt_sub_simd0(t0, t2, p);
    t6 = sp_ntt_add_simd0(t1, t3, p);
    t7 = sp_ntt_sub_simd0(t1, t3, p);

    t0 = sp_ntt_add_simd0(p0, p2, p);
    t1 = sp_ntt_sub_simd0(p0, p2, p);
    t2 = sp_ntt_add_simd0(p1, p3, p);
    t3 = sp_ntt_sub_simd0(p1, p3, p);

    sp_simd_store(t4, x + j01);
    sp_simd_store(t3, x + j06);
    sp_simd_store(t7, x + j11);
    sp_simd_store(t0, x + j16);
    sp_simd_store(t5, x + j21);
    sp_simd_store(t2, x + j26);
    sp_simd_store(t6, x + j31);
    sp_simd_store(t1, x + j36);
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t p0 = a16;
    sp_simd_t p1 = a17;
    sp_simd_t p2 = a18;
    sp_simd_t p3 = a19;
    sp_simd_t p4 = a20;
    sp_simd_t p5 = a21;
    sp_simd_t p6 = a22;
    sp_simd_t p7 = a23;

    t0 = sp_ntt_add_simd0(p4, p5, p);
    t1 = sp_ntt_sub_simd0(p4, p5, p);
    t2 = sp_ntt_add_simd0(p6, p7, p);
    t3 = sp_ntt_sub_simd0(p6, p7, p);
    t4 = sp_ntt_add_simd0(t0, t2, p);
    t5 = sp_ntt_sub_simd0(t0, t2, p);
    t6 = sp_ntt_add_simd0(t1, t3, p);
    t7 = sp_ntt_sub_simd0(t1, t3, p);

    t0 = sp_ntt_add_simd0(p0, p2, p);
    t1 = sp_ntt_sub_simd0(p0, p2, p);
    t2 = sp_ntt_add_simd0(p1, p3, p);
    t3 = sp_ntt_sub_simd0(p1, p3, p);

    sp_simd_store(t2, x + j02);
    sp_simd_store(t6, x + j07);
    sp_simd_store(t1, x + j12);
    sp_simd_store(t4, x + j17);
    sp_simd_store(t3, x + j22);
    sp_simd_store(t7, x + j27);
    sp_simd_store(t0, x + j32);
    sp_simd_store(t5, x + j37);
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t p0 = a24;
    sp_simd_t p1 = a25;
    sp_simd_t p2 = a26;
    sp_simd_t p3 = a27;
    sp_simd_t p4 = a28;
    sp_simd_t p5 = a29;
    sp_simd_t p6 = a30;
    sp_simd_t p7 = a31;

    t0 = sp_ntt_add_simd0(p4, p5, p);
    t1 = sp_ntt_sub_simd0(p4, p5, p);
    t2 = sp_ntt_add_simd0(p6, p7, p);
    t3 = sp_ntt_sub_simd0(p6, p7, p);
    t4 = sp_ntt_add_simd0(t0, t2, p);
    t5 = sp_ntt_sub_simd0(t0, t2, p);
    t6 = sp_ntt_add_simd0(t1, t3, p);
    t7 = sp_ntt_sub_simd0(t1, t3, p);

    t0 = sp_ntt_add_simd0(p0, p2, p);
    t1 = sp_ntt_sub_simd0(p0, p2, p);
    t2 = sp_ntt_add_simd0(p1, p3, p);
    t3 = sp_ntt_sub_simd0(p1, p3, p);

    sp_simd_store(t7, x + j03);
    sp_simd_store(t0, x + j08);
    sp_simd_store(t5, x + j13);
    sp_simd_store(t2, x + j18);
    sp_simd_store(t6, x + j23);
    sp_simd_store(t1, x + j28);
    sp_simd_store(t4, x + j33);
    sp_simd_store(t3, x + j38);
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t p0 = a32;
    sp_simd_t p1 = a33;
    sp_simd_t p2 = a34;
    sp_simd_t p3 = a35;
    sp_simd_t p4 = a36;
    sp_simd_t p5 = a37;
    sp_simd_t p6 = a38;
    sp_simd_t p7 = a39;

    t0 = sp_ntt_add_simd0(p4, p5, p);
    t1 = sp_ntt_sub_simd0(p4, p5, p);
    t2 = sp_ntt_add_simd0(p6, p7, p);
    t3 = sp_ntt_sub_simd0(p6, p7, p);
    t4 = sp_ntt_add_simd0(t0, t2, p);
    t5 = sp_ntt_sub_simd0(t0, t2, p);
    t6 = sp_ntt_add_simd0(t1, t3, p);
    t7 = sp_ntt_sub_simd0(t1, t3, p);

    t0 = sp_ntt_add_simd0(p0, p2, p);
    t1 = sp_ntt_sub_simd0(p0, p2, p);
    t2 = sp_ntt_add_simd0(p1, p3, p);
    t3 = sp_ntt_sub_simd0(p1, p3, p);

    sp_simd_store(t1, x + j04);
    sp_simd_store(t4, x + j09);
    sp_simd_store(t3, x + j14);
    sp_simd_store(t7, x + j19);
    sp_simd_store(t0, x + j24);
    sp_simd_store(t5, x + j29);
    sp_simd_store(t2, x + j34);
    sp_simd_store(t6, x + j39);
  }
}

DECLARE_CORE_ROUTINES_SIMD(40)
