#include "ntt/ntt-impl-simd.h"

#define NC 54

static const uint8_t ntt35_fixed_const[NC] = {1};

extern void X(ntt35_init)(spv_t out, sp_t p, sp_t d, sp_t primroot, 
                            sp_t order, sp_t perm);


static void
ntt35_run_core_simd(spv_t in, spv_size_t istride, spv_size_t idist,
		spv_t out, spv_size_t ostride, spv_size_t odist,
		sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  sp_simd_t a00, a01, a02, a03, a04, a05, a06, 
       a07, a08, a09, a10, a11, a12, a13, 
       a14, a15, a16, a17, a18, a19, a20, 
       a21, a22, a23, a24, a25, a26, a27, 
       a28, a29, a30, a31, a32, a33, a34, 
       a35, a36, a37, a38, a39, a40, a41;

  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x0 = sp_simd_gather(in +  0 * istride, idist, vsize);
    sp_simd_t x1 = sp_simd_gather(in +  7 * istride, idist, vsize);
    sp_simd_t x4 = sp_simd_gather(in + 14 * istride, idist, vsize);
    sp_simd_t x2 = sp_simd_gather(in + 21 * istride, idist, vsize);
    sp_simd_t x3 = sp_simd_gather(in + 28 * istride, idist, vsize);

    t1 = sp_ntt_add_simd(x1, x3, p);
    t3 = sp_ntt_sub_simd(x1, x3, p);
    t2 = sp_ntt_add_simd(x2, x4, p);
    t4 = sp_ntt_sub_simd(x2, x4, p);

    a01 = sp_ntt_add_simd(t1, t2, p);
    a02 = sp_ntt_sub_simd(t1, t2, p);
    a03 = t3;
    a04 = t4;
    a05 = sp_ntt_add_simd(t3, t4, p);

    a00 = sp_ntt_add_simd(x0, a01, p);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x0 = sp_simd_gather(in +  5 * istride, idist, vsize);
    sp_simd_t x1 = sp_simd_gather(in + 12 * istride, idist, vsize);
    sp_simd_t x4 = sp_simd_gather(in + 19 * istride, idist, vsize);
    sp_simd_t x2 = sp_simd_gather(in + 26 * istride, idist, vsize);
    sp_simd_t x3 = sp_simd_gather(in + 33 * istride, idist, vsize);

    t1 = sp_ntt_add_simd(x1, x3, p);
    t3 = sp_ntt_sub_simd(x1, x3, p);
    t2 = sp_ntt_add_simd(x2, x4, p);
    t4 = sp_ntt_sub_simd(x2, x4, p);

    a07 = sp_ntt_add_simd(t1, t2, p);
    a08 = sp_ntt_sub_simd(t1, t2, p);
    a09 = t3;
    a10 = t4;
    a11 = sp_ntt_add_simd(t3, t4, p);

    a06 = sp_ntt_add_simd(x0, a07, p);
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x3 = sp_simd_gather(in +  3 * istride, idist, vsize);
    sp_simd_t x0 = sp_simd_gather(in + 10 * istride, idist, vsize);
    sp_simd_t x1 = sp_simd_gather(in + 17 * istride, idist, vsize);
    sp_simd_t x4 = sp_simd_gather(in + 24 * istride, idist, vsize);
    sp_simd_t x2 = sp_simd_gather(in + 31 * istride, idist, vsize);

    t1 = sp_ntt_add_simd(x1, x3, p);
    t3 = sp_ntt_sub_simd(x1, x3, p);
    t2 = sp_ntt_add_simd(x2, x4, p);
    t4 = sp_ntt_sub_simd(x2, x4, p);

    a13 = sp_ntt_add_simd(t1, t2, p);
    a14 = sp_ntt_sub_simd(t1, t2, p);
    a15 = t3;
    a16 = t4;
    a17 = sp_ntt_add_simd(t3, t4, p);

    a12 = sp_ntt_add_simd(x0, a13, p);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x2 = sp_simd_gather(in +  1 * istride, idist, vsize);
    sp_simd_t x3 = sp_simd_gather(in +  8 * istride, idist, vsize);
    sp_simd_t x0 = sp_simd_gather(in + 15 * istride, idist, vsize);
    sp_simd_t x1 = sp_simd_gather(in + 22 * istride, idist, vsize);
    sp_simd_t x4 = sp_simd_gather(in + 29 * istride, idist, vsize);

    t1 = sp_ntt_add_simd(x1, x3, p);
    t3 = sp_ntt_sub_simd(x1, x3, p);
    t2 = sp_ntt_add_simd(x2, x4, p);
    t4 = sp_ntt_sub_simd(x2, x4, p);

    a19 = sp_ntt_add_simd(t1, t2, p);
    a20 = sp_ntt_sub_simd(t1, t2, p);
    a21 = t3;
    a22 = t4;
    a23 = sp_ntt_add_simd(t3, t4, p);

    a18 = sp_ntt_add_simd(x0, a19, p);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x2 = sp_simd_gather(in +  6 * istride, idist, vsize);
    sp_simd_t x3 = sp_simd_gather(in + 13 * istride, idist, vsize);
    sp_simd_t x0 = sp_simd_gather(in + 20 * istride, idist, vsize);
    sp_simd_t x1 = sp_simd_gather(in + 27 * istride, idist, vsize);
    sp_simd_t x4 = sp_simd_gather(in + 34 * istride, idist, vsize);

    t1 = sp_ntt_add_simd(x1, x3, p);
    t3 = sp_ntt_sub_simd(x1, x3, p);
    t2 = sp_ntt_add_simd(x2, x4, p);
    t4 = sp_ntt_sub_simd(x2, x4, p);

    a25 = sp_ntt_add_simd(t1, t2, p);
    a26 = sp_ntt_sub_simd(t1, t2, p);
    a27 = t3;
    a28 = t4;
    a29 = sp_ntt_add_simd(t3, t4, p);

    a24 = sp_ntt_add_simd(x0, a25, p);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x4 = sp_simd_gather(in +  4 * istride, idist, vsize);
    sp_simd_t x2 = sp_simd_gather(in + 11 * istride, idist, vsize);
    sp_simd_t x3 = sp_simd_gather(in + 18 * istride, idist, vsize);
    sp_simd_t x0 = sp_simd_gather(in + 25 * istride, idist, vsize);
    sp_simd_t x1 = sp_simd_gather(in + 32 * istride, idist, vsize);

    t1 = sp_ntt_add_simd(x1, x3, p);
    t3 = sp_ntt_sub_simd(x1, x3, p);
    t2 = sp_ntt_add_simd(x2, x4, p);
    t4 = sp_ntt_sub_simd(x2, x4, p);

    a31 = sp_ntt_add_simd(t1, t2, p);
    a32 = sp_ntt_sub_simd(t1, t2, p);
    a33 = t3;
    a34 = t4;
    a35 = sp_ntt_add_simd(t3, t4, p);

    a30 = sp_ntt_add_simd(x0, a31, p);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x1 = sp_simd_gather(in +  2 * istride, idist, vsize);
    sp_simd_t x4 = sp_simd_gather(in +  9 * istride, idist, vsize);
    sp_simd_t x2 = sp_simd_gather(in + 16 * istride, idist, vsize);
    sp_simd_t x3 = sp_simd_gather(in + 23 * istride, idist, vsize);
    sp_simd_t x0 = sp_simd_gather(in + 30 * istride, idist, vsize);

    t1 = sp_ntt_add_simd(x1, x3, p);
    t3 = sp_ntt_sub_simd(x1, x3, p);
    t2 = sp_ntt_add_simd(x2, x4, p);
    t4 = sp_ntt_sub_simd(x2, x4, p);

    a37 = sp_ntt_add_simd(t1, t2, p);
    a38 = sp_ntt_sub_simd(t1, t2, p);
    a39 = t3;
    a40 = t4;
    a41 = sp_ntt_add_simd(t3, t4, p);

    a36 = sp_ntt_add_simd(x0, a37, p);
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a00;
    sp_simd_t x1 = a06;
    sp_simd_t x5 = a12;
    sp_simd_t x6 = a18;
    sp_simd_t x3 = a24;
    sp_simd_t x2 = a30;
    sp_simd_t x4 = a36;

    p1 = sp_ntt_add_simd(x1, x4, p);
    p2 = sp_ntt_add_simd(x2, x5, p);
    p3 = sp_ntt_add_simd(x3, x6, p);
    p4 = sp_ntt_sub_simd(x1, x4, p);
    p5 = sp_ntt_sub_simd(x2, x5, p);
    p6 = sp_ntt_sub_simd(x3, x6, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t1 = sp_ntt_add_simd(t1, p3, p);
    t2 = sp_ntt_sub_simd(p1, p3, p);
    t3 = sp_ntt_sub_simd(p2, p3, p);
    t4 = sp_ntt_sub_simd(p4, p5, p);
    t4 = sp_ntt_add_partial_simd(t4, p6, p);
    t5 = sp_ntt_sub_simd(p4, p6, p);
    t6 = sp_ntt_add_simd(p5, p6, p);

    p0 = sp_ntt_add_simd(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd(t5, t6, p);

    p1 = sp_ntt_mul_simd(p1, ntt_const[1], ntt_const[NC+1], p);
    p2 = sp_ntt_mul_simd(p2, ntt_const[2], ntt_const[NC+2], p);
    p3 = sp_ntt_mul_simd(p3, ntt_const[3], ntt_const[NC+3], p);
    p4 = sp_ntt_mul_simd(p4, ntt_const[4], ntt_const[NC+4], p);
    p5 = sp_ntt_mul_simd(p5, ntt_const[5], ntt_const[NC+5], p);
    p6 = sp_ntt_mul_simd(p6, ntt_const[6], ntt_const[NC+6], p);
    p7 = sp_ntt_mul_simd(p7, ntt_const[7], ntt_const[NC+7], p);
    p8 = sp_ntt_mul_simd(p8, ntt_const[8], ntt_const[NC+8], p);

    t1 = sp_ntt_add_simd(p0, p1, p);
    t2 = sp_ntt_add_simd(p2, p4, p);
    t3 = sp_ntt_add_simd(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd(p6, p8, p);
    t6 = sp_ntt_add_simd(p7, p8, p);

    p1 = sp_ntt_add_simd(t1, t2, p);
    p2 = sp_ntt_add_simd(t1, t3, p);
    p3 = sp_ntt_sub_simd(t1, t2, p);
    p3 = sp_ntt_sub_simd(p3, t3, p);
    p4 = sp_ntt_add_simd(t4, t5, p);
    p5 = sp_ntt_sub_simd(t6, t4, p);
    p6 = sp_ntt_sub_simd(t4, t5, p);
    p6 = sp_ntt_add_simd(p6, t6, p);

    t1 = sp_ntt_add_simd(p1, p4, p);
    t2 = sp_ntt_add_simd(p2, p5, p);
    t3 = sp_ntt_add_simd(p3, p6, p);
    t4 = sp_ntt_sub_simd(p1, p4, p);
    t5 = sp_ntt_sub_simd(p2, p5, p);
    t6 = sp_ntt_sub_simd(p3, p6, p);

    a00 = p0;
    a06 = t6;
    a12 = t4;
    a18 = t5;
    a24 = t2;
    a30 = t1;
    a36 = t3;
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a01;
    sp_simd_t x1 = a07;
    sp_simd_t x5 = a13;
    sp_simd_t x6 = a19;
    sp_simd_t x3 = a25;
    sp_simd_t x2 = a31;
    sp_simd_t x4 = a37;

    p1 = sp_ntt_add_simd(x1, x4, p);
    p2 = sp_ntt_add_simd(x2, x5, p);
    p3 = sp_ntt_add_simd(x3, x6, p);
    p4 = sp_ntt_sub_simd(x1, x4, p);
    p5 = sp_ntt_sub_simd(x2, x5, p);
    p6 = sp_ntt_sub_simd(x3, x6, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t1 = sp_ntt_add_simd(t1, p3, p);
    t2 = sp_ntt_sub_simd(p1, p3, p);
    t3 = sp_ntt_sub_simd(p2, p3, p);
    t4 = sp_ntt_sub_simd(p4, p5, p);
    t4 = sp_ntt_add_partial_simd(t4, p6, p);
    t5 = sp_ntt_sub_simd(p4, p6, p);
    t6 = sp_ntt_add_simd(p5, p6, p);

    p0 = sp_ntt_add_partial_simd(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd(t5, t6, p);

    p0 = sp_ntt_mul_simd(p0, ntt_const[9], ntt_const[NC+9], p);
    p1 = sp_ntt_mul_simd(p1, ntt_const[10], ntt_const[NC+10], p);
    p2 = sp_ntt_mul_simd(p2, ntt_const[11], ntt_const[NC+11], p);
    p3 = sp_ntt_mul_simd(p3, ntt_const[12], ntt_const[NC+12], p);
    p4 = sp_ntt_mul_simd(p4, ntt_const[13], ntt_const[NC+13], p);
    p5 = sp_ntt_mul_simd(p5, ntt_const[14], ntt_const[NC+14], p);
    p6 = sp_ntt_mul_simd(p6, ntt_const[15], ntt_const[NC+15], p);
    p7 = sp_ntt_mul_simd(p7, ntt_const[16], ntt_const[NC+16], p);
    p8 = sp_ntt_mul_simd(p8, ntt_const[17], ntt_const[NC+17], p);

    t1 = sp_ntt_add_simd(p0, p1, p);
    t2 = sp_ntt_add_simd(p2, p4, p);
    t3 = sp_ntt_add_simd(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd(p6, p8, p);
    t6 = sp_ntt_add_simd(p7, p8, p);

    p1 = sp_ntt_add_simd(t1, t2, p);
    p2 = sp_ntt_add_simd(t1, t3, p);
    p3 = sp_ntt_sub_simd(t1, t2, p);
    p3 = sp_ntt_sub_simd(p3, t3, p);
    p4 = sp_ntt_add_simd(t4, t5, p);
    p5 = sp_ntt_sub_simd(t6, t4, p);
    p6 = sp_ntt_sub_simd(t4, t5, p);
    p6 = sp_ntt_add_simd(p6, t6, p);

    t1 = sp_ntt_add_simd(p1, p4, p);
    t2 = sp_ntt_add_simd(p2, p5, p);
    t3 = sp_ntt_add_simd(p3, p6, p);
    t4 = sp_ntt_sub_simd(p1, p4, p);
    t5 = sp_ntt_sub_simd(p2, p5, p);
    t6 = sp_ntt_sub_simd(p3, p6, p);

    a01 = p0;
    a07 = t6;
    a13 = t4;
    a19 = t5;
    a25 = t2;
    a31 = t1;
    a37 = t3;
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a02;
    sp_simd_t x1 = a08;
    sp_simd_t x5 = a14;
    sp_simd_t x6 = a20;
    sp_simd_t x3 = a26;
    sp_simd_t x2 = a32;
    sp_simd_t x4 = a38;

    p1 = sp_ntt_add_simd(x1, x4, p);
    p2 = sp_ntt_add_simd(x2, x5, p);
    p3 = sp_ntt_add_simd(x3, x6, p);
    p4 = sp_ntt_sub_simd(x1, x4, p);
    p5 = sp_ntt_sub_simd(x2, x5, p);
    p6 = sp_ntt_sub_simd(x3, x6, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t1 = sp_ntt_add_simd(t1, p3, p);
    t2 = sp_ntt_sub_simd(p1, p3, p);
    t3 = sp_ntt_sub_simd(p2, p3, p);
    t4 = sp_ntt_sub_simd(p4, p5, p);
    t4 = sp_ntt_add_partial_simd(t4, p6, p);
    t5 = sp_ntt_sub_simd(p4, p6, p);
    t6 = sp_ntt_add_simd(p5, p6, p);

    p0 = sp_ntt_add_partial_simd(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd(t5, t6, p);

    p0 = sp_ntt_mul_simd(p0, ntt_const[18], ntt_const[NC+18], p);
    p1 = sp_ntt_mul_simd(p1, ntt_const[19], ntt_const[NC+19], p);
    p2 = sp_ntt_mul_simd(p2, ntt_const[20], ntt_const[NC+20], p);
    p3 = sp_ntt_mul_simd(p3, ntt_const[21], ntt_const[NC+21], p);
    p4 = sp_ntt_mul_simd(p4, ntt_const[22], ntt_const[NC+22], p);
    p5 = sp_ntt_mul_simd(p5, ntt_const[23], ntt_const[NC+23], p);
    p6 = sp_ntt_mul_simd(p6, ntt_const[24], ntt_const[NC+24], p);
    p7 = sp_ntt_mul_simd(p7, ntt_const[25], ntt_const[NC+25], p);
    p8 = sp_ntt_mul_simd(p8, ntt_const[26], ntt_const[NC+26], p);

    t1 = sp_ntt_add_simd(p0, p1, p);
    t2 = sp_ntt_add_simd(p2, p4, p);
    t3 = sp_ntt_add_simd(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd(p6, p8, p);
    t6 = sp_ntt_add_simd(p7, p8, p);

    p1 = sp_ntt_add_simd(t1, t2, p);
    p2 = sp_ntt_add_simd(t1, t3, p);
    p3 = sp_ntt_sub_simd(t1, t2, p);
    p3 = sp_ntt_sub_simd(p3, t3, p);
    p4 = sp_ntt_add_simd(t4, t5, p);
    p5 = sp_ntt_sub_simd(t6, t4, p);
    p6 = sp_ntt_sub_simd(t4, t5, p);
    p6 = sp_ntt_add_simd(p6, t6, p);

    t1 = sp_ntt_add_simd(p1, p4, p);
    t2 = sp_ntt_add_simd(p2, p5, p);
    t3 = sp_ntt_add_simd(p3, p6, p);
    t4 = sp_ntt_sub_simd(p1, p4, p);
    t5 = sp_ntt_sub_simd(p2, p5, p);
    t6 = sp_ntt_sub_simd(p3, p6, p);

    a02 = p0;
    a08 = t6;
    a14 = t4;
    a20 = t5;
    a26 = t2;
    a32 = t1;
    a38 = t3;
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a03;
    sp_simd_t x1 = a09;
    sp_simd_t x5 = a15;
    sp_simd_t x6 = a21;
    sp_simd_t x3 = a27;
    sp_simd_t x2 = a33;
    sp_simd_t x4 = a39;

    p1 = sp_ntt_add_simd(x1, x4, p);
    p2 = sp_ntt_add_simd(x2, x5, p);
    p3 = sp_ntt_add_simd(x3, x6, p);
    p4 = sp_ntt_sub_simd(x1, x4, p);
    p5 = sp_ntt_sub_simd(x2, x5, p);
    p6 = sp_ntt_sub_simd(x3, x6, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t1 = sp_ntt_add_simd(t1, p3, p);
    t2 = sp_ntt_sub_simd(p1, p3, p);
    t3 = sp_ntt_sub_simd(p2, p3, p);
    t4 = sp_ntt_sub_simd(p4, p5, p);
    t4 = sp_ntt_add_partial_simd(t4, p6, p);
    t5 = sp_ntt_sub_simd(p4, p6, p);
    t6 = sp_ntt_add_simd(p5, p6, p);

    p0 = sp_ntt_add_partial_simd(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd(t5, t6, p);

    p0 = sp_ntt_mul_simd(p0, ntt_const[27], ntt_const[NC+27], p);
    p1 = sp_ntt_mul_simd(p1, ntt_const[28], ntt_const[NC+28], p);
    p2 = sp_ntt_mul_simd(p2, ntt_const[29], ntt_const[NC+29], p);
    p3 = sp_ntt_mul_simd(p3, ntt_const[30], ntt_const[NC+30], p);
    p4 = sp_ntt_mul_simd(p4, ntt_const[31], ntt_const[NC+31], p);
    p5 = sp_ntt_mul_simd(p5, ntt_const[32], ntt_const[NC+32], p);
    p6 = sp_ntt_mul_simd(p6, ntt_const[33], ntt_const[NC+33], p);
    p7 = sp_ntt_mul_simd(p7, ntt_const[34], ntt_const[NC+34], p);
    p8 = sp_ntt_mul_simd(p8, ntt_const[35], ntt_const[NC+35], p);

    t1 = sp_ntt_add_simd(p0, p1, p);
    t2 = sp_ntt_add_simd(p2, p4, p);
    t3 = sp_ntt_add_simd(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd(p6, p8, p);
    t6 = sp_ntt_add_simd(p7, p8, p);

    p1 = sp_ntt_add_simd(t1, t2, p);
    p2 = sp_ntt_add_simd(t1, t3, p);
    p3 = sp_ntt_sub_simd(t1, t2, p);
    p3 = sp_ntt_sub_simd(p3, t3, p);
    p4 = sp_ntt_add_simd(t4, t5, p);
    p5 = sp_ntt_sub_simd(t6, t4, p);
    p6 = sp_ntt_sub_simd(t4, t5, p);
    p6 = sp_ntt_add_simd(p6, t6, p);

    t1 = sp_ntt_add_simd(p1, p4, p);
    t2 = sp_ntt_add_simd(p2, p5, p);
    t3 = sp_ntt_add_simd(p3, p6, p);
    t4 = sp_ntt_sub_simd(p1, p4, p);
    t5 = sp_ntt_sub_simd(p2, p5, p);
    t6 = sp_ntt_sub_simd(p3, p6, p);

    a03 = p0;
    a09 = t6;
    a15 = t4;
    a21 = t5;
    a27 = t2;
    a33 = t1;
    a39 = t3;
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a04;
    sp_simd_t x1 = a10;
    sp_simd_t x5 = a16;
    sp_simd_t x6 = a22;
    sp_simd_t x3 = a28;
    sp_simd_t x2 = a34;
    sp_simd_t x4 = a40;

    p1 = sp_ntt_add_simd(x1, x4, p);
    p2 = sp_ntt_add_simd(x2, x5, p);
    p3 = sp_ntt_add_simd(x3, x6, p);
    p4 = sp_ntt_sub_simd(x1, x4, p);
    p5 = sp_ntt_sub_simd(x2, x5, p);
    p6 = sp_ntt_sub_simd(x3, x6, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t1 = sp_ntt_add_simd(t1, p3, p);
    t2 = sp_ntt_sub_simd(p1, p3, p);
    t3 = sp_ntt_sub_simd(p2, p3, p);
    t4 = sp_ntt_sub_simd(p4, p5, p);
    t4 = sp_ntt_add_partial_simd(t4, p6, p);
    t5 = sp_ntt_sub_simd(p4, p6, p);
    t6 = sp_ntt_add_simd(p5, p6, p);

    p0 = sp_ntt_add_partial_simd(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd(t5, t6, p);

    p0 = sp_ntt_mul_simd(p0, ntt_const[36], ntt_const[NC+36], p);
    p1 = sp_ntt_mul_simd(p1, ntt_const[37], ntt_const[NC+37], p);
    p2 = sp_ntt_mul_simd(p2, ntt_const[38], ntt_const[NC+38], p);
    p3 = sp_ntt_mul_simd(p3, ntt_const[39], ntt_const[NC+39], p);
    p4 = sp_ntt_mul_simd(p4, ntt_const[40], ntt_const[NC+40], p);
    p5 = sp_ntt_mul_simd(p5, ntt_const[41], ntt_const[NC+41], p);
    p6 = sp_ntt_mul_simd(p6, ntt_const[42], ntt_const[NC+42], p);
    p7 = sp_ntt_mul_simd(p7, ntt_const[43], ntt_const[NC+43], p);
    p8 = sp_ntt_mul_simd(p8, ntt_const[44], ntt_const[NC+44], p);

    t1 = sp_ntt_add_simd(p0, p1, p);
    t2 = sp_ntt_add_simd(p2, p4, p);
    t3 = sp_ntt_add_simd(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd(p6, p8, p);
    t6 = sp_ntt_add_simd(p7, p8, p);

    p1 = sp_ntt_add_simd(t1, t2, p);
    p2 = sp_ntt_add_simd(t1, t3, p);
    p3 = sp_ntt_sub_simd(t1, t2, p);
    p3 = sp_ntt_sub_simd(p3, t3, p);
    p4 = sp_ntt_add_simd(t4, t5, p);
    p5 = sp_ntt_sub_simd(t6, t4, p);
    p6 = sp_ntt_sub_simd(t4, t5, p);
    p6 = sp_ntt_add_simd(p6, t6, p);

    t1 = sp_ntt_add_simd(p1, p4, p);
    t2 = sp_ntt_add_simd(p2, p5, p);
    t3 = sp_ntt_add_simd(p3, p6, p);
    t4 = sp_ntt_sub_simd(p1, p4, p);
    t5 = sp_ntt_sub_simd(p2, p5, p);
    t6 = sp_ntt_sub_simd(p3, p6, p);

    a04 = p0;
    a10 = t6;
    a16 = t4;
    a22 = t5;
    a28 = t2;
    a34 = t1;
    a40 = t3;
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a05;
    sp_simd_t x1 = a11;
    sp_simd_t x5 = a17;
    sp_simd_t x6 = a23;
    sp_simd_t x3 = a29;
    sp_simd_t x2 = a35;
    sp_simd_t x4 = a41;

    p1 = sp_ntt_add_simd(x1, x4, p);
    p2 = sp_ntt_add_simd(x2, x5, p);
    p3 = sp_ntt_add_simd(x3, x6, p);
    p4 = sp_ntt_sub_simd(x1, x4, p);
    p5 = sp_ntt_sub_simd(x2, x5, p);
    p6 = sp_ntt_sub_simd(x3, x6, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t1 = sp_ntt_add_simd(t1, p3, p);
    t2 = sp_ntt_sub_simd(p1, p3, p);
    t3 = sp_ntt_sub_simd(p2, p3, p);
    t4 = sp_ntt_sub_simd(p4, p5, p);
    t4 = sp_ntt_add_partial_simd(t4, p6, p);
    t5 = sp_ntt_sub_simd(p4, p6, p);
    t6 = sp_ntt_add_simd(p5, p6, p);

    p0 = sp_ntt_add_partial_simd(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd(t5, t6, p);

    p0 = sp_ntt_mul_simd(p0, ntt_const[45], ntt_const[NC+45], p);
    p1 = sp_ntt_mul_simd(p1, ntt_const[46], ntt_const[NC+46], p);
    p2 = sp_ntt_mul_simd(p2, ntt_const[47], ntt_const[NC+47], p);
    p3 = sp_ntt_mul_simd(p3, ntt_const[48], ntt_const[NC+48], p);
    p4 = sp_ntt_mul_simd(p4, ntt_const[49], ntt_const[NC+49], p);
    p5 = sp_ntt_mul_simd(p5, ntt_const[50], ntt_const[NC+50], p);
    p6 = sp_ntt_mul_simd(p6, ntt_const[51], ntt_const[NC+51], p);
    p7 = sp_ntt_mul_simd(p7, ntt_const[52], ntt_const[NC+52], p);
    p8 = sp_ntt_mul_simd(p8, ntt_const[53], ntt_const[NC+53], p);

    t1 = sp_ntt_add_simd(p0, p1, p);
    t2 = sp_ntt_add_simd(p2, p4, p);
    t3 = sp_ntt_add_simd(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd(p6, p8, p);
    t6 = sp_ntt_add_simd(p7, p8, p);

    p1 = sp_ntt_add_simd(t1, t2, p);
    p2 = sp_ntt_add_simd(t1, t3, p);
    p3 = sp_ntt_sub_simd(t1, t2, p);
    p3 = sp_ntt_sub_simd(p3, t3, p);
    p4 = sp_ntt_add_simd(t4, t5, p);
    p5 = sp_ntt_sub_simd(t6, t4, p);
    p6 = sp_ntt_sub_simd(t4, t5, p);
    p6 = sp_ntt_add_simd(p6, t6, p);

    t1 = sp_ntt_add_simd(p1, p4, p);
    t2 = sp_ntt_add_simd(p2, p5, p);
    t3 = sp_ntt_add_simd(p3, p6, p);
    t4 = sp_ntt_sub_simd(p1, p4, p);
    t5 = sp_ntt_sub_simd(p2, p5, p);
    t6 = sp_ntt_sub_simd(p3, p6, p);

    a05 = p0;
    a11 = t6;
    a17 = t4;
    a23 = t5;
    a29 = t2;
    a35 = t1;
    a41 = t3;
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a00;
    sp_simd_t p1 = a01;
    sp_simd_t p2 = a02;
    sp_simd_t p3 = a03;
    sp_simd_t p4 = a04;
    sp_simd_t p5 = a05;

    p1 = sp_ntt_add_simd(p0, p1, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t2 = sp_ntt_sub_simd(p1, p2, p);
    t3 = sp_ntt_add_simd(p3, p5, p);
    t4 = sp_ntt_add_simd(p4, p5, p);

    p1 = sp_ntt_add_simd(t1, t3, p);
    p2 = sp_ntt_add_simd(t2, t4, p);
    p3 = sp_ntt_sub_simd(t1, t3, p);
    p4 = sp_ntt_sub_simd(t2, t4, p);

    sp_simd_scatter(p0, out +  0 * ostride, odist, vsize);
    sp_simd_scatter(p3, out +  7 * ostride, odist, vsize);
    sp_simd_scatter(p2, out + 14 * ostride, odist, vsize);
    sp_simd_scatter(p4, out + 21 * ostride, odist, vsize);
    sp_simd_scatter(p1, out + 28 * ostride, odist, vsize);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a06;
    sp_simd_t p1 = a07;
    sp_simd_t p2 = a08;
    sp_simd_t p3 = a09;
    sp_simd_t p4 = a10;
    sp_simd_t p5 = a11;

    p1 = sp_ntt_add_simd(p0, p1, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t2 = sp_ntt_sub_simd(p1, p2, p);
    t3 = sp_ntt_add_simd(p3, p5, p);
    t4 = sp_ntt_add_simd(p4, p5, p);

    p1 = sp_ntt_add_simd(t1, t3, p);
    p2 = sp_ntt_add_simd(t2, t4, p);
    p3 = sp_ntt_sub_simd(t1, t3, p);
    p4 = sp_ntt_sub_simd(t2, t4, p);

    sp_simd_scatter(p4, out +  1 * ostride, odist, vsize);
    sp_simd_scatter(p1, out +  8 * ostride, odist, vsize);
    sp_simd_scatter(p0, out + 15 * ostride, odist, vsize);
    sp_simd_scatter(p3, out + 22 * ostride, odist, vsize);
    sp_simd_scatter(p2, out + 29 * ostride, odist, vsize);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a12;
    sp_simd_t p1 = a13;
    sp_simd_t p2 = a14;
    sp_simd_t p3 = a15;
    sp_simd_t p4 = a16;
    sp_simd_t p5 = a17;

    p1 = sp_ntt_add_simd(p0, p1, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t2 = sp_ntt_sub_simd(p1, p2, p);
    t3 = sp_ntt_add_simd(p3, p5, p);
    t4 = sp_ntt_add_simd(p4, p5, p);

    p1 = sp_ntt_add_simd(t1, t3, p);
    p2 = sp_ntt_add_simd(t2, t4, p);
    p3 = sp_ntt_sub_simd(t1, t3, p);
    p4 = sp_ntt_sub_simd(t2, t4, p);

    sp_simd_scatter(p3, out +  2 * ostride, odist, vsize);
    sp_simd_scatter(p2, out +  9 * ostride, odist, vsize);
    sp_simd_scatter(p4, out + 16 * ostride, odist, vsize);
    sp_simd_scatter(p1, out + 23 * ostride, odist, vsize);
    sp_simd_scatter(p0, out + 30 * ostride, odist, vsize);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a18;
    sp_simd_t p1 = a19;
    sp_simd_t p2 = a20;
    sp_simd_t p3 = a21;
    sp_simd_t p4 = a22;
    sp_simd_t p5 = a23;

    p1 = sp_ntt_add_simd(p0, p1, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t2 = sp_ntt_sub_simd(p1, p2, p);
    t3 = sp_ntt_add_simd(p3, p5, p);
    t4 = sp_ntt_add_simd(p4, p5, p);

    p1 = sp_ntt_add_simd(t1, t3, p);
    p2 = sp_ntt_add_simd(t2, t4, p);
    p3 = sp_ntt_sub_simd(t1, t3, p);
    p4 = sp_ntt_sub_simd(t2, t4, p);

    sp_simd_scatter(p1, out +  3 * ostride, odist, vsize);
    sp_simd_scatter(p0, out + 10 * ostride, odist, vsize);
    sp_simd_scatter(p3, out + 17 * ostride, odist, vsize);
    sp_simd_scatter(p2, out + 24 * ostride, odist, vsize);
    sp_simd_scatter(p4, out + 31 * ostride, odist, vsize);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a24;
    sp_simd_t p1 = a25;
    sp_simd_t p2 = a26;
    sp_simd_t p3 = a27;
    sp_simd_t p4 = a28;
    sp_simd_t p5 = a29;

    p1 = sp_ntt_add_simd(p0, p1, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t2 = sp_ntt_sub_simd(p1, p2, p);
    t3 = sp_ntt_add_simd(p3, p5, p);
    t4 = sp_ntt_add_simd(p4, p5, p);

    p1 = sp_ntt_add_simd(t1, t3, p);
    p2 = sp_ntt_add_simd(t2, t4, p);
    p3 = sp_ntt_sub_simd(t1, t3, p);
    p4 = sp_ntt_sub_simd(t2, t4, p);

    sp_simd_scatter(p2, out +  4 * ostride, odist, vsize);
    sp_simd_scatter(p4, out + 11 * ostride, odist, vsize);
    sp_simd_scatter(p1, out + 18 * ostride, odist, vsize);
    sp_simd_scatter(p0, out + 25 * ostride, odist, vsize);
    sp_simd_scatter(p3, out + 32 * ostride, odist, vsize);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a30;
    sp_simd_t p1 = a31;
    sp_simd_t p2 = a32;
    sp_simd_t p3 = a33;
    sp_simd_t p4 = a34;
    sp_simd_t p5 = a35;

    p1 = sp_ntt_add_simd(p0, p1, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t2 = sp_ntt_sub_simd(p1, p2, p);
    t3 = sp_ntt_add_simd(p3, p5, p);
    t4 = sp_ntt_add_simd(p4, p5, p);

    p1 = sp_ntt_add_simd(t1, t3, p);
    p2 = sp_ntt_add_simd(t2, t4, p);
    p3 = sp_ntt_sub_simd(t1, t3, p);
    p4 = sp_ntt_sub_simd(t2, t4, p);

    sp_simd_scatter(p0, out +  5 * ostride, odist, vsize);
    sp_simd_scatter(p3, out + 12 * ostride, odist, vsize);
    sp_simd_scatter(p2, out + 19 * ostride, odist, vsize);
    sp_simd_scatter(p4, out + 26 * ostride, odist, vsize);
    sp_simd_scatter(p1, out + 33 * ostride, odist, vsize);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a36;
    sp_simd_t p1 = a37;
    sp_simd_t p2 = a38;
    sp_simd_t p3 = a39;
    sp_simd_t p4 = a40;
    sp_simd_t p5 = a41;

    p1 = sp_ntt_add_simd(p0, p1, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t2 = sp_ntt_sub_simd(p1, p2, p);
    t3 = sp_ntt_add_simd(p3, p5, p);
    t4 = sp_ntt_add_simd(p4, p5, p);

    p1 = sp_ntt_add_simd(t1, t3, p);
    p2 = sp_ntt_add_simd(t2, t4, p);
    p3 = sp_ntt_sub_simd(t1, t3, p);
    p4 = sp_ntt_sub_simd(t2, t4, p);

    sp_simd_scatter(p4, out +  6 * ostride, odist, vsize);
    sp_simd_scatter(p1, out + 13 * ostride, odist, vsize);
    sp_simd_scatter(p0, out + 20 * ostride, odist, vsize);
    sp_simd_scatter(p3, out + 27 * ostride, odist, vsize);
    sp_simd_scatter(p2, out + 34 * ostride, odist, vsize);
  }
}

static void
ntt35_run_core_simd_interleaved(
		spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		sp_simd_t p, sp_simd_t * c)
{
  sp_simd_t a00, a01, a02, a03, a04, a05, a06, 
       a07, a08, a09, a10, a11, a12, a13, 
       a14, a15, a16, a17, a18, a19, a20, 
       a21, a22, a23, a24, a25, a26, a27, 
       a28, a29, a30, a31, a32, a33, a34, 
       a35, a36, a37, a38, a39, a40, a41;

  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x0 = sp_simd_load(in +  0 * istride);
    sp_simd_t x1 = sp_simd_load(in +  7 * istride);
    sp_simd_t x4 = sp_simd_load(in + 14 * istride);
    sp_simd_t x2 = sp_simd_load(in + 21 * istride);
    sp_simd_t x3 = sp_simd_load(in + 28 * istride);

    t1 = sp_ntt_add_simd0(x1, x3, p);
    t3 = sp_ntt_sub_simd0(x1, x3, p);
    t2 = sp_ntt_add_simd0(x2, x4, p);
    t4 = sp_ntt_sub_simd0(x2, x4, p);

    a01 = sp_ntt_add_simd0(t1, t2, p);
    a02 = sp_ntt_sub_simd0(t1, t2, p);
    a03 = t3;
    a04 = t4;
    a05 = sp_ntt_add_simd0(t3, t4, p);

    a00 = sp_ntt_add_simd0(x0, a01, p);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x0 = sp_simd_load(in +  5 * istride);
    sp_simd_t x1 = sp_simd_load(in + 12 * istride);
    sp_simd_t x4 = sp_simd_load(in + 19 * istride);
    sp_simd_t x2 = sp_simd_load(in + 26 * istride);
    sp_simd_t x3 = sp_simd_load(in + 33 * istride);

    t1 = sp_ntt_add_simd0(x1, x3, p);
    t3 = sp_ntt_sub_simd0(x1, x3, p);
    t2 = sp_ntt_add_simd0(x2, x4, p);
    t4 = sp_ntt_sub_simd0(x2, x4, p);

    a07 = sp_ntt_add_simd0(t1, t2, p);
    a08 = sp_ntt_sub_simd0(t1, t2, p);
    a09 = t3;
    a10 = t4;
    a11 = sp_ntt_add_simd0(t3, t4, p);

    a06 = sp_ntt_add_simd0(x0, a07, p);
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x3 = sp_simd_load(in +  3 * istride);
    sp_simd_t x0 = sp_simd_load(in + 10 * istride);
    sp_simd_t x1 = sp_simd_load(in + 17 * istride);
    sp_simd_t x4 = sp_simd_load(in + 24 * istride);
    sp_simd_t x2 = sp_simd_load(in + 31 * istride);

    t1 = sp_ntt_add_simd0(x1, x3, p);
    t3 = sp_ntt_sub_simd0(x1, x3, p);
    t2 = sp_ntt_add_simd0(x2, x4, p);
    t4 = sp_ntt_sub_simd0(x2, x4, p);

    a13 = sp_ntt_add_simd0(t1, t2, p);
    a14 = sp_ntt_sub_simd0(t1, t2, p);
    a15 = t3;
    a16 = t4;
    a17 = sp_ntt_add_simd0(t3, t4, p);

    a12 = sp_ntt_add_simd0(x0, a13, p);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x2 = sp_simd_load(in +  1 * istride);
    sp_simd_t x3 = sp_simd_load(in +  8 * istride);
    sp_simd_t x0 = sp_simd_load(in + 15 * istride);
    sp_simd_t x1 = sp_simd_load(in + 22 * istride);
    sp_simd_t x4 = sp_simd_load(in + 29 * istride);

    t1 = sp_ntt_add_simd0(x1, x3, p);
    t3 = sp_ntt_sub_simd0(x1, x3, p);
    t2 = sp_ntt_add_simd0(x2, x4, p);
    t4 = sp_ntt_sub_simd0(x2, x4, p);

    a19 = sp_ntt_add_simd0(t1, t2, p);
    a20 = sp_ntt_sub_simd0(t1, t2, p);
    a21 = t3;
    a22 = t4;
    a23 = sp_ntt_add_simd0(t3, t4, p);

    a18 = sp_ntt_add_simd0(x0, a19, p);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x2 = sp_simd_load(in +  6 * istride);
    sp_simd_t x3 = sp_simd_load(in + 13 * istride);
    sp_simd_t x0 = sp_simd_load(in + 20 * istride);
    sp_simd_t x1 = sp_simd_load(in + 27 * istride);
    sp_simd_t x4 = sp_simd_load(in + 34 * istride);

    t1 = sp_ntt_add_simd0(x1, x3, p);
    t3 = sp_ntt_sub_simd0(x1, x3, p);
    t2 = sp_ntt_add_simd0(x2, x4, p);
    t4 = sp_ntt_sub_simd0(x2, x4, p);

    a25 = sp_ntt_add_simd0(t1, t2, p);
    a26 = sp_ntt_sub_simd0(t1, t2, p);
    a27 = t3;
    a28 = t4;
    a29 = sp_ntt_add_simd0(t3, t4, p);

    a24 = sp_ntt_add_simd0(x0, a25, p);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x4 = sp_simd_load(in +  4 * istride);
    sp_simd_t x2 = sp_simd_load(in + 11 * istride);
    sp_simd_t x3 = sp_simd_load(in + 18 * istride);
    sp_simd_t x0 = sp_simd_load(in + 25 * istride);
    sp_simd_t x1 = sp_simd_load(in + 32 * istride);

    t1 = sp_ntt_add_simd0(x1, x3, p);
    t3 = sp_ntt_sub_simd0(x1, x3, p);
    t2 = sp_ntt_add_simd0(x2, x4, p);
    t4 = sp_ntt_sub_simd0(x2, x4, p);

    a31 = sp_ntt_add_simd0(t1, t2, p);
    a32 = sp_ntt_sub_simd0(t1, t2, p);
    a33 = t3;
    a34 = t4;
    a35 = sp_ntt_add_simd0(t3, t4, p);

    a30 = sp_ntt_add_simd0(x0, a31, p);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x1 = sp_simd_load(in +  2 * istride);
    sp_simd_t x4 = sp_simd_load(in +  9 * istride);
    sp_simd_t x2 = sp_simd_load(in + 16 * istride);
    sp_simd_t x3 = sp_simd_load(in + 23 * istride);
    sp_simd_t x0 = sp_simd_load(in + 30 * istride);

    t1 = sp_ntt_add_simd0(x1, x3, p);
    t3 = sp_ntt_sub_simd0(x1, x3, p);
    t2 = sp_ntt_add_simd0(x2, x4, p);
    t4 = sp_ntt_sub_simd0(x2, x4, p);

    a37 = sp_ntt_add_simd0(t1, t2, p);
    a38 = sp_ntt_sub_simd0(t1, t2, p);
    a39 = t3;
    a40 = t4;
    a41 = sp_ntt_add_simd0(t3, t4, p);

    a36 = sp_ntt_add_simd0(x0, a37, p);
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a00;
    sp_simd_t x1 = a06;
    sp_simd_t x5 = a12;
    sp_simd_t x6 = a18;
    sp_simd_t x3 = a24;
    sp_simd_t x2 = a30;
    sp_simd_t x4 = a36;

    p1 = sp_ntt_add_simd0(x1, x4, p);
    p2 = sp_ntt_add_simd0(x2, x5, p);
    p3 = sp_ntt_add_simd0(x3, x6, p);
    p4 = sp_ntt_sub_simd0(x1, x4, p);
    p5 = sp_ntt_sub_simd0(x2, x5, p);
    p6 = sp_ntt_sub_simd0(x3, x6, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t1 = sp_ntt_add_simd0(t1, p3, p);
    t2 = sp_ntt_sub_simd0(p1, p3, p);
    t3 = sp_ntt_sub_simd0(p2, p3, p);
    t4 = sp_ntt_sub_simd0(p4, p5, p);
    t4 = sp_ntt_add_partial_simd0(t4, p6, p);
    t5 = sp_ntt_sub_simd0(p4, p6, p);
    t6 = sp_ntt_add_simd0(p5, p6, p);

    p0 = sp_ntt_add_simd0(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd0(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd0(t5, t6, p);

    p1 = sp_ntt_mul_simd0(p1, c+2*1, p);
    p2 = sp_ntt_mul_simd0(p2, c+2*2, p);
    p3 = sp_ntt_mul_simd0(p3, c+2*3, p);
    p4 = sp_ntt_mul_simd0(p4, c+2*4, p);
    p5 = sp_ntt_mul_simd0(p5, c+2*5, p);
    p6 = sp_ntt_mul_simd0(p6, c+2*6, p);
    p7 = sp_ntt_mul_simd0(p7, c+2*7, p);
    p8 = sp_ntt_mul_simd0(p8, c+2*8, p);

    t1 = sp_ntt_add_simd0(p0, p1, p);
    t2 = sp_ntt_add_simd0(p2, p4, p);
    t3 = sp_ntt_add_simd0(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd0(p6, p8, p);
    t6 = sp_ntt_add_simd0(p7, p8, p);

    p1 = sp_ntt_add_simd0(t1, t2, p);
    p2 = sp_ntt_add_simd0(t1, t3, p);
    p3 = sp_ntt_sub_simd0(t1, t2, p);
    p3 = sp_ntt_sub_simd0(p3, t3, p);
    p4 = sp_ntt_add_simd0(t4, t5, p);
    p5 = sp_ntt_sub_simd0(t6, t4, p);
    p6 = sp_ntt_sub_simd0(t4, t5, p);
    p6 = sp_ntt_add_simd0(p6, t6, p);

    t1 = sp_ntt_add_simd0(p1, p4, p);
    t2 = sp_ntt_add_simd0(p2, p5, p);
    t3 = sp_ntt_add_simd0(p3, p6, p);
    t4 = sp_ntt_sub_simd0(p1, p4, p);
    t5 = sp_ntt_sub_simd0(p2, p5, p);
    t6 = sp_ntt_sub_simd0(p3, p6, p);

    a00 = p0;
    a06 = t6;
    a12 = t4;
    a18 = t5;
    a24 = t2;
    a30 = t1;
    a36 = t3;
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a01;
    sp_simd_t x1 = a07;
    sp_simd_t x5 = a13;
    sp_simd_t x6 = a19;
    sp_simd_t x3 = a25;
    sp_simd_t x2 = a31;
    sp_simd_t x4 = a37;

    p1 = sp_ntt_add_simd0(x1, x4, p);
    p2 = sp_ntt_add_simd0(x2, x5, p);
    p3 = sp_ntt_add_simd0(x3, x6, p);
    p4 = sp_ntt_sub_simd0(x1, x4, p);
    p5 = sp_ntt_sub_simd0(x2, x5, p);
    p6 = sp_ntt_sub_simd0(x3, x6, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t1 = sp_ntt_add_simd0(t1, p3, p);
    t2 = sp_ntt_sub_simd0(p1, p3, p);
    t3 = sp_ntt_sub_simd0(p2, p3, p);
    t4 = sp_ntt_sub_simd0(p4, p5, p);
    t4 = sp_ntt_add_partial_simd0(t4, p6, p);
    t5 = sp_ntt_sub_simd0(p4, p6, p);
    t6 = sp_ntt_add_simd0(p5, p6, p);

    p0 = sp_ntt_add_partial_simd0(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd0(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd0(t5, t6, p);

    p0 = sp_ntt_mul_simd0(p0, c+2*9, p);
    p1 = sp_ntt_mul_simd0(p1, c+2*10, p);
    p2 = sp_ntt_mul_simd0(p2, c+2*11, p);
    p3 = sp_ntt_mul_simd0(p3, c+2*12, p);
    p4 = sp_ntt_mul_simd0(p4, c+2*13, p);
    p5 = sp_ntt_mul_simd0(p5, c+2*14, p);
    p6 = sp_ntt_mul_simd0(p6, c+2*15, p);
    p7 = sp_ntt_mul_simd0(p7, c+2*16, p);
    p8 = sp_ntt_mul_simd0(p8, c+2*17, p);

    t1 = sp_ntt_add_simd0(p0, p1, p);
    t2 = sp_ntt_add_simd0(p2, p4, p);
    t3 = sp_ntt_add_simd0(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd0(p6, p8, p);
    t6 = sp_ntt_add_simd0(p7, p8, p);

    p1 = sp_ntt_add_simd0(t1, t2, p);
    p2 = sp_ntt_add_simd0(t1, t3, p);
    p3 = sp_ntt_sub_simd0(t1, t2, p);
    p3 = sp_ntt_sub_simd0(p3, t3, p);
    p4 = sp_ntt_add_simd0(t4, t5, p);
    p5 = sp_ntt_sub_simd0(t6, t4, p);
    p6 = sp_ntt_sub_simd0(t4, t5, p);
    p6 = sp_ntt_add_simd0(p6, t6, p);

    t1 = sp_ntt_add_simd0(p1, p4, p);
    t2 = sp_ntt_add_simd0(p2, p5, p);
    t3 = sp_ntt_add_simd0(p3, p6, p);
    t4 = sp_ntt_sub_simd0(p1, p4, p);
    t5 = sp_ntt_sub_simd0(p2, p5, p);
    t6 = sp_ntt_sub_simd0(p3, p6, p);

    a01 = p0;
    a07 = t6;
    a13 = t4;
    a19 = t5;
    a25 = t2;
    a31 = t1;
    a37 = t3;
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a02;
    sp_simd_t x1 = a08;
    sp_simd_t x5 = a14;
    sp_simd_t x6 = a20;
    sp_simd_t x3 = a26;
    sp_simd_t x2 = a32;
    sp_simd_t x4 = a38;

    p1 = sp_ntt_add_simd0(x1, x4, p);
    p2 = sp_ntt_add_simd0(x2, x5, p);
    p3 = sp_ntt_add_simd0(x3, x6, p);
    p4 = sp_ntt_sub_simd0(x1, x4, p);
    p5 = sp_ntt_sub_simd0(x2, x5, p);
    p6 = sp_ntt_sub_simd0(x3, x6, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t1 = sp_ntt_add_simd0(t1, p3, p);
    t2 = sp_ntt_sub_simd0(p1, p3, p);
    t3 = sp_ntt_sub_simd0(p2, p3, p);
    t4 = sp_ntt_sub_simd0(p4, p5, p);
    t4 = sp_ntt_add_partial_simd0(t4, p6, p);
    t5 = sp_ntt_sub_simd0(p4, p6, p);
    t6 = sp_ntt_add_simd0(p5, p6, p);

    p0 = sp_ntt_add_partial_simd0(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd0(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd0(t5, t6, p);

    p0 = sp_ntt_mul_simd0(p0, c+2*18, p);
    p1 = sp_ntt_mul_simd0(p1, c+2*19, p);
    p2 = sp_ntt_mul_simd0(p2, c+2*20, p);
    p3 = sp_ntt_mul_simd0(p3, c+2*21, p);
    p4 = sp_ntt_mul_simd0(p4, c+2*22, p);
    p5 = sp_ntt_mul_simd0(p5, c+2*23, p);
    p6 = sp_ntt_mul_simd0(p6, c+2*24, p);
    p7 = sp_ntt_mul_simd0(p7, c+2*25, p);
    p8 = sp_ntt_mul_simd0(p8, c+2*26, p);

    t1 = sp_ntt_add_simd0(p0, p1, p);
    t2 = sp_ntt_add_simd0(p2, p4, p);
    t3 = sp_ntt_add_simd0(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd0(p6, p8, p);
    t6 = sp_ntt_add_simd0(p7, p8, p);

    p1 = sp_ntt_add_simd0(t1, t2, p);
    p2 = sp_ntt_add_simd0(t1, t3, p);
    p3 = sp_ntt_sub_simd0(t1, t2, p);
    p3 = sp_ntt_sub_simd0(p3, t3, p);
    p4 = sp_ntt_add_simd0(t4, t5, p);
    p5 = sp_ntt_sub_simd0(t6, t4, p);
    p6 = sp_ntt_sub_simd0(t4, t5, p);
    p6 = sp_ntt_add_simd0(p6, t6, p);

    t1 = sp_ntt_add_simd0(p1, p4, p);
    t2 = sp_ntt_add_simd0(p2, p5, p);
    t3 = sp_ntt_add_simd0(p3, p6, p);
    t4 = sp_ntt_sub_simd0(p1, p4, p);
    t5 = sp_ntt_sub_simd0(p2, p5, p);
    t6 = sp_ntt_sub_simd0(p3, p6, p);

    a02 = p0;
    a08 = t6;
    a14 = t4;
    a20 = t5;
    a26 = t2;
    a32 = t1;
    a38 = t3;
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a03;
    sp_simd_t x1 = a09;
    sp_simd_t x5 = a15;
    sp_simd_t x6 = a21;
    sp_simd_t x3 = a27;
    sp_simd_t x2 = a33;
    sp_simd_t x4 = a39;

    p1 = sp_ntt_add_simd0(x1, x4, p);
    p2 = sp_ntt_add_simd0(x2, x5, p);
    p3 = sp_ntt_add_simd0(x3, x6, p);
    p4 = sp_ntt_sub_simd0(x1, x4, p);
    p5 = sp_ntt_sub_simd0(x2, x5, p);
    p6 = sp_ntt_sub_simd0(x3, x6, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t1 = sp_ntt_add_simd0(t1, p3, p);
    t2 = sp_ntt_sub_simd0(p1, p3, p);
    t3 = sp_ntt_sub_simd0(p2, p3, p);
    t4 = sp_ntt_sub_simd0(p4, p5, p);
    t4 = sp_ntt_add_partial_simd0(t4, p6, p);
    t5 = sp_ntt_sub_simd0(p4, p6, p);
    t6 = sp_ntt_add_simd0(p5, p6, p);

    p0 = sp_ntt_add_partial_simd0(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd0(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd0(t5, t6, p);

    p0 = sp_ntt_mul_simd0(p0, c+2*27, p);
    p1 = sp_ntt_mul_simd0(p1, c+2*28, p);
    p2 = sp_ntt_mul_simd0(p2, c+2*29, p);
    p3 = sp_ntt_mul_simd0(p3, c+2*30, p);
    p4 = sp_ntt_mul_simd0(p4, c+2*31, p);
    p5 = sp_ntt_mul_simd0(p5, c+2*32, p);
    p6 = sp_ntt_mul_simd0(p6, c+2*33, p);
    p7 = sp_ntt_mul_simd0(p7, c+2*34, p);
    p8 = sp_ntt_mul_simd0(p8, c+2*35, p);

    t1 = sp_ntt_add_simd0(p0, p1, p);
    t2 = sp_ntt_add_simd0(p2, p4, p);
    t3 = sp_ntt_add_simd0(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd0(p6, p8, p);
    t6 = sp_ntt_add_simd0(p7, p8, p);

    p1 = sp_ntt_add_simd0(t1, t2, p);
    p2 = sp_ntt_add_simd0(t1, t3, p);
    p3 = sp_ntt_sub_simd0(t1, t2, p);
    p3 = sp_ntt_sub_simd0(p3, t3, p);
    p4 = sp_ntt_add_simd0(t4, t5, p);
    p5 = sp_ntt_sub_simd0(t6, t4, p);
    p6 = sp_ntt_sub_simd0(t4, t5, p);
    p6 = sp_ntt_add_simd0(p6, t6, p);

    t1 = sp_ntt_add_simd0(p1, p4, p);
    t2 = sp_ntt_add_simd0(p2, p5, p);
    t3 = sp_ntt_add_simd0(p3, p6, p);
    t4 = sp_ntt_sub_simd0(p1, p4, p);
    t5 = sp_ntt_sub_simd0(p2, p5, p);
    t6 = sp_ntt_sub_simd0(p3, p6, p);

    a03 = p0;
    a09 = t6;
    a15 = t4;
    a21 = t5;
    a27 = t2;
    a33 = t1;
    a39 = t3;
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a04;
    sp_simd_t x1 = a10;
    sp_simd_t x5 = a16;
    sp_simd_t x6 = a22;
    sp_simd_t x3 = a28;
    sp_simd_t x2 = a34;
    sp_simd_t x4 = a40;

    p1 = sp_ntt_add_simd0(x1, x4, p);
    p2 = sp_ntt_add_simd0(x2, x5, p);
    p3 = sp_ntt_add_simd0(x3, x6, p);
    p4 = sp_ntt_sub_simd0(x1, x4, p);
    p5 = sp_ntt_sub_simd0(x2, x5, p);
    p6 = sp_ntt_sub_simd0(x3, x6, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t1 = sp_ntt_add_simd0(t1, p3, p);
    t2 = sp_ntt_sub_simd0(p1, p3, p);
    t3 = sp_ntt_sub_simd0(p2, p3, p);
    t4 = sp_ntt_sub_simd0(p4, p5, p);
    t4 = sp_ntt_add_partial_simd0(t4, p6, p);
    t5 = sp_ntt_sub_simd0(p4, p6, p);
    t6 = sp_ntt_add_simd0(p5, p6, p);

    p0 = sp_ntt_add_partial_simd0(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd0(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd0(t5, t6, p);

    p0 = sp_ntt_mul_simd0(p0, c+2*36, p);
    p1 = sp_ntt_mul_simd0(p1, c+2*37, p);
    p2 = sp_ntt_mul_simd0(p2, c+2*38, p);
    p3 = sp_ntt_mul_simd0(p3, c+2*39, p);
    p4 = sp_ntt_mul_simd0(p4, c+2*40, p);
    p5 = sp_ntt_mul_simd0(p5, c+2*41, p);
    p6 = sp_ntt_mul_simd0(p6, c+2*42, p);
    p7 = sp_ntt_mul_simd0(p7, c+2*43, p);
    p8 = sp_ntt_mul_simd0(p8, c+2*44, p);

    t1 = sp_ntt_add_simd0(p0, p1, p);
    t2 = sp_ntt_add_simd0(p2, p4, p);
    t3 = sp_ntt_add_simd0(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd0(p6, p8, p);
    t6 = sp_ntt_add_simd0(p7, p8, p);

    p1 = sp_ntt_add_simd0(t1, t2, p);
    p2 = sp_ntt_add_simd0(t1, t3, p);
    p3 = sp_ntt_sub_simd0(t1, t2, p);
    p3 = sp_ntt_sub_simd0(p3, t3, p);
    p4 = sp_ntt_add_simd0(t4, t5, p);
    p5 = sp_ntt_sub_simd0(t6, t4, p);
    p6 = sp_ntt_sub_simd0(t4, t5, p);
    p6 = sp_ntt_add_simd0(p6, t6, p);

    t1 = sp_ntt_add_simd0(p1, p4, p);
    t2 = sp_ntt_add_simd0(p2, p5, p);
    t3 = sp_ntt_add_simd0(p3, p6, p);
    t4 = sp_ntt_sub_simd0(p1, p4, p);
    t5 = sp_ntt_sub_simd0(p2, p5, p);
    t6 = sp_ntt_sub_simd0(p3, p6, p);

    a04 = p0;
    a10 = t6;
    a16 = t4;
    a22 = t5;
    a28 = t2;
    a34 = t1;
    a40 = t3;
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a05;
    sp_simd_t x1 = a11;
    sp_simd_t x5 = a17;
    sp_simd_t x6 = a23;
    sp_simd_t x3 = a29;
    sp_simd_t x2 = a35;
    sp_simd_t x4 = a41;

    p1 = sp_ntt_add_simd0(x1, x4, p);
    p2 = sp_ntt_add_simd0(x2, x5, p);
    p3 = sp_ntt_add_simd0(x3, x6, p);
    p4 = sp_ntt_sub_simd0(x1, x4, p);
    p5 = sp_ntt_sub_simd0(x2, x5, p);
    p6 = sp_ntt_sub_simd0(x3, x6, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t1 = sp_ntt_add_simd0(t1, p3, p);
    t2 = sp_ntt_sub_simd0(p1, p3, p);
    t3 = sp_ntt_sub_simd0(p2, p3, p);
    t4 = sp_ntt_sub_simd0(p4, p5, p);
    t4 = sp_ntt_add_partial_simd0(t4, p6, p);
    t5 = sp_ntt_sub_simd0(p4, p6, p);
    t6 = sp_ntt_add_simd0(p5, p6, p);

    p0 = sp_ntt_add_partial_simd0(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd0(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd0(t5, t6, p);

    p0 = sp_ntt_mul_simd0(p0, c+2*45, p);
    p1 = sp_ntt_mul_simd0(p1, c+2*46, p);
    p2 = sp_ntt_mul_simd0(p2, c+2*47, p);
    p3 = sp_ntt_mul_simd0(p3, c+2*48, p);
    p4 = sp_ntt_mul_simd0(p4, c+2*49, p);
    p5 = sp_ntt_mul_simd0(p5, c+2*50, p);
    p6 = sp_ntt_mul_simd0(p6, c+2*51, p);
    p7 = sp_ntt_mul_simd0(p7, c+2*52, p);
    p8 = sp_ntt_mul_simd0(p8, c+2*53, p);

    t1 = sp_ntt_add_simd0(p0, p1, p);
    t2 = sp_ntt_add_simd0(p2, p4, p);
    t3 = sp_ntt_add_simd0(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd0(p6, p8, p);
    t6 = sp_ntt_add_simd0(p7, p8, p);

    p1 = sp_ntt_add_simd0(t1, t2, p);
    p2 = sp_ntt_add_simd0(t1, t3, p);
    p3 = sp_ntt_sub_simd0(t1, t2, p);
    p3 = sp_ntt_sub_simd0(p3, t3, p);
    p4 = sp_ntt_add_simd0(t4, t5, p);
    p5 = sp_ntt_sub_simd0(t6, t4, p);
    p6 = sp_ntt_sub_simd0(t4, t5, p);
    p6 = sp_ntt_add_simd0(p6, t6, p);

    t1 = sp_ntt_add_simd0(p1, p4, p);
    t2 = sp_ntt_add_simd0(p2, p5, p);
    t3 = sp_ntt_add_simd0(p3, p6, p);
    t4 = sp_ntt_sub_simd0(p1, p4, p);
    t5 = sp_ntt_sub_simd0(p2, p5, p);
    t6 = sp_ntt_sub_simd0(p3, p6, p);

    a05 = p0;
    a11 = t6;
    a17 = t4;
    a23 = t5;
    a29 = t2;
    a35 = t1;
    a41 = t3;
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a00;
    sp_simd_t p1 = a01;
    sp_simd_t p2 = a02;
    sp_simd_t p3 = a03;
    sp_simd_t p4 = a04;
    sp_simd_t p5 = a05;

    p1 = sp_ntt_add_simd0(p0, p1, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t2 = sp_ntt_sub_simd0(p1, p2, p);
    t3 = sp_ntt_add_simd0(p3, p5, p);
    t4 = sp_ntt_add_simd0(p4, p5, p);

    p1 = sp_ntt_add_simd0(t1, t3, p);
    p2 = sp_ntt_add_simd0(t2, t4, p);
    p3 = sp_ntt_sub_simd0(t1, t3, p);
    p4 = sp_ntt_sub_simd0(t2, t4, p);

    sp_simd_store(p0, out +  0 * ostride);
    sp_simd_store(p3, out +  7 * ostride);
    sp_simd_store(p2, out + 14 * ostride);
    sp_simd_store(p4, out + 21 * ostride);
    sp_simd_store(p1, out + 28 * ostride);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a06;
    sp_simd_t p1 = a07;
    sp_simd_t p2 = a08;
    sp_simd_t p3 = a09;
    sp_simd_t p4 = a10;
    sp_simd_t p5 = a11;

    p1 = sp_ntt_add_simd0(p0, p1, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t2 = sp_ntt_sub_simd0(p1, p2, p);
    t3 = sp_ntt_add_simd0(p3, p5, p);
    t4 = sp_ntt_add_simd0(p4, p5, p);

    p1 = sp_ntt_add_simd0(t1, t3, p);
    p2 = sp_ntt_add_simd0(t2, t4, p);
    p3 = sp_ntt_sub_simd0(t1, t3, p);
    p4 = sp_ntt_sub_simd0(t2, t4, p);

    sp_simd_store(p4, out +  1 * ostride);
    sp_simd_store(p1, out +  8 * ostride);
    sp_simd_store(p0, out + 15 * ostride);
    sp_simd_store(p3, out + 22 * ostride);
    sp_simd_store(p2, out + 29 * ostride);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a12;
    sp_simd_t p1 = a13;
    sp_simd_t p2 = a14;
    sp_simd_t p3 = a15;
    sp_simd_t p4 = a16;
    sp_simd_t p5 = a17;

    p1 = sp_ntt_add_simd0(p0, p1, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t2 = sp_ntt_sub_simd0(p1, p2, p);
    t3 = sp_ntt_add_simd0(p3, p5, p);
    t4 = sp_ntt_add_simd0(p4, p5, p);

    p1 = sp_ntt_add_simd0(t1, t3, p);
    p2 = sp_ntt_add_simd0(t2, t4, p);
    p3 = sp_ntt_sub_simd0(t1, t3, p);
    p4 = sp_ntt_sub_simd0(t2, t4, p);

    sp_simd_store(p3, out +  2 * ostride);
    sp_simd_store(p2, out +  9 * ostride);
    sp_simd_store(p4, out + 16 * ostride);
    sp_simd_store(p1, out + 23 * ostride);
    sp_simd_store(p0, out + 30 * ostride);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a18;
    sp_simd_t p1 = a19;
    sp_simd_t p2 = a20;
    sp_simd_t p3 = a21;
    sp_simd_t p4 = a22;
    sp_simd_t p5 = a23;

    p1 = sp_ntt_add_simd0(p0, p1, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t2 = sp_ntt_sub_simd0(p1, p2, p);
    t3 = sp_ntt_add_simd0(p3, p5, p);
    t4 = sp_ntt_add_simd0(p4, p5, p);

    p1 = sp_ntt_add_simd0(t1, t3, p);
    p2 = sp_ntt_add_simd0(t2, t4, p);
    p3 = sp_ntt_sub_simd0(t1, t3, p);
    p4 = sp_ntt_sub_simd0(t2, t4, p);

    sp_simd_store(p1, out +  3 * ostride);
    sp_simd_store(p0, out + 10 * ostride);
    sp_simd_store(p3, out + 17 * ostride);
    sp_simd_store(p2, out + 24 * ostride);
    sp_simd_store(p4, out + 31 * ostride);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a24;
    sp_simd_t p1 = a25;
    sp_simd_t p2 = a26;
    sp_simd_t p3 = a27;
    sp_simd_t p4 = a28;
    sp_simd_t p5 = a29;

    p1 = sp_ntt_add_simd0(p0, p1, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t2 = sp_ntt_sub_simd0(p1, p2, p);
    t3 = sp_ntt_add_simd0(p3, p5, p);
    t4 = sp_ntt_add_simd0(p4, p5, p);

    p1 = sp_ntt_add_simd0(t1, t3, p);
    p2 = sp_ntt_add_simd0(t2, t4, p);
    p3 = sp_ntt_sub_simd0(t1, t3, p);
    p4 = sp_ntt_sub_simd0(t2, t4, p);

    sp_simd_store(p2, out +  4 * ostride);
    sp_simd_store(p4, out + 11 * ostride);
    sp_simd_store(p1, out + 18 * ostride);
    sp_simd_store(p0, out + 25 * ostride);
    sp_simd_store(p3, out + 32 * ostride);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a30;
    sp_simd_t p1 = a31;
    sp_simd_t p2 = a32;
    sp_simd_t p3 = a33;
    sp_simd_t p4 = a34;
    sp_simd_t p5 = a35;

    p1 = sp_ntt_add_simd0(p0, p1, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t2 = sp_ntt_sub_simd0(p1, p2, p);
    t3 = sp_ntt_add_simd0(p3, p5, p);
    t4 = sp_ntt_add_simd0(p4, p5, p);

    p1 = sp_ntt_add_simd0(t1, t3, p);
    p2 = sp_ntt_add_simd0(t2, t4, p);
    p3 = sp_ntt_sub_simd0(t1, t3, p);
    p4 = sp_ntt_sub_simd0(t2, t4, p);

    sp_simd_store(p0, out +  5 * ostride);
    sp_simd_store(p3, out + 12 * ostride);
    sp_simd_store(p2, out + 19 * ostride);
    sp_simd_store(p4, out + 26 * ostride);
    sp_simd_store(p1, out + 33 * ostride);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a36;
    sp_simd_t p1 = a37;
    sp_simd_t p2 = a38;
    sp_simd_t p3 = a39;
    sp_simd_t p4 = a40;
    sp_simd_t p5 = a41;

    p1 = sp_ntt_add_simd0(p0, p1, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t2 = sp_ntt_sub_simd0(p1, p2, p);
    t3 = sp_ntt_add_simd0(p3, p5, p);
    t4 = sp_ntt_add_simd0(p4, p5, p);

    p1 = sp_ntt_add_simd0(t1, t3, p);
    p2 = sp_ntt_add_simd0(t2, t4, p);
    p3 = sp_ntt_sub_simd0(t1, t3, p);
    p4 = sp_ntt_sub_simd0(t2, t4, p);

    sp_simd_store(p4, out +  6 * ostride);
    sp_simd_store(p1, out + 13 * ostride);
    sp_simd_store(p0, out + 20 * ostride);
    sp_simd_store(p3, out + 27 * ostride);
    sp_simd_store(p2, out + 34 * ostride);
  }
}

static void
ntt35_twiddle_run_core_simd(
        spv_t in, spv_size_t istride, spv_size_t idist,
		spv_t out, spv_size_t ostride, spv_size_t odist,
		sp_simd_t *w, sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  sp_simd_t a00, a01, a02, a03, a04, a05, a06, 
       a07, a08, a09, a10, a11, a12, a13, 
       a14, a15, a16, a17, a18, a19, a20, 
       a21, a22, a23, a24, a25, a26, a27, 
       a28, a29, a30, a31, a32, a33, a34, 
       a35, a36, a37, a38, a39, a40, a41;

  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x0 = sp_simd_gather(in +  0 * istride, idist, vsize);
    sp_simd_t x1 = sp_simd_gather(in +  7 * istride, idist, vsize);
    sp_simd_t x4 = sp_simd_gather(in + 14 * istride, idist, vsize);
    sp_simd_t x2 = sp_simd_gather(in + 21 * istride, idist, vsize);
    sp_simd_t x3 = sp_simd_gather(in + 28 * istride, idist, vsize);

    t1 = sp_ntt_add_simd(x1, x3, p);
    t3 = sp_ntt_sub_simd(x1, x3, p);
    t2 = sp_ntt_add_simd(x2, x4, p);
    t4 = sp_ntt_sub_simd(x2, x4, p);

    a01 = sp_ntt_add_simd(t1, t2, p);
    a02 = sp_ntt_sub_simd(t1, t2, p);
    a03 = t3;
    a04 = t4;
    a05 = sp_ntt_add_simd(t3, t4, p);

    a00 = sp_ntt_add_simd(x0, a01, p);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x0 = sp_simd_gather(in +  5 * istride, idist, vsize);
    sp_simd_t x1 = sp_simd_gather(in + 12 * istride, idist, vsize);
    sp_simd_t x4 = sp_simd_gather(in + 19 * istride, idist, vsize);
    sp_simd_t x2 = sp_simd_gather(in + 26 * istride, idist, vsize);
    sp_simd_t x3 = sp_simd_gather(in + 33 * istride, idist, vsize);

    t1 = sp_ntt_add_simd(x1, x3, p);
    t3 = sp_ntt_sub_simd(x1, x3, p);
    t2 = sp_ntt_add_simd(x2, x4, p);
    t4 = sp_ntt_sub_simd(x2, x4, p);

    a07 = sp_ntt_add_simd(t1, t2, p);
    a08 = sp_ntt_sub_simd(t1, t2, p);
    a09 = t3;
    a10 = t4;
    a11 = sp_ntt_add_simd(t3, t4, p);

    a06 = sp_ntt_add_simd(x0, a07, p);
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x3 = sp_simd_gather(in +  3 * istride, idist, vsize);
    sp_simd_t x0 = sp_simd_gather(in + 10 * istride, idist, vsize);
    sp_simd_t x1 = sp_simd_gather(in + 17 * istride, idist, vsize);
    sp_simd_t x4 = sp_simd_gather(in + 24 * istride, idist, vsize);
    sp_simd_t x2 = sp_simd_gather(in + 31 * istride, idist, vsize);

    t1 = sp_ntt_add_simd(x1, x3, p);
    t3 = sp_ntt_sub_simd(x1, x3, p);
    t2 = sp_ntt_add_simd(x2, x4, p);
    t4 = sp_ntt_sub_simd(x2, x4, p);

    a13 = sp_ntt_add_simd(t1, t2, p);
    a14 = sp_ntt_sub_simd(t1, t2, p);
    a15 = t3;
    a16 = t4;
    a17 = sp_ntt_add_simd(t3, t4, p);

    a12 = sp_ntt_add_simd(x0, a13, p);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x2 = sp_simd_gather(in +  1 * istride, idist, vsize);
    sp_simd_t x3 = sp_simd_gather(in +  8 * istride, idist, vsize);
    sp_simd_t x0 = sp_simd_gather(in + 15 * istride, idist, vsize);
    sp_simd_t x1 = sp_simd_gather(in + 22 * istride, idist, vsize);
    sp_simd_t x4 = sp_simd_gather(in + 29 * istride, idist, vsize);

    t1 = sp_ntt_add_simd(x1, x3, p);
    t3 = sp_ntt_sub_simd(x1, x3, p);
    t2 = sp_ntt_add_simd(x2, x4, p);
    t4 = sp_ntt_sub_simd(x2, x4, p);

    a19 = sp_ntt_add_simd(t1, t2, p);
    a20 = sp_ntt_sub_simd(t1, t2, p);
    a21 = t3;
    a22 = t4;
    a23 = sp_ntt_add_simd(t3, t4, p);

    a18 = sp_ntt_add_simd(x0, a19, p);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x2 = sp_simd_gather(in +  6 * istride, idist, vsize);
    sp_simd_t x3 = sp_simd_gather(in + 13 * istride, idist, vsize);
    sp_simd_t x0 = sp_simd_gather(in + 20 * istride, idist, vsize);
    sp_simd_t x1 = sp_simd_gather(in + 27 * istride, idist, vsize);
    sp_simd_t x4 = sp_simd_gather(in + 34 * istride, idist, vsize);

    t1 = sp_ntt_add_simd(x1, x3, p);
    t3 = sp_ntt_sub_simd(x1, x3, p);
    t2 = sp_ntt_add_simd(x2, x4, p);
    t4 = sp_ntt_sub_simd(x2, x4, p);

    a25 = sp_ntt_add_simd(t1, t2, p);
    a26 = sp_ntt_sub_simd(t1, t2, p);
    a27 = t3;
    a28 = t4;
    a29 = sp_ntt_add_simd(t3, t4, p);

    a24 = sp_ntt_add_simd(x0, a25, p);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x4 = sp_simd_gather(in +  4 * istride, idist, vsize);
    sp_simd_t x2 = sp_simd_gather(in + 11 * istride, idist, vsize);
    sp_simd_t x3 = sp_simd_gather(in + 18 * istride, idist, vsize);
    sp_simd_t x0 = sp_simd_gather(in + 25 * istride, idist, vsize);
    sp_simd_t x1 = sp_simd_gather(in + 32 * istride, idist, vsize);

    t1 = sp_ntt_add_simd(x1, x3, p);
    t3 = sp_ntt_sub_simd(x1, x3, p);
    t2 = sp_ntt_add_simd(x2, x4, p);
    t4 = sp_ntt_sub_simd(x2, x4, p);

    a31 = sp_ntt_add_simd(t1, t2, p);
    a32 = sp_ntt_sub_simd(t1, t2, p);
    a33 = t3;
    a34 = t4;
    a35 = sp_ntt_add_simd(t3, t4, p);

    a30 = sp_ntt_add_simd(x0, a31, p);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x1 = sp_simd_gather(in +  2 * istride, idist, vsize);
    sp_simd_t x4 = sp_simd_gather(in +  9 * istride, idist, vsize);
    sp_simd_t x2 = sp_simd_gather(in + 16 * istride, idist, vsize);
    sp_simd_t x3 = sp_simd_gather(in + 23 * istride, idist, vsize);
    sp_simd_t x0 = sp_simd_gather(in + 30 * istride, idist, vsize);

    t1 = sp_ntt_add_simd(x1, x3, p);
    t3 = sp_ntt_sub_simd(x1, x3, p);
    t2 = sp_ntt_add_simd(x2, x4, p);
    t4 = sp_ntt_sub_simd(x2, x4, p);

    a37 = sp_ntt_add_simd(t1, t2, p);
    a38 = sp_ntt_sub_simd(t1, t2, p);
    a39 = t3;
    a40 = t4;
    a41 = sp_ntt_add_simd(t3, t4, p);

    a36 = sp_ntt_add_simd(x0, a37, p);
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a00;
    sp_simd_t x1 = a06;
    sp_simd_t x5 = a12;
    sp_simd_t x6 = a18;
    sp_simd_t x3 = a24;
    sp_simd_t x2 = a30;
    sp_simd_t x4 = a36;

    p1 = sp_ntt_add_simd(x1, x4, p);
    p2 = sp_ntt_add_simd(x2, x5, p);
    p3 = sp_ntt_add_simd(x3, x6, p);
    p4 = sp_ntt_sub_simd(x1, x4, p);
    p5 = sp_ntt_sub_simd(x2, x5, p);
    p6 = sp_ntt_sub_simd(x3, x6, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t1 = sp_ntt_add_simd(t1, p3, p);
    t2 = sp_ntt_sub_simd(p1, p3, p);
    t3 = sp_ntt_sub_simd(p2, p3, p);
    t4 = sp_ntt_sub_simd(p4, p5, p);
    t4 = sp_ntt_add_partial_simd(t4, p6, p);
    t5 = sp_ntt_sub_simd(p4, p6, p);
    t6 = sp_ntt_add_simd(p5, p6, p);

    p0 = sp_ntt_add_simd(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd(t5, t6, p);

    p1 = sp_ntt_mul_simd(p1, ntt_const[1], ntt_const[NC+1], p);
    p2 = sp_ntt_mul_simd(p2, ntt_const[2], ntt_const[NC+2], p);
    p3 = sp_ntt_mul_simd(p3, ntt_const[3], ntt_const[NC+3], p);
    p4 = sp_ntt_mul_simd(p4, ntt_const[4], ntt_const[NC+4], p);
    p5 = sp_ntt_mul_simd(p5, ntt_const[5], ntt_const[NC+5], p);
    p6 = sp_ntt_mul_simd(p6, ntt_const[6], ntt_const[NC+6], p);
    p7 = sp_ntt_mul_simd(p7, ntt_const[7], ntt_const[NC+7], p);
    p8 = sp_ntt_mul_simd(p8, ntt_const[8], ntt_const[NC+8], p);

    t1 = sp_ntt_add_simd(p0, p1, p);
    t2 = sp_ntt_add_simd(p2, p4, p);
    t3 = sp_ntt_add_simd(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd(p6, p8, p);
    t6 = sp_ntt_add_simd(p7, p8, p);

    p1 = sp_ntt_add_simd(t1, t2, p);
    p2 = sp_ntt_add_simd(t1, t3, p);
    p3 = sp_ntt_sub_simd(t1, t2, p);
    p3 = sp_ntt_sub_simd(p3, t3, p);
    p4 = sp_ntt_add_simd(t4, t5, p);
    p5 = sp_ntt_sub_simd(t6, t4, p);
    p6 = sp_ntt_sub_simd(t4, t5, p);
    p6 = sp_ntt_add_simd(p6, t6, p);

    t1 = sp_ntt_add_simd(p1, p4, p);
    t2 = sp_ntt_add_simd(p2, p5, p);
    t3 = sp_ntt_add_simd(p3, p6, p);
    t4 = sp_ntt_sub_simd(p1, p4, p);
    t5 = sp_ntt_sub_simd(p2, p5, p);
    t6 = sp_ntt_sub_simd(p3, p6, p);

    a00 = p0;
    a06 = t6;
    a12 = t4;
    a18 = t5;
    a24 = t2;
    a30 = t1;
    a36 = t3;
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a01;
    sp_simd_t x1 = a07;
    sp_simd_t x5 = a13;
    sp_simd_t x6 = a19;
    sp_simd_t x3 = a25;
    sp_simd_t x2 = a31;
    sp_simd_t x4 = a37;

    p1 = sp_ntt_add_simd(x1, x4, p);
    p2 = sp_ntt_add_simd(x2, x5, p);
    p3 = sp_ntt_add_simd(x3, x6, p);
    p4 = sp_ntt_sub_simd(x1, x4, p);
    p5 = sp_ntt_sub_simd(x2, x5, p);
    p6 = sp_ntt_sub_simd(x3, x6, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t1 = sp_ntt_add_simd(t1, p3, p);
    t2 = sp_ntt_sub_simd(p1, p3, p);
    t3 = sp_ntt_sub_simd(p2, p3, p);
    t4 = sp_ntt_sub_simd(p4, p5, p);
    t4 = sp_ntt_add_partial_simd(t4, p6, p);
    t5 = sp_ntt_sub_simd(p4, p6, p);
    t6 = sp_ntt_add_simd(p5, p6, p);

    p0 = sp_ntt_add_partial_simd(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd(t5, t6, p);

    p0 = sp_ntt_mul_simd(p0, ntt_const[9], ntt_const[NC+9], p);
    p1 = sp_ntt_mul_simd(p1, ntt_const[10], ntt_const[NC+10], p);
    p2 = sp_ntt_mul_simd(p2, ntt_const[11], ntt_const[NC+11], p);
    p3 = sp_ntt_mul_simd(p3, ntt_const[12], ntt_const[NC+12], p);
    p4 = sp_ntt_mul_simd(p4, ntt_const[13], ntt_const[NC+13], p);
    p5 = sp_ntt_mul_simd(p5, ntt_const[14], ntt_const[NC+14], p);
    p6 = sp_ntt_mul_simd(p6, ntt_const[15], ntt_const[NC+15], p);
    p7 = sp_ntt_mul_simd(p7, ntt_const[16], ntt_const[NC+16], p);
    p8 = sp_ntt_mul_simd(p8, ntt_const[17], ntt_const[NC+17], p);

    t1 = sp_ntt_add_simd(p0, p1, p);
    t2 = sp_ntt_add_simd(p2, p4, p);
    t3 = sp_ntt_add_simd(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd(p6, p8, p);
    t6 = sp_ntt_add_simd(p7, p8, p);

    p1 = sp_ntt_add_simd(t1, t2, p);
    p2 = sp_ntt_add_simd(t1, t3, p);
    p3 = sp_ntt_sub_simd(t1, t2, p);
    p3 = sp_ntt_sub_simd(p3, t3, p);
    p4 = sp_ntt_add_simd(t4, t5, p);
    p5 = sp_ntt_sub_simd(t6, t4, p);
    p6 = sp_ntt_sub_simd(t4, t5, p);
    p6 = sp_ntt_add_simd(p6, t6, p);

    t1 = sp_ntt_add_simd(p1, p4, p);
    t2 = sp_ntt_add_simd(p2, p5, p);
    t3 = sp_ntt_add_simd(p3, p6, p);
    t4 = sp_ntt_sub_simd(p1, p4, p);
    t5 = sp_ntt_sub_simd(p2, p5, p);
    t6 = sp_ntt_sub_simd(p3, p6, p);

    a01 = p0;
    a07 = t6;
    a13 = t4;
    a19 = t5;
    a25 = t2;
    a31 = t1;
    a37 = t3;
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a02;
    sp_simd_t x1 = a08;
    sp_simd_t x5 = a14;
    sp_simd_t x6 = a20;
    sp_simd_t x3 = a26;
    sp_simd_t x2 = a32;
    sp_simd_t x4 = a38;

    p1 = sp_ntt_add_simd(x1, x4, p);
    p2 = sp_ntt_add_simd(x2, x5, p);
    p3 = sp_ntt_add_simd(x3, x6, p);
    p4 = sp_ntt_sub_simd(x1, x4, p);
    p5 = sp_ntt_sub_simd(x2, x5, p);
    p6 = sp_ntt_sub_simd(x3, x6, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t1 = sp_ntt_add_simd(t1, p3, p);
    t2 = sp_ntt_sub_simd(p1, p3, p);
    t3 = sp_ntt_sub_simd(p2, p3, p);
    t4 = sp_ntt_sub_simd(p4, p5, p);
    t4 = sp_ntt_add_partial_simd(t4, p6, p);
    t5 = sp_ntt_sub_simd(p4, p6, p);
    t6 = sp_ntt_add_simd(p5, p6, p);

    p0 = sp_ntt_add_partial_simd(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd(t5, t6, p);

    p0 = sp_ntt_mul_simd(p0, ntt_const[18], ntt_const[NC+18], p);
    p1 = sp_ntt_mul_simd(p1, ntt_const[19], ntt_const[NC+19], p);
    p2 = sp_ntt_mul_simd(p2, ntt_const[20], ntt_const[NC+20], p);
    p3 = sp_ntt_mul_simd(p3, ntt_const[21], ntt_const[NC+21], p);
    p4 = sp_ntt_mul_simd(p4, ntt_const[22], ntt_const[NC+22], p);
    p5 = sp_ntt_mul_simd(p5, ntt_const[23], ntt_const[NC+23], p);
    p6 = sp_ntt_mul_simd(p6, ntt_const[24], ntt_const[NC+24], p);
    p7 = sp_ntt_mul_simd(p7, ntt_const[25], ntt_const[NC+25], p);
    p8 = sp_ntt_mul_simd(p8, ntt_const[26], ntt_const[NC+26], p);

    t1 = sp_ntt_add_simd(p0, p1, p);
    t2 = sp_ntt_add_simd(p2, p4, p);
    t3 = sp_ntt_add_simd(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd(p6, p8, p);
    t6 = sp_ntt_add_simd(p7, p8, p);

    p1 = sp_ntt_add_simd(t1, t2, p);
    p2 = sp_ntt_add_simd(t1, t3, p);
    p3 = sp_ntt_sub_simd(t1, t2, p);
    p3 = sp_ntt_sub_simd(p3, t3, p);
    p4 = sp_ntt_add_simd(t4, t5, p);
    p5 = sp_ntt_sub_simd(t6, t4, p);
    p6 = sp_ntt_sub_simd(t4, t5, p);
    p6 = sp_ntt_add_simd(p6, t6, p);

    t1 = sp_ntt_add_simd(p1, p4, p);
    t2 = sp_ntt_add_simd(p2, p5, p);
    t3 = sp_ntt_add_simd(p3, p6, p);
    t4 = sp_ntt_sub_simd(p1, p4, p);
    t5 = sp_ntt_sub_simd(p2, p5, p);
    t6 = sp_ntt_sub_simd(p3, p6, p);

    a02 = p0;
    a08 = t6;
    a14 = t4;
    a20 = t5;
    a26 = t2;
    a32 = t1;
    a38 = t3;
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a03;
    sp_simd_t x1 = a09;
    sp_simd_t x5 = a15;
    sp_simd_t x6 = a21;
    sp_simd_t x3 = a27;
    sp_simd_t x2 = a33;
    sp_simd_t x4 = a39;

    p1 = sp_ntt_add_simd(x1, x4, p);
    p2 = sp_ntt_add_simd(x2, x5, p);
    p3 = sp_ntt_add_simd(x3, x6, p);
    p4 = sp_ntt_sub_simd(x1, x4, p);
    p5 = sp_ntt_sub_simd(x2, x5, p);
    p6 = sp_ntt_sub_simd(x3, x6, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t1 = sp_ntt_add_simd(t1, p3, p);
    t2 = sp_ntt_sub_simd(p1, p3, p);
    t3 = sp_ntt_sub_simd(p2, p3, p);
    t4 = sp_ntt_sub_simd(p4, p5, p);
    t4 = sp_ntt_add_partial_simd(t4, p6, p);
    t5 = sp_ntt_sub_simd(p4, p6, p);
    t6 = sp_ntt_add_simd(p5, p6, p);

    p0 = sp_ntt_add_partial_simd(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd(t5, t6, p);

    p0 = sp_ntt_mul_simd(p0, ntt_const[27], ntt_const[NC+27], p);
    p1 = sp_ntt_mul_simd(p1, ntt_const[28], ntt_const[NC+28], p);
    p2 = sp_ntt_mul_simd(p2, ntt_const[29], ntt_const[NC+29], p);
    p3 = sp_ntt_mul_simd(p3, ntt_const[30], ntt_const[NC+30], p);
    p4 = sp_ntt_mul_simd(p4, ntt_const[31], ntt_const[NC+31], p);
    p5 = sp_ntt_mul_simd(p5, ntt_const[32], ntt_const[NC+32], p);
    p6 = sp_ntt_mul_simd(p6, ntt_const[33], ntt_const[NC+33], p);
    p7 = sp_ntt_mul_simd(p7, ntt_const[34], ntt_const[NC+34], p);
    p8 = sp_ntt_mul_simd(p8, ntt_const[35], ntt_const[NC+35], p);

    t1 = sp_ntt_add_simd(p0, p1, p);
    t2 = sp_ntt_add_simd(p2, p4, p);
    t3 = sp_ntt_add_simd(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd(p6, p8, p);
    t6 = sp_ntt_add_simd(p7, p8, p);

    p1 = sp_ntt_add_simd(t1, t2, p);
    p2 = sp_ntt_add_simd(t1, t3, p);
    p3 = sp_ntt_sub_simd(t1, t2, p);
    p3 = sp_ntt_sub_simd(p3, t3, p);
    p4 = sp_ntt_add_simd(t4, t5, p);
    p5 = sp_ntt_sub_simd(t6, t4, p);
    p6 = sp_ntt_sub_simd(t4, t5, p);
    p6 = sp_ntt_add_simd(p6, t6, p);

    t1 = sp_ntt_add_simd(p1, p4, p);
    t2 = sp_ntt_add_simd(p2, p5, p);
    t3 = sp_ntt_add_simd(p3, p6, p);
    t4 = sp_ntt_sub_simd(p1, p4, p);
    t5 = sp_ntt_sub_simd(p2, p5, p);
    t6 = sp_ntt_sub_simd(p3, p6, p);

    a03 = p0;
    a09 = t6;
    a15 = t4;
    a21 = t5;
    a27 = t2;
    a33 = t1;
    a39 = t3;
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a04;
    sp_simd_t x1 = a10;
    sp_simd_t x5 = a16;
    sp_simd_t x6 = a22;
    sp_simd_t x3 = a28;
    sp_simd_t x2 = a34;
    sp_simd_t x4 = a40;

    p1 = sp_ntt_add_simd(x1, x4, p);
    p2 = sp_ntt_add_simd(x2, x5, p);
    p3 = sp_ntt_add_simd(x3, x6, p);
    p4 = sp_ntt_sub_simd(x1, x4, p);
    p5 = sp_ntt_sub_simd(x2, x5, p);
    p6 = sp_ntt_sub_simd(x3, x6, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t1 = sp_ntt_add_simd(t1, p3, p);
    t2 = sp_ntt_sub_simd(p1, p3, p);
    t3 = sp_ntt_sub_simd(p2, p3, p);
    t4 = sp_ntt_sub_simd(p4, p5, p);
    t4 = sp_ntt_add_partial_simd(t4, p6, p);
    t5 = sp_ntt_sub_simd(p4, p6, p);
    t6 = sp_ntt_add_simd(p5, p6, p);

    p0 = sp_ntt_add_partial_simd(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd(t5, t6, p);

    p0 = sp_ntt_mul_simd(p0, ntt_const[36], ntt_const[NC+36], p);
    p1 = sp_ntt_mul_simd(p1, ntt_const[37], ntt_const[NC+37], p);
    p2 = sp_ntt_mul_simd(p2, ntt_const[38], ntt_const[NC+38], p);
    p3 = sp_ntt_mul_simd(p3, ntt_const[39], ntt_const[NC+39], p);
    p4 = sp_ntt_mul_simd(p4, ntt_const[40], ntt_const[NC+40], p);
    p5 = sp_ntt_mul_simd(p5, ntt_const[41], ntt_const[NC+41], p);
    p6 = sp_ntt_mul_simd(p6, ntt_const[42], ntt_const[NC+42], p);
    p7 = sp_ntt_mul_simd(p7, ntt_const[43], ntt_const[NC+43], p);
    p8 = sp_ntt_mul_simd(p8, ntt_const[44], ntt_const[NC+44], p);

    t1 = sp_ntt_add_simd(p0, p1, p);
    t2 = sp_ntt_add_simd(p2, p4, p);
    t3 = sp_ntt_add_simd(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd(p6, p8, p);
    t6 = sp_ntt_add_simd(p7, p8, p);

    p1 = sp_ntt_add_simd(t1, t2, p);
    p2 = sp_ntt_add_simd(t1, t3, p);
    p3 = sp_ntt_sub_simd(t1, t2, p);
    p3 = sp_ntt_sub_simd(p3, t3, p);
    p4 = sp_ntt_add_simd(t4, t5, p);
    p5 = sp_ntt_sub_simd(t6, t4, p);
    p6 = sp_ntt_sub_simd(t4, t5, p);
    p6 = sp_ntt_add_simd(p6, t6, p);

    t1 = sp_ntt_add_simd(p1, p4, p);
    t2 = sp_ntt_add_simd(p2, p5, p);
    t3 = sp_ntt_add_simd(p3, p6, p);
    t4 = sp_ntt_sub_simd(p1, p4, p);
    t5 = sp_ntt_sub_simd(p2, p5, p);
    t6 = sp_ntt_sub_simd(p3, p6, p);

    a04 = p0;
    a10 = t6;
    a16 = t4;
    a22 = t5;
    a28 = t2;
    a34 = t1;
    a40 = t3;
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a05;
    sp_simd_t x1 = a11;
    sp_simd_t x5 = a17;
    sp_simd_t x6 = a23;
    sp_simd_t x3 = a29;
    sp_simd_t x2 = a35;
    sp_simd_t x4 = a41;

    p1 = sp_ntt_add_simd(x1, x4, p);
    p2 = sp_ntt_add_simd(x2, x5, p);
    p3 = sp_ntt_add_simd(x3, x6, p);
    p4 = sp_ntt_sub_simd(x1, x4, p);
    p5 = sp_ntt_sub_simd(x2, x5, p);
    p6 = sp_ntt_sub_simd(x3, x6, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t1 = sp_ntt_add_simd(t1, p3, p);
    t2 = sp_ntt_sub_simd(p1, p3, p);
    t3 = sp_ntt_sub_simd(p2, p3, p);
    t4 = sp_ntt_sub_simd(p4, p5, p);
    t4 = sp_ntt_add_partial_simd(t4, p6, p);
    t5 = sp_ntt_sub_simd(p4, p6, p);
    t6 = sp_ntt_add_simd(p5, p6, p);

    p0 = sp_ntt_add_partial_simd(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd(t5, t6, p);

    p0 = sp_ntt_mul_simd(p0, ntt_const[45], ntt_const[NC+45], p);
    p1 = sp_ntt_mul_simd(p1, ntt_const[46], ntt_const[NC+46], p);
    p2 = sp_ntt_mul_simd(p2, ntt_const[47], ntt_const[NC+47], p);
    p3 = sp_ntt_mul_simd(p3, ntt_const[48], ntt_const[NC+48], p);
    p4 = sp_ntt_mul_simd(p4, ntt_const[49], ntt_const[NC+49], p);
    p5 = sp_ntt_mul_simd(p5, ntt_const[50], ntt_const[NC+50], p);
    p6 = sp_ntt_mul_simd(p6, ntt_const[51], ntt_const[NC+51], p);
    p7 = sp_ntt_mul_simd(p7, ntt_const[52], ntt_const[NC+52], p);
    p8 = sp_ntt_mul_simd(p8, ntt_const[53], ntt_const[NC+53], p);

    t1 = sp_ntt_add_simd(p0, p1, p);
    t2 = sp_ntt_add_simd(p2, p4, p);
    t3 = sp_ntt_add_simd(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd(p6, p8, p);
    t6 = sp_ntt_add_simd(p7, p8, p);

    p1 = sp_ntt_add_simd(t1, t2, p);
    p2 = sp_ntt_add_simd(t1, t3, p);
    p3 = sp_ntt_sub_simd(t1, t2, p);
    p3 = sp_ntt_sub_simd(p3, t3, p);
    p4 = sp_ntt_add_simd(t4, t5, p);
    p5 = sp_ntt_sub_simd(t6, t4, p);
    p6 = sp_ntt_sub_simd(t4, t5, p);
    p6 = sp_ntt_add_simd(p6, t6, p);

    t1 = sp_ntt_add_simd(p1, p4, p);
    t2 = sp_ntt_add_simd(p2, p5, p);
    t3 = sp_ntt_add_simd(p3, p6, p);
    t4 = sp_ntt_sub_simd(p1, p4, p);
    t5 = sp_ntt_sub_simd(p2, p5, p);
    t6 = sp_ntt_sub_simd(p3, p6, p);

    a05 = p0;
    a11 = t6;
    a17 = t4;
    a23 = t5;
    a29 = t2;
    a35 = t1;
    a41 = t3;
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a00;
    sp_simd_t p1 = a01;
    sp_simd_t p2 = a02;
    sp_simd_t p3 = a03;
    sp_simd_t p4 = a04;
    sp_simd_t p5 = a05;

    p1 = sp_ntt_add_simd(p0, p1, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t2 = sp_ntt_sub_simd(p1, p2, p);
    t3 = sp_ntt_add_simd(p3, p5, p);
    t4 = sp_ntt_add_simd(p4, p5, p);

    p1 = sp_ntt_add_partial_simd(t1, t3, p);
    p2 = sp_ntt_add_partial_simd(t2, t4, p);
    p3 = sp_ntt_sub_partial_simd(t1, t3, p);
    p4 = sp_ntt_sub_partial_simd(t2, t4, p);

    p3 = sp_ntt_twiddle_mul_simd(p3, w + 12, p);
    p2 = sp_ntt_twiddle_mul_simd(p2, w + 26, p);
    p4 = sp_ntt_twiddle_mul_simd(p4, w + 40, p);
    p1 = sp_ntt_twiddle_mul_simd(p1, w + 54, p);

    sp_simd_scatter(p0, out +  0 * ostride, odist, vsize);
    sp_simd_scatter(p3, out +  7 * ostride, odist, vsize);
    sp_simd_scatter(p2, out + 14 * ostride, odist, vsize);
    sp_simd_scatter(p4, out + 21 * ostride, odist, vsize);
    sp_simd_scatter(p1, out + 28 * ostride, odist, vsize);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a06;
    sp_simd_t p1 = a07;
    sp_simd_t p2 = a08;
    sp_simd_t p3 = a09;
    sp_simd_t p4 = a10;
    sp_simd_t p5 = a11;

    p1 = sp_ntt_add_simd(p0, p1, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t2 = sp_ntt_sub_simd(p1, p2, p);
    t3 = sp_ntt_add_simd(p3, p5, p);
    t4 = sp_ntt_add_simd(p4, p5, p);

    p1 = sp_ntt_add_partial_simd(t1, t3, p);
    p2 = sp_ntt_add_partial_simd(t2, t4, p);
    p3 = sp_ntt_sub_partial_simd(t1, t3, p);
    p4 = sp_ntt_sub_partial_simd(t2, t4, p);

    p4 = sp_ntt_twiddle_mul_simd(p4, w + 0, p);
    p1 = sp_ntt_twiddle_mul_simd(p1, w + 14, p);
    p0 = sp_ntt_twiddle_mul_simd(p0, w + 28, p);
    p3 = sp_ntt_twiddle_mul_simd(p3, w + 42, p);
    p2 = sp_ntt_twiddle_mul_simd(p2, w + 56, p);

    sp_simd_scatter(p4, out +  1 * ostride, odist, vsize);
    sp_simd_scatter(p1, out +  8 * ostride, odist, vsize);
    sp_simd_scatter(p0, out + 15 * ostride, odist, vsize);
    sp_simd_scatter(p3, out + 22 * ostride, odist, vsize);
    sp_simd_scatter(p2, out + 29 * ostride, odist, vsize);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a12;
    sp_simd_t p1 = a13;
    sp_simd_t p2 = a14;
    sp_simd_t p3 = a15;
    sp_simd_t p4 = a16;
    sp_simd_t p5 = a17;

    p1 = sp_ntt_add_simd(p0, p1, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t2 = sp_ntt_sub_simd(p1, p2, p);
    t3 = sp_ntt_add_simd(p3, p5, p);
    t4 = sp_ntt_add_simd(p4, p5, p);

    p1 = sp_ntt_add_partial_simd(t1, t3, p);
    p2 = sp_ntt_add_partial_simd(t2, t4, p);
    p3 = sp_ntt_sub_partial_simd(t1, t3, p);
    p4 = sp_ntt_sub_partial_simd(t2, t4, p);

    p3 = sp_ntt_twiddle_mul_simd(p3, w + 2, p);
    p2 = sp_ntt_twiddle_mul_simd(p2, w + 16, p);
    p4 = sp_ntt_twiddle_mul_simd(p4, w + 30, p);
    p1 = sp_ntt_twiddle_mul_simd(p1, w + 44, p);
    p0 = sp_ntt_twiddle_mul_simd(p0, w + 58, p);

    sp_simd_scatter(p3, out +  2 * ostride, odist, vsize);
    sp_simd_scatter(p2, out +  9 * ostride, odist, vsize);
    sp_simd_scatter(p4, out + 16 * ostride, odist, vsize);
    sp_simd_scatter(p1, out + 23 * ostride, odist, vsize);
    sp_simd_scatter(p0, out + 30 * ostride, odist, vsize);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a18;
    sp_simd_t p1 = a19;
    sp_simd_t p2 = a20;
    sp_simd_t p3 = a21;
    sp_simd_t p4 = a22;
    sp_simd_t p5 = a23;

    p1 = sp_ntt_add_simd(p0, p1, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t2 = sp_ntt_sub_simd(p1, p2, p);
    t3 = sp_ntt_add_simd(p3, p5, p);
    t4 = sp_ntt_add_simd(p4, p5, p);

    p1 = sp_ntt_add_partial_simd(t1, t3, p);
    p2 = sp_ntt_add_partial_simd(t2, t4, p);
    p3 = sp_ntt_sub_partial_simd(t1, t3, p);
    p4 = sp_ntt_sub_partial_simd(t2, t4, p);

    p1 = sp_ntt_twiddle_mul_simd(p1, w + 4, p);
    p0 = sp_ntt_twiddle_mul_simd(p0, w + 18, p);
    p3 = sp_ntt_twiddle_mul_simd(p3, w + 32, p);
    p2 = sp_ntt_twiddle_mul_simd(p2, w + 46, p);
    p4 = sp_ntt_twiddle_mul_simd(p4, w + 60, p);

    sp_simd_scatter(p1, out +  3 * ostride, odist, vsize);
    sp_simd_scatter(p0, out + 10 * ostride, odist, vsize);
    sp_simd_scatter(p3, out + 17 * ostride, odist, vsize);
    sp_simd_scatter(p2, out + 24 * ostride, odist, vsize);
    sp_simd_scatter(p4, out + 31 * ostride, odist, vsize);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a24;
    sp_simd_t p1 = a25;
    sp_simd_t p2 = a26;
    sp_simd_t p3 = a27;
    sp_simd_t p4 = a28;
    sp_simd_t p5 = a29;

    p1 = sp_ntt_add_simd(p0, p1, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t2 = sp_ntt_sub_simd(p1, p2, p);
    t3 = sp_ntt_add_simd(p3, p5, p);
    t4 = sp_ntt_add_simd(p4, p5, p);

    p1 = sp_ntt_add_partial_simd(t1, t3, p);
    p2 = sp_ntt_add_partial_simd(t2, t4, p);
    p3 = sp_ntt_sub_partial_simd(t1, t3, p);
    p4 = sp_ntt_sub_partial_simd(t2, t4, p);

    p2 = sp_ntt_twiddle_mul_simd(p2, w + 6, p);
    p4 = sp_ntt_twiddle_mul_simd(p4, w + 20, p);
    p1 = sp_ntt_twiddle_mul_simd(p1, w + 34, p);
    p0 = sp_ntt_twiddle_mul_simd(p0, w + 48, p);
    p3 = sp_ntt_twiddle_mul_simd(p3, w + 62, p);

    sp_simd_scatter(p2, out +  4 * ostride, odist, vsize);
    sp_simd_scatter(p4, out + 11 * ostride, odist, vsize);
    sp_simd_scatter(p1, out + 18 * ostride, odist, vsize);
    sp_simd_scatter(p0, out + 25 * ostride, odist, vsize);
    sp_simd_scatter(p3, out + 32 * ostride, odist, vsize);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a30;
    sp_simd_t p1 = a31;
    sp_simd_t p2 = a32;
    sp_simd_t p3 = a33;
    sp_simd_t p4 = a34;
    sp_simd_t p5 = a35;

    p1 = sp_ntt_add_simd(p0, p1, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t2 = sp_ntt_sub_simd(p1, p2, p);
    t3 = sp_ntt_add_simd(p3, p5, p);
    t4 = sp_ntt_add_simd(p4, p5, p);

    p1 = sp_ntt_add_partial_simd(t1, t3, p);
    p2 = sp_ntt_add_partial_simd(t2, t4, p);
    p3 = sp_ntt_sub_partial_simd(t1, t3, p);
    p4 = sp_ntt_sub_partial_simd(t2, t4, p);

    p0 = sp_ntt_twiddle_mul_simd(p0, w + 8, p);
    p3 = sp_ntt_twiddle_mul_simd(p3, w + 22, p);
    p2 = sp_ntt_twiddle_mul_simd(p2, w + 36, p);
    p4 = sp_ntt_twiddle_mul_simd(p4, w + 50, p);
    p1 = sp_ntt_twiddle_mul_simd(p1, w + 64, p);

    sp_simd_scatter(p0, out +  5 * ostride, odist, vsize);
    sp_simd_scatter(p3, out + 12 * ostride, odist, vsize);
    sp_simd_scatter(p2, out + 19 * ostride, odist, vsize);
    sp_simd_scatter(p4, out + 26 * ostride, odist, vsize);
    sp_simd_scatter(p1, out + 33 * ostride, odist, vsize);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a36;
    sp_simd_t p1 = a37;
    sp_simd_t p2 = a38;
    sp_simd_t p3 = a39;
    sp_simd_t p4 = a40;
    sp_simd_t p5 = a41;

    p1 = sp_ntt_add_simd(p0, p1, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t2 = sp_ntt_sub_simd(p1, p2, p);
    t3 = sp_ntt_add_simd(p3, p5, p);
    t4 = sp_ntt_add_simd(p4, p5, p);

    p1 = sp_ntt_add_partial_simd(t1, t3, p);
    p2 = sp_ntt_add_partial_simd(t2, t4, p);
    p3 = sp_ntt_sub_partial_simd(t1, t3, p);
    p4 = sp_ntt_sub_partial_simd(t2, t4, p);

    p4 = sp_ntt_twiddle_mul_simd(p4, w + 10, p);
    p1 = sp_ntt_twiddle_mul_simd(p1, w + 24, p);
    p0 = sp_ntt_twiddle_mul_simd(p0, w + 38, p);
    p3 = sp_ntt_twiddle_mul_simd(p3, w + 52, p);
    p2 = sp_ntt_twiddle_mul_simd(p2, w + 66, p);

    sp_simd_scatter(p4, out +  6 * ostride, odist, vsize);
    sp_simd_scatter(p1, out + 13 * ostride, odist, vsize);
    sp_simd_scatter(p0, out + 20 * ostride, odist, vsize);
    sp_simd_scatter(p3, out + 27 * ostride, odist, vsize);
    sp_simd_scatter(p2, out + 34 * ostride, odist, vsize);
  }
}

static void
ntt35_twiddle_run_core_simd_interleaved(
		spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		sp_simd_t *w, sp_simd_t p, sp_simd_t * c)
{
  sp_simd_t a00, a01, a02, a03, a04, a05, a06, 
       a07, a08, a09, a10, a11, a12, a13, 
       a14, a15, a16, a17, a18, a19, a20, 
       a21, a22, a23, a24, a25, a26, a27, 
       a28, a29, a30, a31, a32, a33, a34, 
       a35, a36, a37, a38, a39, a40, a41;

  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x0 = sp_simd_load(in +  0 * istride);
    sp_simd_t x1 = sp_simd_load(in +  7 * istride);
    sp_simd_t x4 = sp_simd_load(in + 14 * istride);
    sp_simd_t x2 = sp_simd_load(in + 21 * istride);
    sp_simd_t x3 = sp_simd_load(in + 28 * istride);

    t1 = sp_ntt_add_simd0(x1, x3, p);
    t3 = sp_ntt_sub_simd0(x1, x3, p);
    t2 = sp_ntt_add_simd0(x2, x4, p);
    t4 = sp_ntt_sub_simd0(x2, x4, p);

    a01 = sp_ntt_add_simd0(t1, t2, p);
    a02 = sp_ntt_sub_simd0(t1, t2, p);
    a03 = t3;
    a04 = t4;
    a05 = sp_ntt_add_simd0(t3, t4, p);

    a00 = sp_ntt_add_simd0(x0, a01, p);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x0 = sp_simd_load(in +  5 * istride);
    sp_simd_t x1 = sp_simd_load(in + 12 * istride);
    sp_simd_t x4 = sp_simd_load(in + 19 * istride);
    sp_simd_t x2 = sp_simd_load(in + 26 * istride);
    sp_simd_t x3 = sp_simd_load(in + 33 * istride);

    t1 = sp_ntt_add_simd0(x1, x3, p);
    t3 = sp_ntt_sub_simd0(x1, x3, p);
    t2 = sp_ntt_add_simd0(x2, x4, p);
    t4 = sp_ntt_sub_simd0(x2, x4, p);

    a07 = sp_ntt_add_simd0(t1, t2, p);
    a08 = sp_ntt_sub_simd0(t1, t2, p);
    a09 = t3;
    a10 = t4;
    a11 = sp_ntt_add_simd0(t3, t4, p);

    a06 = sp_ntt_add_simd0(x0, a07, p);
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x3 = sp_simd_load(in +  3 * istride);
    sp_simd_t x0 = sp_simd_load(in + 10 * istride);
    sp_simd_t x1 = sp_simd_load(in + 17 * istride);
    sp_simd_t x4 = sp_simd_load(in + 24 * istride);
    sp_simd_t x2 = sp_simd_load(in + 31 * istride);

    t1 = sp_ntt_add_simd0(x1, x3, p);
    t3 = sp_ntt_sub_simd0(x1, x3, p);
    t2 = sp_ntt_add_simd0(x2, x4, p);
    t4 = sp_ntt_sub_simd0(x2, x4, p);

    a13 = sp_ntt_add_simd0(t1, t2, p);
    a14 = sp_ntt_sub_simd0(t1, t2, p);
    a15 = t3;
    a16 = t4;
    a17 = sp_ntt_add_simd0(t3, t4, p);

    a12 = sp_ntt_add_simd0(x0, a13, p);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x2 = sp_simd_load(in +  1 * istride);
    sp_simd_t x3 = sp_simd_load(in +  8 * istride);
    sp_simd_t x0 = sp_simd_load(in + 15 * istride);
    sp_simd_t x1 = sp_simd_load(in + 22 * istride);
    sp_simd_t x4 = sp_simd_load(in + 29 * istride);

    t1 = sp_ntt_add_simd0(x1, x3, p);
    t3 = sp_ntt_sub_simd0(x1, x3, p);
    t2 = sp_ntt_add_simd0(x2, x4, p);
    t4 = sp_ntt_sub_simd0(x2, x4, p);

    a19 = sp_ntt_add_simd0(t1, t2, p);
    a20 = sp_ntt_sub_simd0(t1, t2, p);
    a21 = t3;
    a22 = t4;
    a23 = sp_ntt_add_simd0(t3, t4, p);

    a18 = sp_ntt_add_simd0(x0, a19, p);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x2 = sp_simd_load(in +  6 * istride);
    sp_simd_t x3 = sp_simd_load(in + 13 * istride);
    sp_simd_t x0 = sp_simd_load(in + 20 * istride);
    sp_simd_t x1 = sp_simd_load(in + 27 * istride);
    sp_simd_t x4 = sp_simd_load(in + 34 * istride);

    t1 = sp_ntt_add_simd0(x1, x3, p);
    t3 = sp_ntt_sub_simd0(x1, x3, p);
    t2 = sp_ntt_add_simd0(x2, x4, p);
    t4 = sp_ntt_sub_simd0(x2, x4, p);

    a25 = sp_ntt_add_simd0(t1, t2, p);
    a26 = sp_ntt_sub_simd0(t1, t2, p);
    a27 = t3;
    a28 = t4;
    a29 = sp_ntt_add_simd0(t3, t4, p);

    a24 = sp_ntt_add_simd0(x0, a25, p);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x4 = sp_simd_load(in +  4 * istride);
    sp_simd_t x2 = sp_simd_load(in + 11 * istride);
    sp_simd_t x3 = sp_simd_load(in + 18 * istride);
    sp_simd_t x0 = sp_simd_load(in + 25 * istride);
    sp_simd_t x1 = sp_simd_load(in + 32 * istride);

    t1 = sp_ntt_add_simd0(x1, x3, p);
    t3 = sp_ntt_sub_simd0(x1, x3, p);
    t2 = sp_ntt_add_simd0(x2, x4, p);
    t4 = sp_ntt_sub_simd0(x2, x4, p);

    a31 = sp_ntt_add_simd0(t1, t2, p);
    a32 = sp_ntt_sub_simd0(t1, t2, p);
    a33 = t3;
    a34 = t4;
    a35 = sp_ntt_add_simd0(t3, t4, p);

    a30 = sp_ntt_add_simd0(x0, a31, p);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x1 = sp_simd_load(in +  2 * istride);
    sp_simd_t x4 = sp_simd_load(in +  9 * istride);
    sp_simd_t x2 = sp_simd_load(in + 16 * istride);
    sp_simd_t x3 = sp_simd_load(in + 23 * istride);
    sp_simd_t x0 = sp_simd_load(in + 30 * istride);

    t1 = sp_ntt_add_simd0(x1, x3, p);
    t3 = sp_ntt_sub_simd0(x1, x3, p);
    t2 = sp_ntt_add_simd0(x2, x4, p);
    t4 = sp_ntt_sub_simd0(x2, x4, p);

    a37 = sp_ntt_add_simd0(t1, t2, p);
    a38 = sp_ntt_sub_simd0(t1, t2, p);
    a39 = t3;
    a40 = t4;
    a41 = sp_ntt_add_simd0(t3, t4, p);

    a36 = sp_ntt_add_simd0(x0, a37, p);
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a00;
    sp_simd_t x1 = a06;
    sp_simd_t x5 = a12;
    sp_simd_t x6 = a18;
    sp_simd_t x3 = a24;
    sp_simd_t x2 = a30;
    sp_simd_t x4 = a36;

    p1 = sp_ntt_add_simd0(x1, x4, p);
    p2 = sp_ntt_add_simd0(x2, x5, p);
    p3 = sp_ntt_add_simd0(x3, x6, p);
    p4 = sp_ntt_sub_simd0(x1, x4, p);
    p5 = sp_ntt_sub_simd0(x2, x5, p);
    p6 = sp_ntt_sub_simd0(x3, x6, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t1 = sp_ntt_add_simd0(t1, p3, p);
    t2 = sp_ntt_sub_simd0(p1, p3, p);
    t3 = sp_ntt_sub_simd0(p2, p3, p);
    t4 = sp_ntt_sub_simd0(p4, p5, p);
    t4 = sp_ntt_add_partial_simd0(t4, p6, p);
    t5 = sp_ntt_sub_simd0(p4, p6, p);
    t6 = sp_ntt_add_simd0(p5, p6, p);

    p0 = sp_ntt_add_simd0(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd0(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd0(t5, t6, p);

    p1 = sp_ntt_mul_simd0(p1, c+2*1, p);
    p2 = sp_ntt_mul_simd0(p2, c+2*2, p);
    p3 = sp_ntt_mul_simd0(p3, c+2*3, p);
    p4 = sp_ntt_mul_simd0(p4, c+2*4, p);
    p5 = sp_ntt_mul_simd0(p5, c+2*5, p);
    p6 = sp_ntt_mul_simd0(p6, c+2*6, p);
    p7 = sp_ntt_mul_simd0(p7, c+2*7, p);
    p8 = sp_ntt_mul_simd0(p8, c+2*8, p);

    t1 = sp_ntt_add_simd0(p0, p1, p);
    t2 = sp_ntt_add_simd0(p2, p4, p);
    t3 = sp_ntt_add_simd0(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd0(p6, p8, p);
    t6 = sp_ntt_add_simd0(p7, p8, p);

    p1 = sp_ntt_add_simd0(t1, t2, p);
    p2 = sp_ntt_add_simd0(t1, t3, p);
    p3 = sp_ntt_sub_simd0(t1, t2, p);
    p3 = sp_ntt_sub_simd0(p3, t3, p);
    p4 = sp_ntt_add_simd0(t4, t5, p);
    p5 = sp_ntt_sub_simd0(t6, t4, p);
    p6 = sp_ntt_sub_simd0(t4, t5, p);
    p6 = sp_ntt_add_simd0(p6, t6, p);

    t1 = sp_ntt_add_simd0(p1, p4, p);
    t2 = sp_ntt_add_simd0(p2, p5, p);
    t3 = sp_ntt_add_simd0(p3, p6, p);
    t4 = sp_ntt_sub_simd0(p1, p4, p);
    t5 = sp_ntt_sub_simd0(p2, p5, p);
    t6 = sp_ntt_sub_simd0(p3, p6, p);

    a00 = p0;
    a06 = t6;
    a12 = t4;
    a18 = t5;
    a24 = t2;
    a30 = t1;
    a36 = t3;
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a01;
    sp_simd_t x1 = a07;
    sp_simd_t x5 = a13;
    sp_simd_t x6 = a19;
    sp_simd_t x3 = a25;
    sp_simd_t x2 = a31;
    sp_simd_t x4 = a37;

    p1 = sp_ntt_add_simd0(x1, x4, p);
    p2 = sp_ntt_add_simd0(x2, x5, p);
    p3 = sp_ntt_add_simd0(x3, x6, p);
    p4 = sp_ntt_sub_simd0(x1, x4, p);
    p5 = sp_ntt_sub_simd0(x2, x5, p);
    p6 = sp_ntt_sub_simd0(x3, x6, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t1 = sp_ntt_add_simd0(t1, p3, p);
    t2 = sp_ntt_sub_simd0(p1, p3, p);
    t3 = sp_ntt_sub_simd0(p2, p3, p);
    t4 = sp_ntt_sub_simd0(p4, p5, p);
    t4 = sp_ntt_add_partial_simd0(t4, p6, p);
    t5 = sp_ntt_sub_simd0(p4, p6, p);
    t6 = sp_ntt_add_simd0(p5, p6, p);

    p0 = sp_ntt_add_partial_simd0(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd0(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd0(t5, t6, p);

    p0 = sp_ntt_mul_simd0(p0, c+2*9, p);
    p1 = sp_ntt_mul_simd0(p1, c+2*10, p);
    p2 = sp_ntt_mul_simd0(p2, c+2*11, p);
    p3 = sp_ntt_mul_simd0(p3, c+2*12, p);
    p4 = sp_ntt_mul_simd0(p4, c+2*13, p);
    p5 = sp_ntt_mul_simd0(p5, c+2*14, p);
    p6 = sp_ntt_mul_simd0(p6, c+2*15, p);
    p7 = sp_ntt_mul_simd0(p7, c+2*16, p);
    p8 = sp_ntt_mul_simd0(p8, c+2*17, p);

    t1 = sp_ntt_add_simd0(p0, p1, p);
    t2 = sp_ntt_add_simd0(p2, p4, p);
    t3 = sp_ntt_add_simd0(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd0(p6, p8, p);
    t6 = sp_ntt_add_simd0(p7, p8, p);

    p1 = sp_ntt_add_simd0(t1, t2, p);
    p2 = sp_ntt_add_simd0(t1, t3, p);
    p3 = sp_ntt_sub_simd0(t1, t2, p);
    p3 = sp_ntt_sub_simd0(p3, t3, p);
    p4 = sp_ntt_add_simd0(t4, t5, p);
    p5 = sp_ntt_sub_simd0(t6, t4, p);
    p6 = sp_ntt_sub_simd0(t4, t5, p);
    p6 = sp_ntt_add_simd0(p6, t6, p);

    t1 = sp_ntt_add_simd0(p1, p4, p);
    t2 = sp_ntt_add_simd0(p2, p5, p);
    t3 = sp_ntt_add_simd0(p3, p6, p);
    t4 = sp_ntt_sub_simd0(p1, p4, p);
    t5 = sp_ntt_sub_simd0(p2, p5, p);
    t6 = sp_ntt_sub_simd0(p3, p6, p);

    a01 = p0;
    a07 = t6;
    a13 = t4;
    a19 = t5;
    a25 = t2;
    a31 = t1;
    a37 = t3;
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a02;
    sp_simd_t x1 = a08;
    sp_simd_t x5 = a14;
    sp_simd_t x6 = a20;
    sp_simd_t x3 = a26;
    sp_simd_t x2 = a32;
    sp_simd_t x4 = a38;

    p1 = sp_ntt_add_simd0(x1, x4, p);
    p2 = sp_ntt_add_simd0(x2, x5, p);
    p3 = sp_ntt_add_simd0(x3, x6, p);
    p4 = sp_ntt_sub_simd0(x1, x4, p);
    p5 = sp_ntt_sub_simd0(x2, x5, p);
    p6 = sp_ntt_sub_simd0(x3, x6, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t1 = sp_ntt_add_simd0(t1, p3, p);
    t2 = sp_ntt_sub_simd0(p1, p3, p);
    t3 = sp_ntt_sub_simd0(p2, p3, p);
    t4 = sp_ntt_sub_simd0(p4, p5, p);
    t4 = sp_ntt_add_partial_simd0(t4, p6, p);
    t5 = sp_ntt_sub_simd0(p4, p6, p);
    t6 = sp_ntt_add_simd0(p5, p6, p);

    p0 = sp_ntt_add_partial_simd0(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd0(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd0(t5, t6, p);

    p0 = sp_ntt_mul_simd0(p0, c+2*18, p);
    p1 = sp_ntt_mul_simd0(p1, c+2*19, p);
    p2 = sp_ntt_mul_simd0(p2, c+2*20, p);
    p3 = sp_ntt_mul_simd0(p3, c+2*21, p);
    p4 = sp_ntt_mul_simd0(p4, c+2*22, p);
    p5 = sp_ntt_mul_simd0(p5, c+2*23, p);
    p6 = sp_ntt_mul_simd0(p6, c+2*24, p);
    p7 = sp_ntt_mul_simd0(p7, c+2*25, p);
    p8 = sp_ntt_mul_simd0(p8, c+2*26, p);

    t1 = sp_ntt_add_simd0(p0, p1, p);
    t2 = sp_ntt_add_simd0(p2, p4, p);
    t3 = sp_ntt_add_simd0(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd0(p6, p8, p);
    t6 = sp_ntt_add_simd0(p7, p8, p);

    p1 = sp_ntt_add_simd0(t1, t2, p);
    p2 = sp_ntt_add_simd0(t1, t3, p);
    p3 = sp_ntt_sub_simd0(t1, t2, p);
    p3 = sp_ntt_sub_simd0(p3, t3, p);
    p4 = sp_ntt_add_simd0(t4, t5, p);
    p5 = sp_ntt_sub_simd0(t6, t4, p);
    p6 = sp_ntt_sub_simd0(t4, t5, p);
    p6 = sp_ntt_add_simd0(p6, t6, p);

    t1 = sp_ntt_add_simd0(p1, p4, p);
    t2 = sp_ntt_add_simd0(p2, p5, p);
    t3 = sp_ntt_add_simd0(p3, p6, p);
    t4 = sp_ntt_sub_simd0(p1, p4, p);
    t5 = sp_ntt_sub_simd0(p2, p5, p);
    t6 = sp_ntt_sub_simd0(p3, p6, p);

    a02 = p0;
    a08 = t6;
    a14 = t4;
    a20 = t5;
    a26 = t2;
    a32 = t1;
    a38 = t3;
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a03;
    sp_simd_t x1 = a09;
    sp_simd_t x5 = a15;
    sp_simd_t x6 = a21;
    sp_simd_t x3 = a27;
    sp_simd_t x2 = a33;
    sp_simd_t x4 = a39;

    p1 = sp_ntt_add_simd0(x1, x4, p);
    p2 = sp_ntt_add_simd0(x2, x5, p);
    p3 = sp_ntt_add_simd0(x3, x6, p);
    p4 = sp_ntt_sub_simd0(x1, x4, p);
    p5 = sp_ntt_sub_simd0(x2, x5, p);
    p6 = sp_ntt_sub_simd0(x3, x6, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t1 = sp_ntt_add_simd0(t1, p3, p);
    t2 = sp_ntt_sub_simd0(p1, p3, p);
    t3 = sp_ntt_sub_simd0(p2, p3, p);
    t4 = sp_ntt_sub_simd0(p4, p5, p);
    t4 = sp_ntt_add_partial_simd0(t4, p6, p);
    t5 = sp_ntt_sub_simd0(p4, p6, p);
    t6 = sp_ntt_add_simd0(p5, p6, p);

    p0 = sp_ntt_add_partial_simd0(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd0(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd0(t5, t6, p);

    p0 = sp_ntt_mul_simd0(p0, c+2*27, p);
    p1 = sp_ntt_mul_simd0(p1, c+2*28, p);
    p2 = sp_ntt_mul_simd0(p2, c+2*29, p);
    p3 = sp_ntt_mul_simd0(p3, c+2*30, p);
    p4 = sp_ntt_mul_simd0(p4, c+2*31, p);
    p5 = sp_ntt_mul_simd0(p5, c+2*32, p);
    p6 = sp_ntt_mul_simd0(p6, c+2*33, p);
    p7 = sp_ntt_mul_simd0(p7, c+2*34, p);
    p8 = sp_ntt_mul_simd0(p8, c+2*35, p);

    t1 = sp_ntt_add_simd0(p0, p1, p);
    t2 = sp_ntt_add_simd0(p2, p4, p);
    t3 = sp_ntt_add_simd0(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd0(p6, p8, p);
    t6 = sp_ntt_add_simd0(p7, p8, p);

    p1 = sp_ntt_add_simd0(t1, t2, p);
    p2 = sp_ntt_add_simd0(t1, t3, p);
    p3 = sp_ntt_sub_simd0(t1, t2, p);
    p3 = sp_ntt_sub_simd0(p3, t3, p);
    p4 = sp_ntt_add_simd0(t4, t5, p);
    p5 = sp_ntt_sub_simd0(t6, t4, p);
    p6 = sp_ntt_sub_simd0(t4, t5, p);
    p6 = sp_ntt_add_simd0(p6, t6, p);

    t1 = sp_ntt_add_simd0(p1, p4, p);
    t2 = sp_ntt_add_simd0(p2, p5, p);
    t3 = sp_ntt_add_simd0(p3, p6, p);
    t4 = sp_ntt_sub_simd0(p1, p4, p);
    t5 = sp_ntt_sub_simd0(p2, p5, p);
    t6 = sp_ntt_sub_simd0(p3, p6, p);

    a03 = p0;
    a09 = t6;
    a15 = t4;
    a21 = t5;
    a27 = t2;
    a33 = t1;
    a39 = t3;
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a04;
    sp_simd_t x1 = a10;
    sp_simd_t x5 = a16;
    sp_simd_t x6 = a22;
    sp_simd_t x3 = a28;
    sp_simd_t x2 = a34;
    sp_simd_t x4 = a40;

    p1 = sp_ntt_add_simd0(x1, x4, p);
    p2 = sp_ntt_add_simd0(x2, x5, p);
    p3 = sp_ntt_add_simd0(x3, x6, p);
    p4 = sp_ntt_sub_simd0(x1, x4, p);
    p5 = sp_ntt_sub_simd0(x2, x5, p);
    p6 = sp_ntt_sub_simd0(x3, x6, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t1 = sp_ntt_add_simd0(t1, p3, p);
    t2 = sp_ntt_sub_simd0(p1, p3, p);
    t3 = sp_ntt_sub_simd0(p2, p3, p);
    t4 = sp_ntt_sub_simd0(p4, p5, p);
    t4 = sp_ntt_add_partial_simd0(t4, p6, p);
    t5 = sp_ntt_sub_simd0(p4, p6, p);
    t6 = sp_ntt_add_simd0(p5, p6, p);

    p0 = sp_ntt_add_partial_simd0(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd0(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd0(t5, t6, p);

    p0 = sp_ntt_mul_simd0(p0, c+2*36, p);
    p1 = sp_ntt_mul_simd0(p1, c+2*37, p);
    p2 = sp_ntt_mul_simd0(p2, c+2*38, p);
    p3 = sp_ntt_mul_simd0(p3, c+2*39, p);
    p4 = sp_ntt_mul_simd0(p4, c+2*40, p);
    p5 = sp_ntt_mul_simd0(p5, c+2*41, p);
    p6 = sp_ntt_mul_simd0(p6, c+2*42, p);
    p7 = sp_ntt_mul_simd0(p7, c+2*43, p);
    p8 = sp_ntt_mul_simd0(p8, c+2*44, p);

    t1 = sp_ntt_add_simd0(p0, p1, p);
    t2 = sp_ntt_add_simd0(p2, p4, p);
    t3 = sp_ntt_add_simd0(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd0(p6, p8, p);
    t6 = sp_ntt_add_simd0(p7, p8, p);

    p1 = sp_ntt_add_simd0(t1, t2, p);
    p2 = sp_ntt_add_simd0(t1, t3, p);
    p3 = sp_ntt_sub_simd0(t1, t2, p);
    p3 = sp_ntt_sub_simd0(p3, t3, p);
    p4 = sp_ntt_add_simd0(t4, t5, p);
    p5 = sp_ntt_sub_simd0(t6, t4, p);
    p6 = sp_ntt_sub_simd0(t4, t5, p);
    p6 = sp_ntt_add_simd0(p6, t6, p);

    t1 = sp_ntt_add_simd0(p1, p4, p);
    t2 = sp_ntt_add_simd0(p2, p5, p);
    t3 = sp_ntt_add_simd0(p3, p6, p);
    t4 = sp_ntt_sub_simd0(p1, p4, p);
    t5 = sp_ntt_sub_simd0(p2, p5, p);
    t6 = sp_ntt_sub_simd0(p3, p6, p);

    a04 = p0;
    a10 = t6;
    a16 = t4;
    a22 = t5;
    a28 = t2;
    a34 = t1;
    a40 = t3;
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a05;
    sp_simd_t x1 = a11;
    sp_simd_t x5 = a17;
    sp_simd_t x6 = a23;
    sp_simd_t x3 = a29;
    sp_simd_t x2 = a35;
    sp_simd_t x4 = a41;

    p1 = sp_ntt_add_simd0(x1, x4, p);
    p2 = sp_ntt_add_simd0(x2, x5, p);
    p3 = sp_ntt_add_simd0(x3, x6, p);
    p4 = sp_ntt_sub_simd0(x1, x4, p);
    p5 = sp_ntt_sub_simd0(x2, x5, p);
    p6 = sp_ntt_sub_simd0(x3, x6, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t1 = sp_ntt_add_simd0(t1, p3, p);
    t2 = sp_ntt_sub_simd0(p1, p3, p);
    t3 = sp_ntt_sub_simd0(p2, p3, p);
    t4 = sp_ntt_sub_simd0(p4, p5, p);
    t4 = sp_ntt_add_partial_simd0(t4, p6, p);
    t5 = sp_ntt_sub_simd0(p4, p6, p);
    t6 = sp_ntt_add_simd0(p5, p6, p);

    p0 = sp_ntt_add_partial_simd0(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd0(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd0(t5, t6, p);

    p0 = sp_ntt_mul_simd0(p0, c+2*45, p);
    p1 = sp_ntt_mul_simd0(p1, c+2*46, p);
    p2 = sp_ntt_mul_simd0(p2, c+2*47, p);
    p3 = sp_ntt_mul_simd0(p3, c+2*48, p);
    p4 = sp_ntt_mul_simd0(p4, c+2*49, p);
    p5 = sp_ntt_mul_simd0(p5, c+2*50, p);
    p6 = sp_ntt_mul_simd0(p6, c+2*51, p);
    p7 = sp_ntt_mul_simd0(p7, c+2*52, p);
    p8 = sp_ntt_mul_simd0(p8, c+2*53, p);

    t1 = sp_ntt_add_simd0(p0, p1, p);
    t2 = sp_ntt_add_simd0(p2, p4, p);
    t3 = sp_ntt_add_simd0(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd0(p6, p8, p);
    t6 = sp_ntt_add_simd0(p7, p8, p);

    p1 = sp_ntt_add_simd0(t1, t2, p);
    p2 = sp_ntt_add_simd0(t1, t3, p);
    p3 = sp_ntt_sub_simd0(t1, t2, p);
    p3 = sp_ntt_sub_simd0(p3, t3, p);
    p4 = sp_ntt_add_simd0(t4, t5, p);
    p5 = sp_ntt_sub_simd0(t6, t4, p);
    p6 = sp_ntt_sub_simd0(t4, t5, p);
    p6 = sp_ntt_add_simd0(p6, t6, p);

    t1 = sp_ntt_add_simd0(p1, p4, p);
    t2 = sp_ntt_add_simd0(p2, p5, p);
    t3 = sp_ntt_add_simd0(p3, p6, p);
    t4 = sp_ntt_sub_simd0(p1, p4, p);
    t5 = sp_ntt_sub_simd0(p2, p5, p);
    t6 = sp_ntt_sub_simd0(p3, p6, p);

    a05 = p0;
    a11 = t6;
    a17 = t4;
    a23 = t5;
    a29 = t2;
    a35 = t1;
    a41 = t3;
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a00;
    sp_simd_t p1 = a01;
    sp_simd_t p2 = a02;
    sp_simd_t p3 = a03;
    sp_simd_t p4 = a04;
    sp_simd_t p5 = a05;

    p1 = sp_ntt_add_simd0(p0, p1, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t2 = sp_ntt_sub_simd0(p1, p2, p);
    t3 = sp_ntt_add_simd0(p3, p5, p);
    t4 = sp_ntt_add_simd0(p4, p5, p);

    p1 = sp_ntt_add_partial_simd0(t1, t3, p);
    p2 = sp_ntt_add_partial_simd0(t2, t4, p);
    p3 = sp_ntt_sub_partial_simd0(t1, t3, p);
    p4 = sp_ntt_sub_partial_simd0(t2, t4, p);

    p3 = sp_ntt_twiddle_mul_simd0(p3, w+12, p);
    p2 = sp_ntt_twiddle_mul_simd0(p2, w+26, p);
    p4 = sp_ntt_twiddle_mul_simd0(p4, w+40, p);
    p1 = sp_ntt_twiddle_mul_simd0(p1, w+54, p);

    sp_simd_store(p0, out +  0 * ostride);
    sp_simd_store(p3, out +  7 * ostride);
    sp_simd_store(p2, out + 14 * ostride);
    sp_simd_store(p4, out + 21 * ostride);
    sp_simd_store(p1, out + 28 * ostride);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a06;
    sp_simd_t p1 = a07;
    sp_simd_t p2 = a08;
    sp_simd_t p3 = a09;
    sp_simd_t p4 = a10;
    sp_simd_t p5 = a11;

    p1 = sp_ntt_add_simd0(p0, p1, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t2 = sp_ntt_sub_simd0(p1, p2, p);
    t3 = sp_ntt_add_simd0(p3, p5, p);
    t4 = sp_ntt_add_simd0(p4, p5, p);

    p1 = sp_ntt_add_partial_simd0(t1, t3, p);
    p2 = sp_ntt_add_partial_simd0(t2, t4, p);
    p3 = sp_ntt_sub_partial_simd0(t1, t3, p);
    p4 = sp_ntt_sub_partial_simd0(t2, t4, p);

    p4 = sp_ntt_twiddle_mul_simd0(p4, w+0, p);
    p1 = sp_ntt_twiddle_mul_simd0(p1, w+14, p);
    p0 = sp_ntt_twiddle_mul_simd0(p0, w+28, p);
    p3 = sp_ntt_twiddle_mul_simd0(p3, w+42, p);
    p2 = sp_ntt_twiddle_mul_simd0(p2, w+56, p);

    sp_simd_store(p4, out +  1 * ostride);
    sp_simd_store(p1, out +  8 * ostride);
    sp_simd_store(p0, out + 15 * ostride);
    sp_simd_store(p3, out + 22 * ostride);
    sp_simd_store(p2, out + 29 * ostride);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a12;
    sp_simd_t p1 = a13;
    sp_simd_t p2 = a14;
    sp_simd_t p3 = a15;
    sp_simd_t p4 = a16;
    sp_simd_t p5 = a17;

    p1 = sp_ntt_add_simd0(p0, p1, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t2 = sp_ntt_sub_simd0(p1, p2, p);
    t3 = sp_ntt_add_simd0(p3, p5, p);
    t4 = sp_ntt_add_simd0(p4, p5, p);

    p1 = sp_ntt_add_partial_simd0(t1, t3, p);
    p2 = sp_ntt_add_partial_simd0(t2, t4, p);
    p3 = sp_ntt_sub_partial_simd0(t1, t3, p);
    p4 = sp_ntt_sub_partial_simd0(t2, t4, p);

    p3 = sp_ntt_twiddle_mul_simd0(p3, w+2, p);
    p2 = sp_ntt_twiddle_mul_simd0(p2, w+16, p);
    p4 = sp_ntt_twiddle_mul_simd0(p4, w+30, p);
    p1 = sp_ntt_twiddle_mul_simd0(p1, w+44, p);
    p0 = sp_ntt_twiddle_mul_simd0(p0, w+58, p);

    sp_simd_store(p3, out +  2 * ostride);
    sp_simd_store(p2, out +  9 * ostride);
    sp_simd_store(p4, out + 16 * ostride);
    sp_simd_store(p1, out + 23 * ostride);
    sp_simd_store(p0, out + 30 * ostride);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a18;
    sp_simd_t p1 = a19;
    sp_simd_t p2 = a20;
    sp_simd_t p3 = a21;
    sp_simd_t p4 = a22;
    sp_simd_t p5 = a23;

    p1 = sp_ntt_add_simd0(p0, p1, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t2 = sp_ntt_sub_simd0(p1, p2, p);
    t3 = sp_ntt_add_simd0(p3, p5, p);
    t4 = sp_ntt_add_simd0(p4, p5, p);

    p1 = sp_ntt_add_partial_simd0(t1, t3, p);
    p2 = sp_ntt_add_partial_simd0(t2, t4, p);
    p3 = sp_ntt_sub_partial_simd0(t1, t3, p);
    p4 = sp_ntt_sub_partial_simd0(t2, t4, p);

    p1 = sp_ntt_twiddle_mul_simd0(p1, w+4, p);
    p0 = sp_ntt_twiddle_mul_simd0(p0, w+18, p);
    p3 = sp_ntt_twiddle_mul_simd0(p3, w+32, p);
    p2 = sp_ntt_twiddle_mul_simd0(p2, w+46, p);
    p4 = sp_ntt_twiddle_mul_simd0(p4, w+60, p);

    sp_simd_store(p1, out +  3 * ostride);
    sp_simd_store(p0, out + 10 * ostride);
    sp_simd_store(p3, out + 17 * ostride);
    sp_simd_store(p2, out + 24 * ostride);
    sp_simd_store(p4, out + 31 * ostride);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a24;
    sp_simd_t p1 = a25;
    sp_simd_t p2 = a26;
    sp_simd_t p3 = a27;
    sp_simd_t p4 = a28;
    sp_simd_t p5 = a29;

    p1 = sp_ntt_add_simd0(p0, p1, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t2 = sp_ntt_sub_simd0(p1, p2, p);
    t3 = sp_ntt_add_simd0(p3, p5, p);
    t4 = sp_ntt_add_simd0(p4, p5, p);

    p1 = sp_ntt_add_partial_simd0(t1, t3, p);
    p2 = sp_ntt_add_partial_simd0(t2, t4, p);
    p3 = sp_ntt_sub_partial_simd0(t1, t3, p);
    p4 = sp_ntt_sub_partial_simd0(t2, t4, p);

    p2 = sp_ntt_twiddle_mul_simd0(p2, w+6, p);
    p4 = sp_ntt_twiddle_mul_simd0(p4, w+20, p);
    p1 = sp_ntt_twiddle_mul_simd0(p1, w+34, p);
    p0 = sp_ntt_twiddle_mul_simd0(p0, w+48, p);
    p3 = sp_ntt_twiddle_mul_simd0(p3, w+62, p);

    sp_simd_store(p2, out +  4 * ostride);
    sp_simd_store(p4, out + 11 * ostride);
    sp_simd_store(p1, out + 18 * ostride);
    sp_simd_store(p0, out + 25 * ostride);
    sp_simd_store(p3, out + 32 * ostride);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a30;
    sp_simd_t p1 = a31;
    sp_simd_t p2 = a32;
    sp_simd_t p3 = a33;
    sp_simd_t p4 = a34;
    sp_simd_t p5 = a35;

    p1 = sp_ntt_add_simd0(p0, p1, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t2 = sp_ntt_sub_simd0(p1, p2, p);
    t3 = sp_ntt_add_simd0(p3, p5, p);
    t4 = sp_ntt_add_simd0(p4, p5, p);

    p1 = sp_ntt_add_partial_simd0(t1, t3, p);
    p2 = sp_ntt_add_partial_simd0(t2, t4, p);
    p3 = sp_ntt_sub_partial_simd0(t1, t3, p);
    p4 = sp_ntt_sub_partial_simd0(t2, t4, p);

    p0 = sp_ntt_twiddle_mul_simd0(p0, w+8, p);
    p3 = sp_ntt_twiddle_mul_simd0(p3, w+22, p);
    p2 = sp_ntt_twiddle_mul_simd0(p2, w+36, p);
    p4 = sp_ntt_twiddle_mul_simd0(p4, w+50, p);
    p1 = sp_ntt_twiddle_mul_simd0(p1, w+64, p);

    sp_simd_store(p0, out +  5 * ostride);
    sp_simd_store(p3, out + 12 * ostride);
    sp_simd_store(p2, out + 19 * ostride);
    sp_simd_store(p4, out + 26 * ostride);
    sp_simd_store(p1, out + 33 * ostride);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a36;
    sp_simd_t p1 = a37;
    sp_simd_t p2 = a38;
    sp_simd_t p3 = a39;
    sp_simd_t p4 = a40;
    sp_simd_t p5 = a41;

    p1 = sp_ntt_add_simd0(p0, p1, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t2 = sp_ntt_sub_simd0(p1, p2, p);
    t3 = sp_ntt_add_simd0(p3, p5, p);
    t4 = sp_ntt_add_simd0(p4, p5, p);

    p1 = sp_ntt_add_partial_simd0(t1, t3, p);
    p2 = sp_ntt_add_partial_simd0(t2, t4, p);
    p3 = sp_ntt_sub_partial_simd0(t1, t3, p);
    p4 = sp_ntt_sub_partial_simd0(t2, t4, p);

    p4 = sp_ntt_twiddle_mul_simd0(p4, w+10, p);
    p1 = sp_ntt_twiddle_mul_simd0(p1, w+24, p);
    p0 = sp_ntt_twiddle_mul_simd0(p0, w+38, p);
    p3 = sp_ntt_twiddle_mul_simd0(p3, w+52, p);
    p2 = sp_ntt_twiddle_mul_simd0(p2, w+66, p);

    sp_simd_store(p4, out +  6 * ostride);
    sp_simd_store(p1, out + 13 * ostride);
    sp_simd_store(p0, out + 20 * ostride);
    sp_simd_store(p3, out + 27 * ostride);
    sp_simd_store(p2, out + 34 * ostride);
  }
}

static void
ntt35_pfa_run_core_simd(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t inc2, spv_size_t n,
	  sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  spv_size_t j00, j01, j02, j03, j04, j05, j06, j07,
	     j08, j09, j10, j11, j12, j13, j14, j15,
	     j16, j17, j18, j19, j20, j21, j22, j23,
	     j24, j25, j26, j27, j28, j29, j30, j31,
	     j32, j33, j34;

  sp_simd_t a00, a01, a02, a03, a04, a05, a06, 
            a07, a08, a09, a10, a11, a12, a13, 
            a14, a15, a16, a17, a18, a19, a20, 
            a21, a22, a23, a24, a25, a26, a27, 
            a28, a29, a30, a31, a32, a33, a34, 
            a35, a36, a37, a38, a39, a40, a41;

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

  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x0 = sp_simd_pfa_gather(x, j00, inc2, n, vsize);
    sp_simd_t x1 = sp_simd_pfa_gather(x, j07, inc2, n, vsize);
    sp_simd_t x4 = sp_simd_pfa_gather(x, j14, inc2, n, vsize);
    sp_simd_t x2 = sp_simd_pfa_gather(x, j21, inc2, n, vsize);
    sp_simd_t x3 = sp_simd_pfa_gather(x, j28, inc2, n, vsize);

    t1 = sp_ntt_add_simd(x1, x3, p);
    t3 = sp_ntt_sub_simd(x1, x3, p);
    t2 = sp_ntt_add_simd(x2, x4, p);
    t4 = sp_ntt_sub_simd(x2, x4, p);

    a01 = sp_ntt_add_simd(t1, t2, p);
    a02 = sp_ntt_sub_simd(t1, t2, p);
    a03 = t3;
    a04 = t4;
    a05 = sp_ntt_add_simd(t3, t4, p);

    a00 = sp_ntt_add_simd(x0, a01, p);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x0 = sp_simd_pfa_gather(x, j05, inc2, n, vsize);
    sp_simd_t x1 = sp_simd_pfa_gather(x, j12, inc2, n, vsize);
    sp_simd_t x4 = sp_simd_pfa_gather(x, j19, inc2, n, vsize);
    sp_simd_t x2 = sp_simd_pfa_gather(x, j26, inc2, n, vsize);
    sp_simd_t x3 = sp_simd_pfa_gather(x, j33, inc2, n, vsize);

    t1 = sp_ntt_add_simd(x1, x3, p);
    t3 = sp_ntt_sub_simd(x1, x3, p);
    t2 = sp_ntt_add_simd(x2, x4, p);
    t4 = sp_ntt_sub_simd(x2, x4, p);

    a07 = sp_ntt_add_simd(t1, t2, p);
    a08 = sp_ntt_sub_simd(t1, t2, p);
    a09 = t3;
    a10 = t4;
    a11 = sp_ntt_add_simd(t3, t4, p);

    a06 = sp_ntt_add_simd(x0, a07, p);
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x3 = sp_simd_pfa_gather(x, j03, inc2, n, vsize);
    sp_simd_t x0 = sp_simd_pfa_gather(x, j10, inc2, n, vsize);
    sp_simd_t x1 = sp_simd_pfa_gather(x, j17, inc2, n, vsize);
    sp_simd_t x4 = sp_simd_pfa_gather(x, j24, inc2, n, vsize);
    sp_simd_t x2 = sp_simd_pfa_gather(x, j31, inc2, n, vsize);

    t1 = sp_ntt_add_simd(x1, x3, p);
    t3 = sp_ntt_sub_simd(x1, x3, p);
    t2 = sp_ntt_add_simd(x2, x4, p);
    t4 = sp_ntt_sub_simd(x2, x4, p);

    a13 = sp_ntt_add_simd(t1, t2, p);
    a14 = sp_ntt_sub_simd(t1, t2, p);
    a15 = t3;
    a16 = t4;
    a17 = sp_ntt_add_simd(t3, t4, p);

    a12 = sp_ntt_add_simd(x0, a13, p);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x2 = sp_simd_pfa_gather(x, j01, inc2, n, vsize);
    sp_simd_t x3 = sp_simd_pfa_gather(x, j08, inc2, n, vsize);
    sp_simd_t x0 = sp_simd_pfa_gather(x, j15, inc2, n, vsize);
    sp_simd_t x1 = sp_simd_pfa_gather(x, j22, inc2, n, vsize);
    sp_simd_t x4 = sp_simd_pfa_gather(x, j29, inc2, n, vsize);

    t1 = sp_ntt_add_simd(x1, x3, p);
    t3 = sp_ntt_sub_simd(x1, x3, p);
    t2 = sp_ntt_add_simd(x2, x4, p);
    t4 = sp_ntt_sub_simd(x2, x4, p);

    a19 = sp_ntt_add_simd(t1, t2, p);
    a20 = sp_ntt_sub_simd(t1, t2, p);
    a21 = t3;
    a22 = t4;
    a23 = sp_ntt_add_simd(t3, t4, p);

    a18 = sp_ntt_add_simd(x0, a19, p);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x2 = sp_simd_pfa_gather(x, j06, inc2, n, vsize);
    sp_simd_t x3 = sp_simd_pfa_gather(x, j13, inc2, n, vsize);
    sp_simd_t x0 = sp_simd_pfa_gather(x, j20, inc2, n, vsize);
    sp_simd_t x1 = sp_simd_pfa_gather(x, j27, inc2, n, vsize);
    sp_simd_t x4 = sp_simd_pfa_gather(x, j34, inc2, n, vsize);

    t1 = sp_ntt_add_simd(x1, x3, p);
    t3 = sp_ntt_sub_simd(x1, x3, p);
    t2 = sp_ntt_add_simd(x2, x4, p);
    t4 = sp_ntt_sub_simd(x2, x4, p);

    a25 = sp_ntt_add_simd(t1, t2, p);
    a26 = sp_ntt_sub_simd(t1, t2, p);
    a27 = t3;
    a28 = t4;
    a29 = sp_ntt_add_simd(t3, t4, p);

    a24 = sp_ntt_add_simd(x0, a25, p);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x4 = sp_simd_pfa_gather(x, j04, inc2, n, vsize);
    sp_simd_t x2 = sp_simd_pfa_gather(x, j11, inc2, n, vsize);
    sp_simd_t x3 = sp_simd_pfa_gather(x, j18, inc2, n, vsize);
    sp_simd_t x0 = sp_simd_pfa_gather(x, j25, inc2, n, vsize);
    sp_simd_t x1 = sp_simd_pfa_gather(x, j32, inc2, n, vsize);

    t1 = sp_ntt_add_simd(x1, x3, p);
    t3 = sp_ntt_sub_simd(x1, x3, p);
    t2 = sp_ntt_add_simd(x2, x4, p);
    t4 = sp_ntt_sub_simd(x2, x4, p);

    a31 = sp_ntt_add_simd(t1, t2, p);
    a32 = sp_ntt_sub_simd(t1, t2, p);
    a33 = t3;
    a34 = t4;
    a35 = sp_ntt_add_simd(t3, t4, p);

    a30 = sp_ntt_add_simd(x0, a31, p);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x1 = sp_simd_pfa_gather(x, j02, inc2, n, vsize);
    sp_simd_t x4 = sp_simd_pfa_gather(x, j09, inc2, n, vsize);
    sp_simd_t x2 = sp_simd_pfa_gather(x, j16, inc2, n, vsize);
    sp_simd_t x3 = sp_simd_pfa_gather(x, j23, inc2, n, vsize);
    sp_simd_t x0 = sp_simd_pfa_gather(x, j30, inc2, n, vsize);

    t1 = sp_ntt_add_simd(x1, x3, p);
    t3 = sp_ntt_sub_simd(x1, x3, p);
    t2 = sp_ntt_add_simd(x2, x4, p);
    t4 = sp_ntt_sub_simd(x2, x4, p);

    a37 = sp_ntt_add_simd(t1, t2, p);
    a38 = sp_ntt_sub_simd(t1, t2, p);
    a39 = t3;
    a40 = t4;
    a41 = sp_ntt_add_simd(t3, t4, p);

    a36 = sp_ntt_add_simd(x0, a37, p);
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a00;
    sp_simd_t x1 = a06;
    sp_simd_t x5 = a12;
    sp_simd_t x6 = a18;
    sp_simd_t x3 = a24;
    sp_simd_t x2 = a30;
    sp_simd_t x4 = a36;

    p1 = sp_ntt_add_simd(x1, x4, p);
    p2 = sp_ntt_add_simd(x2, x5, p);
    p3 = sp_ntt_add_simd(x3, x6, p);
    p4 = sp_ntt_sub_simd(x1, x4, p);
    p5 = sp_ntt_sub_simd(x2, x5, p);
    p6 = sp_ntt_sub_simd(x3, x6, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t1 = sp_ntt_add_simd(t1, p3, p);
    t2 = sp_ntt_sub_simd(p1, p3, p);
    t3 = sp_ntt_sub_simd(p2, p3, p);
    t4 = sp_ntt_sub_simd(p4, p5, p);
    t4 = sp_ntt_add_partial_simd(t4, p6, p);
    t5 = sp_ntt_sub_simd(p4, p6, p);
    t6 = sp_ntt_add_simd(p5, p6, p);

    p0 = sp_ntt_add_simd(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd(t5, t6, p);

    p1 = sp_ntt_mul_simd(p1, ntt_const[1], ntt_const[NC+1], p);
    p2 = sp_ntt_mul_simd(p2, ntt_const[2], ntt_const[NC+2], p);
    p3 = sp_ntt_mul_simd(p3, ntt_const[3], ntt_const[NC+3], p);
    p4 = sp_ntt_mul_simd(p4, ntt_const[4], ntt_const[NC+4], p);
    p5 = sp_ntt_mul_simd(p5, ntt_const[5], ntt_const[NC+5], p);
    p6 = sp_ntt_mul_simd(p6, ntt_const[6], ntt_const[NC+6], p);
    p7 = sp_ntt_mul_simd(p7, ntt_const[7], ntt_const[NC+7], p);
    p8 = sp_ntt_mul_simd(p8, ntt_const[8], ntt_const[NC+8], p);

    t1 = sp_ntt_add_simd(p0, p1, p);
    t2 = sp_ntt_add_simd(p2, p4, p);
    t3 = sp_ntt_add_simd(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd(p6, p8, p);
    t6 = sp_ntt_add_simd(p7, p8, p);

    p1 = sp_ntt_add_simd(t1, t2, p);
    p2 = sp_ntt_add_simd(t1, t3, p);
    p3 = sp_ntt_sub_simd(t1, t2, p);
    p3 = sp_ntt_sub_simd(p3, t3, p);
    p4 = sp_ntt_add_simd(t4, t5, p);
    p5 = sp_ntt_sub_simd(t6, t4, p);
    p6 = sp_ntt_sub_simd(t4, t5, p);
    p6 = sp_ntt_add_simd(p6, t6, p);

    t1 = sp_ntt_add_simd(p1, p4, p);
    t2 = sp_ntt_add_simd(p2, p5, p);
    t3 = sp_ntt_add_simd(p3, p6, p);
    t4 = sp_ntt_sub_simd(p1, p4, p);
    t5 = sp_ntt_sub_simd(p2, p5, p);
    t6 = sp_ntt_sub_simd(p3, p6, p);

    a00 = p0;
    a06 = t6;
    a12 = t4;
    a18 = t5;
    a24 = t2;
    a30 = t1;
    a36 = t3;
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a01;
    sp_simd_t x1 = a07;
    sp_simd_t x5 = a13;
    sp_simd_t x6 = a19;
    sp_simd_t x3 = a25;
    sp_simd_t x2 = a31;
    sp_simd_t x4 = a37;

    p1 = sp_ntt_add_simd(x1, x4, p);
    p2 = sp_ntt_add_simd(x2, x5, p);
    p3 = sp_ntt_add_simd(x3, x6, p);
    p4 = sp_ntt_sub_simd(x1, x4, p);
    p5 = sp_ntt_sub_simd(x2, x5, p);
    p6 = sp_ntt_sub_simd(x3, x6, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t1 = sp_ntt_add_simd(t1, p3, p);
    t2 = sp_ntt_sub_simd(p1, p3, p);
    t3 = sp_ntt_sub_simd(p2, p3, p);
    t4 = sp_ntt_sub_simd(p4, p5, p);
    t4 = sp_ntt_add_partial_simd(t4, p6, p);
    t5 = sp_ntt_sub_simd(p4, p6, p);
    t6 = sp_ntt_add_simd(p5, p6, p);

    p0 = sp_ntt_add_partial_simd(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd(t5, t6, p);

    p0 = sp_ntt_mul_simd(p0, ntt_const[9], ntt_const[NC+9], p);
    p1 = sp_ntt_mul_simd(p1, ntt_const[10], ntt_const[NC+10], p);
    p2 = sp_ntt_mul_simd(p2, ntt_const[11], ntt_const[NC+11], p);
    p3 = sp_ntt_mul_simd(p3, ntt_const[12], ntt_const[NC+12], p);
    p4 = sp_ntt_mul_simd(p4, ntt_const[13], ntt_const[NC+13], p);
    p5 = sp_ntt_mul_simd(p5, ntt_const[14], ntt_const[NC+14], p);
    p6 = sp_ntt_mul_simd(p6, ntt_const[15], ntt_const[NC+15], p);
    p7 = sp_ntt_mul_simd(p7, ntt_const[16], ntt_const[NC+16], p);
    p8 = sp_ntt_mul_simd(p8, ntt_const[17], ntt_const[NC+17], p);

    t1 = sp_ntt_add_simd(p0, p1, p);
    t2 = sp_ntt_add_simd(p2, p4, p);
    t3 = sp_ntt_add_simd(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd(p6, p8, p);
    t6 = sp_ntt_add_simd(p7, p8, p);

    p1 = sp_ntt_add_simd(t1, t2, p);
    p2 = sp_ntt_add_simd(t1, t3, p);
    p3 = sp_ntt_sub_simd(t1, t2, p);
    p3 = sp_ntt_sub_simd(p3, t3, p);
    p4 = sp_ntt_add_simd(t4, t5, p);
    p5 = sp_ntt_sub_simd(t6, t4, p);
    p6 = sp_ntt_sub_simd(t4, t5, p);
    p6 = sp_ntt_add_simd(p6, t6, p);

    t1 = sp_ntt_add_simd(p1, p4, p);
    t2 = sp_ntt_add_simd(p2, p5, p);
    t3 = sp_ntt_add_simd(p3, p6, p);
    t4 = sp_ntt_sub_simd(p1, p4, p);
    t5 = sp_ntt_sub_simd(p2, p5, p);
    t6 = sp_ntt_sub_simd(p3, p6, p);

    a01 = p0;
    a07 = t6;
    a13 = t4;
    a19 = t5;
    a25 = t2;
    a31 = t1;
    a37 = t3;
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a02;
    sp_simd_t x1 = a08;
    sp_simd_t x5 = a14;
    sp_simd_t x6 = a20;
    sp_simd_t x3 = a26;
    sp_simd_t x2 = a32;
    sp_simd_t x4 = a38;

    p1 = sp_ntt_add_simd(x1, x4, p);
    p2 = sp_ntt_add_simd(x2, x5, p);
    p3 = sp_ntt_add_simd(x3, x6, p);
    p4 = sp_ntt_sub_simd(x1, x4, p);
    p5 = sp_ntt_sub_simd(x2, x5, p);
    p6 = sp_ntt_sub_simd(x3, x6, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t1 = sp_ntt_add_simd(t1, p3, p);
    t2 = sp_ntt_sub_simd(p1, p3, p);
    t3 = sp_ntt_sub_simd(p2, p3, p);
    t4 = sp_ntt_sub_simd(p4, p5, p);
    t4 = sp_ntt_add_partial_simd(t4, p6, p);
    t5 = sp_ntt_sub_simd(p4, p6, p);
    t6 = sp_ntt_add_simd(p5, p6, p);

    p0 = sp_ntt_add_partial_simd(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd(t5, t6, p);

    p0 = sp_ntt_mul_simd(p0, ntt_const[18], ntt_const[NC+18], p);
    p1 = sp_ntt_mul_simd(p1, ntt_const[19], ntt_const[NC+19], p);
    p2 = sp_ntt_mul_simd(p2, ntt_const[20], ntt_const[NC+20], p);
    p3 = sp_ntt_mul_simd(p3, ntt_const[21], ntt_const[NC+21], p);
    p4 = sp_ntt_mul_simd(p4, ntt_const[22], ntt_const[NC+22], p);
    p5 = sp_ntt_mul_simd(p5, ntt_const[23], ntt_const[NC+23], p);
    p6 = sp_ntt_mul_simd(p6, ntt_const[24], ntt_const[NC+24], p);
    p7 = sp_ntt_mul_simd(p7, ntt_const[25], ntt_const[NC+25], p);
    p8 = sp_ntt_mul_simd(p8, ntt_const[26], ntt_const[NC+26], p);

    t1 = sp_ntt_add_simd(p0, p1, p);
    t2 = sp_ntt_add_simd(p2, p4, p);
    t3 = sp_ntt_add_simd(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd(p6, p8, p);
    t6 = sp_ntt_add_simd(p7, p8, p);

    p1 = sp_ntt_add_simd(t1, t2, p);
    p2 = sp_ntt_add_simd(t1, t3, p);
    p3 = sp_ntt_sub_simd(t1, t2, p);
    p3 = sp_ntt_sub_simd(p3, t3, p);
    p4 = sp_ntt_add_simd(t4, t5, p);
    p5 = sp_ntt_sub_simd(t6, t4, p);
    p6 = sp_ntt_sub_simd(t4, t5, p);
    p6 = sp_ntt_add_simd(p6, t6, p);

    t1 = sp_ntt_add_simd(p1, p4, p);
    t2 = sp_ntt_add_simd(p2, p5, p);
    t3 = sp_ntt_add_simd(p3, p6, p);
    t4 = sp_ntt_sub_simd(p1, p4, p);
    t5 = sp_ntt_sub_simd(p2, p5, p);
    t6 = sp_ntt_sub_simd(p3, p6, p);

    a02 = p0;
    a08 = t6;
    a14 = t4;
    a20 = t5;
    a26 = t2;
    a32 = t1;
    a38 = t3;
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a03;
    sp_simd_t x1 = a09;
    sp_simd_t x5 = a15;
    sp_simd_t x6 = a21;
    sp_simd_t x3 = a27;
    sp_simd_t x2 = a33;
    sp_simd_t x4 = a39;

    p1 = sp_ntt_add_simd(x1, x4, p);
    p2 = sp_ntt_add_simd(x2, x5, p);
    p3 = sp_ntt_add_simd(x3, x6, p);
    p4 = sp_ntt_sub_simd(x1, x4, p);
    p5 = sp_ntt_sub_simd(x2, x5, p);
    p6 = sp_ntt_sub_simd(x3, x6, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t1 = sp_ntt_add_simd(t1, p3, p);
    t2 = sp_ntt_sub_simd(p1, p3, p);
    t3 = sp_ntt_sub_simd(p2, p3, p);
    t4 = sp_ntt_sub_simd(p4, p5, p);
    t4 = sp_ntt_add_partial_simd(t4, p6, p);
    t5 = sp_ntt_sub_simd(p4, p6, p);
    t6 = sp_ntt_add_simd(p5, p6, p);

    p0 = sp_ntt_add_partial_simd(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd(t5, t6, p);

    p0 = sp_ntt_mul_simd(p0, ntt_const[27], ntt_const[NC+27], p);
    p1 = sp_ntt_mul_simd(p1, ntt_const[28], ntt_const[NC+28], p);
    p2 = sp_ntt_mul_simd(p2, ntt_const[29], ntt_const[NC+29], p);
    p3 = sp_ntt_mul_simd(p3, ntt_const[30], ntt_const[NC+30], p);
    p4 = sp_ntt_mul_simd(p4, ntt_const[31], ntt_const[NC+31], p);
    p5 = sp_ntt_mul_simd(p5, ntt_const[32], ntt_const[NC+32], p);
    p6 = sp_ntt_mul_simd(p6, ntt_const[33], ntt_const[NC+33], p);
    p7 = sp_ntt_mul_simd(p7, ntt_const[34], ntt_const[NC+34], p);
    p8 = sp_ntt_mul_simd(p8, ntt_const[35], ntt_const[NC+35], p);

    t1 = sp_ntt_add_simd(p0, p1, p);
    t2 = sp_ntt_add_simd(p2, p4, p);
    t3 = sp_ntt_add_simd(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd(p6, p8, p);
    t6 = sp_ntt_add_simd(p7, p8, p);

    p1 = sp_ntt_add_simd(t1, t2, p);
    p2 = sp_ntt_add_simd(t1, t3, p);
    p3 = sp_ntt_sub_simd(t1, t2, p);
    p3 = sp_ntt_sub_simd(p3, t3, p);
    p4 = sp_ntt_add_simd(t4, t5, p);
    p5 = sp_ntt_sub_simd(t6, t4, p);
    p6 = sp_ntt_sub_simd(t4, t5, p);
    p6 = sp_ntt_add_simd(p6, t6, p);

    t1 = sp_ntt_add_simd(p1, p4, p);
    t2 = sp_ntt_add_simd(p2, p5, p);
    t3 = sp_ntt_add_simd(p3, p6, p);
    t4 = sp_ntt_sub_simd(p1, p4, p);
    t5 = sp_ntt_sub_simd(p2, p5, p);
    t6 = sp_ntt_sub_simd(p3, p6, p);

    a03 = p0;
    a09 = t6;
    a15 = t4;
    a21 = t5;
    a27 = t2;
    a33 = t1;
    a39 = t3;
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a04;
    sp_simd_t x1 = a10;
    sp_simd_t x5 = a16;
    sp_simd_t x6 = a22;
    sp_simd_t x3 = a28;
    sp_simd_t x2 = a34;
    sp_simd_t x4 = a40;

    p1 = sp_ntt_add_simd(x1, x4, p);
    p2 = sp_ntt_add_simd(x2, x5, p);
    p3 = sp_ntt_add_simd(x3, x6, p);
    p4 = sp_ntt_sub_simd(x1, x4, p);
    p5 = sp_ntt_sub_simd(x2, x5, p);
    p6 = sp_ntt_sub_simd(x3, x6, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t1 = sp_ntt_add_simd(t1, p3, p);
    t2 = sp_ntt_sub_simd(p1, p3, p);
    t3 = sp_ntt_sub_simd(p2, p3, p);
    t4 = sp_ntt_sub_simd(p4, p5, p);
    t4 = sp_ntt_add_partial_simd(t4, p6, p);
    t5 = sp_ntt_sub_simd(p4, p6, p);
    t6 = sp_ntt_add_simd(p5, p6, p);

    p0 = sp_ntt_add_partial_simd(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd(t5, t6, p);

    p0 = sp_ntt_mul_simd(p0, ntt_const[36], ntt_const[NC+36], p);
    p1 = sp_ntt_mul_simd(p1, ntt_const[37], ntt_const[NC+37], p);
    p2 = sp_ntt_mul_simd(p2, ntt_const[38], ntt_const[NC+38], p);
    p3 = sp_ntt_mul_simd(p3, ntt_const[39], ntt_const[NC+39], p);
    p4 = sp_ntt_mul_simd(p4, ntt_const[40], ntt_const[NC+40], p);
    p5 = sp_ntt_mul_simd(p5, ntt_const[41], ntt_const[NC+41], p);
    p6 = sp_ntt_mul_simd(p6, ntt_const[42], ntt_const[NC+42], p);
    p7 = sp_ntt_mul_simd(p7, ntt_const[43], ntt_const[NC+43], p);
    p8 = sp_ntt_mul_simd(p8, ntt_const[44], ntt_const[NC+44], p);

    t1 = sp_ntt_add_simd(p0, p1, p);
    t2 = sp_ntt_add_simd(p2, p4, p);
    t3 = sp_ntt_add_simd(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd(p6, p8, p);
    t6 = sp_ntt_add_simd(p7, p8, p);

    p1 = sp_ntt_add_simd(t1, t2, p);
    p2 = sp_ntt_add_simd(t1, t3, p);
    p3 = sp_ntt_sub_simd(t1, t2, p);
    p3 = sp_ntt_sub_simd(p3, t3, p);
    p4 = sp_ntt_add_simd(t4, t5, p);
    p5 = sp_ntt_sub_simd(t6, t4, p);
    p6 = sp_ntt_sub_simd(t4, t5, p);
    p6 = sp_ntt_add_simd(p6, t6, p);

    t1 = sp_ntt_add_simd(p1, p4, p);
    t2 = sp_ntt_add_simd(p2, p5, p);
    t3 = sp_ntt_add_simd(p3, p6, p);
    t4 = sp_ntt_sub_simd(p1, p4, p);
    t5 = sp_ntt_sub_simd(p2, p5, p);
    t6 = sp_ntt_sub_simd(p3, p6, p);

    a04 = p0;
    a10 = t6;
    a16 = t4;
    a22 = t5;
    a28 = t2;
    a34 = t1;
    a40 = t3;
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a05;
    sp_simd_t x1 = a11;
    sp_simd_t x5 = a17;
    sp_simd_t x6 = a23;
    sp_simd_t x3 = a29;
    sp_simd_t x2 = a35;
    sp_simd_t x4 = a41;

    p1 = sp_ntt_add_simd(x1, x4, p);
    p2 = sp_ntt_add_simd(x2, x5, p);
    p3 = sp_ntt_add_simd(x3, x6, p);
    p4 = sp_ntt_sub_simd(x1, x4, p);
    p5 = sp_ntt_sub_simd(x2, x5, p);
    p6 = sp_ntt_sub_simd(x3, x6, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t1 = sp_ntt_add_simd(t1, p3, p);
    t2 = sp_ntt_sub_simd(p1, p3, p);
    t3 = sp_ntt_sub_simd(p2, p3, p);
    t4 = sp_ntt_sub_simd(p4, p5, p);
    t4 = sp_ntt_add_partial_simd(t4, p6, p);
    t5 = sp_ntt_sub_simd(p4, p6, p);
    t6 = sp_ntt_add_simd(p5, p6, p);

    p0 = sp_ntt_add_partial_simd(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd(t5, t6, p);

    p0 = sp_ntt_mul_simd(p0, ntt_const[45], ntt_const[NC+45], p);
    p1 = sp_ntt_mul_simd(p1, ntt_const[46], ntt_const[NC+46], p);
    p2 = sp_ntt_mul_simd(p2, ntt_const[47], ntt_const[NC+47], p);
    p3 = sp_ntt_mul_simd(p3, ntt_const[48], ntt_const[NC+48], p);
    p4 = sp_ntt_mul_simd(p4, ntt_const[49], ntt_const[NC+49], p);
    p5 = sp_ntt_mul_simd(p5, ntt_const[50], ntt_const[NC+50], p);
    p6 = sp_ntt_mul_simd(p6, ntt_const[51], ntt_const[NC+51], p);
    p7 = sp_ntt_mul_simd(p7, ntt_const[52], ntt_const[NC+52], p);
    p8 = sp_ntt_mul_simd(p8, ntt_const[53], ntt_const[NC+53], p);

    t1 = sp_ntt_add_simd(p0, p1, p);
    t2 = sp_ntt_add_simd(p2, p4, p);
    t3 = sp_ntt_add_simd(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd(p6, p8, p);
    t6 = sp_ntt_add_simd(p7, p8, p);

    p1 = sp_ntt_add_simd(t1, t2, p);
    p2 = sp_ntt_add_simd(t1, t3, p);
    p3 = sp_ntt_sub_simd(t1, t2, p);
    p3 = sp_ntt_sub_simd(p3, t3, p);
    p4 = sp_ntt_add_simd(t4, t5, p);
    p5 = sp_ntt_sub_simd(t6, t4, p);
    p6 = sp_ntt_sub_simd(t4, t5, p);
    p6 = sp_ntt_add_simd(p6, t6, p);

    t1 = sp_ntt_add_simd(p1, p4, p);
    t2 = sp_ntt_add_simd(p2, p5, p);
    t3 = sp_ntt_add_simd(p3, p6, p);
    t4 = sp_ntt_sub_simd(p1, p4, p);
    t5 = sp_ntt_sub_simd(p2, p5, p);
    t6 = sp_ntt_sub_simd(p3, p6, p);

    a05 = p0;
    a11 = t6;
    a17 = t4;
    a23 = t5;
    a29 = t2;
    a35 = t1;
    a41 = t3;
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a00;
    sp_simd_t p1 = a01;
    sp_simd_t p2 = a02;
    sp_simd_t p3 = a03;
    sp_simd_t p4 = a04;
    sp_simd_t p5 = a05;

    p1 = sp_ntt_add_simd(p0, p1, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t2 = sp_ntt_sub_simd(p1, p2, p);
    t3 = sp_ntt_add_simd(p3, p5, p);
    t4 = sp_ntt_add_simd(p4, p5, p);

    p1 = sp_ntt_add_simd(t1, t3, p);
    p2 = sp_ntt_add_simd(t2, t4, p);
    p3 = sp_ntt_sub_simd(t1, t3, p);
    p4 = sp_ntt_sub_simd(t2, t4, p);

    sp_simd_pfa_scatter(p0, x, j00, inc2, n, vsize);
    sp_simd_pfa_scatter(p3, x, j07, inc2, n, vsize);
    sp_simd_pfa_scatter(p2, x, j14, inc2, n, vsize);
    sp_simd_pfa_scatter(p4, x, j21, inc2, n, vsize);
    sp_simd_pfa_scatter(p1, x, j28, inc2, n, vsize);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a06;
    sp_simd_t p1 = a07;
    sp_simd_t p2 = a08;
    sp_simd_t p3 = a09;
    sp_simd_t p4 = a10;
    sp_simd_t p5 = a11;

    p1 = sp_ntt_add_simd(p0, p1, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t2 = sp_ntt_sub_simd(p1, p2, p);
    t3 = sp_ntt_add_simd(p3, p5, p);
    t4 = sp_ntt_add_simd(p4, p5, p);

    p1 = sp_ntt_add_simd(t1, t3, p);
    p2 = sp_ntt_add_simd(t2, t4, p);
    p3 = sp_ntt_sub_simd(t1, t3, p);
    p4 = sp_ntt_sub_simd(t2, t4, p);

    sp_simd_pfa_scatter(p4, x, j01, inc2, n, vsize);
    sp_simd_pfa_scatter(p1, x, j08, inc2, n, vsize);
    sp_simd_pfa_scatter(p0, x, j15, inc2, n, vsize);
    sp_simd_pfa_scatter(p3, x, j22, inc2, n, vsize);
    sp_simd_pfa_scatter(p2, x, j29, inc2, n, vsize);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a12;
    sp_simd_t p1 = a13;
    sp_simd_t p2 = a14;
    sp_simd_t p3 = a15;
    sp_simd_t p4 = a16;
    sp_simd_t p5 = a17;

    p1 = sp_ntt_add_simd(p0, p1, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t2 = sp_ntt_sub_simd(p1, p2, p);
    t3 = sp_ntt_add_simd(p3, p5, p);
    t4 = sp_ntt_add_simd(p4, p5, p);

    p1 = sp_ntt_add_simd(t1, t3, p);
    p2 = sp_ntt_add_simd(t2, t4, p);
    p3 = sp_ntt_sub_simd(t1, t3, p);
    p4 = sp_ntt_sub_simd(t2, t4, p);

    sp_simd_pfa_scatter(p3, x, j02, inc2, n, vsize);
    sp_simd_pfa_scatter(p2, x, j09, inc2, n, vsize);
    sp_simd_pfa_scatter(p4, x, j16, inc2, n, vsize);
    sp_simd_pfa_scatter(p1, x, j23, inc2, n, vsize);
    sp_simd_pfa_scatter(p0, x, j30, inc2, n, vsize);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a18;
    sp_simd_t p1 = a19;
    sp_simd_t p2 = a20;
    sp_simd_t p3 = a21;
    sp_simd_t p4 = a22;
    sp_simd_t p5 = a23;

    p1 = sp_ntt_add_simd(p0, p1, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t2 = sp_ntt_sub_simd(p1, p2, p);
    t3 = sp_ntt_add_simd(p3, p5, p);
    t4 = sp_ntt_add_simd(p4, p5, p);

    p1 = sp_ntt_add_simd(t1, t3, p);
    p2 = sp_ntt_add_simd(t2, t4, p);
    p3 = sp_ntt_sub_simd(t1, t3, p);
    p4 = sp_ntt_sub_simd(t2, t4, p);

    sp_simd_pfa_scatter(p1, x, j03, inc2, n, vsize);
    sp_simd_pfa_scatter(p0, x, j10, inc2, n, vsize);
    sp_simd_pfa_scatter(p3, x, j17, inc2, n, vsize);
    sp_simd_pfa_scatter(p2, x, j24, inc2, n, vsize);
    sp_simd_pfa_scatter(p4, x, j31, inc2, n, vsize);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a24;
    sp_simd_t p1 = a25;
    sp_simd_t p2 = a26;
    sp_simd_t p3 = a27;
    sp_simd_t p4 = a28;
    sp_simd_t p5 = a29;

    p1 = sp_ntt_add_simd(p0, p1, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t2 = sp_ntt_sub_simd(p1, p2, p);
    t3 = sp_ntt_add_simd(p3, p5, p);
    t4 = sp_ntt_add_simd(p4, p5, p);

    p1 = sp_ntt_add_simd(t1, t3, p);
    p2 = sp_ntt_add_simd(t2, t4, p);
    p3 = sp_ntt_sub_simd(t1, t3, p);
    p4 = sp_ntt_sub_simd(t2, t4, p);

    sp_simd_pfa_scatter(p2, x, j04, inc2, n, vsize);
    sp_simd_pfa_scatter(p4, x, j11, inc2, n, vsize);
    sp_simd_pfa_scatter(p1, x, j18, inc2, n, vsize);
    sp_simd_pfa_scatter(p0, x, j25, inc2, n, vsize);
    sp_simd_pfa_scatter(p3, x, j32, inc2, n, vsize);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a30;
    sp_simd_t p1 = a31;
    sp_simd_t p2 = a32;
    sp_simd_t p3 = a33;
    sp_simd_t p4 = a34;
    sp_simd_t p5 = a35;

    p1 = sp_ntt_add_simd(p0, p1, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t2 = sp_ntt_sub_simd(p1, p2, p);
    t3 = sp_ntt_add_simd(p3, p5, p);
    t4 = sp_ntt_add_simd(p4, p5, p);

    p1 = sp_ntt_add_simd(t1, t3, p);
    p2 = sp_ntt_add_simd(t2, t4, p);
    p3 = sp_ntt_sub_simd(t1, t3, p);
    p4 = sp_ntt_sub_simd(t2, t4, p);

    sp_simd_pfa_scatter(p0, x, j05, inc2, n, vsize);
    sp_simd_pfa_scatter(p3, x, j12, inc2, n, vsize);
    sp_simd_pfa_scatter(p2, x, j19, inc2, n, vsize);
    sp_simd_pfa_scatter(p4, x, j26, inc2, n, vsize);
    sp_simd_pfa_scatter(p1, x, j33, inc2, n, vsize);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a36;
    sp_simd_t p1 = a37;
    sp_simd_t p2 = a38;
    sp_simd_t p3 = a39;
    sp_simd_t p4 = a40;
    sp_simd_t p5 = a41;

    p1 = sp_ntt_add_simd(p0, p1, p);

    t1 = sp_ntt_add_simd(p1, p2, p);
    t2 = sp_ntt_sub_simd(p1, p2, p);
    t3 = sp_ntt_add_simd(p3, p5, p);
    t4 = sp_ntt_add_simd(p4, p5, p);

    p1 = sp_ntt_add_simd(t1, t3, p);
    p2 = sp_ntt_add_simd(t2, t4, p);
    p3 = sp_ntt_sub_simd(t1, t3, p);
    p4 = sp_ntt_sub_simd(t2, t4, p);

    sp_simd_pfa_scatter(p4, x, j06, inc2, n, vsize);
    sp_simd_pfa_scatter(p1, x, j13, inc2, n, vsize);
    sp_simd_pfa_scatter(p0, x, j20, inc2, n, vsize);
    sp_simd_pfa_scatter(p3, x, j27, inc2, n, vsize);
    sp_simd_pfa_scatter(p2, x, j34, inc2, n, vsize);
  }
}

static void
ntt35_pfa_run_core_simd_interleaved(
	  spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t n,
	  sp_simd_t p, sp_simd_t * c)
{
  spv_size_t j00, j01, j02, j03, j04, j05, j06, j07,
	     j08, j09, j10, j11, j12, j13, j14, j15,
	     j16, j17, j18, j19, j20, j21, j22, j23,
	     j24, j25, j26, j27, j28, j29, j30, j31,
	     j32, j33, j34;

  sp_simd_t a00, a01, a02, a03, a04, a05, a06, 
            a07, a08, a09, a10, a11, a12, a13, 
            a14, a15, a16, a17, a18, a19, a20, 
            a21, a22, a23, a24, a25, a26, a27, 
            a28, a29, a30, a31, a32, a33, a34, 
            a35, a36, a37, a38, a39, a40, a41;

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

  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x0 = sp_simd_load(x + j00);
    sp_simd_t x1 = sp_simd_load(x + j07);
    sp_simd_t x4 = sp_simd_load(x + j14);
    sp_simd_t x2 = sp_simd_load(x + j21);
    sp_simd_t x3 = sp_simd_load(x + j28);

    t1 = sp_ntt_add_simd0(x1, x3, p);
    t3 = sp_ntt_sub_simd0(x1, x3, p);
    t2 = sp_ntt_add_simd0(x2, x4, p);
    t4 = sp_ntt_sub_simd0(x2, x4, p);

    a01 = sp_ntt_add_simd0(t1, t2, p);
    a02 = sp_ntt_sub_simd0(t1, t2, p);
    a03 = t3;
    a04 = t4;
    a05 = sp_ntt_add_simd0(t3, t4, p);

    a00 = sp_ntt_add_simd0(x0, a01, p);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x0 = sp_simd_load(x + j05);
    sp_simd_t x1 = sp_simd_load(x + j12);
    sp_simd_t x4 = sp_simd_load(x + j19);
    sp_simd_t x2 = sp_simd_load(x + j26);
    sp_simd_t x3 = sp_simd_load(x + j33);

    t1 = sp_ntt_add_simd0(x1, x3, p);
    t3 = sp_ntt_sub_simd0(x1, x3, p);
    t2 = sp_ntt_add_simd0(x2, x4, p);
    t4 = sp_ntt_sub_simd0(x2, x4, p);

    a07 = sp_ntt_add_simd0(t1, t2, p);
    a08 = sp_ntt_sub_simd0(t1, t2, p);
    a09 = t3;
    a10 = t4;
    a11 = sp_ntt_add_simd0(t3, t4, p);

    a06 = sp_ntt_add_simd0(x0, a07, p);
  }
  {
    sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7;

    sp_simd_t x3 = sp_simd_load(x + j03);
    sp_simd_t x0 = sp_simd_load(x + j10);
    sp_simd_t x1 = sp_simd_load(x + j17);
    sp_simd_t x4 = sp_simd_load(x + j24);
    sp_simd_t x2 = sp_simd_load(x + j31);

    t1 = sp_ntt_add_simd0(x1, x3, p);
    t3 = sp_ntt_sub_simd0(x1, x3, p);
    t2 = sp_ntt_add_simd0(x2, x4, p);
    t4 = sp_ntt_sub_simd0(x2, x4, p);

    a13 = sp_ntt_add_simd0(t1, t2, p);
    a14 = sp_ntt_sub_simd0(t1, t2, p);
    a15 = t3;
    a16 = t4;
    a17 = sp_ntt_add_simd0(t3, t4, p);

    a12 = sp_ntt_add_simd0(x0, a13, p);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x2 = sp_simd_load(x + j01);
    sp_simd_t x3 = sp_simd_load(x + j08);
    sp_simd_t x0 = sp_simd_load(x + j15);
    sp_simd_t x1 = sp_simd_load(x + j22);
    sp_simd_t x4 = sp_simd_load(x + j29);

    t1 = sp_ntt_add_simd0(x1, x3, p);
    t3 = sp_ntt_sub_simd0(x1, x3, p);
    t2 = sp_ntt_add_simd0(x2, x4, p);
    t4 = sp_ntt_sub_simd0(x2, x4, p);

    a19 = sp_ntt_add_simd0(t1, t2, p);
    a20 = sp_ntt_sub_simd0(t1, t2, p);
    a21 = t3;
    a22 = t4;
    a23 = sp_ntt_add_simd0(t3, t4, p);

    a18 = sp_ntt_add_simd0(x0, a19, p);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x2 = sp_simd_load(x + j06);
    sp_simd_t x3 = sp_simd_load(x + j13);
    sp_simd_t x0 = sp_simd_load(x + j20);
    sp_simd_t x1 = sp_simd_load(x + j27);
    sp_simd_t x4 = sp_simd_load(x + j34);

    t1 = sp_ntt_add_simd0(x1, x3, p);
    t3 = sp_ntt_sub_simd0(x1, x3, p);
    t2 = sp_ntt_add_simd0(x2, x4, p);
    t4 = sp_ntt_sub_simd0(x2, x4, p);

    a25 = sp_ntt_add_simd0(t1, t2, p);
    a26 = sp_ntt_sub_simd0(t1, t2, p);
    a27 = t3;
    a28 = t4;
    a29 = sp_ntt_add_simd0(t3, t4, p);

    a24 = sp_ntt_add_simd0(x0, a25, p);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x4 = sp_simd_load(x + j04);
    sp_simd_t x2 = sp_simd_load(x + j11);
    sp_simd_t x3 = sp_simd_load(x + j18);
    sp_simd_t x0 = sp_simd_load(x + j25);
    sp_simd_t x1 = sp_simd_load(x + j32);

    t1 = sp_ntt_add_simd0(x1, x3, p);
    t3 = sp_ntt_sub_simd0(x1, x3, p);
    t2 = sp_ntt_add_simd0(x2, x4, p);
    t4 = sp_ntt_sub_simd0(x2, x4, p);

    a31 = sp_ntt_add_simd0(t1, t2, p);
    a32 = sp_ntt_sub_simd0(t1, t2, p);
    a33 = t3;
    a34 = t4;
    a35 = sp_ntt_add_simd0(t3, t4, p);

    a30 = sp_ntt_add_simd0(x0, a31, p);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t x1 = sp_simd_load(x + j02);
    sp_simd_t x4 = sp_simd_load(x + j09);
    sp_simd_t x2 = sp_simd_load(x + j16);
    sp_simd_t x3 = sp_simd_load(x + j23);
    sp_simd_t x0 = sp_simd_load(x + j30);

    t1 = sp_ntt_add_simd0(x1, x3, p);
    t3 = sp_ntt_sub_simd0(x1, x3, p);
    t2 = sp_ntt_add_simd0(x2, x4, p);
    t4 = sp_ntt_sub_simd0(x2, x4, p);

    a37 = sp_ntt_add_simd0(t1, t2, p);
    a38 = sp_ntt_sub_simd0(t1, t2, p);
    a39 = t3;
    a40 = t4;
    a41 = sp_ntt_add_simd0(t3, t4, p);

    a36 = sp_ntt_add_simd0(x0, a37, p);
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a00;
    sp_simd_t x1 = a06;
    sp_simd_t x5 = a12;
    sp_simd_t x6 = a18;
    sp_simd_t x3 = a24;
    sp_simd_t x2 = a30;
    sp_simd_t x4 = a36;

    p1 = sp_ntt_add_simd0(x1, x4, p);
    p2 = sp_ntt_add_simd0(x2, x5, p);
    p3 = sp_ntt_add_simd0(x3, x6, p);
    p4 = sp_ntt_sub_simd0(x1, x4, p);
    p5 = sp_ntt_sub_simd0(x2, x5, p);
    p6 = sp_ntt_sub_simd0(x3, x6, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t1 = sp_ntt_add_simd0(t1, p3, p);
    t2 = sp_ntt_sub_simd0(p1, p3, p);
    t3 = sp_ntt_sub_simd0(p2, p3, p);
    t4 = sp_ntt_sub_simd0(p4, p5, p);
    t4 = sp_ntt_add_partial_simd0(t4, p6, p);
    t5 = sp_ntt_sub_simd0(p4, p6, p);
    t6 = sp_ntt_add_simd0(p5, p6, p);

    p0 = sp_ntt_add_simd0(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd0(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd0(t5, t6, p);

    p1 = sp_ntt_mul_simd0(p1, c+2*1, p);
    p2 = sp_ntt_mul_simd0(p2, c+2*2, p);
    p3 = sp_ntt_mul_simd0(p3, c+2*3, p);
    p4 = sp_ntt_mul_simd0(p4, c+2*4, p);
    p5 = sp_ntt_mul_simd0(p5, c+2*5, p);
    p6 = sp_ntt_mul_simd0(p6, c+2*6, p);
    p7 = sp_ntt_mul_simd0(p7, c+2*7, p);
    p8 = sp_ntt_mul_simd0(p8, c+2*8, p);

    t1 = sp_ntt_add_simd0(p0, p1, p);
    t2 = sp_ntt_add_simd0(p2, p4, p);
    t3 = sp_ntt_add_simd0(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd0(p6, p8, p);
    t6 = sp_ntt_add_simd0(p7, p8, p);

    p1 = sp_ntt_add_simd0(t1, t2, p);
    p2 = sp_ntt_add_simd0(t1, t3, p);
    p3 = sp_ntt_sub_simd0(t1, t2, p);
    p3 = sp_ntt_sub_simd0(p3, t3, p);
    p4 = sp_ntt_add_simd0(t4, t5, p);
    p5 = sp_ntt_sub_simd0(t6, t4, p);
    p6 = sp_ntt_sub_simd0(t4, t5, p);
    p6 = sp_ntt_add_simd0(p6, t6, p);

    t1 = sp_ntt_add_simd0(p1, p4, p);
    t2 = sp_ntt_add_simd0(p2, p5, p);
    t3 = sp_ntt_add_simd0(p3, p6, p);
    t4 = sp_ntt_sub_simd0(p1, p4, p);
    t5 = sp_ntt_sub_simd0(p2, p5, p);
    t6 = sp_ntt_sub_simd0(p3, p6, p);

    a00 = p0;
    a06 = t6;
    a12 = t4;
    a18 = t5;
    a24 = t2;
    a30 = t1;
    a36 = t3;
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a01;
    sp_simd_t x1 = a07;
    sp_simd_t x5 = a13;
    sp_simd_t x6 = a19;
    sp_simd_t x3 = a25;
    sp_simd_t x2 = a31;
    sp_simd_t x4 = a37;

    p1 = sp_ntt_add_simd0(x1, x4, p);
    p2 = sp_ntt_add_simd0(x2, x5, p);
    p3 = sp_ntt_add_simd0(x3, x6, p);
    p4 = sp_ntt_sub_simd0(x1, x4, p);
    p5 = sp_ntt_sub_simd0(x2, x5, p);
    p6 = sp_ntt_sub_simd0(x3, x6, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t1 = sp_ntt_add_simd0(t1, p3, p);
    t2 = sp_ntt_sub_simd0(p1, p3, p);
    t3 = sp_ntt_sub_simd0(p2, p3, p);
    t4 = sp_ntt_sub_simd0(p4, p5, p);
    t4 = sp_ntt_add_partial_simd0(t4, p6, p);
    t5 = sp_ntt_sub_simd0(p4, p6, p);
    t6 = sp_ntt_add_simd0(p5, p6, p);

    p0 = sp_ntt_add_partial_simd0(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd0(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd0(t5, t6, p);

    p0 = sp_ntt_mul_simd0(p0, c+2*9, p);
    p1 = sp_ntt_mul_simd0(p1, c+2*10, p);
    p2 = sp_ntt_mul_simd0(p2, c+2*11, p);
    p3 = sp_ntt_mul_simd0(p3, c+2*12, p);
    p4 = sp_ntt_mul_simd0(p4, c+2*13, p);
    p5 = sp_ntt_mul_simd0(p5, c+2*14, p);
    p6 = sp_ntt_mul_simd0(p6, c+2*15, p);
    p7 = sp_ntt_mul_simd0(p7, c+2*16, p);
    p8 = sp_ntt_mul_simd0(p8, c+2*17, p);

    t1 = sp_ntt_add_simd0(p0, p1, p);
    t2 = sp_ntt_add_simd0(p2, p4, p);
    t3 = sp_ntt_add_simd0(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd0(p6, p8, p);
    t6 = sp_ntt_add_simd0(p7, p8, p);

    p1 = sp_ntt_add_simd0(t1, t2, p);
    p2 = sp_ntt_add_simd0(t1, t3, p);
    p3 = sp_ntt_sub_simd0(t1, t2, p);
    p3 = sp_ntt_sub_simd0(p3, t3, p);
    p4 = sp_ntt_add_simd0(t4, t5, p);
    p5 = sp_ntt_sub_simd0(t6, t4, p);
    p6 = sp_ntt_sub_simd0(t4, t5, p);
    p6 = sp_ntt_add_simd0(p6, t6, p);

    t1 = sp_ntt_add_simd0(p1, p4, p);
    t2 = sp_ntt_add_simd0(p2, p5, p);
    t3 = sp_ntt_add_simd0(p3, p6, p);
    t4 = sp_ntt_sub_simd0(p1, p4, p);
    t5 = sp_ntt_sub_simd0(p2, p5, p);
    t6 = sp_ntt_sub_simd0(p3, p6, p);

    a01 = p0;
    a07 = t6;
    a13 = t4;
    a19 = t5;
    a25 = t2;
    a31 = t1;
    a37 = t3;
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a02;
    sp_simd_t x1 = a08;
    sp_simd_t x5 = a14;
    sp_simd_t x6 = a20;
    sp_simd_t x3 = a26;
    sp_simd_t x2 = a32;
    sp_simd_t x4 = a38;

    p1 = sp_ntt_add_simd0(x1, x4, p);
    p2 = sp_ntt_add_simd0(x2, x5, p);
    p3 = sp_ntt_add_simd0(x3, x6, p);
    p4 = sp_ntt_sub_simd0(x1, x4, p);
    p5 = sp_ntt_sub_simd0(x2, x5, p);
    p6 = sp_ntt_sub_simd0(x3, x6, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t1 = sp_ntt_add_simd0(t1, p3, p);
    t2 = sp_ntt_sub_simd0(p1, p3, p);
    t3 = sp_ntt_sub_simd0(p2, p3, p);
    t4 = sp_ntt_sub_simd0(p4, p5, p);
    t4 = sp_ntt_add_partial_simd0(t4, p6, p);
    t5 = sp_ntt_sub_simd0(p4, p6, p);
    t6 = sp_ntt_add_simd0(p5, p6, p);

    p0 = sp_ntt_add_partial_simd0(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd0(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd0(t5, t6, p);

    p0 = sp_ntt_mul_simd0(p0, c+2*18, p);
    p1 = sp_ntt_mul_simd0(p1, c+2*19, p);
    p2 = sp_ntt_mul_simd0(p2, c+2*20, p);
    p3 = sp_ntt_mul_simd0(p3, c+2*21, p);
    p4 = sp_ntt_mul_simd0(p4, c+2*22, p);
    p5 = sp_ntt_mul_simd0(p5, c+2*23, p);
    p6 = sp_ntt_mul_simd0(p6, c+2*24, p);
    p7 = sp_ntt_mul_simd0(p7, c+2*25, p);
    p8 = sp_ntt_mul_simd0(p8, c+2*26, p);

    t1 = sp_ntt_add_simd0(p0, p1, p);
    t2 = sp_ntt_add_simd0(p2, p4, p);
    t3 = sp_ntt_add_simd0(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd0(p6, p8, p);
    t6 = sp_ntt_add_simd0(p7, p8, p);

    p1 = sp_ntt_add_simd0(t1, t2, p);
    p2 = sp_ntt_add_simd0(t1, t3, p);
    p3 = sp_ntt_sub_simd0(t1, t2, p);
    p3 = sp_ntt_sub_simd0(p3, t3, p);
    p4 = sp_ntt_add_simd0(t4, t5, p);
    p5 = sp_ntt_sub_simd0(t6, t4, p);
    p6 = sp_ntt_sub_simd0(t4, t5, p);
    p6 = sp_ntt_add_simd0(p6, t6, p);

    t1 = sp_ntt_add_simd0(p1, p4, p);
    t2 = sp_ntt_add_simd0(p2, p5, p);
    t3 = sp_ntt_add_simd0(p3, p6, p);
    t4 = sp_ntt_sub_simd0(p1, p4, p);
    t5 = sp_ntt_sub_simd0(p2, p5, p);
    t6 = sp_ntt_sub_simd0(p3, p6, p);

    a02 = p0;
    a08 = t6;
    a14 = t4;
    a20 = t5;
    a26 = t2;
    a32 = t1;
    a38 = t3;
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a03;
    sp_simd_t x1 = a09;
    sp_simd_t x5 = a15;
    sp_simd_t x6 = a21;
    sp_simd_t x3 = a27;
    sp_simd_t x2 = a33;
    sp_simd_t x4 = a39;

    p1 = sp_ntt_add_simd0(x1, x4, p);
    p2 = sp_ntt_add_simd0(x2, x5, p);
    p3 = sp_ntt_add_simd0(x3, x6, p);
    p4 = sp_ntt_sub_simd0(x1, x4, p);
    p5 = sp_ntt_sub_simd0(x2, x5, p);
    p6 = sp_ntt_sub_simd0(x3, x6, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t1 = sp_ntt_add_simd0(t1, p3, p);
    t2 = sp_ntt_sub_simd0(p1, p3, p);
    t3 = sp_ntt_sub_simd0(p2, p3, p);
    t4 = sp_ntt_sub_simd0(p4, p5, p);
    t4 = sp_ntt_add_partial_simd0(t4, p6, p);
    t5 = sp_ntt_sub_simd0(p4, p6, p);
    t6 = sp_ntt_add_simd0(p5, p6, p);

    p0 = sp_ntt_add_partial_simd0(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd0(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd0(t5, t6, p);

    p0 = sp_ntt_mul_simd0(p0, c+2*27, p);
    p1 = sp_ntt_mul_simd0(p1, c+2*28, p);
    p2 = sp_ntt_mul_simd0(p2, c+2*29, p);
    p3 = sp_ntt_mul_simd0(p3, c+2*30, p);
    p4 = sp_ntt_mul_simd0(p4, c+2*31, p);
    p5 = sp_ntt_mul_simd0(p5, c+2*32, p);
    p6 = sp_ntt_mul_simd0(p6, c+2*33, p);
    p7 = sp_ntt_mul_simd0(p7, c+2*34, p);
    p8 = sp_ntt_mul_simd0(p8, c+2*35, p);

    t1 = sp_ntt_add_simd0(p0, p1, p);
    t2 = sp_ntt_add_simd0(p2, p4, p);
    t3 = sp_ntt_add_simd0(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd0(p6, p8, p);
    t6 = sp_ntt_add_simd0(p7, p8, p);

    p1 = sp_ntt_add_simd0(t1, t2, p);
    p2 = sp_ntt_add_simd0(t1, t3, p);
    p3 = sp_ntt_sub_simd0(t1, t2, p);
    p3 = sp_ntt_sub_simd0(p3, t3, p);
    p4 = sp_ntt_add_simd0(t4, t5, p);
    p5 = sp_ntt_sub_simd0(t6, t4, p);
    p6 = sp_ntt_sub_simd0(t4, t5, p);
    p6 = sp_ntt_add_simd0(p6, t6, p);

    t1 = sp_ntt_add_simd0(p1, p4, p);
    t2 = sp_ntt_add_simd0(p2, p5, p);
    t3 = sp_ntt_add_simd0(p3, p6, p);
    t4 = sp_ntt_sub_simd0(p1, p4, p);
    t5 = sp_ntt_sub_simd0(p2, p5, p);
    t6 = sp_ntt_sub_simd0(p3, p6, p);

    a03 = p0;
    a09 = t6;
    a15 = t4;
    a21 = t5;
    a27 = t2;
    a33 = t1;
    a39 = t3;
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a04;
    sp_simd_t x1 = a10;
    sp_simd_t x5 = a16;
    sp_simd_t x6 = a22;
    sp_simd_t x3 = a28;
    sp_simd_t x2 = a34;
    sp_simd_t x4 = a40;

    p1 = sp_ntt_add_simd0(x1, x4, p);
    p2 = sp_ntt_add_simd0(x2, x5, p);
    p3 = sp_ntt_add_simd0(x3, x6, p);
    p4 = sp_ntt_sub_simd0(x1, x4, p);
    p5 = sp_ntt_sub_simd0(x2, x5, p);
    p6 = sp_ntt_sub_simd0(x3, x6, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t1 = sp_ntt_add_simd0(t1, p3, p);
    t2 = sp_ntt_sub_simd0(p1, p3, p);
    t3 = sp_ntt_sub_simd0(p2, p3, p);
    t4 = sp_ntt_sub_simd0(p4, p5, p);
    t4 = sp_ntt_add_partial_simd0(t4, p6, p);
    t5 = sp_ntt_sub_simd0(p4, p6, p);
    t6 = sp_ntt_add_simd0(p5, p6, p);

    p0 = sp_ntt_add_partial_simd0(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd0(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd0(t5, t6, p);

    p0 = sp_ntt_mul_simd0(p0, c+2*36, p);
    p1 = sp_ntt_mul_simd0(p1, c+2*37, p);
    p2 = sp_ntt_mul_simd0(p2, c+2*38, p);
    p3 = sp_ntt_mul_simd0(p3, c+2*39, p);
    p4 = sp_ntt_mul_simd0(p4, c+2*40, p);
    p5 = sp_ntt_mul_simd0(p5, c+2*41, p);
    p6 = sp_ntt_mul_simd0(p6, c+2*42, p);
    p7 = sp_ntt_mul_simd0(p7, c+2*43, p);
    p8 = sp_ntt_mul_simd0(p8, c+2*44, p);

    t1 = sp_ntt_add_simd0(p0, p1, p);
    t2 = sp_ntt_add_simd0(p2, p4, p);
    t3 = sp_ntt_add_simd0(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd0(p6, p8, p);
    t6 = sp_ntt_add_simd0(p7, p8, p);

    p1 = sp_ntt_add_simd0(t1, t2, p);
    p2 = sp_ntt_add_simd0(t1, t3, p);
    p3 = sp_ntt_sub_simd0(t1, t2, p);
    p3 = sp_ntt_sub_simd0(p3, t3, p);
    p4 = sp_ntt_add_simd0(t4, t5, p);
    p5 = sp_ntt_sub_simd0(t6, t4, p);
    p6 = sp_ntt_sub_simd0(t4, t5, p);
    p6 = sp_ntt_add_simd0(p6, t6, p);

    t1 = sp_ntt_add_simd0(p1, p4, p);
    t2 = sp_ntt_add_simd0(p2, p5, p);
    t3 = sp_ntt_add_simd0(p3, p6, p);
    t4 = sp_ntt_sub_simd0(p1, p4, p);
    t5 = sp_ntt_sub_simd0(p2, p5, p);
    t6 = sp_ntt_sub_simd0(p3, p6, p);

    a04 = p0;
    a10 = t6;
    a16 = t4;
    a22 = t5;
    a28 = t2;
    a34 = t1;
    a40 = t3;
  }
  {
    sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
    sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

    sp_simd_t x0 = a05;
    sp_simd_t x1 = a11;
    sp_simd_t x5 = a17;
    sp_simd_t x6 = a23;
    sp_simd_t x3 = a29;
    sp_simd_t x2 = a35;
    sp_simd_t x4 = a41;

    p1 = sp_ntt_add_simd0(x1, x4, p);
    p2 = sp_ntt_add_simd0(x2, x5, p);
    p3 = sp_ntt_add_simd0(x3, x6, p);
    p4 = sp_ntt_sub_simd0(x1, x4, p);
    p5 = sp_ntt_sub_simd0(x2, x5, p);
    p6 = sp_ntt_sub_simd0(x3, x6, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t1 = sp_ntt_add_simd0(t1, p3, p);
    t2 = sp_ntt_sub_simd0(p1, p3, p);
    t3 = sp_ntt_sub_simd0(p2, p3, p);
    t4 = sp_ntt_sub_simd0(p4, p5, p);
    t4 = sp_ntt_add_partial_simd0(t4, p6, p);
    t5 = sp_ntt_sub_simd0(p4, p6, p);
    t6 = sp_ntt_add_simd0(p5, p6, p);

    p0 = sp_ntt_add_partial_simd0(x0, t1, p);
    p1 = t1;
    p2 = t2;
    p3 = t3;
    p4 = sp_ntt_add_partial_simd0(t2, t3, p);
    p5 = t4;
    p6 = t5;
    p7 = t6;
    p8 = sp_ntt_add_partial_simd0(t5, t6, p);

    p0 = sp_ntt_mul_simd0(p0, c+2*45, p);
    p1 = sp_ntt_mul_simd0(p1, c+2*46, p);
    p2 = sp_ntt_mul_simd0(p2, c+2*47, p);
    p3 = sp_ntt_mul_simd0(p3, c+2*48, p);
    p4 = sp_ntt_mul_simd0(p4, c+2*49, p);
    p5 = sp_ntt_mul_simd0(p5, c+2*50, p);
    p6 = sp_ntt_mul_simd0(p6, c+2*51, p);
    p7 = sp_ntt_mul_simd0(p7, c+2*52, p);
    p8 = sp_ntt_mul_simd0(p8, c+2*53, p);

    t1 = sp_ntt_add_simd0(p0, p1, p);
    t2 = sp_ntt_add_simd0(p2, p4, p);
    t3 = sp_ntt_add_simd0(p3, p4, p);
    t4 = p5;
    t5 = sp_ntt_add_simd0(p6, p8, p);
    t6 = sp_ntt_add_simd0(p7, p8, p);

    p1 = sp_ntt_add_simd0(t1, t2, p);
    p2 = sp_ntt_add_simd0(t1, t3, p);
    p3 = sp_ntt_sub_simd0(t1, t2, p);
    p3 = sp_ntt_sub_simd0(p3, t3, p);
    p4 = sp_ntt_add_simd0(t4, t5, p);
    p5 = sp_ntt_sub_simd0(t6, t4, p);
    p6 = sp_ntt_sub_simd0(t4, t5, p);
    p6 = sp_ntt_add_simd0(p6, t6, p);

    t1 = sp_ntt_add_simd0(p1, p4, p);
    t2 = sp_ntt_add_simd0(p2, p5, p);
    t3 = sp_ntt_add_simd0(p3, p6, p);
    t4 = sp_ntt_sub_simd0(p1, p4, p);
    t5 = sp_ntt_sub_simd0(p2, p5, p);
    t6 = sp_ntt_sub_simd0(p3, p6, p);

    a05 = p0;
    a11 = t6;
    a17 = t4;
    a23 = t5;
    a29 = t2;
    a35 = t1;
    a41 = t3;
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a00;
    sp_simd_t p1 = a01;
    sp_simd_t p2 = a02;
    sp_simd_t p3 = a03;
    sp_simd_t p4 = a04;
    sp_simd_t p5 = a05;

    p1 = sp_ntt_add_simd0(p0, p1, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t2 = sp_ntt_sub_simd0(p1, p2, p);
    t3 = sp_ntt_add_simd0(p3, p5, p);
    t4 = sp_ntt_add_simd0(p4, p5, p);

    p1 = sp_ntt_add_simd0(t1, t3, p);
    p2 = sp_ntt_add_simd0(t2, t4, p);
    p3 = sp_ntt_sub_simd0(t1, t3, p);
    p4 = sp_ntt_sub_simd0(t2, t4, p);

    sp_simd_store(p0, x + j00);
    sp_simd_store(p3, x + j07);
    sp_simd_store(p2, x + j14);
    sp_simd_store(p4, x + j21);
    sp_simd_store(p1, x + j28);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a06;
    sp_simd_t p1 = a07;
    sp_simd_t p2 = a08;
    sp_simd_t p3 = a09;
    sp_simd_t p4 = a10;
    sp_simd_t p5 = a11;

    p1 = sp_ntt_add_simd0(p0, p1, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t2 = sp_ntt_sub_simd0(p1, p2, p);
    t3 = sp_ntt_add_simd0(p3, p5, p);
    t4 = sp_ntt_add_simd0(p4, p5, p);

    p1 = sp_ntt_add_simd0(t1, t3, p);
    p2 = sp_ntt_add_simd0(t2, t4, p);
    p3 = sp_ntt_sub_simd0(t1, t3, p);
    p4 = sp_ntt_sub_simd0(t2, t4, p);

    sp_simd_store(p4, x + j01);
    sp_simd_store(p1, x + j08);
    sp_simd_store(p0, x + j15);
    sp_simd_store(p3, x + j22);
    sp_simd_store(p2, x + j29);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a12;
    sp_simd_t p1 = a13;
    sp_simd_t p2 = a14;
    sp_simd_t p3 = a15;
    sp_simd_t p4 = a16;
    sp_simd_t p5 = a17;

    p1 = sp_ntt_add_simd0(p0, p1, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t2 = sp_ntt_sub_simd0(p1, p2, p);
    t3 = sp_ntt_add_simd0(p3, p5, p);
    t4 = sp_ntt_add_simd0(p4, p5, p);

    p1 = sp_ntt_add_simd0(t1, t3, p);
    p2 = sp_ntt_add_simd0(t2, t4, p);
    p3 = sp_ntt_sub_simd0(t1, t3, p);
    p4 = sp_ntt_sub_simd0(t2, t4, p);

    sp_simd_store(p3, x + j02);
    sp_simd_store(p2, x + j09);
    sp_simd_store(p4, x + j16);
    sp_simd_store(p1, x + j23);
    sp_simd_store(p0, x + j30);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a18;
    sp_simd_t p1 = a19;
    sp_simd_t p2 = a20;
    sp_simd_t p3 = a21;
    sp_simd_t p4 = a22;
    sp_simd_t p5 = a23;

    p1 = sp_ntt_add_simd0(p0, p1, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t2 = sp_ntt_sub_simd0(p1, p2, p);
    t3 = sp_ntt_add_simd0(p3, p5, p);
    t4 = sp_ntt_add_simd0(p4, p5, p);

    p1 = sp_ntt_add_simd0(t1, t3, p);
    p2 = sp_ntt_add_simd0(t2, t4, p);
    p3 = sp_ntt_sub_simd0(t1, t3, p);
    p4 = sp_ntt_sub_simd0(t2, t4, p);

    sp_simd_store(p1, x + j03);
    sp_simd_store(p0, x + j10);
    sp_simd_store(p3, x + j17);
    sp_simd_store(p2, x + j24);
    sp_simd_store(p4, x + j31);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a24;
    sp_simd_t p1 = a25;
    sp_simd_t p2 = a26;
    sp_simd_t p3 = a27;
    sp_simd_t p4 = a28;
    sp_simd_t p5 = a29;

    p1 = sp_ntt_add_simd0(p0, p1, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t2 = sp_ntt_sub_simd0(p1, p2, p);
    t3 = sp_ntt_add_simd0(p3, p5, p);
    t4 = sp_ntt_add_simd0(p4, p5, p);

    p1 = sp_ntt_add_simd0(t1, t3, p);
    p2 = sp_ntt_add_simd0(t2, t4, p);
    p3 = sp_ntt_sub_simd0(t1, t3, p);
    p4 = sp_ntt_sub_simd0(t2, t4, p);

    sp_simd_store(p2, x + j04);
    sp_simd_store(p4, x + j11);
    sp_simd_store(p1, x + j18);
    sp_simd_store(p0, x + j25);
    sp_simd_store(p3, x + j32);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a30;
    sp_simd_t p1 = a31;
    sp_simd_t p2 = a32;
    sp_simd_t p3 = a33;
    sp_simd_t p4 = a34;
    sp_simd_t p5 = a35;

    p1 = sp_ntt_add_simd0(p0, p1, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t2 = sp_ntt_sub_simd0(p1, p2, p);
    t3 = sp_ntt_add_simd0(p3, p5, p);
    t4 = sp_ntt_add_simd0(p4, p5, p);

    p1 = sp_ntt_add_simd0(t1, t3, p);
    p2 = sp_ntt_add_simd0(t2, t4, p);
    p3 = sp_ntt_sub_simd0(t1, t3, p);
    p4 = sp_ntt_sub_simd0(t2, t4, p);

    sp_simd_store(p0, x + j05);
    sp_simd_store(p3, x + j12);
    sp_simd_store(p2, x + j19);
    sp_simd_store(p4, x + j26);
    sp_simd_store(p1, x + j33);
  }
  {
    sp_simd_t t1, t2, t3, t4;

    sp_simd_t p0 = a36;
    sp_simd_t p1 = a37;
    sp_simd_t p2 = a38;
    sp_simd_t p3 = a39;
    sp_simd_t p4 = a40;
    sp_simd_t p5 = a41;

    p1 = sp_ntt_add_simd0(p0, p1, p);

    t1 = sp_ntt_add_simd0(p1, p2, p);
    t2 = sp_ntt_sub_simd0(p1, p2, p);
    t3 = sp_ntt_add_simd0(p3, p5, p);
    t4 = sp_ntt_add_simd0(p4, p5, p);

    p1 = sp_ntt_add_simd0(t1, t3, p);
    p2 = sp_ntt_add_simd0(t2, t4, p);
    p3 = sp_ntt_sub_simd0(t1, t3, p);
    p4 = sp_ntt_sub_simd0(t2, t4, p);

    sp_simd_store(p4, x + j06);
    sp_simd_store(p1, x + j13);
    sp_simd_store(p0, x + j20);
    sp_simd_store(p3, x + j27);
    sp_simd_store(p2, x + j34);
  }
}

DECLARE_CORE_ROUTINES_SIMD(35)
