#include "ntt/ntt-impl-simd.h"

#define NC 18

static const uint8_t ntt16_fixed_const[NC] = {1, 1, 1, 0, 1, 0, 0, 0, 1};

extern void X(ntt16_init)(spv_t out, sp_t p, sp_t d, sp_t primroot, 
                        sp_t order, sp_t perm);


static void
ntt16_run_core_simd(spv_t in, spv_size_t istride, spv_size_t idist,
		spv_t out, spv_size_t ostride, spv_size_t odist,
		sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  sp_simd_t x0, x1, x2, x3, x4, x5, x6, x7,
            x8, x9, x10, x11, x12, x13, x14, x15;
  sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7, 
            t8, t9, t10, t11, t12, t13, t14, t15, t16, t17;
  sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7,
            p8, p9, p10, p11, p12, p13, p14, p15, p16, p17;

  x0 = sp_simd_gather(in + 0 * istride, idist, vsize);
  x1 = sp_simd_gather(in + 1 * istride, idist, vsize);
  x2 = sp_simd_gather(in + 2 * istride, idist, vsize);
  x3 = sp_simd_gather(in + 3 * istride, idist, vsize);
  x4 = sp_simd_gather(in + 4 * istride, idist, vsize);
  x5 = sp_simd_gather(in + 5 * istride, idist, vsize);
  x6 = sp_simd_gather(in + 6 * istride, idist, vsize);
  x7 = sp_simd_gather(in + 7 * istride, idist, vsize);
  x8 = sp_simd_gather(in + 8 * istride, idist, vsize);
  x9 = sp_simd_gather(in + 9 * istride, idist, vsize);
  x10 = sp_simd_gather(in + 10 * istride, idist, vsize);
  x11 = sp_simd_gather(in + 11 * istride, idist, vsize);
  x12 = sp_simd_gather(in + 12 * istride, idist, vsize);
  x13 = sp_simd_gather(in + 13 * istride, idist, vsize);
  x14 = sp_simd_gather(in + 14 * istride, idist, vsize);
  x15 = sp_simd_gather(in + 15 * istride, idist, vsize);

  t0 = sp_ntt_add_simd(x0, x8, p);
  t8 = sp_ntt_sub_simd(x0, x8, p);
  t1 = sp_ntt_add_simd(x1, x9, p);
  t9 = sp_ntt_sub_simd(x1, x9, p);
  t2 = sp_ntt_add_simd(x2, x10, p);
  t10 = sp_ntt_sub_simd(x2, x10, p);
  t3 = sp_ntt_add_simd(x3, x11, p);
  t11 = sp_ntt_sub_simd(x3, x11, p);
  t4 = sp_ntt_add_simd(x4, x12, p);
  t12 = sp_ntt_sub_partial_simd(x4, x12, p);
  t5 = sp_ntt_add_simd(x5, x13, p);
  t13 = sp_ntt_sub_simd(x13, x5, p);
  t6 = sp_ntt_add_simd(x6, x14, p);
  t14 = sp_ntt_sub_simd(x14, x6, p);
  t7 = sp_ntt_add_simd(x7, x15, p);
  t15 = sp_ntt_sub_simd(x15, x7, p);

  p0 = sp_ntt_add_simd(t0, t4, p);
  p4 = sp_ntt_sub_simd(t0, t4, p);
  p1 = sp_ntt_add_simd(t1, t5, p);
  p5 = sp_ntt_sub_simd(t1, t5, p);
  p2 = sp_ntt_add_simd(t2, t6, p);
  p6 = sp_ntt_sub_partial_simd(t2, t6, p);
  p3 = sp_ntt_add_simd(t3, t7, p);
  p7 = sp_ntt_sub_simd(t3, t7, p);
  p8 = t8;
  p9 = t12;
  p10 = sp_ntt_add_partial_simd(t10, t14, p);
  p11 = sp_ntt_sub_partial_simd(t10, t14, p);
  p12 = sp_ntt_add_simd(t9, t15, p);
  p13 = sp_ntt_sub_simd(t9, t15, p);
  p14 = sp_ntt_add_simd(t13, t11, p);
  p15 = sp_ntt_sub_simd(t13, t11, p);

  t0 = sp_ntt_add_simd(p0, p2, p);
  t1 = sp_ntt_sub_simd(p0, p2, p);
  t2 = sp_ntt_add_simd(p1, p3, p);
  t3 = sp_ntt_sub_partial_simd(p1, p3, p);
  t4 = p4;
  t5 = p6;
  t6 = sp_ntt_sub_partial_simd(p5, p7, p);
  t7 = sp_ntt_add_partial_simd(p5, p7, p); 
  t8 = p8;
  t9 = p9;
  t10 = p10;
  t11 = p11;
  t12 = p12;
  t13 = p13;
  t14 = p14;
  t15 = p15;
  t16 = sp_ntt_add_partial_simd(p12, p14, p);
  t17 = sp_ntt_add_partial_simd(p13, p15, p);

  t3 = sp_ntt_mul_simd(t3, ntt_const[3], ntt_const[NC+3], p);
  t5 = sp_ntt_mul_simd(t5, ntt_const[5], ntt_const[NC+5], p);
  t6 = sp_ntt_mul_simd(t6, ntt_const[6], ntt_const[NC+6], p);
  t7 = sp_ntt_mul_simd(t7, ntt_const[7], ntt_const[NC+7], p);
  t9 = sp_ntt_mul_simd(t9, ntt_const[9], ntt_const[NC+9], p);
  t10 = sp_ntt_mul_simd(t10, ntt_const[10], ntt_const[NC+10], p);
  t11 = sp_ntt_mul_simd(t11, ntt_const[11], ntt_const[NC+11], p);
  t12 = sp_ntt_mul_simd(t12, ntt_const[12], ntt_const[NC+12], p);
  t13 = sp_ntt_mul_simd(t13, ntt_const[13], ntt_const[NC+13], p);
  t14 = sp_ntt_mul_simd(t14, ntt_const[14], ntt_const[NC+14], p);
  t15 = sp_ntt_mul_simd(t15, ntt_const[15], ntt_const[NC+15], p);
  t16 = sp_ntt_mul_simd(t16, ntt_const[16], ntt_const[NC+16], p);
  t17 = sp_ntt_mul_simd(t17, ntt_const[17], ntt_const[NC+17], p);

  p0 = sp_ntt_add_simd(t4, t5, p);
  p1 = sp_ntt_sub_simd(t4, t5, p);
  p2 = sp_ntt_add_simd(t6, t7, p);
  p3 = sp_ntt_sub_simd(t6, t7, p);
  p4 = sp_ntt_add_simd(p0, p2, p);
  p5 = sp_ntt_sub_simd(p0, p2, p);
  p6 = sp_ntt_add_simd(p1, p3, p);
  p7 = sp_ntt_sub_simd(p1, p3, p);
  p0 = sp_ntt_add_simd(t0, t2, p);
  p1 = sp_ntt_sub_simd(t0, t2, p);
  p2 = sp_ntt_add_simd(t1, t3, p);
  p3 = sp_ntt_sub_simd(t1, t3, p);

  t0 = sp_ntt_add_simd(t12, t16, p);
  t1 = sp_ntt_add_simd(t13, t17, p);
  t2 = sp_ntt_add_simd(t14, t16, p);
  t3 = sp_ntt_add_simd(t15, t17, p);
  t7 = sp_ntt_add_simd(t0, t1, p);
  t6 = sp_ntt_sub_simd(t0, t1, p);
  t5 = sp_ntt_add_simd(t2, t3, p);
  t4 = sp_ntt_sub_simd(t2, t3, p); 

  t2 = sp_ntt_add_simd(t10, t11, p);
  t3 = sp_ntt_sub_simd(t10, t11, p);



  t10 = sp_ntt_add_simd(t8, t9, p);
  t11 = sp_ntt_sub_simd(t8, t9, p);

  t12 = sp_ntt_add_simd(t10, t2, p);
  t13 = sp_ntt_sub_simd(t10, t2, p);
  t14 = sp_ntt_add_simd(t11, t3, p);
  t15 = sp_ntt_sub_simd(t11, t3, p);

  p8 = sp_ntt_add_simd(t12, t4, p);
  p9 = sp_ntt_sub_simd(t12, t4, p);
  p10 = sp_ntt_add_simd(t14, t5, p);
  p11 = sp_ntt_sub_simd(t14, t5, p);
  p12 = sp_ntt_add_simd(t13, t6, p);
  p13 = sp_ntt_sub_simd(t13, t6, p);
  p14 = sp_ntt_add_simd(t15, t7, p);
  p15 = sp_ntt_sub_simd(t15, t7, p);

  sp_simd_scatter(p0,  out +  0 * ostride, odist, vsize);
  sp_simd_scatter(p8,  out +  1 * ostride, odist, vsize);
  sp_simd_scatter(p4,  out +  2 * ostride, odist, vsize);
  sp_simd_scatter(p15, out +  3 * ostride, odist, vsize);
  sp_simd_scatter(p2,  out +  4 * ostride, odist, vsize);
  sp_simd_scatter(p12, out +  5 * ostride, odist, vsize);
  sp_simd_scatter(p7,  out +  6 * ostride, odist, vsize);
  sp_simd_scatter(p11, out +  7 * ostride, odist, vsize);
  sp_simd_scatter(p1,  out +  8 * ostride, odist, vsize);
  sp_simd_scatter(p9,  out +  9 * ostride, odist, vsize);
  sp_simd_scatter(p5,  out + 10 * ostride, odist, vsize);
  sp_simd_scatter(p14, out + 11 * ostride, odist, vsize);
  sp_simd_scatter(p3,  out + 12 * ostride, odist, vsize);
  sp_simd_scatter(p13, out + 13 * ostride, odist, vsize);
  sp_simd_scatter(p6,  out + 14 * ostride, odist, vsize);
  sp_simd_scatter(p10, out + 15 * ostride, odist, vsize);
}

static void
ntt16_run_core_simd_interleaved(
		spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		sp_simd_t p, sp_simd_t * c)
{
  sp_simd_t x0, x1, x2, x3, x4, x5, x6, x7,
            x8, x9, x10, x11, x12, x13, x14, x15;
  sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7, 
            t8, t9, t10, t11, t12, t13, t14, t15, t16, t17;
  sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7,
            p8, p9, p10, p11, p12, p13, p14, p15, p16, p17;

  x0 = sp_simd_load(in + 0 * istride);
  x1 = sp_simd_load(in + 1 * istride);
  x2 = sp_simd_load(in + 2 * istride);
  x3 = sp_simd_load(in + 3 * istride);
  x4 = sp_simd_load(in + 4 * istride);
  x5 = sp_simd_load(in + 5 * istride);
  x6 = sp_simd_load(in + 6 * istride);
  x7 = sp_simd_load(in + 7 * istride);
  x8 = sp_simd_load(in + 8 * istride);
  x9 = sp_simd_load(in + 9 * istride);
  x10 = sp_simd_load(in + 10 * istride);
  x11 = sp_simd_load(in + 11 * istride);
  x12 = sp_simd_load(in + 12 * istride);
  x13 = sp_simd_load(in + 13 * istride);
  x14 = sp_simd_load(in + 14 * istride);
  x15 = sp_simd_load(in + 15 * istride);

  t0 = sp_ntt_add_simd0(x0, x8, p);
  t8 = sp_ntt_sub_simd0(x0, x8, p);
  t1 = sp_ntt_add_simd0(x1, x9, p);
  t9 = sp_ntt_sub_simd0(x1, x9, p);
  t2 = sp_ntt_add_simd0(x2, x10, p);
  t10 = sp_ntt_sub_simd0(x2, x10, p);
  t3 = sp_ntt_add_simd0(x3, x11, p);
  t11 = sp_ntt_sub_simd0(x3, x11, p);
  t4 = sp_ntt_add_simd0(x4, x12, p);
  t12 = sp_ntt_sub_partial_simd0(x4, x12, p);
  t5 = sp_ntt_add_simd0(x5, x13, p);
  t13 = sp_ntt_sub_simd0(x13, x5, p);
  t6 = sp_ntt_add_simd0(x6, x14, p);
  t14 = sp_ntt_sub_simd0(x14, x6, p);
  t7 = sp_ntt_add_simd0(x7, x15, p);
  t15 = sp_ntt_sub_simd0(x15, x7, p);

  p0 = sp_ntt_add_simd0(t0, t4, p);
  p4 = sp_ntt_sub_simd0(t0, t4, p);
  p1 = sp_ntt_add_simd0(t1, t5, p);
  p5 = sp_ntt_sub_simd0(t1, t5, p);
  p2 = sp_ntt_add_simd0(t2, t6, p);
  p6 = sp_ntt_sub_partial_simd0(t2, t6, p);
  p3 = sp_ntt_add_simd0(t3, t7, p);
  p7 = sp_ntt_sub_simd0(t3, t7, p);
  p8 = t8;
  p9 = t12;
  p10 = sp_ntt_add_partial_simd0(t10, t14, p);
  p11 = sp_ntt_sub_partial_simd0(t10, t14, p);
  p12 = sp_ntt_add_simd0(t9, t15, p);
  p13 = sp_ntt_sub_simd0(t9, t15, p);
  p14 = sp_ntt_add_simd0(t13, t11, p);
  p15 = sp_ntt_sub_simd0(t13, t11, p);

  t0 = sp_ntt_add_simd0(p0, p2, p);
  t1 = sp_ntt_sub_simd0(p0, p2, p);
  t2 = sp_ntt_add_simd0(p1, p3, p);
  t3 = sp_ntt_sub_partial_simd0(p1, p3, p);
  t4 = p4;
  t5 = p6;
  t6 = sp_ntt_sub_partial_simd0(p5, p7, p);
  t7 = sp_ntt_add_partial_simd0(p5, p7, p); 
  t8 = p8;
  t9 = p9;
  t10 = p10;
  t11 = p11;
  t12 = p12;
  t13 = p13;
  t14 = p14;
  t15 = p15;
  t16 = sp_ntt_add_partial_simd0(p12, p14, p);
  t17 = sp_ntt_add_partial_simd0(p13, p15, p);

  t3 = sp_ntt_mul_simd0(t3, c+6, p);
  t5 = sp_ntt_mul_simd0(t5, c+10, p);
  t6 = sp_ntt_mul_simd0(t6, c+12, p);
  t7 = sp_ntt_mul_simd0(t7, c+14, p);
  t9 = sp_ntt_mul_simd0(t9, c+18, p);
  t10 = sp_ntt_mul_simd0(t10, c+20, p);
  t11 = sp_ntt_mul_simd0(t11, c+22, p);
  t12 = sp_ntt_mul_simd0(t12, c+24, p);
  t13 = sp_ntt_mul_simd0(t13, c+26, p);
  t14 = sp_ntt_mul_simd0(t14, c+28, p);
  t15 = sp_ntt_mul_simd0(t15, c+30, p);
  t16 = sp_ntt_mul_simd0(t16, c+32, p);
  t17 = sp_ntt_mul_simd0(t17, c+34, p);

  p0 = sp_ntt_add_simd0(t4, t5, p);
  p1 = sp_ntt_sub_simd0(t4, t5, p);
  p2 = sp_ntt_add_simd0(t6, t7, p);
  p3 = sp_ntt_sub_simd0(t6, t7, p);
  p4 = sp_ntt_add_simd0(p0, p2, p);
  p5 = sp_ntt_sub_simd0(p0, p2, p);
  p6 = sp_ntt_add_simd0(p1, p3, p);
  p7 = sp_ntt_sub_simd0(p1, p3, p);
  p0 = sp_ntt_add_simd0(t0, t2, p);
  p1 = sp_ntt_sub_simd0(t0, t2, p);
  p2 = sp_ntt_add_simd0(t1, t3, p);
  p3 = sp_ntt_sub_simd0(t1, t3, p);

  t0 = sp_ntt_add_simd0(t12, t16, p);
  t1 = sp_ntt_add_simd0(t13, t17, p);
  t2 = sp_ntt_add_simd0(t14, t16, p);
  t3 = sp_ntt_add_simd0(t15, t17, p);
  t7 = sp_ntt_add_simd0(t0, t1, p);
  t6 = sp_ntt_sub_simd0(t0, t1, p);
  t5 = sp_ntt_add_simd0(t2, t3, p);
  t4 = sp_ntt_sub_simd0(t2, t3, p); 

  t2 = sp_ntt_add_simd0(t10, t11, p);
  t3 = sp_ntt_sub_simd0(t10, t11, p);



  t10 = sp_ntt_add_simd0(t8, t9, p);
  t11 = sp_ntt_sub_simd0(t8, t9, p);

  t12 = sp_ntt_add_simd0(t10, t2, p);
  t13 = sp_ntt_sub_simd0(t10, t2, p);
  t14 = sp_ntt_add_simd0(t11, t3, p);
  t15 = sp_ntt_sub_simd0(t11, t3, p);

  p8 = sp_ntt_add_simd0(t12, t4, p);
  p9 = sp_ntt_sub_simd0(t12, t4, p);
  p10 = sp_ntt_add_simd0(t14, t5, p);
  p11 = sp_ntt_sub_simd0(t14, t5, p);
  p12 = sp_ntt_add_simd0(t13, t6, p);
  p13 = sp_ntt_sub_simd0(t13, t6, p);
  p14 = sp_ntt_add_simd0(t15, t7, p);
  p15 = sp_ntt_sub_simd0(t15, t7, p);

  sp_simd_store(p0,  out +  0 * ostride);
  sp_simd_store(p8,  out +  1 * ostride);
  sp_simd_store(p4,  out +  2 * ostride);
  sp_simd_store(p15, out +  3 * ostride);
  sp_simd_store(p2,  out +  4 * ostride);
  sp_simd_store(p12, out +  5 * ostride);
  sp_simd_store(p7,  out +  6 * ostride);
  sp_simd_store(p11, out +  7 * ostride);
  sp_simd_store(p1,  out +  8 * ostride);
  sp_simd_store(p9,  out +  9 * ostride);
  sp_simd_store(p5,  out + 10 * ostride);
  sp_simd_store(p14, out + 11 * ostride);
  sp_simd_store(p3,  out + 12 * ostride);
  sp_simd_store(p13, out + 13 * ostride);
  sp_simd_store(p6,  out + 14 * ostride);
  sp_simd_store(p10, out + 15 * ostride);
}


static void
ntt16_twiddle_run_core_simd(
        spv_t in, spv_size_t istride, spv_size_t idist,
		spv_t out, spv_size_t ostride, spv_size_t odist,
		sp_simd_t *w, sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  sp_simd_t x0, x1, x2, x3, x4, x5, x6, x7,
            x8, x9, x10, x11, x12, x13, x14, x15;
  sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7, 
            t8, t9, t10, t11, t12, t13, t14, t15, t16, t17;
  sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7,
            p8, p9, p10, p11, p12, p13, p14, p15, p16, p17;

  x0 = sp_simd_gather(in + 0 * istride, idist, vsize);
  x1 = sp_simd_gather(in + 1 * istride, idist, vsize);
  x2 = sp_simd_gather(in + 2 * istride, idist, vsize);
  x3 = sp_simd_gather(in + 3 * istride, idist, vsize);
  x4 = sp_simd_gather(in + 4 * istride, idist, vsize);
  x5 = sp_simd_gather(in + 5 * istride, idist, vsize);
  x6 = sp_simd_gather(in + 6 * istride, idist, vsize);
  x7 = sp_simd_gather(in + 7 * istride, idist, vsize);
  x8 = sp_simd_gather(in + 8 * istride, idist, vsize);
  x9 = sp_simd_gather(in + 9 * istride, idist, vsize);
  x10 = sp_simd_gather(in + 10 * istride, idist, vsize);
  x11 = sp_simd_gather(in + 11 * istride, idist, vsize);
  x12 = sp_simd_gather(in + 12 * istride, idist, vsize);
  x13 = sp_simd_gather(in + 13 * istride, idist, vsize);
  x14 = sp_simd_gather(in + 14 * istride, idist, vsize);
  x15 = sp_simd_gather(in + 15 * istride, idist, vsize);

  t0 = sp_ntt_add_simd(x0, x8, p);
  t8 = sp_ntt_sub_simd(x0, x8, p);
  t1 = sp_ntt_add_simd(x1, x9, p);
  t9 = sp_ntt_sub_simd(x1, x9, p);
  t2 = sp_ntt_add_simd(x2, x10, p);
  t10 = sp_ntt_sub_simd(x2, x10, p);
  t3 = sp_ntt_add_simd(x3, x11, p);
  t11 = sp_ntt_sub_simd(x3, x11, p);
  t4 = sp_ntt_add_simd(x4, x12, p);
  t12 = sp_ntt_sub_partial_simd(x4, x12, p);
  t5 = sp_ntt_add_simd(x5, x13, p);
  t13 = sp_ntt_sub_simd(x13, x5, p);
  t6 = sp_ntt_add_simd(x6, x14, p);
  t14 = sp_ntt_sub_simd(x14, x6, p);
  t7 = sp_ntt_add_simd(x7, x15, p);
  t15 = sp_ntt_sub_simd(x15, x7, p);

  p0 = sp_ntt_add_simd(t0, t4, p);
  p4 = sp_ntt_sub_simd(t0, t4, p);
  p1 = sp_ntt_add_simd(t1, t5, p);
  p5 = sp_ntt_sub_simd(t1, t5, p);
  p2 = sp_ntt_add_simd(t2, t6, p);
  p6 = sp_ntt_sub_partial_simd(t2, t6, p);
  p3 = sp_ntt_add_simd(t3, t7, p);
  p7 = sp_ntt_sub_simd(t3, t7, p);
  p8 = t8;
  p9 = t12;
  p10 = sp_ntt_add_partial_simd(t10, t14, p);
  p11 = sp_ntt_sub_partial_simd(t10, t14, p);
  p12 = sp_ntt_add_simd(t9, t15, p);
  p13 = sp_ntt_sub_simd(t9, t15, p);
  p14 = sp_ntt_add_simd(t13, t11, p);
  p15 = sp_ntt_sub_simd(t13, t11, p);

  t0 = sp_ntt_add_simd(p0, p2, p);
  t1 = sp_ntt_sub_simd(p0, p2, p);
  t2 = sp_ntt_add_simd(p1, p3, p);
  t3 = sp_ntt_sub_partial_simd(p1, p3, p);
  t4 = p4;
  t5 = p6;
  t6 = sp_ntt_sub_partial_simd(p5, p7, p);
  t7 = sp_ntt_add_partial_simd(p5, p7, p); 
  t8 = p8;
  t9 = p9;
  t10 = p10;
  t11 = p11;
  t12 = p12;
  t13 = p13;
  t14 = p14;
  t15 = p15;
  t16 = sp_ntt_add_partial_simd(p12, p14, p);
  t17 = sp_ntt_add_partial_simd(p13, p15, p);

  t3 = sp_ntt_mul_simd(t3, ntt_const[3], ntt_const[NC+3], p);
  t5 = sp_ntt_mul_simd(t5, ntt_const[5], ntt_const[NC+5], p);
  t6 = sp_ntt_mul_simd(t6, ntt_const[6], ntt_const[NC+6], p);
  t7 = sp_ntt_mul_simd(t7, ntt_const[7], ntt_const[NC+7], p);
  t9 = sp_ntt_mul_simd(t9, ntt_const[9], ntt_const[NC+9], p);
  t10 = sp_ntt_mul_simd(t10, ntt_const[10], ntt_const[NC+10], p);
  t11 = sp_ntt_mul_simd(t11, ntt_const[11], ntt_const[NC+11], p);
  t12 = sp_ntt_mul_simd(t12, ntt_const[12], ntt_const[NC+12], p);
  t13 = sp_ntt_mul_simd(t13, ntt_const[13], ntt_const[NC+13], p);
  t14 = sp_ntt_mul_simd(t14, ntt_const[14], ntt_const[NC+14], p);
  t15 = sp_ntt_mul_simd(t15, ntt_const[15], ntt_const[NC+15], p);
  t16 = sp_ntt_mul_simd(t16, ntt_const[16], ntt_const[NC+16], p);
  t17 = sp_ntt_mul_simd(t17, ntt_const[17], ntt_const[NC+17], p);

  p0 = sp_ntt_add_simd(t4, t5, p);
  p1 = sp_ntt_sub_simd(t4, t5, p);
  p2 = sp_ntt_add_simd(t6, t7, p);
  p3 = sp_ntt_sub_simd(t6, t7, p);
  p4 = sp_ntt_add_partial_simd(p0, p2, p);
  p5 = sp_ntt_sub_partial_simd(p0, p2, p);
  p6 = sp_ntt_add_partial_simd(p1, p3, p);
  p7 = sp_ntt_sub_partial_simd(p1, p3, p);
  p0 = sp_ntt_add_simd(t0, t2, p);
  p1 = sp_ntt_sub_partial_simd(t0, t2, p);
  p2 = sp_ntt_add_partial_simd(t1, t3, p);
  p3 = sp_ntt_sub_partial_simd(t1, t3, p);

  t0 = sp_ntt_add_simd(t12, t16, p);
  t1 = sp_ntt_add_simd(t13, t17, p);
  t2 = sp_ntt_add_simd(t14, t16, p);
  t3 = sp_ntt_add_simd(t15, t17, p);
  t7 = sp_ntt_add_simd(t0, t1, p);
  t6 = sp_ntt_sub_simd(t0, t1, p);
  t5 = sp_ntt_add_simd(t2, t3, p);
  t4 = sp_ntt_sub_simd(t2, t3, p); 

  t2 = sp_ntt_add_simd(t10, t11, p);
  t3 = sp_ntt_sub_simd(t10, t11, p);



  t10 = sp_ntt_add_simd(t8, t9, p);
  t11 = sp_ntt_sub_simd(t8, t9, p);

  t12 = sp_ntt_add_simd(t10, t2, p);
  t13 = sp_ntt_sub_simd(t10, t2, p);
  t14 = sp_ntt_add_simd(t11, t3, p);
  t15 = sp_ntt_sub_simd(t11, t3, p);

  p8 = sp_ntt_add_partial_simd(t12, t4, p);
  p9 = sp_ntt_sub_partial_simd(t12, t4, p);
  p10 = sp_ntt_add_partial_simd(t14, t5, p);
  p11 = sp_ntt_sub_partial_simd(t14, t5, p);
  p12 = sp_ntt_add_partial_simd(t13, t6, p);
  p13 = sp_ntt_sub_partial_simd(t13, t6, p);
  p14 = sp_ntt_add_partial_simd(t15, t7, p);
  p15 = sp_ntt_sub_partial_simd(t15, t7, p);

  p8  = sp_ntt_twiddle_mul_simd(p8, w + 0, p);
  p4  = sp_ntt_twiddle_mul_simd(p4, w + 2, p);
  p15 = sp_ntt_twiddle_mul_simd(p15, w + 4, p);
  p2  = sp_ntt_twiddle_mul_simd(p2, w + 6, p);
  p12 = sp_ntt_twiddle_mul_simd(p12, w + 8, p);
  p7  = sp_ntt_twiddle_mul_simd(p7, w + 10, p);
  p11 = sp_ntt_twiddle_mul_simd(p11, w + 12, p);
  p1  = sp_ntt_twiddle_mul_simd(p1, w + 14, p);
  p9  = sp_ntt_twiddle_mul_simd(p9, w + 16, p);
  p5  = sp_ntt_twiddle_mul_simd(p5, w + 18, p);
  p14 = sp_ntt_twiddle_mul_simd(p14, w + 20, p);
  p3  = sp_ntt_twiddle_mul_simd(p3, w + 22, p);
  p13 = sp_ntt_twiddle_mul_simd(p13, w + 24, p);
  p6  = sp_ntt_twiddle_mul_simd(p6, w + 26, p);
  p10 = sp_ntt_twiddle_mul_simd(p10, w + 28, p);

  sp_simd_scatter(p0, out +  0 * ostride, odist, vsize);
  sp_simd_scatter(p8, out +  1 * ostride, odist, vsize);
  sp_simd_scatter(p4, out +  2 * ostride, odist, vsize);
  sp_simd_scatter(p15, out +  3 * ostride, odist, vsize);
  sp_simd_scatter(p2, out +  4 * ostride, odist, vsize);
  sp_simd_scatter(p12, out +  5 * ostride, odist, vsize);
  sp_simd_scatter(p7, out +  6 * ostride, odist, vsize);
  sp_simd_scatter(p11, out +  7 * ostride, odist, vsize);
  sp_simd_scatter(p1, out +  8 * ostride, odist, vsize);
  sp_simd_scatter(p9, out +  9 * ostride, odist, vsize);
  sp_simd_scatter(p5, out + 10 * ostride, odist, vsize);
  sp_simd_scatter(p14, out + 11 * ostride, odist, vsize);
  sp_simd_scatter(p3, out + 12 * ostride, odist, vsize);
  sp_simd_scatter(p13, out + 13 * ostride, odist, vsize);
  sp_simd_scatter(p6, out + 14 * ostride, odist, vsize);
  sp_simd_scatter(p10, out + 15 * ostride, odist, vsize);
}

static void
ntt16_twiddle_run_core_simd_interleaved(
		spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		sp_simd_t *w, sp_simd_t p, sp_simd_t * c)
{
  sp_simd_t x0, x1, x2, x3, x4, x5, x6, x7,
            x8, x9, x10, x11, x12, x13, x14, x15;
  sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7, 
            t8, t9, t10, t11, t12, t13, t14, t15, t16, t17;
  sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7,
            p8, p9, p10, p11, p12, p13, p14, p15, p16, p17;

  x0 = sp_simd_load(in + 0 * istride);
  x1 = sp_simd_load(in + 1 * istride);
  x2 = sp_simd_load(in + 2 * istride);
  x3 = sp_simd_load(in + 3 * istride);
  x4 = sp_simd_load(in + 4 * istride);
  x5 = sp_simd_load(in + 5 * istride);
  x6 = sp_simd_load(in + 6 * istride);
  x7 = sp_simd_load(in + 7 * istride);
  x8 = sp_simd_load(in + 8 * istride);
  x9 = sp_simd_load(in + 9 * istride);
  x10 = sp_simd_load(in + 10 * istride);
  x11 = sp_simd_load(in + 11 * istride);
  x12 = sp_simd_load(in + 12 * istride);
  x13 = sp_simd_load(in + 13 * istride);
  x14 = sp_simd_load(in + 14 * istride);
  x15 = sp_simd_load(in + 15 * istride);

  t0 = sp_ntt_add_simd0(x0, x8, p);
  t8 = sp_ntt_sub_simd0(x0, x8, p);
  t1 = sp_ntt_add_simd0(x1, x9, p);
  t9 = sp_ntt_sub_simd0(x1, x9, p);
  t2 = sp_ntt_add_simd0(x2, x10, p);
  t10 = sp_ntt_sub_simd0(x2, x10, p);
  t3 = sp_ntt_add_simd0(x3, x11, p);
  t11 = sp_ntt_sub_simd0(x3, x11, p);
  t4 = sp_ntt_add_simd0(x4, x12, p);
  t12 = sp_ntt_sub_partial_simd0(x4, x12, p);
  t5 = sp_ntt_add_simd0(x5, x13, p);
  t13 = sp_ntt_sub_simd0(x13, x5, p);
  t6 = sp_ntt_add_simd0(x6, x14, p);
  t14 = sp_ntt_sub_simd0(x14, x6, p);
  t7 = sp_ntt_add_simd0(x7, x15, p);
  t15 = sp_ntt_sub_simd0(x15, x7, p);

  p0 = sp_ntt_add_simd0(t0, t4, p);
  p4 = sp_ntt_sub_simd0(t0, t4, p);
  p1 = sp_ntt_add_simd0(t1, t5, p);
  p5 = sp_ntt_sub_simd0(t1, t5, p);
  p2 = sp_ntt_add_simd0(t2, t6, p);
  p6 = sp_ntt_sub_partial_simd0(t2, t6, p);
  p3 = sp_ntt_add_simd0(t3, t7, p);
  p7 = sp_ntt_sub_simd0(t3, t7, p);
  p8 = t8;
  p9 = t12;
  p10 = sp_ntt_add_partial_simd0(t10, t14, p);
  p11 = sp_ntt_sub_partial_simd0(t10, t14, p);
  p12 = sp_ntt_add_simd0(t9, t15, p);
  p13 = sp_ntt_sub_simd0(t9, t15, p);
  p14 = sp_ntt_add_simd0(t13, t11, p);
  p15 = sp_ntt_sub_simd0(t13, t11, p);

  t0 = sp_ntt_add_simd0(p0, p2, p);
  t1 = sp_ntt_sub_simd0(p0, p2, p);
  t2 = sp_ntt_add_simd0(p1, p3, p);
  t3 = sp_ntt_sub_partial_simd0(p1, p3, p);
  t4 = p4;
  t5 = p6;
  t6 = sp_ntt_sub_partial_simd0(p5, p7, p);
  t7 = sp_ntt_add_partial_simd0(p5, p7, p); 
  t8 = p8;
  t9 = p9;
  t10 = p10;
  t11 = p11;
  t12 = p12;
  t13 = p13;
  t14 = p14;
  t15 = p15;
  t16 = sp_ntt_add_partial_simd0(p12, p14, p);
  t17 = sp_ntt_add_partial_simd0(p13, p15, p);

  t3 = sp_ntt_mul_simd0(t3, c+6, p);
  t5 = sp_ntt_mul_simd0(t5, c+10, p);
  t6 = sp_ntt_mul_simd0(t6, c+12, p);
  t7 = sp_ntt_mul_simd0(t7, c+14, p);
  t9 = sp_ntt_mul_simd0(t9, c+18, p);
  t10 = sp_ntt_mul_simd0(t10, c+20, p);
  t11 = sp_ntt_mul_simd0(t11, c+22, p);
  t12 = sp_ntt_mul_simd0(t12, c+24, p);
  t13 = sp_ntt_mul_simd0(t13, c+26, p);
  t14 = sp_ntt_mul_simd0(t14, c+28, p);
  t15 = sp_ntt_mul_simd0(t15, c+30, p);
  t16 = sp_ntt_mul_simd0(t16, c+32, p);
  t17 = sp_ntt_mul_simd0(t17, c+34, p);

  p0 = sp_ntt_add_simd0(t4, t5, p);
  p1 = sp_ntt_sub_simd0(t4, t5, p);
  p2 = sp_ntt_add_simd0(t6, t7, p);
  p3 = sp_ntt_sub_simd0(t6, t7, p);
  p4 = sp_ntt_add_partial_simd0(p0, p2, p);
  p5 = sp_ntt_sub_partial_simd0(p0, p2, p);
  p6 = sp_ntt_add_partial_simd0(p1, p3, p);
  p7 = sp_ntt_sub_partial_simd0(p1, p3, p);
  p0 = sp_ntt_add_simd0(t0, t2, p);
  p1 = sp_ntt_sub_partial_simd0(t0, t2, p);
  p2 = sp_ntt_add_partial_simd0(t1, t3, p);
  p3 = sp_ntt_sub_partial_simd0(t1, t3, p);

  t0 = sp_ntt_add_simd0(t12, t16, p);
  t1 = sp_ntt_add_simd0(t13, t17, p);
  t2 = sp_ntt_add_simd0(t14, t16, p);
  t3 = sp_ntt_add_simd0(t15, t17, p);
  t7 = sp_ntt_add_simd0(t0, t1, p);
  t6 = sp_ntt_sub_simd0(t0, t1, p);
  t5 = sp_ntt_add_simd0(t2, t3, p);
  t4 = sp_ntt_sub_simd0(t2, t3, p); 

  t2 = sp_ntt_add_simd0(t10, t11, p);
  t3 = sp_ntt_sub_simd0(t10, t11, p);



  t10 = sp_ntt_add_simd0(t8, t9, p);
  t11 = sp_ntt_sub_simd0(t8, t9, p);

  t12 = sp_ntt_add_simd0(t10, t2, p);
  t13 = sp_ntt_sub_simd0(t10, t2, p);
  t14 = sp_ntt_add_simd0(t11, t3, p);
  t15 = sp_ntt_sub_simd0(t11, t3, p);

  p8 = sp_ntt_add_partial_simd0(t12, t4, p);
  p9 = sp_ntt_sub_partial_simd0(t12, t4, p);
  p10 = sp_ntt_add_partial_simd0(t14, t5, p);
  p11 = sp_ntt_sub_partial_simd0(t14, t5, p);
  p12 = sp_ntt_add_partial_simd0(t13, t6, p);
  p13 = sp_ntt_sub_partial_simd0(t13, t6, p);
  p14 = sp_ntt_add_partial_simd0(t15, t7, p);
  p15 = sp_ntt_sub_partial_simd0(t15, t7, p);

  p8  = sp_ntt_twiddle_mul_simd0(p8, w+0, p);
  p4  = sp_ntt_twiddle_mul_simd0(p4, w+2, p);
  p15 = sp_ntt_twiddle_mul_simd0(p15, w+4, p);
  p2  = sp_ntt_twiddle_mul_simd0(p2, w+6, p);
  p12 = sp_ntt_twiddle_mul_simd0(p12, w+8, p);
  p7  = sp_ntt_twiddle_mul_simd0(p7, w+10, p);
  p11 = sp_ntt_twiddle_mul_simd0(p11, w+12, p);
  p1  = sp_ntt_twiddle_mul_simd0(p1, w+14, p);
  p9  = sp_ntt_twiddle_mul_simd0(p9, w+16, p);
  p5  = sp_ntt_twiddle_mul_simd0(p5, w+18, p);
  p14 = sp_ntt_twiddle_mul_simd0(p14, w+20, p);
  p3  = sp_ntt_twiddle_mul_simd0(p3, w+22, p);
  p13 = sp_ntt_twiddle_mul_simd0(p13, w+24, p);
  p6  = sp_ntt_twiddle_mul_simd0(p6, w+26, p);
  p10 = sp_ntt_twiddle_mul_simd0(p10, w+28, p);

  sp_simd_store(p0, out +  0 * ostride);
  sp_simd_store(p8, out +  1 * ostride);
  sp_simd_store(p4, out +  2 * ostride);
  sp_simd_store(p15, out +  3 * ostride);
  sp_simd_store(p2, out +  4 * ostride);
  sp_simd_store(p12, out +  5 * ostride);
  sp_simd_store(p7, out +  6 * ostride);
  sp_simd_store(p11, out +  7 * ostride);
  sp_simd_store(p1, out +  8 * ostride);
  sp_simd_store(p9, out +  9 * ostride);
  sp_simd_store(p5, out + 10 * ostride);
  sp_simd_store(p14, out + 11 * ostride);
  sp_simd_store(p3, out + 12 * ostride);
  sp_simd_store(p13, out + 13 * ostride);
  sp_simd_store(p6, out + 14 * ostride);
  sp_simd_store(p10, out + 15 * ostride);
}

static void
ntt16_pfa_run_core_simd(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t inc2, spv_size_t n,
	  sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  spv_size_t j0, j1, j2, j3, j4, j5, j6, j7,
             j8, j9, j10, j11, j12, j13, j14, j15;
  sp_simd_t x0, x1, x2, x3, x4, x5, x6, x7,
            x8, x9, x10, x11, x12, x13, x14, x15;
  sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7, 
            t8, t9, t10, t11, t12, t13, t14, t15, t16, t17;
  sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7,
            p8, p9, p10, p11, p12, p13, p14, p15, p16, p17;

  j0 = start;
  j1 = sp_array_inc(j0, inc, n);
  j2 = sp_array_inc(j0, 2 * inc, n);
  j3 = sp_array_inc(j0, 3 * inc, n);
  j4 = sp_array_inc(j0, 4 * inc, n);
  j5 = sp_array_inc(j0, 5 * inc, n);
  j6 = sp_array_inc(j0, 6 * inc, n);
  j7 = sp_array_inc(j0, 7 * inc, n);
  j8 = sp_array_inc(j0, 8 * inc, n);
  j9 = sp_array_inc(j0, 9 * inc, n);
  j10 = sp_array_inc(j0, 10 * inc, n);
  j11 = sp_array_inc(j0, 11 * inc, n);
  j12 = sp_array_inc(j0, 12 * inc, n);
  j13 = sp_array_inc(j0, 13 * inc, n);
  j14 = sp_array_inc(j0, 14 * inc, n);
  j15 = sp_array_inc(j0, 15 * inc, n);

  x0 = sp_simd_pfa_gather(x, j0, inc2, n, vsize);
  x1 = sp_simd_pfa_gather(x, j1, inc2, n, vsize);
  x2 = sp_simd_pfa_gather(x, j2, inc2, n, vsize);
  x3 = sp_simd_pfa_gather(x, j3, inc2, n, vsize);
  x4 = sp_simd_pfa_gather(x, j4, inc2, n, vsize);
  x5 = sp_simd_pfa_gather(x, j5, inc2, n, vsize);
  x6 = sp_simd_pfa_gather(x, j6, inc2, n, vsize);
  x7 = sp_simd_pfa_gather(x, j7, inc2, n, vsize);
  x8 = sp_simd_pfa_gather(x, j8, inc2, n, vsize);
  x9 = sp_simd_pfa_gather(x, j9, inc2, n, vsize);
  x10 = sp_simd_pfa_gather(x, j10, inc2, n, vsize);
  x11 = sp_simd_pfa_gather(x, j11, inc2, n, vsize);
  x12 = sp_simd_pfa_gather(x, j12, inc2, n, vsize);
  x13 = sp_simd_pfa_gather(x, j13, inc2, n, vsize);
  x14 = sp_simd_pfa_gather(x, j14, inc2, n, vsize);
  x15 = sp_simd_pfa_gather(x, j15, inc2, n, vsize);

  t0 = sp_ntt_add_simd(x0, x8, p);
  t8 = sp_ntt_sub_simd(x0, x8, p);
  t1 = sp_ntt_add_simd(x1, x9, p);
  t9 = sp_ntt_sub_simd(x1, x9, p);
  t2 = sp_ntt_add_simd(x2, x10, p);
  t10 = sp_ntt_sub_simd(x2, x10, p);
  t3 = sp_ntt_add_simd(x3, x11, p);
  t11 = sp_ntt_sub_simd(x3, x11, p);
  t4 = sp_ntt_add_simd(x4, x12, p);
  t12 = sp_ntt_sub_partial_simd(x4, x12, p);
  t5 = sp_ntt_add_simd(x5, x13, p);
  t13 = sp_ntt_sub_simd(x13, x5, p);
  t6 = sp_ntt_add_simd(x6, x14, p);
  t14 = sp_ntt_sub_simd(x14, x6, p);
  t7 = sp_ntt_add_simd(x7, x15, p);
  t15 = sp_ntt_sub_simd(x15, x7, p);

  p0 = sp_ntt_add_simd(t0, t4, p);
  p4 = sp_ntt_sub_simd(t0, t4, p);
  p1 = sp_ntt_add_simd(t1, t5, p);
  p5 = sp_ntt_sub_simd(t1, t5, p);
  p2 = sp_ntt_add_simd(t2, t6, p);
  p6 = sp_ntt_sub_partial_simd(t2, t6, p);
  p3 = sp_ntt_add_simd(t3, t7, p);
  p7 = sp_ntt_sub_simd(t3, t7, p);
  p8 = t8;
  p9 = t12;
  p10 = sp_ntt_add_partial_simd(t10, t14, p);
  p11 = sp_ntt_sub_partial_simd(t10, t14, p);
  p12 = sp_ntt_add_simd(t9, t15, p);
  p13 = sp_ntt_sub_simd(t9, t15, p);
  p14 = sp_ntt_add_simd(t13, t11, p);
  p15 = sp_ntt_sub_simd(t13, t11, p);

  t0 = sp_ntt_add_simd(p0, p2, p);
  t1 = sp_ntt_sub_simd(p0, p2, p);
  t2 = sp_ntt_add_simd(p1, p3, p);
  t3 = sp_ntt_sub_partial_simd(p1, p3, p);
  t4 = p4;
  t5 = p6;
  t6 = sp_ntt_sub_partial_simd(p5, p7, p);
  t7 = sp_ntt_add_partial_simd(p5, p7, p); 
  t8 = p8;
  t9 = p9;
  t10 = p10;
  t11 = p11;
  t12 = p12;
  t13 = p13;
  t14 = p14;
  t15 = p15;
  t16 = sp_ntt_add_partial_simd(p12, p14, p);
  t17 = sp_ntt_add_partial_simd(p13, p15, p);

  t3 = sp_ntt_mul_simd(t3, ntt_const[3], ntt_const[NC+3], p);
  t5 = sp_ntt_mul_simd(t5, ntt_const[5], ntt_const[NC+5], p);
  t6 = sp_ntt_mul_simd(t6, ntt_const[6], ntt_const[NC+6], p);
  t7 = sp_ntt_mul_simd(t7, ntt_const[7], ntt_const[NC+7], p);
  t9 = sp_ntt_mul_simd(t9, ntt_const[9], ntt_const[NC+9], p);
  t10 = sp_ntt_mul_simd(t10, ntt_const[10], ntt_const[NC+10], p);
  t11 = sp_ntt_mul_simd(t11, ntt_const[11], ntt_const[NC+11], p);
  t12 = sp_ntt_mul_simd(t12, ntt_const[12], ntt_const[NC+12], p);
  t13 = sp_ntt_mul_simd(t13, ntt_const[13], ntt_const[NC+13], p);
  t14 = sp_ntt_mul_simd(t14, ntt_const[14], ntt_const[NC+14], p);
  t15 = sp_ntt_mul_simd(t15, ntt_const[15], ntt_const[NC+15], p);
  t16 = sp_ntt_mul_simd(t16, ntt_const[16], ntt_const[NC+16], p);
  t17 = sp_ntt_mul_simd(t17, ntt_const[17], ntt_const[NC+17], p);

  p0 = sp_ntt_add_simd(t4, t5, p);
  p1 = sp_ntt_sub_simd(t4, t5, p);
  p2 = sp_ntt_add_simd(t6, t7, p);
  p3 = sp_ntt_sub_simd(t6, t7, p);
  p4 = sp_ntt_add_simd(p0, p2, p);
  p5 = sp_ntt_sub_simd(p0, p2, p);
  p6 = sp_ntt_add_simd(p1, p3, p);
  p7 = sp_ntt_sub_simd(p1, p3, p);
  p0 = sp_ntt_add_simd(t0, t2, p);
  p1 = sp_ntt_sub_simd(t0, t2, p);
  p2 = sp_ntt_add_simd(t1, t3, p);
  p3 = sp_ntt_sub_simd(t1, t3, p);

  t0 = sp_ntt_add_simd(t12, t16, p);
  t1 = sp_ntt_add_simd(t13, t17, p);
  t2 = sp_ntt_add_simd(t14, t16, p);
  t3 = sp_ntt_add_simd(t15, t17, p);
  t7 = sp_ntt_add_simd(t0, t1, p);
  t6 = sp_ntt_sub_simd(t0, t1, p);
  t5 = sp_ntt_add_simd(t2, t3, p);
  t4 = sp_ntt_sub_simd(t2, t3, p); 

  t2 = sp_ntt_add_simd(t10, t11, p);
  t3 = sp_ntt_sub_simd(t10, t11, p);



  t10 = sp_ntt_add_simd(t8, t9, p);
  t11 = sp_ntt_sub_simd(t8, t9, p);

  t12 = sp_ntt_add_simd(t10, t2, p);
  t13 = sp_ntt_sub_simd(t10, t2, p);
  t14 = sp_ntt_add_simd(t11, t3, p);
  t15 = sp_ntt_sub_simd(t11, t3, p);

  p8 = sp_ntt_add_simd(t12, t4, p);
  p9 = sp_ntt_sub_simd(t12, t4, p);
  p10 = sp_ntt_add_simd(t14, t5, p);
  p11 = sp_ntt_sub_simd(t14, t5, p);
  p12 = sp_ntt_add_simd(t13, t6, p);
  p13 = sp_ntt_sub_simd(t13, t6, p);
  p14 = sp_ntt_add_simd(t15, t7, p);
  p15 = sp_ntt_sub_simd(t15, t7, p);

  sp_simd_pfa_scatter(p0, x, j0, inc2, n, vsize);
  sp_simd_pfa_scatter(p8, x, j1, inc2, n, vsize);
  sp_simd_pfa_scatter(p4, x, j2, inc2, n, vsize);
  sp_simd_pfa_scatter(p15,x, j3, inc2, n, vsize);
  sp_simd_pfa_scatter(p2, x, j4, inc2, n, vsize);
  sp_simd_pfa_scatter(p12,x, j5, inc2, n, vsize);
  sp_simd_pfa_scatter(p7, x, j6, inc2, n, vsize);
  sp_simd_pfa_scatter(p11,x, j7, inc2, n, vsize);
  sp_simd_pfa_scatter(p1, x, j8, inc2, n, vsize);
  sp_simd_pfa_scatter(p9, x, j9, inc2, n, vsize);
  sp_simd_pfa_scatter(p5, x, j10, inc2, n, vsize);
  sp_simd_pfa_scatter(p14,x, j11, inc2, n, vsize);
  sp_simd_pfa_scatter(p3, x, j12, inc2, n, vsize);
  sp_simd_pfa_scatter(p13,x, j13, inc2, n, vsize);
  sp_simd_pfa_scatter(p6, x, j14, inc2, n, vsize);
  sp_simd_pfa_scatter(p10,x, j15, inc2, n, vsize);
}

static void
ntt16_pfa_run_core_simd_interleaved(
    	  spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t n,
	  sp_simd_t p, sp_simd_t * c)
{
  spv_size_t j0, j1, j2, j3, j4, j5, j6, j7,
             j8, j9, j10, j11, j12, j13, j14, j15;
  sp_simd_t x0, x1, x2, x3, x4, x5, x6, x7,
            x8, x9, x10, x11, x12, x13, x14, x15;
  sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7, 
            t8, t9, t10, t11, t12, t13, t14, t15, t16, t17;
  sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7,
            p8, p9, p10, p11, p12, p13, p14, p15, p16, p17;

  j0 = start;
  j1 = sp_array_inc(j0, inc, n);
  j2 = sp_array_inc(j0, 2 * inc, n);
  j3 = sp_array_inc(j0, 3 * inc, n);
  j4 = sp_array_inc(j0, 4 * inc, n);
  j5 = sp_array_inc(j0, 5 * inc, n);
  j6 = sp_array_inc(j0, 6 * inc, n);
  j7 = sp_array_inc(j0, 7 * inc, n);
  j8 = sp_array_inc(j0, 8 * inc, n);
  j9 = sp_array_inc(j0, 9 * inc, n);
  j10 = sp_array_inc(j0, 10 * inc, n);
  j11 = sp_array_inc(j0, 11 * inc, n);
  j12 = sp_array_inc(j0, 12 * inc, n);
  j13 = sp_array_inc(j0, 13 * inc, n);
  j14 = sp_array_inc(j0, 14 * inc, n);
  j15 = sp_array_inc(j0, 15 * inc, n);

  x0 = sp_simd_load(x + j0);
  x1 = sp_simd_load(x + j1);
  x2 = sp_simd_load(x + j2);
  x3 = sp_simd_load(x + j3);
  x4 = sp_simd_load(x + j4);
  x5 = sp_simd_load(x + j5);
  x6 = sp_simd_load(x + j6);
  x7 = sp_simd_load(x + j7);
  x8 = sp_simd_load(x + j8);
  x9 = sp_simd_load(x + j9);
  x10 = sp_simd_load(x + j10);
  x11 = sp_simd_load(x + j11);
  x12 = sp_simd_load(x + j12);
  x13 = sp_simd_load(x + j13);
  x14 = sp_simd_load(x + j14);
  x15 = sp_simd_load(x + j15);

  t0 = sp_ntt_add_simd0(x0, x8, p);
  t8 = sp_ntt_sub_simd0(x0, x8, p);
  t1 = sp_ntt_add_simd0(x1, x9, p);
  t9 = sp_ntt_sub_simd0(x1, x9, p);
  t2 = sp_ntt_add_simd0(x2, x10, p);
  t10 = sp_ntt_sub_simd0(x2, x10, p);
  t3 = sp_ntt_add_simd0(x3, x11, p);
  t11 = sp_ntt_sub_simd0(x3, x11, p);
  t4 = sp_ntt_add_simd0(x4, x12, p);
  t12 = sp_ntt_sub_partial_simd0(x4, x12, p);
  t5 = sp_ntt_add_simd0(x5, x13, p);
  t13 = sp_ntt_sub_simd0(x13, x5, p);
  t6 = sp_ntt_add_simd0(x6, x14, p);
  t14 = sp_ntt_sub_simd0(x14, x6, p);
  t7 = sp_ntt_add_simd0(x7, x15, p);
  t15 = sp_ntt_sub_simd0(x15, x7, p);

  p0 = sp_ntt_add_simd0(t0, t4, p);
  p4 = sp_ntt_sub_simd0(t0, t4, p);
  p1 = sp_ntt_add_simd0(t1, t5, p);
  p5 = sp_ntt_sub_simd0(t1, t5, p);
  p2 = sp_ntt_add_simd0(t2, t6, p);
  p6 = sp_ntt_sub_partial_simd0(t2, t6, p);
  p3 = sp_ntt_add_simd0(t3, t7, p);
  p7 = sp_ntt_sub_simd0(t3, t7, p);
  p8 = t8;
  p9 = t12;
  p10 = sp_ntt_add_partial_simd0(t10, t14, p);
  p11 = sp_ntt_sub_partial_simd0(t10, t14, p);
  p12 = sp_ntt_add_simd0(t9, t15, p);
  p13 = sp_ntt_sub_simd0(t9, t15, p);
  p14 = sp_ntt_add_simd0(t13, t11, p);
  p15 = sp_ntt_sub_simd0(t13, t11, p);

  t0 = sp_ntt_add_simd0(p0, p2, p);
  t1 = sp_ntt_sub_simd0(p0, p2, p);
  t2 = sp_ntt_add_simd0(p1, p3, p);
  t3 = sp_ntt_sub_partial_simd0(p1, p3, p);
  t4 = p4;
  t5 = p6;
  t6 = sp_ntt_sub_partial_simd0(p5, p7, p);
  t7 = sp_ntt_add_partial_simd0(p5, p7, p); 
  t8 = p8;
  t9 = p9;
  t10 = p10;
  t11 = p11;
  t12 = p12;
  t13 = p13;
  t14 = p14;
  t15 = p15;
  t16 = sp_ntt_add_partial_simd0(p12, p14, p);
  t17 = sp_ntt_add_partial_simd0(p13, p15, p);

  t3 = sp_ntt_mul_simd0(t3, c+6, p);
  t5 = sp_ntt_mul_simd0(t5, c+10, p);
  t6 = sp_ntt_mul_simd0(t6, c+12, p);
  t7 = sp_ntt_mul_simd0(t7, c+14, p);
  t9 = sp_ntt_mul_simd0(t9, c+18, p);
  t10 = sp_ntt_mul_simd0(t10, c+20, p);
  t11 = sp_ntt_mul_simd0(t11, c+22, p);
  t12 = sp_ntt_mul_simd0(t12, c+24, p);
  t13 = sp_ntt_mul_simd0(t13, c+26, p);
  t14 = sp_ntt_mul_simd0(t14, c+28, p);
  t15 = sp_ntt_mul_simd0(t15, c+30, p);
  t16 = sp_ntt_mul_simd0(t16, c+32, p);
  t17 = sp_ntt_mul_simd0(t17, c+34, p);

  p0 = sp_ntt_add_simd0(t4, t5, p);
  p1 = sp_ntt_sub_simd0(t4, t5, p);
  p2 = sp_ntt_add_simd0(t6, t7, p);
  p3 = sp_ntt_sub_simd0(t6, t7, p);
  p4 = sp_ntt_add_simd0(p0, p2, p);
  p5 = sp_ntt_sub_simd0(p0, p2, p);
  p6 = sp_ntt_add_simd0(p1, p3, p);
  p7 = sp_ntt_sub_simd0(p1, p3, p);
  p0 = sp_ntt_add_simd0(t0, t2, p);
  p1 = sp_ntt_sub_simd0(t0, t2, p);
  p2 = sp_ntt_add_simd0(t1, t3, p);
  p3 = sp_ntt_sub_simd0(t1, t3, p);

  t0 = sp_ntt_add_simd0(t12, t16, p);
  t1 = sp_ntt_add_simd0(t13, t17, p);
  t2 = sp_ntt_add_simd0(t14, t16, p);
  t3 = sp_ntt_add_simd0(t15, t17, p);
  t7 = sp_ntt_add_simd0(t0, t1, p);
  t6 = sp_ntt_sub_simd0(t0, t1, p);
  t5 = sp_ntt_add_simd0(t2, t3, p);
  t4 = sp_ntt_sub_simd0(t2, t3, p); 

  t2 = sp_ntt_add_simd0(t10, t11, p);
  t3 = sp_ntt_sub_simd0(t10, t11, p);



  t10 = sp_ntt_add_simd0(t8, t9, p);
  t11 = sp_ntt_sub_simd0(t8, t9, p);

  t12 = sp_ntt_add_simd0(t10, t2, p);
  t13 = sp_ntt_sub_simd0(t10, t2, p);
  t14 = sp_ntt_add_simd0(t11, t3, p);
  t15 = sp_ntt_sub_simd0(t11, t3, p);

  p8 = sp_ntt_add_simd0(t12, t4, p);
  p9 = sp_ntt_sub_simd0(t12, t4, p);
  p10 = sp_ntt_add_simd0(t14, t5, p);
  p11 = sp_ntt_sub_simd0(t14, t5, p);
  p12 = sp_ntt_add_simd0(t13, t6, p);
  p13 = sp_ntt_sub_simd0(t13, t6, p);
  p14 = sp_ntt_add_simd0(t15, t7, p);
  p15 = sp_ntt_sub_simd0(t15, t7, p);

  sp_simd_store(p0, x + j0);
  sp_simd_store(p8, x + j1);
  sp_simd_store(p4, x + j2);
  sp_simd_store(p15,x + j3);
  sp_simd_store(p2, x + j4);
  sp_simd_store(p12,x + j5);
  sp_simd_store(p7, x + j6);
  sp_simd_store(p11,x + j7);
  sp_simd_store(p1, x + j8);
  sp_simd_store(p9, x + j9);
  sp_simd_store(p5, x + j10);
  sp_simd_store(p14,x + j11);
  sp_simd_store(p3, x + j12);
  sp_simd_store(p13,x + j13);
  sp_simd_store(p6, x + j14);
  sp_simd_store(p10,x + j15);
}

DECLARE_CORE_ROUTINES_SIMD(16)
