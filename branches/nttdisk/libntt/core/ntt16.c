#include "ntt-impl.h"

#define NC 18

static const uint8_t fixed_const[NC] = {1, 1, 1, 0, 1, 0, 0, 0, 1};

static const uint8_t *
ntt16_get_fixed_ntt_const(void)
{
  return fixed_const;
}

void
ntt16_init(spv_t out, sp_t p, sp_t d, 
	  sp_t primroot, sp_t order, sp_t perm)
{
  uint32_t i;
  sp_t w[16];
  sp_t inv2 = sp_inv(2, p, d);
  sp_t inv4 = sp_sqr(inv2, p, d);
  sp_t t0, t1, t2, t3, t4, t5;
  sp_t p0, p1, p2, p3;

  w[1] = sp_pow(primroot, order / 16, p, d);
  for (i = 2; i < 16; i++)
    w[i] = sp_mul(w[i-1], w[1], p, d);

  t3 = sp_sub(w[1], w[9], p);
  t2 = sp_sub(w[15], w[7], p);
  t1 = sp_sub(w[5], w[13], p);
  t0 = sp_sub(w[11], w[3], p);

  p0 = sp_add(t0, t1, p);
  p1 = sp_sub(t0, t1, p);
  p2 = sp_add(t2, t3, p);
  p3 = sp_sub(t2, t3, p);

  t0 = sp_sub(p0, p2, p);
  t1 = sp_sub(p1, p3, p);
  t2 = sp_neg(sp_add(p0, p2, p), p);
  t3 = sp_neg(sp_add(p1, p3, p), p);
  t4 = p2;
  t5 = p3;

  out[0] = 1;
  out[1] = 1;
  out[2] = 1;
  out[3] = w[4];
  out[4] = 1;
  out[5] = w[4];
  out[6] = sp_mul(inv2, sp_sub(w[2], w[6], p), p, d);
  out[7] = sp_mul(inv2, sp_add(w[2], w[6], p), p, d);
  out[8] = 1;
  out[9] = w[4];
  out[10] = out[6];
  out[11] = out[7];
  out[12] = sp_mul(inv4, t0, p, d);
  out[13] = sp_mul(inv4, t1, p, d);
  out[14] = sp_mul(inv4, t2, p, d);
  out[15] = sp_mul(inv4, t3, p, d);
  out[16] = sp_mul(inv4, t4, p, d);
  out[17] = sp_mul(inv4, t5, p, d);
}

static void 
ntt16_run_core(spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		sp_t p, spv_t ntt_const)
{
  sp_t x0, x1, x2, x3, x4, x5, x6, x7,
       x8, x9, x10, x11, x12, x13, x14, x15;
  sp_t t0, t1, t2, t3, t4, t5, t6, t7, 
       t8, t9, t10, t11, t12, t13, t14, t15, t16, t17;
  sp_t p0, p1, p2, p3, p4, p5, p6, p7,
       p8, p9, p10, p11, p12, p13, p14, p15, p16, p17;

  x0 = in[0 * istride];
  x1 = in[1 * istride];
  x2 = in[2 * istride];
  x3 = in[3 * istride];
  x4 = in[4 * istride];
  x5 = in[5 * istride];
  x6 = in[6 * istride];
  x7 = in[7 * istride];
  x8 = in[8 * istride];
  x9 = in[9 * istride];
  x10 = in[10 * istride];
  x11 = in[11 * istride];
  x12 = in[12 * istride];
  x13 = in[13 * istride];
  x14 = in[14 * istride];
  x15 = in[15 * istride];

  t0 = sp_ntt_add(x0, x8, p);
  t8 = sp_ntt_sub(x0, x8, p);
  t1 = sp_ntt_add(x1, x9, p);
  t9 = sp_ntt_sub(x1, x9, p);
  t2 = sp_ntt_add(x2, x10, p);
  t10 = sp_ntt_sub(x2, x10, p);
  t3 = sp_ntt_add(x3, x11, p);
  t11 = sp_ntt_sub(x3, x11, p);
  t4 = sp_ntt_add(x4, x12, p);
  t12 = sp_ntt_sub_partial(x4, x12, p);
  t5 = sp_ntt_add(x5, x13, p);
  t13 = sp_ntt_sub(x13, x5, p);
  t6 = sp_ntt_add(x6, x14, p);
  t14 = sp_ntt_sub(x14, x6, p);
  t7 = sp_ntt_add(x7, x15, p);
  t15 = sp_ntt_sub(x15, x7, p);

  p0 = sp_ntt_add(t0, t4, p);
  p4 = sp_ntt_sub(t0, t4, p);
  p1 = sp_ntt_add(t1, t5, p);
  p5 = sp_ntt_sub(t1, t5, p);
  p2 = sp_ntt_add(t2, t6, p);
  p6 = sp_ntt_sub_partial(t2, t6, p);
  p3 = sp_ntt_add(t3, t7, p);
  p7 = sp_ntt_sub(t3, t7, p);
  p8 = t8;
  p9 = t12;
  p10 = sp_ntt_add_partial(t10, t14, p);
  p11 = sp_ntt_sub_partial(t10, t14, p);
  p12 = sp_ntt_add(t9, t15, p);
  p13 = sp_ntt_sub(t9, t15, p);
  p14 = sp_ntt_add(t13, t11, p);
  p15 = sp_ntt_sub(t13, t11, p);

  t0 = sp_ntt_add(p0, p2, p);
  t1 = sp_ntt_sub(p0, p2, p);
  t2 = sp_ntt_add(p1, p3, p);
  t3 = sp_ntt_sub_partial(p1, p3, p);
  t4 = p4;
  t5 = p6;
  t6 = sp_ntt_sub_partial(p5, p7, p);
  t7 = sp_ntt_add_partial(p5, p7, p); 
  t8 = p8;
  t9 = p9;
  t10 = p10;
  t11 = p11;
  t12 = p12;
  t13 = p13;
  t14 = p14;
  t15 = p15;
  t16 = sp_ntt_add_partial(p12, p14, p);
  t17 = sp_ntt_add_partial(p13, p15, p);

  t3 = sp_ntt_mul(t3, ntt_const[3], ntt_const[NC+3], p);
  t5 = sp_ntt_mul(t5, ntt_const[5], ntt_const[NC+5], p);
  t6 = sp_ntt_mul(t6, ntt_const[6], ntt_const[NC+6], p);
  t7 = sp_ntt_mul(t7, ntt_const[7], ntt_const[NC+7], p);
  t9 = sp_ntt_mul(t9, ntt_const[9], ntt_const[NC+9], p);
  t10 = sp_ntt_mul(t10, ntt_const[10], ntt_const[NC+10], p);
  t11 = sp_ntt_mul(t11, ntt_const[11], ntt_const[NC+11], p);
  t12 = sp_ntt_mul(t12, ntt_const[12], ntt_const[NC+12], p);
  t13 = sp_ntt_mul(t13, ntt_const[13], ntt_const[NC+13], p);
  t14 = sp_ntt_mul(t14, ntt_const[14], ntt_const[NC+14], p);
  t15 = sp_ntt_mul(t15, ntt_const[15], ntt_const[NC+15], p);
  t16 = sp_ntt_mul(t16, ntt_const[16], ntt_const[NC+16], p);
  t17 = sp_ntt_mul(t17, ntt_const[17], ntt_const[NC+17], p);

  p0 = sp_ntt_add(t4, t5, p);
  p1 = sp_ntt_sub(t4, t5, p);
  p2 = sp_ntt_add(t6, t7, p);
  p3 = sp_ntt_sub(t6, t7, p);
  p4 = sp_ntt_add(p0, p2, p);
  p5 = sp_ntt_sub(p0, p2, p);
  p6 = sp_ntt_add(p1, p3, p);
  p7 = sp_ntt_sub(p1, p3, p);
  p0 = sp_ntt_add(t0, t2, p);
  p1 = sp_ntt_sub(t0, t2, p);
  p2 = sp_ntt_add(t1, t3, p);
  p3 = sp_ntt_sub(t1, t3, p);

  t0 = sp_ntt_add(t12, t16, p);
  t1 = sp_ntt_add(t13, t17, p);
  t2 = sp_ntt_add(t14, t16, p);
  t3 = sp_ntt_add(t15, t17, p);
  t7 = sp_ntt_add(t0, t1, p);
  t6 = sp_ntt_sub(t0, t1, p);
  t5 = sp_ntt_add(t2, t3, p);
  t4 = sp_ntt_sub(t2, t3, p); 

  t2 = sp_ntt_add(t10, t11, p);
  t3 = sp_ntt_sub(t10, t11, p);



  t10 = sp_ntt_add(t8, t9, p);
  t11 = sp_ntt_sub(t8, t9, p);

  t12 = sp_ntt_add(t10, t2, p);
  t13 = sp_ntt_sub(t10, t2, p);
  t14 = sp_ntt_add(t11, t3, p);
  t15 = sp_ntt_sub(t11, t3, p);

  p8 = sp_ntt_add(t12, t4, p);
  p9 = sp_ntt_sub(t12, t4, p);
  p10 = sp_ntt_add(t14, t5, p);
  p11 = sp_ntt_sub(t14, t5, p);
  p12 = sp_ntt_add(t13, t6, p);
  p13 = sp_ntt_sub(t13, t6, p);
  p14 = sp_ntt_add(t15, t7, p);
  p15 = sp_ntt_sub(t15, t7, p);

  out[ 0 * ostride] = p0;
  out[ 1 * ostride] = p8;
  out[ 2 * ostride] = p4;
  out[ 3 * ostride] = p15;
  out[ 4 * ostride] = p2;
  out[ 5 * ostride] = p12;
  out[ 6 * ostride] = p7;
  out[ 7 * ostride] = p11;
  out[ 8 * ostride] = p1;
  out[ 9 * ostride] = p9;
  out[10 * ostride] = p5;
  out[11 * ostride] = p14;
  out[12 * ostride] = p3;
  out[13 * ostride] = p13;
  out[14 * ostride] = p6;
  out[15 * ostride] = p10;
}

#ifdef HAVE_SIMD
static void
ntt16_run_core_simd(spv_t in, spv_size_t istride, spv_size_t idist,
		spv_t out, spv_size_t ostride, spv_size_t odist,
		sp_t p, spv_t ntt_const)
{
  sp_simd_t x0, x1, x2, x3, x4, x5, x6, x7,
            x8, x9, x10, x11, x12, x13, x14, x15;
  sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7, 
            t8, t9, t10, t11, t12, t13, t14, t15, t16, t17;
  sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7,
            p8, p9, p10, p11, p12, p13, p14, p15, p16, p17;

  x0 = sp_simd_gather(in + 0 * istride, idist);
  x1 = sp_simd_gather(in + 1 * istride, idist);
  x2 = sp_simd_gather(in + 2 * istride, idist);
  x3 = sp_simd_gather(in + 3 * istride, idist);
  x4 = sp_simd_gather(in + 4 * istride, idist);
  x5 = sp_simd_gather(in + 5 * istride, idist);
  x6 = sp_simd_gather(in + 6 * istride, idist);
  x7 = sp_simd_gather(in + 7 * istride, idist);
  x8 = sp_simd_gather(in + 8 * istride, idist);
  x9 = sp_simd_gather(in + 9 * istride, idist);
  x10 = sp_simd_gather(in + 10 * istride, idist);
  x11 = sp_simd_gather(in + 11 * istride, idist);
  x12 = sp_simd_gather(in + 12 * istride, idist);
  x13 = sp_simd_gather(in + 13 * istride, idist);
  x14 = sp_simd_gather(in + 14 * istride, idist);
  x15 = sp_simd_gather(in + 15 * istride, idist);

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

  sp_simd_scatter(p0,  out +  0 * ostride, odist);
  sp_simd_scatter(p8,  out +  1 * ostride, odist);
  sp_simd_scatter(p4,  out +  2 * ostride, odist);
  sp_simd_scatter(p15, out +  3 * ostride, odist);
  sp_simd_scatter(p2,  out +  4 * ostride, odist);
  sp_simd_scatter(p12, out +  5 * ostride, odist);
  sp_simd_scatter(p7,  out +  6 * ostride, odist);
  sp_simd_scatter(p11, out +  7 * ostride, odist);
  sp_simd_scatter(p1,  out +  8 * ostride, odist);
  sp_simd_scatter(p9,  out +  9 * ostride, odist);
  sp_simd_scatter(p5,  out + 10 * ostride, odist);
  sp_simd_scatter(p14, out + 11 * ostride, odist);
  sp_simd_scatter(p3,  out + 12 * ostride, odist);
  sp_simd_scatter(p13, out + 13 * ostride, odist);
  sp_simd_scatter(p6,  out + 14 * ostride, odist);
  sp_simd_scatter(p10, out + 15 * ostride, odist);
}
#endif


static void
ntt16_run(spv_t in, spv_size_t istride, spv_size_t idist,
    		spv_t out, spv_size_t ostride, spv_size_t odist,
    		spv_size_t num_transforms, sp_t p, spv_t ntt_const)
{
  spv_size_t i = 0;

#ifdef HAVE_SIMD
  spv_size_t num_simd = SP_SIMD_VSIZE * (num_transforms / SP_SIMD_VSIZE);

  for (; i < num_simd; i += SP_SIMD_VSIZE)
    ntt16_run_core_simd(in + i * idist, istride, idist, 
                        out + i * odist, ostride, odist, p, ntt_const);
#endif

  for (; i < num_transforms; i++)
    ntt16_run_core(in + i * idist, istride, 
                out + i * odist, ostride, p, ntt_const);
}


static void
ntt16_twiddle_run_core(spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		spv_t w, sp_t p, spv_t ntt_const)
{
  sp_t x0, x1, x2, x3, x4, x5, x6, x7,
       x8, x9, x10, x11, x12, x13, x14, x15;
  sp_t t0, t1, t2, t3, t4, t5, t6, t7, 
       t8, t9, t10, t11, t12, t13, t14, t15, t16, t17;
  sp_t p0, p1, p2, p3, p4, p5, p6, p7,
       p8, p9, p10, p11, p12, p13, p14, p15, p16, p17;

  x0 = in[0 * istride];
  x1 = in[1 * istride];
  x2 = in[2 * istride];
  x3 = in[3 * istride];
  x4 = in[4 * istride];
  x5 = in[5 * istride];
  x6 = in[6 * istride];
  x7 = in[7 * istride];
  x8 = in[8 * istride];
  x9 = in[9 * istride];
  x10 = in[10 * istride];
  x11 = in[11 * istride];
  x12 = in[12 * istride];
  x13 = in[13 * istride];
  x14 = in[14 * istride];
  x15 = in[15 * istride];

  t0 = sp_ntt_add(x0, x8, p);
  t8 = sp_ntt_sub(x0, x8, p);
  t1 = sp_ntt_add(x1, x9, p);
  t9 = sp_ntt_sub(x1, x9, p);
  t2 = sp_ntt_add(x2, x10, p);
  t10 = sp_ntt_sub(x2, x10, p);
  t3 = sp_ntt_add(x3, x11, p);
  t11 = sp_ntt_sub(x3, x11, p);
  t4 = sp_ntt_add(x4, x12, p);
  t12 = sp_ntt_sub_partial(x4, x12, p);
  t5 = sp_ntt_add(x5, x13, p);
  t13 = sp_ntt_sub(x13, x5, p);
  t6 = sp_ntt_add(x6, x14, p);
  t14 = sp_ntt_sub(x14, x6, p);
  t7 = sp_ntt_add(x7, x15, p);
  t15 = sp_ntt_sub(x15, x7, p);

  p0 = sp_ntt_add(t0, t4, p);
  p4 = sp_ntt_sub(t0, t4, p);
  p1 = sp_ntt_add(t1, t5, p);
  p5 = sp_ntt_sub(t1, t5, p);
  p2 = sp_ntt_add(t2, t6, p);
  p6 = sp_ntt_sub_partial(t2, t6, p);
  p3 = sp_ntt_add(t3, t7, p);
  p7 = sp_ntt_sub(t3, t7, p);
  p8 = t8;
  p9 = t12;
  p10 = sp_ntt_add_partial(t10, t14, p);
  p11 = sp_ntt_sub_partial(t10, t14, p);
  p12 = sp_ntt_add(t9, t15, p);
  p13 = sp_ntt_sub(t9, t15, p);
  p14 = sp_ntt_add(t13, t11, p);
  p15 = sp_ntt_sub(t13, t11, p);

  t0 = sp_ntt_add(p0, p2, p);
  t1 = sp_ntt_sub(p0, p2, p);
  t2 = sp_ntt_add(p1, p3, p);
  t3 = sp_ntt_sub_partial(p1, p3, p);
  t4 = p4;
  t5 = p6;
  t6 = sp_ntt_sub_partial(p5, p7, p);
  t7 = sp_ntt_add_partial(p5, p7, p); 
  t8 = p8;
  t9 = p9;
  t10 = p10;
  t11 = p11;
  t12 = p12;
  t13 = p13;
  t14 = p14;
  t15 = p15;
  t16 = sp_ntt_add_partial(p12, p14, p);
  t17 = sp_ntt_add_partial(p13, p15, p);

  t3 = sp_ntt_mul(t3, ntt_const[3], ntt_const[NC+3], p);
  t5 = sp_ntt_mul(t5, ntt_const[5], ntt_const[NC+5], p);
  t6 = sp_ntt_mul(t6, ntt_const[6], ntt_const[NC+6], p);
  t7 = sp_ntt_mul(t7, ntt_const[7], ntt_const[NC+7], p);
  t9 = sp_ntt_mul(t9, ntt_const[9], ntt_const[NC+9], p);
  t10 = sp_ntt_mul(t10, ntt_const[10], ntt_const[NC+10], p);
  t11 = sp_ntt_mul(t11, ntt_const[11], ntt_const[NC+11], p);
  t12 = sp_ntt_mul(t12, ntt_const[12], ntt_const[NC+12], p);
  t13 = sp_ntt_mul(t13, ntt_const[13], ntt_const[NC+13], p);
  t14 = sp_ntt_mul(t14, ntt_const[14], ntt_const[NC+14], p);
  t15 = sp_ntt_mul(t15, ntt_const[15], ntt_const[NC+15], p);
  t16 = sp_ntt_mul(t16, ntt_const[16], ntt_const[NC+16], p);
  t17 = sp_ntt_mul(t17, ntt_const[17], ntt_const[NC+17], p);

  p0 = sp_ntt_add(t4, t5, p);
  p1 = sp_ntt_sub(t4, t5, p);
  p2 = sp_ntt_add(t6, t7, p);
  p3 = sp_ntt_sub(t6, t7, p);
  p4 = sp_ntt_add_partial(p0, p2, p);
  p5 = sp_ntt_sub_partial(p0, p2, p);
  p6 = sp_ntt_add_partial(p1, p3, p);
  p7 = sp_ntt_sub_partial(p1, p3, p);
  p0 = sp_ntt_add(t0, t2, p);
  p1 = sp_ntt_sub_partial(t0, t2, p);
  p2 = sp_ntt_add_partial(t1, t3, p);
  p3 = sp_ntt_sub_partial(t1, t3, p);

  t0 = sp_ntt_add(t12, t16, p);
  t1 = sp_ntt_add(t13, t17, p);
  t2 = sp_ntt_add(t14, t16, p);
  t3 = sp_ntt_add(t15, t17, p);
  t7 = sp_ntt_add(t0, t1, p);
  t6 = sp_ntt_sub(t0, t1, p);
  t5 = sp_ntt_add(t2, t3, p);
  t4 = sp_ntt_sub(t2, t3, p); 

  t2 = sp_ntt_add(t10, t11, p);
  t3 = sp_ntt_sub(t10, t11, p);



  t10 = sp_ntt_add(t8, t9, p);
  t11 = sp_ntt_sub(t8, t9, p);

  t12 = sp_ntt_add(t10, t2, p);
  t13 = sp_ntt_sub(t10, t2, p);
  t14 = sp_ntt_add(t11, t3, p);
  t15 = sp_ntt_sub(t11, t3, p);

  p8 = sp_ntt_add_partial(t12, t4, p);
  p9 = sp_ntt_sub_partial(t12, t4, p);
  p10 = sp_ntt_add_partial(t14, t5, p);
  p11 = sp_ntt_sub_partial(t14, t5, p);
  p12 = sp_ntt_add_partial(t13, t6, p);
  p13 = sp_ntt_sub_partial(t13, t6, p);
  p14 = sp_ntt_add_partial(t15, t7, p);
  p15 = sp_ntt_sub_partial(t15, t7, p);

  p8  = sp_ntt_mul(p8, w[0], w[1], p);
  p4  = sp_ntt_mul(p4, w[2], w[3], p);
  p15 = sp_ntt_mul(p15, w[4], w[5], p);
  p2  = sp_ntt_mul(p2, w[6], w[7], p);
  p12 = sp_ntt_mul(p12, w[8], w[9], p);
  p7  = sp_ntt_mul(p7, w[10], w[11], p);
  p11 = sp_ntt_mul(p11, w[12], w[13], p);
  p1  = sp_ntt_mul(p1, w[14], w[15], p);
  p9  = sp_ntt_mul(p9, w[16], w[17], p);
  p5  = sp_ntt_mul(p5, w[18], w[19], p);
  p14 = sp_ntt_mul(p14, w[20], w[21], p);
  p3  = sp_ntt_mul(p3, w[22], w[23], p);
  p13 = sp_ntt_mul(p13, w[24], w[25], p);
  p6  = sp_ntt_mul(p6, w[26], w[27], p);
  p10 = sp_ntt_mul(p10, w[28], w[29], p);

  out[ 0 * ostride] = p0;
  out[ 1 * ostride] = p8;
  out[ 2 * ostride] = p4;
  out[ 3 * ostride] = p15;
  out[ 4 * ostride] = p2;
  out[ 5 * ostride] = p12;
  out[ 6 * ostride] = p7;
  out[ 7 * ostride] = p11;
  out[ 8 * ostride] = p1;
  out[ 9 * ostride] = p9;
  out[10 * ostride] = p5;
  out[11 * ostride] = p14;
  out[12 * ostride] = p3;
  out[13 * ostride] = p13;
  out[14 * ostride] = p6;
  out[15 * ostride] = p10;
}


#ifdef HAVE_SIMD
static void
ntt16_twiddle_run_core_simd(
        spv_t in, spv_size_t istride, spv_size_t idist,
		spv_t out, spv_size_t ostride, spv_size_t odist,
		sp_simd_t *w, sp_t p, spv_t ntt_const)
{
  sp_simd_t x0, x1, x2, x3, x4, x5, x6, x7,
            x8, x9, x10, x11, x12, x13, x14, x15;
  sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7, 
            t8, t9, t10, t11, t12, t13, t14, t15, t16, t17;
  sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7,
            p8, p9, p10, p11, p12, p13, p14, p15, p16, p17;

  x0 = sp_simd_gather(in + 0 * istride, idist);
  x1 = sp_simd_gather(in + 1 * istride, idist);
  x2 = sp_simd_gather(in + 2 * istride, idist);
  x3 = sp_simd_gather(in + 3 * istride, idist);
  x4 = sp_simd_gather(in + 4 * istride, idist);
  x5 = sp_simd_gather(in + 5 * istride, idist);
  x6 = sp_simd_gather(in + 6 * istride, idist);
  x7 = sp_simd_gather(in + 7 * istride, idist);
  x8 = sp_simd_gather(in + 8 * istride, idist);
  x9 = sp_simd_gather(in + 9 * istride, idist);
  x10 = sp_simd_gather(in + 10 * istride, idist);
  x11 = sp_simd_gather(in + 11 * istride, idist);
  x12 = sp_simd_gather(in + 12 * istride, idist);
  x13 = sp_simd_gather(in + 13 * istride, idist);
  x14 = sp_simd_gather(in + 14 * istride, idist);
  x15 = sp_simd_gather(in + 15 * istride, idist);

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

  sp_simd_scatter(p0, out +  0 * ostride, odist);
  sp_simd_scatter(p8, out +  1 * ostride, odist);
  sp_simd_scatter(p4, out +  2 * ostride, odist);
  sp_simd_scatter(p15, out +  3 * ostride, odist);
  sp_simd_scatter(p2, out +  4 * ostride, odist);
  sp_simd_scatter(p12, out +  5 * ostride, odist);
  sp_simd_scatter(p7, out +  6 * ostride, odist);
  sp_simd_scatter(p11, out +  7 * ostride, odist);
  sp_simd_scatter(p1, out +  8 * ostride, odist);
  sp_simd_scatter(p9, out +  9 * ostride, odist);
  sp_simd_scatter(p5, out + 10 * ostride, odist);
  sp_simd_scatter(p14, out + 11 * ostride, odist);
  sp_simd_scatter(p3, out + 12 * ostride, odist);
  sp_simd_scatter(p13, out + 13 * ostride, odist);
  sp_simd_scatter(p6, out + 14 * ostride, odist);
  sp_simd_scatter(p10, out + 15 * ostride, odist);
}
#endif

static void
ntt16_twiddle_run(spv_t in, spv_size_t istride, spv_size_t idist,
    			spv_t out, spv_size_t ostride, spv_size_t odist,
    			spv_t w, spv_size_t num_transforms, sp_t p, spv_t ntt_const)
{
  spv_size_t i = 0, j = 0;

#ifdef HAVE_SIMD
  spv_size_t num_simd = SP_SIMD_VSIZE * (num_transforms / SP_SIMD_VSIZE);

  for (; i < num_simd; i += SP_SIMD_VSIZE,
		  	j += 2*(16-1)*SP_SIMD_VSIZE)
    ntt16_twiddle_run_core_simd(
		in + i * idist, istride, idist,
		out + i * odist, ostride, odist,
		(sp_simd_t *)(w + j), p, ntt_const);
#endif

  for (; i < num_transforms; i++, j += 2*(16-1))
    ntt16_twiddle_run_core(in + i * idist, istride, 
			out + i * odist, ostride,
			w + j, p, ntt_const);
}



static void
ntt16_pfa_run_core(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t n,
	  sp_t p, spv_t ntt_const)
{
  spv_size_t j0, j1, j2, j3, j4, j5, j6, j7,
             j8, j9, j10, j11, j12, j13, j14, j15;
  sp_t x0, x1, x2, x3, x4, x5, x6, x7,
       x8, x9, x10, x11, x12, x13, x14, x15;
  sp_t t0, t1, t2, t3, t4, t5, t6, t7, 
       t8, t9, t10, t11, t12, t13, t14, t15, t16, t17;
  sp_t p0, p1, p2, p3, p4, p5, p6, p7,
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

  x0 = x[j0];
  x1 = x[j1];
  x2 = x[j2];
  x3 = x[j3];
  x4 = x[j4];
  x5 = x[j5];
  x6 = x[j6];
  x7 = x[j7];
  x8 = x[j8];
  x9 = x[j9];
  x10 = x[j10];
  x11 = x[j11];
  x12 = x[j12];
  x13 = x[j13];
  x14 = x[j14];
  x15 = x[j15];

  t0 = sp_ntt_add(x0, x8, p);
  t8 = sp_ntt_sub(x0, x8, p);
  t1 = sp_ntt_add(x1, x9, p);
  t9 = sp_ntt_sub(x1, x9, p);
  t2 = sp_ntt_add(x2, x10, p);
  t10 = sp_ntt_sub(x2, x10, p);
  t3 = sp_ntt_add(x3, x11, p);
  t11 = sp_ntt_sub(x3, x11, p);
  t4 = sp_ntt_add(x4, x12, p);
  t12 = sp_ntt_sub_partial(x4, x12, p);
  t5 = sp_ntt_add(x5, x13, p);
  t13 = sp_ntt_sub(x13, x5, p);
  t6 = sp_ntt_add(x6, x14, p);
  t14 = sp_ntt_sub(x14, x6, p);
  t7 = sp_ntt_add(x7, x15, p);
  t15 = sp_ntt_sub(x15, x7, p);

  p0 = sp_ntt_add(t0, t4, p);
  p4 = sp_ntt_sub(t0, t4, p);
  p1 = sp_ntt_add(t1, t5, p);
  p5 = sp_ntt_sub(t1, t5, p);
  p2 = sp_ntt_add(t2, t6, p);
  p6 = sp_ntt_sub_partial(t2, t6, p);
  p3 = sp_ntt_add(t3, t7, p);
  p7 = sp_ntt_sub(t3, t7, p);
  p8 = t8;
  p9 = t12;
  p10 = sp_ntt_add_partial(t10, t14, p);
  p11 = sp_ntt_sub_partial(t10, t14, p);
  p12 = sp_ntt_add(t9, t15, p);
  p13 = sp_ntt_sub(t9, t15, p);
  p14 = sp_ntt_add(t13, t11, p);
  p15 = sp_ntt_sub(t13, t11, p);

  t0 = sp_ntt_add(p0, p2, p);
  t1 = sp_ntt_sub(p0, p2, p);
  t2 = sp_ntt_add(p1, p3, p);
  t3 = sp_ntt_sub_partial(p1, p3, p);
  t4 = p4;
  t5 = p6;
  t6 = sp_ntt_sub_partial(p5, p7, p);
  t7 = sp_ntt_add_partial(p5, p7, p); 
  t8 = p8;
  t9 = p9;
  t10 = p10;
  t11 = p11;
  t12 = p12;
  t13 = p13;
  t14 = p14;
  t15 = p15;
  t16 = sp_ntt_add_partial(p12, p14, p);
  t17 = sp_ntt_add_partial(p13, p15, p);

  t3 = sp_ntt_mul(t3, ntt_const[3], ntt_const[NC+3], p);
  t5 = sp_ntt_mul(t5, ntt_const[5], ntt_const[NC+5], p);
  t6 = sp_ntt_mul(t6, ntt_const[6], ntt_const[NC+6], p);
  t7 = sp_ntt_mul(t7, ntt_const[7], ntt_const[NC+7], p);
  t9 = sp_ntt_mul(t9, ntt_const[9], ntt_const[NC+9], p);
  t10 = sp_ntt_mul(t10, ntt_const[10], ntt_const[NC+10], p);
  t11 = sp_ntt_mul(t11, ntt_const[11], ntt_const[NC+11], p);
  t12 = sp_ntt_mul(t12, ntt_const[12], ntt_const[NC+12], p);
  t13 = sp_ntt_mul(t13, ntt_const[13], ntt_const[NC+13], p);
  t14 = sp_ntt_mul(t14, ntt_const[14], ntt_const[NC+14], p);
  t15 = sp_ntt_mul(t15, ntt_const[15], ntt_const[NC+15], p);
  t16 = sp_ntt_mul(t16, ntt_const[16], ntt_const[NC+16], p);
  t17 = sp_ntt_mul(t17, ntt_const[17], ntt_const[NC+17], p);

  p0 = sp_ntt_add(t4, t5, p);
  p1 = sp_ntt_sub(t4, t5, p);
  p2 = sp_ntt_add(t6, t7, p);
  p3 = sp_ntt_sub(t6, t7, p);
  p4 = sp_ntt_add(p0, p2, p);
  p5 = sp_ntt_sub(p0, p2, p);
  p6 = sp_ntt_add(p1, p3, p);
  p7 = sp_ntt_sub(p1, p3, p);
  p0 = sp_ntt_add(t0, t2, p);
  p1 = sp_ntt_sub(t0, t2, p);
  p2 = sp_ntt_add(t1, t3, p);
  p3 = sp_ntt_sub(t1, t3, p);

  t0 = sp_ntt_add(t12, t16, p);
  t1 = sp_ntt_add(t13, t17, p);
  t2 = sp_ntt_add(t14, t16, p);
  t3 = sp_ntt_add(t15, t17, p);
  t7 = sp_ntt_add(t0, t1, p);
  t6 = sp_ntt_sub(t0, t1, p);
  t5 = sp_ntt_add(t2, t3, p);
  t4 = sp_ntt_sub(t2, t3, p); 

  t2 = sp_ntt_add(t10, t11, p);
  t3 = sp_ntt_sub(t10, t11, p);



  t10 = sp_ntt_add(t8, t9, p);
  t11 = sp_ntt_sub(t8, t9, p);

  t12 = sp_ntt_add(t10, t2, p);
  t13 = sp_ntt_sub(t10, t2, p);
  t14 = sp_ntt_add(t11, t3, p);
  t15 = sp_ntt_sub(t11, t3, p);

  p8 = sp_ntt_add(t12, t4, p);
  p9 = sp_ntt_sub(t12, t4, p);
  p10 = sp_ntt_add(t14, t5, p);
  p11 = sp_ntt_sub(t14, t5, p);
  p12 = sp_ntt_add(t13, t6, p);
  p13 = sp_ntt_sub(t13, t6, p);
  p14 = sp_ntt_add(t15, t7, p);
  p15 = sp_ntt_sub(t15, t7, p);

  x[j0] = p0;
  x[j1] = p8;
  x[j2] = p4;
  x[j3] = p15;
  x[j4] = p2;
  x[j5] = p12;
  x[j6] = p7;
  x[j7] = p11;
  x[j8] = p1;
  x[j9] = p9;
  x[j10] = p5;
  x[j11] = p14;
  x[j12] = p3;
  x[j13] = p13;
  x[j14] = p6;
  x[j15] = p10;
}

#ifdef HAVE_SIMD
static void
ntt16_pfa_run_core_simd(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t inc2, spv_size_t n,
	  sp_t p, spv_t ntt_const)
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

  x0 = sp_simd_pfa_gather(x, j0, inc2, n);
  x1 = sp_simd_pfa_gather(x, j1, inc2, n);
  x2 = sp_simd_pfa_gather(x, j2, inc2, n);
  x3 = sp_simd_pfa_gather(x, j3, inc2, n);
  x4 = sp_simd_pfa_gather(x, j4, inc2, n);
  x5 = sp_simd_pfa_gather(x, j5, inc2, n);
  x6 = sp_simd_pfa_gather(x, j6, inc2, n);
  x7 = sp_simd_pfa_gather(x, j7, inc2, n);
  x8 = sp_simd_pfa_gather(x, j8, inc2, n);
  x9 = sp_simd_pfa_gather(x, j9, inc2, n);
  x10 = sp_simd_pfa_gather(x, j10, inc2, n);
  x11 = sp_simd_pfa_gather(x, j11, inc2, n);
  x12 = sp_simd_pfa_gather(x, j12, inc2, n);
  x13 = sp_simd_pfa_gather(x, j13, inc2, n);
  x14 = sp_simd_pfa_gather(x, j14, inc2, n);
  x15 = sp_simd_pfa_gather(x, j15, inc2, n);

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

  sp_simd_pfa_scatter(p0, x, j0, inc2, n);
  sp_simd_pfa_scatter(p8, x, j1, inc2, n);
  sp_simd_pfa_scatter(p4, x, j2, inc2, n);
  sp_simd_pfa_scatter(p15,x, j3, inc2, n);
  sp_simd_pfa_scatter(p2, x, j4, inc2, n);
  sp_simd_pfa_scatter(p12,x, j5, inc2, n);
  sp_simd_pfa_scatter(p7, x, j6, inc2, n);
  sp_simd_pfa_scatter(p11,x, j7, inc2, n);
  sp_simd_pfa_scatter(p1, x, j8, inc2, n);
  sp_simd_pfa_scatter(p9, x, j9, inc2, n);
  sp_simd_pfa_scatter(p5, x, j10, inc2, n);
  sp_simd_pfa_scatter(p14,x, j11, inc2, n);
  sp_simd_pfa_scatter(p3, x, j12, inc2, n);
  sp_simd_pfa_scatter(p13,x, j13, inc2, n);
  sp_simd_pfa_scatter(p6, x, j14, inc2, n);
  sp_simd_pfa_scatter(p10,x, j15, inc2, n);
}
#endif

static void
ntt16_pfa_run(spv_t x, spv_size_t cofactor,
	  sp_t p, spv_t ntt_const)
{
  spv_size_t i = 0;
  spv_size_t incstart = 0;
  spv_size_t n = 16 * cofactor;
  spv_size_t inc = cofactor;
  spv_size_t inc2 = 16;

#ifdef HAVE_SIMD
  spv_size_t num_simd = SP_SIMD_VSIZE * (cofactor / SP_SIMD_VSIZE);

  for (i = 0; i < num_simd; i += SP_SIMD_VSIZE)
    {
      ntt16_pfa_run_core_simd(x, incstart, inc, inc2, n, p, ntt_const);
      incstart += SP_SIMD_VSIZE * inc2;
    }
#endif

  for (; i < cofactor; i++, incstart += inc2)
    ntt16_pfa_run_core(x, incstart, inc, n, p, ntt_const);

}

const nttconfig_t ntt16_config = 
{
  16,
  NC,
  ntt16_get_fixed_ntt_const,
  ntt16_init,
  ntt16_run,
  ntt16_pfa_run,
  ntt16_twiddle_run
};

