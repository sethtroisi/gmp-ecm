#include "ntt/ntt-impl-scalar.h"

#define NC 18

static const uint8_t ntt16_fixed_const[NC] = {1, 1, 1, 0, 1, 0, 0, 0, 1};

void
X(ntt16_init)(spv_t out, sp_t p, sp_t d, 
	  sp_t primroot, sp_t order, sp_t perm)
{
  X(nttdata_init_generic)(&X(ntt16_config), out, p, d, primroot, order, perm);
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

DECLARE_CORE_ROUTINES(16)
