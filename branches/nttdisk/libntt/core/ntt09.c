#include "ntt-impl.h"

static uint32_t 
ntt9_get_num_const(void)
{
  return 11;
}

void
ntt9_init(spv_t out, sp_t p, sp_t d, 
	  sp_t primroot, sp_t order)
{
  uint32_t i;
  sp_t t1, t2, t3, t4, t5, t6, t7, t8;
  sp_t h1, h2, h3, h4, h5, h6, h7, h8;
  sp_t w1, w2, w3, w4, w5, w6;
  sp_t w[9];
  sp_t inv6 = sp_inv(6, p, d);
  sp_t inv2 = sp_inv(2, p, d);
    
  w[1] = sp_pow(primroot, order / 9, p, d);
  for (i = 2; i < 9; i++)
    w[i] = sp_mul(w[i-1], w[1], p, d);

  w1 = w[2];
  w2 = w[4];
  w3 = w[8];
  w4 = w[7];
  w5 = w[5];
  w6 = w[1];

  t3 = sp_add(w1, w1, p);
  t3 = sp_add(t3, w2, p);
  t3 = sp_sub(t3, w3, p);
  t3 = sp_sub(t3, w4, p);
  t3 = sp_sub(t3, w4, p);
  t3 = sp_sub(t3, w5, p);
  t3 = sp_add(t3, w6, p);

  t4 = sp_add(w1, w2, p);
  t4 = sp_add(t4, w2, p);
  t4 = sp_add(t4, w3, p);
  t4 = sp_sub(t4, w4, p);
  t4 = sp_sub(t4, w5, p);
  t4 = sp_sub(t4, w5, p);
  t4 = sp_sub(t4, w6, p);

  t5 = sp_add(w1, w1, p);
  t5 = sp_sub(t5, w2, p);
  t5 = sp_sub(t5, w3, p);
  t5 = sp_add(t5, w4, p);
  t5 = sp_add(t5, w4, p);
  t5 = sp_sub(t5, w5, p);
  t5 = sp_sub(t5, w6, p);

  t6 = sp_neg(w1, p);
  t6 = sp_add(t6, w2, p);
  t6 = sp_add(t6, w2, p);
  t6 = sp_sub(t6, w3, p);
  t6 = sp_sub(t6, w4, p);
  t6 = sp_add(t6, w5, p);
  t6 = sp_add(t6, w5, p);
  t6 = sp_sub(t6, w6, p);

  h1 = sp_add(w[6], w[3], p);
  h2 = sp_sub(w[6], w[3], p);

  h3 = sp_sub(t3, t4, p);
  h4 = sp_neg(t3, p); 
  h5 = t4;
  h6 = sp_sub(t5, t6, p);
  h7 = sp_neg(sp_add(t5, t6, p), p);
  h7 = sp_sub(h7, t6, p);
  h8 = t6;

  out[0] = 1;

  out[1] = sp_mul(h1, inv2, p, d);
  out[2] = sp_mul(h2, inv2, p, d);

  out[3] = sp_sub(out[1], 1, p);
  out[4] = out[2];

  out[5] = sp_mul(h3, inv6, p, d);
  out[6] = sp_mul(h4, inv6, p, d);
  out[7] = sp_mul(h5, inv6, p, d);
  out[8] = sp_mul(h6, inv6, p, d);
  out[9] = sp_mul(h7, inv6, p, d);
  out[10] = sp_mul(h8, inv6, p, d);
}

static void
ntt9_run(spv_t x, spv_size_t stride,
	  sp_t p, sp_t d, spv_t ntt_const)
{
  sp_t x0, x1, x2, x3, x4, x5, x6;
  sp_t     t1, t2, t3, t4, t5, t6, t7, t8;
  sp_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

  sp_t x0e, x1e, t0e, t1e, p0e, p1e, p2e;


  x0 = x[0 * stride];
  x1 = x[1 * stride];
  x2 = x[2 * stride];
  x0e = x[3 * stride];
  x3 = x[4 * stride];
  x6 = x[5 * stride];
  x1e = x[6 * stride];
  x5 = x[7 * stride];
  x4 = x[8 * stride];

  t0e = sp_add(x0e, x1e, p);
  t1e = sp_sub(x0e, x1e, p);

  p1 = sp_add(x1, x3, p);
  p1 = sp_add(p1, x5, p);
  p2 = sp_add(x2, x4, p);
  p2 = sp_add(p2, x6, p);
  p3 = sp_sub(x1, x5, p);
  p4 = sp_sub(x2, x6, p);
  p5 = sp_sub(x3, x5, p);
  p6 = sp_sub(x4, x6, p);

  t1 = sp_add(p1, p2, p);
  t2 = sp_sub(p1, p2, p);
  t3 = sp_sub(p3, p5, p);
  t5 = sp_add(t3, p6, p);
  t3 = sp_sub(t3, p6, p);
  t4 = sp_add(p4, p5, p);
  t6 = sp_sub(p4, p5, p);

  p0e = sp_add(x0, t0e, p);
  p0 = t1;
  p1 = t1;
  p2 = t2;
  p3 = t3;
  p4 = t4;
  p5 = sp_add(t3, t4, p);
  p6 = t5;
  p7 = t6;
  p8 = sp_add(t5, t6, p);

  p1 = sp_mul(p1, ntt_const[1], p, d);
  p2 = sp_mul(p2, ntt_const[2], p, d);
  t0e = sp_mul(t0e, ntt_const[3], p, d);
  t1e = sp_mul(t1e, ntt_const[4], p, d);
  p3 = sp_mul(p3, ntt_const[5], p, d);
  p4 = sp_mul(p4, ntt_const[6], p, d);
  p5 = sp_mul(p5, ntt_const[7], p, d);
  p6 = sp_mul(p6, ntt_const[8], p, d);
  p7 = sp_mul(p7, ntt_const[9], p, d);
  p8 = sp_mul(p8, ntt_const[10], p, d);

  t0e = sp_add(t0e, p0e, p);
  t1 = sp_add(p1, p2, p);
  t2 = sp_sub(p1, p2, p);
  t3 = sp_add(p3, p5, p);
  t4 = sp_add(p4, p5, p);
  t5 = sp_add(p6, p8, p);
  t6 = sp_add(p7, p8, p);

  p1e = sp_add(t0e, t1e, p);
  p2e = sp_sub(t0e, t1e, p);
  p3 = sp_add(t3, t5, p);
  p4 = sp_add(t4, t6, p);
  p5 = sp_sub(t4, t6, p);
  p5 = sp_sub(p5, p3, p);
  p6 = sp_sub(t5, t3, p);

  p0 = sp_add(p0, p0e, p);
  t1 = sp_add(t1, p0e, p);
  t2 = sp_add(t2, p0e, p);
  t3 = sp_add(p3, p1e, p);
  t4 = sp_add(p4, p2e, p);
  t5 = sp_add(p5, p1e, p);
  t6 = sp_add(p6, p2e, p);
  t7 = sp_add(p3, p5, p);
  t7 = sp_sub(p1e, t7, p);
  t8 = sp_add(p4, p6, p);
  t8 = sp_sub(p2e, t8, p);

  x[0 * stride] = p0;
  x[1 * stride] = t8;
  x[2 * stride] = t3;
  x[3 * stride] = t2;
  x[4 * stride] = t4;
  x[5 * stride] = t7;
  x[6 * stride] = t1;
  x[7 * stride] = t6;
  x[8 * stride] = t5;
}

static void
ntt9_pfa_run_core(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t n,
	  sp_t p, sp_t d, spv_t ntt_const)
{
  spv_size_t j0, j1, j2, j3, j4, j5, j6, j7, j8;
  sp_t x0, x1, x2, x3, x4, x5, x6;
  sp_t     t1, t2, t3, t4, t5, t6, t7, t8;
  sp_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

  sp_t x0e, x1e, t0e, t1e, p0e, p1e, p2e;

  j0 = start;
  j1 = sp_array_inc(j0, inc, n);
  j2 = sp_array_inc(j0, 2 * inc, n);
  j3 = sp_array_inc(j0, 3 * inc, n);
  j4 = sp_array_inc(j0, 4 * inc, n);
  j5 = sp_array_inc(j0, 5 * inc, n);
  j6 = sp_array_inc(j0, 6 * inc, n);
  j7 = sp_array_inc(j0, 7 * inc, n);
  j8 = sp_array_inc(j0, 8 * inc, n);

  x0 = x[j0];
  x1 = x[j1];
  x2 = x[j2];
  x0e = x[j3];
  x3 = x[j4];
  x6 = x[j5];
  x1e = x[j6];
  x5 = x[j7];
  x4 = x[j8];

  t0e = sp_add(x0e, x1e, p);
  t1e = sp_sub(x0e, x1e, p);

  p1 = sp_add(x1, x3, p);
  p1 = sp_add(p1, x5, p);
  p2 = sp_add(x2, x4, p);
  p2 = sp_add(p2, x6, p);
  p3 = sp_sub(x1, x5, p);
  p4 = sp_sub(x2, x6, p);
  p5 = sp_sub(x3, x5, p);
  p6 = sp_sub(x4, x6, p);

  t1 = sp_add(p1, p2, p);
  t2 = sp_sub(p1, p2, p);
  t3 = sp_sub(p3, p5, p);
  t5 = sp_add(t3, p6, p);
  t3 = sp_sub(t3, p6, p);
  t4 = sp_add(p4, p5, p);
  t6 = sp_sub(p4, p5, p);

  p0e = sp_add(x0, t0e, p);
  p0 = t1;
  p1 = t1;
  p2 = t2;
  p3 = t3;
  p4 = t4;
  p5 = sp_add(t3, t4, p);
  p6 = t5;
  p7 = t6;
  p8 = sp_add(t5, t6, p);

  p1 = sp_mul(p1, ntt_const[1], p, d);
  p2 = sp_mul(p2, ntt_const[2], p, d);
  t0e = sp_mul(t0e, ntt_const[3], p, d);
  t1e = sp_mul(t1e, ntt_const[4], p, d);
  p3 = sp_mul(p3, ntt_const[5], p, d);
  p4 = sp_mul(p4, ntt_const[6], p, d);
  p5 = sp_mul(p5, ntt_const[7], p, d);
  p6 = sp_mul(p6, ntt_const[8], p, d);
  p7 = sp_mul(p7, ntt_const[9], p, d);
  p8 = sp_mul(p8, ntt_const[10], p, d);

  t0e = sp_add(t0e, p0e, p);
  t1 = sp_add(p1, p2, p);
  t2 = sp_sub(p1, p2, p);
  t3 = sp_add(p3, p5, p);
  t4 = sp_add(p4, p5, p);
  t5 = sp_add(p6, p8, p);
  t6 = sp_add(p7, p8, p);

  p1e = sp_add(t0e, t1e, p);
  p2e = sp_sub(t0e, t1e, p);
  p3 = sp_add(t3, t5, p);
  p4 = sp_add(t4, t6, p);
  p5 = sp_sub(t4, t6, p);
  p5 = sp_sub(p5, p3, p);
  p6 = sp_sub(t5, t3, p);

  p0 = sp_add(p0, p0e, p);
  t1 = sp_add(t1, p0e, p);
  t2 = sp_add(t2, p0e, p);
  t3 = sp_add(p3, p1e, p);
  t4 = sp_add(p4, p2e, p);
  t5 = sp_add(p5, p1e, p);
  t6 = sp_add(p6, p2e, p);
  t7 = sp_add(p3, p5, p);
  t7 = sp_sub(p1e, t7, p);
  t8 = sp_add(p4, p6, p);
  t8 = sp_sub(p2e, t8, p);

  x[j0] = p0;
  x[j1] = t8;
  x[j2] = t3;
  x[j3] = t2;
  x[j4] = t4;
  x[j5] = t7;
  x[j6] = t1;
  x[j7] = t6;
  x[j8] = t5;
}

#ifdef HAVE_SSE2
static void
ntt9_pfa_run_core_simd(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t inc2, spv_size_t n,
	  sp_t p, sp_t d, spv_t ntt_const)
{
  spv_size_t j0, j1, j2, j3, j4, j5, j6, j7, j8;
  sp_simd_t x0, x1, x2, x3, x4, x5, x6;
  sp_simd_t     t1, t2, t3, t4, t5, t6, t7, t8;
  sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

  sp_simd_t x0e, x1e, t0e, t1e, p0e, p1e, p2e;

  j0 = start;
  j1 = sp_array_inc(j0, inc, n);
  j2 = sp_array_inc(j0, 2 * inc, n);
  j3 = sp_array_inc(j0, 3 * inc, n);
  j4 = sp_array_inc(j0, 4 * inc, n);
  j5 = sp_array_inc(j0, 5 * inc, n);
  j6 = sp_array_inc(j0, 6 * inc, n);
  j7 = sp_array_inc(j0, 7 * inc, n);
  j8 = sp_array_inc(j0, 8 * inc, n);

  x0 = sp_simd_gather(x, j0, inc2, n);
  x1 = sp_simd_gather(x, j1, inc2, n);
  x2 = sp_simd_gather(x, j2, inc2, n);
  x0e = sp_simd_gather(x, j3, inc2, n);
  x3 = sp_simd_gather(x, j4, inc2, n);
  x6 = sp_simd_gather(x, j5, inc2, n);
  x1e = sp_simd_gather(x, j6, inc2, n);
  x5 = sp_simd_gather(x, j7, inc2, n);
  x4 = sp_simd_gather(x, j8, inc2, n);

  t0e = sp_simd_add(x0e, x1e, p);
  t1e = sp_simd_sub(x0e, x1e, p);

  p1 = sp_simd_add(x1, x3, p);
  p1 = sp_simd_add(p1, x5, p);
  p2 = sp_simd_add(x2, x4, p);
  p2 = sp_simd_add(p2, x6, p);
  p3 = sp_simd_sub(x1, x5, p);
  p4 = sp_simd_sub(x2, x6, p);
  p5 = sp_simd_sub(x3, x5, p);
  p6 = sp_simd_sub(x4, x6, p);

  t1 = sp_simd_add(p1, p2, p);
  t2 = sp_simd_sub(p1, p2, p);
  t3 = sp_simd_sub(p3, p5, p);
  t5 = sp_simd_add(t3, p6, p);
  t3 = sp_simd_sub(t3, p6, p);
  t4 = sp_simd_add(p4, p5, p);
  t6 = sp_simd_sub(p4, p5, p);

  p0e = sp_simd_add(x0, t0e, p);
  p0 = t1;
  p1 = t1;
  p2 = t2;
  p3 = t3;
  p4 = t4;
  p5 = sp_simd_add(t3, t4, p);
  p6 = t5;
  p7 = t6;
  p8 = sp_simd_add(t5, t6, p);

  p1 = sp_simd_mul(p1, ntt_const[1], p, d);
  p2 = sp_simd_mul(p2, ntt_const[2], p, d);
  t0e = sp_simd_mul(t0e, ntt_const[3], p, d);
  t1e = sp_simd_mul(t1e, ntt_const[4], p, d);
  p3 = sp_simd_mul(p3, ntt_const[5], p, d);
  p4 = sp_simd_mul(p4, ntt_const[6], p, d);
  p5 = sp_simd_mul(p5, ntt_const[7], p, d);
  p6 = sp_simd_mul(p6, ntt_const[8], p, d);
  p7 = sp_simd_mul(p7, ntt_const[9], p, d);
  p8 = sp_simd_mul(p8, ntt_const[10], p, d);

  t0e = sp_simd_add(t0e, p0e, p);
  t1 = sp_simd_add(p1, p2, p);
  t2 = sp_simd_sub(p1, p2, p);
  t3 = sp_simd_add(p3, p5, p);
  t4 = sp_simd_add(p4, p5, p);
  t5 = sp_simd_add(p6, p8, p);
  t6 = sp_simd_add(p7, p8, p);

  p1e = sp_simd_add(t0e, t1e, p);
  p2e = sp_simd_sub(t0e, t1e, p);
  p3 = sp_simd_add(t3, t5, p);
  p4 = sp_simd_add(t4, t6, p);
  p5 = sp_simd_sub(t4, t6, p);
  p5 = sp_simd_sub(p5, p3, p);
  p6 = sp_simd_sub(t5, t3, p);

  p0 = sp_simd_add(p0, p0e, p);
  t1 = sp_simd_add(t1, p0e, p);
  t2 = sp_simd_add(t2, p0e, p);
  t3 = sp_simd_add(p3, p1e, p);
  t4 = sp_simd_add(p4, p2e, p);
  t5 = sp_simd_add(p5, p1e, p);
  t6 = sp_simd_add(p6, p2e, p);
  t7 = sp_simd_add(p3, p5, p);
  t7 = sp_simd_sub(p1e, t7, p);
  t8 = sp_simd_add(p4, p6, p);
  t8 = sp_simd_sub(p2e, t8, p);

  sp_simd_scatter(p0, x, j0, inc2, n);
  sp_simd_scatter(t8, x, j1, inc2, n);
  sp_simd_scatter(t3, x, j2, inc2, n);
  sp_simd_scatter(t2, x, j3, inc2, n);
  sp_simd_scatter(t4, x, j4, inc2, n);
  sp_simd_scatter(t7, x, j5, inc2, n);
  sp_simd_scatter(t1, x, j6, inc2, n);
  sp_simd_scatter(t6, x, j7, inc2, n);
  sp_simd_scatter(t5, x, j8, inc2, n);
}
#endif

static void
ntt9_pfa_run(spv_t x, spv_size_t stride,
	  spv_size_t cofactor,
	  sp_t p, sp_t d, spv_t ntt_const)
{
  spv_size_t i = 0;
  spv_size_t incstart = 0;
  spv_size_t n = 9 * cofactor * stride;
  spv_size_t inc = cofactor * stride;
  spv_size_t inc2 = 9 * stride;

#ifdef HAVE_SSE2
  spv_size_t num_simd = SP_SIMD_VSIZE * (cofactor / SP_SIMD_VSIZE);

  for (i = 0; i < num_simd; i += SP_SIMD_VSIZE)
    {
      ntt9_pfa_run_core_simd(x, incstart, inc, inc2, n, p, d, ntt_const);
      incstart += SP_SIMD_VSIZE * inc2;
    }
#endif

  for (; i < cofactor; i++, incstart += inc2)
    ntt9_pfa_run_core(x, incstart, inc, n, p, d, ntt_const);

}

const nttconfig_t ntt9_config = 
{
  9,
  ntt9_get_num_const,
  ntt9_init,
  ntt9_run,
  ntt9_pfa_run
};

