#include "ntt-impl.h"

#define NC 8

static const uint8_t fixed_const[NC] = {1, 1, 1, 0, 1, 0, 0, 0};

static const uint8_t *
ntt8_get_fixed_ntt_const(void)
{
  return fixed_const;
}

void
ntt8_init(spv_t out, sp_t p, sp_t d, 
	  sp_t primroot, sp_t order)
{
  sp_t w1 = sp_pow(primroot, order / 8, p, d);
  sp_t w2 = sp_sqr(w1, p, d);
  sp_t w3 = sp_mul(w2, w1, p, d);
  sp_t inv2 = sp_inv(2, p, d);

  out[0] = 1;
  out[1] = 1;
  out[2] = 1;
  out[3] = w2;
  out[4] = 1;
  out[5] = w2;
  out[6] = sp_mul(inv2, sp_sub(w1, w3, p), p, d);
  out[7] = sp_mul(inv2, sp_add(w1, w3, p), p, d);
}

static void
ntt8_run(spv_t x, spv_size_t stride,
	  sp_t p, spv_t ntt_const)
{
  sp_t x0, x1, x2, x3, x4, x5, x6, x7;
  sp_t t0, t1, t2, t3, t4, t5, t6, t7; 
  sp_t p0, p1, p2, p3, p4, p5, p6, p7;

  x0 = x[0 * stride];
  x1 = x[1 * stride];
  x2 = x[2 * stride];
  x3 = x[3 * stride];
  x4 = x[4 * stride];
  x5 = x[5 * stride];
  x6 = x[6 * stride];
  x7 = x[7 * stride];

  t0 = sp_ntt_add(x0, x4, p);
  t4 = sp_ntt_sub(x0, x4, p);
  t1 = sp_ntt_add(x1, x5, p);
  t5 = sp_ntt_sub(x1, x5, p);
  t2 = sp_ntt_add(x2, x6, p);
  t6 = sp_ntt_sub_partial(x2, x6, p);
  t3 = sp_ntt_add(x3, x7, p);
  t7 = sp_ntt_sub(x3, x7, p);

  p0 = sp_ntt_add(t0, t2, p);
  p1 = sp_ntt_sub(t0, t2, p);
  p2 = sp_ntt_add(t1, t3, p);
  p3 = sp_ntt_sub_partial(t1, t3, p);
  p4 = t4;
  p5 = t6;
  p6 = sp_ntt_sub_partial(t5, t7, p);
  p7 = sp_ntt_add_partial(t5, t7, p); 

  p3 = sp_ntt_mul(p3, ntt_const[3], ntt_const[NC+3], p);
  p5 = sp_ntt_mul(p5, ntt_const[5], ntt_const[NC+5], p);
  p6 = sp_ntt_mul(p6, ntt_const[6], ntt_const[NC+6], p);
  p7 = sp_ntt_mul(p7, ntt_const[7], ntt_const[NC+7], p);

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

  x[0 * stride] = t0;
  x[1 * stride] = t4;
  x[2 * stride] = t2;
  x[3 * stride] = t7;
  x[4 * stride] = t1;
  x[5 * stride] = t5;
  x[6 * stride] = t3;
  x[7 * stride] = t6;
}

#ifdef HAVE_SSE2
static void
ntt8_run_simd(spv_t x, spv_size_t stride,
	  sp_t p, spv_t ntt_const)
{
  sp_simd_t x0, x1, x2, x3, x4, x5, x6, x7;
  sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7; 
  sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7;

  x0 = sp_simd_gather(x + 0 * stride);
  x1 = sp_simd_gather(x + 1 * stride);
  x2 = sp_simd_gather(x + 2 * stride);
  x3 = sp_simd_gather(x + 3 * stride);
  x4 = sp_simd_gather(x + 4 * stride);
  x5 = sp_simd_gather(x + 5 * stride);
  x6 = sp_simd_gather(x + 6 * stride);
  x7 = sp_simd_gather(x + 7 * stride);

  t0 = sp_ntt_add_simd(x0, x4, p);
  t4 = sp_ntt_sub_simd(x0, x4, p);
  t1 = sp_ntt_add_simd(x1, x5, p);
  t5 = sp_ntt_sub_simd(x1, x5, p);
  t2 = sp_ntt_add_simd(x2, x6, p);
  t6 = sp_ntt_sub_partial_simd(x2, x6, p);
  t3 = sp_ntt_add_simd(x3, x7, p);
  t7 = sp_ntt_sub_simd(x3, x7, p);

  p0 = sp_ntt_add_simd(t0, t2, p);
  p1 = sp_ntt_sub_simd(t0, t2, p);
  p2 = sp_ntt_add_simd(t1, t3, p);
  p3 = sp_ntt_sub_partial_simd(t1, t3, p);
  p4 = t4;
  p5 = t6;
  p6 = sp_ntt_sub_partial_simd(t5, t7, p);
  p7 = sp_ntt_add_partial_simd(t5, t7, p); 

  p3 = sp_ntt_mul_simd(p3, ntt_const[3], ntt_const[NC+3], p);
  p5 = sp_ntt_mul_simd(p5, ntt_const[5], ntt_const[NC+5], p);
  p6 = sp_ntt_mul_simd(p6, ntt_const[6], ntt_const[NC+6], p);
  p7 = sp_ntt_mul_simd(p7, ntt_const[7], ntt_const[NC+7], p);

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

  sp_simd_scatter(t0, x + 0 * stride);
  sp_simd_scatter(t4, x + 1 * stride);
  sp_simd_scatter(t2, x + 2 * stride);
  sp_simd_scatter(t7, x + 3 * stride);
  sp_simd_scatter(t1, x + 4 * stride);
  sp_simd_scatter(t5, x + 5 * stride);
  sp_simd_scatter(t3, x + 6 * stride);
  sp_simd_scatter(t6, x + 7 * stride);
}
#endif

static void
ntt8_twiddle_run(spv_t x, spv_size_t stride,
	  spv_size_t num_transforms,
	  sp_t p, spv_t ntt_const)
{
  spv_size_t i = 0;

#ifdef HAVE_SSE2
  spv_size_t num_simd = SP_SIMD_VSIZE * (num_transforms / SP_SIMD_VSIZE);

  for (i = 0; i < num_simd; i += SP_SIMD_VSIZE)
      ntt8_run_simd(x + i, stride, p, ntt_const);
#endif

  for (; i < num_transforms; i++)
    ntt8_run(x + i, stride, p, ntt_const);
}

static void
ntt8_pfa_run_core(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t n,
	  sp_t p, spv_t ntt_const)
{
  spv_size_t j0, j1, j2, j3, j4, j5, j6, j7;
  sp_t x0, x1, x2, x3, x4, x5, x6, x7;
  sp_t t0, t1, t2, t3, t4, t5, t6, t7; 
  sp_t p0, p1, p2, p3, p4, p5, p6, p7;

  j0 = start;
  j1 = sp_array_inc(j0, inc, n);
  j2 = sp_array_inc(j0, 2 * inc, n);
  j3 = sp_array_inc(j0, 3 * inc, n);
  j4 = sp_array_inc(j0, 4 * inc, n);
  j5 = sp_array_inc(j0, 5 * inc, n);
  j6 = sp_array_inc(j0, 6 * inc, n);
  j7 = sp_array_inc(j0, 7 * inc, n);

  x0 = x[j0];
  x1 = x[j1];
  x2 = x[j2];
  x3 = x[j3];
  x4 = x[j4];
  x5 = x[j5];
  x6 = x[j6];
  x7 = x[j7];

  t0 = sp_ntt_add(x0, x4, p);
  t4 = sp_ntt_sub(x0, x4, p);
  t1 = sp_ntt_add(x1, x5, p);
  t5 = sp_ntt_sub(x1, x5, p);
  t2 = sp_ntt_add(x2, x6, p);
  t6 = sp_ntt_sub_partial(x2, x6, p);
  t3 = sp_ntt_add(x3, x7, p);
  t7 = sp_ntt_sub(x3, x7, p);

  p0 = sp_ntt_add(t0, t2, p);
  p1 = sp_ntt_sub(t0, t2, p);
  p2 = sp_ntt_add(t1, t3, p);
  p3 = sp_ntt_sub_partial(t1, t3, p);
  p4 = t4;
  p5 = t6;
  p6 = sp_ntt_sub_partial(t5, t7, p);
  p7 = sp_ntt_add_partial(t5, t7, p); 

  p3 = sp_ntt_mul(p3, ntt_const[3], ntt_const[NC+3], p);
  p5 = sp_ntt_mul(p5, ntt_const[5], ntt_const[NC+5], p);
  p6 = sp_ntt_mul(p6, ntt_const[6], ntt_const[NC+6], p);
  p7 = sp_ntt_mul(p7, ntt_const[7], ntt_const[NC+7], p);

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

  x[j0] = t0;
  x[j1] = t4;
  x[j2] = t2;
  x[j3] = t7;
  x[j4] = t1;
  x[j5] = t5;
  x[j6] = t3;
  x[j7] = t6;
}

#ifdef HAVE_SSE2
static void
ntt8_pfa_run_core_simd(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t inc2, spv_size_t n,
	  sp_t p, spv_t ntt_const)
{
  spv_size_t j0, j1, j2, j3, j4, j5, j6, j7;
  sp_simd_t x0, x1, x2, x3, x4, x5, x6, x7;
  sp_simd_t t0, t1, t2, t3, t4, t5, t6, t7; 
  sp_simd_t p0, p1, p2, p3, p4, p5, p6, p7;

  j0 = start;
  j1 = sp_array_inc(j0, inc, n);
  j2 = sp_array_inc(j0, 2 * inc, n);
  j3 = sp_array_inc(j0, 3 * inc, n);
  j4 = sp_array_inc(j0, 4 * inc, n);
  j5 = sp_array_inc(j0, 5 * inc, n);
  j6 = sp_array_inc(j0, 6 * inc, n);
  j7 = sp_array_inc(j0, 7 * inc, n);

  x0 = sp_simd_pfa_gather(x, j0, inc2, n);
  x1 = sp_simd_pfa_gather(x, j1, inc2, n);
  x2 = sp_simd_pfa_gather(x, j2, inc2, n);
  x3 = sp_simd_pfa_gather(x, j3, inc2, n);
  x4 = sp_simd_pfa_gather(x, j4, inc2, n);
  x5 = sp_simd_pfa_gather(x, j5, inc2, n);
  x6 = sp_simd_pfa_gather(x, j6, inc2, n);
  x7 = sp_simd_pfa_gather(x, j7, inc2, n);

  t0 = sp_ntt_add_simd(x0, x4, p);
  t4 = sp_ntt_sub_simd(x0, x4, p);
  t1 = sp_ntt_add_simd(x1, x5, p);
  t5 = sp_ntt_sub_simd(x1, x5, p);
  t2 = sp_ntt_add_simd(x2, x6, p);
  t6 = sp_ntt_sub_partial_simd(x2, x6, p);
  t3 = sp_ntt_add_simd(x3, x7, p);
  t7 = sp_ntt_sub_simd(x3, x7, p);

  p0 = sp_ntt_add_simd(t0, t2, p);
  p1 = sp_ntt_sub_simd(t0, t2, p);
  p2 = sp_ntt_add_simd(t1, t3, p);
  p3 = sp_ntt_sub_partial_simd(t1, t3, p);
  p4 = t4;
  p5 = t6;
  p6 = sp_ntt_sub_partial_simd(t5, t7, p);
  p7 = sp_ntt_add_partial_simd(t5, t7, p); 

  p3 = sp_ntt_mul_simd(p3, ntt_const[3], ntt_const[NC+3], p);
  p5 = sp_ntt_mul_simd(p5, ntt_const[5], ntt_const[NC+5], p);
  p6 = sp_ntt_mul_simd(p6, ntt_const[6], ntt_const[NC+6], p);
  p7 = sp_ntt_mul_simd(p7, ntt_const[7], ntt_const[NC+7], p);

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

  sp_simd_pfa_scatter(t0, x, j0, inc2, n);
  sp_simd_pfa_scatter(t4, x, j1, inc2, n);
  sp_simd_pfa_scatter(t2, x, j2, inc2, n);
  sp_simd_pfa_scatter(t7, x, j3, inc2, n);
  sp_simd_pfa_scatter(t1, x, j4, inc2, n);
  sp_simd_pfa_scatter(t5, x, j5, inc2, n);
  sp_simd_pfa_scatter(t3, x, j6, inc2, n);
  sp_simd_pfa_scatter(t6, x, j7, inc2, n);
}
#endif

static void
ntt8_pfa_run(spv_t x, spv_size_t cofactor,
	  sp_t p, spv_t ntt_const)
{
  spv_size_t i = 0;
  spv_size_t incstart = 0;
  spv_size_t n = 8 * cofactor;
  spv_size_t inc = cofactor;
  spv_size_t inc2 = 8;

#ifdef HAVE_SSE2
  spv_size_t num_simd = SP_SIMD_VSIZE * (cofactor / SP_SIMD_VSIZE);

  for (i = 0; i < num_simd; i += SP_SIMD_VSIZE)
    {
      ntt8_pfa_run_core_simd(x, incstart, inc, inc2, n, p, ntt_const);
      incstart += SP_SIMD_VSIZE * inc2;
    }
#endif

  for (; i < cofactor; i++, incstart += inc2)
    ntt8_pfa_run_core(x, incstart, inc, n, p, ntt_const);

}

const nttconfig_t ntt8_config = 
{
  8,
  NC,
  ntt8_get_fixed_ntt_const,
  ntt8_init,
  ntt8_run,
  ntt8_pfa_run,
  ntt8_twiddle_run
};

