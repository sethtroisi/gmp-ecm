#include "ntt-impl.h"

#define NC 6

static const uint8_t fixed_const[NC] = {1};

static const uint8_t *
ntt5_get_fixed_ntt_const(void)
{
  return fixed_const;
}

void
ntt5_init(spv_t out, sp_t p, sp_t d,
	  sp_t primroot, sp_t order)
{
  nttdata_init_generic(&ntt5_config, out, p, d, primroot, order);
}

static void
ntt5_twiddle_run_core(spv_t x, spv_t w, spv_size_t stride,
			sp_t p, spv_t ntt_const)
{
  sp_t p0, p1, p2, p3, p4, p5;
  sp_t x0, x1, x2, x3, x4;
  sp_t     t1, t2, t3, t4;

  x0 = x[0 * stride];
  x1 = x[1 * stride];
  x4 = x[2 * stride];
  x2 = x[3 * stride];
  x3 = x[4 * stride];

  t1 = sp_ntt_add(x1, x3, p);
  t3 = sp_ntt_sub(x1, x3, p);
  t2 = sp_ntt_add(x2, x4, p);
  t4 = sp_ntt_sub(x2, x4, p);

  p1 = sp_ntt_add(t1, t2, p);
  p2 = sp_ntt_sub_partial(t1, t2, p);
  p3 = t3;
  p4 = t4;
  p5 = sp_ntt_add_partial(t3, t4, p);

  p0 = sp_ntt_add(x0, p1, p);

  p1 = sp_ntt_mul(p1, ntt_const[1], ntt_const[NC+1], p);
  p2 = sp_ntt_mul(p2, ntt_const[2], ntt_const[NC+2], p);
  p3 = sp_ntt_mul(p3, ntt_const[3], ntt_const[NC+3], p);
  p4 = sp_ntt_mul(p4, ntt_const[4], ntt_const[NC+4], p);
  p5 = sp_ntt_mul(p5, ntt_const[5], ntt_const[NC+5], p);

  p1 = sp_ntt_add(p0, p1, p);

  t1 = sp_ntt_add(p1, p2, p);
  t2 = sp_ntt_sub(p1, p2, p);
  t3 = sp_ntt_add(p3, p5, p);
  t4 = sp_ntt_add(p4, p5, p);

  p1 = sp_ntt_add_partial(t1, t3, p);
  p2 = sp_ntt_add_partial(t2, t4, p);
  p3 = sp_ntt_sub_partial(t1, t3, p);
  p4 = sp_ntt_sub_partial(t2, t4, p);

  p4 = sp_ntt_mul(p4, w[0], w[1], p);
  p3 = sp_ntt_mul(p3, w[2], w[3], p);
  p1 = sp_ntt_mul(p1, w[4], w[5], p);
  p2 = sp_ntt_mul(p2, w[6], w[7], p);

  x[0 * stride] = p0;
  x[1 * stride] = p4;
  x[2 * stride] = p3;
  x[3 * stride] = p1;
  x[4 * stride] = p2;
}

#ifdef HAVE_SSE2
static void
ntt5_twiddle_run_core_simd(spv_t x, sp_simd_t *w,
			spv_size_t stride,
			sp_t p, spv_t ntt_const)
{
  sp_simd_t p0, p1, p2, p3, p4, p5;
  sp_simd_t x0, x1, x2, x3, x4;
  sp_simd_t     t1, t2, t3, t4;

  x0 = sp_simd_gather(x + 0 * stride);
  x1 = sp_simd_gather(x + 1 * stride);
  x4 = sp_simd_gather(x + 2 * stride);
  x2 = sp_simd_gather(x + 3 * stride);
  x3 = sp_simd_gather(x + 4 * stride);

  t1 = sp_ntt_add_simd(x1, x3, p);
  t3 = sp_ntt_sub_simd(x1, x3, p);
  t2 = sp_ntt_add_simd(x2, x4, p);
  t4 = sp_ntt_sub_simd(x2, x4, p);

  p1 = sp_ntt_add_simd(t1, t2, p);
  p2 = sp_ntt_sub_partial_simd(t1, t2, p);
  p3 = t3;
  p4 = t4;
  p5 = sp_ntt_add_partial_simd(t3, t4, p);

  p0 = sp_ntt_add_simd(x0, p1, p);

  p1 = sp_ntt_mul_simd(p1, ntt_const[1], ntt_const[NC+1], p);
  p2 = sp_ntt_mul_simd(p2, ntt_const[2], ntt_const[NC+2], p);
  p3 = sp_ntt_mul_simd(p3, ntt_const[3], ntt_const[NC+3], p);
  p4 = sp_ntt_mul_simd(p4, ntt_const[4], ntt_const[NC+4], p);
  p5 = sp_ntt_mul_simd(p5, ntt_const[5], ntt_const[NC+5], p);

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
  p3 = sp_ntt_twiddle_mul_simd(p3, w + 2, p);
  p1 = sp_ntt_twiddle_mul_simd(p1, w + 4, p);
  p2 = sp_ntt_twiddle_mul_simd(p2, w + 6, p);

  sp_simd_scatter(p0, x + 0 * stride);
  sp_simd_scatter(p4, x + 1 * stride);
  sp_simd_scatter(p3, x + 2 * stride);
  sp_simd_scatter(p1, x + 3 * stride);
  sp_simd_scatter(p2, x + 4 * stride);
}
#endif

static void
ntt5_twiddle_run(spv_t x, spv_t w,
	  spv_size_t stride,
	  spv_size_t num_transforms,
	  sp_t p, spv_t ntt_const)
{
  spv_size_t i = 0, j = 0;

#ifdef HAVE_SSE2
  spv_size_t num_simd = SP_SIMD_VSIZE * (num_transforms / SP_SIMD_VSIZE);

  for (i = 0; i < num_simd; i += SP_SIMD_VSIZE,
		  	j += 2*(5-1)*SP_SIMD_VSIZE)
    ntt5_twiddle_run_core_simd(x + i, (sp_simd_t *)(w + j),
				stride, p, ntt_const);
#endif

  for (; i < num_transforms; i++, j += 2*(5-1))
    ntt5_twiddle_run_core(x + i, w + j, stride, p, ntt_const);
}

static void
ntt5_pfa_run_core(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t n,
	  sp_t p, spv_t ntt_const)
{
  spv_size_t j0, j1, j2, j3, j4;
  sp_t p0, p1, p2, p3, p4, p5;
  sp_t x0, x1, x2, x3, x4;
  sp_t     t1, t2, t3, t4;

  j0 = start;
  j1 = sp_array_inc(j0, inc, n);
  j2 = sp_array_inc(j0, 2 * inc, n);
  j3 = sp_array_inc(j0, 3 * inc, n);
  j4 = sp_array_inc(j0, 4 * inc, n);

  x0 = x[j0];
  x1 = x[j1];
  x4 = x[j2];
  x2 = x[j3];
  x3 = x[j4];

  t1 = sp_ntt_add(x1, x3, p);
  t3 = sp_ntt_sub(x1, x3, p);
  t2 = sp_ntt_add(x2, x4, p);
  t4 = sp_ntt_sub(x2, x4, p);

  p1 = sp_ntt_add(t1, t2, p);
  p2 = sp_ntt_sub_partial(t1, t2, p);
  p3 = t3;
  p4 = t4;
  p5 = sp_ntt_add_partial(t3, t4, p);

  p0 = sp_ntt_add(x0, p1, p);

  p1 = sp_ntt_mul(p1, ntt_const[1], ntt_const[NC+1], p);
  p2 = sp_ntt_mul(p2, ntt_const[2], ntt_const[NC+2], p);
  p3 = sp_ntt_mul(p3, ntt_const[3], ntt_const[NC+3], p);
  p4 = sp_ntt_mul(p4, ntt_const[4], ntt_const[NC+4], p);
  p5 = sp_ntt_mul(p5, ntt_const[5], ntt_const[NC+5], p);

  p1 = sp_ntt_add(p0, p1, p);

  t1 = sp_ntt_add(p1, p2, p);
  t2 = sp_ntt_sub(p1, p2, p);
  t3 = sp_ntt_add(p3, p5, p);
  t4 = sp_ntt_add(p4, p5, p);

  p1 = sp_ntt_add(t1, t3, p);
  p2 = sp_ntt_add(t2, t4, p);
  p3 = sp_ntt_sub(t1, t3, p);
  p4 = sp_ntt_sub(t2, t4, p);

  x[j0] = p0;
  x[j1] = p4;
  x[j2] = p3;
  x[j3] = p1;
  x[j4] = p2;
}

#ifdef HAVE_SSE2
static void
ntt5_pfa_run_core_simd(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t inc2, spv_size_t n,
	  sp_t p, spv_t ntt_const)
{
  spv_size_t j0, j1, j2, j3, j4;
  sp_simd_t p0, p1, p2, p3, p4, p5;
  sp_simd_t x0, x1, x2, x3, x4;
  sp_simd_t     t1, t2, t3, t4;

  j0 = start;
  j1 = sp_array_inc(j0, inc, n);
  j2 = sp_array_inc(j0, 2 * inc, n);
  j3 = sp_array_inc(j0, 3 * inc, n);
  j4 = sp_array_inc(j0, 4 * inc, n);

  x0 = sp_simd_pfa_gather(x, j0, inc2, n);
  x1 = sp_simd_pfa_gather(x, j1, inc2, n);
  x4 = sp_simd_pfa_gather(x, j2, inc2, n);
  x2 = sp_simd_pfa_gather(x, j3, inc2, n);
  x3 = sp_simd_pfa_gather(x, j4, inc2, n);

  t1 = sp_ntt_add_simd(x1, x3, p);
  t3 = sp_ntt_sub_simd(x1, x3, p);
  t2 = sp_ntt_add_simd(x2, x4, p);
  t4 = sp_ntt_sub_simd(x2, x4, p);

  p1 = sp_ntt_add_simd(t1, t2, p);
  p2 = sp_ntt_sub_partial_simd(t1, t2, p);
  p3 = t3;
  p4 = t4;
  p5 = sp_ntt_add_partial_simd(t3, t4, p);

  p0 = sp_ntt_add_simd(x0, p1, p);

  p1 = sp_ntt_mul_simd(p1, ntt_const[1], ntt_const[NC+1], p);
  p2 = sp_ntt_mul_simd(p2, ntt_const[2], ntt_const[NC+2], p);
  p3 = sp_ntt_mul_simd(p3, ntt_const[3], ntt_const[NC+3], p);
  p4 = sp_ntt_mul_simd(p4, ntt_const[4], ntt_const[NC+4], p);
  p5 = sp_ntt_mul_simd(p5, ntt_const[5], ntt_const[NC+5], p);

  p1 = sp_ntt_add_simd(p0, p1, p);

  t1 = sp_ntt_add_simd(p1, p2, p);
  t2 = sp_ntt_sub_simd(p1, p2, p);
  t3 = sp_ntt_add_simd(p3, p5, p);
  t4 = sp_ntt_add_simd(p4, p5, p);

  p1 = sp_ntt_add_simd(t1, t3, p);
  p2 = sp_ntt_add_simd(t2, t4, p);
  p3 = sp_ntt_sub_simd(t1, t3, p);
  p4 = sp_ntt_sub_simd(t2, t4, p);

  sp_simd_pfa_scatter(p0, x, j0, inc2, n);
  sp_simd_pfa_scatter(p4, x, j1, inc2, n);
  sp_simd_pfa_scatter(p3, x, j2, inc2, n);
  sp_simd_pfa_scatter(p1, x, j3, inc2, n);
  sp_simd_pfa_scatter(p2, x, j4, inc2, n);
}
#endif

static void
ntt5_pfa_run(spv_t x, spv_size_t cofactor,
	  sp_t p, spv_t ntt_const)
{
  spv_size_t i = 0;
  spv_size_t incstart = 0;
  spv_size_t n = 5 * cofactor;
  spv_size_t inc = cofactor;
  spv_size_t inc2 = 5;

#ifdef HAVE_SSE2
  spv_size_t num_simd = SP_SIMD_VSIZE * (cofactor / SP_SIMD_VSIZE);

  for (i = 0; i < num_simd; i += SP_SIMD_VSIZE)
    {
      ntt5_pfa_run_core_simd(x, incstart, inc, inc2, n, p, ntt_const);
      incstart += SP_SIMD_VSIZE * inc2;
    }
#endif

  for (; i < cofactor; i++, incstart += inc2)
    ntt5_pfa_run_core(x, incstart, inc, n, p, ntt_const);
}

const nttconfig_t ntt5_config = 
{
  5,
  NC,
  ntt5_get_fixed_ntt_const,
  ntt5_init,
  ntt5_pfa_run,
  ntt5_twiddle_run
};

