#include <stdio.h>
#include "sp.h"

typedef uint32_t (*get_num_ntt_const_t)(void);

typedef void (*nttdata_init_t)(spv_t out, 
				sp_t p, sp_t d,
				sp_t primroot, sp_t order);

typedef void (*ntt_run_t)(spv_t out, spv_size_t ostride,
    			  spv_t in, spv_size_t istride,
			  sp_t p, sp_t d, spv_t ntt_const);

typedef struct
{
  uint32_t size;
  get_num_ntt_const_t get_num_ntt_const;
  nttdata_init_t nttdata_init;
  ntt_run_t ntt_run;
} nttconfig_t;

/*---------------------- NTT2 --------------------------------------*/

static uint32_t 
ntt2_get_num_const(void)
{
  return 2;
}

static void
ntt2_init(spv_t out, sp_t p, sp_t d, 
	  sp_t primroot, sp_t order)
{
  out[0] = 1;
  out[1] = 1;
}

static void
ntt2_run(spv_t out, spv_size_t ostride,
    	  spv_t in, spv_size_t istride,
	  sp_t p, sp_t d, spv_t ntt_const)
{
  sp_t x0, x1;

  x0 = in[0 * istride];
  x1 = in[1 * istride];

  out[0 * ostride] = sp_add(x0, x1, p);
  out[1 * ostride] = sp_sub(x0, x1, p);
}

const nttconfig_t ntt2_config = 
{
  2,
  ntt2_get_num_const,
  ntt2_init,
  ntt2_run
};

/*---------------------- NTT3 --------------------------------------*/

static uint32_t 
ntt3_get_num_const(void)
{
  return 3;
}

static void
ntt3_init(spv_t out, sp_t p, sp_t d, 
	  sp_t primroot, sp_t order)
{
  sp_t w1, w2;
  sp_t h1, h2;
  sp_t inv2 = sp_inv(2, p, d);

  w1 = sp_pow(primroot, order / 3, p, d);
  w2 = sp_sqr(w1, p, d);

  h1 = sp_add(w1, w2, p);
  h2 = sp_sub(w1, w2, p);

  out[0] = 1;
  out[1] = sp_sub(sp_mul(h1, inv2, p, d), 1, p);
  out[2] = sp_mul(h2, inv2, p, d);
}

static void
ntt3_run(spv_t out, spv_size_t ostride,
    	  spv_t in, spv_size_t istride,
	  sp_t p, sp_t d, spv_t ntt_const)
{
  sp_t p0, p1, p2;
  sp_t x0, x1, x2;
  sp_t     t1, t2;

  x0 = in[0 * istride];
  x1 = in[1 * istride];
  x2 = in[2 * istride];

  t1 = sp_add(x1, x2, p);
  t2 = sp_sub(x1, x2, p);

  p0 = sp_add(x0, t1, p);

  p1 = sp_mul(t1, ntt_const[2], p, d);
  p2 = sp_mul(t2, ntt_const[3], p, d);

  p1 = sp_add(p0, p1, p);

  t1 = sp_add(p1, p2, p);
  t2 = sp_sub(p1, p2, p);

  out[0 * ostride] = p0;
  out[1 * ostride] = t1;
  out[2 * ostride] = t2;
}

const nttconfig_t ntt3_config = 
{
  3,
  ntt3_get_num_const,
  ntt3_init,
  ntt3_run
};

/*---------------------- NTT4 --------------------------------------*/

static uint32_t 
ntt4_get_num_const(void)
{
  return 4;
}

static void
ntt4_init(spv_t out, sp_t p, sp_t d, 
	  sp_t primroot, sp_t order)
{
  out[0] = 1;
  out[1] = 1;
  out[2] = 1;
  out[3] = sp_pow(primroot, order / 4, p, d);
}

static void
ntt4_run(spv_t out, spv_size_t ostride,
    	  spv_t in, spv_size_t istride,
	  sp_t p, sp_t d, spv_t ntt_const)
{
  sp_t x0, x1, x2, x3;
  sp_t t0, t1, t2, t3;
  sp_t p0, p1, p2, p3;

  x0 = in[0 * istride];
  x1 = in[1 * istride];
  x2 = in[2 * istride];
  x3 = in[3 * istride];

  t0 = sp_add(x0, x2, p);
  t2 = sp_sub(x0, x2, p);
  t1 = sp_add(x1, x3, p);
  t3 = sp_sub(x1, x3, p);

  t3 = sp_mul(t3, ntt_const[3], p, d);

  p0 = sp_add(t0, t1, p);
  p1 = sp_sub(t0, t1, p);
  p2 = sp_add(t2, t3, p);
  p3 = sp_sub(t2, t3, p);

  out[0 * ostride] = p0;
  out[1 * ostride] = p2;
  out[2 * ostride] = p1;
  out[3 * ostride] = p3;
}

const nttconfig_t ntt4_config = 
{
  4,
  ntt4_get_num_const,
  ntt4_init,
  ntt4_run
};

/*---------------------- NTT5 --------------------------------------*/

static uint32_t 
ntt5_get_num_const(void)
{
  return 6;
}

static void
ntt5_init(spv_t out, sp_t p, sp_t d, 
	  sp_t primroot, sp_t order)
{
  sp_t w1, w2, w3, w4;
  sp_t h1, h2, h3, h4, h5;
  sp_t t1, t2, t3, t4;
  sp_t inv4 = sp_inv(4, p, d);
    
  h1 = sp_pow(primroot, order / 5, p, d);
  h2 = sp_mul(h1, h1, p, d);
  h3 = sp_mul(h2, h1, p, d);
  h4 = sp_mul(h3, h1, p, d);

  w1 = h3;
  w2 = h4;
  w3 = h2;
  w4 = h1;

  t1 = sp_add(w1, w2, p);
  t1 = sp_add(t1, w3, p);
  t1 = sp_add(t1, w4, p);

  t2 = sp_sub(w1, w2, p);
  t2 = sp_add(t2, w3, p);
  t2 = sp_sub(t2, w4, p);

  t3 = sp_add(w1, w1, p);
  t3 = sp_sub(t3, w3, p);
  t3 = sp_sub(t3, w3, p);

  t4 = sp_add(w2, w2, p);
  t4 = sp_sub(t4, w4, p);
  t4 = sp_sub(t4, w4, p);

  h1 = t1;
  h2 = t2;
  h3 = sp_sub(t3, t4, p);
  h4 = sp_neg(sp_add(t3, t4, p), p);
  h5 = t4;

  out[0] = 1;
  out[1] = sp_sub(sp_mul(h1, inv4, p, d), 1, p);
  out[2] = sp_mul(h2, inv4, p, d);
  out[3] = sp_mul(h3, inv4, p, d);
  out[4] = sp_mul(h4, inv4, p, d);
  out[5] = sp_mul(h5, inv4, p, d);
}

static void
ntt5_run(spv_t out, spv_size_t ostride,
    	  spv_t in, spv_size_t istride,
	  sp_t p, sp_t d, spv_t ntt_const)
{
  sp_t p0, p1, p2, p3, p4, p5;
  sp_t x0, x1, x2, x3, x4;
  sp_t     t1, t2, t3, t4;

  x0 = in[0 * istride];
  x1 = in[1 * istride];
  x4 = in[2 * istride];
  x2 = in[3 * istride];
  x3 = in[4 * istride];

  t1 = sp_add(x1, x3, p);
  t3 = sp_sub(x1, x3, p);
  t2 = sp_add(x2, x4, p);
  t4 = sp_sub(x2, x4, p);

  p1 = sp_add(t1, t2, p);
  p2 = sp_sub(t1, t2, p);
  p3 = t3;
  p4 = t4;
  p5 = sp_add(t3, t4, p);

  p0 = sp_add(x0, p1, p);

  p1 = sp_mul(p1, ntt_const[1], p, d);
  p2 = sp_mul(p2, ntt_const[2], p, d);
  p3 = sp_mul(p3, ntt_const[3], p, d);
  p4 = sp_mul(p4, ntt_const[4], p, d);
  p5 = sp_mul(p5, ntt_const[5], p, d);

  p1 = sp_add(p0, p1, p);

  t1 = sp_add(p1, p2, p);
  t2 = sp_sub(p1, p2, p);
  t3 = sp_add(p3, p5, p);
  t4 = sp_add(p4, p5, p);

  p1 = sp_add(t1, t3, p);
  p2 = sp_add(t2, t4, p);
  p3 = sp_sub(t1, t3, p);
  p4 = sp_sub(t2, t4, p);

  out[0 * ostride] = p0;
  out[1 * ostride] = p4;
  out[2 * ostride] = p3;
  out[3 * ostride] = p1;
  out[4 * ostride] = p2;
}

const nttconfig_t ntt5_config = 
{
  5,
  ntt5_get_num_const,
  ntt5_init,
  ntt5_run
};

/*---------------------- NTT6 --------------------------------------*/

static uint32_t 
ntt6_get_num_const(void)
{
  return ntt3_get_num_const();
}

static void
ntt6_init(spv_t out, sp_t p, sp_t d, 
	  sp_t primroot, sp_t order)
{
  ntt3_init(out, p, d, primroot, order);
}

static void
ntt6_run(spv_t out, spv_size_t ostride,
    	  spv_t in, spv_size_t istride,
	  sp_t p, sp_t d, spv_t ntt_const)
{
  sp_t t0, t1, t2, t3, t4, t5;

  {
    sp_t x0 = in[0 * istride];
    sp_t x1 = in[2 * istride];
    sp_t x2 = in[4 * istride];
    sp_t b0, b1, b2;

    b1 = sp_add(x1, x2, p);
    b2 = sp_sub(x1, x2, p);

    b0 = sp_add(x0, b1, p);

    b1 = sp_mul(b1, ntt_const[1], p, d);
    b2 = sp_mul(b2, ntt_const[2], p, d);

    b1 = sp_add(b0, b1, p);

    t0 = b0;
    t1 = sp_add(b1, b2, p);
    t2 = sp_sub(b1, b2, p);
  }
  {
    sp_t x5 = in[1 * istride];
    sp_t x3 = in[3 * istride];
    sp_t x4 = in[5 * istride];
    sp_t b0, b1, b2;

    b1 = sp_add(x4, x5, p);
    b2 = sp_sub(x4, x5, p);

    b0 = sp_add(x3, b1, p);

    b1 = sp_mul(b1, ntt_const[1], p, d);
    b2 = sp_mul(b2, ntt_const[2], p, d);

    b1 = sp_add(b0, b1, p);

    t3 = b0;
    t4 = sp_add(b1, b2, p);
    t5 = sp_sub(b1, b2, p);
  }
  {
    sp_t x0 = sp_add(t0, t3, p);
    sp_t x3 = sp_sub(t0, t3, p);

    out[0 * ostride] = x0;
    out[3 * ostride] = x3;
  }
  {
    sp_t x1 = sp_add(t1, t4, p);
    sp_t x4 = sp_sub(t1, t4, p);

    out[1 * ostride] = x4;
    out[4 * ostride] = x1;
  }
  {
    sp_t x2 = sp_add(t2, t5, p);
    sp_t x5 = sp_sub(t2, t5, p);

    out[2 * ostride] = x2;
    out[5 * ostride] = x5;
  }

}

const nttconfig_t ntt6_config = 
{
  6,
  ntt6_get_num_const,
  ntt6_init,
  ntt6_run
};

/*---------------------- NTT7 --------------------------------------*/
static uint32_t 
ntt7_get_num_const(void)
{
  return 9;
}

static void
ntt7_init(spv_t out, sp_t p, sp_t d, 
	  sp_t primroot, sp_t order)
{
  sp_t w1, w2, w3, w4, w5, w6;
  sp_t h1, h2, h3, h4, h5, h6, h7, h8;
  sp_t t1, t2, t3, t4, t5, t6, t7, t8;
  sp_t inv6 = sp_inv(6, p, d);
    
  h1 = sp_pow(primroot, order / 7, p, d);
  h2 = sp_sqr(h1, p, d);
  h3 = sp_mul(h2, h1, p, d);
  h4 = sp_mul(h3, h1, p, d);
  h5 = sp_mul(h4, h1, p, d);
  h6 = sp_mul(h5, h1, p, d);

  w1 = h5;
  w2 = h4;
  w3 = h6;
  w4 = h2;
  w5 = h3;
  w6 = h1;

  t1 = sp_add(w1, w2, p);
  t1 = sp_add(t1, w3, p);
  t1 = sp_add(t1, w4, p);
  t1 = sp_add(t1, w5, p);
  t1 = sp_add(t1, w6, p);

  t2 = sp_add(w1, w1, p);
  t2 = sp_sub(t2, w2, p);
  t2 = sp_sub(t2, w3, p);
  t2 = sp_add(t2, w4, p);
  t2 = sp_add(t2, w4, p);
  t2 = sp_sub(t2, w5, p);
  t2 = sp_sub(t2, w6, p);

  t3 = sp_neg(w1, p);
  t3 = sp_add(t3, w2, p);
  t3 = sp_add(t3, w2, p);
  t3 = sp_sub(t3, w3, p);
  t3 = sp_sub(t3, w4, p);
  t3 = sp_add(t3, w5, p);
  t3 = sp_add(t3, w5, p);
  t3 = sp_sub(t3, w6, p);

  t4 = sp_sub(w1, w2, p);
  t4 = sp_add(t4, w3, p);
  t4 = sp_sub(t4, w4, p);
  t4 = sp_add(t4, w5, p);
  t4 = sp_sub(t4, w6, p);

  t5 = sp_add(w1, w1, p);
  t5 = sp_add(t5, w2, p);
  t5 = sp_sub(t5, w3, p);
  t5 = sp_sub(t5, w4, p);
  t5 = sp_sub(t5, w4, p);
  t5 = sp_sub(t5, w5, p);
  t5 = sp_add(t5, w6, p);

  t6 = sp_add(w1, w2, p);
  t6 = sp_add(t6, w2, p);
  t6 = sp_add(t6, w3, p);
  t6 = sp_sub(t6, w4, p);
  t6 = sp_sub(t6, w5, p);
  t6 = sp_sub(t6, w5, p);
  t6 = sp_sub(t6, w6, p);

  h1 = t1;
  h2 = sp_sub(t2, t3, p);
  h3 = sp_neg(sp_add(t2, t3, p), p);
  h3 = sp_sub(h3, t3, p);
  h4 = t3;
  h5 = t4;
  h6 = sp_sub(t5, t6, p);
  h7 = sp_neg(t5, p); 
  h8 = t6;

  out[0] = 1;
  out[1] = sp_sub(sp_mul(h1, inv6, p, d), 1, p);
  out[2] = sp_mul(h2, inv6, p, d);
  out[3] = sp_mul(h3, inv6, p, d);
  out[4] = sp_mul(h4, inv6, p, d);
  out[5] = sp_mul(h5, inv6, p, d);
  out[6] = sp_mul(h6, inv6, p, d);
  out[7] = sp_mul(h7, inv6, p, d);
  out[8] = sp_mul(h8, inv6, p, d);
}

static void
ntt7_run(spv_t out, spv_size_t ostride,
    	  spv_t in, spv_size_t istride,
	  sp_t p, sp_t d, spv_t ntt_const)
{
  sp_t x0, x1, x2, x3, x4, x5, x6;
  sp_t     t1, t2, t3, t4, t5, t6, t7, t8;
  sp_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

  x0 = in[0 * istride];
  x1 = in[1 * istride];
  x5 = in[2 * istride];
  x6 = in[3 * istride];
  x3 = in[4 * istride];
  x2 = in[5 * istride];
  x4 = in[6 * istride];

  p1 = sp_add(x1, x4, p);
  p2 = sp_add(x2, x5, p);
  p3 = sp_add(x3, x6, p);
  p4 = sp_sub(x1, x4, p);
  p5 = sp_sub(x2, x5, p);
  p6 = sp_sub(x3, x6, p);

  t1 = sp_add(p1, p2, p);
  t1 = sp_add(t1, p3, p);
  t2 = sp_sub(p1, p3, p);
  t3 = sp_sub(p2, p3, p);
  t4 = sp_sub(p4, p5, p);
  t4 = sp_add(t4, p6, p);
  t5 = sp_sub(p4, p6, p);
  t6 = sp_add(p5, p6, p);

  p0 = sp_add(x0, t1, p);
  p1 = t1;
  p2 = t2;
  p3 = t3;
  p4 = sp_add(t2, t3, p);
  p5 = t4;
  p6 = t5;
  p7 = t6;
  p8 = sp_add(t5, t6, p);

  p1 = sp_mul(p1, ntt_const[1], p, d);
  p2 = sp_mul(p2, ntt_const[2], p, d);
  p3 = sp_mul(p3, ntt_const[3], p, d);
  p4 = sp_mul(p4, ntt_const[4], p, d);
  p5 = sp_mul(p5, ntt_const[5], p, d);
  p6 = sp_mul(p6, ntt_const[6], p, d);
  p7 = sp_mul(p7, ntt_const[7], p, d);
  p8 = sp_mul(p8, ntt_const[8], p, d);

  t1 = sp_add(p0, p1, p);
  t2 = sp_add(p2, p4, p);
  t3 = sp_add(p3, p4, p);
  t4 = p5;
  t5 = sp_add(p6, p8, p);
  t6 = sp_add(p7, p8, p);

  p1 = sp_add(t1, t2, p);
  p2 = sp_add(t1, t3, p);
  p3 = sp_sub(t1, t2, p);
  p3 = sp_sub(p3, t3, p);
  p4 = sp_add(t4, t5, p);
  p5 = sp_sub(t6, t4, p);
  p6 = sp_sub(t4, t5, p);
  p6 = sp_add(p6, t6, p);

  t1 = sp_add(p1, p4, p);
  t2 = sp_add(p2, p5, p);
  t3 = sp_add(p3, p6, p);
  t4 = sp_sub(p1, p4, p);
  t5 = sp_sub(p2, p5, p);
  t6 = sp_sub(p3, p6, p);

  out[0 * ostride] = p0;
  out[1 * ostride] = t6;
  out[2 * ostride] = t4;
  out[3 * ostride] = t5;
  out[4 * ostride] = t2;
  out[5 * ostride] = t1;
  out[6 * ostride] = t3;
}

const nttconfig_t ntt7_config = 
{
  7,
  ntt7_get_num_const,
  ntt7_init,
  ntt7_run
};

/*---------------------- NTT8 --------------------------------------*/

static uint32_t 
ntt8_get_num_const(void)
{
  return 8;
}

static void
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
  out[5] = sp_mul(inv2, sp_sub(w1, w3, p), p, d);
  out[6] = w2;
  out[7] = sp_mul(inv2, sp_add(w1, w3, p), p, d);
}

static void
ntt8_run(spv_t out, spv_size_t ostride,
    	  spv_t in, spv_size_t istride,
	  sp_t p, sp_t d, spv_t ntt_const)
{
  sp_t x0, x1, x2, x3, x4, x5, x6, x7;
  sp_t t0, t1, t2, t3, t4, t5, t6, t7; 
  sp_t p0, p1, p2, p3, p4, p5, p6, p7;

  x0 = in[0 * istride];
  x1 = in[1 * istride];
  x2 = in[2 * istride];
  x3 = in[3 * istride];
  x4 = in[4 * istride];
  x5 = in[5 * istride];
  x6 = in[6 * istride];
  x7 = in[7 * istride];

  t0 = sp_add(x0, x4, p);
  t4 = sp_sub(x0, x4, p);
  t1 = sp_add(x1, x5, p);
  t5 = sp_sub(x1, x5, p);
  t2 = sp_add(x2, x6, p);
  t6 = sp_sub(x2, x6, p);
  t3 = sp_add(x3, x7, p);
  t7 = sp_sub(x3, x7, p);

  p0 = sp_add(t0, t2, p);
  p1 = sp_sub(t0, t2, p);
  p2 = sp_add(t1, t3, p);
  p3 = sp_sub(t1, t3, p);
  p4 = t4;
  p5 = sp_sub(t5, t7, p);
  p6 = t6;
  p7 = sp_add(t5, t7, p); 

  p3 = sp_mul(p3, ntt_const[3], p, d);
  p5 = sp_mul(p5, ntt_const[5], p, d);
  p6 = sp_mul(p6, ntt_const[6], p, d);
  p7 = sp_mul(p7, ntt_const[7], p, d);

  t0 = sp_add(p4, p6, p);
  t1 = sp_sub(p4, p6, p);
  t2 = sp_add(p5, p7, p);
  t3 = sp_sub(p5, p7, p);
  t4 = sp_add(t0, t2, p);
  t5 = sp_sub(t0, t2, p);
  t6 = sp_add(t1, t3, p);
  t7 = sp_sub(t1, t3, p);

  t0 = sp_add(p0, p2, p);
  t1 = sp_sub(p0, p2, p);
  t2 = sp_add(p1, p3, p);
  t3 = sp_sub(p1, p3, p);

  out[0 * ostride] = t0;
  out[1 * ostride] = t4;
  out[2 * ostride] = t2;
  out[3 * ostride] = t7;
  out[4 * ostride] = t1;
  out[5 * ostride] = t5;
  out[6 * ostride] = t3;
  out[7 * ostride] = t6;
}

const nttconfig_t ntt8_config = 
{
  8,
  ntt8_get_num_const,
  ntt8_init,
  ntt8_run
};

/*---------------------- NTT9 --------------------------------------*/

static uint32_t 
ntt9_get_num_const(void)
{
  return 11;
}

static void
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
ntt9_run(spv_t out, spv_size_t ostride,
    	  spv_t in, spv_size_t istride,
	  sp_t p, sp_t d, spv_t ntt_const)
{
  sp_t x0, x1, x2, x3, x4, x5, x6;
  sp_t     t1, t2, t3, t4, t5, t6, t7, t8;
  sp_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

  sp_t x0e, x1e, t0e, t1e, p0e, p1e, p2e;


  x0 = in[0 * istride];
  x1 = in[1 * istride];
  x2 = in[2 * istride];
  x0e = in[3 * istride];
  x3 = in[4 * istride];
  x6 = in[5 * istride];
  x1e = in[6 * istride];
  x5 = in[7 * istride];
  x4 = in[8 * istride];

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

  out[0 * ostride] = p0;
  out[1 * ostride] = t8;
  out[2 * ostride] = t3;
  out[3 * ostride] = t2;
  out[4 * ostride] = t4;
  out[5 * ostride] = t7;
  out[6 * ostride] = t1;
  out[7 * ostride] = t6;
  out[8 * ostride] = t5;
}

const nttconfig_t ntt9_config = 
{
  9,
  ntt9_get_num_const,
  ntt9_init,
  ntt9_run
};

/*---------------------- NTT15 --------------------------------------*/

static uint32_t 
ntt15_get_num_const(void)
{
  return ntt3_get_num_const() *
         ntt5_get_num_const();
}

static void
ntt15_init(spv_t out, sp_t p, sp_t d, 
	  sp_t primroot, sp_t order)
{
  uint32_t i, j;
  uint32_t num3 = ntt3_get_num_const();
  uint32_t num5 = ntt5_get_num_const();
  spv_t root3 = (spv_t)alloca(num3 * sizeof(sp_t));
  spv_t root5 = (spv_t)alloca(num5 * sizeof(sp_t));

  ntt3_init(root3, p, d, primroot, order);
  ntt5_init(root5, p, d, primroot, order);

  for (i = 0; i < num3; i++)
    for (j = 0; j < num5; j++)
      out[num5 * i + j] = sp_mul(root3[i], root5[j], p, d);
}

static void
ntt15_run(spv_t out, spv_size_t ostride,
    	  spv_t in, spv_size_t istride,
	  sp_t p, sp_t d, spv_t ntt_const)
{
  sp_t a00, a01, a02, a03, a04, 
       a05, a06, a07, a08, a09, 
       a10, a11, a12, a13, a14;

  {
    sp_t x00 = in[ 0 * istride];
    sp_t x01 = in[ 5 * istride];
    sp_t x02 = in[10 * istride];

    a01 = sp_add(x01, x02, p);
    a02 = sp_sub(x01, x02, p);

    a00 = sp_add(x00, a01, p);
  }
  {
    sp_t x03 = in[ 3 * istride];
    sp_t x04 = in[ 8 * istride];
    sp_t x05 = in[13 * istride];

    a04 = sp_add(x04, x05, p);
    a05 = sp_sub(x04, x05, p);

    a03 = sp_add(x03, a04, p);
  }
  {
    sp_t x08 = in[ 1 * istride];
    sp_t x06 = in[ 6 * istride];
    sp_t x07 = in[11 * istride];

    a07 = sp_add(x07, x08, p);
    a08 = sp_sub(x07, x08, p);

    a06 = sp_add(x06, a07, p);
  }
  {
    sp_t x11 = in[ 4 * istride];
    sp_t x09 = in[ 9 * istride];
    sp_t x10 = in[14 * istride];

    a10 = sp_add(x10, x11, p);
    a11 = sp_sub(x10, x11, p);

    a09 = sp_add(x09, a10, p);
  }
  {
    sp_t x13 = in[ 2 * istride];
    sp_t x14 = in[ 7 * istride];
    sp_t x12 = in[12 * istride];

    a13 = sp_add(x13, x14, p);
    a14 = sp_sub(x13, x14, p);

    a12 = sp_add(x12, a13, p);
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a00;
    b1 = a03;
    b4 = a06;
    b2 = a09;
    b3 = a12;

    c1 = sp_add(b1, b3, p);
    c3 = sp_sub(b1, b3, p);
    c2 = sp_add(b2, b4, p);
    c4 = sp_sub(b2, b4, p);

    b1 = sp_add(c1, c2, p);
    b2 = sp_sub(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_add(c3, c4, p);

    b0 = sp_add(b0, b1, p);

    b1 = sp_mul(b1, ntt_const[1], p, d);
    b2 = sp_mul(b2, ntt_const[2], p, d);
    b3 = sp_mul(b3, ntt_const[3], p, d);
    b4 = sp_mul(b4, ntt_const[4], p, d);
    b5 = sp_mul(b5, ntt_const[5], p, d);

    b1 = sp_add(b0, b1, p);

    c1 = sp_add(b1, b2, p);
    c2 = sp_sub(b1, b2, p);
    c3 = sp_add(b3, b5, p);
    c4 = sp_add(b4, b5, p);

    b1 = sp_add(c1, c3, p);
    b2 = sp_add(c2, c4, p);
    b3 = sp_sub(c1, c3, p);
    b4 = sp_sub(c2, c4, p);

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

    c1 = sp_add(b1, b3, p);
    c3 = sp_sub(b1, b3, p);
    c2 = sp_add(b2, b4, p);
    c4 = sp_sub(b2, b4, p);

    b1 = sp_add(c1, c2, p);
    b2 = sp_sub(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_add(c3, c4, p);

    b0 = sp_add(b0, b1, p);

    b0 = sp_mul(b0, ntt_const[6], p, d);
    b1 = sp_mul(b1, ntt_const[7], p, d);
    b2 = sp_mul(b2, ntt_const[8], p, d);
    b3 = sp_mul(b3, ntt_const[9], p, d);
    b4 = sp_mul(b4, ntt_const[10], p, d);
    b5 = sp_mul(b5, ntt_const[11], p, d);

    b1 = sp_add(b0, b1, p);

    c1 = sp_add(b1, b2, p);
    c2 = sp_sub(b1, b2, p);
    c3 = sp_add(b3, b5, p);
    c4 = sp_add(b4, b5, p);

    b1 = sp_add(c1, c3, p);
    b2 = sp_add(c2, c4, p);
    b3 = sp_sub(c1, c3, p);
    b4 = sp_sub(c2, c4, p);

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

    c1 = sp_add(b1, b3, p);
    c3 = sp_sub(b1, b3, p);
    c2 = sp_add(b2, b4, p);
    c4 = sp_sub(b2, b4, p);

    b1 = sp_add(c1, c2, p);
    b2 = sp_sub(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_add(c3, c4, p);

    b0 = sp_add(b0, b1, p);

    b0 = sp_mul(b0, ntt_const[12], p, d);
    b1 = sp_mul(b1, ntt_const[13], p, d);
    b2 = sp_mul(b2, ntt_const[14], p, d);
    b3 = sp_mul(b3, ntt_const[15], p, d);
    b4 = sp_mul(b4, ntt_const[16], p, d);
    b5 = sp_mul(b5, ntt_const[17], p, d);

    b1 = sp_add(b0, b1, p);

    c1 = sp_add(b1, b2, p);
    c2 = sp_sub(b1, b2, p);
    c3 = sp_add(b3, b5, p);
    c4 = sp_add(b4, b5, p);

    b1 = sp_add(c1, c3, p);
    b2 = sp_add(c2, c4, p);
    b3 = sp_sub(c1, c3, p);
    b4 = sp_sub(c2, c4, p);

    a02 = b0;
    a05 = b4;
    a08 = b3;
    a11 = b1;
    a14 = b2;
  }
  {
    sp_t x00, x01, x02;

    a01 = sp_add(a00, a01, p);

    x00 = a00;
    x01 = sp_add(a01, a02, p);
    x02 = sp_sub(a01, a02, p);

    out[ 0 * ostride] = x00;
    out[ 5 * ostride] = x02;
    out[10 * ostride] = x01;
  }
  {
    sp_t x03, x04, x05;

    a04 = sp_add(a03, a04, p);

    x03 = a03;
    x04 = sp_add(a04, a05, p);
    x05 = sp_sub(a04, a05, p);

    out[ 1 * ostride] = x04;
    out[ 6 * ostride] = x03;
    out[11 * ostride] = x05;
  }
  {
    sp_t x06, x07, x08;

    a07 = sp_add(a06, a07, p);

    x06 = a06;
    x07 = sp_add(a07, a08, p);
    x08 = sp_sub(a07, a08, p);

    out[ 2 * ostride] = x08;
    out[ 7 * ostride] = x07;
    out[12 * ostride] = x06;
  }
  {
    sp_t x09, x10, x11;

    a10 = sp_add(a09, a10, p);

    x09 = a09;
    x10 = sp_add(a10, a11, p);
    x11 = sp_sub(a10, a11, p);

    out[ 3 * ostride] = x09;
    out[ 8 * ostride] = x11;
    out[13 * ostride] = x10;
  }
  {
    sp_t x12, x13, x14;

    a13 = sp_add(a12, a13, p);

    x12 = a12;
    x13 = sp_add(a13, a14, p);
    x14 = sp_sub(a13, a14, p);

    out[ 4 * ostride] = x13;
    out[ 9 * ostride] = x12;
    out[14 * ostride] = x14;
  }

}

const nttconfig_t ntt15_config = 
{
  15,
  ntt15_get_num_const,
  ntt15_init,
  ntt15_run
};

/*---------------------- NTT16 --------------------------------------*/

static uint32_t 
ntt16_get_num_const(void)
{
  return 8;
}

static void
ntt16_init(spv_t out, sp_t p, sp_t d, 
	  sp_t primroot, sp_t order)
{
}

static void
ntt16_run(spv_t out, spv_size_t ostride,
    	  spv_t in, spv_size_t istride,
	  sp_t p, sp_t d, spv_t ntt_const)
{
}

const nttconfig_t ntt16_config = 
{
  16,
  ntt16_get_num_const,
  ntt16_init,
  ntt16_run
};

/*----------------------- generic ----------------------------------*/
/* r = NTT(x) */

static void
bfntt(spv_t r, spv_t x, spv_size_t len, 
     sp_t p, sp_t d, sp_t primroot, sp_t order)
{
  sp_t w0 = primroot;
  sp_t w_inc = 1;
  spv_size_t i, j, k;

  if (order != len)
    w0 = sp_pow(primroot, order / len, p, d);

  for (i = 0; i < len; i++)
    {
      sp_t accum = x[0];
      sp_t curr_w0 = w_inc;
      for (j = 1; j < len; j++)
	{
	  sp_t inc = sp_mul(x[j], curr_w0, p, d);
	  accum = sp_add(accum, inc, p);
	  curr_w0 = sp_mul(curr_w0, w_inc, p, d);
	}

      r[i] = accum;
      w_inc = sp_mul(w_inc, w0, p, d);
    }
}


const nttconfig_t * ntt_config[] = 
{
  &ntt2_config,
  &ntt3_config,
  &ntt4_config,
  &ntt5_config,
  &ntt6_config,
  &ntt7_config,
  &ntt8_config,
  &ntt9_config,
  &ntt15_config,
  &ntt16_config,
};

static void do_test(mpzspm_t mpzspm)
{
  sp_t p = mpzspm->spm[0]->sp;
  sp_t d = mpzspm->spm[0]->mul_c;
  sp_t primroot = mpzspm->spm[0]->prim_root;
  sp_t inv_primroot = mpzspm->spm[0]->inv_prim_root;
  sp_t order = mpzspm->max_ntt_size;
  const nttconfig_t * conf = ntt_config[7];
  spv_size_t len = conf->size;
  spv_size_t i;

  sp_t tmp[20];
  sp_t x[20];
  sp_t r[20];
  sp_t r2[20];

  spv_random(x, len, p);

  for (i = 0; i < len; i++)
    printf("in  %" PRIxsp "\n", x[i]);
  printf("\n");

  bfntt(r, x, len, p, d, primroot, order);

  for (i = 0; i < len; i++)
    printf("ref %" PRIxsp "\n", r[i]);
  printf("\n");

  conf->nttdata_init(tmp, p, d, primroot, order);
  conf->ntt_run(r2, 1, x, 1, p, d, tmp);

  for (i = 0; i < len; i++)
    printf("pfa %" PRIxsp "\n", r2[i]);
  printf("\n");
}

int main(int argc, char **argv)
{
  mpz_t x;
  spv_size_t len = 256*3*3*5*7;
  uint32_t bits = 300;
  mpzspm_t mpzspm;

  mpz_init_set_ui(x, 1);

#if 0
  bits = atol(argv[1]);
  len = atol(argv[2]);
#endif

  mpz_mul_2exp(x, x, bits);
  mpzspm = mpzspm_init((sp_t)len, x);

  if (mpzspm == NULL)
    {
      printf("crap\n");
      return 0;
    }

  do_test(mpzspm);

  mpzspm_clear(mpzspm);
  return 0;
}
