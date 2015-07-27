#include <stdio.h>
#include "ntt-impl.h"

static unsigned long
gcd (unsigned long a, unsigned long b)
{
  unsigned long t;

  while (b != 0UL)
    {
      t = a % b;
      a = b;
      b = t;
    }

  return a;
}

uint64_t
read_clock(void) 
{
#if defined(_MSC_VER)
	LARGE_INTEGER ret;
	QueryPerformanceCounter(&ret);
	return ret.QuadPart;
#else
	uint32_t lo, hi;
	asm volatile("rdtsc":"=d"(hi),"=a"(lo));
	return (uint64_t)hi << 32 | lo;
#endif
}

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

/*------------------------------------------------------------------*/
static void test_core(sp_t p, sp_t d, sp_t primroot, 
    			sp_t order, spv_size_t len,
			nttplan_t *plans, spv_size_t num_plans,
			nttdata_t *data)
{
  uint32_t fail = 0;
  uint32_t m, n;
  spv_t x = (spv_t)alloca(len * sizeof(sp_t));
  spv_t r = (spv_t)alloca(len * sizeof(sp_t));
  uint64_t start, stop, elapsed;

  spv_random(x, len, p);
  
  bfntt(r, x, len, p, d, primroot, order);

  ntt_build_passes(data, plans, num_plans, len, p, primroot, order, d);

  ntt_run(x, p, data);

  for (m = 0; m < len; m++)
     while (x[m] >= p)
     	x[m] -= p;

  for (m = 0; m < len; m++)
    {
      for (n = 0; n < len; n++)
	{
	  if (r[m] == x[n])
	   break;
	}
     if (n == len)
       break;
    }

  fail = (m < len);

  printf("%u ", len);
  for (m = 0; m < num_plans; m++)
    {
      printf("%s%u", m == 0 ? "" : "*", 
		plans[m].codelet_size);
      switch (plans[m].pass_type)
	{
	  case PASS_TYPE_DIRECT: printf("d"); break;
	  case PASS_TYPE_TWIDDLE: printf("t"); break;
	  case PASS_TYPE_PFA: printf("p"); break;
	}
    }

  if (fail)
    printf(" fail");

  elapsed = (uint64_t)(-1);
  for (m = 0; m < 10; m++)
    {
      start = read_clock();
      for (n = 0; n < 10; n++)
	{
	  ntt_run(x, p, data);
	}
      stop = read_clock();
      if (stop - start < elapsed)
		elapsed = stop - start;
    }
  printf(" %lf clocks/point %.2lf usec\n", 
		(double)elapsed / 10 / len,
		(double)elapsed / 10 / 2500
		);

  ntt_reset(data);
}

/*------------------------------------------------------------------*/
static void do_direct_test(mpzspm_t mpzspm)
{
  uint32_t i;
  sp_t p = mpzspm->spm[0]->sp;
  sp_t d = mpzspm->spm[0]->mul_c;
  sp_t primroot = mpzspm->spm[0]->primroot;
  sp_t order = mpzspm->ntt_size;
  nttdata_t *nttdata = (nttdata_t *)mpzspm->spm[0]->ntt_data;
  nttplan_t plans[3];

  plans[0].pass_type = PASS_TYPE_DIRECT;

  for (i = 0; i < nttdata->num_codelets; i++)
    {
      codelet_data_t * c = nttdata->codelets + i;
      spv_size_t len = c->config->size;

      plans[0].codelet_size = len;
      test_core(p, d, primroot, order, len, plans, 1, nttdata);
    }
}

/*------------------------------------------------------------------*/
static void do_twiddle_test(mpzspm_t mpzspm)
{
  uint32_t i, j, k;
  sp_t p = mpzspm->spm[0]->sp;
  sp_t d = mpzspm->spm[0]->mul_c;
  sp_t primroot = mpzspm->spm[0]->primroot;
  sp_t order = mpzspm->ntt_size;
  nttdata_t *nttdata = (nttdata_t *)mpzspm->spm[0]->ntt_data;
  nttplan_t plans[3];

  /* transform pairs */

  plans[0].pass_type = PASS_TYPE_TWIDDLE;
  plans[1].pass_type = PASS_TYPE_DIRECT;

  for (i = 0; i < nttdata->num_codelets; i++)
    {
      for (j = 0; j < nttdata->num_codelets; j++)
	{
	  codelet_data_t * c1 = nttdata->codelets + i;
	  codelet_data_t * c2 = nttdata->codelets + j;
	  spv_size_t len = c1->config->size * c2->config->size;

	  if (order % len != 0)
	    continue;

	  plans[0].codelet_size = c1->config->size;
	  plans[1].codelet_size = c2->config->size;
	  test_core(p, d, primroot, order, len, plans, 2, nttdata);
	}
    }

  /* transform triplets */

  plans[0].pass_type = PASS_TYPE_TWIDDLE;
  plans[1].pass_type = PASS_TYPE_TWIDDLE;
  plans[2].pass_type = PASS_TYPE_DIRECT;

  for (i = 0; i < nttdata->num_codelets; i++)
    {
      for (j = 0; j < nttdata->num_codelets; j++)
	{
      	  for (k = 0; k < nttdata->num_codelets; k++)
    	    {
    	      codelet_data_t * c1 = nttdata->codelets + i;
    	      codelet_data_t * c2 = nttdata->codelets + j;
    	      codelet_data_t * c3 = nttdata->codelets + k;
    	      spv_size_t len = c1->config->size * c2->config->size * c3->config->size;

    	      if (order % len != 0)
    		continue;

    	      plans[0].codelet_size = c1->config->size;
    	      plans[1].codelet_size = c2->config->size;
    	      plans[2].codelet_size = c3->config->size;
    	      test_core(p, d, primroot, order, len, plans, 3, nttdata);
	    }
	}
    }
}

/*------------------------------------------------------------------*/
static void do_pfa_test(mpzspm_t mpzspm)
{
  uint32_t i, j, k;
  sp_t p = mpzspm->spm[0]->sp;
  sp_t d = mpzspm->spm[0]->mul_c;
  sp_t primroot = mpzspm->spm[0]->primroot;
  sp_t order = mpzspm->ntt_size;
  nttdata_t *nttdata = (nttdata_t *)mpzspm->spm[0]->ntt_data;
  nttplan_t plans[3];

  plans[0].pass_type = PASS_TYPE_PFA;
  plans[1].pass_type = PASS_TYPE_PFA;
  plans[2].pass_type = PASS_TYPE_PFA;

  /* transform pairs */

  for (i = 0; i < nttdata->num_codelets - 1; i++)
    {
      for (j = i + 1; j < nttdata->num_codelets; j++)
	{
	  codelet_data_t * c1 = nttdata->codelets + i;
	  codelet_data_t * c2 = nttdata->codelets + j;
	  spv_size_t len = c1->config->size * c2->config->size;

	  if (gcd(c1->config->size, c2->config->size) != 1)
	    continue;

	  plans[0].codelet_size = c1->config->size;
	  plans[1].codelet_size = c2->config->size;
	  test_core(p, d, primroot, order, len, plans, 2, nttdata);
	}
    }

  /* transform triplets */

  for (i = 0; i < nttdata->num_codelets - 2; i++)
    {
      for (j = i + 1; j < nttdata->num_codelets - 1; j++)
	{
      	  for (k = j + 1; k < nttdata->num_codelets; k++)
    	    {
    	      codelet_data_t * c1 = nttdata->codelets + i;
    	      codelet_data_t * c2 = nttdata->codelets + j;
    	      codelet_data_t * c3 = nttdata->codelets + k;
    	      spv_size_t len = c1->config->size * c2->config->size * c3->config->size;

    	      if (gcd(c1->config->size, c2->config->size) != 1 ||
    		  gcd(c1->config->size, c3->config->size) != 1 ||
    		  gcd(c2->config->size, c3->config->size) != 1)
    		continue;

    	      plans[0].codelet_size = c1->config->size;
    	      plans[1].codelet_size = c2->config->size;
    	      plans[2].codelet_size = c3->config->size;
    	      test_core(p, d, primroot, order, len, plans, 3, nttdata);
	    }
	}
    }
}

/*------------------------------------------------------------------*/
int main(int argc, char **argv)
{
  mpz_t x;
  spv_size_t len = 1024*3*3*5*5*7*7;
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

  do_direct_test(mpzspm);
  do_pfa_test(mpzspm);
  do_twiddle_test(mpzspm);

  mpzspm_clear(mpzspm);
  return 0;
}
