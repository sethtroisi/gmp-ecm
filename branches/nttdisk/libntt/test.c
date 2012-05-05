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
	asm("rdtsc":"=d"(hi),"=a"(lo));
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


static void do_test(mpzspm_t mpzspm)
{
  uint32_t i, j, k, m, n;
  sp_t p = mpzspm->spm[0]->sp;
  sp_t d = mpzspm->spm[0]->mul_c;
  sp_t primroot = mpzspm->spm[0]->primroot;
  sp_t order = mpzspm->ntt_size;
  nttdata_t *nttdata = (nttdata_t *)mpzspm->spm[0]->ntt_data;

  /* transform pairs */

  for (i = 0; i < nttdata->num_codelets - 1; i++)
    {
      for (j = i + 1; j < nttdata->num_codelets; j++)
	{
	  codelet_data_t * pfa1 = nttdata->codelets + i;
	  codelet_data_t * pfa2 = nttdata->codelets + j;
	  const nttconfig_t *config1 = pfa1->config;
	  const nttconfig_t *config2 = pfa2->config;
	  spv_size_t len = config1->size * config2->size;
	  spv_t x = (spv_t)alloca(len * sizeof(sp_t));
	  spv_t r = (spv_t)alloca(len * sizeof(sp_t));
	  uint64_t start, stop;

	  if (gcd(config1->size, config2->size) != 1)
	    continue;

	  spv_random(x, len, p);
  
	  bfntt(r, x, len, p, d, primroot, order);

	  config1->ntt_pfa_run(x, 1, len / config1->size, 
	      			p, d, pfa1->ntt_const);
	  config2->ntt_pfa_run(x, 1, len / config2->size, 
	      			p, d, pfa2->ntt_const);

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

	  if (m != len)
	    printf("%u*%u: fail", config1->size, config2->size);
	  else
	    printf("%u*%u:", config1->size, config2->size);

	  start = read_clock();
	  for (m = 0; m < 100; m++)
	    {
	      config1->ntt_pfa_run(x, 1, len / config1->size, 
		    		    p, d, pfa1->ntt_const);
	      config2->ntt_pfa_run(x, 1, len / config2->size, 
			    	    p, d, pfa2->ntt_const);
	    }
	  stop = read_clock();
	  printf("\t%lf clocks/point\n", (double)(stop - start) / 100 / len);
	}
    }

  /* transform triplets */
  for (i = 0; i < nttdata->num_codelets - 2; i++)
    {
      for (j = i + 1; j < nttdata->num_codelets - 1; j++)
	{
      for (k = j + 1; k < nttdata->num_codelets; k++)
	{
	  codelet_data_t * pfa1 = nttdata->codelets + i;
	  codelet_data_t * pfa2 = nttdata->codelets + j;
	  codelet_data_t * pfa3 = nttdata->codelets + k;
	  const nttconfig_t *config1 = pfa1->config;
	  const nttconfig_t *config2 = pfa2->config;
	  const nttconfig_t *config3 = pfa3->config;
	  spv_size_t len = config1->size * config2->size * config3->size;
	  spv_t x = (spv_t)alloca(len * sizeof(sp_t));
	  spv_t r = (spv_t)alloca(len * sizeof(sp_t));
	  uint64_t start, stop;

	  if (gcd(config1->size, config2->size) != 1 ||
	      gcd(config1->size, config3->size) != 1 ||
	      gcd(config2->size, config3->size) != 1)
	    continue;

	  spv_random(x, len, p);
  
	  bfntt(r, x, len, p, d, primroot, order);

	  config1->ntt_pfa_run(x, 1, len / config1->size, 
	      			p, d, pfa1->ntt_const);
	  config2->ntt_pfa_run(x, 1, len / config2->size, 
	      			p, d, pfa2->ntt_const);
	  config3->ntt_pfa_run(x, 1, len / config3->size, 
	      			p, d, pfa3->ntt_const);

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

	  if (m != len)
	    printf("%u*%u*%u: fail", config1->size, config2->size, config3->size);
	  else
	    printf("%u*%u*%u:", config1->size, config2->size, config3->size);

	  start = read_clock();
	  for (m = 0; m < 100; m++)
	    {
	      config1->ntt_pfa_run(x, 1, len / config1->size, 
		    		    p, d, pfa1->ntt_const);
	      config2->ntt_pfa_run(x, 1, len / config2->size, 
			    	    p, d, pfa2->ntt_const);
	      config3->ntt_pfa_run(x, 1, len / config3->size, 
			    	    p, d, pfa3->ntt_const);
	    }
	  stop = read_clock();
	  printf("\t%lf clocks/point\n", (double)(stop - start) / 100 / len);
	}
	}
    }
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
