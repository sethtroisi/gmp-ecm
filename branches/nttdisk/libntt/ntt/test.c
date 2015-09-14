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

static uint64_t
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

  X(spv_random)(x, len, p);
  
  X(ntt_build_passes)(data, plans, num_plans, len, p, primroot, order, d);
#if 0
  bfntt(r, x, len, p, d, primroot, order);

  X(ntt_run)(x, p, data);

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
#endif
  printf("%u ", len);
  for (m = 0; m < num_plans; m++)
    {
      printf("%s%u", m == 0 ? "" : "*", 
		plans[m].codelet_size);
      switch (plans[m].pass_type)
	{
	  case PASS_TYPE_DIRECT: printf("d(%u)", plans[m].group_type); break;
	  case PASS_TYPE_TWIDDLE: printf("t(%u)", plans[m].group_type); break;
	  case PASS_TYPE_PFA: printf("p(%u)", plans[m].group_type); break;
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
	  X(ntt_run)(x, p, data);
	}
      stop = read_clock();
      if (stop - start < elapsed)
		elapsed = stop - start;
    }
  printf(" %lf clocks/point %.2lf usec " SP_NAME_SUFFIX_STR "\n",
		(double)elapsed / 10 / len,
		(double)elapsed / 10 / 2500
		);

  X(ntt_reset)(data);
}

/*------------------------------------------------------------------*/
static void test1(mpzspm_t mpzspm)
{
  uint32_t i, j;
  uint32_t num_groups = X(ntt_master_group_list_size);
  const nttgroup_t **groups = X(ntt_master_group_list);
  sp_t p = mpzspm->spm[0]->sp;
  sp_t d = mpzspm->spm[0]->mul_c;
  sp_t primroot = mpzspm->spm[0]->primroot;
  uint32_t order = mpzspm->ntt_size;
  nttdata_t *nttdata = (nttdata_t *)mpzspm->spm[0]->ntt_data;
  nttplan_t plans[1];

  for (i = 0; i < num_groups; i++)
    {
      uint32_t num_codelets = groups[i]->num_transforms;
      const nttconfig_t **codelets = groups[i]->get_transform_list();

      for (j = 0; j < num_codelets; j++)
	{
	  const nttconfig_t * c = codelets[j];

	  plans[0].codelet_size = c->size;
	  plans[0].group_type = i;
	  plans[0].pass_type = PASS_TYPE_DIRECT;
	  test_core(p, d, primroot, order, c->size, plans, 1, nttdata);
	}
    }
}

/*------------------------------------------------------------------*/
static void test2(mpzspm_t mpzspm)
{
  uint32_t i, j, k, m;
  uint32_t num_groups = X(ntt_master_group_list_size);
  const nttgroup_t **groups = X(ntt_master_group_list);
  sp_t p = mpzspm->spm[0]->sp;
  sp_t d = mpzspm->spm[0]->mul_c;
  sp_t primroot = mpzspm->spm[0]->primroot;
  uint32_t order = mpzspm->ntt_size;
  nttdata_t *nttdata = (nttdata_t *)mpzspm->spm[0]->ntt_data;
  nttplan_t plans[3];

  /* transform pairs */

  for (i = 0; i < num_groups; i++)
    {
      uint32_t num_codelets0 = groups[i]->num_transforms;
      const nttconfig_t **codelets0 = groups[i]->get_transform_list();

      for (j = 0; j < num_groups; j++)
	{
	  uint32_t num_codelets1 = groups[j]->num_transforms;
	  const nttconfig_t **codelets1 = groups[j]->get_transform_list();

	  for (k = 0; k < num_codelets0; k++)
	    {
	      for (m = 0; m < num_codelets1; m++)
		{
		  const nttconfig_t * c1 = codelets0[k];
		  const nttconfig_t * c2 = codelets1[m];
		  spv_size_t len = c1->size * c2->size;

		  if (order % len != 0)
		    continue;

		  plans[0].codelet_size = c1->size;
		  plans[0].group_type = i;
		  plans[0].pass_type = PASS_TYPE_TWIDDLE;
		  plans[1].codelet_size = c2->size;
		  plans[1].group_type = j;
		  plans[1].pass_type = PASS_TYPE_DIRECT;
		  test_core(p, d, primroot, order, len, plans, 2, nttdata);

		  if (gcd(c1->size, c2->size) != 1)
		    continue;

		  plans[0].pass_type = PASS_TYPE_PFA;
		  plans[1].pass_type = PASS_TYPE_PFA;
		  test_core(p, d, primroot, order, len, plans, 2, nttdata);
		}
	    }
	}
    }
}

/*------------------------------------------------------------------*/
static void test3(mpzspm_t mpzspm)
{
  uint32_t i, j, k, m, n, q;
  uint32_t num_groups = X(ntt_master_group_list_size);
  const nttgroup_t **groups = X(ntt_master_group_list);
  sp_t p = mpzspm->spm[0]->sp;
  sp_t d = mpzspm->spm[0]->mul_c;
  sp_t primroot = mpzspm->spm[0]->primroot;
  uint32_t order = mpzspm->ntt_size;
  nttdata_t *nttdata = (nttdata_t *)mpzspm->spm[0]->ntt_data;
  nttplan_t plans[3];

  /* transform triplets */

  for (i = 0; i < num_groups; i++)
    {
      uint32_t num_codelets0 = groups[i]->num_transforms;
      const nttconfig_t **codelets0 = groups[i]->get_transform_list();

      for (j = 0; j < num_groups; j++)
	{
	  uint32_t num_codelets1 = groups[j]->num_transforms;
	  const nttconfig_t **codelets1 = groups[j]->get_transform_list();

	  for (k = 0; k < num_groups; k++)
	    {
	      uint32_t num_codelets2 = groups[k]->num_transforms;
	      const nttconfig_t **codelets2 = groups[k]->get_transform_list();

	      for (m = 0; m < num_codelets0; m++)
		{
		  for (n = 0; n < num_codelets1; n++)
		    {
		      for (q = 0; q < num_codelets2; q++)
			{
			  const nttconfig_t * c1 = codelets0[m];
			  const nttconfig_t * c2 = codelets1[n];
			  const nttconfig_t * c3 = codelets2[q];
			  spv_size_t len = c1->size * c2->size * c3->size;

			  if (order % len != 0)
			    continue;

			  plans[0].codelet_size = c1->size;
			  plans[0].group_type = i;
			  plans[0].pass_type = PASS_TYPE_TWIDDLE;
			  plans[1].codelet_size = c2->size;
			  plans[1].group_type = j;
			  plans[1].pass_type = PASS_TYPE_TWIDDLE;
			  plans[2].codelet_size = c3->size;
			  plans[2].group_type = k;
			  plans[2].pass_type = PASS_TYPE_DIRECT;
			  test_core(p, d, primroot, order, len, plans, 3, nttdata);

	    		  if (gcd(c1->size, c2->size) != 1 ||
		    	      gcd(c1->size, c3->size) != 1 ||
			      gcd(c2->size, c3->size) != 1)
			    continue;

			  plans[0].pass_type = PASS_TYPE_PFA;
			  plans[1].pass_type = PASS_TYPE_PFA;
			  plans[2].pass_type = PASS_TYPE_PFA;
			  test_core(p, d, primroot, order, len, plans, 3, nttdata);
			}
		    }
		}
	    }
	}
    }
}

/*------------------------------------------------------------------*/
int X(test_main)(int argc, char **argv)
{
  mpz_t x;
  uint32_t len = 1024*3*3*5*5*7*7;
  uint32_t bits = 300;
  mpzspm_t mpzspm;

#if SP_NUMB_BITS == 50
  fpu_precision_t old_prec = 0;
  if (!fpu_precision_is_ieee())
    old_prec = fpu_set_precision_ieee();
#endif

  mpz_init_set_ui(x, 1);

#if 0
  bits = atol(argv[1]);
  len = atol(argv[2]);
#endif

  mpz_mul_2exp(x, x, bits);
  mpzspm = X(mpzspm_init)(len, x);

  if (mpzspm == NULL)
    {
      printf("crap\n");
      return 0;
    }

  test1(mpzspm);
  test2(mpzspm);
  test3(mpzspm);

  X(mpzspm_clear)(mpzspm);

#if SP_NUMB_BITS == 50
  if (old_prec != 0)
    fpu_clear_precision(old_prec);
#endif

  return 0;
}
