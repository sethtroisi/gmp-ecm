#include <stdio.h>
#include "libntt.h"

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

#define NUM_CODELETS (sizeof(codelet_sizes) / sizeof(codelet_sizes[0]))

static const uint32_t codelet_sizes[] = {
  2, 3, 4, 5, 7, 8, 9, 15, 16, 35, 40,
};

#define MAX_PLAN_GROUPS 4

nttplangroup_t plan_list[MAX_PLAN_GROUPS];

/*------------------------------------------------------------------*/
static void test_core(nttwork_t nttwork, uint32_t size, uint32_t verify)
{
  uint32_t i, j;
  double elapsed;

  printf("%u ", size);

  for (i = 0; i < nttwork->mpzspm_num; i++)
    {
      printf("w%u:", nttwork->nttinit[i]->sp_bits);

      for (j = 0; j < plan_list[i].num_plans; j++)
	{
	  printf("%u", plan_list[i].plans[j].codelet_size);
	  switch(plan_list[i].plans[j].pass_type)
	    {
	      case PASS_TYPE_DIRECT: printf("d"); break;
	      case PASS_TYPE_PFA: printf("p"); break;
	      case PASS_TYPE_TWIDDLE: printf("t"); break;
	    }
	  printf("(%u)", plan_list[i].plans[j].group_type);
	}
      printf(" ");
    }

  nttwork_ntt_init(nttwork, plan_list);
  nttwork_random(nttwork);
  elapsed = nttwork_ntt_test(nttwork, verify);
  nttwork_ntt_reset(nttwork);
  printf("%.3lf clocks/point\n", elapsed);
}

/*------------------------------------------------------------------*/
static void test1(nttwork_t nttwork, uint32_t verify)
{
  uint32_t imax[MAX_PLAN_GROUPS] = {0};
  uint32_t i0, i1, i2, i3;
  uint32_t j;

  for (j = 0; j < nttwork->mpzspm_num; j++)
    {
      imax[j] = nttwork->nttinit[j]->mpzspm_get_num_groups();
      plan_list[j].num_plans = 1;
    }

  for (j = 0; j < NUM_CODELETS; j++)
    {
      uint32_t size = codelet_sizes[j];

      if (nttwork->max_ntt_size % size != 0)
	continue;

  for (i0 = imax[0]; (int32_t)i0 >= 0; i0--)
    {
      plan_list[0].plans[0].codelet_size = codelet_sizes[j];
      plan_list[0].plans[0].group_type = i0 - 1;
      plan_list[0].plans[0].pass_type = PASS_TYPE_DIRECT;

  for (i1 = imax[1]; (int32_t)i1 >= 0; i1--)
    {
      plan_list[1].plans[0].codelet_size = codelet_sizes[j];
      plan_list[1].plans[0].group_type = i1 - 1;
      plan_list[1].plans[0].pass_type = PASS_TYPE_DIRECT;

  for (i2 = imax[2]; (int32_t)i2 >= 0; i2--)
    {
      plan_list[2].plans[0].codelet_size = codelet_sizes[j];
      plan_list[2].plans[0].group_type = i2 - 1;
      plan_list[2].plans[0].pass_type = PASS_TYPE_DIRECT;

  for (i3 = imax[3]; (int32_t)i3 >= 0; i3--)
    {
      plan_list[3].plans[0].codelet_size = codelet_sizes[j];
      plan_list[3].plans[0].group_type = i3 - 1;
      plan_list[3].plans[0].pass_type = PASS_TYPE_DIRECT;

      test_core(nttwork, size, verify);
    }}}}

    }
}

#if 0
/*------------------------------------------------------------------*/
static void test2(mpzspm_t mpzspm, uint32_t order)
{
  uint32_t i, j, k, m;
  uint32_t num_groups = X(ntt_master_group_list_size);
  const nttgroup_t **groups = X(ntt_master_group_list);
  sp_t p = mpzspm->spm[0]->sp;
  sp_t d = mpzspm->spm[0]->mul_c;
  sp_t primroot = mpzspm->spm[0]->primroot;
  nttdata_t *nttdata = &mpzspm->spm[0]->ntt_data;
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
static void test3(mpzspm_t mpzspm, uint32_t order)
{
  uint32_t i, j, k, m, n, q;
  uint32_t num_groups = X(ntt_master_group_list_size);
  const nttgroup_t **groups = X(ntt_master_group_list);
  sp_t p = mpzspm->spm[0]->sp;
  sp_t d = mpzspm->spm[0]->mul_c;
  sp_t primroot = mpzspm->spm[0]->primroot;
  nttdata_t *nttdata = &mpzspm->spm[0]->ntt_data;
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
#endif

/*------------------------------------------------------------------*/
int main(int argc, char **argv)
{
  nttwork_t nttwork;
  mpz_t modulus;
  uint32_t ntt_size = 1024*3*3*5*5*7*7;
  uint32_t bits = 300;
  uint32_t sp_bits[3];
  uint32_t sp_bits_choices = 0;
  uint32_t interleaved = 0;
  uint32_t verify = 1;

#if GMP_LIMB_BITS == 32
  sp_bits[0] = 30;
  sp_bits[1] = 31;
  sp_bits[2] = 50;
  sp_bits_choices = 3;
#else
  sp_bits[0] = 62;
  sp_bits_choices = 1;
#endif

  mpz_init_set_ui(modulus, 1);
  mpz_mul_2exp(modulus, modulus, bits);

  nttwork = nttwork_init(ntt_size, modulus, interleaved, 
      			sp_bits, sp_bits_choices);

  test1(nttwork, verify);

  nttwork_clear(nttwork);
  mpz_clear(modulus);
  return 0;
}
