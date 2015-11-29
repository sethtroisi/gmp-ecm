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
  printf("  %.3lf clocks/point\n", elapsed);
}

/*------------------------------------------------------------------*/
static void test1(nttwork_t nttwork, uint32_t verify)
{
  uint32_t imax[MAX_PLAN_GROUPS] = {1, 1, 1, 1};
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

  for (i0 = 0; i0 < imax[0]; i0++)
    {
      plan_list[0].plans[0].codelet_size = codelet_sizes[j];
      plan_list[0].plans[0].group_type = i0;
      plan_list[0].plans[0].pass_type = PASS_TYPE_DIRECT;

  for (i1 = 0; i1 < imax[1]; i1++)
    {
      plan_list[1].plans[0].codelet_size = codelet_sizes[j];
      plan_list[1].plans[0].group_type = i1;
      plan_list[1].plans[0].pass_type = PASS_TYPE_DIRECT;

  for (i2 = 0; i2 < imax[2]; i2++)
    {
      plan_list[2].plans[0].codelet_size = codelet_sizes[j];
      plan_list[2].plans[0].group_type = i2;
      plan_list[2].plans[0].pass_type = PASS_TYPE_DIRECT;

  for (i3 = 0; i3 < imax[3]; i3++)
    {
      plan_list[3].plans[0].codelet_size = codelet_sizes[j];
      plan_list[3].plans[0].group_type = i3;
      plan_list[3].plans[0].pass_type = PASS_TYPE_DIRECT;

      test_core(nttwork, size, verify);
    }}}}

    }
}

/*------------------------------------------------------------------*/
static void test2(nttwork_t nttwork, uint32_t verify)
{
  uint32_t imax[MAX_PLAN_GROUPS] = {1, 1, 1, 1};
  uint32_t i0, i1, i2, i3;
  uint32_t j, k, m;

  for (j = 0; j < nttwork->mpzspm_num; j++)
    {
      imax[j] = nttwork->nttinit[j]->mpzspm_get_num_groups();
      plan_list[j].num_plans = 2;
    }

  for (j = 0; j < NUM_CODELETS; j++)
    {
  for (k = 0; k < NUM_CODELETS; k++)
    {
      uint32_t size = codelet_sizes[j] * codelet_sizes[k];

      if (nttwork->max_ntt_size % size != 0)
	continue;

  for (i0 = 0; i0 < imax[0]; i0++)
    {
      plan_list[0].plans[0].codelet_size = codelet_sizes[j];
      plan_list[0].plans[0].group_type = i0;
      plan_list[0].plans[0].pass_type = PASS_TYPE_TWIDDLE;

      plan_list[0].plans[1].codelet_size = codelet_sizes[k];
      plan_list[0].plans[1].group_type = i0;
      plan_list[0].plans[1].pass_type = PASS_TYPE_DIRECT;

  for (i1 = 0; i1 < imax[1]; i1++)
    {
      plan_list[1].plans[0].codelet_size = codelet_sizes[j];
      plan_list[1].plans[0].group_type = i1;
      plan_list[1].plans[0].pass_type = PASS_TYPE_TWIDDLE;

      plan_list[1].plans[1].codelet_size = codelet_sizes[k];
      plan_list[1].plans[1].group_type = i1;
      plan_list[1].plans[1].pass_type = PASS_TYPE_DIRECT;

  for (i2 = 0; i2 < imax[2]; i2++)
    {
      plan_list[2].plans[0].codelet_size = codelet_sizes[j];
      plan_list[2].plans[0].group_type = i2;
      plan_list[2].plans[0].pass_type = PASS_TYPE_TWIDDLE;

      plan_list[2].plans[1].codelet_size = codelet_sizes[k];
      plan_list[2].plans[1].group_type = i2;
      plan_list[2].plans[1].pass_type = PASS_TYPE_DIRECT;

  for (i3 = 0; i3 < imax[3]; i3++)
    {
      plan_list[3].plans[0].codelet_size = codelet_sizes[j];
      plan_list[3].plans[0].group_type = i3;
      plan_list[3].plans[0].pass_type = PASS_TYPE_TWIDDLE;

      plan_list[3].plans[1].codelet_size = codelet_sizes[k];
      plan_list[3].plans[1].group_type = i3;
      plan_list[3].plans[1].pass_type = PASS_TYPE_DIRECT;

      test_core(nttwork, size, verify);

      if (j < k && gcd(codelet_sizes[j], codelet_sizes[k]) == 1)
	{
	  for (m = 0; m < nttwork->mpzspm_num; m++)
	    {
	      plan_list[m].plans[0].pass_type = PASS_TYPE_PFA;
	      plan_list[m].plans[1].pass_type = PASS_TYPE_PFA;
	    }
	  test_core(nttwork, size, verify);
	}

    }}}}

    }
    }
}

/*------------------------------------------------------------------*/
static void test3(nttwork_t nttwork, uint32_t verify)
{
  uint32_t imax[MAX_PLAN_GROUPS] = {1, 1, 1, 1};
  uint32_t i0, i1, i2, i3;
  uint32_t j, k, m, n;

  for (j = 0; j < nttwork->mpzspm_num; j++)
    {
      imax[j] = nttwork->nttinit[j]->mpzspm_get_num_groups();
      plan_list[j].num_plans = 3;
    }

  for (j = 0; j < NUM_CODELETS; j++)
    {
  for (k = 0; k < NUM_CODELETS; k++)
    {
  for (m = 0; m < NUM_CODELETS; m++)
    {
      uint32_t size = codelet_sizes[j] * codelet_sizes[k] *
			codelet_sizes[m];

      if (nttwork->max_ntt_size % size != 0)
	continue;

  for (i0 = 0; i0 < imax[0]; i0++)
    {
      plan_list[0].plans[0].codelet_size = codelet_sizes[j];
      plan_list[0].plans[0].group_type = i0;
      plan_list[0].plans[0].pass_type = PASS_TYPE_TWIDDLE;

      plan_list[0].plans[1].codelet_size = codelet_sizes[k];
      plan_list[0].plans[1].group_type = i0;
      plan_list[0].plans[1].pass_type = PASS_TYPE_TWIDDLE;

      plan_list[0].plans[2].codelet_size = codelet_sizes[m];
      plan_list[0].plans[2].group_type = i0;
      plan_list[0].plans[2].pass_type = PASS_TYPE_DIRECT;

  for (i1 = 0; i1 < imax[1]; i1++)
    {
      plan_list[1].plans[0].codelet_size = codelet_sizes[j];
      plan_list[1].plans[0].group_type = i1;
      plan_list[1].plans[0].pass_type = PASS_TYPE_TWIDDLE;

      plan_list[1].plans[1].codelet_size = codelet_sizes[k];
      plan_list[1].plans[1].group_type = i1;
      plan_list[1].plans[1].pass_type = PASS_TYPE_TWIDDLE;

      plan_list[1].plans[2].codelet_size = codelet_sizes[m];
      plan_list[1].plans[2].group_type = i1;
      plan_list[1].plans[2].pass_type = PASS_TYPE_DIRECT;

  for (i2 = 0; i2 < imax[2]; i2++)
    {
      plan_list[2].plans[0].codelet_size = codelet_sizes[j];
      plan_list[2].plans[0].group_type = i2;
      plan_list[2].plans[0].pass_type = PASS_TYPE_TWIDDLE;

      plan_list[2].plans[1].codelet_size = codelet_sizes[k];
      plan_list[2].plans[1].group_type = i2;
      plan_list[2].plans[1].pass_type = PASS_TYPE_TWIDDLE;

      plan_list[2].plans[2].codelet_size = codelet_sizes[m];
      plan_list[2].plans[2].group_type = i2;
      plan_list[2].plans[2].pass_type = PASS_TYPE_DIRECT;

  for (i3 = 0; i3 < imax[3]; i3++)
    {
      plan_list[3].plans[0].codelet_size = codelet_sizes[j];
      plan_list[3].plans[0].group_type = i3;
      plan_list[3].plans[0].pass_type = PASS_TYPE_TWIDDLE;

      plan_list[3].plans[1].codelet_size = codelet_sizes[k];
      plan_list[3].plans[1].group_type = i3;
      plan_list[3].plans[1].pass_type = PASS_TYPE_TWIDDLE;

      plan_list[3].plans[2].codelet_size = codelet_sizes[m];
      plan_list[3].plans[2].group_type = i3;
      plan_list[3].plans[2].pass_type = PASS_TYPE_DIRECT;

      test_core(nttwork, size, verify);

      if (j < k && k < m && 
	  gcd(codelet_sizes[j], codelet_sizes[k]) == 1 &&
	  gcd(codelet_sizes[k], codelet_sizes[m]) == 1 &&
	  gcd(codelet_sizes[j], codelet_sizes[m]) == 1)
	{
	  for (n = 0; n < nttwork->mpzspm_num; n++)
	    {
	      plan_list[n].plans[0].pass_type = PASS_TYPE_PFA;
	      plan_list[n].plans[1].pass_type = PASS_TYPE_PFA;
	      plan_list[n].plans[2].pass_type = PASS_TYPE_PFA;
	    }
	  test_core(nttwork, size, verify);
	}

      if (k < m && 
	  gcd(codelet_sizes[k], codelet_sizes[m]) == 1)
	{
	  for (n = 0; n < nttwork->mpzspm_num; n++)
	    {
	      plan_list[n].plans[0].pass_type = PASS_TYPE_TWIDDLE;
	      plan_list[n].plans[1].pass_type = PASS_TYPE_PFA;
	      plan_list[n].plans[2].pass_type = PASS_TYPE_PFA;
	    }
	  test_core(nttwork, size, verify);
	}

    }}}}

    }
    }
    }
}

/*------------------------------------------------------------------*/
int main(int argc, char **argv)
{
  nttwork_t nttwork;
  mpz_t modulus;
  uint32_t ntt_size = 1024*3*3*5*5*7*7;
  uint32_t bits = 300;
  uint32_t sp_bits[4];
  uint32_t sp_bits_choices = 0;
  uint32_t interleaved = 1;
  uint32_t verify = 0;

#if GMP_LIMB_BITS == 32
  sp_bits[0] = 30;
  sp_bits[1] = 31;
  sp_bits[2] = 50;
  sp_bits[3] = 62;
  sp_bits_choices = 4;
#else
  sp_bits[0] = 30;
  sp_bits[1] = 31;
  sp_bits[2] = 62;
  sp_bits_choices = 3;
#endif

  mpz_init_set_ui(modulus, 1);
  mpz_mul_2exp(modulus, modulus, bits);

  nttwork = nttwork_init(ntt_size, modulus, interleaved, 
      			sp_bits, sp_bits_choices);

  if (nttwork == NULL)
    {
      printf("failed\n");
      return -1;
    }

  test1(nttwork, verify);
  test2(nttwork, verify);
  test3(nttwork, verify);

  nttwork_clear(nttwork);
  mpz_clear(modulus);
  return 0;
}
