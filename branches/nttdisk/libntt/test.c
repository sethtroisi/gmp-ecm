#include <stdio.h>
#include "ntt-impl.h"

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
  sp_t p = mpzspm->spm[0]->sp;
  sp_t d = mpzspm->spm[0]->mul_c;
  sp_t primroot = mpzspm->spm[0]->primroot;
  sp_t order = mpzspm->ntt_size;
  nttdata_t *nttdata = (nttdata_t *)mpzspm->spm[0]->ntt_data;

  codelet_data_t * pfa1 = nttdata->codelets + 1;
  codelet_data_t * pfa2 = nttdata->codelets + 2;
  const nttconfig_t *config1 = pfa1->config;
  const nttconfig_t *config2 = pfa2->config;
  spv_size_t len = config1->size * config2->size;
  spv_size_t i;

  sp_t tmp1[20];
  sp_t tmp2[20];
  sp_t x[20];
  sp_t r[20];

  spv_random(x, len, p);

  for (i = 0; i < len; i++)
    printf("in  %" PRIxsp "\n", x[i]);
  printf("\n");

  bfntt(r, x, len, p, d, primroot, order);

  for (i = 0; i < len; i++)
    printf("ref %" PRIxsp "\n", r[i]);
  printf("\n");

  config1->ntt_pfa_run(x, 1, config2->size, p, d, pfa1->ntt_const);
  config2->ntt_pfa_run(x, 1, config1->size, p, d, pfa2->ntt_const);

  for (i = 0; i < len; i++)
    printf("pfa %" PRIxsp "\n", x[i]);
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
