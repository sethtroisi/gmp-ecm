#include "ntt-impl.h"

static void 
mpzspm_ntt_reset(void * m)
{
  mpzspm_t mpzspm = (mpzspm_t)m;
  uint32_t i;

  if (mpzspm == NULL)
    return;

  X(ntt_reset)(&mpzspm->nttdata);
  X(ntt_reset)(&mpzspm->inttdata);
  sp_aligned_free (mpzspm->work);
  sp_aligned_free (mpzspm->sp);
  mpzspm->work = NULL;
  mpzspm->sp = NULL;
}

static void 
mpzspm_clear(void * m)
{
  mpzspm_t mpzspm = (mpzspm_t)m;
  uint32_t i;

  if (mpzspm == NULL)
    return;

  mpzspm_ntt_reset(m);

  for (i = 0; i < mpzspm->sp_num; i++)
    {
      X(spm_clear)(mpzspm->spm[i]);
    }

  free (mpzspm->spm);
  free (mpzspm);
}

static void *
mpzspm_init(uint32_t max_ntt_size, mpz_t modulus,
        	mpz_t P, mpz_t S, uint32_t interleaved,
            uint32_t *done)
{
  uint32_t i;
  uint32_t max_spm;
  sp_t a, p;
  mpz_t T, mp, mt;
  mpzspm_t mpzspm;
#if SP_NUMB_BITS == 50
  fpu_precision_t old_prec = 0;
#endif

  mpzspm = (mpzspm_t) calloc (1, sizeof (__mpzspm_struct));
  if (mpzspm == NULL)
    return NULL;
  
#if SP_NUMB_BITS == 50
  if (!fpu_precision_is_ieee())
    old_prec = fpu_set_precision_ieee();
#endif

  /* Upper bound for the number of primes we need.
   * Let minp, maxp denote the min, max permissible prime,
   * S the sum of p_1, p_2, ..., p_ub,
   * P the product of p_1, p_2, ..., p_ub/
   * 
   * Choose ub s.t.
   *
   *     ub * log(minp) >= log(4 * max_ntt_size * modulus^2 * maxp^4)
   * 
   * =>  P >= minp ^ ub >= 4 * max_ntt_size * modulus^2 * maxp^4
   *                    >= 4 * max_ntt_size * modulus^2 * (ub * maxp)^2
   *                    >= 4 * max_ntt_size * modulus^2 * S^2
   * 
   * So we need at most ub primes to satisfy this condition. */
  
  mpz_init_set_ui (mt, max_ntt_size);
  mpz_init (mp);
  mpz_init (T); 

  max_spm = 10;
  mpzspm->sp_num = 0;
  mpzspm->spm = (spm_t *) malloc (max_spm * sizeof (spm_t));
  if (mpzspm->spm == NULL)
    goto cleanup;

  /* find primes congruent to 1 mod max_ntt_size so we can do
   * a ntt of size max_ntt_size */
  /* Find the largest p <= SP_MAX that is p == 1 (mod max_ntt_size) */
  p = ((uint64_t)SP_MAX / max_ntt_size) * max_ntt_size;
  if (p == SP_MAX) /* If max_ntt_size | SP_MAX, the +1 might cause overflow */
    p = p - max_ntt_size + 1;
  else
    p++;
  
  do
    {
      while (p >= SP_MIN && p > max_ntt_size && !X(sp_prime)(p))
        p -= max_ntt_size;

      /* all primes must be in range */
      if (p < SP_MIN || p <= max_ntt_size)
	break;
      
      /* add this p */

      if (mpzspm->sp_num == max_spm)
	{
	  spm_t *tmp;

	  max_spm *= 2;
	  tmp = (spm_t *)realloc(mpzspm->spm, max_spm * sizeof(spm_t));
	  if (tmp == NULL)
	    break;

	  mpzspm->spm = tmp;
	}

      mpzspm->spm[mpzspm->sp_num] = X(spm_init)(max_ntt_size, p);
      if (mpzspm->spm[mpzspm->sp_num] == NULL)
    	break;
      mpzspm->sp_num++;
      
      mpz_set_sp (mp, p);
      mpz_mul (P, P, mp);
      mpz_add (S, S, mp);

      /* we want P > 4 * max_ntt_size * (modulus * S)^2. The S^2 term is due to 
         theorem 3.1 in Bernstein and Sorenson's paper */
      mpz_mul (T, S, modulus);
      mpz_mul (T, T, T);
      mpz_mul (T, T, mt);
      mpz_mul_2exp (T, T, 2UL);
      
      p -= max_ntt_size;
    }
  while (mpz_cmp (P, T) <= 0);

cleanup:
  if (mpzspm->sp_num > 0)
    {
      mpzspm->interleaved = interleaved;
      if (mpz_cmp(P, T) > 0)
	*done = 1;
    }

  mpz_clear (mp);
  mpz_clear (mt);
  mpz_clear (T);
#if SP_NUMB_BITS == 50
  if (old_prec != 0)
    fpu_clear_precision(old_prec);
#endif

  if (mpzspm->sp_num > 0)
    return mpzspm;
  mpzspm_clear(mpzspm);
  return NULL;
}

static void 
mpzspm_ntt_init(void * m, uint32_t ntt_size,
    		uint32_t max_ntt_size, nttplangroup_t *p)
{
  mpzspm_t mpzspm = (mpzspm_t)m;
#if SP_NUMB_BITS == 50
  fpu_precision_t old_prec = 0;
  if (!fpu_precision_is_ieee())
    old_prec = fpu_set_precision_ieee();
#endif

  if (mpzspm->interleaved)
    {
      const nttgroup_t ** groups = X(ntt_master_group_list);
      uint32_t i;
      size_t vsize = 0;

      for (i = 0; i < p->num_plans; i++)
	vsize = MAX(vsize, groups[p->plans[i].group_type]->vsize);

      /* maximum vector length (= stride of work array) */

      mpzspm->work = (spv_t)sp_aligned_malloc(vsize * 
         		 ((mpzspm->sp_num + vsize - 1) / vsize) * 
				ntt_size * sizeof(sp_t));
      /* aligned list of moduli */

      mpzspm->sp = (spv_t)sp_aligned_malloc(vsize * 
         		 ((mpzspm->sp_num + vsize - 1) / vsize) * 
				sizeof(sp_t));

      mpzspm->max_vsize = vsize;
      for (i = 0; i < mpzspm->sp_num; i++)
	mpzspm->sp[i] = mpzspm->spm[i]->sp;
    }
  else
    {
      mpzspm->max_vsize = 1;
      mpzspm->work = (spv_t)sp_aligned_malloc(
	  			mpzspm->sp_num * ntt_size * 
				sizeof(sp_t));
    }

  /* all the other NTT precomputations */

  X(ntt_build_passes) (p, mpzspm, ntt_size, max_ntt_size);

#if SP_NUMB_BITS == 50
  if (old_prec != 0)
    fpu_clear_precision(old_prec);
#endif
}

static void 
mpzspm_random(void * m, uint32_t ntt_size)
{
  mpzspm_t mpzspm = (mpzspm_t)m;
  uint32_t i, j, k;
#if SP_NUMB_BITS == 50
  fpu_precision_t old_prec = 0;
  if (!fpu_precision_is_ieee())
    old_prec = fpu_set_precision_ieee();
#endif

  if (mpzspm->interleaved)
    {
      spv_size_t vsize = mpzspm->max_vsize;
      spv_size_t batches = (mpzspm->sp_num + vsize - 1) / vsize;
      spv_t tmp = (spv_t)sp_aligned_malloc(ntt_size * sizeof(sp_t));

      for (i = 0; i < batches; i++)
	{
	  for (j = 0; j < vsize; j++)
	    {
	      if (i * vsize + j < mpzspm->sp_num)
		{
		  spm_t spm = mpzspm->spm[i * vsize + j];
		  spv_t work = mpzspm->work + i * vsize * ntt_size;

		  X(spv_random)(tmp, ntt_size, spm->sp);

		  for (k = 0; k < ntt_size; k++)
		    work[k * vsize + j] = tmp[k];
		}
	    }
	}

      sp_aligned_free(tmp);
    }
  else
    {
      for (i = 0; i < mpzspm->sp_num; i++)
	X(spv_random)(mpzspm->work + i * ntt_size, 
	      		    ntt_size, mpzspm->spm[i]->sp);
    }

#if SP_NUMB_BITS == 50
  if (old_prec != 0)
    fpu_clear_precision(old_prec);
#endif
}

static void
bfntt(spv_t r, spv_t x, spv_size_t len, 
     sp_t p, sp_t d, sp_t primroot, sp_t order)
{
  sp_t w0 = primroot;
  sp_t w_inc = 1;
  spv_size_t i, j, k;

#if SP_NUMB_BITS == 50
  fpu_precision_t old_prec = 0;
  if (!fpu_precision_is_ieee())
    old_prec = fpu_set_precision_ieee();
#endif

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

#if SP_NUMB_BITS == 50
  if (old_prec != 0)
    fpu_clear_precision(old_prec);
#endif
}

static void compare(spv_t r0, spv_t r1, sp_t p, uint32_t size)
{
  uint32_t j, k;

#ifdef HAVE_PARTIAL_MOD
  for (j = 0; j < size; j++)
    while (r1[j] >= p)
      r1[j] -= p;
#endif

  for (j = 0; j < size; j++)
    {
      for (k = 0; k < size; k++)
	{
	  if (r0[j] == r1[k])
	    break;
	}
      if (k == size)
	break;
    }
  printf("%c", j == size ? '.' : '*');
}

static void 
mpzspm_test(void * m, uint32_t ntt_size, uint32_t max_ntt_size)
{
  mpzspm_t mpzspm = (mpzspm_t)m;
  uint32_t i, j, k;
  spv_t r = (spv_t)malloc(mpzspm->sp_num * ntt_size * sizeof(sp_t));
  spv_t r0 = (spv_t)malloc(ntt_size * sizeof(sp_t));

  if (mpzspm->interleaved)
    {
      spv_size_t vsize = mpzspm->max_vsize;
      spv_size_t batches = (mpzspm->sp_num + vsize - 1) / vsize;

      for (i = 0; i < batches; i++)
	{
	  for (j = 0; j < vsize; j++)
	    {
	      if (i * vsize + j < mpzspm->sp_num)
		{
		  spm_t spm = mpzspm->spm[i * vsize + j];
		  sp_t p = spm->sp;
		  sp_t d = spm->mul_c;
		  sp_t primroot = spm->primroot;
		  spv_t work = mpzspm->work + i * vsize * ntt_size;

		  for (k = 0; k < ntt_size; k++)
		    r0[k] = work[k * vsize + j];

		  bfntt(r + (i * vsize + j) * ntt_size, r0,
			ntt_size, p, d, primroot, max_ntt_size);
		}
	    }
	}
    }
  else
    {
      for (i = 0; i < mpzspm->sp_num; i++)
	{
	  spm_t spm = mpzspm->spm[i];
	  sp_t p = spm->sp;
	  sp_t d = spm->mul_c;
	  sp_t primroot = spm->primroot;

    	  bfntt(r + i * ntt_size, mpzspm->work + i * ntt_size, 
	      	ntt_size, p, d, primroot, max_ntt_size);
	}
    }

  X(ntt_run)(m, NULL, ntt_size);

  if (mpzspm->interleaved)
    {
      spv_size_t vsize = mpzspm->max_vsize;
      spv_size_t batches = (mpzspm->sp_num + vsize - 1) / vsize;

      for (i = 0; i < batches; i++)
	{
	  for (j = 0; j < vsize; j++)
	    {
	      if (i * vsize + j < mpzspm->sp_num)
		{
		  spm_t spm = mpzspm->spm[i * vsize + j];
		  spv_t work = mpzspm->work + i * vsize * ntt_size;

		  for (k = 0; k < ntt_size; k++)
		    r0[k] = work[k * vsize + j];

		  compare(r + (i * vsize + j) * ntt_size, r0,
			spm->sp, ntt_size);
		}
	    }
	}
    }
  else
    {
      for (i = 0; i < mpzspm->sp_num; i++)
	{
	  compare(r + i * ntt_size, 
      		  mpzspm->work + i * ntt_size,
		  mpzspm->spm[i]->sp, ntt_size);
	}
    }
  printf(" ");

  free(r);
  free(r0);
}

const __nttinit_struct
X(nttinit) = {
  SP_NUMB_BITS,
  X(get_master_group_list_size),
  mpzspm_init,
  mpzspm_clear,
  mpzspm_ntt_init,
  X(ntt_run),
  mpzspm_ntt_reset,
  mpzspm_random,
  mpzspm_test
};
