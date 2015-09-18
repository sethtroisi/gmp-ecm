#include "ntt-impl.h"

static void 
mpzspm_clear(void * m)
{
  mpzspm_t mpzspm = (mpzspm_t)m;
  uint32_t i;

  if (mpzspm == NULL)
    return;

  for (i = 0; i < mpzspm->sp_num; i++)
    {
      X(spm_clear)(mpzspm->spm[i]);
    }

  free (mpzspm->spm);
  free (mpzspm);
}

static void *
mpzspm_init(uint32_t max_len, mpz_t modulus,
        	mpz_t P, mpz_t S, uint32_t *done)
{
  uint32_t i;
  uint32_t max_spm;
  sp_t a, p;
  mpz_t T, mp, mt;
  mpzspm_t mpzspm;

  mpzspm = (mpzspm_t) calloc (1, sizeof (__mpzspm_struct));
  if (mpzspm == NULL)
    return NULL;
  
  /* Upper bound for the number of primes we need.
   * Let minp, maxp denote the min, max permissible prime,
   * S the sum of p_1, p_2, ..., p_ub,
   * P the product of p_1, p_2, ..., p_ub/
   * 
   * Choose ub s.t.
   *
   *     ub * log(minp) >= log(4 * max_len * modulus^2 * maxp^4)
   * 
   * =>  P >= minp ^ ub >= 4 * max_len * modulus^2 * maxp^4
   *                    >= 4 * max_len * modulus^2 * (ub * maxp)^2
   *                    >= 4 * max_len * modulus^2 * S^2
   * 
   * So we need at most ub primes to satisfy this condition. */
  
  mpz_init (mp);
  mpz_init (mt);
  mpz_init (T); 

  max_spm = 10;
  mpzspm->sp_num = 0;
  mpzspm->spm = (spm_t *) malloc (max_spm * sizeof (spm_t));
  if (mpzspm->spm == NULL)
    goto cleanup;

  /* find primes congruent to 1 mod max_len so we can do
   * a ntt of size max_len */
  /* Find the largest p <= SP_MAX that is p == 1 (mod max_len) */
  p = ((uint64_t)SP_MAX / max_len) * max_len;
  if (p == SP_MAX) /* If max_len | SP_MAX, the +1 might cause overflow */
    p = p - max_len + 1;
  else
    p++;
  
  do
    {
      while (p >= SP_MIN && p > max_len && !X(sp_prime)(p))
        p -= max_len;

      /* all primes must be in range */
      if (p < SP_MIN || p <= max_len)
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

      mpzspm->spm[mpzspm->sp_num] = X(spm_init)(max_len, p);
      if (mpzspm->spm[mpzspm->sp_num] == NULL)
    	break;
      mpzspm->sp_num++;
      
      mpz_set_sp (mp, p);
      mpz_mul (P, P, mp);
      mpz_add (S, S, mp);

      /* we want P > 4 * max_len * (modulus * S)^2. The S^2 term is due to 
         theorem 3.1 in Bernstein and Sorenson's paper */
      mpz_mul (T, S, modulus);
      mpz_mul (T, T, T);
      mpz_mul (T, T, mt);
      mpz_mul_2exp (T, T, 2UL);
      
      p -= max_len;
    }
  while (mpz_cmp (P, T) <= 0);

cleanup:
  mpz_clear (mp);
  mpz_clear (mt);
  mpz_clear (T);
  if (mpzspm->sp_num > 0)
    {
      if (mpz_cmp(P, T) > 0)
	*done = 1;
      return mpzspm;
    }

  mpzspm_clear(mpzspm);
  return NULL;
}

const __nttinit_struct
X(nttinit) = {
  SP_NUMB_BITS,
  mpzspm_init,
  mpzspm_clear
};
