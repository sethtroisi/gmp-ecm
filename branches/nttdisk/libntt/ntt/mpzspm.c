#include <math.h>
#include "sp.h"

/* This function initializes a mpzspm_t structure which contains the number
   of small primes, the small primes with associated primitive roots and 
   precomputed data for the CRT to allow convolution products of length up 
   to "max_len" with modulus "modulus". 
   Returns NULL in case of an error. */

void *
X(mpzspm_init)(uint32_t max_len_in, mpz_t modulus)
{
  uint32_t ub, i, j;
  uint32_t max_len = max_len_in;
  sp_t a, p;
  mpz_t P, S, T, mp, mt;
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
  
  ub = (2 + 2 * mpz_sizeinbase (modulus, 2) + 
	ceil(log((double)(max_len)) / M_LN2) +
      4 * SP_NUMB_BITS) / (SP_NUMB_BITS - 1);
  
  mpzspm->spm = (spm_t *) malloc (ub * sizeof (spm_t));
  if (mpzspm->spm == NULL)
    goto clear_mpzspm;

  /* product of primes selected so far */
  mpz_init (mp);
  mpz_init_set_ui (P, 1UL);
  /* sum of primes selected so far */
  mpz_init_set_ui (S, 0UL);
  /* T is len*modulus^2, the upper bound on output coefficients of a 
     convolution */
  mpz_init (mt);
  mpz_init (T); 
  mpz_mul (T, modulus, modulus);
  mpz_set_ui (mt, max_len);
  mpz_mul (T, T, mt);
  
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
        {
	  goto clear_mpzspm;
	}
      
      mpzspm->spm[mpzspm->sp_num] = X(spm_init)(max_len, p);
      if (mpzspm->spm[mpzspm->sp_num] == NULL)
        {
          goto clear_mpzspm;
        }
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

  mpz_init_set (mpzspm->modulus, modulus);
  
  mpz_clear (mp);
  mpz_clear (mt);
  mpz_clear (P);
  mpz_clear (S);
  mpz_clear (T);

  return mpzspm;

clear_mpzspm:
  X(mpzspm_clear)(mpzspm);
  return NULL;
}

void X(mpzspm_clear)(void * m)
{
  mpzspm_t mpzspm = (mpzspm_t)m;
  unsigned int i;

  if (mpzspm == NULL)
    return;

  for (i = 0; i < mpzspm->sp_num; i++)
    {
	X(spm_clear)(mpzspm->spm[i]);
    }

  mpz_clear (mpzspm->modulus);
  free (mpzspm->spm);
  free (mpzspm);
}
