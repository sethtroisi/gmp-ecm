/* mpzspm.c - "mpz small prime moduli" - pick a set of small primes large
   enough to represent a mpzv

  Copyright 2005, 2008, 2010 Dave Newman, Jason Papadopoulos, Paul Zimmermann.

  The SP Library is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or (at your
  option) any later version.

  The SP Library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
  License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with the SP Library; see the file COPYING.LIB.  If not, write to
  the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
  MA 02110-1301, USA.
*/

#include <stdlib.h>
#include <math.h>
#include "sp.h"


/* initialize mpzspm->T such that with m[j] := mpzspm->spm[j]->sp
   T[0][0] = m[0], ..., T[0][n-1] = m[n-1]
   ...
   T[d-1][0] = m[0]*...*m[ceil(n/2)-1], T[d-1][1] = m[ceil(n/2)] * ... * m[n-1]
   T[d][0] = m[0] * ... * m[n-1]
   where d = ceil(log(n)/log(2)).
   If n = 5, T[0]: 1, 1, 1, 1, 1
             T[1]: 2, 2, 1
             T[2]: 4, 1
*/
static void
mpzspm_product_tree_init (mpzspm_t mpzspm)
{
  uint32_t d, i, j, oldn;
  uint32_t n = mpzspm->sp_num;
  mpzv_t *T;

  for (i = n, d = 0; i > 1; i = (i + 1) / 2, d ++);
  if (d <= I0_THRESHOLD)
    {
      mpzspm->T = NULL;
      return;
    }
  T = (mpzv_t*) malloc ((d + 1) * sizeof (mpzv_t));
  T[0] = (mpzv_t) malloc (n * sizeof (mpz_t));
  for (j = 0; j < n; j++)
    {
      mpz_init(T[0][j]);
      mpz_set_sp(T[0][j], mpzspm->spm[j]->sp);
    }

  for (i = 1; i <= d; i++)
    {
      oldn = n;
      n = (n + 1) / 2;
      T[i] = (mpzv_t) malloc (n * sizeof (mpz_t));
      for (j = 0; j < n; j++)
        {
          mpz_init (T[i][j]);
          if (2 * j + 1 < oldn)
            mpz_mul (T[i][j], T[i-1][2*j], T[i-1][2*j+1]);
          else /* oldn is odd */
            mpz_set (T[i][j], T[i-1][2*j]);
        }
    }
  mpzspm->T = T;
  mpzspm->d = d;
}

/* This function initializes a mpzspm_t structure which contains the number
   of small primes, the small primes with associated primitive roots and 
   precomputed data for the CRT to allow convolution products of length up 
   to "max_len" with modulus "modulus". 
   Returns NULL in case of an error. */

mpzspm_t
mpzspm_init (sp_t max_len, mpz_t modulus)
{
  uint32_t ub, i, j;
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
  mpz_init (S);
  /* T is len*modulus^2, the upper bound on output coefficients of a 
     convolution */
  mpz_init (mt);
  mpz_init (T); 
  mpz_mul (T, modulus, modulus);
  mpz_set_uint64 (mp, max_len);
  mpz_mul (T, T, mt);
  
  /* find primes congruent to 1 mod max_len so we can do
   * a ntt of size max_len */
  /* Find the largest p <= SP_MAX that is p == 1 (mod max_len) */
  p = (SP_MAX / max_len) * max_len;
  if (p == SP_MAX) /* If max_len | SP_MAX, the +1 might cause overflow */
    p = p - max_len + 1;
  else
    p++;
  
  do
    {
      while (p >= SP_MIN && p > max_len && !sp_prime(p))
        p -= max_len;

      /* all primes must be in range */
      if (p < SP_MIN || p <= max_len)
        {
	  goto clear_mpzspm;
	}
      
      mpzspm->spm[mpzspm->sp_num] = spm_init (max_len, p);
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
  
  mpzspm->max_ntt_size = max_len;
  
  mpzspm->crt1 = (mpzv_t) malloc (mpzspm->sp_num * sizeof (mpz_t));
  mpzspm->crt2 = (mpzv_t) malloc ((mpzspm->sp_num + 2) * sizeof (mpz_t));
  mpzspm->crt3 = (spv_t) malloc (mpzspm->sp_num * sizeof (sp_t));
  mpzspm->crt4 = (spv_t *) calloc (mpzspm->sp_num, sizeof (spv_t));
  mpzspm->crt5 = (spv_t) malloc (mpzspm->sp_num * sizeof (sp_t));
  if (mpzspm->crt1 == NULL || mpzspm->crt2 == NULL || mpzspm->crt3 == NULL ||
      mpzspm->crt4 == NULL || mpzspm->crt5 == NULL)
    {
      goto clear_mpzspm;
    }

  for (i = 0; i < mpzspm->sp_num; i++)
    {
      mpzspm->crt4[i] = (spv_t) malloc (mpzspm->sp_num * sizeof (sp_t));
      if (mpzspm->crt4[i] == NULL)
        goto clear_mpzspm;
    }
  
  for (i = 0; i < mpzspm->sp_num; i++)
    {
      p = mpzspm->spm[i]->sp;
      mpz_set_sp (mp, p);
      
      /* crt3[i] = (P / p)^{-1} mod p */
      mpz_fdiv_q (T, P, mp);
      mpz_fdiv_r (mt, T, mp);
      a = mpz_get_sp (mt);
      mpzspm->crt3[i] = sp_inv (a, p, mpzspm->spm[i]->mul_c);
     
      /* crt1[i] = (P / p) mod modulus */
      mpz_init (mpzspm->crt1[i]);
      mpz_mod (mpzspm->crt1[i], T, modulus);

      /* crt4[i][j] = ((P / p[i]) mod modulus) mod p[j] */
      for (j = 0; j < mpzspm->sp_num; j++)
        {
          mpz_set_sp (mp, mpzspm->spm[j]->sp);
          mpz_fdiv_r (mt, mpzspm->crt1[i], mp);
          mpzspm->crt4[j][i] = mpz_get_sp (mt);
        }
      
      /* crt5[i] = (-P mod modulus) mod p */
      mpz_mod (T, P, modulus);
      mpz_sub (T, modulus, T);
      mpz_set_sp (mp, p);
      mpz_fdiv_r (mt, T, mp);
      mpzspm->crt5[i] = mpz_get_sp (mt);
    }
  
  mpz_set_ui (T, 0);
  
  for (i = 0; i < mpzspm->sp_num + 2; i++)
    {
      mpz_mod (T, T, modulus);
      mpz_init_set (mpzspm->crt2[i], T);
      mpz_sub (T, T, P);
    }
  
  mpz_clear (mp);
  mpz_clear (mt);
  mpz_clear (P);
  mpz_clear (S);
  mpz_clear (T);

  mpzspm_product_tree_init (mpzspm);
  return mpzspm;

clear_mpzspm:
  mpzspm_clear (mpzspm);
  return NULL;
}

/* clear the product tree T */
static void
mpzspm_product_tree_clear (mpzspm_t mpzspm)
{
  uint32_t i, j;
  uint32_t n = mpzspm->sp_num;
  uint32_t d = mpzspm->d;
  mpzv_t *T;

  if (mpzspm == NULL)
    return;

  T = mpzspm->T;
  if (T == NULL) /* use the slow method */
    return;

  for (i = 0; i <= d; i++)
    {
      for (j = 0; j < n; j++)
        mpz_clear (T[i][j]);
      free (T[i]);
      n = (n + 1) / 2;
    }
  free (T);
}

void mpzspm_clear (mpzspm_t mpzspm)
{
  unsigned int i;

  if (mpzspm == NULL)
    return;

  mpzspm_product_tree_clear (mpzspm);

  for (i = 0; i < mpzspm->sp_num; i++)
    {
      if (mpzspm->crt1)
	mpz_clear (mpzspm->crt1[i]);
      if (mpzspm->crt4[i])
	free (mpzspm->crt4[i]);
      if (mpzspm->spm[i])
	spm_clear (mpzspm->spm[i]);
    }

  if (mpzspm->crt2)
    {
      for (i = 0; i < mpzspm->sp_num + 2; i++)
	mpz_clear (mpzspm->crt2[i]);
    }
  
  free (mpzspm->crt1);
  free (mpzspm->crt2);
  free (mpzspm->crt3);
  free (mpzspm->crt4);
  free (mpzspm->crt5);
  
  mpz_clear (mpzspm->modulus);
  free (mpzspm->spm);
  free (mpzspm);
}
