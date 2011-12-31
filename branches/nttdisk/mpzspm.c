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

#include <stdio.h> /* for printf */
#include <stdlib.h>
#include "ecm-impl.h"


/* Tables for the maximum possible modulus (in bit size) for different 
   transform lengths l.
   The modulus is limited by the condition that primes must be 
   p_i == 1 (mod l), and \Prod_i p_i >= 4l (modulus * S)^2, 
   where S=\Sum_i p_i.
   Hence for each l=2^k, we take the product P and sum S of primes p_i,
   SP_MIN <= p_i <= SP_MAX and p_i == 1 (mod l), and store 
   floor (log_2 (sqrt (P / (4l S^2)))) in the table.
   We only consider power-of-two transform lengths <= 2^31 here.

   Table entries generated with
   
   l=2^k;p=1;P=1;S=0;while(p<=SP_MAX, if(p>=SP_MIN && isprime(p), S+=p; P*=p); \
   p+=l);print(floor (log2 (sqrt (P / (4*l * S^2)))))

   in Pari/GP for k=9 ... 24. k<9 simply were doubled and rounded down in 
   each step.

   We curently assume that SP_MIN == 2^(SP_NUMB_BITS-1) and 
   SP_MAX == 2^(SP_NUMB_BITS).
   
*/

#if (SP_NUMB_BITS == 30)
static unsigned long sp_max_modulus_bits[32] = 
  {0, 380000000, 190000000, 95000000, 48000000, 24000000, 12000000, 6000000, 
   3000000, 1512786, 756186, 378624, 188661, 93737, 46252, 23342, 11537, 5791, 
   3070, 1563, 782, 397, 132, 43, 0, 0, 0, 0, 0, 0, 0, 0};
#elif (SP_NUMB_BITS == 31)
static unsigned long sp_max_modulus_bits[32] = 
  {0, 750000000, 380000000, 190000000, 95000000, 48000000, 24000000, 12000000, 
   6000000, 3028766, 1512573, 756200, 379353, 190044, 94870, 47414, 23322, 
   11620, 5891, 2910, 1340, 578, 228, 106, 60, 30, 0, 0, 0, 0, 0, 0};
#elif (SP_NUMB_BITS == 32)
static unsigned long sp_max_modulus_bits[32] = 
  {0, 1520000000, 760000000, 380000000, 190000000, 95000000, 48000000, 
   24000000, 12000000, 6041939, 3022090, 1509176, 752516, 376924, 190107, 
   95348, 47601, 24253, 11971, 6162, 3087, 1557, 833, 345, 172, 78, 46, 15, 
   0, 0, 0, 0};
#elif (SP_NUMB_BITS >= 40)
  /* There are so many primes, we can do pretty much any modulus with 
     any transform length. I didn't bother computing the actual values. */
static unsigned long sp_max_modulus_bits[32] =  
  {0, ULONG_MAX, ULONG_MAX, ULONG_MAX, ULONG_MAX, ULONG_MAX, ULONG_MAX, 
   ULONG_MAX, ULONG_MAX, ULONG_MAX, ULONG_MAX, ULONG_MAX, ULONG_MAX, ULONG_MAX, 
   ULONG_MAX, ULONG_MAX, ULONG_MAX, ULONG_MAX, ULONG_MAX, ULONG_MAX, ULONG_MAX, 
   ULONG_MAX, ULONG_MAX, ULONG_MAX, ULONG_MAX, ULONG_MAX, ULONG_MAX, ULONG_MAX, 
   ULONG_MAX, ULONG_MAX, ULONG_MAX, ULONG_MAX};
#else
#error Table of maximal modulus for transform lengths not defined for this SP_MIN
;
#endif


/* Returns the largest possible transform length we can do for modulus
   without running out of primes */

spv_size_t
mpzspm_max_len (mpz_t modulus)
{
  int i;
  size_t b;

  b = mpz_sizeinbase (modulus, 2); /* b = floor (log_2 (modulus)) + 1 */
  /* Transform length 2^k is ok if log2(modulus) <= sp_max_modulus_bits[k]
     <==> ceil(log2(modulus)) <= sp_max_modulus_bits[k] 
     <==> floor(log_2(modulus)) + 1 <= sp_max_modulus_bits[k] if modulus 
     isn't a power of 2 */
     
  for (i = 0; i < 30; i++)
    {
      if (b > sp_max_modulus_bits[i + 1])
	break;
    }

  return (spv_size_t)1 << i;
}

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
  unsigned int d, i, j, oldn;
  unsigned int n = mpzspm->sp_num;
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
    mpz_init_set_ui (T[0][j], (mp_limb_t) mpzspm->spm[j]->sp);
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
mpzspm_init (spv_size_t max_len, mpz_t modulus)
{
  unsigned int ub, i, j;
  mpz_t P, S, T, mp, mt; /* mp is p as mpz_t, mt is a temp mpz_t */
  sp_t p, a;
  mpzspm_t mpzspm;
  long st;

  st = cputime ();

  mpzspm = (mpzspm_t) malloc (sizeof (__mpzspm_struct));
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
  
  ub = (2 + 2 * mpz_sizeinbase (modulus, 2) + ceil_log_2 (max_len) + \
      4 * SP_NUMB_BITS) / (SP_NUMB_BITS - 1);
  
  mpzspm->spm = (spm_t *) malloc (ub * sizeof (spm_t));
  if (mpzspm->spm == NULL)
    goto error_clear_mpzspm;
  mpzspm->sp_num = 0;

  /* product of primes selected so far */
  mpz_init_set_ui (P, 1UL);
  /* sum of primes selected so far */
  mpz_init (S);
  /* T is len*modulus^2, the upper bound on output coefficients of a 
     convolution */
  mpz_init (T); 
  mpz_mul (T, modulus, modulus);
  mpz_mul_ui (T, T, max_len);
  mpz_init (mp);
  mpz_init (mt);
  
  /* find primes congruent to 1 mod max_len so we can do
   * a ntt of size max_len */
  /* Find the largest p <= SP_MAX that is p == 1 (mod max_len) */
  p = (SP_MAX / (sp_t) max_len) * (sp_t) max_len;
  if (p == SP_MAX) /* If max_len | SP_MAX, the +1 might cause overflow */
    p = p - (sp_t) max_len + (sp_t) 1;
  else
    p++;
  
  do
    {
      while (p >= SP_MIN && p > (sp_t) max_len && !sp_prime(p))
        p -= (sp_t) max_len;

      /* all primes must be in range */
      if (p < SP_MIN || p <= (sp_t) max_len)
        {
	  outputf (OUTPUT_ERROR, 
	           "not enough primes == 1 (mod %lu) in interval\n", 
	           (unsigned long) max_len);
	  goto error_clear_mpzspm_spm;
	}
      
      mpzspm->spm[mpzspm->sp_num] = spm_init (max_len, p);
      if (mpzspm->spm[mpzspm->sp_num] == NULL)
        {
          outputf (OUTPUT_ERROR, "Out of memory in mpzspm_init()\n");
          goto error_clear_mpzspm_spm;
        }
      mpzspm->sp_num++;
      
      mpz_set_sp (mp, p);
      mpz_mul (P, P, mp);
      mpz_add (S, S, mp);

      /* we want P > 4 * max_len * (modulus * S)^2. The S^2 term is due to 
         theorem 3.1 in Bernstein and Sorenson's paper */
      mpz_mul (T, S, modulus);
      mpz_mul (T, T, T);
      mpz_mul_ui (T, T, max_len);
      mpz_mul_2exp (T, T, 2UL);
      
      p -= (sp_t) max_len;
    }
  while (mpz_cmp (P, T) <= 0);

  outputf (OUTPUT_DEVVERBOSE, "mpzspm_init: finding %u primes took %lums\n", 
           mpzspm->sp_num, cputime() - st);

  mpz_init_set (mpzspm->modulus, modulus);
  
  mpzspm->max_ntt_size = max_len;
  
  mpzspm->crt1 = (mpzv_t) malloc (mpzspm->sp_num * sizeof (mpz_t));
  mpzspm->crt2 = (mpzv_t) malloc ((mpzspm->sp_num + 2) * sizeof (mpz_t));
  mpzspm->crt3 = (spv_t) malloc (mpzspm->sp_num * sizeof (sp_t));
  mpzspm->crt4 = (spv_t *) malloc (mpzspm->sp_num * sizeof (spv_t));
  mpzspm->crt5 = (spv_t) malloc (mpzspm->sp_num * sizeof (sp_t));
  if (mpzspm->crt1 == NULL || mpzspm->crt2 == NULL || mpzspm->crt3 == NULL ||
      mpzspm->crt4 == NULL || mpzspm->crt5 == NULL)
    {
      outputf (OUTPUT_ERROR, "Out of memory in mpzspm_init()\n");
      goto error_clear_crt;
    }

  for (i = 0; i < mpzspm->sp_num; i++)
    mpzspm->crt4[i] = NULL;
  for (i = 0; i < mpzspm->sp_num; i++)
    {
      mpzspm->crt4[i] = (spv_t) malloc (mpzspm->sp_num * sizeof (sp_t));
      if (mpzspm->crt4[i] == NULL)
        goto error_clear_crt4;
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

  outputf (OUTPUT_DEVVERBOSE, "mpzspm_init took %lums\n", cputime() - st);

  return mpzspm;
  
  /* Error cases: free memory we allocated so far */

  error_clear_crt4:
  for (i = 0; i < mpzspm->sp_num; i++)
    free (mpzspm->crt4[i]);
  
  error_clear_crt:
  free (mpzspm->crt1);
  free (mpzspm->crt2);
  free (mpzspm->crt3);
  free (mpzspm->crt4);
  free (mpzspm->crt5);
  
  error_clear_mpzspm_spm:
  for (i = 0; i < mpzspm->sp_num; i++)
    free(mpzspm->spm[i]);
  free (mpzspm->spm);

  error_clear_mpzspm:
  free (mpzspm);

  return NULL;
}

/* clear the product tree T */
static void
mpzspm_product_tree_clear (mpzspm_t mpzspm)
{
  unsigned int i, j;
  unsigned int n = mpzspm->sp_num;
  unsigned int d = mpzspm->d;
  mpzv_t *T = mpzspm->T;

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

  mpzspm_product_tree_clear (mpzspm);

  for (i = 0; i < mpzspm->sp_num; i++)
    {
      mpz_clear (mpzspm->crt1[i]);
      free (mpzspm->crt4[i]);
      spm_clear (mpzspm->spm[i]);
    }

  for (i = 0; i < mpzspm->sp_num + 2; i++)
    mpz_clear (mpzspm->crt2[i]);
  
  free (mpzspm->crt1);
  free (mpzspm->crt2);
  free (mpzspm->crt3);
  free (mpzspm->crt4);
  free (mpzspm->crt5);
  
  mpz_clear (mpzspm->modulus);
  free (mpzspm->spm);
  free (mpzspm);
}

