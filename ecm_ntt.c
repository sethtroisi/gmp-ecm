/* ecm_ntt.c - high level poly functions to interface between ecm and sp

  Copyright 2005 Dave Newman.

  This file is part of the ECM Library.
  
  The ECM Library is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or (at your
  option) any later version.

  The ECM Library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
  License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with the ECM Library; see the file COPYING.LIB.  If not, write to
  the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
  MA 02111-1307, USA.
*/

#include <stdio.h>
#include "sp.h"
#include "ecm.h"
#include "ecm-impl.h"

/* external interface to ecm - so we have ntt_mul instead of mpzp_mul */

void
ntt_mul (mpzp_t r, mpzp_t x, mpzp_t y, spv_size_t len, mpzp_t t,
    int monic, mpz_t modulus)
{
  if (len < MUL_NTT_THRESHOLD)
  {
    list_mul (r, x, len, monic, y, len, monic, t);
    return;
  }

  mpzspm_t mpzspm = mpzspm_init (modulus, 1 << ceil_log_2 (len * 2 - 1 + monic));
  mpzspp_t X, Y;
    
  X = mpzspp_init (mpzspm);
  Y = mpzspp_init (mpzspm);
  
  mpzspp_set_mpzp (X, x, len, 0);
  mpzspp_set_mpzp (Y, y, len, 0);
  mpzspp_mul (X, X, Y, monic);
  mpzspp_get_mpzp (X, r, 2 * len - 1 + monic, 0);
  
  mpzspp_clear (X);
  mpzspp_clear (Y);
  mpzspm_clear (mpzspm);
}
 
void
ntt_mul_partial (mpzp_t r, mpzp_t x, spv_size_t x_len, mpzp_t y,
    spv_size_t y_len, spv_size_t k, spv_size_t l, mpzp_t dummy, int monic,
    mpz_t modulus)
{
  mpzspm_t mpzspm = mpzspm_init (modulus, 1 << ceil_log_2 (x_len + y_len - 1 + monic));
  mpzspp_t X, Y;
  int i;
    
  X = mpzspp_init (mpzspm);
  Y = mpzspp_init (mpzspm);
  
  mpzspp_set_mpzp (X, x, x_len, 0);
  mpzspp_set_mpzp (Y, y, y_len, 0);
  mpzspp_mul_partial (X, X, Y, k, l, monic);
  mpzspp_get_mpzp (X, r, l, 0);
  
  mpzspp_clear (X);
  mpzspp_clear (Y);
  mpzspm_clear (mpzspm);
}
  
void
ntt_PolyFromRoots (mpzp_t r, mpzp_t a, spv_size_t len, mpzp_t t, mpz_t n)
{
  ASSERT (len == 1 << ceil_log_2 (len));

  if (len <= MUL_NTT_THRESHOLD)
  {
    PolyFromRoots (r, a, len, t, n);
    return;
  }
  
  mpzspm_t mpzspm;
  mpzspp_t mpzspp;
  spv_t spv;
  spm_t spm;
  sp_t root, inv_root, p, d;
  spv_size_t i, m;
  unsigned int j;
   
  mpzspm = mpzspm_init (n, len);
  mpzspp = mpzspp_init (mpzspm);
  mpzspp_realloc (mpzspp, 2 * len);

  /* FIXME: Are we doing too much memory thrashing in this function? */
  
  for (i = 0; i < len; i += MUL_NTT_THRESHOLD)
    {
      PolyFromRoots (r, a + i, MUL_NTT_THRESHOLD, t, n);
      mpzspp_set_mpzp (mpzspp, r, MUL_NTT_THRESHOLD, 2 * i);
    }
  
  for (m = MUL_NTT_THRESHOLD; m < len; m *= 2)
    {
      for (j = 0; j < mpzspp->mpzspm->sp_num; j++)
        {
          spv = mpzspp->spv[j];
	  spm = mpzspp->mpzspm->spm + j;
	  p = spm->sp;
	  d = spm->mul_c;
	  
	  for (i = 0; i < 2 * len; i += 2 * m)
            {
              /* make monic and zero-pad */
	      spv[m + i] = 1;
	      spv_set_zero (spv + m + i + 1, m - 1);
            }
	  
	  root = sp_pow (spm->prim_root, (p - 1) / (2 * m), p, d); 
	  inv_root = sp_inv (root, p, d);
	  
	  for (i = 0; i < 2 * len; i += 4 * m)
	    {
	      spv_ntt_gfp_dif (spv + i, 2 * m, p, d, root);
	      spv_ntt_gfp_dif (spv + 2 * m + i, 2 * m, p, d, root);
	      spv_pwmul (spv + i, spv + i, spv + 2 * m + i, 2 * m, p, d);
	      spv_ntt_gfp_dit (spv + i, 2 * m, p, d, inv_root);  
	    }

	  root = sp_inv (2 * m, p, d);
	  for (i = 0; i < 2 * len; i += 4 * m)
	    {
	      spv_mul_sp (spv + i, spv + i, root, 2 * m, p, d);
	      
	      /* remove the overflowed x^(2*m) term */
	      spv[i] = sp_sub (spv[i], 1, p);
	    }
	}

      /* theoretically inefficient but seems to be faster
       * (probably due to less mallocing) */
      if (2 * m < len)
        {
	  mpzspp->normalised = 0;
	  mpzspp_normalise (mpzspp, 2 * len - 2 * m, 0);
	}
    }
      
  mpzspp_get_mpzp (mpzspp, r, len, 0);

  mpzspp_clear (mpzspp);
  mpzspm_clear (mpzspm);
}

  
int
ntt_PolyFromRoots_Tree (mpzp_t r, mpzp_t a, spv_size_t len, mpzp_t t, mpz_t n,
    mpzp_t *Tree, FILE *TreeFile)
{
  ASSERT (len == 1 << ceil_log_2 (len));

  if (len <= MUL_NTT_THRESHOLD)
  {
    PolyFromRoots_Tree (r, a, len, t, -1, n, Tree, TreeFile, 0);
    return;
  }
  
  mpzspm_t mpzspm;
  mpzspp_t mpzspp;
  spv_t spv;
  spm_t spm;
  sp_t root, inv_root, p, d;
  spv_size_t i, m;
  unsigned int j;
  unsigned int log_2_len = ceil_log_2 (len);
  mpzp_t src;
  mpzp_t *dst = Tree + log_2_len - 1;
  mpzspm = mpzspm_init (n, len);
  mpzspp = mpzspp_init (mpzspm);
  mpzspp_realloc (mpzspp, 2 * len);

  /* FIXME: Are we doing too much memory thrashing in this function? */
  
  list_set (*dst, a, len);
  
  if (TreeFile && list_out_raw (TreeFile, *dst, len) == ECM_ERROR)
    {
      outputf (OUTPUT_ERROR, "Error writing product tree of F\n");
      return ECM_ERROR;
    }
  
  for (m = 1; m < MUL_NTT_THRESHOLD; m *= 2)
    {
      src = *dst--;
      
      for (i = 0; i < len; i += 2 * m)
	list_mul (*dst + i, src + i, m, 1, src + i + m, m, 1, t);
      list_mod (*dst, *dst, len, n);

      if (TreeFile && list_out_raw (TreeFile, *dst, len) == ECM_ERROR)
        {
          outputf (OUTPUT_ERROR, "Error writing product tree of F\n");
          return ECM_ERROR;
        }
    }
  
  for (m = MUL_NTT_THRESHOLD; m < len; m *= 2)
    {
      src = *dst--;

      if (m == len / 2)
	dst = &r;
      
      for (i = 0; i < len; i += m)
	mpzspp_set_mpzp (mpzspp, src + i, m, 2 * i);
      
      for (j = 0; j < mpzspp->mpzspm->sp_num; j++)
        {
          spv = mpzspp->spv[j];
	  spm = mpzspp->mpzspm->spm + j;
	  p = spm->sp;
	  d = spm->mul_c;
	  
	  for (i = 0; i < 2 * len; i += 2 * m)
            {
              /* make monic and zero-pad */
	      spv[m + i] = 1;
	      spv_set_zero (spv + m + i + 1, m - 1);
            }
	  
	  root = sp_pow (spm->prim_root, (p - 1) / (2 * m), p, d); 
	  inv_root = sp_inv (root, p, d);
	  
	  for (i = 0; i < 2 * len; i += 4 * m)
	    {
	      spv_ntt_gfp_dif (spv + i, 2 * m, p, d, root);
	      spv_ntt_gfp_dif (spv + 2 * m + i, 2 * m, p, d, root);
	      spv_pwmul (spv + i, spv + i, spv + 2 * m + i, 2 * m, p, d);
	      spv_ntt_gfp_dit (spv + i, 2 * m, p, d, inv_root);  
	    }

	  root = sp_inv (2 * m, p, d);
	  for (i = 0; i < 2 * len; i += 4 * m)
	    {
	      spv_mul_sp (spv + i, spv + i, root, 2 * m, p, d);
	      
	      /* remove the overflowed x^(2*m) term */
	      spv[i] = sp_sub (spv[i], 1, p);
	    }
	}

      for (i = 0; i < len; i += 2 * m)
        mpzspp_get_mpzp (mpzspp, *dst + i, 2 * m, 2 * i);
      list_mod (*dst, *dst, len, n);

      if (TreeFile && list_out_raw (TreeFile, *dst, len) == ECM_ERROR)
        {
          outputf (OUTPUT_ERROR, "Error writing product tree of F\n");
          return ECM_ERROR;
        }
    }

  mpzspp_clear (mpzspp);
  mpzspm_clear (mpzspm);

  return 0;
}

void
ntt_PrerevertDivision (mpzp_t a, mpzp_t b, mpzp_t invb, spv_size_t len,
    mpzp_t t, mpz_t n)
{
  mpz_set_ui (a[2 * len - 1], 0);
  mpz_set_ui (invb[len - 1], 0);
  
  if (len < PREREVERT_DIVISION_NTT_THRESHOLD)
    {
      PrerevertDivision (a, b, invb, len, t, n);
      return;
    }
  
  ntt_mul (t, a + len, invb, len, t + 2 * len, 0, n);
  list_mod (t + len - 2, t + len - 2, len, n);
  mpz_set_ui (a[2 * len - 1], 0);

  /* Although our ntt size is only len, each term in the result is the sum
   * of 2 * len coeffs, so we need to make sure we don't overflow */
  mpzspm_t mpzspm = mpzspm_init (n, 2 * len);
  mpzspp_t x, y;
    
  x = mpzspp_init (mpzspm);
  y = mpzspp_init (mpzspm);
  
  mpzspp_set_mpzp (x, t + len - 2, len, 0);
  mpzspp_set_mpzp (y, b, len, 0);
  
  mpzspp_to_ntt (x, len, 0);
  mpzspp_to_ntt (y, len, 1); /* b has leading monomial */
  mpzspp_pwmul (x, x, y);
  mpzspp_get_mpzp (x, t, len, 0);
  
  mpzspp_clear (x);
  mpzspp_clear (y);
  mpzspm_clear (mpzspm);
 
  list_sub (t, t, a + len, len);
  list_sub (a, a, t, len);
  list_mod (a, a, len, n);
}

void ntt_PolyInvert (mpzp_t q, mpzp_t b, spv_size_t len, mpzp_t t, mpz_t n)
{
  if (len < POLYINVERT_NTT_THRESHOLD)
    {
      PolyInvert (q, b, len, t, n);
      return;
    }

  spv_size_t k = POLYINVERT_NTT_THRESHOLD / 2;
  
  PolyInvert (q + len - k, b + len - k, k, t, n);
  
  for (; k < len; k *= 2)
    {
      ntt_mul_partial (t, q + len - k, k, b + len - 2 * k, 2 * k, k - 1,
        2 * k - 1, NULL, 0, n);    
      list_neg (t, t + k - 1, k, n);
      list_mod (t, t, k, n);
  
      ntt_mul (t + k, t, q + len - k, k, t + 3 * k - 1, 0, n);
      list_mod (q + len - 2 * k, t + 2 * k - 1, k, n);
  }

#ifdef DEBUG
  ntt_mul_partial (t, q, len, b, len, 0, 2 * len - 1, NULL, 0, n);
  list_mod (t, t, 2 * len - 1, n);
  spv_size_t i;
  for (i = len - 1; i < 2 * len - 2; i++)
    if (mpz_cmp_ui (t[i], 0))
      printf ("error in ntt_PolyInvert\n");
  if (mpz_cmp_ui (t[2 * len - 2], 1))
    printf ("error in ntt_PolyInvert-\n");
#endif
}
