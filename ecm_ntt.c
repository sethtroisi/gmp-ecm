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
    int monic, mpzspm_t mpzspm)
{
  if (len < MUL_NTT_THRESHOLD)
  {
    list_mul (r, x, len, monic, y, len, monic, t);
    return;
  }

  mpzspp_t X, Y;
    
  X = mpzspp_init2 (mpzspm, 2 * len);
  Y = mpzspp_init2 (mpzspm, 2 * len);
  
  mpzspp_set_mpzp (X, x, len, 0);
  mpzspp_to_ntt (X, 2 * len, monic);
  mpzspp_set_mpzp (Y, y, len, 0);
  mpzspp_to_ntt (Y, 2 * len, monic);
  mpzspp_pwmul (X, X, Y);
  mpzspp_from_ntt (X);
  mpzspp_get_mpzp (X, r, 2 * len - 1 + monic, 0);
  
  mpzspp_clear (X);
  mpzspp_clear (Y);
}

#if 0
void
ntt_mul_partial (mpzp_t r, mpzp_t x, spv_size_t x_len, mpzp_t y,
    spv_size_t y_len, spv_size_t k, spv_size_t l, mpzp_t dummy, int monic,
    mpzspm_t mpzspm)
{
  mpzspp_t X, Y;
    
  X = mpzspp_init (mpzspm);
  Y = mpzspp_init (mpzspm);
  
  mpzspp_set_mpzp (X, x, x_len, 0);
  mpzspp_set_mpzp (Y, y, y_len, 0);
  mpzspp_mul_partial (X, X, Y, k, l, monic);
  mpzspp_get_mpzp (X, r, l, 0);
  
  mpzspp_clear (X);
  mpzspp_clear (Y);
}
#endif

void
ntt_PolyFromRoots (mpzp_t r, mpzp_t a, spv_size_t len, mpzp_t t,
    mpzspm_t mpzspm)
{
  ASSERT (len == 1 << ceil_log_2 (len));

  if (len <= MUL_NTT_THRESHOLD)
  {
    PolyFromRoots (r, a, len, t, mpzspm->modulus);
    return;
  }
  
  mpzspp_t mpzspp;
  spv_t spv;
  spm_t spm;
  sp_t root, inv_root, p, d;
  spv_size_t i, m;
  unsigned int j;
   
  mpzspp = mpzspp_init2 (mpzspm, 2 * len);

  /* FIXME: Are we doing too much memory thrashing in this function? */
  
  for (i = 0; i < len; i += MUL_NTT_THRESHOLD)
    {
      PolyFromRoots (r, a + i, MUL_NTT_THRESHOLD, t, mpzspm->modulus);
      mpzspp_set_mpzp (mpzspp, r, MUL_NTT_THRESHOLD, 2 * i);
    }
  
  for (m = MUL_NTT_THRESHOLD; m < len; m *= 2)
    {
      for (j = 0; j < mpzspm->sp_num; j++)
        {
          spv = mpzspp->spv[j];
	  spm = mpzspm->spm + j;
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
	mpzspp_normalise (mpzspp, 2 * len - 2 * m, 0);
    }
      
  mpzspp_get_mpzp (mpzspp, r, len, 0);

  mpzspp_clear (mpzspp);
}

  
int
ntt_PolyFromRoots_Tree (mpzp_t r, mpzp_t a, spv_size_t len, mpzp_t t,
    mpzspm_t mpzspm, mpzp_t *Tree, FILE *TreeFile)
{
  ASSERT (len == 1 << ceil_log_2 (len));

  if (len <= MUL_NTT_THRESHOLD)
    {
      return PolyFromRoots_Tree
        (r, a, len, t, -1, mpzspm->modulus, Tree, TreeFile, 0);
    }
  
  mpzspp_t mpzspp;
  spv_t spv;
  spm_t spm;
  sp_t root, inv_root, p, d;
  spv_size_t i, m;
  unsigned int j;
  unsigned int log_2_len = ceil_log_2 (len);
  mpzp_t src;
  mpzp_t *dst = Tree + log_2_len - 1;
  
  mpzspp = mpzspp_init2 (mpzspm, 2 * len);

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
      list_mod (*dst, *dst, len, mpzspm->modulus);

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
      
      for (j = 0; j < mpzspm->sp_num; j++)
        {
          spv = mpzspp->spv[j];
	  spm = mpzspm->spm + j;
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
      list_mod (*dst, *dst, len, mpzspm->modulus);

      if (TreeFile && list_out_raw (TreeFile, *dst, len) == ECM_ERROR)
        {
          outputf (OUTPUT_ERROR, "Error writing product tree of F\n");
          return ECM_ERROR;
        }
    }

  mpzspp_clear (mpzspp);

  return 0;
}


/* 2 NTTs of size 2 * len
 * 3 NTTs of size len */
void
ntt_PrerevertDivision (mpzp_t a, mpzp_t b, mpzp_t invb, mpzspp_t sp_invb,
    spv_size_t len, mpzp_t t, mpzspm_t mpzspm)
{
  if (len < PREREVERT_DIVISION_NTT_THRESHOLD)
    {
      PrerevertDivision (a, b, invb, len, t, mpzspm->modulus);
      return;
    }
  
  mpzspp_t x = mpzspp_init2 (mpzspm, 2 * len);
  mpzspp_t y = mpzspp_init2 (mpzspm, 2 * len);

  /* y = TOP (TOP (a) * invb) */
  mpzspp_set_mpzp (x, a + len, len - 1, 0);
  mpzspp_to_ntt (x, 2 * len, 0);
  mpzspp_pwmul (x, x, sp_invb);
  mpzspp_from_ntt (x);
  mpzspp_normalise (x, len, len - 1);
  
  /* can we avoid this (illegal) copy? */
  mpzspp_set (x, x, len, 0, len - 1);
  x->len = len;
  mpzspp_to_ntt (x, len, 0);
  mpzspp_set_mpzp (y, b, len, 0);
  y->len = len;
  mpzspp_to_ntt (y, len, 1); /* b has leading monomial */
  mpzspp_pwmul (x, x, y);
  mpzspp_from_ntt (x);
  mpzspp_get_mpzp (x, t, len, 0);
  
  mpzspp_clear (x);
  mpzspp_clear (y);
 
  list_sub (t, t, a + len, len - 1);
  list_sub (a, a, t, len);
  list_mod (a, a, len, mpzspm->modulus);
}

void ntt_PolyInvert (mpzp_t q, mpzp_t b, spv_size_t len, mpzp_t t,
    mpzspm_t mpzspm)
{
  if (len < POLYINVERT_NTT_THRESHOLD)
    {
      PolyInvert (q, b, len, t, mpzspm->modulus);
      return;
    }

  spv_size_t k = POLYINVERT_NTT_THRESHOLD / 2;
  
  PolyInvert (q + len - k, b + len - k, k, t, mpzspm->modulus);
  
  mpzspp_t w = mpzspp_init2 (mpzspm, len / 2);
  mpzspp_t x = mpzspp_init2 (mpzspm, len + 1);
  mpzspp_t y = mpzspp_init2 (mpzspm, len);
  mpzspp_t z = mpzspp_init2 (mpzspm, len);
  
  mpzspp_set_mpzp (x, q + len - k - 1, k + 1, 0);
  mpzspp_set_mpzp (y, b, len - 1, 0);
  
  for (; k < len; k *= 2)
    {
#if 0
      mpz_set_ui (q[len - k - 1], 0);
      ntt_mul_partial (t, q + len - k - 1, k + 1, b + len - 2 * k, 2 * k - 1, k,
        2 * k, NULL, 0, mpzspm);    
      list_neg (t, t + k, k, mpzspm->modulus);
      list_mod (t, t, k, mpzspm->modulus);
  
      ntt_mul (t + k, t, q + len - k, k, t + 3 * k - 1, 0, mpzspm);
      list_mod (q + len - 2 * k, t + 2 * k - 1, k, mpzspm->modulus);
#else
      mpzspp_set (w, x, k, 0, 1);
      mpzspp_to_ntt (x, 2 * k, 0);
      mpzspp_set (z, y, 2 * k - 1, 0, len - 2 * k);
      z->len = 2 * k - 1;
      mpzspp_to_ntt (z, 2 * k, 0);
      mpzspp_pwmul (z, x, z);
      mpzspp_from_ntt (z);
      mpzspp_normalise (z, k, k);
      mpzspp_neg (z, z, k, k);
      
      mpzspp_to_ntt (z, 2 * k, 0);
      mpzspp_pwmul (x, x, z);
      mpzspp_from_ntt (x);
      mpzspp_normalise (x, k, k);
      mpzspp_set (x, x, k, 1, k); /* illegal copy (overlap) */
      mpzspp_set (x, w, k, k + 1, 0);
#endif
    }

  mpzspp_get_mpzp (x, q, len - POLYINVERT_NTT_THRESHOLD / 2, 1);
 
#if defined DEBUG
  ntt_mul (t, q, b, len, NULL, 0, mpzspm);
  list_mod (t, t, 2 * len - 1, mpzspm->modulus);
  spv_size_t i;
  for (i = len - 1; i < 2 * len - 2; i++)
    if (mpz_cmp_ui (t[i], 0))
      printf ("error in ntt_PolyInvert\n");
  if (mpz_cmp_ui (t[2 * len - 2], 1))
    printf ("error in ntt_PolyInvert-\n");
#endif

  mpzspp_clear (w);
  mpzspp_clear (x);
  mpzspp_clear (y);
  mpzspp_clear (z);
}

void
ntt_polyevalT (mpzp_t b, spv_size_t len, mpzp_t *Tree, mpzp_t T,
                   mpzspp_t sp_invF, mpzspm_t mpzspm, char *TreeFilename)
{
  spv_size_t m, i;
    
  mpzspp_t x = mpzspp_init2 (mpzspm, len);
  mpzspp_t y = mpzspp_init2 (mpzspm, len);
  mpzspp_t z = mpzspp_init2 (mpzspm, 2 * len);

  mpzspp_set_mpzp (z, b, len, 0);
  mpzspp_to_ntt (z, 2 * len, 0);
  mpzspp_pwmul (z, z, sp_invF);
  mpzspp_from_ntt (z);
  mpzspp_normalise (z, len, len - 1);
  mpzspp_set (z, z, len, 0, len - 1);
  mpzspp_reverse (z, len);
    
  for (m = len / 2; m >= POLYEVALT_NTT_THRESHOLD; m /= 2)
    {
      for (i = 0; i < len; i += 2 * m)
        {
          list_revert (*Tree + i, m - 1);
          mpzspp_set_sp (x, 1, 1, 0);   /* x contains reversed monic poly */
          mpzspp_set_mpzp (x, *Tree + i, m, 1);
	  x->len = m + 1;
	  mpzspp_to_ntt (x, 2 * m, 0);
	  /* FIXME: could save this copy by doing transforms in z */
	  mpzspp_set (y, z, 2 * m, 0, i);
	  y->len = 2 * m;
	  mpzspp_to_ntt (y, 2 * m, 0);
	  mpzspp_pwmul (x, x, y);
	  mpzspp_from_ntt (x);
          mpzspp_normalise (x, m, m);
	  mpzspp_set (z, x, m, i + m, m);
	    
	  list_revert (*Tree + i + m, m - 1);
	  mpzspp_set_sp (x, 1, 1, 0);
	  mpzspp_set_mpzp (x, *Tree + i + m, m, 1);
	  x->len = m + 1;
	  mpzspp_to_ntt (x, 2 * m, 0);
	  mpzspp_pwmul (x, x, y);
	  mpzspp_from_ntt (x);
	  mpzspp_normalise (x, m, m);
	  mpzspp_set (z, x, m, i, m);
        }
      Tree++;
    }
    
  mpzspp_get_mpzp (z, T, len, 0);
    
  for (i = 0; i < len; i += 2 * m)
    TUpTree (T + i, Tree, 2 * m, T + len, -1, i, mpzspm->modulus, NULL);
    
  list_swap (b, T, len);
    
  mpzspp_clear (x);
  mpzspp_clear (y);
}
