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

#if 1
void
ntt_mul (mpzv_t r, mpzv_t x, mpzv_t y, spv_size_t len, mpzv_t t,
    int monic, mpzspm_t mpzspm)
{
  if (len < MUL_NTT_THRESHOLD)
  {
    list_mul (r, x, len, monic, y, len, monic, t);
    return;
  }

  mpzspv_t u = mpzspv_init (2 * len, mpzspm);
  mpzspv_t v = mpzspv_init (2 * len, mpzspm);
  
  mpzspv_from_mpzv (u, 0, x, len, mpzspm);
  mpzspv_to_ntt (u, 0, len, 2 * len, monic, mpzspm);
  mpzspv_from_mpzv (v, 0, y, len, mpzspm);
  mpzspv_to_ntt (v, 0, len, 2 * len, monic, mpzspm);
  mpzspv_pwmul (u, 0, u, 0, v, 0, 2 * len, mpzspm);
  mpzspv_from_ntt (u, 0, 2 * len, monic ? 2 * len : 0, mpzspm);
  mpzspv_to_mpzv (u, 0, r, 2 * len - 1 + monic, mpzspm);
  
  mpzspv_clear (u, mpzspm);
  mpzspv_clear (v, mpzspm);
}
#endif

void
ntt_PolyFromRoots (mpzv_t r, mpzv_t a, spv_size_t len, mpzv_t t,
    mpzspm_t mpzspm)
{
  ASSERT (len == 1 << ceil_log_2 (len));

  if (len <= MUL_NTT_THRESHOLD)
  {
    PolyFromRoots (r, a, len, t, mpzspm->modulus);
    return;
  }
  
  mpzspv_t x = mpzspv_init (2 * len, mpzspm);
  spv_size_t i, m;
   
  for (i = 0; i < len; i += MUL_NTT_THRESHOLD)
    {
      PolyFromRoots (r, a + i, MUL_NTT_THRESHOLD, t, mpzspm->modulus);
      mpzspv_from_mpzv (x, 2 * i, r, MUL_NTT_THRESHOLD, mpzspm);
    }
  
  for (m = MUL_NTT_THRESHOLD; m < len; m *= 2)
    {
      for (i = 0; i < 2 * len; i += 4 * m)
        {
	  mpzspv_to_ntt (x, i, m, 2 * m, 1, mpzspm);
	  mpzspv_to_ntt (x, i + 2 * m, m, 2 * m, 1, mpzspm);
	  mpzspv_pwmul (x, i, x, i, x, i + 2 * m, 2 * m, mpzspm);
	  mpzspv_from_ntt (x, i, 2 * m, 2 * m, mpzspm);
	  
	  if (2 * m < len)
	    mpzspv_normalise (x, i, 2 * m, mpzspm);
	}	  
    }
      
  mpzspv_to_mpzv (x, 0, r, len, mpzspm);

  mpzspv_clear (x, mpzspm);
}

  
int
ntt_PolyFromRoots_Tree (mpzv_t r, mpzv_t a, spv_size_t len, mpzv_t t,
    mpzspm_t mpzspm, mpzv_t *Tree, FILE *TreeFile)
{
  ASSERT (len == 1 << ceil_log_2 (len));

  if (len <= MUL_NTT_THRESHOLD)
    {
      return PolyFromRoots_Tree
        (r, a, len, t, -1, mpzspm->modulus, Tree, TreeFile, 0);
    }
  
  mpzspv_t x = mpzspv_init (2 * len, mpzspm);
  spv_size_t i, m;
  mpzv_t src;
  mpzv_t *dst = Tree + ceil_log_2 (len) - 1;
  
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
      
      for (i = 0; i < 2 * len; i += 4 * m)
        {
	  mpzspv_from_mpzv (x, i, src + i / 2, m, mpzspm);
	  mpzspv_to_ntt (x, i, m, 2 * m, 1, mpzspm);
	  mpzspv_from_mpzv (x, i + 2 * m, src + i / 2 + m, m, mpzspm);
	  mpzspv_to_ntt (x, i + 2 * m, m, 2 * m, 1, mpzspm);
	  mpzspv_pwmul (x, i, x, i, x, i + 2 * m, 2 * m, mpzspm);
	  mpzspv_from_ntt (x, i, 2 * m, 2 * m, mpzspm);
	  
          mpzspv_to_mpzv (x, i, *dst + i / 2, 2 * m, mpzspm);

          /* this line can probably be commented out without causing
	   * overflow; will things slow down though? */
#if 0          
	  list_mod (*dst + i / 2, *dst + i / 2, 2 * m, mpzspm->modulus);
#endif
	  if (TreeFile && mpz_out_raw (TreeFile, dst[0][i / 2]) == 0)
	    return ECM_ERROR;
	}
    }

  mpzspv_clear (x, mpzspm);

  return 0;
}


/* 2 NTTs of size 2 * len
 * 3 NTTs of size len */
void
ntt_PrerevertDivision (mpzv_t a, mpzv_t b, mpzv_t invb, mpzspv_t sp_invb,
    spv_size_t len, mpzv_t t, mpzspm_t mpzspm)
{
  if (len < PREREVERT_DIVISION_NTT_THRESHOLD)
    {
      PrerevertDivision (a, b, invb, len, t, mpzspm->modulus);
      return;
    }
  
  mpzspv_t x = mpzspv_init (2 * len, mpzspm);
  mpzspv_t y = mpzspv_init (2 * len, mpzspm);

  /* y = TOP (TOP (a) * invb) */
  mpzspv_set_sp (x, 0, 0, len + 1, mpzspm);
  mpzspv_from_mpzv (x, len + 1, a + len, len - 1, mpzspm);
  mpzspv_to_ntt (x, 0, 2 * len, 2 * len, 0, mpzspm);
  mpzspv_pwmul (x, 0, x, 0, sp_invb, 0, 2 * len, mpzspm);
  mpzspv_from_ntt (x, 0, 2 * len, 0, mpzspm);
  mpzspv_normalise (x, 0, len, mpzspm);
  
  mpzspv_to_ntt (x, 0, len, len, 0, mpzspm);
  mpzspv_from_mpzv (y, 0, b, len, mpzspm);
  mpzspv_to_ntt (y, 0, len, len, 1, mpzspm); /* b has leading monomial */
  mpzspv_pwmul (x, 0, x, 0, y, 0, len, mpzspm);
  mpzspv_from_ntt (x, 0, len, 0, mpzspm);
  mpzspv_to_mpzv (x, 0, t, len, mpzspm);
  
  mpzspv_clear (x, mpzspm);
  mpzspv_clear (y, mpzspm);
 
  list_sub (t, t, a + len, len - 1);
  list_sub (a, a, t, len);
  /* can we avoid this mod without risking overflow later? */
  list_mod (a, a, len, mpzspm->modulus);
}

void ntt_PolyInvert (mpzv_t q, mpzv_t b, spv_size_t len, mpzv_t t,
    mpzspm_t mpzspm)
{
  if (len < POLYINVERT_NTT_THRESHOLD)
    {
      PolyInvert (q, b, len, t, mpzspm->modulus);
      return;
    }

  spv_size_t k = POLYINVERT_NTT_THRESHOLD / 2;
  
  PolyInvert (q + len - k, b + len - k, k, t, mpzspm->modulus);
  
  mpzspv_t w = mpzspv_init (len / 2, mpzspm);
  mpzspv_t x = mpzspv_init (len, mpzspm);
  mpzspv_t y = mpzspv_init (len, mpzspm);
  mpzspv_t z = mpzspv_init (len, mpzspm);
  
  mpzspv_from_mpzv (x, 0, q + len - k - 1, k + 1, mpzspm);
  mpzspv_from_mpzv (y, 0, b, len - 1, mpzspm);
  
  for (; k < len; k *= 2)
    {
      mpzspv_set (w, 0, x, 1, k, mpzspm);
      mpzspv_to_ntt (x, 0, k + 1, 2 * k, 0, mpzspm);
      mpzspv_set (z, 0, y, len - 2 * k, 2 * k - 1, mpzspm);
      mpzspv_to_ntt (z, 0, 2 * k - 1, 2 * k, 0, mpzspm);
      mpzspv_pwmul (z, 0, x, 0, z, 0, 2 * k, mpzspm);
      mpzspv_from_ntt (z, 0, 2 * k, 0, mpzspm);
      mpzspv_normalise (z, k, k, mpzspm);
      mpzspv_neg (z, 0, z, k, k, mpzspm);
      
      mpzspv_to_ntt (z, 0, k, 2 * k, 0, mpzspm);
      mpzspv_pwmul (x, 0, x, 0, z, 0, 2 * k, mpzspm);
      mpzspv_from_ntt (x, 0, 2 * k, 0, mpzspm);
      if (2 * k < len)
	mpzspv_normalise (x, k, k, mpzspm);
      mpzspv_set (x, 1, x, k, k, mpzspm); /* legal overlap */
      mpzspv_set (x, k + 1, w, 0, MIN(k, len / 2), mpzspm);
    }

  mpzspv_to_mpzv (x, 1, q, len - POLYINVERT_NTT_THRESHOLD / 2, mpzspm);
 
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

  mpzspv_clear (w, mpzspm);
  mpzspv_clear (x, mpzspm);
  mpzspv_clear (y, mpzspm);
  mpzspv_clear (z, mpzspm);
}

void
ntt_polyevalT (mpzv_t b, spv_size_t len, mpzv_t *Tree, mpzv_t T,
                   mpzspv_t sp_invF, mpzspm_t mpzspm, char *TreeFilename)
{
  spv_size_t m, i;
    
  mpzspv_t x = mpzspv_init (2 * len, mpzspm);
  mpzspv_t y = mpzspv_init (2 * len, mpzspm);

  mpzspv_from_mpzv (x, 0, b, len, mpzspm);
  mpzspv_to_ntt (x, 0, len, 2 * len, 0, mpzspm);
  mpzspv_pwmul (x, 0, x, 0, sp_invF, 0, 2 * len, mpzspm);
  mpzspv_from_ntt (x, 0, 2 * len, 0, mpzspm);
  mpzspv_normalise (x, len - 1, len, mpzspm);
  mpzspv_set (y, 0, x, len - 1, len, mpzspm);
  mpzspv_reverse (y, 0, len, mpzspm);
    
  for (m = len / 2; m >= POLYEVALT_NTT_THRESHOLD; m /= 2)
    {
      for (i = 0; i < len; i += 2 * m)
        {
       	  mpzspv_to_ntt (y, i, 2 * m, 2 * m, 0, mpzspm);
	  
	  list_revert (*Tree + i, m - 1);
          mpzspv_set_sp (x, 0, 1, 1, mpzspm);
          mpzspv_from_mpzv (x, 1, *Tree + i, m, mpzspm);
	  /* x contains reversed monic poly */
	  mpzspv_to_ntt (x, 0, m + 1, 2 * m, 0, mpzspm);
	  mpzspv_pwmul (x, 0, x, 0, y, i, 2 * m, mpzspm);
	  mpzspv_from_ntt (x, 0, 2 * m, 0, mpzspm);
          if (m > POLYEVALT_NTT_THRESHOLD)
	    mpzspv_normalise (x, m, m, mpzspm);
	    
	  list_revert (*Tree + i + m, m - 1);
	  mpzspv_set_sp (x, 2 * m, 1, 1, mpzspm);
	  mpzspv_from_mpzv (x, 2 * m + 1, *Tree + i + m, m, mpzspm);
	  mpzspv_to_ntt (x, 2 * m, m + 1, 2 * m, 0, mpzspm);
	  mpzspv_pwmul (x, 2 * m, x, 2 * m, y, i, 2 * m, mpzspm);
	  mpzspv_from_ntt (x, 2 * m, 2 * m, 0, mpzspm);
	  if (m > POLYEVALT_NTT_THRESHOLD)
	    mpzspv_normalise (x, 3 * m, m, mpzspm);
	  
	  mpzspv_set (y, i, x, 3 * m, m, mpzspm);
	  mpzspv_set (y, i + m, x, m, m, mpzspm);
        }
      Tree++;
    }
    
  mpzspv_clear (x, mpzspm);
  mpzspv_to_mpzv (y, 0, T, len, mpzspm);
  mpzspv_clear (y, mpzspm);
    
  for (i = 0; i < len; i += 2 * m)
    TUpTree (T + i, Tree, 2 * m, T + len, -1, i, mpzspm->modulus, NULL);
    
  list_swap (b, T, len);
}
