/* ecm_ntt.c - high level poly functions to interface between ecm and sp

Copyright 2005, 2006, 2007, 2008, 2009, 2011, 2012 Dave Newman,
Paul Zimmermann, Alexander Kruppa.

This file is part of the ECM Library.

The ECM Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The ECM Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the ECM Library; see the file COPYING.LIB.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ecm-impl.h"

#ifdef HAVE_UNISTD_H
#include <unistd.h> /* for unlink */
#endif

#define UNUSED 0

/* memory: 4 * len mpspv coeffs */
void
ntt_mul (mpzv_t r, mpzv_t x, mpzv_t y, spv_size_t len, mpzv_t t, 
         mpzspm_t mpzspm)
{
  mpzspv_handle_t u, v;
  
  if (len < MUL_NTT_THRESHOLD)
    {
      list_mul (r, x, len, 0, y, len, 0, t);
      return;
    }

  u = mpzspv_init_handle (NULL, 2 * len, mpzspm);
  v = mpzspv_init_handle (NULL, 2 * len, mpzspm);
  
  mpzspv_fromto_mpzv (v, 0, len, NULL, y, NULL, NULL);
  mpzspv_fromto_mpzv (u, 0, len, NULL, x, NULL, NULL);

  mpzspv_mul_ntt(v, 0, 
                 v, 0, len, 
                 NULL, UNUSED, UNUSED, 
                 2 * len, NTT_MUL_STEP_FFT1);
  mpzspv_mul_ntt(u, 0, 
                 u, 0, len, 
                 v, 0, 2*len, 
                 2 * len, 
                 NTT_MUL_STEP_FFT1 + NTT_MUL_STEP_MUL + NTT_MUL_STEP_IFFT);
  mpzspv_fromto_mpzv (u, 0, 2 * len - 1, NULL, NULL, NULL, r);
  
  mpzspv_clear_handle (u);
  mpzspv_clear_handle (v);
}

/* memory: 2 * len mpzspv coeffs */
void
ntt_PolyFromRoots (mpzspv_handle_t r, spv_size_t offset_r, mpzv_t a, 
    spv_size_t len, mpzv_t t, mpzspm_t mpzspm)
{
  mpzspv_handle_t x;
  spv_size_t i, m;
  
  ASSERT (len == ((spv_size_t)1) << ceil_log2 (len));

  if (len <= MUL_NTT_THRESHOLD)
  {
    /* FIXME: This modifies the input array */
    PolyFromRoots (a, a, len, t, mpzspm->modulus);
    mpzspv_fromto_mpzv (r, offset_r, len, NULL, a, NULL, NULL);
    return;
  }
  
  x = mpzspv_init_handle (NULL, 2 * len, mpzspm);
  
  for (i = 0; i < len; i += MUL_NTT_THRESHOLD)
    {
      /* FIXME: This modifies the input array */
      PolyFromRoots (a + i, a + i, MUL_NTT_THRESHOLD, t, mpzspm->modulus);
      mpzspv_fromto_mpzv (x, 2 * i, MUL_NTT_THRESHOLD, NULL, a + i, NULL, NULL);
    }
  
  for (m = MUL_NTT_THRESHOLD; m < len; m *= 2)
    {
      for (i = 0; i < 2 * len; i += 4 * m)
        {
          mpzspv_set_sp (x, i + 3 * m, (sp_t) 1, 1);
	  mpzspv_mul_ntt (x, i + 2 * m, 
	                  x, i + 2 * m, m + 1, 
	                  NULL, UNUSED, UNUSED, 
	                  2 * m, 
	                  NTT_MUL_STEP_FFT1);
          mpzspv_set_sp (x, i + m, (sp_t) 1, 1);
	  mpzspv_mul_ntt (x, i, 
	                  x, i, m + 1, 
	                  x, i + 2 * m, 2 * m, 
	                  2 * m, 
	                  NTT_MUL_STEP_FFT1 + NTT_MUL_STEP_MUL + NTT_MUL_STEP_IFFT);
          /* Subtract wrapped-around product 1*1 of leading monomials */
	  mpzspv_sub_sp (x, i, x, i, (sp_t) 1, 1);

          mpzspv_fromto_mpzv (x, i, 2 * m, NULL, NULL, NULL, NULL);
	}	  
    }
  
  mpzspv_set (r, offset_r, x, 0, len);
  mpzspv_clear_handle (x);
}

  
/* memory: 2 * len mpzspv coeffs */
int
ntt_PolyFromRoots_Tree (mpzv_t r, mpzv_t a, spv_size_t len, mpzv_t t,
    int dolvl, mpzspm_t mpzspm, mpzv_t *Tree, FILE *TreeFile)
{
  mpzspv_handle_t x;
  spv_size_t i, m, m_max;
  mpzv_t src;
  mpzv_t *dst = Tree + ceil_log2 (len) - 1;

  ASSERT (len == ((spv_size_t)1) << ceil_log2 (len));
  
  x = mpzspv_init_handle (NULL, 2 * len, mpzspm);
  
  if (dolvl >= 0)
    {
      src = a;
      dst = &r;
    }
  else  
    {
      /* Copy the roots into the destination level of the tree (negating
	 if so desired), set the source to this level (which now contains 
	 the possibly negated roots), and advance the destination level 
	 of the tree to the next level */
      src = *dst;
      /* we consider x + root[i], which means we consider negated roots */
      list_set (*dst--, a, len);
    }
  
  m = (dolvl == -1) ? 1 : 1 << (ceil_log2 (len) - 1 - dolvl);
  m_max = (dolvl == -1) ? len : 2 * m;
  
  for (; m < m_max && m < MUL_NTT_THRESHOLD; m *= 2)
    {
      /* dst = &r anyway for dolvl != -1 */
      if (m == len / 2)
	dst = &r;
      
      if (TreeFile && list_out_raw (TreeFile, src, len) == ECM_ERROR)
        {
          outputf (OUTPUT_ERROR, "Error writing product tree of F\n");
          return ECM_ERROR;
        }

      for (i = 0; i < len; i += 2 * m)
	list_mul (t + i, src + i, m, 1, src + i + m, m, 1, t + len);

      list_mod (*dst, t, len, mpzspm->modulus);
      
      src = *dst--;
    }
  
  for (; m < m_max; m *= 2)
    {
      ASSERT (m > 1); /* This code does not do the sign change. Let's assume
			 MUL_NTT_THRESHOLD is always large enough that the
			 degree 1 product are done in the above loop */
      /* dst = &r anyway for dolvl != -1 */
      if (m == len / 2)
        dst = &r;
      
      for (i = 0; i < 2 * len; i += 4 * m)
        {
 	  if (TreeFile &&
	      list_out_raw (TreeFile, src + i / 2, 2 * m) == ECM_ERROR)
	    return ECM_ERROR;
	  
	  mpzspv_fromto_mpzv (x, i, m, NULL, src + i / 2, NULL, NULL);
	  mpzspv_fromto_mpzv (x, i + 2 * m, m, NULL, src + i / 2 + m, NULL, NULL);
	  
	  mpzspv_set_sp (x, i + 3 * m, (sp_t) 1, 1);
          mpzspv_mul_ntt (x, i + 2 * m, 
                          x, i + 2 * m, m + 1, 
                          NULL, UNUSED, UNUSED, 
                          2 * m, NTT_MUL_STEP_FFT1);
	  mpzspv_set_sp (x, i + m, (sp_t) 1, 1);
          mpzspv_mul_ntt (x, i, 
                          x, i, m + 1, 
                          x, i + 2 * m, 2*m, 
                          2 * m, 
                          NTT_MUL_STEP_FFT1 + NTT_MUL_STEP_MUL + NTT_MUL_STEP_IFFT);
          /* Subtract wrapped-around product 1*1 of leading monomials */
          mpzspv_sub_sp (x, i, x, i, (sp_t) 1, 1);
          mpzspv_fromto_mpzv (x, i, 2 * m, NULL, NULL, NULL, *dst + i / 2);

          /* we only do the mod reduction to reduce the file size a bit */
	  if (TreeFile)
	    list_mod (*dst + i / 2, *dst + i / 2, 2 * m, mpzspm->modulus);
	}
    
      src = *dst--;
    }

  mpzspv_clear_handle (x);

  return 0;
}

/* 2 NTTs of size 2 * len
 * 2 NTTs of size len
 *
 * memory: 2 * len mpzspv coeffs */
void
ntt_PrerevertDivision (mpzspv_handle_t a, mpzspv_handle_t sp_b, 
    mpzspv_handle_t sp_invb, spv_size_t len, mpzspm_t mpzspm)
{
  mpzspv_handle_t sum;
  
  sum = mpzspv_init_handle (NULL, len, mpzspm);

  /* sum = a mod (x^len - 1) */
  mpzspv_add (sum, 0, a, 0, a, len, len);

  /* a = TOP (TOP (a) * invb) = trunc(a/f) */
  mpzspv_set (a, 0, a, len, len);
  mpzspv_mul_ntt (a, 0, 
                  a, 0, len, 
                  sp_invb, 0, 2 * len, 
                  2 * len, 
                  NTT_MUL_STEP_FFT1 + NTT_MUL_STEP_MUL + NTT_MUL_STEP_IFFT);
  mpzspv_fromto_mpzv (a, len - 1, len, NULL, NULL, NULL, NULL);
  mpzspv_set (a, 0, a, len - 1, 1);  /* In two pieces to avoid overlap */
  mpzspv_set (a, 1, a, len, len - 1);
  
  /* a := (trunc(a/f) * f) mod (x^len - 1) */
  mpzspv_mul_ntt (a, 0, 
                  a, 0, len, 
                  sp_b, 0, len, 
                  len, 
                  NTT_MUL_STEP_FFT1 + NTT_MUL_STEP_MUL + NTT_MUL_STEP_IFFT);

  /* a := (a - trunc(a/f) * f) mod (x^len - 1) = a mod f */
  mpzspv_sub (a, 0, sum, 0, a, 0, len);
  /* Add modulus to ensure difference is non-negative */
  mpzspv_add_mpz (a, 0, a, 0, mpzspm->modulus, len);
  mpzspv_fromto_mpzv (a, 0, len, NULL, NULL, NULL, NULL);
  
  mpzspv_clear_handle (sum);
  
#if 0
  list_sub (t, t, a + len, len);
  list_sub (a, a, t, len);
  /* can we avoid this mod without risking overflow later? */
  list_mod (a, a, len, mpzspm->modulus);
#endif
}

void
cmp_mpz (void *s, const mpz_t r)
{
  mpz_t **p = s;
  ASSERT (mpz_cmp(**p, r) == 0);
  (*p)++;
}

/* Put in q[0, len - 1] coefficients s.t.
   q[0, len - 1] * b[0, len - 1] \in x^(2k-2) + O(x^(k-2))
   Here O(x^n) means the set of polynomials of degree at most n. */
/* memory: 7/2 * len mpzspv coeffs */
void ntt_PolyInvert (mpzspv_handle_t r, mpzspv_handle_t sp_b, 
    spv_size_t offset_b, spv_size_t len, mpzspm_t mpzspm)
{
  spv_size_t k = 1;
  mpzspv_handle_t w, x, z;
  const int alloc_x = !mpzspv_handle_in_memory(r);
  
  w = mpzspv_init_handle (NULL, len / 2, mpzspm);
  if (alloc_x)
    x = mpzspv_init_handle (NULL, len, mpzspm);
  else
    x = r;
  z = mpzspv_init_handle (NULL, len, mpzspm);
  
  /* We assume that b is monic, thus the degree-0 approximation of the 
     reciprocal is 1 */
  mpzspv_set_sp (x, 0, (sp_t) 1, 1);
  /* This discards the leading 1 of b */
  
  for (; k < len; k *= 2)
    {
      mpzspv_set (w, 0, x, 0, k);
      mpzspv_set (z, 0, sp_b, offset_b + len - 2 * k, 2 * k - 1);
      mpzspv_mul_ntt (x, 0, 
                      x, 0, k, 
                      NULL, UNUSED, UNUSED, 
                      2 * k, NTT_MUL_STEP_FFT1);
      mpzspv_mul_ntt (z, 0, 
                      z, 0, 2 * k - 1, 
                      x, 0, 2 * k, 
                      2 * k, 
                      NTT_MUL_STEP_FFT1 + NTT_MUL_STEP_MUL + NTT_MUL_STEP_IFFT);
      mpzspv_fromto_mpzv (z, k, k, NULL, NULL, NULL, NULL);
      mpzspv_neg (z, 0, z, k - 1, 1); /* In two pieces to avoid overlap */
      mpzspv_neg (z, 1, z, k, k - 1);
      
      mpzspv_mul_ntt (x, 0, 
                      z, 0, k, 
                      x, 0, 2 * k, 
                      2 * k, 
                      NTT_MUL_STEP_FFT1 + NTT_MUL_STEP_MUL + NTT_MUL_STEP_IFFT);
      mpzspv_fromto_mpzv (x, k - 1, k, NULL, NULL, NULL, NULL);
      mpzspv_set (x, 0, x, k - 1, 1); /* In two pieces to avoid overlap */
      if (k > 1)
        mpzspv_set (x, 1, x, k, k - 1);
      mpzspv_set (x, k, w, 0, k);
    }

  if (alloc_x)
    mpzspv_set (r, 0, x, 0, len);
 
  mpzspv_clear_handle (w);
  if (alloc_x)
    mpzspv_clear_handle (x);
  mpzspv_clear_handle (z);

#if defined DEBUG
  ntt_mul (t, q, b, len, NULL, mpzspm);
  list_mod (t, t, 2 * len - 1, mpzspm->modulus);
  
  {
  spv_size_t i;
  for (i = len - 1; i < 2 * len - 2; i++)
    if (mpz_cmp_ui (t[i], 0))
      {
        fprintf (stderr, "error in ntt_PolyInvert %" PRISPVSIZE "\n", i);
        abort();
      }
  }
  if (mpz_cmp_ui (t[2 * len - 2], 1))
    {
      fprintf (stderr, "error in ntt_PolyInvert\n");
      abort();
    }
#endif
}


/* memory: 4 * len mpzspv coeffs */
int
ntt_polyevalT (mpzv_t r, mpzspv_handle_t sp_b, spv_size_t len, mpzv_t *Tree, mpzv_t T,
               mpzspv_handle_t sp_invF, mpzspm_t mpzspm, char *TreeFilenameStem)
{
  spv_size_t m, i;
  FILE *TreeFile = NULL;
  /* assume this "small" malloc will not fail in normal usage */
  char *TreeFilename = NULL;
  mpzv_t *Tree_orig = Tree;
  int level = 0; /* = ceil_log2 (len / m) - 1 */
  mpzspv_handle_t x; 
  mpzspv_handle_t y; 

  x = mpzspv_init_handle (NULL, 2 * len, mpzspm);
  y = mpzspv_init_handle (NULL, 2 * len, mpzspm);

  if (TreeFilenameStem)
    {
      TreeFilename = (char *) malloc (strlen (TreeFilenameStem) + 1 + 2 + 1);
      if (TreeFilename == NULL)
        {
          fprintf (stderr, "Cannot allocate memory in ntt_polyevalT\n");
          exit (1);
        }
    }
  
  mpzspv_set (x, 0, sp_b, 0, len);
  mpzspv_mul_ntt(x, 0, 
                 x, 0, len, 
                 sp_invF, 0, 2 * len, 
                 2 * len, 
                 NTT_MUL_STEP_FFT1 + NTT_MUL_STEP_MUL + NTT_MUL_STEP_IFFT);
  mpzspv_fromto_mpzv (x, len-1, len, NULL, NULL, NULL, NULL);
  mpzspv_set (y, 0, x, len - 1, len); /* y = high (b * invF) */
  mpzspv_reverse (y, 0, y, 0, len); /* y = rev (high (b * invF)) */
  
  for (m = len / 2; m >= POLYEVALT_NTT_THRESHOLD; m /= 2)
    {
      if (TreeFilenameStem)
        {
          Tree = &T;
	  
          sprintf (TreeFilename, "%s.%d", TreeFilenameStem, level);
          
	  TreeFile = fopen (TreeFilename, "rb");
          if (TreeFile == NULL)
            {
              outputf (OUTPUT_ERROR,
		  "Error opening file %s for product tree of F\n",
                        TreeFilename);
              mpzspv_clear_handle (x);
	      mpzspv_clear_handle (y);
	      return ECM_ERROR;
            }

	  list_inp_raw (*Tree, TreeFile, len);

	  fclose (TreeFile);
	  unlink (TreeFilename);
	}

      for (i = 0; i < len; i += 2 * m)
        {
	  list_revert (*Tree + i, m);
          mpzspv_set_sp (x, 0, 1, 1);
          mpzspv_fromto_mpzv (x, 1, m, NULL, *Tree + i, NULL, NULL);
	  /* x contains reversed monic poly */
	  mpzspv_mul_ntt (y, i, 
	                  y, i, 2 * m, 
	                  NULL, UNUSED, UNUSED, 
	                  2 * m, NTT_MUL_STEP_FFT1);
          mpzspv_mul_ntt (x, 0, 
                          x, 0, m + 1, 
                          y, i, 2 * m, 
                          2 * m, 
                          NTT_MUL_STEP_FFT1 + NTT_MUL_STEP_MUL + NTT_MUL_STEP_IFFT);
          if (m > POLYEVALT_NTT_THRESHOLD)
	    mpzspv_fromto_mpzv (x, m, m, NULL, NULL, NULL, NULL);
	    
	  list_revert (*Tree + i + m, m);
	  mpzspv_set_sp (x, 2 * m, 1, 1);
	  mpzspv_fromto_mpzv (x, 2 * m + 1, m, NULL, *Tree + i + m, NULL, NULL);
          mpzspv_mul_ntt(x, 2 * m, 
                         x, 2 * m, m + 1, 
                         y, i, 2 * m, 
                         2 * m, 
                         NTT_MUL_STEP_FFT1 + NTT_MUL_STEP_MUL + NTT_MUL_STEP_IFFT);
	  if (m > POLYEVALT_NTT_THRESHOLD)
	    mpzspv_fromto_mpzv (x, 3 * m, m, NULL, NULL, NULL, NULL);
	  
	  mpzspv_set (y, i, x, 3 * m, m);
	  mpzspv_set (y, i + m, x, m, m);
        }
      
      Tree++;
      level++;
    }
    
  mpzspv_clear_handle (x);
  mpzspv_fromto_mpzv (y, 0, len, NULL, NULL, NULL, T); /* T = rev (high (b * invF)) */
  mpzspv_clear_handle (y);
  for (i = 0; i < len; i++)
    mpz_mod (T[i], T[i], mpzspm->modulus);

  for (; m >= 1; m /= 2)
    {
      if (TreeFilenameStem)
        {
          sprintf (TreeFilename, "%s.%d", TreeFilenameStem, level);

          TreeFile = fopen (TreeFilename, "rb");
          if (TreeFile == NULL)
            {
              outputf (OUTPUT_ERROR,
		  "Error opening file %s for product tree of F\n",
                        TreeFilename);
	      return ECM_ERROR;
            }
	}
      
      TUpTree (T, Tree_orig, len, T + len, level++, 0,
	  mpzspm->modulus, TreeFile);

      if (TreeFilenameStem)
        {
	  fclose (TreeFile);
	  unlink (TreeFilename);
	}
    }
  
  if (TreeFilenameStem)
    free (TreeFilename);
  list_set (r, T, len);
  return 0;
}
