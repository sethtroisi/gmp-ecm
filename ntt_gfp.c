/* ntt_gfp.c - low-level radix-2 dif/dit ntt routines over GF(p)
   
  Copyright 2005 Dave Newman.

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
  the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
  MA 02111-1307, USA.
*/
#include "sp.h"

#define X0 x[i]
#define X1 x[i + m]

#ifdef TUNE
#undef SPV_NTT_GFP_DIF_RECURSIVE_THRESHOLD
#undef SPV_NTT_GFP_DIT_RECURSIVE_THRESHOLD
spv_size_t SPV_NTT_GFP_DIF_RECURSIVE_THRESHOLD = 1;
spv_size_t SPV_NTT_GFP_DIT_RECURSIVE_THRESHOLD = 1;
#endif

void
ntt_scramble (spv_t x, spv_size_t len)
{
  spv_size_t i, j = 0, k, t;

  for (i = 0; i < len - 1; i++)
    {
      if (i < j)
        {
          t = x[i];
	  x[i] = x[j];
	  x[j] = t;
        }
      k = len / 2;
      while (k <= j)
        {
	  j -= k;
	  k /= 2;
	}
      j += k;
    }
}

void
spv_ntt_gfp_dif (spv_t x, spv_size_t len, sp_t p, sp_t d, sp_t root)
{
  sp_t a = 1, t;
  spv_size_t i, j, m;
	
  if (len < SPV_NTT_GFP_DIF_RECURSIVE_THRESHOLD)
    { 
      /* unrolled version for data that
         fits in the L1 cache */
      for (m = len / 2; m >= 1; m >>= 1)
        {
          a = 1;
	  for (j = 0; j < m; j++)
            {
              for (i = j; i < len; i += 2 * m)
                {
	          t = sp_sub (X0, X1, p);
	          X0 = sp_add (X0, X1, p);
	          X1 = sp_mul (t, a, p, d);
 	        }
              a = sp_mul (a, root, p, d);
            }
          root = sp_sqr (root, p, d);
        }
    }
  else
    {
      /* recursive version for data that
         doesn't fit in the L1 cache */
      m = len / 2;
      for (j = 0; j < m; j++)
        {
	  t = sp_sub (x[j], x[j + m], p);
	  x[j] = sp_add (x[j], x[j + m], p);
	  x[j + m] = sp_mul (t, a, p, d);
	  a = sp_mul (a, root, p, d);
	}
      root = sp_sqr (root, p, d);
      spv_ntt_gfp_dif (x, m, p, d, root);
      spv_ntt_gfp_dif (x + m, m, p, d, root);
    }
}

void
spv_ntt_gfp_dit (spv_t x, spv_size_t len, sp_t p, sp_t d, sp_t root)
{
  sp_t a, b, t;
  spv_size_t i, j, m;

  if (len < SPV_NTT_GFP_DIT_RECURSIVE_THRESHOLD)
    {
      for (m = 1; m < len; m <<= 1)
        {
          a = 1;
          b = sp_pow (root, len / (2 * m), p, d);
      
          for (j = 0; j < m; j++)
            {
	      for (i = j; i < len; i += 2 * m)
	        {
	          t = sp_mul (a, x[i + m], p, d);
	          x[i + m] = sp_sub (x[i], t, p);
	          x[i] = sp_add (x[i], t, p);
	        }
	  
	      a = sp_mul (a, b, p, d);
            }
        }
    }
  else
    {
      sp_t root2 = sp_sqr (root, p, d);
      m = len / 2;
      spv_ntt_gfp_dit (x, m, p, d, root2);
      spv_ntt_gfp_dit (x + m, m, p, d, root2);
      a = 1;
      
      for (j = 0; j < m; j++)
        {
	  t = sp_mul (a, x[j + m], p, d);
	  x[j + m] = sp_sub (x[j], t, p);
	  x[j] = sp_add (x[j], t, p);
	  a = sp_mul (a, root, p, d);
	}
    }
}

#if 0
void
spv_mul_ntt_gfp (spv_t r, spv_t x, spv_t y, spv_size_t len, spm_t spm)
{
  sp_t p, d, root;
  
  p = spm->sp;
  d = spm->mul_c;

  /* find primitive len'th root of unity mod p */
  root = sp_pow (spm->prim_root, (p - 1) / len, p, d);
  
  spv_ntt_gfp_dif (x, len, p, d, root);
  spv_ntt_gfp_dif (y, len, p, d, root);
  spv_pwmul (r, x, y, len, p, d);
  spv_ntt_gfp_dit (r, len, p, d, sp_inv (root, p, d));

  spv_mul_sp (r, r, p - (p - 1) / len, len, p, d);
}

void
spv_sqr_ntt_gfp (spv_t r, spv_t x, spv_size_t len, spm_t spm)
{
  sp_t p, d, root, len_inv;
  
  p = spm->sp;
  d = spm->mul_c;

  /* find primitive len'th root of unity mod p */
  root = sp_pow (spm->prim_root, (p - 1) / len, p, d);
  
  spv_ntt_gfp_dif (x, len, p, d, root);
  spv_pwmul (r, x, x, len, p, d);
  spv_ntt_gfp_dit (r, len, p, d, sp_inv (root, p, d));

  len_inv = sp_inv (len, p, d);
  spv_mul_sp (r, r, len_inv, len, p, d);
}
#endif
