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
  the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
  MA 02110-1301, USA.
*/

#include "sp.h"
#include "ecm-impl.h"

#define X0 x[i]
#define X1 x[i + m]

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

/* Take the r-bit integer in "a" and return its bit-reversal. If
   a >= 2^r, only the r lower bits are considered */

static unsigned long ATTRIBUTE_CONST
bitrev (const unsigned int r, unsigned long a)
{
  unsigned long t = 0;
  unsigned int i;

  for (i = 0; i < r; i++)
    {
      t += a & 1UL;
      a >>= 1;
    }
  return t;
}

static unsigned long ATTRIBUTE_CONST
aloc (const unsigned int r, const unsigned long a)
{
  return bitrev (r, (a >> 2) ^ (a >> 1));
}

static sp_t ATTRIBUTE_CONST
sp_V2 (const sp_t a, const sp_t p, const sp_t d)
{
  sp_t r;
  r = sp_mul (a, a, p, d);
  r = sp_sub (r, 2, p);
  return r;
}

void 
spv_ntt_dct2 (spv_t q, spv_size_t len, sp_t p, sp_t d, sp_t root,
		   sp_t invroot)
{
  const unsigned int r = ceil_log2 (len);
  unsigned int s;
  unsigned long a;
  sp_t gamma, v, v_1, vn, t1, t2, t3, t4;
  unsigned long i;

  /* Check that root and invroot are each a 4*len-th primitive root of unity */
  ASSERT (sp_pow (root, 2*len, p, d) == p - 1);
  ASSERT (sp_pow (invroot, 2*len, p, d) == p - 1);
  /* Check that they are inverses of each other */
  ASSERT (sp_mul (root, invroot, p, d) == 1);
  /* Check that transform length is a power of 2 */
  ASSERT (len = 1<<r);
  
  gamma = sp_add (root, invroot, p);
    printf ("/* spv_ntt_dct2 */ p = %lu;\n", p);
  printf ("/* spv_ntt_dct2 */ gam = Mod(%lu, p);\n", gamma);
  
  for (s = 1; s <= r; s++)
    {
      const unsigned long n = 1UL << (r - s);
      /* We have 2^(s-1) pieces, each of length 2n. These will be split into
         2^s pieces, each of length n */

      /* Compute vn = V_n(gamma), v = V_{an}(gamma), v_1 = V_{(a-1)n}(gamma)
	 for a = 1 */
      vn = gamma;
      for (i = 0; i < r - s; i++)
	vn = sp_V2 (vn, p, d);
      v_1 = 2;
      v = vn;
      for (a = 1UL; a < (1UL << s); a += 2)
	{
	  const unsigned long k = aloc (r, a);
	  const spv_t A = q + k, B = A, C = A + n;

	  /* b[0 ... n-1] is at a[0 ... n-1], 
	     c[0 ... n-1] is at a[n ... 2n-1] */
	  
	  t1 = sp_mul (v, A[n], p, d);
	  A[n] = sp_sub (A[0], t1, p); /* c_0 = a_0 - v*a_n */
	  A[0] = sp_add (A[0], t1, p); /* b_0 = a_0 + v*a_n */

	  for (i = 1; i <= n / 2; i++)
	  {
	      t1 = sp_mul (v, A[n + i], p, d);
	      t2 = sp_mul (v, A[2 * n - i], p, d);
	      t3 = sp_sub (A[i], A[2 * n - i], p);
	      t4 = sp_sub (A[n - i], A[n + i], p);
	      B[i] = sp_add (t3, t1, p);
	      B[n - i] = sp_add (t4, t2, p);
	      C[i] = sp_sub (t3, t1, p);
	      C[n - i] = sp_sub (t4, t2, p);
	  }

	  /* Update v_1 and v. V_{n+1} = V_n * V_1 - V_{n-1} */
	  {
	    sp_t t = v;
	    v = sp_mul (v, vn, p, d);
	    v = sp_sub (v, v_1, p);
	    v_1 = t;
	  }
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
