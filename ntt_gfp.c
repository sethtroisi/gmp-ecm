#include "sp.h"

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
	
  if (4 * len <= CACHE_SIZE)
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
      /* iterative version for data that
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

  /* not cache-friendly at all, for large len it's 
   * quicker to use scramble + dif + scramble */
  if (len >= DIT_DIF_THRESHOLD)
    {
      ntt_scramble (x, len);
      spv_ntt_gfp_dif (x, len, p, d, root);
      ntt_scramble (x, len);
      return;
    }
  
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

void
spv_mul_ntt_gfp (spv_t r, spv_t x, spv_t y, spv_size_t len, spm_t spm)
{
  sp_t p, d, root, len_inv;
  
  p = spm->sp;
  d = spm->mul_c;

  /* find primitive len'th root of unity mod p */
  root = sp_pow (spm->prim_root, (p - 1) / len, p, d);
  
  spv_ntt_gfp_dif (x, len, p, d, root);
  spv_ntt_gfp_dif (y, len, p, d, root);
  spv_pwmul (r, x, y, len, p, d);
  spv_ntt_gfp_dit (r, len, p, d, sp_inv (root, p, d));

  len_inv = sp_inv (len, p, d);
  spv_mul_sp (r, r, len_inv, len, p, d);

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

