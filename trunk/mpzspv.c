/* mpzspv.c - "mpz small prime polynomial" functions for arithmetic on mpzv's
   reduced modulo a mpzspm

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

#include <malloc.h>
#include <string.h> /* for memset */
#include "sp.h"

mpzspv_t
mpzspv_init (spv_size_t len, mpzspm_t mpzspm)
{
  unsigned int i;
  mpzspv_t x = (mpzspv_t) malloc (mpzspm->sp_num * sizeof (spv_t));
  
  if (x == NULL)
    return NULL;
  
  for (i = 0; i < mpzspm->sp_num; i++)
    {
      x[i] = (spv_t) valloc (len * sizeof (sp_t));
      
      if (x[i] == NULL)
	{
	  while (i--)
	    free (x[i]);
	  
	  free (x);
	  return NULL;
	}
    }
  
  return x;
}

void
mpzspv_clear (mpzspv_t x, mpzspm_t mpzspm)
{
  unsigned int i;
	
  ASSERT (mpzspv_verify (x, 0, 0, mpzspm));
  
  for (i = 0; i < mpzspm->sp_num; i++)
    free (x[i]);
  
  free (x);
}

/* check that:
 *  - each of the spv's is at least offset + len long
 *  - the data specified by (offset, len) is correctly normalised in the
 *    range [0, sp)
 *
 * return 1 for success, 0 for failure */

int
mpzspv_verify (mpzspv_t x, spv_size_t offset, spv_size_t len, mpzspm_t mpzspm)
{
  unsigned int i;
  spv_size_t j;
  
#ifdef HAVE_malloc_usable_size
  if (malloc_usable_size (x) < mpzspm->sp_num * sizeof (spv_t))
    return 0;
#endif

  for (i = 0; i < mpzspm->sp_num; i++)
    {

#ifdef HAVE_malloc_usable_size
      if (malloc_usable_size (x[i]) < (offset + len) * sizeof (sp_t))
        return 0;
#endif

      for (j = offset; j < offset + len; j++)
	if (x[i][j] >= mpzspm->spm[i].sp)
	  return 0;
    }

  return 1;
}

void
mpzspv_set (mpzspv_t r, spv_size_t r_offset, mpzspv_t x, spv_size_t x_offset,
    spv_size_t len, mpzspm_t mpzspm)
{
  unsigned int i;
  
  ASSERT (mpzspv_verify (r, r_offset + len, 0, mpzspm));
  ASSERT (mpzspv_verify (x, x_offset, len, mpzspm));
  
  for (i = 0; i < mpzspm->sp_num; i++)
    spv_set (r[i] + r_offset, x[i] + x_offset, len);
}

void
mpzspv_set_sp (mpzspv_t r, spv_size_t offset, sp_t c, spv_size_t len,
    mpzspm_t mpzspm)
{
  unsigned int i;
  
  ASSERT (mpzspv_verify (r, offset + len, 0, mpzspm));
  ASSERT (c < SP_MIN); /* not strictly necessary but avoids mod functions */
  
  for (i = 0; i < mpzspm->sp_num; i++)
    spv_set_sp (r[i] + offset, c, len);
}

void
mpzspv_neg (mpzspv_t r, spv_size_t r_offset, mpzspv_t x, spv_size_t x_offset,
    spv_size_t len, mpzspm_t mpzspm)
{
  unsigned int i;
  
  ASSERT (mpzspv_verify (r, r_offset + len, 0, mpzspm));
  ASSERT (mpzspv_verify (x, x_offset, len, mpzspm));
  
  for (i = 0; i < mpzspm->sp_num; i++)
    spv_neg (r[i] + r_offset, x[i] + x_offset, len, mpzspm->spm[i].sp);
}

void
mpzspv_reverse (mpzspv_t x, spv_size_t offset, spv_size_t len, mpzspm_t mpzspm)
{
  unsigned int i;
  spv_size_t j;
  sp_t t;
  spv_t spv;
  
  ASSERT (mpzspv_verify (x, offset, len, mpzspm));
  
  for (i = 0; i < mpzspm->sp_num; i++)
    {
      spv = x[i] + offset;
      for (j = 0; j < len - 1 - j; j++)
        {
	  t = spv[j];
	  spv[j] = spv[len - 1 - j];
	  spv[len - 1 - j] = t;
	}
    }
}

void
mpzspv_from_mpzv (mpzspv_t x, spv_size_t offset, mpzv_t mpzv,
    spv_size_t len, mpzspm_t mpzspm)
{
  unsigned int i, sp_num;
  spv_size_t j;
  
  ASSERT (mpzspv_verify (x, offset + len, 0, mpzspm));
  
  sp_num = mpzspm->sp_num;
  
  /* GMP's comments on mpn_preinv_mod_1:
   *
   * "This function used to be documented, but is now considered obsolete.  It
   * continues to exist for binary compatibility, even when not required
   * internally."
   *
   * It doesn't accept 0 as the dividend so we have to treat this case
   * separately */
  
  for (i = 0; i < len; i++)
    {
      if (SIZ(mpzv[i]) == 0)
	{
	  for (j = 0; j < sp_num; j++)
	    x[j][i + offset] = 0;
	}
      else
        {
	  for (j = 0; j < sp_num; j++)
            x[j][i + offset] =
              mpn_preinv_mod_1 (PTR(mpzv[i]), SIZ(mpzv[i]),
                mpzspm->spm[j].sp, mpzspm->spm[j].mul_c);
	}
    }
}

/* B&S: ecrt mod m */
void
mpzspv_to_mpzv (mpzspv_t x, spv_size_t offset, mpzv_t mpzv,
    spv_size_t len, mpzspm_t mpzspm)
{
  unsigned int i;
  spv_size_t j;
  float *f = (float *) valloc (len * sizeof (float));
  float prime_recip;
  sp_t t;
  spm_t spm;
  
  ASSERT (mpzspv_verify (x, offset, len, mpzspm));
  
  for (j = 0; j < len; j++)
    {
      f[j] = 0.5;
      mpz_set_ui (mpzv[j], 0);
    }
  
  for (i = 0; i < mpzspm->sp_num; i++)
    {
      prime_recip = 1.0 / (float) mpzspm->spm[i].sp;
      spm = mpzspm->spm + i;
      
      for (j = 0; j < len; j++)
        {
	  t = sp_mul (x[i][j + offset], mpzspm->crt3[i], spm->sp, spm->mul_c);
          
	  mpz_addmul_ui (mpzv[j], mpzspm->crt1[i], t);

	  f[j] += (float) t * prime_recip;
        }
    }

  for (j = 0; j < len; j++)
    mpz_add (mpzv[j], mpzv[j], mpzspm->crt2[(unsigned int) f[j]]);
  
  free (f);
}  

void
mpzspv_pwmul (mpzspv_t r, spv_size_t r_offset, mpzspv_t x, spv_size_t x_offset,
    mpzspv_t y, spv_size_t y_offset, spv_size_t len, mpzspm_t mpzspm)
{
  unsigned int i;
  
  ASSERT (mpzspv_verify (r, r_offset + len, 0, mpzspm));
  ASSERT (mpzspv_verify (x, x_offset, len, mpzspm));
  ASSERT (mpzspv_verify (y, y_offset, len, mpzspm));
  
  for (i = 0; i < mpzspm->sp_num; i++)
    spv_pwmul (r[i] + r_offset, x[i] + x_offset, y[i] + y_offset,
	len, mpzspm->spm[i].sp, mpzspm->spm[i].mul_c);
}

#define STRIDE MIN(256,len)

/* B&S: ecrt mod m mod p_j. */
void
mpzspv_normalise (mpzspv_t x, spv_size_t offset, spv_size_t len,
    mpzspm_t mpzspm)
{
  unsigned long i, j, k, l;
  unsigned int sp_num = mpzspm->sp_num;
  sp_t v;
  spv_t s, d, w;
  spm_t spm = mpzspm->spm;
  
  float prime_recip;
  float *f;
  mpzspv_t t;
  
  ASSERT (mpzspv_verify (x, offset, len, mpzspm)); 
  
  f = (float *) valloc (len * sizeof (float));
  t = mpzspv_init (STRIDE, mpzspm);
  
  s = (spv_t) valloc (3 * len * sizeof (sp_t));
  d = (spv_t) valloc (3 * len * sizeof (sp_t));
  memset (s, 0, 3 * len * sizeof (sp_t));

  for (i = 0; i < len; i++)
    f[i] = 0.5;
  
  /* FIXME: use B&S Theorem 2.2 */
  for (i = 0; i < sp_num; i++)
    {
      prime_recip = 1.0 / (float) spm[i].sp;
      
      for (j = 0; j < len; j++)
	{
	  x[i][j + offset] = sp_mul (x[i][j + offset],
	      mpzspm->crt3[i], spm[i].sp, spm[i].mul_c);
	  f[j] += (float) x[i][j + offset] * prime_recip;
	}
    }
  
  for (l = 0; l < len; l += STRIDE)
    {
      for (i = 0; i < sp_num; i++)
        {
          for (j = 0; j < STRIDE; j++)
	    {
	      umul_ppmm (d[3 * j + 1], d[3 * j], mpzspm->crt5[i],
		  (unsigned int) f[j + l]);
              d[3 * j + 2] = 0;
	    }
	
          for (j = 0; j < sp_num; j++)
            {
	      w = x[j] + offset;
	      v = mpzspm->crt4[i][j];
	    
	      for (k = 0; k < STRIDE; k++)
	        umul_ppmm (s[3 * k + 1], s[3 * k], w[k + l], v);
 	      
	      /* this mpn_add_n accounts for about a third of the function's
	       * runtime */
	      mpn_add_n (d, d, s, 3 * STRIDE);
            }      

          /* FIXME: do we need to account for dividend == 0? */
          for (j = 0; j < STRIDE; j++)
	    t[i][j] = mpn_preinv_mod_1 (d + 3 * j, 3, spm[i].sp, spm[i].mul_c);
        }	  
      mpzspv_set (x, l + offset, t, 0, STRIDE, mpzspm);
    }
  
  mpzspv_clear (t, mpzspm);
  
  free (s);
  free (d);
  free (f);
}

void
mpzspv_to_ntt (mpzspv_t x, spv_size_t offset, spv_size_t len,
    spv_size_t ntt_size, int monic, mpzspm_t mpzspm)
{
  unsigned int i;
  spv_size_t j;
  spm_t spm;
  sp_t root;
  spv_t spv;
  
  ASSERT (mpzspv_verify (x, offset, len, mpzspm));
  ASSERT (mpzspv_verify (x, offset + ntt_size, 0, mpzspm));
  
  for (i = 0; i < mpzspm->sp_num; i++)
    {
      spm = mpzspm->spm + i;
      root = sp_pow (spm->prim_root, mpzspm->max_ntt_size / ntt_size,
	  spm->sp, spm->mul_c);
      
      spv = x[i] + offset;
      
      if (ntt_size < len)
        {
	  for (j = ntt_size; j < len; j += ntt_size)
	    spv_add (spv, spv, spv + j, ntt_size, spm->sp);
	}
      if (ntt_size > len)
	spv_set_zero (spv + len, ntt_size - len);

      if (monic)
	spv[len % ntt_size] = sp_add (spv[len % ntt_size], 1, spm->sp);
      
      spv_ntt_gfp_dif (spv, ntt_size, spm->sp, spm->mul_c, root);
    }
}

void mpzspv_from_ntt (mpzspv_t x, spv_size_t offset, spv_size_t ntt_size,
    spv_size_t monic_pos, mpzspm_t mpzspm)
{
  unsigned int i;
  spm_t spm;
  sp_t root;
  spv_t spv;
  
  ASSERT (mpzspv_verify (x, offset, ntt_size, mpzspm));
  
  for (i = 0; i < mpzspm->sp_num; i++)
    {
      spm = mpzspm->spm + i;
      root = sp_pow (spm->inv_prim_root, mpzspm->max_ntt_size / ntt_size,
	  spm->sp, spm->mul_c);
      
      spv = x[i] + offset;
      
      spv_ntt_gfp_dit (spv, ntt_size, spm->sp, spm->mul_c, root);

      /* spm->sp - (spm->sp - 1) / ntt_size is the inverse of ntt_size */
      spv_mul_sp (spv, spv, spm->sp - (spm->sp - 1) / ntt_size,
	  ntt_size, spm->sp, spm->mul_c);
      
      if (monic_pos)
	spv[monic_pos % ntt_size] = sp_sub (spv[monic_pos % ntt_size],
	    1, spm->sp);
    }
}

