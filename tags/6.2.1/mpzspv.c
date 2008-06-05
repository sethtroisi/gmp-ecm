/* mpzspv.c - "mpz small prime polynomial" functions for arithmetic on mpzv's
   reduced modulo a mpzspm

  Copyright 2005, 2008 Dave Newman, Jason Papadopoulos and Alexander Kruppa.

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
#include <string.h> /* for memset */
#include "sp.h"

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

mpzspv_t
mpzspv_init (spv_size_t len, mpzspm_t mpzspm)
{
  unsigned int i;
  mpzspv_t x = (mpzspv_t) malloc (mpzspm->sp_num * sizeof (spv_t));
  
  if (x == NULL)
    return NULL;
  
  for (i = 0; i < mpzspm->sp_num; i++)
    {
      x[i] = (spv_t) sp_aligned_malloc (len * sizeof (sp_t));
      
      if (x[i] == NULL)
	{
	  while (i--)
	    sp_aligned_free (x[i]);
	  
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
    sp_aligned_free (x[i]);
  
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
  
#ifdef HAVE_MALLOC_USABLE_SIZE
  if (malloc_usable_size (x) < mpzspm->sp_num * sizeof (spv_t))
    return 0;
#endif

  for (i = 0; i < mpzspm->sp_num; i++)
    {

#ifdef HAVE_MALLOC_USABLE_SIZE
      if (malloc_usable_size (*( (void **)x[i] - 1 )) < 
          (offset + len) * sizeof (sp_t))
        return 0;
#endif

      for (j = offset; j < offset + len; j++)
	if (x[i][j] >= mpzspm->spm[i]->sp)
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
mpzspv_revcopy (mpzspv_t r, spv_size_t r_offset, mpzspv_t x, 
    spv_size_t x_offset, spv_size_t len, mpzspm_t mpzspm)
{
  unsigned int i;
  
  ASSERT (mpzspv_verify (r, r_offset + len, 0, mpzspm));
  ASSERT (mpzspv_verify (x, x_offset, len, mpzspm));
  
  for (i = 0; i < mpzspm->sp_num; i++)
    spv_rev (r[i] + r_offset, x[i] + x_offset, len);
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
    spv_neg (r[i] + r_offset, x[i] + x_offset, len, mpzspm->spm[i]->sp);
}

void
mpzspv_add (mpzspv_t r, spv_size_t r_offset, mpzspv_t x, spv_size_t x_offset,
            mpzspv_t y, spv_size_t y_offset, spv_size_t len, mpzspm_t mpzspm)
{
  unsigned int i;
  
  ASSERT (mpzspv_verify (r, r_offset + len, 0, mpzspm));
  ASSERT (mpzspv_verify (x, x_offset, len, mpzspm));
  
  for (i = 0; i < mpzspm->sp_num; i++)
    spv_add (r[i] + r_offset, x[i] + x_offset, y[i] + y_offset, len, 
             mpzspm->spm[i]->sp);
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
mpzspv_from_mpzv (mpzspv_t x, const spv_size_t offset, const mpzv_t mpzv,
    const spv_size_t len, mpzspm_t mpzspm)
{
  const unsigned int sp_num = mpzspm->sp_num;
  long i;
  
  ASSERT (mpzspv_verify (x, offset + len, 0, mpzspm));
  ASSERT (sizeof (mp_limb_t) >= sizeof (sp_t));
  
  /* GMP's comments on mpn_preinv_mod_1:
   *
   * "This function used to be documented, but is now considered obsolete.  It
   * continues to exist for binary compatibility, even when not required
   * internally."
   *
   * It doesn't accept 0 as the dividend so we have to treat this case
   * separately */
  
#if defined(_OPENMP)
#pragma omp parallel private(i) if (len > 100)
  {
    /* Multi-threading with dynamic scheduling slows things down */
#pragma omp for schedule(static)
#endif
    for (i = 0; i < (long) len; i++)
    {
      unsigned int j;
      if (mpz_sgn (mpzv[i]) == 0)
	{
	  for (j = 0; j < sp_num; j++)
	    x[j][i + offset] = 0;
	}
      else
        {
	  ASSERT(mpz_sgn (mpzv[i]) > 0); /* We can't handle negative values */
	  for (j = 0; j < sp_num; j++)
            x[j][i + offset] = 
              mpn_mod_1 (PTR(mpzv[i]), SIZ(mpzv[i]), 
                (mp_limb_t) mpzspm->spm[j]->sp);
              /* The typecast to mp_limb_t assumes that mp_limb_t is at least
                 as wide as sp_t */
	}
    }
#if defined(_OPENMP)
  }
#endif
}

/* See: Daniel J. Bernstein and Jonathan P. Sorenson,
 * Modular Exponentiation via the explicit Chinese Remainder Theorem
 *
 * memory: MPZSPV_NORMALISE_STRIDE floats */
void
mpzspv_to_mpzv (mpzspv_t x, spv_size_t offset, mpzv_t mpzv,
    spv_size_t len, mpzspm_t mpzspm)
{
  unsigned int i;
  spv_size_t k, l;
  float *f = (float *) malloc (MPZSPV_NORMALISE_STRIDE * sizeof (float));
  float prime_recip;
  sp_t t;
  spm_t *spm = mpzspm->spm;
  mpz_t mt;
  
  ASSERT (mpzspv_verify (x, offset, len, mpzspm));
  
  mpz_init (mt);
  for (l = 0; l < len; l += MPZSPV_NORMALISE_STRIDE)
    {
      spv_size_t stride = MIN (MPZSPV_NORMALISE_STRIDE, len - l);

      for (k = 0; k < stride; k++)
        {
          f[k] = 0.5;
          mpz_set_ui (mpzv[k + l], 0);
        }
  
    for (i = 0; i < mpzspm->sp_num; i++)
      {
        prime_recip = 1.0f / (float) spm[i]->sp;
      
        for (k = 0; k < stride; k++)
          {
  	    t = sp_mul (x[i][l + k + offset], mpzspm->crt3[i], spm[i]->sp,
                  spm[i]->mul_c);
          
            if (sizeof (sp_t) > sizeof (unsigned long))
              {
                mpz_set_sp (mt, t);
                mpz_addmul (mpzv[l + k], mpzspm->crt1[i], mt);
              }
            else
              {
      	        mpz_addmul_ui (mpzv[l + k], mpzspm->crt1[i], t);
              }

	    f[k] += (float) t * prime_recip;
          }
      }

    for (k = 0; k < stride; k++)
      mpz_add (mpzv[l + k], mpzv[l + k], mpzspm->crt2[(unsigned int) f[k]]);
  }
  
  mpz_clear (mt);
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
	len, mpzspm->spm[i]->sp, mpzspm->spm[i]->mul_c);
}

/* B&S: ecrt mod m mod p_j.
 *
 * memory: MPZSPV_NORMALISE_STRIDE mpzspv coeffs
 *         6 * MPZSPV_NORMALISE_STRIDE sp's
 *         MPZSPV_NORMALISE_STRIDE floats */
void
mpzspv_normalise (mpzspv_t x, spv_size_t offset, spv_size_t len,
    mpzspm_t mpzspm)
{
  unsigned int i, j, sp_num = mpzspm->sp_num;
  spv_size_t k, l;
  sp_t v;
  spv_t s, d, w;
  spm_t *spm = mpzspm->spm;
  
  float prime_recip;
  float *f;
  mpzspv_t t;
  
  ASSERT (mpzspv_verify (x, offset, len, mpzspm)); 
  
  f = (float *) malloc (MPZSPV_NORMALISE_STRIDE * sizeof (float));
  t = mpzspv_init (MPZSPV_NORMALISE_STRIDE, mpzspm);
  
  s = (spv_t) malloc (3 * MPZSPV_NORMALISE_STRIDE * sizeof (sp_t));
  d = (spv_t) malloc (3 * MPZSPV_NORMALISE_STRIDE * sizeof (sp_t));
  memset (s, 0, 3 * MPZSPV_NORMALISE_STRIDE * sizeof (sp_t));

  for (l = 0; l < len; l += MPZSPV_NORMALISE_STRIDE)
    {
      spv_size_t stride = MIN (MPZSPV_NORMALISE_STRIDE, len - l);
      
      /* FIXME: use B&S Theorem 2.2 */
      for (k = 0; k < stride; k++)
	f[k] = 0.5;
      
      for (i = 0; i < sp_num; i++)
        {
          prime_recip = 1.0f / (float) spm[i]->sp;
      
          for (k = 0; k < stride; k++)
	    {
	      x[i][l + k + offset] = sp_mul (x[i][l + k + offset],
	          mpzspm->crt3[i], spm[i]->sp, spm[i]->mul_c);
	      f[k] += (float) x[i][l + k + offset] * prime_recip;
	    }
        }
      
      for (i = 0; i < sp_num; i++)
        {
	  for (k = 0; k < stride; k++)
	    {
	      umul_ppmm (d[3 * k + 1], d[3 * k], mpzspm->crt5[i],
		  (sp_t) f[k]);
              d[3 * k + 2] = 0;
	    }
	
          for (j = 0; j < sp_num; j++)
            {
	      w = x[j] + offset;
	      v = mpzspm->crt4[i][j];
	    
	      for (k = 0; k < stride; k++)
	        umul_ppmm (s[3 * k + 1], s[3 * k], w[k + l], v);
 	      
	      /* this mpn_add_n accounts for about a third of the function's
	       * runtime */
	      mpn_add_n (d, d, s, 3 * stride);
            }      

          for (k = 0; k < stride; k++)
	    t[i][k] = mpn_mod_1 (d + 3 * k, 3, spm[i]->sp);
        }	  
      mpzspv_set (x, l + offset, t, 0, stride, mpzspm);
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
  spv_size_t j, log2_ntt_size;
  spm_t spm;
  spv_t spv;
  
  ASSERT (mpzspv_verify (x, offset, len, mpzspm));
  ASSERT (mpzspv_verify (x, offset + ntt_size, 0, mpzspm));
  
  log2_ntt_size = ceil_log_2 (ntt_size);

  for (i = 0; i < mpzspm->sp_num; i++)
    {
      spm = mpzspm->spm[i];
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
      
      spv_ntt_gfp_dif (spv, log2_ntt_size, spm);
    }
}

void mpzspv_from_ntt (mpzspv_t x, spv_size_t offset, spv_size_t ntt_size,
    spv_size_t monic_pos, mpzspm_t mpzspm)
{
  unsigned int i;
  spv_size_t log2_ntt_size;
  spm_t spm;
  spv_t spv;
  
  ASSERT (mpzspv_verify (x, offset, ntt_size, mpzspm));
  
  log2_ntt_size = ceil_log_2 (ntt_size);

  for (i = 0; i < mpzspm->sp_num; i++)
    {
      spm = mpzspm->spm[i];
      spv = x[i] + offset;
      
      spv_ntt_gfp_dit (spv, log2_ntt_size, spm);

      /* spm->sp - (spm->sp - 1) / ntt_size is the inverse of ntt_size */
      spv_mul_sp (spv, spv, spm->sp - (spm->sp - 1) / ntt_size,
	  ntt_size, spm->sp, spm->mul_c);
      
      if (monic_pos)
	spv[monic_pos % ntt_size] = sp_sub (spv[monic_pos % ntt_size],
	    1, spm->sp);
    }
}

void
mpzspv_random (mpzspv_t x, spv_size_t offset, spv_size_t len, mpzspm_t mpzspm)
{
  unsigned int i;

  ASSERT (mpzspv_verify (x, offset, len, mpzspm));

  for (i = 0; i < mpzspm->sp_num; i++)
    spv_random (x[i] + offset, len, mpzspm->spm[i]->sp);
}

