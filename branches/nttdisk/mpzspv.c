/* mpzspv.c - "mpz small prime polynomial" functions for arithmetic on mpzv's
   reduced modulo a mpzspm

  Copyright 2005, 2008, 2010 Dave Newman, Jason Papadopoulos, Alexander Kruppa
                             and Paul Zimmermann.

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

#include <stdio.h> /* for stderr */
#include <stdlib.h>
#include <string.h> /* for memset */
#include "ecm-impl.h"

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
  
  for (i = 0; i < mpzspm->sp_num; i++)
    {
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

/* convert mpzvi to CRT representation, naive version */
static void
mpzspv_from_mpzv_slow (mpzspv_t x, const spv_size_t offset, mpz_t mpzvi,
                       mpzspm_t mpzspm, 
		       ATTRIBUTE_UNUSED mpz_t rem,
		       unsigned int sp_num)
{
  unsigned int j;

  for (j = 0; j < sp_num; j++)
    { 
#if SP_TYPE_BITS > GMP_LIMB_BITS
      mpz_tdiv_r(rem, mpzvi, mpzspm->spm[j]->mp_sp);
      x[j][offset] = mpz_get_sp(rem);
#else
      x[j][offset] = mpn_mod_1 (PTR(mpzvi), SIZ(mpzvi),
                              (mp_limb_t) mpzspm->spm[j]->sp);
#endif
    }
}

/* convert mpzvi to CRT representation, fast version, assumes
   mpzspm->T has been precomputed (see mpzspm.c) */
static void
mpzspv_from_mpzv_fast (mpzspv_t x, const spv_size_t offset, mpz_t mpzvi,
                       mpzspm_t mpzspm, 
		       ATTRIBUTE_UNUSED mpz_t rem,
		       unsigned int sp_num)
{
  unsigned int i, j, k, i0 = I0_THRESHOLD, I0;
  mpzv_t *T = mpzspm->T;
  unsigned int d = mpzspm->d, ni;

  ASSERT (d > i0);

  /* T[0] serves as vector of temporary mpz_t's, since it contains the small
     primes, which are also in mpzspm->spm[j]->sp */
  /* initially we split mpzvi in two */
  ni = 1 << (d - 1);
  mpz_mod (T[0][0], mpzvi, T[d-1][0]);
  mpz_mod (T[0][ni], mpzvi, T[d-1][1]);
  for (i = d-1; i-- > i0;)
    { /* goes down from depth i+1 to i */
      ni = 1 << i;
      for (j = k = 0; j + ni < sp_num; j += 2*ni, k += 2)
        {
          mpz_mod (T[0][j+ni], T[0][j], T[i][k+1]);
          mpz_mod (T[0][j], T[0][j], T[i][k]);
        }
      /* for the last entry T[0][j] if j < sp_num, there is nothing to do */
    }
  /* last steps */
  I0 = 1 << i0;
  for (j = 0; j < sp_num; j += I0)
    for (k = j; k < j + I0 && k < sp_num; k++)
      {
#if SP_TYPE_BITS > GMP_LIMB_BITS
	mpz_tdiv_r(rem, T[0][j], mpzspm->spm[k]->mp_sp);
	x[k][offset] = mpz_get_sp(rem);
#else
	x[k][offset] = mpn_mod_1 (PTR(T[0][j]), SIZ(T[0][j]),
                                (mp_limb_t) mpzspm->spm[k]->sp);
#endif
      }
}

/* convert an array of len mpz_t numbers to CRT representation modulo
   sp_num moduli */
void
mpzspv_from_mpzv (mpzspv_t x, const spv_size_t offset, const mpzv_t mpzv,
    const spv_size_t len, mpzspm_t mpzspm)
{
  const unsigned int sp_num = mpzspm->sp_num;
  mpz_t rem;
  long i;

  ASSERT (mpzspv_verify (x, offset + len, 0, mpzspm));

  /* GMP's comments on mpn_preinv_mod_1:
   *
   * "This function used to be documented, but is now considered obsolete.  It
   * continues to exist for binary compatibility, even when not required
   * internally."
   *
   * It doesn't accept 0 as the dividend so we have to treat this case
   * separately */
  
#if defined(_OPENMP)
#pragma omp parallel private(i,rem) if (len > 16384)
  {
    /* Multi-threading with dynamic scheduling slows things down */
#pragma omp for schedule(static)
#endif

    mpz_init(rem);
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
          if (mpzspm->T == NULL)
            mpzspv_from_mpzv_slow (x, i + offset, mpzv[i], mpzspm, rem, sp_num);
          else
            mpzspv_from_mpzv_fast (x, i + offset, mpzv[i], mpzspm, rem, sp_num);
	}
    }
    mpz_clear(rem);

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

  if (f == NULL)
    {
      fprintf (stderr, "Cannot allocate memory in mpzspv_to_mpzv\n");
      exit (1);
    }
  
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
          
#if SP_TYPE_BITS > GMP_LIMB_BITS
            mpz_set_sp (mt, t);
            mpz_addmul (mpzv[l + k], mpzspm->crt1[i], mt);
#else
	    mpz_addmul_ui (mpzv[l + k], mpzspm->crt1[i], t);
#endif
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
  s = (spv_t) malloc (3 * MPZSPV_NORMALISE_STRIDE * sizeof (sp_t));
  d = (spv_t) malloc (3 * MPZSPV_NORMALISE_STRIDE * sizeof (sp_t));
  if (f == NULL || s == NULL || d == NULL)
    {
      fprintf (stderr, "Cannot allocate memory in mpzspv_normalise\n");
      exit (1);
    }
  t = mpzspv_init (MPZSPV_NORMALISE_STRIDE, mpzspm);
  
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
	      sp_wide_mul (d[3 * k + 1], d[3 * k], mpzspm->crt5[i],
		  (sp_t) f[k]);
              d[3 * k + 2] = 0;
	    }
	
          for (j = 0; j < sp_num; j++)
            {
	      w = x[j] + offset;
	      v = mpzspm->crt4[i][j];
	    
	      for (k = 0; k < stride; k++)
	        sp_wide_mul (s[3 * k + 1], s[3 * k], w[k + l], v);
 	      
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


/* Do multiplication via NTT. Depending on the value of "steps", does 
   in-place forward transform of x, in-place forward transform of y, 
   pair-wise multiplication of x by y to r, in-place inverse transform of r. 
   Contrary to calling these three operations separately, this function does 
   all three steps on a small-prime vector at a time, resulting in slightly 
   better cache efficiency (also in preparation to storing NTT vectors on disk 
   and reading them in for the multiplication). */

void
mpzspv_mul_ntt (mpzspv_t r, const spv_size_t offsetr, 
    mpzspv_t x, const spv_size_t offsetx, const spv_size_t lenx,
    mpzspv_t y, const spv_size_t offsety, const spv_size_t leny,
    const spv_size_t ntt_size, const int monic, const spv_size_t monic_pos, 
    mpzspm_t mpzspm, const int steps)
{
  spv_size_t log2_ntt_size;
  int i;
  
  ASSERT (mpzspv_verify (x, offsetx, lenx, mpzspm));
  ASSERT (mpzspv_verify (y, offsety, leny, mpzspm));
  ASSERT (mpzspv_verify (x, offsetx + ntt_size, 0, mpzspm));
  ASSERT (mpzspv_verify (y, offsety + ntt_size, 0, mpzspm));
  ASSERT (mpzspv_verify (r, offsetr + ntt_size, 0, mpzspm));
  
  log2_ntt_size = ceil_log_2 (ntt_size);

  /* Need parallelization at higher level (e.g., handling a branch of the 
     product tree in one thread) to make this worthwhile for ECM */
#define MPZSPV_MUL_NTT_OPENMP 0

#if defined(_OPENMP) && MPZSPV_MUL_NTT_OPENMP
#pragma omp parallel if (ntt_size > 16384)
  {
#pragma omp for
#endif
  for (i = 0; i < (int) mpzspm->sp_num; i++)
    {
      spv_size_t j;
      spm_t spm = mpzspm->spm[i];
      spv_t spvr = r[i] + offsetr;
      spv_t spvx = x[i] + offsetx;
      spv_t spvy = y[i] + offsety;

      if ((steps & NTT_MUL_STEP_FFT1) != 0) {
        if (ntt_size < lenx)
          {
            for (j = ntt_size; j < lenx; j += ntt_size)
              spv_add (spvx, spvx, spvx + j, ntt_size, spm->sp);
          }
        if (ntt_size > lenx)
          spv_set_zero (spvx + lenx, ntt_size - lenx);

        if (monic)
          spvx[lenx % ntt_size] = sp_add (spvx[lenx % ntt_size], 1, spm->sp);

        spv_ntt_gfp_dif (spvx, log2_ntt_size, spm);
      }

      if ((steps & NTT_MUL_STEP_FFT2) != 0) {
        if (ntt_size < leny)
          {
            for (j = ntt_size; j < leny; j += ntt_size)
              spv_add (spvy, spvy, spvy + j, ntt_size, spm->sp);
          }
        if (ntt_size > leny)
          spv_set_zero (spvy + leny, ntt_size - leny);

        if (monic)
          spvy[leny % ntt_size] = sp_add (spvy[leny % ntt_size], 1, spm->sp);

        spv_ntt_gfp_dif (spvy, log2_ntt_size, spm);
      }

      if ((steps & NTT_MUL_STEP_MUL) != 0) {
        spv_pwmul (spvr, spvx, spvy, ntt_size, spm->sp, spm->mul_c);
      }

      if ((steps & NTT_MUL_STEP_IFFT) != 0) {
        ASSERT (sizeof (mp_limb_t) >= sizeof (sp_t));

        spv_ntt_gfp_dit (spvr, log2_ntt_size, spm);

        /* spm->sp - (spm->sp - 1) / ntt_size is the inverse of ntt_size */
        spv_mul_sp (spvr, spvr, spm->sp - (spm->sp - 1) / ntt_size,
            ntt_size, spm->sp, spm->mul_c);

        if (monic_pos)
          spvr[monic_pos % ntt_size] = sp_sub (spvr[monic_pos % ntt_size],
              1, spm->sp);
      }
    }
#if defined(_OPENMP) && MPZSPV_MUL_NTT_OPENMP
  }
#endif
}

/* Computes a DCT-I of the length dctlen. Input is the spvlen coefficients
   in spv. tmp is temp space and must have space for 2*dctlen-2 sp_t's */

void
mpzspv_to_dct1 (mpzspv_t dct, const mpzspv_t spv, const spv_size_t spvlen, 
                const spv_size_t dctlen, mpzspv_t tmp, 
		const mpzspm_t mpzspm)
{
  const spv_size_t l = 2 * (dctlen - 1); /* Length for the DFT */
  const spv_size_t log2_l = ceil_log_2 (l);
  int j;

#ifdef _OPENMP
#pragma omp parallel private(j)
  {
#pragma omp for
#endif
  for (j = 0; j < (int) mpzspm->sp_num; j++)
    {
      const spm_t spm = mpzspm->spm[j];
      spv_size_t i;
      
      /* Make a symmetric copy of spv in tmp. I.e. with spv = [3, 2, 1], 
         spvlen = 3, dctlen = 5 (hence l = 8), we want 
         tmp = [3, 2, 1, 0, 0, 0, 1, 2] */
      spv_set (tmp[j], spv[j], spvlen);
      spv_rev (tmp[j] + l - spvlen + 1, spv[j] + 1, spvlen - 1);
      /* Now we have [3, 2, 1, ?, ?, ?, 1, 2]. Fill the ?'s with zeros. */
      spv_set_sp (tmp[j] + spvlen, (sp_t) 0, l - 2 * spvlen + 1);

#if 0
      printf ("mpzspv_to_dct1: tmp[%d] = [", j);
      for (i = 0; i < l; i++)
          printf ("%lu, ", tmp[j][i]);
      printf ("]\n");
#endif
      
      spv_ntt_gfp_dif (tmp[j], log2_l, spm);

#if 0
      printf ("mpzspv_to_dct1: tmp[%d] = [", j);
      for (i = 0; i < l; i++)
          printf ("%lu, ", tmp[j][i]);
      printf ("]\n");
#endif

      /* The forward transform is scrambled. We want elements [0 ... l/2]
         of the unscrabled data, that is all the coefficients with the most 
         significant bit in the index (in log2(l) word size) unset, plus the 
         element at index l/2. By scrambling, these map to the elements with 
         even index, plus the element at index 1. 
         The elements with scrambled index 2*i are stored in h[i], the
         element with scrambled index 1 is stored in h[params->l] */
  
#ifdef WANT_ASSERT
      /* Test that the coefficients are symmetric (if they were unscrambled)
         and that our algorithm for finding identical coefficients in the 
         scrambled data works */
      {
        spv_size_t m = 5;
        for (i = 2; i < l; i += 2L)
          {
            /* This works, but why? */
            if (i + i / 2L > m)
                m = 2L * m + 1L;

            ASSERT (tmp[j][i] == tmp[j][m - i]);
#if 0
            printf ("mpzspv_to_dct1: DFT[%lu] == DFT[%lu]\n", i, m - i);
#endif
          }
      }
#endif

      /* Copy coefficients to dct buffer */
      for (i = 0; i < l / 2; i++)
        dct[j][i] = tmp[j][i * 2];
      dct[j][l / 2] = tmp[j][1];
    }
#ifdef _OPENMP
  }
#endif
}


/* Multiply the polynomial in "dft" by the RLP in "dct", where "dft" 
   contains the polynomial coefficients (not FFT'd yet) and "dct" 
   contains the DCT-I coefficients of the RLP. The latter are 
   assumed to be in the layout produced by mpzspv_to_dct1().
   Output are the coefficients of the product polynomial, stored in dft. 
   The "steps" parameter controls which steps are computed:
   NTT_MUL_STEP_FFT1: do forward transform
   NTT_MUL_STEP_MUL: do point-wise product
   NTT_MUL_STEP_IFFT: do inverse transform 
*/

void
mpzspv_mul_by_dct (mpzspv_t dft, const mpzspv_t dct, const spv_size_t len, 
		   const mpzspm_t mpzspm, const int steps)
{
  int j;
  spv_size_t log2_len = ceil_log_2 (len);
  
#ifdef _OPENMP
#pragma omp parallel private(j)
  {
#pragma omp for
#endif
    for (j = 0; j < (int) (mpzspm->sp_num); j++)
      {
	const spm_t spm = mpzspm->spm[j];
	const spv_t spv = dft[j];
	unsigned long i, m;
	
	/* Forward DFT of dft[j] */
	if ((steps & NTT_MUL_STEP_FFT1) != 0)
	  spv_ntt_gfp_dif (spv, log2_len, spm);
	
	/* Point-wise product */
	if ((steps & NTT_MUL_STEP_MUL) != 0)
	  {
	    m = 5UL;
	    
	    spv[0] = sp_mul (spv[0], dct[j][0], spm->sp, spm->mul_c);
	    spv[1] = sp_mul (spv[1], dct[j][len / 2UL], spm->sp, spm->mul_c);
	    
	    for (i = 2UL; i < len; i += 2UL)
	      {
		/* This works, but why? */
		if (i + i / 2UL > m)
		  m = 2UL * m + 1;
		
		spv[i] = sp_mul (spv[i], dct[j][i / 2UL], spm->sp, spm->mul_c);
		spv[m - i] = sp_mul (spv[m - i], dct[j][i / 2UL], spm->sp, 
				     spm->mul_c);
	      }
	  }
	
	/* Inverse transform of dft[j] */
	if ((steps & NTT_MUL_STEP_IFFT) != 0)
	  {
	    spv_ntt_gfp_dit (spv, log2_len, spm);
	    
	    /* Divide by transform length. FIXME: scale the DCT of h instead */
	    spv_mul_sp (spv, spv, spm->sp - (spm->sp - 1) / len, len, 
			spm->sp, spm->mul_c);
	  }
      }
#ifdef _OPENMP
  }
#endif
}


void 
mpzspv_sqr_reciprocal (mpzspv_t dft, const spv_size_t n, 
                       const mpzspm_t mpzspm)
{
  const spv_size_t log2_n = ceil_log_2 (n);
  const spv_size_t len = ((spv_size_t) 2) << log2_n;
  const spv_size_t log2_len = 1 + log2_n;
  int j;

  ASSERT(mpzspm->max_ntt_size % 3UL == 0UL);
  ASSERT(len % 3UL != 0UL);
  ASSERT(mpzspm->max_ntt_size % len == 0UL);

#ifdef _OPENMP
#pragma omp parallel
  {
#pragma omp for
#endif
    for (j = 0; j < (int) (mpzspm->sp_num); j++)
      {
        const spm_t spm = mpzspm->spm[j];
        const spv_t spv = dft[j];
        sp_t w1, w2, invlen;
        const sp_t sp = spm->sp, mul_c = spm->mul_c;
        spv_size_t i;

        /* Zero out NTT elements [n .. len-n] */
        spv_set_sp (spv + n, (sp_t) 0, len - 2*n + 1);

#ifdef TRACE_ntt_sqr_reciprocal
        if (j == 0)
          {
            printf ("ntt_sqr_reciprocal: NTT vector mod %lu\n", sp);
            ntt_print_vec ("ntt_sqr_reciprocal: before weighting:", spv, len);
          }
#endif

        /* Compute the root for the weight signal, a 3rd primitive root 
           of unity */
        w1 = sp_pow (spm->prim_root, mpzspm->max_ntt_size / 3UL, sp, 
                     mul_c);
        /* Compute iw= 1/w */
        w2 = sp_pow (spm->inv_prim_root, mpzspm->max_ntt_size / 3UL, sp, 
                     mul_c);
#ifdef TRACE_ntt_sqr_reciprocal
        if (j == 0)
          printf ("w1 = %lu ,w2 = %lu\n", w1, w2);
#endif
        ASSERT(sp_mul(w1, w2, sp, mul_c) == (sp_t) 1);
        ASSERT(w1 != (sp_t) 1);
        ASSERT(sp_pow (w1, 3UL, sp, mul_c) == (sp_t) 1);
        ASSERT(w2 != (sp_t) 1);
        ASSERT(sp_pow (w2, 3UL, sp, mul_c) == (sp_t) 1);

        /* Fill NTT elements spv[len-n+1 .. len-1] with coefficients and
           apply weight signal to spv[i] and spv[l-i] for 0 <= i < n
           Use the fact that w^i + w^{-i} = -1 if i != 0 (mod 3). */
        for (i = 0; i + 2 < n; i += 3)
          {
            sp_t t, u;
            
	    if (i > 0)
	      spv[len - i] = spv[i];
            
            t = spv[i + 1];
            u = sp_mul (t, w1, sp, mul_c);
            spv[i + 1] = u;
            spv[len - i - 1] = sp_neg (sp_add (t, u, sp), sp);

            t = spv[i + 2];
            u = sp_mul (t, w2, sp, mul_c);
            spv[i + 2] = u;
            spv[len - i - 2] = sp_neg (sp_add (t, u, sp), sp);
          }
        if (i < n && i > 0)
          {
            spv[len - i] = spv[i];
          }
        if (i + 1 < n)
          {
            sp_t t, u;
            t = spv[i + 1];
            u = sp_mul (t, w1, sp, mul_c);
            spv[i + 1] = u;
            spv[len - i - 1] = sp_neg (sp_add (t, u, sp), sp);
          }

#ifdef TRACE_ntt_sqr_reciprocal
      if (j == 0)
        ntt_print_vec ("ntt_sqr_reciprocal: after weighting:", spv, len);
#endif

        /* Forward DFT of dft[j] */
        spv_ntt_gfp_dif (spv, log2_len, spm);

#ifdef TRACE_ntt_sqr_reciprocal
        if (j == 0)
          ntt_print_vec ("ntt_sqr_reciprocal: after forward transform:", 
			 spv, len);
#endif

        /* Square the transformed vector point-wise */
        spv_pwmul (spv, spv, spv, len, sp, mul_c);
      
#ifdef TRACE_ntt_sqr_reciprocal
        if (j == 0)
          ntt_print_vec ("ntt_sqr_reciprocal: after point-wise squaring:", 
			 spv, len);
#endif

        /* Inverse transform of dft[j] */
        spv_ntt_gfp_dit (spv, log2_len, spm);
      
#ifdef TRACE_ntt_sqr_reciprocal
        if (j == 0)
          ntt_print_vec ("ntt_sqr_reciprocal: after inverse transform:", 
			 spv, len);
#endif

        /* Un-weight and divide by transform length */
        invlen = sp - (sp - (sp_t) 1) / len; /* invlen = 1/len (mod sp) */
        w1 = sp_mul (invlen, w1, sp, mul_c);
        w2 = sp_mul (invlen, w2, sp, mul_c);
        for (i = 0; i < 2 * n - 3; i += 3)
          {
            spv[i] = sp_mul (spv[i], invlen, sp, mul_c);
            spv[i + 1] = sp_mul (spv[i + 1], w2, sp, mul_c);
            spv[i + 2] = sp_mul (spv[i + 2], w1, sp, mul_c);
          }
        if (i < 2 * n - 1)
          spv[i] = sp_mul (spv[i], invlen, sp, mul_c);
        if (i < 2 * n - 2)
          spv[i + 1] = sp_mul (spv[i + 1], w2, sp, mul_c);
        
#ifdef TRACE_ntt_sqr_reciprocal
        if (j == 0)
          ntt_print_vec ("ntt_sqr_reciprocal: after un-weighting:", spv, len);
#endif

        /* Separate the coefficients of R in the wrapped-around product. */

        /* Set w1 = cuberoot(1)^l where cuberoot(1) is the same primitive
           3rd root of unity we used for the weight signal */
        w1 = sp_pow (spm->prim_root, mpzspm->max_ntt_size / 3UL, sp, 
                     mul_c);
        w1 = sp_pow (w1, len % 3UL, sp, mul_c);
        
        /* Set w2 = 1/(w1 - 1/w1). Incidentally, w2 = 1/sqrt(-3) */
        w2 = sp_inv (w1, sp, mul_c);
        w2 = sp_sub (w1, w2, sp);
        w2 = sp_inv (w2, sp, mul_c);
#ifdef TRACE_ntt_sqr_reciprocal
        if (j == 0)
          printf ("For separating: w1 = %lu, w2 = %lu\n", w1, w2);
#endif
        
        for (i = len - (2*n - 2); i <= len / 2; i++)
          {
            sp_t t, u;
            /* spv[i] = s_i + w^{-l} s_{l-i}. 
               spv[l-i] = s_{l-i} + w^{-l} s_i */
            t = sp_mul (spv[i], w1, sp, mul_c); /* t = w^l s_i + s_{l-i} */
            t = sp_sub (t, spv[len - i], sp);   /* t = w^l s_i + w^{-l} s_i */
            t = sp_mul (t, w2, sp, mul_c);      /* t = s_1 */

            u = sp_sub (spv[i], t, sp);         /* u = w^{-l} s_{l-i} */
            u = sp_mul (u, w1, sp, mul_c);      /* u = s_{l-i} */
            spv[i] = t;
            spv[len - i] = u;
            ASSERT(i < len / 2 || t == u);
          }

#ifdef TRACE_ntt_sqr_reciprocal
        if (j == 0)
          ntt_print_vec ("ntt_sqr_reciprocal: after un-wrapping:", spv, len);
#endif
      }
#ifdef _OPENMP
    }
#endif
}
