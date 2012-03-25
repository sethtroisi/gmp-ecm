/* mpzspv.c - "mpz small prime polynomial" functions for arithmetic on mpzv's
   reduced modulo a mpzspm

Copyright 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012 Dave Newman,
Jason Papadopoulos, Alexander Kruppa, Paul Zimmermann.

The SP Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The SP Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the SP Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
MA 02110-1301, USA. */

#define _GNU_SOURCE
#include "config.h"
#include <stdio.h> /* for stderr */
#include <stdlib.h>
#include <errno.h>
#include <string.h> /* for memset */
#ifdef USE_VALGRIND
#include <valgrind/memcheck.h>
#endif
#ifdef HAVE_FCNTL_H
#include <fcntl.h>
#endif
#include "ecm-impl.h"

#define MPZSPV_MUL_NTT_OPENMP 0
#define TRACE_ntt_sqr_reciprocal 0

#define IN_MEMORY(x) ((x) != NULL && (x)->storage == 0)
#define ON_DISK(x) ((x) != NULL && (x)->storage != 0)

static size_t seek_and_read_sp (spv_t, size_t, size_t, FILE *);
static size_t seek_and_write_sp (const spv_t, size_t, size_t, FILE *);
static void  mpzspv_seek_and_read (mpzspv_t, spv_size_t, FILE **, size_t, 
    size_t, mpzspm_t);
static void mpzspv_seek_and_write (mpzspv_t, spv_size_t, FILE **, size_t, 
    size_t, mpzspm_t);


mpzspv_t
mpzspv_init (spv_size_t len, const mpzspm_t mpzspm)
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
mpzspv_clear (mpzspv_t x, const mpzspm_t mpzspm)
{
  unsigned int i;
  
  ASSERT (mpzspv_verify (x, 0, 0, mpzspm));
  
  for (i = 0; i < mpzspm->sp_num; i++)
    sp_aligned_free (x[i]);
  
  free (x);
}


/* If filename == NULL, allocates memory for storage.
   If filename != NULL, opens a file set */

mpzspv_handle_t
mpzspv_init_handle (const char *filename, const spv_size_t len, 
                    const mpzspm_t mpzspm)
{
  mpzspv_handle_t handle;
  
  handle = (mpzspv_handle_t) malloc (sizeof (_mpzspv_handle_t));
  if (handle == NULL)
    return NULL;
  
  handle->mpzspm = mpzspm;

  if (filename == NULL)
    {
      handle->storage = 0;
      handle->mem = mpzspv_init (len, mpzspm);
      handle->files = NULL;
    }
  else
    {
      handle->storage = 1;
      handle->mem = NULL;
      handle->files = mpzspv_open_fileset (filename, len, mpzspm);
    }

  if (handle->mem == NULL && handle->files == NULL)
    {
      free (handle);
      handle = NULL;
    }

  return handle;
}


void
mpzspv_clear_handle (mpzspv_handle_t handle)
{
  if (handle == NULL)
    return;
  
  if (IN_MEMORY(handle))
    {
      mpzspv_clear (handle->mem, handle->mpzspm);
      handle->mem = NULL;
      handle->mpzspm = NULL;
    }
  else 
    {
      mpzspv_close_fileset (handle->files, handle->mpzspm);
      handle->files = NULL;
      handle->mpzspm = NULL;
    }

  free (handle);
}


/* check that:
 *  - each of the spv's is at least offset + len long
 *  - the data specified by (offset, len) is correctly normalised in the
 *    range [0, sp)
 *
 * return 1 for success, 0 for failure */

int
mpzspv_verify (const mpzspv_t x, const spv_size_t offset, 
               const spv_size_t len, const mpzspm_t mpzspm)
{
  unsigned int i;
  
  for (i = 0; i < mpzspm->sp_num; i++)
    {
      if (spv_verify (x[i] + offset, len, mpzspm->spm[i]->sp) == 0)
        return 0;
    }

  return 1;
}

void
mpzspv_set (mpzspv_t r, const spv_size_t r_offset, const mpzspv_t x, 
    const spv_size_t x_offset, const spv_size_t len, const mpzspm_t mpzspm)
{
  unsigned int i;
  
  ASSERT (mpzspv_verify (r, r_offset + len, 0, mpzspm));
  ASSERT (mpzspv_verify (x, x_offset, len, mpzspm));
  
  for (i = 0; i < mpzspm->sp_num; i++)
    spv_set (r[i] + r_offset, x[i] + x_offset, len);
}

void
mpzspv_revcopy (mpzspv_t r, const spv_size_t r_offset, const mpzspv_t x, 
    const spv_size_t x_offset, const spv_size_t len, const mpzspm_t mpzspm)
{
  unsigned int i;
  
  ASSERT (mpzspv_verify (r, r_offset + len, 0, mpzspm));
  ASSERT (mpzspv_verify (x, x_offset, len, mpzspm));
  
  for (i = 0; i < mpzspm->sp_num; i++)
    spv_rev (r[i] + r_offset, x[i] + x_offset, len);
}

void
mpzspv_set_sp (mpzspv_t r, const spv_size_t offset, const sp_t c, 
    const spv_size_t len, const mpzspm_t mpzspm)
{
  unsigned int i;
  
  ASSERT (mpzspv_verify (r, offset + len, 0, mpzspm));
  ASSERT (c < SP_MIN); /* not strictly necessary but avoids mod functions */
  
  for (i = 0; i < mpzspm->sp_num; i++)
    spv_set_sp (r[i] + offset, c, len);
}

void
mpzspv_neg (mpzspv_t r, const spv_size_t r_offset, const mpzspv_t x, 
    const spv_size_t x_offset, const spv_size_t len, const mpzspm_t mpzspm)
{
  unsigned int i;
  
  ASSERT (mpzspv_verify (r, r_offset + len, 0, mpzspm));
  ASSERT (mpzspv_verify (x, x_offset, len, mpzspm));
  
  for (i = 0; i < mpzspm->sp_num; i++)
    spv_neg (r[i] + r_offset, x[i] + x_offset, len, mpzspm->spm[i]->sp);
}

void
mpzspv_add (mpzspv_t r, const spv_size_t r_offset, const mpzspv_t x, 
            const spv_size_t x_offset, const mpzspv_t y, 
            const spv_size_t y_offset, const spv_size_t len, 
            const mpzspm_t mpzspm)
{
  unsigned int i;
  
  ASSERT (mpzspv_verify (r, r_offset + len, 0, mpzspm));
  ASSERT (mpzspv_verify (x, x_offset, len, mpzspm));
  
  for (i = 0; i < mpzspm->sp_num; i++)
    spv_add (r[i] + r_offset, x[i] + x_offset, y[i] + y_offset, len, 
             mpzspm->spm[i]->sp);
}

void
mpzspv_reverse (mpzspv_t x, const spv_size_t offset, const spv_size_t len, 
                const mpzspm_t mpzspm)
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

/* convert mpz to CRT representation, naive version */
static void
mpzspv_from_mpzv_slow (mpzspv_t x, const spv_size_t offset, const mpz_t mpz, 
                       const mpzspm_t mpzspm, ATTRIBUTE_UNUSED mpz_t rem,
		       const unsigned int sp_num)
{
  unsigned int j;

  if (mpz_sgn (mpz) == 0)
    {
      for (j = 0; j < sp_num; j++)
        x[j][offset] = 0;
    }
  else for (j = 0; j < sp_num; j++)
    { 
#if SP_TYPE_BITS > GMP_LIMB_BITS
      mpz_tdiv_r(rem, mpz, mpzspm->spm[j]->mp_sp);
      x[j][offset] = mpz_get_sp(rem);
#else
      x[j][offset] = mpn_mod_1 (PTR(mpz), SIZ(mpz),
                              (mp_limb_t) mpzspm->spm[j]->sp);
#endif
    }
}


/* convert mpzvi to CRT representation, fast version, assumes
   mpzspm->T has been precomputed (see mpzspm.c) */
static void
mpzspv_from_mpzv_fast (mpzspv_t x, const spv_size_t offset, mpz_t mpzvi,
                       const mpzspm_t mpzspm, 
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
    const spv_size_t len, const mpzspm_t mpzspm)
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
#endif
    /* Multi-threading with dynamic scheduling slows things down */

    mpz_init(rem);
#if defined(_OPENMP)
#pragma omp for schedule(static)
#endif
    for (i = 0; i < (long) len; i++)
    {
      ASSERT(mpz_sgn (mpzv[i]) >= 0); /* We can't handle negative values */
      if (mpzspm->T == NULL)
        mpzspv_from_mpzv_slow (x, i + offset, mpzv[i], mpzspm, rem, sp_num);
      else
        mpzspv_from_mpzv_fast (x, i + offset, mpzv[i], mpzspm, rem, sp_num);
    }
    mpz_clear(rem);

#if defined(_OPENMP)
  }
#endif
}


/* See: Daniel J. Bernstein and Jonathan P. Sorenson,
 * Modular Exponentiation via the explicit Chinese Remainder Theorem
 *
 * memory: mpzspm->sp_num floats */

static inline void
mpzspv_to_mpz(mpz_t res, const mpzspv_t x, const spv_size_t offset, 
              const mpzspm_t mpzspm, mpz_t mt)
{
  unsigned int i;
  float f = 0.5;
  mpz_set_ui (res, 0);

  for (i = 0; i < mpzspm->sp_num; i++)
    {
      const sp_t t = sp_mul (x[i][offset], mpzspm->crt3[i], 
          mpzspm->spm[i]->sp, mpzspm->spm[i]->mul_c);

#if SP_TYPE_BITS > GMP_LIMB_BITS
      mpz_set_sp (mt, t);
      mpz_addmul (res, mpzspm->crt1[i], mt);
#else
      mpz_addmul_ui (res, mpzspm->crt1[i], t);
#endif

      f += (float) t * mpzspm->prime_recip[i];
    }

  mpz_add (res, res, mpzspm->crt2[(unsigned int) f]);
}

void
mpzspv_to_mpzv (mpzspv_t x, const spv_size_t offset, mpzv_t mpzv,
                const spv_size_t len, const mpzspm_t mpzspm)
{
  mpzspv_to_mpzv_file (x, offset, NULL, mpzv, NULL, len, len, mpzspm);
}


void
mpzspv_to_mpzv_file (mpzspv_t x, const spv_size_t offset, 
    FILE **sp_files, mpzv_t mpzv, FILE * const mpz_file, 
    const spv_size_t len, const spv_size_t blocklen, const mpzspm_t mpzspm)
{
  spv_size_t len_done = 0;
  mpz_t mt;

  if (!sp_files) {
    ASSERT (mpzspv_verify (x, offset, len, mpzspm));
  } else {
    ASSERT (mpzspv_verify (x, offset + MIN(len, blocklen), 0, mpzspm));
  }
  
#if SP_TYPE_BITS > GMP_LIMB_BITS
  mpz_init (mt);
#endif

  while (len_done < len) 
    {
      const spv_size_t len_now = MIN(len - len_done, blocklen);
      spv_size_t l;
      
      if (sp_files != NULL) {
        mpzspv_seek_and_read (x, offset, sp_files, offset + len_done, len_now, mpzspm);
      }
      
      for (l = 0; l < len_now; l++)
        {
          mpzspv_to_mpz(mpzv[len_done + l], x, l + offset, mpzspm, mt);
          
          /* Write this mpz_t to file, if requested */
          if (mpz_file != NULL) {
            if (mpz_out_raw(mpz_file, mpzv[len_done + l]) == 0) {
              abort();
            }
          }
        }
      len_done += len_now;
    }
#if SP_TYPE_BITS > GMP_LIMB_BITS
  mpz_clear (mt);
#endif
}


void
mpzspv_fromto_mpzv (mpzspv_handle_t x, const spv_size_t offset, 
    const spv_size_t len, 
    mpz_producerfunc_t producer, void * producer_state, 
    mpz_consumerfunc_t consumer, void * consumer_state)
{
  const unsigned int sp_num = x->mpzspm->sp_num;
  spv_size_t block_len = 16384, len_done = 0, buffer_offset = 0;
  mpzspv_t buffer;
  mpz_t mpz1, mpz2, mt;

  ASSERT (sizeof (mp_limb_t) >= sizeof (sp_t));

  mpz_init(mpz1);
  mpz_init(mpz2);
  mpz_init(mt);

  if (IN_MEMORY(x))
    {
      block_len = len; /* Do whole thing at once */
      buffer = x->mem;
      buffer_offset = offset;
    }
  else
    {
      /* Do piecewise, using a temp buffer */
      buffer = mpzspv_init (block_len, x->mpzspm);
      ASSERT_ALWAYS (buffer != NULL);
    }

  while (len_done < len)
  {
    const spv_size_t len_now = MIN(len - len_done, block_len);
    spv_size_t i;

    /* Read x from disk files */
    if (consumer_state != NULL && ON_DISK(x)) {
      mpzspv_seek_and_read (buffer, buffer_offset, x->files, 
                            offset + len_done, len_now, x->mpzspm);
    }

    for (i = 0; i < len_now; i++)
      {
        if (producer_state != NULL)
          {
            if (producer)
              {
                /* Get new mpz1 from producer */
                (*producer)(producer_state, mpz1);
              } else {
                /* Get new mpz1 from listz_t */
                mpz_set (mpz1, ((listz_t)producer_state)[len_done + i]);
              }
          }
        
        if (consumer_state != NULL)
          {
            /* Convert NTT entry to mpz2 */
            mpzspv_to_mpz (mpz2, buffer, buffer_offset + i, x->mpzspm, mt);
            if (consumer != NULL)
              {
                /* Give mpz2 to consumer */
                mpz_mod (mpz2, mpz2, x->mpzspm->modulus);
                (*consumer)(consumer_state, mpz2);
              } else {
                mpz_mod (((listz_t)consumer_state)[len_done + i], mpz2, 
                         x->mpzspm->modulus);
              }
          }
        
        if (producer_state != NULL)
          {
            /* Convert the mpz1 we got from producer to NTT */
            mpzspv_from_mpzv_slow (buffer, buffer_offset + i, mpz1, x->mpzspm, 
                                   mt, sp_num);
          }
      }

    /* Write x to disk files */
    if (producer_state != NULL && ON_DISK(x)) {
      mpzspv_seek_and_write (buffer, buffer_offset, x->files, offset + len_done, 
                             len_now, x->mpzspm);
    }
    len_done += len_now;
    /* If we write NTT data to memory, we need to advance the offset to fill 
       the entire array. If we use a temp buffer, we reuse the same buffer 
       each time */
    if (IN_MEMORY(x))
      buffer_offset += len_now;
  }
  mpz_clear(mpz1);
  mpz_clear(mpz2);
  mpz_clear(mt);

  if (!IN_MEMORY(x))
    mpzspv_clear (buffer, x->mpzspm);
}


void
mpzspv_pwmul (mpzspv_t r, const spv_size_t r_offset, const mpzspv_t x, 
    const spv_size_t x_offset, const mpzspv_t y, const spv_size_t y_offset, 
    const spv_size_t len, const mpzspm_t mpzspm)
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
mpzspv_normalise (mpzspv_t x, const spv_size_t offset, const spv_size_t len,
                  const mpzspm_t mpzspm)
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
mpzspv_random (mpzspv_t x, const spv_size_t offset, const spv_size_t len, 
               const mpzspm_t mpzspm)
{
  unsigned int i;

  ASSERT (mpzspv_verify (x, offset, len, mpzspm));

  for (i = 0; i < mpzspm->sp_num; i++)
    spv_random (x[i] + offset, len, mpzspm->spm[i]->sp);
}


/* Adds or multiplies sp_t's from x and a file and stores result in r. 
   r[i] = x[i] + f[i] (+ x[i + wrap_size] + f[i + wrap_size] ...)
   Adds if add_or_mul == 0 */
static void
add_or_mul_file (spv_t r, const spv_t x, FILE *f, const spv_size_t len, 
    const spv_size_t wrap_size, const spv_size_t block_len, 
    const int add_or_mul, const spm_t spm)
{
  spv_t tmp_block;
  spv_size_t nr_read = 0;

  if (len == 0)
    return; 

  ASSERT(block_len > 0);
  ASSERT(wrap_size > 0);
  ASSERT(block_len <= wrap_size); /* We assume at most 1 wrap per block */

  tmp_block = (spv_t) sp_aligned_malloc (block_len * sizeof (sp_t));
  if (tmp_block == NULL) 
    {
      fprintf (stderr, "add_file(): could not allocate memory\n");
      abort();
    }

  while (nr_read < len)
    {
      const spv_size_t nr_now = MIN(len - nr_read, block_len);
      const spv_size_t offset_within_wrap = nr_read % wrap_size;
      const spv_size_t len_before_wrap = 
        MIN(nr_now, wrap_size - offset_within_wrap);
      const spv_size_t len_after_wrap = nr_now - len_before_wrap;

      fread (tmp_block, sizeof(sp_t), nr_now, f);

      if (add_or_mul == 0)
        spv_add (r + offset_within_wrap, x + nr_read, tmp_block, 
                 len_before_wrap, spm->sp);
      else
        spv_pwmul (r + offset_within_wrap, x + nr_read, tmp_block, 
                   len_before_wrap, spm->sp, spm->mul_c);

      if (len_after_wrap != 0)
        {
          if (add_or_mul == 0)
            spv_add (r, x + nr_read + len_before_wrap, 
                     tmp_block + len_before_wrap, len_after_wrap, 
                     spm->sp);
          else
            spv_pwmul (r, x + nr_read + len_before_wrap, 
                       tmp_block + len_before_wrap, len_after_wrap, 
                       spm->sp, spm->mul_c);
        }
      nr_read += nr_now;
    }
  sp_aligned_free (tmp_block);
}


/* Adds or multiplies sp_t's from two file and stores result in r. 
   Adds if add_or_mul == 0 */
static void
add_or_mul_2file (spv_t r, FILE *f1, FILE *f2, const spv_size_t len, 
    const spv_size_t wrap_size, const spv_size_t block_len, 
    const int add_or_mul, const spm_t spm)
{
  spv_t tmp_block1, tmp_block2;
  spv_size_t nr_read = 0;

  if (len == 0)
    return; 

  ASSERT(block_len > 0);
  ASSERT(wrap_size > 0);
  ASSERT(block_len <= wrap_size); /* We assume at most 1 wrap per block */

  tmp_block1 = (spv_t) sp_aligned_malloc (block_len * sizeof (sp_t));
  tmp_block2 = (spv_t) sp_aligned_malloc (block_len * sizeof (sp_t));
  if (tmp_block1 == NULL || tmp_block2 == NULL) 
    {
      fprintf (stderr, "add_file(): could not allocate memory\n");
      abort();
    }

  while (nr_read < len)
    {
      const spv_size_t nr_now = MIN(len - nr_read, block_len);
      const spv_size_t offset_within_wrap = nr_read % wrap_size;
      const spv_size_t len_before_wrap = 
        MIN(nr_now, wrap_size - offset_within_wrap);
      const spv_size_t len_after_wrap = nr_now - len_before_wrap;

      fread (tmp_block1, sizeof(sp_t), nr_now, f1);
      fread (tmp_block2, sizeof(sp_t), nr_now, f2);

      if (add_or_mul == 0)
        spv_add (r + offset_within_wrap, tmp_block1, tmp_block2, 
                 len_before_wrap, spm->sp);
      else
        spv_pwmul (r + offset_within_wrap, tmp_block1, tmp_block2, 
                   len_before_wrap, spm->sp, spm->mul_c);

      if (len_after_wrap != 0)
        {
          if (add_or_mul == 0)
            spv_add (r, tmp_block1 + len_before_wrap, 
                     tmp_block2 + len_before_wrap, len_after_wrap, spm->sp);
          else
            spv_pwmul (r, tmp_block1 + len_before_wrap, 
                       tmp_block2 + len_before_wrap, len_after_wrap, spm->sp, 
                       spm->mul_c);
        }
      nr_read += nr_now;
    }
  
  sp_aligned_free (tmp_block1);
  sp_aligned_free (tmp_block2);
}


static size_t 
seek_and_read_sp (spv_t ptr, const size_t nread, const size_t offset, FILE *f)
{
  size_t r;
  if (fseek (f, offset * sizeof(sp_t), SEEK_SET) != 0)
    {
      fprintf (stderr, "seek_and_read(): fseek() returned error %d\n", 
               errno);
      abort ();
    }
  
  r = fread(ptr, sizeof(sp_t), nread, f);
  if (r != nread)
    {
      fprintf (stderr, "seek_and_read(): Error reading data, r = %lu, errno = %d\n",
               (unsigned long) r, errno);
      abort();
    }

  return r;
}

static void 
mpzspv_seek_and_read (mpzspv_t dst, spv_size_t offset, FILE **sp_files, 
                      const size_t fileoffset, size_t nread, mpzspm_t mpzspm)
{
  unsigned int j;
  for (j = 0; j < mpzspm->sp_num; j++)
  {
    seek_and_read_sp (dst[j] + offset, nread, fileoffset, sp_files[j]);
  }
}

static size_t 
seek_and_write_sp (const spv_t ptr, const size_t nwrite, const size_t offset, FILE *f)
{
  size_t r;
  if (fseek (f, offset * sizeof(sp_t), SEEK_SET) != 0)
    {
      fprintf (stderr, "seek_and_read(): fseek() returned error %d\n", 
               errno);
      abort ();
    }
  
  r = fwrite(ptr, sizeof(sp_t), nwrite, f);
  if (r != nwrite)
    {
      fprintf (stderr, "seek_and_read(): Error writing data, r = %lu, errno = %d\n",
               (unsigned long) r, errno);
      abort();
    }
  fflush(f);

  return r;
}


static void 
mpzspv_seek_and_write (mpzspv_t src, const spv_size_t offset, FILE **sp_files, 
                       const size_t fileoffset, const size_t nwrite, 
                       const mpzspm_t mpzspm)
{
  unsigned int j;
  for (j = 0; j < mpzspm->sp_num; j++)
  {
    seek_and_write_sp (src[j] + offset, nwrite, fileoffset, sp_files[j]);
  }
}

static void
mul_dct_file (const spv_t r, const spv_t spv, FILE *dct_file, 
              const spv_size_t dftlen, const spv_size_t blocklen, const spm_t spm)
{
  const spv_size_t dctlen = dftlen / 2 + 1;
  spv_size_t nr_read = 0, i;
  unsigned long m = 5UL;
  spv_t tmp;
  
  ASSERT(dftlen % 2 == 0);
  if (dftlen == 0)
    return;
  
  tmp = (spv_t) sp_aligned_malloc (MIN(blocklen, dctlen) * sizeof (sp_t));
  if (tmp == NULL) 
    {
      fprintf (stderr, "mul_dct_file_blocks(): could not allocate memory\n");
      abort();
    }

  while (nr_read < dctlen)
    {
      const spv_size_t read_now = MIN(dctlen - nr_read, blocklen);
      const spv_size_t mul_now = MIN(dctlen - nr_read - 1, blocklen);
      
      seek_and_read_sp (tmp, read_now, nr_read, dct_file);
      
      i = 0;
      if (nr_read == 0)
        {
          r[0] = sp_mul (spv[0], tmp[0], spm->sp, spm->mul_c);
          i = 1;
        }
      
      for ( ; i < mul_now; i++)
        {
          const spv_size_t j = nr_read + i;
          /* This works, but why? */
          if (3*j > m)
            m = 2UL * m + 1;
          
          r[2*j] = sp_mul (spv[2*j], tmp[i], spm->sp, spm->mul_c);
          r[m - 2*j] = sp_mul (spv[m - 2*j], tmp[i], spm->sp, spm->mul_c);
        }
      nr_read += read_now;
      if (nr_read == dctlen)
        {
#ifdef USE_VALGRIND
          VALGRIND_CHECK_VALUE_IS_DEFINED(tmp[i]);
#endif
          r[1] = sp_mul (spv[1], tmp[i], spm->sp, spm->mul_c);
        }
    }
  sp_aligned_free(tmp);
}


/* Multiply the DFT of a polynomial by the DCT-I of a reciprocal Laurent
   polynomial. */
static void
mul_dct(spv_t r, const spv_t spv, const spv_t dct, const spv_size_t len, 
        const spm_t spm)
{
  unsigned long m = 5UL, i;
  
  if (len > 0)
    r[0] = sp_mul (spv[0], dct[0], spm->sp, spm->mul_c);
  if (len > 1)
    r[1] = sp_mul (spv[1], dct[len / 2UL], spm->sp, spm->mul_c);
  
  ASSERT(len % 2 == 0);
  for (i = 2UL; i < len; i += 2UL)
    {
      /* This works, but why? */
      if (i + i / 2UL > m)
        m = 2UL * m + 1;
      
      r[i] = sp_mul (spv[i], dct[i / 2UL], spm->sp, spm->mul_c);
      r[m - i] = sp_mul (spv[m - i], dct[i / 2UL], spm->sp, 
                         spm->mul_c);
    }
}


static void
one_fft_file(spv_t spv, FILE *file, const spv_size_t offset, 
             const spv_size_t len, const spv_size_t ntt_size, const int monic,
             const spv_size_t block_len, const spm_t spm)
{
  spv_size_t log2_ntt_size = ceil_log_2 (ntt_size);
  if (file != NULL)
    {
      seek_and_read_sp (spv, ntt_size, offset, file);
      if (ntt_size < len)
        add_or_mul_file (spv, spv, file, len - ntt_size, ntt_size, block_len, 
                         0, spm);
    } 
  else if (ntt_size < len) 
    {
      spv_size_t j;
      for (j = ntt_size; j < len; j += ntt_size)
        spv_add (spv, spv, spv + j, ntt_size, spm->sp);
    }

  if (ntt_size > len)
    spv_set_zero (spv + len, ntt_size - len);

  if (monic)
    spv[len % ntt_size] = sp_add (spv[len % ntt_size], 1, spm->sp);

  spv_ntt_gfp_dif (spv, log2_ntt_size, spm);
}


/* Do multiplication via NTT. Depending on the value of "steps", does 
   in-place forward transform of x, in-place forward transform of y, 
   pair-wise multiplication of x by y to r, in-place inverse transform of r. 
   Contrary to calling these three operations separately, this function does 
   all steps on a small-prime vector at a time, resulting in slightly 
   better cache efficiency.
   Input and output spv_t's can be stored in files. Files are read or written
   beginning at the current file pointer position. If x_files is non-NULL and
   x is non-NULL, then the data in x is replaced by the data from the files, 
   likewise for y. If r and r_files are non-NULL, then output data is written
   to both r and r_files.
   It is permissible to let any combination of x, y, and r point at the same 
   memory, in which case the pointers x, y, or r, respectively must be 
   identical. Other overlap is not permitted.
*/

void
mpzspv_mul_ntt_file (mpzspv_handle_t r, const spv_size_t offsetr, 
    mpzspv_handle_t x, const spv_size_t offsetx, const spv_size_t lenx, 
    mpzspv_handle_t y, const spv_size_t offsety, const spv_size_t leny, 
    const spv_size_t ntt_size, const int monic, const spv_size_t monic_pos, 
    const int steps)
{
  const spv_size_t block_len = 16384;
  spv_size_t log2_ntt_size;
  int i;
  const int do_fft1 = (steps & NTT_MUL_STEP_FFT1) != 0;
  const int do_fft2 = (steps & NTT_MUL_STEP_FFT2) != 0;
  const int do_pwmul = (steps & NTT_MUL_STEP_MUL) != 0;
  const int do_pwmul_dct = (steps & NTT_MUL_STEP_MULDCT) != 0;
  const int do_ifft = (steps & NTT_MUL_STEP_IFFT) != 0;

  mpzspm_t mpzspm = NULL;

  if (x != NULL) 
    mpzspm = x->mpzspm;
  if (y != NULL) 
    {
      ASSERT_ALWAYS (mpzspm == NULL || y->mpzspm == mpzspm);
      mpzspm = y->mpzspm;
    }
  if (r != NULL)
    {
      ASSERT_ALWAYS (mpzspm == NULL || r->mpzspm == mpzspm);
      mpzspm = r->mpzspm;
    }

  if (x == y && offsetx == offsety && do_fft1 && do_fft2)
    {
      fprintf (stderr, "mpzspv_mul_ntt_file(): Error, x=y and forward "
               "transform requested for both\n");
      abort();
    }
  
  if (do_pwmul && do_pwmul_dct)
    {
      fprintf (stderr, "mpzspv_mul_ntt_file(): Error, both PWMUL "
               "and PWMULDCT requested\n");
      abort();
    }
  
  if (IN_MEMORY(x))
    {
      ASSERT (mpzspv_verify (x->mem, offsetx, lenx, mpzspm));
      ASSERT (mpzspv_verify (x->mem, offsetx + ntt_size, 0, mpzspm));
    }
  if (IN_MEMORY(y)) 
    {
      ASSERT (mpzspv_verify (y->mem, offsety, leny, mpzspm));
      ASSERT (mpzspv_verify (y->mem, offsety + ntt_size, 0, mpzspm));
    }
  if (IN_MEMORY(r))
    {
      ASSERT (mpzspv_verify (r->mem, offsetr + ntt_size, 0, mpzspm));
    }
  
  log2_ntt_size = ceil_log_2 (ntt_size);

  /* Need parallelization at higher level (e.g., handling a branch of the 
     product tree in one thread) to make this worthwhile for ECM */

#if defined(_OPENMP) && MPZSPV_MUL_NTT_OPENMP
#pragma omp parallel if (ntt_size > 32768)
  {
#pragma omp for
#endif
  for (i = 0; i < (int) mpzspm->sp_num; i++)
    {
      const spm_t spm = mpzspm->spm[i];
      spv_t tmp = NULL;
      char tmp_content = ' ';

      /* If we do any transform, we need some memory to do the transform in. 
         The point-wise multiply needa a small block of memory, too, to avoid
         having an fread() for each sp_t.
         If x is in memory, then any forward transform of x is done in-place.
         Same for y. If r is in memory, we can use that as temp storage, so
         long as we don't overwrite input data we still need. */
      
      if ((do_fft1 && ON_DISK(x)) || (do_fft2 && ON_DISK(y)) || 
          (do_ifft && ON_DISK(r)))
        {
          tmp = (spv_t) sp_aligned_malloc (ntt_size * sizeof (sp_t));
          if (tmp == NULL)
            {
              fprintf (stderr, "Cannot allocate tmp memory in "
                       "mpzspv_mul_ntt_file()\n");
              abort();
            }
        }

      spv_t spvx = IN_MEMORY(x) ? x->mem[i] + offsetx : tmp;
      spv_t spvy = IN_MEMORY(y) ? y->mem[i] + offsety : tmp;
      spv_t spvr = IN_MEMORY(r) ? r->mem[i] + offsetr : tmp;

      if (do_fft1) 
        {
          ASSERT_ALWAYS(x != NULL);
          one_fft_file (spvx, ON_DISK(x) ? x->files[i] : NULL, offsetx, 
                        lenx, ntt_size, monic, block_len, spm);
          if (ON_DISK(x))
            {
              tmp_content = 'x';
              seek_and_write_sp (spvx, ntt_size, offsetx, x->files[i]);
            }
        }

      if (do_fft2) 
        {
          ASSERT_ALWAYS(y != NULL);
          one_fft_file(spvy, ON_DISK(y) ? y->files[i] : NULL, offsety, 
                       leny, ntt_size, monic, block_len, spm);
          if (ON_DISK(y))
            {
              tmp_content = 'y';
              seek_and_write_sp (spvy, ntt_size, offsety, y->files[i]);
            }
        }

      if (do_pwmul) 
        {
          ASSERT_ALWAYS(x != NULL && y != NULL && r != NULL);
          if (IN_MEMORY(x) && IN_MEMORY(y))
            {
              spv_pwmul (spvr, spvx, spvy, ntt_size, spm->sp, spm->mul_c);
            }
          else if (ON_DISK(x) && IN_MEMORY(y))
            {
              if (tmp_content == 'x')
                spv_pwmul (spvr, tmp, spvy, ntt_size, spm->sp, spm->mul_c);
              else
                add_or_mul_file (spvr, spvy, x->files[i], ntt_size, ntt_size, 
                                 block_len, 1, spm);
            }
          else if (IN_MEMORY(x) && ON_DISK(y))
            {
              if (tmp_content == 'y')
                spv_pwmul (spvr, spvx, tmp, ntt_size, spm->sp, spm->mul_c);
              else
                add_or_mul_file (spvr, spvx, y->files[i], ntt_size, ntt_size, 
                                 block_len, 1, spm);
            }
          else /* x_files != NULL && y_files != NULL */
            {
              if (tmp_content == 'x')
                add_or_mul_file (spvr, tmp, y->files[i], ntt_size, ntt_size,
                                 block_len, 1, spm);
              else if (tmp_content == 'y')
                add_or_mul_file (spvr, tmp, x->files[i], ntt_size, ntt_size,
                                 block_len, 1, spm);
              else
                add_or_mul_2file (spvr, x->files[i], y->files[i], ntt_size, 
                                  ntt_size, block_len, 1, spm);
            }
          if (ON_DISK(r))
            tmp_content = 'r';
        }
      else if (do_pwmul_dct)
        {
          ASSERT_ALWAYS(x != NULL && y != NULL && r != NULL);
          if (ON_DISK(x) && tmp_content != 'x')
            {
              seek_and_read_sp (spvx, ntt_size, offsetx, x->files[i]);
            }
          if (ON_DISK(y))
            {
              mul_dct_file (spvr, spvx, y->files[i], ntt_size, 
                                   block_len, spm);
            } else { 
              mul_dct (spvr, spvx, spvy, ntt_size, spm);
            }
          if (ON_DISK(r))
            tmp_content = 'r';
        }

      if (do_ifft) 
        {
          ASSERT (sizeof (mp_limb_t) >= sizeof (sp_t));
          ASSERT_ALWAYS(r != NULL);

          if (ON_DISK(r) && tmp_content != 'r')
            {
              seek_and_read_sp (spvr, ntt_size, offsetr, r->files[i]);
            }
          spv_ntt_gfp_dit (spvr, log2_ntt_size, spm);

          /* spm->sp - (spm->sp - 1) / ntt_size is the inverse of ntt_size */
          spv_mul_sp (spvr, spvr, spm->sp - (spm->sp - 1) / ntt_size,
              ntt_size, spm->sp, spm->mul_c);

          if (monic_pos)
            spvr[monic_pos % ntt_size] = sp_sub (spvr[monic_pos % ntt_size],
                1, spm->sp);
        }

      if (ON_DISK(r) && tmp_content == 'r')
        {
          seek_and_write_sp (spvr, ntt_size, offsetr, r->files[i]);
        }
      if (tmp != NULL) 
        {
          sp_aligned_free (tmp);
        }
    }
#if defined(_OPENMP) && MPZSPV_MUL_NTT_OPENMP
  }
#endif
}

void
mpzspv_mul_ntt (mpzspv_t r, const spv_size_t offsetr, 
    mpzspv_t x, const spv_size_t offsetx, const spv_size_t lenx,
    mpzspv_t y, const spv_size_t offsety, const spv_size_t leny,
    const spv_size_t ntt_size, const int monic, const spv_size_t monic_pos, 
    mpzspm_t mpzspm, const int steps)
{
  _mpzspv_handle_t r_handle = {0, mpzspm, r, NULL};
  _mpzspv_handle_t x_handle = {0, mpzspm, x, NULL};
  _mpzspv_handle_t y_handle = {0, mpzspm, y, NULL};

  mpzspv_mul_ntt_file((r != NULL) ? &r_handle : NULL, offsetr, 
                      (x != NULL) ? &x_handle : NULL, offsetx, lenx, 
                      (y != NULL) ? &y_handle : NULL, offsety, leny, 
                      ntt_size, monic, monic_pos, steps);
}


/* Computes a DCT-I of the length dctlen. Input is the spvlen coefficients
   in spv. */

void
mpzspv_to_dct1 (mpzspv_handle_t dct, const mpzspv_handle_t spv,  
                const spv_size_t spvlen, const spv_size_t dctlen)
{
  const spv_size_t ntt_size = 2 * (dctlen - 1); /* Length for the DFT */
  const spv_size_t log2_l = ceil_log_2 (ntt_size);
  int j;

  ASSERT_ALWAYS (dct->mpzspm == spv->mpzspm);

#ifdef _OPENMP
#pragma omp parallel private(j)
  {
#pragma omp for
#endif
  for (j = 0; j < (int) spv->mpzspm->sp_num; j++)
    {
      const spm_t spm = spv->mpzspm->spm[j];
      spv_size_t i;
      
      spv_t tmp = (spv_t) sp_aligned_malloc (ntt_size * sizeof (sp_t));
      if (tmp == NULL)
        {
          fprintf (stderr, "Cannot allocate tmp memory in "
                   "mpzspv_to_dct1()\n");
          abort();
        }

      if (ON_DISK(spv))
        {
          seek_and_read_sp (tmp, spvlen, 0, spv->files[j]);
        } else {
          /* Copy spv to tmp */
          spv_set (tmp, spv->mem[j], spvlen);
        }
      /* Make a symmetric copy of input coefficients in tmp. E.g., 
         with spv = [3, 2, 1], spvlen = 3, dctlen = 5 (hence ntt_size = 8), 
         we want tmp = [3, 2, 1, 0, 0, 0, 1, 2] */
      spv_rev (tmp + ntt_size - spvlen + 1, tmp + 1, spvlen - 1);
      /* Now we have [3, 2, 1, ?, ?, ?, 1, 2]. Fill the ?'s with zeros. */
      spv_set_sp (tmp + spvlen, (sp_t) 0, ntt_size - 2 * spvlen + 1);

#if 0
      printf ("mpzspv_to_dct1: tmp[%d] = [", j);
      for (i = 0; i < ntt_size; i++)
          printf ("%lu, ", tmp[i]);
      printf ("]\n");
#endif
      
      spv_ntt_gfp_dif (tmp, log2_l, spm);

#if 0
      printf ("mpzspv_to_dct1: tmp[%d] = [", j);
      for (i = 0; i < ntt_size; i++)
          printf ("%lu, ", tmp[i]);
      printf ("]\n");
#endif

      /* The forward transform is scrambled. We want elements [0 ... ntt_size/2]
         of the unscrabled data, that is all the coefficients with the most 
         significant bit in the index (in log2(ntt_size) word size) unset, plus the 
         element at index ntt_size/2. By scrambling, these map to the elements with 
         even index, plus the element at index 1. 
         The elements with scrambled index 2*i are stored in h[i], the
         element with scrambled index 1 is stored in h[params->ntt_size] */
  
#ifdef WANT_ASSERT
      /* Test that the coefficients are symmetric (if they were unscrambled)
         and that our algorithm for finding identical coefficients in the 
         scrambled data works */
      {
        spv_size_t m = 5;
        for (i = 2; i < ntt_size; i += 2L)
          {
            /* This works, but why? */
            if (i + i / 2L > m)
                m = 2L * m + 1L;

            ASSERT (tmp[i] == tmp[m - i]);
#if 0
            printf ("mpzspv_to_dct1: DFT[%lu] == DFT[%lu]\n", i, m - i);
#endif
          }
      }
#endif

      /* Copy coefficients to dct buffer */
      {
        spv_t out_buf = IN_MEMORY(dct) ? dct->mem[j] : tmp;
        const sp_t coeff_1 = tmp[1];
        for (i = 0; i < ntt_size / 2; i++)
          out_buf[i] = tmp[i * 2];
        out_buf[ntt_size / 2] = coeff_1;
        if (ON_DISK(dct))
          {
            /* Write data back to file */
            seek_and_write_sp (tmp, dctlen + 1, 0, dct->files[j]);
          }
      }

      sp_aligned_free(tmp);
    }
#ifdef _OPENMP
  }
#endif
}


ATTRIBUTE_UNUSED
static void
spv_print_vec (const char *msg, const spv_t spv, const spv_size_t l)
{
  spv_size_t i;
  printf ("%s [%lu", msg, spv[0]);
  for (i = 1; i < l; i++)
    printf (", %lu", spv[i]);
  printf ("]\n");
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
   Now handled by mpzspv_mul_ntt().
*/

static void
spv_sqr_reciprocal(const spv_size_t n, const spm_t spm, const spv_t spv, 
                   const sp_t max_ntt_size)
{
  const spv_size_t log2_n = ceil_log_2 (n);
  const spv_size_t len = ((spv_size_t) 2) << log2_n;
  const spv_size_t log2_len = 1 + log2_n;
  sp_t w1, w2, invlen;
  const sp_t sp = spm->sp, mul_c = spm->mul_c;
  spv_size_t i;

  /* Zero out NTT elements [n .. len-n] */
  spv_set_sp (spv + n, (sp_t) 0, len - 2*n + 1);

#if TRACE_ntt_sqr_reciprocal
  printf ("ntt_sqr_reciprocal: NTT vector mod %lu\n", sp);
  spv_print_vec ("ntt_sqr_reciprocal: before weighting:", spv, n);
#endif

  /* Compute the root for the weight signal, a 3rd primitive root 
     of unity */
  w1 = sp_pow (spm->prim_root, max_ntt_size / 3UL, sp, mul_c);
  /* Compute iw= 1/w */
  w2 = sp_pow (spm->inv_prim_root, max_ntt_size / 3UL, sp, mul_c);
#if TRACE_ntt_sqr_reciprocal
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

#if TRACE_ntt_sqr_reciprocal
  spv_print_vec ("ntt_sqr_reciprocal: after weighting:", spv, len);
#endif

  /* Forward DFT of dft[j] */
  spv_ntt_gfp_dif (spv, log2_len, spm);

#if TRACE_ntt_sqr_reciprocal
  spv_print_vec ("ntt_sqr_reciprocal: after forward transform:", 
                 spv, len);
#endif

  /* Square the transformed vector point-wise */
  spv_pwmul (spv, spv, spv, len, sp, mul_c);

#if TRACE_ntt_sqr_reciprocal
  spv_print_vec ("ntt_sqr_reciprocal: after point-wise squaring:", 
                 spv, len);
#endif

  /* Inverse transform of dft[j] */
  spv_ntt_gfp_dit (spv, log2_len, spm);

#if TRACE_ntt_sqr_reciprocal
  spv_print_vec ("ntt_sqr_reciprocal: after inverse transform:", 
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
  
#if TRACE_ntt_sqr_reciprocal
  spv_print_vec ("ntt_sqr_reciprocal: after un-weighting:", spv, len);
#endif

  /* Separate the coefficients of R in the wrapped-around product. */

  /* Set w1 = cuberoot(1)^l where cuberoot(1) is the same primitive
     3rd root of unity we used for the weight signal */
  w1 = sp_pow (spm->prim_root, max_ntt_size / 3UL, sp, mul_c);
  w1 = sp_pow (w1, len % 3UL, sp, mul_c);
  
  /* Set w2 = 1/(w1 - 1/w1). Incidentally, w2 = 1/sqrt(-3) */
  w2 = sp_inv (w1, sp, mul_c);
  w2 = sp_sub (w1, w2, sp);
  w2 = sp_inv (w2, sp, mul_c);
#if TRACE_ntt_sqr_reciprocal
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

#if TRACE_ntt_sqr_reciprocal
  spv_print_vec ("ntt_sqr_reciprocal: after un-wrapping:", spv, len);
#endif
}


/* Square an RLP */

void 
mpzspv_sqr_reciprocal (mpzspv_handle_t x, const spv_size_t n)
{
  const spv_size_t log2_n = ceil_log_2 (n);
  const spv_size_t len = ((spv_size_t) 2) << log2_n;
  int j;

  ASSERT(x->mpzspm->max_ntt_size % 3UL == 0UL);
  ASSERT(len % 3UL != 0UL);
  ASSERT(x->mpzspm->max_ntt_size % len == 0UL);

#ifdef _OPENMP
#pragma omp parallel
  {
#pragma omp for
#endif
    for (j = 0; j < (int) (x->mpzspm->sp_num); j++)
      {
        spv_t tmp;
        
        if (ON_DISK(x))
          {
            tmp = (spv_t) sp_aligned_malloc (len * sizeof (sp_t));

            if (tmp == NULL)
              {
                fprintf (stderr, "Cannot allocate tmp memory in "
                         "mpzspv_sqr_reciprocal_file()\n");
                abort();
              }
            seek_and_read_sp (tmp, n, 0, x->files[j]);
          }
        else
          {
            tmp = x->mem[j];
          }
        
        spv_sqr_reciprocal (n, x->mpzspm->spm[j], tmp, x->mpzspm->max_ntt_size);
        
        if (ON_DISK(x))
          {
            seek_and_write_sp (tmp, 2 * n - 1, 0, x->files[j]);
            sp_aligned_free (tmp);
          }
      }
#ifdef _OPENMP
    }
#endif
}

FILE **
mpzspv_open_fileset (const char *file_stem, const spv_size_t len, 
                     const mpzspm_t mpzspm)
{
  FILE **files;
  char *filename;
  const unsigned int spnum = mpzspm->sp_num;
  unsigned int i;
  
#ifdef HAVE_SETVBUF
  files = (FILE **) malloc (2 * spnum * sizeof(FILE *));
#else
  files = (FILE **) malloc (spnum * sizeof(FILE *));
#endif
  filename = (char *) malloc ((strlen(file_stem) + 10) * sizeof(FILE *));
  if (files == NULL || filename == NULL)
    {
      fprintf (stderr, "mpzspv_open_fileset(): could not allocate memory\n");
      abort ();
    }
  
  for (i = 0; i < spnum; i++)
    {
#ifdef HAVE_SETVBUF
      void **buffers = (void**)files + spnum;
#endif
      const size_t bufsize = 1<<22;
      
      sprintf (filename, "%s.%u", file_stem, i);
      files[i] = fopen(filename, "r+");
      if (files[i] == NULL)
        files[i] = fopen(filename, "w+");
      if (files[i] == NULL)
        {
          fprintf (stderr, 
                   "mpzspv_open_fileset(): error opening %s for writing\n", 
                   filename);
          while (i > 0)
            {
              fclose(files[--i]);
            }
          free (filename);
          free (files);
          abort(); /* For now, bomb out for easier debugging */
          return NULL;
        }
#ifdef HAVE_FALLOCATE
      /* Tell the file system to allocate space for the file, avoiding 
         fragementation */
      fallocate (fileno(files[i]), 0, (off_t) 0, len * sizeof(sp_t));
#endif
#ifdef HAVE_SETVBUF
      /* Allocate a bigger buffer so accesses during conversion from/to
         mpz_t's are more sequential. We need to remember the pointers
         to the buffers so we can free them later. For now we allocate
         twice as much memory for the FILE * array and put the pointers
         to the buffers in the second half. This is ugly - FIXME */
      buffers[i] = malloc (bufsize);
      setvbuf (files[i], (char *)buffers[i], _IOFBF, bufsize);
#endif
    }
  
  free (filename);
  
  return files;
}

void
mpzspv_close_fileset (FILE **files, const mpzspm_t mpzspm)
{
  unsigned int i;
#ifdef HAVE_SETVBUF
  void **buffers = (void**)files + mpzspm->sp_num;
#endif
  
  for (i = 0; i < mpzspm->sp_num; i++)
    {
      if (fclose(files[i]) != 0)
        {
          fprintf (stderr, 
                   "mpzspv_close_fileset(): fclose() set error code %d\n", 
                   errno);
          abort();
        }
      files[i] = NULL;
#ifdef HAVE_SETVBUF
      free (buffers[i]);
      buffers[i] = NULL;
#endif
    }
  
  free (files);
}

void
mpzspv_read (mpzspv_t mpzspv, const spv_size_t mpzspv_offset, 
             FILE **files, const spv_size_t file_offset,
             const spv_size_t len, const mpzspm_t mpzspm)
{
  unsigned int i;
  
  for (i = 0; i < mpzspm->sp_num; i++)
    {
      seek_and_read_sp (mpzspv[i] + mpzspv_offset, len, file_offset, files[i]);
    }
}

void
mpzspv_write (const mpzspv_t mpzspv, const spv_size_t mpzspv_offset, 
             FILE **files, const spv_size_t file_offset,
             const spv_size_t len, const mpzspm_t mpzspm)
{
  unsigned int i;
  
  for (i = 0; i < mpzspm->sp_num; i++)
    {
      seek_and_write_sp (mpzspv[i] + mpzspv_offset, len, file_offset, files[i]);
    }
}

static void
mpzspv_print_mem (const mpzspv_t mpzspv, const spv_size_t offset, 
              const spv_size_t len, const char *prefix, 
              const mpzspm_t mpzspm)
{
  unsigned int i;

  if (len == 0)
    {
      printf("%s: Zero length vector\n", prefix);
      return;
    }
  
  for (i = 0; i < mpzspm->sp_num; i++)
    {
      spv_size_t j;
      printf ("%s (%lu", prefix, mpzspv[i][offset]);
      for (j = 1; j < len; j++)
        {
          printf(", %lu", mpzspv[i][offset + j]);
        }
      printf (") (mod %lu)\n", mpzspm->spm[i]->sp);
    }
}

static void
mpzspv_print_file (FILE **files, const spv_size_t offset, 
              const spv_size_t len, const char *prefix, 
              const mpzspm_t mpzspm)
{
  unsigned int i;
  spv_t tmp;

  if (len == 0)
    {
      printf("%s: Zero length vector\n", prefix);
      return;
    }
  
  tmp = (spv_t) sp_aligned_malloc (len * sizeof (sp_t));

  for (i = 0; i < mpzspm->sp_num; i++)
    {
      spv_size_t j;
      seek_and_read_sp (tmp, len, offset, files[i]);
      printf ("%s (%lu", prefix, tmp[0]);
      for (j = 1; j < len; j++)
        {
          printf(", %lu", tmp[j]);
        }
      printf (") (mod %lu)\n", mpzspm->spm[i]->sp);
    }
}

void
mpzspv_print (mpzspv_handle_t handle, const spv_size_t offset, 
              const spv_size_t len, const char *prefix)
{
  if (IN_MEMORY(handle))
    {
      mpzspv_print_mem (handle->mem, offset, len, prefix, handle->mpzspm);
    }
  else
    {
      mpzspv_print_file (handle->files, offset, len, prefix, handle->mpzspm);
    }
}
