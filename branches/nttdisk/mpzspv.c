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
#ifdef HAVE_FCNTL_H
#include <fcntl.h>
#endif
#ifdef HAVE_AIO_H
#include <aio.h>
#endif
#include "ecm-impl.h"

#define TRACE_FUNC 0
#define TRACE_ntt_sqr_reciprocal 0
#define TRACE_ntt_mul 0
#define WANT_PROFILE 0

static void mpzspv_normalise (mpzspv_t, spv_size_t, spv_size_t, mpzspm_t);
static sp_t mpz_mod_spm (const mpz_t, const spm_t, mpz_t);
#ifdef HAVE_AIO_READ
static int mpzspv_lio_rw (struct aiocb *[], mpzspv_t, spv_size_t, FILE *,  
                          spv_size_t, spv_size_t, spv_size_t, const mpzspm_t, 
                          int);
static int mpzspv_lio_suspend (struct aiocb * const[], const mpzspm_t);
#else
static void mpzspv_seek_and_read (mpzspv_t, size_t, FILE *, spv_size_t, 
    spv_size_t, size_t, mpzspm_t);
static void mpzspv_seek_and_write (mpzspv_t, size_t, FILE *, spv_size_t, 
    spv_size_t, size_t, mpzspm_t);
#endif


static mpzspv_t
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

static void
mpzspv_clear (mpzspv_t x, const mpzspm_t mpzspm)
{
  unsigned int i;
  
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
  handle->len = len;

  if (filename == NULL)
    {
      handle->storage = 0;
      handle->mem = mpzspv_init (len, mpzspm);
      if (handle->mem == NULL)
        {
          free (handle);
          return NULL;
        }
      handle->file = NULL;
    }
  else
    {
      handle->storage = 1;
      handle->mem = NULL;

      handle->filename = (char *) malloc ((strlen(filename) + 5) * sizeof(char));
      if (handle->filename == NULL)
        {
          fprintf (stderr, "%s(): could not allocate memory for file name\n", 
                   __func__);
          free (handle);
          return NULL;
        }

      sprintf (handle->filename, "%s.ntt", filename);
      handle->file = fopen(handle->filename, "rb+");
      if (handle->file == NULL)
        handle->file = fopen(handle->filename, "wb+");
      if (handle->file == NULL)
        {
          fprintf (stderr, "%s(): error opening %s for writing\n", 
                   __func__, handle->filename);
          free (handle->filename);
          handle->filename = NULL;
          free (handle);
          return NULL;
        }
#ifdef HAVE_FALLOCATE
      /* Tell the file system to allocate space for the file, avoiding 
         fragementation */
      {
        int r;
        r = fallocate (fileno(handle->file), 0, (off_t) 0, 
                       handle->mpzspm->sp_num * len * sizeof(sp_t));
        if (r != 0)
          {
            int saved_errno = errno;
            perror (NULL);
            fprintf (stderr, "%s(): fallocate() returned %d, errno = %d\n", 
                     __func__, r, saved_errno);
          }
      }
#endif
#if defined(HAVE_SETVBUF) && defined(HAVE_AIO_READ)
      /* Set to unbuffered mode as we use aio_*() functions for reading
         in the background */
      setvbuf (handle->file, NULL, _IONBF, 0);
#endif
    }

  return handle;
}


void
mpzspv_clear_handle (mpzspv_handle_t handle)
{
  if (handle == NULL)
    return;
  
  ASSERT(mpzspv_verify_out (handle, 0, handle->len));
  if (mpzspv_handle_in_memory(handle))
    {
      mpzspv_clear (handle->mem, handle->mpzspm);
      handle->mem = NULL;
      handle->mpzspm = NULL;
    }
  else 
    {
      if (fclose(handle->file) != 0)
        {
         int saved_errno = errno;
          perror (NULL);
          fprintf (stderr, 
                   "%s(): fclose() set error code %d\n", 
                   __func__, saved_errno);
        } 
      else if (remove (handle->filename) != 0)
        {
         int saved_errno = errno;
          perror (NULL);
          fprintf (stderr, 
                   "%s(): remove() set error code %d\n", 
                   __func__, saved_errno);
        }
      free (handle->filename);
      handle->filename = NULL;
      handle->file = NULL;
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
mpzspv_verify_in (const mpzspv_handle_t x, const spv_size_t offset, 
                  const spv_size_t len)
{
  unsigned int i;
  
  if (x->storage == 0)
    {
      for (i = 0; i < x->mpzspm->sp_num; i++)
        {
          if (spv_verify_in (x->mem[i] + offset, len, x->mpzspm->spm[i]->sp) == 0)
            return 0;
        }
    }

  return 1;
}


/* check that:
 *  - each of the spv's is at least offset + len long
 *
 * return 1 for success, 0 for failure */

int
mpzspv_verify_out (const mpzspv_handle_t x, const spv_size_t offset, 
                   const spv_size_t len)
{
  unsigned int i;
  
  if (x->storage == 0)
    {
      for (i = 0; i < x->mpzspm->sp_num; i++)
        {
          if (spv_verify_out (x->mem[i] + offset, len, x->mpzspm->spm[i]->sp) == 0)
            return 0;
        }
    }

  return 1;
}


static inline spv_t
get_mem (mpzspv_handle_t r, const unsigned int i)
{
  return mpzspv_handle_in_memory(r) ? r->mem[i] : NULL;
}


static inline FILE *
get_file (mpzspv_handle_t r, ATTRIBUTE_UNUSED const unsigned int i)
{
  return mpzspv_handle_on_disk(r) ? r->file : NULL;
}


/* Returns the adjusted offset of vector i within the file if the NTT vectors 
   are stored on disk. If they are stored in memory, or if r == NULL, returns 
   the unmodified offset. */
static inline spv_size_t
adjust_offset (mpzspv_handle_t r, const unsigned int i, 
             const spv_size_t offset)
{
  return offset + (mpzspv_handle_on_disk(r) ? i * r->len : 0);
}


static inline void
mpzspv_elementwise (mpzspv_handle_t r, const spv_size_t r_offset,
    const mpzspv_handle_t x, const spv_size_t x_offset,
    const mpzspv_handle_t y, const spv_size_t y_offset,
    const spv_size_t len, const int operation)
{
  unsigned int i;

  ASSERT (r != NULL);
  ASSERT (mpzspv_verify_out (r, r_offset, len));
  ASSERT (x == NULL || r->mpzspm == x->mpzspm);
  ASSERT (x == NULL || mpzspv_verify_in (x, x_offset, len));
  ASSERT (y == NULL || r->mpzspm == y->mpzspm);
  ASSERT (y == NULL || mpzspv_verify_in (y, y_offset, len));

  for (i = 0; i < r->mpzspm->sp_num; i++)
    spv_elementwise (get_mem (r, i), get_file (r, i), adjust_offset (r, i, r_offset), 
                     get_mem (x, i), get_file (x, i), adjust_offset (x, i, x_offset),
                     get_mem (y, i), get_file (y, i), adjust_offset (y, i, y_offset), 
                     r->mpzspm->spm[i]->sp, r->mpzspm->spm[i]->mul_c, len, operation);
}


void
mpzspv_set (mpzspv_handle_t r, const spv_size_t r_offset, 
    const mpzspv_handle_t x, const spv_size_t x_offset, const spv_size_t len)
{
  mpzspv_elementwise (r, r_offset, x, x_offset, NULL, 0, len, 
                      SPV_ELEMENTWISE_SET);
}


void
mpzspv_reverse (mpzspv_handle_t r, const spv_size_t r_offset, 
    const mpzspv_handle_t x, const spv_size_t x_offset, const spv_size_t len)
{
  unsigned int i;
  
  if (r->storage == 1 || x->storage == 1)
    {
      /* Not implemented yet */
      abort();
    }
  
  ASSERT (mpzspv_verify_out (r, r_offset, len));
  ASSERT (mpzspv_verify_in (x, x_offset, len));
  ASSERT_ALWAYS (r->mpzspm == x->mpzspm);

  for (i = 0; i < x->mpzspm->sp_num; i++)
    spv_rev (r->mem[i] + r_offset, x->mem[i] + x_offset, len);
}


void
mpzspv_set_sp (mpzspv_handle_t r, const spv_size_t offset, 
    const sp_t c, const spv_size_t len)
{
  unsigned int i;
  
  for (i = 0; i < r->mpzspm->sp_num; i++)
    spv_elementwise (get_mem (r, i), get_file (r, i), adjust_offset (r, i, offset), 
                     &c, NULL, 0, NULL, NULL, 0, 
                     r->mpzspm->spm[i]->sp, r->mpzspm->spm[i]->mul_c,
                     len, SPV_ELEMENTWISE_SETSP);
}

void
mpzspv_neg (mpzspv_handle_t r, const spv_size_t r_offset, 
    const mpzspv_handle_t x, const spv_size_t x_offset, const spv_size_t len)
{
  mpzspv_elementwise (r, r_offset, x, x_offset, NULL, 0, len, 
                      SPV_ELEMENTWISE_NEG);
}

void
mpzspv_add (mpzspv_handle_t r, const spv_size_t r_offset, 
            const mpzspv_handle_t x, const spv_size_t x_offset, 
            const mpzspv_handle_t y, const spv_size_t y_offset, 
            const spv_size_t len)
{
  mpzspv_elementwise (r, r_offset, x, x_offset, y, y_offset, len, 
                      SPV_ELEMENTWISE_ADD);
}


void
mpzspv_sub (mpzspv_handle_t r, const spv_size_t r_offset, 
            const mpzspv_handle_t x, const spv_size_t x_offset, 
            const mpzspv_handle_t y, const spv_size_t y_offset, 
            const spv_size_t len)
{
  mpzspv_elementwise (r, r_offset, x, x_offset, y, y_offset, len, 
                      SPV_ELEMENTWISE_SUB);
}


void
mpzspv_add_sp (mpzspv_handle_t r, const spv_size_t r_offset, 
               const mpzspv_handle_t x, const spv_size_t x_offset, 
               const sp_t y, const spv_size_t len)
{
  unsigned int i;
  
  for (i = 0; i < r->mpzspm->sp_num; i++)
    spv_elementwise (get_mem (r, i), get_file (r, i), adjust_offset (r, i, r_offset), 
                     get_mem (x, i), get_file (x, i), adjust_offset (x, i, x_offset), 
                     &y, NULL, 0, 
                     r->mpzspm->spm[i]->sp, r->mpzspm->spm[i]->mul_c,
                     len, SPV_ELEMENTWISE_ADDSP);
}


void
mpzspv_sub_sp (mpzspv_handle_t r, const spv_size_t r_offset, 
               const mpzspv_handle_t x, const spv_size_t x_offset, 
               const sp_t y, const spv_size_t len)
{
  unsigned int i;
  
  for (i = 0; i < r->mpzspm->sp_num; i++)
    spv_elementwise (get_mem (r, i), get_file (r, i), adjust_offset (r, i, r_offset), 
                     get_mem (x, i), get_file (x, i), adjust_offset (x, i, x_offset), 
                     &y, NULL, 0, 
                     r->mpzspm->spm[i]->sp, r->mpzspm->spm[i]->mul_c,
                     len, SPV_ELEMENTWISE_SUBSP);
}


void
mpzspv_add_mpz (mpzspv_handle_t r, const spv_size_t r_offset, 
                const mpzspv_handle_t x, const spv_size_t x_offset, 
                const mpz_t y, const spv_size_t len)
{
  unsigned int i;
  mpz_t rem;

  mpz_init (rem);
  for (i = 0; i < r->mpzspm->sp_num; i++)
    {
      const spm_t spm = r->mpzspm->spm[i];
      const sp_t s = mpz_mod_spm (y, spm, rem);
      spv_elementwise (get_mem (r, i), get_file (r, i), adjust_offset (r, i, r_offset), 
                       get_mem (x, i), get_file (x, i), adjust_offset (x, i, x_offset), 
                       &s, NULL, 0, 
                       spm->sp, spm->mul_c, len, SPV_ELEMENTWISE_ADDSP);
    }
  mpz_clear (rem);
}


void
mpzspv_random (mpzspv_handle_t x, const spv_size_t offset, 
               const spv_size_t len)
{
  mpzspv_elementwise (x, offset, NULL, 0, NULL, 0, len, 
                      SPV_ELEMENTWISE_RANDOM);
}


static inline sp_t 
mpz_mod_spm (const mpz_t mpz, const spm_t spm, ATTRIBUTE_UNUSED mpz_t rem)
{
  sp_t r;
#if SP_TYPE_BITS > GMP_LIMB_BITS
  mpz_tdiv_r(rem, mpz, spm->mp_sp);
  r = mpz_get_sp(rem);
#else
#if defined(HAVE___GMPN_MOD_1S_4P_CPS) && defined(HAVE___GMPN_MOD_1S_4P)
#ifdef USE_VALGRIND
  VALGRIND_CHECK_MEM_IS_DEFINED(spm->cps, sizeof(spm->cps));
#endif
  ASSERT (mpz_sgn (mpz) > 0);
  /* __gmpn_mod_1s_4p() does not like input value of zero */
  /* ASSERT (spm->sp << spm->cps[1] >= (mp_limb_t)1 << (8*sizeof(mp_limb_t) - 1)); */ 
  r = __gmpn_mod_1s_4p (PTR(mpz), SIZ(mpz), 
      (mp_limb_t) spm->sp << spm->cps[1], spm->cps);
  /* ASSERT (r == mpn_mod_1 (PTR(mpz), SIZ(mpz), (mp_limb_t) spm->sp)); */
#else
  ASSERT (mpz_sgn (mpz) >= 0);
  r = mpn_mod_1 (PTR(mpz), SIZ(mpz), (mp_limb_t) spm->sp);
#endif
#endif
  return r;
}


/* convert mpz to CRT representation, naive version */
static inline void
spv_from_mpz_slow (spv_t x, const spv_size_t offset, 
                   const size_t sp_per_line, const mpz_t mpz, 
                   const mpzspm_t mpzspm, ATTRIBUTE_UNUSED mpz_t rem)
{
  unsigned int j;

  valgrind_check_mpzin (mpz);
  if (UNLIKELY(mpz_sgn (mpz) == 0))
    {
      for (j = 0; j < mpzspm->sp_num; j++)
        x[j * sp_per_line + offset] = 0;
    } else {
      const spm_t *spm = mpzspm->spm;
      for (j = 0; j < mpzspm->sp_num; j++)
        x[j * sp_per_line + offset] = mpz_mod_spm (mpz, spm[j], rem);
    }
}


/* convert mpzvi to CRT representation, fast version, assumes
   mpzspm->T has been precomputed (see mpzspm.c) */

/* TODO: use fast reduction with invariant modulus. */

static inline void
preinv_mod (mpz_t r, const mpz_t a, const mpz_t m, const mpz_t i,
            mpz_t *mt, const unsigned int b)
{
#if defined(REDC_PREINV_MOD)
  const mp_size_t words = (mp_size_t) b / GMP_NUMB_BITS;
  
  ASSERT (mpz_sizeinbase(i, 2) <= b);
  ASSERT (mpz_sgn (a) >= 0);
  ASSERT (mpz_sgn (i) >= 0);
  ASSERT (mpz_sgn (m) >= 0);
  ASSERT (b % GMP_NUMB_BITS == 0);
  ASSERT (ALLOC(i) >= words);
  ASSERT (ALLOC(m) >= words);
  ASSERT (ALLOC(r) >= words);

#if defined(WANT_ASSERT)
  {
    mp_size_t j;
    for (j = ABSIZ(i); j < words; j++)
      ASSERT (PTR(i)[j] == (mp_limb_t) 0);
    for (j = ABSIZ(m); j < words; j++)
      ASSERT (PTR(m)[j] == (mp_limb_t) 0);
  }
#endif

  mpz_set (mt[0], a);
  {
    mp_size_t j;
    ASSERT (ALLOC(mt[0]) >= 2*words);
    for (j = ABSIZ(mt[0]); j < 2*words; j++)
      PTR(mt[0])[j] = (mp_limb_t) 0;
  }
  __gmpn_redc_n (PTR(r), PTR(mt[0]), PTR(m), words, PTR(i));  
  SIZ(r) = words;
  MPN_NORMALIZE(PTR(r), SIZ(r));
#else
  mpz_tdiv_q_2exp (mt[0], a, b);
  mpz_mul (mt[1], mt[0], i);
  mpz_tdiv_q_2exp (mt[0], mt[1], b);
  mpz_mul (mt[1], mt[0], m);
  mpz_sub (r, a, mt[1]);
  if (0)
    printf ("preinv_mod(): s(a) = %lu, s(m) = %lu, s(i) = %lu, s(r) = %lu\n",
            mpz_sizeinbase (a, 2), mpz_sizeinbase (m, 2), mpz_sizeinbase (i, 2), 
            mpz_sizeinbase (r, 2));
#endif
}

static void
spv_from_mpz_fast (spv_t x, const spv_size_t offset, 
                   const size_t sp_per_line, const mpz_t mpz, 
                   const mpzspm_t mpzspm, ATTRIBUTE_UNUSED mpz_t rem)
{
  const int use_preinv = 0; /* Current code is just a test, needs much better 
                               implementation to be competitive */
  unsigned int i, j;
  const mpzv_t T = mpzspm->T, r = mpzspm->remainders;
  const unsigned int d = mpzspm->d;
  mpz_t mt[2];

#if TRACE_FUNC
  gmp_printf ("%s(%p, %" PRISPVSIZE ", %" PRIuSIZE ", %p, %p)\n", 
              __func__, x, offset, sp_per_line, mpz, mpzspm);
#endif
  ASSERT (d > 0);
  valgrind_check_mpzin (mpz);

  for (j = 0; use_preinv && j < 2; j++)
    mpz_init2 (mt[j], mpzspm->redcbits[0]);

  /* All 0 <= first_i <= d are possible here in principle. We would like to 
     start off with a modulus about half the size of the input number. Since 
     T[0][0] is the product of all NTT primes, it is a little larger than 
     l*N^2, where l is the maximum transform length and N is the modulus for 
     the polynomial (i.e., the number to factor).
     Thus T[1][0 ... 1] are about as large as residues (mod N), 
     and T[2][0 ... 3] are about half as large, so we start with i=2. */

  i = I0_FIRST;
  ASSERT (i <= d);
  /* First pass uses mpz as input */
  for (j = 0; j < 1U << i; j++)
    {
      const size_t n = (1U << i) - 1 + j;
      if (use_preinv)
        preinv_mod (r[j], mpz, T[n], mpzspm->preinv[n], mt, 
                    mpzspm->redcbits[n]);
      else
        mpz_tdiv_r (r[j], mpz, T[n]);
    }
  i++;

  for ( ; i <= d; i++)
    {
      for (j = 1U << i; j-- > 0; )
        {
          const size_t n = (1U << i) - 1 + j;
          if (use_preinv)
            preinv_mod (r[j], r[j/2],  T[n], mpzspm->preinv[n], mt, 
                        mpzspm->redcbits[n]);
          else
            mpz_tdiv_r (r[j], r[j/2], T[n]);
        }
    }

  for (j = 0; use_preinv && j < 2; j++)
    mpz_clear (mt[j]);

  /* Reduce each leaf node modulo individual primes */
  for (i = 0; i < 1U << d; i++)
    {
      const spm_t *spm = mpzspm->spm;
      if (UNLIKELY(mpz_sgn (r[i]) == 0))
        {
          for (j = mpzspm->start_p[i]; j < mpzspm->start_p[i + 1]; j++)
            x[j * sp_per_line + offset] = 0;
        } else {
          if (mpzspm->fixfactors != NULL)
            {
              for (j = mpzspm->start_p[i]; j < mpzspm->start_p[i + 1]; j++)
                {
                  const sp_t t = mpz_mod_spm (r[i], spm[j], rem);
                  x[j * sp_per_line + offset] = sp_mul (t, mpzspm->fixfactors[j], 
                      mpzspm->spm[j]->sp, mpzspm->spm[j]->mul_c);
                }
            }
          else
            {
              for (j = mpzspm->start_p[i]; j < mpzspm->start_p[i + 1]; j++)
                x[j * sp_per_line + offset] = mpz_mod_spm (r[i], spm[j], rem);
            }
            
        }
    }
}

/* See: Daniel J. Bernstein and Jonathan P. Sorenson,
 * Modular Exponentiation via the explicit Chinese Remainder Theorem
 *
 * memory: mpzspm->sp_num floats */


static inline void
spv_to_mpz(mpz_t res, const spv_t x, const spv_size_t offset, 
           const size_t sp_per_line, const mpzspm_t mpzspm, 
           mpz_t mt ATTRIBUTE_UNUSED)
{
  unsigned int i;
  float f = 0.5;
  valgrind_check_mpzout (res);
  mpz_set_ui (res, 0);

  for (i = 0; i < mpzspm->sp_num; i++)
    {
      const sp_t t = sp_mul (x[i * sp_per_line + offset], mpzspm->crt3[i], 
          mpzspm->spm[i]->sp, mpzspm->spm[i]->mul_c);

      if (sizeof (unsigned long) < sizeof (sp_t))
        {
          mpz_set_sp (mt, t);
          mpz_addmul (res, mpzspm->crt1[i], mt);
        } else {
          ASSERT (t <= ULONG_MAX);
          mpz_addmul_ui (res, mpzspm->crt1[i], (unsigned long) t);
        }

      f += (float) t * mpzspm->prime_recip[i];
    }

  mpz_add (res, res, mpzspm->crt2[(unsigned int) f]);
}


void
mpzspv_fromto_mpzv (mpzspv_handle_t x, const spv_size_t offset, 
    const spv_size_t len, 
    mpz_producerfunc_t producer, void * producer_state, 
    mpz_consumerfunc_t consumer, void * consumer_state)
{
  const mpzspm_t mpzspm = x->mpzspm;
  const unsigned int sp_num = mpzspm->sp_num;
  const int have_consumer = consumer != NULL || consumer_state != NULL;
  const int have_producer = producer != NULL || producer_state != NULL;
  const int read_from_x = have_consumer || !have_producer; /* We read if we convert from NTT, or if we normalize */
  const int write_to_x = have_producer || !have_consumer; /* We write if we convert to NTT, or if we normalize */
  static int tested_slow = 0, force_slow = 0, tested_blocklen = 0;
  static spv_size_t force_blocklen = 0;
  int use_slow = (mpzspm->T == NULL);
  spv_size_t block_len = 1<<16, len_done = 0, read_done = 0, buffer_offset;
  mpz_t mpz1, mpz2, mt;
  spv_t cachebuf;
  /* Process entries in chunks of two cache lines */
  const size_t sp_per_line = 2 * CACHE_LINE_SIZE / sizeof(sp_t);
#if defined(HAVE_AIO_READ)
  mpzspv_t buffer[2];
  struct aiocb **aiocb_list = NULL;
#else
  mpzspv_t buffer[1];
#endif
  unsigned int nr_buffers;
  unsigned int work_buffer = 0; /* Which of the two buffers (in case of 
                                   HAVE_AIO_READ) is used for NTT conversion,
                                   the other one is used for disk I/O */

#if TRACE_FUNC
  gmp_printf ("%s(%p, %" PRISPVSIZE ", %" PRISPVSIZE ", %p, %p, %p, %p)\n", 
              __func__, x, offset, len, producer, producer_state, consumer, 
              consumer_state);
#endif

  ASSERT (sizeof (mp_limb_t) >= sizeof (sp_t));

  if (read_from_x) {
    ASSERT(mpzspv_verify_in (x, offset, len));
  }

  if (write_to_x) {
    ASSERT(mpzspv_verify_out (x, offset, len));
  }

  mpz_init(mpz1);
  mpz_init(mpz2);
  mpz_init(mt);
  ASSERT_ALWAYS (sp_per_line > 0 && block_len % sp_per_line == 0);
  cachebuf = (spv_t) sp_aligned_malloc (sp_per_line * sizeof(sp_t) * sp_num);
  ASSERT_ALWAYS (cachebuf != NULL);

  if (have_producer && !tested_slow)
    {
      char *env = getenv ("USE_MPZSPV_FROM_MPZV_SLOW");
      if (env != NULL)
        {
          printf ("%s(): Setting force_slow=1 (use_slow was %d)\n", 
                  __func__, use_slow);
          force_slow = 1;
        }
      tested_slow = 1;
    }
  use_slow |= force_slow;

  if (mpzspv_handle_in_memory(x))
    {
      block_len = len; /* Do whole thing at once */
      buffer[0] = x->mem;
      buffer_offset = offset;
      nr_buffers = 1;
    }
  else
    {
      /* Do piecewise, using a temp buffer */
      unsigned int i;
      if (!tested_blocklen)
        {
          char *env = getenv ("MPZSPV_FROMTO_MPZV_BLOCKLEN");

          if (env != NULL)
            {
              spv_size_t b = strtoul (env, NULL, 10);
              if (b > 0)
                {
                  /* Round up to full cache line */
                  b = ((b - 1) / sp_per_line + 1) * sp_per_line;
                  printf ("%s(): Using block_len = %" PRISPVSIZE " (was %" 
                          PRISPVSIZE ")\n", __func__, b, block_len);
                  force_blocklen = b;
                }
            }
          tested_blocklen = 1;
        }
      if (force_blocklen > 0)
        block_len = force_blocklen;

      block_len = MIN(block_len, len);
      
#if defined(HAVE_AIO_READ)
      nr_buffers = 2;
#else
      nr_buffers = 1;
#endif
      for (i = 0; i < nr_buffers; i++)
        {
          buffer[i] = mpzspv_init (block_len, mpzspm);
          ASSERT_ALWAYS (buffer[i] != NULL);
        }
      buffer_offset = 0;
    }

#if defined(HAVE_AIO_READ)
  if (mpzspv_handle_on_disk(x)) 
    {
      /* Allocate aiocb array and an array of pointers to it */
      unsigned int i;
      aiocb_list = (struct aiocb **) malloc (sp_num * sizeof (struct aiocb *));
      ASSERT_ALWAYS (aiocb_list != NULL);
      /* First of the pointers points to malloc-ed memory */
      aiocb_list[0] = (struct aiocb *) malloc (sp_num * sizeof(struct aiocb));
      ASSERT_ALWAYS (aiocb_list[0] != NULL);
      for (i = 0; i < sp_num; i++)
        aiocb_list[i] = &aiocb_list[0][i];
    }
#endif

#if defined(HAVE_AIO_READ)
  if (read_from_x && mpzspv_handle_on_disk(x)) 
    {
      /* Read first buffer's worth of data from disk files */
      const spv_size_t read_now = MIN(len, block_len);
      if (read_now > 0)
        {
#if WANT_PROFILE
          unsigned long realstart = realtime();
#endif
          int r;
          r = mpzspv_lio_rw (aiocb_list, buffer[0], 0, x->file, x->len, 
                             offset, read_now, mpzspm, LIO_READ);
          ASSERT_ALWAYS (r == 0);
          r = mpzspv_lio_suspend (aiocb_list, mpzspm);
          ASSERT_ALWAYS (r == 0);
          read_done += read_now;
#if WANT_PROFILE
          printf("%s(): read files from position %" PRISPVSIZE 
                 " started at %lu took %lu ms\n", 
                 __func__, offset, realstart, realtime() - realstart);
#endif
        }
    }
#endif

  while (len_done < len)
    {
      const spv_size_t len_now = MIN(len - len_done, block_len);
      const spv_size_t read_now = MIN(len - read_done, block_len);
      spv_size_t i;
#if WANT_PROFILE
      unsigned long realstart;
#endif

      /* Read x from disk files */
      if (read_from_x && mpzspv_handle_on_disk(x) && read_now > 0) 
        {
#if WANT_PROFILE
          unsigned long realstart = realtime();
#endif
#if defined(HAVE_AIO_READ)
          int r;
          r = mpzspv_lio_rw (aiocb_list, buffer[work_buffer ^ 1], 0, x->file, 
                             x->len, offset + read_done, read_now, mpzspm, 
                             LIO_READ);
          ASSERT_ALWAYS (r == 0);
#else
          mpzspv_seek_and_read (buffer[0], 0, x->file, x->len, 
                                offset + read_done, read_now, mpzspm);
#endif
#if WANT_PROFILE
          printf("%s(): scheduling read files from position %" PRISPVSIZE 
                 " started at %lu took %lu ms\n", 
                 __func__, offset + read_done, realstart, 
                 realtime() - realstart);
#endif
        }

      /* Do the conversion */
#if WANT_PROFILE
      realstart = realtime();
#endif
      for (i = 0; i < len_now; i += sp_per_line)
        {
          const unsigned int words_now = MIN(len_now - i, sp_per_line);
          unsigned int j;

          if (have_consumer)
            {
              /* Copy full cache lines from disk buffer/mpzspv_t to cachebuf */
              unsigned int k;
              for (k = 0; k < sp_num; k++)
                {
                  if (words_now == sp_per_line)
                    {
                      /* Hope for good unrolling here */
                      for (j = 0; j < sp_per_line; j++)
                        cachebuf[k * sp_per_line + j] = buffer[work_buffer][k][buffer_offset + i + j];
                    } else {
                      for (j = 0; j < words_now; j++)
                        cachebuf[k * sp_per_line + j] = buffer[work_buffer][k][buffer_offset + i + j];
                    }
                }
            }

          for (j = 0; j < words_now; j++)
            {
              if (have_producer)
                {
                  if (producer != NULL)
                    {
                      /* Get new mpz1 from producer */
                      (*producer)(producer_state, mpz1);
                      valgrind_check_mpzin (mpz1);
                    } else {
                      /* Get new mpz1 from mpzv_t */
                      valgrind_check_mpzin (((mpzv_t)producer_state)[len_done + i + j]);
                      mpz_set (mpz1, ((mpzv_t)producer_state)[len_done + i + j]);
                    }
                }
              
              if (have_consumer)
                {
                  /* Convert NTT entry to mpz2 */
                  spv_to_mpz (mpz2, cachebuf, j, sp_per_line, mpzspm, mt);
                  valgrind_check_mpzin (mpz2);
                  if (consumer != NULL)
                    {
                      /* Give mpz2 to consumer function */
                      mpz_mod (mpz2, mpz2, mpzspm->modulus);
                      (*consumer)(consumer_state, mpz2);
                    } else {
                      mpz_mod (((mpzv_t)consumer_state)[len_done + i + j], mpz2, 
                               mpzspm->modulus);
                    }
                }
              
              if (have_producer)
                {
                  /* Convert the mpz1 we got from producer to NTT */
                  if (use_slow)
                    spv_from_mpz_slow (cachebuf, j, sp_per_line, mpz1, 
                                       mpzspm, mt);
                  else
                    spv_from_mpz_fast (cachebuf, j, sp_per_line, mpz1, 
                                       mpzspm, mt);
                }
            }

          if (have_producer)
            {
              /* Copy full cache lines from cachebuf to disk buffer/mpzspv_t. 
                 TODO: use non-temporal stores that bypass level 1 cache */
              unsigned int k;
              for (k = 0; k < sp_num; k++)
                {
                  if (words_now == sp_per_line)
                    {
                      /* Hope for good unrolling here */
                      for (j = 0; j < sp_per_line; j++)
                        buffer[work_buffer][k][buffer_offset + i + j] = 
                          cachebuf[k * sp_per_line + j];
                    } else {
                      for (j = 0; j < words_now; j++)
                        buffer[work_buffer][k][buffer_offset + i + j] = 
                          cachebuf[k * sp_per_line + j];
                    }
                }
            }
        }

      /* If we have neither producer or consumer, we normalise the NTT vectors */
      if (!have_producer && !have_consumer)
        {
          mpzspv_normalise (buffer[work_buffer], buffer_offset, len_now, mpzspm);
        }

#if WANT_PROFILE
    printf("%s(): processing buffer started at %lu took %lu ms\n", 
           __func__, realstart, realtime() - realstart);
#endif

      if (read_from_x && mpzspv_handle_on_disk(x) && read_now > 0) 
        {
          /* Wait for read to complete */
#if WANT_PROFILE
          unsigned long realstart = realtime();
#endif
#if defined(HAVE_AIO_READ)
          int r;
          r = mpzspv_lio_suspend (aiocb_list, mpzspm);
          ASSERT_ALWAYS (r == 0);
#endif
#if WANT_PROFILE
          printf("%s(): suspend of read files from position %" PRISPVSIZE 
                 " started at %lu took %lu ms\n", 
                 __func__, offset + read_done, realstart, realtime() - realstart);
#endif
          read_done += read_now;
        }

    /* Write current buffer to disk files */
    if (write_to_x && mpzspv_handle_on_disk(x)) {
#if WANT_PROFILE
      unsigned long realstart = realtime();
#endif
#if defined(HAVE_AIO_READ)
      int r;
      r = mpzspv_lio_rw (aiocb_list, buffer[work_buffer], 0, x->file, 
                         x->len, offset + len_done, len_now, mpzspm, 
                         LIO_WRITE);
      ASSERT_ALWAYS (r == 0);
      r = mpzspv_lio_suspend (aiocb_list, mpzspm);
      ASSERT_ALWAYS (r == 0);
#else
      mpzspv_seek_and_write (buffer[work_buffer], 0, x->file, x->len, 
                             offset + len_done, len_now, mpzspm);
#endif
#if WANT_PROFILE
      printf("%s(): write files at position %" PRISPVSIZE 
             " started at %lu took %lu ms\n", 
             __func__, offset + len_done, realstart, realtime() - realstart);
#endif
    }
    len_done += len_now;

    /* Toggle between the two buffers */
    if (nr_buffers == 2) 
      work_buffer ^= 1;

    /* If we write NTT data to memory, we need to advance the offset to fill 
       the entire array. If we use a temp buffer, we reuse the same buffer 
       each time */
    if (mpzspv_handle_in_memory(x))
      buffer_offset += len_now;
  }
  mpz_clear(mpz1);
  mpz_clear(mpz2);
  mpz_clear(mt);
  sp_aligned_free(cachebuf);

  if (!mpzspv_handle_in_memory(x))
    {
      unsigned int i;
      for (i = 0; i < nr_buffers; i++)
        {
          mpzspv_clear (buffer[i], mpzspm);
          buffer[i] = NULL;
        }
    }

#if defined(HAVE_AIO_READ)
  if (mpzspv_handle_on_disk(x)) 
    {
      free (aiocb_list[0]);
      free (aiocb_list);
    }
#endif
}


/* B&S: ecrt mod m mod p_j.
 *
 * memory: MPZSPV_NORMALISE_STRIDE mpzspv coeffs
 *         6 * MPZSPV_NORMALISE_STRIDE sp's
 *         MPZSPV_NORMALISE_STRIDE floats */
static void
mpzspv_normalise (mpzspv_t x, const spv_size_t offset, 
                  const spv_size_t len, mpzspm_t mpzspm)
{
  unsigned int i, j, sp_num = mpzspm->sp_num;
  spv_size_t k, l;
  sp_t v;
  spv_t s, d, w;
  spm_t *spm = mpzspm->spm;
  
  float *f;
  mpzspv_t t;

#ifdef WANT_ASSERT
  for (i = 0; i < mpzspm->sp_num; i++)
    {
      ASSERT(spv_verify_in (x[i] + offset, len, mpzspm->spm[i]->sp));
    }
#endif  
  f = (float *) malloc (MPZSPV_NORMALISE_STRIDE * sizeof (float));
  s = (spv_t) malloc (3 * MPZSPV_NORMALISE_STRIDE * sizeof (sp_t));
  d = (spv_t) malloc (3 * MPZSPV_NORMALISE_STRIDE * sizeof (sp_t));
  if (f == NULL || s == NULL || d == NULL)
    {
      fprintf (stderr, "%s(): Cannot allocate memory\n", __func__);
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
          for (k = 0; k < stride; k++)
	    {
	      x[i][l + k + offset] = sp_mul (x[i][l + k + offset],
	          mpzspm->crt3[i], spm[i]->sp, spm[i]->mul_c);
	      f[k] += (float) x[i][l + k + offset] * mpzspm->prime_recip[i];
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

      for (i = 0; i < sp_num; i++)
        spv_set (x[i] + l + offset, t[i], stride);
    }
  
  mpzspv_clear (t, mpzspm);
  
  free (s);
  free (d);
  free (f);
}

#ifndef HAVE_AIO_READ
static void 
mpzspv_seek_and_read (mpzspv_t dst, const size_t offset, FILE *sp_file, 
                      const spv_size_t veclen, const spv_size_t fileoffset, 
                      const size_t nread, mpzspm_t mpzspm)
{
  unsigned int j;
  for (j = 0; j < mpzspm->sp_num; j++)
  {
    spv_seek_and_read (dst[j] + offset, nread, fileoffset + j * veclen, sp_file);
  }
}


static void  
mpzspv_seek_and_write (mpzspv_t src, const size_t offset, FILE *sp_file, 
                       const spv_size_t veclen, const spv_size_t fileoffset, 
                       const size_t nwrite, const mpzspm_t mpzspm)
{
  unsigned int j;
  for (j = 0; j < mpzspm->sp_num; j++)
  {
    spv_seek_and_write (src[j] + offset, nwrite, fileoffset + j * veclen, sp_file);
  }
}
#endif

static void
mul_dct_file (const spv_t r, const spv_t spv, FILE *dct_file, const spv_size_t offset,
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
      fprintf (stderr, "%s(): could not allocate memory\n", __func__);
      abort();
    }

  while (nr_read < dctlen)
    {
      const spv_size_t read_now = MIN(dctlen - nr_read, blocklen);
      const spv_size_t mul_now = MIN(dctlen - nr_read - 1, blocklen);
      
      spv_seek_and_read (tmp, read_now, nr_read + offset, dct_file);
      
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

static inline void
profile_start (unsigned long *realstart ATTRIBUTE_UNUSED)
{
#if WANT_PROFILE
  *realstart = realtime();
#endif
}

static inline void
profile_end (const unsigned long realstart ATTRIBUTE_UNUSED, 
             const int i ATTRIBUTE_UNUSED, const char *msg ATTRIBUTE_UNUSED)
{
#if WANT_PROFILE
  printf("mpzspv_mul_ntt(): %s %d started at %lu took %lu ms\n", 
         msg, i, realstart, realtime() - realstart);
#endif
}


/* Do multiplication via NTT. Depending on the value of "steps", does 
   forward transform of, pair-wise multiplication, inverse transform. 
   Input and output spv_t's can be stored in files. 
   It is permissible to let any combination of x, y, and r point at the same 
   memory or files. */

void
mpzspv_mul_ntt (mpzspv_handle_t r, const spv_size_t offsetr, 
    mpzspv_handle_t x, const spv_size_t offsetx, const spv_size_t lenx, 
    mpzspv_handle_t y, const spv_size_t offsety, const spv_size_t leny, 
    const spv_size_t ntt_size, const int steps)
{
  const spv_size_t block_len = 16384;
  const spv_size_t log2_ntt_size = ceil_log_2 (ntt_size);
  const int do_fft1 = (steps & NTT_MUL_STEP_FFT1) != 0;
  const int do_pwmul = (steps & NTT_MUL_STEP_MUL) != 0;
  const int do_pwmul_dct = (steps & NTT_MUL_STEP_MULDCT) != 0;
  const int do_ifft = (steps & NTT_MUL_STEP_IFFT) != 0;
  int i;
  mpzspm_t mpzspm = NULL;

  /* Check that the inputs/outputs all use the same NTT definition */
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

  if (do_pwmul && do_pwmul_dct)
    {
      fprintf (stderr, "%s(): Error, both PWMUL and PWMULDCT requested\n",
               __func__);
      abort();
    }
  
  if (mpzspv_handle_in_memory(x))
    {
      ASSERT (mpzspv_verify_in (x, offsetx, lenx));
    }
  if (mpzspv_handle_in_memory(y)) 
    {
      ASSERT (mpzspv_verify_in (y, offsety, leny));
    }
  if (mpzspv_handle_in_memory(r))
    {
      ASSERT (mpzspv_verify_out (r, offsetr, ntt_size));
    }

#if TRACE_ntt_mul
  printf ("%s (r = {%d, %p, %p, %p}, offsetr = %" PRISPVSIZE ", "
          "x = {%d, %p, %p, %p}, offsetx = %" PRISPVSIZE ", lenx = %" PRISPVSIZE ", "
          "y = {%d, %p, %p, %p}, offsety = %" PRISPVSIZE ", leny = %" PRISPVSIZE ", "
          "ntt_size = %" PRISPVSIZE ", steps = %d)\n", __func__, 
          r ? r->storage : 0, r ? r->mpzspm : NULL, r ? r->mem : NULL, r ? r->file : NULL, offsetr, 
          x ? x->storage : 0, x ? x->mpzspm : NULL, x ? x->mem : NULL, x ? x->file : NULL, offsetx, lenx,
          y ? y->storage : 0, y ? y->mpzspm : NULL, y ? y->mem : NULL, y ? y->file : NULL, offsety, leny, 
          ntt_size, steps);
  if (x != NULL)
    mpzspv_print (x, offsetx, lenx, "x");
  if (y != NULL)
    mpzspv_print (y, offsety, leny, "y");
#endif
  
  /* Need parallelization at higher level (e.g., handling a branch of the 
     product tree in one thread) to make this worthwhile for ECM */

#if defined(_OPENMP)
#pragma omp parallel if (ntt_size > 32768)
#endif
  {

#if defined(_OPENMP)
#pragma omp for
#endif
  for (i = 0; i < (int) mpzspm->sp_num; i++)
    {
      const spm_t spm = mpzspm->spm[i];
      const spv_t spvx = mpzspv_handle_in_memory(x) ? x->mem[i] + offsetx : NULL;
      const spv_t spvy = mpzspv_handle_in_memory(y) ? y->mem[i] + offsety : NULL;
      const spv_t spvr = mpzspv_handle_in_memory(r) ? r->mem[i] + offsetr : NULL;
      spv_t tmp = NULL;
      unsigned long realstart;

      /* If we do any arithmetic, we need some memory to do it in. 
         If r is in memory, we can use that as temp storage, so long as we 
         don't overwrite input data we still need. */

      /* This test does not check whether r+offsetr and y+offsety point to
         non-overlapping memory; it simply takes r==y as not allowing r for
         temp space */
      if (mpzspv_handle_in_memory(r) && !(mpzspv_handle_in_memory(y) && r->mem == y->mem))
        tmp = spvr;
      else
        {
          tmp = (spv_t) sp_aligned_malloc (ntt_size * sizeof (sp_t));
          if (tmp == NULL)
            {
              fprintf (stderr, "%s(): Cannot allocate tmp memory\n",
                       __func__);
              abort();
            }
        }

      /* If we do any arithmetic, read the data of x into tmp and do any 
         wrap-around */
      if (do_fft1 || do_pwmul || do_pwmul_dct || do_ifft)
        {
          ASSERT_ALWAYS(x != NULL && r != NULL);
          if (mpzspv_handle_in_memory(x))
            {
              spv_size_t j;
              if (tmp != spvx)
                spv_set (tmp, spvx, MIN(ntt_size, lenx));
              for (j = ntt_size; j < lenx; j += ntt_size)
                {
                  spv_size_t len_now = MIN(lenx - j, ntt_size);
                  spv_add (tmp, tmp, spvx + j, len_now, spm->sp);
                }
            }
          else 
            {
              spv_size_t done = ntt_size;
              profile_start (&realstart);
              spv_seek_and_read (tmp, MIN(ntt_size, lenx), offsetx + i * x->len, 
                                 x->file);
              while (done < lenx)
                {
                  spv_size_t do_now = MIN (ntt_size, lenx - done);
                  spv_elementwise (tmp, NULL, 0, tmp, NULL, 0, NULL, x->file, 
                      adjust_offset(x, i, offsetx + done), spm->sp, spm->mul_c, 
                      do_now, SPV_ELEMENTWISE_ADD);
                  done += do_now;
                }
              profile_end (realstart, i, "read vector");
            } 

          if (ntt_size > lenx)
            spv_set_zero (tmp + lenx, ntt_size - lenx);
        }
      
      if (do_fft1) 
        {
          profile_start (&realstart);

          spv_ntt_gfp_dif (tmp, log2_ntt_size, spm);
          profile_end (realstart, i, "fft on vector");
        }

      if (do_pwmul) 
        {
          profile_start (&realstart);
          ASSERT_ALWAYS(y != NULL);
          ASSERT_ALWAYS(leny == ntt_size);
          spv_elementwise (tmp, NULL, 0, tmp, NULL, 0, get_mem(y, i), 
                           get_file(y, i), adjust_offset(y, i, offsety), 
                           spm->sp, spm->mul_c, ntt_size, SPV_ELEMENTWISE_MUL);
          profile_end (realstart, i, "pwmul on vector");
        }
      else if (do_pwmul_dct)
        {
          ASSERT_ALWAYS(y != NULL);
          ASSERT_ALWAYS(leny == ntt_size / 2 + 1);
          profile_start (&realstart);
          if (mpzspv_handle_in_memory(y))
            mul_dct (tmp, tmp, spvy, ntt_size, spm);
          else
            mul_dct_file (tmp, tmp, y->file, i * y->len, ntt_size, block_len, spm);
          profile_end (realstart, i, "pwmuldct on vector");
        }

      if (do_ifft) 
        {
          ASSERT (sizeof (mp_limb_t) >= sizeof (sp_t));

          profile_start (&realstart);
          spv_ntt_gfp_dit (tmp, log2_ntt_size, spm);

          /* spm->sp - (spm->sp - 1) / ntt_size is the inverse of ntt_size */
          spv_mul_sp (tmp, tmp, spm->sp - (spm->sp - 1) / ntt_size,
                      ntt_size, spm->sp, spm->mul_c);

          profile_end (realstart, i, "ifft on vector");
        }

      if (do_fft1 || do_pwmul || do_pwmul_dct || do_ifft)
        {
          if (mpzspv_handle_in_memory(r))
            {
              if (tmp != spvr)
                spv_set (spvr, tmp, ntt_size);
            }
          else
            {
              profile_start (&realstart);
              spv_seek_and_write (tmp, ntt_size, offsetr + i * r->len, r->file);
              profile_end (realstart, i, "write of vector");
            }

          if (tmp != spvr)
            sp_aligned_free (tmp);
        }
    }
  }

#if TRACE_ntt_mul
  if (r != NULL)
    mpzspv_print (r, offsetr, ntt_size, "r");
#endif
}


/* Computes a DCT-I of length dctlen. Input is the spvlen coefficients
   in spv. FIXME: handle wrap-around in input data */

void
mpzspv_to_dct1 (mpzspv_handle_t dct, const mpzspv_handle_t spv,  
                const spv_size_t spvlen, const spv_size_t dctlen)
{
  const spv_size_t ntt_size = 2 * (dctlen - 1); /* Length for the DFT */
  const spv_size_t log2_l = ceil_log_2 (ntt_size);
  int j;

  ASSERT_ALWAYS (dct->mpzspm == spv->mpzspm);
  ASSERT (mpzspv_verify_out (dct, 0, dctlen));
  ASSERT (mpzspv_verify_in (spv, 0, spvlen));

#ifdef _OPENMP
#pragma omp parallel
#endif
  {

#ifdef _OPENMP
#pragma omp for
#endif
  for (j = 0; j < (int) spv->mpzspm->sp_num; j++)
    {
      const spm_t spm = spv->mpzspm->spm[j];
      spv_size_t i;
      
      spv_t tmp = (spv_t) sp_aligned_malloc (ntt_size * sizeof (sp_t));
      if (tmp == NULL)
        {
          fprintf (stderr, "%s(): Cannot allocate tmp memory in\n", __func__);
          abort();
        }

      if (mpzspv_handle_on_disk(spv))
        {
          spv_seek_and_read (tmp, spvlen, j * spv->len, spv->file);
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
      printf ("%s: tmp[%d] = [", __func__, j);
      for (i = 0; i < ntt_size; i++)
          printf ("%lu, ", tmp[i]);
      printf ("]\n");
#endif
      
      spv_ntt_gfp_dif (tmp, log2_l, spm);

#if 0
      printf ("%s: tmp[%d] = [", __func__, j);
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
            printf ("%s: DFT[%lu] == DFT[%lu]\n", __func__, i, m - i);
#endif
          }
      }
#endif

      /* Copy coefficients to dct buffer */
      {
        spv_t out_buf = mpzspv_handle_in_memory(dct) ? dct->mem[j] : tmp;
        const sp_t coeff_1 = tmp[1];
        for (i = 0; i < dctlen - 1; i++)
          out_buf[i] = tmp[i * 2];
        out_buf[dctlen - 1] = coeff_1;
        if (mpzspv_handle_on_disk(dct))
          {
            /* Write data back to file */
            spv_seek_and_write (tmp, dctlen, j * dct->len, dct->file);
          }
      }

      sp_aligned_free(tmp);
    }
  }
}


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
  printf ("%s: NTT vector mod %lu\n", __func__, sp);
  spv_print_vec (spv, sp, n, "before weighting:", "\n");
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
  spv_print_vec (spv, sp, n, "after weighting:", "\n");
#endif

  /* Forward DFT of dft[j] */
  spv_ntt_gfp_dif (spv, log2_len, spm);

#if TRACE_ntt_sqr_reciprocal
  spv_print_vec (spv, sp, n, "after forward transform:", "\n");
#endif

  /* Square the transformed vector point-wise */
  spv_pwmul (spv, spv, spv, len, sp, mul_c);

#if TRACE_ntt_sqr_reciprocal
  spv_print_vec (spv, sp, n, "after point-wise squaring:", "\n");
#endif

  /* Inverse transform of dft[j] */
  spv_ntt_gfp_dit (spv, log2_len, spm);

#if TRACE_ntt_sqr_reciprocal
  spv_print_vec (spv, sp, n, "after inverse transform:", "\n");
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
  spv_print_vec (spv, sp, n, "after after un-weighting:", "\n");
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
  spv_print_vec (spv, sp, n, "after un-wrapping:", "\n");
#endif
}


/* Square an RLP */

void 
mpzspv_sqr_reciprocal (mpzspv_handle_t x, const spv_size_t n)
{
  const spv_size_t log2_n = ceil_log_2 (n);
  const spv_size_t len = ((spv_size_t) 2) << log2_n;

  ASSERT(x->mpzspm->max_ntt_size % 3UL == 0UL);
  ASSERT(len % 3UL != 0UL);
  ASSERT(x->mpzspm->max_ntt_size % len == 0UL);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    int j;

#ifdef _OPENMP
#pragma omp for
#endif
    for (j = 0; j < (int) (x->mpzspm->sp_num); j++)
      {
        spv_t tmp;
        
        if (mpzspv_handle_on_disk(x))
          {
            tmp = (spv_t) sp_aligned_malloc (len * sizeof (sp_t));

            if (tmp == NULL)
              {
                fprintf (stderr, "%s(): Cannot allocate tmp memory\n", __func__);
                abort();
              }
            spv_seek_and_read (tmp, n, j * x->len, x->file);
          }
        else
          {
            tmp = x->mem[j];
          }
        
        spv_sqr_reciprocal (n, x->mpzspm->spm[j], tmp, x->mpzspm->max_ntt_size);
        
        if (mpzspv_handle_on_disk(x))
          {
            spv_seek_and_write (tmp, 2 * n - 1, j * x->len, x->file);
            sp_aligned_free (tmp);
          }
      }
    }
}


#if defined(HAVE_AIO_READ)
/* If write=0, read data from a set of files and return as soon as the 
   reads are scheduled 
   - OR - 
   if write=1, write data to a set of files and wait until writes are 
   completed (which usually returns as soon as the writes are in the 
   system's disk write cache).
   In both cases, reads/writes "len" entries from/to file position 
   "file_offset" to/from position "mpzspv_offset" in "mpzspv". All 
   lengths/positions use one sp_t as the unit. */
static int 
mpzspv_lio_rw (struct aiocb *aiocb_list[], mpzspv_t mpzspv, 
               const size_t mpzspv_offset, FILE *file, 
               const spv_size_t veclen, const spv_size_t file_offset, 
               const size_t len, const mpzspm_t mpzspm, const int opcode)
{
  unsigned int i;
  struct sigevent sev;
  int r;
  int fn;
  
  ASSERT_ALWAYS (opcode == LIO_READ || opcode == LIO_WRITE || 
                 opcode == LIO_NOP);

  if (0)
    printf("%s(, , %" PRIuSIZE ", , %" PRISPVSIZE ", %" PRIuSIZE ", , %d)\n", 
           __func__, mpzspv_offset, file_offset, len, opcode);
  
  memset (&sev, 0, sizeof(struct sigevent));
  sev.sigev_notify = SIGEV_NONE;
  fn = fileno (file);
  ASSERT_ALWAYS (fn != -1);
  for (i = 0; i < mpzspm->sp_num; i++)
    {
      memset (aiocb_list[i], 0, sizeof (struct aiocb));
      aiocb_list[i]->aio_fildes = fn;
      aiocb_list[i]->aio_offset = (file_offset + i * veclen) * (off_t) sizeof(sp_t);
      aiocb_list[i]->aio_buf = mpzspv[i] + mpzspv_offset;
      aiocb_list[i]->aio_nbytes = len * sizeof(sp_t);
      aiocb_list[i]->aio_reqprio = 0;
      aiocb_list[i]->aio_sigevent = sev;
      aiocb_list[i]->aio_lio_opcode = opcode;
    }
  r = lio_listio (LIO_NOWAIT, aiocb_list, mpzspm->sp_num, NULL);
  return r;
}

/* Wait until all operations in aiocb_list[] have completed */
static int 
mpzspv_lio_suspend (struct aiocb * const aiocb_list[], const mpzspm_t mpzspm)
{
  unsigned int i;

  for (i = 0; i < mpzspm->sp_num; i++)
    {
      /* To avoid warning about mismatch in const-ness */
      const struct aiocb *const p = aiocb_list[i];
      int r, err;
      ssize_t ret;

      r = aio_suspend (&p, 1, NULL);
      while (r != 0)
        {
          if (errno == EINTR)
            {
              err = aio_error (aiocb_list[i]);
              if (err == 0 || err == EINPROGRESS)
                {
                  /* This should not happen if aio_suspend() set 
                     errno == EINTR */
                  fprintf (stderr, "%s(): Error, aio_suspend() returned set "
                           "errno = EINTR but aio_error() returned %d\n", 
                           __func__, err);
                  return r;
                }
              else if (err == ECANCELED)
                {
                  /* This transfer was interrupted. Re-issue it sychroneously.
                     FIXME: do we have to fill in the aiocb struct again? */
                  r = lio_listio (LIO_WAIT, &aiocb_list[i], 1, NULL);
                  /* The lio_listio () call with mode = LIO_WAIT sets r and 
                     errno like aio_suspend() would */
                }
              else 
                {
                  perror (NULL);
                  fprintf (stderr, "mpzspv_lio_suspend(): Error, aio_suspend() "
                           "set errno = EINTR and aio_error() returned %d\n", 
                           err);
                  return r;
                }
            }
          else
            {
              perror (NULL);
              fprintf (stderr, "%s(): aio_suspend() returned %d\n", __func__, r);
              return r;
            }
        }
        
        err = aio_error (aiocb_list[i]);
        if (err != 0)
          {
            fprintf (stderr, "%s(): Error, aio_error() returned %d\n", 
                     __func__, err);
            return -1;
          }
        
        ret = aio_return (aiocb_list[i]);
        if (ret < 0 || (size_t) ret != aiocb_list[i]->aio_nbytes)
          {
            fprintf (stderr, "%s(): Error, transfer was for %" PRIuSIZE 
                     " bytes, but aio_return() returned %" PRIdSIZE "\n", 
                     __func__, aiocb_list[i]->aio_nbytes, ret);
            return -1;
          }
    }

  return 0;
}
#endif


void
mpzspv_print (mpzspv_handle_t handle, const spv_size_t offset, 
              const spv_size_t len, const char *prefix)
{
  unsigned int i;
  spv_t tmp = NULL;
  
  if (len == 0)
    {
      printf("%s: Zero length vector\n", prefix);
      return;
    }
  
  if (mpzspv_handle_on_disk(handle))
    {
      tmp = (spv_t) sp_aligned_malloc (len * sizeof (sp_t));
      if (tmp == NULL)
        {
          fprintf (stderr, "%s(): Could not allocate memory\n", __func__);
          return;
        }
    }

  for (i = 0; i < handle->mpzspm->sp_num; i++)
    {
      if (mpzspv_handle_on_disk(handle))
        spv_seek_and_read (tmp, len, offset + i * handle->len, handle->file);
      else
        tmp = handle->mem[i] + offset;

      spv_print_vec (tmp, handle->mpzspm->spm[i]->sp, len, 
                     prefix, mpzspv_handle_in_memory(handle) ? "(in memory)\n" : "(on disk)\n");
    }

  if (mpzspv_handle_on_disk(handle))
    sp_aligned_free (tmp);
}
