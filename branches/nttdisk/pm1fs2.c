/* Implementation of fast stage 2 for P-1 and P+1 as described in
   "Improved Stage 2 to $P\pm{}1$ Factoring Algorithms" by
   Peter L. Montgomery and Alexander Kruppa, ANTS 2008 (8th Algorithmic 
   Number Theory Symposium).
   
Copyright 2007, 2008, 2009, 2010, 2011, 2012 Alexander Kruppa, Paul Zimmermann.
NTT functions are based on code Copyright 2005 Dave Newman.

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

#define _GNU_SOURCE
#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "ecm-impl.h"
#include "sp.h"
#include <math.h>
#ifdef HAVE_ALLOCA_H
#include <alloca.h>
#endif
#ifdef HAVE_STRING_H
#include <string.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef HAVE_FCNTL_H
#include <fcntl.h>
#endif
#ifdef HAVE_AIO_H
#include <aio.h>
#endif

/* TODO:
   - move functions into their proper files (i.e. NTT functions etc.)
   - later: allow storing NTT vectors on disk
*/

/* Define TEST_ZERO_RESULT to test if any result of the multipoint
   evaluation is equal to zero. If the modulus is composite, this
   happening might indicate a problem in the evalutaion code */
#define TEST_ZERO_RESULT

/* This type is the basis for file I/O of mpz_t */
typedef unsigned long file_word_t;

typedef struct {
  int storage; /* memory = 0, file = 1 */
  uint64_t len;
  size_t words; /* Number of file_word_t in a residue */
  union {
    listz_t mem;
    FILE *file;
  } data;
  char *filename;
} _listz_handle_t;
typedef _listz_handle_t *listz_handle_t;


const int pari = 0;
const int check_eval = 0;


/* Init a listz_handle_t to store up to len residues (modulo m). 
   If filename != NULL, uses disk storage, otherwise memory.
   Returns NULL if something goes wrong (i.e., if a memory allocation
   or opening a file fails) */

listz_handle_t 
listz_handle_init2 (const char *filename, const uint64_t len, const mpz_t m)
{
  listz_handle_t F;
  void *buf;

  F = malloc (sizeof (_listz_handle_t));
  if (F == NULL)
    return NULL;
  
  /* Find out how many file_word_t's  m has */
  buf = (file_word_t *) mpz_export (NULL, &F->words, -1, sizeof(file_word_t), 
                                   -1, 0, m);
  if (buf == NULL)
    {
      free (F);
      return NULL;
    }
  free(buf);

  F->len = len;
  if (filename == NULL)
    {
      F->storage = 0; /* Memory storage */
      F->data.mem = init_list2 (len, mpz_sizeinbase(m, 2));
      if (F->data.mem == NULL)
        {
          free (F);
          F = NULL;
        }
    } else {
      F->storage = 1; /* Disk storage */
      F->filename = (char *) malloc ((strlen (filename) + 1) * sizeof(char));
      if (F->filename == NULL)
        {
          free (F);
          return NULL;
        }
      strcpy (F->filename, filename);
      F->data.file = fopen (F->filename, "w+");
      if (F->data.file == NULL)
        {
          free (F->filename);
          free (F);
          return NULL;
        }
#ifdef HAVE_FALLOCATE
      fallocate (fileno(F->data.file), 0, (off_t) 0, 
                 F->words * sizeof(file_word_t) * len);
#endif
    }

  return F;
}


void 
listz_handle_clear (listz_handle_t F)
{
  if (F->storage == 0)
    {
      clear_list (F->data.mem, F->len);
      F->data.mem = NULL;
    }
  else
    {
      fclose (F->data.file);
      F->data.file = NULL;
      remove (F->filename);
      free (F->filename);
    }
  free (F);
}


static inline void 
export_residue (file_word_t *buf, const size_t bufsize, const mpz_t r)
{
  size_t nr;

  /* Export r to buf */
  mpz_export (buf, &nr, -1, sizeof(file_word_t), 0, 0, r);
  ASSERT_ALWAYS (nr <= bufsize);

  /* Pad buf with zeroes */
  for ( ; nr < bufsize; nr++)
    buf[nr] = 0;
}

static inline void 
write_residue (FILE *f, const mpz_t r, file_word_t *buf, const size_t bufsize)
{
  size_t nr;

  ASSERT_ALWAYS (mpz_sgn (r) >= 0);

  export_residue (buf, bufsize, r);
  nr = fwrite (buf, sizeof(file_word_t), bufsize, f);
  ASSERT_ALWAYS (nr == bufsize);
}

static inline void 
seek_write_residue (FILE *f, const mpz_t r, file_word_t *buf,
              const size_t bufsize, const size_t index)
{
  fseek (f, sizeof(file_word_t) * bufsize * index, SEEK_SET);
  write_residue (f, r, buf, bufsize);
}

static inline void 
read_residue (FILE *f, mpz_t r, file_word_t *buf, const size_t bufsize)
{
  size_t nr;
  
  nr = fread (buf, sizeof(file_word_t), bufsize, f);
  ASSERT_ALWAYS (nr == bufsize);
  
  mpz_import (r, bufsize, -1, sizeof(file_word_t), 0, 0, buf);
}

static inline void 
seek_read_residue (FILE *f, mpz_t r, file_word_t *buf,
              const size_t bufsize, const size_t index)
{
  fseek (f, sizeof(file_word_t) * bufsize * index, SEEK_SET);
  read_residue (f, r, buf, bufsize);
}

static void 
write_residues (FILE *f, const listz_t r, const size_t len, 
                const size_t bufsize)
{
  file_word_t *buf;
  size_t i;

  /* Let GMP allocate a buffer that is large enough for the modulus,
     hence is large enough for any residue */
  buf = (file_word_t *) malloc (bufsize * sizeof(file_word_t));
  ASSERT_ALWAYS (buf != NULL);
  
  for (i = 0; i < len; i++)
    write_residue (f, r[i], buf, bufsize);

  free (buf);
}


/* Fetches one entry from F (either in memory or file) and stores it in r. */

static inline void
listz_handle_get (listz_handle_t F, mpz_t r, file_word_t *buf, 
    const size_t index)
{
  if (F->storage == 0)
    mpz_set (r, F->data.mem[index]);
  else
    seek_read_residue (F->data.file, r, buf, F->words, index);
}

static inline void
listz_handle_get2 (listz_handle_t F, mpz_t r, const size_t index)
{
  file_word_t *buf = NULL;
  if (F->storage == 1)
    buf = malloc (F->words * sizeof (file_word_t));
  listz_handle_get (F, r, buf, index);
  free(buf);
}

/* Stores the value of r in an entry of F (either in memory or file) */

static inline void
listz_handle_set (listz_handle_t F, const mpz_t r, file_word_t *buf,
    const size_t index)
{
  if (F->storage == 0)
    mpz_set (F->data.mem[index], r);
  else
    seek_write_residue (F->data.file, r, buf, F->words, index);
}

typedef struct {
  listz_handle_t l;
  uint64_t index;
  file_word_t *buf;
} listz_handle_state_t;

static void
listz_handle_write (void *statep, const mpz_t m)
{
  listz_handle_state_t *state = statep;
  listz_handle_set (state->l, m, state->buf, state->index++);
}

static void ATTRIBUTE_UNUSED
listz_handle_read (void *statep, mpz_t m)
{
  listz_handle_state_t *state = statep;
  listz_handle_get (state->l, m, state->buf, state->index++);
}


/* Some useful PARI functions:

   V(i,X) = { if (i==0, return(2)); if (i==1, return(X)); if(i%2 == 0, return (V (i/2, X)^2-2)); return (V ((i+1)/2, X) * V ((i-1)/2, X) - X)}

   U(i,X) = { if (i==0, return(0)); if (i==1, return(1)); if(i%2 == 0, return (U (i/2, X) * V(i/2,X))); return (V ((i+1)/2, X) * U( (i-1)/2, X) + 1)}
*/


static void 
get_chunk (uint64_t *chunk_start, uint64_t *chunk_len, const uint64_t len)
{
#ifdef _OPENMP
  if(omp_in_parallel())
    {
      const int nr_chunks = omp_get_num_threads();
      const int thread_nr = omp_get_thread_num();
      uint64_t s, l;

      ASSERT_ALWAYS(nr_chunks > 0);
      if (len == 0)
        {
          *chunk_start = 0;
          *chunk_len = 0;
          return;
        }

      l = (len - 1) / nr_chunks + 1; /* l = ceil(len / nr_chunks) */
      s = thread_nr * l;
      l = MIN(l, (len > s) ? len - s : 0);
      *chunk_start = s;
      *chunk_len = l;
      return;
    }
#endif

  *chunk_start = 0;
  *chunk_len = len;
}

static void 
ntt_sqr_reciprocal (mpzv_t, const mpzv_t, const spv_size_t, const mpzspv_handle_t);

static void
print_elapsed_time (int verbosity, long cpu_start, 
		    ATTRIBUTE_UNUSED long real_start)
{
// #ifdef _OPENMP
  if (real_start != 0L)
    {
      outputf (verbosity, " took %lums (%lums real)\n", 
	       elltime (cpu_start, cputime()), 
	       elltime (real_start, realtime()));
      return;
    }
// #endif
  outputf (verbosity, " took %lums\n", elltime (cpu_start, cputime()));
}


static void
list_output_poly_file (const listz_handle_t l, uint64_t len, int monic, int symmetric,
		  char *prefix, char *suffix, int verbosity)
{
  uint64_t i;
  mpz_t m;
  file_word_t *buf = NULL;

  if (!test_verbose(verbosity))
    return;
  
  if (l->storage != 0)
    buf = (file_word_t *) malloc (l->words * sizeof(file_word_t));

  if (prefix != NULL)
    outputf (verbosity, prefix);

  if (len == 0)
    {
      if (monic)
	outputf (verbosity, "1\n");
      else
	outputf (verbosity, "0\n");
      return;
    }

  mpz_init (m);
  if (monic)
    {
      if (symmetric)
	outputf (verbosity, "(x^%" PRIu64 " + x^-%" PRIu64 ") + ", len, len);
      else
	outputf (verbosity, "x^%" PRIu64 " + ", len);
    }
  for (i = len - 1; i > 0; i--)
    {
      listz_handle_get (l, m, buf, i);
      if (symmetric)
        outputf (verbosity, "Mod(%Zd,N) * (x^%" PRIu64 " + x^-%" PRIu64 ") + ", 
                 m, i, i);
      else
        outputf (verbosity, "Mod(%Zd,N) * x^%" PRIu64 " + ", m, i);
    }
  listz_handle_get (l, m, buf, 0);
  outputf (verbosity, "Mod(%Zd,N)", m);
  if (suffix != NULL)
    outputf (verbosity, suffix);
  free (buf);
  buf = NULL;
  mpz_clear (m);
}

static void
list_output_poly (listz_t l, uint64_t len, int monic, int symmetric,
		  char *prefix, char *suffix, int verbosity)
{
  _listz_handle_t handle;
  handle.storage = 0;
  handle.len = 0;
  handle.data.mem = l;
  /* handle.words is not initialised */
  list_output_poly_file (&handle, len, monic, symmetric, prefix, suffix, 
    verbosity);
}


/* Same, but does squaring which makes things easier */

static void
list_sqr_reciprocal (listz_t R, listz_t S, const uint64_t l, 
		     mpz_t modulus, listz_t tmp, 
		     ATTRIBUTE_UNUSED const uint64_t tmplen)
{
  uint64_t i;
  listz_t Srev, r1 = tmp, r2 = tmp + 2 * l - 1, t = tmp + 4 * l - 2;

  if (l == 0UL)
    return;

  /* FIXME: This modifies the input arguments. */
  /* We have to divide S[0] by 2 */

  ASSERT (tmplen >= 4 * l - 2 + list_mul_mem (l));

#if 0
  list_output_poly (S, l, 0, 1, "/* list_sqr_reciprocal */ S(x) = ", 
                    "\n", OUTPUT_DEVVERBOSE)
#endif

  if (mpz_odd_p (S[0]))
    {
      ASSERT_ALWAYS (mpz_odd_p (modulus));
      mpz_add (S[0], S[0], modulus);
    }
  mpz_tdiv_q_2exp (S[0], S[0], 1UL);
  
  list_mul (r1, S, l, 0, S, l, 0, t);
  /* r1 = f0*g0/4 + (f0*g1 + f1*g0)/2 * x + f1*g1 * x^2 */
#if 0
  for (i = 0; i < 2UL * l - 1UL; i++)
    gmp_printf ("list_sqr_reciprocal: r1[%lu] = %Zd\n", i, r1[i]);
#endif

  Srev = (listz_t) malloc (l * sizeof (mpz_t));
  ASSERT_ALWAYS (Srev != NULL);
  for (i = 0UL; i < l; i++)
      (*Srev)[i] = (*S)[l - 1UL - i];
  list_mul (r2, S, l, 0, Srev, l, 0, t);
  /* r2 is symmetric, r2[i] = r2[2*l - 2 - i]. Check this */
#if 0
  for (i = 0; 0 && i < 2UL * l - 1UL; i++)
    gmp_printf ("list_sqr_reciprocal: r2[%lu] = %Zd\n", i, r2[i]);
#endif
#ifdef WANT_ASSERT
  for (i = 0UL; i < l; i++)
    ASSERT (mpz_cmp (r2[i], r2[2UL * l - 2UL - i]) == 0);
#endif
  free (Srev);
  /* r2 = g1*f0/2 + (g0*f0/4 + g1*f1) * x + g0*f1/2 * x^2 */
#if 0
  for (i = 0; i < 2UL * l - 1UL; i++)
    gmp_printf ("list_sqr_reciprocal: r2[%lu] = %Zd\n", i, r2[i]);
#endif

  mpz_mul_2exp (r1[0], r1[0], 1UL);
  /* r1 = f0*g0/2 + (f0*g1 + f1*g0)/2 * x + f1*g1 * x^2 */
  for (i = 0UL; i < l; i++)
    {
      mpz_mul_2exp (r2[l - i - 1UL], r2[l - i - 1UL], 1UL);
      mpz_add (R[i], r1[i], r2[l - i - 1UL]);
    }
  /* r1 = 3/4*f0*g0 + g1*f1 + (f0*g1 + 2*f1*g0)/2 * x + f1*g1 * x^2 */
  /* r1 = f0*g0 + 2*g1*f1 + (f0*g1 + f1*g0) * x + f1*g1 * x^2 */
  for (i = l; i < 2UL * l - 1UL; i++)
      mpz_set (R[i], r1[i]);

  if (R != S)
    {
      mpz_mul_2exp (S[0], S[0], 1UL);
      if (mpz_cmp(S[0], modulus) >= 0)
        mpz_sub(S[0], S[0], modulus);
      ASSERT_ALWAYS (mpz_cmp(S[0], modulus) < 0);
    }
	
#if 0
  for (i = 0; i < 2UL * l; i++)
    gmp_printf ("list_sqr_reciprocal: R[%lu] = %Zd\n", i, R[i]);
#endif
}

ATTRIBUTE_UNUSED
static void
list_recip_eval1 (mpz_t R, const listz_t S, const uint64_t l)
{
  uint64_t i;

  mpz_set_ui (R, 0UL);
  for (i = 1; i < l; i++)
    mpz_add (R, R, S[i]);
  mpz_mul_2exp (R, R, 1UL);
  if (l > 0UL)
    mpz_add (R, R, S[0]);
}

/* Multiply two reciprocal polynomials of degree 2*l1-2 and 2*l2-2, resp., 
   with coefficients in standard basis

   S_1(x) = S1[0] + sum_{1 \leq i \leq l1 - 1} S1[i] (x^i + x^{-i})
   S_2(x) = S2[0] + sum_{1 \leq i \leq l2 - 1} S2[i] (x^i + x^{-i})

   to the reciprocal polynomial of degree 2*(l1 + l2) - 4

   R(x) = R[0] + sum_{1 \leq i \leq l1 + l2 - 2} R[i] (x^i + x^{-i}) 
        = S_1(x) * S_2(x)

   R == S1 == S2 is permissible, however if S1 == S2, l1 must be equal 
   to l2 (i.e. the multiplication must be a squaring)
*/
  /* FIXME: This modifies the input arguments. */
  /* We have to divide S1[0] and S2[0] by 2 */

static void
list_mul_reciprocal (listz_t R, listz_t S1, uint64_t l1, 
		     listz_t S2, uint64_t l2,
		     mpz_t modulus, listz_t tmp, 
		     ATTRIBUTE_UNUSED const uint64_t tmplen)
{
  uint64_t i;
  const uint64_t lmax = MAX(l1, l2);
  listz_t r1 = tmp, r2 = tmp + 2*lmax - 1, rev = tmp + 4*lmax - 2,
    t = tmp + 6*lmax - 3;
#ifdef WANT_ASSERT
  mpz_t sum1, sum2, prod;
#endif

  ASSERT (S1 < tmp || S1 >= tmp + tmplen);
  ASSERT (S2 < tmp || S2 >= tmp + tmplen);
  ASSERT (R < tmp || R >= tmp + tmplen);

  if (l1 == 0 || l2 == 0)
    return;

  if (S1 == S2)
    {
      ASSERT_ALWAYS (l1 == l2);
      list_sqr_reciprocal (R, S1, l1, modulus, tmp, tmplen);
      return;
    }

  ASSERT (tmplen >= 6*lmax - 3 + list_mul_mem (lmax));
#ifdef WANT_ASSERT
  mpz_init (sum1);
  mpz_init (sum2);
  mpz_init (prod);
  list_recip_eval1 (sum1, S1, l1);
  list_recip_eval1 (sum2, S2, l2);
  mpz_mul (prod, sum1, sum2);
  mpz_mod (prod, prod, modulus);
#endif


  /* Make S1 the longer of the two, i.e. l1 >= l2 */
  if (l2 > l1)
    {
      listz_t St = S1;
      unsigned long lt = l1;
      S1 = S2;
      S2 = St;
      l1 = l2;
      l2 = lt;
    }
  
#if 0
  gmp_printf ("/* list_mul_reciprocal */ S1(x) = %Zd", S1[0]);
  for (i = 1; i < l1; i++)
    gmp_printf (" + %Zd * (x^%lu + 1/x^%lu)", S1[i], i, i);
  gmp_printf ("\n");
  gmp_printf ("/* list_mul_reciprocal */ S2(x) = %Zd", S2[0]);
  for (i = 1; i < l1; i++)
    gmp_printf (" + %Zd * (x^%lu + 1/x^%lu)", S2[i], i, i);
  gmp_printf ("\n");
#endif
  
  /* Divide S1[0] and S2[0] by 2 */
  if (mpz_odd_p (S1[0]))
    {
      ASSERT_ALWAYS (mpz_odd_p (modulus));
      mpz_add (S1[0], S1[0], modulus);
    }
  mpz_tdiv_q_2exp (S1[0], S1[0], 1UL);
  
  if (mpz_odd_p (S2[0]))
    {
      ASSERT_ALWAYS (mpz_odd_p (modulus));
      mpz_add (S2[0], S2[0], modulus);
    }
  mpz_tdiv_q_2exp (S2[0], S2[0], 1UL);

  /* Pad rev with zeros */
  for (i = l2; i < lmax; i++)
    mpz_set_ui (rev[i], 0UL);
  
  for (i = 0UL; i < l2; i++)
    mpz_set (rev[i], S2[l2 - 1UL - i]);
  list_mul (r1, S1, lmax, 0, rev, lmax, 0, t);
  /* r1 = \tilde{f}(x) \rev(\tilde{g}(x)) and has degree l1 + l2 - 2,
     i.e. l1 + l2 - 1 entries. */
#if 0
  for (i = 0; i < 2 * lmax - 1; i++)
    gmp_printf ("list_mul_reciprocal: r1[%lu] = %Zd\n", i, r1[i]);
#endif
  
  for (i = 0UL; i < l2; i++)
    mpz_set(rev[i], S2[i]);
  list_mul (r2, S1, lmax, 0, rev, lmax, 0, t);
  /* \tilde{f}(x) \tilde{g}(x) */
  
#if 0
  for (i = 0; i < 2 * lmax - 1; i++)
    gmp_printf ("list_mul_reciprocal: r2[%lu] = %Zd\n", i, r2[i]);
#endif
  
  /* Add f_0*g_0 by doubling the f_0*g_0 term in r2 */
  mpz_mul_2exp (r2[0], r2[0], 1UL);
  
  /* Add \flloor x^{-d_g} \tilde{f}(x) \rev(\tilde{g}(x)) \rfloor.
     d_g = l2 - 1. */
  for (i = 0; i < l1; i++)
    mpz_add (r2[i], r2[i], r1[i + l2 - 1]);
  
  /* Add \floor x^{-d_f} rev(\tilde{f}(x) \rev(\tilde{g}(x))) \rfloor.
     d_f = l1 - 1. rev(r2)[i] = r2[l1 + l2 - 2 - i]. We want
     rev(r2)[l1 - 1 ... l1 + l2 - 2], hence 
     r2[l2 - 1 ... 0] */
  for (i = 0; i < l2; i++)
    mpz_add (r2[i], r2[i], r1[l2 - 1 - i]);
  
#if 0
  for (i = 0; i < l1 + l2 - 1; i++)
    gmp_printf ("list_mul_reciprocal: r2[%lu] = %Zd\n", i, r2[i]);
#endif
  
  mpz_mul_2exp (S1[0], S1[0], 1UL);
  mpz_mul_2exp (S2[0], S2[0], 1UL);
  
  for (i = 0; i < l1 + l2 - 1; i++)
    mpz_set (R[i], r2[i]);
  
#if 0
  for (i = 0; i < l1 + l2 - 1; i++)
    gmp_printf ("list_mul_reciprocal: R[%lu] = %Zd\n", i, R[i]);
#endif
#ifdef WANT_ASSERT
  list_recip_eval1 (sum1, R, l1 + l2 - 1);
  mpz_mod (sum1, sum1, modulus);
  ASSERT (mpz_cmp (prod, sum1) == 0);
  mpz_clear (sum1);
  mpz_clear (sum2);
  mpz_clear (prod);
#endif
}


/* 
  Computes V_k(S), where the Chebyshev polynomial V_k(X) is defined by 
  V_k(X + 1/X) = X^k + 1/X^k
*/

static void
V (mpres_t R, const mpres_t S, const int64_t k, mpmod_t modulus)
{
  mpres_t V0, Vi, Vi1;
  uint64_t j, uk;
  int po2;

  if (test_verbose(OUTPUT_TRACE))
    {
      mpz_t tz;
      mpz_init (tz);
      mpres_get_z (tz, S, modulus);
      gmp_printf ("\nChebyshev_V(%ld, Mod(%Zd,N)) == ", (long)k, tz);
      mpz_clear (tz);
    }

  if (k == 0)
    {
      mpres_set_ui (R, 2UL, modulus);
      goto exit;
    }

  uk = (k >= 0) ? k : -k;

  for (po2 = 0; uk % 2 == 0; uk >>= 1, po2++);

  mpres_init (V0, modulus);
  mpres_set_ui (V0, 2UL, modulus); /* V0 = V_0(S) = 2 */

  if (uk == 1)
    {
      mpres_set (R, S, modulus);
      while (po2-- > 0)
        {
          mpres_sqr (R, R, modulus);
          mpres_sub (R, R, V0, modulus);
        }
      mpres_clear (V0, modulus);
      goto exit;
    }

  for (j = 1; j <= uk / 2; j <<= 1);

  mpres_init (Vi, modulus);
  mpres_init (Vi1, modulus);

  /* i = 1. Vi = V_i(S), Vi1 = V_{i+1}(S) */
  mpres_set (Vi, S, modulus);
  mpres_sqr (Vi1, S, modulus);
  mpres_sub (Vi1, Vi1, V0, modulus);
  j >>= 1;

  while (j > 1)
    {
      if ((uk & j) != 0)
	{
	  /* i' = 2i + 1.
	     V_{i'} = V_{2i + 1} = V_{i+1 + i} = V_{i+1} * V_{i} - V_1
	     V_{i'+1} = V_{2i + 2} = {V_{i+1}}^2 - V_0. */
	  mpres_mul (Vi, Vi, Vi1, modulus);
	  mpres_sub (Vi, Vi, S, modulus);
	  mpres_sqr (Vi1, Vi1, modulus);
	  mpres_sub (Vi1, Vi1, V0, modulus);
	}
      else
	{
	  /* i' = 2i. 
	     V_{i'} = V_{2i} = {V_i}^2 - V0.
	     V_{i'+1} = V_{2i + 1} = V_{i+1 + i} = V_{i+1} * V_{i} - V_1 */
	  mpres_mul (Vi1, Vi, Vi1, modulus);
	  mpres_sub (Vi1, Vi1, S, modulus);

	  mpres_sqr (Vi, Vi, modulus);
	  mpres_sub (Vi, Vi, V0, modulus);
	}
      j >>= 1;
    }

  /* Least significant bit of uk is always 1 */
  mpres_mul (Vi, Vi, Vi1, modulus);
  mpres_sub (Vi, Vi, S, modulus);

  while (po2-- > 0)
    {
      mpres_sqr (Vi, Vi, modulus);
      mpres_sub (Vi, Vi, V0, modulus);
    }

  mpres_set (R, Vi, modulus);

  mpres_clear (Vi, modulus);
  mpres_clear (Vi1, modulus);
  mpres_clear (V0, modulus);

exit:
  if (test_verbose(OUTPUT_TRACE))
    {
      mpz_t tz;
      mpz_init (tz);
      mpres_get_z (tz, R, modulus);
      gmp_printf ("%Zd /* PARI C V */\n", tz);
      mpz_clear (tz);
    }
}


/* Set R1[i] = V_{i+k}(Q) * F1[i] or U_{i+k}(Q) * F[i], for 0 <= i < len
   (and R2[i] = V_{i+k}(Q) * F2[i] etc, if both R2, F2 are non-NULL)
   We compute V_{i+k+1}(Q) by V_{i+k}(Q)*V_1(Q) - V_{i+k-1}(Q).
   For U, we compute U_{i+k+1}(Q) by U_{i+k}(Q)*V_1(Q) - U_{i+k-1}(Q).
   The values of V_1(Q), V_{k-1}(Q) and V_k(Q) and V_k(Q) are in 
   V1, Vk_1 and Vk, resp. 
   The values of Vk_1 and Vk are clobbered. 
   */
static void
scale_by_chebyshev (listz_t R1, const listz_t F1, 
                    listz_t R2, const listz_t F2, 
                    const uint64_t len,
                    mpmod_t modulus, const mpres_t V1, mpres_t Vk_1, 
                    mpres_t Vk)
{
  mpres_t Vt;
  uint64_t i;

  mpres_init (Vt, modulus);

  for (i = 0; i < len; i++)
    {
      mpres_mul_z_to_z (R1[i], Vk, F1[i], modulus);
      if (R2 != NULL && F2 != NULL)
        mpres_mul_z_to_z (R2[i], Vk, F2[i], modulus);
      mpres_mul (Vt, Vk, V1, modulus);
      mpres_sub (Vt, Vt, Vk_1, modulus);
      mpres_set (Vk_1, Vk, modulus); /* Could be a swap */
      mpres_set (Vk, Vt, modulus); /* Could be a swap */
    }

  mpres_clear (Vt, modulus);
}


/* For a given reciprocal polynomial 
   F(x) = f_0 + sum_{i=1}^{deg} f_i V_i(x+1/x),
   compute F(\gamma x)F(\gamma^{-1} x), with Q = \gamma + 1 / \gamma

   list_scale_V_ntt() needs no temp space.
   
   For list_scale_V(), if NTT is used, needs 4 * deg + 3 entries in tmp.
   If no NTT is used, needs 4 * deg + 2 + (memory use of list_sqr_reciprocal)
*/

const int trace_callbacks = 0;

typedef struct {
  FILE *f;
  void *buf;
  size_t bufsize;
} state_file_t;

static void
file_reader (void * const p, mpz_t r)
{
  state_file_t *state = p;
  
  read_residue (state->f, r, state->buf, state->bufsize);
}


typedef struct {
  mpz_t *mpzv_read, *mpzv_write;
  FILE *file_read, *file_write;
  file_word_t *buf;
  size_t index, bufsize;
  mpres_t V1, Vi, Vi_1, tmp;
  mpmod_t modulus;
  mpz_t mpz;
} stateV_t;

static void 
readerV (void * const p, mpz_t r)
{
  stateV_t * const state = p;

  mpres_mul_z_to_z (r, state->Vi, state->mpzv_read[state->index], 
    state->modulus);
  
  if (trace_callbacks)
    outputf(OUTPUT_TRACE, "Chebyshev_V(%d, Q)*f_%d = %Zd /* readerV */\n",
      (int)state->index, (int)state->index, r);
}

static void 
readerV_file (void * const p, mpz_t r)
{
  stateV_t * const state = p;
  
  seek_read_residue (state->file_read, r, state->buf, state->bufsize, 
    state->index);
  mpres_mul_z_to_z (r, state->Vi, r, state->modulus);

  if (trace_callbacks)
    outputf(OUTPUT_TRACE, "Chebyshev_V(%d, Q)*f_%d = %Zd /* readerV */\n",
      (int)state->index, (int)state->index, r);
}

static void 
writerV (void * const p, const mpz_t r)
{
  stateV_t * const state = p;

  mpres_mul_z_to_z (state->mpzv_write[state->index], state->Vi, r, 
    state->modulus);

  if (trace_callbacks)
    outputf(OUTPUT_TRACE, "r_%d = %Zd, g_%d = %Zd /* writerV */\n",
      (int)state->index, r, (int)state->index, state->mpzv_write[state->index]);
  
  state->index++;
  mpres_mul (state->tmp, state->Vi, state->V1, state->modulus);
  mpres_sub (state->tmp, state->tmp, state->Vi_1, state->modulus);
  mpres_set (state->Vi_1, state->Vi, state->modulus);
  mpres_set (state->Vi, state->tmp, state->modulus);
}

static void 
writerV_file (void * const p, const mpz_t r)
{
  stateV_t * const state = p;
  
  mpres_mul_z_to_z (state->mpz, state->Vi, r, state->modulus);
  seek_write_residue (state->file_write, state->mpz, state->buf, 
    state->bufsize, state->index);

  if (trace_callbacks)
    outputf(OUTPUT_TRACE, "r_%d = %Zd, g_%d = %Zd /* writerV */\n",
      (int)state->index, r, (int)state->index, state->mpz);
  
  state->index++;
  mpres_mul (state->tmp, state->Vi, state->V1, state->modulus);
  mpres_sub (state->tmp, state->tmp, state->Vi_1, state->modulus);
  mpres_set (state->Vi_1, state->Vi, state->modulus);
  mpres_set (state->Vi, state->tmp, state->modulus);
}


typedef struct {
  mpz_t *mpzv;
  FILE *file;
  file_word_t *buf;
  size_t index, bufsize;
  mpz_t mpz, modulus;
} stateD_t;

static void 
writer_diff (void * const p, const mpz_t r)
{
  stateD_t * const state = p;

  mpz_sub (state->mpz, r, state->mpzv[state->index]);
  
  if (mpz_odd_p (state->mpz))
    mpz_add (state->mpz, state->mpz, state->modulus);
  ASSERT (mpz_even_p (state->mpz));
  mpz_tdiv_q_2exp (state->mpz, state->mpz, 1);
  mpz_mod (state->mpzv[state->index], state->mpz, state->modulus);

  if (trace_callbacks)
    outputf(OUTPUT_TRACE, "r_%d = %Zd /* writer_diff */\n",
      (int)state->index, state->mpzv[state->index]);

  state->index++;
}

static void 
writer_diff_file (void * const p, const mpz_t r)
{
  stateD_t * const state = p;

  seek_read_residue (state->file, state->mpz, state->buf, state->bufsize, 
    state->index);
  mpz_sub (state->mpz, r, state->mpz);

  if (mpz_odd_p (state->mpz))
    mpz_add (state->mpz, state->mpz, state->modulus);
  ASSERT (mpz_even_p (state->mpz));
  mpz_tdiv_q_2exp (state->mpz, state->mpz, 1);
  mpz_mod (state->mpz, state->mpz, state->modulus);
  seek_write_residue (state->file, state->mpz, state->buf, state->bufsize, 
    state->index);

  if (trace_callbacks)
    outputf(OUTPUT_TRACE, "r_%d = %Zd /* writer_diff */\n",
      (int)state->index, state->mpz);

  state->index++;
}


static void
list_scale_V_ntt (listz_handle_t R, const listz_handle_t F, 
              const mpres_t Q, const uint64_t deg, mpmod_t modulus, 
	      mpzspv_handle_t ntt_handle)
{
  if (deg == 0)
    {
      ASSERT_ALWAYS (F->storage == 0 && R->storage == 0);
      mpz_t tmp;
      mpz_init (tmp);
      mpz_mul (tmp, F->data.mem[0], F->data.mem[0]);
      mpz_mod (R->data.mem[0], tmp, modulus->orig_modulus);
      mpz_clear (tmp);
      return;
    }
  
  list_output_poly_file (F, deg + 1, 0, 1, "list_scale_V_ntt: F = ", "\n", 
    OUTPUT_TRACE);
  /* Convert F[0, ..., deg] to NTT */
  if (F->storage == 0)
    mpzspv_fromto_mpzv (ntt_handle, (spv_size_t) 0, deg + 1, NULL, F->data.mem, 
        NULL, NULL);
  else
    {
      state_file_t state;
      state.f = F->data.file;
      rewind (state.f);
      state.buf = (file_word_t *) mpz_export (NULL, &state.bufsize, -1, 
          sizeof(file_word_t), -1, 0, modulus->orig_modulus);
      mpzspv_fromto_mpzv (ntt_handle, (spv_size_t) 0, deg + 1, file_reader, 
          &state, NULL, NULL);
      free (state.buf);
    }
  
  if (test_verbose(OUTPUT_TRACE))
    mpzspv_print (ntt_handle, 0, deg + 1, "list_scale_V_ntt: Before squaring ");
  
  /* Compute F^2 in NTT */
  mpzspv_sqr_reciprocal (ntt_handle, deg + 1);
  
  if (test_verbose(OUTPUT_TRACE))
    mpzspv_print (ntt_handle, 0, deg + 1, "list_scale_V_ntt: After squaring ");

#if defined(_OPENMP)
#pragma omp parallel if (deg > 1000)
#endif
  {
    /* Convert F^2 from NTT, add weights and store in F[0, ..., 2*deg], 
       at the same time add weights to F[0, ..., deg] and store that in NTT */
    stateV_t state;
    uint64_t l, start_i;

    mpmod_init_set (state.modulus, modulus);
    mpres_init (state.V1, state.modulus);
    mpres_set (state.V1, Q, state.modulus);
    mpres_init (state.Vi_1, state.modulus);
    mpres_init (state.Vi, state.modulus);
    mpres_init (state.tmp, state.modulus);
    mpz_init (state.mpz);
    state.buf = (file_word_t *) mpz_export (NULL, &state.bufsize, -1, 
        sizeof(file_word_t), -1, 0, state.modulus->orig_modulus);
    
    /* Read and write i = 0, ..., deg */
    get_chunk (&start_i, &l, deg + 1);
    state.index = start_i;
    V (state.Vi_1, state.V1, (int64_t) start_i - 1, state.modulus);
    V (state.Vi, state.V1, start_i, state.modulus);

    if (F->storage == 0)
      {
        state.mpzv_read = F->data.mem;
        state.mpzv_write = R->data.mem;
        state.file_read = NULL;
        state.file_write = NULL;
        mpzspv_fromto_mpzv (ntt_handle, (spv_size_t) 0, l, 
            &readerV, &state, &writerV, &state);
      }
    else
      {
        state.mpzv_read = NULL;
        state.mpzv_write = NULL;
        /* FIXME clone file handles, seek to correct position */
        state.file_read = F->data.file;
        state.file_write = R->data.file;
        rewind (state.file_read);
        rewind (state.file_write);
        mpzspv_fromto_mpzv (ntt_handle, (spv_size_t) 0, l, 
            &readerV_file, &state, &writerV_file, &state);
      }
    
    /* Write the remaining i = deg+1, ..., 2*deg+1 */
    get_chunk (&start_i, &l, deg);
    start_i += deg + 1;
    state.index = start_i;
    V (state.Vi_1, state.V1, (int64_t) start_i - 1, state.modulus);
    V (state.Vi, state.V1, start_i, state.modulus);
    if (F->storage == 0)
      {
        mpzspv_fromto_mpzv (ntt_handle, (spv_size_t) deg + 1, l, 
            NULL, NULL, &writerV, &state);
      } else {
        mpzspv_fromto_mpzv (ntt_handle, (spv_size_t) deg + 1, l, 
            NULL, NULL, &writerV_file, &state);
      }
    
    mpres_clear (state.V1, state.modulus);
    mpres_clear (state.Vi_1, state.modulus);
    mpres_clear (state.Vi, state.modulus);
    mpres_clear (state.tmp, state.modulus);
    mpmod_clear (state.modulus);
    mpz_clear (state.mpz);
    free (state.buf);
  }

  list_output_poly_file (R, 2*deg + 1, 0, 1, "Gw(x) = ", 
      "; /* PARI list_scale_V_ntt */\n", OUTPUT_TRACE);

  /* Square the weighted F in NTT */
  mpzspv_sqr_reciprocal (ntt_handle, deg + 1);
  
  /* Convert from NTT and take half the difference from R */
  {
    stateD_t state;
    
    state.index = 0;
    mpz_init_set (state.modulus, modulus->orig_modulus);
    mpz_init (state.mpz);
    state.buf = (file_word_t *) mpz_export (NULL, &state.bufsize, -1, 
      sizeof(file_word_t), -1, 0, state.modulus);
    if (R->storage == 0)
      {
        state.mpzv = R->data.mem;
        state.file = NULL;
        mpzspv_fromto_mpzv (ntt_handle, (spv_size_t) 0, 2*deg + 1, 
            NULL, NULL, &writer_diff, &state);
      } else {
        state.mpzv = NULL;
        state.file = R->data.file;
        rewind (state.file);
        mpzspv_fromto_mpzv (ntt_handle, (spv_size_t) 0, 2*deg + 1, 
            NULL, NULL, &writer_diff_file, &state);
      }
    mpz_clear (state.modulus);
    mpz_clear (state.mpz);
    free (state.buf);
  }
  list_output_poly_file (R, deg + 1, 0, 1, "list_scale_V_ntt: R = ", "\n", 
    OUTPUT_TRACE);
}


static void
list_scale_V (listz_t R, const listz_t F, const mpres_t Q, 
              const uint64_t deg, mpmod_t modulus, listz_t tmp, 
              const uint64_t tmplen, 
	      mpzspv_handle_t ntt_handle)
{
  uint64_t i;
  const listz_t G = tmp, H = tmp + 2 * deg + 1, newtmp = tmp + 4 * deg + 2;
  const uint64_t newtmplen = tmplen - 4 * deg - 2;
#ifdef WANT_ASSERT
  mpz_t leading;
#endif
  
  if (deg == 0)
    {
      ASSERT(tmplen >= 1);
      mpz_mul (tmp[0], F[0], F[0]);
      mpz_mod (R[0], tmp[0], modulus->orig_modulus);
      return;
    }
  
  outputf (OUTPUT_TRACE, "\nN=%Zd; deg = %lu; /* PARI list_scale_V */\n", 
           modulus->orig_modulus, deg);
  if (test_verbose(OUTPUT_TRACE))
    {
      mpres_t out_t;
      mpres_init (out_t, modulus);
      mpres_get_z (out_t, Q, modulus);
      outputf (OUTPUT_TRACE, "Q = Mod(%Zd,N); /* PARI list_scale_V */\n", out_t);
      mpres_clear (out_t, modulus);
    }
  list_output_poly (F, deg + 1, 0, 1, "F(x) = ", "; /* PARI list_scale_V */\n", 
		    OUTPUT_TRACE);

  /* Make sure newtmplen does not underflow */
  ASSERT_ALWAYS (tmplen >= 4 * deg + 2);
#ifdef WANT_ASSERT
  mpz_init (leading);
  mpz_mul (leading, F[deg], F[deg]);
  mpz_mod (leading, leading, modulus->orig_modulus);
#endif

  for (i = 0; i <= deg; i++)
    {
      ASSERT_ALWAYS (mpz_sgn (F[i]) >= 0 && mpz_cmp (F[i], modulus->orig_modulus) < 0);
    }

  if (ntt_handle != NULL)
    ntt_sqr_reciprocal (G, F, deg + 1, ntt_handle);
  else
    {
      list_sqr_reciprocal (G, F, deg + 1, modulus->orig_modulus, 
                           newtmp, newtmplen);
      list_mod (G, G, 2*deg + 1, modulus->orig_modulus);
    }

  list_output_poly (G, 2 * deg + 1, 0, 1, "G(x) = ", 
		    " /* PARI list_scale_V */\n", OUTPUT_TRACE);
  outputf (OUTPUT_TRACE, "if (G(x) != F(x)^2, print(\"G(x) != F(x)^2 in "
           "list_scale_V()\")); /* PARI list_scale_V */\n");

  /* Compute G[i] = V_i(Q) * G[i] for i = 0, ..., 2*deg
     and H[i] = V_i(Q) * F[i] for i = 0, ..., deg. */

#if defined(_OPENMP)
#pragma omp parallel if (deg > 1000)
#endif
  {
    mpmod_t modulus_local;
    uint64_t l, start_i;
    mpres_t Vi, Vi_1;
    
    mpmod_init_set (modulus_local, modulus);
    mpres_init (Vi_1, modulus_local);
    mpres_init (Vi, modulus_local);
    
    /* Do G[i] and H[i], i = 0, ..., deg */
    get_chunk (&start_i, &l, deg + 1);
    V (Vi_1, Q, (int64_t) start_i - 1, modulus_local);
    V (Vi, Q, start_i, modulus_local);
    scale_by_chebyshev (G + start_i, G + start_i, H + start_i, F + start_i, 
                        l, modulus_local, Q, Vi_1, Vi);
    
    /* Do the remaining entries G[i], i = deg+1, ..., 2*deg+1 */
    get_chunk (&start_i, &l, deg);
    start_i += deg + 1;
    V (Vi, Q, start_i, modulus_local);
    V (Vi_1, Q, (int64_t) start_i - 1, modulus_local);
    scale_by_chebyshev (G + start_i, G + start_i, NULL, NULL, l, modulus_local, 
                        Q, Vi_1, Vi);
    
    mpres_clear (Vi_1, modulus_local);
    mpres_clear (Vi, modulus_local);
    mpmod_clear (modulus_local);
  }

  list_output_poly (G, 2*deg + 1, 0, 1, "Gw(x) = ", 
                    "; /* PARI list_scale_V */\n", OUTPUT_TRACE);

  for (i = 0; i <= deg; i++)
    {
      ASSERT_ALWAYS (mpz_sgn (H[i]) >= 0 && mpz_cmp (H[i], modulus->orig_modulus) < 0);
    }

  list_output_poly (H, deg + 1, 0, 1, "H(x) = ", "; /* PARI list_scale_V */\n", 
		    OUTPUT_TRACE);

  if (ntt_handle != NULL)
    ntt_sqr_reciprocal (H, H, deg + 1, ntt_handle);
  else
    {
      list_sqr_reciprocal (H, H, deg + 1, modulus->orig_modulus, 
		           newtmp, newtmplen);
      list_mod (H, H, 2 * deg + 1, modulus->orig_modulus);
    }

  list_output_poly (H, 2*deg + 1, 0, 1, "H(x)^2 == ", 
                    " /* PARI list_scale_V */\n", OUTPUT_TRACE);

  for (i = 0; i <= 2 * deg; i++)
    {
      mpz_sub (H[i], H[i], G[i]);
      if (mpz_odd_p(H[i]))
        mpz_add (H[i], H[i], modulus->orig_modulus);
      mpz_tdiv_q_2exp(H[i], H[i], 1);
      if (mpz_sgn (H[i]) < 0)
        mpz_add (H[i], H[i], modulus->orig_modulus);
      mpz_set (R[i], H[i]);
      ASSERT_ALWAYS (mpz_sgn (R[i]) >= 0 && mpz_cmp(R[i], modulus->orig_modulus) <= 0);
    }

  list_output_poly (R, 2 * deg, 0, 1, "R(x) = ", " /* PARI list_scale_V */\n", OUTPUT_TRACE);

#ifdef WANT_ASSERT
  mpz_mod (R[2 * deg], R[2 * deg], modulus->orig_modulus);
  ASSERT (mpz_cmp (leading, R[2 * deg]) == 0);
  mpz_clear (leading);
#endif
}


/* Evaluate a polynomial of degree n-1 with all coefficients given in F[],
   or of degree n with an implicit leading 1 monomial not stored in F[],
   at x (mod modulus). Result goes in r. tmp needs 2 entries. */

ATTRIBUTE_UNUSED static void 
list_eval_poly (mpz_t r, const listz_t F, const mpz_t x, 
		const uint64_t n, const int monic, const mpz_t modulus, 
		listz_t tmp)
{
  uint64_t i;

  mpz_set_ui (tmp[0], 1UL);
  mpz_set_ui (r, 0UL);

  for (i = 0UL; i < n; i++)
    {
      /* tmp[0] = x^i */
      mpz_mul (tmp[1], F[i], tmp[0]);
      mpz_mod (tmp[1], tmp[1], modulus);
      mpz_add (r, r, tmp[1]);

      mpz_mul (tmp[1], tmp[0], x);
      mpz_mod (tmp[0], tmp[1], modulus);
    }

  if (monic)
    mpz_add (r, r, tmp[0]);

  mpz_mod (r, r, modulus);
}


/* Build a polynomial with roots r^2i, i in the sumset of the sets in "sets".
   The parameter Q = r + 1/r. This code uses the fact that the polynomials 
   are symmetric. Requires that the first set in "sets" has cardinality 2,
   all sets must be symmetric around 0. The resulting polynomial of degree 
   2*d is F(x) = f_0 + \sum_{1 <= i <= d} f_i (x^i + 1/x^i). The coefficient
   f_i is stored in F[i], which therefore needs d+1 elements. */

static uint64_t
poly_from_sets_V (listz_handle_t F_param, const mpres_t Q, set_list_t *sets, 
		  listz_t tmp, const uint64_t tmplen, mpmod_t modulus,
		  mpzspv_handle_t ntt_handle)
{
  uint32_t nr;
  uint64_t deg;
  mpres_t Qt;
  listz_t F = (F_param->storage == 0) ? F_param->data.mem : tmp;
  
  ASSERT_ALWAYS (sets->num_sets > 0UL);
  /* Check that the cardinality of first set is 2 */
  ASSERT_ALWAYS (sets->sets[0].card == 2UL);
  /* Check that the first set is symmetric around 0 */
  ASSERT_ALWAYS (sets->sets[0].elem[0] == -sets->sets[0].elem[1]);
  /* We allow disk-stored F only with NTT */
  ASSERT_ALWAYS (F_param->storage == 0 || ntt_handle != NULL);

  if (test_verbose (OUTPUT_TRACE))
    {
      mpz_t t;
      mpz_init (t);
      mpres_get_z (t, Q, modulus);
      outputf (OUTPUT_TRACE, "poly_from_sets_V (F, Q = %Zd, sets)\n", t);
      mpz_clear (t);
    }

  mpres_init (Qt, modulus);
  
  outputf (OUTPUT_DEVVERBOSE, " (processing set of size 2");

  V (Qt, Q, sets->sets[0].elem[0], modulus); /* First set in sets is {-k, k}. Qt = V_k(Q) */
  
  ASSERT_ALWAYS (F != tmp || tmplen >= 2);
  mpres_neg (Qt, Qt, modulus);
  mpres_get_z (F[0], Qt, modulus);
  mpz_set_ui (F[1], 1UL);
  deg = 1;
  /* Here, F(x) = (x - r^{2k_1})(x - r^{-2k_1}) / x = 
                  (x^2 - x (r^{2k_1} + r^{-2k_1}) + 1) / x =
		  (x + 1/x) - V_{2k_1}(r + 1/r) */

  for (nr = sets->num_sets - 1; nr > 0; nr--)
    {
      const set_t *curr_set = sets->sets + nr;
      const uint32_t c = curr_set->card;
      unsigned long i;

      /* Assuming the sets are sorted in order of ascending cardinality, 
         we process them back-to-front so the sets of cardinality 2 are 
         processed last, but skipping the first set which we processed 
         already. */
      
      /* Sets of cardinality 2 are processed below */
      if (c == 2)
        break;
      
      /* Process this set of odd cardinality */
      outputf (OUTPUT_DEVVERBOSE, " %lu", c);

      ASSERT_ALWAYS (c % 2UL == 1UL);
      ASSERT_ALWAYS (curr_set->elem[(c - 1UL) / 2UL] == 0);
      /* Generate the F(Q^{2k_i} * X)*F(Q^{-2k_i} * X) polynomials.
         Each is symmetric of degree 2*deg, so each has deg+1 coeffients
         in standard basis. */
      for (i = 0; i < (c - 1) / 2; i++)
        {
          const uint64_t prod_len = 2 * deg + 1;
          const uint64_t prod_offset = deg + 1 + i * prod_len;
          const uint64_t tmpadd = (F == tmp) ? prod_offset + prod_len : 0;
          
          /* Check it's symmetric */
          ASSERT_ALWAYS (curr_set->elem[i] == -curr_set->elem[c - 1L - i]);
          V (Qt, Q, curr_set->elem[i], modulus);
          ASSERT_ALWAYS (mpz_cmp_ui (F[deg], 1) == 0); /* Check it's monic */
          ASSERT_ALWAYS (tmplen >= tmpadd);
          /* Product has degree 2*deg, so has 2*deg+1 coefficients in 
             standard basis, and occupies 
             F[(2 * i + 1) * lenF ... (2 * i + 1) * lenF  + 2 * deg] */
          if (ntt_handle != NULL)
            {
              _listz_handle_t F_handle = {0, deg + 1, F_param->words, 
                  .data.mem = F};
              _listz_handle_t R_handle = {0, 2*deg + 1, F_param->words, 
                  .data.mem = F + prod_offset};
              list_scale_V_ntt (&R_handle, &F_handle, Qt, 
                            deg, modulus, ntt_handle);
            }
          else
            list_scale_V (F + prod_offset, F, Qt, deg, 
                          modulus, tmp + tmpadd, tmplen - tmpadd, NULL);
          ASSERT_ALWAYS (mpz_cmp_ui (F[prod_offset + prod_len - 1], 1) 
                         == 0); /* Check it's monic */
        }
      /* Multiply the polynomials */
      for (i = 0; i < (c - 1) / 2; i++)
        {
          /* So far, we have the product 
             F(X) * F(Q^{2k_j} * X) * F(Q^{-2k_j} * X), 1 <= j <= i,
             at F. This product has degree 2 * deg + i * 4 * deg, that is
             (2 * i + 1) * 2 * deg, which means lenF = (2 * i + 1) * deg + 1
             coefficients in F[0 ... (i * 2 + 1) * deg]. 
             We now multiply by the output of scale_V() which occupies
             prod_len coefficients, starting at prod_offset. The resulting
             RLP will have newlenF coefficients. */
          const uint64_t lenF = (2 * i + 1) * deg + 1; 
          const uint64_t prod_len = 2 * deg + 1; 
          const uint64_t prod_offset = deg + 1 + i * prod_len; 
          const uint64_t newlenF = lenF + prod_len - 1;
          const uint64_t tmpadd = 
            (F == tmp) ? deg + 1 + (c - 1) / 2 * prod_len : 0;

          ASSERT_ALWAYS (mpz_cmp_ui (F[lenF - 1], 1) == 0);
          ASSERT_ALWAYS (mpz_cmp_ui (F[prod_offset + prod_len - 1], 1) == 0);
          ASSERT_ALWAYS (tmplen >= tmpadd);

          list_output_poly (F, lenF, 0, 1, 
                            "f(x) = ", "; /* PARI poly_from_sets_V */\n",
                            OUTPUT_TRACE);
          list_output_poly (F + prod_offset, prod_len, 0, 1, 
                            "g(x) = ", "; /* PARI poly_from_sets_V */\n", 
                            OUTPUT_TRACE);
          list_mul_reciprocal (F, F, lenF, F + prod_offset, prod_len, 
                               modulus->orig_modulus, 
                               tmp + tmpadd, tmplen - tmpadd);
          list_mod (F, F, newlenF, modulus->orig_modulus);
          list_output_poly (F, newlenF, 0, 1, "f(x) * g(x) == ", 
                            " /* PARI poly_from_sets_V */\n", OUTPUT_TRACE);
          ASSERT_ALWAYS (mpz_cmp_ui (F[newlenF - 1], 1) == 0);
        }
      deg *= c;
    }

  /* If we built F in tmp so far (with -treefile) we have to write it to 
     a disk file */
  if (F_param->storage != 0)
    {
      rewind (F_param->data.file);
      write_residues (F_param->data.file, F, deg + 1, F_param->words);
      F = NULL;
    }

  /* Do the sets of cardinality 2 (except first one which we already used) */
  for ( ; nr > 0; nr--)
    {
      const set_t *curr_set = sets->sets + nr;
      const uint32_t c = curr_set->card;

      ASSERT_ALWAYS (c == 2);
      outputf (OUTPUT_DEVVERBOSE, " %lu", c);
      list_output_poly_file (F_param, deg + 1, 0, 1, "poly_from_sets_V: F = ", "\n", 
        OUTPUT_TRACE);
      /* Check it's symmetric */
      ASSERT_ALWAYS (curr_set->elem[0] == -curr_set->elem[1]);
      V (Qt, Q, curr_set->elem[0], modulus);
      if (ntt_handle != NULL)
        {
          list_scale_V_ntt (F_param, F_param, Qt, deg, modulus, ntt_handle);
        }
      else
        list_scale_V (F, F, Qt, deg, modulus, tmp, tmplen, NULL);
      deg *= 2;

      /* Check it's monic */
      listz_handle_get2 (F_param, tmp[0], deg);
      ASSERT_ALWAYS (mpz_cmp_ui (tmp[0], 1UL) == 0);
    }

  mpres_clear (Qt, modulus);
  outputf (OUTPUT_DEVVERBOSE, ")");

  return deg;
}


static listz_handle_t 
build_F_ntt (const mpres_t P_1, set_list_t *S_1, 
	     const faststage2_param_t *params, 
	     mpmod_t modulus)
{
  mpzspm_t F_ntt_context = NULL;
  mpzspv_handle_t ntt_handle = NULL;
  const spv_size_t nttlen = (spv_size_t)1 << ceil_log2 (params->s_1 / 2 + 1);
  uint64_t tmplen;
  listz_t tmp = NULL;
  listz_handle_t F = NULL;
  char *filename_f = NULL, *filename_ntt = NULL;
  long timestart, realstart;
  unsigned long i;

  timestart = cputime ();
  realstart = realtime ();
  
  /* Precompute the small primes, primitive roots and inverses etc. for 
     the NTT. The code to multiply wants a 3*k-th root of unity, where 
     k is the smallest power of 2 with k > s_1/2 */
  
  /* Allocate all the memory we'll need for building f */
  if (params->file_stem != NULL)
    {
      filename_f = malloc ((strlen(params->file_stem) + 6) * sizeof (char));
      filename_ntt = malloc ((strlen(params->file_stem) + 3) * sizeof (char));
      if (filename_f == NULL || filename_ntt == NULL)
        goto clear_and_exit;
      sprintf (filename_f, "%s.fmpz", params->file_stem);
      sprintf (filename_ntt, "%s.f", params->file_stem);
    }
  
  F = listz_handle_init2 (filename_f, params->s_1 / 2 + 1 + 1, 
    modulus->orig_modulus); /* Another +1 in len because 
    poly_from_sets_V stores the leading 1 monomial for each factor */
  free (filename_f);

  F_ntt_context = mpzspm_init (3UL << ceil_log2 (params->s_1 / 2 + 1), 
			       modulus->orig_modulus);
  
  tmplen = params->s_1;
  ASSERT_ALWAYS(tmplen > 0);
  /* All but one factors of 2 are handled by list_scale_V_ntt() 
     which needs no temp space */
  while (tmplen % 4 == 0)
    tmplen /= 2;
  tmplen = 7*tmplen + list_mul_mem (tmplen);
  tmp = init_list2 (tmplen, (unsigned int) abs (modulus->bits));
  ntt_handle = mpzspv_init_handle (filename_ntt, nttlen, F_ntt_context);
  free (filename_ntt);

  if (F == NULL || F_ntt_context == NULL || tmp == NULL || 
      ntt_handle == NULL)
    {
      outputf (OUTPUT_ERROR, "build_F_ntt(): Could not allocate memory\n");
      if (F != NULL)
        {
          listz_handle_clear (F);
          F = NULL;
        }
      goto clear_and_exit;
    }
  
  mpzspm_print_CRT_primes (OUTPUT_DEVVERBOSE, "CRT modulus for building F = ",
		    F_ntt_context);
  
  outputf (OUTPUT_VERBOSE, "Computing F from factored S_1");
  
  i = poly_from_sets_V (F, P_1, S_1, tmp, tmplen, modulus, ntt_handle);
  ASSERT_ALWAYS(2 * i == params->s_1);
  ASSERT_ALWAYS(F->storage != 0 || mpz_cmp_ui (F->data.mem[i], 1UL) == 0);
  
  print_elapsed_time (OUTPUT_VERBOSE, timestart, realstart);
  list_output_poly_file (F, params->s_1 / 2 + 1, 0, 1, "F(x) = ", 
    " /* PARI build_F_ntt */\n", OUTPUT_TRACE);
  
clear_and_exit:
  if (tmp != NULL)
    clear_list (tmp, tmplen);
  tmp = NULL;
  mpzspv_clear_handle (ntt_handle);
  ntt_handle = NULL;
  
  mpzspm_clear (F_ntt_context);
  F_ntt_context = NULL;

  return F;
}

typedef struct {
  mpmod_t modulus;
  mpres_t r[2];
  mpz_t t;
} pm1_g_state_t;

static void
pm1_sequence_g_prod (void *state_p, mpz_t r)
{
  pm1_g_state_t *state = state_p;
  
  mpz_set (r, state->t);
  mpres_mul_z_to_z (state->t, state->r[1], state->t, state->modulus);
  mpres_mul (state->r[1], state->r[1], state->r[0], state->modulus);
}

/* Compute g_i = x_0^{M-i} * r^{(M-i)^2} for 0 <= i < l. 
   x_0 = b_1^{k_2 + m_1 * P}. r = b_1^P. 
   Stores the result in g[0 ... l] and/or in g_ntt[offset ... offset + l] */

static void
pm1_sequence_g (listz_t g_mpz, mpzspv_handle_t g_handle, const mpres_t b_1, 
                const uint64_t P, const uint64_t M_param, 
		const uint64_t l_param, const mpz_t m_1, 
		const int64_t k_2, mpmod_t modulus_param)
{
  long timestart, realstart;

  outputf (OUTPUT_VERBOSE, "Computing g_i");
  outputf (OUTPUT_DEVVERBOSE, "\npm1_sequence_g: P = %" PRIu64
            ", M_param = %" PRIu64 ", l_param = %" PRIu64 
            ", k_2 = %" PRId64 , P, M_param, l_param, k_2);
  outputf (OUTPUT_DEVVERBOSE, ", m_1 = %Zd\n", m_1);
  timestart = cputime ();
  realstart = realtime ();

#ifdef _OPENMP
#pragma omp parallel if (l_param > 100) private(i)
#endif
  {
    uint64_t M, l, offset;
    pm1_g_state_t state;
    mpz_t t, t1;
    mpz_t tM;
    mpres_t r, x_0, x_Mi, r2;
    int want_output = 1;
    uint64_t i;

    /* When multi-threading, we adjust the parameters for each thread */
    get_chunk (&offset, &l, l_param);
    M = M_param - offset;

#ifdef _OPENMP
    outputf (OUTPUT_DEVVERBOSE, 
             "pm1_sequence_g: thread %d has l = %" PRIu64 
             ", offset = %" PRIu64 ", M = %" PRIu64 ".\n", 
             omp_get_thread_num(), l, offset, M);
    
    /* Let only the master thread print stuff */
    want_output = (omp_get_thread_num() == 0);

    if (want_output)
      outputf (OUTPUT_VERBOSE, " using %d threads", omp_get_num_threads());
#endif

    /* Make a private copy of the mpmod_t struct */
    mpmod_init_set (state.modulus, modulus_param);

    mpz_init (t);
    mpz_init (t1);
    mpz_init (tM);
    mpres_init (r, state.modulus);
    mpres_init (r2, state.modulus);
    mpres_init (x_0, state.modulus);
    mpres_init (state.r[0], state.modulus);
    mpres_init (state.r[1], state.modulus);
    mpz_init (state.t);

    if (want_output)
      {
        if (test_verbose (OUTPUT_TRACE))
          { 
            mpres_get_z (t, b_1, state.modulus);
            outputf (OUTPUT_TRACE, "\n/* pm1_sequence_g */ N = %Zd; "
                     "b_1 = Mod(%Zd, N); /* PARI */\n", state.modulus->orig_modulus, t);
            outputf (OUTPUT_TRACE, "/* pm1_sequence_g */ P = %" PRIu64
                     "; M = %"PRIu64"; k_2 = %" PRId64, P, M, k_2);
            outputf (OUTPUT_TRACE, "; m_1 = %Zd; /* PARI */\n", m_1);
            outputf (OUTPUT_TRACE,"/* pm1_sequence_g */ r = b_1^(P/2); /* PARI */\n");
            outputf (OUTPUT_TRACE, "/* pm1_sequence_g */ x_0 = "
                     "b_1^(k_2 + m_1*P); /* PARI */\n");
          }
      }

    /* We use (M-(i+1))^2 = (M-i)^2 + 2*(-M+i) + 1 */
    mpz_set_uint64 (t, P/2);
    mpres_pow (r, b_1, t, state.modulus);     /* r[0] = b_1^(P/2) = r */
    if (test_verbose (OUTPUT_TRACE))
      {
        mpres_get_z (t, r, state.modulus);
        outputf (OUTPUT_TRACE, "/* pm1_sequence_g */ r == %Zd /* PARI C */\n", t);
      }
    
    /* FIXME: This is a huge mess, clean up some time */

    mpz_set_uint64 (tM, M);
    mpres_pow (r2, r, tM, state.modulus);
    mpres_pow (r2, r2, tM, state.modulus);    /* r2 = r^{(M-i)^2}, i = 0 */
    mpz_neg (t, tM);
    mpz_mul_2exp (t, t, 1UL);
    mpz_add_ui (t, t, 1UL);
    mpres_pow (state.r[1], r, t, state.modulus); /* r[1] = r^{2(-M+i)+1}, i = 0 */
    mpres_sqr (state.r[0], r, state.modulus); /* r[0] = r^2 */

    mpz_set_uint64 (t, P);
    mpz_mul (t, t, m_1);
    mpz_set_int64 (t1, k_2);
    mpz_add (t, t, t1);
    if (want_output)
      {
        outputf (OUTPUT_TRACE, "/* pm1_sequence_g */ %" PRId64 , k_2);
        outputf (OUTPUT_TRACE, " + %Zd*P == %Zd /* PARI C */\n", m_1, t);
      }

    mpres_pow (x_0, b_1, t, state.modulus);  /* x_0 = b_1^{k_2 + m_1*P} */
    if (want_output && test_verbose (OUTPUT_TRACE))
      {
        mpres_get_z (t, x_0, state.modulus);
        outputf (OUTPUT_TRACE, "/* pm1_sequence_g */ x_0 == %Zd /* PARI C */\n", 
                 t);
      }
    
    mpres_init (x_Mi, state.modulus);
    mpres_pow (x_Mi, x_0, tM, state.modulus); /* x_Mi = x_0^{M-i}, i = 0 */

    mpres_invert (x_0, x_0, state.modulus);  /* x_0 := x_0^{-1} now */
    mpres_mul (state.r[1], state.r[1], x_0, state.modulus); /* r[1] = x_0^{-1} * r^{-2M+1} */
    
    mpres_mul (r2, r2, x_Mi, state.modulus); /* r2 = x_0^M * r^{M^2} */
    mpres_get_z (state.t, r2, state.modulus);

    mpres_clear (x_Mi, state.modulus);
    mpres_clear (r2, state.modulus);
    mpres_clear (x_0, state.modulus);
    mpres_clear (r, state.modulus);
    mpz_clear (t1);
    mpz_clear (t);
    mpz_clear (tM);

    /* So here we have for i = 0
       t = x_0^(M-i) * r^{(M-i)^2}
       r[1] = x_0^{-1} * r^{2(-M+i)+1}
       r[0] = r^2
    */

    if (g_handle == NULL)
      {
        /* mpz version */
        ASSERT_ALWAYS (g_mpz != NULL);

        for (i = 0; i < l; i++)
          pm1_sequence_g_prod (&state, g_mpz[offset + i]);
      } else {
        /* NTT version */
        ASSERT_ALWAYS (g_mpz == NULL);

        mpzspv_fromto_mpzv (g_handle, offset, l, &pm1_sequence_g_prod, 
                                 &state, NULL, NULL);
      }

    mpres_clear (state.r[0], state.modulus);
    mpres_clear (state.r[1], state.modulus);
    mpz_clear (state.t);
    mpmod_clear (state.modulus); /* Clear our private copy of modulus */
  }

  print_elapsed_time (OUTPUT_VERBOSE, timestart, realstart);

  if (test_verbose (OUTPUT_TRACE))
    {
      uint64_t i;
      mpz_t mt;
      
      mpz_init (mt);
      for (i = 0; i < l_param; i++)
	{
          if (g_handle == NULL)
            mpz_set (mt, g_mpz[i]);
          else
            mpzspv_fromto_mpzv (g_handle, i, 1, NULL, NULL, NULL, &mt);
          outputf (OUTPUT_TRACE, "/* pm1_sequence_g */ g_%" PRIu64 
                   " = %Zd; /* PARI */\n", i, mt);
          outputf (OUTPUT_TRACE, "/* pm1_sequence_g */ g_%" PRIu64
                   " == x_0^(M - %" PRIu64 ") * r^((M - %" PRIu64 ")^2) "
                   "/* PARI C */\n", i, i, i);
	}
      mpz_clear (mt);
      outputf (OUTPUT_TRACE, "/* pm1_sequence_g */ g(x) = g_0");
      for (i = 1; i < l_param; i++)
	outputf (OUTPUT_TRACE, " + g_%" PRIu64 " * x^%" PRIu64 ,  i, i);
      outputf (OUTPUT_TRACE, " /* PARI */\n");
    }
}


typedef struct {
  mpmod_t modulus;
  mpres_t fd[3]; /* finite differences table for r^{-i^2}*/
  _listz_handle_t f;
  uint64_t index;
  file_word_t *buf;
  size_t bufsize;
} pm1_h_state_t;

static void
pm1_sequence_h_prod (void *state_p, mpz_t r)
{
  pm1_h_state_t *state = state_p;

  /* Get a coefficient of F, either from memory or file */
  if (state->f.storage == 0) {
    outputf (OUTPUT_TRACE, 
             "/* pm1_sequence_h */ f_%" PRIu64 " = %Zd; /* PARI */\n", 
             state->index, state->f.data.mem[state->index]);
    mpres_mul_z_to_z (r, state->fd[2], state->f.data.mem[state->index], 
      state->modulus);
  } else {
    seek_read_residue (state->f.data.file, r, state->buf, state->bufsize, 
      state->index);
    outputf (OUTPUT_TRACE, 
             "/* pm1_sequence_h */ f_%" PRIu64 " = %Zd; /* PARI */\n", 
             state->index, r);
    mpres_mul_z_to_z (r, state->fd[2], r, state->modulus);
  }
  outputf (OUTPUT_TRACE, 
           "/* pm1_sequence_h */ h_%" PRIu64 " = %Zd; /* PARI */\n", 
           state->index, r);
  mpres_mul (state->fd[2], state->fd[2], state->fd[1], state->modulus); 
    /* fd[2] = r^{-j^2} */
  mpres_mul (state->fd[1], state->fd[1], state->fd[0], state->modulus); 
    /* fd[1] = r^{-2*j-1} */
  state->index++;
}


/* Compute h_j = r^(-j^2) * f_j for 0 <= j < d as described in section 9 
   of the paper. h == f is ok. */

static void 
pm1_sequence_h (listz_handle_t h, mpzspv_handle_t ntt_handle, listz_handle_t f, 
                const mpres_t r, const uint64_t d, mpmod_t modulus_parm)
{
  mpres_t invr;  /* r^{-1}. Can be shared between threads */
  long timestart, realstart;

  mpres_init (invr, modulus_parm);
  mpres_invert (invr, r, modulus_parm); /* invr = r^{-1}. FIXME: test for 
					   failure, even if theoretically 
					   impossible */

  if (test_verbose (OUTPUT_TRACE))
    {
      mpz_t t;
      mpz_init (t);
      mpres_get_z (t, r, modulus_parm);
      outputf (OUTPUT_TRACE, "\n/* pm1_sequence_h */ N = %Zd; "
	       "r = Mod(%Zd, N); /* PARI */\n", 
	       modulus_parm->orig_modulus, t);
      mpz_clear (t);
    }

  outputf (OUTPUT_VERBOSE, "Computing h");
  outputf (OUTPUT_TRACE, "\n");
  timestart = cputime ();
  realstart = realtime ();

#ifdef _OPENMP
#pragma omp parallel if (d > 100)
#endif
  {
    pm1_h_state_t state;
    mpz_t t;       /* the h_j value as an mpz_t */
    uint64_t i, len, offset;

    /* Adjust offset and length for this thread */
    get_chunk (&offset, &len, d);
#ifdef _OPENMP
    if (omp_get_thread_num() == 0)
      outputf (OUTPUT_VERBOSE, " using %d threads", omp_get_num_threads());
#endif
    
    mpmod_init_set (state.modulus, modulus_parm);
    mpres_init (state.fd[0], state.modulus);
    mpres_init (state.fd[1], state.modulus);
    mpres_init (state.fd[2], state.modulus);
    state.index = offset;
    state.f.storage = f->storage;
    state.f.len = f->len;
    state.buf = (file_word_t *) mpz_export (NULL, &state.bufsize, -1, 
      sizeof(file_word_t), -1, 0, state.modulus->orig_modulus);
    if (f->storage == 0)
      {
        state.f.data.mem = f->data.mem;
      } else {
        /* FIXME Clone the file descriptor */
        state.f.data.file = f->data.file;
      }

    mpz_init (t);
    
    /* We have (n + 1)^2 = n^2 + 2n + 1. For the finite differences we'll 
       need r^{-2}, r^{-(2n+1)}, r^{-n^2}. Init for n = 0. */
    
    /* r^{-2} in fd[0] is constant and could be shared. Computing it 
       separately in each thread has the advantage of putting it in
       local memory. May not make much difference overall */

    mpres_sqr (state.fd[0], invr, state.modulus); /* fd[0] = r^{-2} */
    mpz_set_uint64 (t, offset);
    mpz_mul_2exp (t, t, 1UL);
    mpz_add_ui (t, t, 1UL);                 /* t = 2 * offset + 1 */
    mpres_pow (state.fd[1], invr, t, state.modulus);    /* fd[1] = r^{-(2*offset+1)} */
    mpz_set_uint64 (t, offset);
    mpz_mul (t, t, t);                      /* t = offset^2 */
    mpres_pow (state.fd[2], invr, t, state.modulus);    /* fd[2] = r^{-offset^2} */
    mpz_clear (t);
    
    /* Generate the sequence */
    
    if (ntt_handle == NULL)
      {
        /* mpz version */
        for (i = 0; i < len; i++)
          {
            pm1_sequence_h_prod (&state, h->data.mem[offset + i]);
            outputf (OUTPUT_TRACE, 
                     "/* pm1_sequence_h */ h_%lu = %Zd; /* PARI */\n", 
                     offset + i, h->data.mem[offset + i]);
          }
      } else {
        /* NTT version */
        /* FIXME: Need to print coefficients to make PARI check below work */
        mpzspv_fromto_mpzv (ntt_handle, offset, len, 
                            &pm1_sequence_h_prod, &state, NULL, NULL);
      }
    
    mpres_clear (state.fd[2], state.modulus);
    mpres_clear (state.fd[1], state.modulus);
    mpres_clear (state.fd[0], state.modulus);
    mpmod_clear (state.modulus);
    free (state.buf);
  }

  mpres_clear (invr, modulus_parm);

  print_elapsed_time (OUTPUT_VERBOSE, timestart, realstart);

  if (test_verbose (OUTPUT_TRACE))
    {
      uint64_t j;
      for (j = 0; j < d; j++)
	outputf (OUTPUT_TRACE, "/* pm1_sequence_h */ h_%" PRIu64
		   " == f_%" PRIu64 " * r^(-%" PRIu64 "^2) "
                   "/* PARI C */\n", j, j, j);
      outputf (OUTPUT_TRACE, "/* pm1_sequence_h */ h(x) = h_0");
      for (j = 1; j < d; j++)
        outputf (OUTPUT_TRACE, " + h_%" PRIu64 " * (x^%" PRIu64
                        " + x^(-%" PRIu64 "))", j, j, j);
      outputf (OUTPUT_TRACE, " /* PARI */\n");
    }
}


static int 
make_S_1_S_2 (set_list_t *S_1, int64_t **s2_sumset_out, 
	      uint64_t *s2_sumset_size_out,
              const faststage2_param_t *params)
{
  uint64_t i;
  set_list_t S_2;
  uint64_t s2_sumset_size;
  int64_t *s2_sumset;

  sets_get_factored_sorted (S_1, NULL, params->P);

  {
    mpz_t t1, t2;
    
    mpz_init (t1);
    mpz_init (t2);
    sets_sumset_minmax (t1, S_1, 1);
    sets_get_factored_sorted (NULL, t2, params->P);
    ASSERT_ALWAYS (mpz_cmp (t1, t2) == 0);
    mpz_clear (t1);
    mpz_clear (t2);
  }

  /* Extract sets for S_2 and compute the set of sums */
  
  sets_init(&S_2);
  sets_extract (&S_2, S_1, params->s_2);
  s2_sumset_size = sets_sumset_size(&S_2);
  s2_sumset = (int64_t *)malloc (s2_sumset_size * sizeof(int64_t));
  sets_sumset (s2_sumset, &S_2);
  
  /* Print the sets in devverbose mode */
  if (test_verbose (OUTPUT_DEVVERBOSE))
    {
      outputf (OUTPUT_DEVVERBOSE, "S_1 = ");
      sets_print (OUTPUT_DEVVERBOSE, S_1);
      
      outputf (OUTPUT_DEVVERBOSE, "S_2 = ");
      sets_print (OUTPUT_DEVVERBOSE, &S_2);
      
      outputf (OUTPUT_DEVVERBOSE, "S_2 sums = {");
      for (i = 0UL; i < s2_sumset_size - 1; i++)
	outputf (OUTPUT_DEVVERBOSE, "%" PRId64 ", ", s2_sumset[i]);
      outputf (OUTPUT_DEVVERBOSE, "%" PRId64 "}\n", s2_sumset[i]);
    }

  *s2_sumset_size_out = s2_sumset_size;
  *s2_sumset_out = s2_sumset;
  sets_free(&S_2);
  return 0;
}


/* Square the reciprocal Laurent polynomial S(x) of degree 2*n-2.
   S(x) = s_0 + \sum_{i=1}^{n-1} s_i (x^i + x^{-1}).
   S[i] contains the n coefficients s_i, 0 <= i <= n-1.
   R[i] will contain the 2n-1 coefficients r_i, 0 <= i <= 2*n-2, where 
   R(x) = S(x)^2 = r_0 + \sum_{i=1}^{2n-2} r_i (x^i + x^{-1}).
   dft must have power of 2 length len >= 2n.
   The NTT primes must be == 1 (mod 3*len).
*/

static void
ntt_sqr_reciprocal (mpzv_t R, const mpzv_t S, const spv_size_t n, 
                    const mpzspv_handle_t ntt_handle)
{
  spv_size_t i;
#ifdef WANT_ASSERT
  mpz_t S_eval_1, R_eval_1;
#endif
  
  if (n == 0)
    return;

  if (n == 1)
    {
      mpz_mul (R[0], S[0], S[0]);
      mpz_mod (R[0], R[0], ntt_handle->mpzspm->modulus);
      return;
    }

#ifdef WANT_ASSERT
  mpz_init (S_eval_1);
  list_recip_eval1 (S_eval_1, S, n);
  /* Compute (S(1))^2 */
  mpz_mul (S_eval_1, S_eval_1, S_eval_1);
  mpz_mod (S_eval_1, S_eval_1, ntt_handle->mpzspm->modulus);
#endif

  /* Fill NTT elements [0 .. n-1] with coefficients */
  mpzspv_fromto_mpzv (ntt_handle, (spv_size_t) 0, n, 
      NULL, S, NULL, NULL);
  
  mpzspv_sqr_reciprocal (ntt_handle, n);
  
  mpzspv_fromto_mpzv (ntt_handle, (spv_size_t) 0, 2*n - 1, 
      NULL, NULL, NULL, R);
  
  for (i = 0; i < 2*n - 1; i++)
    mpz_mod (R[i], R[i], ntt_handle->mpzspm->modulus);

#ifdef WANT_ASSERT
  mpz_init (R_eval_1);
  /* Compute (S^2)(1) and compare to (S(1))^2 */
  list_recip_eval1 (R_eval_1, R, 2 * n - 1);
  mpz_mod (R_eval_1, R_eval_1, ntt_handle->mpzspm->modulus);
  if (mpz_cmp (R_eval_1, S_eval_1) != 0)
    {
      gmp_fprintf (stderr, "ntt_sqr_reciprocal: (S(1))^2 = %Zd but "
		   "(S^2)(1) = %Zd\n", S_eval_1, R_eval_1);
#if 0
      list_output_poly (R, 2*n-1, 0, 1, "Output polynomial is ", "\n", 
                        OUTPUT_TRACE);
#endif
      abort ();
    }
  mpz_clear (S_eval_1);
  mpz_clear (R_eval_1);
#endif
}


typedef struct {
  mpmod_t modulus;
  mpres_t prod, tmpres;
  mpz_t sum;
  file_word_t *buf;
  listz_handle_t add;
  uint64_t offset;
} gcd_state_t;

static void
gcd_consumer (void *state_p, const mpz_t s)
{
  gcd_state_t *state = state_p;
  if (state->add != NULL)
    {
      listz_handle_get (state->add, state->sum, state->buf, state->offset);
      mpz_add (state->sum, state->sum, s);
      mpres_set_z_for_gcd (state->tmpres, state->sum, state->modulus);
    } else {
      mpres_set_z_for_gcd (state->tmpres, s, state->modulus);
    }
#define TEST_ZERO_RESULT
#ifdef TEST_ZERO_RESULT
  if (mpres_is_zero (state->tmpres, state->modulus))
    outputf (OUTPUT_VERBOSE, "R_[%" PRIu64 "] = 0\n", state->offset);
#endif
  state->offset++;
  mpres_mul (state->prod, state->prod, state->tmpres, state->modulus); 
}

/* Computes gcd(\prod_{0 <= i < len} (ntt[i + offset] + add[i]), N), 
   the NTT residues are converted to integer residues (mod N) first.
   If add == NULL, add[i] is assumed to be 0. */

static void
ntt_gcd (mpz_t f, mpz_t *product, mpzspv_handle_t ntt, 
         const uint64_t ntt_offset, const listz_handle_t add, 
         const uint64_t len_param, mpmod_t modulus_param)
{
  mpres_t totalprod;
  long timestart, realstart;
  
  outputf (OUTPUT_VERBOSE, "Computing gcd of coefficients and N");
  timestart = cputime ();
  realstart = realtime ();

  /* All the threads will multiply their partial products to this one. */
  mpres_init (totalprod, modulus_param);
  mpres_set_ui (totalprod, 1UL, modulus_param);

#ifdef _OPENMP
#pragma omp parallel if (len_param > 100) shared(totalprod)
  {
#endif
    gcd_state_t state;
    uint64_t len, thread_offset;

#ifdef _OPENMP
#pragma omp master
    {
      outputf (OUTPUT_VERBOSE, " using %d threads", omp_get_num_threads());
    }
#endif

    get_chunk (&thread_offset, &len, len_param);
    /* Make a thread-private copy of the mpmod_t struct */
    mpmod_init_set (state.modulus, modulus_param);
    mpres_init (state.prod, state.modulus);
    mpres_init (state.tmpres, state.modulus);
    mpz_init (state.sum);
    mpres_set_ui (state.prod, 1UL, state.modulus);
    state.add = add;
    state.offset = thread_offset;
    if (add != NULL)
      state.buf = (file_word_t *) malloc (add->words * sizeof(file_word_t));
    else
      state.buf = NULL;
    
    mpzspv_fromto_mpzv (ntt, ntt_offset + thread_offset, len, 
        NULL, NULL, &gcd_consumer, &state);

#ifdef _OPENMP
#pragma omp critical
    {
      mpres_mul (totalprod, totalprod, state.prod, state.modulus);
    }
#else
    mpres_set (totalprod, state.prod, state.modulus);
#endif
    free (state.buf);
    mpres_clear (state.tmpres, state.modulus);
    mpres_clear (state.prod, state.modulus);
    mpmod_clear (state.modulus);
    mpz_clear (state.sum);
#ifdef _OPENMP
  }
#endif

  {
    mpz_t n;
    mpz_init (n);
    mpz_set_ui (n, len_param);
    mpres_set_z_for_gcd_fix (totalprod, totalprod, n, modulus_param);
    mpz_clear (n);
  }

  if (product != NULL)
    mpres_get_z (*product, totalprod, modulus_param);

  mpres_gcd (f, totalprod, modulus_param);
  mpres_clear (totalprod, modulus_param);

  print_elapsed_time (OUTPUT_VERBOSE, timestart, realstart);
}


/* A slow method of evaluating the polynomial on a single point t.
   We iterate through the elements k_1 of S_1 and collect the product 
   \prod_{k_1 \in S_1} (x - P_1^{k_1}) (mod x - t)
   = \prod_{k_1 \in S_1} (t - P_1^{k_1})
   
   = \prod_{k_1 \in S_1, k_1 > 0} (x - P_1^{k_1})(x - P_1^{-k_1}) (mod x - t)
   = \prod_{k_1 \in S_1, k_1 > 0} x (x+1/x - V_{k_1}(P_1 + 1/P_1)) (mod x - t)
   = \prod_{k_1 \in S_1, k_1 > 0} ((t+1/t) - V_{k_1}(P_1 + 1/P_1))
   where both t+1/t and V_{k_1}(P_1 + 1/P_1)) are in the base ring.
*/

static void 
pm1_eval_slow (const mpz_t checkval, const set_list_t *S1, const mpz_t b1, 
               const uint64_t P, const mpz_t m1, const uint64_t m, 
               const int64_t k2, const mpz_t N)
{
  uint32_t *iterator;
  unsigned long timestart, realstart;
  mpz_t B1B, V2B, x0, X, X1X, r, e, product, tmp;
  int extgcd_ok;
  uint64_t degree = 0, maxS, i, evenodd[2]; /* even = [0], odd = [1] */
  char *bitfield;
  int j;
  
  ASSERT_ALWAYS (mpz_sgn (checkval) >= 0);
  ASSERT_ALWAYS (mpz_cmp (checkval, N) < 0);
  
  outputf (OUTPUT_VERBOSE, "Recomputing one polynomial value");
  outputf (OUTPUT_TRACE, "\nN = %Zd; b1 = Mod(%Zd,N); P = %l" PRIu64 
           "; m1 = %Zd; m = %" PRIu64 "; k_2 = %" PRId64 "; /* PARI %s */\n", 
           N, b1, P, m1, m, k2, __func__);

  timestart = cputime ();
  realstart = realtime ();

  mpz_init (e);
  mpz_init (x0);
  mpz_init (r);
  mpz_init (X);
  mpz_init (X1X);
  mpz_init (B1B);
  mpz_init (V2B);
  mpz_init (tmp);
  mpz_init (product);
  mpz_set_ui (product, 1);

  sets_sumset_minmax (tmp, S1, 1);
  maxS = mpz_get_uint64 (tmp);
  bitfield = malloc (maxS / 8 + 1);
  memset (bitfield, 0, maxS / 8 + 1);
  ASSERT_ALWAYS(bitfield != NULL);
  iterator = sets_init_iterator (S1);
  ASSERT_ALWAYS(iterator != NULL);
  evenodd[0] = evenodd[1] = 0;
  while (!sets_end_of_iter (iterator, S1))
    {
      const int64_t k1 = sets_next_iter (iterator, S1);
      if (k1 >= 0)
        {
          evenodd[0] |= k1 ^ 1; /* LSB indicates whether there are any even k1 */
          evenodd[1] |= k1;     /* ... any odd k1 */
          ASSERT_ALWAYS ((uint64_t) k1 <= maxS);
          ASSERT_ALWAYS ((bitfield[k1 / 8] & (1 << (k1 % 8))) == 0);
          /* printf ("pm1_eval_slow(): k1 = %" PRId64 "\n", k1); */
          bitfield[k1 / 8] |= 1 << (k1 % 8);
        }
    }
  free (iterator);
  ASSERT_ALWAYS ((bitfield[maxS / 8] & (1 << (maxS % 8))) != 0);

  evenodd[0] &= 1;
  evenodd[1] &= 1;
  outputf (OUTPUT_DEVVERBOSE, "\nmaxS = %" PRIu64 ", even = %d, odd = %d\n", 
          maxS, evenodd[0] != 0 ? 1 : 0, evenodd[1] != 0 ? 1 : 0);

  /* See section 8 of ANTS paper */
  mpz_set_uint64 (e, P);
  mpz_mul (e, e, m1); /* e = m_1 * P */
  mpz_set_int64 (tmp, k2);
  mpz_add (e, e, tmp); /* e = k_2 + m_1 P */
  mpz_powm (x0, b1, e, N); /* x0 = b_1^{k_2 + m_1 * P} */

  mpz_set_uint64 (e, P/2);
  mpz_powm (r, b1, e, N); /* r = b1^(P/2) */
  mpz_set_uint64 (e, 2*m);
  mpz_powm (X, r, e, N);
  mpz_mul (X, X, x0); /* X = x0 * r^m = b_1^{k_2 + m_1 * P} * b_1^(P*m) = 
                              b_1^{k_2 + (m_1 + m) * P} */
  mpz_mod (X, X, N);

  extgcd_ok = mpz_invert (X1X, X, N);
  ASSERT_ALWAYS (extgcd_ok != 0);
  mpz_add (X1X, X1X, X);
  mpz_mod (X1X, X1X, N); /* X1X = X + 1/X */

  outputf (OUTPUT_TRACE, 
           "x0 = Mod(%Zd,N); X = Mod(%Zd, N); r = Mod(%Zd, N); /* PARI %s */\n", 
           x0, X, r, __func__);
  outputf (OUTPUT_TRACE, "x0 == b1^(k_2 + m1*P) /* PARI C %s */ \n", 
           __func__);
  outputf (OUTPUT_TRACE, "X == b1^(k_2 + (m1+m)*P) /* PARI C %s */ \n", 
           __func__);
  outputf (OUTPUT_TRACE, "r == b1^(P/2) /* PARI C %s */ \n", __func__);
  outputf (OUTPUT_TRACE, "X + 1/X == %Zd /* PARI C %s */ \n", __func__, X1X);

  extgcd_ok = mpz_invert (B1B, b1, N);
  ASSERT_ALWAYS (extgcd_ok != 0);
  mpz_add (B1B, B1B, b1);
  mpz_mod (B1B, B1B, N);

  mpz_mul (V2B, B1B, B1B);
  mpz_sub_ui (V2B, V2B, 2); /* V_2(x) = x^2 - 2, V2B = V_2(B1B) */
  
  for (j = 0; j < 2; j++)
    if (evenodd[j] != 0)
      {
        mpz_t ViB, Vim2B;
        /* We start at index i = j (either 0 or 1) */
        mpz_init (ViB);
        mpz_init (Vim2B);
        if (j == 0)
          {
            mpz_set_ui (ViB, 2); /* ViB = V_i(B1B) = 2 with i=0 */
            mpz_set (Vim2B, V2B);  /* Vim2B = V_{i-2}(B1B) = V_2(B1B) with i=0 */
          } else {
            mpz_set (ViB, B1B); /* ViB = V_i(B1B) = B1B with i=1 */
            mpz_set (Vim2B, B1B); /* Vim2B = V_{i-2}(B1B) = B1B with i=1 */
          }

        for (i = j; i <= maxS; i += 2)
          {
            if ((bitfield[i / 8] & (1 << (i % 8))) != 0)
              {
                /* printf ("pm1_eval_slow(): i = %" PRId64 "\n", i); */
                /* Multiply by (X - b1^i)*(X - b1^{-1})/X = 
                               (X + 1/X - (b1^i + b^1{-i})) =
                               (X + 1/X - V_i(b1 + b1^{-1})) */
                degree+=2;

                mpz_sub (tmp, X1X, ViB);
                mpz_mul (product, product, tmp);
                mpz_mod (product, product, N);
              }
            
            mpz_mul (tmp, ViB, V2B);
            mpz_sub (tmp, tmp, Vim2B);
            mpz_set (Vim2B, ViB);
            mpz_mod (ViB, tmp, N);
          }
        mpz_clear (ViB);
        mpz_clear (Vim2B);
      }

  outputf (OUTPUT_DEVVERBOSE, "degree = %" PRIu64 ", maxS / degree = %f\n", 
           degree, maxS * 1./degree);
  outputf (OUTPUT_TRACE, "F(X) == Mod(%Zd,N) /* PARI C %s */\n", 
           product, __func__);

  /* Compute x_0^m r^(m^2) * F(X) */
  mpz_set_uint64 (e, m);
  mpz_powm (tmp, x0, e, N);
  mpz_mul (product, product, tmp);
  mpz_mod (product, product, N);
  mpz_mul (e, e, e);
  mpz_powm (tmp, r, e, N);
  mpz_mul (product, product, tmp);
  mpz_mod (product, product, N);
  outputf (OUTPUT_TRACE, "x0^m * r^(m^2) * F(X) == Mod(%Zd,N) /* PARI C %s */\n", 
           product, __func__);

  if (mpz_cmp (product, checkval) != 0)
    {
      gmp_printf ("Error, slow result = %Zd, multi-point = %Zd for m = %" PRIu64 "\n", 
                  product, checkval, m);
      abort();
    }

  mpz_clear (product);
  mpz_clear (e);
  mpz_clear (x0);
  mpz_clear (r);
  mpz_clear (X);
  mpz_clear (V2B);
  mpz_clear (B1B);
  mpz_clear (X1X);
  mpz_clear (tmp);
  free(bitfield);
  bitfield = NULL;
  
  print_elapsed_time (OUTPUT_VERBOSE, timestart, realstart);
}


int 
pm1fs2 (mpz_t f, const mpres_t X, mpmod_t modulus, 
	const faststage2_param_t *params)
{
  uint64_t nr;
  uint64_t i, l, lenG, lenR, tmplen;
  set_list_t S_1; /* This is stored as a set of sets (arithmetic 
                     progressions of prime length */
  int64_t *s2_sumset; /* set of sums of S_2 */
  uint64_t s2_sumset_size;
  listz_handle_t F;   /* Polynomial F has roots X^{k_1} for k_1 \in S_1, so has 
		  degree s_1. It is symmetric, so has only s_1 / 2 + 1 
		  distinct coefficients. The sequence h_j will be stored in 
		  the same memory and won't be a monic polynomial, so the 
		  leading 1 monomial of F will be stored explicitly. Hence we 
		  need s_1 / 2 + 1 entries. */
  listz_t g, h, tmp, R;
  mpz_t mt;   /* All-purpose temp mpz_t */
  mpres_t mr; /* All-purpose temp mpres_t */
  int youpi = ECM_NO_FACTOR_FOUND;
  long timetotalstart, realtotalstart, timestart;

  timetotalstart = cputime ();
  realtotalstart = realtime ();

  ASSERT_ALWAYS (eulerphi64 (params->P) == params->s_1 * params->s_2);
  ASSERT_ALWAYS (params->s_1 < params->l);
  nr = params->l - params->s_1; /* Number of points we evaluate */

  sets_init(&S_1);

  if (make_S_1_S_2 (&S_1, &s2_sumset, 
		  &s2_sumset_size, params) == ECM_ERROR)
      return ECM_ERROR;

  /* Allocate all the memory we'll need */
  /* Allocate the correct amount of space for each mpz_t or the 
     reallocations will up to double the time for stage 2! */
  mpz_init (mt);
  F = listz_handle_init2 (NULL, params->s_1 / 2 + 1 + 1, 
    modulus->orig_modulus); /* Another +1 in len because 
    poly_from_sets_V stores the leading 1 monomial for each factor */
  h = malloc ((params->s_1 + 1) * sizeof (mpz_t));
  if (h == NULL)
    {
      fprintf (stderr, "Cannot allocate memory in pm1fs2\n");
      exit (1);
    }
  lenG = params->l;
  g = init_list2 (lenG, (unsigned int) abs (modulus->bits));
  lenR = nr;
  R = init_list2 (lenR, (unsigned int) abs (modulus->bits));    
  tmplen = 3 * params->l + list_mul_mem (params->l / 2);
  outputf (OUTPUT_DEVVERBOSE, "tmplen = %" PRIu64 "\n", tmplen);
  if (TMulGen_space (params->l - 1, params->s_1, lenR) + 12 > tmplen)
    {
      tmplen = TMulGen_space (params->l - 1, params->s_1 - 1, lenR) + 12;
      /* FIXME: It appears TMulGen_space() returns a too small value! */
      outputf (OUTPUT_DEVVERBOSE, "With TMulGen_space, tmplen = %lu\n", 
	       tmplen);
    }
#ifdef SHOW_TMP_USAGE
  tmp = init_list (tmplen);
#else
  tmp = init_list2 (tmplen, (unsigned int) abs (modulus->bits));
#endif
  
  mpres_get_z (mt, X, modulus); /* mpz_t copy of X for printing */
  outputf (OUTPUT_TRACE, 
	   "N = %Zd; X = Mod(%Zd, N); /* PARI */\n", 
	   modulus->orig_modulus, mt);


  /* Compute the polynomial f(x) = \prod_{k_1 in S_1} (x - X^{2k_1}) */
  outputf (OUTPUT_VERBOSE, "Computing F from factored S_1");
  
  timestart = cputime ();
  
  /* First compute X + 1/X */
  mpres_init (mr, modulus);
  mpres_invert (mr, X, modulus);
  mpres_add (mr, mr, X, modulus);
  
  i = poly_from_sets_V (F, mr, &S_1, tmp, tmplen, modulus, NULL);
  ASSERT_ALWAYS(2 * i == params->s_1);
  ASSERT(mpz_cmp_ui (F->data.mem[i], 1UL) == 0);
  
  outputf (OUTPUT_VERBOSE, " took %lums\n", cputime () - timestart);
  list_output_poly_file (F, params->s_1 / 2 + 1, 0, 1, "f(x) = ", 
                         "; /* PARI */ \n", OUTPUT_TRACE);
  
  mpz_set_ui (mt, params->P/2);
  mpres_pow (mr, X, mt, modulus); /* mr = X^(P/2) */
  pm1_sequence_h (F, NULL, F, mr, params->s_1 / 2 + 1, modulus); 

  /* Make a symmetric copy of F in h. It will have length 
     s_1 + 1 = 2*lenF - 1 */
  /* I.e. with F = [3, 2, 1], s_1 = 4, we want h = [1, 2, 3, 2, 1] */
  for (i = 0; i < params->s_1 / 2 + 1; i++)
    *(h[i]) = *(F->data.mem[params->s_1 / 2 - i]); /* Clone the mpz_t. */
  for (i = 0; i < params->s_1 / 2; i++)
    *(h[i + params->s_1 / 2 + 1]) = *(F->data.mem[i + 1]);

  list_output_poly (h, params->s_1, 0, 0, "h(x) = (", 
                    ") / x^((s_1 - 1) / 2) /* PARI */\n", OUTPUT_TRACE);
  
  for (l = 0; l < params->s_2 && youpi == ECM_NO_FACTOR_FOUND; l++)
    {
      const uint64_t M = params->l - 1L - params->s_1 / 2L;
      outputf (OUTPUT_VERBOSE, "Multi-point evaluation %" PRIu64
                        " of %" PRIu64 ":\n", l + 1, params->s_2);
      pm1_sequence_g (g, NULL, X, params->P, M, params->l, 
		      params->m_1, s2_sumset[l], modulus);

      /* Do the convolution */
      /* Use the transposed "Middle Product" algorithm */
      /* TMulGen reverses the first input sequence, but that doesn't matter
	 since h is symmetric. */

      outputf (OUTPUT_VERBOSE, "TMulGen of g and h");
      timestart = cputime ();
      ASSERT(tmplen >= TMulGen_space (nr - 1, params->l - 1, params->s_1));

      /* Computes rev(h)*g, stores coefficients of x^(s_1) to 
	 x^(s_1+nr-1) = x^(len-1) */
      if (TMulGen (R, nr - 1, h, params->s_1, g, params->l - 1, tmp, 
		   modulus->orig_modulus) < 0)
	{
	  outputf (OUTPUT_ERROR, "TMulGen returned error code (probably out "
		   "of memory)\n");
	  youpi = ECM_ERROR;
	  break;
	}
      list_mod (R, R, nr, modulus->orig_modulus);

      outputf (OUTPUT_VERBOSE, " took %lums\n", cputime () - timestart);

      if (test_verbose (OUTPUT_TRACE))
	{
	  for (i = 0; i < nr; i++)
	    outputf (OUTPUT_TRACE, "r_%" PRIu64 " = %Zd; /* PARI */\n", 
                                i, R[i]);
	}

      outputf (OUTPUT_VERBOSE, "Computing product of F(g_i)");
      timestart = cputime ();

      {
	mpres_t tmpres, tmpprod;
	mpres_init (tmpres, modulus);
	mpres_init (tmpprod, modulus);
	mpres_set_z_for_gcd (tmpprod, R[0], modulus);
	for (i = 1; i < nr; i++)
	  {
	    mpres_set_z_for_gcd (tmpres, R[i], modulus);
	    mpres_mul (tmpprod, tmpprod, tmpres, modulus); 
	  }
        mpz_set_uint64 (mt, nr);
        mpres_set_z_for_gcd_fix (tmpprod, tmpprod, mt, modulus);
        mpres_get_z (tmp[1], tmpprod, modulus); /* For printing */
	mpres_gcd (tmp[0], tmpprod, modulus);
	mpres_clear (tmpprod, modulus);
	mpres_clear (tmpres, modulus);
      }

      outputf (OUTPUT_VERBOSE, " took %lums\n", cputime () - timestart);
      outputf (OUTPUT_RESVERBOSE, "Product of R[i] = %Zd\n", tmp[1]);

      if (mpz_cmp_ui (tmp[0], 1UL) > 0)
	{
	  mpz_set (f, tmp[0]);
	  youpi = ECM_FACTOR_FOUND_STEP2;
	}

      /* Evaluate poly for comparison */
      if (check_eval)
        {
          const uint64_t m = 0; /* 0 <= m < nr */
          mpres_get_z (mt, X, modulus);
          pm1_eval_slow (R[nr-1-m], &S_1, mt, params->P, params->m_1, m, 
                         s2_sumset[l], modulus->orig_modulus);
        }
    }

#ifdef SHOW_TMP_USAGE
  for (i = tmplen - 1; i > 0; i--)
    if (tmp[i]->_mp_alloc > 1)
      break;
  outputf (OUTPUT_DEVVERBOSE, "Highest used temp element is tmp[%lu]\n", i);
#endif
  
  sets_free(&S_1);
  free (s2_sumset);
  free (h);
  listz_handle_clear (F);
  clear_list (g, lenG);
  clear_list (R, lenR);    
  clear_list (tmp, tmplen);

  mpz_clear (mt);
  mpres_clear (mr, modulus);

  outputf (OUTPUT_NORMAL, "Step 2");
  /* In normal output mode, print only cpu time as we always have.
     In verbose mode, print real time as well if we used multi-threading */
  if (test_verbose (OUTPUT_VERBOSE))
    print_elapsed_time (OUTPUT_NORMAL, timetotalstart, realtotalstart);
  else
    print_elapsed_time (OUTPUT_NORMAL, timetotalstart, 0L);
  
  return youpi;
}

static void
do_aio_init(ATTRIBUTE_UNUSED const unsigned int sp_num)
{
#ifdef HAVE_AIO_INIT
  /* Set the aio library to use only 1 single thread, as the default 
     of up to 20 threads causes a lot of reads to be executed 
     concurrently, interrupting each other's long sequential reads */
  struct aioinit init;
  memset (&init, 0, sizeof(struct aioinit));
  init.aio_threads = 1;
  init.aio_num = sp_num;
  aio_init (&init);
#endif
}
  

int 
pm1fs2_ntt (mpz_t f, const mpres_t X, mpmod_t modulus, 
	const faststage2_param_t *params)
{
  uint64_t nr;
  uint64_t l;
  set_list_t S_1; /* This is stored as a set of sets (arithmetic 
                       progressions of prime length */
  int64_t *s2_sumset; /* set of sums of S_2 */ 
  uint64_t s2_sumset_size;
  listz_handle_t F;   /* Polynomial F has roots X^{k_1} for k_1 \in S_1, so has 
		  degree s_1. It is symmetric, so has only s_1 / 2 + 1 
		  distinct coefficients. The sequence h_j will be stored in 
		  the same memory and won't be a monic polynomial, so the 
		  leading 1 monomial of F will be stored explicitly. Hence we 
		  need s_1 / 2 + 1 entries. */
  mpzspm_t ntt_context;
  mpzspv_handle_t g_handle, h_handle;
  char *h_filename = NULL, *g_filename = NULL;
  mpz_t mt;   /* All-purpose temp mpz_t */
  mpz_t product; /* Product of each multi-point evaluation */
  mpz_t *product_ptr = NULL;
  mpres_t XP, Q; /* XP = X^P, Q = 1 + 1/X */
  int youpi = ECM_NO_FACTOR_FOUND;
  long timetotalstart, realtotalstart, timestart, realstart;

  timetotalstart = cputime ();
  realtotalstart = realtime ();

  ASSERT_ALWAYS (eulerphi64 (params->P) == params->s_1 * params->s_2);
  ASSERT_ALWAYS (params->s_1 < params->l);
  nr = params->l - params->s_1; /* Number of points we evaluate */

  /* Prepare NTT for computing the h sequence, its DCT-I, and the convolution 
     with g. We need NTT of transform length l. We do it here at the start
     of stage 2 so that in case of a "not enough primes" condition, 
     we don't have to wait until after F is built to get the error. */

  ntt_context = mpzspm_init (params->l, modulus->orig_modulus);
  if (ntt_context == NULL)
    {
      outputf (OUTPUT_ERROR, "Could not initialise ntt_context, "
               "presumably out of memory\n");
      return ECM_ERROR;
    }

  mpzspm_print_CRT_primes (OUTPUT_DEVVERBOSE, "CRT modulus for evaluation = ", 
		    ntt_context);

  sets_init(&S_1);

  if (make_S_1_S_2 (&S_1, &s2_sumset, &s2_sumset_size, params) == ECM_ERROR)
      return ECM_ERROR;

  if (params->file_stem != NULL)
    {
      do_aio_init (ntt_context->sp_num);
      g_filename = malloc ((strlen(params->file_stem) + 3) * sizeof (char));
      h_filename = malloc ((strlen(params->file_stem) + 3) * sizeof (char));
      if (g_filename == NULL || h_filename == NULL)
        {
          fprintf (stderr, 
                   "pm1fs2_ntt(): could not allocate memory for filename\n");
          free (g_filename);
          free (h_filename);
          mpzspm_clear (ntt_context);
          return ECM_ERROR;
        }
      sprintf (g_filename, "%s.g", params->file_stem);
      sprintf (h_filename, "%s.h", params->file_stem);
    }

  /* Compute Q = X + 1/X */
  mpres_init (Q, modulus);
  mpres_invert (Q, X, modulus);
  mpres_add (Q, Q, X, modulus);

  mpz_init (mt);
  mpres_get_z (mt, X, modulus); /* mpz_t copy of X for printing */
  outputf (OUTPUT_TRACE, 
	   "N = %Zd; X = Mod(%Zd, N); /* PARI */\n", 
	   modulus->orig_modulus, mt);

#if 0 && defined (WANT_ASSERT)
  /* For this self test run with a large enough B2 so that enough memory
     is allocated for tmp and F_ntt, otherwise it segfaults. */
  {
    int testlen = 255;
    int i, j;
    /* A test of ntt_sqr_reciprocal() */
    for (j = 1; j <= testlen; j++)
      {
        outputf (OUTPUT_VERBOSE, 
                 "Testing ntt_sqr_reciprocal() for input degree %d\n", 
                 j - 1);
        for (i = 0; i < j; i++)
          mpz_set_ui (tmp[i], 1UL);
        ntt_sqr_reciprocal (tmp, tmp, F_ntt, (spv_size_t) j, ntt_context_F);
        for (i = 0; i < 2 * j - 1; i++)
          {
            ASSERT (mpz_cmp_ui (tmp[i], 2 * j - 1 - i) == 0);
          }
      }
    outputf (OUTPUT_VERBOSE, 
             "Test of ntt_sqr_reciprocal() for input degree 2 ... %d passed\n", 
             testlen - 1);
  }
#endif

  F = build_F_ntt (Q, &S_1, params, modulus);
  if (F == NULL)
    {
      sets_free (&S_1);
      free (s2_sumset);
      mpz_clear (mt);
      mpres_clear (Q, modulus);
      mpzspm_clear (ntt_context);
      return ECM_ERROR;
    }


  mpres_clear (Q, modulus);
  
  h_handle = mpzspv_init_handle (h_filename, params->l / 2 + 1, ntt_context);
  free (h_filename);
  h_filename = NULL;

  /* Compute XP = X^(P/2) */
  mpres_init (XP, modulus);
  mpz_set_ui (mt, params->P/2);
  mpres_pow (XP, X, mt, modulus);

  pm1_sequence_h (NULL, h_handle, F, XP, params->s_1 / 2 + 1, modulus);

  listz_handle_clear (F);
  F = NULL;
  mpres_clear (XP, modulus);

  /* Compute the DCT-I of h */
  outputf (OUTPUT_VERBOSE, "Computing DCT-I of h");
#ifdef _OPENMP
  outputf (OUTPUT_VERBOSE, " using %d threads", omp_get_thread_limit());
#endif
  if (test_verbose (OUTPUT_TRACE))
    {
      mpzspv_print (h_handle, 0, params->s_1 / 2 + 1, "h_ntt");
    }

  timestart = cputime ();
  realstart = realtime ();
  
  mpzspv_to_dct1 (h_handle, h_handle, params->s_1 / 2 + 1, params->l / 2 + 1);

  print_elapsed_time (OUTPUT_VERBOSE, timestart, realstart);
  
  if (test_verbose (OUTPUT_TRACE))
    {
      mpzspv_print (h_handle, 0, params->s_1 / 2 + 1, "DCT-I of h_ntt");
    }

  if (test_verbose (OUTPUT_RESVERBOSE))
    {
      mpz_init (product);
      product_ptr = &product;
    }

  g_handle = mpzspv_init_handle (g_filename, params->l, ntt_context);
  free (g_filename);
  g_filename = NULL;

  for (l = 0; l < params->s_2 && youpi == ECM_NO_FACTOR_FOUND; l++)
    {
      const uint64_t M = params->l - 1L - params->s_1 / 2L;

      outputf (OUTPUT_VERBOSE, "Multi-point evaluation %" PRIu64 
                        " of %" PRIu64 ":\n", l + 1, params->s_2);
      /* Compute the coefficients of the polynomial g(x) */
      pm1_sequence_g (NULL, g_handle, X, params->P, M, params->l, 
		      params->m_1, s2_sumset[l], modulus);

      if (test_verbose (OUTPUT_TRACE))
        mpzspv_print (g_handle, 0, params->s_1 / 2 + 1, "g_ntt");

      /* Do the convolution */
      outputf (OUTPUT_VERBOSE, "Computing g*h");
#ifdef _OPENMP
      outputf (OUTPUT_VERBOSE, " using %d threads", omp_get_thread_limit());
#endif
      timestart = cputime ();
      realstart = realtime ();
      mpzspv_mul_ntt (g_handle, 0, g_handle, 0, params->l, 
          h_handle, 0, params->l / 2 + 1, params->l, 0, 0, 
          NTT_MUL_STEP_FFT1 + NTT_MUL_STEP_MULDCT + NTT_MUL_STEP_IFFT);
      print_elapsed_time (OUTPUT_VERBOSE, timestart, realstart);
      
      if (test_verbose (OUTPUT_TRACE))
        mpzspv_print (g_handle, 0, params->l, "g_ntt * h_ntt");

      if (0)
        mpzspv_fromto_mpzv (g_handle, 0, params->l, NULL, NULL, 
          (mpz_consumerfunc_t) &gmp_printf, "g(x) * h(x) = %Zd\n"); /* hax! */

      /* Compute GCD of N and coefficients of product polynomial */
      ntt_gcd (mt, product_ptr, g_handle, params->s_1 / 2, NULL, nr, modulus);

      /* If we found a factor, stop */
      if (mpz_cmp_ui (mt, 1UL) > 0)
	{
	  mpz_set (f, mt);
	  youpi = ECM_FACTOR_FOUND_STEP2;
	}

      outputf (OUTPUT_RESVERBOSE, "Product of R[i] = %Zd\n", product);

      /* Evaluate poly for comparison */
      if (check_eval)
        {
          mpz_t checkval;
          const uint64_t m = 0; /* 0 <= m < nr */

          mpz_init (checkval);
          mpzspv_fromto_mpzv (g_handle, params->s_1 / 2 + nr - 1 - m, 1, NULL, 
                              NULL, NULL, &checkval);
          mpres_get_z (mt, X, modulus);
          pm1_eval_slow (checkval, &S_1, mt, params->P, params->m_1, m, 
                         s2_sumset[l], modulus->orig_modulus);
          mpz_clear (checkval);
        }

    }

  sets_free (&S_1);

  if (test_verbose (OUTPUT_RESVERBOSE))
    {
      product_ptr = NULL;
      mpz_clear (product);
    }
  mpzspv_clear_handle (h_handle);
  mpzspv_clear_handle (g_handle);
  mpzspm_clear (ntt_context);
  mpz_clear (mt);
  free (s2_sumset);

  outputf (OUTPUT_NORMAL, "Step 2");
  /* In normal output mode, print only cpu time as we always have.
     In verbose mode, print real time as well if we used multi-threading */
  if (test_verbose (OUTPUT_VERBOSE))
    print_elapsed_time (OUTPUT_NORMAL, timetotalstart, realtotalstart);
  else
    print_elapsed_time (OUTPUT_NORMAL, timetotalstart, 0L);
  
  return youpi;
}


/************************************************************************/
/* Functions for arithmetic in a quadratic extension in Z/NZ            */
/************************************************************************/

typedef struct {
  mpres_t x,y;
} _gfp_ext_t;
typedef _gfp_ext_t gfp_ext_t[1];

static void 
gfp_ext_init (gfp_ext_t r, mpmod_t modulus)
{
  mpres_init (r->x, modulus);
  mpres_init (r->y, modulus);
}

static void 
gfp_ext_clear (gfp_ext_t r, mpmod_t modulus)
{
  mpres_clear (r->x, modulus);
  mpres_clear (r->y, modulus);
}

static void 
gfp_ext_print (const gfp_ext_t r, mpmod_t modulus, const int verbose)
{
  mpz_t t1, t2;

  if (!test_verbose (verbose))
    return;

  mpz_init (t1);
  mpz_init (t2);
  mpres_get_z (t1, r->x, modulus);
  mpres_get_z (t2, r->y, modulus);
  outputf (verbose, "Mod(%Zd, N) + Mod(%Zd, N) * w", t1, t2);
  
  mpz_clear (t1);
  mpz_clear (t2);
}



/* Multiplies (a_0 + a_1*sqrt(Delta)) * (b_0 + b_1*sqrt(Delta))
   using four multiplications. Result goes in (r_0 + r_1*sqrt(Delta)). 
   a_0, b_0, r_0 as well as a_1, b_1, r_1 may overlap arbitrarily. t[0], t[1], 
   t[2] and Delta must not overlap with anything. */
/* FIXME: is there a faster multiplication routine if both inputs have 
   norm 1? */

static void 
gfp_ext_mul (gfp_ext_t r, const gfp_ext_t a, const gfp_ext_t b, 
             const mpres_t Delta, 
	     mpmod_t modulus, ATTRIBUTE_UNUSED const uint64_t tmplen, 
	     mpres_t *tmp)
{
  ASSERT (tmplen >= 2UL);
  if (0 && test_verbose (OUTPUT_TRACE))
    {
      mpz_t t;
      mpz_init (t);
      mpres_get_z (t, Delta, modulus);
      outputf (OUTPUT_TRACE, "/* gfp_ext_mul */ w = quadgen (4*%Zd); "
               "N = %Zd; /* PARI */\n", t, modulus->orig_modulus);
      mpz_clear (t);
      outputf (OUTPUT_TRACE, "/* gfp_ext_mul */ (");
      gfp_ext_print (a, modulus, OUTPUT_TRACE);
      outputf (OUTPUT_TRACE, ") * (");
      gfp_ext_print (b, modulus, OUTPUT_TRACE);
    }
  
  mpres_add (tmp[0], a->x, a->y, modulus);
  mpres_add (tmp[1], b->x, b->y, modulus);
  mpres_mul (tmp[1], tmp[0], tmp[1], modulus); /* t[1] = (a_0+a_1)*(b_0+b_1) = 
					    a_0*b_0 + a_0*b_1 + a_1*b_0 + 
					    a_1*b_1 */

  mpres_mul (r->x, a->x, b->x, modulus);    /* r_0 = a_0*b_0. We don't need a_0 
					    or b_0 any more now */
  mpres_sub (tmp[1], tmp[1], r->x, modulus);  /* t[1] = a_0*b_1 + a_1*b_0 + 
						a_1*b_1 */
  
  mpres_mul (tmp[0], a->y, b->y, modulus);   /* t[0] = a_1*b_1. We don't need 
					      a_1 or b_1 any more now */
  mpres_sub (r->y, tmp[1], tmp[0], modulus);  /* r_1 == a_0*b_1 + a_1*b_0 */
  
  mpres_mul (tmp[0], tmp[0], Delta, modulus); /* t[0] = a_1*b_1*Delta */
  mpres_add (r->x, r->x, tmp[0], modulus);   /* r_0 = a_0*b_0 + a_1*b_1*Delta */

  if (0 && test_verbose (OUTPUT_TRACE))
    {
      outputf (OUTPUT_TRACE, ") == ");
      gfp_ext_print (r, modulus, OUTPUT_TRACE);
      outputf (OUTPUT_TRACE, " /* PARI C */\n");
    }
}


/* Computes (a_0 + a_1 * sqrt(Delta))^2, where the norm 
   (a_0^2 - a_1^2*Delta) is assumed to be equal to 1. Hence 
   (a_0 + a_1 * sqrt(Delta))^2 = a_0^2 + 2*a_0*a_1*sqrt(Delta) + a_1^2*Delta
   and a_0^2 + a_1^2*Delta = a_0^2 + a_1^2*Delta + norm - 1 = 2*a_0^2 - 1.
   a_0 and r_0, as well as a_1 and r_1 may overlap */

static void
gfp_ext_sqr_norm1 (gfp_ext_t r, const gfp_ext_t a, mpmod_t modulus)
{
  if (pari)
    gmp_printf ("/* gfp_ext_sqr_norm1 */ (%Zd + %Zd * w)^2 %% N == ", a->x, a->y);
  
  mpres_mul (r->y, a->x, a->y, modulus);
  mpres_add (r->y, r->y, r->y, modulus);       /* r->y = 2*a_0*a_1 */
  
  mpres_sqr (r->x, a->x, modulus);
  mpres_add (r->x, r->x, r->x, modulus);
  mpres_sub_ui (r->x, r->x, 1UL, modulus);    /* r->x = 2*a_0^2 - 1 */

  if (pari)
    gmp_printf ("(%Zd + %Zd * w) %% N /* PARI C */\n", r->x, r->y);
}


/* Raise (a0 + a1*sqrt(Delta)) to the power e which is a 64-bit signed int.
   (a0 + a1*sqrt(Delta)) is assumed to have norm 1, i.e. 
   a0^2 - a1^2*Delta == 1. The result is (r0 * r1*sqrt(Delta)). 
   a0, a1, r0 and r1 must not overlap */

static void 
gfp_ext_pow_norm1_sl (gfp_ext_t r, const gfp_ext_t a, const int64_t e, 
                      const mpres_t Delta, 
                      mpmod_t modulus, const uint64_t tmplen, mpres_t *tmp)
{
  const int64_t abs_e = (e > 0) ? e : -e;
  uint64_t mask = (uint64_t)1 << 63;

  ASSERT (a != r);

  if (e == 0)
    {
      mpres_set_ui (r->x, 1UL, modulus);
      mpres_set_ui (r->y, 0UL, modulus);
      return;
    }

  /* If e < 0, we want 1/(a->x + a->y*sqrt(Delta)). By extending with 
     a->x - a->y*sqrt(Delta), we get 
     (a->x - a->y*sqrt(Delta)) / (a->x^2 - a->y^2 * Delta), but that denomiator
     is the norm which is known to be 1, so the result is 
     a->x - a->y*sqrt(Delta). */

  while ((abs_e & mask) == 0)
    mask >>= 1;

  mpres_set (r->x, a->x, modulus);
  mpres_set (r->y, a->y, modulus);

  while (mask > 1)
    {
      gfp_ext_sqr_norm1 (r, r, modulus);
      mask >>= 1;
      if (abs_e & mask)
	gfp_ext_mul (r, r, a, Delta, modulus, tmplen, tmp);
    }

  if (e < 0)
    mpres_neg (r->y, r->y, modulus);

  if (0 && test_verbose (OUTPUT_TRACE))
    {
      mpz_t t;
      mpz_init (t);
      mpres_get_z (t, Delta, modulus);
      outputf (OUTPUT_TRACE, "/* gfp_ext_pow_norm1_sl */ w = quadgen (4*%Zd); "
               "N = %Zd; /* PARI */\n", t, modulus->orig_modulus);
      mpz_clear (t);
      outputf (OUTPUT_TRACE, "/* gfp_ext_pow_norm1_sl */ (");
      gfp_ext_print (a, modulus, OUTPUT_TRACE);
      outputf (OUTPUT_TRACE, ")^(%" PRId64 ") == ", e);
      gfp_ext_print (r, modulus, OUTPUT_TRACE);
      outputf (OUTPUT_TRACE, " /* PARI C */\n");
    }
}


/* Same, but taking an mpz_t argument for the exponent */

static void 
gfp_ext_pow_norm1 (gfp_ext_t r, const gfp_ext_t a, mpz_t e, 
                   const mpres_t Delta, mpmod_t modulus, 
                   const uint64_t tmplen, mpres_t *tmp)
{
  mpz_t abs_e;
  unsigned long idx;

  ASSERT (a != r);

  if (mpz_sgn (e) == 0)
    {
      mpres_set_ui (r->x, 1UL, modulus);
      mpres_set_ui (r->y, 0UL, modulus);
      return;
    }

  mpz_init (abs_e);
  mpz_abs (abs_e, e);
  idx = mpz_sizeinbase (abs_e, 2) - 1; /* Thus mpz_tstbit (abs_e, idx) == 1 */
  ASSERT (mpz_tstbit (abs_e, idx) == 1);

  mpres_set (r->x, a->x, modulus);
  mpres_set (r->y, a->y, modulus);

  while (idx > 0UL)
    {
      gfp_ext_sqr_norm1 (r, r, modulus);
      idx--;
      if (mpz_tstbit (abs_e, idx))
	gfp_ext_mul (r, r, a, Delta, modulus, tmplen, tmp);
    }

  if (mpz_sgn (e) < 0)
    mpres_neg (r->y, r->y, modulus);

  mpz_clear (abs_e);

  if (test_verbose (OUTPUT_TRACE))
    {
      mpz_t t;
      mpz_init (t);
      mpres_get_z (t, Delta, modulus);
      outputf (OUTPUT_TRACE, "/* gfp_ext_pow_norm1 */ w = quadgen (4*%Zd); "
               "N = %Zd; /* PARI */\n", t, modulus->orig_modulus);
      mpz_clear (t);
      outputf (OUTPUT_TRACE, "/* gfp_ext_pow_norm1 */ (");
      gfp_ext_print (a, modulus, OUTPUT_TRACE);
      outputf (OUTPUT_TRACE, ")^(%Zd) == ", e);
      gfp_ext_print (r, modulus, OUTPUT_TRACE);
      outputf (OUTPUT_TRACE, " /* PARI C */\n");
    }
}


typedef struct {
  gfp_ext_t r1[2], r2[2];
  mpres_t v[2], v2, tmp;
  mpmod_t modulus;
  int imod2, want_y;
} pp1_sequence_g_state_t;

static void
pp1_sequence_g_get_x(void *statep, mpz_t r)
{
  pp1_sequence_g_state_t *state = statep;

  ASSERT (0 <= state->imod2 && state->imod2 < 2);

  /* r1[i] = r2[i-1] * v[i-2] - r2[i-2], with indices of r2 and i taken 
     modulo 2. */
  mpres_mul (state->tmp, state->r2[state->imod2 ^ 1]->x, state->v[state->imod2], state->modulus);
  mpres_sub (state->tmp, state->tmp, state->r2[state->imod2]->x, state->modulus);
  /* r2[i] = r2[i-1] * v[i-1] - r1[i-2] */
  mpres_mul (state->r2[state->imod2]->x, state->r2[state->imod2 ^ 1]->x, state->v[state->imod2 ^ 1], state->modulus);
  mpres_sub (state->r2[state->imod2]->x, state->r2[state->imod2]->x, state->r1[state->imod2]->x, state->modulus);
  mpres_set (state->r1[state->imod2]->x, state->tmp, state->modulus); /* FIXME, avoid this copy */
  mpres_get_z (r, state->tmp, state->modulus);

  /* If we don't call pp1_sequence_g_get_y(), we need to update state.v here */
  if (!state->want_y)
    {
      /* v[i] = v[i - 1] * V_2(a + 1/a) - v[i - 2] */
      mpres_mul (state->tmp, state->v[state->imod2 ^ 1], state->v2, state->modulus);
      mpres_sub (state->v[state->imod2], state->tmp, state->v[state->imod2], state->modulus);
      state->imod2 ^= 1;
    }
}

static void
pp1_sequence_g_get_y(void *statep, mpz_t r)
{
  pp1_sequence_g_state_t *state = statep;

  ASSERT (0 <= state->imod2 && state->imod2 < 2);

  mpres_mul (state->tmp, state->r2[state->imod2 ^ 1]->y, state->v[state->imod2], state->modulus);
  mpres_sub (state->tmp, state->tmp, state->r2[state->imod2]->y, state->modulus);
  mpres_mul (state->r2[state->imod2]->y, state->r2[state->imod2 ^ 1]->y, state->v[state->imod2 ^ 1], state->modulus);
  mpres_sub (state->r2[state->imod2]->y, state->r2[state->imod2]->y,     state->r1[state->imod2]->y, state->modulus);
  mpres_set (state->r1[state->imod2]->y, state->tmp, state->modulus);
  mpres_get_z (r, state->tmp, state->modulus);

  /* v[i] = v[i - 1] * V_2(a + 1/a) - v[i - 2] */
  mpres_mul (state->tmp, state->v[state->imod2 ^ 1], state->v2, state->modulus);
  mpres_sub (state->v[state->imod2], state->tmp, state->v[state->imod2], state->modulus);
  state->imod2 ^= 1;
}

/* Compute g_i = x_0^{M-i} * r^{(M-i)^2} for 0 <= i < l. 
   x_0 = b_1^{k_2 + m_1 * P}. r = b_1^P. */

static void
pp1_sequence_g (listz_t g_x, listz_t g_y, mpzspv_handle_t g_x_ntt, 
                mpzspv_handle_t g_y_ntt,
		const gfp_ext_t b1, const uint64_t P, 
		const mpres_t Delta, const uint64_t M_param, 
		const uint64_t l_param, const mpz_t m_1, 
		const int64_t k_2, const mpmod_t modulus_param)
{
  const int want_x = (g_x != NULL || g_x_ntt != NULL);
  const int want_y = (g_y != NULL || g_y_ntt != NULL);
  long timestart, realstart;

  /* We don't allow filling both NTT and mpz_t vectors */
  {
    const int use_ntt = (g_x_ntt != NULL || g_y_ntt != NULL) ? 1 : 0;
    const int use_mpz = (g_x != NULL || g_y != NULL) ? 1 : 0;
    ASSERT_ALWAYS ((use_ntt ^ use_mpz) == 1);
  }

  outputf (OUTPUT_VERBOSE, "Computing %s%s%s", 
	   (want_x) ? "g_x" : "", 
	   (want_x && want_y) ? " and " : "",
	   (want_y) ? "g_y" : "");
  timestart = cputime ();
  realstart = realtime ();

  /* When multi-threading, we adjust the parameters for each thread */
#ifdef _OPENMP
#pragma omp parallel if (l > 100)
#endif
  {
    pp1_sequence_g_state_t state;
    gfp_ext_t r, x0, tmp2;
    mpres_t tmp[3];
    mpz_t mt, mt1;
    uint64_t i, l, offset, M;
    int want_output = 1;
    const uint64_t tmplen = sizeof(tmp) / sizeof(mpres_t);
  
#ifdef _OPENMP
    want_output = (omp_get_thread_num() == 0);
    if (want_output)
      outputf (OUTPUT_VERBOSE, " using %d threads", omp_get_num_threads());
#endif
    get_chunk (&offset, &l, l_param);
    M = M_param - offset;
    mpmod_init_set (state.modulus, modulus_param);
    gfp_ext_init (r, state.modulus);
    gfp_ext_init (x0, state.modulus);
    gfp_ext_init (tmp2, state.modulus);
    mpres_init (state.v2, state.modulus);
    mpres_init (state.tmp, state.modulus);
    for (i = 0; i < 2UL; i++)
      {
	gfp_ext_init (state.r1[i], state.modulus);
	gfp_ext_init (state.r2[i], state.modulus);
	mpres_init (state.v[i], state.modulus);
      }
    for (i = 0; i < tmplen; i++)
      mpres_init (tmp[i], state.modulus);
    mpz_init (mt);
    mpz_init (mt1);
    
    if (want_output && test_verbose (OUTPUT_TRACE))
      {
	mpres_get_z (mt, Delta, state.modulus);
	outputf (OUTPUT_TRACE, 
		 "\n/* pp1_sequence_g */ w = quadgen (4*%Zd); ", mt);
	outputf (OUTPUT_TRACE, 
                 "P = %" PRIu64 "; M = %" PRIu64 "; k_2 = %" PRId64 ,
		 P, M, k_2);
	outputf (OUTPUT_TRACE, 
                 "; m_1 = %Zd; N = %Zd;/* PARI */\n", 
		 m_1, state.modulus->orig_modulus);
	
	outputf (OUTPUT_TRACE, "/* pp1_sequence_g */ b_1 = ");
	gfp_ext_print (b1, state.modulus, OUTPUT_TRACE);
	outputf (OUTPUT_TRACE, "; /* PARI */\n");
	outputf (OUTPUT_TRACE, 
		 "/* pp1_sequence_g */ r = b_1^(P/2); /* PARI */\n");
	outputf (OUTPUT_TRACE, "/* pp1_sequence_g */ "
		 "x_0 = b_1^(k_2 + m_1 * P); /* PARI */\n");
	outputf (OUTPUT_TRACE, 
		 "/* pp1_sequence_g */ addrec(x) = x + 1/x; /* PARI */\n");
      }
    
    /* Compute r */
    gfp_ext_pow_norm1_sl (r, b1, P/2, Delta, state.modulus, tmplen, tmp);
    if (want_output && test_verbose (OUTPUT_TRACE))
      {
	outputf (OUTPUT_TRACE, "/* pp1_sequence_g */ r == ");
	gfp_ext_print (r, state.modulus, OUTPUT_TRACE);
	outputf (OUTPUT_TRACE, " /* PARI C */\n");
      }
    
    /* Compute x0 = x_0 */
    mpz_set_uint64 (mt1, P);
    mpz_mul (mt, m_1, mt1); /* m_1 * P */
    mpz_set_int64 (mt1, k_2);
    mpz_add (mt, mt, mt1); /* mt = k_2 + m_1 * P */
    gfp_ext_pow_norm1 (x0, b1, mt, Delta, state.modulus, tmplen, tmp);
    if (want_output && test_verbose (OUTPUT_TRACE))
      {
	outputf (OUTPUT_TRACE, "/* pp1_sequence_g */ x_0 == ");
	gfp_ext_print (x0, state.modulus, OUTPUT_TRACE);
	outputf (OUTPUT_TRACE, " /* PARI C */\n");
      }
    
    
    /* Compute g[1] = r1[0] = x0^M * r^(M^2) = (x0 * r^M)^M. */
    gfp_ext_pow_norm1_sl (tmp2, r, M, Delta, state.modulus, tmplen, tmp); /* v[0,1] = r^M */
    gfp_ext_mul (tmp2, tmp2, x0, Delta, state.modulus, tmplen, tmp); /* v[0,1] = r^M * x_0 */
    gfp_ext_pow_norm1_sl (state.r1[0], tmp2, M, Delta, state.modulus, tmplen, tmp); /* r1[0] = (r^M * x_0)^M */
    if (g_x != NULL)
      mpres_get_z (g_x[offset], state.r1[0]->x, state.modulus);
    if (g_y != NULL)
      mpres_get_z (g_y[offset], state.r1[0]->y, state.modulus);
    if (g_x_ntt != NULL)
      {
	mpres_get_z (mt, state.r1[0]->x, state.modulus);
	mpzspv_fromto_mpzv (g_x_ntt, offset, 1, NULL, &mt, NULL, NULL);
      }
    if (g_y_ntt != NULL)
      {
	mpres_get_z (mt, state.r1[0]->y, state.modulus);
	mpzspv_fromto_mpzv (g_y_ntt, offset, 1, NULL, &mt, NULL, NULL);
      }
    
    
    /* Compute g[1] = r1[1] = x0^(M-1) * r^((M-1)^2) = (x0 * r^(M-1))^(M-1). 
       We use v[0,1] as temporary storage. FIXME: simplify, reusing g_0 */
    gfp_ext_pow_norm1_sl (tmp2, r, M - 1, Delta, state.modulus, tmplen, tmp);
    gfp_ext_mul (tmp2, tmp2, x0, Delta, state.modulus, tmplen, tmp);
    gfp_ext_pow_norm1_sl (state.r1[1], tmp2, M - 1, Delta, state.modulus, tmplen, tmp);
    if (g_x != NULL)
      mpres_get_z (g_x[offset + 1], state.r1[1]->x, state.modulus);
    if (g_y != NULL)
      mpres_get_z (g_y[offset + 1], state.r1[1]->y, state.modulus);
    if (g_x_ntt != NULL)
      {
	mpres_get_z (mt, state.r1[1]->x, state.modulus);
	mpzspv_fromto_mpzv (g_x_ntt, offset + 1, 1, NULL, &mt, NULL, NULL);
      }
    if (g_y_ntt != NULL)
      {
	mpres_get_z (mt, state.r1[1]->y, state.modulus);
	mpzspv_fromto_mpzv (g_y_ntt, offset + 1, 1, NULL, &mt, NULL, NULL);
      }
    
    
    /* x0 := $x_0 * r^{2M - 3}$ */
    /* We don't need x0 after this so we overwrite it. We use v[0,1] as 
       temp storage for $r^{2M - 3}$. */
    gfp_ext_pow_norm1_sl (tmp2, r, 2UL*M - 3UL, Delta, state.modulus, tmplen, tmp);
    gfp_ext_mul (x0, x0, tmp2, Delta, state.modulus, tmplen, tmp);
    
    /* Compute r2[0] = r1[0] * r^2 and r2[1] = r1[1] * r^2. */
    /* We only need $r^2$ from here on, so we set r = $r^2$ */
    gfp_ext_sqr_norm1 (r, r, state.modulus);  
    gfp_ext_mul (state.r2[0], state.r1[0], r, Delta, state.modulus, tmplen, tmp);
    gfp_ext_mul (state.r2[1], state.r1[1], r, Delta, state.modulus, tmplen, tmp);
    
    /* v[1] := $x_0 * r^{2*M - 3} + 1/(x_0 * r^{2M - 3}) */
    mpres_add (state.v[1], x0->x, x0->x, state.modulus);
    /* x0 := x0 * r = $x_0 * r^{2M - 1}$ */
    gfp_ext_mul (x0, x0, r, Delta, state.modulus, tmplen, tmp);
    /* v[0] := $x_0 * r^{2M - 1} + 1/(x_0 * r^{2M - 1}) */
    mpres_add (state.v[0], x0->x, x0->x, state.modulus);
    gfp_ext_clear (x0, state.modulus);
    
    /* v2 = V_2 (r + 1/r) = r^2 + 1/r^2 */
    mpres_add (state.v2, r->x, r->x, state.modulus);
    
    gfp_ext_clear (r, state.modulus);

    state.want_y = want_y;
    state.imod2 = 0;
    if (l > 2 && g_x_ntt != NULL && g_y_ntt == NULL)
      {
        mpzspv_fromto_mpzv (g_x_ntt, offset + 2, l - 2, &pp1_sequence_g_get_x, &state, NULL, NULL);
      }
    else if (l > 2 && g_x_ntt == NULL && g_y_ntt != NULL)
      {
        mpzspv_fromto_mpzv (g_y_ntt, offset + 2, l - 2, &pp1_sequence_g_get_y, &state, NULL, NULL);
      }
    else
      {
        if ((g_x_ntt != NULL && g_x_ntt->storage == 1) || (g_y_ntt != NULL && g_y_ntt->storage == 1))
          {
            outputf (OUTPUT_ERROR, "Warning: pp1_sequence_g() uses slow code\n");
          }
        for (i = 2; i < l; i++)
          {
            if (g_x != NULL)
              pp1_sequence_g_get_x (&state, g_x[offset + i]);
            if (g_y != NULL)
              pp1_sequence_g_get_y (&state, g_y[offset + i]);

            /* FIXME check if repeated calls to mpzspv_fromto_mpzv() slow 
               things down. Maybe do blocks with temp buffer */
            if (g_x_ntt != NULL)
              mpzspv_fromto_mpzv (g_x_ntt, offset + i, 1, &pp1_sequence_g_get_x, &state, NULL, NULL);
            if (g_y_ntt != NULL)
              mpzspv_fromto_mpzv (g_y_ntt, offset + i, 1, &pp1_sequence_g_get_y, &state, NULL, NULL);
            
            if (want_output && test_verbose (OUTPUT_TRACE))
              {
                mpres_get_z (mt, state.v[state.imod2], state.modulus);
                outputf (OUTPUT_TRACE, "/* pp1_sequence_g */ "
                         "addrec(x_0 * r^(2*(M-%lu) - 1)) == %Zd /* PARI C */\n", 
                         i, mt);
              }
          }
      }
    
    mpres_clear (state.v2, state.modulus);
    mpres_clear (state.tmp, state.modulus);
    gfp_ext_clear (tmp2, state.modulus);
    for (i = 0; i < 2; i++)
      {
	gfp_ext_clear (state.r1[i], state.modulus);
	gfp_ext_clear (state.r2[i], state.modulus);
	mpres_clear (state.v[i], state.modulus);
      }
    for (i = 0; i < tmplen; i++)
      mpres_clear (tmp[i], state.modulus);
    mpz_clear (mt);
    mpz_clear (mt1);
    mpmod_clear (state.modulus);
  }
  
  print_elapsed_time (OUTPUT_VERBOSE, timestart, realstart);

  if (g_x != NULL && g_y != NULL && test_verbose(OUTPUT_TRACE))
    {
      uint64_t i;
      for (i = 0; i < l_param; i++)
	{
	  outputf (OUTPUT_TRACE, "/* pp1_sequence_g */ g_%" PRIu64 " = "
		   "x_0^(M-%" PRIu64 ") * r^((M-%" PRIu64 ")^2); /* PARI */", 
		   i, i, i);
	  outputf (OUTPUT_TRACE, "/* pp1_sequence_g */ g_%" PRIu64 " == "
		   "%Zd + %Zd*w /* PARI C */\n", 
		   i, g_x[i], g_y[i]);
	}
    }
}


typedef struct {
  gfp_ext_t s[3], s2[2];
  mpres_t v[2], V2, tmp;
  mpmod_t modulus;
  listz_handle_t f;
  file_word_t *buf;
  uint64_t i, offset;
  int want_y;
} pp1_sequence_h_state_t;

static void
pp1_sequence_h_get_x (void *statep, mpz_t r)
{
  pp1_sequence_h_state_t *state = statep;
  const int imod2 = state->i % 2;
  const int imod3 = state->i % 3;

  /* r[i] = r2[i-1] * v[i-2] - r2[i-2], with indices of r2 and i 
     taken modulo 2 */
  mpres_mul (state->s[imod3]->x, state->s2[imod2 ^ 1]->x, state->v[imod2], state->modulus);
  mpres_sub (state->s[imod3]->x, state->s[imod3]->x, state->s2[imod2]->x, state->modulus);
  
  /* r2[i] = r2[i-1] * v[i-1] - r[i-2] */
  mpres_mul (state->s2[imod2]->x, state->s2[imod2 ^ 1]->x, state->v[ imod2 ^ 1], state->modulus);
  mpres_sub (state->s2[imod2]->x, state->s2[imod2]->x, state->s[(imod3 + 1) % 3]->x, state->modulus);

  listz_handle_get (state->f,  state->tmp, state->buf, state->i + state->offset);
  mpres_mul_z_to_z (r, state->s[imod3]->x, state->tmp, state->modulus);
  outputf (OUTPUT_DEVVERBOSE, "h_x,%lu = %Zd\n", state->i, r);

  if (!state->want_y)
    {
      /* v[i] = v[i - 1] * V_2(a + 1/a) - v[i - 2] */
      mpres_mul (state->tmp, state->v[1 - imod2], state->V2, state->modulus);
      mpres_sub (state->v[imod2], state->tmp, state->v[imod2], state->modulus);

      state->i++;
    }
}

static void
pp1_sequence_h_get_y (void *statep, mpz_t r)
{
  pp1_sequence_h_state_t *state = statep;
  const int imod2 = state->i % 2;
  const int imod3 = state->i % 3;

  mpres_mul (state->s[imod3]->y, state->s2[1 - imod2]->y, state->v[imod2], state->modulus);
  mpres_sub (state->s[imod3]->y, state->s[imod3]->y, state->s2[imod2]->y, state->modulus);
  mpres_mul (state->s2[imod2]->y, state->s2[1 - imod2]->y, state->v[1 - imod2], state->modulus);
  mpres_sub (state->s2[imod2]->y, state->s2[imod2]->y, state->s[(imod3 + 1) % 3]->y, state->modulus);

  listz_handle_get (state->f,  state->tmp, state->buf, state->i + state->offset);
  mpres_mul_z_to_z (r, state->s[imod3]->y, state->tmp, state->modulus);
  outputf (OUTPUT_DEVVERBOSE, "h_y,%lu = %Zd\n", state->i, r);

  /* v[i] = v[i - 1] * V_2(a + 1/a) - v[i - 2] */
  mpres_mul (state->tmp, state->v[1 - imod2], state->V2, state->modulus);
  mpres_sub (state->v[imod2], state->tmp, state->v[imod2], state->modulus);

  state->i++;
}

/* Compute r[i] = b1^(-P*(k+i)^2) * f_i for i = 0, 1, ..., l-1, where "b1" is 
   an element of norm 1 in the quadratic extension ring */

static void
pp1_sequence_h (listz_t h_x, listz_t h_y, mpzspv_handle_t h_x_ntt, mpzspv_handle_t h_y_ntt,
		const listz_handle_t f, const gfp_ext_t b1, 
		const int64_t k_param, const uint64_t l_param, 
		const uint64_t P, const mpres_t Delta, 
		mpmod_t modulus_param)
{
  uint64_t i;
  long timestart, realstart;
  const int want_x = (h_x != NULL || h_x_ntt != NULL);
  const int want_y = (h_y != NULL || h_y_ntt != NULL);

  if (l_param == 0UL)
    return;

  ASSERT (f->storage == 0 || f->data.mem != h_x);
  ASSERT (f->storage == 0 || f->data.mem != h_y);

  outputf (OUTPUT_VERBOSE, "Computing %s%s%s", 
	   (want_x) ? "h_x" : "", 
	   (want_x && want_y) ? " and " : "",
	   (want_y) ? "h_y" : "");
  timestart = cputime ();
  realstart = realtime ();

  if (test_verbose (OUTPUT_TRACE))
    {
      mpz_t t;
      mpz_init (t);
      mpres_get_z (t, Delta, modulus_param);
      outputf (OUTPUT_TRACE, "\n/* pp1_sequence_h */ N = %Zd; "
	       "Delta = %Zd; ", modulus_param->orig_modulus, t);
      outputf (OUTPUT_TRACE, "k = %ld; P = %lu; /* PARI */\n", k_param, P);
      outputf (OUTPUT_TRACE, "/* pp1_sequence_h */ b_1 = ");
      gfp_ext_print (b1, modulus_param, OUTPUT_TRACE);
      outputf (OUTPUT_TRACE, "; r = b_1^(P/2); rn = b_1^(-P/2); /* PARI */\n");
      for (i = 0; i < l_param; i++)
        {
          listz_handle_get2 (f, t, i);
          outputf (OUTPUT_TRACE, 
                   "/* pp1_sequence_h */ f_%" PRIu64 " = %Zd; /* PARI */\n", 
                   i, t);
         }
      mpz_clear (t);
    }

#ifdef _OPENMP
#pragma omp parallel if (l_param > 100) private(i)
#endif
  {
    const uint64_t tmplen = 2;
    pp1_sequence_h_state_t state;
    gfp_ext_t rn, tmp2;
    mpres_t tmp[2];
    mpz_t mt;
    uint64_t l, offset;
    int64_t k = k_param;

    /* When multi-threading, we adjust the parameters for each thread */
    get_chunk (&offset, &l, l_param);
#ifdef _OPENMP
    if (omp_get_thread_num() == 0)
      outputf (OUTPUT_VERBOSE, " using %d threads", omp_get_num_threads());
    outputf (OUTPUT_TRACE, "\n");
#endif

    /* Each thread computes r[i + offset] = b1^(-P*(k+i+offset)^2) * f_i 
       for i = 0, 1, ..., l-1, where l is the adjusted length of each thread */

    k += offset;

    mpz_init (mt);
    /* Make thread-local copy of modulus */
    mpmod_init_set (state.modulus, modulus_param);

    /* Init the local mpres_t variables */
    for (i = 0; i < 2; i++)
      {
	gfp_ext_init (state.s[i], state.modulus);
	gfp_ext_init (state.s2[i], state.modulus);
	mpres_init (state.v[i], state.modulus);
      }
    gfp_ext_init (state.s[2], state.modulus);
    mpres_init (state.V2, state.modulus);
    mpres_init (state.tmp, state.modulus);
    gfp_ext_init (rn, state.modulus);
    gfp_ext_init (tmp2, state.modulus);
    for (i = 0; i < tmplen; i++)
      mpres_init (tmp[i], state.modulus);

    /* Compute rn = b_1^{-P}. It has the same value for all threads,
       but we make thread local copies anyway. */
    gfp_ext_pow_norm1_sl (rn, b1, P/2, Delta, state.modulus, tmplen, tmp);
    mpres_neg (rn->y, rn->y, state.modulus);
    
    /* Compute s[0] = rn^(k^2) = r^(-k^2). We do it by two exponentiations 
       by k */
    gfp_ext_pow_norm1_sl (tmp2, rn, k, Delta, state.modulus, tmplen, tmp);
    gfp_ext_pow_norm1_sl (state.s[0], tmp2, k, Delta, state.modulus, tmplen, tmp);
    if (test_verbose (OUTPUT_TRACE))
      {
#ifdef _OPENMP
#pragma omp critical
#endif
	{
	  outputf (OUTPUT_TRACE, 
	  	"/* pp1_sequence_h */ rn^(%" PRId64 "^2) == ", k);
	  gfp_ext_print (state.s[0], state.modulus, OUTPUT_TRACE);
	  outputf (OUTPUT_TRACE, " /* PARI C */\n");
	}
      }
    
    /* Compute s[1] = r^(-(k+1)^2) = r^(-(k^2 + 2k + 1))*/
    if (l > 1)
      {
	/* v[0] + v[1]*sqrt(Delta) still contains rn^k */
	gfp_ext_sqr_norm1 (state.s[1], tmp2, state.modulus);
	/* Now s[1] = r^(-2k) */
	gfp_ext_mul (state.s[1], state.s[1], state.s[0], Delta, state.modulus, tmplen, tmp);
	/* Now s[1] = r^(-(k^2 + 2k)) */
	gfp_ext_mul (state.s[1], state.s[1], rn, Delta, state.modulus, tmplen, tmp);
	/* Now s[1] = r^(-(k^2 + 2k + 1)) = r^(-(k+1)^2) */
	if (test_verbose (OUTPUT_TRACE))
	  {
#ifdef _OPENMP
#pragma omp critical
#endif
	    {
	      outputf (OUTPUT_TRACE, 
		  	"/* pp1_sequence_h */ rn^(%" PRId64 "^2) == ", k + 1);
	      gfp_ext_print (state.s[1], state.modulus, OUTPUT_TRACE);
	      outputf (OUTPUT_TRACE, " /* PARI C */\n");
	    }
	  }
      }
    
    /* Compute s2[0] = r^(k^2+2) = r^(k^2) * r^2 */
    gfp_ext_sqr_norm1 (tmp2, rn, state.modulus);
    gfp_ext_mul (state.s2[0], state.s[0], tmp2, Delta, state.modulus, tmplen, tmp);
    if (test_verbose (OUTPUT_TRACE))
      {
#ifdef _OPENMP
#pragma omp critical
#endif
	{
	  outputf (OUTPUT_TRACE, "/* pp1_sequence_h */ rn^(%ld^2+2) == ", k);
	  gfp_ext_print (state.s2[0], state.modulus, OUTPUT_TRACE);
	  outputf (OUTPUT_TRACE, " /* PARI C */\n");
	}
      }
    
    /* Compute a^((k+1)^2+2) = a^((k+1)^2) * a^2 */
    gfp_ext_mul (state.s2[1], state.s[1], tmp2, Delta, state.modulus, tmplen, tmp);
    if (test_verbose (OUTPUT_TRACE))
      {
#ifdef _OPENMP
#pragma omp critical
#endif
	{
	  outputf (OUTPUT_TRACE, "/* pp1_sequence_h */ rn^(%ld^2+2) == ", 
		   k + 1);
	  gfp_ext_print (state.s2[1], state.modulus, OUTPUT_TRACE);
	  outputf (OUTPUT_TRACE, " /* PARI C */\n");
	}
      }
    
    /* Compute V_2(r + 1/r). Since 1/r = rn.x - rn.y, we have r+1/r = 2*rn.x.
       V_2(x) = x^2 - 2, so we want 4*rn.x^2 - 2. */
    mpres_add (state.V2, rn->x, rn->x, state.modulus); /* V2 = r + 1/r  = 2*rn.x */
    V (state.v[0], state.V2, 2 * k + 1, state.modulus);  /* v[0] = V_{2k+1} (r + 1/r) */
    V (state.v[1], state.V2, 2 * k + 3, state.modulus);  /* v[1] = V_{2k+3} (r + 1/r) */
    mpres_sqr (state.V2, state.V2, state.modulus); /* V2 = 4*a_x^2 */
    mpres_sub_ui (state.V2, state.V2, 2UL, state.modulus); /* V2 = 4*a_x^2 - 2 */
    if (test_verbose (OUTPUT_TRACE))
      {
#ifdef _OPENMP
#pragma omp critical
#endif
	{
	  mpres_get_z (mt, state.V2, state.modulus);
	  outputf (OUTPUT_TRACE, "/* pp1_sequence_h */ r^2 + 1/r^2 == %Zd "
		   "/* PARI C */\n", mt);
	  mpres_get_z (mt, state.v[0], state.modulus);
	  outputf (OUTPUT_TRACE, "/* pp1_sequence_h */ r^(2*%ld+1) + "
		   "1/r^(2*%ld+1) == %Zd /* PARI C */\n", 
		   (long)k, (long)k, mt);
	  mpres_get_z (mt, tmp2->y, state.modulus);
	  outputf (OUTPUT_TRACE, "/* pp1_sequence_h */ r^(2*%ld+3) + "
		   "1/r^(2*%ld+3) == %Zd /* PARI C */\n", 
		   (long)k, (long)k, mt);
	}
      }
    
    for (i = 0; i < 2UL && i < l; i++)
      {
	/* Multiply the 2nd coordinate by Delta, so that after the polynomial
	   multipoint evaluation we get x1 + Delta*x2 */
	mpres_mul (state.s[i]->y, state.s[i]->y, Delta, state.modulus);
	mpres_mul (state.s2[i]->y, state.s2[i]->y, Delta, state.modulus);

        listz_handle_get2 (f, tmp[0], i + offset);

	if (h_x != NULL)
          mpres_mul_z_to_z (h_x[i + offset], state.s[i]->x, tmp[0], state.modulus);
	if (h_y != NULL)
          mpres_mul_z_to_z (h_y[i + offset], state.s[i]->y, tmp[0], state.modulus);
	if (h_x_ntt != NULL)
	  {
	    mpres_mul_z_to_z (mt, state.s[i]->x, tmp[0], state.modulus);
            mpzspv_fromto_mpzv (h_x_ntt, i + offset, 1, NULL, &mt, NULL, NULL);
	  }
	if (h_y_ntt != NULL)
	  {
	    mpres_mul_z_to_z (mt, state.s[i]->y, tmp[0], state.modulus);
            mpzspv_fromto_mpzv (h_y_ntt, i + offset, 1, NULL, &mt, NULL, NULL);
	  }
      }
    
    /* Compute the remaining r^((k+i)^2) values according to Peter's 
       recurrence */
    
    state.want_y = want_y;
    state.f = f;
    state.buf = (file_word_t *) malloc (f->words * sizeof (file_word_t));
    ASSERT_ALWAYS (state.buf != NULL);
    state.i = 2;
    state.offset = offset;

    if (h_x_ntt != NULL && h_y_ntt == NULL && l > 2)
      {
        /* These two cases are for use with -treefile where we handle coordinates separately */
        mpzspv_fromto_mpzv (h_x_ntt, offset + 2, l - 2, &pp1_sequence_h_get_x, &state, NULL, NULL);
      }
    else if (h_x_ntt == NULL && h_y_ntt != NULL && l > 2)
      {
        mpzspv_fromto_mpzv (h_y_ntt, offset + 2, l - 2, &pp1_sequence_h_get_y, &state, NULL, NULL);
      }
    else for (i = 2; i < l; i++)
      {
	if (h_x != NULL || h_x_ntt != NULL)
	  {
	    ASSERT (state.i == i);
	    if (h_x != NULL)
	      {
	        pp1_sequence_h_get_x (&state, h_x[i + offset]); 
              }
	    if (h_x_ntt != NULL)
	      {
		mpzspv_fromto_mpzv (h_x_ntt, i + offset, 1, &pp1_sequence_h_get_x, &state, NULL, NULL);
	      }
	  }
	
	if (h_y != NULL || h_y_ntt != NULL)
	  {
	    ASSERT (state.i == i);
	    /* Same for y coordinate */
	    if (h_y != NULL)
	      {
	        pp1_sequence_h_get_y (&state, h_y[i + offset]);
              }
	    if (h_y_ntt != NULL)
	      {
		mpzspv_fromto_mpzv (h_y_ntt, i + offset, 1, &pp1_sequence_h_get_y, &state, NULL, NULL);
	      }
	  }
      }
    
    free (state.buf);
    /* Clear the local mpres_t variables */
    for (i = 0; i < 2; i++)
      {
	gfp_ext_clear (state.s[i], state.modulus);
	gfp_ext_clear (state.s2[i], state.modulus);
	mpres_clear (state.v[i], state.modulus);
      }
    gfp_ext_clear (state.s[2], state.modulus);
    mpres_clear (state.V2, state.modulus);
    mpres_clear (state.tmp, state.modulus);
    gfp_ext_clear (rn, state.modulus);
    gfp_ext_clear (tmp2, state.modulus);
    for (i = 0; i < tmplen; i++)
      mpres_clear (tmp[i], state.modulus);

    /* Clear the thread-local copy of modulus */
    mpmod_clear (state.modulus);

    mpz_clear (mt);
  }

  print_elapsed_time (OUTPUT_VERBOSE, timestart, realstart);

  if (h_x != NULL && h_y != NULL && test_verbose (OUTPUT_TRACE))
    {
      for (i = 0; i < l_param; i++)
	gmp_printf ("/* pp1_sequence_h */ (rn^((k+%" PRIu64 ")^2) * f_%" 
	            PRIu64 ") == "
		    "(%Zd + Mod(%Zd / Delta, N) * w) /* PARI C */\n", 
		    i, i, h_x[i], h_y[i]);
    }
}


int 
pp1fs2 (mpz_t f, const mpres_t X, mpmod_t modulus, 
	const faststage2_param_t *params)
{
  uint64_t nr;
  uint64_t i, l, lenF, lenH, lenG, lenR, tmplen;
  set_list_t S_1; /* This is stored as a set of sets (arithmetic 
                       progressions of prime length */
  int64_t *s2_sumset; /* set of sums of S_2 */
  uint64_t s2_sumset_size;
  listz_handle_t F;   /* Polynomial F has roots X^{k_1} for k_1 \in S_1, so has 
		  degree s_1. It is symmetric, so has only s_1 / 2 + 1 
		  distinct coefficients. The sequence h_j will be stored in 
		  the same memory and won't be a monic polynomial, so the 
		  leading 1 monomial of F will be stored explicitly. Hence we 
		  need s_1 / 2 + 1 entries. */

  listz_t g_x, g_y, fh_x, fh_y, h_x, h_y, tmp, R_x, R_y; 
  const unsigned long tmpreslen = 2UL;
  gfp_ext_t b1;
  mpres_t Delta, tmpres[2];
  mpz_t mt;   /* All-purpose temp mpz_t */
  int youpi = ECM_NO_FACTOR_FOUND;
  long timetotalstart, realtotalstart, timestart;

  timetotalstart = cputime ();
  realtotalstart = realtime ();

  ASSERT_ALWAYS (eulerphi (params->P) == params->s_1 * params->s_2);
  ASSERT_ALWAYS (params->s_1 < params->l);
  nr = params->l - params->s_1; /* Number of points we evaluate */

  sets_init(&S_1);

  if (make_S_1_S_2 (&S_1, &s2_sumset, 
		  &s2_sumset_size, params) == ECM_ERROR)
      return ECM_ERROR;

  /* Allocate all the memory we'll need */
  /* Allocate the correct amount of space for each mpz_t or the 
     reallocations will up to double the time for stage 2! */
  mpz_init (mt);
  gfp_ext_init (b1, modulus);
  mpres_init (Delta, modulus);
  for (i = 0; i < tmpreslen; i++)
      mpres_init (tmpres[i], modulus);
  lenF = params->s_1 / 2 + 1 + 1; /* Another +1 because poly_from_sets_V stores
				     the leading 1 monomial for each factor */
  lenH = params->s_1 + 1;
  lenG = params->l;
  lenR = nr;
  F = listz_handle_init2 (NULL, lenF, modulus->orig_modulus);
  fh_x = init_list2 (lenF, (unsigned int) abs (modulus->bits));
  fh_y = init_list2 (lenF, (unsigned int) abs (modulus->bits));
  h_x = malloc (lenH * sizeof (mpz_t));
  h_y = malloc (lenH * sizeof (mpz_t));
  if (h_x == NULL || h_y == NULL)
    {
      fprintf (stderr, "Cannot allocate memory in pp1fs2\n");
      exit (1);
    }
  g_x = init_list2 (lenG, (unsigned int) abs (modulus->bits));
  g_y = init_list2 (lenG, (unsigned int) abs (modulus->bits));
  R_x = init_list2 (lenR, (unsigned int) abs (modulus->bits));
  R_y = init_list2 (lenR, (unsigned int) abs (modulus->bits));
  tmplen = 3UL * params->l + list_mul_mem (params->l / 2) + 20;
  outputf (OUTPUT_DEVVERBOSE, "tmplen = %" PRIu64 "\n", tmplen);
  if (TMulGen_space (params->l - 1, params->s_1, lenR) + 12 > tmplen)
    {
      tmplen = TMulGen_space (params->l - 1, params->s_1 - 1, lenR) + 12;
      /* FIXME: It appears TMulGen_space() returns a too small value! */
      outputf (OUTPUT_DEVVERBOSE, "With TMulGen_space, tmplen = %lu\n", 
	       tmplen);
    }

  tmp = init_list2 (tmplen, (unsigned int) abs (modulus->bits));

  if (test_verbose (OUTPUT_TRACE))
    {
      mpres_get_z (mt, X, modulus); /* mpz_t copy of X for printing */
      outputf (OUTPUT_TRACE, 
	       "N = %Zd; X = Mod(%Zd, N); /* PARI */\n", 
	       modulus->orig_modulus, mt);
    }

  /* Compute the polynomial f(x) = \prod_{k_1 in S_1} (x - X^{2 k_1}) */
  outputf (OUTPUT_VERBOSE, "Computing F from factored S_1");
  
  timestart = cputime ();
  i = poly_from_sets_V (F, X, &S_1, tmp, tmplen, modulus, NULL);
  ASSERT_ALWAYS(2 * i == params->s_1);
  ASSERT(mpz_cmp_ui (F->data.mem[i], 1UL) == 0);
  sets_free(&S_1);
  
  outputf (OUTPUT_VERBOSE, " took %lums\n", cputime () - timestart);
  if (test_verbose (OUTPUT_TRACE))
    {
      for (i = 0; i < params->s_1 / 2 + 1; i++)
	outputf (OUTPUT_TRACE, "f_%" PRIu64 " = %Zd; /* PARI */\n", i, F->data.mem[i]);
      outputf (OUTPUT_TRACE, "f(x) = f_0");
      for (i = 1; i < params->s_1 / 2 + 1; i++)
	outputf (OUTPUT_TRACE, "+ f_%" PRIu64 " * (x^%" PRIu64
                        " + x^(-%" PRIu64 "))", i, i, i);
      outputf (OUTPUT_TRACE, "/* PARI */ \n");
    }

  /* Compute Delta and b1_x + b1_y * sqrt(Delta) = X) */
  mpres_sqr (Delta, X, modulus);
  mpres_sub_ui (Delta, Delta, 4UL, modulus);
  mpres_div_2exp (b1->x, X, 1, modulus);
  mpres_set_ui (b1->y, 1UL, modulus);
  mpres_div_2exp (b1->y, b1->y, 1, modulus);
  if (test_verbose (OUTPUT_TRACE))
    {
      mpres_get_z (mt, Delta, modulus);
      outputf (OUTPUT_TRACE, 
	       "Delta = Mod(%Zd, N); w = quadgen (4*lift(Delta)); b_1 = ", mt);
      gfp_ext_print (b1, modulus, OUTPUT_TRACE);
      outputf (OUTPUT_TRACE, "; /* PARI */\n");
      outputf (OUTPUT_TRACE, "X == b_1 + 1/b_1 /* PARI C */\n");
    }

  /* Compute the h sequence h_j = b1^(P*-j^2) * f_j for 0 <= j <= s_1 */
  pp1_sequence_h (fh_x, fh_y, NULL, NULL, F, b1, 0L, 
		  params->s_1 / 2 + 1, params->P, Delta, modulus);
  /* We don't need F(x) any more */
  listz_handle_clear (F);
  F = NULL;

  /* Make a symmetric copy of fh in h. */
  for (i = 0; i < params->s_1 / 2 + 1; i++)
    {
      *(h_x[i]) = *(fh_x[params->s_1 / 2 - i]); /* Clone the mpz_t */
      *(h_y[i]) = *(fh_y[params->s_1 / 2 - i]);
    }
  for (i = 0; i < params->s_1 / 2; i++)
    {
      *(h_x[i + params->s_1 / 2 + 1]) = *(fh_x[i + 1]);
      *(h_y[i + params->s_1 / 2 + 1]) = *(fh_y[i + 1]);
    }
  if (test_verbose (OUTPUT_TRACE))
    {
      for (i = 0; i < params->s_1 + 1; i++)
	outputf (OUTPUT_VERBOSE, "h_%lu = %Zd + %Zd * w; /* PARI */\n", 
		 i, h_x[i], h_y[i]);
    }
  
  for (l = 0; l < params->s_2 && youpi == ECM_NO_FACTOR_FOUND; l++)
    {
      const uint64_t M = params->l - 1 - params->s_1 / 2;
      outputf (OUTPUT_VERBOSE, "Multi-point evaluation %" PRIu64
                        " of %" PRIu64 ":\n", l + 1, params->s_2);
      pp1_sequence_g (g_x, g_y, NULL, NULL, b1, params->P, 
		      Delta, M, params->l, params->m_1, s2_sumset[l], 
		      modulus);
      
      /* Do the two convolution products */
      outputf (OUTPUT_VERBOSE, "TMulGen of g_x and h_x");
      timestart = cputime ();
      if (TMulGen (R_x, nr - 1, h_x, params->s_1, g_x, params->l - 1, tmp,
		   modulus->orig_modulus) < 0)
	{
	  outputf (OUTPUT_ERROR, "TMulGen returned error code (probably out "
		   "of memory)\n");
	  youpi = ECM_ERROR;
	  break;
	}
      outputf (OUTPUT_VERBOSE, " took %lums\n", cputime () - timestart);
      outputf (OUTPUT_VERBOSE, "TMulGen of g_y and h_y");
      timestart = cputime ();
      if (TMulGen (R_y, nr - 1, h_y, params->s_1, g_y, params->l - 1, tmp,
		   modulus->orig_modulus) < 0)
	{
	  outputf (OUTPUT_ERROR, "TMulGen returned error code (probably out "
		   "of memory)\n");
	  youpi = ECM_ERROR;
	  break;
	}
      outputf (OUTPUT_VERBOSE, " took %lums\n", cputime () - timestart);
      
      timestart = cputime ();
      mpres_set_ui (tmpres[1], 1UL, modulus); /* Accumulate product in 
						 tmpres[1] */
      for (i = 0; i < nr; i++)
      {
	  mpz_add (mt, R_x[i], R_y[i]);
	  mpres_set_z_for_gcd (tmpres[0], mt, modulus);
#define TEST_ZERO_RESULT
#ifdef TEST_ZERO_RESULT
	  if (mpres_is_zero (tmpres[0], modulus))
	      outputf (OUTPUT_VERBOSE, "R_[%lu] = 0\n", i);
#endif
	  mpres_mul (tmpres[1], tmpres[1], tmpres[0], modulus); 
      }
      mpz_set_uint64 (mt, nr);
      mpres_set_z_for_gcd_fix (tmpres[1], tmpres[1], mt, modulus);
      outputf (OUTPUT_VERBOSE, "Computing product of F(g_i)^(1) took %lums\n", 
	       cputime () - timestart);
      if (test_verbose(OUTPUT_RESVERBOSE))
      {
	  mpres_get_z (mt, tmpres[1], modulus);
	  outputf (OUTPUT_RESVERBOSE, "Product of R[i] = %Zd\n", mt);
      }
      
      mpres_gcd (mt, tmpres[1], modulus);
      if (mpz_cmp_ui (mt, 1UL) > 0)
	{
	  mpz_set (f, mt);
	  youpi = ECM_FACTOR_FOUND_STEP2;
	  break;
	}
    }

  mpz_clear (mt);
  gfp_ext_clear (b1, modulus);
  mpres_clear (Delta, modulus);
  for (i = 0; i < tmpreslen; i++)
      mpres_clear (tmpres[i], modulus);
  clear_list (fh_x, lenF);
  clear_list (fh_y, lenF);
  free (h_x);
  free (h_y);
  clear_list (g_x, lenG);
  clear_list (g_y, lenG);
  clear_list (R_x, lenR);
  clear_list (R_y, lenR);
  clear_list (tmp, tmplen);
  free (s2_sumset);
 
  outputf (OUTPUT_NORMAL, "Step 2");
  /* In normal output mode, print only cpu time as we always have.
     In verbose mode, print real time as well if we used multi-threading */
  if (test_verbose (OUTPUT_VERBOSE))
    print_elapsed_time (OUTPUT_NORMAL, timetotalstart, realtotalstart);
  else
    print_elapsed_time (OUTPUT_NORMAL, timetotalstart, 0L);

  return youpi;
}


int 
pp1fs2_ntt (mpz_t f, const mpres_t X, mpmod_t modulus,
	    const faststage2_param_t *params, const int twopass_param)
{
  uint64_t nr;
  uint64_t l;
  set_list_t S_1; /* This is stored as a set of sets (arithmetic 
                       progressions of prime length */
  int64_t *s2_sumset; /* set of sums of S_2 */
  uint64_t s2_sumset_size;
  listz_handle_t F;   /* Polynomial F has roots X^{k_1} for k_1 \in S_1, so has 
		  degree s_1. It is symmetric, so has only s_1 / 2 + 1 
		  distinct coefficients. The sequence h_j will be stored in 
		  the same memory and won't be a monic polynomial, so the 
		  leading 1 monomial of F will be stored explicitly. Hence we 
		  need s_1 / 2 + 1 entries. */
  listz_handle_t R = NULL;  /* Is used only for two-pass convolution, has nr 
			entries. R is only ever referenced if twopass == 1,
			but gcc does not realize that and complains about
			uninitialized value, so we set it to NULL. */
  mpzspm_t ntt_context;
  mpzspv_handle_t g_x_ntt, g_y_ntt, h_x_ntt, h_y_ntt;
  gfp_ext_t b1;
  mpres_t Delta;
  mpz_t mt;   /* All-purpose temp mpz_t */
  mpz_t product;
  mpz_t *product_ptr = NULL;
  char *filename = NULL;
  int youpi = ECM_NO_FACTOR_FOUND;
  int twopass = twopass_param;
  long timetotalstart, realtotalstart, timestart, realstart;

  timetotalstart = cputime ();
  realtotalstart = realtime ();

  ASSERT_ALWAYS (eulerphi64 (params->P) == params->s_1 * params->s_2);
  ASSERT_ALWAYS (params->s_1 < params->l);
  nr = params->l - params->s_1; /* Number of points we evaluate */

  sets_init(&S_1);

  if (make_S_1_S_2 (&S_1, &s2_sumset, 
		  &s2_sumset_size, params) == ECM_ERROR)
      return ECM_ERROR;

  
  mpz_init (mt);
  
  /* Prepare NTT for computing the h sequence, its DCT-I, and the convolution 
     with g. We need NTT of transform length l here. If we want to add 
     transformed vectors, we need to double the modulus. */

  /* If we use disk-based NTT vectors, we always use the two-pass variant */
  if (params->file_stem != NULL)
    twopass = 1;

  if (twopass)
    mpz_set (mt, modulus->orig_modulus);
  else
    mpz_mul_2exp (mt, modulus->orig_modulus, 1UL);
  
  ntt_context = mpzspm_init (params->l, mt);

  if (ntt_context == NULL)
    {
      outputf (OUTPUT_ERROR, "Could not initialise ntt_context, "
               "presumably out of memory\n");
      mpz_clear (mt);
      sets_free (&S_1);
      free (s2_sumset);
      return ECM_ERROR;
    }

  mpzspm_print_CRT_primes (OUTPUT_DEVVERBOSE, "CRT modulus for evaluation = ", 
		    ntt_context);

  if (params->file_stem != NULL)
    do_aio_init(ntt_context->sp_num);

  /* Build F */
  F = build_F_ntt (X, &S_1, params, modulus);
  if (F == NULL)
    {
      sets_free (&S_1);
      free (s2_sumset);
      mpz_clear (mt);
      mpzspm_clear (ntt_context);
      return ECM_ERROR;
    }

  sets_free (&S_1);
  
  gfp_ext_init (b1, modulus);
  mpres_init (Delta, modulus);

  /* Compute Delta and b1_x + b1_y * sqrt(Delta) = X) */
  mpres_sqr (Delta, X, modulus);
  mpres_sub_ui (Delta, Delta, 4UL, modulus);
  mpres_div_2exp (b1->x, X, 1, modulus);
  mpres_set_ui (b1->y, 1UL, modulus);
  mpres_div_2exp (b1->y, b1->y, 1, modulus);
  if (test_verbose (OUTPUT_TRACE))
    {
      mpres_get_z (mt, Delta, modulus);
      outputf (OUTPUT_TRACE, 
	       "Delta = Mod(%Zd, N); w = quadgen (4*lift(Delta)); b_1 = ", mt);
      gfp_ext_print (b1, modulus, OUTPUT_TRACE);
      outputf (OUTPUT_TRACE, "; /* PARI */\n");
      outputf (OUTPUT_TRACE, "X == b_1 + 1/b_1 /* PARI C */\n");
    }

  /* Allocate remaining memory for h_ntt */
  if (params->file_stem != NULL)
    {
      filename = malloc(strlen(params->file_stem) + 5);
      sprintf(filename, "%s.h_x", params->file_stem);
    }
  h_x_ntt = mpzspv_init_handle (filename, params->l / 2 + 1, ntt_context);
  if (params->file_stem != NULL)
    {
      sprintf(filename, "%s.h_y", params->file_stem);
    }
  h_y_ntt = mpzspv_init_handle (filename, params->l / 2 + 1, ntt_context);
  free(filename);

  /* Compute the h_j sequence */
  if (h_x_ntt->storage == 0 && h_y_ntt->storage == 0)
    {
      pp1_sequence_h (NULL, NULL, h_x_ntt, h_y_ntt, F, b1, 0L, 
                      params->s_1 / 2 + 1, params->P, Delta, modulus);
    } else {
      /* We don't have any function to fill two NTT vectors with values in 
         lock-step. We'd have to use mpz_t buffers, or fill vectors one entry
         at a time (very slow). For now we generate x and y-coordinates 
         separately, even though this causes some extra computation. */
      pp1_sequence_h (NULL, NULL, h_x_ntt, NULL, F, b1, 0L, 
                      params->s_1 / 2 + 1, params->P, Delta, modulus);
      pp1_sequence_h (NULL, NULL, NULL, h_y_ntt, F, b1, 0L, 
                      params->s_1 / 2 + 1, params->P, Delta, modulus);
    }
  /* We don't need F(x) any more */
  listz_handle_clear (F);
  F = NULL;

  /* compute the forward transform of h and store the distinct coefficients 
     in h_ntt */
  if (twopass)
    {
      char *filename = NULL;
      
      if (params->file_stem != NULL)
        {
          filename = malloc(strlen(params->file_stem) + 5);
          sprintf(filename, "%s.g_x", params->file_stem);
        }
      g_x_ntt = mpzspv_init_handle (filename, params->l, ntt_context);
      g_y_ntt = g_x_ntt;
      if (params->file_stem != NULL)
        sprintf(filename, "%s.R", params->file_stem);
      R = listz_handle_init2 (filename, nr, modulus->orig_modulus);
      free (filename);
    }
  else
    {
      g_x_ntt = mpzspv_init_handle (NULL, params->l, ntt_context);
      g_y_ntt = mpzspv_init_handle (NULL, params->l, ntt_context);
    }
  
  /* Compute DCT-I of h_x and h_y */
  outputf (OUTPUT_VERBOSE, "Computing DCT-I of h_x");
#ifdef _OPENMP
  outputf (OUTPUT_VERBOSE, " using %d threads", omp_get_thread_limit());
#endif
  timestart = cputime ();
  realstart = realtime ();
  mpzspv_to_dct1 (h_x_ntt, h_x_ntt, params->s_1 / 2 + 1, params->l / 2 + 1);
  print_elapsed_time (OUTPUT_VERBOSE, timestart, realstart);

  outputf (OUTPUT_VERBOSE, "Computing DCT-I of h_y");
#ifdef _OPENMP
  outputf (OUTPUT_VERBOSE, " using %d threads", omp_get_thread_limit());
#endif
  timestart = cputime ();
  realstart = realtime ();
  mpzspv_to_dct1 (h_y_ntt, h_y_ntt, params->s_1 / 2 + 1, params->l / 2 + 1);
  print_elapsed_time (OUTPUT_VERBOSE, timestart, realstart);

  if (test_verbose (OUTPUT_RESVERBOSE))
    {
      mpz_init (product);
      product_ptr = &product;
    }

  for (l = 0; l < params->s_2 && youpi == ECM_NO_FACTOR_FOUND; l++)
    {
      const uint64_t M = params->l - 1 - params->s_1 / 2;

      outputf (OUTPUT_VERBOSE, "Multi-point evaluation %" PRIu64 
                        " of %" PRIu64 ":\n", l + 1, params->s_2);
      if (twopass)
	{
	  /* Two-pass variant. Two separate convolutions, 
	     then addition in Z/NZ */
	  pp1_sequence_g (NULL, NULL, g_x_ntt, NULL, b1, params->P, 
			  Delta, M, params->l, params->m_1, s2_sumset[l], 
			  modulus);

	  /* Do the convolution product of g_x * h_x */
	  outputf (OUTPUT_VERBOSE, "Computing g_x*h_x");
#ifdef _OPENMP
          outputf (OUTPUT_VERBOSE, " using %d threads", omp_get_thread_limit());
#endif
	  timestart = cputime ();
	  realstart = realtime ();
          mpzspv_mul_ntt (g_x_ntt, 0, g_x_ntt, 0, params->l, h_x_ntt, 0, params->l / 2 + 1,
              params->l, 0, 0, 
              NTT_MUL_STEP_FFT1 + NTT_MUL_STEP_MULDCT + NTT_MUL_STEP_IFFT);
	  /* Store the product coefficients we want in R */
	  if (R->storage == 0)
	    mpzspv_fromto_mpzv (g_x_ntt, params->s_1 / 2, nr, NULL, NULL, NULL, R->data.mem);
          else
            {
              void *buf = (file_word_t *) malloc (R->words * sizeof(file_word_t));
              ASSERT_ALWAYS(buf != NULL);
              listz_handle_state_t state = {R, 0, buf};
              mpzspv_fromto_mpzv (g_x_ntt, params->s_1 / 2, nr, NULL, NULL, &listz_handle_write, &state);
              free(buf);
            }
	  print_elapsed_time (OUTPUT_VERBOSE, timestart, realstart);

	  /* Compute g_y sequence */
	  pp1_sequence_g (NULL, NULL, NULL, g_y_ntt, b1, params->P, 
			  Delta, M, params->l, params->m_1, s2_sumset[l], 
			  modulus);
          ASSERT (g_y_ntt->storage == 1 || mpzspv_verify_in (g_y_ntt, 0, params->l));
	  
	  /* Do the convolution product of g_y * (Delta * h_y) */
	  outputf (OUTPUT_VERBOSE, "Computing g_y*h_y");
#ifdef _OPENMP
          outputf (OUTPUT_VERBOSE, " using %d threads", omp_get_thread_limit());
#endif
	  timestart = cputime ();
	  realstart = realtime ();
          mpzspv_mul_ntt (g_y_ntt, 0, g_y_ntt, 0, params->l, h_y_ntt, 0, params->l / 2 + 1,
              params->l, 0, 0, 
              NTT_MUL_STEP_FFT1 + NTT_MUL_STEP_MULDCT + NTT_MUL_STEP_IFFT);
	  print_elapsed_time (OUTPUT_VERBOSE, timestart, realstart);
	  
	  /* Compute product of sum of coefficients and gcd with N */
	  ntt_gcd (mt, product_ptr, g_y_ntt, params->s_1 / 2, R, nr, modulus);
	}
      else
	{
	  /* One-pass variant. Two forward transforms and point-wise products,
	     then addition and single inverse transform */
	  pp1_sequence_g (NULL, NULL, g_x_ntt, g_y_ntt, b1, params->P, 
			  Delta, M, params->l, params->m_1, s2_sumset[l], 
			  modulus);

	  outputf (OUTPUT_VERBOSE, "Computing forward NTT of g_x");
#ifdef _OPENMP
          outputf (OUTPUT_VERBOSE, " using %d threads", omp_get_thread_limit());
#endif
	  timestart = cputime ();
	  realstart = realtime ();
          mpzspv_mul_ntt (g_x_ntt, 0, g_x_ntt, 0, params->l, h_x_ntt, 0, params->l / 2 + 1,
              params->l, 0, 0, 
              NTT_MUL_STEP_FFT1 + NTT_MUL_STEP_MULDCT);
	  print_elapsed_time (OUTPUT_VERBOSE, timestart, realstart);
	  
	  outputf (OUTPUT_VERBOSE, "Computing forward NTT of g_y");
#ifdef _OPENMP
          outputf (OUTPUT_VERBOSE, " using %d threads", omp_get_thread_limit());
#endif
	  timestart = cputime ();
	  realstart = realtime ();
          mpzspv_mul_ntt (g_y_ntt, 0, g_y_ntt, 0, params->l, h_y_ntt, 0, params->l / 2 + 1,
              params->l, 0, 0, 
              NTT_MUL_STEP_FFT1 + NTT_MUL_STEP_MULDCT);
	  print_elapsed_time (OUTPUT_VERBOSE, timestart, realstart);
	  
	  outputf (OUTPUT_VERBOSE, "Adding and computing inverse NTT of sum");
#ifdef _OPENMP
          outputf (OUTPUT_VERBOSE, " using %d threads", omp_get_thread_limit());
#endif
	  timestart = cputime ();
	  realstart = realtime ();
	  mpzspv_add (g_x_ntt, (spv_size_t) 0, g_x_ntt, (spv_size_t) 0, 
	              g_y_ntt, (spv_size_t) 0, params->l);
          mpzspv_mul_ntt (g_x_ntt, 0, g_x_ntt, 0, params->l, NULL, 0, 0,
              params->l, 0, 0, NTT_MUL_STEP_IFFT);
	  print_elapsed_time (OUTPUT_VERBOSE, timestart, realstart);
	  
	  ntt_gcd (mt, product_ptr, g_x_ntt, params->s_1 / 2, NULL, nr, modulus);
	}
      
      outputf (OUTPUT_RESVERBOSE, "Product of R[i] = %Zd\n", product);

      if (mpz_cmp_ui (mt, 1UL) > 0)
	{
	  mpz_set (f, mt);
	  youpi = ECM_FACTOR_FOUND_STEP2;
	  break;
	}
    }

  if (test_verbose (OUTPUT_RESVERBOSE))
    {
      product_ptr = NULL;
      mpz_clear (product);
    }
  mpzspv_clear_handle (g_x_ntt);
  if (twopass)
    listz_handle_clear (R);
  else
    mpzspv_clear_handle (g_y_ntt);
  mpzspv_clear_handle (h_x_ntt);
  mpzspv_clear_handle (h_y_ntt);
  mpzspm_clear (ntt_context);
  mpz_clear (mt);
  gfp_ext_clear (b1, modulus);
  mpres_clear (Delta, modulus);
  free (s2_sumset);
 
  outputf (OUTPUT_NORMAL, "Step 2");
  /* In normal output mode, print only cpu time as we always have.
     In verbose mode, print real time as well if we used multi-threading */
  if (test_verbose (OUTPUT_VERBOSE))
    print_elapsed_time (OUTPUT_NORMAL, timetotalstart, realtotalstart);
  else
    print_elapsed_time (OUTPUT_NORMAL, timetotalstart, 0L);

  return youpi;
}

#ifdef TESTDRIVE
static void
testdrive_producer (void * const p, mpz_t r)
{
  mpz_t *s = p;
  mpz_set (r, s[0]);
  /* Multiply s[0] by 2 (mod s[1]) */
  mpz_mul_2exp (s[0], s[0], 1);
  if (mpz_cmp (s[0], s[1]) >= 0)
    mpz_sub (s[0], s[0], s[1]);
}

int main (int argc, char **argv)
{
  uint64_t len = 0;
  listz_handle_t F;  
  mpzspv_handle_t ntt_handle;
  mpz_t N;
  mpmod_t modulus;
  mpzspm_t mpzspm;
  char *filename = NULL;
  int i;
  int do_ntt = 0, do_pwmul = 0, do_intt = 0, do_gcd = 0, do_fill = 0;
  long timestart, realstart, timediff, realdiff;
  
  for (i = 1; i < argc; i++)
    {
      if (argc > i+1 && strcmp (argv[i], "-l") == 0)
        {
          len = strtoul (argv[i + 1], NULL, 10);
          i++;
          continue;
        }
      if (argc > i+1 && strcmp (argv[i], "-treefile") == 0)
        {
          filename = argv[i + 1];
          i++;
          continue;
        }
      if (strcmp (argv[i], "-dontt") == 0)
        {
          do_ntt = 1;
          continue;
        }
      if (strcmp (argv[i], "-dopwmul") == 0)
        {
          do_pwmul = 1;
          continue;
        }
      if (strcmp (argv[i], "-dointt") == 0)
        {
          do_intt = 1;
          continue;
        }
      if (strcmp (argv[i], "-dogcd") == 0)
        {
          do_gcd = 1;
          continue;
        }
      if (strcmp (argv[i], "-dofill") == 0)
        {
          do_fill = 1;
          continue;
        }
      mpz_init (N);
      mpz_set_str (N, argv[i], 10);
      mpmod_init (modulus, N, ECM_MOD_DEFAULT);
    }

    if (do_ntt || do_pwmul || do_intt || do_gcd || do_fill)
      {
        ASSERT_ALWAYS (len != 0);
        mpzspm = mpzspm_init (len, N);
        do_aio_init (mpzspm->sp_num);
        ntt_handle = mpzspv_init_handle (filename, len, mpzspm);
      }

   if (do_gcd)
     {
       mpz_t f;
       uint64_t bytes;
       
       mpz_init (f);
       timestart = cputime ();
       realstart = realtime ();
       ntt_gcd (f, NULL, ntt_handle, 0, NULL, len, modulus);
       timediff = cputime () - timestart;
       realdiff = realtime () - realstart;
       
       printf("GCD took %lu ms, %lu ms elaped\n", timediff, realdiff);
       bytes = len * mpzspm->sp_num * sizeof(sp_t);
       printf("%lu * %u * %lu = %lu bytes, %f MB/s\n", 
               len, mpzspm->sp_num, sizeof(sp_t), bytes, 
               bytes / (realdiff * 0.001) / 1048576.);
       
       mpz_clear (f);
     }

   if (do_fill)
     {
       mpz_t s[2];
       uint64_t bytes;

       mpz_init (s[0]);
       mpz_init (s[1]);
       mpz_set_ui (s[0], 3);
       mpz_set (s[1], N);
       timestart = cputime ();
       realstart = realtime ();
       mpzspv_fromto_mpzv (ntt_handle, (spv_size_t) 0, len, 
                           &testdrive_producer, s, NULL, NULL);
       timediff = cputime () - timestart;
       realdiff = realtime () - realstart;
       printf("Fill took %lu ms, %lu ms elaped\n", timediff, realdiff);
       bytes = len * mpzspm->sp_num * sizeof(sp_t);
       printf("%lu * %u * %lu = %lu bytes, %f MB/s\n", 
               len, mpzspm->sp_num, sizeof(sp_t), bytes, 
               bytes / (realdiff * 0.001) / 1048576.);
       mpz_clear (s[0]);
       mpz_clear (s[1]);
     }

    if (do_ntt || do_pwmul || do_intt || do_gcd || do_fill)
      {
        mpzspv_clear_handle (ntt_handle);
        mpzspm_clear (mpzspm);
      }
   mpz_clear (N);
   mpmod_clear (modulus);
}
#endif
