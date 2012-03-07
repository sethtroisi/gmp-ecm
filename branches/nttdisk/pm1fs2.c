/* 
  Implementation of fast stage 2 for P-1 and P+1 as described in
  "Improved Stage 2 to $P\pm{}1$ Factoring Algorithms" by
  Peter L. Montgomery and Alexander Kruppa, ANTS 2008 (8th Algorithmic 
  Number Theory Symposium).
   
  Copyright 2007, 2008 Alexander Kruppa.
  NTT functions are based on code Copyright 2005 Dave Newman.

  This file is part of the ECM Library.

  The ECM Library is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or (at your
  option) any later version.

  The ECM Library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
  License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with the ECM Library; see the file COPYING.LIB.  If not, write to
  the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
  MA 02110-1301, USA.
*/

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

/* TODO:
   - move functions into their proper files (i.e. NTT functions etc.)
   - later: allow storing NTT vectors on disk
*/

/* Define TEST_ZERO_RESULT to test if any result of the multipoint
   evaluation is equal to zero. If the modulus is composite, this
   happening might indicate a problem in the evalutaion code */
#define TEST_ZERO_RESULT

const int pari = 0;


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
ntt_sqr_reciprocal (mpzv_t, const mpzv_t, mpzspv_t, FILE **, const spv_size_t, 
		    const mpzspm_t);

static void
print_elapsed_time (int verbosity, long cpu_start, 
		    ATTRIBUTE_UNUSED long real_start)
{
#ifdef _OPENMP
  if (real_start != 0L)
    {
      outputf (verbosity, " took %lums (%lums real)\n", 
	       elltime (cpu_start, cputime()), 
	       elltime (real_start, realtime()));
      return;
    }
#endif
  outputf (verbosity, " took %lums\n", elltime (cpu_start, cputime()));
}


static void
list_output_poly (listz_t l, uint64_t len, int monic, int symmetric,
		  char *prefix, char *suffix, int verbosity)
{
  uint64_t i;

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

  if (monic)
    {
      if (symmetric)
	outputf (verbosity, "(x^%" PRIu64 " + x^-%" PRIu64 ") + ", len, len);
      else
	outputf (verbosity, "x^%" PRIu64 " + ", len);
    }
  for (i = len - 1; i > 0; i--)
    if (symmetric)
      outputf (verbosity, "Mod(%Zd,N) * (x^%" PRIu64 " + x^-%" PRIu64 ") + ", 
               l[i], i, i);
    else
      outputf (verbosity, "Mod(%Zd,N) * x^%" PRIu64 " + ", l[i], i);
  outputf (verbosity, "Mod(%Zd,N)", l[0]);
  if (suffix != NULL)
    outputf (verbosity, suffix);
}


/* Multiply P[i] by r^{k(deg-i)}, for 0 <= i <= deg. Needs 3 entries in tmp. */
/* I.e., let P(x) = x^deg + \sum_{i=0}^{deg - 1} P[i] * x^i. The output is 
   R(x) = x^deg + \sum_{i=0}^{deg - 1} R[i] * x^i = r^(k deg) P(r^{-k} x). */
/* The input and output polynomials are monic and have the leading monomial
   implicit, i.e. not actually stored in the array of coefficients. */
/* Returns 0 if a modular inversion failed (in which case R is left 
   unchanged), 1 otherwise */

static int ATTRIBUTE_UNUSED
list_scale_rev (listz_t R, listz_t S, mpz_t r, long k, uint64_t deg, 
		mpz_t modulus, listz_t tmp, 
		ATTRIBUTE_UNUSED const uint64_t tmplen)
{
  uint64_t i;

  ASSERT (tmplen >= 3);
  mpz_powm_ui (tmp[0], r, (unsigned long) labs (k), modulus);
  if (k < 0)
    {
      if (!mpz_invert (tmp[0], tmp[0], modulus)) /* FIXME: get rid of this! */
	return 0;
    }
  /* Here, tmp[0] = r^k */
  mpz_set (tmp[1], tmp[0]);
  /* mpz_set (R[deg], S[deg]); Leading monomial is not stored! */
  for (i = 1; i + 1 <= deg; i++)
    {
      /* Here, tmp[1] = r^(ki) */
      mpz_mul (tmp[2], S[deg-i], tmp[1]);
      mpz_mod (R[deg-i], tmp[2], modulus);
      mpz_mul (tmp[2], tmp[1], tmp[0]);  /* FIXME, avoid unnecessary mul */
      mpz_mod (tmp[1], tmp[2], modulus); /* at end of loop */
    }
  if (i <= deg)
    {
      mpz_mul (tmp[2], S[deg-i], tmp[1]);
      mpz_mod (R[deg-i], tmp[2], modulus);
    }

  return 1;
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
    mpz_mul_2exp (S[0], S[0], 1UL);
	
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


/* Multiply a (possibly monic) polynomial A of length k * len with a 
   (possibly monic) polynomial B of length len. R may be identical to A. */

static void ATTRIBUTE_UNUSED
list_mul_blocks (listz_t R, const listz_t A, int monicA, const listz_t B, 
		 int monicB, const uint64_t len, const uint64_t k,
		 listz_t tmp, ATTRIBUTE_UNUSED const uint64_t tmplen)
{
  uint64_t j;
  
  if (k == 0 || len == 0)
    return;

  ASSERT (R != B);
  ASSERT (tmplen >= 3 * len + list_mul_mem (len));

  /* Do first piece of A */
  list_mul (tmp, A, len, (monicA && k == 1), B, len, monicB, tmp + 2 * len);
  list_set (R, tmp, len); /* May overwrite A[0 ... len-1] */
  list_swap (tmp, tmp + len, len); /* Move high part to tmp[0 ... len-1] */
  
  for (j = 1; j < k; j++) /* Process the remaining k-1 pieces of A */
    {
      list_mul (tmp + len, 
		A + j * len, len, (monicA && j + 1 == k),
		B, len, monicB, tmp + 3 * len);
      /* Add low part of this product and previous product's high part */
      list_add (A + j * len, tmp, tmp + len, len);
      list_swap (tmp, tmp + 2 * len, len); /* Move this product's high 
					      part to beginning of tmp */
    }

  list_set (A + j * len, tmp, len); /* Move the high part of last product */
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
      gmp_printf ("Chebyshev_V(%ld, Mod(%Zd,N)) == ", (long)k, tz);
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
          mpres_mul (R, R, R, modulus);
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
  mpres_mul (Vi1, S, S, modulus);
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
	  mpres_mul (Vi1, Vi1, Vi1, modulus);
	  mpres_sub (Vi1, Vi1, V0, modulus);
	}
      else
	{
	  /* i' = 2i. 
	     V_{i'} = V_{2i} = {V_i}^2 - V0.
	     V_{i'+1} = V_{2i + 1} = V_{i+1 + i} = V_{i+1} * V_{i} - V_1 */
	  mpres_mul (Vi1, Vi, Vi1, modulus);
	  mpres_sub (Vi1, Vi1, S, modulus);

	  mpres_mul (Vi, Vi, Vi, modulus);
	  mpres_sub (Vi, Vi, V0, modulus);
	}
      j >>= 1;
    }

  /* Least significant bit of uk is always 1 */
  mpres_mul (Vi, Vi, Vi1, modulus);
  mpres_sub (Vi, Vi, S, modulus);

  while (po2-- > 0)
    {
      mpres_mul (Vi, Vi, Vi, modulus);
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
      gmp_printf ("%Zd\n", tz);
      mpz_clear (tz);
    }
}

/* 
  Computes U_k(S), where the Chebyshev polynomial U_k(X) is defined by 
  U_k(X + 1/X) = (X^k - 1/X^k) / (X - 1/X)
  If R1 != NULL, stores U_{k+1}(S) there
*/

static void
U (mpres_t R, mpres_t R1, const mpres_t S, const long k, mpmod_t modulus)
{
  mpres_t V0, Vi, Vi1, Ui, Ui1, t;
  unsigned long j, uk;

  if (k == 0L)
    {
      mpres_set_ui (R, 0UL, modulus); /* U_0 = 0 */
      if (R1 != NULL)
	mpres_set_ui (R1, 1UL, modulus); /* U_1 = 1 */
      return;
    }

  uk = labs (k);

  if (uk == 1UL)
    {
      mpres_set_ui (R, 1UL, modulus);
      if (k == -1)
	mpres_neg (R, R, modulus);
      
      if (R1 != NULL)
	{
	  if (k == -1)
	    mpres_set_ui (R1, 0UL, modulus);
	  else
	    mpres_set (R1, S, modulus); /* U_2(S) = S */
	}

      return;
    }

  if (0)
    {
      mpz_t tz;
      mpz_init (tz);
      mpres_get_z (tz, S, modulus);
      gmp_printf ("Chebyshev_U(%ld, Mod(%Zd,N)) == ", k, tz);
      mpz_clear (tz);
    }

  mpres_init (V0, modulus);
  mpres_init (Vi, modulus);
  mpres_init (Vi1, modulus);
  mpres_init (Ui, modulus);
  mpres_init (Ui1, modulus);
  mpres_init (t, modulus);

  for (j = 1UL; j <= uk / 2UL; j <<= 1);

  mpres_set_ui (Ui, 1UL, modulus);   /* Ui = U_1(S) = 1 */
  mpres_set (Ui1, S, modulus);       /* Ui1 = U_2(S) = S */
  mpres_add (V0, Ui, Ui, modulus);   /* V0 = V_0(S) = 2 */
  mpres_set (Vi, S, modulus);        /* Vi = V_1(S) = S */
  mpres_mul (Vi1, Vi, Vi, modulus);
  mpres_sub (Vi1, Vi1, V0, modulus); /* Vi1 = V_2(S) = S^2 - 2 */
  j >>= 1; /* i = 1 */

  while (j != 0)
    {
      if ((uk & j) == 0UL)
	{
	  mpres_mul (Vi1, Vi1, Vi, modulus);
	  mpres_sub (Vi1, Vi1, S, modulus); /* V_{2i+1} = V_{i+1} V_i - V_1 */
	  /* U_{2i+1} = (U_{i+1} + U_i) (U_{i+1} - U_i) */
	  mpres_sub (t, Ui1, Ui, modulus);
	  mpres_add (Ui1, Ui1, Ui, modulus);
	  mpres_mul (Ui1, Ui1, t, modulus); 
	  mpres_mul (Ui, Ui, Vi, modulus); /* U_{2n} = U_n V_n */
	  mpres_mul (Vi, Vi, Vi, modulus);
	  mpres_sub (Vi, Vi, V0, modulus); /* V_{2n} = V_n^2 - 2 */
	}
      else
	{
	  /* U_{2i+1} = (U_{i+1} + U_i) (U_{i+1} - U_i) */
	  mpres_sub (t, Ui1, Ui, modulus);
	  mpres_add (Ui, Ui, Ui1, modulus);
	  mpres_mul (Ui, Ui, t, modulus);
	  mpres_mul (Ui1, Ui1, Vi1, modulus); /* U_{2n+2} = U_{n+1} V_{n+1} */
	  mpres_mul (Vi, Vi, Vi1, modulus);
	  mpres_sub (Vi, Vi, S, modulus); /* V_{2i+1} = V_{i+1} V_i - V_1 */
	  mpres_mul (Vi1, Vi1, Vi1, modulus);
	  mpres_sub (Vi1, Vi1, V0, modulus); /* V_{2n+2} = V_{n+1}^2 - 2 */
	}
      j >>= 1;
    }

  if (k > 0)
    mpres_set (R, Ui, modulus);
  else
    mpres_neg (R, Ui, modulus);

  if (R1 != NULL)
    {
      /* Here k != -1,0,1, so k+1 is negative iff k is */
      if (k > 0)
	mpres_set (R1, Ui1, modulus);
      else
	mpres_neg (R1, Ui1, modulus);
    }

  mpres_clear (V0, modulus);
  mpres_clear (Vi, modulus);
  mpres_clear (Vi1, modulus);
  mpres_clear (Ui, modulus);
  mpres_clear (Ui1, modulus);
  mpres_clear (t, modulus);

  if (0)
    {
      mpz_t tz;
      mpz_init (tz);
      mpres_get_z (tz, R, modulus);
      gmp_printf ("%Zd\n", tz);
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

   If NTT is used, needs 4 * deg + 3 entries in tmp.
   If no NTT is used, needs 4 * deg + 2 + (memory use of list_sqr_reciprocal)
*/

static void ATTRIBUTE_UNUSED
list_scale_V (listz_t R, const listz_t F, const mpres_t Q, 
              const uint64_t deg, mpmod_t modulus, listz_t tmp, 
              const uint64_t tmplen, 
	      mpzspv_t dct, FILE **dct_files, const mpzspm_t ntt_context)
{
  mpres_t Vt;
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
  
  /* Make sure newtmplen does not underflow */
  ASSERT_ALWAYS (tmplen >= 4 * deg + 2);
#ifdef WANT_ASSERT
  mpz_init (leading);
  mpz_mul (leading, F[deg], F[deg]);
  mpz_mod (leading, leading, modulus->orig_modulus);
#endif

  /* Generate V_1(Q)/2 ... V_{deg}(Q)/2, multiply by f_i to form coefficients 
     of G(x). Square the symmetric G(x) polynomial. */

  outputf (OUTPUT_TRACE, "\n/* list_scale_V */ N=%Zd; deg = %lu;\n", 
           modulus->orig_modulus, deg);
  if (test_verbose(OUTPUT_TRACE))
    {
      mpres_t out_t;
      mpres_init (out_t, modulus);
      mpres_get_z (out_t, Q, modulus);
      outputf (OUTPUT_TRACE, "/* list_scale_V */ Q = Mod(%Zd,N);\n", out_t);
      mpres_clear (out_t, modulus);
    }
  list_output_poly (F, deg + 1, 0, 1, "/* list_scale_V */ F(x) = ", ";\n", 
		    OUTPUT_TRACE);

  /* Compute G[i] = V_i(Q)/2 * F[i] for i = 0, ..., deg.
     For i=0, V_0(Q) = 2, so G[0] = F[0], 
     which leaves deg entries to process */

  mpz_mod (G[0], F[0], modulus->orig_modulus);

#if defined(_OPENMP)
#pragma omp parallel if (deg > 1000)
#endif
  {
    mpmod_t modulus_local;
    uint64_t l, start_i;
    mpres_t Vi, Vi_1;
    
    get_chunk (&start_i, &l, deg);
    start_i++;

    mpmod_init_set (modulus_local, modulus);
    mpres_init (Vi_1, modulus_local);
    mpres_init (Vi, modulus_local);
    
    V (Vi, Q, start_i, modulus_local);
    mpres_div_2exp (Vi, Vi, 1, modulus_local);
    V (Vi_1, Q, start_i - 1, modulus_local);
    mpres_div_2exp (Vi_1, Vi_1, 1, modulus_local);
    scale_by_chebyshev (G + start_i, F + start_i, NULL, NULL, l, 
                        modulus_local, Q, Vi_1, Vi);
    
    mpres_clear (Vi_1, modulus_local);
    mpres_clear (Vi, modulus_local);
    mpmod_clear (modulus_local);
  }


  list_output_poly (G, deg + 1, 0, 1, "/* list_scale_V */ G(x) = ", ";\n", 
		    OUTPUT_TRACE);

  /* Now square the G polynomial in G[0 .. deg], put result in
     G[0 .. 2*deg] */

  /* Bugfix: ks_multiply() does not like negative coefficients. FIXME */

  for (i = 0; i <= deg; i++)
    {
      ASSERT (mpz_sgn (G[i]) >= 0 && mpz_cmp(G[i], modulus->orig_modulus) < 0);
    }
  
  if (ntt_context != NULL)
    ntt_sqr_reciprocal (G, G, dct, dct_files, deg + 1, ntt_context);
  else
    list_sqr_reciprocal (G, G, deg + 1, modulus->orig_modulus, 
                         newtmp, newtmplen);

  list_output_poly (G, 2 * deg + 1, 0, 1, "/* list_scale_V */ G(x)^2 == ", 
		    "\n", OUTPUT_TRACE);

  /* Compute H[i-1] = U_i(Q)/2 * F[i] for i = 1, ..., deg */

#if defined(_OPENMP)
#pragma omp parallel if (deg > 1000)
#endif
  {
    mpmod_t modulus_local;
    uint64_t l, start_i;
    mpres_t Ui, Ui_1;
    
    get_chunk (&start_i, &l, deg);
    
    mpmod_init_set (modulus_local, modulus);
    mpres_init (Ui_1, modulus_local);
    mpres_init (Ui, modulus_local);
    
    U (Ui_1, Ui, Q, start_i, modulus_local);
    mpres_div_2exp (Ui, Ui, 1, modulus_local);
    mpres_div_2exp (Ui_1, Ui_1, 1, modulus_local);
    
    scale_by_chebyshev (H + start_i, F + start_i + 1, NULL, NULL, l, 
                        modulus_local, Q, Ui_1, Ui);
    
    mpres_clear (Ui_1, modulus_local);
    mpres_clear (Ui, modulus_local);
    mpmod_clear (modulus_local);
  }

  
  /* Convert H to standard basis */
  /* We can do it in-place with H - 1 = H_U. */

  for (i = deg; i >= 3; i--)
    {
      mpz_add (H[i - 3], H[i - 3], H[i - 1]);
      if (mpz_cmp (H[i - 3], modulus->orig_modulus) >= 0)
        mpz_sub (H[i - 3], H[i - 3], modulus->orig_modulus);
    }
  
  /* U_2(X+1/X) = (X^2 - 1/X^2)/(X-1/X) = X+1/X = V_1(X+1/X),
     so no addition occurs here */
  /* if (deg >= 2)
     mpz_set (H[1], H[1]); Again, a no-op. */
  
  /* U_1(X+1/X) = 1, so this goes to coefficient of index 0 in std. basis */
  /* mpz_set (H[0], H[0]); Another no-op. */
  
  /* Now H[0 ... deg-1] contains the deg coefficients in standard basis
     of symmetric H(X) of degree 2*deg-2. */
  
  list_output_poly (H, deg, 0, 1, "/* list_scale_V */ H(x) = ", ";\n",
		    OUTPUT_TRACE);

  /* Square the symmetric H polynomial of degree 2*deg-2 (i.e. with deg 
     coefficents in standard basis in H[0 ... deg-1]) */

  /* Bugfix: ks_multiply() does not like negative coefficients. */

  for (i = 0; i <= deg; i++)
    if (mpz_sgn (H[i]) < 0)
      {
	mpz_add (H[i], H[i], modulus->orig_modulus);
	if (mpz_sgn (H[i]) < 0)
	  {
	    outputf (OUTPUT_ERROR, "list_scale_V: H[%lu] still negative\n", i);
	    mpz_mod (H[i], H[i], modulus->orig_modulus);
	  }
      }

  if (ntt_context != NULL)
    ntt_sqr_reciprocal (H, H, dct, dct_files, deg, ntt_context);
  else
    list_sqr_reciprocal (H, H, deg, modulus->orig_modulus, 
  		         newtmp, newtmplen);

  /* Now there are the 2*deg-1 coefficients in standard basis of a 
     symmetric polynomial of degree 4*deg - 4 in H[0 ... 2*deg-2] */

  list_output_poly (H, 2*deg - 1, 0, 1, "/* list_scale_V */ H(x)^2 == ", "\n",
		    OUTPUT_TRACE);

  /* Multiply by Q^2-4 */
  mpres_init (Vt, modulus);
  mpres_sqr (Vt, Q, modulus);
  mpres_sub_ui (Vt, Vt, 4, modulus);
  if (test_verbose(OUTPUT_TRACE))
  {
    mpres_t out_t;
    mpres_init (out_t, modulus);
    mpres_get_z (out_t, Vt, modulus);
    outputf (OUTPUT_TRACE, "/* list_scale_V */ Q^2-4 == Mod(%Zd,N)\n", out_t);
    mpres_clear (out_t, modulus);
  }

#if defined(_OPENMP)
#pragma omp parallel if (deg > 1000)
  {
    mpmod_t modulus_local;
    int64_t i; /* OpenMP insists on signed loop iteration var :( */
    
    mpmod_init_set (modulus_local, modulus);
    
#pragma omp for
    for (i = 0; (unsigned long) i <= 2 * deg - 2; i++)
      mpres_mul_z_to_z (H[i], Vt, H[i], modulus_local);
    mpmod_clear (modulus_local);
  }
#else
  for (i = 0; (unsigned long) i <= 2 * deg - 2; i++)
    mpres_mul_z_to_z (H[i], Vt, H[i], modulus);
#endif

  list_output_poly (H, 2 * deg - 1, 0, 1, "/* list_scale_V */ "
		    "H(x)^2*(Q^2-4) == ", "\n", OUTPUT_TRACE);


  /* Multiply by (X - 1/X)^2 = X^2 - 2 + 1/X^2 and subtract from G */
  ASSERT (newtmplen > 0UL);
  if (deg == 1)
    {
      /* H(X) has degree 2*deg-2 = 0, so H(X) = h_0
	 H(X) * (X - 1/X)^2 = -2 h_0 + h_0 V_2(Y)  */
      mpz_mul_2exp (newtmp[0], H[0], 1UL);
      mpz_add (G[0], G[0], newtmp[0]); /* G[0] -= -2*H[0] */
      mpz_sub (G[2], G[2], H[0]);
    }
  else if (deg == 2)
    {
      /* H(X) has degree 2*deg-2 = 2, , so 
	 H(X) = h_0 + h_1 (X+1/X) + h_2 (X^2+1/X^2)

	 H(X) * (X - 1/X)^2 =
	 -2*(h_0 - h_2) - h_1 * V_1(Y) + (h_0 - 2*h_2) * V_2(Y) + 
	 h_1 * V_3(Y) + h_2 * V_4(Y)
      */
      mpz_sub (newtmp[0], H[0], H[2]);          /* h_0 - h_2 */
      mpz_mul_2exp (newtmp[0], newtmp[0], 1UL); /* 2*(h_0 - h_2) */
      mpz_add (G[0], G[0], newtmp[0]);          /* G[0] -= -2*(h_0 - h_2) */

      mpz_add (G[1], G[1], H[1]);               /* G[1] -= -h_1 */
      mpz_sub (newtmp[0], newtmp[0], H[0]);     /* h_0 - 2*h_2 */
      mpz_sub (G[2], G[2], newtmp[0]);          /* G[2] -= h_0 - 2*h_2 */
      mpz_sub (G[3], G[3], H[1]);               /* G[3] -= h_1 */
      mpz_sub (G[4], G[4], H[2]);               /* G[3] -= h_2 */
    }
  else
    {
      /* Let H(X) = h_0 + \sum_{i=1}^{n} h_i V_i(Y), Y = X+1/X. Then
	 (x - 1/x)^2 H(X) = 
	 -2(h_0 - h_2) +
	 (- h_1 + h_3) V_1(Y) +
	 \sum_{i=2}^{n-2} (h_{i-2} - 2h_i + h_{i+2}) V_i(Y) +
	 (h_{n-3} - 2h_{n-1}) V_{n-1}(Y) +
	 (h_{n-2} - 2h_n) V_n(Y) +
	 h_{n-1} V_{n+1}(Y) +
	 h_n V_{n+2}(Y)
	 
	 In our case, n = 2 * deg - 2
      */
      mpz_sub (newtmp[0], H[0], H[2]);
      mpz_mul_2exp (newtmp[0], newtmp[0], 1UL); /* t[0] = 2*(h_0 - h_2) */
      mpz_add (G[0], G[0], newtmp[0]);          /* G[0] -= -2*(h_0 - h_2) */
      
      mpz_add (G[1], G[1], H[1]);
      mpz_sub (G[1], G[1], H[3]); /* G[1] -= -h_1 + h_3 */
      
      for (i = 2; i <= 2 * deg - 4; i++)
	{
	  mpz_mul_2exp (newtmp[0], H[i], 1);
	  mpz_sub (newtmp[0], newtmp[0], H[i - 2]);
	  mpz_sub (newtmp[0], newtmp[0], H[i + 2]); /* 2h_i-h_{i-2}-h_{i+2} */
	  mpz_add (G[i], G[i], newtmp[0]); /* G[i] -= -2h_i+h_{i-2}+h_{i+2} */
	}
      
      for ( ; i <= 2 * deg - 2; i++)
	{
	  mpz_mul_2exp (newtmp[0], H[i], 1UL);
	  mpz_sub (newtmp[0], H[i - 2], newtmp[0]); /* h_{n-3} - 2h_{n-1} */
	  mpz_sub (G[i], G[i], newtmp[0]);
	}
      
      mpz_sub (G[i], G[i], H[i - 2]);
      mpz_sub (G[i + 1], G[i + 1], H[i - 1]);
    }

  for (i = 0; i <= 2 * deg; i++)
    mpz_mod (R[i], G[i], modulus->orig_modulus);

  if (test_verbose (OUTPUT_TRACE))
    for (i = 0; i <= 2 * deg; i++)
      outputf (OUTPUT_TRACE, "list_scale_V: R[%lu] = %Zd\n", i, R[i]);

#ifdef WANT_ASSERT
  mpz_mod (R[2 * deg], R[2 * deg], modulus->orig_modulus);
  ASSERT (mpz_cmp (leading, R[2 * deg]) == 0);
  mpz_clear (leading);
#endif

  mpres_clear (Vt, modulus);
}

/* New, simpler version */
static void
list_scale_V2 (listz_t R, const listz_t F, const mpres_t Q, 
              const uint64_t deg, mpmod_t modulus, listz_t tmp, 
              const uint64_t tmplen, 
	      mpzspv_t dct, FILE **dct_files, const mpzspm_t ntt_context)
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
  
  /* Make sure newtmplen does not underflow */
  ASSERT_ALWAYS (tmplen >= 4 * deg + 2);
#ifdef WANT_ASSERT
  mpz_init (leading);
  mpz_mul (leading, F[deg], F[deg]);
  mpz_mod (leading, leading, modulus->orig_modulus);
#endif

  outputf (OUTPUT_TRACE, "\nN=%Zd; deg = %lu; /* PARI list_scale_V2 */\n", 
           modulus->orig_modulus, deg);
  if (test_verbose(OUTPUT_TRACE))
    {
      mpres_t out_t;
      mpres_init (out_t, modulus);
      mpres_get_z (out_t, Q, modulus);
      outputf (OUTPUT_TRACE, "Q = Mod(%Zd,N); /* PARI list_scale_V2 */\n", out_t);
      mpres_clear (out_t, modulus);
    }
  list_output_poly (F, deg + 1, 0, 1, "F(x) = ", "; /* PARI list_scale_V2 */\n", 
		    OUTPUT_TRACE);

  for (i = 0; i <= deg; i++)
    {
      ASSERT_ALWAYS (mpz_sgn (F[i]) >= 0 && mpz_cmp (F[i], modulus->orig_modulus) < 0);
    }

  if (ntt_context != NULL)
    ntt_sqr_reciprocal (G, F, dct, dct_files, deg + 1, ntt_context);
  else
    list_sqr_reciprocal (G, F, deg + 1, modulus->orig_modulus, 
                         newtmp, newtmplen);

  outputf (OUTPUT_TRACE, "G(x) = F(x)^2;/* PARI list_scale_V2 */\n");
  list_output_poly (G, 2 * deg + 1, 0, 1, "G(x) == ", 
		    " /* PARI list_scale_V2 */\n", OUTPUT_TRACE);

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
                    "; /* PARI list_scale_V2 */\n", OUTPUT_TRACE);

  for (i = 0; i <= deg; i++)
    {
      ASSERT_ALWAYS (mpz_sgn (H[i]) >= 0 && mpz_cmp (H[i], modulus->orig_modulus) < 0);
    }

  list_output_poly (H, deg + 1, 0, 1, "H(x) = ", "; /* PARI list_scale_V2 */\n", 
		    OUTPUT_TRACE);

  if (ntt_context != NULL)
    ntt_sqr_reciprocal (H, H, dct, dct_files, deg + 1, ntt_context);
  else
    list_sqr_reciprocal (H, H, deg + 1, modulus->orig_modulus, 
  		         newtmp, newtmplen);

  list_output_poly (H, 2*deg + 1, 0, 1, "H(x)^2 == ", 
                    " /* PARI list_scale_V2 */\n", OUTPUT_TRACE);

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

  list_output_poly (R, 2 * deg, 0, 1, "R(x) = ", " /* PARI list_scale_V2 */\n", OUTPUT_TRACE);

#ifdef WANT_ASSERT
  mpz_mod (R[2 * deg], R[2 * deg], modulus->orig_modulus);
  ASSERT (mpz_cmp (leading, R[2 * deg]) == 0);
  mpz_clear (leading);
#endif
}


#ifdef WANT_ASSERT
/* Check if l is an (anti-)symmetric, possibly monic, polynomial. 
   Returns -1 if it is (anti-)symmetric, or the smallest index i where 
   l[i] != l[len - 1 + monic - i])
   If anti == 1, the list is checked for symmetry, if it is -1, for
   antisymmetry.
   This function is used only if assertions are enabled.
*/

static long int ATTRIBUTE_UNUSED
list_is_symmetric (listz_t l, uint64_t len, int monic, int anti, 
		   mpz_t modulus, mpz_t tmp)
{
    unsigned long i;

    ASSERT (monic == 0 || monic == 1);
    ASSERT (anti == 1 || anti == -1);

    if (monic && anti == 1 && mpz_cmp_ui (l[0], 1) != 0)
	return 0L;

    if (monic && anti == -1)
      {
	mpz_sub_ui (tmp, modulus, 1);
	if (mpz_cmp (tmp, l[0]) != 0)
	  return 0L;
      }

    for (i = monic; i < len / 2; i++)
      {
	if (anti == -1)
	  {
	    /* Negate (mod modulus) */
	    if (mpz_sgn (l[i]) == 0)
	      {
		if (mpz_sgn (l[len - 1 + monic - i]) != 0)
		  return (long) i;
	      }
	    else
	      {
		mpz_sub (tmp, modulus, l[i]);
		if (mpz_cmp (tmp, l[len - 1 + monic - i]) != 0)
		  return (long) i;
	      }
	  }
	else if (mpz_cmp (l[i], l[len - 1 + monic - i]) != 0)
	    return (long) i;
      }

    return -1L;
}
#endif

/* Evaluate a polynomial of degree n-1 with all coefficients given in F[],
   or of degree n with an implicit leading 1 monomial not stored in F[],
   at x modulo modulus. Result goes in r. tmp needs 2 entries. */

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
poly_from_sets_V (listz_t F, const mpres_t Q, set_list_t *sets, 
		  listz_t tmp, const uint64_t tmplen, mpmod_t modulus,
		  mpzspv_t dct, FILE **dct_files, const mpzspm_t ntt_context)
{
  unsigned long c, i, nr;
  uint64_t deg;
  mpres_t Qt;
  
  ASSERT_ALWAYS (sets->num_sets > 0UL);
  /* Check that the cardinality of first set is 2 */
  ASSERT_ALWAYS (sets->sets[0].card == 2UL);
  /* Check that the first set is symmetric around 0 */
  ASSERT_ALWAYS (sets->sets[0].elem[0] == -sets->sets[0].elem[1]);

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

  V (Qt, Q, sets->sets[0].elem[0], modulus); /* First set in sets is {-k, k} */ 
  V (Qt, Qt, 2UL, modulus);                  /* Qt = V_2k(Q) */
  
  mpres_neg (Qt, Qt, modulus);
  mpres_get_z (F[0], Qt, modulus);
  mpz_set_ui (F[1], 1UL);
  deg = 1UL;
  /* Here, F(x) = (x - r^{2k_1})(x - r^{-2k_1}) / x = 
                  (x^2 - x (r^{2k_1} + r^{-2k_1}) + 1) / x =
		  (x + 1/x) - V_{2k_1}(r + 1/r) */

  for (nr = sets->num_sets - 1; nr > 0; nr--)
    {
      set_t *curr_set = sets->sets + nr;

      /* Assuming the sets are sorted in order of ascending cardinality, 
         we process them back-to-front so the sets of cardinality 2 are 
         processed last, but skipping the first set which we processed 
         already. */
      
      /* Process this set. We assume it is either of cardinality 2, or of 
	 odd cardinality */
      c = curr_set->card;
      outputf (OUTPUT_DEVVERBOSE, " %lu", c);

      if (c == 2UL)
	{
	  /* Check it's symmetric */
	  ASSERT_ALWAYS (curr_set->elem[0] == -curr_set->elem[1]);
	  V (Qt, Q, curr_set->elem[0], modulus);
	  V (Qt, Qt, 2UL, modulus);
	  list_scale_V2 (F, F, Qt, deg, modulus, tmp, tmplen, dct, 
	                dct_files, ntt_context);
	  deg *= 2UL;
	  ASSERT_ALWAYS (mpz_cmp_ui (F[deg], 1UL) == 0); /* Check it's monic */
	}
      else
	{
	  ASSERT_ALWAYS (c % 2UL == 1UL);
	  ASSERT_ALWAYS (curr_set->elem[(c - 1UL) / 2UL] == 0);
	  /* Generate the F(Q^{2k_i} * X)*F(Q^{-2k_i} * X) polynomials.
	     Each is symmetric of degree 2*deg, so each has deg+1 coeffients
	     in standard basis. */
	  for (i = 0UL; i < (c - 1UL) / 2UL; i++)
	    {
              /* Check it's symmetric */
	      ASSERT_ALWAYS (curr_set->elem[i] == -curr_set->elem[c - 1L - i]);
	      V (Qt, Q, curr_set->elem[i], modulus);
	      V (Qt, Qt, 2UL, modulus);
	      ASSERT_ALWAYS (mpz_cmp_ui (F[deg], 1UL) == 0); /* Check it's monic */
	      list_scale_V2 (F + (2UL * i + 1UL) * (deg + 1UL), F, Qt, deg, 
	                    modulus, tmp, tmplen, dct, dct_files, ntt_context);
	      ASSERT_ALWAYS (mpz_cmp_ui (F[(2UL * i + 1UL) * (deg + 1UL) + 2UL * deg], 
	              1UL) == 0); /* Check it's monic */
	    }
	  /* Multiply the polynomials */
	  for (i = 0UL; i < (c - 1UL) / 2UL; i++)
	    {
	      /* So far, we have the product 
		 F(X) * F(Q^{2k_j} * X) * F(Q^{-2k_j} * X), 1 <= j <= i,
		 at F. This product has degree 2 * deg + i * 4 * deg, that is
		 (2 * i + 1) * 2 * deg, which means (2 * i + 1) * deg + 1
		 coefficients in F[0 ... (i * 2 + 1) * deg]. */
	      ASSERT_ALWAYS (mpz_cmp_ui (F[(2UL * i + 1UL) * deg], 1UL) == 0);
	      ASSERT_ALWAYS (mpz_cmp_ui (F[(2UL * i + 1UL) * (deg + 1UL) + 2UL*deg], 
	                          1UL) == 0);
	      list_output_poly (F, (2 * i + 1) * deg + 1, 0, 1, 
				"poly_from_sets_V: Multiplying ", "\n",
				OUTPUT_TRACE);
	      list_output_poly (F + (2 * i + 1) * (deg + 1), 
	                        2UL * deg + 1UL, 0, 1, " and ", "\n", 
	                        OUTPUT_TRACE);
	      list_mul_reciprocal (F, 
		                   F, (2UL * i + 1UL) * deg + 1UL, 
			 	   F + (2UL * i + 1UL) * (deg + 1UL), 
				   2UL * deg + 1UL, modulus->orig_modulus,
				   tmp, tmplen);
	      list_mod (F, F, (2UL * i + 3UL) * deg + 1UL, 
	                modulus->orig_modulus);
	      list_output_poly (F, (2 * i + 3) * deg + 1, 0, 1, 
                                " = ", "\n", OUTPUT_TRACE);
	      ASSERT_ALWAYS (mpz_cmp_ui (F[(2UL * i + 3UL) * deg], 1UL) == 0);
	    }
	  deg *= c;
	}
    }

  mpres_clear (Qt, modulus);
  outputf (OUTPUT_DEVVERBOSE, ")");

  return deg;
}

static int
build_F_ntt (listz_t F, const mpres_t P_1, set_list_t *S_1, 
	     const faststage2_param_t *params, const char *filename, 
	     mpmod_t modulus)
{
  mpzspm_t F_ntt_context;
  mpzspv_t F_ntt;
  FILE **F_files;
  uint64_t tmplen;
  listz_t tmp;
  long timestart, realstart;
  unsigned long i;

  timestart = cputime ();
  realstart = realtime ();
  
  /* Precompute the small primes, primitive roots and inverses etc. for 
     the NTT. The code to multiply wants a 3*k-th root of unity, where 
     k is the smallest power of 2 with k > s_1/2 */
  
  F_ntt_context = mpzspm_init (3UL << ceil_log2 (params->s_1 / 2 + 1), 
			       modulus->orig_modulus);
  if (F_ntt_context == NULL)
    {
      outputf (OUTPUT_ERROR, "Could not initialise F_ntt_context, "
               "presumably out of memory\n");
      return ECM_ERROR;
    }
  
  mpzspm_print_CRT_primes (OUTPUT_DEVVERBOSE, "CRT modulus for building F = ",
		    F_ntt_context);
  
  outputf (OUTPUT_VERBOSE, "Computing F from factored S_1");
  
  tmplen = params->s_1 + 100;
  tmp = init_list2 (tmplen, (unsigned int) abs (modulus->bits));
  if (filename == NULL)
    {
      F_ntt = mpzspv_init (1UL << ceil_log2 (params->s_1 / 2 + 1), F_ntt_context);
      F_files = NULL;
    } else {
      F_ntt = NULL;
      F_files = mpzspv_open_fileset (filename, F_ntt_context);
    }
  
  i = poly_from_sets_V (F, P_1, S_1, tmp, tmplen, modulus, F_ntt, F_files, 
                        F_ntt_context);
  ASSERT_ALWAYS(2 * i == params->s_1);
  ASSERT_ALWAYS(mpz_cmp_ui (F[i], 1UL) == 0);
  
  print_elapsed_time (OUTPUT_VERBOSE, timestart, realstart);
  if (test_verbose (OUTPUT_TRACE))
    {
      for (i = 0; i < params->s_1 / 2 + 1; i++)
	outputf (OUTPUT_TRACE, "f_%lu = %Zd; /* PARI */\n", i, F[i]);
      outputf (OUTPUT_TRACE, "f(x) = f_0");
      for (i = 1; i < params->s_1 / 2 + 1; i++)
	outputf (OUTPUT_TRACE, "+ f_%lu * (x^%lu + x^(-%lu))", i, i, i);
      outputf (OUTPUT_TRACE, "/* PARI */ \n");
    }
  
  clear_list (tmp, tmplen);
  tmp = NULL;
  if (filename == NULL)
    {
      mpzspv_clear (F_ntt, F_ntt_context);
      F_ntt = NULL;
    } else {
      mpzspv_close_fileset (F_files, F_ntt_context);
      F_files = NULL;
    }
  
  mpzspm_clear (F_ntt_context);
  F_ntt_context = NULL;

  return 0;
}

/* Compute g_i = x_0^{M-i} * r^{(M-i)^2} for 0 <= i < l. 
   x_0 = b_1^{2*k_2 + (2*m_1 + 1) * P}. r = b_1^P. 
   Stores the result in g[0 ... l] and/or in g_ntt[offset ... offset + l] */

static void
pm1_sequence_g (listz_t g_mpz, mpzspv_t g_ntt, FILE **ntt_files, const mpres_t b_1, 
                const uint64_t P, const uint64_t M_param, 
		const uint64_t l_param, const mpz_t m_1, 
		const int64_t k_2, mpmod_t modulus_param, 
		const mpzspm_t ntt_context)
{
  mpres_t r[3], x_0, x_Mi;
  mpz_t t, t1;
  uint64_t i;
  long timestart, realstart;
  mpmod_t modulus;
  mpzspv_t tmp_ntt = NULL;
  const size_t buflen = 16384;
  int want_output = 1;

  outputf (OUTPUT_VERBOSE, "Computing g_i");
  outputf (OUTPUT_DEVVERBOSE, "\npm1_sequence_g: P = %" PRIu64
            ", M_param = %" PRIu64 ", l_param = %" PRIu64 
            ", k_2 = %" PRId64 , P, M_param, l_param, k_2);
  outputf (OUTPUT_DEVVERBOSE, ", m_1 = %Zd\n", m_1);
  timestart = cputime ();
  realstart = realtime ();

#ifdef _OPENMP
#pragma omp parallel if (l_param > 100) private(r, x_0, x_Mi, t, t1, i, modulus, want_output) firstprivate(tmp_ntt)
#endif
  {
    uint64_t M, l, offset;
    /* When multi-threading, we adjust the parameters for each thread */
    get_chunk (&offset, &l, l_param);
    M = M_param - offset;

#ifdef _OPENMP
    outputf (OUTPUT_DEVVERBOSE, 
             "pm1_sequence_g: thread %d has l = %" PRIu64 
             ", offset = %" PRIu64 ".\n", omp_get_thread_num(), l, offset);
    
    /* Let only the master thread print stuff */
    want_output = (omp_get_thread_num() == 0);

    if (want_output)
      outputf (OUTPUT_VERBOSE, " using %d threads", omp_get_num_threads());
#endif

  /* Make a private copy of the mpmod_t struct */
  mpmod_init_set (modulus, modulus_param);

  mpz_init (t);
  mpz_init (t1);
  mpres_init (r[0], modulus);
  mpres_init (r[1], modulus);
  mpres_init (r[2], modulus);
  mpres_init (x_0, modulus);
  mpres_init (x_Mi, modulus);

    if (ntt_files != NULL && g_ntt == NULL)
      {
        /* Make an mpzspv buffer from which we write to files */
        tmp_ntt = mpzspv_init (buflen, ntt_context);
        if (tmp_ntt == NULL)
          {
            fprintf (stderr, "pm1_sequence_g(): error, could not initialise tmp_ntt (out of memory?)\n");
            abort();
          }
      }

  if (want_output)
    {
      if (test_verbose (OUTPUT_TRACE))
	{ 
	  mpres_get_z (t, b_1, modulus);
	  outputf (OUTPUT_TRACE, "\n/* pm1_sequence_g */ N = %Zd; "
		   "b_1 = Mod(%Zd, N); /* PARI */\n", modulus->orig_modulus, t);
	  outputf (OUTPUT_TRACE, "/* pm1_sequence_g */ P = %" PRIu64
                   "; M = %" PRIu64 , P, M);
	  outputf (OUTPUT_TRACE, "; m_1 = %Zd; /* PARI */\n", m_1);
	  outputf (OUTPUT_TRACE,"/* pm1_sequence_g */ r = b_1^P; /* PARI */\n");
	  outputf (OUTPUT_TRACE, "/* pm1_sequence_g */ x_0 = "
		   "b_1^(2*% " PRId64 " + (2*m_1 + 1)*P); /* PARI */\n", k_2);
	}
    }

  /* We use (M-(i+1))^2 = (M-i)^2 + 2*(-M+i) + 1 */
  mpz_set_uint64 (t, P);
  mpres_pow (r[0], b_1, t, modulus);     /* r[0] = b_1^P = r */
  if (test_verbose (OUTPUT_TRACE))
    {
      mpres_get_z (t, r[0], modulus);
      outputf (OUTPUT_TRACE, "/* pm1_sequence_g */ r == %Zd /* PARI C */\n", t);
    }
  
  /* FIXME: This is a huge mess, clean up some time */

  mpz_set_uint64 (t, M);
  mpz_neg (t, t);
  mpz_mul_2exp (t, t, 1UL);
  mpz_add_ui (t, t, 1UL);
  mpres_pow (r[1], r[0], t, modulus);    /* r[1] = r^{2(-M+i)+1}, i = 0 */
  mpz_set_uint64 (t, M);
  mpz_mul (t, t, t);                     /* t = M^2 */
  mpres_pow (r[2], r[0], t, modulus);    /* r[2] = r^{(M-i)^2}, i = 0 */
  mpres_mul (r[0], r[0], r[0], modulus); /* r[0] = r^2 */

  mpz_mul_2exp (t, m_1, 1UL);
  mpz_add_ui (t, t, 1UL);
  mpz_set_uint64(t1, P);
  mpz_mul (t, t, t1);
  mpz_set_int64(t1, k_2);
  mpz_addmul_ui (t, t1, 2UL);
  if (want_output)
    {
      outputf (OUTPUT_TRACE, "/* pm1_sequence_g */ 2*%" PRId64 , k_2);
      outputf (OUTPUT_TRACE, " + (2*%Zd + 1)*P == %Zd /* PARI C */\n", m_1, t);
    }

  mpres_pow (x_0, b_1, t, modulus);  /* x_0 = b_1^{2*k_2 + (2*m_1 + 1)*P} */
  if (want_output && test_verbose (OUTPUT_TRACE))
    {
      mpres_get_z (t, x_0, modulus);
      outputf (OUTPUT_TRACE, "/* pm1_sequence_g */ x_0 == %Zd /* PARI C */\n", 
	       t);
    }
  
  mpz_set_uint64 (t, M);
  mpres_pow (x_Mi, x_0, t, modulus); /* x_Mi = x_0^{M-i}, i = 0 */

  mpres_invert (x_0, x_0, modulus);  /* x_0 := x_0^{-1} now */
  mpres_mul (r[1], r[1], x_0, modulus); /* r[1] = x_0^{-1} * r^{-2M+1} */
  
  mpres_mul (r[2], r[2], x_Mi, modulus); /* r[2] = x_0^M * r^{M^2} */
  mpres_get_z (t, r[2], modulus);
  outputf (OUTPUT_TRACE, "/* pm1_sequence_g */ g_%lu = %Zd; /* PARI */\n", 
	   offset, t);
  if (g_mpz != NULL)
    mpz_set (g_mpz[offset], t);
  if (g_ntt != NULL)
    mpzspv_from_mpzv (g_ntt, offset, &t, 1UL, ntt_context);
  if (tmp_ntt != NULL)
    mpzspv_from_mpzv (tmp_ntt, 0, &t, 1UL, ntt_context);

  /* So here we have for i = 0
     r[2] = x_0^(M-i) * r^{(M-i)^2}
     r[1] = x_0^{-1} * r^{2(-M+i)+1}
     r[0] = r^2
     t = r[2]
  */

  for (i = 0; i < l; )
    {
      const unsigned long len_now = MIN(l - i, buflen);
      unsigned long j;
      for (j = (i == 0) ? 1 : 0; j < len_now; j++)
        {
          if (g_mpz != NULL)
            {
              mpres_mul_z_to_z (g_mpz[offset + i + j], r[1], g_mpz[offset + i + j - 1], 
                                modulus);
              outputf (OUTPUT_TRACE, "/* pm1_sequence_g */ g_%lu = %Zd;"
                       " /* PARI */\n", offset + i, g_mpz[offset + i]);
            }
          if (g_ntt != NULL || tmp_ntt != NULL)
            {
              mpres_mul_z_to_z (t, r[1], t, modulus);
              if (g_mpz == NULL) /* Only one should be non-NULL... */
                  outputf (OUTPUT_TRACE, "/* pm1_sequence_g */ g_%lu = %Zd;"
                           " /* PARI */\n", offset + i, t);
              if (g_ntt != NULL)
                mpzspv_from_mpzv (g_ntt, offset + i + j, &t, 1UL, ntt_context);
              if (tmp_ntt != NULL)
                mpzspv_from_mpzv (tmp_ntt, j, &t, 1UL, ntt_context);
            }
          mpres_mul (r[1], r[1], r[0], modulus);
        }
      if (ntt_files != NULL)
        {
          if (tmp_ntt != NULL)
            mpzspv_write (tmp_ntt, 0, ntt_files, offset + i, len_now, ntt_context);
          else
            mpzspv_write (g_ntt, offset + i, ntt_files, offset + i, len_now, ntt_context);
        }
      i += len_now;
    }

  if (g_ntt != NULL)
    mpzspv_verify(g_ntt, offset, l, ntt_context);

  mpres_clear (r[0], modulus);
  mpres_clear (r[1], modulus);
  mpres_clear (r[2], modulus);
  mpres_clear (x_0, modulus);
  mpres_clear (x_Mi, modulus);
  mpz_clear (t);
  mpz_clear (t1);
  mpmod_clear (modulus); /* Clear our private copy of modulus */
  if (ntt_files != NULL && g_ntt == NULL)
    {
      mpzspv_clear (tmp_ntt, ntt_context);
    }
  }

  print_elapsed_time (OUTPUT_VERBOSE, timestart, realstart);
  
  if (test_verbose (OUTPUT_TRACE))
    {
      for (i = 0; i < l_param; i++)
	{
	  outputf (OUTPUT_TRACE, "/* pm1_sequence_g */ g_%" PRIu64
                   " == x_0^(M - %" PRIu64 ") * r^((M - %" PRIu64 ")^2) "
                   "/* PARI C */\n", i, i, i);
	}
      outputf (OUTPUT_TRACE, "/* pm1_sequence_g */ g(x) = g_0");
      for (i = 1; i < l_param; i++)
	outputf (OUTPUT_TRACE, " + g_%" PRIu64 " * x^%" PRIu64 ,  i, i);
      outputf (OUTPUT_TRACE, " /* PARI */\n");
    }
}


/* Compute h_j = r^(-j^2) * f_j for 0 <= j < d as described in section 9 
   of the paper. h == f is ok. */

static void 
pm1_sequence_h (listz_t h, mpzspv_t h_ntt, FILE **ntt_files, mpz_t *f, 
                const mpres_t r, const uint64_t d, mpmod_t modulus_parm, 
		const mpzspm_t ntt_context)
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
  timestart = cputime ();
  realstart = realtime ();

#ifdef _OPENMP
#pragma omp parallel if (d > 100)
#endif
  {
    mpres_t fd[3]; /* finite differences table for r^{-i^2}*/
    mpz_t t;       /* the h_j value as an mpz_t */
    uint64_t i, offset, len;
    const uint64_t blocklen = 16384;
    mpzspv_t tmp_ntt = NULL;
    mpmod_t modulus;

    /* Adjust offset and length for this thread */
    get_chunk (&offset, &len, d);
#ifdef _OPENMP
    if (omp_get_thread_num() == 0)
      outputf (OUTPUT_VERBOSE, " using %d threads", omp_get_num_threads());
#endif
    
    mpmod_init_set (modulus, modulus_parm);
    mpres_init (fd[0], modulus);
    mpres_init (fd[1], modulus);
    mpres_init (fd[2], modulus);
    mpz_init (t);
    if (ntt_files != NULL && h_ntt == NULL)
      {
        /* Make an mpzspv buffer from which we write to files */
        tmp_ntt = mpzspv_init (blocklen, ntt_context);
        if (tmp_ntt == NULL)
          {
            fprintf (stderr, "pm1_sequence_h(): error, could not initialise tmp_ntt (out of memory?)\n");
            abort();
          }
      }
    
    /* We have (n + 1)^2 = n^2 + 2n + 1. For the finite differences we'll 
       need r^{-2}, r^{-(2n+1)}, r^{-n^2}. Init for n = 0. */
    
    /* r^{-2} in fd[0] is constant and could be shared. Computing it 
       separately in each thread has the advantage of putting it in
       local memory. May not make much difference overall */

    mpres_mul (fd[0], invr, invr, modulus); /* fd[0] = r^{-2} */
    mpz_set_uint64 (t, offset);
    mpz_mul_2exp (t, t, 1UL);
    mpz_add_ui (t, t, 1UL);                 /* t = 2 * offset + 1 */
    mpres_pow (fd[1], invr, t, modulus);    /* fd[1] = r^{-(2*offset+1)} */
    mpz_set_uint64 (t, offset);
    mpz_mul (t, t, t);                      /* t = offset^2 */
    mpres_pow (fd[2], invr, t, modulus);    /* fd[2] = r^{-offset^2} */
    
    /* Generate the sequence */
    for (i = 0; i < len; )
      {
        const uint64_t len_now = MIN(len - i, blocklen);
        uint64_t j;
        for (j = 0; j < len_now; j++)
          {
            mpres_mul_z_to_z (t, fd[2], f[offset + i + j], modulus);
            outputf (OUTPUT_TRACE, 
                     "/* pm1_sequence_h */ h_%lu = %Zd; /* PARI */\n", offset + i + j, t);
            
            if (h != NULL)
              mpz_set (h[offset + i + j], t);
            if (h_ntt != NULL)
              mpzspv_from_mpzv (h_ntt, offset + i + j, &t, 1UL, ntt_context);
            if (tmp_ntt != NULL)
              mpzspv_from_mpzv (tmp_ntt, j, &t, 1UL, ntt_context);
            
            mpres_mul (fd[2], fd[2], fd[1], modulus); /* fd[2] = r^{-j^2} */
            mpres_mul (fd[1], fd[1], fd[0], modulus); /* fd[1] = r^{-2*j-1} */
          }
        if (ntt_files != NULL)
          {
            if (tmp_ntt != NULL)
              mpzspv_write (tmp_ntt, 0,        ntt_files, offset + i, len_now, ntt_context);
            else
              mpzspv_write (h_ntt, offset + i, ntt_files, offset + i, len_now, ntt_context);
          }
        i += len_now;
      }    
    
    mpres_clear (fd[2], modulus);
    mpres_clear (fd[1], modulus);
    mpres_clear (fd[0], modulus);
    mpz_clear (t);
    if (tmp_ntt != NULL)
      {
        mpzspv_clear (tmp_ntt, ntt_context);
        tmp_ntt = NULL;
      }
    mpmod_clear (modulus);
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

  sets_get_factored_sorted (S_1, params->P);

  {
    mpz_t t1, t2;
    
    mpz_init (t1);
    mpz_init (t2);
    sets_sumset_minmax (t1, S_1, 1);
    sets_max (t2, params->P);
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


ATTRIBUTE_UNUSED
static mpzspv_t *
mpzspv_init_mt (spv_size_t len, mpzspm_t mpzspm)
{
  int i; /* OpenMP wants the iteration variable a signed type */
  mpzspv_t *x = (mpzspv_t *) malloc (mpzspm->sp_num * sizeof (spv_t *));
  
  if (x == NULL)
    return NULL;

  for (i = 0; i < (int) mpzspm->sp_num; i++)
    x[i] = NULL;
  
#ifdef _OPENMP
#pragma omp parallel private(i) shared(x)
  {
#pragma omp for
#endif
    for (i = 0; i < (int) mpzspm->sp_num; i++)
      x[i] = (spv_t *) sp_aligned_malloc (len * sizeof (sp_t));
	
#ifdef _OPENMP
  }
#endif

  for (i = 0; i < (int) mpzspm->sp_num; i++)
    if (x[i] == NULL)
      break;

  if (i != (int) mpzspm->sp_num) /* There is a NULL pointer */
    {
      for (i = 0; i < (int) mpzspm->sp_num; i++)
	if (x[i] != NULL)
	  sp_aligned_free(x[i]);
      return NULL;
    }

#if 0
  if (test_verbose (OUTPUT_DEVVERBOSE))
    {
      spv_t * last = x[0];
      printf ("mpzspv_init_mt: x[0] = %p\n", x[0]);
      for (i = 1; i < (int) mpzspm->sp_num; i++)
        printf ("mpzspv_init_mt: x[%d] = %p, distance = %ld\n", 
                i, x[i], (long) (x[i] - x[i-1]));
    }
#endif

  return x;
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
ntt_sqr_reciprocal (mpzv_t R, const mpzv_t S, mpzspv_t dft, FILE **dft_files, 
		    const spv_size_t n, const mpzspm_t ntt_context)
{
#ifdef WANT_ASSERT
  mpz_t S_eval_1, R_eval_1;
#endif
  
  if (n == 0)
    return;

  if (n == 1)
    {
      mpz_mul (R[0], S[0], S[0]);
      mpz_mod (R[0], R[0], ntt_context->modulus);
      return;
    }

#ifdef WANT_ASSERT
  mpz_init (S_eval_1);
  list_recip_eval1 (S_eval_1, S, n);
  /* Compute (S(1))^2 */
  mpz_mul (S_eval_1, S_eval_1, S_eval_1);
  mpz_mod (S_eval_1, S_eval_1, ntt_context->modulus);
#endif

  if (dft_files == NULL)
    {
      /* Fill NTT elements [0 .. n-1] with coefficients */
      mpzspv_from_mpzv (dft, (spv_size_t) 0, S, n, ntt_context);
      mpzspv_sqr_reciprocal (dft, n, ntt_context);
    }
  else
    {
      const spv_size_t blocklen = 16384;
      mpzspv_t mpzspv = mpzspv_init (blocklen, ntt_context);
      if (mpzspv == NULL)
        {
          abort();
        }
      mpzspv_from_mpzv_file (mpzspv, (spv_size_t) 0, dft_files, S, NULL, n, 
                             blocklen, ntt_context);
      mpzspv_sqr_reciprocal_file (dft_files, n, ntt_context);
      mpzspv_clear (mpzspv, ntt_context);
    }
  
#if defined(_OPENMP)
#pragma omp parallel if (n > 50)
#endif
  {
    spv_size_t i, offset, chunklen;

    get_chunk (&offset, &chunklen, 2*n - 1);
        
    if (dft_files == NULL)
      {
        mpzspv_to_mpzv (dft, offset, R + offset, chunklen, ntt_context);
      } else {
        const spv_size_t blocklen = 16384;
        mpzspv_t mpzspv = mpzspv_init (blocklen, ntt_context);
        if (mpzspv == NULL)
          {
            abort();
          }
        mpzspv_to_mpzv_file (mpzspv, offset, dft_files, R + offset, NULL, 
                             chunklen, blocklen, ntt_context);
        mpzspv_clear (mpzspv, ntt_context);
      }
    for (i = offset; i < offset + chunklen; i++)
      mpz_mod (R[i], R[i], ntt_context->modulus);
  }

#ifdef WANT_ASSERT
  mpz_init (R_eval_1);
  /* Compute (S^2)(1) and compare to (S(1))^2 */
  list_recip_eval1 (R_eval_1, R, 2 * n - 1);
  mpz_mod (R_eval_1, R_eval_1, ntt_context->modulus);
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


/* Computes gcd(\prod_{0 <= i < len} (ntt[i + offset] + add[i]), N), 
   the NTT residues are converted to integer residues (mod N) first.
   If add == NULL, add[i] is assumed to be 0. */

static void
ntt_gcd (mpz_t f, mpz_t *product, mpzspv_t ntt, FILE **ntt_files, 
         const uint64_t ntt_offset, 
	 const listz_t add, const uint64_t len_param, 
	 const mpzspm_t ntt_context, mpmod_t modulus_param)
{
  uint64_t i, j;
  const unsigned long Rlen = MPZSPV_NORMALISE_STRIDE;
  listz_t R;
  uint64_t len, thread_offset;
  mpres_t tmpres, tmpprod, totalprod;
  mpmod_t modulus;
  long timestart, realstart;
  mpzspv_t tmpspv = NULL;
  
  outputf (OUTPUT_VERBOSE, "Computing gcd of coefficients and N");
  timestart = cputime ();
  realstart = realtime ();

  /* All the threads will multiply their partial products to this one. */
  mpres_init (totalprod, modulus_param);
  mpres_set_ui (totalprod, 1UL, modulus_param);

#ifdef _OPENMP
#pragma omp parallel if (len_param > 100) private(i, j, R, len, thread_offset, tmpres, tmpprod, modulus) shared(totalprod)
  {
#pragma omp master
    {
      outputf (OUTPUT_VERBOSE, " using %d threads", omp_get_num_threads());
    }
#endif

    get_chunk (&thread_offset, &len, len_param);
    /* Make a private copy of the mpmod_t struct */
    mpmod_init_set (modulus, modulus_param);

    MEMORY_TAG;
    R = init_list2 (Rlen, (mpz_size (modulus->orig_modulus) + 2) * 
                           GMP_NUMB_BITS);
    MEMORY_UNTAG;
    if (ntt_files != NULL)
      {
        MEMORY_TAG;
        tmpspv = mpzspv_init (Rlen, ntt_context);
        MEMORY_UNTAG;
      }
    mpres_init (tmpres, modulus);
    mpres_init (tmpprod, modulus);
    mpres_set_ui (tmpprod, 1UL, modulus);
    
    for (i = 0; i < len; i += Rlen)
      {
	const unsigned long blocklen = MIN(len - i, Rlen);

	/* Convert blocklen residues from NTT to integer representatives
	   and store them in R */
        if (ntt_files == NULL)
          {
            mpzspv_to_mpzv (ntt, ntt_offset + thread_offset + i, R, blocklen, 
                            ntt_context);
          } else {
#ifdef _OPENMP
#pragma omp critical
#endif
            {
              mpzspv_read (tmpspv, 0, ntt_files, ntt_offset + thread_offset + i, 
                           blocklen, ntt_context);
            }
            mpzspv_to_mpzv (tmpspv, 0, R, blocklen, ntt_context);
          }

	/* Accumulate product in tmpprod */
	for (j = 0; j < blocklen; j++)
	  {
	    outputf (OUTPUT_TRACE, "r_%" PRIu64 " = %Zd; /* PARI */\n", 
	             i + thread_offset + j, R[j]);
	    if (add != NULL)
	      mpz_add (R[j], R[j], add[i + thread_offset + j]);
	    mpres_set_z_for_gcd (tmpres, R[j], modulus);
#define TEST_ZERO_RESULT
#ifdef TEST_ZERO_RESULT
	    if (mpres_is_zero (tmpres, modulus))
	      outputf (OUTPUT_VERBOSE, "R_[%" PRIu64 "] = 0\n", i + thread_offset + j);
#endif
	    mpres_mul (tmpprod, tmpprod, tmpres, modulus); 
	  }
      }
#ifdef _OPENMP
#pragma omp critical
    {
      mpres_mul (totalprod, totalprod, tmpprod, modulus);
    }
#else
    mpres_set (totalprod, tmpprod, modulus);
#endif
    if (ntt_files != NULL)
      {
        MEMORY_TAG;
        mpzspv_clear (tmpspv, ntt_context);
        MEMORY_UNTAG;
      }
    mpres_clear (tmpres, modulus);
    mpres_clear (tmpprod, modulus);
    mpmod_clear (modulus);
    clear_list (R, Rlen);
#ifdef _OPENMP
  }
#endif

  if (product != NULL)
    mpres_get_z (*product, totalprod, modulus_param);

  mpres_gcd (f, totalprod, modulus_param);
  mpres_clear (totalprod, modulus_param);

  print_elapsed_time (OUTPUT_VERBOSE, timestart, realstart);
}


int 
pm1fs2 (mpz_t f, const mpres_t X, mpmod_t modulus, 
	const faststage2_param_t *params)
{
  uint64_t nr;
  uint64_t i, l, lenF, lenG, lenR, tmplen;
  set_list_t S_1; /* This is stored as a set of sets (arithmetic 
                     progressions of prime length */
  int64_t *s2_sumset; /* set of sums of S_2 */
  uint64_t s2_sumset_size;
  listz_t F;   /* Polynomial F has roots X^{k_1} for k_1 \in S_1, so has 
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
  lenF = params->s_1 / 2 + 1 + 1; /* Another +1 because poly_from_sets_V stores
				     the leading 1 monomial for each factor */
  F = init_list2 (lenF, (unsigned int) abs (modulus->bits));
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
  
  i = poly_from_sets_V (F, mr, &S_1, tmp, tmplen, modulus, NULL, NULL, NULL);
  ASSERT_ALWAYS(2 * i == params->s_1);
  ASSERT(mpz_cmp_ui (F[i], 1UL) == 0);
  sets_free(&S_1);
  
  outputf (OUTPUT_VERBOSE, " took %lums\n", cputime () - timestart);
  if (test_verbose (OUTPUT_TRACE))
    {
      for (i = 0; i < params->s_1 / 2 + 1; i++)
	outputf (OUTPUT_TRACE, "f_%lu = %Zd; /* PARI */\n", i, F[i]);
      outputf (OUTPUT_TRACE, "f(x) = f_0");
      for (i = 1; i < params->s_1 / 2 + 1; i++)
	outputf (OUTPUT_TRACE, "+ f_%lu * (x^%lu + x^(-%lu))", i, i, i);
      outputf (OUTPUT_TRACE, "/* PARI */ \n");
    }
  
  mpz_set_ui (mt, params->P);
  mpres_pow (mr, X, mt, modulus); /* mr = X^P */
  pm1_sequence_h (F, NULL, NULL, F, mr, params->s_1 / 2 + 1, modulus, NULL); 

  /* Make a symmetric copy of F in h. It will have length 
     s_1 + 1 = 2*lenF - 1 */
  /* I.e. with F = [3, 2, 1], s_1 = 4, we want h = [1, 2, 3, 2, 1] */
  for (i = 0; i < params->s_1 / 2 + 1; i++)
    *(h[i]) = *(F[params->s_1 / 2 - i]); /* Clone the mpz_t. */
  for (i = 0; i < params->s_1 / 2; i++)
    *(h[i + params->s_1 / 2 + 1]) = *(F[i + 1]);
  if (test_verbose (OUTPUT_TRACE))
    {
      for (i = 0; i < params->s_1 + 1; i++)
        outputf (OUTPUT_VERBOSE, "h_%" PRIu64 " = %Zd; /* PARI */\n", i, h[i]);
      outputf (OUTPUT_VERBOSE, "h(x) = h_0");
      for (i = 1; i < params->s_1 + 1; i++)
        outputf (OUTPUT_VERBOSE, " + h_%" PRIu64 "* x^%" PRIu64 , i, i);
      outputf (OUTPUT_VERBOSE, " /* PARI */\n");
    }

  for (l = 0; l < params->s_2; l++)
    {
      const uint64_t M = params->l - 1L - params->s_1 / 2L;
      outputf (OUTPUT_VERBOSE, "Multi-point evaluation %" PRIu64
                        " of %" PRIu64 ":\n", l + 1, params->s_2);
      pm1_sequence_g (g, NULL, NULL, X, params->P, M, params->l, 
		      params->m_1, s2_sumset[l], modulus, NULL);

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
        mpres_get_z (tmp[1], tmpprod, modulus); /* For printing */
	mpres_gcd (tmp[0], tmpprod, modulus);
	mpres_clear (tmpprod, modulus);
	mpres_clear (tmpres, modulus);
      }

      outputf (OUTPUT_VERBOSE, " took %lums\n", cputime () - timestart);
      outputf (OUTPUT_RESVERBOSE, "Product of R[i] = %Zd (times some "
	       "power of 2 if REDC was used! Try -mpzmod)\n", tmp[1]);

      if (mpz_cmp_ui (tmp[0], 1UL) > 0)
	{
	  mpz_set (f, tmp[0]);
	  youpi = ECM_FACTOR_FOUND_STEP2;
	  break;
	}
    }

#ifdef SHOW_TMP_USAGE
  for (i = tmplen - 1; i > 0; i--)
    if (tmp[i]->_mp_alloc > 1)
      break;
  outputf (OUTPUT_DEVVERBOSE, "Highest used temp element is tmp[%lu]\n", i);
#endif
  
  free (s2_sumset);
  free (h);
  clear_list (F, lenF);
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


int 
pm1fs2_ntt (mpz_t f, const mpres_t X, mpmod_t modulus, 
	const faststage2_param_t *params)
{
  uint64_t nr;
  uint64_t l, lenF;
  set_list_t S_1; /* This is stored as a set of sets (arithmetic 
                       progressions of prime length */
  int64_t *s2_sumset; /* set of sums of S_2 */ uint64_t s2_sumset_size;
  listz_t F;   /* Polynomial F has roots X^{k_1} for k_1 \in S_1, so has 
		  degree s_1. It is symmetric, so has only s_1 / 2 + 1 
		  distinct coefficients. The sequence h_j will be stored in 
		  the same memory and won't be a monic polynomial, so the 
		  leading 1 monomial of F will be stored explicitly. Hence we 
		  need s_1 / 2 + 1 entries. */
  mpzspm_t ntt_context;
  mpzspv_t g_ntt, h_ntt;
  FILE **g_files = NULL, **h_files = NULL;
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
      char *filename = 
        malloc ((strlen(params->file_stem) + 3) * sizeof (char));
      if (filename == NULL)
        {
          fprintf (stderr, 
                   "pm1fs2_ntt(): could not allocate memory for filename\n");
          mpzspm_clear (ntt_context);
          return ECM_ERROR;
        }
      sprintf (filename, "%s.g", params->file_stem);
      g_files = mpzspv_open_fileset (filename, ntt_context);
      sprintf (filename, "%s.h", params->file_stem);
      h_files = mpzspv_open_fileset (filename, ntt_context);
      free (filename);
      filename = NULL;
    }

  /* Allocate all the memory we'll need for building f */
  lenF = params->s_1 / 2 + 1 + 1; /* Another +1 because poly_from_sets_V stores
				     the leading 1 monomial for each factor */
  F = init_list2 (lenF, (unsigned int) abs (modulus->bits));

  /* Compute Q = X + 1/X */
  mpres_init (Q, modulus);
  mpres_invert (Q, X, modulus);
  mpres_add (Q, Q, X, modulus);

  /* Compute XP = X^P */
  mpres_init (XP, modulus);
  mpz_init (mt);
  mpz_set_ui (mt, params->P);
  mpres_pow (XP, X, mt, modulus);

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

  if (build_F_ntt (F, Q, &S_1, params, params->file_stem, modulus) == ECM_ERROR)
    {
      sets_free (&S_1);
      free (s2_sumset);
      mpz_clear (mt);
      mpres_clear (Q, modulus);
      mpres_clear (XP, modulus);
      mpzspm_clear (ntt_context);
      clear_list (F, lenF);
      return ECM_ERROR;
    }

  sets_free (&S_1);

  mpres_clear (Q, modulus);
  
  if (params->file_stem == NULL)
    h_ntt = mpzspv_init (params->l / 2 + 1, ntt_context);
  else
    h_ntt = NULL;

  pm1_sequence_h (NULL, h_ntt, h_files, F, XP, params->s_1 / 2 + 1, modulus, 
		  ntt_context);

  clear_list (F, lenF);
  mpres_clear (XP, modulus);

  if (params->file_stem == NULL)
    {
      g_ntt = mpzspv_init (params->l, ntt_context);
    } else {
      if (h_ntt != NULL) /* From debugging tests */
        {
          mpzspv_write (h_ntt, 0, h_files, 0, params->l / 2 + 1, ntt_context);
          mpzspv_clear (h_ntt, ntt_context);
          h_ntt = NULL;
        }
      g_ntt = NULL;
    }

  /* Compute the DCT-I of h */
  outputf (OUTPUT_VERBOSE, "Computing DCT-I of h");
#ifdef _OPENMP
  outputf (OUTPUT_VERBOSE, " using %d threads", omp_get_thread_limit());
#endif
  if (test_verbose (OUTPUT_TRACE))
    {
      if (h_ntt != NULL)
        mpzspv_print (h_ntt, 0, params->s_1 / 2 + 1, "h_ntt", ntt_context);
      else
        mpzspv_print_file (h_files, 0, params->s_1 / 2 + 1, "h_ntt", ntt_context);
    }

  timestart = cputime ();
  realstart = realtime ();
  
  mpzspv_to_dct1_file (h_ntt, h_ntt, h_files, params->s_1 / 2 + 1, 
                       params->l / 2 + 1, ntt_context);

  print_elapsed_time (OUTPUT_VERBOSE, timestart, realstart);
  
  if (test_verbose (OUTPUT_TRACE))
    {
      if (h_ntt != NULL)
        mpzspv_print (h_ntt, 0, params->s_1 / 2 + 1, "DCT-I of h_ntt", ntt_context);
      else
        mpzspv_print_file (h_files, 0, params->s_1 / 2 + 1, "DCT-I of h_ntt", ntt_context);
    }

  if (test_verbose (OUTPUT_RESVERBOSE))
    {
      mpz_init (product);
      product_ptr = &product;
    }

  for (l = 0; l < params->s_2; l++)
    {
      const uint64_t M = params->l - 1L - params->s_1 / 2L;

      outputf (OUTPUT_VERBOSE, "Multi-point evaluation %" PRIu64 
                        " of %" PRIu64 ":\n", l + 1, params->s_2);
      /* Compute the coefficients of the polynomial g(x) */
      pm1_sequence_g (NULL, g_ntt, g_files, X, params->P, M, params->l, 
		      params->m_1, s2_sumset[l], modulus, ntt_context);

      if (test_verbose (OUTPUT_TRACE))
        {
          if (g_ntt != NULL)
            mpzspv_print (g_ntt, 0, params->s_1 / 2 + 1, "g_ntt", ntt_context);
          else
            mpzspv_print_file (g_files, 0, params->s_1 / 2 + 1, "g_ntt", ntt_context);
        }

      /* Do the convolution */
      outputf (OUTPUT_VERBOSE, "Computing g*h");
#ifdef _OPENMP
      outputf (OUTPUT_VERBOSE, " using %d threads", omp_get_thread_limit());
#endif
      timestart = cputime ();
      realstart = realtime ();
      mpzspv_mul_ntt_file (g_ntt, 0, g_files, 
          g_ntt, 0, params->l, g_files, 
          h_ntt, 0, params->l / 2 + 1, h_files,
          params->l, 0, 0, ntt_context, 
          NTT_MUL_STEP_FFT1 + NTT_MUL_STEP_MULDCT + NTT_MUL_STEP_IFFT);
      print_elapsed_time (OUTPUT_VERBOSE, timestart, realstart);
      
      if (test_verbose (OUTPUT_TRACE))
        {
          if (g_ntt != NULL)
            mpzspv_print (g_ntt, 0, params->s_1 / 2 + 1, "g_ntt * h_ntt", ntt_context);
          else
            mpzspv_print_file (g_files, 0, params->s_1 / 2 + 1, "g_ntt * h_ntt", ntt_context);
        }

      /* Compute GCD of N and coefficients of product polynomial */
      ntt_gcd (mt, product_ptr, g_ntt, g_files, params->s_1 / 2, NULL, nr, ntt_context, 
	       modulus);

      outputf (OUTPUT_RESVERBOSE, "Product of R[i] = %Zd (times some "
	       "power of 2 if REDC was used! Try -mpzmod)\n", product);

      /* If we found a factor, stop */
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
  if (params->file_stem == NULL)
    {
      mpzspv_clear (g_ntt, ntt_context);
      mpzspv_clear (h_ntt, ntt_context);
    }
  if (params->file_stem != NULL)
    {
      mpzspv_close_fileset (g_files, ntt_context);
      mpzspv_close_fileset (h_files, ntt_context);
    }
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


static void 
gfp_ext_print (const mpres_t r_x, const mpres_t r_y, mpmod_t modulus, 
	       const int verbose)
{
  mpz_t t1, t2;

  if (!test_verbose (verbose))
    return;

  mpz_init (t1);
  mpz_init (t2);
  mpres_get_z (t1, r_x, modulus);
  mpres_get_z (t2, r_y, modulus);
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
gfp_ext_mul (mpres_t r_0, mpres_t r_1, const mpres_t a_0, const mpres_t a_1,
	     const mpres_t b_0, const mpres_t b_1, const mpres_t Delta, 
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
      gfp_ext_print (a_0, a_1, modulus, OUTPUT_TRACE);
      outputf (OUTPUT_TRACE, ") * (");
      gfp_ext_print (b_0, b_1, modulus, OUTPUT_TRACE);
    }
  
  mpres_add (tmp[0], a_0, a_1, modulus);
  mpres_add (tmp[1], b_0, b_1, modulus);
  mpres_mul (tmp[1], tmp[0], tmp[1], modulus); /* t[1] = (a_0+a_1)*(b_0+b_1) = 
					    a_0*b_0 + a_0*b_1 + a_1*b_0 + 
					    a_1*b_1 */

  mpres_mul (r_0, a_0, b_0, modulus);    /* r_0 = a_0*b_0. We don't need a_0 
					    or b_0 any more now */
  mpres_sub (tmp[1], tmp[1], r_0, modulus);  /* t[1] = a_0*b_1 + a_1*b_0 + 
						a_1*b_1 */
  
  mpres_mul (tmp[0], a_1, b_1, modulus);   /* t[0] = a_1*b_1. We don't need 
					      a_1 or b_1 any more now */
  mpres_sub (r_1, tmp[1], tmp[0], modulus);  /* r_1 == a_0*b_1 + a_1*b_0 */
  
  mpres_mul (tmp[0], tmp[0], Delta, modulus); /* t[0] = a_1*b_1*Delta */
  mpres_add (r_0, r_0, tmp[0], modulus);   /* r_0 = a_0*b_0 + a_1*b_1*Delta */

  if (0 && test_verbose (OUTPUT_TRACE))
    {
      outputf (OUTPUT_TRACE, ") == ");
      gfp_ext_print (r_0, r_1, modulus, OUTPUT_TRACE);
      outputf (OUTPUT_TRACE, " /* PARI C */\n");
    }
}


/* Computes (a_0 + a_1 * sqrt(Delta))^2, where the norm 
   (a_0^2 - a_1^2*Delta) is assumed to be equal to 1. Hence 
   (a_0 + a_1 * sqrt(Delta))^2 = a_0^2 + 2*a_0*a_1*sqrt(Delta) + a_1^2*Delta
   and a_0^2 + a_1^2*Delta = a_0^2 + a_1^2*Delta + norm - 1 = 2*a_0^2 - 1.
   a_0 and r_0, as well as a_1 and r_1 may overlap */

static void
gfp_ext_sqr_norm1 (mpres_t r_0, mpres_t r_1, const mpres_t a_0, 
		   const mpres_t a_1, mpmod_t modulus)
{
  ASSERT (a_0 != r_1);  /* a_0 is read after r_1 is written */
  
  if (pari)
    gmp_printf ("/* gfp_ext_sqr_norm1 */ (%Zd + %Zd * w)^2 %% N == ", a_0, a_1);
  
  mpres_mul (r_1, a_0, a_1, modulus);
  mpres_add (r_1, r_1, r_1, modulus);       /* r_1 = 2*a_0*a_1 */
  
  mpres_mul (r_0, a_0, a_0, modulus);
  mpres_add (r_0, r_0, r_0, modulus);
  mpres_sub_ui (r_0, r_0, 1UL, modulus);    /* r_0 = 2*a_0^2 - 1 */

  if (pari)
    gmp_printf ("(%Zd + %Zd * w) %% N /* PARI C */\n", r_0, r_1);
}


/* Raise (a0 + a1*sqrt(Delta)) to the power e which is a 64-bit signed int.
   (a0 + a1*sqrt(Delta)) is assumed to have norm 1, i.e. 
   a0^2 - a1^2*Delta == 1. The result is (r0 * r1*sqrt(Delta)). 
   a0, a1, r0 and r1 must not overlap */

static void 
gfp_ext_pow_norm1_sl (mpres_t r0, mpres_t r1, const mpres_t a0, 
                      const mpres_t a1, const int64_t e, const mpres_t Delta, 
                      mpmod_t modulus, const uint64_t tmplen, mpres_t *tmp)
{
  const int64_t abs_e = (e > 0) ? e : -e;
  uint64_t mask = (uint64_t)1 << 63;

  ASSERT (a0 != r0 && a1 != r0 && a0 != r1 && a1 != r1);

  if (e == 0)
    {
      mpres_set_ui (r0, 1UL, modulus);
      mpres_set_ui (r1, 0UL, modulus);
      return;
    }

  /* If e < 0, we want 1/(a0 + a1*sqrt(Delta)). By extending with 
     a0 - a1*sqrt(Delta), we get 
     (a0 - a1*sqrt(Delta)) / (a0^2 - a1^2 * Delta), but that denomiator
     is the norm which is known to be 1, so the result is 
     a0 - a1*sqrt(Delta). */

  while ((abs_e & mask) == 0)
    mask >>= 1;

  mpres_set (r0, a0, modulus);
  mpres_set (r1, a1, modulus);

  while (mask > 1)
    {
      gfp_ext_sqr_norm1 (r0, r1, r0, r1, modulus);
      mask >>= 1;
      if (abs_e & mask)
	gfp_ext_mul (r0, r1, r0, r1, a0, a1, Delta, modulus, tmplen, tmp);
    }

  if (e < 0)
    mpres_neg (r1, r1, modulus);

  if (0 && test_verbose (OUTPUT_TRACE))
    {
      mpz_t t;
      mpz_init (t);
      mpres_get_z (t, Delta, modulus);
      outputf (OUTPUT_TRACE, "/* gfp_ext_pow_norm1_sl */ w = quadgen (4*%Zd); "
               "N = %Zd; /* PARI */\n", t, modulus->orig_modulus);
      mpz_clear (t);
      outputf (OUTPUT_TRACE, "/* gfp_ext_pow_norm1_sl */ (");
      gfp_ext_print (a0, a1, modulus, OUTPUT_TRACE);
      outputf (OUTPUT_TRACE, ")^(%" PRId64 ") == ", e);
      gfp_ext_print (r0, r1, modulus, OUTPUT_TRACE);
      outputf (OUTPUT_TRACE, " /* PARI C */\n");
    }
}


/* Same, but taking an mpz_t argument for the exponent */

static void 
gfp_ext_pow_norm1 (mpres_t r0, mpres_t r1, const mpres_t a0, 
                   const mpres_t a1, mpz_t e, const mpres_t Delta, 
                   mpmod_t modulus, const uint64_t tmplen, mpres_t *tmp)
{
  mpz_t abs_e;
  unsigned long idx;

  ASSERT (a0 != r0 && a1 != r0 && a0 != r1 && a1 != r1);

  if (mpz_sgn (e) == 0)
    {
      mpres_set_ui (r0, 1UL, modulus);
      mpres_set_ui (r1, 0UL, modulus);
      return;
    }

  mpz_init (abs_e);
  mpz_abs (abs_e, e);
  idx = mpz_sizeinbase (abs_e, 2) - 1; /* Thus mpz_tstbit (abs_e, idx) == 1 */
  ASSERT (mpz_tstbit (abs_e, idx) == 1);

  mpres_set (r0, a0, modulus);
  mpres_set (r1, a1, modulus);

  while (idx > 0UL)
    {
      gfp_ext_sqr_norm1 (r0, r1, r0, r1, modulus);
      idx--;
      if (mpz_tstbit (abs_e, idx))
	gfp_ext_mul (r0, r1, r0, r1, a0, a1, Delta, modulus, tmplen, tmp);
    }

  if (mpz_sgn (e) < 0)
    mpres_neg (r1, r1, modulus);

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
      gfp_ext_print (a0, a1, modulus, OUTPUT_TRACE);
      outputf (OUTPUT_TRACE, ")^(%Zd) == ", e);
      gfp_ext_print (r0, r1, modulus, OUTPUT_TRACE);
      outputf (OUTPUT_TRACE, " /* PARI C */\n");
    }
}


/* Compute r[i] = a^((k+i)^2) for i = 0, 1, ..., l-1, where "a" is an 
   element of norm 1 in the quadratic extension ring */

ATTRIBUTE_UNUSED static void
gfp_ext_rn2 (mpres_t *r_x, mpres_t *r_y, const mpres_t a_x, const mpres_t a_y,
	     const int64_t k, const uint64_t l, const mpres_t Delta, 
	     mpmod_t modulus, const size_t origtmplen, mpres_t *origtmp)
{
  mpres_t *r2_x = origtmp, *r2_y = origtmp + 2, *v = origtmp + 4, 
    *V2 = origtmp + 6;
  const size_t newtmplen = origtmplen - 7;
  mpres_t *newtmp = origtmp + 7;
  uint64_t i;

  if (l == 0UL)
    return;

  ASSERT (origtmplen >= 8UL);

  if (pari)
    gmp_printf ("/* In gfp_ext_rn2 */ ; a = %Zd + %Zd * w; /* PARI */\n",
		a_x, a_y, modulus->orig_modulus);

  /* Compute r[0] = a^(k^2). We do it by two exponentiations by k and use 
     v[0] and v[1] as temp storage */
  gfp_ext_pow_norm1_sl (v[0], v[1], a_x, a_y, k, Delta, modulus, newtmplen, 
		     newtmp);
  gfp_ext_pow_norm1_sl (r_x[0], r_y[0], v[0], v[1], k, Delta, modulus, 
		     newtmplen, newtmp);
  if (pari)
    gmp_printf ("/* In gfp_ext_rn2 */ a^(%ld^2) %% N == (%Zd + %Zd * w) %% N "
		"/* PARI C */\n", k, r_x[0], r_y[0]);

  /* Compute r[1] = a^((k+1)^2) = a^(k^2 + 2k + 1)*/
  if (l > 1)
    {
      /* v[0] + v[1]*sqrt(Delta) still contains a^k */
      gfp_ext_sqr_norm1 (r_x[1], r_y[1], v[0], v[1], modulus);
      /* Now r[1] = a^(2k) */
      gfp_ext_mul (r_x[1], r_y[1], r_x[1], r_y[1], r_x[0], r_y[0], Delta, 
		   modulus, newtmplen, newtmp);
      /* Now r[1] = a^(k^2 + 2k) */
      gfp_ext_mul (r_x[1], r_y[1], r_x[1], r_y[1], a_x, a_y, Delta, modulus, 
		   newtmplen, newtmp);
      /* Now r[1] = a^(k^2 + 2k + 1) = a^((k+1)^2) */
    }
  if (pari)
    {
      gmp_printf ("/* In gfp_ext_rn2 */ a^(%" PRId64 "^2) %% N == ", k + 1);
      gmp_printf ("(%Zd + %Zd * w) %% N /* PARI C */\n", r_x[1], r_y[1]);
    }
  
  /* Compute r2[0] = a^(k^2+2) = a^(k^2) * a^2 */
  gfp_ext_sqr_norm1 (v[0], v[1], a_x, a_y, modulus);
  gfp_ext_mul (r2_x[0], r2_y[0], r_x[0], r_y[0], v[0], v[1], Delta, modulus, 
	       newtmplen, newtmp);
  if (pari)
    gmp_printf ("/* In gfp_ext_rn2 */ a^(%ld^2+2) %% N == (%Zd + %Zd * w) %% N "
		"/* PARI C */\n", k, r2_x[0], r2_y[0]);
  /* Compute a^((k+1)^2+2) = a^((k+1)^2) * a^2 */
  gfp_ext_mul (r2_x[1], r2_y[1], r_x[1], r_y[1], v[0], v[1], Delta, modulus, 
	       newtmplen, newtmp);
  if (pari) 
    { 
      gmp_printf ("/* In gfp_ext_rn2 */ a^(%" PRId64 "^2+2) %% N == ", k + 1);
      gmp_printf ("(%Zd + %Zd * w) %% N /* PARI C */\n", r2_x[1], r2_y[1]);
    }
  
  /* Compute V_2(a + 1/a). Since 1/a = a_x - a_y, we have a+1/a = 2*a_x.
     V_2(x) = x^2 - 2, so we want 4*a_x^2 - 2. */
  mpres_add (*V2, a_x, a_x, modulus); /* V2 = a + 1/a  = 2*a_x*/
  V (v[0], *V2, 2 * k + 1, modulus);  /* v[0] = V_{2k+1} (a + 1/a) */
  V (v[1], *V2, 2 * k + 3, modulus);  /* v[0] = V_{2k+3} (a + 1/a) */
  mpres_mul (*V2, *V2, *V2, modulus); /* V2 = 4*a_x^2 */
  mpres_sub_ui (*V2, *V2, 2UL, modulus); /* V2 = 4*a_x^2 - 2 */
  if (pari)
    {
      gmp_printf ("/* In gfp_ext_rn2 */ ((a + 1/a)^2 - 2) %% N == "
		  "%Zd %% N /* PARI C */\n", *V2);
      gmp_printf ("/* In gfp_ext_rn2 */ V(%lu, a + 1/a) %% N == %Zd %% N "
		  "/* PARI C */\n", 2 * k + 1, v[0]);
      gmp_printf ("/* In gfp_ext_rn2 */ V(%lu, a + 1/a) %% N == %Zd %% N "
		  "/* PARI C */\n", 2 * k + 3, v[1]);
    }
  
  /* Compute the remaining a^((k+i)^2) values according to Peter's 
     recurrence */
  
  for (i = 2; i < l; i++)
    {
      /* r[i] = r2[i-1] * v[i-2] - r2[i-2], with indices of r2 and i taken
	 modulo 2 */
      mpres_mul (r_x[i], r2_x[1 - i % 2], v[i % 2], modulus);
      mpres_sub (r_x[i], r_x[i], r2_x[i % 2], modulus);
      mpres_mul (r_y[i], r2_y[1 - i % 2], v[i % 2], modulus);
      mpres_sub (r_y[i], r_y[i], r2_y[i % 2], modulus);
      
      /* r2[i] = r2[i-1] * v[i-1] - r[i-2] */
      mpres_mul (r2_x[i % 2], r2_x[1 - i % 2], v[1 - i % 2], modulus);
      mpres_sub (r2_x[i % 2], r2_x[i % 2], r_x[i - 2], modulus);
      mpres_mul (r2_y[i % 2], r2_y[1 - i % 2], v[1 - i % 2], modulus);
      mpres_sub (r2_y[i % 2], r2_y[i % 2], r_y[i - 2], modulus);
      
      /* v[i] = v[i - 1] * V_2(a + 1/a) - v[i - 2] */
      mpres_mul (newtmp[0], v[1 - i % 2], *V2, modulus);
      mpres_sub (v[i % 2], newtmp[0], v[i % 2], modulus);
      if (pari) 
	{
	  gmp_printf ("/* In gfp_ext_rn2 */ V(%" PRId64 ", a + 1/a) ",
			      2 * (k + i) + 1);
	  gmp_printf ("%% N == %Zd %% N /* PARI C */\n", v[i % 2]);
	}
    }
}


/* Compute g_i = x_0^{M-i} * r^{(M-i)^2} for 0 <= i < l. 
   x_0 = b_1^{2*k_2 + (2*m_1 + 1) * P}. r = b_1^P. */

static void
pp1_sequence_g (listz_t g_x, listz_t g_y, mpzspv_t g_x_ntt, mpzspv_t g_y_ntt,
		const mpres_t b1_x, const mpres_t b1_y, const uint64_t P, 
		const mpres_t Delta, const uint64_t M_param, 
		const uint64_t l_param, const mpz_t m_1, 
		const int64_t k_2, const mpmod_t modulus_param, 
		const mpzspm_t ntt_context)
{
  const uint64_t tmplen = 3;
  const int want_x = (g_x != NULL || g_x_ntt != NULL);
  const int want_y = (g_y != NULL || g_y_ntt != NULL);
  mpres_t r_x, r_y, x0_x, x0_y, v2,
      r1_x[2], r1_y[2], r2_x[2], r2_y[2], 
      v[2], tmp[3];
  mpz_t mt, mt1, mt2;
  mpmod_t modulus; /* Thread-local copy of modulus_param */
  uint64_t i, l = l_param, offset = 0;
  uint64_t M = M_param;
  long timestart, realstart;
  int want_output = 1;

  outputf (OUTPUT_VERBOSE, "Computing %s%s%s", 
	   (want_x) ? "g_x" : "", 
	   (want_x && want_y) ? " and " : "",
	   (want_y) ? "g_y" : "");
  timestart = cputime ();
  realstart = realtime ();

  /* When multi-threading, we adjust the parameters for each thread */
#ifdef _OPENMP
#pragma omp parallel if (l > 100) private(r_x, r_y, x0_x, x0_y, v2, r1_x, r1_y, r2_x, r2_y, v, tmp, mt, mt1, mt2, modulus, i, l, offset, M, want_output)
  {
    want_output = (omp_get_thread_num() == 0);
    if (want_output)
      outputf (OUTPUT_VERBOSE, " using %d threads", omp_get_num_threads());
#endif
    get_chunk (&offset, &l, l_param);
    M = M_param - offset;
    mpmod_init_set (modulus, modulus_param);
    mpres_init (r_x, modulus);
    mpres_init (r_y, modulus);
    mpres_init (x0_x, modulus);
    mpres_init (x0_y, modulus);
    mpres_init (v2, modulus);
    for (i = 0; i < 2UL; i++)
      {
	mpres_init (r1_x[i], modulus);
	mpres_init (r1_y[i], modulus);
	mpres_init (r2_x[i], modulus);
	mpres_init (r2_y[i], modulus);
	mpres_init (v[i], modulus);
      }
    for (i = 0; i < tmplen; i++)
      mpres_init (tmp[i], modulus);
    mpz_init (mt);
    mpz_init (mt1);
    mpz_init (mt2);
    
    if (want_output && test_verbose (OUTPUT_TRACE))
      {
	mpres_get_z (mt, Delta, modulus);
	outputf (OUTPUT_TRACE, 
		 "\n/* pp1_sequence_g */ w = quadgen (4*%Zd); ", mt);
	outputf (OUTPUT_TRACE, 
                 "P = %" PRIu64 "; M = %" PRIu64 "; k_2 = %" PRId64 ,
		 P, M, k_2);
	outputf (OUTPUT_TRACE, 
                 "; m_1 = %Zd; N = %Zd;/* PARI */\n", 
		 m_1, modulus->orig_modulus);
	
	outputf (OUTPUT_TRACE, "/* pp1_sequence_g */ b_1 = ");
	gfp_ext_print (b1_x, b1_y, modulus, OUTPUT_TRACE);
	outputf (OUTPUT_TRACE, "; /* PARI */\n");
	outputf (OUTPUT_TRACE, 
		 "/* pp1_sequence_g */ r = b_1^P; /* PARI */\n");
	outputf (OUTPUT_TRACE, "/* pp1_sequence_g */ "
		 "x_0 = b_1^(2*k_2 + (2*m_1 + 1) * P); /* PARI */\n");
	outputf (OUTPUT_TRACE, 
		 "/* pp1_sequence_g */ addrec(x) = x + 1/x; /* PARI */\n");
      }
    
    /* Compute r */
    gfp_ext_pow_norm1_sl (r_x, r_y, b1_x, b1_y, P, Delta, modulus, 
			  tmplen, tmp);
    if (want_output && test_verbose (OUTPUT_TRACE))
      {
	outputf (OUTPUT_TRACE, "/* pp1_sequence_g */ r == ");
	gfp_ext_print (r_x, r_y, modulus, OUTPUT_TRACE);
	outputf (OUTPUT_TRACE, " /* PARI C */\n");
      }
    
    /* Compute x0 = x_0 */
    mpz_set_int64(mt1, k_2);
    mpz_set_uint64(mt2, P);
    mpz_mul_2exp (mt, m_1, 1UL);
    mpz_add_ui (mt, mt, 1UL);
    mpz_mul (mt, mt, mt2);
    mpz_addmul_ui (mt, mt1, 2UL); /* mt = 2*k_2 + (2*m_1 + 1) * P */
    gfp_ext_pow_norm1 (x0_x, x0_y, b1_x, b1_y, mt, Delta, modulus, 
		       tmplen, tmp);
    if (want_output && test_verbose (OUTPUT_TRACE))
      {
	outputf (OUTPUT_TRACE, "/* pp1_sequence_g */ x_0 == ");
	gfp_ext_print (x0_x, x0_y, modulus, OUTPUT_TRACE);
	outputf (OUTPUT_TRACE, " /* PARI C */\n");
      }
    
    
    /* Compute g[1] = r1[0] = x0^M * r^(M^2) = (x0 * r^M)^M.
       We use v[0,1] as temporary storage */
    gfp_ext_pow_norm1_sl (v[0], v[1], r_x, r_y, M, Delta, modulus, 
			  tmplen, tmp); /* v[0,1] = r^M */
    gfp_ext_mul (v[0], v[1], v[0], v[1], x0_x, x0_y, Delta, modulus, 
		 tmplen, tmp); /* v[0,1] = r^M * x_0 */
    gfp_ext_pow_norm1_sl (r1_x[0], r1_y[0], v[0], v[1], M, Delta, modulus, 
			  tmplen, tmp); /* r1[0] = (r^M * x_0)^M */
    if (g_x != NULL)
      mpres_get_z (g_x[offset], r1_x[0], modulus);
    if (g_y != NULL)
      mpres_get_z (g_y[offset], r1_y[0], modulus);
    if (g_x_ntt != NULL)
      {
	mpres_get_z (mt, r1_x[0], modulus);
	mpzspv_from_mpzv (g_x_ntt, offset, &mt, 1UL, ntt_context);
      }
    if (g_y_ntt != NULL)
      {
	mpres_get_z (mt, r1_y[0], modulus);
	mpzspv_from_mpzv (g_y_ntt, offset, &mt, 1UL, ntt_context);
      }
    
    
    /* Compute g[1] = r1[1] = x0^(M-1) * r^((M-1)^2) = (x0 * r^(M-1))^(M-1). 
       We use v[0,1] as temporary storage. FIXME: simplify, reusing g_0 */
    gfp_ext_pow_norm1_sl (v[0], v[1], r_x, r_y, M - 1, Delta, modulus, 
			  tmplen, tmp);
    gfp_ext_mul (v[0], v[1], v[0], v[1], x0_x, x0_y, Delta, modulus, 
		 tmplen, tmp);
    gfp_ext_pow_norm1_sl (r1_x[1], r1_y[1], v[0], v[1], M - 1, Delta, 
			  modulus, tmplen, tmp);
    if (g_x != NULL)
      mpres_get_z (g_x[offset + 1], r1_x[1], modulus);
    if (g_y != NULL)
      mpres_get_z (g_y[offset + 1], r1_y[1], modulus);
    if (g_x_ntt != NULL)
      {
	mpres_get_z (mt, r1_x[1], modulus);
	mpzspv_from_mpzv (g_x_ntt, offset + 1, &mt, 1UL, ntt_context);
      }
    if (g_y_ntt != NULL)
      {
	mpres_get_z (mt, r1_y[1], modulus);
	mpzspv_from_mpzv (g_y_ntt, offset + 1, &mt, 1UL, ntt_context);
      }
    
    
    /* x0 := $x_0 * r^{2M - 3}$ */
    /* We don't need x0 after this so we overwrite it. We use v[0,1] as 
       temp storage for $r^{2M - 3}$. */
    gfp_ext_pow_norm1_sl (v[0], v[1], r_x, r_y, 2UL*M - 3UL, Delta, modulus,
			  tmplen, tmp);
    gfp_ext_mul (x0_x, x0_y, x0_x, x0_y, v[0], v[1], Delta, modulus,
		 tmplen, tmp);
    
    /* Compute r2[0] = r1[0] * r^2 and r2[1] = r1[1] * r^2. */
    /* We only need $r^2$ from here on, so we set r = $r^2$ */
    gfp_ext_sqr_norm1 (r_x, r_y, r_x, r_y, modulus);  
    gfp_ext_mul (r2_x[0], r2_y[0], r1_x[0], r1_y[0], r_x, r_y, Delta, 
		 modulus, tmplen, tmp);
    gfp_ext_mul (r2_x[1], r2_y[1], r1_x[1], r1_y[1], r_x, r_y, Delta, 
		 modulus, tmplen, tmp);
    
    /* v[1] := $x_0 * r^{2*M - 3} + 1/(x_0 * r^{2M - 3}) */
    mpres_add (v[1], x0_x, x0_x, modulus);
    /* x0 := x0 * r = $x_0 * r^{2M - 1}$ */
    gfp_ext_mul (x0_x, x0_y, x0_x, x0_y, r_x, r_y, Delta, modulus,
		 tmplen, tmp);
    /* v[0] := $x_0 * r^{2M - 1} + 1/(x_0 * r^{2M - 1}) */
    mpres_add (v[0], x0_x, x0_x, modulus);
    
    /* v2 = V_2 (r + 1/r) = r^2 + 1/r^2 */
    mpres_add (v2, r_x, r_x, modulus);
    
    /* We don't need the contents of r any more and use it as a temp var */
    
    for (i = 2; i < l; i++)
      {
	if (want_x)
	  {
	    /* r1[i] = r2[i-1] * v[i-2] - r2[i-2], with indices of r2 and i 
	       taken modulo 2. We store the new r1_x[i] in r_x for now */
	    mpres_mul (r_x, r2_x[1 - i % 2], v[i % 2], modulus);
	    mpres_sub (r_x, r_x,            r2_x[i % 2], modulus);
	    /* r2[i] = r2[i-1] * v[i-1] - r1[i-2] */
	    mpres_mul (r2_x[i % 2], r2_x[1 - i % 2], v[1 - i % 2], modulus);
	    mpres_sub (r2_x[i % 2], r2_x[i % 2],     r1_x[i % 2], modulus);
	    mpres_set (r1_x[i % 2], r_x, modulus); /* FIXME, avoid this copy */
	    if (g_x != NULL)
	      mpres_get_z (g_x[offset + i], r_x, modulus); /* FIXME, avoid these REDC */
	    if (g_x_ntt != NULL)
	      {
		mpres_get_z (mt, r_x, modulus);
		mpzspv_from_mpzv (g_x_ntt, offset + i, &mt, 1UL, ntt_context);
	      }
	  }
	
	if (want_y)
	  {
	    /* Same for y coordinate */
	    mpres_mul (r_y, r2_y[1 - i % 2], v[i % 2], modulus);
	    mpres_sub (r_y, r_y,             r2_y[i % 2], modulus);
	    mpres_mul (r2_y[i % 2], r2_y[1 - i % 2], v[1 - i % 2], modulus);
	    mpres_sub (r2_y[i % 2], r2_y[i % 2],     r1_y[i % 2], modulus);
	    mpres_set (r1_y[i % 2], r_y, modulus);
	    if (g_y != NULL)
	      mpres_get_z (g_y[offset + i], r_y, modulus); /* Keep r1, r2 in mpz_t ? */
	    if (g_y_ntt != NULL)
	      {
		mpres_get_z (mt, r_y, modulus);
		mpzspv_from_mpzv (g_y_ntt, offset + i, &mt, 1UL, ntt_context);
	      }
	  }
	
	/* v[i] = v[i - 1] * V_2(a + 1/a) - v[i - 2] */
	mpres_mul (r_x, v[1 - i % 2], v2, modulus);
	mpres_sub (v[i % 2], r_x, v[i % 2], modulus);
	if (want_output && test_verbose (OUTPUT_TRACE))
	  {
	    mpz_t t;
	    mpz_init (t);
	    mpres_get_z (t, v[i % 2], modulus);
	    outputf (OUTPUT_TRACE, "/* pp1_sequence_g */ "
		     "addrec(x_0 * r^(2*(M-%lu) - 1)) == %Zd /* PARI C */\n", 
		     i, t);
	    mpz_clear (t);
	  }
      }
    
    mpres_clear (r_x, modulus);
    mpres_clear (r_y, modulus);
    mpres_clear (x0_x, modulus);
    mpres_clear (x0_y, modulus);
    mpres_clear (v2, modulus);
    for (i = 0; i < 2; i++)
      {
	mpres_clear (r1_x[i], modulus);
	mpres_clear (r1_y[i], modulus);
	mpres_clear (r2_x[i], modulus);
	mpres_clear (r2_y[i], modulus);
	mpres_clear (v[i], modulus);
      }
    for (i = 0; i < tmplen; i++)
      mpres_clear (tmp[i], modulus);
    mpz_clear (mt);
    mpz_clear (mt1);
    mpz_clear (mt2);
    mpmod_clear (modulus);
#ifdef _OPENMP
  }
#endif
  
  print_elapsed_time (OUTPUT_VERBOSE, timestart, realstart);

  if (g_x != NULL && g_y != NULL && test_verbose(OUTPUT_TRACE))
    {
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


/* Compute r[i] = b1^(-P*(k+i)^2) * f_i for i = 0, 1, ..., l-1, where "b1" is 
   an element of norm 1 in the quadratic extension ring */

static void
pp1_sequence_h (listz_t h_x, listz_t h_y, mpzspv_t h_x_ntt, mpzspv_t h_y_ntt,
		const listz_t f, const mpres_t b1_x, const mpres_t b1_y, 
		const int64_t k_param, const uint64_t l_param, 
		const uint64_t P, const mpres_t Delta, 
		mpmod_t modulus_param, const mpzspm_t ntt_context)
{
  uint64_t i;
  long timestart, realstart;

  if (l_param == 0UL)
    return;

  ASSERT (f != h_x);
  ASSERT (f != h_y);

  outputf (OUTPUT_VERBOSE, "Computing h_x and h_y");
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
      gfp_ext_print (b1_x, b1_y, modulus_param, OUTPUT_TRACE);
      outputf (OUTPUT_TRACE, "; r = b_1^P; rn = b_1^(-P); /* PARI */\n");
      for (i = 0; i < l_param; i++)
	outputf (OUTPUT_TRACE, 
		 "/* pp1_sequence_h */ f_%" PRIu64 " = %Zd; /* PARI */\n", 
		 i, f[i]);
      mpz_clear (t);
    }

#ifdef _OPENMP
#pragma omp parallel if (l_param > 100) private(i)
#endif
  {
    const uint64_t tmplen = 2;
    mpres_t s_x[3], s_y[3], s2_x[2], s2_y[2], v[2], V2, rn_x, rn_y, 
      tmp[2];
    mpmod_t modulus; /* Thread-local copy of modulus_param */
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
    mpmod_init_set (modulus, modulus_param);

    /* Init the local mpres_t variables */
    for (i = 0; i < 2; i++)
      {
	mpres_init (s_x[i], modulus);
	mpres_init (s_y[i], modulus);
	mpres_init (s2_x[i], modulus);
	mpres_init (s2_y[i], modulus);
	mpres_init (v[i], modulus);
      }
    mpres_init (s_x[2], modulus);
    mpres_init (s_y[2], modulus);
    mpres_init (V2, modulus);
    mpres_init (rn_x, modulus);
    mpres_init (rn_y, modulus);
    for (i = 0; i < tmplen; i++)
      mpres_init (tmp[i], modulus);

    /* Compute rn = b_1^{-P}. It has the same value for all threads,
       but we make thread local copies anyway. */
    gfp_ext_pow_norm1_sl (rn_x, rn_y, b1_x, b1_y, P, Delta, modulus, tmplen, 
			  tmp);
    mpres_neg (rn_y, rn_y, modulus);
    
    /* Compute s[0] = rn^(k^2) = r^(-k^2). We do it by two 
       exponentiations by k and use v[0] and v[1] as temp storage */
    gfp_ext_pow_norm1_sl (v[0], v[1], rn_x, rn_y, k, Delta, modulus, 
			  tmplen, tmp);
    gfp_ext_pow_norm1_sl (s_x[0], s_y[0], v[0], v[1], k, Delta, modulus, 
			  tmplen, tmp);
    if (test_verbose (OUTPUT_TRACE))
      {
#ifdef _OPENMP
#pragma omp critical
#endif
	{
	  outputf (OUTPUT_TRACE, 
	  	"/* pp1_sequence_h */ rn^(%" PRId64 "^2) == ", k);
	  gfp_ext_print (s_x[0], s_y[0], modulus, OUTPUT_TRACE);
	  outputf (OUTPUT_TRACE, " /* PARI C */\n");
	}
      }
    
    /* Compute s[1] = r^(-(k+1)^2) = r^(-(k^2 + 2k + 1))*/
    if (l > 1)
      {
	/* v[0] + v[1]*sqrt(Delta) still contains rn^k */
	gfp_ext_sqr_norm1 (s_x[1], s_y[1], v[0], v[1], modulus);
	/* Now s[1] = r^(-2k) */
	gfp_ext_mul (s_x[1], s_y[1], s_x[1], s_y[1], s_x[0], s_y[0], Delta, 
		     modulus, tmplen, tmp);
	/* Now s[1] = r^(-(k^2 + 2k)) */
	gfp_ext_mul (s_x[1], s_y[1], s_x[1], s_y[1], rn_x, rn_y, Delta, 
		     modulus, tmplen, tmp);
	/* Now s[1] = r^(-(k^2 + 2k + 1)) = r^(-(k+1)^2) */
	if (test_verbose (OUTPUT_TRACE))
	  {
#ifdef _OPENMP
#pragma omp critical
#endif
	    {
	      outputf (OUTPUT_TRACE, 
		  	"/* pp1_sequence_h */ rn^(%" PRId64 "^2) == ", k + 1);
	      gfp_ext_print (s_x[1], s_y[1], modulus, OUTPUT_TRACE);
	      outputf (OUTPUT_TRACE, " /* PARI C */\n");
	    }
	  }
      }
    
    /* Compute s2[0] = r^(k^2+2) = r^(k^2) * r^2 */
    gfp_ext_sqr_norm1 (v[0], v[1], rn_x, rn_y, modulus);
    gfp_ext_mul (s2_x[0], s2_y[0], s_x[0], s_y[0], v[0], v[1], Delta, modulus, 
		 tmplen, tmp);
    if (test_verbose (OUTPUT_TRACE))
      {
#ifdef _OPENMP
#pragma omp critical
#endif
	{
	  outputf (OUTPUT_TRACE, "/* pp1_sequence_h */ rn^(%ld^2+2) == ", k);
	  gfp_ext_print (s2_x[0], s2_y[0], modulus, OUTPUT_TRACE);
	  outputf (OUTPUT_TRACE, " /* PARI C */\n");
	}
      }
    
    /* Compute a^((k+1)^2+2) = a^((k+1)^2) * a^2 */
    gfp_ext_mul (s2_x[1], s2_y[1], s_x[1], s_y[1], v[0], v[1], Delta, modulus, 
		 tmplen, tmp);
    if (test_verbose (OUTPUT_TRACE))
      {
#ifdef _OPENMP
#pragma omp critical
#endif
	{
	  outputf (OUTPUT_TRACE, "/* pp1_sequence_h */ rn^(%ld^2+2) == ", 
		   k + 1);
	  gfp_ext_print (s2_x[1], s2_y[1], modulus, OUTPUT_TRACE);
	  outputf (OUTPUT_TRACE, " /* PARI C */\n");
	}
      }
    
    /* Compute V_2(r + 1/r). Since 1/r = rn_x - rn_y, we have r+1/r = 2*rn_x.
       V_2(x) = x^2 - 2, so we want 4*rn_x^2 - 2. */
    mpres_add (V2, rn_x, rn_x, modulus); /* V2 = r + 1/r  = 2*rn_x */
    V (v[0], V2, 2 * k + 1, modulus);  /* v[0] = V_{2k+1} (r + 1/r) */
    V (v[1], V2, 2 * k + 3, modulus);  /* v[1] = V_{2k+3} (r + 1/r) */
    mpres_mul (V2, V2, V2, modulus); /* V2 = 4*a_x^2 */
    mpres_sub_ui (V2, V2, 2UL, modulus); /* V2 = 4*a_x^2 - 2 */
    if (test_verbose (OUTPUT_TRACE))
      {
#ifdef _OPENMP
#pragma omp critical
#endif
	{
	  mpres_get_z (mt, V2, modulus);
	  outputf (OUTPUT_TRACE, "/* pp1_sequence_h */ r^2 + 1/r^2 == %Zd "
		   "/* PARI C */\n", mt);
	  mpres_get_z (mt, v[0], modulus);
	  outputf (OUTPUT_TRACE, "/* pp1_sequence_h */ r^(2*%ld+1) + "
		   "1/r^(2*%ld+1) == %Zd /* PARI C */\n", 
		   (long)k, (long)k, mt);
	  mpres_get_z (mt, v[1], modulus);
	  outputf (OUTPUT_TRACE, "/* pp1_sequence_h */ r^(2*%ld+3) + "
		   "1/r^(2*%ld+3) == %Zd /* PARI C */\n", 
		   (long)k, (long)k, mt);
	}
      }
    
    for (i = 0; i < 2UL && i < l; i++)
      {
	/* Multiply the 2nd coordinate by Delta, so that after the polynomial
	   multipoint evaluation we get x1 + Delta*x2 */
	mpres_mul (s_y[i], s_y[i], Delta, modulus);
	mpres_mul (s2_y[i], s2_y[i], Delta, modulus);
	
	if (h_x != NULL)
	  mpres_mul_z_to_z (h_x[i + offset], s_x[i], f[i + offset], modulus);
	if (h_y != NULL)
	  mpres_mul_z_to_z (h_y[i + offset], s_y[i], f[i + offset], modulus);
	if (h_x_ntt != NULL)
	  {
	    mpres_mul_z_to_z (mt, s_x[i], f[i + offset], modulus);
	    mpzspv_from_mpzv (h_x_ntt, i + offset, &mt, 1UL, ntt_context);
	  }
	if (h_y_ntt != NULL)
	  {
	    mpres_mul_z_to_z (mt, s_y[i], f[i + offset], modulus);
	    mpzspv_from_mpzv (h_y_ntt, i + offset, &mt, 1UL, ntt_context);
	  }
      }
    
    /* Compute the remaining r^((k+i)^2) values according to Peter's 
       recurrence */
    
    for (i = 2; i < l; i++)
      {
	if (h_x != NULL || h_x_ntt != NULL)
	  {
	    /* r[i] = r2[i-1] * v[i-2] - r2[i-2], with indices of r2 and i 
	       taken modulo 2 */
	    mpres_mul (s_x[i % 3], s2_x[1 - i % 2], v[i % 2], modulus);
	    mpres_sub (s_x[i % 3], s_x[i % 3], s2_x[i % 2], modulus);
	    
	    /* r2[i] = r2[i-1] * v[i-1] - r[i-2] */
	    mpres_mul (s2_x[i % 2], s2_x[1 - i % 2], v[1 - i % 2], modulus);
	    mpres_sub (s2_x[i % 2], s2_x[i % 2], s_x[(i - 2) % 3], modulus);
	    if (h_x != NULL)
	      mpres_mul_z_to_z (h_x[i + offset], s_x[i % 3], f[i + offset], 
				modulus);
	    if (h_x_ntt != NULL)
	      {
		mpres_mul_z_to_z (mt, s_x[i % 3], f[i + offset], modulus);
		mpzspv_from_mpzv (h_x_ntt, i + offset, &mt, 1UL, ntt_context);
	      }
	  }
	
	if (h_y != NULL || h_y_ntt != NULL)
	  {
	    /* Same for y coordinate */
	    mpres_mul (s_y[i % 3], s2_y[1 - i % 2], v[i % 2], modulus);
	    mpres_sub (s_y[i % 3], s_y[i % 3], s2_y[i % 2], modulus);
	    mpres_mul (s2_y[i % 2], s2_y[1 - i % 2], v[1 - i % 2], modulus);
	    mpres_sub (s2_y[i % 2], s2_y[i % 2], s_y[(i - 2) % 3], modulus);
	    if (h_y != NULL)
	      mpres_mul_z_to_z (h_y[i + offset], s_y[i % 3], f[i + offset], 
				modulus);
	    if (h_y_ntt != NULL)
	      {
		mpres_mul_z_to_z (mt, s_y[i % 3], f[i + offset], modulus);
		mpzspv_from_mpzv (h_y_ntt, i + offset, &mt, 1UL, ntt_context);
	      }
	  }
	
	/* v[i] = v[i - 1] * V_2(a + 1/a) - v[i - 2] */
	mpres_mul (tmp[0], v[1 - i % 2], V2, modulus);
	mpres_sub (v[i % 2], tmp[0], v[i % 2], modulus);
      }
    
    /* Clear the local mpres_t variables */
    for (i = 0; i < 2; i++)
      {
	mpres_clear (s_x[i], modulus);
	mpres_clear (s_y[i], modulus);
	mpres_clear (s2_x[i], modulus);
	mpres_clear (s2_y[i], modulus);
	mpres_clear (v[i], modulus);
      }
    mpres_clear (s_x[2], modulus);
    mpres_clear (s_y[2], modulus);
    mpres_clear (V2, modulus);
    mpres_clear (rn_x, modulus);
    mpres_clear (rn_y, modulus);
    for (i = 0; i < tmplen; i++)
      mpres_clear (tmp[i], modulus);

    /* Clear the thread-local copy of modulus */
    mpmod_clear (modulus);

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
  listz_t F;   /* Polynomial F has roots X^{k_1} for k_1 \in S_1, so has 
		  degree s_1. It is symmetric, so has only s_1 / 2 + 1 
		  distinct coefficients. The sequence h_j will be stored in 
		  the same memory and won't be a monic polynomial, so the 
		  leading 1 monomial of F will be stored explicitly. Hence we 
		  need s_1 / 2 + 1 entries. */

  listz_t g_x, g_y, fh_x, fh_y, h_x, h_y, tmp, R_x, R_y; 
  const unsigned long tmpreslen = 2UL;
  mpres_t b1_x, b1_y, Delta, tmpres[2];
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
  mpres_init (b1_x, modulus);
  mpres_init (b1_y, modulus);
  mpres_init (Delta, modulus);
  for (i = 0; i < tmpreslen; i++)
      mpres_init (tmpres[i], modulus);
  lenF = params->s_1 / 2 + 1 + 1; /* Another +1 because poly_from_sets_V stores
				     the leading 1 monomial for each factor */
  lenH = params->s_1 + 1;
  lenG = params->l;
  lenR = nr;
  F = init_list2 (lenF, (unsigned int) abs (modulus->bits));
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
  i = poly_from_sets_V (F, X, &S_1, tmp, tmplen, modulus, NULL, NULL, NULL);
  ASSERT_ALWAYS(2 * i == params->s_1);
  ASSERT(mpz_cmp_ui (F[i], 1UL) == 0);
  sets_free(&S_1);
  
  outputf (OUTPUT_VERBOSE, " took %lums\n", cputime () - timestart);
  if (test_verbose (OUTPUT_TRACE))
    {
      for (i = 0; i < params->s_1 / 2 + 1; i++)
	outputf (OUTPUT_TRACE, "f_%" PRIu64 " = %Zd; /* PARI */\n", i, F[i]);
      outputf (OUTPUT_TRACE, "f(x) = f_0");
      for (i = 1; i < params->s_1 / 2 + 1; i++)
	outputf (OUTPUT_TRACE, "+ f_%" PRIu64 " * (x^%" PRIu64
                        " + x^(-%" PRIu64 "))", i, i, i);
      outputf (OUTPUT_TRACE, "/* PARI */ \n");
    }

  /* Compute Delta and b1_x + b1_y * sqrt(Delta) = X) */
  mpres_mul (Delta, X, X, modulus);
  mpres_sub_ui (Delta, Delta, 4UL, modulus);
  mpres_div_2exp (b1_x, X, 1, modulus);
  mpres_set_ui (b1_y, 1UL, modulus);
  mpres_div_2exp (b1_y, b1_y, 1, modulus);
  if (test_verbose (OUTPUT_TRACE))
    {
      mpres_get_z (mt, Delta, modulus);
      outputf (OUTPUT_TRACE, 
	       "Delta = Mod(%Zd, N); w = quadgen (4*lift(Delta)); b_1 = ", mt);
      gfp_ext_print (b1_x, b1_y, modulus, OUTPUT_TRACE);
      outputf (OUTPUT_TRACE, "; /* PARI */\n");
      outputf (OUTPUT_TRACE, "X == b_1 + 1/b_1 /* PARI C */\n");
    }

  /* Compute the h sequence h_j = b1^(P*-j^2) * f_j for 0 <= j <= s_1 */
  pp1_sequence_h (fh_x, fh_y, NULL, NULL, F, b1_x, b1_y, 0L, 
		  params->s_1 / 2 + 1, params->P, Delta, modulus, NULL);
  /* We don't need F(x) any more */
  clear_list (F, lenF);

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
  
  for (l = 0; l < params->s_2; l++)
    {
      const uint64_t M = params->l - 1 - params->s_1 / 2;
      outputf (OUTPUT_VERBOSE, "Multi-point evaluation %" PRIu64
                        " of %" PRIu64 ":\n", l + 1, params->s_2);
      pp1_sequence_g (g_x, g_y, NULL, NULL, b1_x, b1_y, params->P, 
		      Delta, M, params->l, params->m_1, s2_sumset[l], 
		      modulus, NULL);
      
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
      for (i = 0; i < nr; i++)
	  mpz_add (R_x[i], R_x[i], R_y[i]);
      
      timestart = cputime ();
      mpres_set_ui (tmpres[1], 1UL, modulus); /* Accumulate product in 
						 tmpres[1] */
      for (i = 0; i < nr; i++)
      {
	  mpres_set_z_for_gcd (tmpres[0], R_x[i], modulus);
#define TEST_ZERO_RESULT
#ifdef TEST_ZERO_RESULT
	  if (mpres_is_zero (tmpres[0], modulus))
	      outputf (OUTPUT_VERBOSE, "R_[%lu] = 0\n", i);
#endif
	  mpres_mul (tmpres[1], tmpres[1], tmpres[0], modulus); 
      }
      outputf (OUTPUT_VERBOSE, "Computing product of F(g_i)^(1) took %lums\n", 
	       cputime () - timestart);
      if (test_verbose(OUTPUT_RESVERBOSE))
      {
	  mpres_get_z (mt, tmpres[1], modulus);
	  outputf (OUTPUT_RESVERBOSE, "Product of R[i] = %Zd (times some "
		   "power of 2 if REDC was used! Try -mpzmod)\n", mt);
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
  mpres_clear (b1_x, modulus);
  mpres_clear (b1_y, modulus);
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
	    const faststage2_param_t *params, const int twopass)
{
  uint64_t nr;
  uint64_t l, lenF;
  set_list_t S_1; /* This is stored as a set of sets (arithmetic 
                       progressions of prime length */
  int64_t *s2_sumset; /* set of sums of S_2 */
  uint64_t s2_sumset_size;
  listz_t F;   /* Polynomial F has roots X^{k_1} for k_1 \in S_1, so has 
		  degree s_1. It is symmetric, so has only s_1 / 2 + 1 
		  distinct coefficients. The sequence h_j will be stored in 
		  the same memory and won't be a monic polynomial, so the 
		  leading 1 monomial of F will be stored explicitly. Hence we 
		  need s_1 / 2 + 1 entries. */
  listz_t R = NULL;  /* Is used only for two-pass convolution, has nr 
			entries. R is only ever referenced if twopass == 1,
			but gcc does not realize that and complains about
			uninitialized value, so we set it to NULL. */
  mpzspm_t ntt_context;
  mpzspv_t g_x_ntt, g_y_ntt, h_x_ntt, h_y_ntt;
  mpres_t b1_x, b1_y, Delta;
  mpz_t mt;   /* All-purpose temp mpz_t */
  mpz_t product;
  mpz_t *product_ptr = NULL;
  int youpi = ECM_NO_FACTOR_FOUND;
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

  /* Allocate memory for F with correct amount of space for each mpz_t */
  lenF = params->s_1 / 2 + 1 + 1; /* Another +1 because poly_from_sets_V stores
				     the leading 1 monomial for each factor */
  MEMORY_TAG;
  F = init_list2 (lenF, (unsigned int) abs (modulus->bits) + GMP_NUMB_BITS);
  MEMORY_UNTAG;
  
  /* Build F */
  if (build_F_ntt (F, X, &S_1, params, NULL, modulus) == ECM_ERROR)
    {
      sets_free (&S_1);
      free (s2_sumset);
      mpz_clear (mt);
      mpzspm_clear (ntt_context);
      clear_list (F, lenF);
      return ECM_ERROR;
    }

  sets_free (&S_1);
  
  mpres_init (b1_x, modulus);
  mpres_init (b1_y, modulus);
  mpres_init (Delta, modulus);

  /* Compute Delta and b1_x + b1_y * sqrt(Delta) = X) */
  mpres_mul (Delta, X, X, modulus);
  mpres_sub_ui (Delta, Delta, 4UL, modulus);
  mpres_div_2exp (b1_x, X, 1, modulus);
  mpres_set_ui (b1_y, 1UL, modulus);
  mpres_div_2exp (b1_y, b1_y, 1, modulus);
  if (test_verbose (OUTPUT_TRACE))
    {
      mpres_get_z (mt, Delta, modulus);
      outputf (OUTPUT_TRACE, 
	       "Delta = Mod(%Zd, N); w = quadgen (4*lift(Delta)); b_1 = ", mt);
      gfp_ext_print (b1_x, b1_y, modulus, OUTPUT_TRACE);
      outputf (OUTPUT_TRACE, "; /* PARI */\n");
      outputf (OUTPUT_TRACE, "X == b_1 + 1/b_1 /* PARI C */\n");
    }

  /* Allocate remaining memory for h_ntt */
  h_x_ntt = mpzspv_init (params->l / 2 + 1, ntt_context);
  h_y_ntt = mpzspv_init (params->l / 2 + 1, ntt_context);
  /* Compute the h_j sequence */
  pp1_sequence_h (NULL, NULL, h_x_ntt, h_y_ntt, F, b1_x, b1_y, 0L, 
		  params->s_1 / 2 + 1, params->P, Delta, modulus, 
		  ntt_context);
  /* We don't need F(x) any more */
  clear_list (F, lenF);

  /* compute the forward transform of h and store the distinct coefficients 
     in h_ntt */
  g_x_ntt = mpzspv_init (params->l, ntt_context);
  if (twopass)
    {
      g_y_ntt = g_x_ntt;
      MEMORY_TAG;
      R = init_list2 (nr, (mpz_size (modulus->orig_modulus) + 2) *  
                          GMP_NUMB_BITS);
      MEMORY_UNTAG;
    }
  else
    g_y_ntt = mpzspv_init (params->l, ntt_context);
  
  /* Compute DCT-I of h_x and h_y */
  outputf (OUTPUT_VERBOSE, "Computing DCT-I of h_x");
#ifdef _OPENMP
  outputf (OUTPUT_VERBOSE, " using %d threads", omp_get_thread_limit());
#endif
  timestart = cputime ();
  realstart = realtime ();
  mpzspv_to_dct1 (h_x_ntt, h_x_ntt, params->s_1 / 2 + 1, params->l / 2 + 1,
		  ntt_context);
  print_elapsed_time (OUTPUT_VERBOSE, timestart, realstart);

  outputf (OUTPUT_VERBOSE, "Computing DCT-I of h_y");
#ifdef _OPENMP
  outputf (OUTPUT_VERBOSE, " using %d threads", omp_get_thread_limit());
#endif
  timestart = cputime ();
  realstart = realtime ();
  mpzspv_to_dct1 (h_y_ntt, h_y_ntt, params->s_1 / 2 + 1, params->l / 2 + 1,
		  ntt_context);
  print_elapsed_time (OUTPUT_VERBOSE, timestart, realstart);

  if (test_verbose (OUTPUT_RESVERBOSE))
    {
      mpz_init (product);
      product_ptr = &product;
    }

  for (l = 0; l < params->s_2; l++)
    {
      const uint64_t M = params->l - 1 - params->s_1 / 2;

      outputf (OUTPUT_VERBOSE, "Multi-point evaluation %" PRIu64 
                        " of %" PRIu64 ":\n", l + 1, params->s_2);
      if (twopass)
	{
	  /* Two-pass variant. Two separate convolutions, 
	     then addition in Z/NZ */
	  pp1_sequence_g (NULL, NULL, g_x_ntt, NULL, b1_x, b1_y, params->P, 
			  Delta, M, params->l, params->m_1, s2_sumset[l], 
			  modulus, ntt_context);

	  /* Do the convolution product of g_x * h_x */
	  outputf (OUTPUT_VERBOSE, "Computing g_x*h_x");
#ifdef _OPENMP
          outputf (OUTPUT_VERBOSE, " using %d threads", omp_get_thread_limit());
#endif
	  timestart = cputime ();
	  realstart = realtime ();
          mpzspv_mul_ntt (g_x_ntt, 0, g_x_ntt, 0, params->l, h_x_ntt, 0, params->l / 2,
              params->l, 0, 0, ntt_context, 
              NTT_MUL_STEP_FFT1 + NTT_MUL_STEP_MULDCT + NTT_MUL_STEP_IFFT);
	  /* Store the product coefficients we want in R */
	  mpzspv_to_mpzv (g_x_ntt, params->s_1 / 2, R, nr, ntt_context);
	  print_elapsed_time (OUTPUT_VERBOSE, timestart, realstart);

	  /* Compute g_y sequence */
	  pp1_sequence_g (NULL, NULL, NULL, g_y_ntt, b1_x, b1_y, params->P, 
			  Delta, M, params->l, params->m_1, s2_sumset[l], 
			  modulus, ntt_context);
	  
	  /* Do the convolution product of g_y * (Delta * h_y) */
	  outputf (OUTPUT_VERBOSE, "Computing g_y*h_y");
#ifdef _OPENMP
          outputf (OUTPUT_VERBOSE, " using %d threads", omp_get_thread_limit());
#endif
	  timestart = cputime ();
	  realstart = realtime ();
          mpzspv_mul_ntt (g_y_ntt, 0, g_y_ntt, 0, params->l, h_y_ntt, 0, params->l / 2,
              params->l, 0, 0, ntt_context, 
              NTT_MUL_STEP_FFT1 + NTT_MUL_STEP_MULDCT + NTT_MUL_STEP_IFFT);
	  print_elapsed_time (OUTPUT_VERBOSE, timestart, realstart);
	  
	  /* Compute product of sum of coefficients and gcd with N */
	  ntt_gcd (mt, product_ptr, g_y_ntt, NULL, params->s_1 / 2, R, nr, 
		   ntt_context, modulus);
	}
      else
	{
	  /* One-pass variant. Two forward transforms and point-wise products,
	     then addition and single inverse transform */
	  pp1_sequence_g (NULL, NULL, g_x_ntt, g_y_ntt, b1_x, b1_y, params->P, 
			  Delta, M, params->l, params->m_1, s2_sumset[l], 
			  modulus, ntt_context);

	  outputf (OUTPUT_VERBOSE, "Computing forward NTT of g_x");
#ifdef _OPENMP
          outputf (OUTPUT_VERBOSE, " using %d threads", omp_get_thread_limit());
#endif
	  timestart = cputime ();
	  realstart = realtime ();
          mpzspv_mul_ntt (g_x_ntt, 0, g_x_ntt, 0, params->l, h_x_ntt, 0, params->l / 2,
              params->l, 0, 0, ntt_context, 
              NTT_MUL_STEP_FFT1 + NTT_MUL_STEP_MULDCT);
	  print_elapsed_time (OUTPUT_VERBOSE, timestart, realstart);
	  
	  outputf (OUTPUT_VERBOSE, "Computing forward NTT of g_y");
#ifdef _OPENMP
          outputf (OUTPUT_VERBOSE, " using %d threads", omp_get_thread_limit());
#endif
	  timestart = cputime ();
	  realstart = realtime ();
          mpzspv_mul_ntt (g_y_ntt, 0, g_y_ntt, 0, params->l, h_y_ntt, 0, params->l / 2,
              params->l, 0, 0, ntt_context, 
              NTT_MUL_STEP_FFT1 + NTT_MUL_STEP_MULDCT);
	  print_elapsed_time (OUTPUT_VERBOSE, timestart, realstart);
	  
	  outputf (OUTPUT_VERBOSE, "Adding and computing inverse NTT of sum");
#ifdef _OPENMP
          outputf (OUTPUT_VERBOSE, " using %d threads", omp_get_thread_limit());
#endif
	  timestart = cputime ();
	  realstart = realtime ();
	  mpzspv_add (g_x_ntt, (spv_size_t) 0, g_x_ntt, (spv_size_t) 0, 
	              g_y_ntt, (spv_size_t) 0, params->l, ntt_context);
          mpzspv_mul_ntt (g_x_ntt, 0, g_x_ntt, 0, params->l, NULL, 0, 0,
              params->l, 0, 0, ntt_context, NTT_MUL_STEP_IFFT);
	  print_elapsed_time (OUTPUT_VERBOSE, timestart, realstart);
	  
	  ntt_gcd (mt, product_ptr, g_x_ntt, NULL, params->s_1 / 2, NULL, nr, 
		   ntt_context, modulus);
	}
      
      outputf (OUTPUT_RESVERBOSE, "Product of R[i] = %Zd (times some "
	       "power of 2 if REDC was used! Try -mpzmod)\n", product);

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
  mpzspv_clear (g_x_ntt, ntt_context);
  if (twopass)
    clear_list (R, nr);
  else
    mpzspv_clear (g_y_ntt, ntt_context);
  mpzspv_clear (h_x_ntt, ntt_context);
  mpzspv_clear (h_y_ntt, ntt_context);
  mpzspm_clear (ntt_context);
  mpz_clear (mt);
  mpres_clear (b1_x, modulus);
  mpres_clear (b1_y, modulus);
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