/* 
  Functions for sets of long ints, to factor (Z/NZ)* into a set of sums
  as described in section 5 of "Improved Stage 2 to $P\pm{}1$ Factoring 
  Algorithms" by Peter L. Montgomery and Alexander Kruppa, ANTS 2008
  (8th Algorithmic Number Theory Symposium).
  
  Copyright 2007, 2008 Alexander Kruppa.

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

#include "ecm-impl.h"
#include <stdlib.h>
#include <string.h>
#ifdef HAVE_ALLOCA_H
#include <alloca.h>
#endif
#ifdef TESTDRIVE
#include <stdio.h>
FILE *ECM_STDOUT, *ECM_STDERR;
#endif

/*****************************************************************

          Functions for processing sets

*****************************************************************/
void
sets_init (set_list_t *L)
{
  memset(L, 0, sizeof(set_list_t));
  L->num_sets_alloc = 10;
  L->sets = (set_t *)malloc (L->num_sets_alloc * sizeof(set_t));
}


void
sets_free (set_list_t *L)
{
  uint32_t i;

  for (i = 0; i < L->num_sets; i++)
    {
      free (L->sets[i].elem);
    }

  free (L->sets);
  memset(L, 0, sizeof(set_list_t));
}


/* Returns the smallest prime factor of N. If N == 1, return 1. */
static int64_t
smallest_factor (const uint64_t N)
{
  unsigned long i;

  ASSERT_ALWAYS (N != 0);

  if (N == 1)
    return 1;

  if (N % 2 == 0)
    return 2;

  for (i = 3; i*i <= N; i += 2)
    if (N % i == 0)
      return i;

  return N;
}

/* Returns max(S), where S == (Z/\beta Z)* as chosen by
   sets_get_factored_sorted() */
/* Assumes that S == 0 at recursion entry */
static void
sets_max_recurse (mpz_t S, mpz_t tmp, const uint64_t beta)
{
  uint64_t P = beta; 
  uint64_t inc = 0;
  uint64_t p, pk;
  uint32_t k;
  
  if (beta == 1UL)
    return;
  
  p = smallest_factor (P);
  k = 1; pk = p; P /= p;
  while (P % p == 0)
    {
      k++;
      pk *= p;
      P /= p; /* P*pk=beta is invariant */
    }
  sets_max_recurse (S, tmp, P);

  mpz_set_uint64 (tmp, pk);
  mpz_mul (S, S, tmp);

  if (p == 2 && k == 1)
    inc = P;
  else if (p == 2)
    inc = P * (pk / 2 - 1);
  else if (p % 4 == 1)
    inc = P * ((pk + p) / 2 - 2);
  else if (p % 4 == 3)
    inc = P * ((pk - 1) / 2);
  else
    abort();

  mpz_set_uint64 (tmp, inc);
  mpz_add (S, S, tmp);
}

void
sets_max (mpz_t S, const uint64_t beta)
{
  mpz_t tmp;

  mpz_set_ui (S, 0UL);
  mpz_init (tmp);
  sets_max_recurse (S, tmp, beta);
  mpz_clear (tmp);
}


/* Compute the set of sums over the different sets in "*sets".
   The value of "add" is added to each element of the set of sums. 
   "*sum" will have {\prod_{S \in "*sets"} #S} entries and must have
   enough memory allocated. Returns number of elements in the 
   set of sums */

static uint32_t 
sets_sumset_recurse (int64_t *sum, const set_list_t *sets, 
                    const uint32_t which_set, const int64_t add)
{
  uint32_t i, j = 0;
  set_t * curr_set;

  if (which_set == sets->num_sets)
    {
      sum[0] = add;
      return 1UL;
    }

  curr_set = sets->sets + which_set;
  for (i = 0UL; i < curr_set->card; i++)
    {
      int64_t elem = curr_set->elem[i];
      /* Test for overflow */
      ASSERT_ALWAYS (add <= 0 || add + elem > elem);
      ASSERT_ALWAYS (add >= 0 || add + elem < elem);
      j += sets_sumset_recurse (sum + j, sets, 
		      	which_set + 1, add + elem);
    }

  return j;
}


uint32_t
sets_sumset (int64_t *sum, const set_list_t *sets)
{
  return sets_sumset_recurse (sum, sets, 0, 0);
}


uint32_t
sets_sumset_size (const set_list_t *L)
{
  uint32_t i;
  uint32_t size = 1;

  for (i = 0; i < L->num_sets; i++)
    size *= L->sets[i].card;

  return size;
}


/* Returns the minimal (if minmax == -1) or maximal (minmax == 1) value
   in the set of sums over the sets in "*sets". */

void
sets_sumset_minmax (mpz_t sum, const set_list_t *sets, const int minmax)
{
  uint32_t i, j;
  int64_t extremum;
  mpz_t tmp;

  ASSERT (minmax == 1 || minmax == -1);
  mpz_set_ui (sum, 0UL);
  mpz_init (tmp);

  for (i = 0; i < sets->num_sets; i++)
    {
      set_t *curr_set = sets->sets + i;

      ASSERT (curr_set->card > 0UL);

      extremum = curr_set->elem[0];
      for (j = 1UL; j < curr_set->card; j++)
      	{
	  int64_t elem = curr_set->elem[j];

	  if ((minmax == -1 && elem < extremum) ||
	      (minmax == 1 && elem > extremum))
	    extremum = elem;
	}
      
      if (extremum >= 0)
	{
	  mpz_set_uint64 (tmp, extremum);
	  mpz_add (sum, sum, tmp);
	}
      else
	{
	  mpz_set_uint64 (tmp, -extremum);
	  mpz_sub (sum, sum, tmp);
	}
    }

  mpz_clear (tmp);
}


static void
sets_add_new (set_list_t *L, const uint32_t card)
{
  set_t *curr_set;

  if (L->num_sets == L->num_sets_alloc)
  {
    L->num_sets_alloc *= 2;
    L->sets = (set_t *)realloc (L->sets, L->num_sets_alloc * sizeof(set_t));
  }

  curr_set = L->sets + L->num_sets++;
  curr_set->card = card;
  curr_set->elem = (int64_t *)malloc (card * sizeof(int64_t));
}


/* Store in *L arithmetic progressions of prime length whose sumset is 
   k/2*R_n, an arithmetic progression centered at 0 of common difference k 
   and cardinality n. If n is even, k must be as well to ensure integer
   results.
   I.e. n = 1: k/2*R_n = {0}, 
        n = 2: k/2*R_n = k/2 * {1, -1}, 
        n = 3: k/2*R_n = k * {-1, 0, 1}, 
        n = 4: k/2*R_n = k/2 * {-3, -1, 1, 3}, 
        n = 5: k/2*R_n = k * {-2, -1, 0, 1, 2} etc. 
*/

static void
sets_factored_Rn2 (set_list_t *L, const uint32_t n, const uint64_t k)
{
  uint32_t i, q;
  uint64_t m, r;

  /* n must be odd, or n and k both even */
  ASSERT_ALWAYS(n % 2 == 1 || k % 2 == 0);

  m = k; /* The multiplier accumulated so far, init to k */
  r = n; /* The remaining cofactor of n */
  for (q = 2; r > 1; q = (q + 1) | 1) /* Find prime factors of n */
    {
      ASSERT (q <= r);
      while (r % q == 0L)
	{
	  /* Add m*R_q/2 to list */

	  set_t *curr_set;

	  sets_add_new(L, q);
	  curr_set = L->sets + L->num_sets - 1;
	  for (i = 0L; i < q; i++)
	    {
	      const int64_t t = (int64_t)m * ((int64_t)2 * i - q + 1);
	      ASSERT(t % 2L == 0L);
	      curr_set->elem[i] = t / 2;
	    }

	  /* Multiply this t to multiplier and treat remaining
	     factors of the set */
	  m *= q;
	  r /= q;
	}
    }
}


static int
sort_set_ascending(const void *x, const void *y)
{
  set_t *xx = (set_t *)x;
  set_t *yy = (set_t *)y;
  return (int)xx->card - (int)yy->card;
}


/* Return a set L of sets M_i so that M_1 + ... + M_k is congruent to 
   (Z/nZ)*, which is the set of residue classes coprime to n. The M_i all
   have prime cardinality and are sorted by increasing cardinality
*/

void
sets_get_factored_sorted (set_list_t *L, const uint64_t n)
{
  uint64_t r, np;
  uint32_t k, p;
  
  ASSERT (n > 0UL);

  r = n;
  while (r > 1)
    {
      p = smallest_factor (r);
      for (k = 0; r % p == 0; k++, r /= p); /* Find p^k || r */
      np = n/p;

      if (p == 2 && k == 1) /* Case 2^1. Deal with it before the */
        {		       /* while loop below decreases k. */
	  sets_add_new (L, 1);
	  L->sets[L->num_sets - 1].elem[0] = np;
        }

      /* If k > 1, do the \sum_{i=1}^{k-1} p^i (Z/pZ) part here.
	 (Z/pZ) is represented by an arithmetic progression of
	 common difference 1 and length p. */
		
      while (k-- > 1)
        {
	  sets_factored_Rn2 (L, p, np);
	  np /= p;
        }

      if (p % 4 == 3)
        {
	  /* We can use \hat{S}_p. Factor as 
	     {-(p+1)/4, (p+1)/4} + C_{(p-1)/2} */
	  
	  /* Add the {-(p+1)/4, (p+1)/4} set to L */
	  sets_factored_Rn2 (L, 2, (p + 1) / 2 * np);

	  /* Add the np / 2 * R_{(p-1)/2} set to L */
	  sets_factored_Rn2 (L, (p - 1) / 2, np);
        }
      else if (p % 4 == 1)
        {
	  /* Factor into arithmetic progressions of prime length.
	     R_{p} = {-p+1, -p+3, ..., p-3, p+1}, i.e.
	     R_2 = {-1, 1}, R_3 = {-2, 0, 2}, R_4 = {-3, -1, 1, 3}
	     We have R_{sq} = R_q + q*R_s */
	  
	  sets_factored_Rn2 (L, p - 1, 2 * np);
        }
    }

  qsort (L->sets, L->num_sets, sizeof(set_t), sort_set_ascending);
}


/* Print all the sets in *L, formatted as a sum of sets */

void
sets_print (const int verbosity, set_list_t *L)
{
  unsigned long i, j;

  for (i = 0; i < L->num_sets; i++)
    {
      set_t *curr_set = L->sets + i;

      if (i == 0UL)
        outputf (verbosity, "{");
      else
	outputf (verbosity, " + {");

      ASSERT(curr_set->card > 0);
      outputf (verbosity, "%" PRId64, curr_set->elem[0]);
      for (j = 1; j < curr_set->card; j++)
        outputf (verbosity, ", %" PRId64, curr_set->elem[j]);
      outputf (verbosity, "}");
  }
  outputf (verbosity, "\n");
}


/* Extract sets whose set of sums has cardinality "d". We expect that
   "d" divides the cardinality of the set of sums of L and that
   the cardinalities of the sets in L are all prime. */

void
sets_extract (set_list_t *extracted, set_list_t *L, 
              const uint64_t d)
{
  uint64_t remaining_d = d;
  uint32_t i, j;

  ASSERT_ALWAYS (d > 0UL);

  if (d == 1UL)
    {
      /* d == 1 means we need to extract a set of cardinality 1, which we
         most likely don't have in L. (FIXME: check for set of 
         cardinality 1?) We return the set containing only zero, which
         can be added to any set of sets without changing the set of sums */

      sets_add_new (extracted, 1);
      extracted->sets[extracted->num_sets - 1].elem[0] = 0;
      return;
    }

  for (i = j = 0; i < L->num_sets; i++)
    {
      uint32_t c = L->sets[i].card;

      if (remaining_d % c == 0)
	{
	  if (extracted->num_sets == extracted->num_sets_alloc)
	    {
	      extracted->num_sets_alloc *= 2;
	      extracted->sets = (set_t *)realloc (extracted->sets, 
		  			extracted->num_sets_alloc *
					sizeof(set_t));
	    }

	  extracted->sets[extracted->num_sets++] = L->sets[i];
	  remaining_d /= c;
	}
      else
	{
	  L->sets[j++] = L->sets[i];
	}
    }

  L->num_sets = j;
  ASSERT_ALWAYS (remaining_d == 1);
}


#ifdef TESTDRIVE

static int
sort_int64_ascending(const void *x, const void *y)
{
  int64_t *xx = (int64_t *)x;
  int64_t *yy = (int64_t *)y;
  return (int)(*xx - *yy);
}

static void
selftest (const uint64_t beta)
{
  set_list_t L;
  int64_t *sumset;
  uint32_t sumset_size;
  uint32_t i, j, phibeta;
  mpz_t max, tmp;

  ASSERT_ALWAYS (beta > 0);

  sets_init (&L);
  sets_get_factored_sorted (&L, beta);
  
#if 0
  printf("sets = \n");
  sets_print (OUTPUT_TRACE, &L);
#endif

  /* Test that the sumset % beta is equal to (Z/betaZ)* % beta */
  phibeta = eulerphi (beta);
  sumset_size = sets_sumset_size(&L);
  ASSERT_ALWAYS(sumset_size > 0);

  sumset = malloc (sumset_size * sizeof(int64_t));
  if (sumset == NULL)
    {
      fprintf (stderr, "Cannot allocate memory in selftest\n");
      exit (1);
    }
  sets_sumset (sumset, &L);

  /* Also test that max (sumset) == sets_max (beta) */
  mpz_init (max);
  mpz_init (tmp);
  sets_max (max, beta);
  if (phibeta > 0)
    {
      int64_t maxelem;

      ASSERT_ALWAYS (sumset_size == phibeta);
      maxelem = sumset[0];
      for (i = 1; i < phibeta; i++)
	if (maxelem < sumset[i])
	  maxelem = sumset[i];

      mpz_set_int64(tmp, maxelem);
      ASSERT_ALWAYS (mpz_cmp (max, tmp) == 0);
    }
  else
    {
      ASSERT_ALWAYS (mpz_cmp_ui (max, 0UL) == 0);
    }
  mpz_clear (max);
  mpz_clear (tmp);

#if 0
  printf ("sumset, before reduction: ");
  for (i = 0; i < phibeta; i++)
    printf ("%" PRId64 "%s", sumset[i], i < phibeta-1 ? ", " : "\n");
#endif

  for (i = 0; i < phibeta; i++)
    {
      sumset[i] = (sumset[i] < 0) ?  beta - (-sumset[i]) % beta : 
				sumset[i] % beta;
      ASSERT_ALWAYS (sumset[i] >= 0L);
      ASSERT_ALWAYS (sumset[i] < beta);
    }
#if 0
  printf ("sumset, after reduction: ");
  for (i = 0; i < phibeta; i++)
    printf ("%" PRId64 "%s", sumset[i], i < phibeta-1 ? ", " : "\n");
#endif

  qsort(sumset, sumset_size, sizeof(int64_t), sort_int64_ascending);

#if 0
  printf ("sumset, after sorting: ");
  for (i = 0; i < phibeta; i++)
    printf ("%" PRId64 "%s", sumset[i], i < phibeta-1 ? ", " : "\n");
#endif

  j = 0;
  for (i = 1; i < beta; i++)
    {
      if (gcd (i, beta) == 1)
        {
          if (sumset[j] != i)
            {
              printf ("sumset->elem[%u] = %" PRId64 " != %u\n", 
                      j, sumset[j], i);
              abort();
            }
          j++;
        }
    }
  free (sumset);
  sets_free (&L);
}

int
main (int argc, char **argv)
{
  uint64_t beta;
  const uint64_t selftest_max = 1000;
  int loop = 1;

  ECM_STDOUT = stdout;
  ECM_STDERR = stderr;
  
  if (argc > 1)
    {
      beta = atol (argv[1]);
      loop = 0;
    }
  
  if (!loop)
    set_verbose (OUTPUT_TRACE);

  if (!loop)
    selftest (beta);
  else
    {
      printf ("Testing beta = 1, ..., %" PRId64 "\n", selftest_max);
      for (beta = 1; beta < selftest_max; beta++)
        selftest (beta);
    }

  return 0;
}
#endif
