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

#include "config.h"
#include "ecm-impl.h"
#include <stdlib.h>
#ifdef HAVE_ALLOCA_H
#include <alloca.h>
#endif
#ifdef TESTDRIVE
#include <stdio.h>
FILE *ECM_STDOUT, *ECM_STDERR;
#endif

/*****************************************************************

          Functions for processing sets

  A set is a cardinality of unsigned long type and an array of 
  long ints. 
  A set of sets is an unsigned long telling the number of sets,
  an array that has several sets stored back-to-back.

*****************************************************************/


/* Copy a set from "*S" to "*T". Assumes that the sets do not overlap,
   or that T < S. */

static void
set_copy (set_long_t *T, set_long_t *S)
{
  unsigned long i;
  const unsigned long c = S->card; /* We might overwrite S->card */
  
  T->card = c;
  for (i = 0UL; i < c; i++)
    T->elem[i] = S->elem[i];
}


/* Exchange two adjacent sets in memory. */

static void 
set_swap (set_long_t *T)
{
  set_long_t *next, *tmp;
  
  next = sets_nextset (T);
  tmp = alloca (set_sizeof (T->card));
  ASSERT(tmp != NULL);
  set_copy (tmp, T);
  set_copy (T, next);
  set_copy (sets_nextset(T), tmp);  
}


/* Functions for sorting an array of longs */

static inline void
swap_long (long *a, long *b)
{
  long t;
  t = *a;
  *a = *b;
  *b = t;
}

static inline void 
swapsort_long (long *a, long *b)
{
  if (*a > *b)
    swap_long (a, b);
}

void 
quicksort_long (long *a, unsigned long l)
{
  unsigned long i, j;
  long pivot;

  if (l < 2)
    return;

  j = l - 1;
  swapsort_long (a, a+j);
  if (l == 2)
    return;

  i = j / 2;
  swapsort_long (a, a+i);
  swapsort_long (a+i, a+j);
  if (l == 3)
    return;

  pivot = a[i]; /* Median of three */

  /* Stuff <= pivot goes in first list */

  /* Invariant: a[0 ... i-1] <= pivot, a[j+1 ... l-1] > pivot */
  for (i = 1; i < j;)
    if (a[i] > pivot)
      {
	for (; a[j] > pivot; j--);
	if (i < j)
	  swap_long (a+(i++), a+j);
      }
    else
      i++;

#ifdef WANT_ASSERT
  for (j = 0; j < i; j++)
    ASSERT (a[j] <= pivot);
  for (j = i; j < l; j++)
    ASSERT(a[j] > pivot);
#endif

  quicksort_long (a, i);
  quicksort_long (a + i, l - i);

#ifdef WANT_ASSERT
  for (j = 0; i < l - 1; i++)
    ASSERT (a[j] <= a[j + 1]);
#endif
}


/* Returns max(S), where S == (Z/\beta Z)* as chosen by
   sets_get_factored_sorted() */
/* Assumes that S == 0 at recursion entry */
static void
sets_max_recurse (mpz_t S, const unsigned long beta)
{
  unsigned long P = beta, p, pk;
  unsigned int k;
  
  if (beta == 1UL)
    return;
  
  p = find_factor (P);
  k = 1; pk = p; P /= p;
  while (P % p == 0)
    {
      k++;
      pk *= p;
      P /= p; /* P*pk=beta is invariant */
    }
  sets_max_recurse (S, P);
  mpz_mul_ui (S, S, pk);
  if (p == 2UL && k == 1)
    mpz_add_ui (S, S, P);
  else if (p == 2UL)
    mpz_add_ui (S, S, P * (pk / 2UL - 1UL));
  else if (p % 4UL == 1UL)
    mpz_add_ui (S, S, P * ((pk + p) / 2UL - 2UL));
  else if (p % 4UL == 3UL)
    mpz_add_ui (S, S, P * ((pk - 1UL) / 2UL));
  else
    abort();
}

void
sets_max (mpz_t S, const unsigned long beta)
{
  mpz_set_ui (S, 0UL);
  sets_max_recurse (S, beta);
}


/* Compute the set of sums over the "nr_sets" different sets in "*sets".
   The value of "add" is added to each element of the set of sums. 
   "*sum" will have {\prod_{S \in "*sets"} #S} entries and must have
   enough memory allocated. This number of elements in the set of sums 
   is the return value. In case of nr_sets == 0, "add" is written to *sets 
   and 1 is returned. The sets in "*sets" are assumed to be non-empty.
   If "*sum" is NULL, nothing is written, but the return value is computed
   correctly. */

static unsigned long 
sets_sumset_recurse (long *sum, const set_long_t *sets, 
                    const unsigned long nr_sets, const long add)
{
  unsigned long i, j = 0UL;

  if (nr_sets == 0UL)
    {
      if (sum != NULL)
        sum[0] = add;
      return 1UL;
    }

  ASSERT (sets->card > 0UL);
  for (i = 0UL; i < sets->card; i++)
    {
      /* Test for overflow */
      ASSERT_ALWAYS (add <= 0 || add + sets->elem[i] > sets->elem[i]);
      ASSERT_ALWAYS (add >= 0 || add + sets->elem[i] < sets->elem[i]);
      j += sets_sumset_recurse (sum + j, sets_nextset(sets), nr_sets - 1UL, 
				add + sets->elem[i]);
    }

  return j;
}


void 
sets_sumset (set_long_t *sum, const sets_long_t *sets)
{
  sum->card = sets_sumset_recurse (sum->elem, sets->sets, sets->nr, 0L);
}


/* Returns the minimal (if minmax == -1) or maximal (minmax == 1) value
   in the set of sums over the sets in "*sets". */

void
sets_sumset_minmax (mpz_t sum, const sets_long_t *sets, const int minmax)
{
  unsigned long i, nr;
  const set_long_t *set = sets->sets;
  long extremum;

  ASSERT (minmax == 1 || minmax == -1);
  mpz_set_ui (sum, 0UL);

  for (nr = 0; nr < sets->nr; nr++)
    {
      ASSERT (set->card > 0UL);
      extremum = set->elem[0];

      for (i = 1UL; i < set->card; i++)
	if ((minmax == -1 && set->elem[i] < extremum) ||
	    (minmax == 1 && set->elem[i] > extremum))
	  extremum = set->elem[i];
      
      if (extremum >= 0)
	mpz_add_ui (sum, sum, extremum);
      else
	mpz_sub_ui (sum, sum, -extremum);
      set = sets_nextset (set);
    }

  return;
}


/* Store in (**L) arithmetic progressions of prime length whose sumset is 
   k/2*R_n, an arithmetic progression centered at 0 of common difference k 
   and cardinality n. If n is even, k must be as well to ensure integer
   results.
   I.e. n = 1: k/2*R_n = {0}, 
        n = 2: k/2*R_n = k/2 * {1, -1}, 
        n = 3: k/2*R_n = k * {-1, 0, 1}, 
        n = 4: k/2*R_n = k/2 * {-3, -1, 1, 3}, 
        n = 5: k/2*R_n = k * {-2, -1, 0, 1, 2} etc. 
  _ADDS_ the size in bytes of the set to "*sets_size"
*/

static unsigned long
sets_factored_Rn2 (set_long_t **L, size_t *sets_size, const long n, 
                   const long k)
{
  unsigned long nr = 0UL;
  long i, m, q, r;
  size_t size = 0;

  /* n must be odd, or n and k both even */
  ASSERT_ALWAYS(n % 2L == 1L || k % 2L == 0L);
  ASSERT(L != NULL);
  m = k; /* The multiplier accumulated so far, init to k */
  r = n; /* The remaining cofactor of n */
  for (q = 2L; r > 1L; q = (q + 1L) | 1L) /* Find prime factors of n */
    {
      ASSERT (q <= r);
      while (r % q == 0L)
	{
	  if (*L != NULL)
	    {
	      /* Add m*R_q/2 to list */
	      (*L)->card = q;
	      for (i = 0L; i < q; i++)
	        {
	          const long t = m * (2L * i - q + 1L);
	          ASSERT(t % 2L == 0L);
		  (*L)->elem[i] = t / 2L;
		}
	      *L = sets_nextset (*L);
	      nr++;
	    }
	  size += set_sizeof ((unsigned long) q);
	  /* Multiply this t to multiplier and treat remaining
	     factors of the set */
	  m *= q;
	  r /= q;
	}
    }
  if (sets_size != NULL)
    *sets_size += size;
  return nr;
}


/* Return a set L of sets M_i so that M_1 + ... + M_k is congruent to 
   (Z/nZ)*, which is the set of residue classes coprime to n. The M_i all
   have prime cardinality.
   The size of the set of sets "*L" in bytes is computed and stored in  
   "*sets_size" unless "*sets_size" is NULL.
   Return the number of sets in L. 
   If L is the NULL pointer, nothing will be stored in L. The correct
   return value (number of set in L) and "*sets_size" value will still 
   be computed, for example so that the correct amount of space can be 
   allocated and factor_coprimeset() be called again.
*/

static unsigned long 
sets_factor_coprime (sets_long_t *sets, size_t *sets_size, 
                     const unsigned long n)
{
  unsigned long r, k, nr = 0UL;
  long p, np;
  size_t size = sizeof (unsigned long);
  set_long_t *set = NULL;
  
  ASSERT (n > 0UL);
  if (sets != NULL)
    set = sets->sets;

  r = n;
  while (r > 1UL)
    {
      for (p = 2L; r % p > 0L; p++); /* Find smallest prime p that divides r */
      for (k = 0UL; r % p == 0UL; k++, r /= p); /* Find p^k || r */
      np = n/p;

      if (p == 2L && k == 1UL) /* Case 2^1. Deal with it before the */
        {		       /* while loop below decreases k. */
	  if (set != NULL)
	    {
	      set->card = 1UL;
	      set->elem[0] = np;
	      set = sets_nextset (set);
	    }
	  size += set_sizeof (1UL);
	  nr++;
        }

      /* If k > 1, do the \sum_{i=1}^{k-1} p^i (Z/pZ) part here.
	 (Z/pZ) is represented by an arithmetic progression of
	 common difference 1 and length p. */
		
      while (k-- > 1UL)
        {
	  nr += sets_factored_Rn2 (&set, &size, p, np);
	  np /= p;
        }

      if (p % 4L == 3L)
        {
	  /* We can use \hat{S}_p. Factor as 
	     {-(p+1)/4, (p+1)/4} + C_{(p-1)/2} */
	  
	  /* Add the {-(p+1)/4, (p+1)/4} set to L */
	  nr += sets_factored_Rn2 (&set, &size, 2L, (p + 1L) / 2L * np);

	  /* Add the np / 2 * R_{(p-1)/2} set to L */
	  nr += sets_factored_Rn2 (&set, &size, (p - 1L) / 2L, np);
        }
      else if (p % 4L == 1L)
        {
	  /* Factor into arithmetic progressions of prime length.
	     R_{p} = {-p+1, -p+3, ..., p-3, p+1}, i.e.
	     R_2 = {-1, 1}, R_3 = {-2, 0, 2}, R_4 = {-3, -1, 1, 3}
	     We have R_{sq} = R_q + q*R_s */
	  
	  nr += sets_factored_Rn2 (&set, &size, p - 1L, 2L * np);
        }
    }
  
  if (sets_size != NULL)
    *sets_size = size;

  if (sets != NULL)
    sets->nr = nr;

  return nr;
}


/* Sort the sets in F into order of ascending cardinality. Uses a simple
   Bubble sort. */

static void 
sets_sort (sets_long_t *sets)
{
  unsigned long i, nr_unsorted, highest_swap;
  set_long_t *set;

  /* The last sets->nr - nr_unsorted sets in "*sets" are known to be
     sorted and each one larger than any of the first nr_unsorted sets 
     in "*sets". */
  nr_unsorted = sets->nr;
  while (nr_unsorted > 1UL)
    {
      outputf (OUTPUT_TRACE, "nr_unsorted = %lu. ", nr_unsorted);
      sets_print (OUTPUT_TRACE, sets);
      set = sets->sets;
      highest_swap = 1UL;
      for (i = 1UL; i < nr_unsorted; i++)
	{
	  if (set->card > sets_nextset(set)->card)
	    {
	      outputf (OUTPUT_TRACE, "sets_sort: swapping %lu and %lu\n", 
                       i - 1, i);
	      set_swap (set);
	      highest_swap = i;
	    }
          set = sets_nextset (set);
	}
      nr_unsorted = highest_swap;
    }

#ifdef WANT_ASSERT
  set = sets->sets;
  for (i = 0UL; i + 1UL < sets->nr; i++)
    {
      ASSERT(set->card <= sets_nextset (set)->card);
      set = sets_nextset (set);
    }
#endif
}

/* Print all the sets in "*sets", formatted as a sum of sets */

void
sets_print (const int verbosity, sets_long_t *sets)
{
  unsigned long i, j;
  set_long_t *set = sets->sets;

  for (i = 0UL; i < sets->nr; i++)
  {
      if (i == 0UL)
        outputf (verbosity, "{");
      else
	outputf (verbosity, " + {");

      ASSERT(set->card > 0UL);
      outputf (verbosity, "%ld", set->elem[0]);
      for (j = 1UL; j < set->card; j++)
        outputf (verbosity, ", %ld", set->elem[j]);
      outputf (verbosity, "}");
      set = sets_nextset (set);
  }
  outputf (verbosity, "\n");
}


/* Extract sets whose set of sums has cardinality "d". We expect that
   "d" divides the cardinality of the set of sums of "sets" and that
   the cardinalities of the sets in "sets" are all prime. 
   The amount of memory in bytes needed to store the extracted sets in 
   "*extracted" is stored at "*extr_size". The number of sets extracted
   is returned. (If d = p_1 * ... * p_k, the return value is k and 
   "*extr_size" is set_sizeof(p_1) + ... + set_sizeof(p_k).)
   If "*extracted" is NULL, nothing is written and no sets are removed
   from "*sets", but "*extr_size" is computed as if they were.  */

void
sets_extract (sets_long_t *extracted, size_t *extr_size, sets_long_t *sets, 
              const unsigned long d)
{
  unsigned long i, c, remaining_d = d;
  set_long_t *readfrom, *readnext, *moveto, *extractto = NULL;
  size_t extracted_size = sizeof (unsigned long);

  ASSERT_ALWAYS (d > 0UL);

  if (d == 1UL)
    {
      /* d == 1 means we need to extract a set of cardinality 1, which we
         most likely don't have in "*sets". (FIXME: check for set of 
         cardinality 1?) We return the set containing only zero, which
         can be added to any set of sets without changing the set of sums */
      if (extracted != NULL)
        {
          extracted->nr = 1;
          extractto = extracted->sets;
          extractto->card = 1UL;
          extractto->elem[0] = 0L;
        }
      if (extr_size != NULL)
        *extr_size = sizeof (unsigned long) + set_sizeof (1UL);
      return;
    }

  if (extracted != NULL)
    {
      extracted->nr = 0UL;
      extractto = extracted->sets;
    }
  /* All sets from *sets are read via *readfrom, and (assuming we actually
     extract them) are either copied to *extractto to *moveto */
  readfrom = moveto = sets->sets;
  for (i = 0UL; i < sets->nr; i++)
    {
      c = readfrom->card; /* readfrom->card may get garbled */
      readnext = sets_nextset (readfrom);
      if (remaining_d % c == 0UL)
        {
          if (extracted != NULL)
            {
              /* Copy this set to extractto */
              set_copy (extractto, readfrom);
              extractto = sets_nextset (extractto);
              extracted->nr++;
            }
          remaining_d /= c;
          extracted_size += set_sizeof (c);
        } else {
          if (extracted != NULL)
            {
              /* Move this set within "*sets", filling the gaps left by 
                 extracted sets */
              set_copy (moveto, readfrom);
              moveto = sets_nextset (moveto);
            }
        }
      readfrom = readnext;
    }

  ASSERT_ALWAYS (remaining_d == 1UL);
  if (extr_size != NULL)
    *extr_size = extracted_size;
  if (extracted != NULL)
    sets->nr -= extracted->nr;
}

sets_long_t *
sets_get_factored_sorted (const unsigned long beta)
{
  sets_long_t *sets;
  size_t size;

  sets_factor_coprime (NULL, &size, beta);
  sets = malloc (size);
  if (sets == NULL)
    return NULL;
  sets_factor_coprime (sets, NULL, beta);
  
  if (test_verbose (OUTPUT_TRACE))
    {
      outputf (OUTPUT_TRACE, 
              "sets_get_factored_sorted: Factored sets before sorting are ");
      sets_print (OUTPUT_TRACE, sets);
    }
  
  sets_sort (sets);
  
  if (test_verbose (OUTPUT_TRACE))
    {
      outputf (OUTPUT_TRACE, "Factored sets after sorting are ");
      sets_print (OUTPUT_TRACE, sets);
    }
  
  return sets;
}

#ifdef TESTDRIVE
static void
selftest (const unsigned long beta)
{
  sets_long_t *sets;
  set_long_t *sumset;
  unsigned long i, j, phibeta;
  mpz_t max;

  ASSERT_ALWAYS (beta > 0);
  sets = sets_get_factored_sorted (beta);
  
  /* Test that the sumset % beta is equal to (Z/betaZ)* % beta */
  phibeta = eulerphi (beta);
  sumset = malloc (set_sizeof (phibeta));
  sets_sumset (sumset, sets);
  ASSERT_ALWAYS (sumset->card = phibeta);
  /* Also test that max (sumset) == sets_max (beta) */
  mpz_init (max);
  sets_max (max, beta);
  if (phibeta > 0)
    {
      long maxelem;
      maxelem = sumset->elem[0];
      for (i = 1; i < phibeta; i++)
	if (maxelem < sumset->elem[i])
	  maxelem = sumset->elem[i];
      ASSERT_ALWAYS (mpz_cmp_si (max, maxelem) == 0);
    }
  else
    {
      ASSERT_ALWAYS (mpz_cmp_ui (max, 0UL) == 0);
    }
  mpz_clear (max);

 /*  printf ("sumset, before reduction: ");
  for (i = 0; i < phibeta; i++)
    printf ("%ld%s", sumset->elem[i], i < phibeta-1 ? ", " : "\n"); */
  for (i = 0; i < phibeta; i++)
    {
      sumset->elem[i] = (sumset->elem[i] < 0L) ? 
          beta - (long) ((unsigned long) (-sumset->elem[i]) % beta)
        : (unsigned long) sumset->elem[i] % beta;
      ASSERT_ALWAYS (sumset->elem[i] >= 0L);
      ASSERT_ALWAYS (sumset->elem[i] < (long) beta);
    }
  /* printf ("sumset, after reduction: ");
  for (i = 0; i < phibeta; i++)
    printf ("%ld%s", sumset->elem[i], i < phibeta-1 ? ", " : "\n"); */

  quicksort_long (sumset->elem, sumset->card);
  /* printf ("sumset, after sorting: ");
  for (i = 0; i < phibeta; i++)
    printf ("%ld%s", sumset->elem[i], i < phibeta-1 ? ", " : "\n"); */

  j = 0;
  for (i = 1; i < beta; i++)
    {
      if (gcd (i, beta) == 1)
        {
          if (sumset->elem[j] != (long) i)
            {
              printf ("sumset->elem[%ld] = %ld != %ld\n", 
                      j, sumset->elem[j], i);
              abort();
            }
          j++;
        }
    }
  free (sumset);
  free (sets);
}

int
main (int argc, char **argv)
{
  unsigned long beta;
  const unsigned long selftest_max = 1000;
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
      printf ("Testing beta = 1, ..., %lu\n", selftest_max);
      for (beta = 1; beta < selftest_max; beta++)
        selftest (beta);
    }

  return 0;
}
#endif
