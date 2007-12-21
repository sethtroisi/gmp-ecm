#include "config.h"
#include "ecm-impl.h"
#include <stdlib.h>
#ifdef TESTDRIVE
#include <stdio.h>
FILE *ECM_STDOUT, *ECM_STDERR;
#endif

/*****************************************************************

          Functions for processing sets

  A set is a cardinality and an array of long ints. 
  A set of sets is an array that has several sets stored 
  back-to-back. The number of sets in the set is not stored but
  must be supplied separately.

*****************************************************************/

/* Return the size in bytes of a set of cardinality c */

static unsigned long 
set_sizeof (const unsigned long c)
{
  return sizeof (long) + c * sizeof (unsigned long);
}

/* Return pointer to the next set in "*sets" */

static set_long_t *
set_nextset (const set_long_t *sets)
{
  return (set_long_t *) ((void *)sets + sizeof(unsigned long) + 
                         sets->card * sizeof(long));
}


/* Copy a set from "*S" to "*T" */

static void
set_copy (set_long_t *T, set_long_t *S)
{
  unsigned long i;
  
  T->card = S->card;
  for (i = 0UL; i < S->card; i++)
    T->elem[i] = S->elem[i];
}


/* Exchange two adjacent sets in memory. */

static void 
set_swap (set_long_t *T)
{
  set_long_t *next, *tmp;
  
  next = set_nextset (T);
  tmp = alloca (set_sizeof (T->card));
  ASSERT(tmp != NULL);
  set_copy (tmp, T);
  set_copy (T, next);
  set_copy (set_nextset(T), tmp);  
}


/* Functions for sorting an array of longs */

static inline void
swap_long (long *a, long *b, long *t)
{
  *t = *a;
  *a = *b;
  *b = *t;
}

static inline void 
swapsort_long (long *a, long *b, long *t)
{
  if (*a > *b)
    swap_long (a, b, t);
}

static void 
quicksort_long (long *a, unsigned long l)
{
  unsigned long i, j;
  long pivot, t;

  if (l < 2)
    return;

  j = l - 1;
  swapsort_long (a, a+j, &t);
  if (l == 2)
    return;

  i = j / 2;
  swapsort_long (a, a+i, &t);
  swapsort_long (a+i, a+j, &t);
  if (a[i] > a[j])
    swap_long (a+i, a+j, &t);
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
	  swap_long (a+(i++), a+j, &t);
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


/* Compute the set of sums over the "nr_sets" different sets in "*sets".
   The value of "add" is added to each element of the set of sums, it is
   used mostly for the recursion (i.e. the top level call should use 
   add = 0). "*sum" will have {\prod_{S \in "*sets"} #S} entries and must 
   have enough memory allocated. This number of elements in the set of sums 
   is the return value. In case of nr_sets == 0, "add" written to *sets 
   and 1 is returned. The sets in "*sets" are assumed to be non-empty. */

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

  ASSERT (sets->card > 0);
  for (i = 0UL; i < sets->card; i++)
    j += sets_sumset_recurse (sum + j, set_nextset(sets), nr_sets - 1UL, 
	              add + sets->elem[i]);

  return j;
}


unsigned long 
sets_sumset (set_long_t *sum, const set_long_t *sets, 
                    const unsigned long nr_sets)
{
  unsigned long c;
  c = sets_sumset_recurse (sum->elem, sets, nr_sets, 0L);
  sum->card = c;
  return c;                      
}


/* Returns the minimal (if minmax == -1) or maximal (minmax == 1) value
   in the set of sums over the sets in "*sets". */

long 
sets_sumset_minmax (const set_long_t *sets, unsigned long nr_sets, 
                    const int minmax)
{
  unsigned long i;
  long sum = 0, extremum;

  ASSERT (minmax == 1 || minmax == -1);

  while (nr_sets > 0UL)
    {
      ASSERT (sets->card > 0UL);
      extremum = sets->elem[0];

      for (i = 1; i < sets->card; i++)
	if ((minmax == -1 && sets->elem[i] < extremum) ||
	    (minmax == 1 && sets->elem[i] > extremum))
	  extremum = sets->elem[i];
      
      sum += extremum;
      sets = set_nextset(sets);
      nr_sets--;
    }

  return sum;
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
*/

static unsigned long
sets_factored_Rn2 (set_long_t **L, unsigned long *nr_sets, 
                   const unsigned long n, const unsigned long k)
{
  unsigned long q, m, size = 0UL, r, nr = 0;
  long i;

  ASSERT_ALWAYS(n % 2 == 1 || k % 2 == 0); /* n must be odd, or both even */
  ASSERT(L != NULL);
  m = k; /* The multiplier accumulated so far, init to k */
  r = n;
  for (q = 2UL; r > 1UL; q = (q + 1UL) | 1UL) /* Find prime factors of n */
    {
      ASSERT (q <= r);
      while (r % q == 0UL)
	{
	  if (*L != NULL)
	    {
	      /* Add m*R_q/2 to list */
	      (*L)->card = q;
	      for (i = 0L; (unsigned long) i < q; i++)
	        {
	          long t = (signed long) m * (2L * i - (signed long) q + 1L);
	          ASSERT(t % 2L == 0L);
		  (*L)->elem[i] = t / 2L;
		}
	      *L = set_nextset (*L);
	      nr++;
	    }
	  size += set_sizeof (q);
	  /* Multiply this t to multiplier and treat remaining
	     factors of the set */
	  m *= q;
	  r /= q;
	}
    }
  if (nr_sets != NULL)
    *nr_sets = nr;
  return size;
}


/* Return a set L of sets M_i so that M_1 + ... + M_k is congruent to 
   (Z/nZ)*, which is the set of residue classes coprime to n. The M_i all
   have prime cardinality.
   L is stored as #M_1, M_1, #M_2, M_2, ..., #M_k, M_k. I.e. for n=15,
   S_n = 5*S_3 + 3*S_5 = 5*{-1,1} + 3*{-3,-1,1,3} = 
         5*{-1,1} + 3*{-2, 2} + 3*{-1,1}
   L = [2, -5, 5, 2, -6, 6, 2, -3, 3]
   Return the space (in longs) needed in L. 
   If L is the NULL pointer, nothing will be stored in L. The correct
   return value (amount of space needed in L) will still be computed, for
   example so that the correct amount of space can be allocated and 
   factor_coprimeset() be called again.
*/

static int 
sets_factor_coprime (set_long_t *L, unsigned long *nr_sets,
                     const unsigned long n)
{
  unsigned long r, k, size = 0, nr = 0, Rn2nr;
  long p, q, np;
  
  ASSERT (n > 0UL);

  r = n;
  while (r > 1UL)
    {
      for (p = 2; r % p > 0; p++); /* Find smallest prime p that divides r */
      for (k = 0; r % p == 0; k++, r /= p); /* Find p^k || r */
      np = n/p;
      /* Choose \hat{S}_p or \tilde{S}_p */
      if (p == 2)
	{
	  if (k == 1) /* Case 2^1 */
	    {
	      if (L != NULL)
		{
		  L->card = 1UL;
		  L->elem[0] = np;
		  L = set_nextset (L);
		}
	      size += set_sizeof (1UL);
	      nr++;
	    }
	  else /* Case 2^k, k > 1 */
	    while (k-- > 1)
	      {
		np /= 2L;
		if (L != NULL)
		  {
		    L->card = 2UL;
		    L->elem[0] = -np;
		    L->elem[1] = np;
		    L = set_nextset (L);
		  }
		size += set_sizeof (2UL);
		nr++;
	      }
	}
      else if (p % 4L == 3L)
	{
	  /* Can use \hat{S}_p. Factor as {(p+1)/4, (p+1)/4} + C_{(p-1)/2} */
	  
	  /* If k > 1, so the \sum_{i=1}^{k-1} p^i (Z/pZ) part here.
	     (Z/pZ) is represented by an arithmetic progression of
	     common difference 1 and length p which is prime, so it
	     can't be factored any further. */
	  while (k-- > 1)
	    {
	      if (L != NULL)
		{
		  L->card = p;
		  for (q = 0; q <= p - 1; q++)
		    L->elem[q] = (q - (p-1UL)/2UL) * np;
                  L = set_nextset (L);
		}
	      size += set_sizeof (p);
	      nr++;
	      np /= p;
	    }

	  /* Add the {(p+1)/4, (p+1)/4} set to L */
	  if (L != NULL)
	    {
	      L->card = 2UL;
	      L->elem[0] = -((p+1UL)/4L) * np;
	      L->elem[1] = (p+1UL)/4L * np;
	      L = set_nextset (L);
	    }
	  size += set_sizeof (2UL);
	  nr++;

	  /* Add np / 2 * R_{(p-1)/2}/2 = np / 2 * {-(p-3)/4, ..., (p-3)/4}.*/
	  size += sets_factored_Rn2 (&L, &Rn2nr, (p - 1UL) / 2UL, np);
	  nr += Rn2nr;
	}
      else /* case p%4 == 1 */
	{
	  /* Factor into arithmetic progressions of prime length.
	     R_{p} = {-p+1, -p+3, ..., p-3, p+1}, i.e.
	     R_2 = {-1, 1}, R_3 = {-2, 0, 2}, R_4 = {-3, -1, 1, 3}
	     We have R_{sq} = R_q + q*R_s */

	  /* If k > 1, so the \sum_{i=1}^{k-1} p^i (Z/pZ) part here.
	     (Z/pZ) is represented by an arithmetic progression of
	     common difference 1 and length p which is prime, so it
	     can't be factored any further. */
	  while (k-- > 1)
	    {
	      if (L != NULL)
		{
		  L->card = p;
		  for (q = 0; q <= p - 1; q++)
		    L->elem[q] = (q - (p-1UL)/2UL) * np;
		  L = set_nextset (L);
		}
	      size += set_sizeof (p);
	      nr++;
	      np /= p;
	    }

	  /* Add np * R_{p-1} */
	  size += sets_factored_Rn2 (&L, &Rn2nr, p - 1UL, 2UL * np);
	  nr += Rn2nr;
	}
    }
  
  if (nr_sets != NULL)
    *nr_sets = nr;

  return size;
}


/* Sort the sets in F into order of ascending cardinality. Uses a simple
   Bubble sort. */

static void 
sets_sort (set_long_t *sets, unsigned long nr_sets)
{
  unsigned long i, nr_unsorted, highest_swap;
  set_long_t *this;

  /* The last nr_sets - nr_unsorted sets in "*sets" are known to be
     sorted and each one larger than any of the first nr_unsorted sets 
     in "*sets". */
  nr_unsorted = nr_sets;
  while (nr_unsorted > 1UL)
    {
      this = sets;
      for (i = 1UL; i < nr_unsorted; i++)
	{
	  highest_swap = 1UL;
	  if (this->card > set_nextset(this)->card)
	    {
	      set_swap (this);
	      highest_swap = i;
	    }
          this = set_nextset (this);
	}
      nr_unsorted = highest_swap;
    }
}

/* Print all the sets in "*sets", formatted as a sum of sets */

void
sets_print (const int verbosity, set_long_t *sets, 
            const unsigned long nr_sets)
{
  unsigned long i, j;
  for (i = 0UL; i < nr_sets; i++)
  {
      if (i == 0UL)
        outputf (verbosity, "{");
      else
	outputf (verbosity, " + {");

      ASSERT(sets->card > 0UL);
      outputf (verbosity, "%ld", sets->elem[0]);
      for (j = 1UL; j < sets->card; j++)
        outputf (verbosity, ", %ld", sets->elem[j]);
      outputf (verbosity, "}");
      sets = set_nextset (sets);
  }
  outputf (verbosity, "\n");
}


/* Extract sets whose set of sums has cardinality d. We expect that
   d divides the cardinality of the set of sums of "sets" and that
   the cardinalities of the sets in "sets" are all prime. 
   If "*extracted" is NULL, nothing is written and no sets are removed
   from "*sets".
   Returns the amount of memory needed in "*extracted" in bytes. */

unsigned long
sets_extract (set_long_t *extracted, set_long_t *sets, 
              const unsigned long nr_sets, const unsigned long d)
{
  unsigned long i, s, extracted_size = 0UL, remaining_d = d;
  set_long_t *readfrom, *moveto;

  /* All sets from *sets are read via *readfrom, and (assuming we actually
     extract them) are either copied to *extracted to *moveto */
  readfrom = moveto = sets;
  for (i = 0UL; i < nr_sets; i++)
    {
      s = readfrom->card;
      if (remaining_d % s == 0)
        {
          if (extracted != NULL)
            {
              /* Copy this set to extracted */
              set_copy (extracted, readfrom);
              extracted = set_nextset (extracted);
            }
          remaining_d /= s;
          extracted_size += set_sizeof (s);
        } else {
          if (extracted != NULL)
            {
              set_copy (moveto, readfrom);
              moveto = set_nextset (moveto);
            }
        }
      readfrom = set_nextset (readfrom);
    }

  ASSERT_ALWAYS (remaining_d == 1);
  return extracted_size;
}

set_long_t *
sets_get_factored_sorted (unsigned long *nr_sets, const unsigned long beta)
{
  set_long_t *sets;
  unsigned long size, i, nr;

  size = sets_factor_coprime (NULL, &nr, beta);
  sets = malloc (size);
  if (sets == NULL)
    return NULL;
  i = sets_factor_coprime (sets, &nr, beta);
  ASSERT_ALWAYS (i == size);
  
  if (test_verbose (OUTPUT_TRACE))
    {
      outputf (OUTPUT_TRACE, "Factored sets before sorting are ");
      sets_print (OUTPUT_TRACE, sets, nr);
    }
  
  sets_sort (sets, nr);
  
  if (test_verbose (OUTPUT_TRACE))
    {
      outputf (OUTPUT_TRACE, "Factored sets after sorting are ");
      sets_print (OUTPUT_TRACE, sets, nr);
    }
  
  if (nr_sets != NULL)
    *nr_sets = nr;
  
  return sets;
}

#ifdef TESTDRIVE
static void
selftest (const unsigned long beta)
{
  set_long_t *sets, *sumset;
  unsigned long i, j, nr_sets, phibeta;

  sets = sets_get_factored_sorted (&nr_sets, beta);

  /* Test that the sumset % beta is equal to (Z/betaZ)* % beta */
  phibeta = eulerphi (beta);
  sumset = malloc (set_sizeof (phibeta));
  sets_sumset (sumset, sets, nr_sets);
  ASSERT_ALWAYS (sumset->card = phibeta);

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
      printf ("Testing beta = 1, ..., 10000\n");
      for (beta = 1; beta < 10000; beta++)
        selftest (beta);
    }

  return 0;
}
#endif
