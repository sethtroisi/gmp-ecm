#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include "ecm-impl.h"
#include "sp.h"
#include <math.h>

#ifdef TESTDRIVE
#include <string.h>
static int verbose = 0;
static int pari = 0;
#else
const int pari = 0;
#endif

const unsigned long Pvalues[] = {
    3UL, 5UL, 9UL, 15UL, 21UL, 17UL, 27UL, 33UL, 45UL, 51UL, 63UL, 75UL, 
    105UL, 99UL, 135UL, 165UL, 195UL, 189UL, 231UL, 255UL, 315UL, 345UL, 
    357UL, 375UL, 405UL, 435UL, 525UL, 585UL, 615UL, 735UL, 765UL, 825UL, 
    945UL, 1155UL, 1065UL, 1365UL, 1305UL, 1335UL, 1575UL, 1785UL, 1995UL, 
    2145UL, 2205UL, 2415UL, 2625UL, 2805UL, 3045UL, 3465UL, 3675UL, 4095UL, 
    4305UL, 4515UL, 4725UL, 4785UL, 5355UL, 5775UL, 5985UL, 5865UL, 6825UL, 
    7245UL, 8085UL, 8925UL, 9555UL, 10395UL, 10725UL, 11025UL, 12285UL, 
    12705UL, 15015UL, 14175UL, 15225UL, 16065UL, 17325UL, 19635UL, 21945UL, 
    23205UL, 24255UL, 25935UL, 26775UL, 28875UL, 31395UL, 33495UL, 35805UL, 
    36465UL, 38115UL, 39585UL, 40425UL, 45045UL, 45885UL, 49665UL, 51765UL, 
    58905UL, 65835UL, 69615UL, 75075UL, 77805UL, 82005UL, 84315UL, 86625UL, 
    88935UL, 94185UL, 98175UL, 105105UL, 109725UL, 116025UL, 118755UL, 
    121275UL, 135135UL, 137445UL, 137655UL, 144375UL, 153615UL, 165165UL, 
    167475UL, 176715UL, 179025UL, 185955UL, 197505UL, 208845UL, 215985UL, 
    225225UL, 255255UL, 250635UL, 285285UL, 277095UL, 294525UL, 315315UL, 
    345345UL, 373065UL, 368445UL, 405405UL, 435435UL, 451605UL, 465465UL, 
    454545UL, 504735UL, 525525UL, 555555UL, 569415UL, 596505UL, 645645UL, 
    647955UL, 672945UL, 687225UL, 765765UL, 770385UL, 805035UL, 855855UL, 
    858585UL, 915915UL, 945945UL, 962115UL, 1036035UL, 1066065UL, 1119195UL, 
    1156155UL, 1276275UL, 1306305UL, 1354815UL, 1426425UL, 1456455UL, 
    1514205UL, 1576575UL, 1666665UL, 1726725UL, 1786785UL, 1789515UL, 
    1865325UL, 1996995UL, 1983135UL, 2177175UL, 2297295UL, 2327325UL, 
    2417415UL, 2567565UL, 2611455UL, 2807805UL, 2847075UL, 2878785UL, 
    3048045UL, 3161235UL, 3258255UL, 3357585UL, 3401475UL, 3533145UL, 
    3828825UL, 3918915UL, 3985905UL, 4279275UL, 4849845UL, 4789785UL, 
    4967655UL, 5180175UL, 5360355UL, 5870865UL, 5990985UL, 6561555UL, 
    6531525UL, 6891885UL, 7402395UL, 7912905UL, 8273265UL, 8580495UL, 
    8843835UL, 9444435UL, 10015005UL, 10465455UL, 10705695UL, 10885875UL, 
    11696685UL, 12267255UL, 12507495UL, 12785955UL, 13498485UL, 14549535UL, 
    14849835UL, 15570555UL, 16111095UL, 16291275UL, 17612595UL, 18123105UL, 
    18633615UL, 19684665UL, 20255235UL, 20825805UL, 22207185UL, 22717695UL, 
    24249225UL, 24819795UL, 25741485UL, 26531505UL, 28333305UL, 29354325UL, 
    30045015UL, 31396365UL, 32807775UL, 33948915UL, 33528495UL, 34879845UL, 
    37011975UL, 37522485UL, 39564525UL, 41096055UL, 43648605UL, 44219175UL, 
    45930885UL, 47222175UL, 48333285UL, 50075025UL, 51816765UL, 52777725UL, 
    55390335UL, 55547415UL, 59053995UL, 60063465UL, 61906845UL, 64579515UL, 
    66621555UL, 67492425UL, 70105035UL, 73258185UL, 74939865UL, 77224455UL, 
    79594515UL, 81876795UL, 84999915UL, 88062975UL, 91005915UL, 94189095UL, 
    98423325UL, 101846745UL, 111546435UL, 111035925UL, 115120005UL, 
    121246125UL, 124098975UL, 130945815UL, 140645505UL, 150345195UL, 
    150225075UL, 155450295UL, 158333175UL, 170255085UL, 179444265UL, 
    190285095UL, 198843645UL, 203408205UL, 206831625UL, 217222005UL, 
    229474245UL, 240705465UL, 252447195UL, 254999745UL, 269023755UL, 
    282146865UL, 287672385UL, 294076965UL, 306110805UL, 318302985UL, 
    334639305UL, 344338995UL, 354038685UL, 363738375UL, 373438065UL,
    387221835UL, 400254855UL, 421936515UL, 431636205UL, 451035585UL,
    451035585UL, 470434965UL, 480134655UL};


/* Some useful PARI functions:
   sumset(a,b) = {local(i, j, l); l = listcreate (length(a) * length(b)); for (i = 1, length(a), for (j = 1, length(b), listput(l, a[i] + b[j]))); listsort (l, 1); l}

   V(i,X) = { if (i==0, return(2)); if (i==1, return(X)); if(i%2 == 0, return (V (i/2, X)^2-2)); return (V ((i+1)/2, X) * V ((i-1)/2, X) - X)}

   U(i,X) = { if (i==0, return(0)); if (i==1, return(1)); if(i%2 == 0, return (U (i/2, X) * V(i/2,X))); return (V ((i+1)/2, X)  *U( (i-1)/2, X) + 1)}

*/

static int pariline = 0;


/*****************************************************************

          Functions for processing sets

  A set is an array of long ints. The first element of the array 
  is the length of the set including the length value itself, i.e.
  equal to cardinality + 1. It is followed by the length - 1 
  elements of the set.
  A set of sets is an array that has several sets stored 
  back-to-back.

*****************************************************************/


/* Compute the set of sums over the sets in "*sets", whose total memory is
   "setsize" words. The value of "add" is added to each element of the set 
   of sums, it is used mostly for the recursion. "*sum" will have 
   {\prod_{S \in "*sets"} #S} entries and must have enough memory allocated. 
   This number of elements in the set of sums is the return value. In case 
   of setsize == 0, nothing is written to *sets and 0 is returned. */

static unsigned long 
sum_sets (long *sum, const long *sets, const unsigned long setsize, 
	  const long add)
{
  unsigned long i, j = 0, l;

  if (setsize == 0)
    return 0;

  ASSERT (sets[0] > 1); /* Zero or negative cardinality means bug */
  l = (unsigned long) sets[0];
  ASSERT (l <= setsize);

  for (i = 1; i < l; i++)
    {
      if (setsize - l > 0)
	j += sum_sets (sum + j, sets + l, setsize - l, add + sets[i]);
      else
	sum[j++] = add + sets[i];
    }

  return j;
}


/* Returns the minimal (if minmax == -1) or maximal (minmax == 1) value
   in the set of sums over the sets in "*sets". */

static long 
sum_sets_minmax (long *sets, unsigned long setsize, int minmax)
{
  unsigned long i, l;
  long sum = 0, extremum;

  ASSERT (minmax == 1 || minmax == -1);

  while (setsize > 0)
    {
      ASSERT (sets[0] > 1); /* Zero or negative cardinality means bug */
      l = (unsigned long) sets[0];
      ASSERT (l <= setsize);

      extremum = sets[1];

      for (i = 2; i < l; i++)
	if ((minmax == -1 && sets[i] < extremum) ||
	    (minmax == 1 && sets[i] > extremum))
	  extremum = sets[i];
      
      sum += extremum;
      sets += l;
      setsize -= l;
    }

  return sum;
}


/* Return a set L of sets M_i so that M_1 + ... + M_k is congruent to 
   (Z/nZ)*, which is the set of residue classes coprime to n. The M_i all
   have prime cardinality.
   L is stored as #M_1+1, M_1, #M_2+1, M_2, ..., #M_k+1, M_k. I.e. for n=15,
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
factor_coprimeset (long *L, const unsigned long n)
{
  unsigned long r, i = 0, j, m;
  long p, q, s, np;
  
  ASSERT (n > 0);

  r = n;
  while (r > 1)
    {
      for (p = 2; r % p > 0; p++);
      r /= p;
      np = n/p;
      /* Choose \hat{S}_p or \tilde{S}_p */
      if (p == 2)
	{
	  if (r % 2 == 1) /* Case 2^1 */
	    {
	      if (L != NULL)
		{
		  L[i++] = 2L;
		  L[i++] = np;
		}
	      else
		i += 2;
	    }
	  else /* Case 2^k, k > 1 */
	    while (r % 2 == 0)
	      {
		r /= 2;
		np /= 2;
		if (L != NULL)
		  {
		    L[i++] = 3L;
		    L[i++] = -np;
		    L[i++] = np;
		  }
		else
		  i += 3;
	      }
	}
      else if (p % 4 == 3)
	{
	  /* Can use \hat{S}_p. Factor as {(p+1)/4, (p+1)/4} + C_{(p-1)/2} */
	  
	  while (r % p == 0)
	    {
	      if (L != NULL)
		{
		  L[i++] = p + 1;
		  for (q = -(p-1)/2; q <= (p-1)/2; q++)
		    L[i++] = q * np;
		}
	      else
		i += p + 1;
	      
	      r /= p;
	      np /= p;
	    }

	  /* Add the {(p+1)/4, (p+1)/4} set to L */
	  if (L != NULL)
	    {
	      L[i++] = 3L;
	      L[i++] = -((p+1)/4L) * np;
	      L[i++] = (p+1)/4L * np;
	    }
	  else
	    i += 3;

	  /* Now deal with C_{(p-1)/2} = {-(p-3)/4, ..., (p-3)/4}. 
	     We have C_{sq} = C_q + q*C_s */
	  s = (p-1) / 2; /* s is odd */
	  m = 1; /* The multiplier accumulated so far */
	  for (q = 3; s > 1; q += 2)
	    {
	      ASSERT (q <= s);
	      while (s % q == 0)
		{
		  if (L != NULL)
		    {
		      /* Add m*C_q to list */
		      L[i++] = q + 1;
		      for (j = 0; j < (unsigned long) q; j++)
			L[i++] = m * (-(q - 1) / 2 + j) * np;
		    }
		  else
		    i += q + 1;

		  /* Multiply this t to multiplier and treat remaining
		     factors of the set */
		  m *= q;
		  s /= q;
		}
	    }
	}
      else
	{
	  /* Factor into arithmetic progressions of prime length.
	     D_{p} = {-p+1, -p+3, ..., p-3, p+1}, i.e.
	     D_2 = {-1, 1}, D_3 = {-2, 0, 2}, D_4 = {-3, -1, 1, 3}
	     We have D_{sq} = D_q + q*D_s */

	  while (r % p == 0)
	    {
	      if (L != NULL)
		{
		  L[i++] = p + 1;
		  for (q = -(p-1)/2; q <= (p-1)/2; q++)
		    L[i++] = q * np;
		}
	      else
		i += p + 1;
	      
	      r /= p;
	      np /= p;
	    }

	  s = p - 1;
	  m = 1;
	  for (q = 2; s > 1; q += 1)
	    {
	      ASSERT (q <= s);
	      while (s % q == 0)
		{
		  if (L != NULL)
		    {
		      /* Add m*C_q to list */
		      L[i++] = q + 1;
		      for (j = 0; j < (unsigned long) q; j++)
			L[i++] = m * (-q + 1 + 2 * j) * np;
		    }
		  else
		    i += q + 1;

		  m *= q;
		  s /= q;
		}
	    }
	}
    }

  return i;
}


/* Exchange two adjacent sets in memory. The sets have l and m elements of
   unsigned long type, respectively. In our case, that includes the size 
   value at the start of each set. */

static void 
swap_sets (long *T, const unsigned long l, const unsigned long m)
{
  const unsigned long d = lgcd (l, m);
  unsigned long i, j;
  long t1, t2;
  
  /* Swapping two sets of length l and m, respectively, is a cyclic 
     permutation by m of a length l+m sequence. We can split that 
     into d = gcd(l, m) cyclic permutations of length (l+m)/d sequences. 
     That way we can do with 2 instead of l+m temp vars, so we don't 
     need to allocate memory */
  
  for (i = 0; i < d; i++)
    {
      t1 = T[i];
      for (j = 1; j < (l + m) / d; j += 1)
	{
	  t2 = T[(m * j + i) % (l + m)];
	  T[(m * j + i) % (l + m)] = t1;
	  t1 = t2;
	}
      /* Here, j == (l + m) / d and d|m, so
	 (m*j + i) % (l + m) == i and were're back at the start */
      ASSERT((m*j + i) % (l + m) == i);
      T[i] = t1;
    }
}


/* Sort the sets in F into order of ascending cardinality. Uses a simple
   Bubble sort. */

static void 
sort_sets (long *F, unsigned long size)
{
  unsigned long a, l, m;
  int more = 1;

  while (more)
    {
      more = 0;
      a = 0;
      while (a + F[a] < size) /* Are there at least 2 more sets? */
	{
	  l = F[a];
	  m = F[a + l];
	  ASSERT (l > 0 && m > 0);
	  ASSERT (a + l + m <= size); /* Second set must be within array */
	  if (l > m) /* Is the first set bigger than the second? */
	    {
	      swap_sets (F + a, l, m);
	      a += m; /* The bigger set moved forward m positions */
	      more = 1;
	    }
	  else
	    a += l;
	}
    }
}

static void
print_sets (int verbosity, long *sets, unsigned long setsize)
{
  unsigned long i, j;
  for (i = 0; i < setsize; )
  {
      if (i > 0UL)
	outputf (verbosity, " + ");
      outputf (verbosity, "{");
      for (j = 1; j < (unsigned long) sets[i]; j++)
      {
	  if (j > 1)
	      outputf (verbosity, ", ");
	  outputf (verbosity, "%ld", sets[i + j]);
      }
      outputf (verbosity, "}");
      i += sets[i];
  }
  outputf (verbosity, "\n");
}


/* Extract sets whose set of sums has cardinality d. We expect that
   d divides the cardinality of the set of sums of "sets" and that
   the cardinalities of the sets in "sets" are all prime. */

static unsigned long
extract_sets (long *extracted, long *sets, const unsigned long setsize, 
	      const unsigned long d)
{
    unsigned long i, j, s, extracted_size = 0UL, remaining_d = d, 
	remaining_size = setsize;

    j = 0;
    for (i = 0; i < remaining_size; )
    {
	s = sets[i];
	if (remaining_d % (s - 1) == 0)
	{
	    remaining_d /= (s - 1);
	    if (extracted == NULL) /* Only count the size */
	    {
		i += s;
	    }
	    else
	    {
		/* Copy this set to extracted */
		for (j = 0; j < s; j++)
		    extracted[extracted_size + j] = sets[i + j];
		/* Move remaining sets in "sets" forward */
		remaining_size -= s;
		for (j = i; j < remaining_size; j++)
		    sets[j] = sets[j + s];
	    }
	    extracted_size += s;
	} else {
	    i += s;
	}
    }

    ASSERT_ALWAYS (i == remaining_size);
    ASSERT_ALWAYS (remaining_d == 1);
    return extracted_size;
}

static long *
get_factored_sorted_sets (unsigned long *setsize, const unsigned long beta)
{
  long *sets;
  unsigned long size, i;

  size = factor_coprimeset (NULL, beta);
  sets = malloc (size * sizeof (long));
  ASSERT_ALWAYS (sets != NULL); /* FIXME, do error handling */
  i = factor_coprimeset (sets, beta);
  ASSERT_ALWAYS (i == size);
  
  if (test_verbose (OUTPUT_TRACE))
    {
      outputf (OUTPUT_TRACE, "Factored sets before sorting are ");
      print_sets (OUTPUT_TRACE, sets, size);
    }
  
  sort_sets (sets, size);
  
  if (setsize != NULL)
    *setsize = size;
  
  return sets;
}


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


/* Simple, slow methods for testing / finding primes */
int 
is_prime_ul (unsigned long n)
{
  unsigned long i;

  if (n < 2UL)
    return 0;

  if (n % 2UL == 0UL)
    return n == 2UL;

  for (i = 3UL; i*i <= n; i += 2UL)
    if (n % i == 0UL)
      return 0;
  
  return 1;
}

unsigned long
next_prime_ul (const unsigned long n)
{
  unsigned long m;

  if (n < 2UL)
    return 2UL;
  
  if (n == 2UL)
    return 3UL;

  m = n + 2;
  while (! is_prime_ul (m))
    m += 2UL;

  return m;
}


/* Find a divisor of n close to l. FIXME: should find only even divisors
   to ensure that s_1 is even */
unsigned long
find_nearby_divisor (const unsigned long n, const unsigned long l)
{
  unsigned long i;

  ASSERT_ALWAYS (l > 0);
  
  if (n <= l)
    return n;
  
  for (i = 0UL; i <= l - 1UL; i++)
    {
      if (n % (l - i) == 0UL)
	return l - i;
      if (n % (l + i) == 0UL)
	return l + i;
    }

  /* Unreachable, because we had l - (l - 1) = 1 as candidate divisor */
  abort ();
  return 0;
}


/* Find an even divisor of n close to l. */

unsigned long
find_nearby_even_divisor (const unsigned long n, const unsigned long l)
{
  unsigned long i;

  ASSERT_ALWAYS (l > 0);
  ASSERT (n % 2 == 0);
  
  if (n <= l)
    return n;
  
  for (i = l % 2UL; i <= l - 1UL; i += 2UL)
    {
      if (n % (l - i) == 0UL)
	return l - i;
      if (n % (l + i) == 0UL)
	return l + i;
    }

  /* Unreachable, because if l is odd, we had l - (l - 1) = 1 as candidate 
     divisor, and if l is even, we had 2 as candidate divisor */
  abort ();
  return 0;
}


/* Find the smallest prime factor of N. If N == 1, return 1. */

static unsigned long
find_factor_ul (unsigned long N)
{
  unsigned long i;

  ASSERT_ALWAYS (N != 0UL);

  if (N == 1UL)
    return 1UL;

  if (N % 2UL == 0)
    return 2UL;

  for (i = 3; i*i <= N; i += 2)
    if (N % i == 0)
      return i;

  return N;
}


static unsigned long
maxS (unsigned long P)
{
  unsigned long p, pk;
  unsigned int k;

  if (P == 1UL)
    return 0L;

  p = find_factor_ul (P);
  k = 1; pk = p; P /= p;
  while (P % p == 0)
    {
      k++;
      pk *= p;
      P /= p;
    }

  if (p % 4UL == 1UL)
    return (P * ((pk + p) / 2UL - 2UL) + pk * maxS(P));
  if (p % 4UL == 3UL)
    return (P * ((pk - 1UL) / 2UL) + pk * maxS(P));

  abort();
}

int
test_P (const mpz_t B2min, const mpz_t B2, mpz_t m_1, const unsigned long P, 
	const unsigned long nr, mpz_t effB2min, mpz_t effB2)
{
  unsigned long m = maxS(P);
  /* We need B2min >= 2 * max(S_1 + S_2) + (2*m_1 - 1)*P + 1, or
     B2min - 2 * max(S_1 + S_2) - 1 >= (2*m_1)*P - P, or
     (B2min - 2*max(S_1 + S_2) + P - 1)/(2P) >= m_1
     Choose m_1 accordingly */
  
  mpz_sub_ui (m_1, B2min, 2UL * m + 1UL);
  mpz_add_ui (m_1, m_1, P);
  mpz_fdiv_q_ui (m_1, m_1, 2UL * P);
  
  /* Compute effB2min = 2 * max(S_1 + S_2) + (2*(m_1 - 1) + 1)*P + 1 */
  
  mpz_mul_2exp (effB2min, m_1, 1UL);
  mpz_sub_ui (effB2min, effB2min, 1UL);
  mpz_mul_ui (effB2min, effB2min, P);
  mpz_add_ui (effB2min, effB2min, 2UL * m + 1UL);
  ASSERT_ALWAYS (mpz_cmp (effB2min, B2min) <= 0);

  /* Compute the smallest value coprime to P at the high end of the stage 2
     interval that will not be covered: 
     2*(min(S_1 + S_2)) + (2*(m_1 + nr) + 1)*P. 
     We assume min(S_1 + S_2) = -max(S_1 + S_2) */
  mpz_add_ui (effB2, m_1, nr);
  mpz_mul_2exp (effB2, effB2, 1UL);
  mpz_add_ui (effB2, effB2, 1UL);
  mpz_mul_ui (effB2, effB2, P);
  mpz_sub_ui (effB2, effB2, 2UL*m);

  /* The effective B2 values is that value, minus 1 */
  mpz_sub_ui (effB2, effB2, 1UL);

  return (mpz_cmp (B2, effB2) <= 0);
}


/* Choose s_1 so that s_1 * s_2 = phiP, s_1 is positive and even, 
   s_2 > min_s2 and s_2 is minimal under those conditions. 
   Returns 0 if no such choice is possible */

unsigned long 
choose_s_1 (const unsigned long phiP, const unsigned long min_s2,
	    const unsigned long lmax)
{
  unsigned long s_1, s_2;

  s_2 = MAX(1UL, min_s2);
  if (s_2 > phiP / 2) /* s_1 >= 2, so s_1 *  min_s2 > phiP */
    return 0;

  for (; s_2 < phiP / 2 && phiP % s_2 != 0UL && (phiP / s_2) % 2 == 0 ; s_2++);
  if (phiP % s_2 != 0 || (phiP / s_2) % 2 != 0)
    return 0;
  s_1 = find_nearby_even_divisor (phiP / s_2, lmax / 2);
  ASSERT (s_1 % 2 == 0);
  return s_1;
}

/* Choose P so that a stage 2 range of length B2len can be covered with
   multipoint evaluations, each using a convolution of length lmax. 
   The parameters for stage 2 are stored in finalparams, the final effective
   B2min and B2 values in final_B2min and final_B2, respecively. Each of these
   may be NULL, in which case the value is not stored. It is permissible
   to let B2min and final_B2min, or B2 and final_B2 point at the same mpz_t. */

long
choose_P (const mpz_t B2min, const mpz_t B2, const unsigned long lmax,
	  const unsigned long min_s2, faststage2_param_t *finalparams, 
	  mpz_t final_B2min, mpz_t final_B2)
{
  /* Let S_1 + S_2 == (Z/PZ)* (mod P).

     Let F(x) = \prod_{k_1 \in S_1} (x - b_1^{2 k_1}).

     If we evaluate F(b_1^{2 k_2 + (2m + 1)P}) for all k_2 \in S_2 with 
     m_1 <= m < m_1+nr, the largest value coprime to P at the 
     low end of the stage 2 interval *not* covered will be 
       2*max(S_2) + (2*m_1 - 1)*P - 2*min(S_1).
     The smallest value at the high end not covered will be
       2*min(S_2) + (2*m_1 + 2*nr + 1)*P - 2*max(S_1).
     Assume S_1 and S_2 are symmetric around 0, so i.e. max(S_1) = -min(S_1).
     Then the largest ... is:
       2*(max(S_1) + max(S_2)) + (2*m_1 - 1)*P
     The smallest ... is:
       -2*(max(S_1) + max(S_2)) + (2*m_1 + 2*nr + 1)*P

     Then the difference is
       -4*(max(S_1) + max(S_2)) + 2P*(nr + 1)

     max(S_1) + max(S_2) = max(S_1 + S_2) <= P-1, let's use P tho.
     Then highest not covered is at least P*(2*m_1 +3), so we need
     P*(2*m_1 +3) < B2min.
     Smallest not covered is at most P*(2*m_1 + 2*nr - 3), so we need
     B2 < P*(2*m_1 + 2*nr - 3).
     Subtracting both yields
     2*P*(nr - 3) > B2 - B2min
     
     Hence we are looking for an odd P with s_1 * s_2 = eulerphi(P) so that
     s_1 ~= lmax / 2 and the whole stage 2 interval is covered. s_2 should 
     be small, as long as s_1 is small enough.

  */

  mpz_t m_1, t, B2l, effB2min, effB2;
  unsigned long P = 0, s_1 = 0, s_2 = 0, l;
  unsigned int i;
  const unsigned int Pvalues_len = sizeof (Pvalues) / sizeof (unsigned long);

  if (mpz_cmp (B2, B2min) < 0)
    return 0L;

  mpz_init (m_1);
  mpz_init (t);
  mpz_init (B2l);
  mpz_init (effB2min);
  mpz_init (effB2);

  /* Find the smallest P that can cover the B2 - B2min interval */
  /* We have nr < lmax, so we want 2*P*(lmax - 3) > B2l,
     or P >= B2l / (2*lmax - 6) */
  mpz_sub (B2l, B2, B2min);
  mpz_tdiv_q_ui (t, B2l, 2*lmax - 6UL);
  outputf (OUTPUT_DEVVERBOSE, "choose_P: We need P >= %Zd\n", t);
  for (i = 0; i < Pvalues_len; i++)
    if (mpz_cmp_ui (t, Pvalues[i]) <= 0)
      break;

  for ( ; i <  Pvalues_len; i++)
    {
      unsigned long phiP;
      /* Now a careful check to see if this P is large enough */
      P = Pvalues[i];
      phiP = eulerphi (P);
      s_1 = choose_s_1 (phiP, min_s2, lmax);
      if (s_1 == 0)
	continue;
      s_2 = phiP / s_1;
      outputf (OUTPUT_DEVVERBOSE, 
	       "Testing P = %lu, phiP = %lu, s_1 = %lu, s_2 = %lu, nr = %lu\n", 
	       P, phiP, s_1, s_2, lmax - s_1);
      if (test_P (B2min, B2, m_1, P, lmax - s_1, effB2min, effB2))
	{
	    outputf (OUTPUT_DEVVERBOSE, 
		     "This P is acceptable, B2 = %Zd\n", effB2);
	    break;
	}
      else
	outputf (OUTPUT_DEVVERBOSE, "Not good enough, trying next P\n");
  }
  
  if (i == Pvalues_len)
    return ECM_ERROR; /* Could not find suitable P */

  /* s_2 cannot be reduced further, but the transform length l could. */
  l = lmax;
  while (l / 2 > s_1 && test_P (B2min, B2, m_1, P, l / 2 - s_1, effB2min, t))
  {
      l /= 2;
      outputf (OUTPUT_DEVVERBOSE, 
	       "Reducing transform length to %ld, effB2 = %Zd\n", l, t);
      mpz_set (effB2, t);
  }

  for ( ; i + 1 < Pvalues_len; i++)
    {
      unsigned long tryP, tryphiP, trys_1, trys_2;

      /* We only found the smallest P that works so far. Maybe a larger one
	 works as well, and better */
      tryP = Pvalues[i + 1];
      tryphiP = eulerphi (tryP);
      /* tryphiP is strictly increasing and trys_2 >= tryphiP / l. Stop if
	 we can't possibly find a trys_2 <= s_2 any more */
      if (s_2 < tryphiP / l)
	break;
      trys_1 = choose_s_1 (tryphiP, min_s2, l);
      if (trys_1 == 0)
	continue;
      trys_2 = tryphiP / trys_1;
      outputf (OUTPUT_DEVVERBOSE, "Trying if P = %lu, phiP = %lu, s_1 = %lu, "
	       "s_2 = %lu works as well\n", tryP, tryphiP, trys_1, trys_2);
      if (trys_2 > s_2) /* We want to keep the minimal */
      {                 /* number of multipoint evaluations */
	  outputf (OUTPUT_DEVVERBOSE, "No, s_2 would become %lu\n", trys_2);
	  /* break; */
	  continue;
      }
      if (!test_P (B2min, B2, m_1, tryP, l - trys_1, effB2min, t))
      {
	  outputf (OUTPUT_DEVVERBOSE, 
		   "No, does not cover B2min - B2 range, effB2 = %Zd\n", t);
      }
      else
      {
	  if (mpz_cmp (t, effB2) >= 0)
	  {
	      outputf (OUTPUT_DEVVERBOSE, 
		       "Yes, works and gives higher B2 = %Zd\n", t);
	      P = tryP;
	      s_1 = trys_1;
	      s_2 = trys_2;
	      mpz_set (effB2, t);
	      while (l / 2 > s_1 && 
		     test_P (B2min, B2, m_1, P, l / 2 - s_1, effB2min, t))
	      {
		  l /= 2;
		  outputf (OUTPUT_DEVVERBOSE, 
			   "Reducing transform length to %ld, B2 = %Zd\n", 
			   l, t);
		  mpz_set (effB2, t);
	      }
	  }
	  else
	  {
	      outputf (OUTPUT_DEVVERBOSE, 
		       "Works, but does not give higher B2, %Zd <= %Zd\n",
		       t, effB2);
	  }
      }
    }

  /* Compute the correct values again */
  test_P (B2min, B2, m_1, P, l - s_1, effB2min, effB2);
  outputf (OUTPUT_DEVVERBOSE, "choose_P: final choice is: P = %lu, s_1 = %lu, "
	   "s_2 = %lu, l = %lu, m_1 = %Zd, effB2min = %Zd, effB2 = %Zd\n", 
	   P, s_1, s_2, l, m_1, effB2min, effB2);

  if (finalparams != NULL)
    {
      finalparams->P = P;
      finalparams->s_1 = s_1;
      finalparams->s_2 = s_2;
      finalparams->l = l;
      mpz_set (finalparams->m_1, m_1);
    }
  if (final_B2min != NULL)
    mpz_set (final_B2min, effB2min);
  if (final_B2 != NULL)
    mpz_set (final_B2, effB2);

  mpz_clear (m_1);
  mpz_clear (t);
  mpz_clear (B2l);
  mpz_clear (effB2);
  mpz_clear (effB2min);

  return P;
}



static void
list_output_poly (listz_t l, unsigned long len, int monic, int symmetric,
		  char *prefix, int verbosity)
{
  unsigned long i;

  if (prefix != NULL)
    outputf (verbosity, prefix);

  if (len == 0)
    {
      if (monic)
	outputf (verbosity, "1\n", len, len);
      else
	outputf (verbosity, "0\n", len);
      return;
    }

  if (monic)
    {
      if (symmetric)
	outputf (verbosity, "(x^%lu + x^-%lu) + ", len, len);
      else
	outputf (verbosity, "x^%lu + ", len);
    }
  for (i = len - 1; i > 0; i--)
    if (symmetric)
      outputf (verbosity, "%Zd * (x^%lu + x^-%lu) + ", l[i], i, i);
    else
      outputf (verbosity, "%Zd * x^%lu + ", l[i], i);
  outputf (verbosity, "%Zd\n", l[0]);
}


/* Multiply P[i] by r^{k(deg-i)}, for 0 <= i <= deg. Needs 3 entries in tmp. */
/* I.e., let P(x) = x^deg + \sum_{i=0}^{deg - 1} P[i] * x^i. The output is 
   R(x) = x^deg + \sum_{i=0}^{deg - 1} R[i] * x^i = r^(k deg) P(r^{-k} x). */
/* The input and output polynomials are monic and have the leading monomial
   implicit, i.e. not actually stored in the array of coefficients. */
/* Returns 0 if a modular inversion failed (in which case R is left 
   unchanged), 1 otherwise */

static int
list_scale_rev (listz_t R, listz_t S, mpz_t r, long k, unsigned long deg, 
		mpz_t modulus, listz_t tmp, 
		ATTRIBUTE_UNUSED const unsigned long tmplen)
{
  unsigned long i;

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


/* Multiply two reciprocal polynomials with coefficients in standard basis

   S_1(x) = S1[0] + sum_{1 \leq i \leq l1 - 1} S1[i] (x^i + x^{-i})
   S_2(x) = S2[0] + sum_{1 \leq i \leq l2 - 1} S2[i] (x^i + x^{-i})

   to the reciprocal polynomial

   R(x) = R[0] + sum_{1 \leq i \leq l1 + l2 - 2} R[i] (x^i + x^{-i}) 
        = S_1(x) * S_2(x)
*/

static void
list_mul_symmetric (listz_t R, const listz_t S1, const unsigned long l1, 
		    const listz_t S2, const unsigned long l2,
		    listz_t tmp, ATTRIBUTE_UNUSED const unsigned long tmplen)
{
  unsigned long i, lmax;
  const unsigned long dsum = l1 + l2 - 2; /* Half the degree of product */
  listz_t t1, t2, r;

  lmax = (l1 > l2) ? l1 : l2;

  ASSERT (tmplen >= 8 * lmax - 2 + list_mul_mem (2 * lmax - 1));

  if (l1 == 0 || l2 == 0)
    return;

  t1 = tmp;
  t2 = tmp + 2 * lmax - 1;
  r = tmp + 4 * lmax - 2;
  
  /* Full copy of S_1(x). S1 = [1,2,3,4] => t1 = [4,3,2,1,2,3,4]
     There are 2*l1 - 1 coefficients in monomial basis, which go in 
     t1[0 ... 2*l1-2]. We pad the high end with zeroes up to t1[2*lmax-2] */
  for (i = 0; i < l1; i++)
    mpz_set (t1[i], S1[l1 - 1 - i]); /* S[l1-1 ... 0] -> t1[0 ... l1-1] */
  for (i = 1; i < l1; i++)
    mpz_set (t1[l1 - 1 + i], S1[i]); /* S[1 ... l1-1] -> t1[l1 ... 2*l1-2] */
  for (i = 2 * l1 - 1; i <= 2 * lmax - 2; i++)
    mpz_set_ui (t1[i], 0UL);         /* t1[2*l1-1 ... 2*lmax-2] = 0 */
  
  /* Same for S_2(x) */
  for (i = 0; i < l2; i++)
    mpz_set (t2[i], S2[l2 - 1 - i]);
  for (i = 1; i < l2; i++)
    mpz_set (t2[l2 - 1 + i], S2[i]);
  for (i = 2 * l2 - 1; i <= 2 * lmax - 2; i++)
    mpz_set_ui (t2[i], 0UL);
  
  if (l1 == l2)
    {
      /* Simple case, the two polynomials are the same length. We can use
	 list_mul_high () */

      list_mul_high (r, t1, t2, 2 * lmax - 1, tmp + 8 * lmax - 2);
    }
  else
    {
      /* More difficult case, the lengths are different. We just do a full
	 multiply and take the coefficients we want from that. This is not
	 very efficient, but it'll only happen during building the polynomial
	 and for sets of odd cardinality, i.e. when the polynomials to be
	 multiplied are still quite small. The inefficiency of the code here 
	 does not realy matter too much. */

      list_mul (r, t1, 2 * lmax - 1, 0, t2, 2 * lmax - 1, 0, 
		tmp + 8 * lmax - 2);

      /* Now r/x^(dsum) is the product polynomial. It has degree 2*dsum and 
	 so has 2 * dsum + 1 coefficients in monomial basis, which reside in
	 r[0 ... 2 * sum] */

#ifdef WANT_ASSERT
      /* Check the lower terms are symmetric */
      for (i = 1; i <= dsum; i++)
	ASSERT (mpz_cmp (r[dsum - i], r[dsum + i]) == 0);
      
      /* Check the high terms are zero */
      for (i = 2 * dsum + 1; i <= 2 * lmax - 2; i++)
	ASSERT (mpz_sgn (r[i]) == 0);
#endif
    }
      
  for (i = 0; i <= dsum; i++)
    mpz_set (R[i], r[dsum + i]);
}


/* Multiply a (possibly monic) polynomial A of length k * len with a 
   (possibly monic) polynomial B of length len. R may be identical to A. */

static void
list_mul_blocks (listz_t R, const listz_t A, int monicA, const listz_t B, 
		 int monicB, const unsigned long len, const unsigned int k,
		 listz_t tmp, ATTRIBUTE_UNUSED const unsigned long tmplen)
{
  unsigned int j;
  
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


/* compute V_k(S), where V(x) is defined by V_k(X + 1/X) = X^k + 1/X^k */

static void
V (mpres_t R, const mpres_t S, const long k, mpmod_t modulus)
{
  mpres_t Vi, Vi1;
  unsigned long i, j, uk;

  if (k == 0)
    {
      mpres_set_ui (R, 2UL, modulus);
      return;
    }

  uk = labs (k);

  if (uk == 1)
    {
      mpres_set (R, S, modulus);
      return;
    }

  mpres_init (Vi, modulus);
  mpres_init (Vi1, modulus);

  for (j = 1; j <= uk / 2; j <<= 1);
  ASSERT ((uk & j) > 0);

  j >>= 1;
  i = 1;
  mpres_set (Vi, S, modulus);
  mpres_mul (Vi1, S, S, modulus);
  mpres_sub_ui (Vi1, Vi1, 2, modulus);

  while (j)
    {
      if ((uk & j) != 0)
	{
	  /* i' = 2i + 1.
	     V_{i'} = V_{2i + 1} = V_{i+1 + i} = V_{i+1} * V_{i} - V_1
	     V_{i'+1} = V_{2i + 2} = {V_{i+1}}^2 - 2. */
	  mpres_mul (Vi, Vi, Vi1, modulus);
	  mpres_sub (Vi, Vi, S, modulus);
	  mpres_mul (Vi1, Vi1, Vi1, modulus);
	  mpres_sub_ui (Vi1, Vi1, 2, modulus);
	  i = 2*i + 1;
	}
      else
	{
	  /* i' = 2i. 
	     V_{i'} = V_{2i} = {V_i}^2 - 2.
	     V_{i'+1} = V_{2i + 1} = V_{i+1 + i} = V_{i+1} * V_{i} - V_1 */
	  mpres_mul (Vi1, Vi, Vi1, modulus);
	  mpres_sub (Vi1, Vi1, S, modulus);
	  mpres_mul (Vi, Vi, Vi, modulus);
	  mpres_sub_ui (Vi, Vi, 2, modulus);
	  i = 2*i;
	}
      j >>= 1;
    }

  ASSERT (i == uk);
  mpres_set (R, Vi, modulus);

  mpres_clear (Vi, modulus);
  mpres_clear (Vi1, modulus);
}

/* For a given reciprocal polynomial 
   F(x) = f_0 + sum_{i=1}^{deg} f_i V_i(x+1/x),
   compute F(\gamma x)F(\gamma^{-1} x), with Q = \gamma + 1 / \gamma

   Hence, the result is 
   f_i * (gamma*x)^deg * f_i * (1/gamma*x)^deg + ... 
   + f_i * (gamma*x)^-deg * f_i * (1/gamma*x)^-deg
   = f_i * gamma^deg * x^deg * f_i * 1/gamma^deg * x^deg + ...
     f_i * gamma^-deg * x^-deg * f_i * 1/gamma^-deg * x^-deg
   = f_i * x^deg * f_i * x^deg + ... + f_i * x^-deg * f_i * x^-deg
   = f_i^2 * x^(2*deg) + ... + f_i^2 * x^(-2*deg)
*/

static void
list_scale_V (listz_t R, listz_t F, mpres_t Q, unsigned long deg,
	      mpmod_t modulus, listz_t tmp, const unsigned long tmplen)
{
  mpres_t Vi_1, Vi, Vt;
  unsigned long i;
  const listz_t G = tmp, H = tmp + 2 * deg + 1, newtmp = tmp + 4 * deg + 2;
  listz_t H_U;
  const unsigned long newtmplen = tmplen - 4 * deg - 2;
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
  
#ifdef WANT_ASSERT
  mpz_init (leading);
  mpz_mul (leading, F[deg], F[deg]);
  mpz_mod (leading, leading, modulus->orig_modulus);
#endif

  /* Generate V_1(Q)/2 ... V_{deg}(Q)/2, multiply by f_i to form coefficients 
     of G(x). Square the symmetric G(x) polynomial. */

  outputf (OUTPUT_TRACE, "list_scale_V: Q=%Zd, deg = %lu\n", Q, deg);
  list_output_poly (F, deg + 1, 0, 1, "list_scale_V: F(x) = ", 
		    OUTPUT_TRACE);

  mpres_init (Vi_1, modulus);
  mpres_init (Vi, modulus);
  mpres_init (Vt, modulus);
  mpres_set_ui (Vi_1, 1UL, modulus); /* Vi_1 = V_0(Q) / 2 = 1*/
  mpres_div_2exp (Vi, Q, 1, modulus); /* Vi = V_1(Q) = Q/2 */

  mpz_set (G[0], F[0]); /* G_0 = S_0 * V_0(Q)/2 = S_0 * 1 */
  outputf (OUTPUT_TRACE, "list_scale_V: G_%lu = %Zd\n", 0, G[0]);
  for (i = 1; i <= deg; i++)
    {
      /* Here, Vi = V_i(Q)/2, Vi_1 = V_{i-1}(Q)/2. */
      mpres_mul_z_to_z (G[i], Vi, F[i], modulus); /* G[i] = S_i * V_i(Q)/2 */
      outputf (OUTPUT_TRACE, 
	       "list_scale_V: G_%lu = F_%lu * V_%lu(Q)/2 = %Zd * %Zd = %Zd\n", 
	       i, i, i, F[i], Vi, G[i]);
      
      mpres_mul (Vt, Vi, Q, modulus);
      mpres_sub (Vt, Vt, Vi_1, modulus);
      mpres_set (Vi_1, Vi, modulus); /* Could be a swap */
      mpres_set (Vi, Vt, modulus); /* Could be a swap */
    }

  list_output_poly (G, deg + 1, 0, 1, "list_scale_V: G(x) = ", OUTPUT_TRACE);

  /* Now square the G polynomial in G[0 .. deg], put result in
     G[0 .. 2*deg] */

  /* Bugfix: ks_multiply() does not like negative coefficients. */

  for (i = 0; i <= deg; i++)
    if (mpz_sgn (G[i]) < 0)
      {
	mpz_add (G[i], G[i], modulus->orig_modulus);
	/* FIXME: make sure the absolute size does not "run away" */
	if (mpz_sgn (G[i]) < 0)
	  {
	    outputf (OUTPUT_ERROR, "list_scale_V: G[%lu] still negative\n", i);
	    mpz_mod (G[i], G[i], modulus->orig_modulus);
	  }
      }

  list_mul_symmetric (G, G, deg + 1, G, deg + 1, newtmp, newtmplen);
  
  list_output_poly (G, 2 * deg + 1, 0, 1, "list_scale_V: G(x)^2 == ", 
		    OUTPUT_TRACE);

  /* Generate U_1(Q)/2 ... U_deg(Q)/2, multpliy by S_i to form H. Convert H 
     to standard basis. Square the symmetic H polynomial. Multiply H^2 by
     (X + 1/X)^2 = X^2 + 2 + 1/X^2. Multiply that by (Q^2 - 4). */

  /* We'll reuse the Vi and Vi_1 variables here, but now they hold the 
     U_i(Q)/2 and U_{i-1}(Q)/2 values, respectively. */
  mpres_set_ui (Vi_1, 0UL, modulus); /* Vi_1 = U_0(Q) / 2 = 0 */
  mpres_set_ui (Vi, 1UL, modulus);
  mpres_div_2exp (Vi, Vi, 1, modulus); /* V_i = U_1(Q) / 2 = 1/2 */

  /* We later want to convert H in U_i basis to monomial basis. To do so,
     we'll need one list element below H_U[0], so H_U gets stored shifted
     up by one index */
  H_U = H - 1;

  /* H_U[i] = h_i =  F[i] * U_i(Q) / 2, for 1 <= i <= deg. H[0] is undefined
     and has no storage allocated (H_U[0] = H[-1]) */
  for (i = 1; i <= deg; i++)
    {
      /* Here, Vi = U_i(Q) / 2, Vi_1 = U_{i-1}(Q) / 2. */
      /* h_i = S_i * U_i(Q)/2 */
      mpres_mul_z_to_z (H_U[i], Vi, F[i], modulus);
      outputf (OUTPUT_TRACE, 
	       "list_scale_V: H_%lu (in U_i basis) = F_%lu * U_%lu(Q)/2 = %Zd * %Zd = %Zd\n", 
	       i, i, i, F[i], Vi, H_U[i]);
      
      mpres_mul (Vt, Vi, Q, modulus);
      mpres_sub (Vt, Vt, Vi_1, modulus);
      mpres_set (Vi_1, Vi, modulus); /* Could be a swap */
      mpres_set (Vi, Vt, modulus); /* Could be a swap */
    }

  /* Convert H_U to standard basis */

  /* We can do it in-place with H - 1 = H_U. */

  for (i = deg; i >= 3; i--)
    {
      mpz_add (H_U[i - 2], H_U[i - 2], H_U[i]);
      /* mpz_set (H[i - 1], H_U[i]); A no-op, since H - 1 = H_U. */
    }
  
  /* U_2(X+1/X) = (X^2 - 1/X^2)/(X-1/X) = X+1/X = V_1(X+1/X),
     so no addition occures here */
  /* if (deg >= 2)
     mpz_set (H[1], H_U[2]); Again, a no-op. */
  
  /* U_1(X+1/X) = 1, so this goes to coefficient of index 0 in std. basis */
  /* mpz_set (H[0], H_U[1]); Another no-op. */
  
  /* Now H[0 ... deg-1] contains the deg coefficients in standard basis
     of symmetric H(X) of degree 2*deg-2. */
  
  list_output_poly (H, deg, 0, 1, "list_scale_V: H(x) = ", OUTPUT_TRACE);

  /* Square the symmetric H polynomial of degree 2*deg-2 (i.e. with deg 
     coefficents in standard basis in H[0 ... deg-1]) */

  /* Bugfix: ks_multiply() does not like negative coefficients. */

  for (i = 0; i <= deg; i++)
    if (mpz_sgn (H[i]) < 0)
      mpz_mod (H[i], H[i], modulus->orig_modulus);

  list_mul_symmetric (H, H, deg, H, deg, newtmp, newtmplen);

  /* Now there are the 2*deg-1 coefficients in standard basis of a 
     symmetric polynomial of degree 4*deg - 4 in H[0 ... 2*deg-2] */

  list_output_poly (H, 2*deg - 1, 0, 1, "list_scale_V: H(x)^2 == ", 
		    OUTPUT_TRACE);

  /* Multiply by Q^2-4 */
  ASSERT (newtmplen >= 2);
  mpres_mul (Vt, Q, Q, modulus);
  mpres_sub_ui (Vt, Vt, 4, modulus);
  for (i = 0; i <= 2 * deg - 2; i++)
    mpres_mul_z_to_z (H[i], Vt, H[i], modulus);
  list_output_poly (H, 2 * deg - 1, 0, 1, "list_scale_V: H(x)^2*(Q^2-4) == ", 
		    OUTPUT_TRACE);

  /* Multiplying by (X - 1/X)^2 = X^2 - 2 + 1/X^2 */

  if (deg == 1)
    {
      /* H(X) has degree 2*deg-2 = 0, so H(X) = h_0
	 H(X) * (X - 1/X)^2 = -2 h_0 + h_0 V_2(Y)  */
      mpz_mul_2exp (newtmp[0], H[0], 1UL);
      mpz_neg (newtmp[0], newtmp[0]);
      mpz_set_ui (newtmp[1], 0UL);
      mpz_set (newtmp[2], H[0]);
    }
  else if (deg == 2)
    {
      /* H(X) has degree 2*deg-2 = 2, , so 
	 H(X) = h_0 + h_1 (X+1/X) + h_2 (X^2+1/X^2)

	 H(X) * (X - 1/X)^2 =
	 2*(h_2 - h_0) - h_1 * V_1(Y) + (h_0 - 2*h_2) * V_2(Y) + 
	 h_1 * V_3(Y) + h_2 * V_4(Y)
      */
      mpz_sub (newtmp[0], H[2], H[0]);
      mpz_mul_2exp (newtmp[0], newtmp[0], 1UL); /* 2*(h_2 - h_0) */
      mpz_neg (newtmp[1], H[1]);                /* -h_1 */
      mpz_mul_2exp (newtmp[2], H[2], 1UL);
      mpz_sub (newtmp[2], H[0], newtmp[2]);     /* h_0 - 2*h_2 */
      mpz_set (newtmp[3], H[1]);
      mpz_set (newtmp[4], H[2]);
    }
  else
    {
      /* Let H(X) = h_0 + \sum_{i=1}^{n} h_i V_i(Y), Y = X+1/X. Then
	 (x - 1/x)^2 H(X) = 
	 2(-h_0 + h_2) +
	 (- h_1 + h_3) V_1(Y) +
	 \sum_{i=2}^{n-2} (h_{i-2} - 2h_i + h_{i+2}) V_i(Y) +
	 (h_{n-3} - 2h_{n-1}) V_{n-1}(Y) +
	 (h_{n-2} - 2h_n) V_n(Y) +
	 h_{n-1} V_{n+1}(Y) +
	 h_n V_{n+2}(Y)
	 
	 In our case, n = 2 * deg - 2
      */
      mpz_sub (newtmp[0], H[2], H[0]);
      mpz_mul_2exp (newtmp[0], newtmp[0], 1UL); /* t[0] = 2*(h_0 + h_2) */
      
      mpz_sub (newtmp[1], H[3], H[1]); /* t[1] = -h_1 + h_3 */
      
      for (i = 2; i <= 2 * deg - 4; i++)
	{
	  mpz_add (newtmp[i], H[i-2], H[i+2]);
	  mpz_sub (newtmp[i], newtmp[i], H[i]); /* t[i] = h_{i-2}-2h_i+h_{i+2} */
	  mpz_sub (newtmp[i], newtmp[i], H[i]); /* for 2 <= i <= n-2 */
	}
      
      mpz_mul_2exp (newtmp[2 * deg - 3], H[2 * deg - 3], 1UL);
      mpz_sub (newtmp[2 * deg - 3], H[2 * deg - 5], newtmp[2 * deg - 3]); 
      /* t[n-1] = h_{n-3} - 2h_{n-1} */
      
      mpz_mul_2exp (newtmp[2 * deg - 2], H[2 * deg - 2], 1UL);
      mpz_sub (newtmp[2 * deg - 2], H[2 * deg - 4], newtmp[2 * deg - 2]);
      /* t[n] = h_{n-2} - 2h_n */
      
      mpz_set (newtmp[2 * deg - 1], H[2 * deg - 3]); /* t[n+1] = h_{n-1} */
      mpz_set (newtmp[2 * deg], H[2 * deg - 2]);  /* t[n+2] = h_n */
    }
  
  for (i = 0; i <= 2 * deg; i++)
    mpz_set (H[i], newtmp[i]);

  /* Now H[0 ... 2*deg] contains the 2*deg+1 coefficients in standard basis
     of a degree 4*deg symmetric polynomial */

  /* Subtract the two polynomials, reduce mod modulus and store in R */

  for (i = 0; i <= 2 * deg; i++)
    {
      mpz_sub (G[i], G[i], H[i]);
      mpz_mod (R[i], G[i], modulus->orig_modulus);
    }

#ifdef WANT_ASSERT
  mpz_mod (R[2 * deg], R[2 * deg], modulus->orig_modulus);
  ASSERT (mpz_cmp (leading, R[2 * deg]) == 0);
  mpz_clear (leading);
#endif

  mpres_clear (Vt, modulus);
  mpres_clear (Vi, modulus);
  mpres_clear (Vi_1, modulus);
}


#ifdef WANT_ASSERT
/* Check if l is an (anti-)symmetric, possibly monic, polynomial. 
   Returns -1 if it is (anti-)symmetric, or the smallest index i where 
   l[i] != l[len - 1 + monic - i])
   If anti == 1, the list is checked for symmetry, if it is -1, for
   antisymmetry.
   This function is used only if assertions are enabled.
*/

static long int
list_is_symmetric (listz_t l, unsigned long len, int monic, int anti, 
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
		const unsigned long n, const int monic, const mpz_t modulus, 
		listz_t tmp)
{
  unsigned long i;

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

/* Build a polynomial F from the arithmetic progressions in sets.
   The roots of F will be r^i for all i in the sumset of the sets.
   Returns the degree of the resulting polynomial. */

static int
poly_from_sets (mpz_t *F, mpz_t r, long *sets, unsigned long setsize,
		mpz_t *tmp, const unsigned long tmplen, mpz_t modulus)
{

  unsigned long l, c; /* Cardinality of this set */
  unsigned long i, j, deg;
  long k;
  
  if (setsize == 0)
    {
      mpz_set_si (F[0], -1L);
      return 1;
    }

  ASSERT (sets[0] >= 0);
  l = (unsigned long) sets[0];
  ASSERT (l > 1); /* Empty sets indicate a bug somewhere */
  ASSERT (l <= setsize);
  c = l - 1;
  deg = poly_from_sets (F, r, sets + l, setsize - l, tmp, tmplen, modulus);

  if (test_verbose (OUTPUT_TRACE))
    list_output_poly (F, deg, 1, 0, "poly_from_sets: F = ", OUTPUT_TRACE);

  j = c;
  for (i = 0; i < c; i++)
    {
      k = sets[i + 1];
      if (k != 0)
	{
	  outputf (OUTPUT_TRACE, 
		   "list_scale_rev(F+%lu, F, %Zd, %ld, %lu) = ", 
		   (j-1)*deg, r, k, deg);
	  if (list_scale_rev (F + (--j) * deg, F, r, k, deg, modulus, tmp, 
			      tmplen) == 0)
	    abort (); /* FIXME, deal with factor! */
	  if (test_verbose (OUTPUT_TRACE))
	    list_output_poly (F + j, deg, 1, 0, "poly_from_sets: ", 
			      OUTPUT_TRACE);
	}
    }
  /* The following assert assumes that the sets are symmetric around 0,
     that is, they contain 0 iff the cardinality is odd.
     The set expanded so far is not symmetric around 0 if a set of 
     cardinality 1 appeared somewhere. Assuming that the sets are sorted,
     this can only happen at the top of the recursion, so if this one
     doesn't have cardinality one, none of the already processed ones did. */

  ASSERT(c == 1 || j == (c & 1)); /* Cardinality odd <=> one less list_scale */

  ASSERT (c == 1 || (c == 2  && tmplen >= 2 * deg + list_mul_mem (deg))
	  || (c > 2  && tmplen >= 3 * deg + list_mul_mem (deg)));

  for (i = 1; i < c; i++)
    {
      list_mul_blocks (F, F, 1, F + i * deg, 1, deg, i, tmp, tmplen);
      list_mod (F, F, (i + 1) * deg, modulus);
      if (test_verbose (OUTPUT_TRACE))
	list_output_poly (F, (i + 1) * deg, 1, 0, "poly_from_sets: ", 
			  OUTPUT_TRACE);
    }

#if defined(WANT_ASSERT)
  /* Test that the polynomial is symmetric if the degree is even, or anti-
     symmetric if it is odd */
  if (c != 1)
    ASSERT (list_is_symmetric (F, c * deg, 1, (c * deg) % 2 == 0 ? 1 : -1, 
			       modulus, tmp[0]) == -1);
#endif
  
  return c * deg;
}

/* Build a polynomial with roots r^i, i in the sumset of the sets in "sets".
   The parameter Q = r + 1/r. This code uses the fact that the polynomials 
   are symmetric. Requires that the first set in sets has cardinality 2,
   all sets must be symmetric around 0. The resulting polynomial of degree 
   2*d is F(x) = f_0 + \sum_{1 <= i <= d} f_i (x^i + 1/x^i). The coefficient
   f_i is stored in F[i], which therefore needs d+1 elements. */

static int
poly_from_sets_V (listz_t F, const mpres_t Q, const long *sets, 
		  unsigned long setsize, listz_t tmp, 
		  const unsigned long tmplen, mpmod_t modulus)
{
  unsigned long c, deg, lastidx, i;
  const long *curset;
  mpres_t Qt;
  
  ASSERT_ALWAYS (sets[0] == 3); /* Check that the cardinality of first set 
				   is 2 */
  ASSERT_ALWAYS (sets[1] == -sets[2]); /* Check that first set is symmetric 
					  around 0 */

  if (test_verbose (OUTPUT_TRACE))
    {
      mpz_t t;
      mpz_init (t);
      mpres_get_z (t, Q, modulus);
      outputf (OUTPUT_TRACE, 
	       "poly_from_sets_V (F, Q = %Zd, sets, setsize = %lu)\n", 
	       t, setsize);
      mpz_clear (t);
    }

  mpres_init (Qt, modulus);
  
  V (Qt, Q, sets[1], modulus); /* Qt = V_k(Q) */

  mpres_get_z (F[0], Qt, modulus);
  mpz_neg (F[0], F[0]);
  mpz_set_ui (F[1], 1UL);
  deg = 1;
  /* Here, F(x) = (x - r^{k_1})(x - r^{-k_1}) / x = 
                  (x^2 - x (r^{k_1} + r^{-k_1}) + 1) / x =
		  (x + 1/x) - V_{k_1}(r + 1/r) */

  sets += 3;
  setsize -= 3;

  outputf (OUTPUT_DEVVERBOSE, " (processing set of size");

  while (setsize > 0)
    {
      /* Assuming the sets are sorted in ascending order, we process them
	 back-to-front, so the sets of cardinality 2 are processed last. */
      for (lastidx = 0; 
	   lastidx + sets[lastidx] < setsize; /* More sets after this one? */
	   lastidx += sets[lastidx]); /* If yes, skip this one */
      
      /* Process this set. We assume it is either of cardinality 2, or of 
	 odd cardinality */
      c = sets[lastidx] - 1;
      curset = sets + lastidx + 1;
      outputf (OUTPUT_DEVVERBOSE, " %lu", c);

      if (c == 2)
	{
	  ASSERT_ALWAYS (curset[0] == -curset[1]); /* Check it's symmetric */
	  V (Qt, Q, curset[0], modulus);
	  list_scale_V (F, F, Qt, deg, modulus, tmp, tmplen);
	  deg *= 2;
	  ASSERT (mpz_cmp_ui (F[deg], 1UL) == 0); /* Check it's monic */
	}
      else
	{
	  ASSERT_ALWAYS (c % 2 == 1);
	  /* Generate the F(Q^{k_i} * X)*F(Q^{-k_i} * X) polynomials.
	     Each is symmetric of degree 2*deg, so each has deg+1 coeffients
	     in standard basis. */
	  for (i = 0; i < (c - 1) / 2; i++)
	    {
	      ASSERT (curset[i] == -curset[c-1-i]); /* Check it's symmetric */
	      V (Qt, Q, curset[i], modulus);
	      ASSERT (mpz_cmp_ui (F[deg], 1UL) == 0); /* Check it's monic */
	      list_scale_V (F + (2 * i + 1) * (deg + 1), F, Qt, deg, modulus, 
			    tmp, tmplen);
	      ASSERT (mpz_cmp_ui (F[(2 * i + 1) * (deg + 1) + 2 * deg], 1UL) == 0); /* Check it's monic */
	    }
	  /* Multiply the polynomials together */
	  for (i = 0; i < (c - 1) / 2; i++)
	    {
	      /* So far, we have the product 
		 F(X) * F(Q^{k_j} * X) * F(Q^{-k_j} * X), 1 <= j <= i,
		 at F. This product has degree 2 * deg + i * 4 * deg, that is
		 (2 * i + 1) * 2 * deg, which means (2 * i + 1) * deg + 1
		 coefficients in F[0 ... (i * 2 + 1) * deg]. */
	      ASSERT (mpz_cmp_ui (F[(2 * i + 1) * deg], 1UL) == 0);
	      ASSERT (mpz_cmp_ui (F[(2 * i + 1) * (deg + 1) + 2*deg], 1UL) == 0);
	      list_output_poly (F, (2 * i + 1) * deg + 1, 0, 1, 
				"poly_from_sets_V: Multiplying ", 
				OUTPUT_TRACE);
	      list_output_poly (F + (2 * i + 1) * (deg + 1), 2 * deg + 1, 0, 1,
				" and ", OUTPUT_TRACE);
	      list_mul_symmetric (F, 
				  F, (2 * i + 1) * deg + 1, 
				  F + (2 * i + 1) * (deg + 1), 2 * deg + 1, 
				  tmp, tmplen);
	      list_mod (F, F, (2 * i + 3) * deg + 1, modulus->orig_modulus);
	      list_output_poly (F, (2 * i + 3) * deg + 1, 0, 1, " = ",
				OUTPUT_TRACE);
	      ASSERT (mpz_cmp_ui (F[(2 * i + 3) * deg], 1UL) == 0);
	    }
	  deg *= c;
	}

      ASSERT(setsize >= (unsigned long) sets[lastidx]);
      setsize -= sets[lastidx];
    }

  mpres_clear (Qt, modulus);
  outputf (OUTPUT_DEVVERBOSE, ")");

  return deg;
}


/* Compute g_i = x_0^{M-i} * r^{(M-i)^2} for 0 <= i < l. 
   x_0 = b_1^{2*k_2 + (2*m_1 + 1) * P}. r = b_1^P. */

static void
pm1_sequence_g (mpzspv_t g_ntt, const mpres_t b_1, const unsigned long P, 
		const long M, const unsigned long l, const mpz_t m_1, 
		const long k_2, mpmod_t modulus, mpzspm_t ntt_context)
{
  mpres_t r[3], x_0, x_Mi;
  mpz_t t;
  unsigned long i;
  long timestart, timestop;

  mpz_init (t);
  mpres_init (r[0], modulus);
  mpres_init (r[1], modulus);
  mpres_init (r[2], modulus);
  mpres_init (x_0, modulus);
  mpres_init (x_Mi, modulus);

  outputf (OUTPUT_DEVVERBOSE, "sequence_g (g, b_1, P = %lu, M = %ld, l = %lu, "
	   "m_1 = %Zd, k_2 = %ld, modulus)\n", P, M, l, m_1, k_2);
  outputf (OUTPUT_VERBOSE, "Computing g_i");
  if (test_verbose (OUTPUT_TRACE))
    { 
      mpres_get_z (t, b_1, modulus);
      outputf (OUTPUT_TRACE, "\n/* pm1_sequence_g */ N = %Zd; "
	       "b_1 = Mod(%Zd, N); /* PARI */\n", modulus->orig_modulus, t);
      outputf (OUTPUT_TRACE, "/* pm1_sequence_g */ P = %lu; M = %lu; "
	       "m_1 = %Zd; /* PARI */\n", P, M, m_1);
      outputf (OUTPUT_TRACE, 
	       "/* pm1_sequence_g */ r = b_1^P; /* PARI */\n");
      outputf (OUTPUT_TRACE, "/* pm1_sequence_g */ x_0 = "
	       "b_1^(2*%ld + (2*m_1 + 1)*P); /* PARI */\n", k_2);
    }
  timestart = cputime ();

  /* We use (M-(i+1))^2 = (M-i)^2 + 2*(-M+i) + 1 */
  mpz_set_ui (t, P);
  mpres_pow (r[0], b_1, t, modulus);     /* r[0] = b_1^P = r */
  if (test_verbose (OUTPUT_TRACE))
    {
      mpres_get_z (t, r[0], modulus);
      outputf (OUTPUT_TRACE, "/* pm1_sequence_g */ r == %Zd /* PARI C */\n", t);
    }
  
  /* FIXME: This is a huge mess, clean up some time */

  mpz_set_ui (t, M);
  mpz_neg (t, t);
  mpz_mul_2exp (t, t, 1UL);
  mpz_add_ui (t, t, 1UL);
  mpres_pow (r[1], r[0], t, modulus);    /* r[1] = r^{2(-M+i)+1}, i = 0 */
  mpz_set_ui (t, M);
  mpz_mul_ui (t, t, M);
  mpres_pow (r[2], r[0], t, modulus);    /* r[2] = r^{(M-i)^2}, i = 0 */
  mpres_mul (r[0], r[0], r[0], modulus); /* r[0] = r^2 */

  mpz_mul_2exp (t, m_1, 1UL);
  mpz_add_ui (t, t, 1UL);
  mpz_mul_ui (t, t, P);
  mpz_add_si (t, t, k_2);
  mpz_add_si (t, t, k_2);
  outputf (OUTPUT_TRACE, 
	   "/* pm1_sequence_g */ 2*%ld + (2*%Zd + 1)*P == %Zd /* PARI C */\n", 
           k_2, m_1, t);

  mpres_pow (x_0, b_1, t, modulus);  /* x_0 = b_1^{2*k_2 + (2*m_1 + 1)*P} */
  if (test_verbose (OUTPUT_TRACE))
    {
      mpres_get_z (t, x_0, modulus);
      outputf (OUTPUT_TRACE, "/* pm1_sequence_g */ x_0 == %Zd /* PARI C */\n", 
	       t);
    }
  
  mpz_set_ui (t, M);
  mpres_pow (x_Mi, x_0, t, modulus); /* x_Mi = x_0^{M-i}, i = 0 */

  mpres_invert (x_0, x_0, modulus);  /* x_0 := x_0^{-1} now */
  mpres_mul (r[1], r[1], x_0, modulus); /* r[1] = x_0^{-1} * r^{-2M+1} */
  
  mpres_mul (r[2], r[2], x_Mi, modulus); /* r[2] = x_0^M * r^{M^2} */
  mpres_get_z (t, r[2], modulus);
  outputf (OUTPUT_TRACE, "/* pm1_sequence_g */ g_0 = %Zd; /* PARI */\n", t);
  mpzspv_from_mpzv (g_ntt, 0, &t, 1, ntt_context);

  /* So here we have for i = 0
     r[2] = x_0^(M-i) * r^{(M-i)^2}
     r[1] = x_0^{-1} * r^{2(-M+i)+1}
     r[0] = r^2
     t = r[2]
  */

  for (i = 1; i < l; i++)
    {
      mpres_mul_z_to_z (t, r[1], t, modulus);
      /* Hence t = x_0^(M-i) * r^{(M-i)^2} * x_0^{-1} * r^{2*(-M+i)+1} 
                 = x_0^(M-(i+1)) * r^{(M-(i+1))^2} */
      mpres_mul (r[1], r[1], r[0], modulus);
      /* Hence r[1] = x_0^{-1} * r^{2(-M+i)+1} * r^2 
                    = x_0^{-1} * r^{2(-M+(i+1))+1} */
      outputf (OUTPUT_TRACE, "/* pm1_sequence_g */ g_%lu = %Zd; /* PARI */\n", 
	       i, t);
      mpzspv_from_mpzv (g_ntt, i, &t, 1, ntt_context);
    }

  mpres_clear (r[0], modulus);
  mpres_clear (r[1], modulus);
  mpres_clear (r[2], modulus);
  mpres_clear (x_0, modulus);
  mpres_clear (x_Mi, modulus);
  mpz_clear (t);

  timestop = cputime ();
  outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);
  
  if (test_verbose (OUTPUT_TRACE))
    {
      for (i = 0; i < l; i++)
        {
          outputf (OUTPUT_TRACE, "/* pm1_sequence_g */ g_%lu == "
		   "x_0^(M - %lu) * r^((M - %lu)^2) /* PARI C */\n", i, i, i);
        }
      outputf (OUTPUT_TRACE, "/* pm1_sequence_g */ g(x) = g_0");
      for (i = 1; i < l; i++)
        outputf (OUTPUT_TRACE, " + g_%lu * x^%lu", i, i);
      outputf (OUTPUT_TRACE, " /* PARI */\n");
    }
}


/* Compute h_j = r^(-j^2) * f_j for 0 <= j < d as described in section 9 
   of the paper. h == f is ok. */

static void 
pm1_sequence_h (mpzspv_t h, mpz_t *f, const mpres_t r, const unsigned long d, 
		mpmod_t modulus, mpzspm_t ntt_context)
{
  mpres_t invr;  /* r^{-1} */
  mpres_t fd[3]; /* finite differences table for r^{-i^2}*/
  mpz_t t[1]; /* the h_j value as an mpz_t */
  unsigned long j;
  long timestart, timestop;

  outputf (OUTPUT_VERBOSE, "Computing h_i = r^(-i^2) * f_i, 0 <= i < %lu", d);
  if (test_verbose (OUTPUT_TRACE))
    {
      mpz_t tmp;
      mpz_init (tmp);
      mpres_get_z (tmp, r, modulus);
      outputf (OUTPUT_TRACE, 
	       "\n/* pm1_sequence_h */ N = %Zd; r = Mod(%Zd, N); /* PARI */\n", 
	       modulus->orig_modulus, tmp);
      mpz_clear (tmp);
    }

  timestart = cputime ();

  mpres_init (invr, modulus);
  mpres_init (fd[0], modulus);
  mpres_init (fd[1], modulus);
  mpres_init (fd[2], modulus);
  mpz_init (t[0]);

  mpres_invert (invr, r, modulus); /* invr = r^{-1}. FIXME: test for failure,
				      even if it is theoretically impossible */

  /* We have (n + 1)^2 = n^2 + 2n + 1. For the finite differences we'll need
     2, 2n+1, n^2. Init for n = 0. */

  mpres_mul (fd[0], invr, invr, modulus); /* fd[0] = r^{-2} */
  mpres_set (fd[1], invr, modulus);       /* fd[1] = r^{-1} */
  /* mpres_set_ui (fd[2], 1UL, modulus); fd[2] = r^0. We set fd[2] below */

  /* For j = 0, we have h_j = f_j */
  if (d > 0)
    {
      outputf (OUTPUT_TRACE, "/* pm1_sequence_h */ h_0 = f_0; /* PARI */\n");
      mpzspv_from_mpzv (h, 0, f, 1, ntt_context);
    }

  /* Do j = 1 */
  if (d > 1)
    {
      mpres_set (fd[2], fd[1], modulus);        /* fd[2] = r^{-1} */
      mpres_mul (fd[1], fd[1], fd[0], modulus); /* fd[1] = r^{-3} */
      mpres_mul_z_to_z (t[0], fd[2], f[1], modulus);
      outputf (OUTPUT_TRACE,
	       "/* pm1_sequence_h */ h_1 = %Zd; /* PARI */\n", t);
      mpzspv_from_mpzv (h, 1, t, 1, ntt_context);
    }
  
  /* Do the rest */
  for (j = 2; j < d; j++)
    {
      mpres_mul (fd[2], fd[2], fd[1], modulus); /* fd[2] = r^{-j^2} */
      mpres_mul (fd[1], fd[1], fd[0], modulus); /* fd[1] = r^{-2*j-1} */
      mpres_mul_z_to_z (t[0], fd[2], f[j], modulus);
      outputf (OUTPUT_TRACE, 
	       "/* pm1_sequence_h */ h_%lu = %Zd; /* PARI */\n", j, t);
      mpzspv_from_mpzv (h, j, t, 1, ntt_context);
    }
  
  timestop = cputime ();
  outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);
  
  mpres_clear (fd[2], modulus);
  mpres_clear (fd[1], modulus);
  mpres_clear (fd[0], modulus);
  mpres_clear (invr, modulus);
  mpz_clear (t[0]);

  if (test_verbose (OUTPUT_TRACE))
    {
      for (j = 0; j < d; j++)
	outputf (OUTPUT_TRACE, "/* pm1_sequence_h */ h_%lu == "
		   "f_%lu * r^(-%lu^2) /* PARI C */\n", j, j, j);
      outputf (OUTPUT_TRACE, "/* pm1_sequence_h */ h(x) = h_0");
      for (j = 1; j < d; j++)
        outputf (OUTPUT_TRACE, " + h_%lu * (x^%lu + x^(-%lu))", j, j, j);
      outputf (OUTPUT_TRACE, " /* PARI */\n");
    }
}



/* Build polynomial F(x) with roots X^i for i covering all the residue classes
   coprime to beta. F must have space for eulerphi(beta) coefficients.
   method can be 0, 1 or 2, which mean: 0 = old way of computing all the 
   roots and doing a product tree, 1 = using recursive expansion of polynomial
   *without* Chebychev polynomials to utilize symmetry, 2 = using recursive 
   expansion of polynomial *with* Chebychev polynomials */

ATTRIBUTE_UNUSED static void 
pm1_build_poly_F (mpz_t *F, const mpres_t X, mpmod_t modulus, 
		  const unsigned long beta, const int method, 
		  const mpz_t i0, const unsigned long nr, 
		  const unsigned long blocks, const unsigned long tmplen, 
		  mpz_t *tmp)
{
  long timestart, timestop;
  long *sets = NULL;
  unsigned long setsize = 0UL, i;
  mpz_t mt;
  const unsigned long dF = eulerphi (beta);
  
  ASSERT (0 <= method && method <= 2);
  
  mpz_init (mt);
  
  if (method == 1 || method == 2)
    {
      sets = get_factored_sorted_sets (&setsize, beta);
      
      mpz_mul_ui (mt, i0, beta);
      mpz_add_si (mt, mt, sum_sets_minmax (sets, setsize, 1));
      outputf (OUTPUT_VERBOSE, "Effective B2min = %Zd\n", mt);
      mpz_add_ui (mt, i0, blocks * nr);
      mpz_mul_ui (mt, mt, beta);
      mpz_sub_si (mt, mt, sum_sets_minmax (sets, setsize, -1));
      outputf (OUTPUT_VERBOSE, "Effective B2max = %Zd\n", mt);
    }
  else
    {
      mpz_mul_ui (mt, i0, beta);
      outputf (OUTPUT_VERBOSE, "Effective B2min ~= %Zd\n", mt);
      mpz_add_ui (mt, i0, blocks * nr);
      mpz_mul_ui (mt, mt, beta);
      outputf (OUTPUT_VERBOSE, "Effective B2max ~= %Zd\n", mt);
    }

  if (method == 1 || (method == 2 && sets[0] == 2))
    {
      /* poly_from_sets_V () can't handle set of cardinality 1, so use 
	 poly_from_sets () in that case */
      
      outputf (OUTPUT_VERBOSE, 
	       "Computing F from factored set of units mod %lu", beta);
      outputf (OUTPUT_DEVVERBOSE, "\n");
      
      timestart = cputime ();
      
      mpres_get_z (mt, X, modulus); /* poly_from_sets() expects mpz_t atm */
      i = poly_from_sets (F, mt, sets, setsize, tmp, tmplen, 
			modulus->orig_modulus);
      ASSERT_ALWAYS (i == dF);

      timestop = cputime ();
      outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);

#if defined(WANT_ASSERT)
      if (sets[0] != 2) /* Unless the first set has one element, 
			   the polynomial should be symmetric */
	{
	  long i = list_is_symmetric (F, dF, 1, 1, modulus->orig_modulus, mt);
	  if (i != -1)
	    {
	      outputf (OUTPUT_ERROR, 
		       "Polynomial not symmetric! F[%ld] != F[%ld]\n", 
		       i, dF - i);
	      list_output_poly (F, dF, 1, 0, "F(x) = ", OUTPUT_ERROR);
	      outputf (OUTPUT_ERROR, "Factored sets: ");
	      print_sets (OUTPUT_ERROR, sets, setsize);
	      abort ();
	    }
	}
#endif
    }
  else if (method == 2)
    {
      /* Use poly_from_sets_V () */
      mpres_t X1X;

      outputf (OUTPUT_VERBOSE, "Computing F from symmetric factored set "
	       "of units mod %lu", beta);
      if (test_verbose (OUTPUT_DEVVERBOSE))
	outputf (OUTPUT_VERBOSE, "\n");

      timestart = cputime ();

      /* First compute X + 1/X */
      mpres_init (X1X, modulus);
      mpres_invert (X1X, X, modulus);
      mpres_add (X1X, X1X, X, modulus);

      i = poly_from_sets_V (F, X1X, sets, setsize, tmp, tmplen, 
			    modulus);
      ASSERT(2 * i == dF);
      ASSERT(mpz_cmp_ui (F[i], 1UL) == 0);
      
      mpres_clear (X1X, modulus);

      /* Make symmetric copy. The leading 1 monomial will not get stored,
         but will be implicit from here on. */

      for (i = 0; i < dF / 2; i++)
	mpz_set (F[dF / 2 + i], F[i]); /* F[dF/2]:=f_0, F[dF-1]:=f_{dF/2-1} */
      
      for (i = 1; i < dF / 2; i++)
	mpz_set (F[i], F[dF - i]); /* F[1] := F[dF - 1] = f_{dF / 2 - 1}, 
				      F[dF / 2 - 1] := F[dF / 2 + 1] = f_1 */
      mpz_set_ui (F[0], 1UL);
      
      timestop = cputime ();
      outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);
    }
  else /* method == 0 */
    {
      /* Build F the old way, computing all the roots and doing a 
	 product tree */

      mpres_t Xs, Xi;
      const unsigned long s = 2 - beta % 2; /* 2 if beta is even, 1 if odd */
      unsigned long j;

      outputf (OUTPUT_VERBOSE, "Computing roots of F");
      outputf (OUTPUT_TRACE, "\n"); /* So the "Computing" does not interfere */
      
      timestart = cputime ();
      mpres_init (Xs, modulus);
      mpres_init (Xi, modulus);

      if (s == 1) /* Xs = X^s */
	mpres_set (Xs, X, modulus);
      else
	mpres_mul (Xs, X, X, modulus);

      mpres_set (Xi, X, modulus);    /* Xi = X^i for i = 1 */

      /* Prepare polynomial F(x), which is monic of degree dF. The leading
	 monomial is not stored. */
      /* Put in F[0 .. dF-1] the values of X^i, 1<=i<beta, gcd(i, beta) == 1 */
      for (i = 1, j = 0; i < beta; i += s)
	{
	  if (lgcd (i, beta) == 1)
	    mpres_get_z (F[j++], Xi, modulus);
	  mpres_mul (Xi, Xi, Xs, modulus);
	}
      
      ASSERT_ALWAYS (j == dF);
      mpres_clear (Xs, modulus);
      mpres_clear (Xi, modulus);
      timestop = cputime ();
      outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);
      
      /* Print the polynomial in linear factors form */
      outputf (OUTPUT_TRACE, "F(x) = ");
      for (j = 0; j < dF - 1; j++)
	outputf (OUTPUT_TRACE, "(x - %Zd) * ", F[j]);
      outputf (OUTPUT_TRACE, "(x - %Zd); /* PARI %ld */\n", F[dF - 1], 
	       pariline++);
      
      /* Multiply all the (x - f_i) to form F(x) in monomial basis */
      outputf (OUTPUT_VERBOSE, "Building F from its roots");
      timestart = cputime ();
      PolyFromRoots (F, F, dF, tmp, modulus->orig_modulus);
      timestop = cputime ();
      outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);
    }

  mpz_clear (mt);
  if (method == 1 || method == 2)
    free (sets);

  /* Print the final polynomial in monomial form */
  if (test_verbose (OUTPUT_TRACE))
    list_output_poly (F, dF, 1, 0, "F(x) == ", OUTPUT_TRACE);
}


static int 
make_S_1_S_2 (long **S_1, unsigned long *S_1_size, long **S_2, 
	      const faststage2_param_t *params)
{
  unsigned long i;

  *S_1 = get_factored_sorted_sets (S_1_size, params->P);
  if (*S_1 == NULL)
    return ECM_ERROR;
  ASSERT (sum_sets_minmax (*S_1, *S_1_size, 1) == (long) maxS (params->P));
  *S_2 = malloc ((params->s_2 + 1) * sizeof (long));
  if (*S_2 == NULL)
    {
      free (*S_1);
      return ECM_ERROR;
    }
  if (params->s_2 == 1)
    {
      (*S_2)[0] = 2L;
      (*S_2)[1] = 0L; /* A set containing only 0 */
    }
  else
    {
      /* Expand the sum of sets for S_2 into a single set */
      long *factored_S_2;
      unsigned long factored_S_2_size;
      
      factored_S_2_size = extract_sets (NULL, *S_1, *S_1_size, params->s_2);
      factored_S_2 = malloc (factored_S_2_size * sizeof (long));
      if (factored_S_2 == NULL)
	{
	  free (*S_1);
	  free (*S_2);
	  return ECM_ERROR;
	}
      extract_sets (factored_S_2, *S_1, *S_1_size, params->s_2);
      *S_1_size -= factored_S_2_size;
      (*S_2)[0] = params->s_2 + 1;
      sum_sets ((*S_2) + 1, factored_S_2, factored_S_2_size, 0L);
      free (factored_S_2);
      quicksort_long ((*S_2) + 1, params->s_2);
    }
  
  /* Print the sets in devverbose mode */
  if (test_verbose (OUTPUT_DEVVERBOSE))
    {
      outputf (OUTPUT_DEVVERBOSE, "S_1 = ");
      print_sets (OUTPUT_DEVVERBOSE, *S_1, *S_1_size);
      
      outputf (OUTPUT_DEVVERBOSE, "S_2 = {");
      for (i = 0UL; i + 1UL < params->s_2; i++)
	outputf (OUTPUT_DEVVERBOSE, "%ld, ", (*S_2)[i + 1UL]);
      if (i < params->s_2)
	outputf (OUTPUT_DEVVERBOSE, "%ld", (*S_2)[i + 1UL]); 
      outputf (OUTPUT_DEVVERBOSE, "}\n");
    }

  return 0;
}


int 
pm1fs2 (mpz_t f, const mpres_t X, mpmod_t modulus, 
	const faststage2_param_t *params)
{
  unsigned long nr;
  unsigned long i, j, l, lenF, tmplen;
  long *S_1; /* This is stored as a set of sets (arithmetic progressions of
		prime length */
  long *S_2; /* This is stored as a regular set */
  unsigned long S_1_size; /* The size in longs of S_1 */
  listz_t F;   /* Polynomial F has roots X^{k_1} for k_1 \in S_1, so has 
		  degree s_1. It is symmetric, so has only s_1 / 2 + 1 
		  distinct coefficients. The sequence h_j will be stored in 
		  the same memory and won't be a monic polynomial, so the 
		  leading 1 monomial of F will be stored explicitly. Hence we 
		  need s_1 / 2 + 1 entries. */
  listz_t tmp;
  mpzspm_t ntt_context;
  mpzspv_t g_ntt, h_ntt;
  mpz_t mt;   /* All-purpose temp mpz_t */
  mpres_t mr, product; /* All-purpose temp mpres_t */
  int youpi = ECM_NO_FACTOR_FOUND;
  long timetotalstart, timestart, timestop;

  timetotalstart = cputime ();

  ASSERT_ALWAYS (eulerphi (params->P) == params->s_1 * params->s_2);
  ASSERT_ALWAYS (params->s_1 < params->l);
  nr = params->l - params->s_1; /* Number of points we evaluate */

  outputf (OUTPUT_TRACE, "compare(a, b, n) = if (a != b,print(\"In PARI \", n, \" line, \", a \" != \" b)); /* PARI %ld */\n", pariline++);

  /* Precompute the small primes, primitive roots and inverses etc. for 
     the NTT. mpzspm_init() chooses the NTT primes large enough for 
     residues up to 4*l*modulus^2, so adding in Fourier space is ok. */
  ntt_context = mpzspm_init (params->l, modulus->orig_modulus);
  if (test_verbose (OUTPUT_DEVVERBOSE))
    {
      double modbits = 0.;
      outputf (OUTPUT_DEVVERBOSE, "CRT modulus = %lu", ntt_context->spm[0]->sp);
      modbits += log (ntt_context->spm[0]->sp);
      for (i = 1; i < ntt_context->sp_num; i++)
	{
	  outputf (OUTPUT_DEVVERBOSE, " * %lu", ntt_context->spm[i]->sp);
	  modbits += log (ntt_context->spm[i]->sp);
	}
      outputf (OUTPUT_DEVVERBOSE, ", has %f bits\n", modbits / log (2.));
    }


  /* Allocate memory for h_ntt. We allocate it before the memory for f 
     because that will get freed again, and we want to avoid heap 
     fragmentation. */
  h_ntt = mpzspv_init (params->l / 2 + 1, ntt_context);

  make_S_1_S_2 (&S_1, &S_1_size, &S_2, params);

  /* Allocate all the memory we'll need */
  /* Allocate the correct amount of space for each mpz_t or the 
     reallocations will up to double the time for stage 2! */
  mpz_init (mt);
  mpres_init (mr, modulus);
  mpres_init (product, modulus);
  lenF = params->s_1 / 2 + 1 + 1; /* Another +1 because poly_from_sets_V stores
				     the leading 1 monomial for each factor */
  tmplen = 4 * params->s_1 + 1000;
  F = init_list2 (lenF, abs (modulus->bits));
  tmp = init_list2 (tmplen, abs (modulus->bits));

  mpres_get_z (mt, X, modulus); /* mpz_t copy of X for printing */
  outputf (OUTPUT_TRACE, 
	   "N = %Zd; X = Mod(%Zd, N); /* PARI */\n", 
	   modulus->orig_modulus, mt);


  /* Compute the polynomial f(x) = \prod_{k_1 in S_1} (x - X^{2k_1}) */
  outputf (OUTPUT_VERBOSE, "Computing F from factored S_1");
  
  timestart = cputime ();
  
  /* First compute X^2 + 1/X^2 */
  mpres_init (mr, modulus);
  mpres_invert (mr, X, modulus);
  mpres_add (mr, mr, X, modulus);
  V (mr, mr, 2UL, modulus);
  
  i = poly_from_sets_V (F, mr, S_1, S_1_size, tmp, tmplen, modulus);
  ASSERT(2 * i == params->s_1);
  ASSERT(mpz_cmp_ui (F[i], 1UL) == 0);
  
  timestop = cputime ();
  outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);
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
  pm1_sequence_h (h_ntt, F, mr, params->s_1 / 2 + 1, modulus, ntt_context);

  clear_list (tmp, tmplen);
  clear_list (F, lenF);
  g_ntt = mpzspv_init (params->l, ntt_context);

  /* Make a symmetric copy of h in g_ntt. I.e. with h = [3, 2, 1], s_1 = 4, 
     l = 8, we want g_ntt = [3, 2, 1, 0, 0, 0, 1, 2] */
  mpzspv_set (g_ntt, 0, h_ntt, 0, params->s_1 / 2 + 1, ntt_context);
  mpzspv_revcopy (g_ntt, params->l - params->s_1 / 2, h_ntt, 1, 
		  params->s_1 / 2, ntt_context);
  /* Now we have [3, 2, 1, ?, ?, ?, 1, 2]. Fill the ?'s with zeros. */
  mpzspv_set_sp (g_ntt, params->s_1 / 2 + 1, 0, params->l - params->s_1 - 1,
		 ntt_context);
#if 0
  for (i = 0; i < params->l; i++)
      printf ("%lu ", g_ntt[0][i]);
  printf ("\n");
#endif
  /* Compute forward transform */
  outputf (OUTPUT_VERBOSE, "Computing forward NTT of h");
  timestart = cputime ();
  mpzspv_to_ntt (g_ntt, 0, params->l, params->l, 0, ntt_context);
  timestop = cputime ();
  outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);

#if 0
  for (i = 0; i < params->l; i++)
      printf ("%lu ", g_ntt[0][i]);
  printf ("\n");
#endif
  /* The forward transform is scrambled. We want elements [0 ... l/2]
     of the unscrabled data, that is all the coefficients with the most 
     significant bit in the index (in log2(l) word size) unset, plus the 
     element at index l/2. By scrambling, these map to the elements with 
     even index, plus the element at index 1. 
     The elements with scrambled index 2*i are stored in h[i], the
     element with scrambled index 1 is stored in h[params->l] */
  
#ifdef WANT_ASSERT
  /* Test that the coefficients are symmetic (if they were unscambled) and 
     that our algorithm for finding identical coefficients in the scrambled 
     data works */
  l = 5;
  for (i = 2; i < params->l; i += 2)
  {
      unsigned long t;
      /* This works, but why? */
      if (i + i / 2 > l)
	  l = l * 2 + 1;

      t = l - i;

      for (j = 0; j < ntt_context->sp_num; j++)
	  ASSERT (g_ntt[j][i] == g_ntt[j][t]);
  }
#endif
  /* Copy coefficients to h_ntt */
  for (j = 0; j < ntt_context->sp_num; j++)
  {
      for (i = 0; i < params->l / 2; i++)
	  h_ntt[j][i] = g_ntt[j][i * 2];
      h_ntt[j][params->l / 2] = g_ntt[j][1];
  }

  for (l = 0; l < params->s_2; l++)
    {
      pm1_sequence_g (g_ntt, X, params->P, params->l - 1 - params->s_1 / 2, 
		      params->l, params->m_1, S_2[l + 1], modulus, 
		      ntt_context);

      /* Do the convolution */
      outputf (OUTPUT_VERBOSE, "Computing forward NTT of g");
      timestart = cputime ();
      mpzspv_to_ntt (g_ntt, 0, params->l, params->l, 0, ntt_context);
      timestop = cputime ();
      outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);
      
      outputf (OUTPUT_VERBOSE, "Computing point-wise product");
      timestart = cputime ();
      for (j = 0; j < ntt_context->sp_num; j++)
      {
	  unsigned long m = 5UL;
	  const sp_t sp = ntt_context->spm[j]->sp; 
	  const sp_t mul_c = ntt_context->spm[j]->mul_c;
	  
	  g_ntt[j][0] = sp_mul (g_ntt[j][0], h_ntt[j][0], sp, mul_c);

	  for (i = 2; i < params->l; i += 2)
	  {
	      unsigned long t;

	      /* This works, but why? */
	      if (i + i / 2 > m)
		  m = m * 2 + 1;
	      
	      t = m - i;
	      ASSERT (t < params->l);
	      g_ntt[j][i] = sp_mul (g_ntt[j][i], h_ntt[j][i / 2], sp, mul_c);
	      g_ntt[j][t] = sp_mul (g_ntt[j][t], h_ntt[j][i / 2], sp, mul_c);
	  }
	  g_ntt[j][1] = sp_mul (g_ntt[j][1], h_ntt[j][params->l / 2], sp, 
				mul_c);
      }
      timestop = cputime ();
      outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);

      outputf (OUTPUT_VERBOSE, "Computing inverse NTT of g*h");
      timestart = cputime ();
      mpzspv_from_ntt (g_ntt, 0, params->l, 0, ntt_context);
      timestop = cputime ();
      outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);
     
      outputf (OUTPUT_VERBOSE, "Computing product of g*h coefficients");
      outputf (OUTPUT_TRACE, "\n");
      timestart = cputime ();
      mpres_set_ui (product, 1UL, modulus);
      for (i = 0; i < nr; i++)
      {
          /* FIXME: do in pieces of MPZSPV_NORMALISE_STRIDE length */
	  mpzspv_to_mpzv (g_ntt, params->s_1 / 2  + i, &mt, 1, ntt_context);
	  outputf (OUTPUT_TRACE, "r_%lu = %Zd; /* PARI */\n", i, mt);
	  if (mpz_sgn (mt) == 0)
	      outputf (OUTPUT_VERBOSE, "r_%lu = 0\n", i);
	  mpres_set_z_for_gcd (mr, mt, modulus);
	  mpres_mul (product, product, mr, modulus);
      }
      timestop = cputime ();
      outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);

      mpres_gcd (mt, product, modulus);
      if (mpz_cmp_ui (mt, 1UL) > 0)
	{
	  outputf (OUTPUT_NORMAL, "Found factor in stage 2\n");
	  mpz_set (f, mt);
	  youpi = ECM_FACTOR_FOUND_STEP2;
	  break;
	}
    }

  
  mpzspv_clear (g_ntt, ntt_context);
  mpzspv_clear (h_ntt, ntt_context);
  mpzspm_clear (ntt_context);
  mpz_clear (mt);
  mpres_clear (mr, modulus);
  mpres_clear (product, modulus);

  timestop = cputime ();
  outputf (OUTPUT_NORMAL, "Step 2 took %ldms\n", 
           timestop - timetotalstart);
  
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
	     mpmod_t modulus, ATTRIBUTE_UNUSED const unsigned long tmplen, 
	     mpz_t *t)
{
  ASSERT (tmplen >= 2);
  if (test_verbose (OUTPUT_TRACE))
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
  
  mpres_add (t[0], a_0, a_1, modulus);
  mpres_add (t[1], b_0, b_1, modulus);
  mpres_mul (t[1], t[0], t[1], modulus); /* t[1] = (a_0 + a_1) * (b_0 + b_1) = 
					    a_0*b_0 + a_0*b_1 + a_1*b_0 + 
					    a_1*b_1 */

  mpres_mul (r_0, a_0, b_0, modulus);    /* r_0 = a_0*b_0. We don't need a_0 
					    or b_0 any more now */
  mpres_sub (t[1], t[1], r_0, modulus);  /* t[1] = a_0*b_1 + a_1*b_0 + 
					    a_1*b_1 */
  
  mpres_mul (t[0], a_1, b_1, modulus);   /* t[0] = a_1*b_1. We don't need a_1 
					    or b_1 any more now */
  mpres_sub (r_1, t[1], t[0], modulus);  /* r_1 == a_0*b_1 + a_1*b_0 */
  
  mpres_mul (t[0], t[0], Delta, modulus); /* t[0] = a_1*b_1*Delta */
  mpres_add (r_0, r_0, t[0], modulus);   /* r_0 = a_0*b_0 + a_1*b_1*Delta */

  if (test_verbose (OUTPUT_TRACE))
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
  
  if (0 && pari)
    gmp_printf ("/* gfp_ext_sqr_norm1 */ (%Zd + %Zd * w)^2 %% N == ", a_0, a_1);
  
  mpres_mul (r_1, a_0, a_1, modulus);
  mpres_add (r_1, r_1, r_1, modulus);       /* r_1 = 2*a_0*a_1 */
  
  mpres_mul (r_0, a_0, a_0, modulus);
  mpres_add (r_0, r_0, r_0, modulus);
  mpres_sub_ui (r_0, r_0, 1UL, modulus);    /* r_0 = 2*a_0^2 - 1 */

  if (0 && pari)
    gmp_printf ("(%Zd + %Zd * w) %% N /* PARI C */\n", r_0, r_1);
}


/* Raise (a0 + a1*sqrt(Delta)) to the power e. (a0 + a1*sqrt(Delta)) is 
   assumed to have norm 1, i.e. a0^2 - a1^2*Delta == 1. The result is 
   (r0 * r1*sqrt(Delta)). a0, a1, r0 and r1 must not overlap */

static void 
gfp_ext_pow_norm1_ul (mpres_t r0, mpres_t r1, const mpres_t a0, 
                      const mpres_t a1, const long e, const mpres_t Delta, 
                      mpmod_t modulus, unsigned long tmplen, mpres_t *tmp)
{
  const unsigned long abs_e = labs (e);
  unsigned long mask = ~0UL - (~0UL >> 1);

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

  while ((abs_e & mask) == 0UL)
    mask >>= 1;

  mpres_set (r0, a0, modulus);
  mpres_set (r1, a1, modulus);

  while (mask > 1UL)
    {
      gfp_ext_sqr_norm1 (r0, r1, r0, r1, modulus);
      mask >>= 1;
      if (abs_e & mask)
	gfp_ext_mul (r0, r1, r0, r1, a0, a1, Delta, modulus, tmplen, tmp);
    }

  if (e < 0)
    mpres_neg (r1, r1, modulus);

  if (test_verbose (OUTPUT_TRACE))
    {
      mpz_t t;
      mpz_init (t);
      mpres_get_z (t, Delta, modulus);
      outputf (OUTPUT_TRACE, "/* gfp_ext_pow_norm1_ul */ w = quadgen (4*%Zd); "
               "N = %Zd; /* PARI */\n", t, modulus->orig_modulus);
      mpz_clear (t);
      outputf (OUTPUT_TRACE, "/* gfp_ext_pow_norm1_ul */ (");
      gfp_ext_print (a0, a1, modulus, OUTPUT_TRACE);
      outputf (OUTPUT_TRACE, ")^(%ld) == ", e);
      gfp_ext_print (r0, r1, modulus, OUTPUT_TRACE);
      outputf (OUTPUT_TRACE, " /* PARI C */\n");
    }
}


/* Same, but taking an mpz_t argument for the exponent */

static void 
gfp_ext_pow_norm1 (mpres_t r0, mpres_t r1, const mpres_t a0, 
                   const mpres_t a1, mpz_t e, const mpres_t Delta, 
                   mpmod_t modulus, unsigned long tmplen, mpres_t *tmp)
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
	     const long k, const unsigned long l, const mpres_t Delta, 
	     mpmod_t modulus, const unsigned long origtmplen, mpres_t *origtmp)
{
  mpres_t *r2_x = origtmp, *r2_y = origtmp + 2, *v = origtmp + 4, 
    *V2 = origtmp + 6;
  const unsigned long newtmplen = origtmplen - 7;
  mpres_t *newtmp = origtmp + 7;
  unsigned long i;

  if (l == 0UL)
    return;

  ASSERT (origtmplen >= 8UL);

  if (pari)
    gmp_printf ("/* In gfp_ext_rn2 */ ; a = %Zd + %Zd * w; /* PARI */\n",
		a_x, a_y, modulus->orig_modulus);

  /* Compute r[0] = a^(k^2). We do it by two exponentiations by k and use 
     v[0] and v[1] as temp storage */
  gfp_ext_pow_norm1_ul (v[0], v[1], a_x, a_y, k, Delta, modulus, newtmplen, 
		     newtmp);
  gfp_ext_pow_norm1_ul (r_x[0], r_y[0], v[0], v[1], k, Delta, modulus, 
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
    gmp_printf ("/* In gfp_ext_rn2 */ a^(%ld^2) %% N == (%Zd + %Zd * w) %% N "
		"/* PARI C */\n", k + 1, r_x[1], r_y[1]);
  
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
    gmp_printf ("/* In gfp_ext_rn2 */ a^(%ld^2+2) %% N == (%Zd + %Zd * w) %% N "
		"/* PARI C */\n", k + 1, r2_x[1], r2_y[1]);
  
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
      if (0 && pari)
	gmp_printf ("/* In gfp_ext_rn2 */ V(%lu, a + 1/a) %% N == %Zd %% N "
		    "/* PARI C */\n", 2 * (k + i) + 1, v[i % 2]);
    }
}


/* Compute g_i = x_0^{M-i} * r^{(M-i)^2} for 0 <= i < l. 
   x_0 = b_1^{2*k_2 + (2*m_1 + 1) * P}. r = b_1^P. */

static void
pp1_sequence_g (listz_t g_x, listz_t g_y, const mpres_t b1_x, 
		const mpres_t b1_y, const unsigned long P, 
		const mpres_t Delta, const long M, const unsigned long l, 
		const mpz_t m_1, const long k_2, mpmod_t modulus, 
		const unsigned long tmplen, mpz_t *tmp)
{
  mpres_t r_x, r_y, x0_x, x0_y, r1_x[2], r1_y[2], r2_x[2], r2_y[2], v2, v[2];
  unsigned long i;
  long timestart, timestop;

  outputf (OUTPUT_VERBOSE, "Computing g_i");
  timestart = cputime ();
  
  mpres_init (r_x, modulus);
  mpres_init (r_y, modulus);
  mpres_init (x0_x, modulus);
  mpres_init (x0_y, modulus);
  mpres_init (r1_x[0], modulus);
  mpres_init (r1_y[0], modulus);
  mpres_init (r1_x[1], modulus);
  mpres_init (r1_y[1], modulus);
  mpres_init (r2_x[0], modulus);
  mpres_init (r2_y[0], modulus);
  mpres_init (r2_x[1], modulus);
  mpres_init (r2_y[1], modulus);
  mpres_init (v2, modulus);
  mpres_init (v[0], modulus);
  mpres_init (v[1], modulus);

  ASSERT (tmplen >= 1);

  if (test_verbose (OUTPUT_TRACE))
    {
      mpz_t t;
      mpz_init (t);
      mpres_get_z (t, Delta, modulus);
      outputf (OUTPUT_TRACE, 
	       "\n/* pp1_sequence_g */ w = quadgen (4*%Zd); P = %lu; M = %ld; "
	       "k_2 = %ld; m_1 = %Zd; N = %Zd; /* PARI */\n", 
	       t, P, M, k_2, m_1, modulus->orig_modulus);
      
      outputf (OUTPUT_TRACE, "/* pp1_sequence_g */ b_1 = ");
      gfp_ext_print (b1_x, b1_y, modulus, OUTPUT_TRACE);
      outputf (OUTPUT_TRACE, "; /* PARI */\n", t);
      outputf (OUTPUT_TRACE, 
	       "/* pp1_sequence_g */ r = b_1^P; /* PARI */\n");
      outputf (OUTPUT_TRACE, "/* pp1_sequence_g */ "
	       "x_0 = b_1^(2*k_2 + (2*m_1 + 1) * P); /* PARI */\n");
      outputf (OUTPUT_TRACE, 
	       "/* pp1_sequence_g */ addrec(x) = x + 1/x; /* PARI */\n");
      mpz_clear (t);
    }

  /* Compute r */
  gfp_ext_pow_norm1_ul (r_x, r_y, b1_x, b1_y, P, Delta, modulus, tmplen, tmp);
  if (test_verbose (OUTPUT_TRACE))
    {
      outputf (OUTPUT_TRACE, "/* pp1_sequence_g */ r == ");
      gfp_ext_print (r_x, r_y, modulus, OUTPUT_TRACE);
      outputf (OUTPUT_TRACE, " /* PARI C */\n");
    }
  
  /* Compute x0 = x_0 */
  mpz_mul_2exp (tmp[0], m_1, 1UL);
  mpz_add_ui (tmp[0], tmp[0], 1UL);
  mpz_mul_ui (tmp[0], tmp[0], P);
  mpz_add_si (tmp[0], tmp[0], k_2);
  mpz_add_si (tmp[0], tmp[0], k_2); /* tmp[0] = 2*k_2 + (2*m_1 + 1) * P */
  outputf (OUTPUT_TRACE, "/* pp1_sequence_g */ 2*k_2 + (2*m_1 + 1) * P == "
           "%Zd /* PARI C */\n", tmp[0]);
  gfp_ext_pow_norm1 (x0_x, x0_y, b1_x, b1_y, tmp[0], Delta, modulus, 
                     tmplen - 1, tmp + 1);
  if (test_verbose (OUTPUT_TRACE))
    {
      outputf (OUTPUT_TRACE, "/* pp1_sequence_g */ x_0 == ");
      gfp_ext_print (x0_x, x0_y, modulus, OUTPUT_TRACE);
      outputf (OUTPUT_TRACE, " /* PARI C */\n");
    }

  /* Compute g[1] = r1[0] = x0^M * r^(M^2) = (x0 * r^M)^M.
     We use v[0,1] as temporary storage */
  gfp_ext_pow_norm1_ul (v[0], v[1], r_x, r_y, M, Delta, modulus, 
			tmplen, tmp); /* v[0,1] = r^M */
  gfp_ext_mul (v[0], v[1], v[0], v[1], x0_x, x0_y, Delta, modulus, 
	       tmplen, tmp); /* v[0,1] = r^M * x_0 */
  gfp_ext_pow_norm1_ul (r1_x[0], r1_y[0], v[0], v[1], M, Delta, modulus, 
			tmplen, tmp); /* r1[0] = (r^M * x_0)^M */
  mpres_get_z (g_x[0], r1_x[0], modulus);
  mpres_get_z (g_y[0], r1_y[0], modulus);
  outputf (OUTPUT_TRACE, 
	   "/* pp1_sequence_g */ g_0 = x_0^M * r^(M^2); /* PARI */\n");
  outputf (OUTPUT_TRACE, 
	   "/* pp1_sequence_g */ g_0 == %Zd + %Zd*w /* PARI C */\n", 
	   g_x[0], g_y[0]);

  /* Compute g[1] = r1[1] = x0^(M-1) * r^((M-1)^2) = (x0 * r^(M-1))^(M-1). 
     We use v[0,1] as temporary storage. FIXME: simplify, reusing g_0 */
  gfp_ext_pow_norm1_ul (v[0], v[1], r_x, r_y, M - 1, Delta, modulus, 
			tmplen, tmp);
  gfp_ext_mul (v[0], v[1], v[0], v[1], x0_x, x0_y, Delta, modulus, 
	       tmplen, tmp);
  gfp_ext_pow_norm1_ul (r1_x[1], r1_y[1], v[0], v[1], M - 1, Delta, 
			modulus, tmplen, tmp);
  mpres_get_z (g_x[1], r1_x[1], modulus);
  mpres_get_z (g_y[1], r1_y[1], modulus);
  outputf (OUTPUT_TRACE, 
	   "/* pp1_sequence_g */ g_1 = x_0^(M-1) * r^((M-1)^2); /* PARI */\n");
  outputf (OUTPUT_TRACE, 
	   "/* pp1_sequence_g */ g_1 == %Zd + %Zd*w /* PARI C */\n", 
	   g_x[1], g_y[1]);

  /* x0 := $x_0 * r^{2M - 3}$ */
  /* We don't need x0 after this so we overwrite it. We use v[0,1] as 
     temp storage for $r^{2M - 3}$. */
  gfp_ext_pow_norm1_ul (v[0], v[1], r_x, r_y, 2UL*M - 3UL, Delta, modulus,
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
      /* r1[i] = r2[i-1] * v[i-2] - r2[i-2], with indices of r2 and i taken
	 modulo 2. We store the new r1[i] in r for now */
      mpres_mul (r_x, r2_x[1 - i % 2], v[i % 2], modulus);
      mpres_sub (r_x, r_x,             r2_x[i % 2], modulus);
      mpres_mul (r_y, r2_y[1 - i % 2], v[i % 2], modulus);
      mpres_sub (r_y, r_y,             r2_y[i % 2], modulus);
      
      /* r2[i] = r2[i-1] * v[i-1] - r1[i-2] */
      mpres_mul (r2_x[i % 2], r2_x[1 - i % 2], v[1 - i % 2], modulus);
      mpres_sub (r2_x[i % 2], r2_x[i % 2],     r1_x[i % 2], modulus);
      mpres_mul (r2_y[i % 2], r2_y[1 - i % 2], v[1 - i % 2], modulus);
      mpres_sub (r2_y[i % 2], r2_y[i % 2],     r1_y[i % 2], modulus);
      mpres_set (r1_x[i % 2], r_x, modulus); /* FIXME, avoid this copy */
      mpres_set (r1_y[i % 2], r_y, modulus);
      mpres_get_z (g_x[i], r_x, modulus); /* FIXME, avoid these REDC */
      mpres_get_z (g_y[i], r_y, modulus); /* Keep r1, r2 in mpz_t ? */
      outputf (OUTPUT_TRACE, "/* pp1_sequence_g */ "
	       "x_0^(M-%lu) * r^((M-%lu)^2) == %Zd + %Zd*w /* PARI C */\n", 
	       i, i, g_x[i], g_y[i]);
     
      /* v[i] = v[i - 1] * V_2(a + 1/a) - v[i - 2] */
      mpres_mul (r_x, v[1 - i % 2], v2, modulus);
      mpres_sub (v[i % 2], r_x, v[i % 2], modulus);
      if (test_verbose (OUTPUT_TRACE))
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
  mpres_clear (r2_x[0], modulus);
  mpres_clear (r2_y[0], modulus);
  mpres_clear (r2_x[1], modulus);
  mpres_clear (r2_y[1], modulus);
  mpres_clear (v2, modulus);
  mpres_clear (v[0], modulus);
  mpres_clear (v[1], modulus);
  
  timestop = cputime ();
  outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);
}


/* Compute r[i] = b1^(-P*(k+i)^2) * f_i for i = 0, 1, ..., l-1, where "b1" is 
   an element of norm 1 in the quadratic extension ring */

static void
pp1_sequence_h (listz_t h_x, listz_t h_y, const listz_t f, const mpres_t b1_x, 
		const mpres_t b1_y, const long k, const unsigned long l, 
		const unsigned long P, const mpres_t Delta, mpmod_t modulus, 
		const unsigned long origtmplen, mpres_t *origtmp)
{
  mpres_t *s_x = origtmp, *s_y = origtmp + 3, *s2_x = origtmp + 6, 
    *s2_y = origtmp + 8, *v = origtmp + 10, *V2 = origtmp + 12;
  const unsigned long newtmplen = origtmplen - 13;
  mpres_t *newtmp = origtmp + 13;
  mpres_t rn_x, rn_y;
  unsigned long i;
  long timestart, timestop;

  if (l == 0UL)
    return;

  ASSERT (origtmplen >= 13);
  ASSERT (f != h_x);
  ASSERT (f != h_y);

  outputf (OUTPUT_VERBOSE, "Computing h_i");
  timestart = cputime ();

  mpres_init (rn_x, modulus);
  mpres_init (rn_y, modulus);

  if (test_verbose (OUTPUT_TRACE))
    {
      mpz_t t;
      mpz_init (t);
      mpres_get_z (t, Delta, modulus);
      outputf (OUTPUT_TRACE, "\n/* pp1_sequence_h */ w = quadgen (4*%Zd); "
	       "k = %ld; P = %lu; N = %Zd; /* PARI */\n", 
	       t, k, P, modulus->orig_modulus);
      outputf (OUTPUT_TRACE, "/* pp1_sequence_h */ b_1 = ");
      gfp_ext_print (b1_x, b1_y, modulus, OUTPUT_TRACE);
      outputf (OUTPUT_TRACE, "; rn = b_1^(-P); /* PARI */\n", t);
      mpz_clear (t);
      for (i = 0; i < l; i++)
	outputf (OUTPUT_TRACE, 
		 "/* pp1_sequence_h */ f_%lu = %Zd; /* PARI */\n", i, f[i]);
    }

  ASSERT (origtmplen >= 13UL);

  /* Compute r = b_1^P */
  gfp_ext_pow_norm1_ul (rn_x, rn_y, b1_x, b1_y, P, Delta, modulus, newtmplen,
		     newtmp);
  mpres_neg (rn_y, rn_y, modulus);

  /* Compute s[0] = r^(k^2). We do it by two exponentiations by k and use 
     v[0] and v[1] as temp storage */
  gfp_ext_pow_norm1_ul (v[0], v[1], rn_x, rn_y, k, Delta, modulus, newtmplen, 
		     newtmp);
  gfp_ext_pow_norm1_ul (s_x[0], s_y[0], v[0], v[1], k, Delta, modulus, 
		     newtmplen, newtmp);
  mpres_mul_z_to_z (h_x[0], s_x[0], f[0], modulus);
  mpres_mul_z_to_z (h_y[0], s_y[0], f[0], modulus);
  if (test_verbose (OUTPUT_TRACE))
    {
      outputf (OUTPUT_TRACE, "/* pp1_sequence_h */ rn^(k^2) == ");
      gfp_ext_print (s_x[0], s_y[0], modulus, OUTPUT_TRACE);
      outputf (OUTPUT_TRACE, " /* PARI C */\n");
    }

  /* Compute s[1] = r^((k+1)^2) = r^(k^2 + 2k + 1)*/
  if (l > 1)
    {
      /* v[0] + v[1]*sqrt(Delta) still contains r^k */
      gfp_ext_sqr_norm1 (s_x[1], s_y[1], v[0], v[1], modulus);
      /* Now s[1] = r^(2k) */
      gfp_ext_mul (s_x[1], s_y[1], s_x[1], s_y[1], s_x[0], s_y[0], Delta, 
		   modulus, newtmplen, newtmp);
      /* Now s[1] = r^(k^2 + 2k) */
      gfp_ext_mul (s_x[1], s_y[1], s_x[1], s_y[1], rn_x, rn_y, Delta, modulus, 
		   newtmplen, newtmp);
      /* Now s[1] = r^(k^2 + 2k + 1) = r^((k+1)^2) */
      if (test_verbose (OUTPUT_TRACE))
	{
	  outputf (OUTPUT_TRACE, "/* pp1_sequence_h */ rn^((k+1)^2) == ");
	  gfp_ext_print (s_x[1], s_y[1], modulus, OUTPUT_TRACE);
	  outputf (OUTPUT_TRACE, " /* PARI C */\n");
	}
      mpres_mul_z_to_z (h_x[1], s_x[1], f[1], modulus);
      mpres_mul_z_to_z (h_y[1], s_y[1], f[1], modulus);
    }
  
  /* Compute s2[0] = r^(k^2+2) = r^(k^2) * r^2 */
  gfp_ext_sqr_norm1 (v[0], v[1], rn_x, rn_y, modulus);
  gfp_ext_mul (s2_x[0], s2_y[0], s_x[0], s_y[0], v[0], v[1], Delta, modulus, 
	       newtmplen, newtmp);
  if (test_verbose (OUTPUT_TRACE))
    {
      outputf (OUTPUT_TRACE, "/* pp1_sequence_h */ rn^(k^2+2) == ");
      gfp_ext_print (s2_x[0], s2_y[0], modulus, OUTPUT_TRACE);
      outputf (OUTPUT_TRACE, " /* PARI C */\n");
    }

  /* Compute a^((k+1)^2+2) = a^((k+1)^2) * a^2 */
  gfp_ext_mul (s2_x[1], s2_y[1], s_x[1], s_y[1], v[0], v[1], Delta, modulus, 
	       newtmplen, newtmp);
  if (test_verbose (OUTPUT_TRACE))
    {
      outputf (OUTPUT_TRACE, "/* pp1_sequence_h */ rn^((k+1)^2+2) == ");
      gfp_ext_print (s2_x[1], s2_y[1], modulus, OUTPUT_TRACE);
      outputf (OUTPUT_TRACE, " /* PARI C */\n");
    }
  
  /* Compute V_2(r + 1/r). Since 1/r = rn_x - rn_y, we have r+1/r = 2*rn_x.
     V_2(x) = x^2 - 2, so we want 4*rn_x^2 - 2. */
  mpres_add (*V2, rn_x, rn_x, modulus); /* V2 = r + 1/r  = 2*rn_x */
  V (v[0], *V2, 2 * k + 1, modulus);  /* v[0] = V_{2k+1} (r + 1/r) */
  V (v[1], *V2, 2 * k + 3, modulus);  /* v[0] = V_{2k+3} (r + 1/r) */
  mpres_mul (*V2, *V2, *V2, modulus); /* V2 = 4*a_x^2 */
  mpres_sub_ui (*V2, *V2, 2UL, modulus); /* V2 = 4*a_x^2 - 2 */

  
  /* Compute the remaining r^((k+i)^2) values according to Peter's 
     recurrence */
  
  for (i = 2; i < l; i++)
    {
      /* r[i] = r2[i-1] * v[i-2] - r2[i-2], with indices of r2 and i taken
	 modulo 2 */
      mpres_mul (s_x[i % 3], s2_x[1 - i % 2], v[i % 2], modulus);
      mpres_mul (s_y[i % 3], s2_y[1 - i % 2], v[i % 2], modulus);
      mpres_sub (s_x[i % 3], s_x[i % 3], s2_x[i % 2], modulus);
      mpres_sub (s_y[i % 3], s_y[i % 3], s2_y[i % 2], modulus);
      
      /* r2[i] = r2[i-1] * v[i-1] - r[i-2] */
      mpres_mul (s2_x[i % 2], s2_x[1 - i % 2], v[1 - i % 2], modulus);
      mpres_mul (s2_y[i % 2], s2_y[1 - i % 2], v[1 - i % 2], modulus);
      mpres_sub (s2_x[i % 2], s2_x[i % 2], s_x[(i - 2) % 3], modulus);
      mpres_sub (s2_y[i % 2], s2_y[i % 2], s_y[(i - 2) % 3], modulus);
      
      /* v[i] = v[i - 1] * V_2(a + 1/a) - v[i - 2] */
      mpres_mul (newtmp[0], v[1 - i % 2], *V2, modulus);
      mpres_sub (v[i % 2], newtmp[0], v[i % 2], modulus);

      mpres_mul_z_to_z (h_x[i], s_x[i % 3], f[i], modulus);
      mpres_mul_z_to_z (h_y[i], s_y[i % 3], f[i], modulus);
    }

  if (test_verbose (OUTPUT_TRACE))
    {
      for (i = 0; i < l; i++)
	gmp_printf ("/* pp1_sequence_h */ (rn^((k+%lu)^2) * f_%lu) == (%Zd"
		    " + %Zd * w) /* PARI C */\n", i, i, h_x[i], h_y[i]);
    }

  timestop = cputime ();
  outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);
}


int 
pp1fs2 (mpz_t f, const mpres_t X, mpmod_t modulus, 
	const faststage2_param_t *params)
{
  unsigned long nr;
  unsigned long i, l, lenF, lenH, lenG, lenR, tmplen;
  long *S_1; /* This is stored as a set of sets (arithmetic progressions of
		prime length */
  long *S_2; /* This is stored as a regular set */
  unsigned long S_1_size; /* The size in longs of S_1 */
  listz_t F;   /* Polynomial F has roots X^{k_1} for k_1 \in S_1, so has 
		  degree s_1. It is symmetric, so has only s_1 / 2 + 1 
		  distinct coefficients. The sequence h_j will be stored in 
		  the same memory and won't be a monic polynomial, so the 
		  leading 1 monomial of F will be stored explicitly. Hence we 
		  need s_1 / 2 + 1 entries. */
  listz_t g_x, g_y, fh_x, fh_y, h_x, h_y, tmp, R_x, R_y; 
  /* g, h and R are in extension ring, g = g_1 + g_2 * sqrt(Delta) (likewise 
     for h, R) */
  mpres_t b1_x, b1_y, Delta;
  mpz_t mt;   /* All-purpose temp mpz_t */
  mpres_t mr; /* All-purpose temp mpres_t */
  int youpi = ECM_NO_FACTOR_FOUND;
  long timetotalstart, timestart, timestop;

  timetotalstart = cputime ();

  ASSERT_ALWAYS (eulerphi (params->P) == params->s_1 * params->s_2);
  ASSERT_ALWAYS (params->s_1 < params->l);
  nr = params->l - params->s_1; /* Number of points we evaluate */

  outputf (OUTPUT_TRACE, "compare(a, b, n) = if (a != b,print(\"In PARI \", "
	   "n, \" line, \", a \" != \" b)); /* PARI %ld */\n", pariline++);

  make_S_1_S_2 (&S_1, &S_1_size, &S_2, params);

  /* Allocate all the memory we'll need */
  /* Allocate the correct amount of space for each mpz_t or the 
     reallocations will up to double the time for stage 2! */
  mpz_init (mt);
  mpres_init (mr, modulus);
  mpres_init (b1_x, modulus);
  mpres_init (b1_y, modulus);
  mpres_init (Delta, modulus);
  lenF = params->s_1 / 2 + 1 + 1; /* Another +1 because poly_from_sets_V stores
				     the leading 1 monomial for each factor */
  lenH = params->s_1 + 1;
  lenG = params->l;
  lenR = nr;
  F = init_list2 (lenF, abs (modulus->bits));
  fh_x = init_list2 (lenF, abs (modulus->bits));
  fh_y = init_list2 (lenF, abs (modulus->bits));
  h_x = malloc (lenH * sizeof (mpz_t));
  h_y = malloc (lenH * sizeof (mpz_t));
  g_x = init_list2 (lenG, abs (modulus->bits));
  g_y = init_list2 (lenG, abs (modulus->bits));
  R_x = init_list2 (lenR, abs (modulus->bits));    
  R_y = init_list2 (lenR, abs (modulus->bits));    
  tmplen = 3UL * params->l + list_mul_mem (params->l / 2) + 20;
  outputf (OUTPUT_DEVVERBOSE, "tmplen = %lu\n", tmplen);
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
  tmp = init_list2 (tmplen, abs (modulus->bits));
#endif

  mpres_get_z (mt, X, modulus); /* mpz_t copy of X for printing */
  outputf (OUTPUT_TRACE, 
	   "N = %Zd; X = Mod(%Zd, N); /* PARI */\n", 
	   modulus->orig_modulus, mt);

  /* Compute the polynomial f(x) = \prod_{k_1 in S_1} (x - X^{2 k_1}) */
  outputf (OUTPUT_VERBOSE, "Computing F from factored S_1");
  
  timestart = cputime ();
  V (mr, X, 2UL, modulus);
  i = poly_from_sets_V (F, mr, S_1, S_1_size, tmp, tmplen, modulus);
  ASSERT(2 * i == params->s_1);
  ASSERT(mpz_cmp_ui (F[i], 1UL) == 0);
  timestop = cputime ();
  outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);
  if (test_verbose (OUTPUT_TRACE))
    {
      for (i = 0; i < params->s_1 / 2 + 1; i++)
	outputf (OUTPUT_TRACE, "f_%lu = %Zd; /* PARI */\n", i, F[i]);
      outputf (OUTPUT_TRACE, "f(x) = f_0");
      for (i = 1; i < params->s_1 / 2 + 1; i++)
	outputf (OUTPUT_TRACE, "+ f_%lu * (x^%lu + x^(-%lu))", i, i, i);
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
      mpz_t t;
      mpz_init (t);
      mpres_get_z (t, Delta, modulus);
      outputf (OUTPUT_TRACE, 
	       "Delta = Mod(%Zd, N); w = quadgen (4*lift(Delta)); b_1 = ", t);
      gfp_ext_print (b1_x, b1_y, modulus, OUTPUT_TRACE);
      outputf (OUTPUT_TRACE, "; /* PARI */\n");
      outputf (OUTPUT_TRACE, "X == b_1 + 1/b_1 /* PARI C */\n");
      mpz_clear (t);
    }

  /* Compute the h sequence h_j = b1^(P*-j^2) * f_j for 0 <= j <= s_1 */
  pp1_sequence_h (fh_x, fh_y, F, b1_x, b1_y, 0L, params->s_1 / 2 + 1, 
		  params->P, Delta, modulus, tmplen, tmp);
  /* We don't need F(x) any more */
  clear_list (F, lenF);

  /* Make a symmetric copy of fh in h. */
  for (i = 0; i < params->s_1 / 2 + 1; i++)
    {
      *(h_x[i]) = *(fh_x[params->s_1 / 2 - i]); /* Clone the mpz_t. 
						   Don't tell Torbjörn */
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
      const long M = params->l - 1 - params->s_1 / 2;
      pp1_sequence_g (g_x, g_y, b1_x, b1_y, params->P, Delta, M, params->l, 
		      params->m_1, S_2[l + 1], modulus, tmplen, tmp);
      
      /* Do the two convolution products */
      outputf (OUTPUT_VERBOSE, "TMulGen of g_x and h_x");
      timestart = cputime ();
      TMulGen (R_x, nr - 1, h_x, params->s_1, g_x, params->l - 1, tmp,
	       modulus->orig_modulus);
      timestop = cputime ();
      outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);
      outputf (OUTPUT_VERBOSE, "TMulGen of g_y and h_y");
      timestart = cputime ();
      TMulGen (R_y, nr - 1, h_y, params->s_1, g_y, params->l - 1, tmp,
	       modulus->orig_modulus);
      timestop = cputime ();
      outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);
      for (i = 0; i < nr; i++)
	mpres_mul_z_to_z (R_y[i], Delta, R_y[i], modulus);
      list_add (R_x, R_x, R_y, nr);

      timestart = cputime ();
      {
        mpres_t tmpres, tmpprod;
        mpres_init (tmpres, modulus);
        mpres_init (tmpprod, modulus);
#define TEST_ZERO_RESULT
#ifdef TEST_ZERO_RESULT
        mpz_mod (R_x[0], R_x[0], modulus->orig_modulus);
        if (mpz_sgn (R_x[0]) == 0)
          outputf (OUTPUT_VERBOSE, "R[0] == 0\n");
#endif
        mpres_set_z_for_gcd (tmpprod, R_x[0], modulus);
        for (i = 1; i < nr; i++)
          {
#ifdef TEST_ZERO_RESULT
            mpz_mod (R_x[i], R_x[i], modulus->orig_modulus);
            if (mpz_sgn (R_x[i]) == 0)
              outputf (OUTPUT_VERBOSE, "R[%lu] == 0\n", i);
#endif
            mpres_set_z_for_gcd (tmpres, R_x[i], modulus);
            mpres_mul (tmpprod, tmpprod, tmpres, modulus); 
          }
        if (test_verbose(OUTPUT_RESVERBOSE))
          {
            mpres_get_z (tmp[0], tmpprod, modulus);
            outputf (OUTPUT_RESVERBOSE, "Product of R[i] = %Zd (times some "
                     "power of 2 if REDC was used! Try -mpzmod)\n", tmp[0]);
          }
        mpres_gcd (tmp[0], tmpprod, modulus);
        mpres_clear (tmpprod, modulus);
        mpres_clear (tmpres, modulus);
      }
      timestop = cputime ();
      outputf (OUTPUT_VERBOSE, "Computing product of F(g_i)^(1) took %lu ms\n", 
	       timestop - timestart);
      
      if (mpz_cmp_ui (tmp[0], 1UL) > 0)
	{
	  outputf (OUTPUT_NORMAL, "Found factor in stage 2\n");
	  mpz_set (f, tmp[0]);
	  youpi = ECM_FACTOR_FOUND_STEP2;
	  break;
	}
    }

  mpz_clear (mt);
  mpres_clear (mr, modulus);
  mpres_clear (b1_x, modulus);
  mpres_clear (b1_y, modulus);
  mpres_clear (Delta, modulus);
  clear_list (fh_x, lenF);
  clear_list (fh_y, lenF);
  free (h_x);
  free (h_y);
  clear_list (g_x, lenG);
  clear_list (g_y, lenG);
  clear_list (R_x, lenR);
  clear_list (R_y, lenR);
  clear_list (tmp, tmplen);

  timestop = cputime ();
  outputf (OUTPUT_NORMAL, "Step 2 took %ld ms\n", 
           timestop - timetotalstart);

  return youpi;
}


#ifdef TESTDRIVE


int main (int argc, char **argv)
{
  unsigned long pn, d, i, j, tmplen, setsize, lmax = 1024;
  listz_t F, tmp;
  mpz_t r, N, B2min, B2;
  long *L;
  int selftest = 0;
  mpmod_t modulus;
  faststage2_param_t params;

  ECM_STDOUT = stdout;
  ECM_STDERR = stderr;
  set_verbose (OUTPUT_DEVVERBOSE);

  mpz_init (N);
  mpz_init (B2min);
  mpz_init (B2);

  pn = 1;
  while (argc > 1)
    {
      if (strcmp (argv[1], "-v") == 0)
	{
	  verbose++;
	  inc_verbose ();
	}
      else if (strcmp (argv[1], "-p") == 0)
	pari = 1;
      else if (strcmp (argv[1], "-t") == 0)
	selftest = 1;
      else if (argc > 2 && strcmp (argv[1], "-N") == 0)
	{
	  mpz_set_str (N, argv[2], 0);
	  argc--;
	  argv++;
	}
      else if (argc > 2 && strcmp (argv[1], "-l") == 0)
        {
	  lmax = strtoul (argv[2], NULL, 10);
	  argc--;
	  argv++;
        }
      else
	break;
      argc--;
      argv++;
    }
  if (argc > 1)
    pn = strtoul (argv[1], NULL, 10);
  if (argc > 2)
    mpz_set_str (B2, argv[2], 0);
  if (argc > 3)
  {
      mpz_set (B2min, B2);
      mpz_set_str (B2, argv[3], 0);
  }
  gmp_printf ("B2min = %Zd, B2 = %Zd\n", B2min, B2);

  if (mpz_cmp (B2, B2min) > 0)
    {
      mpz_init (params.m_1);
      pn = choose_P (B2min, B2, lmax, &params, NULL, NULL);
    }

  d = eulerphi (pn);

  L = get_factored_sorted_sets (&setsize, pn);

  F = init_list (d);
  tmplen = 10*d+10;
  if (tmplen < 2000)
    tmplen = 2000;
  tmp = init_list (tmplen);
  mpz_init (r);
  mpz_set_ui (r, 3UL);
  if (mpz_sgn (N) == 0)
    {
      /* By default, use the Mersenne prime 2^31-1 as the modulus */
      mpz_set_ui (N, 1UL);
      mpz_mul_2exp (N, N, 31UL);
      mpz_sub_ui (N, N, 1UL);
    }
  mpmod_init (modulus, N, 0);
  /* We don't need N anymore now */
  mpz_clear (N);
  if (pari)
    gmp_printf ("N = %Zd; r = Mod(%Zd, N); /* PARI */\n", 
		modulus->orig_modulus, r);

  /************************************************************
      Simple check of list_mul_symmetric() 
  ************************************************************/

  for (i = 0; i < 4; i++)
    mpz_set_ui (tmp[i], i + 1);
  /* tmp[0 .. 3] = [1, 2, 3, 4] = 4*x^3 + 3*x^2 + 2*x + 1 */

  /* Compute (4*(x^3+x^{-3}) + 3*(x^2+x^{-2}) + 2*(x+x^{-1}) + 1)^2 =
     16*(x^6+x^{-6}) + 24*(x^5+x^{-5}) + 25*(x^4+x^{-4}) + 20*(x^3+x^{-3}) + 
     26*(x^2+x^{-2}) + 40*(x+x^{-1}) + 59, so we expect in 
     tmp[0 .. 6] = [59, 40, 26, 20, 25, 24, 16] */
  list_mul_symmetric (tmp, tmp, 4, tmp, 4, tmp + 7, tmplen - 7);

  if (mpz_cmp_ui (tmp[6], 16UL) != 0 || mpz_cmp_ui (tmp[5], 24UL) != 0 ||
      mpz_cmp_ui (tmp[4], 25UL) != 0 || mpz_cmp_ui (tmp[3], 20UL) != 0 ||
      mpz_cmp_ui (tmp[2], 26UL) != 0 || mpz_cmp_ui (tmp[1], 40UL) != 0 ||
      mpz_cmp_ui (tmp[0], 59UL) != 0)
    {
      list_output_poly (tmp, 7, 0, 0, "Error, list_mul_symmetric produced ", 
			OUTPUT_ERROR);
      abort ();
    }

  for (i = 0; i < 4; i++)
    mpz_set_ui (tmp[i], i + 1);
  /* tmp[0 .. 3] = [1, 2, 3, 4] = 4*x^3 + 3*x^2 + 2*x + 1 */

  /* Compute (4*(x^3+x^{-3}) + 3*(x^2+x^{-2}) + 2*(x+x^{-1}) + 1) *
     (3*(x^2+x^-2) + 2*(x+x^-1) + 1), so we expect in 
     tmp[0 .. 5] = [27, 28, 18, 16, 17, 12] */
  list_mul_symmetric (tmp, tmp, 4, tmp, 3, tmp + 6, tmplen - 6);

  if (mpz_cmp_ui (tmp[0], 27UL) != 0 || mpz_cmp_ui (tmp[1], 28UL) != 0 ||
      mpz_cmp_ui (tmp[2], 18UL) != 0 || mpz_cmp_ui (tmp[3], 16UL) != 0 ||
      mpz_cmp_ui (tmp[4], 17UL) != 0 || mpz_cmp_ui (tmp[5], 12UL) != 0)
    {
      list_output_poly (tmp, 6, 0, 0, "Error, list_mul_symmetric produced ", 
			OUTPUT_ERROR);
      abort ();
    }

  /* Simple test of list_scale_V(). Set F(x) = \sum_{i=0}^{7} (i+1) x^i
     Compute F(2 x) F(1/2 x) */

  {
    mpres_t Q;
    unsigned long len;

    mpres_init (Q, modulus);
    mpres_set_ui (Q, 2UL, modulus); /* This corresponds to Q = 1 + 1/1, 
				       or gamma = 1 */
    for (len = 1; len <= 100; len++)
      {
	for (i = 0; i < len; i++)
	  mpz_set_ui (tmp[i], i + 1);
	list_output_poly (tmp, len, 0, 1, "Input to list_scale_V: ", 
			  OUTPUT_TRACE);
	
	list_scale_V (tmp + len, tmp, Q, len - 1, modulus, tmp + 3*len, 
		      tmplen - 3*len);
	
	list_mod (tmp + len, tmp + len, 2 * len - 1, modulus->orig_modulus);
	list_output_poly (tmp + len, 2*len-1, 0, 1, 
			  "Output of list_scale_V: ", OUTPUT_TRACE);
	
	/* With Q = 2 = 1 + 1/1, gamma = 1 and F(gamma*X)*F(1/gamma *X) =F(X)^2
	   Compare with a simple symmetic multiply */
	list_mul_symmetric (tmp + 3 * len, tmp, len, tmp, len, tmp + 5*len,
			    tmplen - 5*len);
	
	list_mod (tmp + 3*len, tmp + 3*len, 2*len - 1, modulus->orig_modulus);
	
	for (i = 0; i <= 2 * len - 2; i++)
	  ASSERT(mpz_cmp (tmp[len + i], tmp[3*len + i]) == 0);
      }
    mpres_clear (Q, modulus);
  }

  /* Build the polynomial */
  
  poly_from_sets (F, r, L, setsize, tmp, tmplen, modulus->orig_modulus);

  if (pn % 4 != 2)
    {
      /* The leading monomial of F in implicit */
      if (mpz_cmp_ui (F[0], 1UL) != 0)
	printf ("Error, F[0] != 1, F is not symmetric!\n");
      for (i = 1; i < d / 2; i++)
	{
	  if (mpz_cmp (F[i], F[d - i]) != 0)
	    {
	      printf ("Error, F[%lu] != F[%lu - %lu], F is not symmetric!\n",
		      i, d, i);
	    }
	}
    }
  
  if (pari)
    {
      printf ("F(x) = x^%lu", d);
      for (i = d - 1; i > 0; i--)
	if (mpz_sgn (F[i]) != 0)
	  gmp_printf (" + %Zd * x^%lu", F[i], i);
      gmp_printf(" + %Zd /* PARI */\n", F[0]);
    }

  {
    mpres_t Q, mr;
    mpres_init (Q, modulus);
    mpres_init (mr, modulus);
    mpres_set_z (mr, r, modulus);
    mpres_invert (Q, mr, modulus);
    mpres_add (Q, Q, mr, modulus);
    poly_from_sets_V (tmp, Q, L, setsize, tmp + d, tmplen - d, modulus);
    mpres_clear (Q, modulus);
    mpres_clear (mr, modulus);
  }
  list_mod (tmp, tmp, d/2, modulus->orig_modulus);
  /* Check that the polynomials produced by poly_from_sets() and by 
     poly_from_sets_V() are identical */

  for (i = 0; i < d / 2; i++)
    {
      ASSERT(mpz_cmp (tmp[i], F[i + d / 2]) == 0);
      if (mpz_cmp (tmp[i], F[i + d / 2]) != 0)
	break; /* In case we don't have ASSERT on */
    }

  if (i == d / 2)
    outputf (OUTPUT_DEVVERBOSE, "Polynomials produced by poly_from_sets() "
	     "and poly_from_sets_V() agree.\n");
  
  /* Test gfp_ext_pow () */


  /* Test gfp_ext_rn2() */

  {
    mpres_t a_x, a_y, Delta;
    const long k = 3;

    mpres_init (a_x, modulus);
    mpres_init (a_y, modulus);
    mpres_init (Delta, modulus);

    mpres_set_ui (a_x, 9UL, modulus);
    mpres_set_ui (a_y, 4UL, modulus);
    mpres_set_ui (Delta, 5UL, modulus); /* norm = 9^2 - 4^2*5 = 1 */
    printf ("w = quadgen(20); /* PARI */\n");

    gfp_ext_rn2 (tmp, tmp+d, a_x, a_y, k, d, Delta, modulus,
		 tmplen - 2*d, tmp + 2*d);

    for (i = 0; i < d; i++)
      {
	gmp_printf ("(%Zd + %Zd*w)^%d %% N == (%Zd + %Zd*w) %% N "
		    "/* PARI from gfp_ext_rn2 */\n",
		    a_x, a_y, (k+i)*(k+i), tmp[i], tmp[d+i]);
	gfp_ext_pow_norm1_ul (tmp[i], tmp[d+i], a_x, a_y, 
			   (k+(long)i)*(k+(long)i), Delta, modulus, 
			   tmplen - 2*d, tmp + 2*d);
	gmp_printf ("(%Zd + %Zd*w)^%d %% N == (%Zd + %Zd*w) %% N "
		    "/* PARI from gfp_ext_pow_norm1_ul */\n",
		    a_x, a_y, (k+i)*(k+i), tmp[i], tmp[d+i]);
      }

    mpz_clear (a_x);
    mpz_clear (a_y);
    mpz_clear (Delta);
  }


  if (selftest) /* Do some self-tests */
    {
      long *sumset = malloc (d *sizeof (long));
      unsigned long t;
      
      t = sum_sets (sumset, L, setsize, 0L);
      ASSERT (t == d);
      
      if (pari)
	{
	  printf ("exponents = [");
	  for (i = 0; i < d - 1; i++)
	    printf ("%ld, ", sumset[i]);
	  printf ("%ld]; /* PARI */\n", sumset[d - 1]);
	  printf ("lift(prod(k=1,length(exponents),(x-r^exponents[k]))) "
		  "== F(x) /* PARI C */\n");
	}
      
      mpz_invert (tmp[0], r, modulus->orig_modulus);
      
      /* Check that all the elements in the sumset are really exponents of the
	 roots of F */
      gmp_printf ("Selftest: checking that %Zd^((Z/%luZ)*) are roots of F(x)\n", 
	      r, pn);
      for (i = 0; i < d; i++)
	{
	  if (sumset[i] < 0)
	    mpz_powm_ui (tmp[1], tmp[0], (unsigned long) (-sumset[i]), 
			 modulus->orig_modulus);
	  else
	    mpz_powm_ui (tmp[1], r, (unsigned long) sumset[i], 
			 modulus->orig_modulus);
	  list_eval_poly (tmp[2], F, tmp[1], d, 1, modulus->orig_modulus, 
			  tmp + 3);
	  if (mpz_sgn (tmp[2]) != 0)
	    printf ("Error, r^%ld is not a root of F\n", sumset[i]);
	}
      
      /* Check that the set of sums is really a complete set of representatives
	 of residue classes coprime to pn */

      printf ("Selftest: checking that set of sums is congruent to (Z/%luZ)*\n",
	      pn);

      for (i = 0; i < d; i++)
	if (sumset[i] >= 0)
	  sumset[i] = sumset[i] % pn;
	else
	  sumset[i] = pn - (-sumset[i]) % pn;
      quicksort_long (sumset, d);
      for (i = 1, j = 0; i < pn; i++)
	if (lgcd (i, pn) == 1)
	  {
	    ASSERT((unsigned long) sumset[j] == i);
	    j++;
	  }
      
      free (sumset);

      printf ("Selftest finished\n");
    }

  mpz_clear (r);
  
  return 0;
}
#endif
