#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include "ecm-impl.h"

#ifdef TESTDRIVE
#include <string.h>
static int verbose = 0;
#endif

/* Some useful PARI functions:
   sumset(a,b) = {local(i, j, l); l = listcreate (length(a) * length(b)); for (i = 1, length(a), for (j = 1, length(b), listput(l, a[i] + b[j]))); listsort (l, 1); l}

*/

static int pariline = 0;

static int 
sum_sets (long *sum, long *sets, unsigned long setsize, long add)
{
  unsigned long i, j = 0, l;

  if (setsize == 0)
    return 0;

  ASSERT (sets[0] >= 0);
  l = (unsigned long) sets[0];
  ASSERT (l > 1);
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

static int 
sum_sets_minmax (long *sets, unsigned long setsize, int minmax)
{
  unsigned long i, l, sum = 0;
  long extremum;

  ASSERT (minmax == 1 || minmax == -1);

  while (setsize > 0)
    {
      ASSERT (sets[0] >= 0);
      l = (unsigned long) sets[0];
      ASSERT (l > 1);
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

  swapsort_long (a, a+l-1, &t);
  if (l == 2)
    return;

  i = (l-1) / 2;
  swapsort_long (a, a+i, &t);
  swapsort_long (a+i, a+l-1, &t);
  if (l == 3)
    return;

  pivot = a[i]; /* Median of three */

  /* Stuff <= pivot goes in first list */
  swap_long (a+1, a+i, &t); /* a[i] == pivot, so swap to start */

  /* Invariant: a[0 ... i-1] <= pivot, a[j+1 ... l-1] > pivot */
  for (i = 2, j = l - 1; i < j;)
    if (a[i] > pivot)
      {
	for (; a[j] > pivot; j--);
	if (i < j)
	  swap_long (a+(i++), a+j, &t);
      }
    else
      i++;

  ASSERT(a[i-1] <= pivot);
  ASSERT(a[i] > pivot);

  swap_long (a+1, a+i-1, &t); /* Now a[i-1] == pivot */

  quicksort_long (a, i - 1);
  quicksort_long (a + i, l - i);

  for (i = 1; i < l; i++)
    ASSERT (a[i - 1] <= a[i]);
}


/* Return a set L of sets M_i so that M_1 + ... + M_k is congruent to 
   S_n, which is the set of residue classes coprime to n.
   L is stored as #M_1+1, M_1, #M_2+1, M_2, ..., #M_k+1, M_k. I.e. for n=15,
   S_n = 5*S_3 + 3*S_5 = 5*{-1,1} + 3*{-3,-1,1,3} = 
         5*{-1,1} + 3*{-2, 2} + 3*{-1,1}
   L = [2, -5, 5, 2, -6, 6, 2, -3, 3]
   Return the space (in longs) needed in L. 
   If L is the NULL pointer, nothing will be stored in L. The correct
   return value (amount of space needed in L) will still be computed.
   Currently assumes that n is squarefree. */
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
      ASSERT ((unsigned long) p <= n);
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


static void 
swap_sets (long *T, const unsigned long l, const unsigned long m)
{
  const unsigned long d = gcd (l, m);
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


/* Sort the sets in F into order of ascending cardinality */

static void 
sort_sets (long *F, unsigned long size)
{
  unsigned long a, l, m;
  int more = 1;

  /* Simple bubble sort */
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
  unsigned int i, j;
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


static void
list_output_poly (listz_t l, unsigned long len, int monic, char *prefix,
		  int verbosity)
{
  unsigned long i;

  if (prefix != NULL)
    outputf (verbosity, prefix);
  if (monic)
    outputf (verbosity, "x^%lu + ", len);
  for (i = len - 1; i > 0; i--)
    outputf (verbosity, "%Zd * x^%lu + ", l[i], i);
  if (len > 0)
    outputf (verbosity, "%Zd\n", l[0]);
  else
    outputf (verbosity, "\n");
}

/* Multiply P[deg-i] by r^(ki), for 1 <= i <= deg. Needs 3 entries in tmp. */
/* Returns 0 if a modular inversion failed (in which case R is left 
   unchanged), 1 otherwise */

static int
list_scale_rev (listz_t R, listz_t S, mpz_t r, long k, unsigned long deg, 
	    mpz_t modulus, listz_t tmp, unsigned long tmplen)
{
  unsigned long i;

  ASSERT (tmplen >= 3);
  mpz_powm_ui (tmp[0], r, labs (k), modulus);
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
      mpz_mul (tmp[2], tmp[1], tmp[0]);
      mpz_mod (tmp[1], tmp[2], modulus);
    }
  if (i <= deg)
    {
      mpz_mul (tmp[2], S[deg-i], tmp[1]);
      mpz_mod (R[deg-i], tmp[2], modulus);
    }

  return 1;
}

/* Check if l is an (anti-) symmetric, possibly monic, polynomial. 
   Returns -1 if it is (anti-)symmetric, or the smallest index i where 
   l[i] != l[len - 1 + monic - i])
   If anti == 1, the list is checked for symmetry, if it is -1, for
   antisymmetry.
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

/* Evaluate a polynomial of degree n-1 with all coefficients given in F[],
   or of degree n with an implicit leading 1 monomial not stored in F[],
   at x modulo modulus. Result goes in r. tmp needs 2 entries. */

static void 
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
		mpz_t *tmp, unsigned long tmplen, mpz_t modulus)
{

  unsigned long l, c; /* Cardinality of this set */
  unsigned long i, j, deg;
  long k;
  
  if (setsize == 0)
    {
      mpz_set_si (F[0], -1);
      return 1;
    }

  ASSERT (sets[0] >= 0);
  l = (unsigned long) sets[0];
  ASSERT (l > 1); /* Empty sets indicate a bug somewhere */
  ASSERT (l <= setsize);
  c = l - 1;
  deg = poly_from_sets (F, r, sets + l, setsize - l, tmp, tmplen, modulus);

  if (test_verbose (OUTPUT_TRACE))
    list_output_poly (F, deg, 1, "F = ", OUTPUT_TRACE);

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
	    list_output_poly (F + j, deg, 1, NULL, OUTPUT_TRACE);
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
      /* Multiply a monic polynomial A of degree i*deg at F
	 with a monic polynomial B of degree deg at F + i*deg */
      const listz_t B = F + i*deg;
      
      /* First piece of A */
      list_mul (tmp, F, deg, (i == 1), B, deg, 1, tmp + 2 * deg);
      list_mod (F, tmp, deg, modulus);
      list_swap (tmp, tmp + deg, deg);
      
      for (j = 1; j < i; j++) /* Process the remaining i-1 pieces of A */
	{
	  list_mul (tmp + deg, 
		    F + j * deg, deg, (j + 1 == i),
		    B, deg, 1, /* B */
		    tmp + 3 * deg);
	  /* Add low part of this product to previous product's high part */
	  list_add (tmp, tmp, tmp + deg, deg);
	  list_mod (F + j * deg, tmp, deg, modulus);
	  list_swap (tmp, tmp + 2 * deg, deg); /* Move this product's high 
						  part to beginning of tmp */
	}
      list_mod (F + j * deg, tmp, deg, modulus);
      if (test_verbose (OUTPUT_TRACE))
	list_output_poly (F, (i + 1) * deg, 1, NULL, OUTPUT_TRACE);
    }

#if 0 && defined(WANT_ASSERT)
  /* Test that the polynomial is symmetric if the degree is even, or anti-
     symmetric if it is odd */
  if (c != 1)
    ASSERT (list_is_symmetric (F, c * deg, 1, (c * deg) % 2 == 0 ? 1 : -1, 
			       modulus) == -1);
#endif
  
  return c * deg;
}


/* Build polynomial F(x) = \prod_{0<k<n, (k,n)=1} (x - r^k). Returns
   the number of multiplications used */

/* For Pari: F(n,r)=prod(k=1,n-1,if(gcd(k,n)==1,(x-r^k),1)) */

static int
poly_from_roots_coprime (mpz_t *F, mpz_t r, const unsigned long pn, 
			 mpz_t *tmp, unsigned long tmplen, mpz_t modulus)
{
  unsigned long p = 1, n, i, k, muls = 0, deg;

  if (pn == 1)
    {
#if NEGATED_ROOTS == 0
      mpz_set_si (F[0], -1L);
#else
      mpz_set_ui (F[0], 1UL);
#endif
      return 0;
    }

  /* Find smallest prime factor p of pn */
  for (p = 2; pn % p != 0; p++);
  n = pn / p;

  deg = eulerphi(n);

  ASSERT(tmplen > 0);
  mpz_powm_ui (tmp[0], r, p, modulus);
#ifdef TESTDRIVE
  if (verbose > 1)
    gmp_printf ("Calling poly_from_roots_coprime(F, %Zd, %lu, tmp, %Zd)\n",
		tmp[0], n, modulus);
#endif
  poly_from_roots_coprime (F, tmp[0], n, tmp + 1, tmplen - 1, modulus);
#if defined(TESTDRIVE)
  for (i = 0; verbose > 1 && i < deg; i++)
    gmp_printf ("poly_from_roots_coprime(F, %Zd, %lu, tmp, %Zd) returned"
		" f_%ld = %Zd\n", tmp[0], n, modulus, i, F[i]);
#endif
  for (i = p - 1; i >= 1; i--)
    {
      ASSERT(tmplen > 2);
      mpz_powm_ui (tmp[0], r, i*n, modulus);
      mpz_set (tmp[1], tmp[0]); /* tmp[1] = r^(i*n*k) */
      for (k = 1; k <= deg; k++)
	{
#ifdef TESTDRIVE
	  if (verbose > 1)
	    gmp_printf ("Setting f_%lu = f_%lu * r^%lu = %Zd * %Zd = ", 
			i * deg - k, deg - k, i*n*k, F[deg - k], tmp[1]);
#endif
	  mpz_mul (tmp[2], F[deg - k], tmp[1]);
	  mpz_mod (F[i * deg - k], tmp[2], modulus);
#ifdef TESTDRIVE
	  if (verbose > 1)
	    gmp_printf ("%Zd\n", F[i * deg - k]);
#endif
	  mpz_mul (tmp[2], tmp[1], tmp[0]);
	  mpz_mod (tmp[1], tmp[2], modulus);
	  muls += 2;
	}
    }

#ifdef TESTDRIVE
  for (i = 1; verbose > 1 && i < p; i++)
    {
      gmp_printf ("%Zd^(%lu*%lu) * F_{%lu, %Zd}(x/%Zd^%lu) = ", 
		  r, i, deg, n, r, r, i);
      for (k = 0; k < deg; k++)
	gmp_printf ("+ %Zd *x^%lu", F[(i - 1)*deg + k], k);
      printf ("\n");
    }
#endif

  i = eulerphi(p); /* The number of pieces we have, hopefully a power of 2 */
  while (i > 1)
    {
      ASSERT (i % 2 == 0);
      for (k = 0; k < i; k += 2) /* Multiply pieces in pairs */
	{
	  ASSERT(tmplen > 2 * deg + list_mul_mem (deg));
	  list_mul (tmp, F + deg * k, deg, 1, 
		         F + deg * (k + 1), deg, 1,
		    tmp + 2 * deg);
	  list_mod (F + deg * k, tmp, 2 * deg, modulus);
	}
      /* Now there are half as many pieces, each twice as long */
      i /= 2;
      deg *= 2;
    }

#ifdef TESTDRIVE
  if (verbose > 1)
    printf ("Final polynomial at this stage is ");
  for (i = 0; verbose > 1 && i < eulerphi(p) * eulerphi(n); i++)
    gmp_printf ("f_%lu = %Zd\n", i, F[i]);
#endif

#ifdef TESTDRIVE
  if (verbose > 1)
    gmp_printf ("Leaving poly_from_roots_coprime(F, %Zd, %lu, tmp, %Zd)\n", 
		r, pn, modulus);
#endif

  return muls;
}


/* Generate sequence B for the evaluation of F(X^(alpha + beta*i)), 
   for 0 <= i < d */

static int 
sequence_B (mpz_t *B, mpres_t X, long d, long beta, mpmod_t modulus)
{
  /* Prepare the polynomial B of degree 2*d-1, but not necessarily monic. 
     Since this is invariant over the different blocks, we need to 
     compute it only once */
  
  /* We want b_j = X^(-beta*(j-d)^2/2) = X^(-2P*(j-d)^2/2)
     for j = 0 ... 2*d-1. This sequence is symmetric around j = d,
     so only compute one half, and copy the other. This b_j sequence
     is the same for all blocks, so we need to compute it only once.
     We can compute it with the finite differences of the exponents
     -(j+1-d)^2 - -(j-d)^2 = 2(d-j)-1,
     2(d-(j+1))-1 - 2(d-j)-1 = -2j
  */

  mpres_t XP,    /* X^{\beta / 2} */
          invXP; /* X^{-\beta / 2} */
  mpres_t findiff[3];
  mpz_t tmp;
  long j, timestart, timestop, P;

  ASSERT (beta % 2 == 0);
  P = beta / 2;
  outputf (OUTPUT_VERBOSE, "Computing coefficients of B");
  timestart = cputime ();
  mpz_init (tmp);
  mpres_init (XP, modulus);
  mpz_set_ui (tmp, P);
  mpres_pow (XP, X, tmp, modulus);  /* XP = X^P */
  mpres_init (invXP, modulus);
  mpres_invert (invXP, XP, modulus); /* invXP = X^{-P} */

  mpres_init (findiff[0], modulus);
  mpres_init (findiff[1], modulus);
  mpres_init (findiff[2], modulus);
  mpres_mul (findiff[0], invXP, invXP, modulus); 
    /* findiff[0] = X^{-2P} */
  
  mpz_set_ui (tmp, 2 * d - 1);
  mpres_pow (findiff[1], XP, tmp, modulus); 
    /* findiff[1] = X^{P(2*d-1)} */
  
  mpz_set_ui (tmp, d);
  mpz_mul (tmp, tmp, tmp); /* mt = d^2 (may exceed 2^32, so square in mpz_t) */
  mpres_pow (findiff[2], invXP, tmp, modulus); /* findiff[2] = X^(-Pd^2) */
  
  for (j = 0; j <= d; j++)
    {
      mpres_get_z (B[j], findiff[2], modulus);
      mpres_mul (findiff[2], findiff[2], findiff[1], modulus);
      mpres_mul (findiff[1], findiff[1], findiff[0], modulus);
    }
  
  /* B[d] = XP^(-(d-d)^2) = 1. Check that it is so */
  ASSERT(mpz_cmp_ui (B[d], 1) == 0);
  timestop = cputime ();
  outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);
  
  outputf (OUTPUT_VERBOSE, "Copying high half of B to low half");
  timestart = cputime ();
  /* Now mirror-copy the low half into the high half */
  for (j = 1; j < d; j++)
    mpz_set (B[d + j], B[d - j]);
  timestop = cputime ();
  outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);

  if (test_verbose (OUTPUT_TRACE))
    {
      for (j = 0; j < 2*d; j++)
	{
	  outputf (OUTPUT_TRACE, 
		   "B_%lu = XP^(-(%lu-d)^2); ", j, j);
	  outputf (OUTPUT_TRACE, 
		   "compare(B_%lu, Mod(%Zd, N), %ld); /* PARI %ld */\n", 
		   j, B[j], pariline, pariline); pariline++;
	}
  
      outputf (OUTPUT_TRACE, "B(x) = ");
      for (j = 0; j < 2 * d; j++)
	outputf (OUTPUT_TRACE, " + B_%lu*x^%lu", j, j);
      outputf (OUTPUT_TRACE, "; /* PARI %ld */\n", pariline++);
    }
  
  mpres_clear (findiff[2], modulus);
  mpres_clear (findiff[1], modulus);
  mpres_clear (findiff[0], modulus);
  mpres_clear (invXP, modulus);
  mpres_clear (XP, modulus);
  mpz_clear (tmp);

  return 0;
}


static int 
sequence_C (mpz_t *C, mpz_t *F, const mpres_t X, long d, const mpz_t alpha, 
	    const long beta, const long l, mpmod_t modulus)
{
      /* Prepare polynomial C. We want 
	 c_j = f_j * X^(beta * j^2/2 + alpha * j), j = 0 ... d
	 
	 We have beta*j + alpha + beta/2 and beta for the finite differences 
	 of the exponents of X.
      */

  mpres_t findiff[3];
  mpz_t tmp;
  long j, timestart, timestop;

  ASSERT (beta % 2 == 0);

  mpz_init (tmp);
  mpres_init (findiff[0], modulus);
  mpres_init (findiff[1], modulus);
  mpres_init (findiff[2], modulus);

  outputf (OUTPUT_VERBOSE, "Computing coefficients of C, block %lu", l);
  timestart = cputime ();
  mpz_set_ui (tmp, beta);
  mpres_pow (findiff[0], X, tmp, modulus); /* fd[0] = X^beta */
  mpz_add_ui (tmp, alpha, beta / 2);
  mpres_pow (findiff[1], X, tmp, modulus); /* fd[1] = X^(alpha + beta/2) */
  mpres_set_ui (findiff[2], 1, modulus); /* j=0, X^(beta*j^2/2 + alpha*j)=1 */
  /* Can we just init this once and let the findiff stuff continue
     over all the blocks? */
  
  mpz_set (C[0], F[0]); /* fd[2] = 1, so c_0 = f_0 */
  
  for (j = 1; j < d; j++)
    {
      mpres_mul (findiff[2], findiff[2], findiff[1], modulus);
      mpres_mul (findiff[1], findiff[1], findiff[0], modulus);
      mpres_get_z (tmp, findiff[2], modulus);
      mpz_mul (tmp, tmp, F[j]);
      mpz_mod (C[j], tmp, modulus->orig_modulus);
    }
  /* F[d] = 1 is not actually stored anywhere. Treat it separately */
  mpres_mul (findiff[2], findiff[2], findiff[1], modulus);
  mpres_get_z (C[j], findiff[2], modulus);
  timestop = cputime ();
  outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);
  
  if (test_verbose (OUTPUT_TRACE))
    {
      for (j = 0; j <= d; j++)
	{
	  outputf (OUTPUT_TRACE, 
		   "C_%lu = F_%lu * XP^(%lu^2 + 2*(%lu*d+%Zd)*%lu); ", 
		   j, j, j, l, alpha, j);
	  outputf (OUTPUT_TRACE, "compare(C_%lu, %Zd, %ld); /* PARI %ld */\n", 
		   j, C[j], pariline, pariline); pariline++;
	}
      
      outputf (OUTPUT_TRACE, "C(x) = ");
      for (j = 0; j <= d; j++)
	outputf (OUTPUT_TRACE, " + C_%lu*x^%lu", j, j);
      outputf (OUTPUT_TRACE, "; /* PARI %ld */\n", pariline++);
    }

  mpres_init (findiff[2], modulus);
  mpres_init (findiff[1], modulus);
  mpres_init (findiff[0], modulus);
  mpz_clear (tmp);

  return 0;
}

int 
pm1fs2(mpz_t f, mpres_t X, mpmod_t modulus, root_params_t *root_params, 
       unsigned long dF_param, unsigned long blocks)
{
  /* Polynomial F has roots X^i for one i from each class coprime to 
     root_params->d1, so has degree dF = eulerphi(d1) */
  /* We want to evaluate the polynomial F of degree dF at nr points, 
     where dF is the actual degree of the polynomial (not rounded
     up to a power of 2). We probably want len = dF + nr to be a power 
     of 2, so we use dF_param, which probably is a power of 2 that's >= dF/2,
     and quadruple that to get len. */

  const unsigned long beta = root_params->d1;
  const unsigned long dF = eulerphi (beta);
  const unsigned long len = 4 * dF_param;
#if 1
  /* For now, we require that the number of points of evaluation and the
     degree of the polynomial (rounded up if necessary) are the same, 
     each len / 2 */
  const unsigned long nr = len / 2;
#else
  /* In the future, we'll want as many points of evaluation as we can
     fit in a convolution of length len */
  const unsigned long nr = len - dF;
#endif

  unsigned long i, j, l, lenF, lenB, lenC, lenR, lentmp, setsize;
  listz_t F, B, C, tmp, R;
  mpres_t Xi, X2;
  mpz_t mt;
  int youpi = 0;
  extern unsigned int Fermat;
  long timetotalstart, timestart, timestop, muls;
  long *sets;

  ASSERT (dF <= len / 2);
  ASSERT  (beta % 2 == 0);
  
  if (modulus->repr == ECM_MOD_BASE2 && modulus->Fermat > 0)
    Fermat = modulus->Fermat;

  timetotalstart = cputime ();

  outputf (OUTPUT_VERBOSE, "P = %lu; dF = %lu; /* PARI %ld */\n", 
	   beta / 2, dF, pariline++);
  outputf (OUTPUT_TRACE, "compare(a, b, n) = if (a != b,print(\"In PARI \", n, \" line, \", a \" != \" b));"
	   "/* PARI %ld */\n", pariline++);

  /* Allocate all the memory we'll need */
  /* Allocate the correct amount of space for each mpz_t or the 
     reallocations will up to double the time for stage 2! */
  mpz_init (mt);    /* All-purpose temp mpz_t */
  lenF = len / 2;   /* This should be dF sometime later */
  F = init_list2 (lenF, labs (modulus->bits));
  lenB = len;
  B = init_list2 (lenB, labs (modulus->bits));
  lenC = lenF + 1; /* Like F, but with leading monomial actually stored */
  C = init_list2 (lenC, labs (modulus->bits));
  lenR = nr;
  R = init_list2 (lenR, labs (modulus->bits));    
  lentmp = len + list_mul_mem (len / 2);
  outputf (OUTPUT_DEVVERBOSE, "lentmp = %lu\n", lentmp);
  if (TMulGen_space (lenC - 1, lenB - 1, lenR) + 12 > lentmp)
    lentmp = TMulGen_space (lenC - 1, lenB - 1, lenR) + 12;
  /* FIXME: It appears TMulGen_space() returns a too small value! */
  outputf (OUTPUT_DEVVERBOSE, "lentmp = %lu\n", lentmp);
  tmp = init_list2 (lentmp, labs (modulus->bits));
  
  mpres_get_z (mt, X, modulus); /* mpz_t copy of X for printing */
  outputf (OUTPUT_RESVERBOSE, 
	   "N = %Zd; X = Mod(%Zd, N); XP=X^P; /* PARI %ld */\n", 
	   modulus->orig_modulus, mt, pariline++);
  
#if 1
  /* Build polynomial for fully factored set of residues coprime to beta */
  outputf (OUTPUT_VERBOSE, "Factoring and sorting sets");
  timestart = cputime ();
  setsize = factor_coprimeset (NULL, beta);
  sets = malloc (setsize * sizeof (long));
  ASSERT(sets != NULL); /* FIXME, error handling */
  i = factor_coprimeset (sets, beta);
  ASSERT(i == setsize);
  sort_sets (sets, setsize);
  timestop = cputime ();
  outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);
  
  if (test_verbose (OUTPUT_DEVVERBOSE))
  {
      outputf (OUTPUT_DEVVERBOSE, "The fully factored, sorted sets are ");
      print_sets (OUTPUT_DEVVERBOSE, sets, setsize);
  }

  mpz_mul_ui (mt, root_params->i0, beta);
  mpz_add_si (mt, mt, sum_sets_minmax (sets, setsize, 1));
  outputf (OUTPUT_VERBOSE, "Effective B2min = %Zd\n", mt);
  mpz_add_ui (mt, root_params->i0, blocks * nr);
  mpz_mul_ui (mt, mt, beta);
  mpz_sub_si (mt, mt, sum_sets_minmax (sets, setsize, -1));
  outputf (OUTPUT_VERBOSE, "Effective B2max = %Zd\n", mt);

  outputf (OUTPUT_VERBOSE, "Computing F from factored set of units mod %lu",
	   beta);
  if (test_verbose (OUTPUT_DEVVERBOSE))
      outputf (OUTPUT_VERBOSE, "\n");

  mpres_get_z (mt, X, modulus);
  i = poly_from_sets (F, mt, sets, setsize, tmp, lentmp, 
		      modulus->orig_modulus);
  ASSERT (i == dF);

#if defined(WANT_ASSERT)
  if (sets[0] != 2) /* Unless the first set has one element, 
		       the polynomial should be symmetric */
    {
      long i = list_is_symmetric (F, dF, 1, 1, modulus->orig_modulus, mt);
      if (i != -1)
	{

	  outputf (OUTPUT_ERROR, 
		   "Polynomial not symmetric! F[%ld] != F[%ld]\n", i, dF - i);
	  list_output_poly (F, dF, 1, "F(x) = ", OUTPUT_ERROR);
	  outputf (OUTPUT_ERROR, "Factored sets: ");
	  print_sets (OUTPUT_ERROR, sets, setsize);
	  abort ();
	}
      else
	outputf (OUTPUT_DEVVERBOSE, "Polynomial is symmetric\n");
    }
#endif
  
  free (sets);
  sets = NULL;

  timestop = cputime ();
  outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);

  /*
    outputf (OUTPUT_VERBOSE, "Computing F via product of coprime sets");
    timestart = cputime ();

    poly_from_roots_coprime (F, X, beta, tmp, lentmp, modulus->orig_modulus);
    
    timestop = cputime ();
    outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);
  */

#else
  /* Build F the old way, computing all the roots and doing a product tree */
  mpz_mul_ui (mt, root_params->i0, beta);
  outputf (OUTPUT_VERBOSE, "Effective B2min ~= %Zd\n", mt);
  mpz_add_ui (mt, root_params->i0, blocks * nr);
  mpz_mul_ui (mt, mt, beta);
  outputf (OUTPUT_VERBOSE, "Effective B2max ~= %Zd\n", mt);

  outputf (OUTPUT_VERBOSE, "Computing roots of F");
  outputf (OUTPUT_TRACE, "\n"); /* So the "Computing" does not interfere */
  timestart = cputime ();
  mpres_init (X2, modulus);
  mpres_init (Xi, modulus);
  mpres_mul (X2, X, X, modulus); /* X2 = X^2 */
  mpres_set (Xi, X, modulus);    /* Xi = X^i */
  /* Prepare polynomial F(x), which is monic of degree dF. The leading
     monomial is not stored. */
  /* Put in F[0 .. dF-1] the values of X^i, 1<=i<beta, gcd(i, beta) == 1 */
  for (i = 1, j = 0; i < beta; i += 2) /* "+= 2" assumes beta is even */
    {
      if (gcd (i, beta) == 1)
	{
	  mpres_get_z (F[j], Xi, modulus);
	  outputf (OUTPUT_TRACE, "f_%lu = X^%lu;"
		   "compare(f_%lu, Mod(%Zd, N), %ld); /* PARI %ld */\n", 
		   j, i, j, F[j], pariline, pariline); pariline++;
	  j++;
	}
      mpres_mul (Xi, Xi, X2, modulus);
    }
  
  ASSERT(j == dF);
  mpres_clear (X2, modulus);
  mpres_clear (Xi, modulus);
  timestop = cputime ();
  outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);

  outputf (OUTPUT_TRACE, "F(x) = ");
  for (j = 0; j < dF - 1; j++)
    outputf (OUTPUT_TRACE, "(x - f_%lu) * ", j);
  outputf (OUTPUT_TRACE, "(x - f_%lu); /* PARI %ld */\n", dF - 1, pariline++);
  
  /* Multiply all the (x - f_i) to form F(x) in monomial basis */
  outputf (OUTPUT_VERBOSE, "Building F from its roots");
  timestart = cputime ();
  PolyFromRoots (F, F, dF, tmp, modulus->orig_modulus);
  timestop = cputime ();
  outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);
#endif
  
  if (test_verbose (OUTPUT_TRACE))
    list_output_poly (F, dF, 1, "F(x) == ", OUTPUT_TRACE);
  
  sequence_B (B, X, len / 2, beta, modulus);
  
  for (l = 0; l < blocks; l++)
    {
      /* Now the multipoint evaluation. We want to evaluate F(x) on
	 X^(beta*(i0 + l*len/2 + i)), for len/2 successive values of i */

      mpz_set_ui (mt, l);
      mpz_mul_ui (mt, mt, nr);   /* mt = l * nr */
      mpz_add (mt, mt, root_params->i0);     /* mt = i0 + l*nr */
      mpz_mul_ui (mt, mt, beta); /* mt = alpha = beta*(i0 + l*nr) */
      sequence_C (C, F, X, dF, mt, beta, l, modulus);

      /* Do the convolution */
#if 1
      /* Use the transposed "Middle Product" algorithm */
      /* TMulGen reverses the first input sequence, which we don't want.
	 We could fill C[] in reverse order, for now reverse it 
	 separately here. */
      outputf (OUTPUT_VERBOSE, "Swapping C\n");
      /* Remember that the C array has length len/2+1 ! */
      for (j = 0; j < lenC - 1 - j; j++)
	mpres_swap (C[j], C[lenC - 1 - j], modulus);

      outputf (OUTPUT_VERBOSE, "TMulGen of B and C");
      ASSERT(lentmp >= TMulGen_space (nr - 1, lenC - 1, lenB - 1));
      muls = TMulGen (R, nr - 1, C, lenC - 1, B, lenB - 1, tmp, 
		      modulus->orig_modulus);
      list_mod (R, R, nr, modulus->orig_modulus);

      timestop = cputime ();
      outputf (OUTPUT_VERBOSE, " took %lu muls, %lu ms\n", 
	       muls, timestop - timestart);

      /* Undo swap. Not doing so causes incorrect results!! */
      /* This is because we can have lenC > dF, so by swapping, some
	 non-zero values get swapped to the highest coefficients of C
	 that will not be overwritten by the next sequence_C() call. 
	 This bug took a while to find. :( */
	 
      outputf (OUTPUT_VERBOSE, "Un-Swapping C\n");
      for (j = 0; j < lenC - 1 - j; j++)
	mpres_swap (C[j], C[lenC - 1 - j], modulus);
#else
      /* Use two non-transposed multiplications */
      outputf (OUTPUT_VERBOSE, "Computing B * C");
      timestart = cputime ();
      list_mul (tmp, C, lenF + 1, 0, B, len/2, 0, tmp + lenF + 1 + len/2);
      list_set (R, tmp + len/2, nr);
      
      /* len/2 is enough for the length of C here, as the coefficient in
	 C[len/2] will only affect R[i], i >= len, which we don't need */
      list_mul (tmp, C, lenF, 0, B + len/2, len/2, 0, tmp + lenF + len/2);
      list_add (R, R, tmp, nr);
      list_mod (R, R, nr, modulus->orig_modulus);
      timestop = cputime ();
      outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);
#endif
#if defined(WANT_ASSERT)
      if (l == 0 && mpz_cmp_ui (root_params->i0, 0) == 0)
	{
	  /* i0 == 0 so we evaluated F(X^0) = F(1). We can check that for
	     correctness easily. */
	  mpz_set_ui (mt, 0);
	  for (i = 0; i < dF; i++)
	    mpz_add (mt, mt, F[i]);
	  /* Add the leading 1 term which is not stored anywhere */
	  mpz_add_ui (mt, mt, 1);
	  mpz_mod (mt, mt, modulus->orig_modulus);
	  outputf (OUTPUT_DEVVERBOSE, "Testing R[0] = F(1)\n");
	  outputf (OUTPUT_TRACE, "%Zd == %Zd /* PARI */\n", R[0], mt);
	  ASSERT (mpz_cmp (R[0], mt) == 0);
	}
#endif

#if 0 && defined(WANT_ASSERT)

      /* See if R[i] is correct, with a test that works even if i0 != 0 */
      /* More expensive self-test */
      /* alpha = beta*(i0 + l*nr) */

      outputf (OUTPUT_VERBOSE, "Verifying all results (slow)");
      for (i = 0; i < nr; i++)
	{
	  mpz_set_ui (mt, nr * l);
	  mpz_add (mt, mt, root_params->i0);
	  mpz_add_ui (mt, mt, i);
	  mpz_mul_ui (mt, mt, beta);
	  mpres_get_z (tmp[0], X, modulus);
	  mpz_powm (tmp[0], tmp[0], mt, modulus->orig_modulus);
	  /* Hence, tmp[0] = X^(alpha + i * beta) */
	  list_eval_poly (tmp[1], F, tmp[0], dF, 1, modulus->orig_modulus, 
			  tmp + 2);

	  mpz_set_ui (mt, i);
	  mpz_mul_ui (mt, mt, i);
	  mpz_mul_ui (mt, mt, beta / 2); /* h(i) = beta*i^2/2 */
	  mpres_get_z (tmp[0], X, modulus);
	  mpz_powm (tmp[0], tmp[0], mt, modulus->orig_modulus); /* X^h(1) */
	  mpz_mul (tmp[0], tmp[0], R[i]);
	  mpz_mod (tmp[0], tmp[0], modulus->orig_modulus);
	  if (mpz_cmp (tmp[0], tmp[1]) != 0)
	    {
	      outputf (OUTPUT_ERROR, "Result in R[%ld] incorrect.\n", i);
	      outputf (OUTPUT_ERROR, "R[%ld] = %Zd\n", i, R[i]);
	      abort ();
	    }
	}
      outputf (OUTPUT_VERBOSE, " - everything's correct! :-D\n");
#endif

      outputf (OUTPUT_VERBOSE, "Computing product of F(g_i)");
      timestart = cputime ();
#if 1
      /* Try a faster way of multiplying up the R[i] for the gcd */
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
	mpres_gcd (tmp[0], tmpprod, modulus);
	mpres_clear (tmpprod, modulus);
	mpres_clear (tmpres, modulus);
      }
#else
      list_mulup (R, nr, modulus->orig_modulus, tmp[0]); 
      mpz_gcd (tmp[0], R[nr - 1], modulus->orig_modulus);
#endif
      timestop = cputime ();
      outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);

      if (mpz_cmp_ui (tmp[0], 1) > 0)
	{
	  outputf (OUTPUT_NORMAL, "Found factor in stage 2\n");
	  mpz_set (f, tmp[0]);
	  youpi = ECM_FACTOR_FOUND_STEP2;
	  break;
	}
      outputf (OUTPUT_RESVERBOSE, "Product of F(g_i) = %Zd\n", R[nr - 1]);
    }
  
  clear_list (F, lenF);
  clear_list (B, lenB);
  clear_list (C, lenC);
  clear_list (R, lenR);    
  clear_list (tmp, lentmp);

  timestop = cputime ();
  outputf (OUTPUT_NORMAL, "Step 2 took %ld ms\n", 
           timestop - timetotalstart);
  
  return youpi;
}

#ifdef TESTDRIVE
int main (int argc, char **argv)
{
  unsigned long pn, d, i, j, tmplen, setsize;
  listz_t F, tmp;
  mpz_t r, modulus;
  long *L;
  int selftest = 0, pari = 0;
  

  ECM_STDOUT = stdout;
  ECM_STDERR = stderr;
  set_verbose (OUTPUT_DEVVERBOSE);

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
      else
	break;
      argc--;
      argv++;
    }
  if (argc > 1)
    pn = strtoul (argv[1], NULL, 10);

  d = eulerphi (pn);

  setsize = factor_coprimeset (NULL, pn);
  L = malloc (setsize * sizeof (long));
  i = factor_coprimeset (L, pn);
  ASSERT (i == setsize);
  
#if 1
  if (verbose)
    {
      printf ("Factored sets before sorting\n");
      print_sets (OUTPUT_DEVVERBOSE, L, setsize);
    }
#endif

  sort_sets (L, setsize);

  if (verbose)
    {
      printf ("Factored sets after sorting\n");
      print_sets (OUTPUT_DEVVERBOSE, L, setsize);
    }

  F = init_list (d);
  tmplen = 5*d;
  tmp = init_list (tmplen);
  mpz_init (r);
  mpz_init (modulus);
  mpz_set_ui (r, 3);
  mpz_ui_pow_ui (modulus, 2, 31);
  mpz_sub_ui (modulus, modulus, 1);
  if (pari)
    gmp_printf ("N = %Zd; r = Mod(%Zd, N);\n", modulus, r);

#if 0
  poly_from_roots_coprime (F, r, pn, tmp, tmplen, modulus);
#else
  poly_from_sets (F, r, L, setsize, tmp, tmplen, modulus);
#endif

  if (pn % 4 != 2)
    for (i = 1; i < d / 2; i++)
      {
	if (mpz_cmp (F[i], F[d - i]) != 0)
	  {
	    printf ("Error, F[%lu] != F[%lu - %lu], F is not symmetric!\n",
		    i, d, i);
	  }
      }
  
  if (pari)
    {
      printf ("F(x) = x^%lu", d);
      for (i = d - 1; i > 0; i--)
	if (mpz_sgn (F[i]) != 0)
	  gmp_printf (" + %Zd * x^%lu", F[i], i);
      gmp_printf(" + %Zd\n", F[0]);
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
	  printf ("%ld];\n", sumset[d - 1]);
	  printf ("lift(prod(k=1,length(exponents),(x-r^exponents[k]))) "
		  "== F(x)\n");
	}
      
      mpz_invert (tmp[0], r, modulus);
      
      /* Check that all the elements in the sumset are really exponents of the
	 roots of F */
      for (i = 0; i < d; i++)
	{
	  if (sumset[i] < 0)
	    mpz_powm_ui (tmp[1], tmp[0], -sumset[i], modulus);
	  else
	    mpz_powm_ui (tmp[1], r, sumset[i], modulus);
	  list_eval_poly (tmp[2], F, tmp[1], d, 1, modulus, tmp + 3);
	  if (mpz_sgn (tmp[2]) != 0)
	    printf ("Error, r^%ld is not a root of F\n", sumset[i]);
	}
      
      /* Check that the set of sums is really a complete set of representatives
	 of residue classes coprime to pn */
      for (i = 0; i < d; i++)
	if (sumset[i] >= 0)
	  sumset[i] = sumset[i] % pn;
	else
	  sumset[i] = pn - (-sumset[i]) % pn;
      quicksort_long (sumset, d);
      for (i = 1, j = 0; i < pn; i++)
	if (gcd (i, pn) == 1)
	  {
	    ASSERT((unsigned long) sumset[j] == i);
	    j++;
	  }
      
      free (sumset);
    }
  
  return 0;
}
#endif
