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

   V(i,X)={ if (i==0, return(2)); if (i==1, return(X)); if(i%2 == 0, return (V (i/2, X)^2-2)); return (V ((i+1)/2, X) * V ((i-1)/2, X) - X)}

   U(i,X)={if(i==0,return(0));if(i==1,return(1));if(i%2==0,return(U(i/2,X)*V(i/2,X)));return(V((i+1)/2,X)*U((i-1)/2,X)+1)}

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
sum_sets (long *sum, long *sets, unsigned long setsize, long add)
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


/* Exchange two adjacent sets in memory. The sets have l and m elements of
   unsigned long type, respectively. In our case, that includes the size 
   value at the start of each set. */

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
  
  if (test_verbose (OUTPUT_DEVVERBOSE))
  {
    outputf (OUTPUT_DEVVERBOSE, "The fully factored, sorted sets are for "
	     "beta = %lu are ", beta);
      print_sets (OUTPUT_DEVVERBOSE, sets, size);
  }
  
  if (setsize != NULL)
    *setsize = size;
  
  return sets;
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
		mpz_t modulus, listz_t tmp, const unsigned long tmplen)
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
		    listz_t tmp, const unsigned long tmplen)
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
    mpz_set_ui (t1[i], 0);           /* t1[2*l1-1 ... 2*lmax-2] = 0 */
  
  /* Same for S_2(x) */
  for (i = 0; i < l2; i++)
    mpz_set (t2[i], S2[l2 - 1 - i]);
  for (i = 1; i < l2; i++)
    mpz_set (t2[l2 - 1 + i], S2[i]);
  for (i = 2 * l2 - 1; i <= 2 * lmax - 2; i++)
    mpz_set_ui (t2[i], 0);
  
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
		 listz_t tmp, const unsigned long tmplen)
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
      mpres_set_ui (R, 2, modulus);
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
	ASSERT(mpz_sgn (G[i]) >= 0);
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
      mpz_mul_2exp (newtmp[0], H[0], 1);
      mpz_neg (newtmp[0], newtmp[0]);
      mpz_set_ui (newtmp[1], 0);
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
      mpz_mul_2exp (newtmp[0], newtmp[0], 1); /* 2*(h_2 - h_0) */
      mpz_neg (newtmp[1], H[1]);              /* -h_1 */
      mpz_mul_2exp (newtmp[2], H[2], 1);
      mpz_sub (newtmp[2], H[0], newtmp[2]);   /* h_0 - 2*h_2 */
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
      mpz_mul_2exp (newtmp[0], newtmp[0], 1); /* t[0] = 2*(h_0 + h_2) */
      
      mpz_sub (newtmp[1], H[3], H[1]); /* t[1] = -h_1 + h_3 */
      
      for (i = 2; i <= 2 * deg - 4; i++)
	{
	  mpz_add (newtmp[i], H[i-2], H[i+2]);
	  mpz_sub (newtmp[i], newtmp[i], H[i]); /* t[i] = h_{i-2}-2h_i+h_{i+2} */
	  mpz_sub (newtmp[i], newtmp[i], H[i]); /* for 2 <= i <= n-2 */
	}
      
      mpz_mul_2exp (newtmp[2 * deg - 3], H[2 * deg - 3], 1);
      mpz_sub (newtmp[2 * deg - 3], H[2 * deg - 5], newtmp[2 * deg - 3]); 
      /* t[n-1] = h_{n-3} - 2h_{n-1} */
      
      mpz_mul_2exp (newtmp[2 * deg - 2], H[2 * deg - 2], 1);
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


/* Check if l is an (anti-)symmetric, possibly monic, polynomial. 
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
		mpz_t *tmp, const unsigned long tmplen, mpz_t modulus)
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
   all sets must be symmetric around 0. */

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

  if (test_verbose (OUTPUT_DEVVERBOSE))
    {
      mpz_t t;
      mpz_init (t);
      mpres_get_z (t, Q, modulus);
      outputf (OUTPUT_DEVVERBOSE, 
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

  sets += 3;
  setsize -= 3;

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

  return deg;
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


/* Build polynomial F(x) with roots X^i for i covering all the residue classes
   coprime to beta. F must have space for eulerphi(beta) coefficients.
   method can be 0, 1 or 2, which mean: 0 = old way of computing all the 
   roots and doing a product tree, 1 = using recursive expansion of polynomial
   *without* Chebychev polynomials to utilize symmetry, 2 = using recursive 
   expansion of polynomial *with* Chebychev polynomials */

static void 
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
      if (test_verbose (OUTPUT_DEVVERBOSE))
	outputf (OUTPUT_VERBOSE, "\n");
      
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
	  if (gcd (i, beta) == 1)
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


int 
pm1fs2 (mpz_t f, mpres_t X, mpmod_t modulus, root_params_t *root_params, 
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

  unsigned long i, j, l, lenF, lenB, lenC, lenR, tmplen;
  listz_t F, B, C, tmp, R;
  mpz_t mt;
  int youpi = 0;
  extern unsigned int Fermat;
  long timetotalstart, timestart, timestop, muls;

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
  lenF = len / 2;   /* FIXME: This should be dF sometime later */
  F = init_list2 (lenF, labs (modulus->bits));
  lenB = len;
  B = init_list2 (lenB, labs (modulus->bits));
  lenC = lenF + 1; /* Like F, but with leading coefficient actually stored */
  C = init_list2 (lenC, labs (modulus->bits));
  lenR = nr;
  R = init_list2 (lenR, labs (modulus->bits));    
  tmplen = 4 * len + list_mul_mem (len / 2);
  outputf (OUTPUT_DEVVERBOSE, "tmplen = %lu\n", tmplen);
  if (TMulGen_space (lenC - 1, lenB - 1, lenR) + 12 > tmplen)
    {
      tmplen = TMulGen_space (lenC - 1, lenB - 1, lenR) + 12;
      /* FIXME: It appears TMulGen_space() returns a too small value! */
      outputf (OUTPUT_DEVVERBOSE, "With TMulGen_space, tmplen = %lu\n", 
	       tmplen);
    }
  tmp = init_list2 (tmplen, labs (modulus->bits));
  
  mpres_get_z (mt, X, modulus); /* mpz_t copy of X for printing */
  outputf (OUTPUT_RESVERBOSE, 
	   "N = %Zd; X = Mod(%Zd, N); XP=X^P; /* PARI %ld */\n", 
	   modulus->orig_modulus, mt, pariline++);
  
  pm1_build_poly_F (F, X, modulus, beta, 0, root_params->i0, nr, 
		    blocks, tmplen, tmp);
  
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
      timestart = cputime ();
      ASSERT(tmplen >= TMulGen_space (nr - 1, lenC - 1, lenB - 1));
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
  clear_list (tmp, tmplen);

  timestop = cputime ();
  outputf (OUTPUT_NORMAL, "Step 2 took %ld ms\n", 
           timestop - timetotalstart);
  
  return youpi;
}


/* Multiplies (a0 + a1*sqrt(Delta)) * (b0 + b1*sqrt(Delta))
   using four multiplications. Result goes is (r0 + r1*sqrt(Delta)). 
   a0, a1, b0, b1, r0, r1 may overlap arbitrarily. t[0], t[1], t[2] and Delta
   must not overlap with anything. */
/* FIXME: is there a faster multiplication routine if both inputs have 
   norm 1? */

static void 
gfp_ext_mul (mpz_t r0, mpz_t r1, const mpz_t a0, const mpz_t a1,
	     const mpz_t b0, const mpz_t b1, const mpz_t Delta, 
	     const mpz_t N, const unsigned long tmplen, mpz_t *t)
{
  ASSERT (tmplen >= 3);
  mpz_add (t[0], a0, a1);
  mpz_add (t[1], b0, b1);
  mpz_mul (t[2], t[0], t[1]); /* t2 = (a0 + a1) * (b0 + b1) = a0*b0 + a0*b1 + 
				 a1*b0 + a1*b1 */

  mpz_mul (t[0], a0, b0);     /* t[0] = a0*b0. We don't need a0 or b0 now */
  mpz_sub (t[2], t[2], t[0]); /* r1 = a0*b1 + a1*b0 + a1*b1 */
  
  mpz_mul (t[1], a1, b1);     /* t[1] = a1*b1. We don't need a1 or b1 now */
  mpz_sub (t[2], t[2], t[1]);
  mpz_mod (r1, t[2], N);      /* r1 == a0*b1 + a1*b0 */

  mpz_mod (t[2], t[1], N);    /* t[2] = a1*b1 % N */
  mpz_mul (t[1], t[2], Delta); /* t[1] = a1*b1*Delta */
  mpz_add (t[0], t[0], t[1]);  /* t[0] = a0*b0 + a1*b1*Delta */
  mpz_mod (r0, t[0], N);       /* r0 = (a0*b0 + a1*b1*Delta) % N */
}

static void
gfp_ext_sqr_norm1 (mpz_t r0, mpz_t r1, const mpz_t a0, const mpz_t a1,
		   const mpz_t N, const unsigned long tmplen, mpz_t *t)
{
  ASSERT (tmplen >= 1);
  mpz_mul (t[0], a0, a1);
  mpz_mod (r1, t[0], N);
  mpz_mul_2exp (r1, r1, 1);
  
  mpz_mul (t[0], a0, a0);
  mpz_mod (r0, t[0], N);
  mpz_mul_2exp (r0, r0, 1);
  mpz_sub_ui (r0, r0, 1);
}


/* Raise (a0 + a1*sqrt(Delta)) to the power e. (a0 + a1*sqrt(Delta)) is 
   assumed to have norm 1, i.e. a0^2 - a1^2*Delta == 1. The result is 
   (r0 * r1*sqrt(Delta)). */

static void 
gfp_ext_pow_norm1 (mpz_t r0, mpz_t r1, const mpz_t a0, const mpz_t a1, 
		   const mpz_t Delta, const long e, const mpz_t N, 
		   unsigned long tmplen, mpz_t *tmp)
{
  const unsigned long abs_e = labs (e);
  unsigned long mask = ~0UL - (~0UL >> 1);

  if (e == 0)
    {
      mpz_set_ui (r0, 1UL);
      mpz_set_ui (r1, 0UL);
      return;
    }

  /* If e < 0, we want 1/(a0 + a1*sqrt(Delta)). By extending with 
     a0 - a1*sqrt(Delta), we get 
     (a0 - a1*sqrt(Delta)) / (a0^2 - a1^2 * Delta), but that denomiator
     is the norm which is known to be 1, so the result is 
     a0 - a1*sqrt(Delta). */

  while ((abs_e & mask) == 0UL)
    mask >>= 1;

  mpz_set (r0, a0);
  mpz_set (r1, a1);
  if (e < 0)
    mpz_neg (r1, r1);

  while (mask > 1UL)
    {
      gfp_ext_sqr_norm1 (r0, r1, r0, r1, N, tmplen, tmp);
      mask >>= 1;
      if (abs_e & mask)
	{
	  if (e < 0)
	    mpz_neg (r1, r1); /* We compute 1/(1/r * a) for r*1/a */
	  gfp_ext_mul (r0, r1, r0, r1, a0, a1, Delta, N, tmplen, tmp);
	  if (e < 0)
	    mpz_neg (r1, r1);
	}
    }
}



int 
pp1fs2 ()
{
  return 0;
}



#ifdef TESTDRIVE

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

int main (int argc, char **argv)
{
  unsigned long pn, d, i, j, tmplen, setsize;
  listz_t F, tmp;
  mpz_t r, N;
  long *L;
  int selftest = 0, pari = 0;
  mpmod_t modulus;
  

  ECM_STDOUT = stdout;
  ECM_STDERR = stderr;
  set_verbose (OUTPUT_DEVVERBOSE);

  mpz_init (N);

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
      else
	break;
      argc--;
      argv++;
    }
  if (argc > 1)
    pn = strtoul (argv[1], NULL, 10);

  d = eulerphi (pn);

  L = get_factored_sorted_sets (&setsize, pn);

  F = init_list (d);
  tmplen = 10*d+10;
  if (tmplen < 2000)
    tmplen = 2000;
  tmp = init_list (tmplen);
  mpz_init (r);
  mpz_set_ui (r, 3);
  if (mpz_sgn (N) == 0)
    {
      /* By default, use the Mersenne prime 2^31-1 as the modulus */
      mpz_set_ui (N, 1);
      mpz_mul_2exp (N, N, 31);
      mpz_sub_ui (N, N, 1);
    }
  mpmod_init (modulus, N, 0);
  /* We don't need N anymore now */
  mpz_clear (N);
  if (pari)
    gmp_printf ("N = %Zd; r = Mod(%Zd, N);\n", modulus->orig_modulus, r);

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
    mpres_set_ui (Q, 2, modulus); /* This corresponds to Q = 1 + 1/1, 
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
      if (mpz_cmp_ui (F[0], 1) != 0)
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
      gmp_printf(" + %Zd\n", F[0]);
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
      
      mpz_invert (tmp[0], r, modulus->orig_modulus);
      
      /* Check that all the elements in the sumset are really exponents of the
	 roots of F */
      for (i = 0; i < d; i++)
	{
	  if (sumset[i] < 0)
	    mpz_powm_ui (tmp[1], tmp[0], -sumset[i], modulus->orig_modulus);
	  else
	    mpz_powm_ui (tmp[1], r, sumset[i], modulus->orig_modulus);
	  list_eval_poly (tmp[2], F, tmp[1], d, 1, modulus->orig_modulus, 
			  tmp + 3);
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

  mpz_clear (r);
  
  return 0;
}
#endif
