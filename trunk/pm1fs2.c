
#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "ecm-impl.h"
#include "sp.h"
#include <math.h>
#if HAVE_ALLOCA_H
#include <alloca.h>
#endif
#if HAVE_STRING_H
#include <string.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

/* Define TEST_ZERO_RESULT to test if any result of the multipoint
   evaluation is equal to zero. If the modulus is composite, this
   happening might indicate a problem in the evalutaion code */
#define TEST_ZERO_RESULT

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
    453888435UL, 470434965UL, 480134655UL, 510765255UL, 522506985UL,
    557732175UL, 570855285UL, 596530935UL, 610224615UL, 627912285UL,
    654729075UL, 703227525UL, 722116395UL, 751725975UL, 780825045UL,
    790524735UL, 821665845UL, 851275425UL, 863017155UL, 909984075UL,
    936020085UL, 984518535UL, 1017041025UL, 1052416365UL, 1086110025UL,
    1120314195UL, 1191785595UL, 1213887675UL, 1265809545UL, 1282356075UL,
    1331995665UL, 1391905515UL, 1423857435UL, 1479202725UL, 1547100555UL,
    1567520955UL, 1673196525UL, 1712565855UL, 1767130365UL, 1830673845UL,
    1883166285UL, 1973886915UL, 2001274275UL, 2109682575UL, 2187280095UL,
    2255177925UL, 2342475135UL, 2389442055UL, 2449171725UL, 2553826275UL,
    2647760115UL, 2788660875UL, 2808060255UL, 2953555605UL, 3234846615UL,
    3457939485UL, 3516137625UL, 3681032355UL, 3758629875UL, 3904125225UL,
    4127218095UL
#if (ULONG_MAX > 4294967295)
  , 4360010655UL, 4573403835UL, 4796496705UL, 4844995155UL,
    5019589575UL, 5203883685UL, 5262081825UL, 5465775315UL, 5766465705UL,
    5898837945UL, 6164152995UL, 6464843385UL, 6804332535UL, 6980458485UL,
    7320968655UL, 7526103585UL, 7919796885UL, 8142889755UL, 8589075495UL,
    8906102205UL, 9171056895UL, 9704539845UL, 9927632715UL, 10373818455UL
#endif
};

/* All the prime factors that can appear in eulerphi(P) */
const unsigned long phiPfactors[] = {2UL, 3UL, 5UL, 7UL, 11UL, 13UL};


/* Some useful PARI functions:
   sumset(a,b) = {local(i, j, l); l = listcreate (length(a) * length(b)); for (i = 1, length(a), for (j = 1, length(b), listput(l, a[i] + b[j]))); listsort (l, 1); l}

   V(i,X) = { if (i==0, return(2)); if (i==1, return(X)); if(i%2 == 0, return (V (i/2, X)^2-2)); return (V ((i+1)/2, X) * V ((i-1)/2, X) - X)}

   U(i,X) = { if (i==0, return(0)); if (i==1, return(1)); if(i%2 == 0, return (U (i/2, X) * V(i/2,X))); return (V ((i+1)/2, X)  *U( (i-1)/2, X) + 1)}
*/

void 
ntt_sqr_recip (mpzv_t, const mpzv_t, mpzspv_t, const spv_size_t, 
               const mpzspm_t);

static unsigned long
maxS (unsigned long P)
{
  unsigned long p, pk;
  unsigned int k;

  if (P == 1UL)
    return 0L;

  p = find_factor (P);
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


static void
factor_phiP (int *exponents, const unsigned long phiP)
{
    const int nrprimes = sizeof (phiPfactors) / sizeof (unsigned long);
    unsigned long cofactor = phiP;
    int i;
    
    ASSERT_ALWAYS (phiP > 0UL);

    for (i = 0; i < nrprimes; i++)
	for (exponents[i] = 0; cofactor % phiPfactors[i] == 0UL; exponents[i]++)
	    cofactor /= phiPfactors[i];

    ASSERT_ALWAYS (cofactor == 1UL);
}


static unsigned long 
pow_ul (const unsigned long b, const unsigned int e)
{
    unsigned long r = 1UL;
    unsigned int i;

    for (i = 0; i < e; i++)
	r *= b;

    return r;
}

static unsigned long
absdiff_ul (unsigned long a, unsigned long b)
{
    return (a > b) ? a - b : b - a;
}

/* Choose s_1 so that s_1 * s_2 = phiP, s_1 is positive and even, 
   s_2 >= min_s2 and s_2 is minimal and abs(s_1 - l) is minimal 
   under those conditions. 
   Returns 0 if no such choice is possible */

static unsigned long 
choose_s_1 (const unsigned long phiP, const unsigned long min_s2,
	    const unsigned long l)
{
  const int nrprimes = sizeof (phiPfactors) / sizeof (unsigned long);
  int phiPexponents[nrprimes], exponents[nrprimes];
  unsigned long s_1 = 0UL, s_2 = 0UL, trys_1;
  int i;

  ASSERT_ALWAYS (phiP > 0 && phiP % 2 == 0);

  /* We want only even s_1. We divide one 2 out of phiP here... */
  factor_phiP (phiPexponents, phiP / 2);
  for (i = 0; i < nrprimes; i++)
      exponents[i] = 0;

  do {
      trys_1 = 2; /* ... and add a 2 here */
      for (i = 0; i < nrprimes; i++)
	  trys_1 *= pow_ul (phiPfactors[i], exponents[i]);
#if 0
      printf ("choose_s_1: Trying trys_1 = %lu\n", trys_1);
#endif
      /* See if it satisfies all the required conditions and is an 
	 improvement over the previous choice */
      if (phiP / trys_1 >= min_s2 && 
	  (s_2 == 0UL || phiP / trys_1 < s_2) && 
	  absdiff_ul (trys_1, l) < absdiff_ul (s_1, l))
      {
#if 0
	  printf ("choose_s_1: New best s_1 for phiP = %lu, min_s2 = %lu, "
		  "l = %lu : %lu\n", phiP, min_s2, l, trys_1);
#endif
	  s_1 = trys_1;
      }
      for (i = 0; i < nrprimes; i++)
      {
	  if (++(exponents[i]) <= phiPexponents[i])
	      break;
	  exponents[i] = 0;
      }
  } while (i < nrprimes);

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
  mpz_tdiv_q_ui (t, B2l, (lmax > 3UL) ? 2UL * lmax - 6UL : 1UL);
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
      s_1 = choose_s_1 (phiP, min_s2, lmax / 2);
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
      trys_1 = choose_s_1 (tryphiP, min_s2, l / 2);
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
		  char *prefix, char *suffix, int verbosity)
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
  outputf (verbosity, "%Zd", l[0]);
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


/* Multiply two reciprocal polynomials of degree 2*l1-2 and 2*l2-2, resp., 
   with coefficients in standard basis

   S_1(x) = S1[0] + sum_{1 \leq i \leq l1 - 1} S1[i] (x^i + x^{-i})
   S_2(x) = S2[0] + sum_{1 \leq i \leq l2 - 1} S2[i] (x^i + x^{-i})

   to the reciprocal polynomial of degree 2*(l1 + l2) - 4

   R(x) = R[0] + sum_{1 \leq i \leq l1 + l2 - 2} R[i] (x^i + x^{-i}) 
        = S_1(x) * S_2(x)

   R == S1 == S2 is permissible.
*/

static void
list_mul_reciprocal (listz_t R, listz_t S1, const unsigned long l1, 
		     listz_t S2, const unsigned long l2,
		     mpz_t modulus, listz_t tmp, 
		     ATTRIBUTE_UNUSED const unsigned long tmplen)
{
  unsigned long i, lmax;

  lmax = (l1 > l2) ? l1 : l2;

  if (l1 == 0UL || l2 == 0UL)
    return;

  if (l1 == l2)
    {
	/* FIXME: This modifies the input arguments. */
	/* We have to divide S1[0] and S2[0] by 2 */

	/* Assume l1 = l2 = 2, S1 = f0 + f1 * x, S2 = g0 + g1 * x */
	/* We want the coefficients at non-negative powers of x in 
	   (f0 + f1 * (x + 1/x)) * (g0 + g1 * (x + 1/x)), that means
	   (g0*f0 + 2*g1*f1) + (g1*f0 + g0*f1) * x + g1*f1 * x^2 */
	listz_t S2rev, r1 = tmp, r2 = tmp + 2 * l1 - 1, t = tmp + 4 * l1 - 2;

	ASSERT (tmplen >= 4 * lmax - 2 + list_mul_mem (lmax));

#if 0
	gmp_printf ("/* list_mul_reciprocal */ S1(x) = %Zd", S1[0]);
	for (i = 1; i < l1; i++)
	    gmp_printf (" + %Zd * (x^%lu + 1/x^%lu)", S1[i], i, i);
	gmp_printf ("\n");
	gmp_printf ("/* list_mul_reciprocal */ S2(x) = %Zd", S2[0]);
	for (i = 1; i < l2; i++)
	    gmp_printf (" + %Zd * (x^%lu + 1/x^%lu)", S2[i], i, i);
	gmp_printf ("\n");
#endif

	if (mpz_odd_p (S1[0]))
	  {
	    ASSERT_ALWAYS (mpz_odd_p (modulus));
	    mpz_add (S1[0], S1[0], modulus);
	  }
	mpz_tdiv_q_2exp (S1[0], S1[0], 1UL);
	if (S1 != S2)
	  {
	    if (mpz_odd_p (S2[0]))
	      {
		ASSERT_ALWAYS (mpz_odd_p (modulus));
		mpz_add (S2[0], S2[0], modulus);
	      }
	    mpz_tdiv_q_2exp (S2[0], S2[0], 1UL);
	  }
	
	list_mul (r1, S1, l1, 0, S2, l2, 0, t);
	/* r1 = f0*g0/4 + (f0*g1 + f1*g0)/2 * x + f1*g1 * x^2 */
#if 0
	for (i = 0; i < 2 * l1 - 1; i++)
	    gmp_printf ("list_mul_reciprocal: r1[%lu] = %Zd\n", i, r1[i]);
#endif

	/* Remember that S1 == S2 is possible */
	S2rev = (listz_t) malloc (l2 * sizeof (mpz_t));
	ASSERT_ALWAYS (S2rev != NULL);
	for (i = 0UL; i < l2; i++)
	    (*S2rev)[i] = (*S2)[l2 - 1UL - i];
	list_mul (r2, S1, l1, 0, S2rev, l2, 0, t);
	free (S2rev);
	/* r2 = g1*f0/2 + (g0*f0/4 + g1*f1) * x + g0*f1/2 * x^2 */
#if 0
	for (i = 0; i < 2 * l1 - 1; i++)
	    gmp_printf ("list_mul_reciprocal: r2[%lu] = %Zd\n", i, r2[i]);
#endif

	mpz_mul_2exp (r1[0], r1[0], 1UL);
	/* r1 = f0*g0/2 + (f0*g1 + f1*g0)/2 * x + f1*g1 * x^2 */
	for (i = 0; i < l1; i++)
	    mpz_add (r1[i], r1[i], r2[i + l1 - 1]);
	/* r1 = 3/4*f0*g0 + g1*f1 + (f0*g1 + 2*f1*g0)/2 * x + f1*g1 * x^2 */
	for (i = 0; i < l1; i++)
	    mpz_add (r1[i], r1[i], r2[l1 - i - 1]);
	/* r1 = f0*g0 + 2*g1*f1 + (f0*g1 + f1*g0) * x + f1*g1 * x^2 */
	for (i = 0; i < 2*l1 - 1; i++)
	    mpz_set (R[i], r1[i]);

	if (R != S1)
	    mpz_mul_2exp (S1[0], S1[0], 1UL);
	if (S1 != S2 && R != S2)
	    mpz_mul_2exp (S2[0], S2[0], 1UL);
	
#if 0
	for (i = 0; i < 2*l1; i++)
	    gmp_printf ("list_mul_reciprocal: R[%lu] = %Zd\n", i, R[i]);
#endif
    }
  else
    {
      /* More difficult case, the lengths are different. We just do a full
	 multiply and take the coefficients we want from that. This is not
	 very efficient, but it'll only happen during building the polynomial
	 and for sets of odd cardinality, i.e. when the polynomials to be
	 multiplied are still quite small. The inefficiency of the code here 
	 does not really matter too much. */
      const unsigned long dsum = l1 + l2 - 2; /* Half the degree of prod. */
      listz_t t1, t2, r;

      ASSERT (tmplen >= 8 * lmax - 2 + list_mul_mem (2 * lmax - 1));

      t1 = tmp;
      t2 = tmp + 2 * lmax - 1;
      r = tmp + 4 * lmax - 2;
      
      /* Full copy of S_1(x). S1 = [1,2,3,4] => t1 = [4,3,2,1,2,3,4]
	 There are 2*l1 - 1 coefficients in monomial basis, which go in 
	 t1[0 ... 2*l1-2]. We pad the high end with zeros up to t1[2*lmax-2] */
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
      
      for (i = 0; i <= dsum; i++)
	  mpz_set (R[i], r[dsum + i]);
    }
}


/* Same, but does squaring which makes things easier */

static void
list_sqr_reciprocal (listz_t R, listz_t S, const unsigned long l, 
		     mpz_t modulus, listz_t tmp, 
		     ATTRIBUTE_UNUSED const unsigned long tmplen)
{
  unsigned long i;
  listz_t Srev, r1 = tmp, r2 = tmp + 2 * l - 1, t = tmp + 4 * l - 2;

  if (l == 0UL)
    return;

  /* FIXME: This modifies the input arguments. */
  /* We have to divide S[0] by 2 */

  ASSERT (tmplen >= 4 * l - 2 + list_mul_mem (l));

#if 0
  gmp_printf ("/* list_sqr_reciprocal */ S(x) = %Zd", S[0]);
  for (i = 1; i < l1; i++)
    gmp_printf (" + %Zd * (x^%lu + 1/x^%lu)", S[i], i, i);
  gmp_printf ("\n");
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

  if (k == 0L)
    {
      mpres_set_ui (R, 2UL, modulus);
      return;
    }

  uk = labs (k);

  if (uk == 1UL)
    {
      mpres_set (R, S, modulus);
      return;
    }

  mpres_init (Vi, modulus);
  mpres_init (Vi1, modulus);

  for (j = 1UL; j <= uk / 2UL; j <<= 1);
  ASSERT ((uk & j) > 0UL);

  j >>= 1;
  i = 1UL;
  mpres_set (Vi, S, modulus);
  mpres_mul (Vi1, S, S, modulus);
  mpres_sub_ui (Vi1, Vi1, 2, modulus);

  while (j)
    {
      if ((uk & j) != 0UL)
	{
	  /* i' = 2i + 1.
	     V_{i'} = V_{2i + 1} = V_{i+1 + i} = V_{i+1} * V_{i} - V_1
	     V_{i'+1} = V_{2i + 2} = {V_{i+1}}^2 - 2. */
	  mpres_mul (Vi, Vi, Vi1, modulus);
	  mpres_sub (Vi, Vi, S, modulus);
	  mpres_mul (Vi1, Vi1, Vi1, modulus);
	  mpres_sub_ui (Vi1, Vi1, 2, modulus);
	  i = 2UL*i + 1UL;
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
	  i = 2UL*i;
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
	      mpmod_t modulus, listz_t tmp, const unsigned long tmplen,
	      mpzspv_t dct, const mpzspm_t ntt_context)
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
  
  ASSERT (tmplen >= 4 * deg + 2); /* Make sure newtmplen does not underflow */
#ifdef WANT_ASSERT
  mpz_init (leading);
  mpz_mul (leading, F[deg], F[deg]);
  mpz_mod (leading, leading, modulus->orig_modulus);
#endif

  /* Generate V_1(Q)/2 ... V_{deg}(Q)/2, multiply by f_i to form coefficients 
     of G(x). Square the symmetric G(x) polynomial. */

  outputf (OUTPUT_TRACE, "list_scale_V: Q=%Zd, deg = %lu\n", Q, deg);
  list_output_poly (F, deg + 1, 0, 1, "/* list_scale_V */ F(x) = ", "\n", 
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

  list_output_poly (G, deg + 1, 0, 1, "/* list_scale_V */ G(x) = ", "\n", 
		    OUTPUT_TRACE);

  /* Now square the G polynomial in G[0 .. deg], put result in
     G[0 .. 2*deg] */

  /* Bugfix: ks_multiply() does not like negative coefficients. FIXME */

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

  if (dct != NULL && ntt_context != NULL)
    ntt_sqr_recip (G, G, dct, deg + 1, ntt_context);
  else
    list_sqr_reciprocal (G, G, deg + 1, modulus->orig_modulus, 
                         newtmp + 2*(deg+2), newtmplen - 2*(deg+2));

  list_output_poly (G, 2 * deg + 1, 0, 1, "/* list_scale_V */ G(x)^2 == ", 
		    "\n", OUTPUT_TRACE);

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
  H_U = H - 1; /* H_U[0] is undefined, no need to store it */

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
  
  list_output_poly (H, deg, 0, 1, "/* list_scale_V */ H(x) = ", "\n",
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

  if (dct != NULL && ntt_context != NULL)
    ntt_sqr_recip (H, H, dct, deg, ntt_context);
  else
    list_sqr_reciprocal (H, H, deg, modulus->orig_modulus, 
  		         newtmp, newtmplen);

  /* Now there are the 2*deg-1 coefficients in standard basis of a 
     symmetric polynomial of degree 4*deg - 4 in H[0 ... 2*deg-2] */

  list_output_poly (H, 2*deg - 1, 0, 1, "/* list_scale_V */ H(x)^2 == ", "\n",
		    OUTPUT_TRACE);

  /* Multiply by Q^2-4 */
  ASSERT (newtmplen >= 2);
  mpres_mul (Vt, Q, Q, modulus);
  mpres_sub_ui (Vt, Vt, 4, modulus);
  for (i = 0; i <= 2 * deg - 2; i++)
    mpres_mul_z_to_z (H[i], Vt, H[i], modulus);
  list_output_poly (H, 2 * deg - 1, 0, 1, "/* list_scale_V */ "
		    "H(x)^2*(Q^2-4) == ", "\n", OUTPUT_TRACE);

  /* Multiply by (X - 1/X)^2 = X^2 - 2 + 1/X^2 */

  if (deg == 1)
    {
      /* H(X) has degree 2*deg-2 = 0, so H(X) = h_0
	 H(X) * (X - 1/X)^2 = -2 h_0 + h_0 V_2(Y)  */
      ASSERT (newtmplen > 2UL);
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
      ASSERT (newtmplen > 4UL);
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
      ASSERT (newtmplen > 2 * deg);
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
    list_output_poly (F, deg, 1, 0, "poly_from_sets: F = ", "\n", OUTPUT_TRACE);

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
	    list_output_poly (F + j, deg, 1, 0, "poly_from_sets: ", "\n", 
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
	list_output_poly (F, (i + 1) * deg, 1, 0, "poly_from_sets: ", "\n",
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
   are symmetric. Requires that the first set in "sets" has cardinality 2,
   all sets must be symmetric around 0. The resulting polynomial of degree 
   2*d is F(x) = f_0 + \sum_{1 <= i <= d} f_i (x^i + 1/x^i). The coefficient
   f_i is stored in F[i], which therefore needs d+1 elements. */

static int
poly_from_sets_V (listz_t F, const mpres_t Q, sets_long_t *sets, 
		  listz_t tmp, const unsigned long tmplen, mpmod_t modulus,
		  mpzspv_t dct, const mpzspm_t ntt_context)
{
  unsigned long c, deg, i, nr;
  set_long_t *set = sets->sets;
  mpres_t Qt;
  
  ASSERT_ALWAYS (sets->nr > 0UL);
  ASSERT_ALWAYS (set->card == 2UL); /* Check that the cardinality of 
                                       first set is 2 */
  /* Check that first set is symmetric around 0 */
  ASSERT_ALWAYS (set->elem[0] == -set->elem[1]);

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

  V (Qt, Q, set->elem[0], modulus); /* First set in sets is {-k, k}, 
                                       Qt = V_k(Q) */
  mpres_neg (Qt, Qt, modulus);
  mpres_get_z (F[0], Qt, modulus);
  mpz_set_ui (F[1], 1UL);
  deg = 1UL;
  /* Here, F(x) = (x - r^{k_1})(x - r^{-k_1}) / x = 
                  (x^2 - x (r^{k_1} + r^{-k_1}) + 1) / x =
		  (x + 1/x) - V_{k_1}(r + 1/r) */

  for (nr = sets->nr - 1UL; nr > 0UL; nr--)
    {
      /* Assuming the sets are sorted in order of ascending cardinality, 
         we process them back-to-front so the sets of cardinality 2 are 
         processed last, but skipping the first set which we processed 
         already. */
      
      set = sets_nextset (sets->sets); /* Skip first set */
      for (i = 1UL; i < nr; i++) /* Skip over remaining sets but one */
        set = sets_nextset (set);
        
      /* Process this set. We assume it is either of cardinality 2, or of 
	 odd cardinality */
      c = set->card;
      outputf (OUTPUT_DEVVERBOSE, " %lu", c);

      if (c == 2UL)
	{
	  /* Check it's symmetric */
	  ASSERT_ALWAYS (set->elem[0] == -set->elem[1]);
	  V (Qt, Q, set->elem[0], modulus);
	  list_scale_V (F, F, Qt, deg, modulus, tmp, tmplen, dct, 
	                ntt_context);
	  deg *= 2UL;
	  ASSERT (mpz_cmp_ui (F[deg], 1UL) == 0); /* Check it's monic */
	}
      else
	{
	  ASSERT_ALWAYS (c % 2UL == 1UL);
	  /* Generate the F(Q^{k_i} * X)*F(Q^{-k_i} * X) polynomials.
	     Each is symmetric of degree 2*deg, so each has deg+1 coeffients
	     in standard basis. */
	  for (i = 0UL; i < (c - 1UL) / 2UL; i++)
	    {
              /* Check it's symmetric */
	      ASSERT (set->elem[i] == -set->elem[c - 1L - i]);
	      V (Qt, Q, set->elem[i], modulus);
	      ASSERT (mpz_cmp_ui (F[deg], 1UL) == 0); /* Check it's monic */
	      list_scale_V (F + (2UL * i + 1UL) * (deg + 1UL), F, Qt, deg, 
	                    modulus, tmp, tmplen, dct, ntt_context);
	      ASSERT (mpz_cmp_ui (F[(2UL * i + 1UL) * (deg + 1UL) + 2UL * deg], 
	              1UL) == 0); /* Check it's monic */
	    }
	  /* Multiply the polynomials together */
	  for (i = 0UL; i < (c - 1UL) / 2UL; i++)
	    {
	      /* So far, we have the product 
		 F(X) * F(Q^{k_j} * X) * F(Q^{-k_j} * X), 1 <= j <= i,
		 at F. This product has degree 2 * deg + i * 4 * deg, that is
		 (2 * i + 1) * 2 * deg, which means (2 * i + 1) * deg + 1
		 coefficients in F[0 ... (i * 2 + 1) * deg]. */
	      ASSERT (mpz_cmp_ui (F[(2UL * i + 1UL) * deg], 1UL) == 0);
	      ASSERT (mpz_cmp_ui (F[(2UL * i + 1UL) * (deg + 1UL) + 2UL*deg], 
	                          1UL) == 0);
	      list_output_poly (F, (2UL * i + 1UL) * deg + 1, 0, 1, 
				"poly_from_sets_V: Multiplying ", "\n",
				OUTPUT_TRACE);
	      list_output_poly (F + (2UL * i + 1UL) * (deg + 1UL), 
	                        2UL * deg + 1UL, 0, 1, " and ", "\n", 
	                        OUTPUT_TRACE);
	      list_mul_reciprocal (F, 
		                   F, (2UL * i + 1UL) * deg + 1UL, 
			 	   F + (2UL * i + 1UL) * (deg + 1UL), 
				   2UL * deg + 1UL, modulus->orig_modulus,
				   tmp, tmplen);
	      list_mod (F, F, (2UL * i + 3UL) * deg + 1UL, 
	                modulus->orig_modulus);
	      list_output_poly (F, (2UL * i + 3UL) * deg + 1UL, 0, 1, 
                                " = ", "\n", OUTPUT_TRACE);
	      ASSERT (mpz_cmp_ui (F[(2UL * i + 3UL) * deg], 1UL) == 0);
	    }
	  deg *= c;
	}
    }

  mpres_clear (Qt, modulus);
  outputf (OUTPUT_DEVVERBOSE, ")");

  return deg;
}


/* Compute g_i = x_0^{M-i} * r^{(M-i)^2} for 0 <= i < l. 
   x_0 = b_1^{2*k_2 + (2*m_1 + 1) * P}. r = b_1^P. 
   Stores the result in g[0 ... l] and/or in g_ntt[offset ... offset + l] */

static void
pm1_sequence_g (listz_t g_mpz, mpzspv_t g_ntt, const mpres_t b_1, 
                const unsigned long P, const unsigned long M_param, 
		const unsigned long l_param, const mpz_t m_1, const long k_2, 
		mpmod_t modulus_param, const mpzspm_t ntt_context)
{
  mpres_t r[3], x_0, x_Mi;
  mpz_t t;
  unsigned long i;
  long timestart, timestop;
  unsigned long M = M_param;
  unsigned long l = l_param, offset = 0UL;
  mpmod_t modulus;
  int want_output = 1;

  outputf (OUTPUT_VERBOSE, "Computing g_i");

#ifdef _OPENMP
#pragma omp parallel private(r, x_0, x_Mi, t, i, M, l, offset, modulus, want_output) shared(timestart, timestop)
  {
    /* When multi-threading, we adjust the parameters for each thread */

    const int nr_chunks = omp_get_num_threads();
    const int thread_nr = omp_get_thread_num();
    
    l = (l_param - 1) / nr_chunks + 1;
    offset = thread_nr * l;
    l = MIN(l, l_param - offset);
    ASSERT_ALWAYS (M_param >= offset);
    M = M_param - (unsigned long) offset;
    
    /* Let only the master thread print stuff */
    want_output = (omp_get_thread_num() == 0);

    if (want_output)
      outputf (OUTPUT_VERBOSE, " using %d threads", nr_chunks);
#endif

  /* Make a private copy of the mpmod_t struct */
  mpmod_copy (modulus, modulus_param);

  mpz_init (t);
  mpres_init (r[0], modulus);
  mpres_init (r[1], modulus);
  mpres_init (r[2], modulus);
  mpres_init (x_0, modulus);
  mpres_init (x_Mi, modulus);

  if (want_output)
    {
      outputf (OUTPUT_DEVVERBOSE, "sequence_g (g, b_1, P = %lu, M = %lu, "
	       "l = %lu, m_1 = %Zd, k_2 = %ld, modulus)\n", P, M, l, m_1, k_2);
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
    }

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
  if (want_output)
    outputf (OUTPUT_TRACE, "/* pm1_sequence_g */ 2*%ld + (2*%Zd + 1)*P == "
	     "%Zd /* PARI C */\n", k_2, m_1, t);

  mpres_pow (x_0, b_1, t, modulus);  /* x_0 = b_1^{2*k_2 + (2*m_1 + 1)*P} */
  if (want_output && test_verbose (OUTPUT_TRACE))
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
  outputf (OUTPUT_TRACE, "/* pm1_sequence_g */ g_%lu = %Zd; /* PARI */\n", 
	   offset, t);
  if (g_mpz != NULL)
    mpz_set (g_mpz[offset], t);
  if (g_ntt != NULL)
    mpzspv_from_mpzv (g_ntt, offset, &t, 1UL, ntt_context);

  /* So here we have for i = 0
     r[2] = x_0^(M-i) * r^{(M-i)^2}
     r[1] = x_0^{-1} * r^{2(-M+i)+1}
     r[0] = r^2
     t = r[2]
  */

  for (i = 1; i < l; i++)
    {
      if (g_mpz != NULL)
        {
	  mpres_mul_z_to_z (g_mpz[offset + i], r[1], g_mpz[offset + i - 1], 
			    modulus);
	  outputf (OUTPUT_TRACE, "/* pm1_sequence_g */ g_%lu = %Zd;"
		   " /* PARI */\n", offset + i, g_mpz[offset + i]);
        }
      if (g_ntt != NULL)
      {
	  mpres_mul_z_to_z (t, r[1], t, modulus);
	  if (g_mpz == NULL) /* Only one should be non-NULL... */
	      outputf (OUTPUT_TRACE, "/* pm1_sequence_g */ g_%lu = %Zd;"
		       " /* PARI */\n", offset + i, t);
	  mpzspv_from_mpzv (g_ntt, offset + i, &t, 1UL, ntt_context);
      }
      mpres_mul (r[1], r[1], r[0], modulus);
    }

  mpres_clear (r[0], modulus);
  mpres_clear (r[1], modulus);
  mpres_clear (r[2], modulus);
  mpres_clear (x_0, modulus);
  mpres_clear (x_Mi, modulus);
  mpz_clear (t);
  mpmod_clear (modulus); /* Clear our private copy of modulus */

#ifdef _OPENMP
  }
#endif

  timestop = cputime ();
  outputf (OUTPUT_VERBOSE, " took %lums\n", timestop - timestart);
  
  if (test_verbose (OUTPUT_TRACE))
    {
      for (i = 0; i < l_param; i++)
	{
	  outputf (OUTPUT_TRACE, "/* pm1_sequence_g */ g_%lu == x_0^"
		   "(M - %lu) * r^((M - %lu)^2) /* PARI C */\n", i, i, i);
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
pm1_sequence_h (listz_t h, mpzspv_t h_ntt, mpz_t *f, const mpres_t r, 
		const unsigned long d, mpmod_t modulus, 
		const mpzspm_t ntt_context)
{
  mpres_t invr;  /* r^{-1} */
  mpres_t fd[3]; /* finite differences table for r^{-i^2}*/
  mpz_t t;       /* the h_j value as an mpz_t */
  unsigned long j;
  long timestart, timestop;

  outputf (OUTPUT_VERBOSE, "Computing h");
  timestart = cputime ();

  mpres_init (invr, modulus);
  mpres_init (fd[0], modulus);
  mpres_init (fd[1], modulus);
  mpres_init (fd[2], modulus);
  mpz_init (t);

  if (test_verbose (OUTPUT_TRACE))
    {
      mpres_get_z (t, r, modulus);
      outputf (OUTPUT_TRACE, "\n/* pm1_sequence_h */ N = %Zd; "
	       "r = Mod(%Zd, N); /* PARI */\n", 
	       modulus->orig_modulus, t);
    }

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
      mpz_mod (t, f[0], modulus->orig_modulus);
      outputf (OUTPUT_TRACE, "/* pm1_sequence_h */ h_0 = %Zd; /* PARI */\n", t);
      if (h != NULL)
	mpz_set (h[0], t);
      if (h_ntt != NULL)
	mpzspv_from_mpzv (h_ntt, 0UL, &t, 1UL, ntt_context);
    }

  /* Do j = 1 */
  if (d > 1)
    {
      mpres_set (fd[2], fd[1], modulus);        /* fd[2] = r^{-1} */
      mpres_mul (fd[1], fd[1], fd[0], modulus); /* fd[1] = r^{-3} */
      mpres_mul_z_to_z (t, fd[2], f[1], modulus);
      outputf (OUTPUT_TRACE,
	       "/* pm1_sequence_h */ h_1 = %Zd; /* PARI */\n", t);
      if (h != NULL)
	  mpz_set (h[1], t);
      if (h_ntt != NULL)
	  mpzspv_from_mpzv (h_ntt, 1UL, &t, 1UL, ntt_context);
    }
  
  /* Do the rest */
  for (j = 2; j < d; j++)
    {
      mpres_mul (fd[2], fd[2], fd[1], modulus); /* fd[2] = r^{-j^2} */
      mpres_mul (fd[1], fd[1], fd[0], modulus); /* fd[1] = r^{-2*j-1} */
      mpres_mul_z_to_z (t, fd[2], f[j], modulus);
      outputf (OUTPUT_TRACE, "/* pm1_sequence_h */ h_%lu = %Zd; /* PARI */\n", 
	       j, t);

      if (h != NULL)
	mpz_set (h[j], t);
      if (h_ntt != NULL)
	mpzspv_from_mpzv (h_ntt, j, &t, 1UL, ntt_context);
    }
  
  timestop = cputime ();
  outputf (OUTPUT_VERBOSE, " took %lums\n", timestop - timestart);
  
  mpres_clear (fd[2], modulus);
  mpres_clear (fd[1], modulus);
  mpres_clear (fd[0], modulus);
  mpres_clear (invr, modulus);
  mpz_clear (t);

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


#if 0
/* Build polynomial F(x) with roots X^i for i covering all the residue classes
   coprime to beta. F must have space for eulerphi(beta) coefficients.
   method can be 0, 1 or 2, which mean: 0 = old way of computing all the 
   roots and doing a product tree, 1 = using recursive expansion of polynomial
   *without* Chebychev polynomials to utilize symmetry, 2 = using recursive 
   expansion of polynomial *with* Chebychev polynomials */

ATTRIBUTE_UNUSED static int 
pm1_build_poly_F (mpz_t *F, const mpres_t X, mpmod_t modulus, 
		  const unsigned long beta, const int method, 
		  const mpz_t i0, const unsigned long nr, 
		  const unsigned long blocks, const unsigned long tmplen, 
		  mpz_t *tmp)
{
  long timestart, timestop;
  sets_long_t *sets = NULL;
  unsigned long setsize = 0UL, i;
  mpz_t mt;
  const unsigned long dF = eulerphi (beta);
  
  ASSERT (0 <= method && method <= 2);
  
  mpz_init (mt);
  
  if (method == 1 || method == 2)
    {
      sets = get_factored_sorted_sets (&setsize, beta);
      if (sets == NULL)
	  return ECM_ERROR;
      mpz_mul_ui (mt, i0, beta);
      mpz_add_si (mt, mt, set_of_sums_minmax (sets, setsize, 1));
      outputf (OUTPUT_VERBOSE, "Effective B2min = %Zd\n", mt);
      mpz_add_ui (mt, i0, blocks * nr);
      mpz_mul_ui (mt, mt, beta);
      mpz_sub_si (mt, mt, set_of_sums_minmax (sets, setsize, -1));
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
      outputf (OUTPUT_VERBOSE, " took %lums\n", timestop - timestart);

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
	      list_output_poly (F, dF, 1, 0, "F(x) = ", "\n", OUTPUT_ERROR);
	      outputf (OUTPUT_ERROR, "Factored sets: ");
	      sets_print (OUTPUT_ERROR, sets, setsize);
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
      ASSERT_ALWAYS(2 * i == dF);
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
      outputf (OUTPUT_VERBOSE, " took %lums\n", timestop - timestart);
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
	  if (gcd (i, beta) == 1UL)
	    mpres_get_z (F[j++], Xi, modulus);
	  mpres_mul (Xi, Xi, Xs, modulus);
	}
      
      ASSERT_ALWAYS (j == dF);
      mpres_clear (Xs, modulus);
      mpres_clear (Xi, modulus);
      timestop = cputime ();
      outputf (OUTPUT_VERBOSE, " took %lums\n", timestop - timestart);
      
      /* Print the polynomial in linear factors form */
      outputf (OUTPUT_TRACE, "F(x) = ");
      for (j = 0; j < dF - 1; j++)
	outputf (OUTPUT_TRACE, "(x - %Zd) * ", F[j]);
      outputf (OUTPUT_TRACE, "(x - %Zd); /* PARI */\n", F[dF - 1]);
      
      /* Multiply all the (x - f_i) to form F(x) in monomial basis */
      outputf (OUTPUT_VERBOSE, "Building F from its roots");
      timestart = cputime ();
      PolyFromRoots (F, F, dF, tmp, modulus->orig_modulus);
      timestop = cputime ();
      outputf (OUTPUT_VERBOSE, " took %lums\n", timestop - timestart);
    }

  mpz_clear (mt);
  if (method == 1 || method == 2)
    free (sets);

  /* Print the final polynomial in monomial form */
  if (test_verbose (OUTPUT_TRACE))
    list_output_poly (F, dF, 1, 0, "F(x) == ", "\n", OUTPUT_TRACE);

  return 0;
}
#endif

static int 
make_S_1_S_2 (sets_long_t **S_1, set_long_t **S_2, 
              const faststage2_param_t *params)
{
  unsigned long i;
  sets_long_t *facS_2;
  size_t facS_2_size;

  *S_1 = sets_get_factored_sorted (params->P);
  if (*S_1 == NULL)
    return ECM_ERROR;
  ASSERT (sets_sumset_minmax (*S_1, 1) == (long) maxS (params->P));

  *S_2 = malloc (set_sizeof(params->s_2));
  if (*S_2 == NULL)
    {
      free (*S_1);
      return ECM_ERROR;
    }

  /* Extract sets for S_2 and compute the set of sums */
  
  sets_extract (NULL, &facS_2_size, *S_1, params->s_2);
  facS_2 = malloc (facS_2_size);
  if (facS_2 == NULL)
    {
      free (*S_1);
      free (*S_2);
      return ECM_ERROR;
    }
  sets_extract (facS_2, NULL, *S_1, params->s_2);
  sets_sumset (*S_2, facS_2);
  ASSERT_ALWAYS ((*S_2)->card == params->s_2);
  free (facS_2);
  quicksort_long ((*S_2)->elem, (*S_2)->card);
  
  /* Print the sets in devverbose mode */
  if (test_verbose (OUTPUT_DEVVERBOSE))
    {
      outputf (OUTPUT_DEVVERBOSE, "S_1 = ");
      sets_print (OUTPUT_DEVVERBOSE, *S_1);
      
      outputf (OUTPUT_DEVVERBOSE, "S_2 = {");
      for (i = 0UL; i + 1UL < params->s_2; i++)
	outputf (OUTPUT_DEVVERBOSE, "%ld, ", (*S_2)->elem[i]);
      if (i < params->s_2)
	outputf (OUTPUT_DEVVERBOSE, "%ld", (*S_2)->elem[i]); 
      outputf (OUTPUT_DEVVERBOSE, "}\n");
    }

  return 0;
}



static mpzspv_t
mpzspv_init_mt (spv_size_t len, mpzspm_t mpzspm)
{
  int i; /* OpenMP wants the iteration variable a signed type */
  mpzspv_t x = (mpzspv_t) malloc (mpzspm->sp_num * sizeof (spv_t));
  
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
      x[i] = (spv_t) malloc (len * sizeof (sp_t));
	
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
	  free(x[i]);
      return NULL;
    }
  else
    return x;
}


/* Copy the len/2 + 1 distinct coefficients of the DFT of a symmetric signal,
   kinda as if a DCT had been used. The coefficients in the DFT are 
   assumed to be scrambled. The output is not arranged like a regular 
   scrambled DCT would be. */

ATTRIBUTE_UNUSED
static void
ntt_dft_to_dct (mpzspv_t dct, mpzspv_t dft, unsigned long len, 
		const mpzspm_t ntt_context)
{
  unsigned long i, j;
#ifdef WANT_ASSERT
  unsigned long m;
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
  m = 5UL;
  for (i = 2UL; i < len; i += 2UL)
    {
      /* This works, but why? */
      if (i + i / 2UL > m)
	  m = 2UL * m + 1UL;

      for (j = 0UL; j < ntt_context->sp_num; j++)
	{
	  ASSERT (dft[j][i] == dft[j][m - i]);
	}
      outputf (OUTPUT_TRACE, "ntt_dft_to_dct: DFT[%lu] == DFT[%lu]\n",
	       i, m - i);
    }
#endif
  /* Copy coefficients to dct buffer */
  for (j = 0UL; j < ntt_context->sp_num; j++)
    {
      for (i = 0UL; i < len / 2UL; i++)
	  dct[j][i] = dft[j][i * 2UL];
      dct[j][len / 2UL] = dft[j][1];
    }
}


/* Computes a DCT-I of the length dctlen. Input is the spvlen coefficients
   in spv. tmp is temp space and must have space for 2*dctlen-2 sp_t's */

static void
ntt_spv_to_dct (mpzspv_t dct, const mpzspv_t spv, const spv_size_t spvlen, 
                const spv_size_t dctlen, mpzspv_t tmp, 
		const mpzspm_t ntt_context)
{
  const spv_size_t l = 2 * (dctlen - 1); /* Length for the DFT */
  int j;

#ifdef _OPENMP
#pragma omp parallel private(j)
  {
#pragma omp master 
    outputf (OUTPUT_VERBOSE, " using %d threads", omp_get_num_threads());
#pragma omp for
#endif
  for (j = 0; j < (int) ntt_context->sp_num; j++)
    {
      const spm_t spm = ntt_context->spm[j];
      const sp_t root = sp_pow (spm->prim_root, ntt_context->max_ntt_size / l,
                                spm->sp, spm->mul_c);
      spv_size_t i;
      
      /* Make a symmetric copy of spv in tmp. I.e. with spv = [3, 2, 1], 
         spvlen = 3, dctlen = 5 (hence l = 8), we want 
         tmp = [3, 2, 1, 0, 0, 0, 1, 2] */
      spv_set (tmp[j], spv[j], spvlen);
      spv_rev (tmp[j] + l - spvlen + 1, spv[j] + 1, spvlen - 1);
      /* Now we have [3, 2, 1, ?, ?, ?, 1, 2]. Fill the ?'s with zeros. */
      spv_set_sp (tmp[j] + spvlen, (sp_t) 0, l - 2 * spvlen + 1);

#if 0
      printf ("ntt_spv_to_dct: tmp[%d] = [", j);
      for (i = 0; i < l; i++)
          printf ("%lu, ", tmp[j][i]);
      printf ("]\n");
#endif
      
      spv_ntt_gfp_dif (tmp[j], l, spm->sp, spm->mul_c, root);

#if 0
      printf ("ntt_spv_to_dct: tmp[%d] = [", j);
      for (i = 0; i < l; i++)
          printf ("%lu, ", tmp[j][i]);
      printf ("]\n");
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
      {
        spv_size_t m = 5;
        for (i = 2; i < l; i += 2L)
          {
            /* This works, but why? */
            if (i + i / 2L > m)
                m = 2L * m + 1L;

            ASSERT (tmp[j][i] == tmp[j][m - i]);
#if 0
            outputf (OUTPUT_TRACE, "ntt_spv_to_dct: DFT[%lu] == DFT[%lu]\n",
                     i, m - i);
#endif
          }
      }
#endif

      /* Copy coefficients to dct buffer */
      for (i = 0; i < l / 2; i++)
        dct[j][i] = tmp[j][i * 2];
      dct[j][l / 2] = tmp[j][1];
    }
#ifdef _OPENMP
  }
#endif
}


/* Multiply the NTT coefficients in dft, assumed to be a scrambled DFT,
   by the coefficients in dct, assumed to be in for layout produced by 
   ntt_dft_to_dct(); */

static void
ntt_dft_mul_dct (mpzspv_t dft, const mpzspv_t dct, const unsigned long len, 
		 const mpzspm_t ntt_context)
{
  unsigned long i, j, m;
  
  for (j = 0; j < ntt_context->sp_num; j++)
  {
      const sp_t sp = ntt_context->spm[j]->sp; 
      const sp_t mul_c = ntt_context->spm[j]->mul_c;
      m = 5UL;
      
      dft[j][0] = sp_mul (dft[j][0], dct[j][0], sp, mul_c);
      
      for (i = 2UL; i < len; i += 2UL)
      {
	  /* This works, but why? */
	  if (i + i / 2UL > m)
	      m = 2UL * m + 1UL;
	  
	  dft[j][i] = sp_mul (dft[j][i], dct[j][i / 2UL], sp, mul_c);
	  dft[j][m - i] = sp_mul (dft[j][m - i], dct[j][i / 2UL], sp, mul_c);
      }
      dft[j][1] = sp_mul (dft[j][1], dct[j][len / 2UL], sp, mul_c);
  }
}


/* Multiply the polynomial in "dft" by the RLP in "dct", where "dft" 
   contains the polynomial coefficients (not FFT'd yet) and "dct" 
   contains the DCT-I coefficients of the RLP. The latter are 
   assumed to be in the layout produced by ntt_dft_to_dct().
   Output are the coefficients of the product polynomial, stored in dft. */

static void
ntt_mul_by_dct (mpzspv_t dft, const mpzspv_t dct, const unsigned long len, 
		const mpzspm_t ntt_context)
{
  int j;
  
#ifdef _OPENMP
#pragma omp parallel private(j)
  {
#pragma omp master 
    {
      outputf (OUTPUT_VERBOSE, " using %d threads", omp_get_num_threads());
    }
#pragma omp for
#endif
    for (j = 0; j < (int) (ntt_context->sp_num); j++)
      {
	const spm_t spm = ntt_context->spm[j];
	const spv_t spv = dft[j];
	sp_t root;
	unsigned long i, m;
	
	root = sp_pow (spm->prim_root, ntt_context->max_ntt_size / len,
		       spm->sp, spm->mul_c);
	
	/* Forward DFT of dft[j] */
	spv_ntt_gfp_dif (spv, len, spm->sp, spm->mul_c, root);
	
	m = 5UL;
	
	/* Point-wise product */
	spv[0] = sp_mul (spv[0], dct[j][0], spm->sp, spm->mul_c);
	spv[1] = sp_mul (spv[1], dct[j][len / 2UL], spm->sp, spm->mul_c);
	
	for (i = 2UL; i < len; i += 2UL)
	  {
	    /* This works, but why? */
	    if (i + i / 2UL > m)
	      m = 2UL * m + 1;
	    
	    spv[i] = sp_mul (spv[i], dct[j][i / 2UL], spm->sp, spm->mul_c);
	    spv[m - i] = sp_mul (spv[m - i], dct[j][i / 2UL], spm->sp, 
				 spm->mul_c);
	  }
	
	root = sp_pow (spm->inv_prim_root, ntt_context->max_ntt_size / len,
		       spm->sp, spm->mul_c);
	
	/* Inverse transform of dft[j] */
	spv_ntt_gfp_dit (spv, len, spm->sp, spm->mul_c, root);
	
	/* Divide by transform length. FIXME: scale the DCT of h instead */
	spv_mul_sp (spv, spv, spm->sp - (spm->sp - 1) / len, len, spm->sp, 
		    spm->mul_c);
      }
#ifdef _OPENMP
  }
#endif
}



/* Square the reciprocal Laurent polynomial S(x) of degree 2*n-2.
   S(x) = s_0 + \sum_{i=1}^{n-1} s_i (x^i + x^{-1}).
   S[i] contains the n coefficients s_i, 0 <= i <= n-1.
   R[i] will contain the 2n-1 coefficients r_i, 0 <= i <= 2*n-2, where 
   R(x) = S(x)^2 = r_0 + \sum_{i=1}^{2n-2} r_i (x^i + x^{-1}).
   dft must have power of 2 length len >= 2n.
   The NTT primes must be == 1 (mod 4*len).
*/
#undef TRACE_ntt_sqr_recip
void
ntt_sqr_recip (mpzv_t R, const mpzv_t S, mpzspv_t dft, 
               const spv_size_t n, const mpzspm_t ntt_context)
{
  const spv_size_t len = ((spv_size_t) 2) << ceil_log2 (n);
  const spv_size_t d_h = 4*n - 4;
  int j;
  
  outputf (OUTPUT_DEVVERBOSE, 
           "ntt_sqr_recip: called with n = %lu. len = %lu, d_h = %lu\n",
           n, len, d_h);
  
  if (n == 1)
    {
      mpz_mul (R[0], S[0], S[0]);
      return;
    }

  /* We choose an 4*len-th primitive root of unity w.
     Let a = 1/w^len, so that a^2 = -1.
     We want the convolution product (mod x^len - 1/a) = (mod x^len - w^len)
     So we weight input vectors with 1, w, w^2, ... so that wrap around
     leaves excess of w^len = 1/a. */

  ASSERT_ALWAYS (len > 0);
#ifdef TRACE_ntt_sqr_recip
  printf ("ntt_sqr_recip: n %lu, length %lu, d_h %lu\n", n, len, d_h);
  gmp_printf ("Input polynomial is %Zd", S[0]);
  for (j = 1; (spv_size_t) j < n; j++)
    gmp_printf (" + %Zd * (x^%lu + x^(-%lu))", S[j], j, j);
  printf ("\n");
#endif
  
  ASSERT_ALWAYS (ntt_context->max_ntt_size % (4 * len) == 0UL);
  /* Fill NTT elements [n-1 .. 2n-2] with coefficients */
  mpzspv_from_mpzv (dft, n - 1, S, (spv_size_t) n, ntt_context);
  /* Fill NTT elements [0 .. n-2] with coefficients */
  mpzspv_revcopy (dft, 0, dft, n, n - 1, ntt_context);
  /* Zero out NTT elements [2n-1 .. len-1] */
  mpzspv_set_sp (dft, 2*n - 1, (sp_t) 0, len - 2*n + 1, ntt_context);

  for (j = 0UL; j < (int) (ntt_context->sp_num); j++)
    {
      const spm_t spm = ntt_context->spm[j];
      const spv_t spv = dft[j];
      sp_t root, iw, w, weight;
      const sp_t sp = spm->sp, mul_c = spm->mul_c;
      spv_size_t i;

#ifdef TRACE_ntt_sqr_recip
      if (j == 0UL)
        {
	  printf ("ntt_sqr_recip: NTT vector mod %lu\n", sp);
	  printf ("ntt_sqr_recip: before weighting: [%lu", spv[0]);
	  for (i = 1; i < len; i++)
	    printf (", %lu", spv[i]);
	  printf ("]\n");
        }
#endif

      /* Compute the root for the weight signal, a 4*len-th primitive root 
	 of unity */
      w = sp_pow (spm->prim_root, ntt_context->max_ntt_size / (4 * len), 
                  sp, mul_c);
      /* Compute the root for the un-weight signal, the reciprocal of
         the previous w */
      iw = sp_pow (spm->inv_prim_root, ntt_context->max_ntt_size / (4 * len), 
                   sp, mul_c);
#ifdef TRACE_ntt_sqr_recip
      if (j == 0)
        printf ("w = %lu ,iw = %lu\n", w, iw);
#endif
      ASSERT(sp_mul(w, iw, sp, mul_c) == (sp_t) 1);

      /* Apply weight signal */
      weight = w;
      for (i = 1; i < len; i++)
        {
          spv[i] = sp_mul (spv[i], weight, sp, mul_c);
          weight = sp_mul (weight, w, sp, mul_c);
        }
#ifdef TRACE_ntt_sqr_recip
      if (j == 0UL)
        {
	  printf ("ntt_sqr_recip: after weighting: [%lu", spv[0]);
	  for (i = 1; i < len; i++)
	    printf (", %lu", spv[i]);
	  printf ("]\n");
        }
#endif

      /* Compute root for the transform */
      root = sp_pow (w, (sp_t) 4, sp, mul_c);

      /* Forward DFT of dft[j] */
      spv_ntt_gfp_dif (spv, len, sp, mul_c, root);

#if 0 && defined(TRACE_ntt_sqr_recip)
      if (j == 0UL)
        {
	  printf ("ntt_sqr_recip: after forward transform: [%lu", spv[0]);
	  for (i = 1; i < len; i++)
	    printf (", %lu", spv[i]);
	  printf ("]\n");
        }
#endif

      /* Square the transformed vector point-wise */
      spv_pwmul (spv, spv, spv, len, sp, mul_c);
      
#if 0 && defined(TRACE_ntt_sqr_recip)
      if (j == 0UL)
        {
	  printf ("ntt_sqr_recip: after point-wise squaring: [%lu", spv[0]);
	  for (i = 1; i < len; i++)
	    printf (", %lu", spv[i]);
	  printf ("]\n");
        }
#endif

      /* Compute root for inverse transform */
      root = sp_pow (iw, (sp_t) 4, sp, mul_c);

      /* Inverse transform of dft[j] */
      spv_ntt_gfp_dit (spv, len, sp, mul_c, root);
      
#if 0 && defined(TRACE_ntt_sqr_recip)
      if (j == 0UL)
        {
	  printf ("ntt_sqr_recip: after inverse transform: [%lu", spv[0]);
	  for (i = 1; i < len; i++)
	    printf (", %lu", spv[i]);
	  printf ("]\n");
        }
#endif

      /* Divide by transform length */
      weight = sp - (sp - (sp_t) 1) / len; /* weight = 1/len (mod sp) */
      for (i = 0; i < len; i++)
	spv[i] = sp_mul (spv[i], weight, sp, mul_c);
#ifdef TRACE_ntt_sqr_recip
      if (j == 0UL)
        {
	  printf ("ntt_sqr_recip: after dividing by transform length: [%lu", 
		  spv[0]);
	  for (i = 1; i < len; i++)
	    printf (", %lu", spv[i]);
	  printf ("]\n");
        }
#endif

      /* Un-weight */
      weight = iw;
      for (i = 1; i < len; i++)
        {
          spv[i] = sp_mul (spv[i], weight, sp, mul_c);
          weight = sp_mul (weight, iw, sp, mul_c);
        }
#ifdef TRACE_ntt_sqr_recip
      if (j == 0UL)
        {
      printf ("ntt_sqr_recip: after un-weighting: [%lu", spv[0]);
      for (i = 1; i < len; i++)
        printf (", %lu", spv[i]);
      printf ("]\n");
        }
#endif

      /* Separate the coefficients of R in the wrapped-around product. */
      /* The product h(x) = \sum_{i=0}^{d_h} h_i x^i % (x^len - a)
	 has h_i + a*h_{i + len} at position 0 <= i <= d_h - len
	 and h_i at position d_h - len < i < len.
	 The coefficients h_i and h_{d_h - i} are identical.
	 
	 Hence for 0 <= i <= d_h - len we have 
	 at position i:             h_i + a*h_{d_h - len - i}
	 at position d_h - len - i: a*h_i + h_{d_h - len - i}
	 This implies at position i = (d_h - len) / 2: (A+1) * h_i
      */

      /* Put a = w^len = sqrt(-1) in weight */
      weight = sp_pow (w, len, sp, mul_c);
      
      for (i = 0; i <= (d_h - len) / 2; i++)
        {
          sp_t t, u;
	  const spv_size_t i2 = d_h - len - i;
          t = sp_mul (spv[i], weight, sp, mul_c); 
          /* t = a*(h_i + a*h_{d_h-len-i}) */
          t = sp_sub (t, spv[i2], sp);
          /* t = a*(h_i + a*h_{d_h-len-i}) - (a*h_i + h_{d_h-len-i}) = 
                 (a^2-1) h_{d_h-len-i} = -2*h_{d_h-len-i} */
          if (t & (sp_t) 1)
            t = (sp - t) >> 1;
          else
            {
              if ((t >>= 1) != (sp_t) 0)
                t = sp - t;
            }
          /* t = h_{d_h-len-i} */
          u = t;
          t = sp_mul (t, weight, sp, mul_c);
          /* t = a*h_{d_h-len-i} */
          t = sp_sub (spv[i], t, sp);
	  /* t = h_i + a*h_{d_h-len- i} - a*h_{d_h-len-i} = h_i */
          spv[i] = t;
          spv[i2] = u;
        }

#ifdef TRACE_ntt_sqr_recip
      if (j == 0UL)
        {
      printf ("ntt_sqr_recip: after un-wrapping: [%lu", spv[0]);
      for (i = 1; i < len; i++)
        printf (", %lu", spv[i]);
      printf ("]\n");
        }
#endif

      /* Revert the entries in spv[0 ... d_h / 2 + 1] */
      for (i = 0; i < d_h / 4; i++)
	{
	  sp_t t, u;
	  t = spv[i];
	  u = spv[d_h / 2 - i];
	  spv[i] = u;
	  spv[d_h / 2 - i] = t;
	}
#ifdef TRACE_ntt_sqr_recip
      if (j == 0UL)
        {
      printf ("ntt_sqr_recip: after mirroring: [%lu", spv[0]);
      for (i = 1; i < len; i++)
        printf (", %lu", spv[i]);
      printf ("]\n");
        }
#endif
    }

  mpzspv_to_mpzv (dft, (spv_size_t) 0, R, d_h / 2 + 1, ntt_context);
  for (j = 0; (spv_size_t) j < d_h / 2 + 1; j++)
    mpz_mod (R[j], R[j], ntt_context->modulus);
#ifdef TRACE_ntt_sqr_recip
  gmp_printf ("ntt_sqr_recip: Output polynomial is %Zd", R[0]);
  for (j = 1; (spv_size_t) j < d_h / 2 + 1; j++)
    gmp_printf (" + %Zd * (x^%lu + x^(-%lu))", R[j], j, j);
  printf ("\n");
#endif
}


/* Computes gcd(\prod_{0 <= i < len} (ntt[i + offset] + add[i]), N), 
   the NTT residues are converted to integer residues (mod N) first.
   If add == NULL, add[i] is assumed to be 0. */

static void
ntt_gcd (mpz_t f, mpzspv_t ntt, const unsigned long ntt_offset,
	 const listz_t add, const unsigned long len_param, 
	 const mpzspm_t ntt_context, mpmod_t modulus_param)
{
  unsigned long i, j;
  const unsigned long Rlen = MPZSPV_NORMALISE_STRIDE;
  listz_t R;
  unsigned long len = len_param, thread_offset = 0;
  mpres_t tmpres, tmpprod, totalprod;
  mpmod_t modulus;
  long timestart, timestop;
  
  outputf (OUTPUT_VERBOSE, "Computing gcd of coefficients and N");
  timestart = cputime ();

  /* All the threads will multiply their partial products to this one. */
  mpres_init (totalprod, modulus_param);
  mpres_set_ui (totalprod, 1UL, modulus_param);

#ifdef _OPENMP
#pragma omp parallel private(i, j, R, len, thread_offset, tmpres, tmpprod, modulus) shared(totalprod)
  {
    const int nr_chunks = omp_get_num_threads();
    const int thread_nr = omp_get_thread_num();
    len = (len_param - 1) / nr_chunks + 1;
    thread_offset = thread_nr * len;
    len = MIN(len, len_param - thread_offset);
#pragma omp master
    {
      outputf (OUTPUT_VERBOSE, " using %d threads", nr_chunks);
    }
#endif

    /* Make a private copy of the mpmod_t struct */
    mpmod_copy (modulus, modulus_param);

    R = init_list (Rlen);
    mpres_init (tmpres, modulus);
    mpres_init (tmpprod, modulus);
    mpres_set_ui (tmpprod, 1UL, modulus);
    
    for (i = 0; i < len; i += Rlen)
      {
	const unsigned long blocklen = MIN(len - i, Rlen);

	/* Convert blocklen residues from NTT to integer representatives
	   and store them in R */
	mpzspv_to_mpzv (ntt, ntt_offset + thread_offset + i, R, blocklen, 
			ntt_context);

	/* Accumulate product in tmpprod */
	for (j = 0; j < blocklen; j++)
	  {
	    outputf (OUTPUT_TRACE, "r_%lu = %Zd; /* PARI */\n", i, R[j]);
	    if (add != NULL)
	      mpz_add (R[j], R[j], add[i + thread_offset + j]);
	    mpres_set_z_for_gcd (tmpres, R[j], modulus);
#define TEST_ZERO_RESULT
#ifdef TEST_ZERO_RESULT
	    if (mpres_is_zero (tmpres, modulus))
	      outputf (OUTPUT_VERBOSE, "R_[%lu] = 0\n", i);
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
    mpres_clear (tmpres, modulus);
    mpres_clear (tmpprod, modulus);
    mpmod_clear (modulus);
    clear_list (R, Rlen);
#ifdef _OPENMP
  }
#endif

  if (test_verbose(OUTPUT_RESVERBOSE))
    {
      mpres_t t;
      mpres_init (t, modulus_param);
      mpres_get_z (t, totalprod, modulus_param);
      outputf (OUTPUT_RESVERBOSE, "\nProduct of R[i] = %Zd (times some "
	       "power of 2 if REDC was used! Try -mpzmod)\n", t);
      mpres_clear (t, modulus_param);
    }

  mpres_gcd (f, totalprod, modulus_param);
  mpres_clear (totalprod, modulus_param);

  timestop = cputime ();
  outputf (OUTPUT_VERBOSE, " took %lums\n", timestop - timestart);
}



int 
pm1fs2 (mpz_t f, const mpres_t X, mpmod_t modulus, 
	const faststage2_param_t *params)
{
  unsigned long phiP, nr;
  unsigned long i, l, lenF, lenG, lenR, tmplen;
  sets_long_t *S_1; /* This is stored as a set of sets (arithmetic 
                       progressions of prime length */
  set_long_t *S_2; /* This is stored as a regular set */
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
  long timetotalstart, timestart, timestop;

  timetotalstart = cputime ();

  phiP = eulerphi (params->P);
  ASSERT_ALWAYS (phiP == params->s_1 * params->s_2);
  ASSERT_ALWAYS (params->s_1 < params->l);
  nr = params->l - params->s_1; /* Number of points we evaluate */

  outputf (OUTPUT_TRACE, "compare(a, b, n) = if (a != b,print(\"In PARI \", n, \" line, \", a \" != \" b)); /* PARI */\n");

  if (make_S_1_S_2 (&S_1, &S_2, params) == ECM_ERROR)
      return ECM_ERROR;

  /* Allocate all the memory we'll need */
  /* Allocate the correct amount of space for each mpz_t or the 
     reallocations will up to double the time for stage 2! */
  mpz_init (mt);
  mpres_init (mr, modulus);
  lenF = params->s_1 / 2 + 1 + 1; /* Another +1 because poly_from_sets_V stores
				     the leading 1 monomial for each factor */
  F = init_list2 (lenF, (unsigned int) abs (modulus->bits));
  h = malloc ((params->s_1 + 1) * sizeof (mpz_t));
  lenG = params->l;
  g = init_list2 (lenG, (unsigned int) abs (modulus->bits));
  lenR = nr;
  R = init_list2 (lenR, (unsigned int) abs (modulus->bits));    
  tmplen = 3UL * params->l + list_mul_mem (params->l / 2);
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
  tmp = init_list2 (tmplen, (unsigned int) abs (modulus->bits));
#endif
  
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
  
  i = poly_from_sets_V (F, mr, S_1, tmp, tmplen, modulus, NULL, NULL);
  ASSERT_ALWAYS(2 * i == params->s_1);
  ASSERT(mpz_cmp_ui (F[i], 1UL) == 0);
  free (S_1);
  S_1 = NULL;
  
  timestop = cputime ();
  outputf (OUTPUT_VERBOSE, " took %lums\n", timestop - timestart);
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
  pm1_sequence_h (F, NULL, F, mr, params->s_1 / 2 + 1, modulus, NULL); 

  /* Make a symmetric copy of F in h. It will have length 
     s_1 + 1 = 2*lenF - 1 */
  /* I.e. with F = [3, 2, 1], s_1 = 4, we want h = [1, 2, 3, 2, 1] */
  for (i = 0; i < params->s_1 / 2 + 1; i++)
    *(h[i]) = *(F[params->s_1 / 2 - i]); /* Clone the mpz_t. 
					    Don't tell Torbjorn */
  for (i = 0; i < params->s_1 / 2; i++)
    *(h[i + params->s_1 / 2 + 1]) = *(F[i + 1]);
  if (test_verbose (OUTPUT_TRACE))
    {
      for (i = 0; i < params->s_1 + 1; i++)
        outputf (OUTPUT_VERBOSE, "h_%lu = %Zd; /* PARI */\n", i, h[i]);
      outputf (OUTPUT_VERBOSE, "h(x) = h_0");
      for (i = 1; i < params->s_1 + 1; i++)
        outputf (OUTPUT_VERBOSE, " + h_%lu * x^%lu", i, i);
      outputf (OUTPUT_VERBOSE, " /* PARI */\n");
    }

  for (l = 0; l < params->s_2; l++)
    {
      const unsigned long M = params->l - 1L - params->s_1 / 2L;
      pm1_sequence_g (g, NULL, X, params->P, M, params->l, 
		      params->m_1, S_2->elem[l], modulus, NULL);

      /* Do the convolution */
      /* Use the transposed "Middle Product" algorithm */
      /* TMulGen reverses the first input sequence, but that doesn't matter
	 since h is symmetric. */

      outputf (OUTPUT_VERBOSE, "TMulGen of g and h");
      timestart = cputime ();
      ASSERT(tmplen >= TMulGen_space (nr - 1, params->l - 1, params->s_1));

      /* Computes rev(h)*g, stores coefficients of x^(s_1) to 
	 x^(s_1+nr-1) = x^(len-1) */
      TMulGen (R, nr - 1, h, params->s_1, g, params->l - 1, tmp, 
	       modulus->orig_modulus);
      list_mod (R, R, nr, modulus->orig_modulus);

      timestop = cputime ();
      outputf (OUTPUT_VERBOSE, " took %lums\n", timestop - timestart);

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

      if (test_verbose (OUTPUT_TRACE))
	{
	  for (i = 0; i < nr; i++)
	    outputf (OUTPUT_TRACE, "r_%lu = %Zd; /* PARI */\n", i, R[i]);
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

      timestop = cputime ();
      outputf (OUTPUT_VERBOSE, " took %lums\n", timestop - timestart);
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
  
  free (h);
  clear_list (F, lenF);
  clear_list (g, lenG);
  clear_list (R, lenR);    
  clear_list (tmp, tmplen);

  mpz_clear (mt);
  mpres_clear (mr, modulus);

  timestop = cputime ();
  outputf (OUTPUT_NORMAL, "Step 2 took %ldms\n", 
           timestop - timetotalstart);
  
  return youpi;
}


int 
pm1fs2_ntt (mpz_t f, const mpres_t X, mpmod_t modulus, 
	const faststage2_param_t *params)
{
  unsigned long nr;
  unsigned long i, l, lenF, tmplen;
  sets_long_t *S_1; /* This is stored as a set of sets (arithmetic 
                       progressions of prime length */
  set_long_t *S_2; /* This is stored as a regular set */
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
  mpres_t tmpres; /* All-purpose temp mpres_t */
  int youpi = ECM_NO_FACTOR_FOUND;
  long timetotalstart, timestart, timestop;

  timetotalstart = cputime ();

  ASSERT_ALWAYS (eulerphi (params->P) == params->s_1 * params->s_2);
  ASSERT_ALWAYS (params->s_1 < params->l);
  nr = params->l - params->s_1; /* Number of points we evaluate */

  outputf (OUTPUT_TRACE, "compare(a, b, n) = if (a != b,print(\"In PARI \", n, \" line, \", a \" != \" b)); /* PARI */\n");

  if (make_S_1_S_2 (&S_1, &S_2, params) == ECM_ERROR)
      return ECM_ERROR;

  /* Precompute the small primes, primitive roots and inverses etc. for 
     the NTT. mpzspm_init() chooses the NTT primes large enough for 
     residues up to 4*l*modulus^2, so adding in Fourier space is ok. 
     The code to multiply wants a 4*len-th root of unity, where len is
     the smallest power of 2 > s_1 */
  ntt_context = mpzspm_init (MAX(params->l, 2UL << ceil_log2 (params->s_1)), 
                             modulus->orig_modulus);
  if (test_verbose (OUTPUT_DEVVERBOSE))
    {
      double modbits = 0.;
      outputf (OUTPUT_DEVVERBOSE, "CRT modulus = %lu", ntt_context->spm[0]->sp);
      modbits += log ((double) ntt_context->spm[0]->sp);
      for (i = 1; i < ntt_context->sp_num; i++)
	{
	  outputf (OUTPUT_DEVVERBOSE, " * %lu", ntt_context->spm[i]->sp);
	  modbits += log ((double) ntt_context->spm[i]->sp);
	}
      outputf (OUTPUT_DEVVERBOSE, ", has %f bits\n", modbits / log (2.));
    }

  /* Allocate all the memory we'll need for building f */
  mpz_init (mt);
  mpres_init (tmpres, modulus);
  lenF = params->s_1 / 2 + 1 + 1; /* Another +1 because poly_from_sets_V stores
				     the leading 1 monomial for each factor */
  tmplen = 3 * params->s_1 + 1000;
  F = init_list2 (lenF, (unsigned int) abs (modulus->bits));
  tmp = init_list2 (tmplen, (unsigned int) abs (modulus->bits));
  /* Allocate memory for h_ntt */
  h_ntt = mpzspv_init_mt (params->l / 2 + 1, ntt_context);

  mpres_get_z (mt, X, modulus); /* mpz_t copy of X for printing */
  outputf (OUTPUT_TRACE, 
	   "N = %Zd; X = Mod(%Zd, N); /* PARI */\n", 
	   modulus->orig_modulus, mt);


  /* Compute the polynomial f(x) = \prod_{k_1 in S_1} (x - X^{2k_1}) */
  outputf (OUTPUT_VERBOSE, "Computing F from factored S_1");
  
  timestart = cputime ();
  
  /* First compute X^2 + 1/X^2 */
  mpres_invert (tmpres, X, modulus);
  mpres_add (tmpres, tmpres, X, modulus);
  V (tmpres, tmpres, 2UL, modulus);
  
  i = poly_from_sets_V (F, tmpres, S_1, tmp, tmplen, modulus, h_ntt, 
                        ntt_context);
  ASSERT_ALWAYS(2 * i == params->s_1);
  ASSERT(mpz_cmp_ui (F[i], 1UL) == 0);
  free (S_1);
  S_1 = NULL;
  
  timestop = cputime ();
  outputf (OUTPUT_VERBOSE, " took %lums\n", timestop - timestart);
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

  mpz_set_ui (mt, params->P);
  mpres_pow (tmpres, X, mt, modulus); /* tmpres = X^P */
  pm1_sequence_h (NULL, h_ntt, F, tmpres, params->s_1 / 2 + 1, modulus, 
		  ntt_context);

  clear_list (F, lenF);
  g_ntt = mpzspv_init_mt (params->l, ntt_context);

  /* Compute the DCT-I of h */
  outputf (OUTPUT_VERBOSE, "Computing DCT-I of h");
  timestart = cputime ();
  ntt_spv_to_dct (h_ntt, h_ntt, params->s_1 / 2 + 1, params->l / 2 + 1, 
                  g_ntt, ntt_context);
  timestop = cputime ();
  outputf (OUTPUT_VERBOSE, " took %lums\n", timestop - timestart);
                
  for (l = 0; l < params->s_2; l++)
    {
      const unsigned long M = params->l - 1L - params->s_1 / 2L;

      /* Compute the coefficients of the polynomial g(x) */
      pm1_sequence_g (NULL, g_ntt, X, params->P, M, params->l, 
		      params->m_1, S_2->elem[l], modulus, ntt_context);

      /* Do the convolution */
      outputf (OUTPUT_VERBOSE, "Computing g*h");
      timestart = cputime ();
      ntt_mul_by_dct (g_ntt, h_ntt, params->l, ntt_context);
      timestop = cputime ();
      outputf (OUTPUT_VERBOSE, " took %lums\n", timestop - timestart);      
      
      /* Compute GCD of N and coefficients of product polynomial */
      ntt_gcd (mt, g_ntt, params->s_1 / 2, NULL, nr, ntt_context, modulus);

      /* If we found a factor, stop */
      if (mpz_cmp_ui (mt, 1UL) > 0)
	{
	  mpz_set (f, mt);
	  youpi = ECM_FACTOR_FOUND_STEP2;
	  break;
	}
    }

  mpzspv_clear (g_ntt, ntt_context);
  mpzspv_clear (h_ntt, ntt_context);
  mpzspm_clear (ntt_context);
  mpres_clear (tmpres, modulus);
  mpz_clear (mt);

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
	     mpres_t *tmp)
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
pp1_sequence_g (listz_t g_x, listz_t g_y, mpzspv_t g_x_ntt, mpzspv_t g_y_ntt,
		const mpres_t b1_x, const mpres_t b1_y, const unsigned long P, 
		const mpres_t Delta, const long M_param, 
		const unsigned long l_param, const mpz_t m_1, const long k_2, 
		const mpmod_t modulus_param, const mpzspm_t ntt_context)
{
  const unsigned long tmplen = 3;
  const int want_x = (g_x != NULL || g_x_ntt != NULL);
  const int want_y = (g_y != NULL || g_y_ntt != NULL);
  mpres_t r_x, r_y, x0_x, x0_y, v2,
      r1_x[2], r1_y[2], r2_x[2], r2_y[2], 
      v[2], tmp[tmplen];
  mpz_t mt;
  mpmod_t modulus; /* Thread-local copy of modulus_param */
  unsigned long i, l = l_param, offset = 0;
  long M = M_param;
  long timestart, timestop;
  int want_output = 1;

  outputf (OUTPUT_VERBOSE, "Computing %s%s%s", 
	   (want_x) ? "g_x" : "", 
	   (want_x && want_y) ? " and " : "",
	   (want_y) ? "g_y" : "");
  timestart = cputime ();
  
#ifdef _OPENMP
#pragma omp parallel private(r_x, r_y, x0_x, x0_y, v2, r1_x, r1_y, r2_x, r2_y, v, tmp, mt, modulus, i, l, offset, M, want_output)
  {
    /* When multi-threading, we adjust the parameters for each thread */

    const int nr_chunks = omp_get_num_threads();
    const int thread_nr = omp_get_thread_num();

    l = (l_param - 1) / nr_chunks + 1;
    offset = thread_nr * l;
    l = MIN(l, l_param - offset);
    M = M_param - (long) offset;

    want_output = (omp_get_thread_num() == 0);
    if (want_output)
      outputf (OUTPUT_VERBOSE, " using %d threads", nr_chunks);
#endif
    mpmod_copy (modulus, modulus_param);
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
    
    if (want_output && test_verbose (OUTPUT_TRACE))
      {
	mpres_get_z (mt, Delta, modulus);
	outputf (OUTPUT_TRACE, 
		 "\n/* pp1_sequence_g */ w = quadgen (4*%Zd); P = %lu; "
		 "M = %ld; k_2 = %ld; m_1 = %Zd; N = %Zd; /* PARI */\n", 
		 mt, P, M, k_2, m_1, modulus->orig_modulus);
	
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
    gfp_ext_pow_norm1_ul (r_x, r_y, b1_x, b1_y, P, Delta, modulus, 
			  tmplen, tmp);
    if (want_output && test_verbose (OUTPUT_TRACE))
      {
	outputf (OUTPUT_TRACE, "/* pp1_sequence_g */ r == ");
	gfp_ext_print (r_x, r_y, modulus, OUTPUT_TRACE);
	outputf (OUTPUT_TRACE, " /* PARI C */\n");
      }
    
    /* Compute x0 = x_0 */
    mpz_mul_2exp (mt, m_1, 1UL);
    mpz_add_ui (mt, mt, 1UL);
    mpz_mul_ui (mt, mt, P);
    mpz_add_si (mt, mt, k_2);
    mpz_add_si (mt, mt, k_2); /* mt = 2*k_2 + (2*m_1 + 1) * P */
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
    gfp_ext_pow_norm1_ul (v[0], v[1], r_x, r_y, M, Delta, modulus, 
			  tmplen, tmp); /* v[0,1] = r^M */
    gfp_ext_mul (v[0], v[1], v[0], v[1], x0_x, x0_y, Delta, modulus, 
		 tmplen, tmp); /* v[0,1] = r^M * x_0 */
    gfp_ext_pow_norm1_ul (r1_x[0], r1_y[0], v[0], v[1], M, Delta, modulus, 
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
    gfp_ext_pow_norm1_ul (v[0], v[1], r_x, r_y, M - 1, Delta, modulus, 
			  tmplen, tmp);
    gfp_ext_mul (v[0], v[1], v[0], v[1], x0_x, x0_y, Delta, modulus, 
		 tmplen, tmp);
    gfp_ext_pow_norm1_ul (r1_x[1], r1_y[1], v[0], v[1], M - 1, Delta, 
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
    mpmod_clear (modulus);
#ifdef _OPENMP
  }
#endif
  
  timestop = cputime ();
  outputf (OUTPUT_VERBOSE, " took %lums\n", timestop - timestart);

  if (g_x != NULL && g_y != NULL && test_verbose(OUTPUT_TRACE))
    {
      for (i = 0; i < l_param; i++)
	{
	  outputf (OUTPUT_TRACE, "/* pp1_sequence_g */ g_%lu = "
		   "x_0^(M-%lu) * r^((M-%lu)^2); /* PARI */", i, i, i);
	  outputf (OUTPUT_TRACE, "/* pp1_sequence_g */ g_%lu == "
		   "%Zd + %Zd*w /* PARI C */\n", i, g_x[i], g_y[i]);
	}
    }
}


/* Compute r[i] = b1^(-P*(k+i)^2) * f_i for i = 0, 1, ..., l-1, where "b1" is 
   an element of norm 1 in the quadratic extension ring */

static void
pp1_sequence_h (listz_t h_x, listz_t h_y, mpzspv_t h_x_ntt, mpzspv_t h_y_ntt,
		const listz_t f, const mpres_t b1_x, const mpres_t b1_y, 
		const long k, const unsigned long l, const unsigned long P, 
		const mpres_t Delta, mpmod_t modulus, 
		const mpzspm_t ntt_context, const unsigned long origtmplen, 
		mpres_t *origtmp)
{
  mpres_t *s_x = origtmp, *s_y = origtmp + 3, *s2_x = origtmp + 6, 
    *s2_y = origtmp + 8, *v = origtmp + 10, *V2 = origtmp + 12,
    *rn_x = origtmp + 13, *rn_y = origtmp + 14;
  const unsigned long newtmplen = origtmplen - 15;
  mpres_t *newtmp = origtmp + 15;
  mpz_t mt;
  unsigned long i;
  long timestart, timestop;

  if (l == 0UL)
    return;

  ASSERT (origtmplen >= 15UL);
  ASSERT (f != h_x);
  ASSERT (f != h_y);

  outputf (OUTPUT_VERBOSE, "Computing h_x and h_y");
  timestart = cputime ();

  mpz_init (mt);

  if (test_verbose (OUTPUT_TRACE))
    {
      mpres_get_z (mt, Delta, modulus);
      outputf (OUTPUT_TRACE, "\n/* pp1_sequence_h */ w = quadgen (4*%Zd); "
	       "k = %ld; P = %lu; N = %Zd; /* PARI */\n", 
	       mt, k, P, modulus->orig_modulus);
      outputf (OUTPUT_TRACE, "/* pp1_sequence_h */ b_1 = ");
      gfp_ext_print (b1_x, b1_y, modulus, OUTPUT_TRACE);
      outputf (OUTPUT_TRACE, "; rn = b_1^(-P); /* PARI */\n");
      for (i = 0; i < l; i++)
	outputf (OUTPUT_TRACE, 
		 "/* pp1_sequence_h */ f_%lu = %Zd; /* PARI */\n", i, f[i]);
    }

  ASSERT (origtmplen >= 13UL);

  /* Compute rn = b_1^{-P} */
  gfp_ext_pow_norm1_ul (*rn_x, *rn_y, b1_x, b1_y, P, Delta, modulus, newtmplen,
			newtmp);
  mpres_neg (*rn_y, *rn_y, modulus);

  /* Compute s[0] = rn^(k^2) = r^(-k^2). We do it by two exponentiations by 
     k and use v[0] and v[1] as temp storage */
  gfp_ext_pow_norm1_ul (v[0], v[1], *rn_x, *rn_y, k, Delta, modulus, 
			newtmplen, newtmp);
  gfp_ext_pow_norm1_ul (s_x[0], s_y[0], v[0], v[1], k, Delta, modulus, 
			newtmplen, newtmp);
  if (test_verbose (OUTPUT_TRACE))
    {
      outputf (OUTPUT_TRACE, "/* pp1_sequence_h */ rn^(k^2) == ");
      gfp_ext_print (s_x[0], s_y[0], modulus, OUTPUT_TRACE);
      outputf (OUTPUT_TRACE, " /* PARI C */\n");
    }

  /* Compute s[1] = r^(-(k+1)^2) = r^(-(k^2 + 2k + 1))*/
  if (l > 1)
    {
      /* v[0] + v[1]*sqrt(Delta) still contains rn^k */
      gfp_ext_sqr_norm1 (s_x[1], s_y[1], v[0], v[1], modulus);
      /* Now s[1] = r^(-2k) */
      gfp_ext_mul (s_x[1], s_y[1], s_x[1], s_y[1], s_x[0], s_y[0], Delta, 
		   modulus, newtmplen, newtmp);
      /* Now s[1] = r^(-(k^2 + 2k)) */
      gfp_ext_mul (s_x[1], s_y[1], s_x[1], s_y[1], *rn_x, *rn_y, Delta, 
		   modulus, newtmplen, newtmp);
      /* Now s[1] = r^(-(k^2 + 2k + 1)) = r^(-(k+1)^2) */
      if (test_verbose (OUTPUT_TRACE))
	{
	  outputf (OUTPUT_TRACE, "/* pp1_sequence_h */ rn^((k+1)^2) == ");
	  gfp_ext_print (s_x[1], s_y[1], modulus, OUTPUT_TRACE);
	  outputf (OUTPUT_TRACE, " /* PARI C */\n");
	}
    }

  /* Compute s2[0] = r^(k^2+2) = r^(k^2) * r^2 */
  gfp_ext_sqr_norm1 (v[0], v[1], *rn_x, *rn_y, modulus);
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
  mpres_add (*V2, *rn_x, *rn_x, modulus); /* V2 = r + 1/r  = 2*rn_x */
  V (v[0], *V2, 2 * k + 1, modulus);  /* v[0] = V_{2k+1} (r + 1/r) */
  V (v[1], *V2, 2 * k + 3, modulus);  /* v[0] = V_{2k+3} (r + 1/r) */
  mpres_mul (*V2, *V2, *V2, modulus); /* V2 = 4*a_x^2 */
  mpres_sub_ui (*V2, *V2, 2UL, modulus); /* V2 = 4*a_x^2 - 2 */

  for (i = 0; i < 2 && i < l; i++)
  {
      mpres_mul (s_y[i], s_y[i], Delta, modulus);
      mpres_mul (s2_y[i], s2_y[i], Delta, modulus);
      if (h_x != NULL)
	  mpres_mul_z_to_z (h_x[i], s_x[i], f[i], modulus);
      if (h_y != NULL)
	  mpres_mul_z_to_z (h_y[i], s_y[i], f[i], modulus);
      if (h_x_ntt != NULL)
      {
	  mpres_mul_z_to_z (mt, s_x[i], f[i], modulus);
	  mpzspv_from_mpzv (h_x_ntt, i, &mt, 1UL, ntt_context);
      }
      if (h_y_ntt != NULL)
      {
	  mpres_mul_z_to_z (mt, s_y[i], f[i], modulus);
	  mpzspv_from_mpzv (h_y_ntt, i, &mt, 1UL, ntt_context);
      }
  }
  
  
  /* Compute the remaining r^((k+i)^2) values according to Peter's 
     recurrence */
  
  for (i = 2; i < l; i++)
    {
      if (h_x != NULL || h_x_ntt != NULL)
	{
	  /* r[i] = r2[i-1] * v[i-2] - r2[i-2], with indices of r2 and i taken
	     modulo 2 */
	  mpres_mul (s_x[i % 3], s2_x[1 - i % 2], v[i % 2], modulus);
	  mpres_sub (s_x[i % 3], s_x[i % 3], s2_x[i % 2], modulus);
	  
	  /* r2[i] = r2[i-1] * v[i-1] - r[i-2] */
	  mpres_mul (s2_x[i % 2], s2_x[1 - i % 2], v[1 - i % 2], modulus);
	  mpres_sub (s2_x[i % 2], s2_x[i % 2], s_x[(i - 2) % 3], modulus);
	  if (h_x != NULL)
	    mpres_mul_z_to_z (h_x[i], s_x[i % 3], f[i], modulus);
	  if (h_x_ntt != NULL)
	    {
	      mpres_mul_z_to_z (mt, s_x[i % 3], f[i], modulus);
	      mpzspv_from_mpzv (h_x_ntt, i, &mt, 1UL, ntt_context);
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
	    mpres_mul_z_to_z (h_y[i], s_y[i % 3], f[i], modulus);
	  if (h_y_ntt != NULL)
	    {
	      mpres_mul_z_to_z (mt, s_y[i % 3], f[i], modulus);
	      mpzspv_from_mpzv (h_y_ntt, i, &mt, 1UL, ntt_context);
	    }
	}
      
      /* v[i] = v[i - 1] * V_2(a + 1/a) - v[i - 2] */
      mpres_mul (newtmp[0], v[1 - i % 2], *V2, modulus);
      mpres_sub (v[i % 2], newtmp[0], v[i % 2], modulus);
    }

  timestop = cputime ();
  outputf (OUTPUT_VERBOSE, " took %lums\n", timestop - timestart);

  if (h_x != NULL && h_y != NULL && test_verbose (OUTPUT_TRACE))
    {
      for (i = 0; i < l; i++)
	gmp_printf ("/* pp1_sequence_h */ (rn^((k+%lu)^2) * f_%lu) == (%Zd"
		    " + %Zd * w) /* PARI C */\n", i, i, h_x[i], h_y[i]);
    }

  mpz_clear (mt);
}


int 
pp1fs2 (mpz_t f, const mpres_t X, mpmod_t modulus, 
	const faststage2_param_t *params)
{
  unsigned long nr;
  unsigned long i, l, lenF, lenH, lenG, lenR, tmplen;
  sets_long_t *S_1; /* This is stored as a set of sets (arithmetic 
                       progressions of prime length */
  set_long_t *S_2; /* This is stored as a regular set */
  listz_t F;   /* Polynomial F has roots X^{k_1} for k_1 \in S_1, so has 
		  degree s_1. It is symmetric, so has only s_1 / 2 + 1 
		  distinct coefficients. The sequence h_j will be stored in 
		  the same memory and won't be a monic polynomial, so the 
		  leading 1 monomial of F will be stored explicitly. Hence we 
		  need s_1 / 2 + 1 entries. */

  listz_t g_x, g_y, fh_x, fh_y, h_x, h_y, tmp, R_x, R_y; 
  const unsigned long tmpreslen = 20UL;
  mpres_t b1_x, b1_y, Delta, tmpres[tmpreslen];
  mpz_t mt;   /* All-purpose temp mpz_t */
  int youpi = ECM_NO_FACTOR_FOUND;
  long timetotalstart, timestart, timestop;

  timetotalstart = cputime ();

  ASSERT_ALWAYS (eulerphi (params->P) == params->s_1 * params->s_2);
  ASSERT_ALWAYS (params->s_1 < params->l);
  nr = params->l - params->s_1; /* Number of points we evaluate */

  outputf (OUTPUT_TRACE, "compare(a, b, n) = if (a != b,print(\"In PARI \", "
	   "n, \" line, \", a \" != \" b)); /* PARI */\n");

  if (make_S_1_S_2 (&S_1, &S_2, params) == ECM_ERROR)
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
  g_x = init_list2 (lenG, (unsigned int) abs (modulus->bits));
  g_y = init_list2 (lenG, (unsigned int) abs (modulus->bits));
  R_x = init_list2 (lenR, (unsigned int) abs (modulus->bits));
  R_y = init_list2 (lenR, (unsigned int) abs (modulus->bits));
  tmplen = 3UL * params->l + list_mul_mem (params->l / 2) + 20;
  outputf (OUTPUT_DEVVERBOSE, "tmplen = %lu\n", tmplen);
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
  V (tmpres[0], X, 2UL, modulus);
  i = poly_from_sets_V (F, tmpres[0], S_1, tmp, tmplen, modulus, NULL, NULL);
  ASSERT_ALWAYS(2 * i == params->s_1);
  ASSERT(mpz_cmp_ui (F[i], 1UL) == 0);
  free (S_1);
  S_1 = NULL;
  
  timestop = cputime ();
  outputf (OUTPUT_VERBOSE, " took %lums\n", timestop - timestart);
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
      mpres_get_z (mt, Delta, modulus);
      outputf (OUTPUT_TRACE, 
	       "Delta = Mod(%Zd, N); w = quadgen (4*lift(Delta)); b_1 = ", mt);
      gfp_ext_print (b1_x, b1_y, modulus, OUTPUT_TRACE);
      outputf (OUTPUT_TRACE, "; /* PARI */\n");
      outputf (OUTPUT_TRACE, "X == b_1 + 1/b_1 /* PARI C */\n");
    }

  /* Compute the h sequence h_j = b1^(P*-j^2) * f_j for 0 <= j <= s_1 */
  pp1_sequence_h (fh_x, fh_y, NULL, NULL, F, b1_x, b1_y, 0L, 
		  params->s_1 / 2 + 1, params->P, Delta, modulus, 
		  NULL, tmpreslen, tmpres);
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
      const long M = params->l - 1 - params->s_1 / 2;
      pp1_sequence_g (g_x, g_y, NULL, NULL, b1_x, b1_y, params->P, 
		      Delta, M, params->l, params->m_1, S_2->elem[l], 
		      modulus, NULL);
      
      /* Do the two convolution products */
      outputf (OUTPUT_VERBOSE, "TMulGen of g_x and h_x");
      timestart = cputime ();
      TMulGen (R_x, nr - 1, h_x, params->s_1, g_x, params->l - 1, tmp,
	       modulus->orig_modulus);
      timestop = cputime ();
      outputf (OUTPUT_VERBOSE, " took %lums\n", timestop - timestart);
      outputf (OUTPUT_VERBOSE, "TMulGen of g_y and h_y");
      timestart = cputime ();
      TMulGen (R_y, nr - 1, h_y, params->s_1, g_y, params->l - 1, tmp,
	       modulus->orig_modulus);
      timestop = cputime ();
      outputf (OUTPUT_VERBOSE, " took %lums\n", timestop - timestart);
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
      if (test_verbose(OUTPUT_RESVERBOSE))
      {
	  mpres_get_z (mt, tmpres[1], modulus);
	  outputf (OUTPUT_RESVERBOSE, "Product of R[i] = %Zd (times some "
		   "power of 2 if REDC was used! Try -mpzmod)\n", mt);
      }
      mpres_gcd (mt, tmpres[1], modulus);

      timestop = cputime ();
      outputf (OUTPUT_VERBOSE, "Computing product of F(g_i)^(1) took %lums\n", 
	       timestop - timestart);
      
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
 
  timestop = cputime ();
  outputf (OUTPUT_NORMAL, "Step 2 took %ldms\n", 
           timestop - timetotalstart);

  return youpi;
}


int 
pp1fs2_ntt (mpz_t f, const mpres_t X, mpmod_t modulus,
	    const faststage2_param_t *params, const int twopass)
{
  unsigned long nr;
  unsigned long i, l, lenF, tmplen;
  sets_long_t *S_1; /* This is stored as a set of sets (arithmetic 
                       progressions of prime length */
  set_long_t *S_2; /* This is stored as a regular set */
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
  listz_t tmp;
  mpzspm_t ntt_context;
  mpzspv_t g_x_ntt, g_y_ntt, h_x_ntt, h_y_ntt;
  const unsigned long tmpreslen = 20;
  mpres_t b1_x, b1_y, Delta, tmpres[tmpreslen];
  mpz_t mt;   /* All-purpose temp mpz_t */
  int youpi = ECM_NO_FACTOR_FOUND;
  long timetotalstart, timestart, timestop;

  timetotalstart = cputime ();

  ASSERT_ALWAYS (eulerphi (params->P) == params->s_1 * params->s_2);
  ASSERT_ALWAYS (params->s_1 < params->l);
  nr = params->l - params->s_1; /* Number of points we evaluate */

  outputf (OUTPUT_TRACE, "compare(a, b, n) = if (a != b,print(\"In PARI \", "
	   "n, \" line, \", a \" != \" b)); /* PARI */\n");

  if (make_S_1_S_2 (&S_1, &S_2, params) == ECM_ERROR)
      return ECM_ERROR;

  /* Init the NTT context */
  ntt_context = mpzspm_init (MAX(params->l, 2 << ceil_log2 (params->s_1)), 
                             modulus->orig_modulus);

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
  F = init_list2 (lenF, (unsigned int) abs (modulus->bits));
  tmplen = 3UL * params->l + list_mul_mem (params->l / 2) + 20;
  outputf (OUTPUT_DEVVERBOSE, "tmplen = %lu\n", tmplen);
  tmp = init_list2 (tmplen, (unsigned int) abs (modulus->bits));
  /* Allocate memory for h_ntt */
  h_x_ntt = mpzspv_init_mt (params->l / 2 + 1, ntt_context);

  if (test_verbose (OUTPUT_TRACE))
    {
      mpres_get_z (mt, X, modulus); /* mpz_t copy of X for printing */
      outputf (OUTPUT_TRACE, "N = %Zd; X = Mod(%Zd, N); /* PARI */\n", 
	       modulus->orig_modulus, mt);
    }
  
  /* Compute the polynomial f(x) = \prod_{k_1 in S_1} (x - X^{2 k_1}) */
  outputf (OUTPUT_VERBOSE, "Computing F from factored S_1");
  
  timestart = cputime ();
  V (tmpres[0], X, 2UL, modulus);
  i = poly_from_sets_V (F, tmpres[0], S_1, tmp, tmplen, modulus, h_x_ntt,
                        ntt_context);
  ASSERT_ALWAYS(2 * i == params->s_1);
  ASSERT(mpz_cmp_ui (F[i], 1UL) == 0);
  free (S_1);
  S_1 = NULL;
  
  timestop = cputime ();
  outputf (OUTPUT_VERBOSE, " took %lums\n", timestop - timestart);
  if (test_verbose (OUTPUT_TRACE))
    {
      for (i = 0; i < params->s_1 / 2 + 1; i++)
	outputf (OUTPUT_TRACE, "f_%lu = %Zd; /* PARI */\n", i, F[i]);
      outputf (OUTPUT_TRACE, "f(x) = f_0");
      for (i = 1; i < params->s_1 / 2 + 1; i++)
	outputf (OUTPUT_TRACE, "+ f_%lu * (x^%lu + x^(-%lu))", i, i, i);
      outputf (OUTPUT_TRACE, "/* PARI */ \n");
    }
  clear_list (tmp, tmplen); /* The NTT variant won't use tmp anymore */

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
  h_y_ntt = mpzspv_init_mt (params->l / 2 + 1, ntt_context);
  /* Compute the h_j sequence */
  pp1_sequence_h (NULL, NULL, h_x_ntt, h_y_ntt, F, b1_x, b1_y, 0L, 
		  params->s_1 / 2 + 1, params->P, Delta, modulus, 
		  ntt_context, tmpreslen, tmpres);
  /* We don't need F(x) any more */
  clear_list (F, lenF);

  /* If we use NTT, compute the forward transform of h and store the distinct
     coefficients in h_ntt */

  /* Compute forward transform of h */
  g_x_ntt = mpzspv_init_mt (params->l, ntt_context);
  if (twopass)
    {
      g_y_ntt = g_x_ntt;
      R = init_list (nr);
    }
  else
    g_y_ntt = mpzspv_init_mt (params->l, ntt_context);
  
  /* Compute DCT-I of h_x and h_y */
  outputf (OUTPUT_VERBOSE, "Computing DCT-I of h_x");
  timestart = cputime ();
  ntt_spv_to_dct (h_x_ntt, h_x_ntt, params->s_1 / 2 + 1, params->l / 2 + 1,
		  g_x_ntt, ntt_context);
  timestop = cputime ();
  outputf (OUTPUT_VERBOSE, " took %lums\n", timestop - timestart);

  outputf (OUTPUT_VERBOSE, "Computing DCT-I of h_y");
  timestart = cputime ();
  ntt_spv_to_dct (h_y_ntt, h_y_ntt, params->s_1 / 2 + 1, params->l / 2 + 1,
		  g_x_ntt, ntt_context);
  timestop = cputime ();
  outputf (OUTPUT_VERBOSE, " took %lums\n", timestop - timestart);


  for (l = 0; l < params->s_2; l++)
    {
      const long M = params->l - 1 - params->s_1 / 2;

      if (twopass)
	{
	  pp1_sequence_g (NULL, NULL, g_x_ntt, NULL, b1_x, b1_y, params->P, 
			  Delta, M, params->l, params->m_1, S_2->elem[l], 
			  modulus, ntt_context);

	  /* Do the convolution product of g_x * h_x */
	  outputf (OUTPUT_VERBOSE, "Computing g_x*h_x");
	  timestart = cputime ();
	  ntt_mul_by_dct (g_x_ntt, h_x_ntt, params->l, ntt_context);
	  /* Store the product coefficients we want in R */
	  mpzspv_to_mpzv (g_x_ntt, params->s_1 / 2, R, nr, ntt_context);
	  timestop = cputime ();
	  outputf (OUTPUT_VERBOSE, " took %lums\n", timestop - timestart);

	  /* Compute g_y sequence */
	  pp1_sequence_g (NULL, NULL, NULL, g_y_ntt, b1_x, b1_y, params->P, 
			  Delta, M, params->l, params->m_1, S_2->elem[l], 
			  modulus, ntt_context);
	  
	  /* Do the convolution product of g_y * (Delta * h_y) */
	  outputf (OUTPUT_VERBOSE, "Computing g_y*h_y");
	  timestart = cputime ();
	  ntt_mul_by_dct (g_y_ntt, h_y_ntt, params->l, ntt_context);
	  timestop = cputime ();
	  outputf (OUTPUT_VERBOSE, " took %lums\n", timestop - timestart);
	  
	  /* Compute product of sum of coefficients and gcd with N */
	  ntt_gcd (mt, g_y_ntt, params->s_1 / 2, R, nr, ntt_context, modulus);
	}
      else
	{
	  /* Currently does not work */
	  abort();
	  
	  pp1_sequence_g (NULL, NULL, g_x_ntt, g_y_ntt, b1_x, b1_y, params->P, 
			  Delta, M, params->l, params->m_1, S_2->elem[l], 
			  modulus, ntt_context);

	  
	  outputf (OUTPUT_VERBOSE, "Computing forward NTT of g_x");
	  timestart = cputime ();
	  mpzspv_to_ntt (g_x_ntt, (spv_size_t) 0, (spv_size_t) params->l, 
	                 (spv_size_t) params->l, 0, ntt_context);
	  timestop = cputime ();
	  outputf (OUTPUT_VERBOSE, " took %lums\n", timestop - timestart);
	  
	  outputf (OUTPUT_VERBOSE, "Computing point-wise product of g_x and h_x");
	  timestart = cputime ();
	  ntt_dft_mul_dct (g_x_ntt, h_x_ntt, params->l, ntt_context);
	  timestop = cputime ();
	  outputf (OUTPUT_VERBOSE, " took %lums\n", timestop - timestart);
	  
	  outputf (OUTPUT_VERBOSE, "Computing forward NTT of g_y");
	  timestart = cputime ();
	  mpzspv_to_ntt (g_y_ntt, (spv_size_t) 0, (spv_size_t) params->l, 
	                 (spv_size_t) params->l, 0, ntt_context);
	  timestop = cputime ();
	  outputf (OUTPUT_VERBOSE, " took %lums\n", timestop - timestart);
	  
	  outputf (OUTPUT_VERBOSE, "Computing point-wise product of g_y and h_y");
	  timestart = cputime ();
	  ntt_dft_mul_dct (g_y_ntt, h_y_ntt, params->l, ntt_context);
	  timestop = cputime ();
	  outputf (OUTPUT_VERBOSE, " took %lums\n", timestop - timestart);
	  
	  outputf (OUTPUT_VERBOSE, "Adding and computing inverse NTT of sum");
	  timestart = cputime ();
	  mpzspv_add (g_x_ntt, (spv_size_t) 0, g_x_ntt, (spv_size_t) 0, 
	              g_y_ntt, (spv_size_t) 0, params->l, ntt_context);
	  mpzspv_from_ntt (g_x_ntt, (spv_size_t) 0, params->l, (spv_size_t) 0,
	                   ntt_context);
	  timestop = cputime ();
	  outputf (OUTPUT_VERBOSE, " took %lums\n", timestop - timestart);
	  
	  ntt_gcd (mt, g_x_ntt, params->s_1 / 2, twopass ? R : NULL, nr, 
		   ntt_context, modulus);
	}
      
      if (mpz_cmp_ui (mt, 1UL) > 0)
	{
	  mpz_set (f, mt);
	  youpi = ECM_FACTOR_FOUND_STEP2;
	  break;
	}
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
  for (i = 0; i < tmpreslen; i++)
      mpres_clear (tmpres[i], modulus);
 
  timestop = cputime ();
  outputf (OUTPUT_NORMAL, "Step 2 took %ldms\n", 
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
      pn = choose_P (B2min, B2, lmax, 1UL, &params, NULL, NULL);
    }

  d = eulerphi (pn);

  L = get_factored_sorted_sets (&setsize, pn);
  if (L == NULL)
  {
      printf ("Error, get_factored_sorted_sets() returned NULL pointer\n");
      exit (EXIT_FAILURE);
  }

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
      Simple check of list_mul_reciprocal() 
  ************************************************************/

  for (i = 0; i < 4; i++)
    mpz_set_ui (tmp[i], i + 2);
  /* tmp[0 .. 3] = [2, 3, 4, 5] */

  /* Compute (5*(x^3+x^{-3}) + 4*(x^2+x^{-2}) + 3*(x+x^{-1}) + 2)^2 =
     25*(x^6+x^{-6}) + 40*(x^5+x^{-5}) + 46*(x^4+x^{-4}) + 44*(x^3+x^{-3}) + 
     55*(x^2+x^{-2}) + 76*(x+x^{-1}) + 104, so we expect in 
     tmp[0 .. 6] = [104, 76, 55, 44, 46, 40, 25] */
  list_mul_reciprocal (tmp, tmp, 4, tmp, 4, modulus->orig_modulus, 
		       tmp + 7, tmplen - 7);

  if (mpz_cmp_ui (tmp[6], 25UL) != 0 || mpz_cmp_ui (tmp[5], 40UL) != 0 ||
      mpz_cmp_ui (tmp[4], 46UL) != 0 || mpz_cmp_ui (tmp[3], 44UL) != 0 ||
      mpz_cmp_ui (tmp[2], 55UL) != 0 || mpz_cmp_ui (tmp[1], 76UL) != 0 ||
      mpz_cmp_ui (tmp[0], 104UL) != 0)
    {
      list_output_poly (tmp, 7, 0, 0, "Error, list_mul_reciprocal produced ", 
			"\n", OUTPUT_ERROR);
      abort ();
    }

  for (i = 0; i < 4; i++)
    mpz_set_ui (tmp[i], i + 1);
  /* tmp[0 .. 3] = [1, 2, 3, 4] = 4*x^3 + 3*x^2 + 2*x + 1 */

  /* Compute (4*(x^3+x^{-3}) + 3*(x^2+x^{-2}) + 2*(x+x^{-1}) + 1) *
     (3*(x^2+x^-2) + 2*(x+x^-1) + 1), so we expect in 
     tmp[0 .. 5] = [27, 28, 18, 16, 17, 12] */
  list_mul_reciprocal (tmp, tmp, 4, tmp, 3, modulus->orig_modulus, 
		       tmp + 6, tmplen - 6);

  if (mpz_cmp_ui (tmp[0], 27UL) != 0 || mpz_cmp_ui (tmp[1], 28UL) != 0 ||
      mpz_cmp_ui (tmp[2], 18UL) != 0 || mpz_cmp_ui (tmp[3], 16UL) != 0 ||
      mpz_cmp_ui (tmp[4], 17UL) != 0 || mpz_cmp_ui (tmp[5], 12UL) != 0)
    {
      list_output_poly (tmp, 6, 0, 0, "Error, list_mul_reciprocal produced ", 
			"\n", OUTPUT_ERROR);
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
	list_output_poly (tmp, len, 0, 1, "Input to list_scale_V: ", "\n", 
			  OUTPUT_TRACE);
	
	list_scale_V (tmp + len, tmp, Q, len - 1, modulus, tmp + 3*len, 
		      tmplen - 3*len);
	
	list_mod (tmp + len, tmp + len, 2 * len - 1, modulus->orig_modulus);
	list_output_poly (tmp + len, 2*len-1, 0, 1, 
			  "Output of list_scale_V: ", "\n", OUTPUT_TRACE);
	
	/* With Q = 2 = 1 + 1/1, gamma = 1 and F(gamma*X)*F(1/gamma *X) =F(X)^2
	   Compare with a simple symmetic multiply */
	list_mul_reciprocal (tmp + 3 * len, tmp, len, tmp, len, 
		             modulus->orig_modulus, tmp + 5*len, tmplen - 5*len);
	
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
      
      t = set_of_sums (sumset, L, setsize, 0L);
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
	if (gcd (i, pn) == 1UL)
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
