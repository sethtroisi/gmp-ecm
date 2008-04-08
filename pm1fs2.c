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
#if HAVE_ALLOCA_H
#include <alloca.h>
#endif
#if HAVE_STRING_H
#include <string.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

/* TODO:
   - fix parameter selection
   - allow choosing parameters according to -maxmem
   - check one pass P+1 stage 2 carefully, allow its use
   - move functions into their proper files (i.e. NTT functions etc.)
   - later: allow storing NTT vectors on disk
*/

/* Define TEST_ZERO_RESULT to test if any result of the multipoint
   evaluation is equal to zero. If the modulus is composite, this
   happening might indicate a problem in the evalutaion code */
#define TEST_ZERO_RESULT

const int pari = 0;

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
ntt_sqr_reciprocal (mpzv_t, const mpzv_t, mpzspv_t, const spv_size_t, 
		    const mpzspm_t);

static void
print_elapsed_time (long cpu, ATTRIBUTE_UNUSED long real)
{
#ifdef _OPENMP
  if (real != 0L)
    outputf (OUTPUT_VERBOSE, " took %lums (%lums real)\n", cpu, real);
  else
    outputf (OUTPUT_VERBOSE, " took %lums\n", cpu);
#else
  outputf (OUTPUT_VERBOSE, " took %lums\n", cpu);
#endif
}


static void
print_CRT_primes (const int verbosity, const char *prefix, 
		   const mpzspm_t ntt_context)
{
  double modbits = 0.;
  unsigned int i;
  
  if (test_verbose (verbosity))
    {
      outputf (verbosity, "%s = %lu", prefix, ntt_context->spm[0]->sp);
      modbits += log ((double) ntt_context->spm[0]->sp);
      for (i = 1; i < ntt_context->sp_num; i++)
	{
	  outputf (verbosity, " * %lu", ntt_context->spm[i]->sp);
	  modbits += log ((double) ntt_context->spm[i]->sp);
	}
      outputf (verbosity, ", has %d primes, %f bits\n", 
               ntt_context->sp_num, modbits / log (2.));
    }
}

/* Approximate amount of memory in bytes each coefficient in an NTT takes 
   so that NTT can do transforms up to length lmax with modulus, or
   with 2*modulus if twice != 0 */
static size_t
ntt_coeff_mem (const unsigned long lmax, const mpz_t modulus, const int twice)
{
  mpz_t t;
  size_t n;
  
  mpz_init (t);
  mpz_mul (t, modulus, modulus);
  mpz_mul_ui (t, t, lmax);
  if (twice)
    mpz_mul_2exp (t, t, 1UL);
  /* +4: +1 for rounding up, +3 for extra words due to ECRT */
  n = (mpz_sizeinbase (t, 2) - 1) / SP_NUMB_BITS + 4;
  mpz_clear (t);
  return n * sizeof (sp_t);
}

size_t
pm1fs2_memory_use (const unsigned long lmax, const mpz_t modulus, 
		   const int use_ntt)
{
  if (use_ntt)
    {
      /* We store lmax / 2 + 1 coefficients for the DCT-I of F and lmax 
	 coefficients for G in NTT ready format. Each coefficient in 
	 NTT-ready format occupies approx. 
	 ceil(log(lmax*modulus^2)/log(bits per sp_t)) + 3 words. */
      
      size_t n;
  
      
      n = ntt_coeff_mem (lmax, modulus, 0) * (size_t) (3 * lmax / 2 + 1);
      outputf (OUTPUT_DEVVERBOSE, "pm1fs2_memory_use: Estimated memory use "
	       "with lmax = %lu NTT is %lu bytes\n", lmax, n);
      return n;
    }
  else
    {
      /* F stores s_1/2 residues,
	 h stores s_1 mpz_t structs (residues get cloned from F)
	 g stores lmax residues, 
	 R stores lmax-s_1 residues, 
	 and tmp stores 3*lmax+list_mul_mem (lmax / 2) residues.
	 Assume s_1 is close to lmax/2.
	 Then we have 
	 lmax/4 + lmax/2 + lmax + lmax/2 + 3*lmax + list_mul_mem (lmax / 2)
	 = (5+1/4)*lmax + list_mul_mem (lmax / 2) residues, plus s_1 mpz_t.
      */
      
      size_t n;
      
      n = mpz_size (modulus) * sizeof (mp_limb_t) + sizeof (mpz_t);
      n *= 5 * lmax + lmax / 4 + list_mul_mem (lmax / 2);
      n += lmax / 2 * sizeof (mpz_t);
      /* Memory use due to temp space allocation in TMulKS appears to 
	 approximately triple the estimated memory use. This is hard to
	 estimate precisely, so let's got with the fudge factor of 3 here */
      n *= 3;
      outputf (OUTPUT_DEVVERBOSE, "pm1fs2_memory_use: Estimated memory use "
	       "with lmax = %lu is %lu bytes\n", lmax, n);
      return n;
    }
}

/* return the possible lmax for given memory use and modulus */

unsigned long
pm1fs2_maxlen (const size_t memory, const mpz_t modulus, const int use_ntt)
{
  if (use_ntt)
    {
      size_t n, lmax = 1;
  
      n = ntt_coeff_mem (lmax, modulus, 0);
      lmax = 1UL << (ceil_log2 (2 * memory / n / 3) - 1);
      return lmax;
    }
  else
    {
      size_t lmax, n;
      
      n = mpz_size (modulus) * sizeof (mp_limb_t) + sizeof (mpz_t);

      /* Guess an initial value of lmax for list_mul_mem (lmax / 2) */
      /* memory = n * 25/4 * lmax + lmax / 2 * sizeof (mpz_t); */
      /* Fudge factor of 3 for TMulKS as above */
      lmax = memory / (3 * n * 25 / 4 + 1 / 2 * sizeof (mpz_t));
      return lmax;
    }
}

size_t 
pp1fs2_memory_use (const unsigned long lmax, const mpz_t modulus, 
		   const int use_ntt, const int twopass)
{
  size_t n, m;
  
  m = mpz_size (modulus) * sizeof (mp_limb_t) + sizeof (mpz_t);
  if (use_ntt)
    {
      /* In one pass mode, we store h_x_ntt and h_y_ntt, each of length 
	 lmax/2(+1), and g_x_ntt and g_y_ntt, each of length lmax, all in 
	 NTT ready format. In two pass mode, we store h_x_ntt, h_y_ntt and 
	 g_x_ntt as before, plus R which is lmax - s_1 mpz_t. 
	 We assume s_1 ~= lmax/2.
      */

      n = ntt_coeff_mem (lmax, modulus, !twopass);
      if (twopass)
	return lmax * (2 * n + m / 2);
      else
	return lmax * 3 * n;
    }
  else
    {
      /* We allocate:
	 F: s_1/2 coefficients
	 fh_x, fh_y: s_1/2 coefficients
	 h_x, h_y: s_1 mpz_t's (cloned from fh_x and fh_y)
	 g_x, g_y: lmax coefficients
	 R_x, R_y: lmax - s_1 coefficients
	 tmp: 3UL * lmax + list_mul_mem (lmax / 2)
	 Assuming s_1 ~ lmax/2, that's
	 lmax/2 + 2*lmax/4 + 2*lmax + 2*lmax/2 * 3*lmax + 
           list_mul_mem (lmax / 2) =
	 7 + list_mul_mem (lmax / 2) coefficients and lmax mpz_t.
       */
      
      n = m * (7 * lmax + list_mul_mem (lmax / 2));
      n += lmax * sizeof (mpz_t);
      n = 5 * n / 2; /* A fudge factor again */
      return n;
    }
}

unsigned long 
pp1fs2_maxlen (const size_t memory, const mpz_t modulus, const int use_ntt, 
	       const int twopass)
{
  size_t n, m;
  
  m = mpz_size (modulus) * sizeof (mp_limb_t) + sizeof (mpz_t);
  if (use_ntt)
    {
      n = ntt_coeff_mem (1, modulus, !twopass);
      if (twopass)
	n = memory / (2 * n + m / 2);
      else
	n = memory / (3 * n);
      return 1UL << (ceil_log2 (n) - 1); /* Rounded down to power of 2 */
    }
  else
    {
      return 2 * memory / 5 / (m * 8 + sizeof (mpz_t));
    }
}


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
  unsigned long m = maxS (P);
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
   under those conditions. If use_ntt == 1, we require s_1 < l.
   Returns 0 if no such choice is possible */

static unsigned long 
choose_s_1 (const unsigned long phiP, const unsigned long min_s2,
	    const unsigned long l, const int use_ntt)
{
  const int nrprimes = sizeof (phiPfactors) / sizeof (unsigned long);
  /* Using [nrprimes] here makes the compiler complain about variable-sized
     arrays */
  int phiPexponents[sizeof (phiPfactors) / sizeof (unsigned long)], 
    exponents[sizeof (phiPfactors) / sizeof (unsigned long)];
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
	  absdiff_ul (trys_1, l) < absdiff_ul (s_1, l) &&
	  (use_ntt == 0 || trys_1 < l))
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
	  mpz_t final_B2min, mpz_t final_B2, const int use_ntt)
{
  /* Let S_1 + S_2 == (Z/PZ)* (mod P).

     Let F(x) = \prod_{k_1 \in S_1} (x - b_1^{2 k_1}).

     If we evaluate F(b_1^{2 k_2 + (2m + 1)P}) for all k_2 \in S_2 with 
     m_1 <= m < m_1+nr, we test all exponents 2 k_2 + (2m + 1)P - 2 k_1.
     The largest value coprime to P at the low end of the stage 2 interval 
     *not* covered will be 
       2*max(S_2) + (2*(m_1-1) + 1)*P - 2*min(S_1).
     The smallest value at the high end not covered will be
       2*min(S_2) + (2*(m_1 + nr) + 1)*P - 2*max(S_1).
     Assume S_1 and S_2 are symmetric around 0, so that max(S_1) = -min(S_1).
     Then the largest ... is:
       2*(max(S_1) + max(S_2)) + (2*m_1 - 1)*P
     The smallest ... is:
       -2*(max(S_1) + max(S_2)) + (2*m_1 + 2*nr + 1)*P
     The effective B2min = 2*(max(S_1) + max(S_2)) + (2*m_1 - 1)*P + 1
     The effective B2max = -2*(max(S_1) + max(S_2)) + (2*m_1 + 2*nr + 1)*P - 1

     Then the difference effB2max - effB2min =
       -4*(max(S_1) + max(S_2)) + 2P*(nr + 1) - 2

     We obviously require B2max - B2min <= 2*nr*P
     Since nr < lmax, B2max - B2min <= 2*lmax*P or
     P >= ceil((B2max - B2min)/(2*lmax))

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

  /* Find the smallest P that can cover an interval of length B2 - B2min */
  mpz_sub (B2l, B2, B2min);
  mpz_cdiv_q_ui (t, B2l, 2UL*lmax);
  outputf (OUTPUT_TRACE, "choose_P: We need P >= %Zd\n", t);
  for (i = 0; i < Pvalues_len; i++)
    if (mpz_cmp_ui (t, Pvalues[i]) <= 0)
      break;

  for ( ; i <  Pvalues_len; i++)
    {
      unsigned long phiP;
      /* Now a careful check to see if this P is large enough */
      P = Pvalues[i];
      phiP = eulerphi (P);
      s_1 = choose_s_1 (phiP, min_s2, lmax / 2, use_ntt);
      if (s_1 == 0)
	continue;
      s_2 = phiP / s_1;
      outputf (OUTPUT_TRACE, "choose_P: Testing P = %lu, phiP = %lu, "
	       "s_1 = %lu, s_2 = %lu, nr = %lu\n", 
	       P, phiP, s_1, s_2, lmax - s_1);
      if (test_P (B2min, B2, m_1, P, lmax - s_1, effB2min, effB2))
	{
	    outputf (OUTPUT_TRACE, 
		     "choose_P: This P is acceptable, B2 = %Zd\n", effB2);
	    break;
	}
      else
	outputf (OUTPUT_TRACE, "choose_P: Not good enough, trying next P\n");
  }
  
  if (i == Pvalues_len)
    return ECM_ERROR; /* Could not find suitable P */

  /* For each P, determine the smallest s_2 * lmax that lets us cover
     B2min, B2max */


  /* s_2 cannot be reduced further, but the transform length l could. */
  l = lmax;
  while (l / 2 > s_1 && test_P (B2min, B2, m_1, P, l / 2 - s_1, effB2min, t))
  {
      l /= 2;
      outputf (OUTPUT_TRACE, "choose_P: Reducing transform length to "
	       "%ld, effB2 = %Zd\n", l, t);
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
      trys_1 = choose_s_1 (tryphiP, min_s2, l / 2, use_ntt);
      if (trys_1 == 0)
	continue;
      trys_2 = tryphiP / trys_1;
      outputf (OUTPUT_TRACE, "choose_P: Trying if P = %lu, phiP = %lu, "
	       "s_1 = %lu, s_2 = %lu works as well\n", 
	       tryP, tryphiP, trys_1, trys_2);
      if (trys_2 > s_2) /* We want to keep the minimal */
      {                 /* number of multipoint evaluations */
	  outputf (OUTPUT_TRACE, "choose_P: No, s_2 would become %lu\n", 
		   trys_2);
	  /* break; */
	  continue;
      }
      if (!test_P (B2min, B2, m_1, tryP, l - trys_1, effB2min, t))
      {
	  outputf (OUTPUT_TRACE, "choose_P: No, does not cover B2min - B2 "
		   "range, effB2 = %Zd\n", t);
      }
      else
      {
	  if (mpz_cmp (t, effB2) >= 0)
	  {
	      outputf (OUTPUT_TRACE, 
		       "choose_P: Yes, works and gives higher B2 = %Zd\n", t);
	      P = tryP;
	      s_1 = trys_1;
	      s_2 = trys_2;
	      mpz_set (effB2, t);
	      while (l / 2 > s_1 && 
		     test_P (B2min, B2, m_1, P, l / 2 - s_1, effB2min, t))
	      {
		  l /= 2;
		  outputf (OUTPUT_TRACE, "choose_P: Reducing transform length "
			   "to %ld, B2 = %Zd\n", l, t);
		  mpz_set (effB2, t);
	      }
	  }
	  else
	  {
	      outputf (OUTPUT_TRACE, "choose_P: Works, but does not give "
		       "higher B2, %Zd <= %Zd\n", t, effB2);
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

ATTRIBUTE_UNUSED
static void
list_recip_eval1 (mpz_t R, const listz_t S, const unsigned long l)
{
  unsigned long i;

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
list_mul_reciprocal (listz_t R, listz_t S1, unsigned long l1, 
		     listz_t S2, unsigned long l2,
		     mpz_t modulus, listz_t tmp, 
		     ATTRIBUTE_UNUSED const unsigned long tmplen)
{
  unsigned long i;
  const unsigned long lmax = MAX(l1, l2);
  listz_t r1 = tmp, r2 = tmp + 2*lmax - 1, rev = tmp + 4*lmax - 2,
    t = tmp + 6*lmax - 3;
#ifdef WANT_ASSERT
  mpz_t sum1, sum2, prod;
#endif

  ASSERT (S1 < tmp || S1 >= tmp + tmplen);
  ASSERT (S2 < tmp || S2 >= tmp + tmplen);
  ASSERT (R < tmp || R >= tmp + tmplen);

  if (l1 == 0UL || l2 == 0UL)
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

   If NTT is used, needs 4 * deg + 3 entries in tmp.
   If no NTT is used, needs 4 * deg + 2 + (memory use of list_sqr_reciprocal)
*/

static void
list_scale_V (listz_t R, listz_t F, mpres_t Q, unsigned long deg,
	      mpmod_t modulus, listz_t tmp, const unsigned long tmplen,
	      mpzspv_t dct, const mpzspm_t ntt_context)
{
  mpres_t Vi_1, Vi, Vt;
  unsigned long i;
  const listz_t G = tmp, H = tmp + 2 * deg + 1, newtmp = tmp + 4 * deg + 2;
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
  
  /* Make sure newtmplen does not underflow */
  ASSERT_ALWAYS (tmplen >= 4 * deg + 2);
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

  mpz_mod (G[0], F[0], modulus->orig_modulus); /* G_0 = S_0 * V_0(Q)/2 = S_0 * 1 */
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
    ntt_sqr_reciprocal (G, G, dct, deg + 1, ntt_context);
  else
    list_sqr_reciprocal (G, G, deg + 1, modulus->orig_modulus, 
                         newtmp, newtmplen);

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

  /* H[i-1] = h_i =  F[i] * U_i(Q) / 2, for 1 <= i <= deg. h_0 is undefined
     and has no storage allocated */
  for (i = 1; i <= deg; i++)
    {
      /* Here, Vi = U_i(Q) / 2, Vi_1 = U_{i-1}(Q) / 2. */
      /* h_i = S_i * U_i(Q)/2 */
      mpres_mul_z_to_z (H[i - 1], Vi, F[i], modulus);
      outputf (OUTPUT_TRACE, 
	       "list_scale_V: H_%lu (in U_i basis) = F_%lu * U_%lu(Q)/2 = %Zd * %Zd = %Zd\n", 
	       i, i, i, F[i], Vi, H[i - 1]);
      
      mpres_mul (Vt, Vi, Q, modulus);
      mpres_sub (Vt, Vt, Vi_1, modulus);
      mpres_set (Vi_1, Vi, modulus); /* Could be a swap */
      mpres_set (Vi, Vt, modulus); /* Could be a swap */
    }

  /* Convert H to standard basis */
  /* We can do it in-place with H - 1 = H_U. */

  for (i = deg; i >= 3; i--)
    mpz_add (H[i - 3], H[i - 3], H[i - 1]);
  
  /* U_2(X+1/X) = (X^2 - 1/X^2)/(X-1/X) = X+1/X = V_1(X+1/X),
     so no addition occures here */
  /* if (deg >= 2)
     mpz_set (H[1], H[1]); Again, a no-op. */
  
  /* U_1(X+1/X) = 1, so this goes to coefficient of index 0 in std. basis */
  /* mpz_set (H[0], H[0]); Another no-op. */
  
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
    ntt_sqr_reciprocal (H, H, dct, deg, ntt_context);
  else
    list_sqr_reciprocal (H, H, deg, modulus->orig_modulus, 
  		         newtmp, newtmplen);

  /* Now there are the 2*deg-1 coefficients in standard basis of a 
     symmetric polynomial of degree 4*deg - 4 in H[0 ... 2*deg-2] */

  list_output_poly (H, 2*deg - 1, 0, 1, "/* list_scale_V */ H(x)^2 == ", "\n",
		    OUTPUT_TRACE);

  /* Multiply by Q^2-4 */
  mpres_mul (Vt, Q, Q, modulus);
  mpres_sub_ui (Vt, Vt, 4, modulus);
  for (i = 0; i <= 2 * deg - 2; i++)
    mpres_mul_z_to_z (H[i], Vt, H[i], modulus);
  list_output_poly (H, 2 * deg - 1, 0, 1, "/* list_scale_V */ "
		    "H(x)^2*(Q^2-4) == ", "\n", OUTPUT_TRACE);

  /* Multiply by (X - 1/X)^2 = X^2 - 2 + 1/X^2 and subtract from G*/
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


/* Build a polynomial with roots r^2i, i in the sumset of the sets in "sets".
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
                                        Qt = V_2k(Q) */
  V (Qt, Qt, 2UL, modulus);
  
  mpres_neg (Qt, Qt, modulus);
  mpres_get_z (F[0], Qt, modulus);
  mpz_set_ui (F[1], 1UL);
  deg = 1UL;
  /* Here, F(x) = (x - r^{2k_1})(x - r^{-2k_1}) / x = 
                  (x^2 - x (r^{2k_1} + r^{-2k_1}) + 1) / x =
		  (x + 1/x) - V_{2k_1}(r + 1/r) */

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
	  V (Qt, Qt, 2UL, modulus);
	  list_scale_V (F, F, Qt, deg, modulus, tmp, tmplen, dct, 
	                ntt_context);
	  deg *= 2UL;
	  ASSERT (mpz_cmp_ui (F[deg], 1UL) == 0); /* Check it's monic */
	}
      else
	{
	  ASSERT_ALWAYS (c % 2UL == 1UL);
	  /* Generate the F(Q^{2k_i} * X)*F(Q^{-2k_i} * X) polynomials.
	     Each is symmetric of degree 2*deg, so each has deg+1 coeffients
	     in standard basis. */
	  for (i = 0UL; i < (c - 1UL) / 2UL; i++)
	    {
              /* Check it's symmetric */
	      ASSERT (set->elem[i] == -set->elem[c - 1L - i]);
	      V (Qt, Q, set->elem[i], modulus);
	      V (Qt, Qt, 2UL, modulus);
	      ASSERT (mpz_cmp_ui (F[deg], 1UL) == 0); /* Check it's monic */
	      list_scale_V (F + (2UL * i + 1UL) * (deg + 1UL), F, Qt, deg, 
	                    modulus, tmp, tmplen, dct, ntt_context);
	      ASSERT (mpz_cmp_ui (F[(2UL * i + 1UL) * (deg + 1UL) + 2UL * deg], 
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

static int
build_F_ntt (listz_t F, const mpres_t P_1, sets_long_t *S_1, 
	     const faststage2_param_t *params, mpmod_t modulus)
{
  mpzspm_t F_ntt_context;
  mpzspv_t F_ntt;
  unsigned long tmplen;
  listz_t tmp;
  long timestart, realstart;
  int i;

  timestart = cputime ();
  realstart = realtime ();
  
  /* Precompute the small primes, primitive roots and inverses etc. for 
     the NTT. The code to multiply wants a 3*k-th root of unity, where 
     k is the smallest power of 2 with k > s_1/2 */
  
  F_ntt_context = mpzspm_init (3UL << ceil_log2 (params->s_1 / 2 + 1), 
			       modulus->orig_modulus);
  if (F_ntt_context == NULL)
    return ECM_ERROR;
  
  print_CRT_primes (OUTPUT_DEVVERBOSE, "CRT modulus for building F = ",
		    F_ntt_context);
  
  outputf (OUTPUT_VERBOSE, "Computing F from factored S_1");
  
  tmplen = params->s_1 + 100;
  tmp = init_list2 (tmplen, (unsigned int) abs (modulus->bits));
  F_ntt = mpzspv_init (1UL << ceil_log2 (params->s_1 / 2 + 1), F_ntt_context);
  
  i = poly_from_sets_V (F, P_1, S_1, tmp, tmplen, modulus, F_ntt, 
                        F_ntt_context);
  ASSERT_ALWAYS(2 * i == params->s_1);
  ASSERT_ALWAYS(mpz_cmp_ui (F[i], 1UL) == 0);
  
  print_elapsed_time (cputime() - timestart, realtime() - realstart);
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
  mpzspv_clear (F_ntt, F_ntt_context);
  F_ntt = NULL;
  mpzspm_clear (F_ntt_context);
  F_ntt_context = NULL;

  return 0;
}

/* Compute g_i = x_0^{M-i} * r^{(M-i)^2} for 0 <= i < l. 
   x_0 = b_1^{2*k_2 + (2*m_1 + 1) * P}. r = b_1^P. 
   Stores the result in g[0 ... l] and/or in g_ntt[offset ... offset + l] */

static void
pm1_sequence_g (listz_t g_mpz, mpzspv_t g_ntt, const mpres_t b_1, 
                const unsigned long P, const long M_param, 
		const unsigned long l_param, const mpz_t m_1, const long k_2, 
		mpmod_t modulus_param, const mpzspm_t ntt_context)
{
  mpres_t r[3], x_0, x_Mi;
  mpz_t t;
  unsigned long i;
  long timestart, realstart;
  long M = M_param;
  unsigned long l = l_param, offset = 0UL;
  mpmod_t modulus;
  int want_output = 1;

  outputf (OUTPUT_VERBOSE, "Computing g_i");
  outputf (OUTPUT_DEVVERBOSE, "\npm1_sequence_g: P = %lu, M_param = %lu, "
           "l_param = %lu, m_1 = %Zd, k_2 = %lu\n", 
	   P, M_param, l_param, m_1, k_2);
  timestart = cputime ();
  realstart = realtime ();

#ifdef _OPENMP
#pragma omp parallel if (l > 100) private(r, x_0, x_Mi, t, i, M, l, offset, modulus, want_output)
  {
    /* When multi-threading, we adjust the parameters for each thread */

    const int nr_chunks = omp_get_num_threads();
    const int thread_nr = omp_get_thread_num();
    
    l = (l_param - 1) / nr_chunks + 1; /* = ceil(l_param / nr_chunks) */
    offset = thread_nr * l;
    outputf (OUTPUT_DEVVERBOSE, 
             "pm1_sequence_g: thread %d has l = %lu, offset = %lu.\n", 
             thread_nr, l, offset);
    ASSERT_ALWAYS (l_param >= offset);
    l = MIN(l, l_param - offset);
    M = M_param - (long) offset;
    
    /* Let only the master thread print stuff */
    want_output = (thread_nr == 0);

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
      if (test_verbose (OUTPUT_TRACE))
	{ 
	  mpres_get_z (t, b_1, modulus);
	  outputf (OUTPUT_TRACE, "\n/* pm1_sequence_g */ N = %Zd; "
		   "b_1 = Mod(%Zd, N); /* PARI */\n", modulus->orig_modulus, t);
	  outputf (OUTPUT_TRACE, "/* pm1_sequence_g */ P = %lu; M = %ld; "
		   "m_1 = %Zd; /* PARI */\n", P, M, m_1);
	  outputf (OUTPUT_TRACE, 
		   "/* pm1_sequence_g */ r = b_1^P; /* PARI */\n");
	  outputf (OUTPUT_TRACE, "/* pm1_sequence_g */ x_0 = "
		   "b_1^(2*%ld + (2*m_1 + 1)*P); /* PARI */\n", k_2);
	}
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

  mpz_set_si (t, M);
  mpz_neg (t, t);
  mpz_mul_2exp (t, t, 1UL);
  mpz_add_ui (t, t, 1UL);
  mpres_pow (r[1], r[0], t, modulus);    /* r[1] = r^{2(-M+i)+1}, i = 0 */
  mpz_set_si (t, M);
  mpz_mul (t, t, t);                     /* t = M^2 */
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
  
  mpz_set_si (t, M);
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

  print_elapsed_time (cputime() - timestart, realtime() - realstart);
  
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
		const unsigned long d, mpmod_t modulus_parm, 
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
    unsigned long j;
    unsigned long offset = 0UL, len = d;
    mpmod_t modulus;

    /* Adjust offset and length for this thread */
#ifdef _OPENMP
    {
      const int nr_chunks = omp_get_num_threads();
      const int thread_nr = omp_get_thread_num();
      unsigned long chunklen;
      
      if (thread_nr == 0)
	outputf (OUTPUT_VERBOSE, " using %d threads", nr_chunks);

      chunklen = (len - 1UL) / (unsigned long) nr_chunks + 1UL;
      offset = chunklen * (unsigned long) thread_nr;
      len = MIN(chunklen, len - offset);
    }
#endif
    
    mpmod_copy (modulus, modulus_parm);
    mpres_init (fd[0], modulus);
    mpres_init (fd[1], modulus);
    mpres_init (fd[2], modulus);
    mpz_init (t);
    
    /* We have (n + 1)^2 = n^2 + 2n + 1. For the finite differences we'll 
       need r^{-2}, r^{-(2n+1)}, r^{-n^2}. Init for n = 0. */
    
    /* r^{-2} in fd[0] is constant and could be shared. Computing it 
       separately in each thread has the advantage of putting it in
       local memory. May not make much difference overall */

    mpres_mul (fd[0], invr, invr, modulus); /* fd[0] = r^{-2} */
    mpz_set_ui (t, offset);
    mpz_mul_2exp (t, t, 1UL);
    mpz_add_ui (t, t, 1UL);                 /* t = 2 * offset + 1 */
    mpres_pow (fd[1], invr, t, modulus);    /* fd[1] = r^{-(2*offset+1)} */
    mpz_set_ui (t, offset);
    mpz_mul (t, t, t);                      /* t = offset^2 */
    mpres_pow (fd[2], invr, t, modulus);    /* fd[2] = r^{-offset^2} */
    
    /* Generate the sequence */
    for (j = offset; j < offset + len; j++)
      {
	mpres_mul_z_to_z (t, fd[2], f[j], modulus);
	outputf (OUTPUT_TRACE, 
		 "/* pm1_sequence_h */ h_%lu = %Zd; /* PARI */\n", j, t);
	
	if (h != NULL)
	  mpz_set (h[j], t);
	if (h_ntt != NULL)
	  mpzspv_from_mpzv (h_ntt, j, &t, 1UL, ntt_context);
	
	mpres_mul (fd[2], fd[2], fd[1], modulus); /* fd[2] = r^{-j^2} */
	mpres_mul (fd[1], fd[1], fd[0], modulus); /* fd[1] = r^{-2*j-1} */
      }
    
    mpres_clear (fd[2], modulus);
    mpres_clear (fd[1], modulus);
    mpres_clear (fd[0], modulus);
    mpz_clear (t);
    mpmod_clear (modulus);
  }

  mpres_clear (invr, modulus_parm);

  print_elapsed_time (cputime() - timestart, realtime() - realstart);

  if (test_verbose (OUTPUT_TRACE))
    {
      unsigned long j;
      for (j = 0; j < d; j++)
	outputf (OUTPUT_TRACE, "/* pm1_sequence_h */ h_%lu == "
		   "f_%lu * r^(-%lu^2) /* PARI C */\n", j, j, j);
      outputf (OUTPUT_TRACE, "/* pm1_sequence_h */ h(x) = h_0");
      for (j = 1; j < d; j++)
        outputf (OUTPUT_TRACE, " + h_%lu * (x^%lu + x^(-%lu))", j, j, j);
      outputf (OUTPUT_TRACE, " /* PARI */\n");
    }
}


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
  const spv_size_t log2_l = ceil_log_2 (l);
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
  for (j = 0; j < (int) ntt_context->sp_num; j++)
    {
      const spm_t spm = ntt_context->spm[j];
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
      
      spv_ntt_gfp_dif (tmp[j], log2_l, spm);

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
  spv_size_t log2_len = ceil_log_2 (len);
  
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
	unsigned long i, m;
	
	/* Forward DFT of dft[j] */
	spv_ntt_gfp_dif (spv, log2_len, spm);
	
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
	
	/* Inverse transform of dft[j] */
	spv_ntt_gfp_dit (spv, log2_len, spm);
	
	/* Divide by transform length. FIXME: scale the DCT of h instead */
	spv_mul_sp (spv, spv, spm->sp - (spm->sp - 1) / len, len, spm->sp, 
		    spm->mul_c);
      }
#ifdef _OPENMP
  }
#endif
}

ATTRIBUTE_UNUSED
static void
ntt_print_vec (const char *msg, const spv_t spv, const spv_size_t l)
{
  spv_size_t i;
  printf ("%s [%lu", msg, spv[0]);
  for (i = 1; i < l; i++)
    printf (", %lu", spv[i]);
  printf ("]\n");
}


/* Square the reciprocal Laurent polynomial S(x) of degree 2*n-2.
   S(x) = s_0 + \sum_{i=1}^{n-1} s_i (x^i + x^{-1}).
   S[i] contains the n coefficients s_i, 0 <= i <= n-1.
   R[i] will contain the 2n-1 coefficients r_i, 0 <= i <= 2*n-2, where 
   R(x) = S(x)^2 = r_0 + \sum_{i=1}^{2n-2} r_i (x^i + x^{-1}).
   dft must have power of 2 length len >= 2n.
   The NTT primes must be == 1 (mod 4*len).
*/

#undef TRACE_ntt_sqr_reciprocal
void
ntt_sqr_reciprocal (mpzv_t R, const mpzv_t S, mpzspv_t dft, 
		    const spv_size_t n, const mpzspm_t ntt_context)
{
  const spv_size_t log2_n = ceil_log2 (n);
  const spv_size_t len = ((spv_size_t) 2) << log2_n;
  const spv_size_t log2_len = 1 + log2_n;
#ifdef WANT_ASSERT
  mpz_t S_eval_1, R_eval_1;
#endif
  int j;
  
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

#ifdef TRACE_ntt_sqr_reciprocal
  printf ("ntt_sqr_reciprocal: n %lu, length %lu\n", n, len);
  gmp_printf ("Input polynomial is %Zd", S[0]);
  for (j = 1; (spv_size_t) j < n; j++)
    gmp_printf (" + %Zd * (x^%lu + x^(-%lu))", S[j], j, j);
  printf ("\n");
#endif
  
  ASSERT(ntt_context->max_ntt_size % 3UL == 0UL);
  ASSERT(len % 3UL != 0UL);
  ASSERT(ntt_context->max_ntt_size % len == 0UL);
  /* Fill NTT elements [0 .. n-1] with coefficients */
  mpzspv_from_mpzv (dft, (spv_size_t) 0, S, (spv_size_t) n, ntt_context);

#ifdef _OPENMP
#pragma omp parallel
  {
#pragma omp for
#endif
    for (j = 0UL; j < (int) (ntt_context->sp_num); j++)
      {
        const spm_t spm = ntt_context->spm[j];
        const spv_t spv = dft[j];
        sp_t w1, w2, invlen;
        const sp_t sp = spm->sp, mul_c = spm->mul_c;
        spv_size_t i;

        /* Zero out NTT elements [n .. len-n] */
        spv_set_sp (spv + n, (sp_t) 0, len - 2*n + 1);

#ifdef TRACE_ntt_sqr_reciprocal
        if (j == 0UL)
          {
            printf ("ntt_sqr_reciprocal: NTT vector mod %lu\n", sp);
            ntt_print_vec ("ntt_sqr_reciprocal: before weighting:", spv, len);
          }
#endif

        /* Compute the root for the weight signal, a 3rd primitive root 
           of unity */
        w1 = sp_pow (spm->prim_root, ntt_context->max_ntt_size / 3UL, sp, 
                     mul_c);
        /* Compute iw= 1/w */
        w2 = sp_pow (spm->inv_prim_root, ntt_context->max_ntt_size / 3UL, sp, 
                     mul_c);
#ifdef TRACE_ntt_sqr_reciprocal
        if (j == 0)
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
        for (i = 0; i < n - 2UL; i += 3)
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
        if (i < n - 1UL)
          {
            sp_t t, u;
            t = spv[i + 1];
            u = sp_mul (t, w1, sp, mul_c);
            spv[i + 1] = u;
            spv[len - i - 1] = sp_neg (sp_add (t, u, sp), sp);
          }

#ifdef TRACE_ntt_sqr_reciprocal
      if (j == 0UL)
        ntt_print_vec ("ntt_sqr_reciprocal: after weighting:", spv, len);
#endif

        /* Forward DFT of dft[j] */
        spv_ntt_gfp_dif (spv, log2_len, spm);

#if 0 && defined(TRACE_ntt_sqr_reciprocal)
        if (j == 0UL)
          ntt_print_vec ("ntt_sqr_reciprocal: after forward transform:", 
			 spv, len);
#endif

        /* Square the transformed vector point-wise */
        spv_pwmul (spv, spv, spv, len, sp, mul_c);
      
#if 0 && defined(TRACE_ntt_sqr_reciprocal)
        if (j == 0UL)
          ntt_print_vec ("ntt_sqr_reciprocal: after point-wise squaring:", 
			 spv, len);
#endif

        /* Inverse transform of dft[j] */
        spv_ntt_gfp_dit (spv, log2_len, spm);
      
#if 0 && defined(TRACE_ntt_sqr_reciprocal)
        if (j == 0UL)
          ntt_print_vec ("ntt_sqr_reciprocal: after inverse transform:", 
			 spv, len);
#endif

        /* Un-weight and divide by transform length */
        invlen = sp - (sp - (sp_t) 1) / len; /* weight = 1/len (mod sp) */
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
        
#ifdef TRACE_ntt_sqr_reciprocal
        if (j == 0UL)
          ntt_print_vec ("ntt_sqr_reciprocal: after un-weighting:", spv, len);
#endif

        /* Separate the coefficients of R in the wrapped-around product. */

        /* Set w1 = cuberoot(1)^l where cuberoot(1) is the same primitive
           3rd root of unity we used for the weight signal */
        w1 = sp_pow (spm->prim_root, ntt_context->max_ntt_size / 3UL, sp, 
                     mul_c);
        w1 = sp_pow (w1, len % 3UL, sp, mul_c);
        
        /* Set w2 = 1/(w1 - 1/w1). Incidentally, w2 = 1/sqrt(-3) */
        w2 = sp_inv (w1, sp, mul_c);
        w2 = sp_sub (w1, w2, sp);
        w2 = sp_inv (w2, sp, mul_c);
#ifdef TRACE_ntt_sqr_reciprocal
        if (j == 0UL)
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

#ifdef TRACE_ntt_sqr_reciprocal
        if (j == 0UL)
          ntt_print_vec ("ntt_sqr_reciprocal: after un-wrapping:", spv, len);
#endif
      }
#ifdef _OPENMP
    }
#endif

#if defined(_OPENMP)
#pragma omp parallel if (n > 50)
#endif
  {
    spv_size_t i, offset = 0, chunklen = 2*n - 1;

#if defined(_OPENMP)
    {
      const int nr_chunks = omp_get_num_threads();
      const int thread_nr = omp_get_thread_num();
      
      chunklen = (chunklen - 1) / (spv_size_t) nr_chunks + 1;
      offset = (spv_size_t) thread_nr * chunklen;
      if (2*n - 1 > offset)
        chunklen = MIN(chunklen, (2*n - 1) - offset);
      else
        chunklen = 0UL;
    }
#endif
    
    mpzspv_to_mpzv (dft, offset, R + offset, chunklen, ntt_context);
    for (i = offset; i < offset + chunklen; i++)
      mpz_mod (R[i], R[i], ntt_context->modulus);
  }

#ifdef TRACE_ntt_sqr_reciprocal
  gmp_printf ("ntt_sqr_reciprocal: Output polynomial is %Zd", R[0]);
  for (j = 1; (spv_size_t) j < 2*n - 1; j++)
    gmp_printf (" + %Zd * (x^%lu + x^(-%lu))", R[j], j, j);
  printf ("\n");
#endif

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
      gmp_printf ("Output polynomial is %Zd", R[0]);
      for (j = 1; (spv_size_t) j < 2*n - 1; j++)
	gmp_printf (" + %Zd * (x^%lu + x^(-%lu))", R[j], j, j);
      printf ("\n");
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
ntt_gcd (mpz_t f, mpz_t *product, mpzspv_t ntt, const unsigned long ntt_offset,
	 const listz_t add, const unsigned long len_param, 
	 const mpzspm_t ntt_context, mpmod_t modulus_param)
{
  unsigned long i, j;
  const unsigned long Rlen = MPZSPV_NORMALISE_STRIDE;
  listz_t R;
  unsigned long len = len_param, thread_offset = 0;
  mpres_t tmpres, tmpprod, totalprod;
  mpmod_t modulus;
  long timestart, realstart;
  
  outputf (OUTPUT_VERBOSE, "Computing gcd of coefficients and N");
  timestart = cputime ();
  realstart = realtime ();

  /* All the threads will multiply their partial products to this one. */
  mpres_init (totalprod, modulus_param);
  mpres_set_ui (totalprod, 1UL, modulus_param);

#ifdef _OPENMP
#pragma omp parallel if (len > 100) private(i, j, R, len, thread_offset, tmpres, tmpprod, modulus) shared(totalprod)
  {
    const int nr_chunks = omp_get_num_threads();
    const int thread_nr = omp_get_thread_num();

    len = (len_param - 1) / nr_chunks + 1;
    thread_offset = thread_nr * len;
    ASSERT (len_param >= thread_offset);
    len = MIN(len, len_param - thread_offset);
#pragma omp master
    {
      outputf (OUTPUT_VERBOSE, " using %d threads", nr_chunks);
    }
#endif

    /* Make a private copy of the mpmod_t struct */
    mpmod_copy (modulus, modulus_param);

    MEMORY_TAG;
    R = init_list2 (Rlen, (mpz_size (modulus->orig_modulus) + 2) * 
                           mp_bits_per_limb);
    MEMORY_UNTAG;
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

  if (product != NULL)
    mpres_get_z (*product, totalprod, modulus_param);

  mpres_gcd (f, totalprod, modulus_param);
  mpres_clear (totalprod, modulus_param);

  print_elapsed_time (cputime() - timestart, realtime() - realstart);
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
  long timetotalstart, timestart;

  timetotalstart = cputime ();

  phiP = eulerphi (params->P);
  ASSERT_ALWAYS (phiP == params->s_1 * params->s_2);
  ASSERT_ALWAYS (params->s_1 < params->l);
  nr = params->l - params->s_1; /* Number of points we evaluate */

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
  
  /* First compute X + 1/X */
  mpres_init (mr, modulus);
  mpres_invert (mr, X, modulus);
  mpres_add (mr, mr, X, modulus);
  
  i = poly_from_sets_V (F, mr, S_1, tmp, tmplen, modulus, NULL, NULL);
  ASSERT_ALWAYS(2 * i == params->s_1);
  ASSERT(mpz_cmp_ui (F[i], 1UL) == 0);
  free (S_1);
  S_1 = NULL;
  
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
  pm1_sequence_h (F, NULL, F, mr, params->s_1 / 2 + 1, modulus, NULL); 

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

#if 0 && defined(WANT_ASSERT)

      /* See if R[i] is correct, with a test that works even if i0 != 0 */
      /* More expensive self-test */
      /* alpha = beta*(i0 + l*nr) */
      /* This code is old and probably does not work. */

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
  
  free (S_2);
  free (h);
  clear_list (F, lenF);
  clear_list (g, lenG);
  clear_list (R, lenR);    
  clear_list (tmp, tmplen);

  mpz_clear (mt);
  mpres_clear (mr, modulus);

  if (youpi == ECM_NO_FACTOR_FOUND)
    outputf (OUTPUT_NORMAL, "Step 2 took %ldms\n", 
	     cputime () - timetotalstart);
  
  return youpi;
}


int 
pm1fs2_ntt (mpz_t f, const mpres_t X, mpmod_t modulus, 
	const faststage2_param_t *params)
{
  unsigned long nr;
  unsigned long i, l, lenF;
  sets_long_t *S_1; /* This is stored as a set of sets (arithmetic 
                       progressions of prime length */
  set_long_t *S_2; /* This is stored as a regular set */
  listz_t F;   /* Polynomial F has roots X^{k_1} for k_1 \in S_1, so has 
		  degree s_1. It is symmetric, so has only s_1 / 2 + 1 
		  distinct coefficients. The sequence h_j will be stored in 
		  the same memory and won't be a monic polynomial, so the 
		  leading 1 monomial of F will be stored explicitly. Hence we 
		  need s_1 / 2 + 1 entries. */
  mpzspm_t ntt_context;
  mpzspv_t g_ntt, h_ntt;
  mpz_t mt;   /* All-purpose temp mpz_t */
  mpz_t product; /* Product of each multi-point evaluation */
  mpz_t *product_ptr = NULL;
  mpres_t tmpres; /* All-purpose temp mpres_t */
  int youpi = ECM_NO_FACTOR_FOUND;
  long timetotalstart, timestart, realstart;

  timetotalstart = cputime ();

  ASSERT_ALWAYS (eulerphi (params->P) == params->s_1 * params->s_2);
  ASSERT_ALWAYS (params->s_1 < params->l);
  nr = params->l - params->s_1; /* Number of points we evaluate */

  /* Prepare NTT for computing the h sequence, its DCT-I, and the convolution 
     with g. We need NTT of transform length l. We do it here at the start
     of stage 2 so that in case of a "not enough primes" condition, 
     we don't have to wait until after F is built to get the error. */

  ntt_context = mpzspm_init (params->l, modulus->orig_modulus);
  if (ntt_context == NULL)
    return ECM_ERROR;

  print_CRT_primes (OUTPUT_DEVVERBOSE, "CRT modulus for evaluation = ", 
		    ntt_context);

  if (make_S_1_S_2 (&S_1, &S_2, params) == ECM_ERROR)
      return ECM_ERROR;

  /* Allocate all the memory we'll need for building f */
  mpz_init (mt);
  mpres_init (tmpres, modulus);
  lenF = params->s_1 / 2 + 1 + 1; /* Another +1 because poly_from_sets_V stores
				     the leading 1 monomial for each factor */
  F = init_list2 (lenF, (unsigned int) abs (modulus->bits));

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


  /* First compute X + 1/X */
  mpres_invert (tmpres, X, modulus);
  mpres_add (tmpres, tmpres, X, modulus);

  if (build_F_ntt (F, tmpres, S_1, params, modulus) == ECM_ERROR)
    {
      free (S_1);
      free (S_2);
      mpz_clear (mt);
      mpres_clear (tmpres, modulus);
      mpzspm_clear (ntt_context);
      clear_list (F, lenF);
      return ECM_ERROR;
    }

  free (S_1);
  S_1 = NULL;
  
  h_ntt = mpzspv_init (params->l / 2 + 1, ntt_context);

  mpz_set_ui (mt, params->P);
  mpres_pow (tmpres, X, mt, modulus); /* tmpres = X^P */
  pm1_sequence_h (NULL, h_ntt, F, tmpres, params->s_1 / 2 + 1, modulus, 
		  ntt_context);

  clear_list (F, lenF);
  g_ntt = mpzspv_init (params->l, ntt_context);

  /* Compute the DCT-I of h */
  outputf (OUTPUT_VERBOSE, "Computing DCT-I of h");
  timestart = cputime ();
  realstart = realtime ();
  
  ntt_spv_to_dct (h_ntt, h_ntt, params->s_1 / 2 + 1, params->l / 2 + 1, 
                  g_ntt, ntt_context);
  print_elapsed_time (cputime() - timestart, realtime() - realstart);
  
  if (test_verbose (OUTPUT_RESVERBOSE))
    {
      mpz_init (product);
      product_ptr = &product;
    }

  for (l = 0; l < params->s_2; l++)
    {
      const unsigned long M = params->l - 1L - params->s_1 / 2L;

      /* Compute the coefficients of the polynomial g(x) */
      pm1_sequence_g (NULL, g_ntt, X, params->P, M, params->l, 
		      params->m_1, S_2->elem[l], modulus, ntt_context);

      /* Do the convolution */
      outputf (OUTPUT_VERBOSE, "Computing g*h");
      timestart = cputime ();
      realstart = realtime ();
      ntt_mul_by_dct (g_ntt, h_ntt, params->l, ntt_context);
      print_elapsed_time (cputime() - timestart, realtime() - realstart);
      
      /* Compute GCD of N and coefficients of product polynomial */
      ntt_gcd (mt, product_ptr, g_ntt, params->s_1 / 2, NULL, nr, ntt_context, 
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
  mpzspv_clear (g_ntt, ntt_context);
  mpzspv_clear (h_ntt, ntt_context);
  mpzspm_clear (ntt_context);
  mpres_clear (tmpres, modulus);
  mpz_clear (mt);
  free (S_2);

  outputf (OUTPUT_NORMAL, "Step 2 took %ldms\n", cputime () - timetotalstart);
  
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


/* Raise (a0 + a1*sqrt(Delta)) to the power e which is a signed long int. 
   (a0 + a1*sqrt(Delta)) is assumed to have norm 1, i.e. 
   a0^2 - a1^2*Delta == 1. The result is (r0 * r1*sqrt(Delta)). 
   a0, a1, r0 and r1 must not overlap */

static void 
gfp_ext_pow_norm1_sl (mpres_t r0, mpres_t r1, const mpres_t a0, 
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
      if (pari)
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
  long timestart, realstart;
  int want_output = 1;

  outputf (OUTPUT_VERBOSE, "Computing %s%s%s", 
	   (want_x) ? "g_x" : "", 
	   (want_x && want_y) ? " and " : "",
	   (want_y) ? "g_y" : "");
  timestart = cputime ();
  realstart = realtime ();

#ifdef _OPENMP
#pragma omp parallel if (l > 100) private(r_x, r_y, x0_x, x0_y, v2, r1_x, r1_y, r2_x, r2_y, v, tmp, mt, modulus, i, l, offset, M, want_output)
  {
    /* When multi-threading, we adjust the parameters for each thread */

    const int nr_chunks = omp_get_num_threads();
    const int thread_nr = omp_get_thread_num();

    l = (l_param - 1) / nr_chunks + 1;
    offset = thread_nr * l;
    ASSERT_ALWAYS (l_param >= offset);
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
    gfp_ext_pow_norm1_sl (r_x, r_y, b1_x, b1_y, P, Delta, modulus, 
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
    mpmod_clear (modulus);
#ifdef _OPENMP
  }
#endif
  
  print_elapsed_time (cputime() - timestart, realtime() - realstart);

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
		const long k_param, const unsigned long l_param, 
		const unsigned long P, const mpres_t Delta, 
		mpmod_t modulus_param, const mpzspm_t ntt_context)
{
  unsigned long i;
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
	       "Delta = %Zd; w = quadgen (4*Delta); k = %ld; P = %lu; "
	       "/* PARI */\n", modulus_param->orig_modulus, t, k_param, P);
      outputf (OUTPUT_TRACE, "/* pp1_sequence_h */ b_1 = ");
      gfp_ext_print (b1_x, b1_y, modulus_param, OUTPUT_TRACE);
      outputf (OUTPUT_TRACE, "; r = b_1^P; rn = b_1^(-P); /* PARI */\n");
      for (i = 0; i < l_param; i++)
	outputf (OUTPUT_TRACE, 
		 "/* pp1_sequence_h */ f_%lu = %Zd; /* PARI */\n", i, f[i]);
      mpz_clear (t);
    }

#ifdef _OPENMP
#pragma omp parallel if (l_param > 100) private(i)
#endif
  {
    const size_t tmplen = 2;
    mpres_t s_x[3], s_y[3], s2_x[2], s2_y[2], v[2], V2, rn_x, rn_y, 
      tmp[tmplen];
    mpmod_t modulus; /* Thread-local copy of modulus_param */
    mpz_t mt;
    unsigned long l = l_param, offset = 0;
    long k = k_param;

#ifdef _OPENMP
    /* When multi-threading, we adjust the parameters for each thread */

    const int nr_chunks = omp_get_num_threads();
    const int thread_nr = omp_get_thread_num();

    l = (l_param - 1) / nr_chunks + 1;
    offset = thread_nr * l;
    ASSERT_ALWAYS (l_param >= offset);
    l = MIN(l, l_param - offset);

    if (thread_nr == 0)
      outputf (OUTPUT_VERBOSE, " using %d threads", nr_chunks);
    outputf (OUTPUT_TRACE, "\n");
#endif

    /* Each thread computes r[i + offset] = b1^(-P*(k+i+offset)^2) * f_i 
       for i = 0, 1, ..., l-1, where l is the adjusted length of each thread */

    /* Test that k+offset does not overflow */
    ASSERT_ALWAYS (offset <= (unsigned long) LONG_MAX && 
		   k <= LONG_MAX - (long) offset);
    k += (long) offset;

    mpz_init (mt);
    /* Make thread-local copy of modulus */
    mpmod_copy (modulus, modulus_param);

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
    for (i = 0; i < (unsigned long) tmplen; i++)
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
	  outputf (OUTPUT_TRACE, "/* pp1_sequence_h */ rn^(%ld^2) == ", k);
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
	      outputf (OUTPUT_TRACE, "/* pp1_sequence_h */ rn^(%ld^2) == ", 
		       k + 1);
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
		   "1/r^(2*%ld+1) == %Zd /* PARI C */\n", k, k, mt);
	  mpres_get_z (mt, v[1], modulus);
	  outputf (OUTPUT_TRACE, "/* pp1_sequence_h */ r^(2*%ld+3) + "
		   "1/r^(2*%ld+3) == %Zd /* PARI C */\n", k, k, mt);
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

  print_elapsed_time (cputime() - timestart, realtime() - realstart);

  if (h_x != NULL && h_y != NULL && test_verbose (OUTPUT_TRACE))
    {
      for (i = 0; i < l_param; i++)
	gmp_printf ("/* pp1_sequence_h */ (rn^((k+%lu)^2) * f_%lu) == "
		    "(%Zd + Mod(%Zd / Delta, N) * w) /* PARI C */\n", 
		    i, i, h_x[i], h_y[i]);
    }
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
  const unsigned long tmpreslen = 2UL;
  mpres_t b1_x, b1_y, Delta, tmpres[tmpreslen];
  mpz_t mt;   /* All-purpose temp mpz_t */
  int youpi = ECM_NO_FACTOR_FOUND;
  long timetotalstart, timestart;

  timetotalstart = cputime ();

  ASSERT_ALWAYS (eulerphi (params->P) == params->s_1 * params->s_2);
  ASSERT_ALWAYS (params->s_1 < params->l);
  nr = params->l - params->s_1; /* Number of points we evaluate */

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
  i = poly_from_sets_V (F, X, S_1, tmp, tmplen, modulus, NULL, NULL);
  ASSERT_ALWAYS(2 * i == params->s_1);
  ASSERT(mpz_cmp_ui (F[i], 1UL) == 0);
  free (S_1);
  S_1 = NULL;
  
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
      const long M = params->l - 1 - params->s_1 / 2;
      pp1_sequence_g (g_x, g_y, NULL, NULL, b1_x, b1_y, params->P, 
		      Delta, M, params->l, params->m_1, S_2->elem[l], 
		      modulus, NULL);
      
      /* Do the two convolution products */
      outputf (OUTPUT_VERBOSE, "TMulGen of g_x and h_x");
      timestart = cputime ();
      if (TMulGen (R_x, nr - 1, h_x, params->s_1, g_x, params->l - 1, tmp,
		   modulus->orig_modulus) < 0)
	{
	  outputf (OUTPUT_ERROR, "TMulGen returned error code (probably out "
		   "of memory)\n");
	  youpi == ECM_ERROR;
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
	  youpi == ECM_ERROR;
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
  free (S_2);
 
  outputf (OUTPUT_NORMAL, "Step 2 took %ldms\n", cputime () - timetotalstart);

  return youpi;
}


int 
pp1fs2_ntt (mpz_t f, const mpres_t X, mpmod_t modulus,
	    const faststage2_param_t *params, const int twopass)
{
  unsigned long nr;
  unsigned long i, l, lenF;
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
  mpzspm_t ntt_context;
  mpzspv_t g_x_ntt, g_y_ntt, h_x_ntt, h_y_ntt;
  mpres_t b1_x, b1_y, Delta;
  mpz_t mt;   /* All-purpose temp mpz_t */
  mpz_t product;
  mpz_t *product_ptr = NULL;
  int youpi = ECM_NO_FACTOR_FOUND;
  long timetotalstart, timestart, realstart;

  timetotalstart = cputime ();

  ASSERT_ALWAYS (eulerphi (params->P) == params->s_1 * params->s_2);
  ASSERT_ALWAYS (params->s_1 < params->l);
  nr = params->l - params->s_1; /* Number of points we evaluate */

  if (make_S_1_S_2 (&S_1, &S_2, params) == ECM_ERROR)
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
      mpz_clear (mt);
      free (S_1);
      S_1 = NULL;
      free (S_2);
      S_2 = NULL;
      return ECM_ERROR;
    }

  print_CRT_primes (OUTPUT_DEVVERBOSE, "CRT modulus for evaluation = ", 
		    ntt_context);

  /* Allocate memory for F with correct amount of space for each mpz_t */
  lenF = params->s_1 / 2 + 1 + 1; /* Another +1 because poly_from_sets_V stores
				     the leading 1 monomial for each factor */
  MEMORY_TAG;
  F = init_list2 (lenF, (unsigned int) abs (modulus->bits) + mp_bits_per_limb);
  MEMORY_UNTAG;
  
  /* Build F */
  if (build_F_ntt (F, X, S_1, params, modulus) == ECM_ERROR)
    {
      free (S_1);
      free (S_2);
      mpz_clear (mt);
      mpzspm_clear (ntt_context);
      clear_list (F, lenF);
      return ECM_ERROR;
    }

  free (S_1);
  S_1 = NULL;
  
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
                          mp_bits_per_limb);
      MEMORY_UNTAG;
    }
  else
    g_y_ntt = mpzspv_init (params->l, ntt_context);
  
  /* Compute DCT-I of h_x and h_y */
  outputf (OUTPUT_VERBOSE, "Computing DCT-I of h_x");
  timestart = cputime ();
  realstart = realtime ();
  ntt_spv_to_dct (h_x_ntt, h_x_ntt, params->s_1 / 2 + 1, params->l / 2 + 1,
		  g_x_ntt, ntt_context);
  print_elapsed_time (cputime() - timestart, realtime() - realstart);

  outputf (OUTPUT_VERBOSE, "Computing DCT-I of h_y");
  timestart = cputime ();
  realstart = realtime ();
  ntt_spv_to_dct (h_y_ntt, h_y_ntt, params->s_1 / 2 + 1, params->l / 2 + 1,
		  g_x_ntt, ntt_context);
  print_elapsed_time (cputime() - timestart, realtime() - realstart);

  if (test_verbose (OUTPUT_RESVERBOSE))
    {
      mpz_init (product);
      product_ptr = &product;
    }

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
	  realstart = realtime ();
	  ntt_mul_by_dct (g_x_ntt, h_x_ntt, params->l, ntt_context);
	  /* Store the product coefficients we want in R */
	  mpzspv_to_mpzv (g_x_ntt, params->s_1 / 2, R, nr, ntt_context);
	  print_elapsed_time (cputime() - timestart, realtime() - realstart);

	  /* Compute g_y sequence */
	  pp1_sequence_g (NULL, NULL, NULL, g_y_ntt, b1_x, b1_y, params->P, 
			  Delta, M, params->l, params->m_1, S_2->elem[l], 
			  modulus, ntt_context);
	  
	  /* Do the convolution product of g_y * (Delta * h_y) */
	  outputf (OUTPUT_VERBOSE, "Computing g_y*h_y");
	  timestart = cputime ();
	  realstart = realtime ();
	  ntt_mul_by_dct (g_y_ntt, h_y_ntt, params->l, ntt_context);
	  print_elapsed_time (cputime() - timestart, realtime() - realstart);
	  
	  /* Compute product of sum of coefficients and gcd with N */
	  ntt_gcd (mt, product_ptr, g_y_ntt, params->s_1 / 2, R, nr, 
		   ntt_context, modulus);
	}
      else
	{
	  pp1_sequence_g (NULL, NULL, g_x_ntt, g_y_ntt, b1_x, b1_y, params->P, 
			  Delta, M, params->l, params->m_1, S_2->elem[l], 
			  modulus, ntt_context);

	  outputf (OUTPUT_VERBOSE, "Computing forward NTT of g_x");
	  timestart = cputime ();
	  realstart = realtime ();
	  mpzspv_to_ntt (g_x_ntt, (spv_size_t) 0, (spv_size_t) params->l, 
	                 (spv_size_t) params->l, 0, ntt_context);
	  print_elapsed_time (cputime() - timestart, realtime() - realstart);
	  
	  outputf (OUTPUT_VERBOSE, "Computing point-wise product of g_x and h_x");
	  timestart = cputime ();
	  realstart = realtime ();
	  ntt_dft_mul_dct (g_x_ntt, h_x_ntt, params->l, ntt_context);
	  print_elapsed_time (cputime() - timestart, realtime() - realstart);
	  
	  outputf (OUTPUT_VERBOSE, "Computing forward NTT of g_y");
	  timestart = cputime ();
	  realstart = realtime ();
	  mpzspv_to_ntt (g_y_ntt, (spv_size_t) 0, (spv_size_t) params->l, 
	                 (spv_size_t) params->l, 0, ntt_context);
	  print_elapsed_time (cputime() - timestart, realtime() - realstart);
	  
	  outputf (OUTPUT_VERBOSE, "Computing point-wise product of g_y and h_y");
	  timestart = cputime ();
	  realstart = realtime ();
	  ntt_dft_mul_dct (g_y_ntt, h_y_ntt, params->l, ntt_context);
	  print_elapsed_time (cputime() - timestart, realtime() - realstart);
	  
	  outputf (OUTPUT_VERBOSE, "Adding and computing inverse NTT of sum");
	  timestart = cputime ();
	  realstart = realtime ();
	  mpzspv_add (g_x_ntt, (spv_size_t) 0, g_x_ntt, (spv_size_t) 0, 
	              g_y_ntt, (spv_size_t) 0, params->l, ntt_context);
	  mpzspv_from_ntt (g_x_ntt, (spv_size_t) 0, params->l, (spv_size_t) 0,
	                   ntt_context);
	  print_elapsed_time (cputime() - timestart, realtime() - realstart);
	  
	  ntt_gcd (mt, product_ptr, g_x_ntt, params->s_1 / 2, NULL, nr, 
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
  free (S_2);
 
  outputf (OUTPUT_NORMAL, "Step 2 took %ldms\n", cputime () - timetotalstart);

  return youpi;
}
