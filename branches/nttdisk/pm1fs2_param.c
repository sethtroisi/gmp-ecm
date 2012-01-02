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

/* Some useful PARI functions:
   sumset(a,b) = {local(i, j, l); l = listcreate (length(a) * length(b)); for (i = 1, length(a), for (j = 1, length(b), listput(l, a[i] + b[j]))); listsort (l, 1); l}
*/

const uint64_t Pvalues[] = {
3ULL, 5ULL, 9ULL, 15ULL, 21ULL, 17ULL, 27ULL, 33ULL, 45ULL, 
51ULL, 63ULL, 75ULL, 105ULL, 99ULL, 135ULL, 165ULL, 195ULL, 189ULL, 
231ULL, 255ULL, 315ULL, 345ULL, 357ULL, 375ULL, 405ULL, 435ULL, 
525ULL, 585ULL, 615ULL, 735ULL, 765ULL, 825ULL, 945ULL, 1155ULL, 
1065ULL, 1365ULL, 1305ULL, 1335ULL, 1575ULL, 1785ULL, 1995ULL, 
2145ULL, 2205ULL, 2415ULL, 2625ULL, 2805ULL, 3045ULL, 3465ULL, 
3675ULL, 4095ULL, 4305ULL, 4515ULL, 4725ULL, 4785ULL, 5355ULL, 
5775ULL, 5985ULL, 5865ULL, 6825ULL, 7245ULL, 8085ULL, 8925ULL, 
9555ULL, 10395ULL, 10725ULL, 11025ULL, 12285ULL, 12705ULL, 
15015ULL, 14175ULL, 15225ULL, 16065ULL, 17325ULL, 19635ULL, 
21945ULL, 23205ULL, 24255ULL, 25935ULL, 26775ULL, 28875ULL, 
31395ULL, 33495ULL, 35805ULL, 36465ULL, 38115ULL, 39585ULL, 
40425ULL, 45045ULL, 45885ULL, 49665ULL, 51765ULL, 58905ULL, 
65835ULL, 69615ULL, 75075ULL, 77805ULL, 82005ULL, 84315ULL, 
86625ULL, 88935ULL, 94185ULL, 98175ULL, 105105ULL, 109725ULL, 
116025ULL, 118755ULL, 121275ULL, 135135ULL, 137445ULL, 137655ULL, 
144375ULL, 153615ULL, 165165ULL, 167475ULL, 176715ULL, 179025ULL, 
185955ULL, 197505ULL, 208845ULL, 215985ULL, 225225ULL, 255255ULL, 
250635ULL, 285285ULL, 277095ULL, 294525ULL, 315315ULL, 345345ULL, 
373065ULL, 368445ULL, 405405ULL, 435435ULL, 451605ULL, 465465ULL, 
454545ULL, 504735ULL, 525525ULL, 555555ULL, 569415ULL, 596505ULL, 
645645ULL, 647955ULL, 672945ULL, 687225ULL, 765765ULL, 770385ULL, 
805035ULL, 855855ULL, 858585ULL, 915915ULL, 945945ULL, 962115ULL, 
1036035ULL, 1066065ULL, 1119195ULL, 1156155ULL, 1276275ULL, 
1306305ULL, 1354815ULL, 1426425ULL, 1456455ULL, 1514205ULL, 
1576575ULL, 1666665ULL, 1726725ULL, 1786785ULL, 1789515ULL, 
1865325ULL, 1996995ULL, 1983135ULL, 2177175ULL, 2297295ULL, 
2327325ULL, 2417415ULL, 2567565ULL, 2611455ULL, 2807805ULL, 2847075ULL, 
2878785ULL, 3048045ULL, 3161235ULL, 3258255ULL, 3357585ULL, 3401475ULL, 
3533145ULL, 3828825ULL, 3918915ULL, 3985905ULL, 4279275ULL, 4849845ULL, 
4789785ULL, 4967655ULL, 5180175ULL, 5360355ULL, 5870865ULL, 5990985ULL, 
6561555ULL, 6531525ULL, 6891885ULL, 7402395ULL, 7912905ULL, 8273265ULL, 
8580495ULL, 8843835ULL, 9444435ULL, 10015005ULL, 10465455ULL, 10705695ULL, 
10885875ULL, 11696685ULL, 12267255ULL, 12507495ULL, 12785955ULL, 
13498485ULL, 14549535ULL, 14849835ULL, 15570555ULL, 16111095ULL, 
16291275ULL, 17612595ULL, 18123105ULL, 18633615ULL, 19684665ULL, 
20255235ULL, 20825805ULL, 22207185ULL, 22717695ULL, 24249225ULL, 
24819795ULL, 25741485ULL, 26531505ULL, 28333305ULL, 29354325ULL, 
30045015ULL, 31396365ULL, 32807775ULL, 33948915ULL, 33528495ULL, 
34879845ULL, 37011975ULL, 37522485ULL, 39564525ULL, 41096055ULL, 
43648605ULL, 44219175ULL, 45930885ULL, 47222175ULL, 48333285ULL, 
50075025ULL, 51816765ULL, 52777725ULL, 55390335ULL, 55547415ULL, 
59053995ULL, 60063465ULL, 61906845ULL, 64579515ULL, 66621555ULL, 
67492425ULL, 70105035ULL, 73258185ULL, 74939865ULL, 77224455ULL, 
79594515ULL, 81876795ULL, 84999915ULL, 88062975ULL, 91005915ULL, 
94189095ULL, 98423325ULL, 101846745ULL, 111546435ULL, 111035925ULL, 
115120005ULL, 121246125ULL, 124098975ULL, 130945815ULL, 140645505ULL, 
150345195ULL, 150225075ULL, 155450295ULL, 158333175ULL, 170255085ULL, 
179444265ULL, 190285095ULL, 198843645ULL, 203408205ULL, 206831625ULL, 
217222005ULL, 229474245ULL, 240705465ULL, 252447195ULL, 254999745ULL, 
269023755ULL, 282146865ULL, 287672385ULL, 294076965ULL, 306110805ULL, 
318302985ULL, 334639305ULL, 344338995ULL, 354038685ULL, 363738375ULL, 
373438065ULL, 387221835ULL, 400254855ULL, 421936515ULL, 431636205ULL, 
451035585ULL, 453888435ULL, 470434965ULL, 480134655ULL, 510765255ULL, 
522506985ULL, 557732175ULL, 570855285ULL, 596530935ULL, 610224615ULL, 
627912285ULL, 654729075ULL, 703227525ULL, 722116395ULL, 751725975ULL, 
780825045ULL, 790524735ULL, 821665845ULL, 851275425ULL, 863017155ULL, 
909984075ULL, 936020085ULL, 984518535ULL, 1017041025ULL, 1052416365ULL,
1086110025ULL, 1110614505ULL, 1147371225ULL, 1191785595ULL, 1213887675ULL, 
1265809545ULL, 1282356075ULL, 1331995665ULL, 1391905515ULL, 1450103655ULL, 
1479202725ULL, 1547100555ULL, 1555088535ULL, 1673196525ULL, 1712565855ULL, 
1767130365ULL, 1830673845ULL, 1883166285ULL, 1954487535ULL, 2001964965ULL, 
2119382265ULL, 2187280095ULL, 2255177925ULL, 2342475135ULL, 2390973585ULL, 
2421213795ULL, 2555868315ULL, 2672264595ULL, 2788660875ULL, 2856558705ULL, 
2953555605ULL, 3050552505ULL, 3234846615ULL, 3457939485ULL, 3516137625ULL, 
3681032355ULL, 3758629875ULL, 3904125225ULL, 4127218095ULL, 4360010655ULL, 
4573403835ULL, 4796496705ULL, 4844995155ULL, 5019589575ULL, 5203883685ULL, 
5262081825ULL, 5465775315ULL, 5766465705ULL, 5898837945ULL, 6164152995ULL, 
6358146795ULL, 6411780375ULL, 6804332535ULL, 6980458485ULL, 7172920755ULL, 
7473611145ULL, 7716103395ULL, 7968295335ULL, 8182259085ULL, 8342499165ULL, 
8812168365ULL, 9023519505ULL, 9704539845ULL, 9927632715ULL, 10373818455ULL, 
10439434005ULL, 10820004195ULL, 11043097065ULL, 11489282805ULL, 11877270405ULL, 
12381654285ULL, 12604747155ULL, 13080031965ULL, 13274025765ULL, 13642613985ULL,
14389490115ULL, 14583483915ULL, 15058768725ULL, 15611651055ULL, 16174233075ULL,
16397325945ULL, 17289697425ULL, 17735883165ULL, 18143270145ULL, 18381678315ULL,
19074440385ULL, 19559424885ULL, 20636090475ULL, 20941375455ULL, 21800053275ULL,
22643926305ULL, 23148310185ULL, 24205576395ULL, 24546777255ULL, 25544133615ULL,
26389538175ULL, 26863291455ULL, 27813861075ULL, 29113619535ULL, 29494189725ULL,
30520074585ULL, 30684969315ULL, 31790733975ULL, 33575476935ULL, 34467848415ULL,
35202742575ULL, 36427185795ULL, 38037334335ULL, 39240095895ULL, 40365259935ULL,
42053005995ULL, 43168470345ULL, 44953213305ULL, 45845584785ULL, 48522699225ULL,
50307442185ULL, 51869092275ULL, 53653835235ULL, 54546206715ULL, 56680138515ULL,
58784971245ULL, 59386352025ULL, 61908271425ULL, 63431122755ULL, 65700850215ULL,
67931778915ULL, 70162707615ULL, 72616729185ULL, 74120181135ULL, 75740029365ULL,
78417143805ULL, 80871165375ULL, 82840202445ULL, 86448487125ULL, 88466022645ULL, 
91133437395ULL, 92918180355ULL, 100280245065ULL, 100726430805ULL, 
102811864155ULL, 106749938295ULL, 109000266375ULL, 113219631525ULL, 
119689324755ULL, 121027881975ULL, 127943760945ULL, 132628711215ULL, 
134859639915ULL, 141775518885ULL, 148691397855ULL, 150922326555ULL, 
155607276825ULL, 161320394235ULL, 164977177365ULL, 171446870595ULL, 
177470378085ULL, 183270792705ULL,
};

/* All the prime factors that can appear in eulerphi(P) */
const unsigned long phiPfactors[] = {2UL, 3UL, 5UL, 7UL, 11UL, 13UL, 
				     17UL, 19UL};

/* returns Euler's totient phi function */
static uint64_t
eulerphi64 (uint64_t n)
{
  uint64_t phi = 1, p;

  for (p = 2UL; p * p <= n; p += 2)
    {
      if (n % p == 0)
	{
	  phi *= p - 1;
	  n /= p;
	  while (n % p == 0)
	    {
	      phi *= p;
	      n /= p;
	    }
	}

      if (p == 2UL)
	p--;
    }

  /* now n is prime or 1 */
  return (n == 1) ? phi : phi * (n - 1);
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
	 estimate precisely, so let's go with the fudge factor of 3 here */
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
      lmax = 1UL << ceil_log2 (memory / n / 3);
      return lmax;
    }
  else
    {
      size_t lmax, n;
      
      n = mpz_size (modulus) * sizeof (mp_limb_t) + sizeof (mpz_t);

      /* Guess an initial value of lmax for list_mul_mem (lmax / 2) */
      /* memory = n * 25/4 * lmax + lmax / 2 * sizeof (mpz_t); */
      /* Fudge factor of 3 for TMulKS as above */
      lmax = memory / (3 * 25 * n / 4 + 3 * sizeof (mpz_t) / 2);
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
      return 1UL << (ceil_log2 (n / 2)); /* Rounded down to power of 2 */
    }
  else
    {
      return memory / 5 / (m * 8 + sizeof (mpz_t)) * 2;
    }
}


/* Test if for given P, nr, B2min and B2 we can choose an m_1 so that the 
   stage 2 interval [B2min, B2] is covered. The effective B2min and B2
   are stored in effB2min and effB2 */

static int
test_P (const mpz_t B2min, const mpz_t B2, mpz_t m_1, const uint64_t P, 
	const unsigned long nr, mpz_t effB2min, mpz_t effB2)
{
  mpz_t m, Ptmp;
  /* We need B2min >= 2 * max(S_1 + S_2) + (2*m_1 - 1)*P + 1, or
     B2min - 2 * max(S_1 + S_2) - 1 >= (2*m_1)*P - P, or
     (B2min - 2*max(S_1 + S_2) + P - 1)/(2P) >= m_1
     Choose m_1 accordingly */
  
  mpz_init (m);
  mpz_init (Ptmp);
  sets_max (m, P);
  mpz_mul_2exp (m, m, 1UL); /* m = 2*max(S_1 + S_2) */

  mpz_sub (m_1, B2min, m);
  mpz_sub_ui (m_1, m_1, 1UL); /* m_1 = B2min - 2*max(S_1 + S_2) - 1 */
  mpz_set_uint64(Ptmp, P);
  mpz_add (m_1, m_1, Ptmp);
  mpz_fdiv_q_2exp (m_1, m_1, 1UL);
  mpz_fdiv_q (m_1, m_1, Ptmp);    /* 2UL*P may overflow */
  
  /* Compute effB2min = 2 * max(S_1 + S_2) + (2*(m_1 - 1) + 1)*P + 1 */
  
  mpz_mul_2exp (effB2min, m_1, 1UL);
  mpz_sub_ui (effB2min, effB2min, 1UL);
  mpz_mul (effB2min, effB2min, Ptmp);
  mpz_add (effB2min, effB2min, m);
  mpz_add_ui (effB2min, effB2min, 1UL);
  ASSERT_ALWAYS (mpz_cmp (effB2min, B2min) <= 0);

  /* Compute the smallest value coprime to P at the high end of the stage 2
     interval that will not be covered: 
     2*(min(S_1 + S_2)) + (2*(m_1 + nr) + 1)*P. 
     We assume min(S_1 + S_2) = -max(S_1 + S_2) */
  mpz_add_ui (effB2, m_1, nr);
  mpz_mul_2exp (effB2, effB2, 1UL);
  mpz_add_ui (effB2, effB2, 1UL);
  mpz_mul (effB2, effB2, Ptmp);
  mpz_sub (effB2, effB2, m);

  /* The effective B2 values is that value, minus 1 */
  mpz_sub_ui (effB2, effB2, 1UL);

  mpz_clear (m);
  mpz_clear (Ptmp);
  return (mpz_cmp (B2, effB2) <= 0);
}


static void
factor_phiP (int *exponents, const uint64_t phiP)
{
    const int nrprimes = sizeof (phiPfactors) / sizeof (unsigned long);
    uint64_t cofactor = phiP;
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

static uint64_t
absdiff_ul (uint64_t a, uint64_t b)
{
    return (a > b) ? a - b : b - a;
}

/* Choose s_1 so that s_1 * s_2 = phiP, s_1 is positive and even, 
   s_2 >= min_s2 and s_2 is minimal and abs(s_1 - l) is minimal 
   under those conditions. If use_ntt == 1, we require s_1 < l.
   Returns 0 if no such choice is possible */

static unsigned long 
choose_s_1 (const uint64_t phiP, const uint64_t min_s2,
	    const unsigned long l, const int use_ntt)
{
  const int nrprimes = sizeof (phiPfactors) / sizeof (unsigned long);
  /* Using [nrprimes] here makes the compiler complain about variable-sized
     arrays */
  int phiPexponents[sizeof (phiPfactors) / sizeof (unsigned long)], 
    exponents[sizeof (phiPfactors) / sizeof (unsigned long)];
  uint64_t s_1 = 0UL, s_2 = 0UL, trys_1;
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
      printf ("choose_s_1: Trying trys_1 = %" PRId64 "\n", trys_1);
#endif
      /* See if it satisfies all the required conditions and is an 
	 improvement over the previous choice */
      if (phiP / trys_1 >= min_s2 && 
	  (s_2 == 0 || phiP / trys_1 < s_2) && 
	  absdiff_ul (trys_1, l) < absdiff_ul (s_1, l) &&
	  (use_ntt == 0 || trys_1 < l))
      {
#if 0
	  printf ("choose_s_1: New best s_1 for "
	          "phiP = %" PRId64 ", min_s2 = %" PRId64 
		  ", l = %lu : %" PRId64\n", phiP, min_s2, l, trys_1);
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


/* Approximate cost of stage 2. Cost with and without ntt are not 
   comparable. We have l > s_1 and s_1 * s_2 = eulerphi(P), hence
   s_2*l > eulerphi(P) and so cost (s_2, l) > eulerphi(P) for all P */
static uint64_t
est_cost (const uint64_t s_2, const unsigned long l, const int use_ntt,
          const int method)
{
  if (method == ECM_PM1)
    {
      /* The time for building f, h and DCT-I of h seems to be about 
         7/6 of the time of computing g, h*g and gcd with NTT, and 
         3/2 of the time of computing g, h*g and gcd without NTT */

      if (use_ntt)
        return (7 * l) / 6 + s_2 * l;
      else
        return (3 * l) / 2 + s_2 * l;
    }
  else if (method == ECM_PP1)
    {
      /* Building f is the same, building h and its forward transform is
         twice about as expensive as for P-1. Each multi-point evaluation
         is twice as expensive as for P-1.
         FIXME: The estimate for NTT assumes the "one-pass" variant, in 
         "two-pass" the multipoint evaluations are slower, so the optimum 
         shifts towards smaller s_2 some more */
      if (use_ntt)
        return (4 * l) / 5 + s_2 * l;
      else
        return (3 * l) / 4 + s_2 * l;
    }
  else
    abort (); /* Invalid value for method */
}

/* Choose P so that a stage 2 range from B2min to B2 can be covered with
   multipoint evaluations, each using a convolution of length at most lmax. 
   The parameters for stage 2 are stored in finalparams, the final effective
   B2min and B2 values in final_B2min and final_B2, respecively. Each of these
   may be NULL, in which case the value is not stored. It is permissible
   to let B2min and final_B2min, or B2 and final_B2 point at the same mpz_t. */

long
choose_P (const mpz_t B2min, const mpz_t B2, const unsigned long lmax,
	  const uint64_t min_s2, faststage2_param_t *finalparams, 
	  mpz_t final_B2min, mpz_t final_B2, const int use_ntt, 
	  const int method)
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

  mpz_t B2l, m_1, effB2min, tryeffB2, effB2, lmin;
  /* The best parameters found so far, P == 0 means that no suitable P
     has been found yet: */
  uint64_t P = 0, s_1 = 0, s_2 = 0, l = 0, cost = 0;
  unsigned int i;
  const unsigned int Pvalues_len = sizeof (Pvalues) / sizeof (uint64_t);
  int r;

  outputf (OUTPUT_TRACE, 
           "choose_P(B2min = %Zd, B2 = %Zd, lmax = %lu, min_s2 = %" PRId64 ", "
           "use_ntt = %d, method = %d\n", 
           B2min, B2, lmax, min_s2, use_ntt, method);

  if (mpz_cmp (B2, B2min) < 0)
    return 0L;

  /* If we use the NTT, we allow only power-of-two transform lengths.
     In that case, the code below assumes that lmax is a power of two.
     If that is not the case, print error and return. */
  if (use_ntt && (lmax & (lmax - 1UL)) != 0)
    {
      outputf (OUTPUT_ERROR, 
               "choose_P: Error, lmax = %lu is not a power of two\n", lmax);
      return ECM_ERROR;
    }
  
  mpz_init (effB2);
  mpz_init (tryeffB2);
  mpz_init (effB2min);
  mpz_init (B2l);
  mpz_init (m_1);
  mpz_init (lmin);
  
  mpz_sub (B2l, B2, B2min);
  mpz_add_ui (B2l, B2l, 1UL); /* +1 due to closed interval */
  
  /* For each candidate P, check if [B2min, B2] can be covered at all,
     and if so, what the best parameters (minimizing the cost, maximizing 
     effB2) are. If they are better than the best parameters for the best P 
     so far, remember them. */

  for (i = 0 ; i < Pvalues_len; i++)
    {
      uint64_t tryP, tryphiP, trys_1, trys_2, trycost;
      unsigned long tryl;
      
      tryP = Pvalues[i];
      tryphiP = eulerphi64 (tryP);
      
      outputf (OUTPUT_TRACE, 
	       "choose_P: trying P = %" PRId64 ", eulerphi(P) = "
	       "%" PRId64 "\n", tryP, tryphiP);
      
      /* If we have a good P already and this tryphiP >= cost, then 
	 there's no hope for this tryP, since cost(s_2, l) > eulerphi(P) */
      if (P != 0 && tryphiP >= cost)
	{
	  outputf (OUTPUT_TRACE, 
		   "choose_P: tryphiP > cost = %" PRId64 
		   ", this P is too large\n", cost);
	  continue;
	}
      
      /* We have nr < l and effB2-effB2min <= 2*nr*P. Hence we need 
	 l >= B2l/P/2 */
      mpz_cdiv_q_ui (lmin, B2l, tryP);
      mpz_cdiv_q_2exp (lmin, lmin, 1UL);
      outputf (OUTPUT_TRACE, "choose_P: lmin = %Zd for "
	  			"P = %" PRId64 "\n", lmin, tryP);
      if (mpz_cmp_ui (lmin, lmax) > 0)
	{
	  outputf (OUTPUT_TRACE, 
		   "choose_P: lmin > lmax, this P is too small\n");
	  continue;
	}
      
      /* Try all possible transform lengths and store parameters in 
	 P, s_1, s_2, l if they are better than the previously best ones */
       
      /* Keep reducing tryl to find best parameters. For NTT, we only have 
	 power of 2 lengths so far, so we can simply divide by 2. 
	 For non-NTT, we have arbitrary transform lengths so we can decrease 
	 in smaller steps... let's say by, umm, 25% each time? */
      for (tryl = lmax; mpz_cmp_ui (lmin, tryl) <= 0;
	   tryl = (use_ntt) ? tryl / 2 : 3 * tryl / 4)
	{
	  trys_1 = choose_s_1 (tryphiP, min_s2, tryl / 2, use_ntt);
	  if (trys_1 == 0)
	    {
	      outputf (OUTPUT_TRACE, 
		       "choose_P: could not choose s_1 for "
		       "P = %" PRId64 ", l = %lu\n", tryP, tryl);
	      continue;
	    }
	  ASSERT (tryphiP % trys_1 == 0UL);
	  trys_2 = tryphiP / trys_1;
	  outputf (OUTPUT_TRACE, 
	      	"choose_P: chose s_1 = %" PRId64 ", k = s_2 = %" PRId64 
		   "for P = %" PRId64 ", l = %lu\n", 
		   trys_1, trys_2, tryP, tryl);
	  
	  if (test_P (B2min, B2, m_1, tryP, tryl - trys_1, effB2min, tryeffB2))
	    {
	      outputf (OUTPUT_TRACE, 
		       "choose_P: P = %" PRId64 ", l = %" PRId64 
		       ", s_1 = %" PRId64 ", k = s_2 = %" PRId64
		       "works, m_1 = %Zd, effB2min = %Zd, effB2 = %zZd\n",
		       tryP, tryl, trys_1, trys_2, m_1, effB2min, tryeffB2);
	      /* We use these parameters if we 
		 1. didn't have any suitable ones yet, or 
		 2. these cover [B2min, B2] and are cheaper than the best 
                    ones so far, or 
		 3. they are as expensive but reach greater effB2. */
	      trycost = est_cost (trys_2, tryl, use_ntt, method);
	      ASSERT (tryphiP < trycost);
	      if (P == 0 || trycost < cost ||
		  (trycost == cost && mpz_cmp (tryeffB2, effB2) > 0))
		{
		  outputf (OUTPUT_TRACE, 
			   "choose_P: and is the new "
			   "optimum (cost = %" PRId64 ")\n", trycost);
		  P = tryP;
		  s_1 = trys_1;
		  s_2 = trys_2;
		  l = tryl;
		  cost = trycost;
		  mpz_set (effB2, tryeffB2);
		}
	    }
	}
  }
  
  if (P != 0) /* If we found a suitable P */
    {
      /* Compute m_1, effB2min, effB2 again */
      r = test_P (B2min, B2, m_1, P, l - s_1, effB2min, effB2);
      ASSERT_ALWAYS(r != 0);
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
    }
  
  mpz_clear (effB2);
  mpz_clear (tryeffB2);
  mpz_clear (effB2min);
  mpz_clear (B2l);
  mpz_clear (m_1);
  mpz_clear (lmin);

  return (P != 0) ? (long) P : ECM_ERROR;
}
