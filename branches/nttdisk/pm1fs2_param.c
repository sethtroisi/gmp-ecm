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
    936020085UL, 984518535UL, 1017041025UL, 1052416365UL
#if (ULONG_MAX > 4294967295)
   ,1086110025UL, 1110614505UL, 1147371225UL, 1191785595UL, 1213887675UL,
    1265809545UL, 1282356075UL, 1331995665UL, 1391905515UL, 1450103655UL,
    1479202725UL, 1547100555UL, 1555088535UL, 1673196525UL, 1712565855UL,
    1767130365UL, 1830673845UL, 1883166285UL, 1954487535UL, 2001964965UL,
    2119382265UL, 2187280095UL, 2255177925UL, 2342475135UL, 2390973585UL,
    2421213795UL, 2555868315UL, 2672264595UL, 2788660875UL, 2856558705UL,
    2953555605UL, 3050552505UL, 3234846615UL, 3457939485UL, 3516137625UL,
    3681032355UL, 3758629875UL, 3904125225UL, 4127218095UL, 4360010655UL,
    4573403835UL, 4796496705UL, 4844995155UL, 5019589575UL, 5203883685UL,
    5262081825UL, 5465775315UL, 5766465705UL, 5898837945UL, 6164152995UL,
    6358146795UL, 6411780375UL, 6804332535UL, 6980458485UL, 7172920755UL,
    7473611145UL, 7716103395UL, 7968295335UL, 8182259085UL, 8342499165UL,
    8812168365UL, 9023519505UL, 9704539845UL, 9927632715UL, 10373818455UL,
    10439434005UL, 10820004195UL, 11043097065UL, 11489282805UL,
    11877270405UL, 12381654285UL, 12604747155UL, 13080031965UL,
    13274025765UL, 13642613985UL, 14389490115UL, 14583483915UL,
    15058768725UL, 15611651055UL, 16174233075UL, 16397325945UL,
    17289697425UL, 17735883165UL, 18143270145UL, 18381678315UL,
    19074440385UL, 19559424885UL, 20636090475UL, 20941375455UL,
    21800053275UL, 22643926305UL, 23148310185UL, 24205576395UL,
    24546777255UL, 25544133615UL, 26389538175UL, 26863291455UL,
    27813861075UL, 29113619535UL, 29494189725UL, 30520074585UL,
    30684969315UL, 31790733975UL, 33575476935UL, 34467848415UL,
    35202742575UL, 36427185795UL, 38037334335UL, 39240095895UL,
    40365259935UL, 42053005995UL, 43168470345UL, 44953213305UL,
    45845584785UL, 48522699225UL, 50307442185UL, 51869092275UL,
    53653835235UL, 54546206715UL, 56680138515UL, 58784971245UL,
    59386352025UL, 61908271425UL, 63431122755UL, 65700850215UL,
    67931778915UL, 70162707615UL, 72616729185UL, 74120181135UL,
    75740029365UL, 78417143805UL, 80871165375UL, 82840202445UL,
    86448487125UL, 88466022645UL, 91133437395UL, 92918180355UL,
    100280245065UL, 100726430805UL, 102811864155UL, 106749938295UL,
    109000266375UL, 113219631525UL, 119689324755UL, 121027881975UL,
    127943760945UL, 132628711215UL, 134859639915UL, 141775518885UL,
    148691397855UL, 150922326555UL, 155607276825UL, 161320394235UL,
    164977177365UL, 171446870595UL, 177470378085UL, 183270792705UL
#endif
};

/* All the prime factors that can appear in eulerphi(P) */
const unsigned long phiPfactors[] = {2UL, 3UL, 5UL, 7UL, 11UL, 13UL, 
				     17UL, 19UL};


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
test_P (const mpz_t B2min, const mpz_t B2, mpz_t m_1, const unsigned long P, 
	const unsigned long nr, mpz_t effB2min, mpz_t effB2)
{
  mpz_t m;
  /* We need B2min >= 2 * max(S_1 + S_2) + (2*m_1 - 1)*P + 1, or
     B2min - 2 * max(S_1 + S_2) - 1 >= (2*m_1)*P - P, or
     (B2min - 2*max(S_1 + S_2) + P - 1)/(2P) >= m_1
     Choose m_1 accordingly */
  
  mpz_init (m);
  sets_max (m, P);
  mpz_mul_2exp (m, m, 1UL); /* m = 2*max(S_1 + S_2) */

  mpz_sub (m_1, B2min, m);
  mpz_sub_ui (m_1, m_1, 1UL); /* m_1 = B2min - 2*max(S_1 + S_2) - 1 */
  mpz_add_ui (m_1, m_1, P);
  mpz_fdiv_q_2exp (m_1, m_1, 1UL);
  mpz_fdiv_q_ui (m_1, m_1, P);    /* 2UL*P may overflow */
  
  /* Compute effB2min = 2 * max(S_1 + S_2) + (2*(m_1 - 1) + 1)*P + 1 */
  
  mpz_mul_2exp (effB2min, m_1, 1UL);
  mpz_sub_ui (effB2min, effB2min, 1UL);
  mpz_mul_ui (effB2min, effB2min, P);
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
  mpz_mul_ui (effB2, effB2, P);
  mpz_sub (effB2, effB2, m);

  /* The effective B2 values is that value, minus 1 */
  mpz_sub_ui (effB2, effB2, 1UL);

  mpz_clear (m);
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


/* Approximate cost of stage 2. Cost with and without ntt are not 
   comparable. We have l > s_1 and s_1 * s_2 = eulerphi(P), hence
   s_2*l > eulerphi(P) and so cost (s_2, l) > eulerphi(P) for all P */
static unsigned long 
est_cost (const unsigned long s_2, const unsigned long l, const int use_ntt,
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
	  const unsigned long min_s2, faststage2_param_t *finalparams, 
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
  unsigned long P = 0, s_1 = 0, s_2 = 0, l = 0, cost = 0;
  unsigned int i;
  const unsigned int Pvalues_len = sizeof (Pvalues) / sizeof (unsigned long);
  int r;

  outputf (OUTPUT_TRACE, 
           "choose_P(B2min = %Zd, B2 = %Zd, lmax = %lu, min_s2 = %ld, "
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
      unsigned long tryP, tryphiP, trys_1, trys_2, tryl, trycost;
      
      tryP = Pvalues[i];
      tryphiP = eulerphi (tryP);
      
      outputf (OUTPUT_TRACE, 
	       "choose_P: trying P = %lu, eulerphi(P) = %lu\n", tryP, tryphiP);
      
      /* If we have a good P already and this tryphiP >= cost, then 
	 there's no hope for this tryP, since cost(s_2, l) > eulerphi(P) */
      if (P != 0 && tryphiP >= cost)
	{
	  outputf (OUTPUT_TRACE, 
		   "choose_P: tryphiP > cost = %lu, this P is too large\n",
		   cost);
	  continue;
	}
      
      /* We have nr < l and effB2-effB2min <= 2*nr*P. Hence we need 
	 l >= B2l/P/2 */
      mpz_cdiv_q_ui (lmin, B2l, tryP);
      mpz_cdiv_q_2exp (lmin, lmin, 1UL);
      outputf (OUTPUT_TRACE, "choose_P: lmin = %Zd for P = %lu\n", lmin, tryP);
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
		       "choose_P: could not choose s_1 for P = %lu, l = %lu\n",
		       tryP, tryl);
	      continue;
	    }
	  ASSERT (tryphiP % trys_1 == 0UL);
	  trys_2 = tryphiP / trys_1;
	  outputf (OUTPUT_TRACE, "choose_P: chose s_1 = %lu, k = s_2 = %lu "
		   "for P = %lu, l = %lu\n", trys_1, trys_2, tryP, tryl);
	  
	  if (test_P (B2min, B2, m_1, tryP, tryl - trys_1, effB2min, tryeffB2))
	    {
	      outputf (OUTPUT_TRACE, 
		       "choose_P: P = %lu, l = %lu, s_1 = %lu, k = s_2 = %lu "
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
			   "choose_P: and is the new optimum (cost = %lu)\n",
			   trycost);
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
