/* Choice of best parameters for stage 2.

  Copyright 2001, 2002, 2003 Paul Zimmermann and Alexander Kruppa.

  This program is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 2 of the License, or (at your
  option) any later version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
  more details.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h> /* for DBL_MAX */
#include "gmp.h"
#include "ecm.h"

/* number of multiplication of Karatsuba */
static unsigned long
muls_kara (unsigned int K)
{
  unsigned int k, l;

  if (K <= 1)
    return K;

  if (K == 3)
    return 6; /* Weimerskirch/Paar trick */
  
  k = K / 2;
  l = K - k;
  return (k == l) ? 3 * muls_kara (k) : 2 * muls_kara (l) + muls_kara (k);
}

/* number of multiplication of toomcook3 */
unsigned long
muls_toom3 (unsigned int n)
{
  unsigned int l, k;
  if (n <= 2 || n == 4)
    return muls_kara (n);

  l = (n + 2) / 3;
  k = n - 2 * l;
  return (l == k) ? 5 * muls_toom3 (l) : 4 * muls_toom3 (l) + muls_toom3 (k);
}

/* number of multiplication of toomcook4 */
static unsigned long
muls_toom4 (unsigned int n)
{
  unsigned int l, k;

  if (n <= 2)
    return muls_kara (n);

  if (n == 3 || n == 5 || n == 6 || n == 9 || n == 17 || n == 18 ||
      (25 <= n && n <= 27) || (77 <= n && n <= 81))
    return muls_toom3 (n);

  l = (n + 3) / 4;
  k = n - 3 * l;
  return (l == k) ? 7 * muls_toom4 (l) : 6 * muls_toom4 (l) + muls_toom4 (k);
}

/* number of multiplication of toomcook4, short division */
static unsigned long
muls_toom4_short (unsigned int n)
{
  unsigned int p, q;
  unsigned long muls;

  switch (n)
    {
    case 1:
      return 1;
    case 2:
      return 3;
    case 3:
      return 5;
    default:
      for (p = 1; 4 * p <= n; p *= 4);
      p = (n / p) * p;
      muls = muls_toom4 (p);
      if ((q = n - p))
        muls += 2 * muls_toom4_short (q);
      return muls;
    }
}

unsigned long
muls_gen (unsigned int n)
{
#if (MULT == TOOM4)
  return muls_toom4 (n);
#elif (MULT == TOOM3)
  return muls_toom3 (n);
#elif (MULT == KARA)
  return muls_kara (n);
#elif (MULT == KS)
  return muls_toom4 (n); /* approximate */
#else
#error "MULT is neither KS, nor TOOM4, nor TOOM3, nor KARA"
#endif
}

unsigned long
muls_gen_short (unsigned int n)
{
  unsigned int p, q, muls;
  switch (n)
    {
    case 1:
      return 1;
    case 2:
      return 3;
    case 3:
      return 5;
    default:
      for (p = 1; MULT * p <= n; p *= MULT);
      p = (n / p) * p;
      muls = muls_gen (p);
      if ((q = n - p))
        muls += 2 * muls_gen_short (q);
      return muls;
    }
}

static unsigned long
muls_polyfromroots (unsigned int k)
{
  unsigned int l, m;
  unsigned long muls;

  if (k == 1)
    return 0;
  
   m = k / 2;
   l = k - m;
   muls = (m == l) ? 2 * muls_polyfromroots (l)
     : muls_polyfromroots (l) + muls_polyfromroots (m);
   muls += muls_gen (m);
   if (l > m)
     muls += m;
   return muls;
}

static unsigned long
muls_polyinvert (unsigned int K)
{
  unsigned int k, l;
  unsigned long muls;

  if (K == 1)
    return 0;

  k = K / 2;
  l = K - k;

  muls = muls_polyinvert (l) + muls_gen (l);
  if (k > 1)
    {
      muls += muls_gen (k - 1);
      if (l > k)
        muls += k - 1;
    }
  muls += muls_gen (k);
  return muls;
}

static unsigned long
muls_prerevertdiv (unsigned int n)
{
  return muls_gen_short (n - 1) + muls_gen_short (n);
}

static unsigned long
muls_recdiv (unsigned int K)
{
  unsigned int k, l;
  unsigned long muls;

  if (K == 1)
    return 1;

  k = K / 2;
  l = K - k;

  muls = (k < l) ? 2 * k + muls_recdiv (l) + muls_recdiv (k)
    : 2 * muls_recdiv (k);
  muls += 2 * muls_gen (k);
  return muls;
}

static unsigned long
muls_polyeval (unsigned int k)
{
  unsigned int l, m;
  unsigned long muls;

  if (k == 1)
    return 0;

  m = k / 2;
  l = k - m;
  muls = (m < l) ? m + muls_recdiv (m) + muls_recdiv (l) : 2 * muls_recdiv (l);
  muls += (m < l) ? muls_polyeval (l) + muls_polyeval (m)
    : 2 * muls_polyeval (m);
  return muls;
}

/* estimated number of multiplications for stage 2 [ecm] */
static double
muls_stage2 (unsigned int dF, unsigned int d, unsigned int S, unsigned int k)
{
  double muls;

  muls = (double) (d - 6) * S; /* ecm_rootsF */
  muls += (double) (k + 1) * muls_polyfromroots (dF);
  muls += (double) muls_polyinvert (dF);
  muls += (double) k * 6 * S * dF;
  muls += (k - 1) * (double) muls_gen (dF);
  muls += (k - 1) * (double) muls_prerevertdiv (dF);
#ifdef POLYEVALTELLEGEN
  muls += (double) muls_polyeval_tellegen (dF);
#else
  muls += (double) muls_polyeval (dF);
#endif
  return muls;
}

/* returns Euler's totient phi function */
unsigned long
phi (unsigned long n)
{
  unsigned long phi = 1, p;

  for (p = 2; p * p <= n; p += 2)
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

      if (p == 2)
	p--;
    }

  /* now n is prime */

  return (n == 1) ? phi : phi * (n - 1);
}

/* return the maximal block size corresponding to d */
double
block_size (unsigned long d)
{
  return (double) d * ((double) phi (d) / 2.0);
}

/*
  Return (d,k) such that:
  (0) k >= k0
  (1) d is a multiple of 6
  (2) k * d * (phi(d)/2) >= B2
  (3) the cost of step 2 is minimal
 */
void
bestD (double B2min, double B2, unsigned int k0, unsigned int *finald, 
       unsigned int *finald2, unsigned int *k)
{
/* the following list contains successive values of b with
   increasing values of phi(b). It was generated by the following Maple
   program:
l := [[1,1]]:
for b from 12 by 6 do
   d:=numtheory[phi](b)/2;
   while d <= l[nops(l)][2] do l:=subsop(nops(l)=NULL, l) od;
   n := nops(l);
   if b>1.1*l[n][1] then l := [op(l), [b,d]]; lprint(l) fi;
od:
*/
#define N 100
  static double l[N] = {12, 18, 30, 42, 60, 90, 120, 150, 210, 240, 270, 330, 420, 510, 630, 840, 1050, 1260, 1470, 1680, 1890, 2310, 2730, 3150, 3570, 3990, 4620, 5460, 6090, 6930, 8190, 9240, 10920, 12180, 13860, 16170, 18480, 20790, 23100, 30030, 34650, 39270, 43890, 48510, 60060, 66990, 78540, 90090, 99330, 120120, 133980, 150150, 180180, 210210, 240240, 270270, 300300, 334950, 371280, 420420, 510510, 570570, 630630, 746130, 870870, 1021020, 1141140, 1291290, 1531530, 1711710, 1891890, 2081310, 2312310, 2552550, 2852850, 3183180, 3573570, 3993990, 4594590, 5105100, 5705700, 6322470, 7147140, 7987980, 8978970, 10210200, 11741730, 13123110, 14804790, 16546530, 19399380, 21411390, 23993970, 26816790, 29609580, 33093060, 36606570, 40330290, 44414370, 49639590};

  unsigned int i;
  unsigned int d, d2, dmin = 0, d2min = 1, dF;
  double j, dd2, jmin = DBL_MAX;
  int found = 0; /* non-zero when some j >= k0 was found */

  for (i = 0; i < N; i++)
    {
      d = (unsigned long) l[N - 1 - i];
      dF = phi (d) / 2;
      /* Look for smallest prime < 25 that does not divide d */
      for (d2 = 5; d2 < 25; d2 += 2)
        {
          if (d2 % 3 == 0)
            continue;
          if (d % d2 > 0)
            break;
        }
      if (d2 >= 25 || d2 - 1 > dF)
        d2 = 1;
      dd2 = (double) d * (double) d2;
      j = ceil (B2 / dd2) - floor (B2min / dd2);
      j = ceil (j / (double) dF * (double) phi (d2));
      /* warning: if B2min=B2, we may always have j=1 here */
      if (j >= (double) k0 || found == 0)
	{
	  /* if found=0 and j>=k0: forgot previous and keep that one */
	  if ((found == 0 && j >= (double) k0) || j < jmin
	      || (j == jmin && d < dmin))
	    {
	      dmin = d;
	      d2min = d2;
	      jmin = j;
	    }
	  found = (j >= (double) k0);
	}
    }

  if (dmin == 0)
    {
      fprintf (stderr, "Error, too large step2 range, please increase k\n");
      exit (1);
    }

  *k = (unsigned int) jmin;
  *finald = dmin;
  *finald2 = d2min;
  
  return;
}

/*
  Return (d,d2,k) such that:
  (0) k == k0 if k0 > 0, 2 <= k <= 9 otherwise
  (1) d is a multiple of 6
  (2) dF = 2^ceil(log_2(phi(d)/2)), a power of 2
  (3) k * d * d2 / phi(d2) * dF + floor(B2min / d / d2) * d * d2 >= B2
  (4) d*d2/phi(d2) maximal
*/
void
bestD_po2 (double B2min, double B2, unsigned int *finald, 
           unsigned int *finald2, unsigned int *k, unsigned long *est_muls)
{
/* List of d values where phi(d)/2 is just below or equal a power of 2 */
#define Npo2 23
  static unsigned int l[Npo2] = {12, 30, 60, 120, 240, 510, 1020, 2310, 4620, 
                                 9240, 19110, 39270, 79170, 158340, 324870, 
                                 690690, 1345890, 2852850, 5705700, 11741730, 
                                 23130030, 48498450, 96996900};
  /* l2[i] is always the smallest prime (or 1) not present in l[i] and < dF */
  static unsigned int l2[Npo2] = {1, 1, 7, 7, 7, 7, 7, 13, 13, 
                                 13, 11, 13, 11, 11, 11, 
                                 17, 11, 17, 17, 19, 
                                 13, 23, 23};

  unsigned int i, j, d, d2, dF;
  double dd2;

  /* Find the smallest d so that the required number of blocks to cover
     a stage 2 interval of length B2-B2min is no greater than the given k, 
     or no greater than 9 if *k == 0 */
  for (i = 0; i < Npo2; i++)
    {
      d = l[i];
      d2 = l2[i];
      j = phi (d) / 2;
      for (dF = 1; dF < j; dF <<= 1); /* dF = 2^ceil(log_2(phi(d))) */
      dd2 = (double) d * (double) d2;
      /* How many blocks will we need */
      j = ceil ( (ceil (B2 / dd2) - floor (B2min / dd2)) * 
                 (double) phi (d2) / (double) dF ); 
      if (j <= *k || (*k == 0 && j <= 9))
        break;
    }

  if (i == Npo2)
    {
      fprintf (stderr, "Error, too large step2 range, please increase k\n");
      exit (1);
    }
  
  /* If the user specified a number of blocks, we'll use that no matter what.
     Since j may be smaller than k, this may increase the B2 limit */
  if (*k == 0)
    *k = j;
  *finald = d;
  *finald2 = d2;

  /* Q&D cost estimation. Doesn't account for Karatsuba and Toom-Cook in F_mul
     and gets the Polyeval cost rather wrong */
  *est_muls = 6 * (1 + *k) * dF +   /* Roots of F and G (for ECM) */
              (*k + 1) * dF * (unsigned long)(log(dF) / log(2) + 0.5) + /* Building from roots */
              6 * dF - 16 +         /* Inverting F */
              (*k - 1) * 5 * dF +   /* G = G*H (mod F) */
              2 * dF * (unsigned long)(log(dF) / log(2.) + 0.5) + dF; /* Polyeval */
  
  return;
}
