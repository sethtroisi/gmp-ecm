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
#if defined (__MINGW32__) || defined (_MSC_VER)
#include <float.h> /* for DBL_MAX, in MinGW and VC */
#else
#include <values.h> /* for DBL_MAX */
#endif
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
#else
#error "MULT is neither TOOM4, nor TOOM3, nor KARA"
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
unsigned long
bestD (double B2, unsigned int k0, unsigned int *k, unsigned int S,
       unsigned long *estimated_cost)
{
/* the following list contains successive values of (d,a,b) with
   increasing values of a and decreasing values of b/(d*phi(d)/2)
   where the cost of step2 is a + b * ceil(B2 / (d*phi(d)/2)).
*/
#define N 154  
  static double l[N] = {12.0, 18.0, 30.0, 60.0, 90.0, 120.0, 210.0, 240.0, 270.0, 330.0, 420.0, 510.0, 840.0, 1020.0, 1260.0, 1470.0, 1680.0, 1890.0, 2040.0, 2310.0, 3570.0, 3990.0, 4080.0, 4620.0, 5460.0, 6090.0, 7140.0, 8190.0, 9240.0, 14280.0, 16170.0, 16380.0, 18480.0, 19110.0, 21840.0, 24570.0, 24990.0, 25410.0, 28560.0, 30030.0, 32760.0, 34650.0, 39270.0, 60060.0, 66990.0, 67830.0, 71610.0, 78540.0, 79170.0, 90090.0, 92820.0, 103740.0, 103530.0, 106260.0, 106470.0, 108570.0, 120120.0, 131670.0, 133980.0, 139230.0, 150150.0, 157080.0, 158340.0, 240240.0, 237510.0, 237930.0, 270270.0, 274890.0, 278460.0, 300300.0, 314160.0, 324870.0, 330330.0, 371280.0, 420420.0, 417690.0, 431970.0, 450450.0, 510510.0, 570570.0, 565110.0, 600600.0, 628320.0, 649740.0, 690690.0, 1021020.0, 1141140.0, 1138830.0, 1231230.0, 1256640.0, 1291290.0, 1345890.0, 1381380.0, 1711710.0, 1741740.0, 1763580.0, 1806420.0, 1861860.0, 1845690.0, 2042040.0, 2072070.0, 2282280.0, 2277660.0, 2312310.0, 2316930.0, 2372370.0, 2386020.0, 2552550.0, 2513280.0, 2612610.0, 2645370.0, 2691780.0, 2762760.0, 2732730.0, 2852850.0, 4084080.0, 4037670.0, 4144140.0, 4114110.0, 4594590.0, 4555320.0, 4654650.0, 4834830.0, 5105100.0, 5135130.0, 5026560.0, 5225220.0, 5290740.0, 5383560.0, 5615610.0, 5705700.0, 7147140.0, 7054320.0, 7087080.0, 7225680.0, 7417410.0, 7657650.0, 7537530.0, 7987980.0, 8168160.0, 8288280.0, 8198190.0, 8678670.0, 9699690.0, 9669660.0, 9606870.0, 9639630.0, 10210200.0, 10270260.0, 10360350.0, 10330320.0, 10307220.0, 10292100.0, 10341870.0};

  unsigned int i, j;
  unsigned long d, dmin, dF;
  double cost, mincost = DBL_MAX;

  for (i = 0; i < N; i++)
    {
      d = (unsigned long) l[i];
      dF = phi (d) / 2;
      j = (unsigned int) ceil ((B2 / (double) d + 2.0) / (double) dF);
      cost = muls_stage2 (dF, d, S, j);
      if (j >= k0 && cost < mincost)
        {
          mincost = cost;
          dmin = d;
          *k = j;
        }
    }

  if (mincost == DBL_MAX)
    {
      fprintf (stderr, "Error, too large step2 range, please increase k\n");
      exit (1);
    }

  *estimated_cost = (unsigned long) mincost;

  return dmin;
}
