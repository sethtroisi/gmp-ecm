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
#include "gmp.h"
#include "ecm.h"

#define N 391

/* returns Euler's totient phi function */
unsigned int
phi (unsigned int n)
{
  unsigned int phi = 1, p;

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
double block_size (unsigned int d)
{
  return (double) d * ((double) phi (d) / 2.0 - 1.0);
}

/*
  Return d such that:
  (1) d is a multiple of 6
  (2) d * (phi(d)/2 - 1) >= B2
  (3) phi(d) is minimal over all d satisfying (1) and (2)
 */
unsigned int
bestD (double B2)
{
/* the following list contains successive values of d with strictly increasing
   values of phi(d) and d*(phi(d)-1) */
   static unsigned int l[N] = 
        {6, 12, 18, 30, 42, 60, 54, 66, 90, 120, 126, 150, 138, 210, 198, 240,
         270, 330, 420, 378, 354, 462, 510, 630, 660, 690, 840, 810, 870, 1050,
        1020, 1260, 1320, 1470, 1386, 1410, 1680, 1650, 1590, 1890, 1770, 2310,
        2130, 2730, 2640, 2940, 3150, 3570, 3990, 4620, 4410, 4830, 5460, 5250,
        5070, 5610, 5670, 6090, 6930, 7140, 7350, 8190, 9240, 9030, 9660, 9450,
        9870, 10920, 11550, 11130, 11220, 11970, 12180, 12390, 13860, 13230,
        14280, 14490, 16170, 16380, 15750, 15510, 16590, 18480, 19110, 19320,
        20790, 21840, 23100, 22260, 22050, 22470, 22440, 24570, 25410, 24990,
        24780, 30030, 28560, 28980, 28350, 28770, 29610, 32340, 32760, 34650,
        33810, 35490, 39270, 38220, 37170, 38640, 43890, 46410, 46200, 44520,
        48510, 51870, 53130, 50190, 60060, 57750, 57330, 57120, 62790, 62370,
        66990, 67830, 71610, 70980, 78540, 76230, 79170, 80850, 82110, 90090,
        92820, 91770, 94710, 99330, 103740, 106260, 103950, 108570, 120120,
        115500, 117390, 122430, 125580, 131670, 133980, 139230, 150150, 144690,
        145530, 157080, 155610, 159390, 158340, 161070, 164010, 164220, 180180,
        173250, 171990, 185640, 188370, 196350, 191730, 189210, 210210, 207480,
        212520, 219450, 217140, 215670, 240240, 233310, 237510, 237930, 244860,
        242550, 251160, 270270, 265650, 274890, 278460, 275310, 300300, 289380,
        307230, 314160, 311220, 330330, 324870, 325710, 334950, 360360, 358050,
        363090, 371280, 371910, 390390, 392700, 395010, 385770, 420420, 417690,
        431970, 450450, 439530, 434070, 436590, 441210, 510510, 489720, 501270,
        502320, 570570, 542850, 556920, 554190, 565110, 600600, 603330, 630630,
        628320, 622440, 690690, 669900, 674310, 746130, 750750, 742560, 743820,
        780780, 785400, 810810, 796950, 870870, 881790, 903210, 930930, 921690,
        1021020, 1009470, 985530, 977130, 979440, 1051050, 1067430, 1141140,
        1138830, 1113840, 1115730, 1171170, 1193010, 1231230, 1206660, 1291290,
        1381380, 1360590, 1411410, 1531530, 1501500, 1504230, 1485120, 1540770,
        1591590, 1610070, 1711710, 1741740, 1763580, 1771770, 1806420, 1861860,
        1845690, 1891890, 2042040, 2072070, 2062830, 2081310, 2132130, 2134860,
        2282280, 2277660, 2227680, 2312310, 2316930, 2372370, 2386020, 2552550,
        2492490, 2459730, 2612610, 2645370, 2762760, 2852850, 2822820, 2788170,
        2765070, 2815890, 3063060, 3028410, 3033030, 3093090, 3081540, 3183180,
        3153150, 3202290, 3213210, 3220140, 3423420, 3453450, 3573570, 3527160,
        3543540, 3513510, 3612840, 3579030, 3730650, 3691380, 3666390, 3753750,
        3993990, 4084080, 4144140, 4081770, 4114110, 4234230, 4354350, 4294290,
        4594590, 4516050, 4555320, 4654650, 4834830, 4774770, 4772040, 5105100,
        5135130, 5047350, 5225220, 5195190, 5290740, 5315310, 5615610, 5705700,
        5645640, 5694150, 5675670, 5735730, 6126120, 6096090, 6276270, 6172530,
        6186180, 6322470, 6636630, 6846840, 6906900, 7147140, 7057050, 7054320,
        7066290, 7087080, 7417410, 7657650, 7597590, 7477470, 7472010, 7379190,
        7537530, 7987980, 8168160, 8288280, 8558550, 8678670, 8708700, 8978970,
        9699690, 9279270, 9669660, 9549540, 9524130, 9519510, 9606870, 9639630,
        10210200, 10270260, 10360350};
   double B;
   int i, j, k;

   B = block_size (l[N-1]);
   if (B2 > B)
     {
       fprintf (stderr, "Error: too large B2, maximal value is k*%1.0f where k is the number of step 2 blocks (see option -k)\n", B);
       exit (EXIT_FAILURE);
     }

   /* perform dichotomy search through l */
   i = 0;
   j = N - 1;
   /* invariant block_size(l[i]) < B2 <= block_size(l[j]) */
   while (i + 1 < j)
     {
      k = (i + j) / 2;
      B = block_size (l[k]);
      if (B < B2)
	i = k;
      else
	j = k;
     }

   return l[j];
}
