/* Code to compute "Automatic calculated" B1 incrementation

  Copyright 2003, 2005 Jim Fougeron, Paul Zimmermann.

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
#include <gmp.h>
#include "ecm.h"
#include "ecm-ecm.h"

/* 
 * Version #2 function is the one we are using with a const
 * adjustment of 1.33
 */

/* Version #1 and Version #3 are not being used, but they have been kept in 
   the source, so that they can be refered to if and when changes are made */
double calc_B1_AutoIncrement_v3 (double cur_B1, double incB1val, int calcInc);
double calc_B1_AutoIncrement_v1 (double cur_B1, double incB1val, int calcInc);

/* Here is my "first" attempt at a B1 adjustment function.
 * Parameters:
 *  cur_B1    the current B1 level
 *  incB1val  This is ether a constant or an "adjustment factor"
 *  calcInc   Tells whether incB1val is a constant, or whether we
 *            should compute the optimal B1 adjustment, and then
 *            "adjust" that optimal value based up incB1val, which
 *            is then treadted as a scaling factor.
 *
 * Returns:   the new B1 value to use.
 *
 * Assumption The return value is based upon providing the recommended
 *            optimal number of curves for a "range" of B1's, then 
 *            computing the amount of adjustment needed to push B1 to 
 *            the next level.  NOTE this may be too slow of a push.
 *            If it optimal curves at 250000 is 500, using 500 curves
 *            in an ever advancing B1 level from 250000 to the next
 *            level (1e6) is considerably more work than the simple
 *            500 curves at B1=250000.  It might be prudent to make
 *            the adjustment deal with the low bound value, the high
 *            bound value, and how far the current B1 is from the low
 *            and the high B1 boundary.
 */
double
calc_B1_AutoIncrement_v1 (double cur_B1, double incB1val, int calcInc)
{
  double B1Mod;
  if (!calcInc)
    return cur_B1 + incB1val;  /* incB1val is a constant to add to B1 */

  /* This simple table was "created" based upon the "Optimal B1 table"
     in the README file */
  if (cur_B1 < 2000.)
    B1Mod = 200.;
  else if (cur_B1 < 11000.)  /* 30 curves from B1=2000 to 11000 */
    B1Mod = 300.;
  else if (cur_B1 < 50000.)  /* 90 curves from B1=11000 to 50000 */
    B1Mod = 433.3334;
  else if (cur_B1 < 250000.)  /* 240 curves from B1=50000 to 250000 */
    B1Mod = 833.3334;
  else if (cur_B1 < 1000000.)  /* 500 curves from B1=250000 to 1e6 */
    B1Mod = 1500.;
  else if (cur_B1 < 3000000.)  /* 1100 curves from B1=1e6 to 3e6 */
    B1Mod = 1818.18182;
  else if (cur_B1 < 11000000.)  /* 2900 curves from B1=3e6 to 11e6 */
    B1Mod = 2758.621;
  else if (cur_B1 < 43000000.)  /* 5500 curves from B1=11e6 to 43e6 */
    B1Mod = 5818.18182;
  else if (cur_B1 < 110000000.)  /* 9000 curves from B1=43e6 to 11e7 */
    B1Mod = 7444.44445;
  else if (cur_B1 < 260000000.)  /* 22000 curves from B1=11e7 to 26e7 */
    B1Mod = 6818.18182;
  else if (cur_B1 < 850000000.)  /* 52000 curves from B1=26e7 to 85e7 */
    B1Mod = 11346.1539;
  else if (cur_B1 < 2900000000.)  /* 83000 curves from B1=85e7 to 29e8 */
    B1Mod = 24698.8;
  else
    B1Mod = 35000.;

  return floor (cur_B1 + (B1Mod*incB1val) + 0.5);
}

/* Here is my "second" attempt at a B1 adjustment function.
 * this version looks pretty good 
 *
 * THIS is the version being used.
 */
double
calc_B1_AutoIncrement (double cur_B1, double incB1val, int calcInc)
{
  const double const_adj = 1.33;
  double B1Mod;
  if (!calcInc)
    return cur_B1 + incB1val;  /* incB1val is a constant to add to B1 */

  /* This simple table was "created" based upon the "Optimal B1 table"
     in the README file */
  if (cur_B1 < 2000.)
    B1Mod = 200.;
  else if (cur_B1 < 11000.)  /* 30 curves from B1=2000 to 11000 */
    {
      B1Mod = 300.    * (1. - ((cur_B1 - 2000.) / 9000.));
      B1Mod +=433.334 * (1. - ((11000. - cur_B1) / 9000.));
    }
  else if (cur_B1 < 50000.)  /* 90 curves from B1=11000 to 50000 */
    {
      B1Mod = 433.334 * (1. - ((cur_B1 - 11000.) / 39000.));
      B1Mod +=833.334 * (1. - ((50000. - cur_B1) / 39000.));
    }
  else if (cur_B1 < 250000.)  /* 240 curves from B1=50000 to 250000 */
    {
      B1Mod = 833.334 * (1. - ((cur_B1 - 50000.) / 200000.));
      B1Mod +=1500.   * (1. - ((250000. - cur_B1) / 200000.));
    }
  else if (cur_B1 < 1000000.)  /* 500 curves from B1=250000 to 1e6 */
    {
      B1Mod = 1500.        * (1. - ((cur_B1 - 250000.) / 750000.));
      B1Mod +=1818.18182   * (1. - ((1000000. - cur_B1) / 750000.));
    }
  else if (cur_B1 < 3000000.)  /* 1100 curves from B1=1e6 to 3e6 */
    {
      B1Mod = 1818.18182   * (1. - ((cur_B1 - 1000000.) / 2000000.));
      B1Mod +=2758.621     * (1. - ((3000000. - cur_B1) / 2000000.));
    }
  else if (cur_B1 < 11000000.)  /* 2900 curves from B1=3e6 to 11e6 */
    {
      B1Mod = 2758.621     * (1. - ((cur_B1 - 3000000.) / 8000000.));
      B1Mod +=5818.18182   * (1. - ((11000000. - cur_B1) / 8000000.));
    }
  else if (cur_B1 < 43000000.)  /* 5500 curves from B1=11e6 to 43e6 */
    {
      B1Mod = 5818.18182   * (1. - ((cur_B1 - 11000000.) / 32000000.));
      B1Mod +=7444.44445   * (1. - ((43000000. - cur_B1) / 32000000.));
    }
  else if (cur_B1 < 110000000.)  /* 9000 curves from B1=43e6 to 11e7 */
    {
      B1Mod = 7444.44445   * (1. - ((cur_B1 - 43000000.)  / 67000000.));
      B1Mod +=6818.18182   * (1. - ((110000000. - cur_B1) / 67000000.));
    }
  else if (cur_B1 < 260000000.)  /* 22000 curves from B1=11e7 to 26e7 */
    {
      B1Mod = 6818.18182   * (1. - ((cur_B1 - 110000000.) / 150000000.));
      B1Mod +=11346.1539   * (1. - ((260000000. - cur_B1) / 150000000.));
    }
  else if (cur_B1 < 850000000.)  /* 52000 curves from B1=26e7 to 85e7 */
    {
      B1Mod = 11346.1539   * (1. - ((cur_B1 - 260000000.) / 590000000.));
      B1Mod +=24698.8      * (1. - ((850000000. - cur_B1) / 590000000.));
    }
  else if (cur_B1 < 2900000000.)  /* 83000 curves from B1=85e7 to 29e8 */
    {
      B1Mod = 24698.8      * (1. - ((cur_B1 - 850000000.)  / 2050000000.));
      B1Mod +=50000.0      * (1. - ((2900000000. - cur_B1) / 2050000000.));
    }
  else
    B1Mod = 50000.;

  return floor (cur_B1 + const_adj*(B1Mod*incB1val) + 0.5);
}


/* Here is my "third" attempt at a B1 adjustment function.
 * It seems to adjust too quickly
 */
double B1Min[12] = 
{ 2000.0,  11000.0, 50000.0,  250000.0,  1000000.0, 3000000.0,  11000000.0, 43000000.0,  110000000.0, 260000000.0, 850000000.0,  2900000000.0 };
double B1Max[12] = 
{ 11000.0, 50000.0, 250000.0, 1000000.0, 3000000.0, 11000000.0, 43000000.0, 110000000.0, 260000000.0, 850000000.0, 2900000000.0, 9000000000.0 };
double B1Inc[12] = 
{ 300.0,   433.334, 833.334,  1500.0,    1818.1819, 2758.621,   5818.1819,  7444.4445,   6818.1819,   11346.1539,  24698.8,      50000.0 };
/*B1Table_t B1Table[12] =
    {300,0       ,2000.0       ,11000.0 },
    {433.334,    ,11000.0      ,50000.0 },
    {833.334,    ,50000.0      ,250000.0 },
    {1500.0      ,250000.0     ,1000000.0 },
    {1818.1819,  ,1000000.0    ,3000000.0 },
    {2758.621,   ,3000000.0    ,11000000.0 },
    {5818.1819,  ,11000000.0   ,43000000.0 },
    {7444.4445,  ,43000000.0   ,110000000.0 },
    {6818.1819,  ,110000000.0  ,260000000.0 }, NOTE the increment does not look larger enough here!!
    {11346.1539, ,260000000.0  ,850000000.0 },
    {24698.8,    ,850000000.0  ,2900000000.0 },
    {50000.0,    ,2900000000.0 ,9000000000.0 };
*/

double
calc_B1_AutoIncrement_v3 (double cur_B1, double incB1val, int calcInc)
{
  double B1Mod;
  if (!calcInc)
    return cur_B1 + incB1val;  /* incB1val is a constant to add to B1 */

  /* This simple table was "created" based upon the "Optimal B1 table"
     in the README file */
  if (cur_B1 < 2000.)
    B1Mod = 200.;
  else if (cur_B1 > 2900000000.)
    B1Mod = 50000;
  else
    {
      double OrigMin;
      int i = 0;
      while (i < 11 && B1Max[i] < cur_B1)
	  ++i;
      B1Mod =  B1Inc[i] * (1. - ((cur_B1 - B1Min[i]) / (B1Max[i] - B1Min[i])));
      OrigMin = B1Min[i];
      while (++i < 12)
	{
	  B1Mod +=  B1Inc[i] * (1. - ((B1Min[i] - cur_B1) / (B1Min[i] - OrigMin)));
	}
    }
  return floor (cur_B1 + (B1Mod*incB1val) + 0.5);
}
