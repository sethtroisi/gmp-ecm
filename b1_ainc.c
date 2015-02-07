/* Code to compute "Automatic calculated" B1 incrementation

Copyright 2003, 2005, 2006, 2015 Jim Fougeron, Paul Zimmermann.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
more details.

You should have received a copy of the GNU General Public License
along with this program; see the file COPYING.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#include "ecm-ecm.h"
#include <math.h>

/* return a new value of B1: cur_B1 + incB1val * sqrt(cur_B1) */
double
calc_B1_AutoIncrement (double cur_B1, double incB1val)
{
  return cur_B1 + incB1val * sqrt (cur_B1);
}
