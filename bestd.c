/* Choice of best parameters for stage 2.

  Copyright 2001, 2002, 2003, 2004, 2005 Paul Zimmermann and Alexander Kruppa.

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
  the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
  MA 02111-1307, USA.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h> /* for UINT_MAX */
#include "gmp.h"
#include "ecm.h"
#include "ecm-impl.h"

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

  /* now n is prime or 1 */

  return (n == 1) ? phi : phi * (n - 1);
}

/*
  Compute (d, d2, k) such that:
  (0) k >= k0
  (1) d is a multiple of 6
  (2) k * d * (phi(d)/2) * d2 / phi(d2) >= B2 - B2min
  (3) gcd(d, d2) == 1
  (4) k is minimal, subject to previous conditions
  (5) if parameter po2 is != 0, rounds dF up to a power of 2
  Return non-zero iff an error occurred (too large step 2 bound).
 */
int
bestD (mpz_t B2min, mpz_t B2, int po2, unsigned int *finald, 
       unsigned int *finald2, unsigned int *k, unsigned int *finaldF,
       mpz_t finali0)
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
  static unsigned int l[N] = {12, 18, 30, 42, 60, 90, 120, 150, 210, 240, 270, 330, 420, 510, 630, 840, 1050, 1260, 1470, 1680, 1890, 2310, 2730, 3150, 3570, 3990, 4620, 5460, 6090, 6930, 8190, 9240, 10920, 12180, 13860, 16170, 18480, 20790, 23100, 30030, 34650, 39270, 43890, 48510, 60060, 66990, 78540, 90090, 99330, 120120, 133980, 150150, 180180, 210210, 240240, 270270, 300300, 334950, 371280, 420420, 510510, 570570, 630630, 746130, 870870, 1021020, 1141140, 1291290, 1531530, 1711710, 1891890, 2081310, 2312310, 2552550, 2852850, 3183180, 3573570, 3993990, 4594590, 5105100, 5705700, 6322470, 7147140, 7987980, 8978970, 10210200, 11741730, 13123110, 14804790, 16546530, 19399380, 21411390, 23993970, 26816790, 29609580, 33093060, 36606570, 40330290, 44414370, 49639590};
#define Npo2 23
  static unsigned int lpo2[Npo2] = {12, 30, 60, 120, 240, 510, 1020, 2310, 
                                 4620, 9240, 19110, 39270, 79170, 158340, 
                                 324870, 690690, 1345890, 2852850, 5705700, 
                                 11741730, 23130030, 48498450, 96996900};

  unsigned int i, d, d2, dF, phid;
  mpz_t j, t, i0, i1;
  int r = 0;
  
  mpz_init (i0);
  mpz_init (i1);
  mpz_init (j);
  mpz_init (t);
  
  for (i = 0; i < ((po2) ? Npo2 : N); i++)
    {
      d = (po2) ? lpo2[i] : l[i];
      phid = phi (d) / 2;
      if (po2)
        for (dF = 1; dF < phid; dF <<= 1); /* dF = 2^ceil(log_2(phi(d))) */
      else
        dF = phid;
      /* Look for smallest prime < 25 that does not divide d */
      for (d2 = 5; d2 < 25; d2 += 2)
        {
          if (d2 % 3 == 0)
            continue;
          if (d % d2 > 0)
            break;
        }
      if (d2 >= 25 || d2 - 1 > dF || (d2 > 1 &&
          mpz_cmp_ui (B2min, (d - 1) * d2 - d) <= 0)) /* Would make i0 < 0 */
        d2 = 1;
      
      mpz_set_ui (i0, d - 1);
      mpz_mul_ui (i0, i0, d2);
      mpz_set (j, B2);
      mpz_add (i1, j, i0); /* i1 = B2 + (d - 1) * d2 */
      mpz_set (j, B2min);
      mpz_sub (i0, j, i0); /* i0 = B2min - (d - 1) * d2 */
      mpz_cdiv_q_ui (i0, i0, d); /* i0 = ceil ((B2min - (d - 1) * d2) / d) */
      mpz_fdiv_q_ui (i1, i1, d); /* i1 = floor ((B2 + (d - 1) * d2) / d) */
      
      /* How many roots of G will we need ? */
      mpz_sub (j, i1, i0);
      mpz_add_ui (j, j, 1);

      /* Integer multiples of d2 are skipped (if d2 > 1) */
      if (d2 > 1)
        {
          mpz_fdiv_q_ui (t, i1, d2);
          mpz_sub (j, j, t);
          mpz_fdiv_q_ui (t, i0, d2);
          mpz_add (j, j, t); /* j -= floor (i1 / d2) - floor (i0 / d2) */
        }
      
      /* How many blocks will we need ? Divide lines by dF, rounding up */
      mpz_cdiv_q_ui (j, j, dF);
      
      if (mpz_cmp_ui (j, *k) <= 0 || 
          (*k == ECM_DEFAULT_K && mpz_cmp_ui (j, (po2) ? 9 : 2) <= 0))
        break;
    }

  if (i == N)
    {
      if (*k != 0)
        {
          /* The user asked for a specific k and we couldn't satisfy the
             condition. Nothing we can do ... */
          outputf (OUTPUT_ERROR, "Error: too large step 2 bound, increase -k\n");
          r = ECM_ERROR;
          goto clear_and_exit;
        }
      else if (!mpz_fits_uint_p (j))
        {
          /* Can't fit the number of blocks in an unsigned int. Nothing we
             can do ... */
          outputf (OUTPUT_ERROR, "Error: stage 2 interval too large, cannot "
                   "generate suitable parameters.\nTry a smaller B2 value.\n");
          r = ECM_ERROR;
          goto clear_and_exit;
        }
      /* else: We can fit the number of blocks into an unsigned int, so we use 
         it. This may be a very large value for huge B2-B2min, the user
         is going to notice sooner or later */
    }
  
  /* If the user specified a number of blocks, we'll use that no matter what.
     Since j may be smaller than k, this may increase the B2 limit */
  if (*k == ECM_DEFAULT_K)
    *k = mpz_get_ui (j);

  /* Now that we have the number of blocks, compute real i1. There will be
     k * dF roots of G computed, starting at i0, skipping all that are not
     coprime to d2. While d2 is prime, that means: are not multiples of d2.
     Hence we want i1 so that 
       i1 - floor(i1 / d2) - i0 + ceil((i0 / d2) == k * dF
       i1 - floor(i1 / d2) == k * dF + i0 - ceil((i0 / d2)
  */
  
  mpz_set_ui (j, *k);
  mpz_mul_ui (j, j, dF);
  if (d2 == 1)
    mpz_add (i1, i0, j);
  else
    {
      mpz_add (j, j, i0);
      mpz_cdiv_q_ui (t, i0, d2);
      mpz_sub (j, j, t); /* j = k * dF + i0 - ceil((i0 / d2) */
      mpz_fdiv_qr_ui (j, t, j, d2 - 1);
      mpz_mul_ui (j, j, d2);
      mpz_add (i1, j, t);
    }

  *finald = d;
  *finald2 = d2;
  *finaldF = dF;
  mpz_set (finali0, i0);
  mpz_mul_ui (B2, i1, d);
  
clear_and_exit:
  mpz_clear (t);
  mpz_clear (j);
  mpz_clear (i1);
  mpz_clear (i0);

  return r;
}
