/* 
  Interface code for George Woltman's gwnum library
  
  Copyright 2004-2006 Paul Zimmermann and Alexander Kruppa.
  
  Contains code based on the GWNUM library, 
    copyright 2002-2005 George Woltman, Just For Fun Software, Inc.
  
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
  Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
  02111-1307, USA.
*/

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "ecm-gmp.h"
#include "ecm.h"
#include "ecm-impl.h"
#define ADD_UNDERSCORES
#include "gwdbldbl.h"
#include "gwnum.h"
#include "cpuid.h"

void __gxx_personality_v0()
{
  exit (EXIT_FAILURE);
}

void __cxa_guard_acquire ()
{
  return;
}

void __cxa_guard_release ()
{
  return;
}

static int 
sgn (const int i)
{
  if (i == 0)
    return 0;
  return i > 0 ? 1 : -1;
}


int 
gw_ecm_stage1 (mpz_t f, curve *P, mpmod_t modulus, 
	       double B1, double *B1done, mpz_t go)
{
  const double gw_k = 1.;
  const unsigned long gw_b = 2;
  const unsigned long gw_n = abs (modulus->bits);
  const signed long gw_c = sgn (modulus->bits);
  unsigned long gw_B1done = *B1done;
  unsigned long siz_x, siz_z; /* Size of gw_x and gw_y as longs */
  mpz_t gw_x, gw_z, gw_A;
  int youpi;

  if (mpz_cmp_ui (go, 1) > 0)
    {
      mpres_t b;
      mpres_init (b, modulus);
      mpres_add_ui (b, P->A, 2, modulus);
      mpres_div_2exp (b, b, 2, modulus); /* b == (A+2)/4 */
      ecm_mul (P->x, P->y, go, modulus, b);
      mpres_clear (b, modulus);
    }
  
  outputf (OUTPUT_VERBOSE, 
           "Using gwnum_ecmStage1(%f, %d, %d, %d, , %.0f, %ld)\n",
           gw_k, gw_b, gw_n, gw_c, B1, gw_B1done);

  /* Copy x, z and A values from modular representation to 
     plain integers */

  /* Allocate enough memory for any residue (mod 2^n+-1) for x, z */
  mpz_init2 (gw_x, gw_n + 64);
  mpz_init2 (gw_z, gw_n + 64);
  mpz_init (gw_A);

  /* mpres_get_z always produces non-negative integers */
  mpres_get_z (gw_x, P->x, modulus);
  mpres_get_z (gw_z, P->y, modulus);
  mpres_get_z (gw_A, P->A, modulus);

  /* gwnum_ecmStage1() wants long int pointers for size_x, size_z, 
     so copy them into long int vars */
  siz_x = SIZ(gw_x);
  siz_z = SIZ(gw_z);
  
  youpi = gwnum_ecmStage1 (gw_k, gw_b, gw_n, gw_c, 
      PTR(modulus->orig_modulus), ABSIZ(modulus->orig_modulus), 
      B1, &gw_B1done, PTR(gw_A), ABSIZ(gw_A), 
      PTR(gw_x), &siz_x, PTR(gw_z), &siz_z, NULL, 0);
  
  /* Test that not more was written to gw_x and gw_z than we had space for */
  ASSERT_ALWAYS (siz_x <= (unsigned long) ALLOC(gw_x));
  ASSERT_ALWAYS (siz_z <= (unsigned long) ALLOC(gw_z));
  
  SIZ(gw_x) = siz_x;
  SIZ(gw_z) = siz_z;

  outputf (OUTPUT_DEVVERBOSE, 
           "gw_ecm_stage1: after gwnum_ecmStage1, \n"
           "B1done = %lu, x = %Zd\nz = %Zd\n",
           gw_B1done, gw_x, gw_z);
  
  /* Copy x, z back to P and clean up the temp vars */
  mpres_set_z (P->x, gw_x, modulus);
  mpres_set_z (P->y, gw_z, modulus);
  mpz_clear (gw_A);
  mpz_clear (gw_z);
  mpz_clear (gw_x);

  *B1done = gw_B1done;

  if (youpi > 1)
    {
      outputf (OUTPUT_ERROR, "GW stage 1 returned code %d\n", youpi);
      youpi = ECM_ERROR;
      goto end_of_gwecm;
    }

  if (youpi == 1)
    {
      /* How did that happen? Since we passed z, GWNUM should not do
         an extgcd and so not find factors... but if it did anyways, 
         we deal with it. Who's going to turn down a factor? */
      outputf (OUTPUT_DEVVERBOSE, 
               "gw_ecm_stage1: Strange, gwnum_ecmStage1 reports a factor\n");
      mpres_get_z (f, P->x, modulus);
      youpi = ECM_FACTOR_FOUND_STEP1;
      goto end_of_gwecm;
    }

  /* Normalize z (in P->y) to 1 */
  youpi = ECM_NO_FACTOR_FOUND;
  if (!mpres_invert (P->y, P->y, modulus)) /* Factor found? */
    {
      mpres_gcd (f, P->y, modulus);
      youpi = ECM_FACTOR_FOUND_STEP1;
    } 
  else
    {
      mpres_mul (P->x, P->x, P->y, modulus);
      mpres_set_ui (P->y, 1UL, modulus);
    }

end_of_gwecm:

  return youpi;
}
