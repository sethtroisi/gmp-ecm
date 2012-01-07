/* 
  Interface code for George Woltman's gwnum library
  
  Copyright 2004, 2005, 2006, 2008, 2011 Paul Zimmermann and Alexander Kruppa.
  
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
#include <math.h> /* for rint */
#include <gmp.h>
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

/* With the following 2 functions, we try to find a representation of an
   input number in the form of z = k*b^n+c. If such a representation was
   found, set the the appropriate values and return 1. Otherwise, set b to
   zero and return 0. */

/* This function searches for a representation of z of the form k*b^n+c */
int
kbnc_z (double *k, unsigned long *b, unsigned long *n, signed long *c, mpz_t z)
{
  int i = 0;
  int j = 0;
  int exp = 1;
  int check_it_out = 0;
  int ret = 0;
  mpz_t diff;
  mpz_t abs_diff;
  mpz_t b_n; /* this will = base^exp */
  mpz_t k_b_n; /* this will = k*b^n */
  mpz_t test_k;
  mpz_t max_k;
  mpz_t lhs; /* used for finding the k value */
  mpz_t rhs; /* used for finding the k value */
  mpz_t base;
  mpz_t base_min;
  mpz_t base_max;
  unsigned long test_k_ui = 0;

  /* this puts a bound on how large our C value can be */
  int max_diff = 8388607;

  /* make sure we have a place to put our results */
  if (k == NULL || b == NULL || n == NULL || c == NULL)
    return 0;

  /* make sure the input meets some sort of minimum size requirement.
     The gwnum library reports ES1_CANNOT_DO_QUICKLY for number < 2^350 */
  if (mpz_sizeinbase(z, 2) < 350)
  {
    *b = 0;
    return 0;
  }

  mpz_init (diff);
  mpz_init (abs_diff);
  mpz_init (b_n);
  mpz_init (k_b_n);
  mpz_init (lhs);
  mpz_init (rhs);
  mpz_init (test_k);
  mpz_init (base);
  mpz_init_set_ui (base_min, 2);
  mpz_init_set_ui (base_max, 10000);
  /* this puts a bound on how large of a k value we want to find */
  mpz_init_set_str (max_k, "562949953421312", 10);

  /* when dividing: z/(base^exp) this will give us a possible k value */
  /* we want a quick test to see if this might be a viable k value */
  /* so, we want this k value to be close to an integer */
  /* ie, test_k = 13.99999, is pretty close to the integer 14 */
  /* since it is "pretty close", we can test this k value. */
  /* whereas test_k = 13.5689, is not "pretty close" to an integer */
  /* so, we will not run extra tests with this k value */
  /* should we change this based on the size of z? */
  /* for now, the code checks to see whether test_k is with 1/1000 of an integer */

  for (mpz_set (base, base_min); mpz_cmp (base, base_max) <= 0;
       mpz_add_ui (base, base, 1))
  {
    exp = (mpz_sizeinbase (z, 2) - 1) / (mpz_sizeinbase (base, 2) - 1) + 1;
    mpz_pow_ui (b_n, base, exp); /* base^exp should be > z here */

    while (1)
      {
        check_it_out = 0; /* 0 */

        mpz_tdiv_q (test_k, z, b_n);

        if (mpz_cmp(test_k, max_k) > 0)
          break;

        /* check to see if test_k is "pretty close" to the next smallest
           integer: z/b_n - test_k <= 1/1000  # z/b_n should be > test_k here
           z/b_n <= 1/1000 + test_k
           1000*z/b_n <= 1 + 1000*test_k
           if (1000*z <= b_n + 1000*b_n*test_k) */
        mpz_mul_ui (lhs, z, 1000);
        mpz_mul (rhs, b_n, test_k);
        mpz_mul_ui (rhs, rhs, 1000);
        mpz_add (rhs, rhs, b_n);
        if (mpz_cmp (lhs, rhs) <= 0)
          check_it_out = 1;

        /* check to see if test_k is "pretty close" to the next largest
           integer */
        if (!check_it_out)
          {
            mpz_add_ui (test_k, test_k, 1);
            /* test_k - z/b_n <= 1/1000  # test_k should be > z/b_n here */
            /* test_k <= 1/1000 + z/b_n */
            /* test_k - 1/1000 <= z/b_n */
            /* 1000*test_k - 1 <= 1000*z/b_n */
            /* if (1000*b_n*test_k - b_n <= 1000*z) */
            mpz_mul (lhs, b_n, test_k);
            mpz_mul_ui (lhs, lhs, 1000);
            mpz_sub (lhs, lhs, b_n);
            mpz_mul_ui (rhs, z, 1000);
            if (mpz_cmp (lhs, rhs) <= 0)
              check_it_out = 1;
          }

        if (check_it_out)
          {
            mpz_mul (k_b_n, b_n, test_k);
            mpz_sub (diff, z, k_b_n);
            mpz_abs (abs_diff, diff);

            if (mpz_cmp_ui (abs_diff, max_diff) <= 0)
        {
          /* make sure k and c are relatively prime */
          if (mpz_gcd_ui (NULL, test_k, mpz_get_ui (diff)) == 1)
            {
              /* we are done!!! */
              *k = mpz_get_d (test_k);
              *b = mpz_get_ui (base);
              *n = exp;
              *c = mpz_get_si (diff);

              ret = 1;
              goto end_kbnc;
            }
          else
            {
              *b = 0;
              ret = 0;
              goto end_kbnc;
            }
        }
          }

        mpz_divexact (b_n, b_n, base);
        exp--;
    }
  }

  /* if we get down here, then we couldn't find a representation k*b^n + c */
 end_kbnc:
  mpz_clear (diff);
  mpz_clear (abs_diff);
  mpz_clear (b_n);
  mpz_clear (k_b_n);
  mpz_clear (lhs);
  mpz_clear (rhs);
  mpz_clear (test_k);
  mpz_clear (max_k);
  mpz_clear (base);
  mpz_clear (base_min);
  mpz_clear (base_max);

  return ret;
}

/* This function searches for a nice representation of z
   We are trying to see if z = k*b^n + c
   Some examples that we can find:
   "3^218+5123"
   "(199*3^218+5123)/(2*17*587*1187)"
   "(199*3^218 + 5123)/2/17/587/1187" */
int
kbnc_str (double *k, unsigned long *b, unsigned long *n, signed long *c,
          char *z, mpz_t num)
{
  int i = 0;
  int total = 0;
  char strk[11];
  char strb[11];
  char strn[11];
  char strc[11];
  mpz_t tmp;

  /* make sure we have a place to put our results */
  if (k == NULL || b == NULL || n == NULL || c == NULL || z == NULL)
    return 0;

  *b = 0;

  for (i = 0; i < strlen(z); i++)
    {
      if (z[i] == '(' || z[i] == '{' || z[i] == '[')
        continue;

      /* check to see if the input is k*b^n+c */
      total = sscanf (z+i, "%10[0-9]*%10[0-9]^%10[0-9]%*[ +]%10[0-9]",
                      strk, strb, strn, strc);
      if (total == 4)
        {
          *k = (double) strtoul (strk, NULL, 10);
          *b = strtoul (strb, NULL, 10);
          *n = strtoul (strn, NULL, 10);
          *c = strtol (strc, NULL, 10);
          break;
        }

      /* check to see if the input is k*b^n-c */
      total = sscanf (z+i, "%10[0-9]*%10[0-9]^%10[0-9]%*[ -]%10[0-9]",
                      strk, strb, strn, strc);
      if (total == 4)
        {
          *k = (double) strtoul (strk, NULL, 10);
          *b = strtoul (strb, NULL, 10);
          *n = strtoul (strn, NULL, 10);
          *c = strtol (strc, NULL, 10);
          *c *= -1;
          break;
        }

      /* check to see if the input is b^n+c (k = 1) */
      total = sscanf (z+i, "%10[0-9]^%10[0-9]%*[ +]%10[0-9]",
                      strb, strn, strc);
      if (total == 3)
        {
          *k = 1.0;
          *b = strtoul (strb, NULL, 10);
          *n = strtoul (strn, NULL, 10);
          *c = strtol (strc, NULL, 10);
          break;
        }

    /* check to see if the input is b^n-c (k = 1) */
      total = sscanf (z+i, "%10[0-9]^%10[0-9]%*[ -]%10[0-9]",
                      strb, strn, strc);
      if (total == 3)
        {
          *k = 1.0;
          *b = strtoul (strb, NULL, 10);
          *n = strtoul (strn, NULL, 10);
          *c = strtol (strc, NULL, 10);
          *c *= -1;
          break;
        }
      break;
    }

  /* first, check to see if we found a k*b^n+c */
  if (*b)
    {
      /* if we did, make sure that (k*b^n+c) is divisible by num */
      mpz_init_set_ui (tmp, *b);
      mpz_pow_ui (tmp, tmp, *n);
      mpz_mul_ui (tmp, tmp, (unsigned long) *k);
      if (*c >= 0)
        mpz_add_ui (tmp, tmp, *c);
      else
        mpz_sub_ui (tmp, tmp, (*c * -1));

      if (mpz_divisible_p (tmp, num))
        return 1;
    }

  /* set b to zero so users have a second way to know we didn't find k,b,n,c */
  *b = 0;

  /* if we get here, we didn't find a formula k*b^n+c for z */
  return 0;
}

int 
gw_ecm_stage1 (mpz_t f, curve *P, mpmod_t modulus, 
	       double B1, double *B1done, mpz_t go, double gw_k,
               unsigned long gw_b, unsigned long gw_n, signed long gw_c)
{
  ecm_uint gw_B1done = *B1done;
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
           "Using gwnum_ecmStage1(%.0f, %d, %d, %d, %.0f, %ld)\n",
           gw_k, gw_b, gw_n, gw_c, B1, gw_B1done);

  /* Copy x, z and A values from modular representation to 
     plain integers */

  /* Allocate enough memory for any residue (mod k*b^n+c) for x, z */
  mpz_init2 (gw_x, mpz_sizeinbase (modulus->orig_modulus, 2) + 64);
  mpz_init2 (gw_z, mpz_sizeinbase (modulus->orig_modulus, 2) + 64);
  mpz_init (gw_A);

  /* mpres_get_z always produces non-negative integers */
  mpres_get_z (gw_x, P->x, modulus);
  mpres_get_z (gw_z, P->y, modulus);
  mpres_get_z (gw_A, P->A, modulus);

  /* gwnum_ecmStage1() wants long int pointers for size_x, size_z, 
     so copy them into long int vars */
  siz_x = SIZ(gw_x);
  siz_z = SIZ(gw_z);

  /* George Woltman says that the gwnum library can handle k values up to 49
     or 50 bits long, and the maximum c value is +/-8388607 */
  ASSERT_ALWAYS (gw_k == rint (gw_k)); /* check that k is an integer */
  ASSERT_ALWAYS (1.0 <= gw_k && gw_k <= 562949953421312.0);
  ASSERT_ALWAYS (-8388607 <= gw_c && gw_c <= 8388607);
#if GMP_NUMB_BITS <= 32
  youpi = gwnum_ecmStage1_u32 (gw_k, gw_b, gw_n, gw_c, 
      PTR(modulus->orig_modulus), ABSIZ(modulus->orig_modulus), 
      B1, &gw_B1done, PTR(gw_A), ABSIZ(gw_A), 
      PTR(gw_x), &siz_x, PTR(gw_z), &siz_z, NULL, 0);
#else /* contributed by David Cleaver */
  youpi = gwnum_ecmStage1_u64 (gw_k, gw_b, gw_n, gw_c,
      PTR(modulus->orig_modulus), ABSIZ(modulus->orig_modulus),
      B1, &gw_B1done, PTR(gw_A), ABSIZ(gw_A),
      PTR(gw_x), &siz_x, PTR(gw_z), &siz_z, NULL, 0);
#endif
  
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

/* Here is a list of gwnum return codes. */
/* In the case of 2 or 5, we should continue on and let gmp-ecm */
/* do stage 1, instead of throwing an error and quitting */
/* #define ES1_SUCCESS           0 *//* Success, but no factor */
/* #define ES1_FACTOR_FOUND      1 *//* Success, factor found */
/* #define ES1_CANNOT_DO_IT      2 *//* This k,b,n,c cannot be handled */
/* #define ES1_MEMORY            3 *//* Out of memory */
/* #define ES1_INTERRUPT         4 *//* Execution interrupted */
/* #define ES1_CANNOT_DO_QUICKLY 5 *//* Requires 3-multiply reduction */
/* #define ES1_HARDWARE_ERROR    6 *//* An error was detected, most */
  if (youpi == ES1_CANNOT_DO_IT || youpi == ES1_CANNOT_DO_QUICKLY)
  {
    outputf (OUTPUT_VERBOSE, 
           "Notice: Did not use gwnum_ecmStage1(%.0f, %d, %d, %d, %.0f, %ld)\n",
           gw_k, gw_b, gw_n, gw_c, B1, gw_B1done);
    youpi = ECM_NO_FACTOR_FOUND;
    goto end_of_gwecm;
  }

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
