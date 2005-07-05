/* 
  Interface code for George Woltman's gwnum library
  
  Copyright 2004 Paul Zimmermann and Alexander Kruppa.
  
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
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
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

#define USE_FCONV
/* #define DEBUG */

/* This prototype should be in gwnum.h */
int gwnum_ecmStage1 (
        double,                         /* K in K*B^N+C */
        unsigned long,                  /* B in K*B^N+C */
        unsigned long,                  /* N in K*B^N+C */
        signed long,                    /* C in K*B^N+C */
        unsigned long *,                /* Number to factor */
        unsigned long,
        unsigned long,                  /* Stage 1 bound */
        unsigned long *,                /* Stage 1 that is already done */
        unsigned long *,                /* A - caller derives it from sigma */
        unsigned long,
        unsigned long *,                /* X value of point */
        unsigned long *,
        unsigned long *,                /* Z value of point */
        unsigned long *,
        int     (*)(void),              /* Ptr to proc that returns TRUE */
                                        /* if user interrupts processing */
        unsigned long options);


typedef struct {
  char *bigword;
  double *weights;
  double *invweights;
  unsigned long *offsets;
} gwconvplan_struct;
typedef gwconvplan_struct gwconvplan_t;

extern unsigned long FFTLEN;	/* The FFT size we are using */
extern double KARG;		/* K in K*2^N+C */
extern unsigned long BITS_PER_WORD;
extern unsigned long GW_ZEROWORDSLOW;
extern int ZERO_PADDED_FFT;

static gwnum g1 = NULL, g2 = NULL;
static gwconvplan_t *gwplan = NULL;

void __gxx_personality_v0()
{
  exit (EXIT_FAILURE);
}

static int 
sgn (int i)
{
  if (i == 0)
    return 0;
  return i > 0 ? 1 : -1;
}

static gwconvplan_t *
make_conv_plan ()
{
  gwconvplan_t *plan;
  unsigned int i;
  unsigned long offset, lastoff;
  const unsigned int planlen = FFTLEN + 1; 
  /* One extra or preloading values in loop runs off end of array */

  plan = (gwconvplan_t *) malloc (sizeof (gwconvplan_t));
  plan->bigword = (char *) malloc (planlen);
  plan->weights = (double *) malloc (planlen * sizeof (double));
  plan->invweights = (double *) malloc (planlen * sizeof (double));
  plan->offsets = (unsigned long *) malloc (planlen * sizeof (long));

  lastoff = 0;

  for (i = 0; i < FFTLEN; i++)
    {
      plan->bigword[i] = is_big_word (i) ? 1 : 0;
      offset = addr_offset (FFTLEN, i);
      plan->offsets[i] = offset - lastoff;
      plan->weights[i] = gwfft_weight (i);;
      plan->invweights[i] = gwfft_weight_inverse (i);
      lastoff = offset;
    }
  plan->bigword[i] = 0;
  plan->offsets[i] = 0;
  plan->weights[i] = 0.;
  plan->invweights[i] = 0.;

  return plan;
}

#ifdef DEBUG
static void 
gwdump (gwnum g)
{
  unsigned long i, offset;
  for (i = 0; i < FFTLEN; i++)
    {
      offset = addr_offset (FFTLEN, i);
      printf ("%ld: [%ld] = %.1f  ", i, offset, *(double *)((unsigned long)g + offset));
    }
  printf ("\n");
}
#endif

/* Convert a number stored as an array of dwords to the gwnum FFT format. */
/* Assumes 16 bits per FFT element and all DWT weights == 1.0 */

void 
Fdwordstogw (
	unsigned long *e1,
	long e1len,
	gwnum	g,
	gwconvplan_t *plan)
{
	unsigned long binval;
	unsigned long nextoff;
	long i, lim;
	long value, carry;
	unsigned long *op;

/* Convert the dword array to FFT format */

	if (FFTLEN == 0)
		return;

	binval = 0;
	if (e1len)
	  {
            /* e1 has e1len dwords of valid data left */
            binval = *e1++;
            e1len--;
          }
        lim = e1len * 2 - 2;
        if ((signed) FFTLEN < lim)
          lim = FFTLEN;

	carry = 0;
	op = plan->offsets;
	nextoff = (unsigned long) g;

	/* Unrolled, does 2 input dwords -> 4 FFT elements per loop.
	   A pass decreases e1len by 2, so needs lim+2 <= e1len. */
	for (i = 0; i < lim; i += 4) {
		nextoff += *op++;
		value = (long)((unsigned short) binval);
		binval >>= 16;
		
#ifdef BALANCED
		value += carry;
		if (value > 0x7FFF) {
			value -= 0x10000;
			binval++;
		}
#endif
		
		*(double *) nextoff = (double) value;
		nextoff += *op++;

		value = binval;
#ifdef BALANCED
		carry = (value > 0x7FFF) ? 1 : 0;
		if (value > 0x7FFF)
			value -= 0x10000;
#endif
		binval = *e1++;
		
		*(double *) nextoff = (double) value;
		nextoff += *op++;


		value = (long)((unsigned short) binval);
		binval >>= 16;
		
#ifdef BALANCED
		value += carry; 
		if (value > 0x7FFF) {
			value -= 0x10000;
			binval++;
		}
#endif
		
		*(double *) nextoff = (double) value;
		nextoff += *op++;

		value = binval;
#ifdef BALANCED
		carry = (value > 0x7FFF) ? 1 : 0;
		if (value > 0x7FFF)
			value -= 0x10000;
#endif
		binval = *e1++;
		
		*(double *) nextoff = (double) value;
	}
	
	e1len -= i / 2;
	
	for ( ; i < (signed) FFTLEN - 2; i += 2) {
		nextoff += *op++;
		value = (long)((unsigned short) binval);
		binval >>= 16;

#ifdef BALANCED
		value += carry; 
		if (value > 0x7FFF) {
			value -= 0x10000;
			binval++;
		}
#endif
		
		*(double *) nextoff = (double) value;
		nextoff += *op++;

		value = binval;
#ifdef BALANCED
		carry = (value > 0x7FFF) ? 1 : 0;
		if (value > 0x7FFF)
		  value -= 0x10000;
#endif
		if (e1len) {
			binval = *e1++;
			e1len--;
		} else
			binval = 0;
		
		*(double *) nextoff = (double) value;
	}

	nextoff += *op++;
	value = (long)(unsigned short) binval;
	binval >>= 16;

	value += carry; 
#ifdef BALANCED
	if (value > 0x7FFF) {
		value -= 0x10000;
		binval++;
	}
#endif
	
	*(double *) nextoff = (double) value;
/*	printf ("Fdwordstogw: next-to-last value written: %d\n", value); */
	nextoff += *op++;
	
	/* This is the last value written. Don't generate a carry out */
	*(double *) nextoff = (double) binval;
/*	printf ("Fdwordstogw: Last value %d written to offset %d\n", 
		binval, nextoff - (unsigned long) g); */

	((long *) g)[-1] = 0;	/* Clear needs-normalize counter */
	((long *) g)[-7] = 0;	/* Clear has been FFTed flag */
}

/* Convert a gwnum value to an array of dwords. Assumes a DWT with 16 bits
   per FFT element */

void 
Fgwtodwords (
        gwnum   g,
        unsigned long *outptrparam,
        int *outlen,
        int outalloc,
        gwconvplan_t *plan)
{
	long val, binval;
	double dval;
	int i;
	unsigned long nextoff;
	unsigned long *op, *outptr;
	
	val = 0;
	
	if ((signed) FFTLEN > outalloc * 2)
	  {
	    fprintf (stderr, "Fgwtodwords: output array has only %d dwords allocated, but need %ld\n",
	             outalloc, FFTLEN / 2);
	    exit (EXIT_FAILURE);
	  }
	
	op = plan->offsets;
	nextoff = (unsigned long) g;
	nextoff += *op++;
	outptr = outptrparam;
	
	for (i = 0; i < (signed) FFTLEN; i += 2) {
		dval = *(double *) nextoff;
		nextoff += *op++;

#if 0
		if (dval < -0.5)
			val += (long) (dval - 0.5);
		else
			val += (long) (dval + 0.5);
#else
		val += (long) dval;
#endif		
		binval = (long)(unsigned short) val;
		val >>= 16; /* -1 if val was negative, 0 otherwise */
		
		dval = *(double *) nextoff;
		nextoff += *op++;

#if 0
                if (dval < -0.5)
                        val += (long) (dval - 0.5);
                else
                        val += (long) (dval + 0.5);
#else
		val += (long) dval;
#endif                
                binval += ((unsigned long)(unsigned short) val) << 16;
                val >>= 16;
		*outptr++ = binval;
	}
	
	*outlen = FFTLEN / 2;
	/* Output carry */
	if (val != 0 && val != -1) {
		if ((*outlen) + 1 >= outalloc) {
			fprintf (stderr, "Fgwtodwords: carry, but out of allocated space\n");
			exit (EXIT_FAILURE);
		}
		*outptr++ = val;
		outlen++;
	}
	
	/* Set length, normalize */
	while (*outlen && outptrparam[(*outlen) - 1] == 0) 
		(*outlen)--;
	
	if (val < 0) {
		int len = *outlen;
		/* Negative. :(   Flip bits */
		for (i = 0; i < len; i++) 
			outptrparam[i] = ~outptrparam[i];
		/* Normalize again */
		while (*outlen && outptrparam[(*outlen)-1] == 0) 
			(*outlen)--;
		/* Add 1 */
		for (i = 0; i < *outlen; i++)
			if (++outptrparam[i])
				break;
		if (i == *outlen && outptrparam[i] == 0) {
			(*outlen)++;
			outptrparam[i] = 1;
		}
		*outlen = -(*outlen);
	}
}

void 
Fgwinit (int Fermat)
{
  /* Init GW routines with 16 bits per FFT element */
  guessCpuType ();
  guessCpuSpeed ();
  gwsetup (1., 2, Fermat, 1, Fermat / 16);
  g1 = gwalloc();
  g2 = gwalloc();
  gwplan = make_conv_plan ();
}

void 
Fgwclear ()
{
  free (gwplan);
  gwplan = NULL;
  gwfree (g2);
  g2 = NULL;
  gwfree (g1);
  g1 = NULL;
  gwdone ();
}

void
Fgwmul (mpz_t R, mpz_t S1, mpz_t S2)
{
  int sgnS1, sgnS2;
#ifdef DEBUG
  mpz_t t;
#endif

  if (g1 == NULL || g2 == NULL || gwplan == NULL)
    {
      fprintf (stderr, "Fgwmul: g1, g2 or gwplan not initialised\n");
      exit (EXIT_FAILURE);
    }

  /* Make both input arguments nonnegative, remember original sign */
  sgnS1 = mpz_sgn (S1);
  sgnS2 = mpz_sgn (S2);
  if (sgnS1 < 0)
    mpz_neg (S1, S1);
  /* Warning: S1 == S2 is possible! */
  if (sgnS2 < 0 && S1 != S2)
    mpz_neg (S2, S2);

#ifdef DEBUG
  /* Test if converting to gwnum and back is identity */
  mpz_init (t);
  _mpz_realloc (t, FFTLEN / 2);

  /* Test for S1 */
#ifdef USE_FCONV
  Fdwordstogw (PTR(S1), ABSIZ(S1), g1, gwplan);
  Fgwtodwords (g1, PTR(t), &SIZ(t), ALLOC(t), gwplan);
#else
  binarytogw (PTR(S1), ABSIZ(S1), g1);
  SIZ(t) = gwtobinary (g1, PTR(t), ALLOC(t));
#endif
  if (mpz_cmp (S1, t) != 0)
    {
      gmp_printf (
#ifdef USE_FCONV
                  "Fdwordstogw/Fgwtodwords"
#else
                  "binarytogw/gwtobinary"
#endif
                  " is not identity for S1 = %Zd\nResult: %Zd\ng = ", S1, t);
      gwdump (g1);
    }

  /* Test for S2 */
#ifdef USE_FCONV
  Fdwordstogw (PTR(S2), ABSIZ(S2), g1, gwplan);
  Fgwtodwords (g1, PTR(t), &SIZ(t), ALLOC(t), gwplan);
#else
  binarytogw (PTR(S2), ABSIZ(S2), g1);
  SIZ(t) = gwtobinary (g1, PTR(t), ALLOC(t));
#endif
  if (mpz_cmp (S2, t) != 0)
    {
      gmp_printf (
#ifdef USE_FCONV
                  "Fdwordstogw/Fgwtodwords"
#else
                  "binarytogw/gwtobinary"
#endif
                  "is not identity for S2 = %Zd\nResult: %Zd\ng = ", S2, t);
      gwdump (g1);
    }
/*  printf ("Identitiy test passed\n"); */
  mpz_clear (t);
#endif

  if (ALLOC(R) < (signed) FFTLEN / 2)
    _mpz_realloc (R, FFTLEN / 2);

#ifdef USE_FCONV
  Fdwordstogw (PTR(S1), ABSIZ(S1), g1, gwplan);
#else
  binarytogw (PTR(S1), ABSIZ(S1), g1);
#endif
  if (S1 == S2)
    gwsquare (g1);
  else
    {
#ifdef USE_FCONV
      Fdwordstogw (PTR(S2), ABSIZ(S2), g2, gwplan);
#else
      binarytogw (PTR(S2), ABSIZ(S2), g2);
#endif
      gwmul (g2, g1); /* g1 = g1 * g2, g2 left FFT'd */
    }
  if (gw_test_for_error () != 0)
    {
      fprintf (stderr, "Fgwmul: gwsquare/gwmul reports error %d\n", 
               gw_test_for_error ());
      exit (EXIT_FAILURE);
    }
#ifdef USE_FCONV
  Fgwtodwords (g1, PTR(R), &SIZ(R), ALLOC(R), gwplan);
#else
  SIZ(R) = gwtobinary (g1, PTR(R), ALLOC(R));
#endif
  ASSERT (SIZ(R) >= 0);

  /* Undo sign change */
  /* S1 == R or S2 == R is possible! Don't change sign of R here */
  if (sgnS1 < 0 && S1 != R)
    mpz_neg (S1, S1);
  if (sgnS2 < 0 && S1 != S2 && S2 != R)
    mpz_neg (S2, S2);
  
  /* If one source operand was negative and the other positive, 
     change sign of R (also happens if both were 0, but changing sign has 
     no effect then */
  if (sgnS1 + sgnS2 == 0)
    mpz_neg (R, R);
}

int 
gw_ecm_stage1 (mpz_t f, curve *P, mpmod_t modulus, 
	       double B1, double *B1done, mpz_t go)
{
  const double gw_k = 1.;
  const unsigned long gw_b = 2;
  unsigned long gw_n = abs (modulus->bits);
  signed long gw_c = sgn (modulus->bits);
  unsigned long gw_B1done = *B1done;
  unsigned long siz_x, siz_z; /* Size of gw_x and gw_y as longs */
  mpz_t gw_x, gw_z, gw_A;
  int youpi = ECM_NO_FACTOR_FOUND;

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

  /* Temporarily pause the use of GWNUM by mpmod.c */
  mpmod_pausegw (modulus);

  /* gwnum_ecmStage1() wants long int pointers for size_x, size_z, 
     so copy them into long int vars */
  siz_x = SIZ(gw_x);
  siz_z = SIZ(gw_z);
  
  youpi = gwnum_ecmStage1 (gw_k, gw_b, gw_n, gw_c, 
      PTR(modulus->orig_modulus), ABSIZ(modulus->orig_modulus), 
      B1, &gw_B1done, PTR(gw_A), ABSIZ(gw_A), 
      PTR(gw_x), &siz_x, PTR(gw_z), &siz_z, NULL, 0);
  
  SIZ(gw_x) = siz_x;
  SIZ(gw_z) = siz_z;

  /* Resume use of GWNUM by mpmod.c */
  mpmod_contgw (modulus);

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
      mpres_set_ui (P->y, 1, modulus);
    }

end_of_gwecm:

  return youpi;
}

