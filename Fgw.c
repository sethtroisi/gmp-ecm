/* 
  Interface code for George Woltman's gwnum library
  
  Copyright 2004 Paul Zimmermann and Alexander Kruppa.
  
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

#ifndef HAVE_GWNUM

/* This file does nothing at all if HAVE_GWNUM is not defined */

#else

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "ecm-gmp.h"
#include "ecm-impl.h"
#define ADD_UNDERSCORES
#include "gwnum.h"
#include "cpuid.h"

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

static gwconvplan_t *
make_conv_plan ()
{
  gwconvplan_t *plan;
  unsigned int i;
  unsigned long offset, lastoff;

  plan = (gwconvplan_t *) malloc (sizeof (gwconvplan_t));
  plan->bigword = (char *) malloc (FFTLEN);
  plan->weights = (double *) malloc (FFTLEN * sizeof (double));
  plan->invweights = (double *) malloc (FFTLEN * sizeof (double));
  plan->offsets = (unsigned long *) malloc (FFTLEN * sizeof (long));

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

  return plan;
}

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

/* Convert a number stored as an array of dwords to the gwnum FFT format. */
/* Only for b^n+c numbers, i.e. k==1 */

static void 
dwordstogw (
	unsigned long *e1,
	long e1len,
	gwnum	g,
	gwconvplan_t *plan)
{
	unsigned long i, mask1, mask2;
	int	bits1, bits2, bits_in_next_binval = 0;
	unsigned long binval, carry, nextoff;

	ASSERTG (KARG == 1);

/* Convert the dword array to FFT format */

	bits1 = BITS_PER_WORD;
	bits2 = bits1 + 1;
	mask1 = (1L << bits1) - 1;
	mask2 = (1L << bits2) - 1;
	if (e1len) {binval = *e1++; e1len--; bits_in_next_binval = 32;}
	else binval = 0;
	carry = 0;
	nextoff = (unsigned long) g;
	for (i = 0; i < FFTLEN - 1; i++) {
		int	big_word, bits;
		long	value, mask;
		big_word = plan->bigword[i];
		nextoff += plan->offsets[i];
		mask = big_word ? mask2 : mask1;
		bits = bits1 + big_word;
		value = (binval & mask) + carry;
		carry = 0;
		if (value > (mask >> 1)) {
			value = value - mask - 1;
			carry = 1;
		}
		*(double *) nextoff = (double) value * plan->weights[i];

		binval >>= bits;
		if (e1len == 0) continue;
		if (bits_in_next_binval == 0) {
			e1++; 
			e1len--; 
			bits_in_next_binval = 32;
			if (e1len == 0) continue;
		} else if (bits_in_next_binval < bits) {
				binval |= (*e1 >> (32 - bits_in_next_binval)) << (32 - bits);
				bits -= bits_in_next_binval;
				e1++; 
				e1len--; 
				bits_in_next_binval = 32;
				if (e1len == 0) continue;
			}
		if (bits) {
			binval |= (*e1 >> (32 - bits_in_next_binval)) << (32 - bits);
			bits_in_next_binval -= bits;
		}
	}
	nextoff += plan->offsets[i];
	*(double *) nextoff = (double) (binval + carry) * plan->weights[i];
	((long *) g)[-1] = 0;	/* Clear needs-normalize counter */
	((long *) g)[-7] = 0;	/* Clear has been FFTed flag */

}


/* Convert a gwnum value to an array of dwords */

static void 
gwtodwords (
	gwnum	gg,
	unsigned long *outptrparam,
	int *outlen,
	int outalloc,
	gwconvplan_t *plan)
{
	long	val, len;
	int	j, bits, bitsout, carry;
	unsigned long i, limit, *outptr = outptrparam, nextoff;

	ASSERTG (((long *) gg)[-7] == 0);	// Number not FFTed?
	ASSERTG (KARG == 1);

/* If this is a general-purpose mod, then only convert the needed words */
/* which will be less than half the FFT length.  If this is a zero padded */
/* FFT, then only convert a little more than half of the FFT data words. */
/* For a DWT, convert all the FFT data. */

	if (GENERAL_MOD) limit = GW_ZEROWORDSLOW + 3;
	else if (ZERO_PADDED_FFT) limit = FFTLEN / 2 + 4;
	else limit = FFTLEN;

/* Collect bits until we have all of them */

	carry = 0;
	bitsout = 0;
	*outptr = 0;
	nextoff = (unsigned long) gg;
	for (i = 0; i < limit; i++) {
		double dval;
		/* get_fft_value (gg, i, &val); */
		nextoff += plan->offsets[i];
		dval = *((double *) nextoff) * plan->invweights[i];
		if (dval < -0.5)
			val = (long) (dval - 0.5);
		else
			val = (long) (dval + 0.5);
		
		bits = BITS_PER_WORD;
		bits += plan->bigword[i];
		val += carry;
		for (j = 0; j < bits; j++) {
			*outptr >>= 1;
			if (val & 1) *outptr += 0x80000000;
			val >>= 1;
			bitsout++;
			if (bitsout == 32) {
				outptr++;
				bitsout = 0;
				outalloc--;
				if (outalloc == 0 && i + 1 < limit) {
					fprintf (stderr, "gwtodwords: not enough allocated memory in target, i = %lu, limit = %lu\n", i, limit);
					exit (EXIT_FAILURE);
				}
			}
		}
		carry = val;
	}

/* Finish outputting the last word and any carry data */

	while (bitsout || (carry != -1 && carry != 0)) {
		*outptr >>= 1;
		if (carry & 1) *outptr += 0x80000000;
		carry >>= 1;
		bitsout++;
		if (bitsout == 32) {
			if (outalloc == 0) {
				fprintf (stderr, "gwtodwords: not enough allocated memory in target\n");
				exit (EXIT_FAILURE);
			}
			outptr++;
			bitsout = 0;
			outalloc--;
		}
	}

/* Set the length */

	len = (outptr - outptrparam);
	while (len && outptrparam[len-1] == 0) len--;

/* If carry is -1, the gwnum is negative.  Ugh.  Flip the bits and sign. */

	if (carry == -1) {
	        /* printf ("gwtodwords: negative\n"); */
		for (j = 0; j < len; j++) outptrparam[j] = ~outptrparam[j];
		while (len && outptrparam[len-1] == 0) len--;
		/* iaddg (1, v); */
		outptr = outptrparam;
		do {
		  (*outptr)++;			/* Add 1 */
		} while (*outptr++ == 0);	/* and propagate carry */
		len = -len;
	}

	*outlen = len;

/* The gwnum is not guaranteed to be smaller than k*b^n+c.  Handle this */
/* possibility.  This also converts negative values to positive. */

/*	specialmodg (v); */

}


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
//	printf ("Fdwordstogw: next-to-last value written: %d\n", value);
	nextoff += *op++;
	
	/* This is the last value written. Don't generate a carry out */
	*(double *) nextoff = (double) binval;
//	printf ("Fdwordstogw: Last value %d written to offset %d\n", 
//		binval, nextoff - (unsigned long) g);

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
  int sgnchg = 0;
#ifdef DEBUG
  mpz_t t;
#endif

  if (g1 == NULL || g2 == NULL || gwplan == NULL)
    {
      fprintf (stderr, "Fgwmul: g1, g2 or gwplan not initialised\n");
      exit (EXIT_FAILURE);
    }

#ifdef DEBUG
  /* Test if converting to gwnum and back is identity */
  mpz_init (t);
  _mpz_realloc (t, FFTLEN / 2);
  Fdwordstogw (PTR(S1), ABSIZ(S1), g1, gwplan);
  Fgwtodwords (g1, PTR(t), &SIZ(t), ALLOC(t), gwplan);
  if (mpz_sgn (S1) < 0)
    mpz_neg (t, t);
  if (mpz_cmp (S1, t) != 0)
    {
      gmp_printf ("Fdwordstogw/Fgwtodwords is not identity for %Zd\nResult: %Zd\ng = ", 
                  S1, t);
      gwdump (g1);
    }
  Fdwordstogw (PTR(S2), ABSIZ(S2), g1, gwplan);
  Fgwtodwords (g1, PTR(t), &SIZ(t), ALLOC(t), gwplan);
  if (mpz_sgn (S2) < 0)
    mpz_neg (t, t);
  if (mpz_cmp (S2, t) != 0)
    {
      gmp_printf ("Fdwordstogw/Fgwtodwords is not identity for %Zd\nResult: %Zd\ng = ", 
                  S2, t);
      gwdump (g1);
    }
  mpz_clear (t);
#endif

  if (ALLOC(R) < (signed) FFTLEN / 2)
    _mpz_realloc (R, FFTLEN / 2);

  Fdwordstogw (PTR(S1), ABSIZ(S1), g1, gwplan);
  if (S1 == S2)
    gwsquare (g1);
  else
    {
      sgnchg = mpz_sgn (S1) ^ mpz_sgn (S2);
      Fdwordstogw (PTR(S2), ABSIZ(S2), g2, gwplan);
      gwmul (g2, g1); /* g1 = g1 * g2, g2 left FFT'd */
    }
  if (gw_test_for_error () != 0)
    {
      fprintf (stderr, "Fgwmul: gwsquare/gwmul reports error %ld\n", 
               gw_test_for_error ());
      exit (EXIT_FAILURE);
    }
  Fgwtodwords (g1, PTR(R), &SIZ(R), ALLOC(R), gwplan);
  if (sgnchg < 0)
    mpz_neg (R, R);
}


#endif /* endelse of #ifndef HAVE_GWNUM */
