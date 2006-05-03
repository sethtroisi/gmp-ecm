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

#ifdef TESTDRIVE
/* For rdtscll() */
#include <asm/msr.h>
#endif

#if 0 || 0 && defined (TESTDRIVE)
  #define DEBUG
#endif

typedef struct {
  unsigned long offset1, offset2, offset3;
  unsigned long *offsets;
  int layout;
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

static inline int
dtoi (double d)
{
    int i;
    __asm__ ("fistpl %0"
	 : "=m" (i) : "t" (d) : "st");
    return i;
}


/*

>  You can reduce the offsets table size by a factor of 2,4,or 8 if you want 
>  to do some loop unrolling.  The distance between elements 0 and 1 is the 
>  same as 2 and 3, 4 and 5, etc.  So, to unroll by 4:
>  
>  Remember the offsets for 1,2,3  and 0,4,8,12,...  In your code, output 
>  to 0,0+1,0+2,0+3,loop,4,4+1,4+2,4+3,etc.  You do more additions but the 
>  offsets table is 1/4 of the size.
>  
>  Regards,
>  George 

*/

static gwconvplan_t *
make_conv_plan ()
{
  gwconvplan_t *plan;
  unsigned int i;
  long int offset, lastoff;
  const int stride = 2; /* We only store every stride-th offset */
  const char *nosimplemsg = "Fgw.c, make_conv_plan(): No simple layout";

  plan = (gwconvplan_t *) malloc (sizeof (gwconvplan_t));
  if (plan == NULL)
    {
      outputf (OUTPUT_ERROR, 
	       "make_conv_plan: could not allocate memory for plan\n");
      return NULL;
    }
  
  /* One extra or preloading values in loop runs off end of array */
  plan->offsets = (unsigned long *) malloc ((FFTLEN / stride + 1) 
					    * sizeof (long));
  if (plan->offsets == NULL)
    {
      outputf (OUTPUT_ERROR, 
	       "make_conv_plan: could not allocate memory for offsets\n");
      return NULL;
    }

  plan->offset1 = addr_offset (FFTLEN, 1);

  lastoff = 0;
  for (i = 0; i < FFTLEN / stride; i++)
    {
      offset = addr_offset (FFTLEN, stride * i);
      outputf (OUTPUT_TRACE, "Fgw.c, make_conv_plan(): offset[%u] = %ld\n", 
	       i, offset);
      plan->offsets[i] = offset - lastoff;
      lastoff = offset;
    }

  plan->offsets[i] = 0;

  plan->layout = 2; /* Standard layout, uses the offsets table */

  /* See if this is a simple layout */

  if (plan->offset1 != 16)
  {
      outputf (OUTPUT_DEVVERBOSE, "%s, plan->offset1 = %ld != 16\n", 
	       nosimplemsg, plan->offset1);
      goto exitplan;
  }

  lastoff = plan->offsets[0];
  if (plan->offsets[0] != 0)
  {
      outputf (OUTPUT_DEVVERBOSE, "%s, plan->offsets[0] = %ld != 0\n", 
	       nosimplemsg, plan->offsets[0]);
      goto exitplan;
  }

  for (i = 1; i < FFTLEN / stride / 2; i++)
      if (plan->offsets[i] != 16 * stride)
      {
	  outputf (OUTPUT_DEVVERBOSE, "%s, plan->offsets[%d] = %ld\n", 
		   nosimplemsg, i, plan->offsets[i]);
	  goto exitplan;
      }

  if ((lastoff = addr_offset (FFTLEN, stride * i)) != 8)
  {
      outputf (OUTPUT_DEVVERBOSE, "%s, addr_offset (%d) = %ld != 8\n", 
	       nosimplemsg, i, lastoff);
      goto exitplan;
  }

  for (i = FFTLEN / stride / 2 + 1; i < FFTLEN / stride; i++)
      if (plan->offsets[i] != 16 * stride)
      {
	  outputf (OUTPUT_DEVVERBOSE, "%s, plan->offsets[%d] = %ld\n", 
		   nosimplemsg, i, plan->offsets[i]);
	  goto exitplan;
      }

  /* addr_offset (FFTLEN, i) = 
     (i < FFTLEN / 2) ? 16 * i : 
                        8 + 16 * (i - FFTLEN / 2) */
  plan->layout = 1;
  free (plan->offsets);
  outputf (OUTPUT_DEVVERBOSE, 
	   "Fgw.c, make_conv_plan(): Using simple layout\n");

 exitplan:

  return plan;
}

#ifdef DEBUG
static void 
gwdump (char *s, gwnum g)
{
  unsigned long i, j, offset;
  const unsigned long gl = (unsigned long) g;
  
  if (s != NULL)
      printf (s);

  for (i = 0; i < FFTLEN; i++)
    {
      offset = addr_offset (FFTLEN, i);
      if (*(double *)(gl + offset) == 0.)
      {
	  for (j = i + 1; j < FFTLEN; j++)
	      if (*(double *)(gl + addr_offset (FFTLEN, j)) != 0.)
		  break;
	  if (j - i > 2)
	  {
	      printf ("%ld ... %ld: 0.  ", i, j);
	      i = j;
	      continue;
	  }
      }
      printf ("%ld: [%ld] = %.1f  ", i, offset, *(double *)(gl + offset));
    }
  printf ("\n");
  fflush (stdout);
}
#endif

/* Convert a number stored as an array of dwords to the gwnum FFT format.
   Assumes 16 bits per FFT element, all DWT weights == 1.0 and simple
   FFT data layout (0, 16, 32, ..., 8*FFTLEN - 16, 8, 24, ..., 8*FFTLEN - 8)
*/

static void 
Fdwordstogw_simple (
	unsigned long *e1param,
	long e1len,
	gwnum	g)
{
	double *nextoff;
	unsigned long i, lim, maxlen;
	double dval1, dval2;
	unsigned short *e1;
	unsigned short *e2;

/* Convert the dword array to FFT format */

	if (FFTLEN == 0)
		return;

	nextoff = (double *) g;
	e1 = (unsigned short *) e1param;
	e2 = ((unsigned short *) e1param) + FFTLEN / 2;

	/* maxlen: how many words can be read from e1 at most */
	ASSERT (e1len >= 0);
	maxlen = 2 * (unsigned long) e1len;
	if (maxlen > FFTLEN)
	    maxlen = FFTLEN;

	((long *) g)[-1] = 0;	/* Clear needs-normalize counter */
	((long *) g)[-7] = 0;	/* Clear has been FFTed flag */

	/* In the first pass, we write FFT elements 
	   0, FFTLEN/2, 1, FFTLEN/2+1, ..., where FFT elem. FFTLEN/2+i 
	   contains data from e1[FFTLEN/4+i/2]. So we can write at most
	   i < maxlen - FFTLEN/2 <= FFTLEN/2 */
	if (maxlen > FFTLEN / 2)
	{
	    lim = maxlen - FFTLEN / 2;
	    /* Fill low and high halves with data from e1. This results in
	       sequential write accesses to the gwnum */
	    for (i = lim; i > 0; i--)
	    {
		dval1 = (double) *e1++;
		dval2 = (double) *e2++;
		*(nextoff++) = dval1;
		*(nextoff++) = dval2;
	    }

	    if (lim == FFTLEN / 2)
		return;
	    i = lim;
	}
	else 
	    i = 0;

	/* Now lim words from the low half of e1 and lim (=all) words
	   from the high half of e1 have been used. Use the remaining
	   words from the low half. There are min(maxlen, FFTLEN/2) words
	   in the low half. */

	/* Fill up the low half (addresses == 0 (mod 16) of the FFT array 
	   with data from e1, the high half (== 8 (mod 16) with 0. */
	lim = maxlen;
	if (FFTLEN / 2 < lim)
	    lim = FFTLEN / 2;
	dval2 = 0.;

	for ( ; i < lim; i++)
	{
	    dval1 = (double) *e1++;
	    *(nextoff++) = dval1;
	    *(nextoff++) = dval2;
	}

	/* Fill any rest with 0. */
	for ( ; i < FFTLEN / 2; i++) 
	{
	    *(nextoff++) = dval2;
	    *(nextoff++) = dval2;
	}

#ifdef DEBUG
	e1 = (unsigned short *) e1param;
	for (i = 0; i < FFTLEN; i++)
	{
	    unsigned short ival = (i < maxlen) ? *(e1++) : 0;
	    dval1 = *(double *) ((unsigned int) g + addr_offset (FFTLEN, i));
	    if ((unsigned long) dval1 != ival)
		printf ("Fdwordstogw_simple: g[%ld] = %f, e1[%ld] = %u\n",
			i, dval1, i, ival);
	}
#endif
}

static void 
Fdwordstogw2 (
	unsigned long *e1param,
	long e1len,
	gwnum	g,
	gwconvplan_t *plan)
{
	unsigned long nextoff, off1;
	long i, lim;
	unsigned long *op;
	double dval1, dval2;
	unsigned short *e1 = (unsigned short *) e1param;

/* Convert the dword array to FFT format */

	if (FFTLEN == 0)
		return;

	op = plan->offsets;
	nextoff = (unsigned long) g + *op++;
        lim = e1len / 2;
        if ((signed) FFTLEN / 4 < lim)
	    lim = FFTLEN;    /* We can process 4*lim FFT elements in the 
				first loop */
	off1 = plan->offset1;

	/* Unrolled, does 2 input dwords -> 4 FFT elements per loop. */
	for (i = lim; i > 0; i--)
	{
		dval1 = (double) *e1++;
		dval2 = (double) *e1++;
		*(double *) nextoff = dval1;
		*(double *) (nextoff + off1) = dval2;
		nextoff += *op++;

		dval1 = (double) *e1++;
		dval2 = (double) *e1++;
		*(double *) nextoff = dval1;
		*(double *) (nextoff + off1) = dval2;
		nextoff += *op++;
	}

	i = lim * 4; /* We have done i FFT elements so far */
	lim = e1len - lim * 2; /* e1 has lim dwords worth of data left */
	((long *) g)[-1] = 0;	/* Clear needs-normalize counter */
	((long *) g)[-7] = 0;	/* Clear has been FFTed flag */
	
	for ( ; i < (signed) FFTLEN; i += 2) 
	{
	    if (lim > 0) 
	    {
		dval1 = (double) *e1++;
		dval2 = (double) *e1++;
		lim--;
	    } else {
		dval1 = 0.;
		dval2 = 0.;
	    }
	    
	    *(double *) nextoff = dval1;
	    *(double *) (nextoff + off1) = dval2;
	    nextoff += *op++;
	}
}

/* Convert a gwnum value to an array of dwords. Assumes a DWT with 16 bits
   per FFT element */

static int 
Fgwtodwords_simple (
        gwnum   g,
        unsigned long *outptrparam,
        int *outlen,
        const int outalloc)
{
    long vals[2];
#ifndef HAVE_SSE2
    double dval1, dval2;
#endif
    int i;
    double *nextoff;
    unsigned long binval1, binval2, *outptr1, *outptr2;
    
    if ((signed) FFTLEN / 2 > outalloc)
    {
	outputf (OUTPUT_ERROR, "Fgwtodwords_simple: output array has only %d "
		 "dwords allocated, but needs %ld\n", outalloc, FFTLEN / 2);
	return ECM_ERROR;
    }

#if 0
    gwdump ("Fgwtodwords_simple: g = ", g);
#endif
	
    nextoff = (double *) g;
    outptr1 = outptrparam;
    outptr2 = outptrparam + FFTLEN / 4;
    vals[0] = 0;
    vals[1] = 0;
    
    for (i = FFTLEN / 4; i > 0; i--) 
    {
#ifdef HAVE_SSE2
	long long lldummy;

	asm ("cvtpd2pi (%2), %1 \n"
	     "paddd %0, %1 \n"
             "movq %1, %0" : "+m" (vals[0]), "=y" (lldummy) : "r" (nextoff));
	nextoff += 2;
#else
	dval1 = *(nextoff++);
	dval2 = *(nextoff++);
	vals[0] += dtoi (dval1);
	vals[1] += dtoi (dval2);
#endif

	binval1 = (unsigned long)(unsigned short) vals[0];
	binval2 = (unsigned long)(unsigned short) vals[1];
	vals[0] >>= 16;
	vals[1] >>= 16;
	
#ifdef HAVE_SSE2
	asm ("cvtpd2pi %2, %1 \n"
	     "paddd %0, %1 \n"
             "movq %1, %0" : "+m" (vals[0]), "=y" (lldummy) : "m" (*nextoff));
	nextoff += 2;
#else
	dval1 = *(nextoff++);
	dval2 = *(nextoff++);
	vals[0] += dtoi (dval1);
	vals[1] += dtoi (dval2);
#endif
	binval1 |= ((unsigned long)(unsigned short) vals[0]) << 16;
	binval2 |= ((unsigned long)(unsigned short) vals[1]) << 16;
	*outptr1++ = binval1;
	*outptr2++ = binval2;
	vals[0] >>= 16;
	vals[1] >>= 16;
    }

    /* Add val1 to first element of high half of output */
    binval1 = *outptr1;
    binval1 += vals[0];
    *outptr1++ = binval1;

    /* Do we have a carry? */
    if (((vals[0] < 0) && (binval1 < (unsigned long) -vals[0])) || 
	((vals[0] > 0) && (binval1 < (unsigned long) vals[0])))
    {
	/* Add to e1 */
	do
	{
	    binval1 = *outptr1;
	    binval1++;
	    *outptr1++ = binval1;
	} while  (binval1 == 0 && outptr1 < outptr2);

	/* Still carry? Add to val2 */
	vals[1] += (binval1 == 0) ? 1 : 0;
    }

    *outlen = FFTLEN / 2;
    /* Output carry */
    if (vals[1] != 0 && vals[1] != -1) 
    {
	if ((*outlen) + 1 > outalloc) 
	{
	    outputf (OUTPUT_ERROR, "Fgwtodwords_simple: carry = %lu, but out "
	             "of allocated space (%d allocated)\n", vals[1], outalloc);
	    return ECM_ERROR;
	}
	*outptr2++ = vals[1];
	(*outlen)++;
    }
    
#if 0
    printf ("Fgwtodwords_simple: ");
    for (i = 0; i < FFTLEN / 2; i++)
	printf ("%lx ", outptrparam[i]);
    printf ("\n");
    fflush (stdout);
#endif

    /* Set length, normalize */
    while (*outlen && outptrparam[(*outlen) - 1] == 0) 
	(*outlen)--;
    
    if (vals[1] < 0) 
    {
	/* Negative. :(   Flip bits */
	int len = *outlen;
#ifdef DEBUG
	printf ("Fgwtodwords_simple: negating\n");
#endif
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
    return 0;
}

static void 
Fgwtodwords2 (
        gwnum   g,
        unsigned long *outptrparam,
        int *outlen,
        int outalloc,
        gwconvplan_t *plan)
{
	long val, binval;
	double dval1, dval2;
	int i;
	unsigned long nextoff, off1;
	unsigned long *op, *outptr;
	
	if ((signed) FFTLEN / 2 > outalloc)
	  {
	    outputf (OUTPUT_ERROR, "Fgwtodwords: output array has only %d "
	             "dwords allocated, but needs %ld\n", 
	             outalloc, FFTLEN / 2);
	    exit (EXIT_FAILURE);
	  }
	
	op = plan->offsets;
	nextoff = (unsigned long) g + *op++;
	off1 = plan->offset1;
	outptr = outptrparam;
	val = 0;
	
	/* Not unrolled, does 2 FFT elements per loop */
	for (i = FFTLEN / 2; i > 0; i--) 
	{
	    dval1 = *(double *) nextoff;
	    dval2 = *(double *) (nextoff + off1);
	    val += (long) dval1;
	    binval = (long)(unsigned short) val;
	    nextoff += *op++;
	    val >>= 16;
	    val += (long) dval2;
	    binval |= ((unsigned long)(unsigned short) val) << 16;
	    val >>= 16;
	    *outptr++ = binval;
	}
	
	*outlen = FFTLEN / 2;
	/* Output carry */
	if (val != 0 && val != -1) {
		if ((*outlen) + 1 > outalloc) {
			outputf (OUTPUT_ERROR, "Fgwtodwords: carry, but out "
			         "of allocated space\n");
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

static inline void
Fmpztogw (gwnum g, mpz_t S, gwconvplan_t *plan)

{
    ASSERT(mpz_sgn (S) >= 0);
    if (plan->layout == 1)
	Fdwordstogw_simple (PTR(S), ABSIZ(S), g);
    else if (plan->layout == 2)
	Fdwordstogw2 (PTR(S), ABSIZ(S), g, plan);
    else
	binarytogw (PTR(S), ABSIZ(S), g2);
}

static inline int 
Fgwtompz (mpz_t R, gwnum g, gwconvplan_t *plan)
{
    if (plan->layout == 1)
	return Fgwtodwords_simple (g, PTR(R), &SIZ(R), ALLOC(R));
    else if (plan->layout == 2)
    {
	Fgwtodwords2 (g, PTR(R), &SIZ(R), ALLOC(R), plan);
	return ECM_NO_FACTOR_FOUND;
    } else {
	SIZ(R) = gwtobinary (g, PTR(R), ALLOC(R));
	return ECM_NO_FACTOR_FOUND;
    }
}

void 
Fgwinit (int Fermat)
{
  char buf[128];
  /* Init GW routines with 16 bits per FFT element */
  guessCpuType ();
  guessCpuSpeed ();
  gwsetup (1., 2, Fermat, 1, (Fermat + 15) / 16);
  gwfft_description (buf);
  outputf (OUTPUT_VERBOSE, "Using %s\n", buf);
  outputf (OUTPUT_VERBOSE, "Modulus is %s\n", gwmodulo_as_string ());
  g1 = gwalloc ();
  g2 = gwalloc ();
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

/* Compute R = S1 * S2 (mod F_m) using George Woltman's GWNUM code
   R == S1 || R == S2 is permissible */

void
Fgwmul (mpz_t R, mpz_t S1, mpz_t S2)
{
  int sgnS1, sgnS2;
#ifdef DEBUG
  mpz_t t;
#endif

  if (g1 == NULL || g2 == NULL || gwplan == NULL)
    {
      outputf (OUTPUT_ERROR, "Fgwmul: g1, g2 or gwplan not initialised\n");
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
  mpz_init2 (t, (FFTLEN + 1) * 16);

  /* Test identity for S1 */

  Fmpztogw (g1, S1, gwplan);
  Fgwtompz (t, g1, gwplan);

  if (mpz_cmp (S1, t) != 0)
    {
      gmp_printf ("Fmpztogw/Fgwtompz is not identity for "
		  "S1 = %ZX\nResult: %ZX\n", S1, t);
      gwdump ("g1 = ", g1);
    }

  /* Test identity for S2 */
  Fmpztogw (g1, S2, gwplan);
  Fgwtompz (t, g1, gwplan);
  
  if (mpz_cmp (S2, t) != 0)
    {
      gmp_printf ("Fmpztogw/Fgwtompz is not identity for "
		  "S2 = %ZX\nResult: %ZX\n", S2, t);
      gwdump ("g1 = ", g1);
    }
  mpz_clear (t);
#endif /* DEBUG */

  if (ALLOC(R) < (signed) FFTLEN / 2 + 1)
    {
      /* WARNING: _mpz_realloc() does not keep the value! Since we allow 
         R == S1 || R == S2, we cannot use _mpz_realloc() here */
      mpz_realloc2 (R, (FFTLEN * 16) + mp_bits_per_limb);
    }

  Fmpztogw (g1, S1, gwplan);
  
  if (S1 == S2)
  {
#ifdef DEBUG
      if (test_verbose (OUTPUT_TRACE))
	  gwdump ("Fgwmul: before gwsquare, g1 = ", g1);
#endif
      
      gwsquare (g1);
      
#ifdef DEBUG
      if (test_verbose (OUTPUT_TRACE))
	  gwdump ("Fgwmul: after gwsquare, g1 = ", g1);
#endif
  } else {
      Fmpztogw (g2, S2, gwplan);

#ifdef DEBUG
      if (test_verbose (OUTPUT_TRACE))
      {
	  gwdump ("Fgwmul: before gwmul,\ng1 = ", g1);
	  gwdump ("g2 = ", g2);
      }
#endif

      gwmul (g2, g1); /* g1 = g1 * g2, g2 left FFT'd */

#ifdef DEBUG
      if (test_verbose (OUTPUT_TRACE))
	  gwdump ("Fgwmul: after gwmul, g1 = ", g1);
#endif
    }
  if (gw_test_for_error () != 0)
    {
      outputf (OUTPUT_ERROR, "Fgwmul: gwsquare/gwmul reports error %d\n", 
               gw_test_for_error ());
      exit (EXIT_FAILURE);
    }

  Fgwtompz (R, g1, gwplan);

#ifdef DEBUG
  if (mpz_sgn (R) < 0)
      outputf (OUTPUT_DEVVERBOSE, "Fgwmul: R is negative\n");
#endif

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

#ifdef TESTDRIVE
int
main (int argc, char **argv)
{
  unsigned int Fermat = 1024, i;
  int error = 0, debug = 0;
  mpz_t t1, t2, t3, Fm;
  unsigned long long tsc1, tsc2;
  unsigned long besttime1, besttime2, besttime3;

  ECM_STDOUT = stdout;

  while (argc > 1)
  {
      if (strcmp (argv[1], "-d") == 0)
	  debug = 1;      
      else if (strcmp (argv[1], "-v") == 0)
	  set_verbose (get_verbose () + 1);
      else
	  break;

      argc--;
      argv++;
  }

  if (argc > 1)
  {
      i = atoi (argv[1]);
      if (9 <= i && i <= 20)
	  Fermat = 1 << i;
  }

  Fgwinit (Fermat);
  mpz_init2 (Fm, Fermat + 1);
  mpz_init2 (t1, 2 * Fermat);
  mpz_init2 (t2, 2 * Fermat);
  mpz_init2 (t3, 2 * Fermat);
  mpz_set_ui (Fm, 1);
  mpz_mul_2exp (Fm, Fm, Fermat);
  mpz_add_ui (Fm, Fm, 1); /* Fm = 2^Fermat + 1 */
  mpz_set_ui (t1, 3);
  
  for (i = 1; i <= 10; i++)
  {
      mpz_mul (t3, t1, t1);
      mpz_mod (t3, t3, Fm);
      
      if (debug)
      {
	  Fgwmul (t2, t1, t1);
	  mpz_mod (t2, t2, Fm);
	  if (mpz_cmp (t2, t3) != 0)
	  {
	      fprintf (stderr, "Error: 3^(2*3^%u) %% (2^%u+1) t2 != t3\n", 
		       i - 1, Fermat);
	      if (1)
		  gmp_fprintf (stderr, "t2 = %ZX,\nt3 = %ZX\n", t2, t3);
	      error = 1;
	      mpz_set (t2, t3);
	  }
      }
      
      mpz_mul (t3, t3, t1);
      mpz_mod (t3, t3, Fm);
      
      if (debug)
      {
	  Fgwmul (t2, t2, t1);
	  mpz_mod (t2, t2, Fm);
	  if (mpz_cmp (t2, t3) != 0)
	  {
	      fprintf (stderr, "Error: 3^(3^%u) %% (2^%u+1) t2 != t3\n", 
		       i, Fermat);
	      if (1)
		  gmp_fprintf (stderr, "t1 = %ZX,\nt2 = %ZX\n", t2, t3);
	      error = 1;
	  }
      }
      
      mpz_set (t1, t3);
  }
  
  /* Do timing test. We'll take the best of 100 runs  */
  mpz_set (t1, t2);
  besttime1 = -1;
  besttime2 = -1;
  besttime3 = -1;
  for (i = 0; i < 100; i++)
  {
      rdtscll (tsc1);
      Fmpztogw (g1, t1, gwplan);
      rdtscll (tsc2);
      tsc2 -= tsc1;
      if ((unsigned long) tsc2 < besttime1)
	  besttime1 = (unsigned long) tsc2;
      
      rdtscll (tsc1);
      gwsquare (g1);          
      rdtscll (tsc2);
      tsc2 -= tsc1;
      if ((unsigned long) tsc2 < besttime3)
	  besttime3 = (unsigned long) tsc2;
      
      rdtscll (tsc1);
      Fgwtompz (t1, g1, gwplan);
      rdtscll (tsc2);
      tsc2 -= tsc1;
      if ((unsigned long) tsc2 < besttime2)
	  besttime2 = (unsigned long) tsc2;

      mpz_mod (t1, t1, Fm);
  }

  printf ("Best time Fmpztogw() for 2^%d+1: %lu ticks, "
	  "%f per FFT element\n", 
	  Fermat, besttime1, (double) besttime1 / FFTLEN);
  
  printf ("Best time Fgwtompz() for 2^%d+1: %lu ticks, "
	  "%f per FFT element\n", 
	  Fermat, besttime2, (double) besttime2 / FFTLEN);
  
  printf ("Best time gwsqr() for 2^%d+1: %lu ticks, "
	  "%f per FFT element\n", 
	  Fermat, besttime3, (double) besttime3 / FFTLEN);
  
  mpz_clear (t3);
  mpz_clear (t2);
  mpz_clear (t1);
  mpz_clear (Fm);
  Fgwclear ();

  return error;
}
#endif
