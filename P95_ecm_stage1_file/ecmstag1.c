/* ECM stage 1 using GWNUM -- for use by GMP-ECM

  Copyright 1996-2024 Mersenne Research, Inc.

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

/**************************************************************
 *
 *	ecm.c
 *
 *	ECM stage 1 factoring program
 *
 *	Original author:  Richard Crandall - www.perfsci.com
 *	Adapted to Mersenne numbers and optimized by George Woltman
 *	Further optimizations from Paul Zimmerman's GMP-ECM program
 *	Other important ideas courtesy of Peter Montgomery.
 *
 *************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "cpuid.h"
#include "gmp.h"		// GMP library
#include "gwnum.h"
#include "math.h"
#include "memory.h"

/* PBMcL additions  for Lucas chain codes */

/* 3-bit chain start sequences. Covers all possible Lucas chains */
#define CHAIN_START_5_8_13 0x7
#define CHAIN_START_5_8_11 0x6
#define CHAIN_START_5_8_10 0x5
#define CHAIN_START_5_7    0x4
#define CHAIN_START_5_6    0x3
#define CHAIN_START_4_7    0x2
#define CHAIN_START_4_5    0x1
#define CHAIN_START_4_6    0x0 /* precludes a completely zero code */

typedef struct
{
	/* note: "parent" is the immediately previous chain element */
	uint64_t	value;			/* integer value of this chain element */
	uint8_t	comp_offset_1;	/* larger summand (summand_1) component index counting back from parent (parent = 0) */
	uint8_t	comp_offset_2;	/* smaller summand (summand_2) component index counting back from parent */
	uint8_t	dif_offset;		/* component index of (summand_1 - summand_2) counting back from parent */
								/* note: dif_offset = 0 will indicate that this is a doubled element */
} chain_element;

/* prototypes */
static void gw_max_continuation( chain_element *, uint8_t *, uint8_t );
static uint8_t gw_generate_Lucas_chain( uint64_t, uint64_t, chain_element * );

/* end PBMcL additions */

gwhandle gwdata;

/* from ecm.cpp */

struct xz {
	gwnum	x;		/* x or FFT(x) */
	gwnum	z;		/* z or FFT(z) */
};

/* This routine initializes an xz pair with two allocated gwnums */

__inline int alloc_xz ( 
        gwhandle *gwdata,			/* Returns TRUE if successful */
		struct xz *arg)
{
	arg->x = gwalloc (gwdata);
	if (arg->x == NULL) return (FALSE);
	arg->z = gwalloc (gwdata);
	if (arg->z == NULL) return (FALSE);
	return (TRUE);
}

/* This routine cleans up an xz pair with two allocated gwnums */

void gwfree_xz (
    gwhandle *gwdata,
    struct xz *arg)
{
	gwfree (gwdata, arg->z); arg->z = NULL;
	gwfree (gwdata, arg->x); arg->x = NULL;
}

/* Macro to swap to xz structs */

#define xzswap(a,b)	{ struct xz t; t = a; a = b; b = t; }

/* This routine gwcopies an xz pair */

__inline void gwcopy_xz (
    gwhandle *gwdata,
	struct xz *src,
	struct xz *dst)
{
	gwcopy (gwdata, src->x, dst->x);
	gwcopy (gwdata, src->z, dst->z);
}

/* Global variables */

giant	N = NULL;		/* Number being factored */
giant	FAC = NULL;		/* Found factor */
int	PRAC_SEARCH = 7;

gwnum	Ad4 = NULL;

struct xz A = {NULL, NULL};
struct xz B = {NULL, NULL};
struct xz C = {NULL, NULL};
struct xz T = {NULL, NULL};
struct xz scr = {NULL, NULL};

gwnum t1_es1 = NULL;
gwnum t2_es1 = NULL;
gwnum t3_es1 = NULL;

gwnum	xA_es1 = NULL;
gwnum	zA_es1 = NULL;
gwnum	xB_es1 = NULL;
gwnum	zB_es1 = NULL;
gwnum	xC_es1 = NULL;
gwnum	zC_es1 = NULL;
gwnum	xt_es1 = NULL;
gwnum	zt_es1 = NULL;

/* Bit manipulation macros */

#define bitset(a,i)	{ a[(i) >> 3] |= (1 << ((i) & 7)); }
#define bitclr(a,i)	{ a[(i) >> 3] &= ~(1 << ((i) & 7)); }
#define bittst(a,i)	(a[(i) >> 3] & (1 << ((i) & 7)))


/* Perform cleanup functions. */

void ecm_cleanup ()
{
	free (N);
	N = NULL;
	free (FAC);
	FAC = NULL;
	gwdone (&gwdata);
}


/* Determine if a number is prime */

int isPrime (
	unsigned long p)
{
	unsigned long i;
	for (i = 2; i * i <= p; i = (i + 1) | 1)
		if (p % i == 0) return (FALSE);
	return (TRUE);
}

/* Use a simple sieve to find prime numbers */

#define MAX_PRIMES	6542
static	unsigned int *primes = NULL;
static	struct sieve_info {
	uint64_t first_number;
	unsigned int bit_number;
	unsigned int num_primes;
	uint64_t start;
	char	array[4096];
} si = {0};

/* Fill up the sieve array */

void fill_sieve (void)
{
	unsigned int i, fmax;

/* Determine the first bit to clear */

	fmax = (unsigned int)
		sqrt ((double) (si.first_number + sizeof (si.array) * 8 * 2));
	for (i = si.num_primes; i < MAX_PRIMES * 2; i += 2) {
		unsigned long f, r, bit;
		f = primes[i];
		if (f > fmax) break;
		if (si.first_number == 3) {
			bit = (f * f - 3) >> 1;
		} else {
			r = (unsigned long) (si.first_number % f);
			if (r == 0) bit = 0;
			else if (r & 1) bit = (f - r) / 2;
			else bit = (f + f - r) / 2;
			if (f == si.first_number + 2 * bit) bit += f;
		}
		primes[i+1] = bit;
	}
	si.num_primes = i;

/* Fill the sieve with ones, then zero out the composites */

	memset (si.array, 0xFF, sizeof (si.array));
	for (i = 0; i < si.num_primes; i += 2) {
		unsigned int f, bit;
		f = primes[i];
		for (bit = primes[i+1]; bit < sizeof (si.array) * 8; bit += f)
			bitclr (si.array, bit);
		primes[i+1] = bit - sizeof (si.array) * 8;
	}
	si.bit_number = 0;
}

/* Start sieve by allocate a sieve info structure */

void start_sieve (
	uint64_t start)
{
	unsigned int i;

/* Remember starting point (in case its 2) and make real start odd */

	if (start < 2) start = 2;
	si.start = start;
	start |= 1;

/* See if we can just reuse the existing sieve */

	if (si.first_number &&
	    start >= si.first_number &&
	    start < si.first_number + sizeof (si.array) * 8 * 2) {
		si.bit_number = (unsigned int) (start - si.first_number) / 2;
		return;
	}

/* Initialize sieve */

	if (primes == NULL) {
		unsigned int f;
		primes = (unsigned int *)
			malloc (MAX_PRIMES * 2 * sizeof (unsigned int));
		for (i = 0, f = 3; i < MAX_PRIMES * 2; f += 2)
			if (isPrime (f)) primes[i] = f, i += 2;
	}

	si.first_number = start;
	si.num_primes = 0;
	fill_sieve ();
}

/* Return next prime from the sieve */

uint64_t sieve (void)
{
	if (si.start == 2) {
		si.start = 3;
		return (2);
	}
	for ( ; ; ) {
		unsigned int bit;
		if (si.bit_number == sizeof (si.array) * 8) {
			si.first_number += 2 * sizeof (si.array) * 8;
			fill_sieve ();
		}
		bit = si.bit_number++;
		if (bittst (si.array, bit))
			return (si.first_number + 2 * bit);
	}
}

/**************************************************************
 *
 *	Functions
 *
 **************************************************************/

uint8_t gw_generate_Lucas_chain( uint64_t prime, uint64_t chain_code, chain_element *Lchain )
{
	uint8_t code_fragment, chain_length, i, k;
	uint64_t dif;

	/* the code file starts at p = 11, so handle 2, 3, 5, or 7 separately */
	if( prime < 11 )
	{
		if(prime == 2)
		{
			chain_length = 1;
			return chain_length;
		}
		else if( prime == 3 )
		{
			chain_length = 2;
			return chain_length;
		}
		else if( prime == 5 )
		{
			Lchain[3].value = 5;
			Lchain[3].comp_offset_1 = 0;
			Lchain[3].comp_offset_2 = 1;
			Lchain[3].dif_offset = 2;

			chain_length = 3;
			return chain_length;
		}
		else if( prime == 7 )
		{
			Lchain[3].value = 4;
			Lchain[3].comp_offset_1 = 1;
			Lchain[3].comp_offset_2 = 1;
			Lchain[3].dif_offset = 0;

			Lchain[4].value = 7;
			Lchain[4].comp_offset_1 = 0;
			Lchain[4].comp_offset_2 = 1;
			Lchain[4].dif_offset = 3;

			chain_length = 4;
			return chain_length;
		}
		else
		{
			printf("ERROR: generate_Lucas_chain entered with prime = %lu < 11 but != 2, 3, 5 or 7\n", prime);
			return 0;
		}
	}

	/* first 3 bits of code give the next two or three chain components */
	code_fragment = (uint8_t)(chain_code & 0x7);
	chain_code >>= 3;
	switch( code_fragment )
	{
	case CHAIN_START_5_8_13:
		Lchain[3].value = 5;
		Lchain[3].comp_offset_1 = 0;
		Lchain[3].comp_offset_2 = 1;
		Lchain[3].dif_offset = 2;

		Lchain[4].value = 8;
		Lchain[4].comp_offset_1 = 0;
		Lchain[4].comp_offset_2 = 1;
		Lchain[4].dif_offset = 2;

		Lchain[5].value = 13;
		Lchain[5].comp_offset_1 = 0;
		Lchain[5].comp_offset_2 = 1;
		Lchain[5].dif_offset = 2;

		chain_length = 5;
		break;

	case CHAIN_START_5_8_11:
		Lchain[3].value = 5;
		Lchain[3].comp_offset_1 = 0;
		Lchain[3].comp_offset_2 = 1;
		Lchain[3].dif_offset = 2;

		Lchain[4].value = 8;
		Lchain[4].comp_offset_1 = 0;
		Lchain[4].comp_offset_2 = 1;
		Lchain[4].dif_offset = 2;

		Lchain[5].value = 11;
		Lchain[5].comp_offset_1 = 0;
		Lchain[5].comp_offset_2 = 2;
		Lchain[5].dif_offset = 1;

		chain_length = 5;
		break;

	case CHAIN_START_5_8_10:
		Lchain[3].value = 5;
		Lchain[3].comp_offset_1 = 0;
		Lchain[3].comp_offset_2 = 1;
		Lchain[3].dif_offset = 2;

		Lchain[4].value = 8;
		Lchain[4].comp_offset_1 = 0;
		Lchain[4].comp_offset_2 = 1;
		Lchain[4].dif_offset = 2;

		Lchain[5].value = 10;
		Lchain[5].comp_offset_1 = 1;
		Lchain[5].comp_offset_2 = 1;
		Lchain[5].dif_offset = 0;

		chain_length = 5;
		break;

	case CHAIN_START_5_7:
		Lchain[3].value = 5;
		Lchain[3].comp_offset_1 = 0;
		Lchain[3].comp_offset_2 = 1;
		Lchain[3].dif_offset = 2;

		Lchain[4].value = 7;
		Lchain[4].comp_offset_1 = 0;
		Lchain[4].comp_offset_2 = 2;
		Lchain[4].dif_offset = 1;

		chain_length = 4;
		break;

	case CHAIN_START_5_6:
		Lchain[3].value = 5;
		Lchain[3].comp_offset_1 = 0;
		Lchain[3].comp_offset_2 = 1;
		Lchain[3].dif_offset = 2;

		Lchain[4].value = 6;
		Lchain[4].comp_offset_1 = 1;
		Lchain[4].comp_offset_2 = 1;
		Lchain[4].dif_offset = 0;

		chain_length = 4;
		break;

	case CHAIN_START_4_7:
		Lchain[3].value = 4;
		Lchain[3].comp_offset_1 = 1;
		Lchain[3].comp_offset_2 = 1;
		Lchain[3].dif_offset = 0;

		Lchain[4].value = 7;
		Lchain[4].comp_offset_1 = 0;
		Lchain[4].comp_offset_2 = 1;
		Lchain[4].dif_offset = 3;

		chain_length = 4;
		break;

	case CHAIN_START_4_6:
		Lchain[3].value = 4;
		Lchain[3].comp_offset_1 = 1;
		Lchain[3].comp_offset_2 = 1;
		Lchain[3].dif_offset = 0;

		Lchain[4].value = 6;
		Lchain[4].comp_offset_1 = 1;
		Lchain[4].comp_offset_2 = 1;
		Lchain[4].dif_offset = 0;

		chain_length = 4;
		break;

	case CHAIN_START_4_5:
		Lchain[3].value = 4;
		Lchain[3].comp_offset_1 = 1;
		Lchain[3].comp_offset_2 = 1;
		Lchain[3].dif_offset = 0;

		Lchain[4].value = 5;
		Lchain[4].comp_offset_1 = 0;
		Lchain[4].comp_offset_2 = 3;
		Lchain[4].dif_offset = 1;

		chain_length = 4;
		break;

	default: /* should never happen */
		printf("ERROR: bad chain code start value in gw_generate_Lucas_chain = %u\n", code_fragment);
		return 0;
	}

	/* rebuild chain from code fragments */
	while( chain_code != 0 )
	{
		code_fragment = (uint8_t)( chain_code & 0xF );
		chain_code >>= 4;
		switch( code_fragment )
		{
			case 0: /* step type 1 or 4 */
			{
				i = (uint8_t)( chain_code & 0x3 );
				chain_code >>= 2;
				if( i == 3 )
				{
					 /* a rare double, but does occur */
					Lchain[ chain_length+1 ].value = 2*Lchain[ chain_length-2 ].value;
					Lchain[ chain_length+1 ].comp_offset_1 = 2;
					Lchain[ chain_length+1 ].comp_offset_2 = 2;
					Lchain[ chain_length+1 ].dif_offset = 0;
					chain_length++;
				}
				else
				{
					i = 12*(i + 1);
					gw_max_continuation( Lchain, &chain_length, i );
				}
				break;
			}
			case 1:
			{
				Lchain[ chain_length+1 ].value = Lchain[ chain_length ].value + Lchain[ chain_length-2 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 0;
				Lchain[ chain_length+1 ].comp_offset_2 = 2;
				dif = Lchain[ chain_length ].value - Lchain[ chain_length-2 ].value;
				k = 1;
				while( dif < Lchain[ chain_length-k ].value )
					k++;
//				ASSERT( dif == Lchain[ chain_length-k ].value );

				Lchain[ chain_length+1 ].dif_offset = k;
				chain_length++;
				break;
			}
			case 2:
			{
				Lchain[ chain_length+1 ].value = 2*Lchain[ chain_length-1 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 1;
				Lchain[ chain_length+1 ].comp_offset_2 = 1;
				Lchain[ chain_length+1 ].dif_offset = 0;
				chain_length++;
				break;
			}
			case 3:
			{
				i = 1;
				gw_max_continuation( Lchain, &chain_length, i );
				Lchain[ chain_length+1 ].value = Lchain[ chain_length ].value + Lchain[ chain_length-2 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 0;
				Lchain[ chain_length+1 ].comp_offset_2 = 2;
				Lchain[ chain_length+1 ].dif_offset = 1;
				chain_length++;
				break;
			}
			case 4:
			{
				i = 1;
				gw_max_continuation( Lchain, &chain_length, i );
				Lchain[ chain_length+1 ].value = 2*Lchain[ chain_length-1 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 1;
				Lchain[ chain_length+1 ].comp_offset_2 = 1;
				Lchain[ chain_length+1 ].dif_offset = 0;
				chain_length++;
				break;
			}
			case 5:
			{
				i = 2;
				gw_max_continuation( Lchain, &chain_length, i );
				Lchain[ chain_length+1 ].value = Lchain[ chain_length ].value + Lchain[ chain_length-2 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 0;
				Lchain[ chain_length+1 ].comp_offset_2 = 2;
				Lchain[ chain_length+1 ].dif_offset = 1;
				chain_length++;
				break;
			}
			case 6:
			{
				i = 2;
				gw_max_continuation( Lchain, &chain_length, i );
				Lchain[ chain_length+1 ].value = 2*Lchain[ chain_length-1 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 1;
				Lchain[ chain_length+1 ].comp_offset_2 = 1;
				Lchain[ chain_length+1 ].dif_offset = 0;
				chain_length++;
				break;
			}
			case 7:
			{
				i = 3;
				gw_max_continuation( Lchain, &chain_length, i );
				Lchain[ chain_length+1 ].value = Lchain[ chain_length ].value + Lchain[ chain_length-2 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 0;
				Lchain[ chain_length+1 ].comp_offset_2 = 2;
				Lchain[ chain_length+1 ].dif_offset = 1;
				chain_length++;
				break;
			}
			case 8:
			{
				i = 3;
				gw_max_continuation( Lchain, &chain_length, i );
				Lchain[ chain_length+1 ].value = 2*Lchain[ chain_length-1 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 1;
				Lchain[ chain_length+1 ].comp_offset_2 = 1;
				Lchain[ chain_length+1 ].dif_offset = 0;
				chain_length++;
				break;
			}
			case 9:
			{
				i = (uint8_t)( chain_code & 0x3 );
				chain_code >>= 2;
				i += 4;
				gw_max_continuation( Lchain, &chain_length, i );
				Lchain[ chain_length+1 ].value = Lchain[ chain_length ].value + Lchain[ chain_length-2 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 0;
				Lchain[ chain_length+1 ].comp_offset_2 = 2;
				Lchain[ chain_length+1 ].dif_offset = 1;
				chain_length++;
				break;
			}
			case 10:
			{
				i = (uint8_t)( chain_code & 0x3 );
				chain_code >>= 2;
				i += 4;
				gw_max_continuation( Lchain, &chain_length, i );
				Lchain[ chain_length+1 ].value = 2*Lchain[ chain_length-1 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 1;
				Lchain[ chain_length+1 ].comp_offset_2 = 1;
				Lchain[ chain_length+1 ].dif_offset = 0;
				chain_length++;
				break;
			}
			case 11:
			{
				i = (uint8_t)( chain_code & 0x3 );
				chain_code >>= 2;
				i += 8;
				gw_max_continuation( Lchain, &chain_length, i );
				Lchain[ chain_length+1 ].value = Lchain[ chain_length ].value + Lchain[ chain_length-2 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 0;
				Lchain[ chain_length+1 ].comp_offset_2 = 2;
				Lchain[ chain_length+1 ].dif_offset = 1;
				chain_length++;
				break;
			}
			case 12:
			{
				i = (uint8_t)( chain_code & 0x3 );
				chain_code >>= 2;
				i += 8;
				gw_max_continuation( Lchain, &chain_length, i );
				Lchain[ chain_length+1 ].value = 2*Lchain[ chain_length-1 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 1;
				Lchain[ chain_length+1 ].comp_offset_2 = 1;
				Lchain[ chain_length+1 ].dif_offset = 0;
				chain_length++;
				break;
			}
			case 13:
			{
				i = (uint8_t)( chain_code & 0x3 );
				chain_code >>= 2;

				switch( i )
				{
				case 0: /* chain encoding did not require this step until chain length = 44 */
					Lchain[ chain_length+1 ].value = Lchain[ chain_length ].value + Lchain[ chain_length-7 ].value;
					Lchain[ chain_length+1 ].comp_offset_1 = 0;
					Lchain[ chain_length+1 ].comp_offset_2 = 7;
					dif = Lchain[ chain_length ].value - Lchain[ chain_length-7 ].value;
					k = 1;
					while( dif < Lchain[ chain_length-k ].value )
						k++;
//					ASSERT( dif == Lchain[ chain_length-k ].value );

					Lchain[ chain_length+1 ].dif_offset = k;
					chain_length++;
					break;

				case 1:
					Lchain[ chain_length+1 ].value = Lchain[ chain_length ].value + Lchain[ chain_length-8 ].value;
					Lchain[ chain_length+1 ].comp_offset_1 = 0;
					Lchain[ chain_length+1 ].comp_offset_2 = 8;
					dif = Lchain[ chain_length ].value - Lchain[ chain_length-8 ].value;
					k = 1;
					while( dif < Lchain[ chain_length-k ].value )
						k++;
//					ASSERT( dif == Lchain[ chain_length-k ].value );

					Lchain[ chain_length+1 ].dif_offset = k;
					chain_length++;
					break;

				case 2: /* example chain 1 2 3 5 8 10 12 13 */
					Lchain[ chain_length+1 ].value = Lchain[ chain_length - 2 ].value + Lchain[ chain_length - 3 ].value;
					Lchain[ chain_length+1 ].comp_offset_1 = 2;
					Lchain[ chain_length+1 ].comp_offset_2 = 3;
					dif = Lchain[ chain_length - 2 ].value - Lchain[ chain_length - 3 ].value;
					k = 1;
					while( dif < Lchain[ chain_length-k ].value )
						k++;
//					ASSERT( dif == Lchain[ chain_length-k ].value );

					Lchain[ chain_length+1 ].dif_offset = k;
					chain_length++;
					break;

				case 3: /* can occur but don't have an example (yet) */
					Lchain[ chain_length+1 ].value = Lchain[ chain_length - 2 ].value + Lchain[ chain_length - 4 ].value;
					Lchain[ chain_length+1 ].comp_offset_1 = 2;
					Lchain[ chain_length+1 ].comp_offset_2 = 4;
					dif = Lchain[ chain_length - 2 ].value - Lchain[ chain_length - 4 ].value;
					k = 1;
					while( dif < Lchain[ chain_length-k ].value )
						k++;
//					ASSERT( dif == Lchain[ chain_length-k ].value );

					Lchain[ chain_length+1 ].dif_offset = k;
					chain_length++;
					break;

				default: /* should never happen */
					printf("ERROR: unimplemented code fragment 0xiD, i = %u\n", i);
					return 0;
				}
				break;
			}
			case 14:
			{
				i = (uint8_t)( chain_code & 0x3 );
				chain_code >>= 2;

				Lchain[ chain_length+1 ].value = Lchain[ chain_length ].value + Lchain[ chain_length-3-i ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 0;
				Lchain[ chain_length+1 ].comp_offset_2 = 3+i;
				dif = Lchain[ chain_length ].value - Lchain[ chain_length-3-i ].value;
				k = 1;
				while( dif < Lchain[ chain_length-k ].value )
					k++;
//				ASSERT( dif == Lchain[ chain_length-k ].value );

				Lchain[ chain_length+1 ].dif_offset = k;
				chain_length++;
				break;
			}
			case 15:
			{
				i = (uint8_t)( chain_code & 0x3 );
				chain_code >>= 2;

				Lchain[ chain_length+1 ].value = Lchain[ chain_length-1 ].value + Lchain[ chain_length-2-i ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 1;
				Lchain[ chain_length+1 ].comp_offset_2 = 2+i;
				dif = Lchain[ chain_length-1 ].value - Lchain[ chain_length-2-i ].value;
				k = 2;
				while( dif < Lchain[ chain_length-k ].value )
					k++;
//				ASSERT( dif == Lchain[ chain_length-k ].value );

				Lchain[ chain_length+1 ].dif_offset = k;
				chain_length++;
				break;
			}
		}
	}

	/* finish the chain with maximum continuations until prime is reached */
	i = 1;
	while( Lchain[chain_length].value < prime )
	{
		gw_max_continuation( Lchain, &chain_length, i );
	}

	if( Lchain[ chain_length ].value != prime )
		printf("ERROR: prime/prime code mismatch for p = %lu\n", prime);

	return chain_length;
}

/* extend the chain with maximum elements for i steps */

void gw_max_continuation( chain_element *Lchain, uint8_t *chain_length, uint8_t i )
{
	uint8_t k;
	uint64_t dif;

	Lchain[ *chain_length+1 ].value = Lchain[ *chain_length ].value + Lchain[ *chain_length-1 ].value;
	Lchain[ *chain_length+1 ].comp_offset_1 = 0;
	Lchain[ *chain_length+1 ].comp_offset_2 = 1;
	dif = Lchain[ *chain_length ].value - Lchain[ *chain_length-1 ].value;
	k = 2;
	while( dif < Lchain[ *chain_length-k ].value )
		k++;
//	ASSERT( dif == Lchain[ *chain_length-k ].value );

	Lchain[ *chain_length+1 ].dif_offset = k;
	(*chain_length)++;
	if( i > 1 )
	{
		for( k = 1; k < i; k++)
		{
			Lchain[ *chain_length+1 ].value = Lchain[ *chain_length ].value + Lchain[ *chain_length-1 ].value;
			Lchain[ *chain_length+1 ].comp_offset_1 = 0;
			Lchain[ *chain_length+1 ].comp_offset_2 = 1;
			Lchain[ *chain_length+1 ].dif_offset = 2;
			(*chain_length)++;
		}
	}
}



/* old ecmstag1.c routines */
#if 0


/* end of old routines */
#endif

/* updated routines copied from ecm.cpp */
/* Computes 2P=(out.x:out.z) from P=(in.x:in.z), uses the global variable Ad4. */
/* Input arguments may be in FFTed state.  Out argument can be same as input argument. */
/* Scratch xz argument can equal in but cannot equal out. */

void ell_dbl_xz_scr (
	struct xz *in,		/* Input value to double */
	struct xz *out,		/* Output value */
	struct xz *scr)		// Scratch registers (only the .x gwnum is used)
{				/* 10 or 11 FFTs, 4 adds */
	gwnum	t2, t3, t4;

	ASSERTG (scr != out);

	/* If we have extra_bits to square the output of a gwadd, then we can use an algorithm that minimizes FFTs. */
	if (square_safe (&gwdata, 1)) {
		t3 = scr->x;
		t2 = out->x;
		t4 = out->z;

		gwsetmulbyconst (&gwdata, 4);
		gwmul3 (&gwdata, in->x, in->z, t3, GWMUL_FFT_S12 | GWMUL_MULBYCONST | GWMUL_STARTNEXTFFT); /* t3 = 4*x*z */
		gwsub3o (&gwdata, in->x, in->z, t2, GWADD_DELAY_NORMALIZE);			/* Compute x - z */
		gwsquare2 (&gwdata, t2, t2, GWMUL_STARTNEXTFFT);				/* t2 = (x - z)^2 */
		gwmul3 (&gwdata, t2, Ad4, t4, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);	/* t4 = t2 * Ad4 */
		gwaddmul4 (&gwdata, t2, t3, t4, out->x, GWMUL_FFT_S23 | GWMUL_STARTNEXTFFT);	/* outx = (t1 = t2 + t3) * t4 */
		gwaddmul4 (&gwdata, t4, t3, t3, out->z, GWMUL_STARTNEXTFFT);			/* outz = (t4 + t3) * t3 */
	}

	/* This algorithm uses an extra FFT to assure that it will work even if there are no extra bits.  We've made a conscious decision to */
	/* spend an extra FFT to make ell_add as clean as possible.  Ell_adds are ten times more common than ell_dbls.  This version lets */
	/* ell_add FFT x and z without any worries. */
	else {
		int	t1_safe = addmul_safe (&gwdata, 2,0,0);

		t2 = out->x;
		t3 = out->z;
		t4 = scr->x;

		gwsquare2 (&gwdata, in->x, scr->x, GWMUL_FFT_S1);					/* x^2 */
		gwsquare2 (&gwdata, in->z, scr->z, GWMUL_FFT_S1);					/* z^2 */
		gwsetmulbyconst (&gwdata, 2);
		gwmul3 (&gwdata, in->x, in->z, t3, GWMUL_MULBYCONST);					/* t3 = 2*x*z */

		gwadd3o (&gwdata, scr->x, scr->z, t2, GWADD_DELAY_NORMALIZE);				/* x^2 + z^2 */
		gwsub3o (&gwdata, t2, t3, t2, GWADD_DELAYNORM_IF(t1_safe));				/* t2 = x^2 - 2xz + z^2 */
		gwadd3o (&gwdata, t3, t3, t3, GWADD_FORCE_NORMALIZE);					/* t3 = 4*x*z */

		gwmul3 (&gwdata, t2, Ad4, t4, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);		/* t4 = t2 * Ad4 */
		gwaddmul4 (&gwdata, t2, t3, t4, out->x, GWMUL_FFT_S23 | GWMUL_STARTNEXTFFT);		/* outx = (t1 = t2 + t3) * t4 */
		gwaddmul4 (&gwdata, t4, t3, t3, out->z, GWMUL_STARTNEXTFFT);				/* outz = (t4 + t3) * t3 */
	}
}

/* Like ell_dbl_xz_scr, but the output arguments are not partially FFTed.  The input argument is assumed to never be used again. */

void ell_dbl_xz_scr_last (
	struct xz *in,		/* Input value to double */
	struct xz *out,		/* Output value */
	struct xz *scr)		// Scratch registers (only the .x gwnum is used)
{				/* 10 or 11 FFTs, 4 adds */
	gwnum	t2, t3, t4;

	ASSERTG (scr != out);

	/* If we have extra_bits to square the output of a gwadd, then we can use an algorithm that minimizes FFTs. */
	if (square_safe (&gwdata, 1)) {
		t3 = scr->x;
		t2 = out->x;
		t4 = out->z;

		gwsetmulbyconst (&gwdata, 4);
		gwmul3 (&gwdata, in->x, in->z, t3, GWMUL_FFT_S12 | GWMUL_MULBYCONST | GWMUL_STARTNEXTFFT); /* t3 = 4*x*z */
		gwsub3o (&gwdata, in->x, in->z, t2, GWADD_DELAY_NORMALIZE);			/* Compute x - z */
		gwsquare2 (&gwdata, t2, t2, GWMUL_STARTNEXTFFT);				/* t2 = (x - z)^2 */
		gwmul3 (&gwdata, t2, Ad4, t4, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);	/* t4 = t2 * Ad4 */
		gwaddmul4 (&gwdata, t2, t3, t4, out->x, GWMUL_FFT_S23);			/* outx = (t1 = t2 + t3) * t4 */
		gwaddmul4 (&gwdata, t4, t3, t3, out->z, 0);					/* outz = (t4 + t3) * t3 */
	}

	/* This algorithm uses an extra FFT to assure that it will work even if there are no extra bits.  We've made a conscious decision to */
	/* spend an extra FFT to make ell_add as clean as possible.  Ell_adds are ten times more common than ell_dbls.  This version lets */
	/* ell_add FFT x and z without any worries. */
	else {
		int	t1_safe = addmul_safe (&gwdata, 2,0,0);

		t2 = out->x;
		t3 = out->z;
		t4 = scr->x;

		gwsquare2 (&gwdata, in->x, scr->x, GWMUL_FFT_S1);					/* x^2 */
		gwsquare2 (&gwdata, in->z, scr->z, GWMUL_FFT_S1);					/* z^2 */
		gwsetmulbyconst (&gwdata, 2);
		gwmul3 (&gwdata, in->x, in->z, t3, GWMUL_MULBYCONST);					/* t3 = 2*x*z */

		gwadd3o (&gwdata, scr->x, scr->z, t2, GWADD_DELAY_NORMALIZE);				/* x^2 + z^2 */
		gwsub3o (&gwdata, t2, t3, t2, GWADD_DELAYNORM_IF(t1_safe));				/* t2 = x^2 - 2xz + z^2 */
		gwadd3o (&gwdata, t3, t3, t3, GWADD_FORCE_NORMALIZE);					/* t3 = 4*x*z */

		gwmul3 (&gwdata, t2, Ad4, t4, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);		/* t4 = t2 * Ad4 */
		gwaddmul4 (&gwdata, t2, t3, t4, out->x, GWMUL_FFT_S23);				/* outx = (t1 = t2 + t3) * t4 */
		gwaddmul4 (&gwdata, t4, t3, t3, out->z, 0);						/* outz = (t4 + t3) * t3 */
	}
}

/* Adds Q=(in2.x:in2.z) and R=(in1.x:in1.z) and puts the result in (out.x:out.z).  Assumes that Q-R=P or R-Q=P where P=(diff.x:diff.z). */
/* Input arguments may be in FFTed format.  Out argument can be same as any of the 3 input arguments. */
/* Scratch xz argument cannot equal in1, in2, or diff. */

void ell_add_xz_scr (
	struct xz *in1,
	struct xz *in2,
	struct xz *diff,
	struct xz *out,
	struct xz *scr)		// Scratch registers
{				/* 12 FFTs, 6 adds */
	gwnum	t1, t2;
	int	options;

	ASSERTG (scr != in1 && scr != in2 && scr != diff);
	t1 = scr->z;
	t2 = scr->x;
	gwmulmulsub5 (&gwdata, in1->x, in2->z, in1->z, in2->x, t1, GWMUL_FFT_S12 | GWMUL_STARTNEXTFFT);	/* t1 = x1z2 - z1x2 */
	gwmulmulsub5 (&gwdata, in1->x, in2->x, in1->z, in2->z, t2, GWMUL_STARTNEXTFFT);			/* t2 = x1x2 - z1z2 */
	gwsquare2 (&gwdata, t1, t1, GWMUL_STARTNEXTFFT);				/* t1 = t1^2 */
	gwsquare2 (&gwdata, t2, t2, GWMUL_STARTNEXTFFT);				/* t2 = t2^2 */
	options = (diff == out) ? GWMUL_STARTNEXTFFT : GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT;
	gwmul3 (&gwdata, t2, diff->z, t2, options);					/* t2 = t2 * zdiff (will become outx) */
	gwmul3 (&gwdata, t1, diff->x, t1, options);					/* t1 = t1 * xdiff (will become outz) */
	if (out != scr) xzswap (*scr, *out);
}

// Like ell_add_xz_scr except that out is not equal any input (including diff) so that scratch register can be inferred.
// This is simply a shortcut for better readability.

#define ell_add_xz(i1,i2,dif,o)	ell_add_xz_scr(i1,i2,dif,o,o);

/* Like ell_add_xz_scr but in1, in2 and diff are assumed to not be used again.  Out argument never has its forward FFT begun. */

void ell_add_xz_last (
	struct xz *in1,
	struct xz *in2,
	struct xz *diff,
	struct xz *out,
	struct xz *scr)
{				/* 12 FFTs, 6 adds */
	gwnum	t1, t2;

	ASSERTG (scr != in1 && scr != in2 && scr != diff);
	t1 = scr->z;
	t2 = scr->x;
	gwmulmulsub5 (&gwdata, in1->x, in2->z, in1->z, in2->x, t1, GWMUL_FFT_S12 | GWMUL_STARTNEXTFFT);	/* t1 = x1z2 - z1x2 */
	gwmulmulsub5 (&gwdata, in1->x, in2->x, in1->z, in2->z, t2, GWMUL_STARTNEXTFFT);			/* t2 = x1x2 - z1z2 */
	gwsquare2 (&gwdata, t1, t1, GWMUL_STARTNEXTFFT);				/* t1 = t1^2 */
	gwsquare2 (&gwdata, t2, t2, GWMUL_STARTNEXTFFT);				/* t2 = t2^2 */
	gwmul3 (&gwdata, t2, diff->z, t2, 0);						/* t2 = t2 * zdiff (will become outx) */
	gwmul3 (&gwdata, t1, diff->x, t1, 0);						/* t1 = t1 * xdiff (will become outz) */
	if (out != scr) xzswap (*scr, *out);
}


/* Perform an elliptic multiply using an algorithm developed by Peter Montgomery.  Basically, we try to find a near optimal Lucas */
/* chain of additions that generates the number we are multiplying by.  This minimizes the number of calls to ell_dbl and ell_add. */

/* The costing function assigns an ell_dbl call a cost of 10 and an ell_add call a cost of 12. */
/* This cost estimates the number of forward and inverse transforms performed. */

#define swap(a,b)	{uint64_t t=a;a=b;b=t;}

int lucas_cost (
	uint64_t n,
	uint64_t d)
{
	uint64_t e;//, dmod3, emod3;
	unsigned long c;

	if (d >= n || d <= n/2) return (999999999);		/* Catch invalid costings */

	c = 0;
	while (n != 1) {
	    e = n - d;
	    d = d - e;

	    c += 12;

	    while (d != e) {
		if (d < e) {
			swap (d,e);
		}
//		if (d <= e + (e >> 2)) {
//			if ((dmod3 = d%3) == 3 - (emod3 = e%3)) {
//				t = d;
//				d = (d+d-e)/3;
//				e = (e+e-t)/3;
//				c += 36;
//				continue;
//			}
//			if (dmod3 == emod3 && (d&1) == (e&1)) {
//				d = (d-e) >> 1;
//				c += 22;
//				continue;
//			}
//		}
		if (100 * d <= 296 * e) {
			d = d-e;
			c += 12;
		} else if ((d&1) == (e&1)) {
			d = (d-e) >> 1;
			c += 22;
		} else if ((d&1) == 0) {
			d = d >> 1;
			c += 22;
//		} else if ((dmod3 = d%3) == 0) {
//			d = d/3-e;
//			c += 46;
//		} else if (dmod3 == 3 - (emod3 = e%3)) {
//			d = (d-e-e)/3;
//			c += 46;
//		} else if (dmod3 == emod3) {
//			d = (d-e)/3;
//			c += 46;
		} else {
			e = e >> 1;
			c += 22;
		}
	    }
	    c += 10;
	    if (d == 1) break;
	    n = d;
	    d = (uint64_t) ((double) n * 0.6180339887498948);
	}

	return (c);
}

int lucas_mul (
	struct xz *A,
	uint64_t n,
	uint64_t d,
	int	last_mul)	// TRUE if the last multiply should not use the GWMUL_STARTNEXTFFT option
{
	uint64_t e;//, dmod3, emod3;

	while (n != 1) {
	    ell_dbl_xz_scr (A, &B, &C);				/* B = 2*A, scratch reg = C */
										/* C = A (but we delay setting that up) */

	    e = n - d;
	    d = d - e;

	    // To save two gwcopies setting C=A, we handle the most common case for the first iteration of the following loop.
	    // I've only seen three cases that end up doing the gwcopies, n=3, n=11 and n=17.  With change to 2.96 there are a couple more.
	    if (e > d && 100 * e <= 296 * d) {
		swap (d, e);
		xzswap (*A, B);							/* swap A & B, thus diff C = B */
		ell_add_xz (A, &B, &B, &C);				/* B = A+B */
		xzswap (B, C);							/* C = B */
		d = d-e;
	    }
	    else if (d > e && 100 * d <= 296 * e) {
		ell_add_xz (A, &B, A, &C);				/* B = A+B */
		xzswap (B, C);							/* C = B */
		d = d-e;
	    } else {
		gwcopy (&gwdata, A->x, C.x);
		gwcopy (&gwdata, A->z, C.z);				/* C = A */
	    }

	    while (d != e) {
		if (d < e) {
			swap (d, e);
			xzswap (*A, B);
		}

		if (100 * d <= 296 * e) {					/* d <= 2.96 * e (Montgomery used 4.00) */
			ell_add_xz_scr (A, &B, &C, &C, &T);		/* B = A+B, scratch reg = T */
			xzswap (B, C);						/* C = B */
			d = d-e;
		} else if ((d&1) == (e&1)) {
			ell_add_xz_scr (A, &B, &C, &B, &T);		/* B = A+B, scratch reg = T */
			ell_dbl_xz_scr (A, A, &T);			/* A = 2*A, scratch reg = T */
			d = (d-e) >> 1;
		} else if ((d&1) == 0) {
			ell_add_xz_scr (A, &C, &B, &C, &T);		/* C = A+C, scratch reg = T */
			ell_dbl_xz_scr (A, A, &T);			/* A = 2*A, scratch reg = T */
			d = d >> 1;
		}
		else {
			ell_add_xz_scr (&B, &C, A, &C, &T);		/* C = C-B, scratch reg = T */
			ell_dbl_xz_scr (&B, &B, &T);			/* B = 2*B, scratch reg = T */
			e = e >> 1;
		}
	    }

	    if (d == 1 && last_mul) {
		ell_add_xz_last (&B, A, &C, A, &T);			/* A = A+B, scratch reg = T */
		break;
	    } else {
		ell_add_xz_scr (&B, A, &C, A, &T);			/* A = A+B, scratch reg = T */
	    }

	    n = d;
	    d = (uint64_t) ((double) n * 0.6180339887498948);
	}
	return (0);
}
#undef swap

/* Try a series of Lucas chains to find the cheapest. */
/* First try v = (1+sqrt(5))/2, then (2+v)/(1+v), then (3+2*v)/(2+v), */
/* then (5+3*v)/(3+2*v), etc.  Finally, execute the cheapest. */

__inline int lucas_cost_several (uint64_t n, uint64_t *d) {
	int	i, c, min;
	uint64_t testd;
	for (i = 0, testd = *d - PRAC_SEARCH / 2; i < PRAC_SEARCH; i++, testd++) {
		c = lucas_cost (n, testd);
		if (i == 0 || c < min) min = c, *d = testd;
	}
	return (min);
}

int ell_mul (
	struct xz *arg,
	uint64_t n,
	int	last_mul)	// TRUE if the this is the last mul in a series and we do not want to apply the GWMUL_STARTNEXTFFT option to the result
{
	unsigned long zeros;
	int	stop_reason;

	for (zeros = 0; (n & 1) == 0; zeros++) n >>= 1;

	if (n > 1) {
		int	c, min;
		uint64_t d, mind;

		mind = (uint64_t) ceil((double) 0.6180339887498948 * n);		/*v=(1+sqrt(5))/2*/
		min = lucas_cost_several (n, &mind);

		d = (uint64_t) ceil ((double) 0.7236067977499790 * n);			/*(2+v)/(1+v)*/
		c = lucas_cost_several (n, &d);
		if (c < min) min = c, mind = d;

		d = (uint64_t) ceil ((double) 0.5801787282954641 * n);			/*(3+2*v)/(2+v)*/
		c = lucas_cost_several (n, &d);
		if (c < min) min = c, mind = d;

		d = (uint64_t) ceil ((double) 0.6328398060887063 * n);			/*(5+3*v)/(3+2*v)*/
		c = lucas_cost_several (n, &d);
		if (c < min) min = c, mind = d;

		d = (uint64_t) ceil ((double) 0.6124299495094950 * n);			/*(8+5*v)/(5+3*v)*/
		c = lucas_cost_several (n, &d);
		if (c < min) min = c, mind = d;

		d = (uint64_t) ceil ((double) 0.6201819808074158 * n);			/*(13+8*v)/(8+5*v)*/
		c = lucas_cost_several (n, &d);
		if (c < min) min = c, mind = d;

		d = (uint64_t) ceil ((double) 0.6172146165344039 * n);			/*(21+13*v)/(13+8*v)*/
		c = lucas_cost_several (n, &d);
		if (c < min) min = c, mind = d;

		d = (uint64_t) ceil ((double) 0.6183471196562281 * n);			/*(34+21*v)/(21+13*v)*/
		c = lucas_cost_several (n, &d);
		if (c < min) min = c, mind = d;

		d = (uint64_t) ceil ((double) 0.6179144065288179 * n);			/*(55+34*v)/(34+21*v)*/
		c = lucas_cost_several (n, &d);
		if (c < min) min = c, mind = d;

		d = (uint64_t) ceil ((double) 0.6180796684698958 * n);			/*(89+55*v)/(55+34*v)*/
		c = lucas_cost_several (n, &d);
		if (c < min) min = c, mind = d;

		stop_reason = lucas_mul (arg, n, mind, zeros == 0 && last_mul);
		if (stop_reason) return (stop_reason);
	}
	while (zeros--) {
		if (zeros || !last_mul) ell_dbl_xz_scr (arg, arg, &scr);
		else ell_dbl_xz_scr_last (arg, arg, &scr);
	}
	return (0);
}



/* end of updated routines */

/* Test if factor divides N, return TRUE if it does.  Destroys N. */

int testFactor (
	giant	f)
{
	modg (f, N);
	return (isZero (N));
}

/* Computes the modular inverse of a number */
/* This is done using the extended GCD algorithm */
/* The GCD is returned in FAC.  Function returns FALSE */
/* if it was interrupted by an escape. */

int modinv (
	gwnum b)
{
	giant	v;

/* Convert input number to binary */

	v = popg (&gwdata.gdata, ((unsigned long) gwdata.bit_length >> 5) + 5);
	gwtogiant (&gwdata, b, v);

#ifdef MODINV_USING_GIANTS

	int	stop_reason;

/* Let the invg code use gwnum b's memory. */
/* Compute 1/v mod N */

	gwfree_temporarily (&gwdata, b);
	stop_reason = invgi (&gwdata.gdata, 0, N, v);
	gwrealloc_temporarily (&gwdata, b);
	if (stop_reason) {
		pushg (&gwdata.gdata, 1);
		return (FALSE);
	}

/* If a factor was found, save it in FAC */

	if (v->sign < 0) {
		negg (v);
		FAC = allocgiant (v->sign);
		gtog (v, FAC);
	}

/* Otherwise, convert the inverse to FFT-ready form */

	else {
		gianttogw (&gwdata, v, b);
	}

/* Use the faster GMP library to do an extended GCD which gives us 1/v mod N */

#else
	{
	mpz_t	__v, __N, __gcd, __inv;

/* Do the extended GCD */

	mpz_init (__v);
	mpz_init (__N);
	mpz_init (__gcd);
	mpz_init (__inv);
	gtompz (v, __v);
	gtompz (N, __N);
	mpz_gcdext (__gcd, __inv, NULL, __v, __N);
	mpz_clear (__v);

/* If a factor was found (gcd != 1 && gcd != N), save it in FAC */

	if (mpz_cmp_ui (__gcd, 1) && mpz_cmp (__gcd, __N)) {
		FAC = allocgiant ((int) mpz_sizeinbase (__gcd, 32));
		mpztog (__gcd, FAC);
	}

/* Otherwise, convert the inverse to FFT-ready form */

	else {
		if (mpz_sgn (__inv) < 0) mpz_add (__inv, __inv, __N);
		mpztog (__inv, v);
		gianttogw (&gwdata, v, b);
	}

/* Cleanup and return */

	mpz_clear (__gcd);
	mpz_clear (__inv);
	mpz_clear (__N);
	}
#endif

/* Clean up */

	pushg (&gwdata.gdata, 1);

/* Increment count and return */

	return (TRUE);
}

/* Takes a point (a,b) and multiplies it by a value such that b will be one */
/* If we find a factor it is returned in FAC.  Function returns FALSE if it */
/* was interrupted. */

int normalize (
	gwnum	a,
	gwnum	b)
{
	giant	g;

/* Compute the modular inverse and scale up the first input value */

	if (!modinv (b)) return (FALSE);
	if (FAC != NULL) return (TRUE);
	gwmul (&gwdata, b, a);

/* Now make sure value is less than N */

	g = popg (&gwdata.gdata, ((unsigned long) gwdata.bit_length >> 5) + 5);
	gwtogiant (&gwdata, a, g);
	modg (N, g);
	gianttogw (&gwdata, g, a);
	pushg (&gwdata.gdata, 1);

/* All done */

	return (TRUE);
}

/* common macros for 32-bit & 64-bit subroutines */

#define ES1_ALLOCATE_AND_OPEN_CODE_FILE \
\
	/* temp vars now global; speedup ~2% */ \
\
	if (!alloc_xz (&gwdata, &current_xz)) goto no_mem;\
	if (!alloc_xz (&gwdata, &A)) goto no_mem;\
	if (!alloc_xz (&gwdata, &B)) goto no_mem;\
	if (!alloc_xz (&gwdata, &C)) goto no_mem;\
	if (!alloc_xz (&gwdata, &T)) goto no_mem;\
	if (!alloc_xz (&gwdata, &scr)) goto no_mem;\
\
	Ad4 = gwalloc (&gwdata);\
	if (Ad4 == NULL) goto no_mem;\
	x = gwalloc (&gwdata);\
	if (x == NULL) goto no_mem;\
	z = gwalloc (&gwdata);\
	if (z == NULL) goto no_mem;\
\
	/* Init Lucas chain states */\
	for( i = 0; i < 16; i++ )\
	{\
    	if (!alloc_xz (&gwdata, &LCS[i])) goto no_mem;\
	}\
	using_code_file = 0;\
\
	if( B1 > *B1_done )\
	{\
		chain_code_file = fopen("Lchain_codes.dat", "r");\
		if(chain_code_file != (FILE *)NULL )\
		{\
			using_code_file = 1;\
			base_indx = 0;\
			printf("Using Lucas chain codes in (updated) gwnum\n");\
\
			/* the first 3 chain elements are always value = 1, 2, and 3, respectively */\
			Lchain[0].value = 1;\
\
			Lchain[1].value = 2;\
			Lchain[1].comp_offset_1 = 0;\
			Lchain[1].comp_offset_2 = 0;\
			Lchain[1].dif_offset = 0;\
\
			Lchain[2].value = 3;\
			Lchain[2].comp_offset_1 = 0;\
			Lchain[2].comp_offset_2 = 1;\
			Lchain[2].dif_offset = 1;\
		}\
		else\
		{\
			printf("Lchain_codes.dat file failed to open, using prac\n");\
		}\
	}\
	chain_length = 0; /* not necessary but stops a compiler warning */


#define GENERATE_LCHAIN \
\
		if(using_code_file)\
		{\
			if (prime >= 11) /* code file starts at prime = 11 */\
				dum = fread((int8_t *)&chain_code, sizeof(uint64_t), 1, chain_code_file);\
\
			if(dum || (prime < 11))\
			{\
				/* the following call doesn't need a chain code for prime < 11 */\
				chain_length = gw_generate_Lucas_chain( prime, chain_code, Lchain );\
			}\
			else /* reached end of code file; revert to use prac */\
			{\
				fclose(chain_code_file);\
				using_code_file = 0;\
\
                /* copy LCS[base_indx] x,z states over to current_xz */\
                gwcopy_xz ( &gwdata, &LCS[base_indx], &current_xz);\
\
				printf ("Reached Lchain_codes.dat EOF at p = %lu, reverting to use prac\n", prime);\
			}\
		} /* end if( using_code_file) */


#define PROCESS_PRIME \
				/* compute P = prime*P */\
				if (using_code_file)\
				{\
					/* process Lchain */\
					next_indx = (base_indx + 1) & 0xF;\
                    ell_dbl_xz_scr (&LCS[base_indx], &LCS[next_indx], &scr);\
					base_indx = next_indx;\
\
					for(i=2; i<=chain_length; i++)\
					{\
						next_indx = (base_indx + 1) & 0xF;\
						s1_indx = (base_indx - Lchain[i].comp_offset_1) & 0xF;\
						if(Lchain[i].dif_offset == 0)\
						{\
                            ell_dbl_xz_scr (&LCS[s1_indx], &LCS[next_indx], &scr);\
						}\
						else\
						{\
							s2_indx = (base_indx - Lchain[i].comp_offset_2) & 0xF;\
							dif_indx = (base_indx - Lchain[i].dif_offset) & 0xF;\
							if( (i < chain_length) || (!last_mul) )\
                            {\
	                        	ell_add_xz (&LCS[s1_indx], &LCS[s2_indx], &LCS[dif_indx], &LCS[next_indx]);\
                            }\
							else\
                            {\
								ell_add_xz_last (&LCS[s1_indx], &LCS[s2_indx], &LCS[dif_indx], &current_xz, &current_xz);\
                            }\
						}\
						base_indx = next_indx;\
					}\
				}\
				else /* use prac */\
				    res = ell_mul (&current_xz, prime, last_mul);


#define FREE_GWNUMS \
	if(using_code_file)\
		fclose(chain_code_file);\
	for( i = 15; i >= 0; i-- )\
	{\
        gwfree_xz(&gwdata, &LCS[i]);\
	}\
	gwfree (&gwdata, z);\
	gwfree (&gwdata, x);\
	gwfree (&gwdata, Ad4);\
\
    gwfree_xz(&gwdata, &scr);\
    gwfree_xz(&gwdata, &T);\
    gwfree_xz(&gwdata, &C);\
    gwfree_xz(&gwdata, &B);\
    gwfree_xz(&gwdata, &A);\
    gwfree_xz(&gwdata, &current_xz);

 

/**************************************************************
 *
 *	ECM Function for 32-bit inputs...
 *
 **************************************************************/

/* Do ECM stage 1 for GMP-ECM using gwnum library.  See gwnum.h for */
/* a detailed explanation of inputs and outputs. */

int gwnum_ecmStage1_u32 (
	double	k,			/* K in K*B^N+C */
	unsigned long b,		/* B in K*B^N+C */
	unsigned long n,		/* N in K*B^N+C */
	signed long c,			/* C in K*B^N+C */
	uint32_t *num_being_factored_array, /* Number to factor */
	unsigned long num_being_factored_array_len,
	uint64_t B1,			/* Stage 1 bound */
	uint64_t *B1_done,		/* Stage 1 that is already done */
	uint32_t *A_array,		/* A - caller derives it from sigma */
	unsigned long A_array_len,
	uint32_t *x_array,		/* X value of point */
	unsigned long *x_array_len,
	uint32_t *z_array,		/* Z value of point */
	unsigned long *z_array_len,
	int	(*stop_check_proc)(int),/* Ptr to proc that returns TRUE */
					/* if user interrupts processing */
	unsigned long options)
{
	unsigned long bits;
	uint64_t prime, next_prime;
	uint32_t word_count_32;
	int	res;
	long	reslong;
	gwnum	x, z;
    struct xz current_xz;
	uint64_t mult;
    int last_mul;

/* Lucas chain code file mods */
	uint64_t chain_code;
	uint8_t dum, chain_length;
	int32_t i;
	FILE *chain_code_file = NULL;
	chain_element Lchain[64];
	uint8_t using_code_file; /* logical */

	/* Elliptic curve states as we follow the Lucas chain */
    struct xz LCS[16];

//	gwnum LCS_x[16];
//	gwnum LCS_z[16];
	uint8_t base_indx, next_indx, s1_indx, s2_indx, dif_indx;
/* end mods */

/* Calculate an upper bound on the number of bits in the numbers we will be */
/* FFTing.  Note: We allocate 60 extra bits to handle any possible k value. */

	if (b) 
		bits = (unsigned long) (n * log ((double) b) / log ((double) 2.0)) + 60;
	else
		bits = num_being_factored_array_len * 8 * sizeof (uint32_t);

	word_count_32 = (uint32_t)((bits >> 5) + 1);

/* Setup the assembly code */

	guessCpuType ();
	gwinit (&gwdata);
	if (b)
		res = gwsetup (&gwdata, k, b, n, c);
	else if (sizeof (unsigned long) == sizeof (uint32_t))
		res = gwsetup_general_mod (&gwdata,
					   (uint32_t *) num_being_factored_array,
					   num_being_factored_array_len);
	else
		res = gwsetup_general_mod (&gwdata,
					   (uint32_t *) num_being_factored_array,
					   num_being_factored_array_len * 2);
	if (res == GWERROR_MALLOC) return (ES1_MEMORY);
	if (res) return (ES1_CANNOT_DO_IT);
	StopCheckRoutine = stop_check_proc;

/* If we cannot handle this very efficiently, let caller know it */

	if (gwdata.GENERAL_MOD && ! (options & ES1_DO_SLOW_CASE)) {
		ecm_cleanup ();
		return (ES1_CANNOT_DO_QUICKLY);
	}

/* Allocate memory; open Lchain_codes.dat file (if it exists) */

	ES1_ALLOCATE_AND_OPEN_CODE_FILE

/* Turn the input number we are factoring into a giant.  Either use the */
/* number we were passed in or calculate k*b^n+c */

	N = allocgiant ((bits >> 5) + 1);
	if (N == NULL) goto no_mem;
	if (num_being_factored_array != NULL && num_being_factored_array_len) {
		giantstruct tmp;
		tmp.sign = num_being_factored_array_len;
		tmp.n = (uint32_t *) num_being_factored_array;
		while (tmp.sign && tmp.n[tmp.sign-1] == 0) tmp.sign--;
		gtog (&tmp, N);
	} else {
		ultog (b, N);
		power (N, n);
		dblmulg (k, N);
		iaddg (c, N);
	}

/* Convert the input A value to a gwnum.  For extra speed we precompute */
/* A * 4 and FFT that value. */

	binarytogw (&gwdata, A_array, A_array_len, Ad4);
	gwaddsmall (&gwdata, Ad4, 2);	/* Compute A+2 */
	modinv (Ad4);
	if (FAC != NULL) goto bingo;

	dbltogw (&gwdata, 4.0, x);	/* For extra speed, precompute 4 / (A+2) */
	gwmul (&gwdata, x, Ad4);
	gwfft (&gwdata, Ad4, Ad4);	/* Even more speed, save FFT of Ad4 */

/* Convert the input x value to a gwnum */

	binarytogw (&gwdata, x_array, *x_array_len, current_xz.x);

/* Convert the input z value to a gwnum.  If the input z value was not */
/* given, then assume z is one. */

	if (z_array != NULL && z_array_len != NULL && *z_array_len)
		binarytogw (&gwdata, z_array, *z_array_len, current_xz.z);
	else
		dbltogw (&gwdata, 1.0, current_xz.z);

/* Output a startup message */

//	{
//		char	fft_desc[100];
//		gwfft_description (fft_desc);
//		sprintf (buf, "Using %s\n", fft_desc);
//		OutputStr (buf);
//	}

/* Do ECM stage 1 */

	if( B1_done == NULL ) goto error; /* should never happen */

    gwfft (&gwdata, current_xz.x, current_xz.x);
    gwfft (&gwdata, current_xz.z, current_xz.z);

    last_mul = 0;	

	/* do 2's separately */
	for (mult = 2; mult <= B1; mult *= 2)
		if (mult > *B1_done)
			ell_dbl_xz_scr (&current_xz, &current_xz, &scr);

    if (using_code_file)
        gwcopy_xz ( &gwdata, &current_xz, &LCS[base_indx]);

	start_sieve (3);

	for (prime = sieve (); prime <= B1; prime = next_prime) {

		next_prime = sieve ();

		GENERATE_LCHAIN

/* Apply as many powers of prime as long as prime^n <= B1 */
/* MEMUSED: 3 gwnums (x, z, AD4) + 10 for ell_mul, + 16 for Lucas chain states */

		for (mult = prime; mult <= B1; mult *= prime)
		{
            if (next_prime > B1)
                if( mult*prime > B1) last_mul = 1;

			if (mult > *B1_done)
			{
				PROCESS_PRIME
			}
		}

/* Check for errors */

		if (gw_test_for_error (&gwdata)) goto error;

/* Check for interrupt.  If one occurs return normalized x OR x,z pair. */

		if (stop_check_proc != NULL && (*stop_check_proc)(0)) {
			*B1_done = prime;

			if (z_array == NULL) {
				StopCheckRoutine = NULL;
				normalize (current_xz.x, current_xz.z);
				if (FAC != NULL) goto bingo;
				reslong = gwtobinary (&gwdata, current_xz.x, x_array, word_count_32);
				if (reslong < 0) goto error;
				*x_array_len = reslong;
			}

			else {
				reslong = gwtobinary (&gwdata, current_xz.x, x_array, word_count_32);
				if (reslong < 0) goto error;
				*x_array_len = reslong;

				reslong = gwtobinary (&gwdata, current_xz.z, z_array, word_count_32);
				if (reslong < 0) goto error;
				*z_array_len = reslong;
			}

			FREE_GWNUMS
			ecm_cleanup ();
			return (ES1_INTERRUPT);
		}
	}
	*B1_done = B1;

/* Normalize the x value OR return the x,z pair */

	if (z_array == NULL) {
		StopCheckRoutine = NULL;
		normalize (current_xz.x, current_xz.z);
		if (FAC != NULL) goto bingo;
		reslong = gwtobinary (&gwdata, current_xz.x, x_array, word_count_32);
		if (reslong < 0) goto error;
		*x_array_len = reslong;
	} else {
		reslong = gwtobinary (&gwdata, current_xz.x, x_array, word_count_32);
		if (reslong < 0) goto error;
		*x_array_len = reslong;

		reslong = gwtobinary (&gwdata, current_xz.z, z_array, word_count_32);
		if (reslong < 0) goto error;
		*z_array_len = reslong;
	}

/* Free memory and return */

	FREE_GWNUMS
	ecm_cleanup ();
	return (ES1_SUCCESS);

/* Print a message if we found a factor! */

bingo:	//printf ("ECM found a factor\n");
	if (!testFactor (FAC)) goto error;
	gianttogw (&gwdata, FAC, x);
	reslong = gwtobinary (&gwdata, x, x_array, word_count_32);
	if (reslong < 0) goto error;
	*x_array_len = reslong;
	if (z_array != NULL) {
		z_array[0] = 1;
		*z_array_len = 1;
	}
	FREE_GWNUMS
	ecm_cleanup ();
	return (ES1_FACTOR_FOUND);

/* Return a hardware error occurred code */

error:	FREE_GWNUMS
	ecm_cleanup ();
	return (ES1_HARDWARE_ERROR);

/* Return out-of-memory error */

no_mem:	ecm_cleanup ();
	return (ES1_MEMORY);
}

/**************************************************************
 *
 *	ECM Function for 64-bit inputs...
 *
 **************************************************************/

/* Do ECM stage 1 for GMP-ECM using gwnum library.  See gwnum.h for */
/* a detailed explanation of inputs and outputs. */

int gwnum_ecmStage1_u64 (
	double	k,			/* K in K*B^N+C */
	unsigned long b,		/* B in K*B^N+C */
	unsigned long n,		/* N in K*B^N+C */
	signed long c,			/* C in K*B^N+C */
	uint64_t *num_being_factored_array, /* Number to factor */
	unsigned long num_being_factored_array_len,
	uint64_t B1,			/* Stage 1 bound */
	uint64_t *B1_done,		/* Stage 1 that is already done */
	uint64_t *A_array,	/* A - caller derives it from sigma */
	unsigned long A_array_len,
	uint64_t *x_array,	/* X value of point */
	unsigned long *x_array_len,
	uint64_t *z_array,	/* Z value of point */
	unsigned long *z_array_len,
	int	(*stop_check_proc)(int),/* Ptr to proc that returns TRUE */
					/* if user interrupts processing */
	unsigned long options)
{
	unsigned long bits;
	uint64_t prime, next_prime;
	uint32_t word_count_64;
	int	res;
	long	reslong;
	gwnum	x, z;
    struct xz current_xz;
	uint64_t mult;
    int last_mul;

	/* Lucas chain code file mods */
	uint64_t chain_code;
	uint8_t dum, chain_length;
	int32_t i;
	FILE *chain_code_file = NULL;
	chain_element Lchain[64];
	uint8_t using_code_file; /* logical */

	/* Elliptic curve states as we follow the Lucas chain */
    struct xz LCS[16];

//	gwnum LCS_x[16];
//	gwnum LCS_z[16];
	uint8_t base_indx, next_indx, s1_indx, s2_indx, dif_indx;
	/* end mods */

/* Calculate an upper bound on the number of bits in the numbers we will be */
/* FFTing.  Note: We allocate 60 extra bits to handle any possible k value. */

	if (b) 
		bits = (unsigned long) (n * log ((double) b) / log ((double) 2.0)) + 60;
	else
		bits = num_being_factored_array_len * 8 * sizeof (uint64_t);

	word_count_64 = (uint32_t)((bits >> 6) + 1);

/* Setup the assembly code */

	guessCpuType ();
	gwinit (&gwdata);

	if (b)
		res = gwsetup (&gwdata, k, b, n, c);
	else
		res = gwsetup_general_mod_64 (&gwdata,
					      (uint64_t *) num_being_factored_array,
					      num_being_factored_array_len);

	if (res == GWERROR_MALLOC) return (ES1_MEMORY);
	if (res) return (ES1_CANNOT_DO_IT);
	StopCheckRoutine = stop_check_proc;

/* If we cannot handle this very efficiently, let caller know it */

	if (gwdata.GENERAL_MOD && ! (options & ES1_DO_SLOW_CASE)) {
		ecm_cleanup ();
		return (ES1_CANNOT_DO_QUICKLY);
	}

/* Allocate memory; open Lchain_codes.dat file (if it exists) */

	ES1_ALLOCATE_AND_OPEN_CODE_FILE

/* Turn the input number we are factoring into a giant.  Either use the */
/* number we were passed in or calculate k*b^n+c */

	N = allocgiant ((bits >> 5) + 1);
	if (N == NULL) goto no_mem;
	if (num_being_factored_array != NULL && num_being_factored_array_len) {
		int	i;
		for (i = 0; i < (int) (num_being_factored_array_len * 2); i += 2) {
			N->n[i] = (uint32_t) num_being_factored_array[i/2];  /* bottom half of the 64-bit value */
			N->n[i+1] = (uint32_t) (num_being_factored_array[i/2] >> 32); /* top half of the 64-bit value */
		}
		N->sign = num_being_factored_array_len * 2;
		while (N->sign && N->n[N->sign-1] == 0) N->sign--;
	} else {
		ultog (b, N);
		power (N, n);
		dblmulg (k, N);
		iaddg (c, N);
	}

/* Convert the input A value to a gwnum.  For extra speed we precompute */
/* A * 4 and FFT that value. */

	binary64togw (&gwdata, A_array, A_array_len, Ad4);
	gwaddsmall (&gwdata, Ad4, 2);	/* Compute A+2 */
	modinv (Ad4);
	if (FAC != NULL) goto bingo;

	dbltogw (&gwdata, 4.0, x);	/* For extra speed, precompute 4 / (A+2) */
	gwmul (&gwdata, x, Ad4);
	gwfft (&gwdata, Ad4, Ad4);	/* Even more speed, save FFT of Ad4 */

/* Convert the input x value to a gwnum */

	binary64togw (&gwdata, x_array, *x_array_len, current_xz.x);

/* Convert the input z value to a gwnum.  If the input z value was not */
/* given, then assume z is one. */

	if (z_array != NULL && z_array_len != NULL && *z_array_len)
		binary64togw (&gwdata, z_array, *z_array_len, current_xz.z);
	else
		dbltogw (&gwdata, 1.0, current_xz.z);

/* Output a startup message */

//	{
//		char	fft_desc[100];
//		gwfft_description (fft_desc);
//		sprintf (buf, "Using %s\n", fft_desc);
//		OutputStr (buf);
//	}

/* Do ECM stage 1 */

	if( B1_done == NULL ) goto error; /* should never happen */

    gwfft (&gwdata, current_xz.x, current_xz.x);
    gwfft (&gwdata, current_xz.z, current_xz.z);
	
    last_mul = 0;	

	/* do 2's separately */
	for (mult = 2; mult <= B1; mult *= 2)
		if (mult > *B1_done)
			ell_dbl_xz_scr (&current_xz, &current_xz, &scr);

    if (using_code_file)
        gwcopy_xz ( &gwdata, &current_xz, &LCS[base_indx]);

	start_sieve (3);

	for (prime = sieve (); prime <= B1; prime = next_prime) {

		next_prime = sieve ();

		GENERATE_LCHAIN

/* Apply as many powers of prime as long as prime^n <= B1 */
/* MEMUSED: 3 gwnums (x, z, AD4) + 10 for ell_mul, + 16 for Lucas chain states */

		for (mult = prime; mult <= B1; mult *= prime)
		{
            if (next_prime > B1)
                if( mult*prime > B1) last_mul = 1;

			if (mult > *B1_done)
			{
				PROCESS_PRIME
			}
		}

/* Check for errors */

		if (gw_test_for_error (&gwdata)) goto error;

/* Check for interrupt.  If one occurs return normalized x OR x,z pair. */

		if (stop_check_proc != NULL && (*stop_check_proc)(0)) {
			*B1_done = prime;

			if (z_array == NULL) {
				StopCheckRoutine = NULL;
				normalize (current_xz.x, current_xz.z);
				if (FAC != NULL) goto bingo;
				reslong = gwtobinary64 (&gwdata, current_xz.x, x_array, word_count_64);
				if (reslong < 0) goto error;
				*x_array_len = reslong;
			}

			else {
				reslong = gwtobinary64 (&gwdata, current_xz.x, x_array, word_count_64);
				if (reslong < 0) goto error;
				*x_array_len = reslong;

				reslong = gwtobinary64 (&gwdata, current_xz.z, z_array, word_count_64);
				if (reslong < 0) goto error;
				*z_array_len = reslong;
			}

			FREE_GWNUMS
			ecm_cleanup ();
			return (ES1_INTERRUPT);
		}
	}
	*B1_done = B1;

/* Normalize the x value OR return the x,z pair */

	if (z_array == NULL) {
		StopCheckRoutine = NULL;
		normalize (current_xz.x, current_xz.z);
		if (FAC != NULL) goto bingo;
		reslong = gwtobinary64 (&gwdata, current_xz.x, x_array, word_count_64);
		if (reslong < 0) goto error;
		*x_array_len = reslong;
	} else {
		reslong = gwtobinary64 (&gwdata, current_xz.x, x_array, word_count_64);
		if (reslong < 0) goto error;
		*x_array_len = reslong;

		reslong = gwtobinary64 (&gwdata, current_xz.z, z_array, word_count_64);
		if (reslong < 0) goto error;
		*z_array_len = reslong;
	}

/* Free memory and return */

	FREE_GWNUMS
	ecm_cleanup ();
	return (ES1_SUCCESS);

/* Print a message if we found a factor! */

bingo:	//printf ("ECM found a factor\n");
	if (!testFactor (FAC)) goto error;
	gianttogw (&gwdata, FAC, x);
	reslong = gwtobinary64 (&gwdata, x, x_array, word_count_64);
	if (reslong < 0) goto error;
	*x_array_len = reslong;
	if (z_array != NULL) {
		z_array[0] = 1;
		*z_array_len = 1;
	}
	FREE_GWNUMS
	ecm_cleanup ();
	return (ES1_FACTOR_FOUND);

/* Return a hardware error occurred code */

error:	FREE_GWNUMS
	ecm_cleanup ();
	return (ES1_HARDWARE_ERROR);

/* Return out-of-memory error */

no_mem:	ecm_cleanup ();
	return (ES1_MEMORY);
}
