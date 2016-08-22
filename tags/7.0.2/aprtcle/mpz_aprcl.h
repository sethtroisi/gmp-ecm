/* Copyright 2011-2015 David Cleaver
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __MPZ_APRCL__
#define __MPZ_APRCL__

#ifndef HAVE_U64_T
#define HAVE_U64_T
typedef long long s64_t;
typedef unsigned long long u64_t;
#endif

#include "jacobi_sum.h"

/*******************************************************/
/*******************************************************/
/* These are the definitions for the APRT-CLE routines */
/*******************************************************/
/*******************************************************/
/* verbose = 0 means to only return the status  */
/*          it will not print out progress info */
/* verbose = 1 means to print out progress info */
/* verbose = 2 means to print out progress/failure/failover info */
#define APRTCLE_VERBOSE0 0
#define APRTCLE_VERBOSE1 1
#define APRTCLE_VERBOSE2 2

#define APRTCLE_ERROR -1
#define APRTCLE_COMPOSITE 0
#define APRTCLE_PRP 1
#define APRTCLE_PRIME 2

/* **********************************************************************************
 * APR-CL (also known as APRT-CLE) is a prime proving algorithm developed by:
 * L. Adleman, C. Pomerance, R. Rumely, H. Cohen, and H. W. Lenstra
 * APRT-CLE = Adleman-Pomerance-Rumely Test Cohen-Lenstra Extended version
 * You can find all the details of this implementation in the Cohen & Lenstra paper:
 *    H. Cohen and A. K. Lenstra, "Implementation of a new primality test",
 *    Math. Comp., 48 (1987) 103--121
 *
 * ----------------------------------------------------------------------------------
 *
 * This C/GMP version is a conversion of Dario Alpern's Java based APRT-CLE code
 * His code was based on Yuji Kida's UBASIC code
 *
 * Based on APRT-CLE Written by Dario Alejandro Alpern (Buenos Aires - Argentina)
 * From: Last updated September 10th, 2011. See http://www.alpertron.com.ar/ECM.HTM
 *
 * On 2012/11/12 Dario Alpern has approved the conversion, from Java to C/GMP, of 
 * his implementation of the APR-CL algorithm, and that it be licensed under the LGPL.
 *
 * ----------------------------------------------------------------------------------
 *
 * With improvements based on Jason Moxham's APRCL v1.15 code, from 2003/01/01
 *
 * On 2013/04/14 Toby Moxham has approved the APR-CL code and data tables, 
 * originally written by his brother Jason Moxham on 2003/01/01, to be released 
 * under the LGPL.
 *
 * ----------------------------------------------------------------------------------
 *
 * v1.0 to v1.1 improvements:
 *
 * [The following fix was recommended by Dana Jacobsen and verified by Jon Grantham]
 *      - Bug fix: Removed unnecessary vl==0 check in mpz_extrastronglucas_prp
 * [The following improvements/fixes were recommended by Laurent Desnogues in 2013/08]
 *      - Speed improvement 1: Removed extraneous NormalizeJS calls in ARPCL
 *      - Speed improvement 2: Removed/consolidated calls to mpz_mod in APRCL
 *        (these improvements make the APRCL code about 1.5-2.2x faster)
 *      - Bug fix: Final test in APRCL routine is now correct
 *
 *
 * *********************************************************************************/

int mpz_aprcl(mpz_t N); /* Just return the status of the input, no progress is printed out */
int mpz_aprtcle(mpz_t N, int verbose);

#endif
