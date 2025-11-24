/* defines and prototypes for ecmstag1.c

Copyright 2025 Paul Zimmermann, Philip McLaughlin.

This file is part of the ECM Library.

The ECM Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The ECM Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the ECM Library; see the file COPYING.LIB.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#ifndef _GW_ECMSTAG1_H
#define _GW_ECMSTAG1_H 1

#ifdef __cplusplus
extern "C" {
#endif

/* Return codes */
#define S1_SUCCESS		0	/* Success, but no factor */
#define S1_FACTOR_FOUND	1	/* Success, factor found */
#define S1_CANNOT_DO_IT	2	/* This k,b,n,c cannot be handled */
#define S1_MEMORY		3	/* Out of memory */
#define S1_INTERRUPT		4	/* Execution interrupted */
#define S1_CANNOT_DO_QUICKLY	5	/* Requires 3-multiply reduction */
#define S1_HARDWARE_ERROR	6	/* An error was detected, most likely a hardware error. */

/* Option codes */
#define S1_DO_SLOW_CASE	0x1	/* Set this if ecmStage1 should do slow 3-multiply reduction cases. */

/* INPUTS:

Input number (3 possibilities):

1) k,b,n,c set and num_being_factored_array = NULL.  k*b^n+c is factored.
2) k,b,n,c zero and num_being_factored_array set.  num_being_factored is
   worked on using generic 3-multiply reduction
3) k,b,n,c set and num_being_factored_array set.  num_being_factored is
   worked on - it must be a factor of k*b^n+c.

A_array, B1 are required

B1_done is optional.  Use it to resume a stage 1 calculation.

x_array, z_array is the starting point.  If z_array is not given, then
z is assumed to be one.

stop_check_proc is optional

options are defined above


   OUTPUTS:

On success:

   if z_array is NULL, then x_array is set to normalized point
   else x_array, z_array is set to the unnormalized point

On factor found:

   x_array is set to the factor found

On interrupt:

   B1_done is set to the last prime below B1 that was processed.
   If z_array != NULL (preferred) then x_array and z_array are set to the
current point.  The (x,z) point is not normalized because it will
be slow for large numbers.  This is unacceptable during system shutdown.
Caller must allocate x and z arrays large enough to hold any k*b^n+c value.
   If z_array == NULL, then a normalized x is returned. Caller must allocate
x array large enough to hold any value less than num_being_factored.

*/

int gw_ecmStage1_u32 (
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
	unsigned long options);

int gw_ecmStage1_u64 (
	double	k,			/* K in K*B^N+C */
	unsigned long b,		/* B in K*B^N+C */
	unsigned long n,		/* N in K*B^N+C */
	signed long c,			/* C in K*B^N+C */
	uint64_t *num_being_factored_array, /* Number to factor */
	unsigned long num_being_factored_array_len,
	uint64_t B1,			/* Stage 1 bound */
	uint64_t *B1_done,		/* Stage 1 that is already done */
	uint64_t *A_array,		/* A - caller derives it from sigma */
	unsigned long A_array_len,
	uint64_t *x_array,		/* X value of point */
	unsigned long *x_array_len,
	uint64_t *z_array,		/* Z value of point */
	unsigned long *z_array_len,
	int	(*stop_check_proc)(int),/* Ptr to proc that returns TRUE */
					/* if user interrupts processing */
	unsigned long options);

#ifdef __cplusplus
}
#endif

#endif /* _GW_ECMSTAG1_H */

