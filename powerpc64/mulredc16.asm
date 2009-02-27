dnl ******************************************************************************
dnl   Copyright 2009 Paul Zimmermann and Alexander Kruppa.
dnl 
dnl   This file is part of the ECM Library.
dnl 
dnl   The ECM Library is free software; you can redistribute it and/or modify
dnl   it under the terms of the GNU Lesser General Public License as published by
dnl   the Free Software Foundation; either version 2.1 of the License, or (at your
dnl   option) any later version.
dnl 
dnl   The ECM Library is distributed in the hope that it will be useful, but
dnl   WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
dnl   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
dnl   License for more details.
dnl 
dnl   You should have received a copy of the GNU Lesser General Public License
dnl   along with the ECM Library; see the file COPYING.LIB.  If not, write to
dnl   the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
dnl   MA 02110-1301, USA.
dnl ******************************************************************************

// mp_limb_t mulredc16(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y,
//                 const mp_limb_t *m, mp_limb_t inv_m);
//
// arguments:
// r3 = ptr to result z least significant limb
// r4 = ptr to input x least significant limb
// r5 = ptr to input y least significant limb
// r6 = ptr to modulus m least significant limb
// r7 = -1/m mod 2^64
//
// final carry returned in r3



include(`config.m4')

	TEXT
.align 5 // powerPC 32 byte alignment
	GLOBL GSYM_PREFIX`'mulredc16
	TYPE(GSYM_PREFIX`'mulredc`'16,`function')

/* Implements multiplication and REDC for two input numbers of 16 words */

// The algorithm:
//   (Notation: a:b:c == a * 2^128 + b * 2^64 + c)
//
// T1:T0 = x[i]*y[0] ;
// u = (T0*invm) % 2^64 ;
// cy:T1 = (m[0]*u + T1:T0) / 2^64 ; /* cy:T1 <= 2*2^64 - 4 (see note 1) */
// for (j = 1; j < len; j++)
//   {
//     cy:T1:T0 = x[i]*y[j] + m[j]*u + cy:T1 ;
//        /* for all j result cy:T1 <= 2*2^64 - 3 (see note 2) */
//     tmp[j-1] = T0;
//   }
// tmp[len-1] = T1 ;
// tmp[len] = cy ; /* cy <= 1 (see note 2) */
// for (i = 1; i < len; i++)
//   {
//     cy:T1:T0 = x[i]*y[0] + tmp[1]:tmp[0] ;
//     u = (T0*invm) % 2^64 ;
//     cy:T1 = (m[0]*u + cy:T1:T0) / 2^64 ; /* cy:T1 <= 3*2^64 - 4 (see note 3) */
//     for (j = 1; j < len; j++)
//       {
//         cy:T1:T0 = x[i]*y[j] + m[j]*u + (tmp[j+1] + cy):T1 ;
//         /* for all j < (len-1), result cy:T1 <= 3*2^64 - 3
//            for j = (len-1), result cy:T1 <= 2*2^64 - 1  (see note 4) */
//         tmp[j-1] = T0;
//       }
//     tmp[len-1] = T1 ;
//     tmp[len] = cy ; /* cy <= 1 for all i (see note 4) */
//   }
// z[0 ... len-1] = tmp[0 ... len-1] ;
// return (tmp[len]) ;
//
// notes:
//
// 1:  m[0]*u + T1:T0 <= 2*(2^64 - 1)^2 <= 2*2^128 - 4*2^64 + 2,
//     so cy:T1 <= 2*2^64 - 4.
// 2:  For j = 1, x[i]*y[j] + m[j]*u + cy:T1 <= 2*(2^64 - 1)^2 + 2*2^64 - 4
//                 <= 2*2^128 - 2*2^64 - 2 = 1:(2^64-3):(2^64-2),
//     so cy:T1 <= 2*2^64 - 3. For j > 1,
//     x[i]*y[j] + m[j]*u + cy:T1 <= 2*2^128 - 2*2^64 - 1 = 1:(2^64-3):(2^64-1),
//     so cy:T1 <= 2*2^64 - 3 = 1:(2^64-3) holds for all j.
// 3:  m[0]*u + cy:T1:T0 <= 2*(2^64 - 1)^2 + 2^128 - 1 = 3*2^128 - 4*2^64 + 1,
//     so cy:T1 <= 3*2^64 - 4 = 2:(2^64-4)
// 4:  For j = 1, x[i]*y[j] + m[j]*u + (tmp[j+1] + cy):T1
//                  <= 2*(2^64 - 1)^2 + (3*2^64 - 4) + (2^64-1)*2^64
//                  <= 3*2^128 - 2*2^64 - 2 = 2:(2^64-3):(2^64-2),
//     so cy:T1 <= 3*2^64 - 3. For j > 1,
//     x[i]*y[j] + m[j]*u + (tmp[j+1] + cy):T1 <= 2:(2^64-3):(2^64-1),
//     so cy:T1 <= 3*2^64 - 3 = 2:(2^64-3) holds for all j < len - 1.
//     For j = len - 1, we know from note 2 that tmp(len) <= 1 for i = 0.
//     Assume this is true for index i-1, Then
//                x[i]*y[len-1] + m[len-1]*u + (tmp[len] + cy):T1
//                  <= 2*(2^64 - 1)^2 + (3*2^64 - 3) + 2^64
//                  <= 2*2^128 - 1 = 1:(2^64-1):(2^64-1),
//     so cy:T1 <= 1:(2^64-1) and tmp[len] <= 1 for all i by induction.
//
// Register vars: T0 = r13, T1 = r14, CY = r10, XI = r12, U = r11
//                YP = r5, MP = r6, TP = r1 (stack ptr)
//

// local variables: tmp[0 ... 16] array, having 16+1 8-byte words
// The tmp array needs 16+1 entries, but tmp[16] is stored in
// r15, so only 16 entries are used in the stack.



GSYM_PREFIX`'mulredc16:

//########################################################################
// i = 0 pass
//########################################################################

// Pass for j = 0. We need to fetch x[i] from memory and compute the new u

	ld      r12, 0(r4)		// XI = x[0]
	ld      r0, 0(r5)		// y[0]
	stdu    r13, -8(r1)		// save r13
	mulld   r8, r0, r12		// x[0]*y[0] low half
	stdu    r14, -8(r1)		// save r14
	mulhdu  r9, r0, r12		// x[0]*y[0] high half
	ld      r0, 0(r6)		// m[0]
	mulld   r11, r7, r8		// U = T0*invm mod 2^64
	stdu    r15, -8(r1)		// save r15
	mulld   r13, r0, r11		// T0 = U*m[0] low
	stdu    r16, -8(r1)		// save r16
	li      r16, 0			// set r16 to zero for carry propagation
	subi    r1, r1, 128		// set tmp stack space
	mulhdu  r14, r0, r11		// T1 = U*m[0] high
	ld      r0, 8(r5)		// y[1]
	addc    r8, r8, r13		//
	adde    r13, r9, r14		// T0 = initial tmp(0)
	addze   r10, r16		// carry to CY
	// CY:T1:T0 <= 2*(2^64-1)^2 <= 2^2*128 - 4*2^64 + 2, hence
	// CY:T1 <= 2*2^64 - 4

// Pass for j = 1

	mulld   r8, r0, r12		// x[i]*y[j] low half
	mulhdu  r9, r0, r12		// x[i]*y[j] high half
	ld      r0, 8(r6)		// m[j]
	addc    r13, r8, r13		// add low word to T0
	adde    r14, r9, r10		// add high word with carry + CY to T1
	// T1:T0 <= 2^128 - 2*2^64 + 1 + 2*2^64 - 3 <= 2^128 - 2, no carry!

	mulld   r8, r0, r11		// U*m[j] low
	mulhdu  r9, r0, r11		// U*m[j] high
	addc    r8, r8, r13		// add T0 and low word
	ld      r0, 16(r5)		// y[j+1]
	adde    r13, r9, r14		// add high word with carry to T1
	addze   r10, r16		// carry to CY
	std     r8, 0(r1)		// store tmp[j-1]
	// CY:T1:T0 <= 2^128 - 2 + 2^128 - 2*2^64 + 1 <=
	//             2 * 2^128 - 2*2^64 - 1 ==> CY:T1 <= 2 * 2^64 - 3

// Pass for j = 2

	mulld   r8, r0, r12		// x[i]*y[j] low half
	mulhdu  r9, r0, r12		// x[i]*y[j] high half
	ld      r0, 16(r6)		// m[j]
	addc    r13, r8, r13		// add low word to T0
	adde    r14, r9, r10		// add high word with carry + CY to T1
	// T1:T0 <= 2^128 - 2*2^64 + 1 + 2*2^64 - 3 <= 2^128 - 2, no carry!

	mulld   r8, r0, r11		// U*m[j] low
	mulhdu  r9, r0, r11		// U*m[j] high
	addc    r8, r8, r13		// add T0 and low word
	ld      r0, 24(r5)		// y[j+1]
	adde    r13, r9, r14		// add high word with carry to T1
	addze   r10, r16		// carry to CY
	std     r8, 8(r1)		// store tmp[j-1]
	// CY:T1:T0 <= 2^128 - 2 + 2^128 - 2*2^64 + 1 <=
	//             2 * 2^128 - 2*2^64 - 1 ==> CY:T1 <= 2 * 2^64 - 3

// Pass for j = 3

	mulld   r8, r0, r12		// x[i]*y[j] low half
	mulhdu  r9, r0, r12		// x[i]*y[j] high half
	ld      r0, 24(r6)		// m[j]
	addc    r13, r8, r13		// add low word to T0
	adde    r14, r9, r10		// add high word with carry + CY to T1
	// T1:T0 <= 2^128 - 2*2^64 + 1 + 2*2^64 - 3 <= 2^128 - 2, no carry!

	mulld   r8, r0, r11		// U*m[j] low
	mulhdu  r9, r0, r11		// U*m[j] high
	addc    r8, r8, r13		// add T0 and low word
	ld      r0, 32(r5)		// y[j+1]
	adde    r13, r9, r14		// add high word with carry to T1
	addze   r10, r16		// carry to CY
	std     r8, 16(r1)		// store tmp[j-1]
	// CY:T1:T0 <= 2^128 - 2 + 2^128 - 2*2^64 + 1 <=
	//             2 * 2^128 - 2*2^64 - 1 ==> CY:T1 <= 2 * 2^64 - 3

// Pass for j = 4

	mulld   r8, r0, r12		// x[i]*y[j] low half
	mulhdu  r9, r0, r12		// x[i]*y[j] high half
	ld      r0, 32(r6)		// m[j]
	addc    r13, r8, r13		// add low word to T0
	adde    r14, r9, r10		// add high word with carry + CY to T1
	// T1:T0 <= 2^128 - 2*2^64 + 1 + 2*2^64 - 3 <= 2^128 - 2, no carry!

	mulld   r8, r0, r11		// U*m[j] low
	mulhdu  r9, r0, r11		// U*m[j] high
	addc    r8, r8, r13		// add T0 and low word
	ld      r0, 40(r5)		// y[j+1]
	adde    r13, r9, r14		// add high word with carry to T1
	addze   r10, r16		// carry to CY
	std     r8, 24(r1)		// store tmp[j-1]
	// CY:T1:T0 <= 2^128 - 2 + 2^128 - 2*2^64 + 1 <=
	//             2 * 2^128 - 2*2^64 - 1 ==> CY:T1 <= 2 * 2^64 - 3

// Pass for j = 5

	mulld   r8, r0, r12		// x[i]*y[j] low half
	mulhdu  r9, r0, r12		// x[i]*y[j] high half
	ld      r0, 40(r6)		// m[j]
	addc    r13, r8, r13		// add low word to T0
	adde    r14, r9, r10		// add high word with carry + CY to T1
	// T1:T0 <= 2^128 - 2*2^64 + 1 + 2*2^64 - 3 <= 2^128 - 2, no carry!

	mulld   r8, r0, r11		// U*m[j] low
	mulhdu  r9, r0, r11		// U*m[j] high
	addc    r8, r8, r13		// add T0 and low word
	ld      r0, 48(r5)		// y[j+1]
	adde    r13, r9, r14		// add high word with carry to T1
	addze   r10, r16		// carry to CY
	std     r8, 32(r1)		// store tmp[j-1]
	// CY:T1:T0 <= 2^128 - 2 + 2^128 - 2*2^64 + 1 <=
	//             2 * 2^128 - 2*2^64 - 1 ==> CY:T1 <= 2 * 2^64 - 3

// Pass for j = 6

	mulld   r8, r0, r12		// x[i]*y[j] low half
	mulhdu  r9, r0, r12		// x[i]*y[j] high half
	ld      r0, 48(r6)		// m[j]
	addc    r13, r8, r13		// add low word to T0
	adde    r14, r9, r10		// add high word with carry + CY to T1
	// T1:T0 <= 2^128 - 2*2^64 + 1 + 2*2^64 - 3 <= 2^128 - 2, no carry!

	mulld   r8, r0, r11		// U*m[j] low
	mulhdu  r9, r0, r11		// U*m[j] high
	addc    r8, r8, r13		// add T0 and low word
	ld      r0, 56(r5)		// y[j+1]
	adde    r13, r9, r14		// add high word with carry to T1
	addze   r10, r16		// carry to CY
	std     r8, 40(r1)		// store tmp[j-1]
	// CY:T1:T0 <= 2^128 - 2 + 2^128 - 2*2^64 + 1 <=
	//             2 * 2^128 - 2*2^64 - 1 ==> CY:T1 <= 2 * 2^64 - 3

// Pass for j = 7

	mulld   r8, r0, r12		// x[i]*y[j] low half
	mulhdu  r9, r0, r12		// x[i]*y[j] high half
	ld      r0, 56(r6)		// m[j]
	addc    r13, r8, r13		// add low word to T0
	adde    r14, r9, r10		// add high word with carry + CY to T1
	// T1:T0 <= 2^128 - 2*2^64 + 1 + 2*2^64 - 3 <= 2^128 - 2, no carry!

	mulld   r8, r0, r11		// U*m[j] low
	mulhdu  r9, r0, r11		// U*m[j] high
	addc    r8, r8, r13		// add T0 and low word
	ld      r0, 64(r5)		// y[j+1]
	adde    r13, r9, r14		// add high word with carry to T1
	addze   r10, r16		// carry to CY
	std     r8, 48(r1)		// store tmp[j-1]
	// CY:T1:T0 <= 2^128 - 2 + 2^128 - 2*2^64 + 1 <=
	//             2 * 2^128 - 2*2^64 - 1 ==> CY:T1 <= 2 * 2^64 - 3

// Pass for j = 8

	mulld   r8, r0, r12		// x[i]*y[j] low half
	mulhdu  r9, r0, r12		// x[i]*y[j] high half
	ld      r0, 64(r6)		// m[j]
	addc    r13, r8, r13		// add low word to T0
	adde    r14, r9, r10		// add high word with carry + CY to T1
	// T1:T0 <= 2^128 - 2*2^64 + 1 + 2*2^64 - 3 <= 2^128 - 2, no carry!

	mulld   r8, r0, r11		// U*m[j] low
	mulhdu  r9, r0, r11		// U*m[j] high
	addc    r8, r8, r13		// add T0 and low word
	ld      r0, 72(r5)		// y[j+1]
	adde    r13, r9, r14		// add high word with carry to T1
	addze   r10, r16		// carry to CY
	std     r8, 56(r1)		// store tmp[j-1]
	// CY:T1:T0 <= 2^128 - 2 + 2^128 - 2*2^64 + 1 <=
	//             2 * 2^128 - 2*2^64 - 1 ==> CY:T1 <= 2 * 2^64 - 3

// Pass for j = 9

	mulld   r8, r0, r12		// x[i]*y[j] low half
	mulhdu  r9, r0, r12		// x[i]*y[j] high half
	ld      r0, 72(r6)		// m[j]
	addc    r13, r8, r13		// add low word to T0
	adde    r14, r9, r10		// add high word with carry + CY to T1
	// T1:T0 <= 2^128 - 2*2^64 + 1 + 2*2^64 - 3 <= 2^128 - 2, no carry!

	mulld   r8, r0, r11		// U*m[j] low
	mulhdu  r9, r0, r11		// U*m[j] high
	addc    r8, r8, r13		// add T0 and low word
	ld      r0, 80(r5)		// y[j+1]
	adde    r13, r9, r14		// add high word with carry to T1
	addze   r10, r16		// carry to CY
	std     r8, 64(r1)		// store tmp[j-1]
	// CY:T1:T0 <= 2^128 - 2 + 2^128 - 2*2^64 + 1 <=
	//             2 * 2^128 - 2*2^64 - 1 ==> CY:T1 <= 2 * 2^64 - 3

// Pass for j = 10

	mulld   r8, r0, r12		// x[i]*y[j] low half
	mulhdu  r9, r0, r12		// x[i]*y[j] high half
	ld      r0, 80(r6)		// m[j]
	addc    r13, r8, r13		// add low word to T0
	adde    r14, r9, r10		// add high word with carry + CY to T1
	// T1:T0 <= 2^128 - 2*2^64 + 1 + 2*2^64 - 3 <= 2^128 - 2, no carry!

	mulld   r8, r0, r11		// U*m[j] low
	mulhdu  r9, r0, r11		// U*m[j] high
	addc    r8, r8, r13		// add T0 and low word
	ld      r0, 88(r5)		// y[j+1]
	adde    r13, r9, r14		// add high word with carry to T1
	addze   r10, r16		// carry to CY
	std     r8, 72(r1)		// store tmp[j-1]
	// CY:T1:T0 <= 2^128 - 2 + 2^128 - 2*2^64 + 1 <=
	//             2 * 2^128 - 2*2^64 - 1 ==> CY:T1 <= 2 * 2^64 - 3

// Pass for j = 11

	mulld   r8, r0, r12		// x[i]*y[j] low half
	mulhdu  r9, r0, r12		// x[i]*y[j] high half
	ld      r0, 88(r6)		// m[j]
	addc    r13, r8, r13		// add low word to T0
	adde    r14, r9, r10		// add high word with carry + CY to T1
	// T1:T0 <= 2^128 - 2*2^64 + 1 + 2*2^64 - 3 <= 2^128 - 2, no carry!

	mulld   r8, r0, r11		// U*m[j] low
	mulhdu  r9, r0, r11		// U*m[j] high
	addc    r8, r8, r13		// add T0 and low word
	ld      r0, 96(r5)		// y[j+1]
	adde    r13, r9, r14		// add high word with carry to T1
	addze   r10, r16		// carry to CY
	std     r8, 80(r1)		// store tmp[j-1]
	// CY:T1:T0 <= 2^128 - 2 + 2^128 - 2*2^64 + 1 <=
	//             2 * 2^128 - 2*2^64 - 1 ==> CY:T1 <= 2 * 2^64 - 3

// Pass for j = 12

	mulld   r8, r0, r12		// x[i]*y[j] low half
	mulhdu  r9, r0, r12		// x[i]*y[j] high half
	ld      r0, 96(r6)		// m[j]
	addc    r13, r8, r13		// add low word to T0
	adde    r14, r9, r10		// add high word with carry + CY to T1
	// T1:T0 <= 2^128 - 2*2^64 + 1 + 2*2^64 - 3 <= 2^128 - 2, no carry!

	mulld   r8, r0, r11		// U*m[j] low
	mulhdu  r9, r0, r11		// U*m[j] high
	addc    r8, r8, r13		// add T0 and low word
	ld      r0, 104(r5)		// y[j+1]
	adde    r13, r9, r14		// add high word with carry to T1
	addze   r10, r16		// carry to CY
	std     r8, 88(r1)		// store tmp[j-1]
	// CY:T1:T0 <= 2^128 - 2 + 2^128 - 2*2^64 + 1 <=
	//             2 * 2^128 - 2*2^64 - 1 ==> CY:T1 <= 2 * 2^64 - 3

// Pass for j = 13

	mulld   r8, r0, r12		// x[i]*y[j] low half
	mulhdu  r9, r0, r12		// x[i]*y[j] high half
	ld      r0, 104(r6)		// m[j]
	addc    r13, r8, r13		// add low word to T0
	adde    r14, r9, r10		// add high word with carry + CY to T1
	// T1:T0 <= 2^128 - 2*2^64 + 1 + 2*2^64 - 3 <= 2^128 - 2, no carry!

	mulld   r8, r0, r11		// U*m[j] low
	mulhdu  r9, r0, r11		// U*m[j] high
	addc    r8, r8, r13		// add T0 and low word
	ld      r0, 112(r5)		// y[j+1]
	adde    r13, r9, r14		// add high word with carry to T1
	addze   r10, r16		// carry to CY
	std     r8, 96(r1)		// store tmp[j-1]
	// CY:T1:T0 <= 2^128 - 2 + 2^128 - 2*2^64 + 1 <=
	//             2 * 2^128 - 2*2^64 - 1 ==> CY:T1 <= 2 * 2^64 - 3

// Pass for j = 14

	mulld   r8, r0, r12		// x[i]*y[j] low half
	mulhdu  r9, r0, r12		// x[i]*y[j] high half
	ld      r0, 112(r6)		// m[j]
	addc    r13, r8, r13		// add low word to T0
	adde    r14, r9, r10		// add high word with carry + CY to T1
	// T1:T0 <= 2^128 - 2*2^64 + 1 + 2*2^64 - 3 <= 2^128 - 2, no carry!

	mulld   r8, r0, r11		// U*m[j] low
	mulhdu  r9, r0, r11		// U*m[j] high
	addc    r8, r8, r13		// add T0 and low word
	ld      r0, 120(r5)		// y[j+1]
	adde    r13, r9, r14		// add high word with carry to T1
	addze   r10, r16		// carry to CY
	std     r8, 104(r1)		// store tmp[j-1]
	// CY:T1:T0 <= 2^128 - 2 + 2^128 - 2*2^64 + 1 <=
	//             2 * 2^128 - 2*2^64 - 1 ==> CY:T1 <= 2 * 2^64 - 3

// Pass for j = 15. Don't fetch new data from y[j+1].

	mulld   r8, r0, r12		// x[i]*y[j] low half
	mulhdu  r9, r0, r12		// x[i]*y[j] high half
	ld      r0, 120(r6)		// m[j]
	addc    r13, r8, r13		// add low word to T0
	adde    r14, r9, r10		// add high word with carry + CY to T1
	// T1:T0 <= 2^128 - 2*2^64 + 1 + 2*2^64 - 3 <= 2^128 - 2, no carry!

	mulld   r8, r0, r11		// U*m[j] low
	mulhdu  r9, r0, r11		// U*m[j] high
	addc    r8, r8, r13		// add T0 and low word
	adde    r13, r9, r14		// add high word with carry to T1
	std     r8, 112(r1)		// store tmp[len-2]
	addze   r15, r16		// put carry in r15 (tmp[len] <= 1)
	std     r13, 120(r1)		// store tmp[len-1]


//########################################################################
// i > 0 passes
//########################################################################


	li      r9, 15			// outer loop count
	mtctr   r9

1:

// Pass for j = 0. We need to fetch x[i], tmp[i] and tmp[i+1] from memory
// and compute the new u

	ldu     r12, 8(r4)		// x[i]
	ld      r0, 0(r5)		// y[0]
	ld      r13, 0(r1)		// tmp[0]
	mulld   r8, r0, r12		// x[i]*y[0] low half
	ld      r14, 8(r1)		// tmp[1]
	mulhdu  r9, r0, r12		// x[i]*y[0] high half
	addc    r13, r8, r13		// T0
	ld      r0, 0(r6)		// m[0]
	mulld   r11, r7, r13		// U = T0*invm mod 2^64
	adde    r14, r9, r14		// T1
	mulld   r8, r0, r11		// U*m[0] low
	addze   r10, r16		// CY
	mulhdu  r9, r0, r11		// U*m[0] high
	ld      r0, 8(r5)		// y[1]
	addc    r8, r8, r13		// result = 0
	adde    r13, r9, r14		// T0, carry pending
	// cy:T1:T0 <= 2*(2^64 - 1)^2 + 2^128 - 1 = 3*2^128 - 4*2^64 + 1,
	// so cy:T1 <= 3*2^64 - 4

// Pass for j = 1

	ld      r14, 16(r1)		// tmp[j+1]
	mulld   r8, r0, r12		// x[i]*y[j] low half
	adde    r14, r14, r10		// tmp[j+1] + CY + pending carry
	addze   r10, r16		// carry to CY
	mulhdu  r9, r0, r12		// x[i]*y[j] high half
	ld      r0, 8(r6)		// m[j]
	addc    r13, r8, r13		// add low word to T0
	mulld   r8, r0, r11		// U*m[j] low
	adde    r14, r9, r14		// add high to T1
	addze   r10, r10		// add carry to CY
	mulhdu  r9, r0, r11		// U*m[j] high
	addc    r8, r8, r13		// add T0 and low word
	ld      r0, 16(r5)		// y[j+1]
	adde    r13, r9, r14		// T1, carry pending
	std     r8, 0(r1)		// store tmp[j-1]
	// CY:T1:T0 <= 2*(2^64 - 1)^2 + (3*2^64 - 3) + (2^64-1)*2^64
	//          <= 3*2^128 - 2*2^64 - 1 ==> CY:T1 <= 3*2^64 - 3

// Pass for j = 2

	ld      r14, 24(r1)		// tmp[j+1]
	mulld   r8, r0, r12		// x[i]*y[j] low half
	adde    r14, r14, r10		// tmp[j+1] + CY + pending carry
	addze   r10, r16		// carry to CY
	mulhdu  r9, r0, r12		// x[i]*y[j] high half
	ld      r0, 16(r6)		// m[j]
	addc    r13, r8, r13		// add low word to T0
	mulld   r8, r0, r11		// U*m[j] low
	adde    r14, r9, r14		// add high to T1
	addze   r10, r10		// add carry to CY
	mulhdu  r9, r0, r11		// U*m[j] high
	addc    r8, r8, r13		// add T0 and low word
	ld      r0, 24(r5)		// y[j+1]
	adde    r13, r9, r14		// T1, carry pending
	std     r8, 8(r1)		// store tmp[j-1]
	// CY:T1:T0 <= 2*(2^64 - 1)^2 + (3*2^64 - 3) + (2^64-1)*2^64
	//          <= 3*2^128 - 2*2^64 - 1 ==> CY:T1 <= 3*2^64 - 3

// Pass for j = 3

	ld      r14, 32(r1)		// tmp[j+1]
	mulld   r8, r0, r12		// x[i]*y[j] low half
	adde    r14, r14, r10		// tmp[j+1] + CY + pending carry
	addze   r10, r16		// carry to CY
	mulhdu  r9, r0, r12		// x[i]*y[j] high half
	ld      r0, 24(r6)		// m[j]
	addc    r13, r8, r13		// add low word to T0
	mulld   r8, r0, r11		// U*m[j] low
	adde    r14, r9, r14		// add high to T1
	addze   r10, r10		// add carry to CY
	mulhdu  r9, r0, r11		// U*m[j] high
	addc    r8, r8, r13		// add T0 and low word
	ld      r0, 32(r5)		// y[j+1]
	adde    r13, r9, r14		// T1, carry pending
	std     r8, 16(r1)		// store tmp[j-1]
	// CY:T1:T0 <= 2*(2^64 - 1)^2 + (3*2^64 - 3) + (2^64-1)*2^64
	//          <= 3*2^128 - 2*2^64 - 1 ==> CY:T1 <= 3*2^64 - 3

// Pass for j = 4

	ld      r14, 40(r1)		// tmp[j+1]
	mulld   r8, r0, r12		// x[i]*y[j] low half
	adde    r14, r14, r10		// tmp[j+1] + CY + pending carry
	addze   r10, r16		// carry to CY
	mulhdu  r9, r0, r12		// x[i]*y[j] high half
	ld      r0, 32(r6)		// m[j]
	addc    r13, r8, r13		// add low word to T0
	mulld   r8, r0, r11		// U*m[j] low
	adde    r14, r9, r14		// add high to T1
	addze   r10, r10		// add carry to CY
	mulhdu  r9, r0, r11		// U*m[j] high
	addc    r8, r8, r13		// add T0 and low word
	ld      r0, 40(r5)		// y[j+1]
	adde    r13, r9, r14		// T1, carry pending
	std     r8, 24(r1)		// store tmp[j-1]
	// CY:T1:T0 <= 2*(2^64 - 1)^2 + (3*2^64 - 3) + (2^64-1)*2^64
	//          <= 3*2^128 - 2*2^64 - 1 ==> CY:T1 <= 3*2^64 - 3

// Pass for j = 5

	ld      r14, 48(r1)		// tmp[j+1]
	mulld   r8, r0, r12		// x[i]*y[j] low half
	adde    r14, r14, r10		// tmp[j+1] + CY + pending carry
	addze   r10, r16		// carry to CY
	mulhdu  r9, r0, r12		// x[i]*y[j] high half
	ld      r0, 40(r6)		// m[j]
	addc    r13, r8, r13		// add low word to T0
	mulld   r8, r0, r11		// U*m[j] low
	adde    r14, r9, r14		// add high to T1
	addze   r10, r10		// add carry to CY
	mulhdu  r9, r0, r11		// U*m[j] high
	addc    r8, r8, r13		// add T0 and low word
	ld      r0, 48(r5)		// y[j+1]
	adde    r13, r9, r14		// T1, carry pending
	std     r8, 32(r1)		// store tmp[j-1]
	// CY:T1:T0 <= 2*(2^64 - 1)^2 + (3*2^64 - 3) + (2^64-1)*2^64
	//          <= 3*2^128 - 2*2^64 - 1 ==> CY:T1 <= 3*2^64 - 3

// Pass for j = 6

	ld      r14, 56(r1)		// tmp[j+1]
	mulld   r8, r0, r12		// x[i]*y[j] low half
	adde    r14, r14, r10		// tmp[j+1] + CY + pending carry
	addze   r10, r16		// carry to CY
	mulhdu  r9, r0, r12		// x[i]*y[j] high half
	ld      r0, 48(r6)		// m[j]
	addc    r13, r8, r13		// add low word to T0
	mulld   r8, r0, r11		// U*m[j] low
	adde    r14, r9, r14		// add high to T1
	addze   r10, r10		// add carry to CY
	mulhdu  r9, r0, r11		// U*m[j] high
	addc    r8, r8, r13		// add T0 and low word
	ld      r0, 56(r5)		// y[j+1]
	adde    r13, r9, r14		// T1, carry pending
	std     r8, 40(r1)		// store tmp[j-1]
	// CY:T1:T0 <= 2*(2^64 - 1)^2 + (3*2^64 - 3) + (2^64-1)*2^64
	//          <= 3*2^128 - 2*2^64 - 1 ==> CY:T1 <= 3*2^64 - 3

// Pass for j = 7

	ld      r14, 64(r1)		// tmp[j+1]
	mulld   r8, r0, r12		// x[i]*y[j] low half
	adde    r14, r14, r10		// tmp[j+1] + CY + pending carry
	addze   r10, r16		// carry to CY
	mulhdu  r9, r0, r12		// x[i]*y[j] high half
	ld      r0, 56(r6)		// m[j]
	addc    r13, r8, r13		// add low word to T0
	mulld   r8, r0, r11		// U*m[j] low
	adde    r14, r9, r14		// add high to T1
	addze   r10, r10		// add carry to CY
	mulhdu  r9, r0, r11		// U*m[j] high
	addc    r8, r8, r13		// add T0 and low word
	ld      r0, 64(r5)		// y[j+1]
	adde    r13, r9, r14		// T1, carry pending
	std     r8, 48(r1)		// store tmp[j-1]
	// CY:T1:T0 <= 2*(2^64 - 1)^2 + (3*2^64 - 3) + (2^64-1)*2^64
	//          <= 3*2^128 - 2*2^64 - 1 ==> CY:T1 <= 3*2^64 - 3

// Pass for j = 8

	ld      r14, 72(r1)		// tmp[j+1]
	mulld   r8, r0, r12		// x[i]*y[j] low half
	adde    r14, r14, r10		// tmp[j+1] + CY + pending carry
	addze   r10, r16		// carry to CY
	mulhdu  r9, r0, r12		// x[i]*y[j] high half
	ld      r0, 64(r6)		// m[j]
	addc    r13, r8, r13		// add low word to T0
	mulld   r8, r0, r11		// U*m[j] low
	adde    r14, r9, r14		// add high to T1
	addze   r10, r10		// add carry to CY
	mulhdu  r9, r0, r11		// U*m[j] high
	addc    r8, r8, r13		// add T0 and low word
	ld      r0, 72(r5)		// y[j+1]
	adde    r13, r9, r14		// T1, carry pending
	std     r8, 56(r1)		// store tmp[j-1]
	// CY:T1:T0 <= 2*(2^64 - 1)^2 + (3*2^64 - 3) + (2^64-1)*2^64
	//          <= 3*2^128 - 2*2^64 - 1 ==> CY:T1 <= 3*2^64 - 3

// Pass for j = 9

	ld      r14, 80(r1)		// tmp[j+1]
	mulld   r8, r0, r12		// x[i]*y[j] low half
	adde    r14, r14, r10		// tmp[j+1] + CY + pending carry
	addze   r10, r16		// carry to CY
	mulhdu  r9, r0, r12		// x[i]*y[j] high half
	ld      r0, 72(r6)		// m[j]
	addc    r13, r8, r13		// add low word to T0
	mulld   r8, r0, r11		// U*m[j] low
	adde    r14, r9, r14		// add high to T1
	addze   r10, r10		// add carry to CY
	mulhdu  r9, r0, r11		// U*m[j] high
	addc    r8, r8, r13		// add T0 and low word
	ld      r0, 80(r5)		// y[j+1]
	adde    r13, r9, r14		// T1, carry pending
	std     r8, 64(r1)		// store tmp[j-1]
	// CY:T1:T0 <= 2*(2^64 - 1)^2 + (3*2^64 - 3) + (2^64-1)*2^64
	//          <= 3*2^128 - 2*2^64 - 1 ==> CY:T1 <= 3*2^64 - 3

// Pass for j = 10

	ld      r14, 88(r1)		// tmp[j+1]
	mulld   r8, r0, r12		// x[i]*y[j] low half
	adde    r14, r14, r10		// tmp[j+1] + CY + pending carry
	addze   r10, r16		// carry to CY
	mulhdu  r9, r0, r12		// x[i]*y[j] high half
	ld      r0, 80(r6)		// m[j]
	addc    r13, r8, r13		// add low word to T0
	mulld   r8, r0, r11		// U*m[j] low
	adde    r14, r9, r14		// add high to T1
	addze   r10, r10		// add carry to CY
	mulhdu  r9, r0, r11		// U*m[j] high
	addc    r8, r8, r13		// add T0 and low word
	ld      r0, 88(r5)		// y[j+1]
	adde    r13, r9, r14		// T1, carry pending
	std     r8, 72(r1)		// store tmp[j-1]
	// CY:T1:T0 <= 2*(2^64 - 1)^2 + (3*2^64 - 3) + (2^64-1)*2^64
	//          <= 3*2^128 - 2*2^64 - 1 ==> CY:T1 <= 3*2^64 - 3

// Pass for j = 11

	ld      r14, 96(r1)		// tmp[j+1]
	mulld   r8, r0, r12		// x[i]*y[j] low half
	adde    r14, r14, r10		// tmp[j+1] + CY + pending carry
	addze   r10, r16		// carry to CY
	mulhdu  r9, r0, r12		// x[i]*y[j] high half
	ld      r0, 88(r6)		// m[j]
	addc    r13, r8, r13		// add low word to T0
	mulld   r8, r0, r11		// U*m[j] low
	adde    r14, r9, r14		// add high to T1
	addze   r10, r10		// add carry to CY
	mulhdu  r9, r0, r11		// U*m[j] high
	addc    r8, r8, r13		// add T0 and low word
	ld      r0, 96(r5)		// y[j+1]
	adde    r13, r9, r14		// T1, carry pending
	std     r8, 80(r1)		// store tmp[j-1]
	// CY:T1:T0 <= 2*(2^64 - 1)^2 + (3*2^64 - 3) + (2^64-1)*2^64
	//          <= 3*2^128 - 2*2^64 - 1 ==> CY:T1 <= 3*2^64 - 3

// Pass for j = 12

	ld      r14, 104(r1)		// tmp[j+1]
	mulld   r8, r0, r12		// x[i]*y[j] low half
	adde    r14, r14, r10		// tmp[j+1] + CY + pending carry
	addze   r10, r16		// carry to CY
	mulhdu  r9, r0, r12		// x[i]*y[j] high half
	ld      r0, 96(r6)		// m[j]
	addc    r13, r8, r13		// add low word to T0
	mulld   r8, r0, r11		// U*m[j] low
	adde    r14, r9, r14		// add high to T1
	addze   r10, r10		// add carry to CY
	mulhdu  r9, r0, r11		// U*m[j] high
	addc    r8, r8, r13		// add T0 and low word
	ld      r0, 104(r5)		// y[j+1]
	adde    r13, r9, r14		// T1, carry pending
	std     r8, 88(r1)		// store tmp[j-1]
	// CY:T1:T0 <= 2*(2^64 - 1)^2 + (3*2^64 - 3) + (2^64-1)*2^64
	//          <= 3*2^128 - 2*2^64 - 1 ==> CY:T1 <= 3*2^64 - 3

// Pass for j = 13

	ld      r14, 112(r1)		// tmp[j+1]
	mulld   r8, r0, r12		// x[i]*y[j] low half
	adde    r14, r14, r10		// tmp[j+1] + CY + pending carry
	addze   r10, r16		// carry to CY
	mulhdu  r9, r0, r12		// x[i]*y[j] high half
	ld      r0, 104(r6)		// m[j]
	addc    r13, r8, r13		// add low word to T0
	mulld   r8, r0, r11		// U*m[j] low
	adde    r14, r9, r14		// add high to T1
	addze   r10, r10		// add carry to CY
	mulhdu  r9, r0, r11		// U*m[j] high
	addc    r8, r8, r13		// add T0 and low word
	ld      r0, 112(r5)		// y[j+1]
	adde    r13, r9, r14		// T1, carry pending
	std     r8, 96(r1)		// store tmp[j-1]
	// CY:T1:T0 <= 2*(2^64 - 1)^2 + (3*2^64 - 3) + (2^64-1)*2^64
	//          <= 3*2^128 - 2*2^64 - 1 ==> CY:T1 <= 3*2^64 - 3

// Pass for j = 14

	ld      r14, 120(r1)		// tmp[j+1]
	mulld   r8, r0, r12		// x[i]*y[j] low half
	adde    r14, r14, r10		// tmp[j+1] + CY + pending carry
	addze   r10, r16		// carry to CY
	mulhdu  r9, r0, r12		// x[i]*y[j] high half
	ld      r0, 112(r6)		// m[j]
	addc    r13, r8, r13		// add low word to T0
	mulld   r8, r0, r11		// U*m[j] low
	adde    r14, r9, r14		// add high to T1
	addze   r10, r10		// add carry to CY
	mulhdu  r9, r0, r11		// U*m[j] high
	addc    r8, r8, r13		// add T0 and low word
	ld      r0, 120(r5)		// y[j+1]
	adde    r13, r9, r14		// T1, carry pending
	std     r8, 104(r1)		// store tmp[j-1]
	// CY:T1:T0 <= 2*(2^64 - 1)^2 + (3*2^64 - 3) + (2^64-1)*2^64
	//          <= 3*2^128 - 2*2^64 - 1 ==> CY:T1 <= 3*2^64 - 3

// Pass for j = 15. Don't fetch new data from y[j+1].

	mulld   r8, r0, r12		// x[i]*y[j] low half
	adde    r14, r15, r10		// T1 = tmp[len] + CY + pending carry
	// since tmp[len] <= 1, T1 <= 3 and carry is zero
	mulhdu  r9, r0, r12		// x[i]*y[j] high half
	ld      r0, 120(r6)		// m[j]
	addc    r13, r8, r13		// add low word to T0
	mulld   r8, r0, r11		// U*m[j] low
	adde    r14, r9, r14		// add high to T1
	addze   r10, r16		// CY
	mulhdu  r9, r0, r11		// U*m[j] high
	addc    r8, r8, r13		// add T0 and low word
	adde    r13, r9, r14		// T1, carry pending
	std     r8, 112(r1)		// store tmp[len-2]
	addze   r15, r10		// store tmp[len] <= 1
	std     r13, 120(r1)		// store tmp[len-1]
	// CY:T1:T0 <= 2*(2^64 - 1)^2 + (3*2^64 - 3) + 2^64
	//          <= 2*2^128 - 1 ==> CY:T1 <= 2*2^64 - 1 = 1:(2^64-1)

	bdnz 1b

// Copy result from tmp memory to z

	ld      r8, 0(r1)
	ldu     r9, 8(r1)
	std     r8, 0(r3)
	stdu    r9, 8(r3)
	ldu     r8, 8(r1)
	ldu     r9, 8(r1)
	stdu    r8, 8(r3)
	stdu    r9, 8(r3)
	ldu     r8, 8(r1)
	ldu     r9, 8(r1)
	stdu    r8, 8(r3)
	stdu    r9, 8(r3)
	ldu     r8, 8(r1)
	ldu     r9, 8(r1)
	stdu    r8, 8(r3)
	stdu    r9, 8(r3)
	ldu     r8, 8(r1)
	ldu     r9, 8(r1)
	stdu    r8, 8(r3)
	stdu    r9, 8(r3)
	ldu     r8, 8(r1)
	ldu     r9, 8(r1)
	stdu    r8, 8(r3)
	stdu    r9, 8(r3)
	ldu     r8, 8(r1)
	ldu     r9, 8(r1)
	stdu    r8, 8(r3)
	stdu    r9, 8(r3)
	ldu     r8, 8(r1)
	ldu     r9, 8(r1)
	stdu    r8, 8(r3)
	stdu    r9, 8(r3)

	mr      r3, r15         // return tmp(len)
	ldu     r16, 8(r1)
	ldu     r15, 8(r1)
	ldu     r14, 8(r1)
	ldu     r13, 8(r1)
	addi    r1, r1, 8
	blr
