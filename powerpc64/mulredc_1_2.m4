`dnl ******************************************************************************
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
dnl ******************************************************************************'

ifelse(eval(LENGTH),1,
`//  mp_limb_t mulredc1(mp_limb_t * z, const mp_limb_t x, const mp_limb_t y,
//                 const mp_limb_t m, mp_limb_t inv_m);
//
// arguments:
// r3  : ptr to result z
// r4  : input x
// r5  : input y
// r6  : modulus m',
`// mp_limb_t mulredc'LENGTH`(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y,
//                 const mp_limb_t *m, mp_limb_t inv_m);
//
// arguments:
// r3 = ptr to result z least significant limb
// r4 = ptr to input x least significant limb
// r5 = ptr to input y least significant limb
// r6 = ptr to modulus m least significant limb')
// r7 = -1/m mod 2^64
//
// final carry returned in r3

divert(-1)
dnl  forloop(i, from, to, stmt)

define(`forloop', `pushdef(`$1', `$2')_forloop(`$1', `$2', `$3', `$4')popdef(`$1')')
define(`_forloop',
       `$4`'ifelse($1, `$3', ,
                         `define(`$1', incr($1))_forloop(`$1', `$2', `$3', `$4')')')
divert

`include(`config.m4')'

	TEXT
.align 5 // powerPC 32 byte alignment
	GLOBL GSYM_PREFIX``''mulredc`'LENGTH
	TYPE(GSYM_PREFIX``''mulredc``''LENGTH,``function'')

GSYM_PREFIX``''mulredc`'LENGTH:
ifelse(eval(LENGTH),1,
`		mulld   r8, r4, r5			// x*y low half T0
		mulhdu  r9, r4, r5			// x*y high half T1
		mulld   r0, r7, r8			// u = t0 * invm
		mulld   r10, r0, r6			// u*m low
		mulhdu  r11, r0, r6			// u*m high
		addc    r8, r8, r10			// x*y + u*m low (= zero)
		adde    r9, r9, r11			// result
		std     r9, 0(r3)			// store in z
		addze   r3, r8				// return carry
		blr',
eval(LENGTH),2,
`		ld      r12, 0(r4)          // XI = x[0]
		ld      r0, 0(r5)           // y[0]
		stdu    r13, -8(r1)         // save r13
		mulld   r8, r0, r12         // x[0]*y[0] low half
		stdu    r14, -8(r1)         // save r14
		mulhdu  r9, r0, r12         // x[0]*y[0] high half
		ld      r0, 0(r6)           // m[0]
		mulld   r11, r7, r8         // U = T0*invm mod 2^64
		stdu    r15, -8(r1)         // save r15
		mulld   r13, r0, r11        // T0 = U*m[0] low
		stdu    r16, -8(r1)         // save r16
		li      r16, 0              // set r16 to zero for carry propagation
		mulhdu  r14, r0, r11        // T1 = U*m[0] high
		ld      r0, 8(r5)           // y[1]
		addc    r8, r8, r13         // result zero
		mulld   r8, r0, r12         // x[0]*y[1] low half
		adde    r13, r9, r14        // T0 = initial tmp(0)
		addze   r10, r16            // carry to CY

		mulhdu  r9, r0, r12         // x[0]*y[1] high half
		ld      r0, 8(r6)           // m[1]
		addc    r13, r8, r13        // add low word to T0
		mulld   r8, r0, r11         // U*m[1] low
		adde    r14, r9, r10        // add high word with carry + CY to T1
		// T1:T0 <= 2^128 - 2*2^64 + 1 + 2*2^64 - 3 <= 2^128 - 2, no carry!

		mulhdu  r9, r0, r11         // U*m[1] high
		ldu     r12, 8(r4)          // x[1]
		ld      r0, 0(r5)           // y[0]
		addc    r13, r8, r13        // add T0 and low word
		mulld   r8, r0, r12         // x[1]*y[0] low half
		adde    r14, r9, r14        // add high word with carry to T1
		addze   r15, r16            // put carry in r15 (tmp[len] <= 1)
		mulhdu  r9, r0, r12         // x[1]*y[0] high half
		addc    r13, r8, r13        // T0
		ld      r0, 0(r6)           // m[0]
		mulld   r11, r7, r13        // U = T0*invm mod 2^64
		adde    r14, r9, r14        // T1
		mulld   r8, r0, r11         // U*m[0] low
		addze   r10, r16            // CY
		mulhdu  r9, r0, r11         // T1 = U*m[0] high
		ld      r0, 8(r5)           // y[1]
		addc    r8, r8, r13         // result = 0
		adde    r13, r9, r14        // T0, carry pending

		mulld   r8, r0, r12         // x[1]*y[1] low half
		adde    r14, r15, r10       // T1 = tmp[len] + CY + pending carry
		// since tmp[len] <= 1, T1 <= 3 and carry is zero
		mulhdu  r9, r0, r12         // x[1]*y[1] high half
		ld      r0, 8(r6)           // m[1]
		addc    r13, r8, r13        // add low word to T0
		mulld   r8, r0, r11         // U*m[1] low
		adde    r14, r9, r14        // add high to T1
		addze   r10, r16            // CY
		mulhdu  r9, r0, r11         // U*m[1] high
		addc    r8, r8, r13         // add T0 and low word
		adde    r9, r9, r14         // T1, carry pending
		std     r8, 0(r3)           // copy result to z
		stdu    r9, 8(r3)

		addze   r3, r10             // return tmp(len)
		ld      r16, 0(r1)
		ldu     r15, 8(r1)
		ldu     r14, 8(r1)
		ldu     r13, 8(r1)
		addi    r1, r1, 8
		blr')
