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

//  mp_limb_t mulredc1(mp_limb_t * z, const mp_limb_t x, const mp_limb_t y,
//                 const mp_limb_t m, mp_limb_t inv_m);
//
// arguments:
// r3  : ptr to result z
// r4  : input x
// r5  : input y
// r6  : modulus m
// r7 = -1/m mod 2^64
//
// final carry returned in r3



include(`config.m4')

	TEXT
.align 5 // powerPC 32 byte alignment
	GLOBL GSYM_PREFIX`'mulredc1
	TYPE(GSYM_PREFIX`'mulredc`'1,`function')

GSYM_PREFIX`'mulredc1:
		mulld   r8, r4, r5			// x*y low half T0
		mulhdu  r9, r4, r5			// x*y high half T1
		mulld   r0, r7, r8			// u = t0 * invm
		mulld   r10, r0, r6			// u*m low
		mulhdu  r11, r0, r6			// u*m high
		addc    r8, r8, r10			// x*y + u*m low (= zero)
		adde    r9, r9, r11			// result
		std     r9, 0(r3)			// store in z
		addze   r3, r8				// return carry
		blr
