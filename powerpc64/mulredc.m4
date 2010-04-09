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

dnl Use `C' to remove comments in .asm -> .s conversion.
dnl Copied from GMP 4.2.
`define(C, `
dnl')'

C mp_limb_t mulredc`'LENGTH`'(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y,
C                 const mp_limb_t *m, mp_limb_t inv_m);
C
C arguments:
C r3 = ptr to result z least significant limb
C r4 = ptr to input x least significant limb
C r5 = ptr to input y least significant limb
C r6 = ptr to modulus m least significant limb
C r7 = -1/m mod 2^64
C
C final carry returned in r3

divert(-1)
dnl  forloop(i, from, to, stmt)

define(`forloop', `pushdef(`$1', `$2')_forloop(`$1', `$2', `$3', `$4')popdef(`$1')')
define(`_forloop',
       `$4`'ifelse($1, `$3', ,
                         `define(`$1', incr($1))_forloop(`$1', `$2', `$3', `$4')')')
divert

`include(`config.m4')'



	TEXT
.align 5 C powerPC 32 byte alignment
	GLOBL GSYM_PREFIX``''mulredc`'LENGTH
	TYPE(GSYM_PREFIX``''mulredc``''LENGTH,``function'')

C Implements multiplication and REDC for two input numbers of LENGTH words

C The algorithm:
C   (Notation: a:b:c == a * 2^128 + b * 2^64 + c)
C
C T1:T0 = x[i]*y[0] ;
C u = (T0*invm) % 2^64 ;
C cy:T1 = (m[0]*u + T1:T0) / 2^64 ; /* cy:T1 <= 2*2^64 - 4 (see note 1) */
C for (j = 1; j < len; j++)
C   {
C     cy:T1:T0 = x[i]*y[j] + m[j]*u + cy:T1 ;
C        /* for all j result cy:T1 <= 2*2^64 - 3 (see note 2) */
C     tmp[j-1] = T0;
C   }
C tmp[len-1] = T1 ;
C tmp[len] = cy ; /* cy <= 1 (see note 2) */
C for (i = 1; i < len; i++)
C   {
C     cy:T1:T0 = x[i]*y[0] + tmp[1]:tmp[0] ;
C     u = (T0*invm) % 2^64 ;
C     cy:T1 = (m[0]*u + cy:T1:T0) / 2^64 ; /* cy:T1 <= 3*2^64 - 4 (see note 3) */
C     for (j = 1; j < len; j++)
C       {
C         cy:T1:T0 = x[i]*y[j] + m[j]*u + (tmp[j+1] + cy):T1 ;
C         /* for all j < (len-1), result cy:T1 <= 3*2^64 - 3
C            for j = (len-1), result cy:T1 <= 2*2^64 - 1  (see note 4) */
C         tmp[j-1] = T0;
C       }
C     tmp[len-1] = T1 ;
C     tmp[len] = cy ; /* cy <= 1 for all i (see note 4) */
C   }
C z[0 ... len-1] = tmp[0 ... len-1] ;
C return (tmp[len]) ;
C
C notes:
C
C 1:  m[0]*u + T1:T0 <= 2*(2^64 - 1)^2 <= 2*2^128 - 4*2^64 + 2,
C     so cy:T1 <= 2*2^64 - 4.
C 2:  For j = 1, x[i]*y[j] + m[j]*u + cy:T1 <= 2*(2^64 - 1)^2 + 2*2^64 - 4
C                 <= 2*2^128 - 2*2^64 - 2 = 1:(2^64-3):(2^64-2),
C     so cy:T1 <= 2*2^64 - 3. For j > 1,
C     x[i]*y[j] + m[j]*u + cy:T1 <= 2*2^128 - 2*2^64 - 1 = 1:(2^64-3):(2^64-1),
C     so cy:T1 <= 2*2^64 - 3 = 1:(2^64-3) holds for all j.
C 3:  m[0]*u + cy:T1:T0 <= 2*(2^64 - 1)^2 + 2^128 - 1 = 3*2^128 - 4*2^64 + 1,
C     so cy:T1 <= 3*2^64 - 4 = 2:(2^64-4)
C 4:  For j = 1, x[i]*y[j] + m[j]*u + (tmp[j+1] + cy):T1
C                  <= 2*(2^64 - 1)^2 + (3*2^64 - 4) + (2^64-1)*2^64
C                  <= 3*2^128 - 2*2^64 - 2 = 2:(2^64-3):(2^64-2),
C     so cy:T1 <= 3*2^64 - 3. For j > 1,
C     x[i]*y[j] + m[j]*u + (tmp[j+1] + cy):T1 <= 2:(2^64-3):(2^64-1),
C     so cy:T1 <= 3*2^64 - 3 = 2:(2^64-3) holds for all j < len - 1.
C     For j = len - 1, we know from note 2 that tmp(len) <= 1 for i = 0.
C     Assume this is true for index i-1, Then
C                x[i]*y[len-1] + m[len-1]*u + (tmp[len] + cy):T1
C                  <= 2*(2^64 - 1)^2 + (3*2^64 - 3) + 2^64
C                  <= 2*2^128 - 1 = 1:(2^64-1):(2^64-1),
C     so cy:T1 <= 1:(2^64-1) and tmp[len] <= 1 for all i by induction.
C
C Register vars: T0 = r13, T1 = r14, CY = r10, XI = r12, U = r11
C                YP = r5, MP = r6, TP = r1 (stack ptr)
C

C local variables: tmp[0 ... LENGTH] array, having LENGTH+1 8-byte words
C The tmp array needs LENGTH+1 entries, but tmp[LENGTH] is stored in
C r15, so only LENGTH entries are used in the stack.



GSYM_PREFIX``''mulredc`'LENGTH:
define(`S', `eval(8 * LENGTH)')dnl

C ########################################################################
C # i = 0 pass
C #########################################################################

C Pass for j = 0. We need to fetch x[i] from memory and compute the new u

	ld      r12, 0(r4)		C XI = x[0]
	ld      r0, 0(r5)		C y[0]
	stdu    r13, -8(r1)		C save r13
	mulld   r8, r0, r12		C x[0]*y[0] low half
	stdu    r14, -8(r1)		C save r14
	mulhdu  r9, r0, r12		C x[0]*y[0] high half
	ld      r0, 0(r6)		C m[0]
	mulld   r11, r7, r8		C U = T0*invm mod 2^64
	stdu    r15, -8(r1)		C save r15
	mulld   r13, r0, r11		C T0 = U*m[0] low
	stdu    r16, -8(r1)		C save r16
	li      r16, 0			C set r16 to zero for carry propagation
	subi    r1, r1, S		C set tmp stack space
	mulhdu  r14, r0, r11		C T1 = U*m[0] high
	ld      r0, 8(r5)		C y[1]
	addc    r8, r8, r13		C
	adde    r13, r9, r14		C T0 = initial tmp(0)
	addze   r10, r16		C carry to CY
	C CY:T1:T0 <= 2*(2^64-1)^2 <= 2^2*128 - 4*2^64 + 2, hence
	C CY:T1 <= 2*2^64 - 4

forloop(`UNROLL', 1, eval(LENGTH - 2), `dnl
define(`J', `eval(8 * UNROLL)')dnl
define(`J8', `eval(J + 8)')dnl
define(`JM8', `eval(J - 8)')dnl
C Pass for j = UNROLL

	mulld   r8, r0, r12		C x[i]*y[j] low half
	mulhdu  r9, r0, r12		C x[i]*y[j] high half
	ld      r0, J`'(r6)		C m[j]
	addc    r13, r8, r13		C add low word to T0
	adde    r14, r9, r10		C add high word with carry + CY to T1
	C T1:T0 <= 2^128 - 2*2^64 + 1 + 2*2^64 - 3 <= 2^128 - 2, no carry!

	mulld   r8, r0, r11		C U*m[j] low
	mulhdu  r9, r0, r11		C U*m[j] high
	addc    r8, r8, r13		C add T0 and low word
	ld      r0, J8`'(r5)		C y[j+1]
	adde    r13, r9, r14		C add high word with carry to T1
	addze   r10, r16		C carry to CY
	std     r8, JM8`'(r1)		C store tmp[j-1]
	C CY:T1:T0 <= 2^128 - 2 + 2^128 - 2*2^64 + 1 <=
	C             2 * 2^128 - 2*2^64 - 1 ==> CY:T1 <= 2 * 2^64 - 3

')dnl end forloop
C Pass for j = eval(LENGTH - 1). Don't fetch new data from y[j+1].
define(`J', `eval(8*LENGTH - 8)')dnl
define(`JM8', `eval(J - 8)')dnl

	mulld   r8, r0, r12		C x[i]*y[j] low half
	mulhdu  r9, r0, r12		C x[i]*y[j] high half
	ld      r0, J`'(r6)		C m[j]
	addc    r13, r8, r13		C add low word to T0
	adde    r14, r9, r10		C add high word with carry + CY to T1
	C T1:T0 <= 2^128 - 2*2^64 + 1 + 2*2^64 - 3 <= 2^128 - 2, no carry!

	mulld   r8, r0, r11		C U*m[j] low
	mulhdu  r9, r0, r11		C U*m[j] high
	addc    r8, r8, r13		C add T0 and low word
	adde    r13, r9, r14		C add high word with carry to T1
	std     r8, JM8`'(r1)		C store tmp[len-2]
	addze   r15, r16		C put carry in r15 (tmp[len] <= 1)
	std     r13, J`'(r1)		C store tmp[len-1]


C #########################################################################
C # i > 0 passes
C #########################################################################

define(`LM1', `eval(LENGTH - 1)')dnl

	li      r9, LM1			C outer loop count
	mtctr   r9

1:

C Pass for j = 0. We need to fetch x[i], tmp[i] and tmp[i+1] from memory
C and compute the new u

	ldu     r12, 8(r4)		C x[i]
	ld      r0, 0(r5)		C y[0]
	ld      r13, 0(r1)		C tmp[0]
	mulld   r8, r0, r12		C x[i]*y[0] low half
	ld      r14, 8(r1)		C tmp[1]
	mulhdu  r9, r0, r12		C x[i]*y[0] high half
	addc    r13, r8, r13		C T0
	ld      r0, 0(r6)		C m[0]
	mulld   r11, r7, r13		C U = T0*invm mod 2^64
	adde    r14, r9, r14		C T1
	mulld   r8, r0, r11		C U*m[0] low
	addze   r10, r16		C CY
	mulhdu  r9, r0, r11		C U*m[0] high
	ld      r0, 8(r5)		C y[1]
	addc    r8, r8, r13		C result = 0
	adde    r13, r9, r14		C T0, carry pending
	C cy:T1:T0 <= 2*(2^64 - 1)^2 + 2^128 - 1 = 3*2^128 - 4*2^64 + 1,
	C so cy:T1 <= 3*2^64 - 4

forloop(`UNROLL', 1, eval(LENGTH - 2), `dnl
define(`J', `eval(8 * UNROLL)')dnl
define(`J8', `eval(J + 8)')dnl
define(`JM8', `eval(J - 8)')dnl
C Pass for j = UNROLL

	ld      r14, J8`'(r1)		C tmp[j+1]
	mulld   r8, r0, r12		C x[i]*y[j] low half
	adde    r14, r14, r10		C tmp[j+1] + CY + pending carry
	addze   r10, r16		C carry to CY
	mulhdu  r9, r0, r12		C x[i]*y[j] high half
	ld      r0, J`'(r6)		C m[j]
	addc    r13, r8, r13		C add low word to T0
	mulld   r8, r0, r11		C U*m[j] low
	adde    r14, r9, r14		C add high to T1
	addze   r10, r10		C add carry to CY
	mulhdu  r9, r0, r11		C U*m[j] high
	addc    r8, r8, r13		C add T0 and low word
	ld      r0, J8`'(r5)		C y[j+1]
	adde    r13, r9, r14		C T1, carry pending
	std     r8, JM8`'(r1)		C store tmp[j-1]
	C CY:T1:T0 <= 2*(2^64 - 1)^2 + (3*2^64 - 3) + (2^64-1)*2^64
	C          <= 3*2^128 - 2*2^64 - 1 ==> CY:T1 <= 3*2^64 - 3

')dnl end forloop
C Pass for j = eval(LENGTH - 1). Don't fetch new data from y[j+1].
define(`J', `eval(8*LENGTH - 8)')dnl
define(`JM8', `eval(J - 8)')dnl

	mulld   r8, r0, r12		C x[i]*y[j] low half
	adde    r14, r15, r10		C T1 = tmp[len] + CY + pending carry
	C since tmp[len] <= 1, T1 <= 3 and carry is zero
	mulhdu  r9, r0, r12		C x[i]*y[j] high half
	ld      r0, J`'(r6)		C m[j]
	addc    r13, r8, r13		C add low word to T0
	mulld   r8, r0, r11		C U*m[j] low
	adde    r14, r9, r14		C add high to T1
	addze   r10, r16		C CY
	mulhdu  r9, r0, r11		C U*m[j] high
	addc    r8, r8, r13		C add T0 and low word
	adde    r13, r9, r14		C T1, carry pending
	std     r8, JM8`'(r1)		C store tmp[len-2]
	addze   r15, r10		C store tmp[len] <= 1
	std     r13, J`'(r1)		C store tmp[len-1]
	C CY:T1:T0 <= 2*(2^64 - 1)^2 + (3*2^64 - 3) + 2^64
	C          <= 2*2^128 - 1 ==> CY:T1 <= 2*2^64 - 1 = 1:(2^64-1)

	bdnz 1b

C Copy result from tmp memory to z

dnl ==== THIS LOOP WILL NOT WORK FOR LENGTH <= 1 ====
forloop(`UNROLL', 0, eval(LENGTH / 2 - 1), `dnl
define(`J', `eval(UNROLL)')dnl
ifelse(J, `0', dnl
`	ld      r8, 0(r1)', dnl
`	ldu     r8, 8(r1)')
	ldu     r9, 8(r1)
ifelse(J, `0', dnl
`	std     r8, 0(r3)', dnl
`	stdu    r8, 8(r3)')
	stdu    r9, 8(r3)
')dnl
ifelse(eval(LENGTH % 2), 1, `dnl
	ldu     r8, 8(r1)
	stdu    r8, 8(r3)
')dnl

	mr      r3, r15         C return tmp(len)
	ldu     r16, 8(r1)
	ldu     r15, 8(r1)
	ldu     r14, 8(r1)
	ldu     r13, 8(r1)
	addi    r1, r1, 8
	blr
