dnl ******************************************************************************
dnl   Copyright 2009 Paul Zimmermann and Alexander Kruppa.
dnl 
dnl   This file is part of the ECM Library.
dnl
dnl  The ECM Library is free software; you can redistribute it and/or modify
dnl  it under the terms of the GNU Lesser General Public License as published by
dnl  the Free Software Foundation; either version 2.1 of the License, or (at your
dnl  option) any later version.
dnl
dnl  The ECM Library is distributed in the hope that it will be useful, but
dnl  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
dnl  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
dnl  License for more details.
dnl
dnl  You should have received a copy of the GNU Lesser General Public License
dnl  along with the ECM Library; see the file COPYING.LIB.  If not, write to
dnl  the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
dnl  MA 02110-1301, USA.
dnl ******************************************************************************
dnl 
dnl  void ecm_redc3(mp_limb_t * c, const mp_limb_t * m, size_t n, mp_limb_t m_inv)
dnl 
dnl  input arguments:
dnl 
dnl  r3: ptr to c[0], the least significant word of the number to be reduced
dnl      c[0 ... 2*n-1] is of length 2*n words
dnl  r4: ptr to m[0], the least significant word of the modulus m of length n
dnl  r5: the length n
dnl  r6: m_inv = -1/m mod 2^64
dnl 
dnl  the residue (before adding the word carries) will be in c[n ... 2*n-1].
dnl  c[0 ... n-1] will contain the high word carries from each inner loop pass.
dnl  These carry words are added by the calling routine to obtain the final
dnl  residue.

dnl Use `C' to remove comments in .asm -> .s conversion.
dnl Copied from GMP 4.2.
define(C, `
dnl')

include(`config.m4')
        TEXT
        GLOBL GSYM_PREFIX`'ecm_redc3
        TYPE(GSYM_PREFIX`'ecm_redc3,`function')

GSYM_PREFIX`'ecm_redc3:

	cmpdi	r5, 1				C length = 1?
	bne		1f

    ld      r12, 0(r3)          C c[0]
    ld      r0, 0(r4)           C m[0]
    mulld   r7, r6, r12         C u = c[0] * m_inv mod 2^64
    mulld   r11, r0, r7         C m[0]*u low
    mulhdu  r10, r0, r7         C m[0]*u high
	addc	r11, r11, r12		C c[0] + m[0]*u low = 0
	addze	r10, r10			C carry to high half
	std		r10, 0(r3)			C store the "carry" word
	blr

    nop
	nop
	nop
	nop
	nop
1:
    mflr    r0                  C save return addr
    stdu    r0, -8(r1)          C on the stack
    stdu    r13, -8(r1)         C save r13
dnl 
dnl        get inner loop count and jump offset
dnl 
    subi    r7, r5, 2           C r7 = n - 2
    andi.   r8, r7, 15          C r8 =  (n - 2) mod 16
    sldi    r8, r8, 5           C r8 * 32 = byte offset
    srdi    r7, r7, 4           C int((n - 2)/16)
dnl 
dnl  compute the address of inner loop end and subtract the offset
dnl 
    bl      nxt                 C put the address of the next instruction
    							C into the link register
nxt:							C
    mflr    r9                  C r9 = address of this instruction
    addi    r9, r9, 640  	C add offset to v_1 from nxt
	                        C WARNING: any changes to the code between
	                        C the labels "nxt" and "v_1" may require
	                        C recomputation of the offset above.
    sub     r9, r9, r8          C offset back to desired starting point
    mtlr    r9                  C and now we can branch directly to our target
    mtctr   r5                  C outer loop count n
    addi    r13, r7, 1          C inner loop counter

	nop
	nop

OuterLoop:                      C execute n times

dnl  compute u, set addr's
    
    ld      r12, 0(r3)          C c[0]
    mr      r8, r4              C r8 = working copy of m address
    ld      r0, 0(r8)           C m[0]
    mulld   r7, r6, r12         C u = c[0] * m_inv mod 2^64
    mfctr   r5                  C save current outer loop count

dnl  start inner
    mulld   r11, r0, r7         C m[0]*u low
    mtctr   r13                 C inner loop count
    mulhdu  r10, r0, r7         C m[0]*u high
    ldu     r0, 8(r8)           C m[1]
    addc    r11, r11, r12       C m[0]*u low + c[0] (don't bother storing zero)
    mulld   r11, r0, r7         C m[1]*u low
    ldu     r12, 8(r3)          C c[1], update c address
    mr      r9, r3              C r9 = working copy of c addr
    mulhdu  r0, r0, r7          C m[1]*u high
    adde    r11, r10, r11       C m[1]*u low + m[0]*u high + cy
    addze   r10, r0             C m[1]*u high + cy

    blr                         C jump to start of the (n-2) mod 16 section
								C (or to v_1, if (n-2) mod 16 = 0)
    nop
	nop
	nop
	nop
	nop
	nop
	nop

ILoop:

    ldu     r0, 8(r8)           C m[i]
    addc    r11, r11, r12       C m[i-1]*u low + m[i-2]*u high + c[i-1]
    std     r11, 0(r9)          C store it in c[i-1]
    mulld   r11, r0, r7         C m[i]*u low
    ldu     r12, 8(r9)          C c[i]
    mulhdu  r0, r0, r7          C m[i]*u high
    adde    r11, r10, r11       C m[i]*u low + m[i-1]*u high + cy
    addze   r10, r0             C r10 =  m[i]*u + cy
dnl 
dnl  start (n-2) mod 16 = 15
dnl 
    ldu     r0, 8(r8)           C m[i]
    addc    r11, r11, r12       C m[i-1]*u low + m[i-2]*u high + c[i-1]
    std     r11, 0(r9)          C store it in c[i-1]
    mulld   r11, r0, r7         C m[i]*u low
    ldu     r12, 8(r9)          C c[i]
    mulhdu  r0, r0, r7          C m[i]*u high
    adde    r11, r10, r11       C m[i]*u low + m[i-1]*u high + cy
    addze   r10, r0             C r10 =  m[i]*u + cy
dnl 
dnl  start (n-2) mod 16 = 14
dnl 
    ldu     r0, 8(r8)           C m[i]
    addc    r11, r11, r12       C m[i-1]*u low + m[i-2]*u high + c[i-1]
    std     r11, 0(r9)          C store it in c[i-1]
    mulld   r11, r0, r7         C m[i]*u low
    ldu     r12, 8(r9)          C c[i]
    mulhdu  r0, r0, r7          C m[i]*u high
    adde    r11, r10, r11       C m[i]*u low + m[i-1]*u high + cy
    addze   r10, r0             C r10 =  m[i]*u + cy
dnl 
dnl  start (n-2) mod 16 = 13
dnl 
    ldu     r0, 8(r8)           C m[i]
    addc    r11, r11, r12       C m[i-1]*u low + m[i-2]*u high + c[i-1]
    std     r11, 0(r9)          C store it in c[i-1]
    mulld   r11, r0, r7         C m[i]*u low
    ldu     r12, 8(r9)          C c[i]
    mulhdu  r0, r0, r7          C m[i]*u high
    adde    r11, r10, r11       C m[i]*u low + m[i-1]*u high + cy
    addze   r10, r0             C r10 =  m[i]*u + cy
dnl 
dnl  start (n-2) mod 16 = 12
dnl 
    ldu     r0, 8(r8)           C m[i]
    addc    r11, r11, r12       C m[i-1]*u low + m[i-2]*u high + c[i-1]
    std     r11, 0(r9)          C store it in c[i-1]
    mulld   r11, r0, r7         C m[i]*u low
    ldu     r12, 8(r9)          C c[i]
    mulhdu  r0, r0, r7          C m[i]*u high
    adde    r11, r10, r11       C m[i]*u low + m[i-1]*u high + cy
    addze   r10, r0             C r10 =  m[i]*u + cy
dnl 
dnl  start (n-2) mod 16 = 11
dnl 
    ldu     r0, 8(r8)           C m[i]
    addc    r11, r11, r12       C m[i-1]*u low + m[i-2]*u high + c[i-1]
    std     r11, 0(r9)          C store it in c[i-1]
    mulld   r11, r0, r7         C m[i]*u low
    ldu     r12, 8(r9)          C c[i]
    mulhdu  r0, r0, r7          C m[i]*u high
    adde    r11, r10, r11       C m[i]*u low + m[i-1]*u high + cy
    addze   r10, r0             C r10 =  m[i]*u + cy
dnl 
dnl  start (n-2) mod 16 = 10
dnl 
    ldu     r0, 8(r8)           C m[i]
    addc    r11, r11, r12       C m[i-1]*u low + m[i-2]*u high + c[i-1]
    std     r11, 0(r9)          C store it in c[i-1]
    mulld   r11, r0, r7         C m[i]*u low
    ldu     r12, 8(r9)          C c[i]
    mulhdu  r0, r0, r7          C m[i]*u high
    adde    r11, r10, r11       C m[i]*u low + m[i-1]*u high + cy
    addze   r10, r0             C r10 =  m[i]*u + cy
dnl 
dnl  start (n-2) mod 16 = 9
dnl 
    ldu     r0, 8(r8)           C m[i]
    addc    r11, r11, r12       C m[i-1]*u low + m[i-2]*u high + c[i-1]
    std     r11, 0(r9)          C store it in c[i-1]
    mulld   r11, r0, r7         C m[i]*u low
    ldu     r12, 8(r9)          C c[i]
    mulhdu  r0, r0, r7          C m[i]*u high
    adde    r11, r10, r11       C m[i]*u low + m[i-1]*u high + cy
    addze   r10, r0             C r10 =  m[i]*u + cy
dnl 
dnl  start (n-2) mod 16 = 8
dnl 
    ldu     r0, 8(r8)           C m[i]
    addc    r11, r11, r12       C m[i-1]*u low + m[i-2]*u high + c[i-1]
    std     r11, 0(r9)          C store it in c[i-1]
    mulld   r11, r0, r7         C m[i]*u low
    ldu     r12, 8(r9)          C c[i]
    mulhdu  r0, r0, r7          C m[i]*u high
    adde    r11, r10, r11       C m[i]*u low + m[i-1]*u high + cy
    addze   r10, r0             C r10 =  m[i]*u + cy
dnl 
dnl  start (n-2) mod 16 = 7
dnl 
    ldu     r0, 8(r8)           C m[i]
    addc    r11, r11, r12       C m[i-1]*u low + m[i-2]*u high + c[i-1]
    std     r11, 0(r9)          C store it in c[i-1]
    mulld   r11, r0, r7         C m[i]*u low
    ldu     r12, 8(r9)          C c[i]
    mulhdu  r0, r0, r7          C m[i]*u high
    adde    r11, r10, r11       C m[i]*u low + m[i-1]*u high + cy
    addze   r10, r0             C r10 =  m[i]*u + cy
dnl 
dnl  start (n-2) mod 16 = 6
dnl 
    ldu     r0, 8(r8)           C m[i]
    addc    r11, r11, r12       C m[i-1]*u low + m[i-2]*u high + c[i-1]
    std     r11, 0(r9)          C store it in c[i-1]
    mulld   r11, r0, r7         C m[i]*u low
    ldu     r12, 8(r9)          C c[i]
    mulhdu  r0, r0, r7          C m[i]*u high
    adde    r11, r10, r11       C m[i]*u low + m[i-1]*u high + cy
    addze   r10, r0             C r10 =  m[i]*u + cy
dnl 
dnl  start (n-2) mod 16 = 5
dnl 
    ldu     r0, 8(r8)           C m[i]
    addc    r11, r11, r12       C m[i-1]*u low + m[i-2]*u high + c[i-1]
    std     r11, 0(r9)          C store it in c[i-1]
    mulld   r11, r0, r7         C m[i]*u low
    ldu     r12, 8(r9)          C c[i]
    mulhdu  r0, r0, r7          C m[i]*u high
    adde    r11, r10, r11       C m[i]*u low + m[i-1]*u high + cy
    addze   r10, r0             C r10 =  m[i]*u + cy
dnl 
dnl  start (n-2) mod 16 = 4
dnl 
    ldu     r0, 8(r8)           C m[i]
    addc    r11, r11, r12       C m[i-1]*u low + m[i-2]*u high + c[i-1]
    std     r11, 0(r9)          C store it in c[i-1]
    mulld   r11, r0, r7         C m[i]*u low
    ldu     r12, 8(r9)          C c[i]
    mulhdu  r0, r0, r7          C m[i]*u high
    adde    r11, r10, r11       C m[i]*u low + m[i-1]*u high + cy
    addze   r10, r0             C r10 =  m[i]*u + cy
dnl 
dnl  start (n-2) mod 16 = 3
dnl 
    ldu     r0, 8(r8)           C m[i]
    addc    r11, r11, r12       C m[i-1]*u low + m[i-2]*u high + c[i-1]
    std     r11, 0(r9)          C store it in c[i-1]
    mulld   r11, r0, r7         C m[i]*u low
    ldu     r12, 8(r9)          C c[i]
    mulhdu  r0, r0, r7          C m[i]*u high
    adde    r11, r10, r11       C m[i]*u low + m[i-1]*u high + cy
    addze   r10, r0             C r10 =  m[i]*u + cy
dnl 
dnl  start (n-2) mod 16 = 2
dnl 
	ldu     r0, 8(r8)           C m[i]
    addc    r11, r11, r12       C m[i-1]*u low + m[i-2]*u high + c[i-1]
    std     r11, 0(r9)          C store it in c[i-1]
    mulld   r11, r0, r7         C m[i]*u low
    ldu     r12, 8(r9)          C c[i]
    mulhdu  r0, r0, r7          C m[i]*u high
    adde    r11, r10, r11       C m[i]*u low + m[i-1]*u high + cy
    addze   r10, r0             C r10 =  m[i]*u + cy
dnl 
dnl  start (n-2) mod 16 = 1
dnl 
    ldu     r0, 8(r8)           C m[i]
    addc    r11, r11, r12       C m[i-1]*u low + m[i-2]*u high + c[i-1]
    std     r11, 0(r9)          C store it in c[i-1]
    mulld   r11, r0, r7         C m[i]*u low
    ldu     r12, 8(r9)          C c[i]
    mulhdu  r0, r0, r7          C m[i]*u high
    adde    r11, r10, r11       C m[i]*u low + m[i-1]*u high + cy
    addze   r10, r0             C r10 =  m[i]*u + cy
v_1:
    bdnz    ILoop               C blr above jumps directly to this bdnz instruction
                                C when (n-2) mod 16 = 0
dnl  finish inner
    addc    r11, r11, r12       C m[i-1]*u low + m[i-2]*u high + c[i-1]
    std     r11, 0(r9)          C store it in c[i-1]
    addze   r10, r10            C result cy = 0 always
    std     r10, -8(r3)         C store the "carry" word
    mtctr   r5                  C restore outer loop count
    bdnz    OuterLoop

    ld      r13, 0(r1)          C restore r13
    ld      r0, 8(r1)           C original return address
    addi    r1, r1, 16          C restore stack ptr
    mtlr    r0
    blr
