dnl  Copyright 1999, 2000, 2001, 2002, 2005 Free Software Foundation, Inc.
dnl
dnl  This file is a modified part of the GNU MP Library.
dnl
dnl  The GNU MP Library is free software; you can redistribute it and/or
dnl  modify it under the terms of the GNU Lesser General Public License as
dnl  published by the Free Software Foundation; either version 2.1 of the
dnl  License, or (at your option) any later version.
dnl
dnl  The GNU MP Library is distributed in the hope that it will be useful,
dnl  but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
dnl  Lesser General Public License for more details.
dnl
dnl  You should have received a copy of the GNU Lesser General Public
dnl  License along with the GNU MP Library; see the file COPYING.LIB.  If
dnl  not, write to the Free Software Foundation, Inc., 59 Temple Place -
dnl  Suite 330, Boston, MA 02111-1307, USA.

include(`config.m4')
        TEXT
	GLOBL GSYM_PREFIX`'ecm_redc3
	TYPE(GSYM_PREFIX`'ecm_redc3,`function')

GSYM_PREFIX`'ecm_redc3:
	push	%ebp					# Push registers
	push	%edi
	push	%esi
	push	%ebx
	subl	$16, %esp				# SF: 2 Cpt + Jump +1

	movl	44(%esp), %ecx                          # Read size
	movl	36(%esp), %edi				# Read Dest Ptr
	movl	%ecx, (%esp)				# Save counter
        cmpl    $5, %ecx
        jae     Unroll		
Loop:	
		movl	48(%esp), %ebp			# Read invm
	        movl    40(%esp), %esi                  # Read Source Ptr
		imull	(%edi), %ebp			# Dest[0] * invm
		movl	%edi, 36(%esp)			# Save new Dest
		movl	44(%esp), %ecx			# Read Size (2)
		xorl	%ebx, %ebx			# Initial Carry
InnerLoop:
		        # esi:	  Source
		        # edi:	  Dest
			# ebp:	  Multiplier
			# ecx:	  Counter
		        movl    (%esi), %eax		# U1
			addl    $4, %edi		# V1
			mull    %ebp			# U2
			addl    $4, %esi		# V2
			addl    %ebx, %eax		# U3
		        adcl    $0, %edx		# U4
			addl    %eax, -4(%edi)		# V4
			adcl    $0, %edx		# U5
			decl    %ecx			# V5
			movl    %edx, %ebx		# U6
			jnz     InnerLoop		# V6
		movl	36(%esp), %edi
		movl    %ebx, (%edi)                    # Save final carry
		decl	(%esp)
		lea	4(%edi), %edi			# Advance Dest
		jnz     Loop				# Loop
End:
	addl	$16, %esp
	pop	%ebx
	pop	%esi
	pop	%edi
	pop	%ebp
	ret
	
Unroll:
# %ecx Read size // %edi Dest Ptr
	# Precalcul du saut 
	movl    %ecx, %edx
        decl    %ecx
	subl    $2, %edx
	negl    %ecx	
	shrl    $4, %edx
	andl    $15, %ecx
	movl    %edx, 8(%esp)				# Org Cpt of 4(%esp)
	movl    %ecx, %edx
	shll    $4, %edx
	negl    %ecx
        leal    UnrollEntry (%edx, %ecx,1), %edx
	movl	%ecx, 44(%esp)				# (-size)%16
	movl	%edx, 12(%esp)				# Org PC inside	

UnrollLoop:	
                movl    48(%esp), %ebp                  # Read invm
                movl    40(%esp), %esi                  # Read Source Ptr
                imull    (%edi), %ebp                   # Dest[0] * invm
                movl    %edi, 36(%esp)                  # Save new Dest
                movl    44(%esp), %ecx                  # Read Size %16
		movl    8(%esp), %edx			# Read InnerLoop Cpt
		movl	%edx, 4(%esp)			# Set InnerLoop Cpt
	
		# First mull and set initial carry
	        movl    (%esi), %eax
	        leal    4(%esi,%ecx,4), %esi
	        mull    %ebp
		leal    (%edi,%ecx,4), %edi
	        movl    %edx, %ebx
	
		# Do the Jump inside the unrolling loop
		# And set up the registers differently if odd
	        movl    12(%esp), %edx
	        testl   $1, %ecx
	        movl    %eax, %ecx
		cmovnz  %ebx, %ecx
	        cmovnz  %eax, %ebx	
	        jmp     *%edx
		
		        # eax   scratch
			# ebx   carry hi
			# ecx   carry lo
			# edx   scratch
			# esi   src
			# edi   dst
			# ebp   multiplier

	       .align  32, 0x90
UnrollInnerLoop:	
		addl    $64, %edi
UnrollEntry:	
#	        movl    0(%esi), %eax # Can't use this instruction
	        .byte   0x8b,0x46,0x00
	        mull    %ebp
#	        addl    %ecx, 0(%edi) # Can't use this instruction
	        .byte   0x01,0x4f,0x00
	        adcl    %eax, %ebx
	        movl    %edx, %ecx
	        adcl    $0, %ecx

	        movl    4(%esi), %eax
	        mull    %ebp
	        addl    %ebx, 4(%edi)
	        adcl    %eax, %ecx
	        movl    %edx, %ebx
	        adcl    $0, %ebx

	        movl    8(%esi), %eax
	        mull    %ebp
	        addl    %ecx, 8(%edi)
	        adcl    %eax, %ebx
	        movl    %edx, %ecx
	        adcl    $0, %ecx

	        movl    12(%esi), %eax
	        mull    %ebp
	        addl    %ebx, 12(%edi)
	        adcl    %eax, %ecx
	        movl    %edx, %ebx
	        adcl    $0, %ebx

	        movl    16(%esi), %eax
	        mull    %ebp
	        addl    %ecx, 16(%edi)
	        adcl    %eax, %ebx
	        movl    %edx, %ecx
	        adcl    $0, %ecx

	        movl    20(%esi), %eax
	        mull    %ebp
	        addl    %ebx, 20(%edi)
	        adcl    %eax, %ecx
	        movl    %edx, %ebx
	        adcl    $0, %ebx

	        movl    24(%esi), %eax
	        mull    %ebp
	        addl    %ecx, 24(%edi)
	        adcl    %eax, %ebx
	        movl    %edx, %ecx
	        adcl    $0, %ecx
	
	        movl    28(%esi), %eax
	        mull    %ebp
	        addl    %ebx, 28(%edi)
	        adcl    %eax, %ecx
	        movl    %edx, %ebx
	        adcl    $0, %ebx

	        movl    32(%esi), %eax
	        mull    %ebp
	        addl    %ecx, 32(%edi)
	        adcl    %eax, %ebx
	        movl    %edx, %ecx
	        adcl    $0, %ecx

	        movl    36(%esi), %eax
	        mull    %ebp
	        addl    %ebx, 36(%edi)
	        adcl    %eax, %ecx
	        movl    %edx, %ebx
	        adcl    $0, %ebx

	        movl    40(%esi), %eax
	        mull    %ebp
	        addl    %ecx, 40(%edi)
	        adcl    %eax, %ebx
	        movl    %edx, %ecx
	        adcl    $0, %ecx

	        movl    44(%esi), %eax
	        mull    %ebp
	        addl    %ebx, 44(%edi)
	        adcl    %eax, %ecx
	        movl    %edx, %ebx
	        adcl    $0, %ebx

	        movl    48(%esi), %eax
	        mull    %ebp
	        addl    %ecx, 48(%edi)
	        adcl    %eax, %ebx
	        movl    %edx, %ecx
	        adcl    $0, %ecx

	        movl    52(%esi), %eax
	        mull    %ebp
	        addl    %ebx, 52(%edi)
	        adcl    %eax, %ecx
	        movl    %edx, %ebx
	        adcl    $0, %ebx

	        movl    56(%esi), %eax
	        mull    %ebp
	        addl    %ecx, 56(%edi)
	        adcl    %eax, %ebx
	        movl    %edx, %ecx
	        adcl    $0, %ecx

	        movl    60(%esi), %eax
	        mull    %ebp
	        addl    %ebx, 60(%edi)
	        adcl    %eax, %ecx
	        movl    %edx, %ebx
	        adcl    $0, %ebx

	        decl    4(%esp)
	        leal    64(%esi), %esi
	        jns     UnrollInnerLoop

	        addl    %ecx, 64(%edi)
	        movl    36(%esp), %edi
	        adcl    $0, %ebx
                movl    %ebx, (%edi)                    # Save final carry
                decl    (%esp)
                lea     4(%edi), %edi                   # Advance Dest
                jnz     UnrollLoop                      # Loop
End2:	
        addl    $16, %esp
        pop     %ebx
        pop     %esi
        pop     %edi
        pop     %ebp
        ret
