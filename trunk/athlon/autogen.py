#!/usr/bin/python

import re
import sys


def offaddr(addr, offset):
	if offset == 0:
		return "("+addr+")"
	else:
		return str(offset)+"("+addr+")"

# Generate asm for addmul1_k
# src and dst are pointers (stored in regs) + offsets
# multiplier is in a register
# rax, rbx, rcx, rdx are free for use.

def addmul1_k(src, off_src, dst, off_dst, mult, k):
	init = "### addmul1: src[0] is " + offaddr(src, off_src) + "\n"
	init = init + "###          dst[0] is " + offaddr(dst, off_dst) + "\n"
	init = init + "###          mult is " + mult + "\n"
	init = init + "###          k is " + str(k) + "\n"
	init = init + "###          kills %eax, %ebx, %ecx, %edx\n"
	init = init + "###   dst[0,k[ += mult*src[0,k[  plus carry put in ecx or ebx\n"
	init = init + "	movl	" + offaddr(src, off_src) + ", %eax\n"
	init = init + "	mull	" + mult + "\n"
	init = init + "	movl	%eax, %ebx\n"
	init = init + "	movl	%edx, %ecx\n"
	init = init + "	movl	" + offaddr(src, off_src+4) + ", %eax\n"

	block = """
	mull	__mult__
	addl	__cylo__, __zi__
	movl	$0, __cylo__
	adcl	%eax, __cyhi__
	movl	__xi2__, %eax
	adcl	%edx, __cylo__
"""
	
	code = init
	
	cylo = "%ebx"
	cyhi = "%ecx"
	for i in range(0,k-2):
		blocki = re.sub('__cylo__', cylo, block)
		blocki = re.sub('__cyhi__', cyhi, blocki)
		blocki = re.sub('__xi2__', offaddr(src, off_src+(i+2)*4), blocki)
		blocki = re.sub('__zi__', offaddr(dst, off_dst+i*4), blocki)
		blocki = re.sub('__mult__', mult, blocki)
		code = code + blocki
		tmp = cylo
		cylo = cyhi
		cyhi = tmp
	
	final = "	mull	" + mult + "\n"
	final = final + "	addl	" + cylo + ", " + offaddr(dst, off_dst+(k-2)*4) + "\n"
	final = final + "	adcl	" + cyhi + ", %eax\n"
	final = final + "	adcl	$0, %edx\n"
	final = final + "	addl	%eax, " + offaddr(dst, off_dst+4*(k-1)) + "\n"
	final = final + "	adcl	$0, %edx\n"
	final = final + "### carry limb is in %edx\n"

	code = code + final
	return code, "%edx"

### Try mmx/sse2 addmul_1, copying the one of GMP for Pentium4

def addmul1_k_var(src, off_src, dst, off_dst, mult, k):
	init = "### addmul1: src[0] is " + offaddr(src, off_src) + "\n"
	init = init + "###          dst[0] is " + offaddr(dst, off_dst) + "\n"
	init = init + "###          mult is " + mult + "\n"
	init = init + "###          k is " + str(k) + "\n"
	init = init + "###          kills %eax, %edx and mmx regs \n"
	init = init + "###   dst[0,k[ += mult*src[0,k[  plus carry put in ecx\n"
	init = init + "	pxor	%mm0, %mm0\n"
	init = init + "	movd	" + mult + ", %mm7\n"

	block = """
	movd	__xi__, %mm1
	movd	__zi__, %mm2
	pmuludq	%mm7, %mm1
	paddq	%mm1, %mm2
	paddq	%mm2, %mm0
	movd	%mm0, __zi__
	psrlq	$32, %mm0
"""
	
	code = init
	
	for i in range(0,k):
		blocki = re.sub('__xi__', offaddr(src, off_src+i*4), block)
		blocki = re.sub('__zi__', offaddr(dst, off_dst+i*4), blocki)
		code = code + blocki
	
	final = "	movd	%mm0, %ecx\n"
	final = final + "### carry limb is in %ecx\n"

	code = code + final
	return code, "%ecx"


def mulredc_k_rolled(k):
	header = """# mp_limb_t mulredc__k(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y,
#                 const mp_limb_t *m, mp_limb_t inv_m);
#
#  Stack:
#    inv_m    ## parameters
#    m
#    y
#    x
#    z							(4*(2k+7))%esp
#    ???   (1 limb???)
#    ebp      ## pushed registers                  (4*(2k+5))%esp
#    edi
#    esi
#    ebx
#    ...      ## counter (1 mp_limb_t)             (4*(2k+1))%esp
#    ...      ## tmp space (2*k+1 mp_limb_t)

include(`config.m4')
	TEXT
	GLOBL GSYM_PREFIX`'mulredc__k
	TYPE(GSYM_PREFIX`'mulredc__k,`function')

GSYM_PREFIX`'mulredc__k:
"""
	init = re.sub("__k", str(k), header)

	INV_M  = offaddr("%esp", 4*(2*k+1) + 40)
	ADDR_M = offaddr("%esp", 4*(2*k+1) + 36)
	ADDR_Y = offaddr("%esp", 4*(2*k+1) + 32)
	ADDR_X = offaddr("%esp", 4*(2*k+1) + 28)
	ADDR_Z = offaddr("%esp", 4*(2*k+1) + 24)

	init = init + """	pushl	%ebp
	pushl	%edi
	pushl	%esi
	pushl	%ebx
"""
	init = init + "	subl	$" + str(4*(2*k+2)) + ", %esp\n"
	init = init + "	movl	%esp, %edi\n"
	init = init + "### set tmp[0..2k+1[ to 0\n"
	for i in range(0,2*k+1):
		init = init + "	movl	$0, " + offaddr("%edi", 4*i) + "\n"
	
	code = init
	
	middle_code = "###########################################\n"
	middle_code = middle_code + "	movl	$" + str(k) + ", " + offaddr("%esp", 4*(2*k+1)) + "\n"
	middle_code = middle_code + """
	.align 32
Loop:
	## compute u and store in %ebp
"""
	middle_code = middle_code + "	movl	" + ADDR_X + ", %eax\n"
	middle_code = middle_code + "	movl	" + ADDR_Y + ", %esi\n"
	middle_code = middle_code + """	movl	(%eax), %eax
	mull	(%esi)
	addl	(%edi), %eax
"""
	middle_code = middle_code + "	mull	" + INV_M + "\n"
	middle_code = middle_code + "	movl    %eax, %ebp\n"

	middle_code = middle_code + "	movl	" + ADDR_M + ", %esi\n"
	codeaddmul, carry = addmul1_k("%esi", 0, "%edi", 0, "%ebp", k)
	middle_code = middle_code + codeaddmul
	middle_code = middle_code + "	addl	" + carry + ", " + offaddr("%edi", 4*k) + "\n"
	middle_code = middle_code + "	adcl	$0, " + offaddr("%edi", 4*(k+1)) + "\n"
	middle_code = middle_code + "	movl	" + ADDR_X + ", %eax\n"
	middle_code = middle_code + "	movl	(%eax), %ebp\n"
	middle_code = middle_code + "	movl	" + ADDR_Y + ", %esi\n"
	codeaddmul, carry = addmul1_k("%esi", 0, "%edi", 0, "%ebp", k)
	middle_code = middle_code + codeaddmul
	middle_code = middle_code + "   addl    " + carry + ", " + offaddr("%edi", 4*k) + "\n"
	middle_code = middle_code + "   adcl    $0, " + offaddr("%edi", 4*(k+1)) + "\n\n"
	middle_code = middle_code + "	addl	$4, " + ADDR_X + "\n	addl	$4, %edi\n"
	middle_code = middle_code + "	decl	" + offaddr("%esp", 4*(2*k+1)) + "\n	jnz	Loop\n"
	code = code + middle_code

	final = "###########################################\n"
	final = final + "### Copy result in z\n"
	final = final + "	movl	" + ADDR_Z + ", %ebx\n"
	for i in range(0,k):
		final = final + "	movl	" + offaddr("%edi", 4*i) + ", %eax\n"
		final = final + "	movl	%eax, " + offaddr("%ebx", 4*i) + "\n"
	final = final + "	movl	" + offaddr("%edi", 4*k) + ", %eax	# carry\n"
	final = final + "	addl    $" + str(4*(2*k+2)) + ", %esp\n"
	final = final + "	popl	%ebx\n"
	final = final + "	popl	%esi\n"
	final = final + "	popl	%edi\n"
	final = final + "	popl	%ebp\n"
#	final = final + "	emms\n"
	final = final + "	ret\n"

	code = code + final
	
	return code

	
k = int(sys.argv[1])
if k == 1:
	print """#
#  mp_limb_t mulredc1(mp_limb_t *z, const mp_limb_t x, const mp_limb_t y,
#                 const mp_limb_t m, mp_limb_t inv_m)
#
#  Compute z := x*y mod m, in Montgomery representation, where x, y < m
#  and m is n limb wide.  inv_m is the less significant limb of the
#  inverse of m modulo 2^(n*GMP_LIMB_BITS)
#
#  The result might be unreduced (larger than m) but becomes reduced
#  after subtracting m. The calling function should take care of that.
#
#  We use a temporary space for unreduced product on the stack.
#  Therefore, this can not be used for large integers (anyway, the
#  algorithm is quadratic).
#
#  WARNING: z is only n limbs but since it might be unreduced, there
#  could be a carry that does not fit in z. This carry is returned.

include(`config.m4')
	TEXT
	GLOBL GSYM_PREFIX`'mulredc1
	TYPE(GSYM_PREFIX`'mulredc1,`function')

GSYM_PREFIX`'mulredc1:
# Stack:
#    inv_m  20(%esp)
#    m      16
#    y      12(%esp)
#    x      8
#    z      4(%esp)

	movl	12(%esp), %eax
	mull	8(%esp)
	movl	%edx, 12(%esp)
	movl	%eax, 8(%esp)   # store xy in [8(%esp):12(%esp)]
	mull	20(%esp)          # compute u
	mull	16(%esp)         # compute u*m
	addl	8(%esp), %eax       # eax is 0, now (carry is important)
	adcl	12(%esp), %edx
	movl	4(%esp), %ecx
	movl    %edx, (%ecx)
	adcl	$0, %eax
	ret
"""
else:
	print mulredc_k_rolled(k)

