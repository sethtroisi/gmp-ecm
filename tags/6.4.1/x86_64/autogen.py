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
	init = init + "###          kills %rax, %rbx, %rcx, %rdx\n"
	init = init + "###   dst[0,k[ += mult*src[0,k[  plus carry put in rcx or rbx\n"
	init = init + "	movq	" + offaddr(src, off_src) + ", %rax\n"
	init = init + "	mulq	" + mult + "\n"
	init = init + "	movq	%rax, %rbx\n"
	init = init + "	movq	%rdx, %rcx\n"

	block = """
	movq	__xii__, %rax
	mulq	__mult__
	addq	__cylo__, __zi__
	adcq	%rax, __cyhi__
	movq	%rdx, __cylo__
	adcq	$0, __cylo__
"""
	
	code = init
	
	cylo = "%rbx"
	cyhi = "%rcx"
	for i in range(0,k-1):
		blocki = re.sub('__cylo__', cylo, block)
		blocki = re.sub('__cyhi__', cyhi, blocki)
		blocki = re.sub('__xii__', offaddr(src, off_src+(i+1)*8), blocki)
		blocki = re.sub('__zi__', offaddr(dst, off_dst+i*8), blocki)
		blocki = re.sub('__mult__', mult, blocki)
		code = code + blocki
		tmp = cylo
		cylo = cyhi
		cyhi = tmp
	
	final = "	addq	" + cylo + ", " + offaddr(dst, off_dst+8*(k-1)) + "\n"
	final = final + "	adcq	$0, " + cyhi + "\n"
	final = final + "### carry limb is in " + cyhi + "\n"

	code = code + final
	return code, cyhi


######## TODO: improve this code!!!!

def mul1_k(src, off_src, dst, off_dst, mult, k):
	init = "### mul1: src[0] is " + offaddr(src, off_src) + "\n"
	init = init + "###          dst[0] is " + offaddr(dst, off_dst) + "\n"
	init = init + "###          mult is " + mult + "\n"
	init = init + "###          k is " + str(k) + "\n"
	init = init + "###          kills %rax, %rbx, %rcx, %rdx\n"
	init = init + "###   dst[0,k[ = mult*src[0,k[  plus carry put in rcx or rbx\n"
	init = init + "	movq	" + offaddr(src, off_src) + ", %rax\n"
	init = init + "	mulq	" + mult + "\n"
	init = init + "	movq	%rax, %rbx\n"
	init = init + "	movq	%rdx, %rcx\n"

	block = """
	movq	__xii__, %rax
	mulq	__mult__
	movq	__cylo__, __zi__
	addq	%rax, __cyhi__
	movq	%rdx, __cylo__
	adcq	$0, __cylo__
"""
	
	code = init
	
	cylo = "%rbx"
	cyhi = "%rcx"
	for i in range(0,k-1):
		blocki = re.sub('__cylo__', cylo, block)
		blocki = re.sub('__cyhi__', cyhi, blocki)
		blocki = re.sub('__xii__', offaddr(src, off_src+(i+1)*8), blocki)
		blocki = re.sub('__zi__', offaddr(dst, off_dst+i*8), blocki)
		blocki = re.sub('__mult__', mult, blocki)
		code = code + blocki
		tmp = cylo
		cylo = cyhi
		cyhi = tmp
	
	final = "	movq	" + cylo + ", " + offaddr(dst, off_dst+8*(k-1)) + "\n"
	final = final + "### carry limb is in " + cyhi + "\n"

	code = code + final
	return code


def mulredc_k_rolled(k):
	header = """# mp_limb_t mulredc__k(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y,
#                 const mp_limb_t *m, mp_limb_t inv_m);
#

include(`config.m4')
	TEXT
	GLOBL GSYM_PREFIX`'mulredc__k
	TYPE(GSYM_PREFIX`'mulredc__k,`function')

GSYM_PREFIX`'mulredc__k:
"""
	init = re.sub("__k", str(k), header)
  
	init = init + """	movq	%rdx, %r11
	movq	%rcx, %r10
	pushq	%rbx
	pushq	%rbp
"""
	init = init + "	subq	$" + str(8*(2*k+1)) + ", %rsp\n"
	init = init + """#      %r8  : inv_m
#     %r10 : m
#     %r11 : y
#     %rsi : x
#     %rdi : z
#     %rsp : tmp
# Free registers
#     %rax, %rbx, %rcx, %rdx, %r9

### set tmp[0..2k+1[ to 0
"""
	for i in range(0,2*k+1):
		init = init + "	movq	$0, " + offaddr("%rsp", 8*i) + "\n"
	
	code = init
	
	middle_code = "###########################################\n"
	middle_code = middle_code + "	movq	$" + str(k) + ", %rbp\n"
	middle_code = middle_code + """
	.align 64
Loop:
	## compute u and store in %r9
	movq	(%rsi), %rax
	mulq	(%r11)
	addq	(%rsp), %rax
	mulq	%r8
	movq    %rax, %r9
"""
	codeaddmul, carry = addmul1_k("%r10", 0, "%rsp", 0, "%r9", k)
	middle_code = middle_code + codeaddmul
	middle_code = middle_code + "	addq	" + carry + ", " + offaddr("%rsp", 8*k) + "\n"
	middle_code = middle_code + "	adcq	$0, " + offaddr("%rsp", 8*(k+1)) + "\n"
	middle_code = middle_code + "	movq	(%rsi), %r9\n"
	codeaddmul, carry = addmul1_k("%r11", 0, "%rsp", 0, "%r9", k)
	middle_code = middle_code + codeaddmul
	middle_code = middle_code + "   addq    " + carry + ", " + offaddr("%rsp", 8*k) + "\n"
	middle_code = middle_code + "   adcq    $0, " + offaddr("%rsp", 8*(k+1)) + "\n\n"
	middle_code = middle_code + """
	addq	$8, %rsi
	addq	$8, %rsp
	decq	%rbp
	jnz	Loop
"""
	code = code + middle_code

	final = "###########################################\n"
	final = final + "### Copy result in z\n"
	for i in range(0,k):
		final = final + "	movq	" + offaddr("%rsp", 8*i) + ", %rax\n"
		final = final + "	movq	%rax, " + offaddr("%rdi", 8*i) + "\n"
	final = final + "	movq	" + offaddr("%rsp", 8*k) + ", %rax	# carry\n"
	final = final + "	addq    $" + str(8*(k+1)) + ", %rsp\n"
	final = final + "	popq	%rbp\n"
	final = final + "	popq	%rbx\n"
	final = final + "	ret\n"

	code = code + final
	
	return code



def mulredc_k(k):
	header = """# mp_limb_t mulredc__k(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y,
#                 const mp_limb_t *m, mp_limb_t inv_m);
#

include(`config.m4')
	TEXT
	GLOBL GSYM_PREFIX`'mulredc__k
	TYPE(GSYM_PREFIX`'mulredc__k,`function')

GSYM_PREFIX`'mulredc__k:
"""
	init = re.sub("__k", str(k), header)
  
	init = init + """	movq	%rdx, %r11
	movq	%rcx, %r10
	pushq	%rbx
"""
	init = init + "	subq	$" + str(8*(2*k+1)) + ", %rsp\n"
	init = init + """#      %r8  : inv_m
#     %r10 : m
#     %r11 : y
#     %rsi : x
#     %rdi : z
#     %rsp : tmp
# Free registers
#     %rax, %rbx, %rcx, %rdx, %r9

### set tmp[0..2k+1[ to 0
"""
	for i in range(0,2*k+1):
		init = init + "	movq	$0, " + offaddr("%rsp", 8*i) + "\n"
	
	code = init
	
	for i in range(0,k):
		blocki = "###########################################\n"
		blocki = blocki + "### Step " + str(i) + "\n"
		blocki = blocki + "### Compute u and store in %r9\n"
		blocki = blocki + "	movq	" + offaddr("%rsi", 8*i) + ", %rax\n"
		blocki = blocki + "	mulq	(%r11)\n"
		blocki = blocki + "	addq	" + offaddr("%rsp", 8*i) + ", %rax\n"
		blocki = blocki + "	mulq	%r8\n"
		blocki = blocki + "	movq	%rax, %r9\n"
		blocki = blocki + "### tmp[i,i+k] += x[i]*y + u*m\n"
		codeaddmul, carry = addmul1_k("%r10", 0, "%rsp", 8*i, "%r9", k)
		blocki = blocki + codeaddmul
		blocki = blocki + "	addq	" + carry + ", " + offaddr("%rsp", 8*(k+i)) + "\n"
		blocki = blocki + "	adcq	$0, " + offaddr("%rsp", 8*(k+i+1)) + "\n"
		blocki = blocki + "	movq	" + offaddr("%rsi", 8*i) + ", %r9\n"
		codeaddmul, carry = addmul1_k("%r11", 0, "%rsp", 8*i, "%r9", k)
		blocki = blocki + codeaddmul
		blocki = blocki + "	addq	" + carry + ", " + offaddr("%rsp", 8*(k+i)) + "\n"
		blocki = blocki + "	adcq	$0, " + offaddr("%rsp", 8*(k+i+1)) + "\n"
		code = code + blocki
	
	final = "###########################################\n"
	final = final + "### Copy result in z\n"
	for i in range(0,k):
		final = final + "	movq	" + offaddr("%rsp", 8*(k+i)) + ", %rax\n"
		final = final + "	movq	%rax, " + offaddr("%rdi", 8*i) + "\n"
	final = final + "	movq	" + offaddr("%rsp", 16*k) + ", %rax	# carry\n"
	final = final + "	addq    $" + str(8*(2*k+1)) + ", %rsp\n"
	final = final + "	popq	%rbx\n"
	final = final + "	ret\n"

	code = code + final
	
	return code

	
##print addmul1_k("%rsi", 0, "%dsi", 0, "%r9", 3)

k = int(sys.argv[1])
if k == 1:
	print """#
#  mp_limb_t mulredc1(mp_limb_t * z, const mp_limb_t x, const mp_limb_t y,
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
#     %r8  : inv_m
#     %rcx : m
#     %rdx : y
#     %rsi : x
#     %rdi : z
	movq	%rdx, %rax
	mulq	%rsi
	movq	%rdx, %r10
	movq	%rax, %r9       # store xy in [r9:r10]
	mulq	%r8             # compute u
	mulq	%rcx          # compute u*m
	addq	%r9, %rax       # rax is 0, now (carry is important)
	adcq	%r10, %rdx
	movq	%rdx, (%rdi)
	adcq	$0, %rax
	ret
"""
else:
	print mulredc_k_rolled(k)

