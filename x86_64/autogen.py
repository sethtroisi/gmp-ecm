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

ifdef(`WINDOWS64_ABI',
# stack: inv_m, %r9: m, %r8: y, %rdx: x, %rcx: *z
`define(`INV_M', `0x28(%rsp)')
define(`M', `%r9')
define(`Y', `%r8')
define(`X', `%rdx')
define(`Z', `%rcx')
define(`TMP2', `%r10')
define(`TMP1', `%r8')',
# %r8: inv_m, %rcx: m, %rdx: y, %rsi : x, %rdi : *z
`define(`INV_M', `%r8')
define(`M', `%rcx')
define(`Y', `%rdx')
define(`X', `%rsi')
define(`Z', `%rdi')
define(`TMP2', `%r10')
define(`TMP1', `%r9')')

GSYM_PREFIX`'mulredc1:
	movq	Y, %rax
	mulq	X
	movq	%rdx, TMP2
	movq	%rax, TMP1      # store xy in [r9:r10]
	mulq	INV_M           # compute u
	mulq	M               # compute u*m
	addq	TMP1, %rax      # rax is 0, now (carry is important)
	adcq	TMP2, %rdx
	movq	%rdx, (Z)
	adcq	$0, %rax
	ret

ifdef(`WINDOWS64_ABI',
,
`
`#'if defined(__linux__) && defined(__ELF__)
.section .note.GNU-stack,"",%progbits
`#'endif
') dnl
"""
elif k == 2:
  print """#
# mp_limb_t mulredc2(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y,
#                 const mp_limb_t *m, mp_limb_t inv_m);
#
# Linux:   z: %rdi, x: %rsi, y: %rdx, m: %rcx, inv_m: %r8
#          Needs %rbx, %rsp, %rbp, %r12-%r15 restored
# Windows: z: %rcx, x: %rdx, y: %r8,  m: %r9, inv_m: 28(%rsp)
#          Needs %rbx, %rbp, %rdi, %rsi, %r12...%15 restored

# This stuff is run through M4 twice, first when generating the
# mulredc*.asm files from the mulredc.m4 file (when preparing the distro)
# and again when generating the mulredc*.s files from the mulredc*.asm files
# when the user compiles the program.
# We used to substitute XP etc. by register names in the first pass,
# but now with switching between Linux and Windows ABI, we do it in
# the second pass instead when we know which ABI we have, as that 
# allows us to assign registers differently for the two ABIs.
# That means that the defines for XP etc., need to be quoted once to be 
# protected in the first M4 pass, so that they are processed and 
# occurrences of XP etc. happen only in the second pass.



include(`config.m4')

	TEXT
.p2align 6 # x86_64 L1 code cache line is 64 bytes long
	GLOBL GSYM_PREFIX`'mulredc2
	TYPE(GSYM_PREFIX`'mulredc`'2,`function')

# Implements multiplication and REDC for two input numbers of LENGTH words
ifdef(`WINDOWS64_ABI', `# Uses Windows ABI', `# Uses Linux ABI')
# tmp[0 ... len+1] = 0
# for (i = 0; i < len; i++)
#   {
#     t = x[i] * y[0]; /* Keep and reuse this product */
#     u = ((t + tmp[0]) * invm) % 2^64
#     tmp[0] += (t + m[0]*u) / 2^64; /* put carry in cy. */
#     for (j = 1; j < len; j++)
#       {
#         tmp[j-1 ... j] += x[i]*y[j] + m[j]*u + (cy << BITS_PER_WORD);
#         /* put new carry in cy */
#       }
#     tmp[len] = cy;
#   }
# z[0 ... len-1] = tmp[0 ... len-1]
# return (tmp[len])


# Values that are referenced only once in the loop over j go into r8 .. r14,
# In the inner loop (over j), tmp, x[i], y, m, and u are constant.
# tmp[j], tmp[j+1], tmp[j+2] are updated frequently. These 8 values
# stay in registers and are referenced as
# TP = tmp, YP = y, MP = m, 
# XI = x[i], T0 = tmp[j], T1 = tmp[j+1], CY = carry

define(`T0', `rsi')dnl
define(`T0l', `esi')dnl
define(`T1', `rbx')dnl
define(`T1l', `ebx')dnl
define(`CY', `rcx')dnl
define(`CYl', `ecx')dnl
define(`CYb', `cl')dnl
define(`XI', `r14')dnl		# register that holds x[i] value
define(`U', `r11')dnl
define(`XP', `r13')dnl		# register that points to the x arraz
define(`TP', `rbp')dnl		# register that points to t + i
define(`I', `r12')dnl		# register that holds loop counter i
define(`Il', `r12d')dnl		# register that holds loop counter i
define(`ZP', `rdi')dnl		# register that holds z. Same as passed in
ifdef(`WINDOWS64_ABI',
`define(`YP', `r8')dnl		# points to y array, same as passed in
define(`MP', `r9')dnl		# points to m array, same as passed in
define(`INVM', `r10')dnl	# register that holds invm. Same as passed in'
,
`define(`YP', `r9')dnl		# register that points to the y array
define(`MP', `r10')dnl		# register that points to the m array
define(`INVM', `r8')dnl		# register that holds invm. Same as passed in'
)dnl

`#' Register vars: `T0' = T0, `T1' = T1, `CY' = CY, `XI' = XI, `U' = U
`#'                `YP' = YP, `MP' = MP, `TP' = TP

# local variables: tmp[0 ... LENGTH] array, having LENGTH+1 8-byte words
# The tmp array needs LENGTH+1 entries, the last one is so that we can 
# store CY at tmp[j+1] for j == len-1



GSYM_PREFIX`'mulredc2:
	pushq	%rbx
	pushq	%rbp
	pushq	%r12
	pushq	%r13
	pushq	%r14
ifdef(`WINDOWS64_ABI',
`	pushq	%rsi
	pushq	%rdi
') dnl
ifdef(`WINDOWS64_ABI',
`	movq	%rdx, %XP
	movq	%rcx, %ZP
	movq	96(%rsp), %INVM # 7 push, ret addr, 4 reg vars = 96 bytes'
,
`	movq	%rsi, %XP		# store x in XP
	movq	%rdx, %YP		# store y in YP
	movq	%rcx, %MP		# store m in MP'
) dnl
	subq	$24, %rsp	# subtract size of local vars


#########################################################################
# i = 0 pass
#########################################################################

# register values at loop entry: %TP = tmp, %I = i, %YP = y, %MP = m
# %CY < 255 (i.e. only low byte may be != 0)

# Pass for j = 0. We need to fetch x[i] from memory and compute the new u

	movq	(%XP), %XI		# XI = x[0]
	movq	(%YP), %rax		# rax = y[0]

	xorl	%CYl, %CYl		# set %CY to 0
	lea	(%rsp), %TP		# store addr of tmp array in TP
	movl	%CYl, %Il		# Set %I to 0

	mulq	%XI			# rdx:rax = y[0] * x[i]
	addq	$1, %I

	movq 	%rax, %T0		# Move low word of product to T0
	movq	%rdx, %T1		# Move high word of product to T1

ifdef(`MULREDC_SVOBODA',
, `'
`	imulq	%INVM, %rax		# %rax = ((x[i]*y[0]+tmp[0])*invm)%2^64'
) 	movq	%rax, %U		# this is the new u value

	mulq	(%MP)			# multipy u*m[0]
	addq	%rax, %T0		# Now %T0 = 0, need not be stored
	movq	8(%YP), %rax		# Fetch y[1]
	adcq	%rdx, %T1		# 
	setc	%CYb
	# CY:T1:T0 <= 2*(2^64-1)^2 <= 2^2*128 - 4*2^64 + 2, hence
	# CY:T1 <= 2*2^64 - 4

define(`TT', defn(`T0'))dnl
define(`TTl', defn(`T0l'))dnl
define(`T0', defn(`T1'))dnl
define(`T0l', defn(`T1l'))dnl
define(`T1', defn(`TT'))dnl
define(`T1l', defn(`TTl'))dnl
undefine(`TT')dnl
undefine(`TTl')dnl
`#' Now `T0' = T0, `T1' = T1


`#' Pass for j = 1. Don't fetch new data from y[j+1].

	movl	%CYl, %T1l	# T1 = CY <= 1
	
	mulq	%XI		# y[j] * x[i]
	addq	%rax, %T0	# Add low word to T0
	movq	8(%MP), %rax	# Fetch m[j] into %rax
	adcq	%rdx, %T1 	# Add high word with carry to T1
	mulq    %U		# m[j]*u
	addq	%rax, %T0	# Add low word to T0
	movq	%T0, 0(%TP)	# Store T0 in tmp[j-1]
	adcq	%rdx, %T1	# Add high word with carry to T1
	movq	%T1, 8(%TP)	# Store T1 in tmp[j]
	setc	%CYb		# %CY <= 1
	movq	%CY, 16(%TP)	# Store CY in tmp[j+1]

#########################################################################
# i > 0 passes
#########################################################################

.p2align 5,,4
LABEL_SUFFIX(1)

# register values at loop entry: %TP = tmp, %I = i, %YP = y, %MP = m
# %CY < 255 (i.e. only low byte may be > 0)

# Pass for j = 0. We need to fetch x[i], tmp[i] and tmp[i+1] from memory
# and compute the new u

	movq	(%XP,%I,8), %XI		# XI = x[i]
	movq	(%YP), %rax		# rax = y[0]
#init the register tmp ring buffer
        movq	(%TP), %T0		# Load tmp[0] into T0
	movq	8(%TP), %T1		# Load tmp[1] into T1

	mulq	%XI			# rdx:rax = y[0] * x[i]
	addq	$1, %I

	addq	%T0, %rax		# Add T0 to low word
	adcq	%rdx, %T1		# Add high word with carry to T1
	setc	%CYb			# %CY <= 1

	movq 	%rax, %T0		# Save sum of low words in T0
	imulq	%INVM, %rax		# %rax = ((x[i]*y[0]+tmp[0])*invm)%2^64
	movq	%rax, %U		# this is the new u value

	mulq	(%MP)			# multipy u*m[0]
	addq	%rax, %T0		# Now %T0 = 0, need not be stored
	adcq	%rdx, %T1		# 

	movq	8(%YP), %rax		# Fetch y[1]

define(`TT', defn(`T0'))dnl
define(`TTl', defn(`T0l'))dnl
define(`T0', defn(`T1'))dnl
define(`T0l', defn(`T1l'))dnl
define(`T1', defn(`TT'))dnl
define(`T1l', defn(`TTl'))dnl
undefine(`TT')dnl
undefine(`TTl')dnl
`#' Now `T0' = T0, `T1' = T1


`#' Pass for j = 1. Don't fetch new data from y[j+1].

	movq	16(%TP), %T1
	adcq	%CY, %T1	# T1 = CY + tmp[j+1]
	
	mulq	%XI		# y[j] * x[i]
	addq	%rax, %T0	# Add low word to T0
	movq	8(%MP), %rax	# Fetch m[j] into %rax
	adcq	%rdx, %T1 	# Add high word with carry to T1
	setc	%CYb	    	# %CY <= 1
	mulq    %U		# m[j]*u
	addq	%rax, %T0	# Add low word to T0
	movq	%T0, 0(%TP)	# Store T0 in tmp[j-1]
	adcq	%rdx, %T1	# Add high word with carry to T1
	movq	%T1, 8(%TP)	# Store T1 in tmp[j]
	adcb	$0, %CYb	# %CY <= 2
	movq	%CY, 16(%TP)	# Store CY in tmp[j+1]

	cmpq	$2, %I
	jb	1b

# Copy result from tmp memory to z
	movq	(%TP), %rax
	movq	8(%TP), %rdx
	movq	%rax, (%ZP)
	movq	%rdx, 8(%ZP)

	movl	%CYl, %eax	# use carry as return value
	addq	$24, %rsp
ifdef(`WINDOWS64_ABI',
`	popq	%rdi
	popq	%rsi
') dnl
	popq	%r14
	popq	%r13
	popq	%r12
	popq	%rbp
	popq	%rbx
	ret

ifdef(`WINDOWS64_ABI',
,
`
`#'if defined(__linux__) && defined(__ELF__)
.section .note.GNU-stack,"",%progbits
`#'endif
') dnl
"""
else:
	print mulredc_k_rolled(k)

