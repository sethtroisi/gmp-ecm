`# mp_limb_t mulredc'LENGTH`(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y,'
#                 const mp_limb_t *m, mp_limb_t inv_m);
#
# z: %rdi, x: %rsi, y: %rdx, m: %rcx, inv_m: %r8

divert(-1)
include(`../config.m4')dnl
divert`'dnl
include(`forloop.m4')dnl

	TEXT
	GLOBL GSYM_PREFIX`'mulredc`'LENGTH
	TYPE(GSYM_PREFIX`'mulredc`'LENGTH,`function')

/* Implements multiplication and REDC for two input numbers of LENGTH words */

# tmp[0 ... 2*len] = 0
# for (i = 0; i < len; i++)
#   {
#     t = x[i] * y[0]; /* Keep and reuse this product */
#     u = ((t + tmp[i]) * invm) % 2^64
#     tmp[i ... i+1] += t + m[0]*u; /* put carry in cy. tmp[i] become 0 here */
#     for (j = 1; j < len; j++)
#       {
#         tmp[i+j ... i+j+1] += x[i]*y[j] + m[j]*u + (cy << BITS_PER_WORD);
#         /* put new carry in cy */
#       }
#     tmp[i+len+1] = cy;
#   }
# z[0 ... len-1] = tmp[len ... 2*len-1]
# return (tmp[2*len])

# Optimization hints: adcq %rxx, (%rxx), that is an add with carry to a
# memory location, has latency 4. The carry propagation in one long
# dependency chain. Therefore, adc to memory should be avoided - work on
# registers instead, except maybe for the last adcq in a carry chain.
# Athlon64 and Opteron decoders read aligned 16 byte packets.

# Values that are referenced only once in the loop over j go into r8 .. r14,

# In the inner loop (over j), tmp + i, x[i], y, m, and u are constant.
# tmp[i+j], tmp[i+j+1], tmp[i+j+2] are updated frequently. These 8 values
# stay in registers and are referenced as
# TP = tmp, YP = y, MP = m, 
# XI = x[i], T0 = tmp[i+j], T1 = tmp[i+j+1], CY = carry

define(`T0', `rsi')dnl
define(`T1', `rbx')dnl
define(`CY', `rcx')dnl
define(`CYb', `cl')dnl
define(`XI', `r14')dnl		# register that holds x[i] value
define(`U', `r11')dnl
define(`YP', `r9')dnl		# register that points to the y array
define(`MP', `r10')dnl		# register that points to the m array
define(`XP', `r13')dnl		# register that points to the x arraz
define(`TP', `rbp')dnl		# register that points to t + i
define(`I', `r12')dnl		# register that holds loop counter i
define(`INVM', `r8')		# register that holds invm. Same as passed in
define(`ZP', `rdi')		# register that holds z. Same as passed in

dnl Put overview of register allocation into generated code
`#' Register vars: `T0' = T0, `T1' = T1, `CY' = CY, `XI' = XI, `U' = U
`#'                `YP' = YP, `MP' = MP, `TP' = TP

# The tmp array need 2*LENGTH+1 entries, the last one is so that we can 
# store CY at tmp[i+j+2] for i == j == len-1

# local variables: tmp[0 ... LENGTH] array, having LENGTH+1 8-byte words

define(`LOCALSPACE', `eval(8*(LENGTH + 1))')dnl
define(`LOCALTMP', `(%rsp)')dnl

GSYM_PREFIX`'mulredc`'LENGTH:
	pushq	%rbx
	pushq	%rbp
	pushq	%r12
	pushq	%r13
	pushq	%r14
	subq	$LOCALSPACE, %rsp	# subtract size of local vars
	movq	%rsi, %XP		# store x in XP
	movq	%rdx, %YP		# store y in YP
	movq	%rcx, %MP		# store m in MP
	xorq	%rax, %rax		# set %rax to 0
	movq	%rax, %CY		# Set CY to 0

# Clear tmp memory
	lea	LOCALTMP, %TP		# store addr of tmp array in TP
	movq	%rax, %I 		# set I to 0
forloop(`UNROLL', 0, LENGTH, `dnl	# tmp[0 ... 2*len] = 0
ifelse(UNROLL, `0', dnl
`	movq	%rax, (%TP)', dnl
`	movq	%rax, eval(UNROLL * 8)(%TP)')
')

.align 32
1:

# register values at loop entry: %TP = tmp + i, %I = i, %YP = y, %MP = m

# Pass for j = 0. We need to fetch x[i], tmp[i] and tmp[i+1] from memory
# and compute the new u

	movq	(%XP,%I,8), %XI		# XI = x[i]
	movq	(%YP), %rax		# rax = y[0]
#init the register tmp ring buffer
        movq	(%TP), %T0		# Load tmp[i] into T0
	movq	8(%TP), %T1		# Load tmp[i+1] into T1

	mulq	%XI			# rdx:rax = y[0] * x[i]
	addq	$1, %I

	addq	%T0, %rax		# Add T0 to low word
	adcq	%rdx, %T1		# Add high word with carry to T1
	setc	%CYb			# %CY <= 1

	movq 	%rax, %T0		# Save sum of low words in T0
	imulq	%INVM, %rax		# %rax = ((x[i]*y[0]+tmp[i])*invm)%2^64
	movq	%rax, %U		# this is the new u value

	mulq	(%MP)			# multipy u*m[0]
	addq	%rax, %T0		# Now %T0 = 0, need not be stored
	adcq	%rdx, %T1		# 

	movq	8(%YP), %rax		# Fetch y[1]

ifdef(`WANT_ASSERT', `
        pushf
	testq	%T0, %T0
	jz	assert1
	call	abort
assert1:
	popf
',`')
dnl Cycle ring buffer. Only mappings of T0 and T1 to regs change, no MOVs!
define(`TT', defn(`T0'))dnl
define(`T0', defn(`T1'))dnl
define(`T1', defn(`TT'))dnl
undefine(`TT')dnl
undefine(`TTl')dnl
`#' Now `T0' = T0, `T1' = T1

forloop(`UNROLL', 1, eval(LENGTH - 2), `dnl
define(`J', `eval(8 * UNROLL)')dnl
define(`J8', `eval(J + 8)')dnl
define(`JM8', `eval(J - 8)')dnl

`#' Pass for j = UNROLL
`#' Register values at entry: 
`#' %rax = y[j], %XI = x[i], %U = u
`#' %TP = tmp + i, %T0 = value to store in tmp[i+j], %T1 value to store in 
`#' tmp[i+j+1], %CY = carry into T1, carry flag: also carry into T1

	movq	%CY, %T1	# T1 = CY
	adcq	J8`'(%TP), %T1	# T1 += tmp[j+1]
	setc	%CYb		# %CY <= 1

	mulq	%XI		# y[j] * x[i]
	addq	%rax, %T0	# Add low word to T0
	movq	J`'(%MP), %rax	`#' Fetch m[UNROLL] into %rax
	adcq	%rdx, %T1	# Add high word with carry to T1
	adcb	$0, %CYb	# %CY <= 2
	
	mulq	%U		`#' m[UNROLL]*u
	addq	%T0, %rax	# Add T0 to low word
	movq	%rax, JM8`'(%TP)	`#' Store T0 in tmp[UNROLL-1]
	adcq	%rdx, %T1	# Add high word with carry to T1
	movq	J8`'(%YP), %rax	`#' Fetch y[j+1] = y[eval(UNROLL+1)]

dnl Cycle ring buffer. Only mappings of T0 and T1 to regs change, no MOVs!
define(`TT', defn(`T0'))dnl
define(`T0', defn(`T1'))dnl
define(`T1', defn(`TT'))dnl
undefine(`TT')dnl
undefine(`TTl')dnl
`#' Now `T0' = T0, `T1' = T1

')dnl # end forloop

`#' Pass for j = eval(LENGTH - 1). Don't fetch new data from y[j+1].
define(`J', `eval(8*LENGTH - 8)')dnl
define(`J8', `eval(J + 8)')dnl
define(`JM8', `eval(J - 8)')dnl

	movq	%CY, %T1	# T1 = CY
	adcq	J8`'(%TP), %T1	# T1 += tmp[j + 1]
	setc	%CYb	    	# %CY <= 1
	
	mulq	%XI		# y[j] * x[i]
	addq	%rax, %T0	# Add low word to T0
	movq	J`'(%MP), %rax	`#' Fetch m[eval(LENGTH-1)] into %rax
	adcq	%rdx, %T1 	# Add high word with carry to T1
	adcb	$0, %CYb	# %CY <= 2
	mulq    %U		`#' m[eval(LENGTH-1)]*u
	addq	%rax, %T0	# Add low word to T0
	movq	%T0, JM8`'(%TP)	`#' Store `T0' in tmp[i+eval(LENGTH-1)-1]
	adcq	%rdx, %T1	# Add high word with carry to T1
	movq	%T1, J`'(%TP)	`#' Store `T1' in tmp[i+LENGTH]
	adcb	$0, %CYb	# %CY <= 3
	movq	%CY, J8`'(%TP)

	cmpq	$LENGTH, %I	# increase i by 1
	jb	1b

# Copy result from tmp memory to z
dnl ==== THIS LOOP WILL NOT WORK FOR LENGTH <= 1 ====
forloop(`UNROLL', 0, eval(LENGTH / 2 - 1), `dnl
define(`J', `eval(2 * UNROLL * 8)')dnl
define(`J8', `eval(J + 8)')dnl
ifelse(J, `0', dnl
`	movq	(%TP), %rax', dnl
`	movq	J`'(%TP), %rax')
	movq	J8`'(%TP), %rdx
ifelse(J, `0', dnl
`	movq	%rax, (%ZP)', dnl
`	movq	%rax, J`'(%ZP)')
	movq	%rdx, J8`'(%ZP)
')dnl
ifelse(eval(LENGTH % 2), 1, `dnl
define(`J', `eval(LENGTH * 8 - 8)')dnl
	movq	J`'(%TP), %rax
	movq	%rax, J`'(%ZP)
')dnl

	movq	%CY, %rax	# use carry as return value
	addq	$LOCALSPACE, %rsp
	popq	%r14
	popq	%r13
	popq	%r12
	popq	%rbp
	popq	%rbx
	ret
