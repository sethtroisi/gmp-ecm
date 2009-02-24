`# mp_limb_t mulredc1_'LENGTH`(mp_limb_t * z, const mp_limb_t x, const mp_limb_t * y,'
#                 const mp_limb_t *m, mp_limb_t inv_m);
#
# z: %rdi, x: %rsi, y: %rdx, m: %rcx, inv_m: %r8

divert(-1)
# forloop(i, from, to, stmt)

define(`forloop', `pushdef(`$1', `$2')_forloop(`$1', `$2', `$3', `$4')popdef(`$1')')
define(`_forloop',
       `ifelse(eval($1 <= `$3'), 1, 
            `$4'`define(`$1', incr($1))_forloop(`$1', `$2', `$3', `$4')')')
divert

`include(`config.m4')'

	TEXT
.align 64 # Opteron L1 code cache line is 64 bytes long
	GLOBL GSYM_PREFIX``''mulredc1_`'LENGTH
	TYPE(GSYM_PREFIX``''mulredc1_``''LENGTH,``function'')

/* Implements multiplication and REDC for one input numbers of LENGTH words
   and a multiplier of one word */

# Values that are referenced only once in the loop over j go into r8 .. r14,
# In the inner loop (over j), tmp, x[i], y, m, and u are constant.
# tmp[j], tmp[j+1], tmp[j+2] are updated frequently. These 8 values
# stay in registers and are referenced as
# YP = y, MP = m, 
# X = x, T0 = tmp[j], T1 = tmp[j+1], CY = carry

define(`T0', `rsi')dnl
define(`T1', `rbx')dnl
define(`CY', `rcx')dnl
define(`CYl', `ecx')dnl
define(`CYb', `cl')dnl
define(`X', `r14')dnl		# register that holds x[i] value
define(`U', `r11')dnl
define(`YP', `r9')dnl		# register that points to the y array
define(`MP', `r10')dnl		# register that points to the m array
define(`INVM', `r8')dnl		# register that holds invm. Same as passed in
define(`ZP', `rdi')dnl		# register that holds z. Same as passed in

dnl Put overview of register allocation into generated code
`#' Register vars: `T0' = T0, `T1' = T1, `CY' = CY, `X' = X, `U' = U
`#'                `YP' = YP, `MP' = MP

GSYM_PREFIX``''mulredc1_`'LENGTH:


#########################################################################
# i = 0 pass
#########################################################################

# register values at loop entry: %YP = y, %MP = m

# We need to compute u

	movq	(%rdx), %rax		# rax = y[0] (time critical, do first)
	pushq	%rbx
	pushq	%r14
	movq	%rdx, %YP		# store y in YP
	movq	%rcx, %MP		# store m in MP
	movq	%rsi, %X		# store x in X

	xorl	%CYl, %CYl		# set %CY to 0

	mulq	%X			# rdx:rax = y[0] * x

	movq 	%rax, %T0		# Move low word of product to T0
	movq	%rdx, %T1		# Move high word of product to T1

	imulq	%INVM, %rax		# %rax = ((x[i]*y[0]+tmp[0])*invm)%2^64
	movq	%rax, %U		# this is the new u value

	mulq	(%MP)			# multipy u*m[0]
	addq	%rax, %T0		# Now %T0 = 0, need not be stored
	movq	8(%YP), %rax		# Fetch y[1]
	adcq	%rdx, %T1		# 
	setc	%CYb
	# CY:T1:T0 <= 2*(2^64-1)^2 <= 2^2*128 - 4*2^64 + 2, hence
	# CY:T1 <= 2*2^64 - 4

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
`#' Now `T0' = T0, `T1' = T1

forloop(`UNROLL', 1, eval(LENGTH - 2), `dnl
define(`J', `eval(8 * UNROLL)')dnl
define(`J8', `eval(J + 8)')dnl
define(`JM8', `eval(J - 8)')dnl

`#' Pass for j = UNROLL
`#' Register values at entry: 
`#' %rax = y[j], %X = x, %U = u
`#' %T0 = value to store in tmp[j], %T1 undefined 
`#' %CY = carry into T1 (is <= 2)
# We have %CY:%T1 <= 2 * 2^64 - 2

	movq	%CY, %T1	# T1 = CY <= 1

	# Here, T1:T0 <= 2*2^64 - 2
	mulq	%X		# y[j] * x
	# rdx:rax <= (2^64-1)^2 <= 2^128 - 2*2^64 + 1
	addq	%rax, %T0	# Add low word to T0
	movq	J`'(%MP), %rax	# Fetch m[j] into %rax
	adcq	%rdx, %T1	# Add high word with carry to T1
	# T1:T0 <= 2^128 - 2*2^64 + 1 + 2*2^64 - 2 <= 2^128 - 1, no carry!
ifdef(`WANT_ASSERT', `
	jnc	1f
	call	abort
1:
',`')
	
	mulq	%U		# m[j]*u
	# rdx:rax <= 2^128 - 2*2^64 + 1, T1:T0 <= 2^128 - 1
	addq	%T0, %rax	# Add T0 and low word
	movq	%rax, JM8`'(%ZP)	`#' Store T0 in z[UNROLL-1]
	movq	J8`'(%YP), %rax	`#' Fetch y[j+1] = y[eval(UNROLL+1)] into %rax
	adcq	%rdx, %T1	# Add high word with carry to T1
	setc	%CYb		# %CY <= 1
	# CY:T1:T0 <= 2^128 - 1 + 2^128 - 2*2^64 + 1 <=
	#             2 * 2^128 - 2*2^64 ==> CY:T1 <= 2 * 2^64 - 2

dnl Cycle ring buffer. Only mappings of T0 and T1 to regs change, no MOVs!
define(`TT', defn(`T0'))dnl
define(`T0', defn(`T1'))dnl
define(`T1', defn(`TT'))dnl
undefine(`TT')dnl
`#' Now `T0' = T0, `T1' = T1

')dnl # end forloop

`#' Pass for j = eval(LENGTH - 1). Don't fetch new data from y[j+1].
define(`J', `eval(8*LENGTH - 8)')dnl
define(`J8', `eval(J + 8)')dnl
define(`JM8', `eval(J - 8)')dnl

	movq	%CY, %T1	# T1 = CY <= 1
	
	mulq	%X		# y[j] * x[i]
	addq	%rax, %T0	# Add low word to T0
	movq	J`'(%MP), %rax	# Fetch m[j] into %rax
	adcq	%rdx, %T1 	# Add high word with carry to T1
	mulq    %U		# m[j]*u
	addq	%rax, %T0	# Add low word to T0
	movq	%T0, JM8`'(%ZP)	# Store T0 in z[j-1]
	adcq	%rdx, %T1	# Add high word with carry to T1
	movq	%T1, J`'(%ZP)	# Store T1 in tmp[j]
	setc	%CYb		# %CY <= 1

	movq	%CY, %rax	# use carry as return value
	popq	%r14
	popq	%rbx
	ret
