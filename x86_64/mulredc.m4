`# mp_limb_t mulredc'LENGTH`(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y,'
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

divert(-1)
# forloop(i, from, to, stmt)

define(`forloop', `pushdef(`$1', `$2')_forloop(`$1', `$2', `$3', `$4')popdef(`$1')')
define(`_forloop',
       `ifelse(eval($1 <= `$3'), 1, 
            `$4'`define(`$1', incr($1))_forloop(`$1', `$2', `$3', `$4')')')
divert

`include(`config.m4')'

	TEXT
.p2align 6 # x86_64 L1 code cache line is 64 bytes long
	GLOBL GSYM_PREFIX``''mulredc`'LENGTH
	TYPE(GSYM_PREFIX``''mulredc``''LENGTH,``function'')

# Implements multiplication and REDC for two input numbers of LENGTH words
`ifdef(`WINDOWS64_ABI', `# Uses Windows ABI', `# Uses Linux ABI')'
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

`define(`T0', `rsi')dnl
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
)dnl'

dnl Put overview of register allocation into .s file
``#'' `Register vars: `T0' = T0, `T1' = T1, `CY' = CY, `XI' = XI, `U' = U'
``#'' `               `YP' = YP, `MP' = MP, `TP' = TP'

# local variables: tmp[0 ... LENGTH] array, having LENGTH+1 8-byte words
# The tmp array needs LENGTH+1 entries, the last one is so that we can 
# store CY at tmp[j+1] for j == len-1


define(`LOCALSPACE', `eval(8*(LENGTH + 1))')dnl
define(`LOCALTMP', `(%rsp)')dnl

GSYM_PREFIX``''mulredc`'LENGTH:
	pushq	%rbx
	pushq	%rbp
	pushq	%r12
	pushq	%r13
	pushq	%r14
`ifdef(`WINDOWS64_ABI',
`	pushq	%rsi
	pushq	%rdi
') dnl'
`ifdef(`WINDOWS64_ABI',
`	movq	%rdx, %XP
	movq	%rcx, %ZP
	movq	96(%rsp), %INVM # 7 push, ret addr, 4 reg vars = 96 bytes'
,
`	movq	%rsi, %XP		# store x in XP
	movq	%rdx, %YP		# store y in YP
	movq	%rcx, %MP		# store m in MP'
) dnl'
	subq	$LOCALSPACE, %rsp	# subtract size of local vars


#########################################################################
# i = 0 pass
#########################################################################

# register values at loop entry: %TP = tmp, %I = i, %YP = y, %MP = m
# %CY < 255 (i.e. only low byte may be != 0)

# Pass for j = 0. We need to fetch x[i] from memory and compute the new u

	movq	(%XP), %XI		# XI = x[0]
	movq	(%YP), %rax		# rax = y[0]

	xorl	%CYl, %CYl		# set %CY to 0
	lea	LOCALTMP, %TP		# store addr of tmp array in TP
	movl	%CYl, %Il		# Set %I to 0

	mulq	%XI			# rdx:rax = y[0] * x[i]
	addq	$1, %I

	movq 	%rax, %T0		# Move low word of product to T0
	movq	%rdx, %T1		# Move high word of product to T1

`ifdef(`MULREDC_SVOBODA',
, `'
`	imulq	%INVM, %rax		# %rax = ((x[i]*y[0]+tmp[0])*invm)%2^64'
') dnl'
	movq	%rax, %U		# this is the new u value

	mulq	(%MP)			# multipy u*m[0]
	addq	%rax, %T0		# Now %T0 = 0, need not be stored
	movq	8(%YP), %rax		# Fetch y[1]
	adcq	%rdx, %T1		# 
	setc	%CYb
	# CY:T1:T0 <= 2*(2^64-1)^2 <= 2^2*128 - 4*2^64 + 2, hence
	# CY:T1 <= 2*2^64 - 4

dnl Cycle ring buffer. Only mappings of T0 and T1 to regs change, no MOVs!
`define(`TT', defn(`T0'))dnl
define(`TTl', defn(`T0l'))dnl
define(`T0', defn(`T1'))dnl
define(`T0l', defn(`T1l'))dnl
define(`T1', defn(`TT'))dnl
define(`T1l', defn(`TTl'))dnl
undefine(`TT')dnl
undefine(`TTl')dnl'
``#'' Now ``T0'' = T0, ``T1'' = T1

forloop(`UNROLL', 1, eval(LENGTH - 2), `dnl
define(`J', `eval(8 * UNROLL)')dnl
define(`J8', `eval(J + 8)')dnl
define(`JM8', `eval(J - 8)')dnl

``#'' Pass for j = UNROLL
``#'' Register values at entry: 
``#'' %rax = y[j], %XI = x[i], %U = u
``#'' %TP = tmp, %T0 = value to store in tmp[j], %T1 undefined 
``#'' %CY = carry into T1 (is <= 2)
# We have %CY:%T1 <= 2 * 2^64 - 2

	movl	%CYl, %T1l	# T1 = CY <= 1

	# Here, T1:T0 <= 2*2^64 - 2
	mulq	%XI		# y[j] * x[i]
	# rdx:rax <= (2^64-1)^2 <= 2^128 - 2*2^64 + 1
	addq	%rax, %T0	# Add low word to T0
	movq	J`'(%MP), %rax	# Fetch m[j] into %rax
	adcq	%rdx, %T1	# Add high word with carry to T1
	# T1:T0 <= 2^128 - 2*2^64 + 1 + 2*2^64 - 2 <= 2^128 - 1, no carry!
	
	mulq	%U		# m[j]*u
	# rdx:rax <= 2^128 - 2*2^64 + 1, T1:T0 <= 2^128 - 1
	addq	%T0, %rax	# Add T0 and low word
	movq	%rax, JM8`'(%TP)	``#'' Store T0 in tmp[UNROLL-1]
	movq	J8`'(%YP), %rax	``#'' Fetch y[j+1] = y[eval(UNROLL+1)] into %rax
	adcq	%rdx, %T1	# Add high word with carry to T1
	setc	%CYb		# %CY <= 1
	# CY:T1:T0 <= 2^128 - 1 + 2^128 - 2*2^64 + 1 <=
	#             2 * 2^128 - 2*2^64 ==> CY:T1 <= 2 * 2^64 - 2

dnl Cycle ring buffer. Only mappings of T0 and T1 to regs change, no MOVs!
`define(`TT', defn(`T0'))dnl
define(`TTl', defn(`T0l'))dnl
define(`T0', defn(`T1'))dnl
define(`T0l', defn(`T1l'))dnl
define(`T1', defn(`TT'))dnl
define(`T1l', defn(`TTl'))dnl
undefine(`TT')dnl
undefine(`TTl')dnl'
``#'' Now ``T0'' = T0, ``T1'' = T1

')dnl # end forloop

``#'' Pass for j = eval(LENGTH - 1). Don't fetch new data from y[j+1].
define(`J', `eval(8*LENGTH - 8)')dnl
define(`J8', `eval(J + 8)')dnl
define(`JM8', `eval(J - 8)')dnl

	movl	%CYl, %T1l	# T1 = CY <= 1
	
	mulq	%XI		# y[j] * x[i]
	addq	%rax, %T0	# Add low word to T0
	movq	J`'(%MP), %rax	# Fetch m[j] into %rax
	adcq	%rdx, %T1 	# Add high word with carry to T1
	mulq    %U		# m[j]*u
	addq	%rax, %T0	# Add low word to T0
	movq	%T0, JM8`'(%TP)	# Store T0 in tmp[j-1]
	adcq	%rdx, %T1	# Add high word with carry to T1
	movq	%T1, J`'(%TP)	# Store T1 in tmp[j]
	setc	%CYb		# %CY <= 1
	movq	%CY, J8`'(%TP)	# Store CY in tmp[j+1]

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

dnl Cycle ring buffer. Only mappings of T0 and T1 to regs change, no MOVs!
`define(`TT', defn(`T0'))dnl
define(`TTl', defn(`T0l'))dnl
define(`T0', defn(`T1'))dnl
define(`T0l', defn(`T1l'))dnl
define(`T1', defn(`TT'))dnl
define(`T1l', defn(`TTl'))dnl
undefine(`TT')dnl
undefine(`TTl')dnl'
``#'' Now ``T0'' = T0, ``T1'' = T1

forloop(`UNROLL', 1, eval(LENGTH - 2), `dnl
define(`J', `eval(8 * UNROLL)')dnl
define(`J8', `eval(J + 8)')dnl
define(`JM8', `eval(J - 8)')dnl

``#'' Pass for j = UNROLL
``#'' Register values at entry: 
``#'' %rax = y[j], %XI = x[i], %U = u
``#'' %TP = tmp, %T0 = value to store in tmp[j], %T1 value to store in 
``#'' tmp[j+1], %CY = carry into T1, carry flag: also carry into T1

	movq	J8`'(%TP), %T1
	adcq	%CY, %T1	# T1 = CY + tmp[j+1]
	setc	%CYb		# %CY <= 1

	mulq	%XI		# y[j] * x[i]
	addq	%rax, %T0	# Add low word to T0

	movq	%U, %rax
	adcq	%rdx, %T1	# Add high word with carry to T1
	adcb	$0, %CYb	# %CY <= 2
	
	mulq	J`'(%MP)	# m[j]*u
	addq	%rax, %T0	# Add T0 and low word

	movq	J8`'(%YP), %rax	``#'' Fetch y[j+1] = y[eval(UNROLL+1)] into %rax
	adcq	%rdx, %T1	# Add high word with carry to T1
	movq	%T0, JM8`'(%TP)	``#'' Store T0 in tmp[UNROLL-1]

dnl Cycle ring buffer. Only mappings of T0 and T1 to regs change, no MOVs!
`define(`TT', defn(`T0'))dnl
define(`TTl', defn(`T0l'))dnl
define(`T0', defn(`T1'))dnl
define(`T0l', defn(`T1l'))dnl
define(`T1', defn(`TT'))dnl
define(`T1l', defn(`TTl'))dnl
undefine(`TT')dnl
undefine(`TTl')dnl'
``#'' Now ``T0'' = T0, ``T1'' = T1

')dnl # end forloop

``#'' Pass for j = eval(LENGTH - 1). Don't fetch new data from y[j+1].
define(`J', `eval(8*LENGTH - 8)')dnl
define(`J8', `eval(J + 8)')dnl
define(`JM8', `eval(J - 8)')dnl

	movq	J8`'(%TP), %T1
	adcq	%CY, %T1	# T1 = CY + tmp[j+1]
	
	mulq	%XI		# y[j] * x[i]
	addq	%rax, %T0	# Add low word to T0
	movq	J`'(%MP), %rax	# Fetch m[j] into %rax
	adcq	%rdx, %T1 	# Add high word with carry to T1
	setc	%CYb	    	# %CY <= 1
	mulq    %U		# m[j]*u
	addq	%rax, %T0	# Add low word to T0
	movq	%T0, JM8`'(%TP)	# Store T0 in tmp[j-1]
	adcq	%rdx, %T1	# Add high word with carry to T1
	movq	%T1, J`'(%TP)	# Store T1 in tmp[j]
	adcb	$0, %CYb	# %CY <= 2
	movq	%CY, J8`'(%TP)	# Store CY in tmp[j+1]

	cmpq	$LENGTH, %I
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

	movl	%CYl, %eax	# use carry as return value
	addq	$LOCALSPACE, %rsp
`ifdef(`WINDOWS64_ABI',
`	popq	%rdi
	popq	%rsi
') dnl'
	popq	%r14
	popq	%r13
	popq	%r12
	popq	%rbp
	popq	%rbx
	ret

`ifdef(`WINDOWS64_ABI',
,
`
`#'if defined(__linux__) && defined(__ELF__)
.section .note.GNU-stack,"",%progbits
`#'endif
') dnl'
