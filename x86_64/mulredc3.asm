# mp_limb_t mulredc3(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y,
#                 const mp_limb_t *m, mp_limb_t inv_m);
#
# z: %rdi, x: %rsi, y: %rdx, m: %rcx, inv_m: %r8



include(`config.m4')

	TEXT
	GLOBL GSYM_PREFIX`'mulredc3
	TYPE(GSYM_PREFIX`'mulredc`'3,`function')

/* Implements multiplication and REDC for two input numbers of 3 words */

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

		# register that holds invm. Same as passed in
		# register that holds z. Same as passed in

# Register vars: T0 = rsi, T1 = rbx, CY = rcx, XI = r14, U = r11
#                YP = r9, MP = r10, TP = rbp

# The tmp array need 2*LENGTH+1 entries, the last one is so that we can 
# store CY at tmp[i+j+2] for i == j == len-1

# local variables: tmp[0 ... LENGTH] array, having LENGTH+1 8-byte words


GSYM_PREFIX`'mulredc3:
	pushq	%rbx
	pushq	%rbp
	pushq	%r12
	pushq	%r13
	pushq	%r14
	subq	$32, %rsp	# subtract size of local vars
	movq	%rsi, %r13		# store x in XP
	movq	%rdx, %r9		# store y in YP
	movq	%rcx, %r10		# store m in MP
	xorq	%rax, %rax		# set %rax to 0
	movq	%rax, %rcx		# Set CY to 0

# Clear tmp memory
	lea	(%rsp), %rbp		# store addr of tmp array in TP
	movq	%rax, %r12 		# set I to 0
	movq	%rax, (%rbp)
	movq	%rax, 8(%rbp)
	movq	%rax, 16(%rbp)
	movq	%rax, 24(%rbp)


.align 32
1:

# register values at loop entry: %TP = tmp + i, %I = i, %YP = y, %MP = m

# Pass for j = 0. We need to fetch x[i], tmp[i] and tmp[i+1] from memory
# and compute the new u

	movq	(%r13,%r12,8), %r14		# XI = x[i]
	movq	(%r9), %rax		# rax = y[0]
#init the register tmp ring buffer
        movq	(%rbp), %rsi		# Load tmp[i] into T0
	movq	8(%rbp), %rbx		# Load tmp[i+1] into T1

	mulq	%r14			# rdx:rax = y[0] * x[i]
	addq	$1, %r12

	addq	%rsi, %rax		# Add T0 to low word
	adcq	%rdx, %rbx		# Add high word with carry to T1
	setc	%cl			# %CY <= 1

	movq 	%rax, %rsi		# Save sum of low words in T0
	imulq	%r8, %rax		# %rax = ((x[i]*y[0]+tmp[i])*invm)%2^64
	movq	%rax, %r11		# this is the new u value

	mulq	(%r10)			# multipy u*m[0]
	addq	%rax, %rsi		# Now %T0 = 0, need not be stored
	adcq	%rdx, %rbx		# 

	movq	8(%r9), %rax		# Fetch y[1]


# Now T0 = rbx, T1 = rsi


# Pass for j = 1
# Register values at entry: 
# %rax = y[j], %r14 = x[i], %r11 = u
# %rbp = tmp + i, %rbx = value to store in tmp[i+j], %rsi value to store in 
# tmp[i+j+1], %rcx = carry into rsi, carry flag: also carry into rsi

	movq	%rcx, %rsi	# T1 = CY
	adcq	16(%rbp), %rsi	# T1 += tmp[j+1]
	setc	%cl		# %CY <= 1

	mulq	%r14		# y[j] * x[i]
	addq	%rax, %rbx	# Add low word to T0
	movq	8(%r10), %rax	# Fetch m[1] into %rax
	adcq	%rdx, %rsi	# Add high word with carry to T1
	adcb	$0, %cl	# %CY <= 2
	
	mulq	%r11		# m[1]*u
	addq	%rbx, %rax	# Add T0 to low word
	movq	%rax, 0(%rbp)	# Store rbx in tmp[1-1]
	adcq	%rdx, %rsi	# Add high word with carry to T1
	movq	16(%r9), %rax	# Fetch y[j+1] = y[2]

# Now T0 = rsi, T1 = rbx


# Pass for j = 2. Don't fetch new data from y[j+1].

	movq	%rcx, %rbx	# T1 = CY
	adcq	24(%rbp), %rbx	# T1 += tmp[j + 1]
	setc	%cl	    	# %CY <= 1
	
	mulq	%r14		# y[j] * x[i]
	addq	%rax, %rsi	# Add low word to T0
	movq	16(%r10), %rax	# Fetch m[2] into %rax
	adcq	%rdx, %rbx 	# Add high word with carry to T1
	adcb	$0, %cl	# %CY <= 2
	mulq    %r11		# m[2]*u
	addq	%rax, %rsi	# Add low word to T0
	movq	%rsi, 8(%rbp)	# Store T0 in tmp[i+2-1]
	adcq	%rdx, %rbx	# Add high word with carry to T1
	movq	%rbx, 16(%rbp)	# Store T1 in tmp[i+3]
	adcb	$0, %cl	# %CY <= 3
	movq	%rcx, 24(%rbp)

	cmpq	$3, %r12	# increase i by 1
	jb	1b

# Copy result from tmp memory to z
	movq	(%rbp), %rax
	movq	8(%rbp), %rdx
	movq	%rax, (%rdi)
	movq	%rdx, 8(%rdi)
	movq	16(%rbp), %rax
	movq	%rax, 16(%rdi)

	movq	%rcx, %rax	# use carry as return value
	addq	$32, %rsp
	popq	%r14
	popq	%r13
	popq	%r12
	popq	%rbp
	popq	%rbx
	ret
