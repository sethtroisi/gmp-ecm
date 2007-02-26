# mp_limb_t mulredc2(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y,
#                 const mp_limb_t *m, mp_limb_t inv_m);
#

include(`config.m4')
	TEXT
	GLOBL GSYM_PREFIX`'mulredc2
	TYPE(GSYM_PREFIX`'mulredc2,`function')

GSYM_PREFIX`'mulredc2:
	movq	%rdx, %r11
	movq	%rcx, %r10
	pushq	%rbx
	pushq	%r12
	pushq	%rdi
	subq	$40, %rsp
#     %r8  : inv_m
#     %r10 : m
#     %r11 : y
#     %rsi : x
#     %rdi : x[i]
#     %rsp : tmp
# Free registers
#     %rax, %rbx, %rcx, %rdx, %r9

# Optimization hints: adcq %rxx, (%rxx), that is an add with carry to a
# memory location, has latency 4. The carry propagation in one long
# dependency chain. Therefore, adc to memory should be avoided - work on
# registers instead, except maybe for the last adcq in a carry chain.

### set tmp[0..2k+1[ to 0
#	movq	$0, (%rsp)
#	movq	$0, 8(%rsp)
#	movq	$0, 16(%rsp)
#	movq	$0, 24(%rsp)
#	movq	$0, 32(%rsp)
###########################################

# Athlon64 and Opteron decoders read aligned 16 byte packets

# First pass
	## compute u=((x[i]*y[0]+tmp[i])*invm)%2^64 and store in %r9
	movq	(%r11), %rax	# 1
	movq	(%rsi), %rdi	# 1
	xorq	%rbx, %rbx

	addq	$8, %rsi
	mulq	%rdi		# 2

	movq	%rax, %r12	# %r12 to be added to (%rsp)
	movq	%rdx, %rcx	# %rcx to be added to 8(%rsp)

	mulq	%r8		# 3

	movq    %rax, %r9	# %r9 = u
# Compute x*y + m*u (%r11[0,8]*%rsi[0] + %r10[0,8] * %r9) and 
# add to tmp[0,8,16,24]
	mulq	(%r10)		# 4 %rax still contains u

	addq	%rax, %r12	# %r12 to be moved to (%rsp), %rax free
	movq	8(%r10), %rax	# Fetch data for next mul

#	movq	%r12, (%rsp)	# Is always zero! %r12 free
	adcq	%rdx, %rcx	# %rcx to be added to 8(%rsp), %rdx free

	setc	%bl		# %bl <= 1 is the carry into 16(%rsp)
	addq	$8, %rsp

	mulq	%r9		# get the next mul on the way

	addq	%rax, %rcx	# %rcx to be added to 8(%rsp), %rax free
	movq	8(%r11), %rax	# Fetch data for next mul

	adcq	%rbx, %rdx	# propagate carry for 16(%rsp)

	movq	%rdx, %r12	# %r12 to be added to 16(%rsp), %rdx free
	setc	%bl		# %ebx <= 1 is carry into 24(%rsp)

	mulq	%rdi
	movq	(%rsi), %rdi

	addq	%rax, %rcx	# %rcx to be added to 8(%rsp), %rax free
	movq	(%r11), %rax

	movq	%rcx, (%rsp)	# %rcx free
	adcq    %rdx, %r12	# %r12 to be added to 16(%rsp), %rdx free

	adcl	$0, %ebx	# %ebx <= 2 is carry into 24(%rsp)

	movq	%r12, 8(%rsp)	# %r12 free   These two memory transfers
	movq	%rbx, 16(%rsp)	# %rbx free   are paired


# Last pass
	## compute u=((x[i]*y[0]+tmp[i])*invm)%2^64 and store in %r9

	mulq	%rdi
	movq	(%r10), %rbx	# Preload m

	movq	%rax, %r12	# %r12 to be added to (%rsp)
	movq	%rdx, %rcx	# %rcx to be added to 8(%rsp)
	addq	(%rsp), %rax

	mulq	%r8

	movq    %rax, %r9
# Compute m*u (%r10[0,8] * %r9) and add to tmp[0,8,16,24]
	movq	%rbx, %rax

	mulq	%r9

	addq	%rax, %r12	# %r12 to be added to (%rsp), %rax free
	movq	8(%r10), %rax	# Fetch data for next mul

	adcq	%rdx, %rcx	# %rcx to be added to 8(%rsp), %rdx free

	setb	%bl		# %bl <= 1 is the carry into 16(%rsp)

	mulq	%r9		# get the next mul on the way
	movzx	%bl, %rbx	# clear bytes 1-7 of %rbx
	addq	%r12, (%rsp)	# %r12 free

	adcq	%rax, %rcx	# %rcx to be added to 8(%rsp), %rax free
	movq	8(%r11), %rax	# Fetch data for next mul

	adcq	%rbx, %rdx	# %rbx <= 2 is the carry into 16(%rsp)

	movq	%rdx, %r12	# %r12 to be added to 16(%rsp), %rdx free
	setc	%bl		# %ebx <= 1 is carry into 24(%rsp)

	mulq	%rdi

	addq	%rax, %rcx	# %rcx to be added to 8(%rsp), %rax free

	adcq    %rdx, %r12	# %r12 to be added to 16(%rsp), %rdx free

	adcl	$0, %ebx	# %ebx <= 2 is carry into 24(%rsp)

	addq	%rcx, 8(%rsp)	# %rcx free

	adcq	%r12, 16(%rsp)	# %r12 free

	adcl	$0, %ebx	# %ebx <= 3 is carry into 24(%rsp)

###########################################
### Copy result in z
	movq	8(%rsp), %rax
	movq	16(%rsp), %rdx
	addq    $32, %rsp
	popq	%rdi
	movq	%rax, (%rdi)
	movq	%rdx, 8(%rdi)
	movq	%rbx, %rax	# carry
	popq	%r12
	popq	%rbx
	ret

# Used per outer loop (over i): invm
# Used per inner loop (over j): u, x[i], y[j], tmp[j], m[j]
# In reg.: x[i] (%rdi), u (%r9), *y (%r11), y[j] (%rax), y[j+1] (%r12), 
#          *m (%r10), m[j] (%r13), m[j+1] (%r14), 
#          *tmp (%rbp), tmp[j] (%rsi), tmp[j+1] (%rbx), tmp[j+2] (%rcx)
# Compute x[i]*y[j] + u*m[j], add to tmp[j]:tmp[j+1]:tmp[j+2]
# Compute x[i]*y[j+1] + u*m[j+1], add to tmp[j+1]:tmp[j+2]:tmp[j+3]
# j+=2
	mulq	%rdi		# x[i]*y[j]
	movq	j(%r10), %r13
	movq	j+8(%r10), %(r14)
	
	addq	%rax, %rsi
	movq	%rdx, %rcx
	movq	%r9, %rax

	mulq	%r13		# u*m[j]
