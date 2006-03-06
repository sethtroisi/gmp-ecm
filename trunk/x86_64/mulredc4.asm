# mp_limb_t mulredc4(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y,
#                 const mp_limb_t *m, mp_limb_t inv_m);
#

include(`config.m4')
	TEXT
	GLOBL GSYM_PREFIX`'mulredc4
	TYPE(GSYM_PREFIX`'mulredc4,`function')

GSYM_PREFIX`'mulredc4:
	movq	%rdx, %r11
	movq	%rcx, %r10
	pushq	%rbx
	pushq	%rbp
	subq	$72, %rsp
#      %r8  : inv_m
#     %r10 : m
#     %r11 : y
#     %rsi : x
#     %rdi : z
#     %rsp : tmp
# Free registers
#     %rax, %rbx, %rcx, %rdx, %r9

### set tmp[0..2k+1[ to 0
	movq	$0, (%rsp)
	movq	$0, 8(%rsp)
	movq	$0, 16(%rsp)
	movq	$0, 24(%rsp)
	movq	$0, 32(%rsp)
	movq	$0, 40(%rsp)
	movq	$0, 48(%rsp)
	movq	$0, 56(%rsp)
	movq	$0, 64(%rsp)
###########################################
	movq	$4, %rbp

	.align 64
Loop:
	## compute u and store in %r9
	movq	(%rsi), %rax
	mulq	(%r11)
	addq	(%rsp), %rax
	mulq	%r8
	movq    %rax, %r9
### addmul1: src[0] is (%r10)
###          dst[0] is (%rsp)
###          mult is %r9
###          k is 4
###          kills %rax, %rbx, %rcx, %rdx
###   dst[0,k[ += mult*src[0,k[  plus carry put in rcx or rbx
	movq	(%r10), %rax
	mulq	%r9
	movq	%rax, %rbx
	movq	%rdx, %rcx

	movq	8(%r10), %rax
	mulq	%r9
	addq	%rbx, (%rsp)
	adcq	%rax, %rcx
	movq	%rdx, %rbx
	adcq	$0, %rbx

	movq	16(%r10), %rax
	mulq	%r9
	addq	%rcx, 8(%rsp)
	adcq	%rax, %rbx
	movq	%rdx, %rcx
	adcq	$0, %rcx

	movq	24(%r10), %rax
	mulq	%r9
	addq	%rbx, 16(%rsp)
	adcq	%rax, %rcx
	movq	%rdx, %rbx
	adcq	$0, %rbx
	addq	%rcx, 24(%rsp)
	adcq	$0, %rbx
### carry limb is in %rbx
	addq	%rbx, 32(%rsp)
	adcq	$0, 40(%rsp)
	movq	(%rsi), %r9
### addmul1: src[0] is (%r11)
###          dst[0] is (%rsp)
###          mult is %r9
###          k is 4
###          kills %rax, %rbx, %rcx, %rdx
###   dst[0,k[ += mult*src[0,k[  plus carry put in rcx or rbx
	movq	(%r11), %rax
	mulq	%r9
	movq	%rax, %rbx
	movq	%rdx, %rcx

	movq	8(%r11), %rax
	mulq	%r9
	addq	%rbx, (%rsp)
	adcq	%rax, %rcx
	movq	%rdx, %rbx
	adcq	$0, %rbx

	movq	16(%r11), %rax
	mulq	%r9
	addq	%rcx, 8(%rsp)
	adcq	%rax, %rbx
	movq	%rdx, %rcx
	adcq	$0, %rcx

	movq	24(%r11), %rax
	mulq	%r9
	addq	%rbx, 16(%rsp)
	adcq	%rax, %rcx
	movq	%rdx, %rbx
	adcq	$0, %rbx
	addq	%rcx, 24(%rsp)
	adcq	$0, %rbx
### carry limb is in %rbx
   addq    %rbx, 32(%rsp)
   adcq    $0, 40(%rsp)


	addq	$8, %rsi
	addq	$8, %rsp
	decq	%rbp
	jnz	Loop
###########################################
### Copy result in z
	movq	(%rsp), %rax
	movq	%rax, (%rdi)
	movq	8(%rsp), %rax
	movq	%rax, 8(%rdi)
	movq	16(%rsp), %rax
	movq	%rax, 16(%rdi)
	movq	24(%rsp), %rax
	movq	%rax, 24(%rdi)
	movq	32(%rsp), %rax	# carry
	addq    $40, %rsp
	popq	%rbp
	popq	%rbx
	ret

