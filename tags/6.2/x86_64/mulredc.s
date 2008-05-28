#
#  mp_limb_t mulredc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y,
#                 const mp_limb_t *m, mp_size_t n, mp_limb_t inv_m)
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


	.text
	.globl mulredc
	.type   mulredc,@function
mulredc:
	pushq	%rbx
	pushq	%r12				
	pushq	%r13				
	pushq	%r14				
	pushq	%r15				# push registers of caller
	movq	%r8, %rax
	shlq	$4, %rax
	addq	$8, %rax
	subq	%rax, %rsp			# allocate TMP on stack
	movq	%rsp, %r10			# r10 <- TMP
	pushq	%r10				# pointer to TMP
	pushq	%rdi				# push parameters
	pushq	%rsi
	pushq	%rdx
	pushq	%rcx

#  Stack:			Other parameters:
#    (%rsp) : m				%r8 : n
#   8(%rsp) : y				%r9 : inv_m
#  16(%rsp) : x
#  24(%rsp) : z
#  32(%rsp)  : TMP  (pointer to)
#  40+%rsp   : TMP  (first limb)

### Init: put 0 in TMP[0..2n]        
#############################
	movq	%r8, %rcx
	shlq	$1, %rcx
	incq	%rcx
	movq	32(%rsp), %rdi
InitLoop:
	movq	$0, (%rdi)
	addq	$8, %rdi
	decq	%rcx
	jnz	InitLoop

### Outer:
##############################
###   for (i = 0; i < n; ++i) { 
###     u = (TMP[0]+x[i]*y[0])*inv_m   mod 2^GMP_LIMB_BITS;
###     TMP[0..n+1] += x[i]*y + u*m;   // this is "Inner"
###     TMP++;
###   }
##############################
	movq	%r8, %r15
OuterLoop:
	## compute u
	movq	16(%rsp), %rbx
	movq	(%rbx), %rax
	movq	8(%rsp), %rbx
	mulq	(%rbx)
	movq	32(%rsp), %rbx
	addq	(%rbx), %rax
	mulq	%r9
	movq	%rax, %r10		# %r10 gets the multiplier u
	
	### Inner:
	######################
	### TMP[0..n+1] += x[i]*y + u*m
	###
	###  ax: ...		r8 : n (for outer)
	###  bx: m		r9 : inv_m (for outer)
	###  cx: inner cpt	r10: u
	###  dx: ...		r11: x[i]
	###  di: TMP		r12, r13, r14: ...
	###  si: y		r15: outer cpt
	######################
	movq	(%rsp), %rbx
	movq	32(%rsp), %rdi
	movq	8(%rsp), %rsi
	movq	16(%rsp), %rdx
	movq	(%rdx), %r11
	movq	%r8, %rcx
	xorq	%r12, %r12
	xorq	%r13, %r13
	InnerLoop:    ### r12: carry lo, r13: carry hi (can be at most 2)
		movq	(%rsi), %rax
		mulq	%r11
		addq	%rax, %r12
		adcq	$0, %rdx
		addq	%r12, (%rdi)
		adcq	$0, %rdx     ## carry flag is clean
		movq	%rdx, %r14
		movq    (%rbx), %rax
		mulq	%r10
		movq	%r13, %r12   ## carry hi becomes carry low
		xorq	%r13, %r13
		addq	%rax, (%rdi)
		adcq	%rdx, %r12
		adcq	$0, %r13
		addq	%r14, %r12
		adcq	$0, %r13
		
		addq	$8, %rdi
		addq	$8, %rsi
		addq	$8, %rbx
		decq	%rcx
		jnz InnerLoop
	addq	%r12, (%rdi)
	adcq	%r13, 8(%rdi)
	######################

	## advance TMP and x
	movq	32(%rsp), %rdx
	addq	$8, %rdx
	movq	%rdx, 32(%rsp)
	movq	16(%rsp), %rdx
	addq	$8, %rdx
	movq	%rdx, 16(%rsp)
	decq	%r15
	jnz	OuterLoop

### Finish:
##############################
### Copy TMP into z
##############################
	movq	%r8, %rcx
	movq	32(%rsp), %rsi
	movq	24(%rsp), %rdi
FinishLoop:
	movq	(%rsi), %rax
	movq	%rax, (%rdi)
	addq	$8, %rdi
	addq	$8, %rsi
	decq	%rcx
	jnz	FinishLoop


	addq	$40, %rsp	# clean parameters
	movq    %r8, %rax
	shlq    $4, %rax
	addq    $8, %rax 
	addq    %rax, %rsp      # free TMP on stack

	popq	%r15
	popq	%r14
	popq	%r13
	popq	%r12
	popq	%rbx

	movq	(%rsi), %rax   ## returned value
	ret
