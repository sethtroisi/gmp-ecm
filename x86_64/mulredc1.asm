#
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
	movq	%rdx, %rax	# %rax = y

	mulq	%rsi		# [rdx:rax] = x*y

	movq	%rdx, %r10
	movq	%rax, %r9	# store xy in [r9:r10]

	mulq	%r8		# compute u = (xy % 2^32) * inv_m

	mulq	%rcx		# compute u*m

	addq	%r9, %rax	# rax is 0 now (carry is important)

	adcq	%r10, %rdx

	movq	%rdx, (%rdi)
	adcq	$0, %rax	# return carry in %rax

	ret

