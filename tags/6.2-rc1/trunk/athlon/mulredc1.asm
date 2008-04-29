#
#  mp_limb_t mulredc1(mp_limb_t *z, const mp_limb_t x, const mp_limb_t y,
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
# Stack:
#    inv_m  20(%esp)
#    m      16
#    y      12(%esp)
#    x      8
#    z      4(%esp)

	movl	12(%esp), %eax
	mull	8(%esp)
	movl	%edx, 12(%esp)
	movl	%eax, 8(%esp)   # store xy in [8(%esp):12(%esp)]
	mull	20(%esp)          # compute u
	mull	16(%esp)         # compute u*m
	addl	8(%esp), %eax       # eax is 0, now (carry is important)
	adcl	12(%esp), %edx
	movl	4(%esp), %ecx
	movl    %edx, (%ecx)
	adcl	$0, %eax
	ret

