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

ifdef(`MINGW64_ABI',
# stack: inv_m, %r9: m, %r8: y, %rdx: x, %rcx: *z
define(`INV_M', `(%rsp)')
define(`M', `%r9')
define(`Y', `%r8')
define(`X', `%rdx')
define(`Z', `%rcx')
define(`TMP2', `%r10')
define(`TMP1', `%r11'),
# %r8: inv_m, %rcx: m, %rdx: y, %rsi : x, %rdi : *z
define(`INV_M', `%r8')
define(`M', `%rcx')
define(`Y', `%rdx')
define(`X', `%rsi')
define(`Z', `%rdi')
define(`TMP2', `%r10')
define(`TMP1', `%r9'))

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

