# mp_limb_t mulredc15(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y,
#                 const mp_limb_t *m, mp_limb_t inv_m);
#
#  Stack:
#    inv_m    ## parameters
#    m
#    y
#    x
#    z							(4*(2k+7))%esp
#    ???   (1 limb???)
#    ebp      ## pushed registers                  (4*(2k+5))%esp
#    edi
#    esi
#    ebx
#    ...      ## counter (1 mp_limb_t)             (4*(2k+1))%esp
#    ...      ## tmp space (2*k+1 mp_limb_t)

include(`config.m4')
	TEXT
	GLOBL GSYM_PREFIX`'mulredc15
	TYPE(GSYM_PREFIX`'mulredc15,`function')

GSYM_PREFIX`'mulredc15:
	pushl	%ebp
	pushl	%edi
	pushl	%esi
	pushl	%ebx
	subl	$128, %esp
	movl	%esp, %edi
### set tmp[0..2k+1[ to 0
	movl	$0, (%edi)
	movl	$0, 4(%edi)
	movl	$0, 8(%edi)
	movl	$0, 12(%edi)
	movl	$0, 16(%edi)
	movl	$0, 20(%edi)
	movl	$0, 24(%edi)
	movl	$0, 28(%edi)
	movl	$0, 32(%edi)
	movl	$0, 36(%edi)
	movl	$0, 40(%edi)
	movl	$0, 44(%edi)
	movl	$0, 48(%edi)
	movl	$0, 52(%edi)
	movl	$0, 56(%edi)
	movl	$0, 60(%edi)
	movl	$0, 64(%edi)
	movl	$0, 68(%edi)
	movl	$0, 72(%edi)
	movl	$0, 76(%edi)
	movl	$0, 80(%edi)
	movl	$0, 84(%edi)
	movl	$0, 88(%edi)
	movl	$0, 92(%edi)
	movl	$0, 96(%edi)
	movl	$0, 100(%edi)
	movl	$0, 104(%edi)
	movl	$0, 108(%edi)
	movl	$0, 112(%edi)
	movl	$0, 116(%edi)
	movl	$0, 120(%edi)
###########################################
	movl	$15, 124(%esp)

	.align 32
Loop:
	## compute u and store in %ebp
	movl	152(%esp), %eax
	movl	156(%esp), %esi
	movl	(%eax), %eax
	mull	(%esi)
	addl	(%edi), %eax
	mull	164(%esp)
	movl    %eax, %ebp
	movl	160(%esp), %esi
### addmul1: src[0] is (%esi)
###          dst[0] is (%edi)
###          mult is %ebp
###          k is 15
###          kills %eax, %ebx, %ecx, %edx
###   dst[0,k[ += mult*src[0,k[  plus carry put in ecx or ebx
	movl	(%esi), %eax
	mull	%ebp
	movl	%eax, %ebx
	movl	%edx, %ecx
	movl	4(%esi), %eax

	mull	%ebp
	addl	%ebx, (%edi)
	movl	$0, %ebx
	adcl	%eax, %ecx
	movl	8(%esi), %eax
	adcl	%edx, %ebx

	mull	%ebp
	addl	%ecx, 4(%edi)
	movl	$0, %ecx
	adcl	%eax, %ebx
	movl	12(%esi), %eax
	adcl	%edx, %ecx

	mull	%ebp
	addl	%ebx, 8(%edi)
	movl	$0, %ebx
	adcl	%eax, %ecx
	movl	16(%esi), %eax
	adcl	%edx, %ebx

	mull	%ebp
	addl	%ecx, 12(%edi)
	movl	$0, %ecx
	adcl	%eax, %ebx
	movl	20(%esi), %eax
	adcl	%edx, %ecx

	mull	%ebp
	addl	%ebx, 16(%edi)
	movl	$0, %ebx
	adcl	%eax, %ecx
	movl	24(%esi), %eax
	adcl	%edx, %ebx

	mull	%ebp
	addl	%ecx, 20(%edi)
	movl	$0, %ecx
	adcl	%eax, %ebx
	movl	28(%esi), %eax
	adcl	%edx, %ecx

	mull	%ebp
	addl	%ebx, 24(%edi)
	movl	$0, %ebx
	adcl	%eax, %ecx
	movl	32(%esi), %eax
	adcl	%edx, %ebx

	mull	%ebp
	addl	%ecx, 28(%edi)
	movl	$0, %ecx
	adcl	%eax, %ebx
	movl	36(%esi), %eax
	adcl	%edx, %ecx

	mull	%ebp
	addl	%ebx, 32(%edi)
	movl	$0, %ebx
	adcl	%eax, %ecx
	movl	40(%esi), %eax
	adcl	%edx, %ebx

	mull	%ebp
	addl	%ecx, 36(%edi)
	movl	$0, %ecx
	adcl	%eax, %ebx
	movl	44(%esi), %eax
	adcl	%edx, %ecx

	mull	%ebp
	addl	%ebx, 40(%edi)
	movl	$0, %ebx
	adcl	%eax, %ecx
	movl	48(%esi), %eax
	adcl	%edx, %ebx

	mull	%ebp
	addl	%ecx, 44(%edi)
	movl	$0, %ecx
	adcl	%eax, %ebx
	movl	52(%esi), %eax
	adcl	%edx, %ecx

	mull	%ebp
	addl	%ebx, 48(%edi)
	movl	$0, %ebx
	adcl	%eax, %ecx
	movl	56(%esi), %eax
	adcl	%edx, %ebx
	mull	%ebp
	addl	%ecx, 52(%edi)
	adcl	%ebx, %eax
	adcl	$0, %edx
	addl	%eax, 56(%edi)
	adcl	$0, %edx
### carry limb is in %edx
	addl	%edx, 60(%edi)
	adcl	$0, 64(%edi)
	movl	152(%esp), %eax
	movl	(%eax), %ebp
	movl	156(%esp), %esi
### addmul1: src[0] is (%esi)
###          dst[0] is (%edi)
###          mult is %ebp
###          k is 15
###          kills %eax, %ebx, %ecx, %edx
###   dst[0,k[ += mult*src[0,k[  plus carry put in ecx or ebx
	movl	(%esi), %eax
	mull	%ebp
	movl	%eax, %ebx
	movl	%edx, %ecx
	movl	4(%esi), %eax

	mull	%ebp
	addl	%ebx, (%edi)
	movl	$0, %ebx
	adcl	%eax, %ecx
	movl	8(%esi), %eax
	adcl	%edx, %ebx

	mull	%ebp
	addl	%ecx, 4(%edi)
	movl	$0, %ecx
	adcl	%eax, %ebx
	movl	12(%esi), %eax
	adcl	%edx, %ecx

	mull	%ebp
	addl	%ebx, 8(%edi)
	movl	$0, %ebx
	adcl	%eax, %ecx
	movl	16(%esi), %eax
	adcl	%edx, %ebx

	mull	%ebp
	addl	%ecx, 12(%edi)
	movl	$0, %ecx
	adcl	%eax, %ebx
	movl	20(%esi), %eax
	adcl	%edx, %ecx

	mull	%ebp
	addl	%ebx, 16(%edi)
	movl	$0, %ebx
	adcl	%eax, %ecx
	movl	24(%esi), %eax
	adcl	%edx, %ebx

	mull	%ebp
	addl	%ecx, 20(%edi)
	movl	$0, %ecx
	adcl	%eax, %ebx
	movl	28(%esi), %eax
	adcl	%edx, %ecx

	mull	%ebp
	addl	%ebx, 24(%edi)
	movl	$0, %ebx
	adcl	%eax, %ecx
	movl	32(%esi), %eax
	adcl	%edx, %ebx

	mull	%ebp
	addl	%ecx, 28(%edi)
	movl	$0, %ecx
	adcl	%eax, %ebx
	movl	36(%esi), %eax
	adcl	%edx, %ecx

	mull	%ebp
	addl	%ebx, 32(%edi)
	movl	$0, %ebx
	adcl	%eax, %ecx
	movl	40(%esi), %eax
	adcl	%edx, %ebx

	mull	%ebp
	addl	%ecx, 36(%edi)
	movl	$0, %ecx
	adcl	%eax, %ebx
	movl	44(%esi), %eax
	adcl	%edx, %ecx

	mull	%ebp
	addl	%ebx, 40(%edi)
	movl	$0, %ebx
	adcl	%eax, %ecx
	movl	48(%esi), %eax
	adcl	%edx, %ebx

	mull	%ebp
	addl	%ecx, 44(%edi)
	movl	$0, %ecx
	adcl	%eax, %ebx
	movl	52(%esi), %eax
	adcl	%edx, %ecx

	mull	%ebp
	addl	%ebx, 48(%edi)
	movl	$0, %ebx
	adcl	%eax, %ecx
	movl	56(%esi), %eax
	adcl	%edx, %ebx
	mull	%ebp
	addl	%ecx, 52(%edi)
	adcl	%ebx, %eax
	adcl	$0, %edx
	addl	%eax, 56(%edi)
	adcl	$0, %edx
### carry limb is in %edx
   addl    %edx, 60(%edi)
   adcl    $0, 64(%edi)

	addl	$4, 152(%esp)
	addl	$4, %edi
	decl	124(%esp)
	jnz	Loop
###########################################
### Copy result in z
	movl	148(%esp), %ebx
	movl	(%edi), %eax
	movl	%eax, (%ebx)
	movl	4(%edi), %eax
	movl	%eax, 4(%ebx)
	movl	8(%edi), %eax
	movl	%eax, 8(%ebx)
	movl	12(%edi), %eax
	movl	%eax, 12(%ebx)
	movl	16(%edi), %eax
	movl	%eax, 16(%ebx)
	movl	20(%edi), %eax
	movl	%eax, 20(%ebx)
	movl	24(%edi), %eax
	movl	%eax, 24(%ebx)
	movl	28(%edi), %eax
	movl	%eax, 28(%ebx)
	movl	32(%edi), %eax
	movl	%eax, 32(%ebx)
	movl	36(%edi), %eax
	movl	%eax, 36(%ebx)
	movl	40(%edi), %eax
	movl	%eax, 40(%ebx)
	movl	44(%edi), %eax
	movl	%eax, 44(%ebx)
	movl	48(%edi), %eax
	movl	%eax, 48(%ebx)
	movl	52(%edi), %eax
	movl	%eax, 52(%ebx)
	movl	56(%edi), %eax
	movl	%eax, 56(%ebx)
	movl	60(%edi), %eax	# carry
	addl    $128, %esp
	popl	%ebx
	popl	%esi
	popl	%edi
	popl	%ebp
	ret

