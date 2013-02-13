# mp_limb_t mulredc7(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y,
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
	GLOBL GSYM_PREFIX`'mulredc7
	TYPE(GSYM_PREFIX`'mulredc7,`function')

GSYM_PREFIX`'mulredc7:
	pushl	%ebp
	pushl	%edi
	pushl	%esi
	pushl	%ebx
	subl	$64, %esp
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
###########################################
	movl	$7, 60(%esp)

	.align 32
Loop:
	## compute u and store in %ebp
	movl	88(%esp), %eax
	movl	92(%esp), %esi
	movl	(%eax), %eax
	mull	(%esi)
	addl	(%edi), %eax
	mull	100(%esp)
	movl    %eax, %ebp
	movl	96(%esp), %esi
### addmul1: src[0] is (%esi)
###          dst[0] is (%edi)
###          mult is %ebp
###          k is 7
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
	adcl	%ebx, %eax
	adcl	$0, %edx
	addl	%eax, 24(%edi)
	adcl	$0, %edx
### carry limb is in %edx
	addl	%edx, 28(%edi)
	adcl	$0, 32(%edi)
	movl	88(%esp), %eax
	movl	(%eax), %ebp
	movl	92(%esp), %esi
### addmul1: src[0] is (%esi)
###          dst[0] is (%edi)
###          mult is %ebp
###          k is 7
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
	adcl	%ebx, %eax
	adcl	$0, %edx
	addl	%eax, 24(%edi)
	adcl	$0, %edx
### carry limb is in %edx
   addl    %edx, 28(%edi)
   adcl    $0, 32(%edi)

	addl	$4, 88(%esp)
	addl	$4, %edi
	decl	60(%esp)
	jnz	Loop
###########################################
### Copy result in z
	movl	84(%esp), %ebx
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
	movl	28(%edi), %eax	# carry
	addl    $64, %esp
	popl	%ebx
	popl	%esi
	popl	%edi
	popl	%ebp
	ret

