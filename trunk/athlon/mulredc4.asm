# mp_limb_t mulredc4(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y,
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
	GLOBL GSYM_PREFIX`'mulredc4
	TYPE(GSYM_PREFIX`'mulredc4,`function')

GSYM_PREFIX`'mulredc4:
	pushl	%ebp
	pushl	%edi
	pushl	%esi
	pushl	%ebx
	subl	$40, %esp
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
###########################################
	movl	$4, 36(%esp)

	.align 32
Loop:
	## compute u and store in %ebp
	movl	64(%esp), %eax
	movl	68(%esp), %esi
	movl	(%eax), %eax
	mull	(%esi)
	addl	(%edi), %eax
	mull	76(%esp)
	movl    %eax, %ebp
	movl	72(%esp), %esi
### addmul1: src[0] is (%esi)
###          dst[0] is (%edi)
###          mult is %ebp
###          k is 4
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
	adcl	%ecx, %eax
	adcl	$0, %edx
	addl	%eax, 12(%edi)
	adcl	$0, %edx
### carry limb is in %edx
	addl	%edx, 16(%edi)
	adcl	$0, 20(%edi)
	movl	64(%esp), %eax
	movl	(%eax), %ebp
	movl	68(%esp), %esi
### addmul1: src[0] is (%esi)
###          dst[0] is (%edi)
###          mult is %ebp
###          k is 4
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
	adcl	%ecx, %eax
	adcl	$0, %edx
	addl	%eax, 12(%edi)
	adcl	$0, %edx
### carry limb is in %edx
   addl    %edx, 16(%edi)
   adcl    $0, 20(%edi)

	addl	$4, 64(%esp)
	addl	$4, %edi
	decl	36(%esp)
	jnz	Loop
###########################################
### Copy result in z
	movl	60(%esp), %ebx
	movl	(%edi), %eax
	movl	%eax, (%ebx)
	movl	4(%edi), %eax
	movl	%eax, 4(%ebx)
	movl	8(%edi), %eax
	movl	%eax, 8(%ebx)
	movl	12(%edi), %eax
	movl	%eax, 12(%ebx)
	movl	16(%edi), %eax	# carry
	addl    $40, %esp
	popl	%ebx
	popl	%esi
	popl	%edi
	popl	%ebp
	ret

