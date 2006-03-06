#
# void ecm_redc3(mp_limb_t * z, const mp_limb_t * x, size_t n, mp_limb_t m)
#                     %rdi           %rsi              %rdx      %rcx
#
#  save them in       %r8            %r9               %r10      %r11


include(`config.m4')
        TEXT
        GLOBL GSYM_PREFIX`'ecm_redc3
        TYPE(GSYM_PREFIX`'ecm_redc3,`function')

GSYM_PREFIX`'ecm_redc3:
	push	%rbp					# Push registers
	push	%rbx
	subq	$32, %rsp				# SF: 2 Cpt + Jump +1

	movq    %rdi, %r8
	movq    %rsi, %r9
	movq    %rdx, %r10
	movq    %rcx, %r11

	movq	%r10, %rcx                          # Read size
	movq	%rcx, (%rsp)				# Save counter
        cmpq    $3, %rcx
        jae     Unroll		
Loop:	
		movq	%r11, %rax			# Read invm
	        movq    %r9, %rsi                  # Read Source Ptr
		mulq	(%rdi)				# Dest[0] * invm
		movq	%rdi, %r8			# Save new Dest
		movq	%r10, %rcx			# Read Size (2)
		xorq	%rbx, %rbx			# Initial Carry
		movq	%rax, %rbp			# Multiplier
InnerLoop:
		        # rsi:	  Source
		        # rdi:	  Dest
			# rbp:	  Multiplier
			# rcx:	  Counter
		        movq    (%rsi), %rax		# U1
			addq    $8, %rdi		# V1
			mulq    %rbp			# U2
			addq    $8, %rsi		# V2
			addq    %rbx, %rax		# U3
		        adcq    $0, %rdx		# U4
			addq    %rax, -8(%rdi)		# V4
			adcq    $0, %rdx		# U5
			decq    %rcx			# V5
			movq    %rdx, %rbx		# U6
			jnz     InnerLoop		# V6
		movq	%r8, %rdi
		movq    %rbx, (%rdi)                    # Save final carry
		decq	(%rsp)
		lea	8(%rdi), %rdi			# Advance Dest
		jnz     Loop				# Loop
End:
	addq	$32, %rsp
	pop	%rbx
	pop	%rbp
	ret
	
Unroll:
# %rcx Read size // %rdi Dest Ptr
	# Precalcul du saut.   21 bytes per (was 15 on x86)
	movq    %rcx, %rdx
        decq    %rcx
	subq    $2, %rdx
	negq    %rcx	
	shrq    $4, %rdx
	andq    $15, %rcx
	movq    %rdx, 16(%rsp)				# Org Cpt of 8(%rsp)
	movq    %rcx, %rdx
	shlq    $4, %rdx
        leaq    UnrollEntry (%rdx, %rcx,4), %rdx
	addq	%rcx, %rdx
	negq    %rcx
	movq	%rcx, %r10				# (-size)%16
	movq	%rdx, 24(%rsp)				# Org PC inside	

UnrollLoop:	
                movq    %r11, %rax                  # Read invm
                movq    %r9, %rsi                  # Read Source Ptr
                mulq    (%rdi)                          # Dest[0] * invm
                movq    %rdi, %r8                  # Save new Dest
                movq    %r10, %rcx                  # Read Size %16
		movq    16(%rsp), %rdx			# Read InnerLoop Cpt
	        movq    %rax, %rbp                      # Set Multiplier
		movq	%rdx, 8(%rsp)			# Set InnerLoop Cpt
	
		# First mull and set initial carry
	        movq    (%rsi), %rax
	        leaq    8(%rsi,%rcx,8), %rsi
	        mulq    %rbp
		leaq    (%rdi,%rcx,8), %rdi
	        movq    %rdx, %rbx
	
		# Do the Jump inside the unrolling loop
		# And set up the registers differently if odd
	        movq    24(%rsp), %rdx
	        testq   $1, %rcx
	        movq    %rax, %rcx
		cmovnz  %rbx, %rcx
	        cmovnz  %rax, %rbx	
	        jmp     *%rdx
		
		        # rax   scratch
			# rbx   carry hi
			# rcx   carry lo
			# rdx   scratch
			# rsi   src
			# rdi   dst
			# rbp   multiplier

	       .align  64, 0x90
UnrollInnerLoop:	
		addq    $128, %rdi   
UnrollEntry:	
#	        movq    0(%rsi), %rax # Can't use this instruction
	        .byte   0x48,0x8b,0x46,0x00
	        mulq    %rbp
#	        addq    %rcx, 0(%rdi) # Can't use this instruction
	        .byte   0x48,0x01,0x4f,0x00
	        adcq    %rax, %rbx
	        movq    %rdx, %rcx
	        adcq    $0, %rcx

	        movq    8(%rsi), %rax
	        mulq    %rbp
	        addq    %rbx, 8(%rdi)
	        adcq    %rax, %rcx
	        movq    %rdx, %rbx
	        adcq    $0, %rbx

	        movq    16(%rsi), %rax
	        mulq    %rbp
	        addq    %rcx, 16(%rdi)
	        adcq    %rax, %rbx
	        movq    %rdx, %rcx
	        adcq    $0, %rcx

	        movq    24(%rsi), %rax
	        mulq    %rbp
	        addq    %rbx, 24(%rdi)
	        adcq    %rax, %rcx
	        movq    %rdx, %rbx
	        adcq    $0, %rbx

	        movq    32(%rsi), %rax
	        mulq    %rbp
	        addq    %rcx, 32(%rdi)
	        adcq    %rax, %rbx
	        movq    %rdx, %rcx
	        adcq    $0, %rcx

	        movq    40(%rsi), %rax
	        mulq    %rbp
	        addq    %rbx, 40(%rdi)
	        adcq    %rax, %rcx
	        movq    %rdx, %rbx
	        adcq    $0, %rbx

	        movq    48(%rsi), %rax
	        mulq    %rbp
	        addq    %rcx, 48(%rdi)
	        adcq    %rax, %rbx
	        movq    %rdx, %rcx
	        adcq    $0, %rcx
	
	        movq    56(%rsi), %rax
	        mulq    %rbp
	        addq    %rbx, 56(%rdi)
	        adcq    %rax, %rcx
	        movq    %rdx, %rbx
	        adcq    $0, %rbx
		
	        movq    64(%rsi), %rax
	        mulq    %rbp
	        addq    %rcx, 64(%rdi)
	        adcq    %rax, %rbx
	        movq    %rdx, %rcx
	        adcq    $0, %rcx

	        movq    72(%rsi), %rax
	        mulq    %rbp
	        addq    %rbx, 72(%rdi)
	        adcq    %rax, %rcx
	        movq    %rdx, %rbx
	        adcq    $0, %rbx

	        movq    80(%rsi), %rax
	        mulq    %rbp
	        addq    %rcx, 80(%rdi)
	        adcq    %rax, %rbx
	        movq    %rdx, %rcx
	        adcq    $0, %rcx

	        movq    88(%rsi), %rax
	        mulq    %rbp
	        addq    %rbx, 88(%rdi)
	        adcq    %rax, %rcx
	        movq    %rdx, %rbx
	        adcq    $0, %rbx

	        movq    96(%rsi), %rax
	        mulq    %rbp
	        addq    %rcx, 96(%rdi)
	        adcq    %rax, %rbx
	        movq    %rdx, %rcx
	        adcq    $0, %rcx
	
	        movq    104(%rsi), %rax
	        mulq    %rbp
	        addq    %rbx, 104(%rdi)
	        adcq    %rax, %rcx
	        movq    %rdx, %rbx
	        adcq    $0, %rbx
		
	        movq    112(%rsi), %rax
	        mulq    %rbp
	        addq    %rcx, 112(%rdi)
	        adcq    %rax, %rbx
	        movq    %rdx, %rcx
	        adcq    $0, %rcx
	
		movq    120(%rsi), %rax
	        mulq    %rbp
	        addq    %rbx, 120(%rdi)
	        adcq    %rax, %rcx
	        movq    %rdx, %rbx
	        adcq    $0, %rbx
	
	        decq    8(%rsp)
	        leaq    128(%rsi), %rsi
	        jns     UnrollInnerLoop

	        addq    %rcx, 128(%rdi)
	        movq    %r8, %rdi
	        adcq    $0, %rbx
                movq    %rbx, (%rdi)                    # Save final carry
                decq    (%rsp)
                lea     8(%rdi), %rdi                   # Advance Dest
                jnz     UnrollLoop                      # Loop
End2:	
        addq    $32, %rsp
        pop     %rbx
        pop     %rbp
        ret
