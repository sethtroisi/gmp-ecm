
; mp_limb_t mulredc<limbs>(       MSVC
;       mp_limb_t       *z,        rcx
;       const mp_limb_t *x,        rdx
;       const mp_limb_t *y,         r8
;       const mp_limb_t *m,         r9
;       mp_limb_t inv_m     [rsp+0x28]
;   )

%macro  mseq_1 4
	mov	    %2, rcx
	mul	    r14
	add	    %1, rax
	mov	    rax, [r9+8*%3]
	adc	    %2, rdx
	mul	    r11
%if %3 < %4 - 1
	    add	    rax, %1
	    mov     [rbp+8*(%3-1)], rax
	    mov	    rax, [r8+8*(%3+1)]
	    adc	    %2, rdx
	    setc	cl
%else
        add     %1, rax
	    mov	    [rbp+8*(%3-1)], %1
	    adc     %2, rdx
	    mov     [rbp+8*%3], %2
	    setc	cl
	    mov	    [rbp+8*(%3+1)], rcx
%endif
%endmacro

%macro mseq_20 2
	mov	    r14, [r13+r12*8]
	mov	    rax, [r8]
    mov	    %1, [rbp]
	mov	    %2, [rbp+8]
	mul	    r14
	add	    r12, 1
	add	    rax, %1
	adc	    %2, rdx
	setc	cl
	mov 	%1, rax
	imul	rax, r10
	mov	    r11, rax
	mul	    qword [r9]
	add	    %1, rax
	adc	    %2, rdx
	mov	    rax, [r8+8]
%endmacro

%macro mseq_2 4
	mov     %2, rcx
	adc	    %2, [rbp+8*(%3+1)]
	setc	cl
	mul	    r14
	add	    %1, rax
	mov	    rax, [r9+8*%3]
	adc	    %2, rdx
	adc	    cl, 0
	mul	    r11
%if %3 < %4 - 1
	    add	    rax, %1
	    mov	    [rbp+8*(%3-1)], rax
	    adc	    %2, rdx
	    mov	    rax, [r8+8*(%3+1)]
%else
        add     %1, rax
        mov     [rbp+8*(%3-1)], %1
        adc     %2, rdx
	    mov	    [rbp+8*%3],%2
	    adc	    cl, 0
	    mov	    [rbp+8*(%3+1)], rcx
%endif
%endmacro

%macro store 1
%assign i 0
%rep %1
    %if i == %1 - 1 && (%1 & 1)
	    mov	    rax, [rbp + 8 * i]
	    mov 	[rdi + 8 * i], rax
    %elif (i & 1)
	    mov 	[rdi + 8 * (i - 1)], rax
	    mov 	[rdi + 8 * i], rdx
    %else
	    mov	    rax, [rbp + 8 * i]
    	mov	    rdx, [rbp + 8 * (i + 1)]
    %endif
    %assign i i + 1
%endrep
%endmacro

%macro mulredc 1

%assign limbs       %1
%define f_name(x)   mulredc %+ x
%define stack_space 8 * (limbs + 2 - (limbs & 1))

	global	f_name(limbs)
%ifdef DLL
	export	f_name(limbs)
%endif

    align   64

PROC_FRAME	f_name(limbs)
    push_reg    rbp
    push_reg    rbx
    push_reg    rsi
    push_reg    rdi
    push_reg    r12
    push_reg    r13
    push_reg    r14
    alloc_stack stack_space
END_PROLOGUE
                                        ;   *y ->  r8
	mov     rdi, rcx                    ;   *z -> rdi
	mov	    r13, rdx                    ;   *x -> r13
    mov     r10, [rsp+0x28+stack_space] ; invm -> r10
                                        ;   *m in  r9    
	mov     r14, [r13]
	mov	    rax, [r8]
	xor     rcx, rcx
	lea	    rbp, [rsp]
	mov     r12, rcx
	mul     qword r14
	add	    r12, 1
	mov 	rsi, rax
	mov	    rbx, rdx
	imul	rax, r10
	mov	    r11, rax
	mul	    qword [r9]
	add	    rsi, rax
	mov	    rax, [r8+8]
	adc	    rbx, rdx
	setc	cl

%assign j 1
%rep    limbs - 1
%if (j & 1)
    mseq_1  rbx, rsi, j, limbs
%else
    mseq_1  rsi, rbx, j, limbs
%endif
   %assign j j + 1
%endrep

    align 32
.1:

%assign j 1
%if (limbs & 1)
    mseq_20 rsi, rbx
    %rep    limbs - 1
        %if (j & 1)
            mseq_2  rbx, rsi, j, limbs
        %else
            mseq_2  rsi, rbx, j, limbs
        %endif
        %assign j j + 1
    %endrep
%else
    mseq_20 rbx, rsi
    %rep    limbs - 1
        %if (j & 1)
            mseq_2  rsi, rbx, j, limbs
        %else
            mseq_2  rbx, rsi, j, limbs
        %endif
        %assign j j + 1
    %endrep
%endif

	cmp     r12, limbs
	jb	    .1

    store   limbs

	mov	    rax, rcx
    add     rsp, stack_space
    pop     r14
    pop     r13
    pop     r12
    pop     rdi
    pop     rsi
    pop     rbx
    pop     rbp
    ret
ENDPROC_FRAME
%endmacro

	bits    64
	section .text
	
	global	mulredc1
%ifdef DLL
	export	mulredc1
%endif

    align   64
mulredc1:
	mov	    rax, r10
	mul	    rdx
	mov	    r9, rax
	mov	    r11, rdx
	mul	    qword [rsp+0x28]
	mul	    r8
	add	    rax, r9
	adc	    rdx, r11
	mov	    [rcx], rdx
	adc	    rax, 0
	ret
	
%assign i 2
%rep    19      ; 3..20 inclusive
    mulredc i
    %assign i i + 1
%endrep

    end
