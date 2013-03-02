;
; Part of GMP-ECM
;
; void ecm_redc3(
;       mp_limb_t       *z,     rdi  r8  <- rcx
;       const mp_limb_t *x,     rsi  r9  <- rdx
;       size_t           n,     rdx r10  <-  r8
;       mp_limb_t        m      rcx r11  <-  r9
;   )

%macro rloop 3
    mov     rax,[byte rsi+8*%3]
    mul     rbp
    add     [byte rdi+8*%3], %1
    adc     %2, rax
    mov     %1, rdx
    adc     %1, 0
%endmacro

	bits 64
	section .text

	global	ecm_redc3
%ifdef DLL
	export	ecm_redc3
%endif

PROC_FRAME	ecm_redc3
    push_reg    rbp
    push_reg    rbx
    push_reg    rsi
    push_reg    rdi
	alloc_stack	5*8
END_PROLOGUE
	mov     rdi, rcx
	mov     rsi, rdx
	mov     rdx, r8
	mov     rcx, r9

	mov     r8, rdi
	mov     r9, rsi
	mov     r10, rdx
	mov     r11, rcx

	mov	    rcx, r10
	mov	    [rsp], rcx
    cmp     rcx, 3
    jae     .unroll

.1: mov	    rbp, r11
	mov     rsi, r9
	imul    rbp, [rdi]
	mov	    r8, rdi
	mov	    rcx, r10
	xor	    rbx, rbx

.2: mov     rax, [rsi]
	add     rdi, 8
	mul     rbp
	add     rsi, 8
	add     rax, rbx
    adc     rdx, 0
	add     [rdi-8], rax
	adc     rdx, 0
	dec     rcx
	mov     rbx, rdx
	jnz     .2
	mov	    rdi, r8
	mov     [rdi], rbx
	dec	    qword [rsp]
	lea	    rdi, [rdi+8]
	jnz     .1

    add     rsp, 5*8
    pop     rdi
    pop     rsi
    pop     rbx
    pop     rbp
    ret

.unroll:
	mov     rdx, rcx
    dec     rcx
	sub     rdx, 2
	neg     rcx
	shr     rdx, 4
	and     rcx, 15
	mov     [rsp+16], rdx
	mov     rdx, rcx
	shl     rdx, 4
    lea     r10, [.loop_base wrt rip]
    add     rdx, r10
    lea     rdx, [rdx+rcx*4]
	add	    rdx, rcx
	neg     rcx
	mov	    r10, rcx
	mov	    [rsp+24], rdx

.4:	mov     rbp, r11
    mov     rsi, r9
    imul    rbp, [rdi]
    mov     r8, rdi
    mov     rcx, r10
	mov     rdx, [rsp+16]
	mov	    [rsp+8], rdx

    mov     rax, [rsi]
    lea     rsi, [rsi+rcx*8+8]
    mul     rbp
	lea     rdi, [rdi+rcx*8]
	mov     rbx, rdx

    mov     rdx, [rsp+24]
    test    rcx, 1
    mov     rcx, rax
	cmovnz  rcx, rbx
	cmovnz  rbx, rax
	jmp     rdx

    align   64

.5:	add     rdi, 128
.loop_base:
    rloop rcx, rbx,  0
    rloop rbx, rcx,  1
    rloop rcx, rbx,  2
    rloop rbx, rcx,  3
    rloop rcx, rbx,  4
    rloop rbx, rcx,  5
    rloop rcx, rbx,  6
    rloop rbx, rcx,  7
    rloop rcx, rbx,  8
    rloop rbx, rcx,  9
    rloop rcx, rbx, 10
    rloop rbx, rcx, 11
    rloop rcx, rbx, 12
    rloop rbx, rcx, 13
    rloop rcx, rbx, 14
    rloop rbx, rcx, 15

    dec     qword [rsp+8]
    lea     rsi, [rsi+128]
    jns     .5

    add     [rdi+128], rcx
    mov     rdi, r8
    adc     rbx, 0
    mov     [rdi], rbx
    dec     qword [rsp]
    lea     rdi, [rdi+8]
    jnz     .4

    add     rsp, 5*8
    pop     rdi
    pop     rsi
    pop     rbx
    pop     rbp
    ret
ENDPROC_FRAME

    end
