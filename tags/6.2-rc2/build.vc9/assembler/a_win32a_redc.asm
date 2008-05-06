;
; Part of GMP-ECM
;
; void ecm_redc3(
;       mp_limb_t       *z,     rdi  r8  <- rcx
;       const mp_limb_t *x,     rsi  r9  <- rdx
;       size_t           n,     rdx r10  <-  r8
;       mp_limb_t        m      rcx r11  <-  r9
;   )

%macro  seq 3
    mov     eax, [byte esi+4*%3]
    mul     ebp
    add     [byte edi+4*%3], %2
    adc     %1, eax
    mov     %2, edx
    adc     %2, 0
%endmacro

    text
    global _ecm_redc3
     
_ecm_redc3:
	push	ebp					
	push	edi
	push	esi
	push	ebx
	sub	    esp, 16
	mov	    ecx, [esp+44]
	mov	    edi, [esp+36]
	mov	    [esp], ecx
    cmp     ecx, 5
    jae     .3		

.1:	mov	    eax,[esp+48]
    mov     esi, [esp+40]
	mul	    dword [edi]				
	mov	    [esp+36], edi
	mov	    ecx, [esp+44]
	xor	    ebx, ebx			
	mov	    ebp, eax

.2: mov     eax, [esi]
    add     edi, 4
    mul     ebp			
    add     esi, 4
    add     eax, ebx
    adc     edx, 0		
    add     [edi-4], eax
    adc     edx, 0		
    dec     ecx			
    mov     ebx, edx
    jnz     .2		
    mov	    edi, [esp+36]
    mov     [edi], ebx
    dec	    dword [esp]
    lea	    edi, [edi+4]			
    jnz     .1				

	add	    esp, 16
	pop	    ebx
	pop	    esi
	pop	    edi
	pop 	ebp
	ret
	
.3: mov     edx, ecx
    dec     ecx
    sub     edx, 2
    neg     ecx	
    shr     edx, 4
    and     ecx, 15
    mov     [esp+8], edx
    mov     edx, ecx
    shl     edx, 4
    neg     ecx
    lea    edx, [edx+ecx+.6]
    mov	    [esp+44], ecx
    mov	    [esp+12], edx

.4:	mov     eax, [esp+48]
    mov     esi, [esp+40]
    mul     dword [edi]                          
    mov     [esp+36], edi                  
    mov     ecx, [esp+44]
    mov     edx, [esp+8]
    mov     ebp, eax
    mov	    [esp+4], edx
    mov     eax, [esi]
    lea     esi, [esi+ecx*4+4]
    mul     ebp
    lea     edi, [edi+ecx*4]
    mov     ebx, edx
    mov     edx, [esp+12]
    test    ecx, 1
    mov     ecx, eax
    cmovnz  ecx, ebx
    cmovnz  ebx, eax
    jmp     edx

    align   32
.5:	add     edi, 64
.6:	

%assign i 0
%rep 16
    %if (i & 1)
        seq     ecx, ebx,  i
    %else
        seq     ebx, ecx,  i
    %endif
    %assign i i + 1
%endrep

    dec     dword [esp+4]
    lea     esi, [esi+64]
    jns     .5

    add     [edi+64], ecx
    mov     edi, [esp+36]
    adc     ebx, 0
    mov     [edi], ebx
    dec     dword [esp]
    lea     edi, [edi+4]
    jnz     .4                      

    add     esp, 16
    pop     ebx
    pop     esi
    pop     edi
    pop     ebp
    ret

    end
    