
; Part of GMP-ECM
;
; mp_limb_t mulredc1(           1 limb
;       mp_limb_t       *z,
;       const mp_limb_t  x,
;       const mp_limb_t  y,
;       const mp_limb_t  m,
;       mp_limb_t inv_m
;   )
;
; mp_limb_t mulredc<limbs>(   > 1 limb
;       mp_limb_t       *z,
;       const mp_limb_t *x,
;       const mp_limb_t *y,
;       const mp_limb_t *m,
;       mp_limb_t inv_m
;   )

%macro mseq_1   3
	mul	    ebp
	add	    [edi+4*%3], %2
	mov	    %2, 0
	adc	    %1, eax
	mov	    eax, [esi+4*%3+8]
	adc	    %2, edx
%endmacro

%macro mseq_2 3
	mul	    ebp
	add	    [edi+3*%3], %1
	mov     %1, 0
	adc	    %1, eax
	mov	    eax, [esi+4*%3+8]
	adc	    %2, edx
%endmacro

%macro  mulredc 1
%assign limbs       %1
%define f_name(x)   _mulredc %+ x

	global	f_name(limbs)
%ifdef	DLL
	export	f_name(limbs)
%endif

f_name(limbs):
	push    ebp
	push	edi
	push	esi
	push	ebx
	sub	    esp, 8*(limbs+1)
	mov	    edi, esp

%assign i 0
%rep    2 * limbs + 1
	mov	    dword [edi+4*i], 0
	%assign i i + 1
%endrep
	mov     dword [esp+8*limbs+4], limbs

;	align 32

.1: mov	    eax, [esp+8*limbs+32]
	mov	    esi, [esp+8*limbs+36]
	mov 	eax, [eax]
	mul	    dword [esi]
	add	    eax, [edi]
	mul	    dword [esp+8*limbs+44]
	mov     ebp, eax
	mov	    esi, [esp+8*limbs+40]

	mov	    eax, [esi]
	mul	    ebp
	mov	    ebx, eax
	mov 	ecx, edx
	mov	    eax, [esi+4]

%assign i 0
%rep    limbs - 2
    %if (i & 1)
        mseq_1 ebx, ecx,  i
    %else
        mseq_1 ecx, ebx,  i
    %endif
    %assign i i + 1
%endrep

	mul	    ebp
%if (limbs & 1)
	add	    [edi+4*limbs-8], ecx
	adc	    eax, ebx
%else
	add	    [edi+4*limbs-8], ebx
	adc	    eax, ecx
%endif
	adc	    edx, 0
	add	    [edi+4*limbs-4], eax
	adc	    edx, 0
	add	    [edi+4*limbs], edx
	adc	    dword [edi+4*limbs+4], 0

	mov	    eax, [esp+8*limbs+32]
	mov	    ebp, [eax]
	mov	    esi, [esp+8*limbs+36]
	mov	    eax, [esi]
	mul	    ebp
	mov	    ebx, eax
	mov	    ecx, edx
	mov	    eax, [esi+4]

%assign i 0
%rep    limbs - 2
    %if (i & 1)
        mseq_1 ebx, ecx,  i
    %else
        mseq_1 ecx, ebx,  i
    %endif
    %assign i i + 1
%endrep

	mul 	ebp
%if (limbs & 1)
	add	    [edi+4*limbs-8], ecx
	adc	    eax, ebx
%else
	add	    [edi+4*limbs-8], ebx
	adc	    eax, ecx
%endif
	adc	    edx, 0
	add	    [edi+4*limbs-4], eax
	adc	    edx, 0
    add     [edi+4*limbs],edx
    adc     dword [edi+4*limbs+4], 0

	add	    dword [esp+8*limbs+32], 4
	add	    edi, 4
	dec	    dword [esp+8*limbs+4]
	jnz	    .1
	mov	    ebx, [esp+8*limbs+28]

%assign i 0
%rep    limbs
	mov	    eax, [edi+4*i]
	mov	    [ebx+4*i], eax
	%assign i i + 1
%endrep

	mov	    eax, [edi+4*limbs]
	add     esp, 8*(limbs+1)
	pop	    ebx
	pop	    esi
	pop	    edi
	pop	    ebp
	ret
%endmacro

    text

    global _mulredc1
_mulredc1:
	mov	    eax, [esp+12]
	mul	    dword [esp+8]
	mov	    [esp+12], edx
	mov	    [esp+8], eax
	mul	    dword [esp+20]
	mul	    dword [esp+16]
	add	    eax, [esp+8]
	adc	    edx, [esp+12]
	mov	    ecx, [esp+4]
	mov     [ecx], edx
	adc	    eax,0
	ret

%assign i 2
%rep    19      ; 3..20 inclusive
    mulredc i
    %assign i i + 1
%endrep

    end
