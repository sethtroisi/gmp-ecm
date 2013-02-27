
%ifdef _WIN64
%include "a_x64_mulredc.asm"
%elif AMD_ASM
%include "a_win32a_mulredc.asm"
%else
%include "a_win32p_mulredc.asm"
%endif
