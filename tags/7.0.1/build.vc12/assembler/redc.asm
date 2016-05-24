%ifdef _WIN64
%include "a_x64_redc.asm"
%elif AMD_ASM
%include "a_win32a_redc.asm"
%else
%include "a_win32p_redc.asm"
%endif
