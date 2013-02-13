Architecture-specifc assembly code is stored in different subdirectories.
Currently (March 2006), these are
  pentium4
  x86_64
  athlon

The code for pentium4 uses MMX/SSE2 instructions, and therefore can not
run on old x86. The code in the 'athlon' subdir is pure i486 and can
therefore be used for any x86 but the asm is optimized for athlon.

In these subdirs, there is size-specific asm code for combined
multiplication and redc. The sizes are currently 1 to 20 limbs. If
needed, one could go to higher sizes, but is there a need?  There is also
a redc function coded in asm (without mul).

The files are automatically generated using a Python script. This
generation is not done at configure or compile time, to avoid a
dependency to Python. However, the script is given, for developpers to
play with.

At configure, if asm-redc is enabled, symbolic links are done to the
.asm files in the appropriate directory. 
