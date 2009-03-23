
Building GMP-ECM with Microsoft Visual C++ 2008 (version 9)
===========================================================

If you wish to build the assembler code support you will need to 
install the YASM assembler that is available at:

  http://www.tortall.net/projects/yasm/

You should ensure that the binary is named yasm.exe and put it in 
the same directory as your Visual C++ compiler, which is typically:

C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\bin

You also need to install the enclosed yasm.rules file in a suitable 
path, for example:

C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\VCProjectDefaults

so that it is recognised by the Visual Studio 2008 IDE.

GMP
===

To build and run GMP-ECM using Visual Studio you first need to obtain 
and build GMP. GMP is available here:

  http://gmplib.org/

and the build files for GMP using Visual Studio 2008 are provided here:

  http://fp.gladman.plus.com/computing/gmp4win.htm 

The build files for GMP-ECM assume that the GMP and ECM build directories
are in a common parent directory as follows:

  Parent Directory
    GMP
      build.vc9    -- GMP build files
      ...
    GMP-ECM
      buid.vc9	 -- ECM build files 
      
The root directories for GMP and GMP-ECM are assumed to have these names
irreespective of which version is being used (they used to be followed by 
version numbers but this meant that the build projects had to be updated
too frequently). 

The ECM build project includes source files for 32-bit and 64-bit assembler
code builds.  

The default for 64-bit builds is to use the assembler code provided in the 
files:

    a_x64_mulredc.asm
    a_x64_redc.asm

For 32-bit builds there are two alternative assembler file sets:

AMD64:
    a_win32a_mulredc.asm
    a_win32a_redc.asm

Intel Core2:
    a_win32p_mulredc.asm
    a_win32p_redc.asm

the first being the default. If you wish to build with 32-bit core2 
assembler support, you will need to open the build project in the Visual
Studio IDE, exclude the 32-bit AMMD64 assembler code files and include 
those for Intel core2.  

If you don't want assembler support you need to change the define:      

#define NATIVE_REDC   1         

in config.h (in the build.vc9 subdirectory) to: 

#undef NATIVE_REDC

Tests
=====

The file tests.py is a python script that runs the GMP ECM tests. It runs 
the x64/release version by default but can be edited to test other builds.

GMP Version
===========

The previous VC++ build project for GMP-ECM assumed that the GMP root 
directory was GMP-4.2.x but this meant that the build project had to be
updated as the current GMP version changed.  The VC++ build projects now
assume that the GMP root directory is just GMP alone so GMP-ECM can be
linked without change to any GMP version provided only that its root
directory is named (or renamed) to 'GMP'. 

    Brian Gladman, 24th November 2008
