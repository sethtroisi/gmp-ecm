
Building GMP-ECM with Microsoft Visual C++ 2008 (version 9)
===========================================================

If you wish to build the assembler code support you will need to install 
the YASM assembler that is available at:

  http://www.tortall.net/projects/yasm/

You should ensure that the binary is named yasm.exe and put it in the 
same directory as your Visual C++ compiler, which is typically:

C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\bin

You also need to install the enclosed yasm.rules file in a suitable path, 
for example:

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
    GMP-4.2.2
      build.vc9    -- GMP build files
      ...
    ECM-6.2
      buid.vc9	 -- ECM build files 

The ECM build project includes source files for 32-bit and 64-bit assembler
code builds.  It is set up to use these files optimised for AMD64. If you 
want 32-bit code optimised for Intel core2 processors you will need to open 
the build project in the Visual Studio IDE, exclude the 32-bit athlon assembler 
files and include those for Intel core2.  

If you don't want assembler support you need to change the define      

#define NATIVE_REDC   1         

in config.h in the build.vc9 subdirectory to: 

#undef NATIVE_REDC         

Tests
=====

The file tests.py is a python script that runs the GMP ECM tests. It runs 
the x64/release version by default but can be edited to test other builds.

    Brian Gladman, 1st May 2008
