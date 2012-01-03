
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

The Multi-Precision Library - GMP and MPIR
==========================================

GMP-ECM works with either GMP or MPIR, a fork of GMP.  To build and run
GMP-ECM using Visual Studio you first need to obtain and build either 
GMP or MPIR.   MPIR has a fully integrated Visual Studio build system
for Windows but GMP does not.  

The VC++ build of GMP-ECM now defaults to MPIR but the property sheet 
mp_lib.vsprops can be edited to set the macro mp_lib to 'gmp' instead 
of 'mpir' to build ECM using GMP.

GMP
===

GMP can be built from the GMP source code available here:

  http://gmplib.org/
  
using the Visual Studio build files I provide here:

  http://fp.gladman.plus.com/computing/gmp4win.htm 
  
But these are based on GMP 4.2.x and are no longer being maintained.

GMP 4.3.x can be built using cygwin or mingw for win32 and it is reported 
that the resulting libraries work with Visual Studio when appropriately 
renamed. It may also be possible to build the generic C version of GMP for 
64-bit Windows systems using mingw64. But this version will be fairly slow
because it cannot use the fast assembler normally used by GMP because this
is not available in Windows format.

MPIR
====

MPIR is available here:

  http://www.mpir.org
  
It has full support for building MPIR for 32 and 64 bit Windows systems 
with x86 assembler support using the YASM assembler.  In particular it  
includes fast assembler code for modern AMD and Intel architectures 
running in 64-bit mode on Windows (not available in GMP).

Building GMP-ECM
================

The build files for GMP-ECM assume that the GMP and ECM build directories
are in a common parent directory as follows:

  Parent Directory
    MPIR (or GMP)
      build.vc9    -- MPIR (or GMP) build files
      ...
    GMP-ECM
      buid.vc9	   -- ECM build files 
      
The root directories for GMP and GMP-ECM are assumed to have these names
irrespective of which version is being used (they used to be followed by 
version numbers but this meant that the build projects had to be updated
too frequently). 

There are three build projects in build.vc9:

    ecm     - the ECM application 
    ecmlib  - the ECM library
    tune    - a program for tuning 
    
and each of these has the following configurations:

    win32\release-amd
    win32\release-intel
    win32\debug-amd     (not tune)
    win32\debug-intel   (not tune)
    x64\release-amd
    x64\release-intel
    x64\debug-amd       (not tune)
    x64\debug-intel     (not tune)

When a version of ecm and ecmlib are built the library and the application
are put in the directory matching the configuation that has been built:

   build.vc9\win32\release
   build.vc9\win32\debug
   build.vc9\x64\release
   build.vc9\x64\debug

If you don't want assembler support you need to change the define:      

#define NATIVE_REDC   1         

in config.h (in the build.vc9 subdirectory) to: 

#undef NATIVE_REDC

Tune
====

If tune is compiled and run for a particular configuration it will output a
file with appropriate parameters for this configuration with a name suuch as:

    ecm-params.h.win32.amd.new

To use this file when building ecm and ecmlib, remove the '.new' extension.

Tests
=====

The file tests.py is a python script that runs the ECM tests. It runs the
x64/release-amd version by default but can be edited to test other builds.

    Brian Gladman, 11th August 2009
