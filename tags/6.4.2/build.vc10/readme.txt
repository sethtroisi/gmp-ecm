
Building GMP-ECM with Microsoft Visual C++ 2010 (version 10)
===========================================================

If you wish to build the assembler code support you will need to 
install the YASM assembler that is available at:

  http://www.tortall.net/projects/yasm/

THe version you need is vsyasm, which should be put it in the same
directory as your Visual C++ compiler, which is typically:

C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\bin


The Multi-Precision Library - GMP and MPIR
==========================================

GMP-ECM works with either GMP or MPIR, a fork of GMP. To build and run
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

  http://www.gladman.me.uk/computing/gmp4win.php
  
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
      build.vc10    -- MPIR (or GMP) build files
      ...
    GMP-ECM
      buid.vc10     -- ECM build files 
      
The root directories for GMP and GMP-ECM are assumed to have these names
irrespective of which version is being used (they used to be followed by 
version numbers but this meant that the build projects had to be updated
too frequently). 

There are three build projects in build.vc10:

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
are put in the directory matching the configuration that has been built:

    GMP-ECM
      build.vc10    -- ECM build files 
      lib           -- ECM static library files
      dll           -- ECM dynamic library files
      bin           -- ECM executable files
      
within these lib, dll and bin directories, the outputs are located in
sub-directories determined by the platform and configuration:
 
   win32\release
   win32\debug
   x64\release
   x64\debug

If you don't want assembler support you need to change the define:      

#define NATIVE_REDC   1         

in config.h (in the build.vc10 subdirectory) to: 

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

    Brian Gladman, 3rd January 2012

