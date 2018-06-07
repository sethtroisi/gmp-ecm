
Building GMP-ECM with Microsoft Visual C++ 2017
===============================================

If you wish to build the assembler code support you will need to 
install the YASM assembler that is available at:

    http://www.tortall.net/projects/yasm/

THe version you need is vsyasm, which should be put it in the directory

    C:\Program Files\yasm

Alternatively vsyasm can be installed anywhere provided the environment
variable YASMPATH is set to its absolute file path.

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
  
GMP can be built with mingw for 32-bit Windows and mingw64 for Windows x64.
It is reported that the resulting libraries work with Visual Studio when
appropriately renamed. 

MPIR
====

MPIR is available here:

  http://www.mpir.org
  
It has full support for building MPIR for 32 and 64 bit Windows systems 
with x86 assembler support using the YASM assembler.  

Building GMP-ECM
================

The build files for GMP-ECM assume that the GMP and ECM build directories
are in a common parent directory as follows:

  Parent Directory
    MPIR (or GMP)
      build.vc15    -- MPIR (or GMP) build files
      ...
    GMP-ECM
      buid.vc15     -- ECM build files 
      
The root directories for GMP and GMP-ECM are assumed to have these names
irrespective of which version is being used (they used to be followed by 
version numbers but this meant that the build projects had to be updated
too frequently). 

The normal (non GPU) build is opened by loading the file ecm.sln (from 
the build.vc14 directory) into Visual Studio. This provides these build
projects in build.vc15 for the non GPU build:

    ecm           - the ECM application 
    ecmlib        - the ECM library
    tune          - a program for tuning 
    bench_mulredc - for benchmarking mulredc 
    multiecm      - work in progress (not working)

The GPU build is opened by loading the file ecm.sln (from  the build.vc14
directory) into Visual Studio. This provides two build projects in 
build.vc15:

    ecm_gpu     - the ECM application 
    ecmlib_gpu  - the ECM library

In all cases you have to choose either a win32 or x64 build and either a
Release or Debug configuration (however the win32 builds are no longer 
actively supported and may not work). 

The non GPU Build
-----------------

Before starting a build, there are a number of configuration options
that need to be set:

1. If you wish to compile GMP-ECM for use on a particular processor,
   select the appropriate define from the file 'ecm-params.h' in the
   GMP-ECM root directory and decide which of the defines suit your
   needs (e.g. __tune_corei7__).  Then replace the existing define:
   
   /* define Windows tuning here */
   #  define __tune_corei7__
 
   towards the end of the file config.h file in the 'build.vc14'
   directory (build.vc14\config.h) with the chosen define.

2. The file at 'build.vc14\mul_fft-params.h' allows the FFT code to
   be tuned to 32 or 64-bit systems by selecting an option by 
   changing the appropriate '#elif 0' to #elif 1'.   If you wish to 
   use the win32 AMD assembler files, you also have to use the 
   Visual Studio property page to define AMD_ASM (alternatively
   you can edit the two files mulredc.asm and redc.asm in the 
   build.vc14\assembler\ directory to include the AMD assembler).

The GPU Build
-------------

1. If you wish to build with a GPU capability you will need to 
   install Nvidia Nsight for Visual Studio version 5.4 and the
   CUDA Toolkit v9.0.  You can then build the libecm_gpu and
   ecm_gpu projects
   
2. The choices above for the non GPU build aslo apply when
   building for a GPU based system. 
   
   By default, the GPU  configuration is "compute_50,sm_50". If
   you need to change this, select libecm_gpu and ecm_gpu and 
   set the propertiesfor "CUDA C/C++|Device|Code Generation" for
   your GPU capability. 
   
   Also under "C/C++|Preprocessor|Preprocessor Definitions" for 
   both these projects, change the current definition GPU_CC50 to 
   that for your GPU capability

Build Configurations
--------------------

When a version of ecm and ecmlib are built, the library and the application
are put in the directory matching the configuration that has been built:

    GMP-ECM
      build.vc15    -- ECM build files 
      lib           -- ECM static library files
      bin           -- ECM executable files
      
within these lib, dll and bin directories, the outputs are located in
sub-directories determined by the platform and configuration:
 
   win32\release
   win32\debug
   x64\release
   x64\debug

If you don't want assembler support you need to change the define:      

#define NATIVE_REDC   1         

in config.h (in the build.vc14 subdirectory) to: 

#undef NATIVE_REDC

Tune
====

If tune is compiled and run for a particular configuration it will output
suitable values for optimising GMP-ECM to the console window.  To optimise
GMP-ECM these values should be put in a suitably named file whcih then has
to be integrated in ecm-params.h.

Tests
=====

The file test.py is a python script that runs the ECM tests. It runs the
x64/release-amd (non GPU) version by default but can be edited to test other
builds.  It cannot run some tests as a result of the diifficulty in the
conversion of the Unix shell scripts for the tests for use on Windows. 

    Brian Gladman, November 2017
