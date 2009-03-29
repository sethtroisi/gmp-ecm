
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
      buid.vc9	   -- ECM build files 
      
The root directories for GMP and GMP-ECM are assumed to have these names
irreespective of which version is being used (they used to be followed by 
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
are put in the directory matching the configuation that has been built.

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

The file tests.py is a python script that runs the GMP ECM tests. It runs 
the x64/release-amd version by default but can be edited to test other builds.

GMP Version
===========

The previous VC++ build project for GMP-ECM assumed that the GMP root 
directory was GMP-4.2.x but this meant that the build project had to be
updated as the current GMP version changed.  The VC++ build projects now
assume that the GMP root directory is just GMP alone so GMP-ECM can be
linked without change to any GMP version provided only that its root
directory is named (or renamed) to 'GMP'. 

    Brian Gladman, 28th March 2009
