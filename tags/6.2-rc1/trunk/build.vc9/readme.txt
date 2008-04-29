
Building GMP-ECM with Microsft Visual C++ (version 9)
=====================================================

If you wish to build the assembler code support you will need to install
the YASM assembler that is available at:

  http://www.tortall.net/projects/yasm/

You should ensure that the binary is named yasm.exe and put it in the same
directory as your Visual C++ compiler, which is typically:

C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\bin

You also need to install the enclosed yasm.rules file in a suitable path, 
for example:

C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\VCProjectDefaults

so that it is recognised by the Visual Studio 2008 IDE.

Tests
=====

The file tests.py is a python script that runs the GMP ECM tests. It runs 
the x64/release version by default but can be edited to test other builds.

    Brian Gladman, 16th April 2008
