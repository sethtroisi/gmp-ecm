
GWNUM bug:

There is a (very rare) bug buried in the gwnum fft code used by the P95/gwnum/ecmstag1.c file.
The bug shows itself (repeatably, on multiple machines, for PBMcL, anyway) when the following
GMP-ECM command is executed:

echo "(2^4242+1)/((2^2121-2^1061+1)*((2^707-2^354+1)*(2^303+2^152+1)*(2^21-2^11+1)/((2^101-2^51+1)*(2^7+2^4+1)))*45898441*21364316209019049762889*110142421356349014645013*16545613673048228577304069968130717)" | ./ecm -sigma 0:1666967755823658690 43e6 3e12

The run (using gwnum) completes with no factor found.

The same command with "-force-no-gwnum" added ends (as it should) with:

********** Factor found in step 2: 4395390156887311588965295007732214030665341
Found prime factor of 43 digits: 4395390156887311588965295007732214030665341
Composite cofactor ((2^4242+1)/((2^2121-2^1061+1)*((2^707-2^354+1)*(2^303+2^152+1)*(2^21-2^11+1)/((2^101-2^51+1)*(2^7+2^4+1)))*45898441*21364316209019049762889*110142421356349014645013*16545613673048228577304069968130717))/4395390156887311588965295007732214030665341 has 231 digits

BUG FIX:

George Woltman's file P95/ecm.cpp includes updated addition and duplication routines for
points on classic Montgomery curves. These routines have been adapted for use in the
P95/gwnum/ecmstag1.c file. The above bug disappears when the new file, included here,
is used to build the gwnum library linked with GMP-ECM.

To use the file, add the following (between steps 0) and 1)) to the build instructions
in INSTALL-gwnum:

0.5) Replace the file <P95 parent folder>/gwnum/ecmstag1.c with the file
   <ECM parent folder>/P95_ecm_stage1_file/ecmstag1.c. This file includes
   Woltman's updated Montgomery fft code and also code to use the file
   "Lchain_codes.dat" for generating Lucas chains.




















