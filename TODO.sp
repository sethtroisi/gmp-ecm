TODOs:

General:

1) Check for malloc errors.

2) grep FIXME *.c *.h

asm level:

All the SP code is *heavily* dependent on the C compiler. The performance gap
between unoptimized and -O3 code is huge. Maybe a couple of asm routines would
be appropriate if some kind soul were willing to write them.

1) spv_dotproduct. Given spvs a and b, compute 
   a[0] * b[0] + ... a[n-1] * b[n-1] mod p

   Useful in mpzspp_normalize.

2) spv_ntt_gfp_dif, spv_ntt_gfp_dit. Most of stage 2 is spent on these
   functions; maybe a hand-written version will outperform gcc's efforts.


SP level:

1) Maybe use 64->32 bit Montgomery reduction instead of GMP's
   udiv_qrnnd_preinv2norm. Preliminary tests indicate this gives a 10% speedup
   in ntt_gfp_dif. See sp_montmul in sp.h.

   UPDATE: This code has been written (but is not in cvs). It turns out that
           sp_montmul as presented is incorrect as we have to check for
	   overflows. Unfortunately this cancels out the 10% speedup.

SPV level:

1) toom3 and toom4 have been translated directly from toomcook.c but seem to be
   very slow. Try lowering the operation count (i.e. increase the number of
   muls) and see if we can fit a toom3/4 in-between karatsuba and ntt.

2) Better thresholds for when to start using ntt, especially for the "generic"
   poly mul routine.

4) Try out the NTT over GP(p^2). Old code doing this exists (DN), I'll see if
   I can get it to work at some point. This is why the function is called
   spv_ntt_gfp_dif, not just spv_ntt_dif.

5) Split-radix NTT. Old code doing this exists (DN), etc... I read somewhere
   that "A well-implemented radix-2 NTT can outperform a badly-implemented
   split-radix NTT"

MPZSPP level:

1) Find a better way of doing mpzspp_normalize.

2) Use a division tree for mpzspp_set_mpzp. Note however that these aren't
   too time-critical if we rewrite the high-level poly functions (see below)
   so they don't convert to/from mpzspps for every poly mul.

