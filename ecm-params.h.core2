/* produced 11 Feb 2012 on coing.loria.fr with GMP 5.0.4 for ecm svn 1723
   (Intel(R) Core(TM)2 Duo CPU     U9600  @ 1.60GHz) */

#ifndef HAVE_MPIR
/* 0:mulredc 1:mul+redc_1 2:mul+redc_2 */
#define TUNE_MULREDC_TABLE {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
/* 0:mulredc 1:sqr+redc_1 2:sqr+redc_2 */
#define TUNE_SQRREDC_TABLE {0,0,0,0,0,0,0,0,1,1,2,2,2,2,2,2,2,2,2,2,2}
#else /* tuning parameters for MPIR 2.5.0 */
/* 0:mulredc 1:mul+redc_1 2:mul+redc_2 3:mul+redc_n */
#define TUNE_MULREDC_TABLE {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
/* 0:mulredc 1:sqr+redc_1 2:sqr+redc_2 3:sqr+redc_n */
#define TUNE_SQRREDC_TABLE {0,0,0,0,0,0,0,0,1,1,0,1,1,1,1,1,1,1,1,1,1}
#endif

#define MPZMOD_THRESHOLD 81
#define REDC_THRESHOLD 491
#define MPN_MUL_LO_THRESHOLD_TABLE {0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 1, 1, 10, 10, 1, 12, 12, 12, 16, 12, 13, 14, 15, 16, 17, 18, 19}
#define NTT_GFP_TWIDDLE_DIF_BREAKOVER 16
#define NTT_GFP_TWIDDLE_DIT_BREAKOVER 17
#define MUL_NTT_THRESHOLD 512
#define PREREVERTDIVISION_NTT_THRESHOLD 16
#define POLYINVERT_NTT_THRESHOLD 128
#define POLYEVALT_NTT_THRESHOLD 256
#define MPZSPV_NORMALISE_STRIDE 64
