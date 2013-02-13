/* produced on pasta.loria.fr (Intel(R) Core(TM)2 CPU 6700  @ 2.66GHz) */

#ifndef HAVE_MPIR /* tuning parameters for GMP, tuned for GMP 5.0.4 */

/* 0:mulredc 1:mul+redc_1 2:mul+redc_2 3:mul+redc_n */
#define TUNE_MULREDC_TABLE {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0}
/* 0:mulredc 1:sqr+redc_1 2:sqr+redc_2 3:sqr+redc_n */
#define TUNE_SQRREDC_TABLE {0,0,0,0,0,0,0,2,0,2,0,2,1,1,1,1,2,2,1,2,2}
#define MPZMOD_THRESHOLD 21
#define REDC_THRESHOLD 512
#define MPN_MUL_LO_THRESHOLD_TABLE {0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 8, 8, 9, 10, 12, 11, 12, 13, 12, 12, 14, 16, 16, 16, 18, 18, 18}
#define NTT_GFP_TWIDDLE_DIF_BREAKOVER 17
#define NTT_GFP_TWIDDLE_DIT_BREAKOVER 17
#define MUL_NTT_THRESHOLD 256
#define PREREVERTDIVISION_NTT_THRESHOLD 8
#define POLYINVERT_NTT_THRESHOLD 128
#define POLYEVALT_NTT_THRESHOLD 256
#define MPZSPV_NORMALISE_STRIDE 128

#else /* tuning parameters for MPIR, tuned for MPIR 2.5.1 */

/* 0:mulredc 1:mul+redc_1 2:mul+redc_2 3:mul+redc_n */
#define TUNE_MULREDC_TABLE {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0}
/* 0:mulredc 1:sqr+redc_1 2:sqr+redc_2 3:sqr+redc_n */
#define TUNE_SQRREDC_TABLE {0,0,0,0,0,0,0,0,1,1,2,2,1,1,1,1,1,1,2,1,2}
#define MPZMOD_THRESHOLD 21
#define REDC_THRESHOLD 512
#define MPN_MUL_LO_THRESHOLD_TABLE {0, 0, 1, 1, 1, 1, 1, 0, 6, 6, 7, 8, 9, 9, 11, 10, 10, 11, 12, 13, 14, 14, 11, 13, 18, 18, 14, 20, 16, 18, 18, 20}
#define NTT_GFP_TWIDDLE_DIF_BREAKOVER 17
#define NTT_GFP_TWIDDLE_DIT_BREAKOVER 17
#define MUL_NTT_THRESHOLD 256
#define PREREVERTDIVISION_NTT_THRESHOLD 16
#define POLYINVERT_NTT_THRESHOLD 128
#define POLYEVALT_NTT_THRESHOLD 256
#define MPZSPV_NORMALISE_STRIDE 32

#endif
