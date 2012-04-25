/* tuned on confit.loria.fr (Intel(R) Core(TM) i5-2500 CPU) */

#ifndef HAVE_MPIR /* tuning parameters for GMP, tuned for GMP 5.0.4 */

/* 0:mulredc 1:mul+redc_1 2:mul+redc_2 3:mul+redc_n */
#define TUNE_MULREDC_TABLE {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1}
/* 0:mulredc 1:sqr+redc_1 2:sqr+redc_2 3:sqr+redc_n */
#define TUNE_SQRREDC_TABLE {0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1}
#define MPZMOD_THRESHOLD 21
#define REDC_THRESHOLD 512
#define MPN_MUL_LO_THRESHOLD_TABLE {0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 12, 13, 13, 13, 14, 14, 15, 16, 16, 17, 20, 22}
#define NTT_GFP_TWIDDLE_DIF_BREAKOVER 17
#define NTT_GFP_TWIDDLE_DIT_BREAKOVER 17
#define MUL_NTT_THRESHOLD 256
#define PREREVERTDIVISION_NTT_THRESHOLD 8
#define POLYINVERT_NTT_THRESHOLD 128
#define POLYEVALT_NTT_THRESHOLD 128
#define MPZSPV_NORMALISE_STRIDE 512

#else /* tuning parameters for MPIR, tuned for MPIR 2.5.1 */

/* 0:mulredc 1:mul+redc_1 2:mul+redc_2 3:mul+redc_n */
#define TUNE_MULREDC_TABLE {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,2,2,2}
/* 0:mulredc 1:sqr+redc_1 2:sqr+redc_2 3:sqr+redc_n */
#define TUNE_SQRREDC_TABLE {0,0,0,0,0,0,0,0,0,0,0,2,2,1,2,2,2,2,2,2,2}
#define MPZMOD_THRESHOLD 21
#define REDC_THRESHOLD 512
#define MPN_MUL_LO_THRESHOLD_TABLE {0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 11, 12, 13, 14, 15, 14, 16, 18, 18, 20, 18, 20}
#define NTT_GFP_TWIDDLE_DIF_BREAKOVER 17
#define NTT_GFP_TWIDDLE_DIT_BREAKOVER 17
#define MUL_NTT_THRESHOLD 256
#define PREREVERTDIVISION_NTT_THRESHOLD 16
#define POLYINVERT_NTT_THRESHOLD 128
#define POLYEVALT_NTT_THRESHOLD 256
#define MPZSPV_NORMALISE_STRIDE 256

#endif
