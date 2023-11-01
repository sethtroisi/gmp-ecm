/* coriandre.loria.fr, Intel(R) Core(TM) i7-8700, gcc 13.2.0 */

/* 0:mulredc 1:mul+redc_1 2:mul+redc_2 3:mul+redc_n */
#define TUNE_MULREDC_TABLE {0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,1,1,1}
/* 0:mulredc 1:sqr+redc_1 2:sqr+redc_2 3:sqr+redc_n */
#define TUNE_SQRREDC_TABLE {0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}
#define LIST_MUL_TABLE {0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,3,1,1,1,1,1,1,1,1,1,3,1,1,1}
#define MPZMOD_THRESHOLD 54
#define REDC_THRESHOLD 512
#define MPN_MUL_LO_THRESHOLD_TABLE {0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
#define NTT_GFP_TWIDDLE_DIF_BREAKOVER 17
#define NTT_GFP_TWIDDLE_DIT_BREAKOVER 17
#define MUL_NTT_THRESHOLD 16384
#define PREREVERTDIVISION_NTT_THRESHOLD 32
#define POLYINVERT_NTT_THRESHOLD 512
#define POLYEVALT_NTT_THRESHOLD 512
#define MPZSPV_NORMALISE_STRIDE 512
