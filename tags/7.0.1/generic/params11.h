#ifndef MPZMOD_THRESHOLD
#define MPZMOD_THRESHOLD 170
#endif

#ifndef REDC_THRESHOLD
#define REDC_THRESHOLD 294
#endif

#ifndef MPN_MUL_LO_THRESHOLD_TABLE
#define MPN_MUL_LO_THRESHOLD_TABLE {0, 0, 0, 0, 0, 0, 0, 1, 7, 8, 1, 1, 8, 1, 1, 10, 1, 1, 1, 1, 1, 1, 1, 16, 1, 1, 16, 16, 1, 1, 16, 1}
#endif

#ifndef NTT_GFP_TWIDDLE_DIF_BREAKOVER
#define NTT_GFP_TWIDDLE_DIF_BREAKOVER 11
#endif

#ifndef NTT_GFP_TWIDDLE_DIT_BREAKOVER
#define NTT_GFP_TWIDDLE_DIT_BREAKOVER 11
#endif

#ifndef MUL_NTT_THRESHOLD
#define MUL_NTT_THRESHOLD 1024
#endif

#ifndef PREREVERTDIVISION_NTT_THRESHOLD
#define PREREVERTDIVISION_NTT_THRESHOLD 64
#endif

#ifndef POLYINVERT_NTT_THRESHOLD
#define POLYINVERT_NTT_THRESHOLD 512
#endif

#ifndef POLYEVALT_NTT_THRESHOLD
#define POLYEVALT_NTT_THRESHOLD 512
#endif

#ifndef MPZSPV_NORMALISE_STRIDE
#define MPZSPV_NORMALISE_STRIDE 512
#endif

/* 0:mulredc 1:mul+redc_1 2:mul+redc_2 3:mul+redc_n */
#ifndef TUNE_MULREDC_TABLE
#define TUNE_MULREDC_TABLE {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}
#endif

/* 0:mulredc 1:sqr+redc_1 2:sqr+redc_2 3:sqr+redc_n */
#ifndef TUNE_SQRREDC_TABLE
#define TUNE_SQRREDC_TABLE {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}
#endif

#ifndef LIST_MUL_TABLE
#define LIST_MUL_TABLE {0,0,0,0,0,0,0,0,0,0,3,3,3,3,3,3}
#endif
