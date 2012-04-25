/* those parameters were obtained on gcc42.fsffrance.org with ecm-6.4.1-rc3
   gmp-5.0.2, and gcc 4.3.1 -O2 -mabi=n32 (mips64el-unknown-linux-gnu) */

/* 0:mulredc 1:mul+redc_1 2:mul+redc_2 3:mul+redc_n */
#define TUNE_MULREDC_TABLE {0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}
/* 0:mulredc 1:sqr+redc_1 2:sqr+redc_2 3:sqr+redc_n */
#define TUNE_SQRREDC_TABLE {0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}

#define MPZMOD_THRESHOLD 23
#define REDC_THRESHOLD 512
#define MPN_MUL_LO_THRESHOLD_TABLE {0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
#define NTT_GFP_TWIDDLE_DIF_BREAKOVER 12
#define NTT_GFP_TWIDDLE_DIT_BREAKOVER 13
#define MUL_NTT_THRESHOLD 1024
#define PREREVERTDIVISION_NTT_THRESHOLD 16
#define POLYINVERT_NTT_THRESHOLD 256
#define POLYEVALT_NTT_THRESHOLD 256
#define MPZSPV_NORMALISE_STRIDE 128
