/* those parameters were obtained on gcc61.fsffrance.org with ecm-6.4.1-rc3
   gmp-5.0.2, and gcc 4.4.1 -O2 -pedantic -mpa-risc-1-1
   (note that GMP must be configured with ABI=1.0, see
   http://gmplib.org/list-archives/gmp-bugs/2009-August/001585.html */

/* 0:mulredc 1:mul+redc_1 2:mul+redc_2 3:mul+redc_n */
#define TUNE_MULREDC_TABLE {0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}
/* 0:mulredc 1:sqr+redc_1 2:sqr+redc_2 3:sqr+redc_n */
#define TUNE_SQRREDC_TABLE {0,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1}

#define MPZMOD_THRESHOLD 49
#define REDC_THRESHOLD 512
#define MPN_MUL_LO_THRESHOLD_TABLE {0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
#define NTT_GFP_TWIDDLE_DIF_BREAKOVER 17
#define NTT_GFP_TWIDDLE_DIT_BREAKOVER 17
#define MUL_NTT_THRESHOLD 262144
#define PREREVERTDIVISION_NTT_THRESHOLD 262144
#define POLYINVERT_NTT_THRESHOLD 262144
#define POLYEVALT_NTT_THRESHOLD 262144
#define MPZSPV_NORMALISE_STRIDE 256
