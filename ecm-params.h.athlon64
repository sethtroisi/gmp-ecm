/* updated 09 Mar 2012 on frite.loria.fr (AMD Phenom(tm) II X2 B55 Processor)
   for ecm svn 1724 with GMP 5.0.4 */

#ifndef HAVE_MPIR
/* 0:mulredc 1:mul+redc_1 2:mul+redc_2 3:mul+redc_n */
#define TUNE_MULREDC_TABLE {0,0,0,0,0,0,0,0,0,0,0,0,2,0,2,1,2,1,1,1,2}
/* 0:mulredc 1:sqr+redc_1 2:sqr+redc_2 3:sqr+redc_n */
#define TUNE_SQRREDC_TABLE {0,0,0,0,0,0,0,0,2,2,1,2,2,1,2,1,2,1,1,1,2}
#else /* tuning parameters for MPIR 2.5.0 */
/* 0:mulredc 1:mul+redc_1 2:mul+redc_2 3:mul+redc_n */
#define TUNE_MULREDC_TABLE {0,0,0,0,0,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1}
/* 0:mulredc 1:sqr+redc_1 2:sqr+redc_2 3:sqr+redc_n */
#define TUNE_SQRREDC_TABLE {0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}
#endif

#define MPZMOD_THRESHOLD 103
#define REDC_THRESHOLD 512
#define MPN_MUL_LO_THRESHOLD_TABLE {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 7, 0, 8, 10, 11, 11, 12, 12, 12, 13, 14, 15, 16, 17, 18, 19, 16, 18, 18, 18, 20}
#define NTT_GFP_TWIDDLE_DIF_BREAKOVER 12
#define NTT_GFP_TWIDDLE_DIT_BREAKOVER 17
#define MUL_NTT_THRESHOLD 256
#define PREREVERTDIVISION_NTT_THRESHOLD 16
#define POLYINVERT_NTT_THRESHOLD 512
#define POLYEVALT_NTT_THRESHOLD 128
#define MPZSPV_NORMALISE_STRIDE 512
