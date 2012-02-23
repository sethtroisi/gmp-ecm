#ifdef CC20
  #define MAJOR 2
  #define MINOR 0
  #define CURVES_BY_BLOCK 32//16
  #define BLOCKS_BY_MP 1//3
  #define CURVES_BY_MP CURVES_BY_BLOCK*BLOCKS_BY_MP
#endif

#ifdef CC13
  #define MAJOR 1
  #define MINOR 3
  #define CURVES_BY_BLOCK 16
  #define BLOCKS_BY_MP 1
  #define CURVES_BY_MP CURVES_BY_BLOCK*BLOCKS_BY_MP
#endif

#define errCheck(err) cuda_errCheck (err, __FILE__, __LINE__)

__global__ void Cuda_Ell_DblAdd (biguint_t *xarg, biguint_t *zarg, 
                  biguint_t *x2arg, biguint_t *z2arg, unsigned int firstinvd);
