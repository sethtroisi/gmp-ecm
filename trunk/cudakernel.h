#ifndef _CUDAKERNEL_H
#define _CUDAKERNEL_H 1

#ifdef GPU_CC20
  #define MAJOR 2
  #define MINOR 0
  #define CURVES_BY_BLOCK 32//16
  #define BLOCKS_BY_MP 1//3
  #define CURVES_BY_MP CURVES_BY_BLOCK*BLOCKS_BY_MP
#endif

__global__ void Cuda_Ell_DblAdd (biguint_t *xarg, biguint_t *zarg, 
                  biguint_t *x2arg, biguint_t *z2arg, unsigned int firstinvd);


#endif /* _CUDAKERNEL_H */
