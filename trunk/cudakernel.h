#ifndef _CUDAKERNEL_H
#define _CUDAKERNEL_H 1

#ifdef GPU_CC20
  #define ECM_GPU_MAJOR 2
  #define ECM_GPU_MINOR 0
  #define ECM_GPU_CURVES_BY_BLOCK 32
  #define ECM_GPU_BLOCKS_BY_MP 1
  #define ECM_GPU_CURVES_BY_MP ECM_GPU_CURVES_BY_BLOCK * ECM_GPU_BLOCKS_BY_MP
#endif

#ifdef GPU_CC21
  #define ECM_GPU_MAJOR 2
  #define ECM_GPU_MINOR 1
  #define ECM_GPU_CURVES_BY_BLOCK 32
  #define ECM_GPU_BLOCKS_BY_MP 1
  #define ECM_GPU_CURVES_BY_MP ECM_GPU_CURVES_BY_BLOCK * ECM_GPU_BLOCKS_BY_MP
#endif

#ifdef GPU_CC30
  #define ECM_GPU_MAJOR 3
  #define ECM_GPU_MINOR 0
  #define ECM_GPU_CURVES_BY_BLOCK 32
  #define ECM_GPU_BLOCKS_BY_MP 1
  #define ECM_GPU_CURVES_BY_MP ECM_GPU_CURVES_BY_BLOCK * ECM_GPU_BLOCKS_BY_MP
#endif

__global__ void Cuda_Ell_DblAdd (biguint_t *xarg, biguint_t *zarg, 
                  biguint_t *x2arg, biguint_t *z2arg, unsigned int firstinvd);


#endif /* _CUDAKERNEL_H */
