#ifndef _CUDAKERNEL_H
#define _CUDAKERNEL_H 1

__global__ void Cuda_Ell_DblAdd (biguint_t *xarg, biguint_t *zarg, 
                  biguint_t *x2arg, biguint_t *z2arg, unsigned int firstinvd);


#endif /* _CUDAKERNEL_H */
