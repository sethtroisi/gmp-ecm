#include "utils.h"

__constant__ __device__ biguint_t invmodcst;
__device__ biguint_t Ncst;

#ifdef TEST
__global__ void Cuda_Test(biguint_t *A, biguint_t *B);
#else
__global__ void Cuda_Ell_DblAdd(biguint_t *xarg, biguint_t *zarg, biguint_t *x2arg, biguint_t *z2arg, unsigned int firstd);
#endif

__host__ void cuda_copy_cst(biguint_t h_N, biguint_t h_invmod);
