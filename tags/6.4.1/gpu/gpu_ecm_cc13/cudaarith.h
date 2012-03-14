#include "utils.h"

__constant__ __device__ biguint_t Ncst;
//__constant__ __device__ biguint_t N2cst; //2*N
//__constant__ __device__ biguint_t N3cst; //3*N
//__constant__ __device__ biguint_t N4cst; //4*N
__constant__ __device__ biguint_t invmodcst;
//__constant__ __device__ biguint_t dcst[MAX_CURVES];

#ifdef CC13
__global__ void Cuda_Ell_Mul_2_3(biguint_t *xAarg, biguint_t *zAarg, biguint_t *xBarg, biguint_t *zBarg, biguint_t *xCarg, biguint_t *zCarg, const biguint_t *darg, const unsigned int power2, const unsigned int power3);
__global__ void Cuda_Ell_Add_Dbl(biguint_t *xAarg, biguint_t *zAarg, biguint_t *xBarg, biguint_t *zBarg, biguint_t *xCarg, biguint_t *zCarg, const biguint_t *darg);
__global__ void Cuda_Ell_3Add(biguint_t *xAarg, biguint_t *zAarg, biguint_t *xBarg, biguint_t *zBarg, biguint_t *xCarg, biguint_t *zCarg);
__global__ void Cuda_Ell_Add_Perm(biguint_t *xAarg, biguint_t *zAarg, biguint_t *xBarg, biguint_t *zBarg, biguint_t *xCarg, biguint_t *zCarg);
__global__ void Cuda_Ell_Dbl_3Add(biguint_t *xAarg, biguint_t *zAarg, biguint_t *xBarg, biguint_t *zBarg, biguint_t *xCarg, biguint_t *zCarg, const biguint_t *darg);
__global__ void Cuda_Ell_2Add_Dbl_Add(biguint_t *xAarg, biguint_t *zAarg, biguint_t *xBarg, biguint_t *zBarg, biguint_t *xCarg, biguint_t *zCarg, const biguint_t *darg, const int version);
__global__ void Cuda_Ell_Final(biguint_t *xAarg, biguint_t *zAarg, biguint_t *xBarg, biguint_t *zBarg, biguint_t *xCarg, biguint_t *zCarg, const biguint_t *darg);
#endif

#ifdef CC20
__global__ void Cuda_Ell_Mul_2_3(biguint_t *xAarg, biguint_t *zAarg, const biguint_t *darg, const unsigned int power2, const unsigned int power3);
__global__ void Cuda_PRAC(biguint_t *xAarg, biguint_t *zAarg, const biguint_t *darg, unsigned int PI, unsigned int B1, double val);
#endif

#ifdef TEST
__global__ void Cuda_Test(biguint_t *Aarg, biguint_t *Barg);
#endif

__host__ void cuda_copy_cst(biguint_t h_N, biguint_t h_invmod); 

__global__ void cpt(biguint_t *zBarg);
