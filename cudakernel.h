#ifndef _CUDAKERNEL_H
#define _CUDAKERNEL_H 1

#ifdef __cplusplus
__global__ void Cuda_Ell_DblAdd (biguint_t *xarg, biguint_t *zarg, 
                  biguint_t *x2arg, biguint_t *z2arg, unsigned int firstinvd);
#endif


#ifdef __cplusplus
extern "C" {
#endif

int select_and_init_GPU (int, unsigned int*, int, int);
float cuda_Main (biguint_t, biguint_t, biguint_t, digit_t, biguint_t*,
                        biguint_t*, biguint_t*, biguint_t*, mpz_t, unsigned int,
                        unsigned int, int);
#ifdef __cplusplus
}
#endif


#endif /* _CUDAKERNEL_H */
