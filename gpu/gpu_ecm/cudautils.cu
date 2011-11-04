#include "main.h"
#include "utils.h"
#include "cudaarith.h"

#define errCheck(err) cuda_errCheck (err, __FILE__, __LINE__)

inline void cuda_errCheck (cudaError err, const char *file, const int line)
{
  if( err != cudaSuccess ) 
  {
    fprintf(stderr, "%s(%i) : Error cuda : %s.\n",
              file, line, cudaGetErrorString( err) );
    exit(EXIT_FAILURE);
  }
}


void cuda_Main (biguint_t h_N, biguint_t h_invmod, biguint_t *h_xarray, 
                    biguint_t *h_zarray, biguint_t *h_x2array, 
                    biguint_t *h_z2array, mpz_t s, unsigned int firstinvd, 
                    unsigned int number_of_curves, FILE *OUTPUT_VERBOSE,
                    FILE *OUTPUT_VVERBOSE)
{
  cudaError_t err;
  
  biguint_t *d_xA;
  biguint_t *d_zA;
  biguint_t *d_xB;
  biguint_t *d_zB;

  size_t array_size = sizeof(biguint_t) * number_of_curves;

  dim3 dimBlock(NB_DIGITS,CURVES_BY_BLOCK);
  dim3 dimGrid(number_of_curves/CURVES_BY_BLOCK);

  fprintf(OUTPUT_VVERBOSE, "Block: %ux%ux%u Grid: %ux%ux%u\n", dimBlock.x, 
                      dimBlock.y, dimBlock.z, dimGrid.x, dimGrid.y, dimGrid.z);

  errCheck( cudaMalloc(&d_xA, array_size) );
  errCheck( cudaMalloc(&d_zA, array_size) );
  errCheck( cudaMalloc(&d_xB, array_size) );
  errCheck( cudaMalloc(&d_zB, array_size) );

#ifndef TEST
  //Copy into the gpu memory
  cuda_copy_cst(h_N, h_invmod);

  errCheck( cudaMemcpy((void *) d_xA, (void *) h_xarray, array_size, 
                                                    cudaMemcpyHostToDevice) );
  errCheck( cudaMemcpy((void *) d_zA, (void *) h_zarray, array_size, 
                                                    cudaMemcpyHostToDevice) );
  errCheck( cudaMemcpy((void *) d_xB, (void *) h_x2array, array_size, 
                                                    cudaMemcpyHostToDevice) );
  errCheck( cudaMemcpy((void *) d_zB, (void *) h_z2array, array_size, 
                                                    cudaMemcpyHostToDevice) );

  size_t j;

  for (j=mpz_sizeinbase(s,2)-1; j>0; j-- )
  {
    if (mpz_tstbit(s, j-1) ==1 )
      Cuda_Ell_DblAdd<<<dimGrid,dimBlock>>>(d_xB, d_zB, d_xA, d_zA, firstinvd);
    else
      Cuda_Ell_DblAdd<<<dimGrid,dimBlock>>>(d_xA, d_zA, d_xB, d_zB, firstinvd);

    //cudaThreadSynchronize();
    
    //maybe only for debug mode??
    err = cudaGetLastError(); 
    if (err != cudaSuccess )
    {
      fprintf(stderr, "%s(%i) : Error cuda : %s.\n",
              __FILE__, __LINE__, cudaGetErrorString( err) );
      exit(EXIT_FAILURE);
    }
    
  }

  errCheck( cudaMemcpy((void *) h_xarray, (void *) d_xA, array_size, 
                                                   cudaMemcpyDeviceToHost) );
  errCheck( cudaMemcpy((void *) h_zarray, (void *) d_zA, array_size, 
                                                   cudaMemcpyDeviceToHost) );

#else
  cuda_copy_cst(h_N, h_invmod);
  errCheck( cudaMemcpy(d_xB, h_xarray, array_size, cudaMemcpyHostToDevice) );
  errCheck( cudaMemcpy(d_zB, h_zarray, array_size, cudaMemcpyHostToDevice) );
  
  Cuda_Test<<<dimGrid,dimBlock>>>(d_xB,d_zB);
  errCheck( cudaMemcpy(h_xarray, d_xB, array_size, cudaMemcpyDeviceToHost) );
  errCheck( cudaMemcpy(h_zarray, d_zB, array_size, cudaMemcpyDeviceToHost) );

#endif


  cudaFree((void *) d_xA);
  cudaFree((void *) d_zA);
  cudaFree((void *) d_xB);
  cudaFree((void *) d_zB);

}
