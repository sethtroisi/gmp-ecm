#include "main.h"
#include "utils.h"
#include "cudaarith.h"

clock2_t cuda_Main(biguint_t h_N, biguint_t h_invmod, biguint_t *h_xarray, biguint_t *h_zarray, biguint_t *h_x2array, biguint_t *h_z2array, unsigned int B1, unsigned int firstinvd, unsigned int number_of_curves)
{
	clock_t begingpu,temp1,endgpu;
	clock2_t gputime;
	begingpu=clock();
	
	biguint_t *d_xA;
	biguint_t *d_zA;
	biguint_t *d_xB;
	biguint_t *d_zB;
	mpz_t s;
	mpz_init_set_ui(s,1);

	dim3 dimBlock(SIZE_NUMBER,CURVES_BY_BLOCK);
	dim3 dimGrid(number_of_curves/CURVES_BY_BLOCK);

	fprintf(stdout,"%u %u %u %u %u %u\n",dimBlock.x,dimBlock.y,dimBlock.z,dimGrid.x,dimGrid.y,dimGrid.z);

	cudaMalloc(&d_xA,number_of_curves*sizeof(biguint_t));
	cudaMalloc(&d_zA,number_of_curves*sizeof(biguint_t));
	cudaMalloc(&d_xB,number_of_curves*sizeof(biguint_t));
	cudaMalloc(&d_zB,number_of_curves*sizeof(biguint_t));
	//Copy into the gpu memory

#ifndef TEST
	int j;
	cuda_copy_cst(h_N, h_invmod);
	cudaMemcpy(d_xA, h_xarray, sizeof(biguint_t)*number_of_curves, cudaMemcpyHostToDevice) ;
	cudaMemcpy(d_zA, h_zarray, sizeof(biguint_t)*number_of_curves, cudaMemcpyHostToDevice) ;
	cudaMemcpy(d_xB, h_x2array, sizeof(biguint_t)*number_of_curves, cudaMemcpyHostToDevice) ;
	cudaMemcpy(d_zB, h_z2array, sizeof(biguint_t)*number_of_curves, cudaMemcpyHostToDevice) ;

	endgpu=clock();
	compute_s(s,B1);

	temp1=clock();
	//cudaSetDeviceFlags(cudaDeviceScheduleSpin);	
	//cudaSetDeviceFlags(cudaDeviceBlockingSync);	
	gmp_fprintf(stdout,"#s has %lu bits\nCompute in %.3f sec\n", 
		mpz_sizeinbase(s,2), (double) (temp1-endgpu)/CLOCKS_PER_SEC);


	for (j=mpz_sizeinbase(s,2)-2;j>=0;j--)
	{
		if (mpz_tstbit(s,j)==1)
			Cuda_Ell_DblAdd<<<dimGrid,dimBlock>>>(d_xB,d_zB,d_xA,d_zA,firstinvd);
		else
			Cuda_Ell_DblAdd<<<dimGrid,dimBlock>>>(d_xA,d_zA,d_xB,d_zB,firstinvd);
	}

	fprintf(stdout,"#All kernels launched, waiting for results...\n");
	cudaMemcpy(h_xarray, d_xA, sizeof(biguint_t)*number_of_curves, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_zarray, d_zA, sizeof(biguint_t)*number_of_curves, cudaMemcpyDeviceToHost);

#else
	cuda_copy_cst(h_N, h_invmod);
	cudaMemcpy(d_xB, h_xarray, sizeof(biguint_t)*number_of_curves, cudaMemcpyHostToDevice) ;
	cudaMemcpy(d_zB, h_zarray, sizeof(biguint_t)*number_of_curves, cudaMemcpyHostToDevice) ;
	
	temp1=clock();
	Cuda_Test<<<dimGrid,dimBlock>>>(d_xB,d_zB);
	cudaMemcpy(h_xarray, d_xB, sizeof(biguint_t)*number_of_curves, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_zarray, d_zB, sizeof(biguint_t)*number_of_curves, cudaMemcpyDeviceToHost);

#endif


	cudaFree(d_xA);
	cudaFree(d_zA);
	cudaFree(d_xB);
	cudaFree(d_zB);
	mpz_clear(s);

	cudaThreadExit();
	endgpu=clock();

	gputime.init=temp1-begingpu;
	gputime.computation=endgpu-temp1;
	return gputime;
}
