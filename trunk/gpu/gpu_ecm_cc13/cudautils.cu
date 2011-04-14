#include "cudaarith.h"

clock2_t cuda_Main(biguint_t h_invmod, biguint_t h_N, biguint_t *h_xarray, biguint_t *h_zarray,biguint_t *h_darray, unsigned int B1, unsigned int number_of_curves,int device)
{
	clock_t begingpu,temp1,endgpu;
	clock2_t gputime;
	begingpu=clock();
	unsigned int PI=3,pp;

	//variables for prac
  unsigned int d,i;
  double cmin,c;
#ifdef CC13
  unsigned int e, r;
#endif

  static double val[NV] =
    { 0.61803398874989485, 0.72360679774997897, 0.58017872829546410,
      0.63283980608870629, 0.61242994950949500, 0.62018198080741576,
      0.61721461653440386, 0.61834711965622806, 0.61791440652881789,
      0.61807966846989581};

	biguint_t *d_xA;
	biguint_t *d_zA;
#ifdef CC13
	biguint_t *d_xB;
	biguint_t *d_zB;
	biguint_t *d_xC;
	biguint_t *d_zC;
	biguint_t *tmp;//for swapping
#endif

	dim3 dimBlock(SIZE_NUMBER);
	dim3 dimGrid(number_of_curves);

	//int deviceCount;
	//int bestDevice=-1,bestMajor=-1,bestMinor=-1;
	cudaDeviceProp deviceProp;
	//cudaGetDeviceCount(&deviceCount);
				
	printf("#Compiled for a NVIDIA GPU with compute capability at least %d.%d.\n",MAJOR,MINOR);
/*
	if (device==-1)
	{
		for (device = 0; device < deviceCount; ++device) 
		{
			cudaGetDeviceProperties(&deviceProp, device);
			if (device == 0) 
			{
				if (deviceProp.major == 9999 && deviceProp.minor == 9999)
					printf("#There is no device supporting CUDA.\n");
				else if (deviceCount == 1)
					printf("#There is 1 device supporting CUDA :\n");
				else
				printf("#There are %d devices supporting CUDA :\n",deviceCount);
			}
	
			printf("#%d) %s (compute capability %d.%d)\n",device,deviceProp.name,deviceProp.major,deviceProp.minor);

			if (deviceProp.major > MAJOR || (deviceProp.major==MAJOR && deviceProp.minor >= MINOR))
			{
				if (deviceProp.major > bestMajor || (deviceProp.major==bestMajor && deviceProp.minor > bestMinor))
				{
				bestDevice=device;
				bestMajor=deviceProp.major;
				bestMinor=deviceProp.minor;
				}
			}
		}

		if (bestDevice<0)
		{
			printf("#Error : there is no device with compute capability at least %d.%d.\n",MAJOR,MINOR);
			exit(1);
		}
		cudaSetDevice(bestDevice);
	}
*/	
	if (device!=-1)
	{
		printf("#Device %d was required.\n",device);
		cudaGetDeviceProperties(&deviceProp,device);
		if (deviceProp.major < MAJOR || (deviceProp.major== MAJOR && deviceProp.minor< MINOR))
		{
			printf("#Error : Device %d have a compute capability of %d.%d (required %d.%d).\n",device,deviceProp.major,deviceProp.minor,MAJOR,MINOR);
			exit(1);
		}
		else
			cudaSetDevice(device);
	}
	
	cudaGetDevice(&device);
	cudaGetDeviceProperties(&deviceProp,device);
	printf("#Will use device %d : %s, compute capability %d.%d, %d MPs.\n",device,deviceProp.name,deviceProp.major,deviceProp.minor,deviceProp.multiProcessorCount);



	cudaMalloc(&d_xA,number_of_curves*sizeof(biguint_t));
	cudaMalloc(&d_zA,number_of_curves*sizeof(biguint_t));
#ifdef CC13
	cudaMalloc(&d_xB,number_of_curves*sizeof(biguint_t));
	cudaMalloc(&d_zB,number_of_curves*sizeof(biguint_t));
	cudaMalloc(&d_xC,number_of_curves*sizeof(biguint_t));
	cudaMalloc(&d_zC,number_of_curves*sizeof(biguint_t));
#endif

	//printf("#size %lu\n",sizeof(unsigned int));	

	//Copy into the gpu memory
	cuda_copy_cst(h_N,h_invmod,h_darray,number_of_curves);
	cudaMemcpy(d_xA, h_xarray, sizeof(biguint_t)*number_of_curves, cudaMemcpyHostToDevice) ;
	cudaMemcpy(d_zA, h_zarray, sizeof(biguint_t)*number_of_curves, cudaMemcpyHostToDevice) ;

	temp1=clock();
#ifndef TEST
	//cudaSetDeviceFlags(cudaDeviceScheduleSpin);	
	//cudaSetDeviceFlags(cudaDeviceBlockingSync);	
	unsigned int power2,power3;

#ifdef CC13
	//P=2 and P=3 are treated separately
	power2=0;
	pp=2;
	while(pp<=B1)
	{
		power2++;
		pp*=2;
	}
	power3=0;
	pp=3;
	while(pp<=B1)
	{
		power3++;
		pp*=3;
	}
	//printf("#Calling gpu for prime 2 with power %u 3 with power %u...\n",power2,power3);	
	Cuda_Ell_Mul_2_3<<<dimGrid,dimBlock>>>(d_xA,d_zA,d_xB,d_zB,d_xC,d_zC,power2,power3);
	//cudaThreadSynchronize();


	//for prime >=5
	PI=getprime(PI);
	PI=getprime(PI);
	do
	{
		pp=PI;
		i=0;
			
		for (d=0, cmin=ADD* (double) PI; d < NV; d++)
		{
			c = lucas_cost (PI,val[d]);
			if (c < cmin)
			{
				cmin=c;
				i = d;
			}
		}

		while(pp<=B1)
		{
  		d = PI;
  		r = (unsigned int) ((double) d * val[i] + 0.5);
  		d = PI - r;
  		e = 2 * r - PI;

			//printf("#Calling gpu for prime %u d=%u e=%u...\n",PI,d,e);	

			while (d != e)
    	{
 				if (d < e)
  			{
	  			r = d;
    			d = e;
   	 			e = r;
					tmp=d_xB;
					d_xB=d_xA;
					d_xA=tmp;
					tmp=d_zB;
					d_zB=d_zA;
					d_zA=tmp;
  			}
      
				// do the first line of Table 4 whose condition qualifies 
  			if (d - e <= e / 4 && ((d + e) % 3) == 0)
  			{ // condition 1 
  				d = (2 * d - e) / 3;
    			e = (e - d) / 2;
					Cuda_Ell_3Add<<<dimGrid,dimBlock>>>(d_xA,d_zA,d_xB,d_zB,d_xC,d_zC);
				}
  
				else if (d - e <= e / 4 && (d - e) % 6 == 0)
  			{ // condition 2 
  				d = (d - e) / 2;
					Cuda_Ell_Add_Dbl<<<dimGrid,dimBlock>>>(d_xA,d_zA,d_xB,d_zB,d_xC,d_zC);
  			}
  
				else if ((d + 3) / 4 <= e)
  			{ // condition 3 
  				d -= e;
					Cuda_Ell_Add_Perm<<<dimGrid,dimBlock>>>(d_xA,d_zA,d_xB,d_zB,d_xC,d_zC);
				}
  
				else if ((d + e) % 2 == 0)
				{ // condition 4 
  				d = (d - e) / 2;
					Cuda_Ell_Add_Dbl<<<dimGrid,dimBlock>>>(d_xA,d_zA,d_xB,d_zB,d_xC,d_zC);
				}
  
				// now d+e is odd 
  			else if (d % 2 == 0)
				{ // condition 5 
  				d /= 2;
					Cuda_Ell_Add_Dbl<<<dimGrid,dimBlock>>>(d_xA,d_zA,d_xC,d_zC,d_xB,d_zB);
				}
  
				// now d is odd, e is even 
  			else if (d % 3 == 0)
				{ // condition 6 
  				d = d / 3 - e;
					Cuda_Ell_Dbl_3Add<<<dimGrid,dimBlock>>>(d_xA,d_zA,d_xB,d_zB,d_xC,d_zC);
				}
  
				else if ((d + e) % 3 == 0)
				{ // condition 7 
  				d = (d - 2 * e) / 3;
					Cuda_Ell_2Add_Dbl_Add<<<dimGrid,dimBlock>>>(d_xA,d_zA,d_xB,d_zB,d_xC,d_zC,1);
				}
  
				else if ((d - e) % 3 == 0)
				{ // condition 8 
  				d = (d - e) / 3;
					Cuda_Ell_2Add_Dbl_Add<<<dimGrid,dimBlock>>>(d_xA,d_zA,d_xB,d_zB,d_xC,d_zC,2);
				}
					
				else // necessarily e is even here 
				{ // condition 9 
  				e /= 2;
					Cuda_Ell_Add_Dbl<<<dimGrid,dimBlock>>>(d_xB,d_zB,d_xC,d_zC,d_xA,d_zA);
				}
				//cudaThreadSynchronize();
			}
			Cuda_Ell_Final<<<dimGrid,dimBlock>>>(d_xA,d_zA,d_xB,d_zB,d_xC,d_zC);
			//cudaThreadSynchronize();
			pp*=PI;
		}
		
		PI=getprime(PI);
	}	while(PI<=B1);

#endif

#ifdef CC20
	//P=2 and P=3 are treated separately
	power2=0;
	pp=2;
	while(pp<=B1)
	{
		power2++;
		pp*=2;
	}
	power3=0;
	pp=3;
	while(pp<=B1)
	{
		power3++;
		pp*=3;
	}
	printf("#Calling gpu for prime 2 with power %u 3 with power %u...\n",power2,power3);	
	Cuda_Ell_Mul_2_3<<<dimGrid,dimBlock>>>(d_xA,d_zA,power2,power3);
	cudaThreadSynchronize();
	//cudaError_t error=cudaGetLastError();
	//if (error != cudaSuccess)
	//	fprintf(stderr, "%s(%i) : cudaSafeCall() Runtime API error : %s.\n",__FILE__,__LINE__,cudaGetErrorString(error));

	//for prime >=5
	PI=getprime(PI);
	PI=getprime(PI);
	while(PI<=B1)
	{
	
		i=0;
			
		for (d=0, cmin=ADD* (double) PI; d < NV; d++)
		{
			c = lucas_cost (PI,val[d]);
			if (c < cmin)
			{
				cmin=c;
				i = d;
			}
		}
				
		printf("#Calling gpu for prime %u ...\n",PI);	
		Cuda_PRAC<<<dimGrid,dimBlock>>>(d_xA,d_zA,PI,B1,val[i]);
		cudaThreadSynchronize();
		
		PI=getprime(PI);
	}	

#endif

	getprime(0);
	printf("#All kernels launched, waiting for results...\n");
#ifdef CC13
	cudaMemcpy(h_xarray, d_xB, sizeof(biguint_t)*number_of_curves, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_zarray, d_zB, sizeof(biguint_t)*number_of_curves, cudaMemcpyDeviceToHost);
#endif
#ifdef CC20
	cudaMemcpy(h_xarray, d_xA, sizeof(biguint_t)*number_of_curves, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_zarray, d_zA, sizeof(biguint_t)*number_of_curves, cudaMemcpyDeviceToHost);
#endif

#endif

#ifdef TEST
	Cuda_Test<<<dimGrid,dimBlock>>>(d_xA,d_zA);
	cudaMemcpy(h_xarray, d_xA, sizeof(biguint_t)*number_of_curves, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_zarray, d_zA, sizeof(biguint_t)*number_of_curves, cudaMemcpyDeviceToHost);
#endif


	cudaFree(d_xA);
	cudaFree(d_zA);
#ifdef CC13
	cudaFree(d_xB);
	cudaFree(d_zB);
	cudaFree(d_xC);
	cudaFree(d_zC);
#endif

	cudaThreadExit();
	endgpu=clock();

	gputime.init=temp1-begingpu;
	gputime.computation=endgpu-temp1;
	return gputime;
}

// returns the number of modular multiplications for computing
// V_n from V_r * V_{n-r} - V_{n-2r}.
double lucas_cost (unsigned long n, double v)
{
  unsigned long d, e, r;
  double c; // cost

  d = n;
  r = (unsigned long) ((double) d * v + 0.5);
  if (r >= n)
    return (ADD * (double) n);
  d = n - r;
  e = 2 * r - n;
  c = DUP + ADD; // initial duplicate and final addition 
	while (d != e)
    {
      if (d < e)
        {
          r = d;
          d = e;
          e = r;
        }
      if (d - e <= e / 4 && ((d + e) % 3) == 0)
        { //condition 1 
          d = (2 * d - e) / 3;
          e = (e - d) / 2;
          c += 3.0 * ADD; //3 additions 
        }
      else if (d - e <= e / 4 && (d - e) % 6 == 0)
        { //condition 2 
          d = (d - e) / 2;
          c += ADD + DUP; //one addition, one duplicate 
        }
      else if ((d + 3) / 4 <= e)
        { //condition 3 
          d -= e;
          c += ADD; //one addition 
        }
      else if ((d + e) % 2 == 0)
        { //condition 4 
          d = (d - e) / 2;
          c += ADD + DUP; //one addition, one duplicate 
        }
      //now d+e is odd 
      else if (d % 2 == 0)
        { //condition 5 
          d /= 2;
          c += ADD + DUP; //one addition, one duplicate 
        }
      //now d is odd and e is even 
      else if (d % 3 == 0)
        { //condition 6 
          d = d / 3 - e;
          c += 3.0 * ADD + DUP; //three additions, one duplicate 
        }
      else if ((d + e) % 3 == 0)
        { //condition 7 
          d = (d - 2 * e) / 3;
          c += 3.0 * ADD + DUP; //three additions, one duplicate 
        }
      else if ((d - e) % 3 == 0)
        { //condition 8 
          d = (d - e) / 3;
          c += 3.0 * ADD + DUP; //three additions, one duplicate 
        }
      else //necessarily e is even: catches all cases 
        { //condition 9 
          e /= 2;
          c += ADD + DUP; //one addition, one duplicate 
        }
    }
  
  return c;
}
