#include "main.h"
#include "utils.h"
#include "cudautils.h"

int main (int argc, char * argv[]) 
{
	unsigned int i;
#ifdef TEST
	unsigned int j;
#endif
	unsigned int B1;
	unsigned int number_of_curves;
	unsigned int nbfactor=0;

	clock2_t time_gpu;
	int device;

	char usage[]="./gpu_ecm N B1 [-s firstsigma] [-n number_of_curves] [-d device]";

	biguint_t *h_xarray;
	biguint_t *h_zarray;
	biguint_t h_invmod;
	biguint_t h_N;
	biguint_t h_N2;
	biguint_t *h_x2array;
	biguint_t *h_z2array;
 
  mpz_t N;
  mpz_t N2;
	mpz_t mpz_max; //2^(32*SIZE_NUMBER)
	mpz_t mpz_invmod; //N^-1 mod (2^(32*SIZE_NUMBER))
	mpz_t mpz_Rinv; // 2^-(32*SIZE_NUMBER) mod N
  mpz_t sigma;
  mpz_t mpz_B1;
  mpz_t mpz_d;
  mpz_t xp;
  mpz_t zp;
	mpz_t x2p;
	mpz_t z2p;
  mpz_t xfin;
  mpz_t zfin;

	gmp_randstate_t state;
	gmp_randinit_default(state);
	unsigned long seed=0;

  mpz_init (N);
  mpz_init (N2);
  mpz_init (mpz_max);
  mpz_init (sigma);
  mpz_init (mpz_B1);
  mpz_init (mpz_d);
  mpz_init_set_ui (xp,2);
  mpz_init_set_ui (zp,1);
	mpz_init_set_ui (x2p,9);
	mpz_init (z2p);
	mpz_init (xfin);
  mpz_init (zfin);
  mpz_init (mpz_invmod);
  mpz_init (mpz_Rinv);

  if (argc < 3)
  {  
		printf ("Error in call function : not enough arguments.\n%s\n",usage);
		exit(1);
	}
	
  mpz_set_str (N, argv[1], 10); // in base 10 
  mpz_set_str (mpz_B1, argv[2], 10);

	argc-=2;
	argv+=2;

	//default values
	number_of_curves=0;
	device=-1;
	mpz_set_ui(sigma,0);
	
	while (argc > 2 && argv[1][0]=='-')
	{
		if (argv[1][1] == 's')
		{
			mpz_set_str(sigma, argv[2], 10);
			argc-=2;
			argv+=2;
		}
		else if (argv[1][1] == 'n')
		{
			sscanf(argv[2],"%u",&number_of_curves);
			argc-=2;
			argv+=2;
		}
		else if (argv[1][1] == 'd')
		{
			sscanf(argv[2],"%d",&device);
			argc-=2;
			argv+=2;
		}
	}
	
	if (argc!=1)
	{
			printf ("Error in call function : wrong number of arguments.\n%s\n",usage);
			exit(1);
	}


	//init data according to the arguments
	if (mpz_cmp_ui(sigma,0)==0)
	{
		seed=time(NULL);
		gmp_randseed_ui(state,seed);
  	mpz_urandomb(sigma,state,3);//between 0 and 2^3
		mpz_add_ui(sigma,sigma,6);//add 6
	}
	gmp_randclear (state);
	
  // check N is odd 
  if (mpz_divisible_ui_p (N, 2))
  {
  	fprintf (stderr, "Error, N should be odd\n");
   	exit (1);
  }

	mpz_ui_pow_ui(mpz_max,2,32*SIZE_NUMBER-1);	
	if (mpz_cmp(N,mpz_max) >=0)
  {
  	fprintf (stderr, "Error, N should be stricly lower than 2^%d\n",32*SIZE_NUMBER-1);
   	exit (1);
  }
	mpz_ui_pow_ui(mpz_max,2,32*SIZE_NUMBER);	

	if (mpz_cmp_ui(mpz_B1,TWO32) >=0)
  {
  	fprintf (stderr, "Error, B1 should be stricly lower than 2^%d\n",32);
   	exit (1);
  }
	B1=(unsigned int)mpz_get_ui(mpz_B1);

	cudaDeviceProp deviceProp;
				
	printf("#Compiled for a NVIDIA GPU with compute capability %d.%d.\n",MAJOR,MINOR);
	if (device!=-1)
	{
		printf("#Device %d is required.\n",device);
		cudaError_t err= cudaSetDevice(device);
		if (err != cudaSuccess)
		{
			printf("#Error : Could not use device %d\n",device);
			exit(1);
		}
	}
	
	cudaGetDevice(&device);
	cudaGetDeviceProperties(&deviceProp,device);
	if (deviceProp.major < MAJOR || (deviceProp.major== MAJOR && deviceProp.minor< MINOR))
	{
		printf("#Error : Device %d have a compute capability of %d.%d (required %d.%d).\n",device,deviceProp.major,deviceProp.minor,MAJOR,MINOR);
		exit(1);
	}
	else if (deviceProp.major==MAJOR && deviceProp.minor==MINOR)
		printf("#Will use device %d : %s, compute capability %d.%d, %d MPs.\n",device,deviceProp.name,deviceProp.major,deviceProp.minor,deviceProp.multiProcessorCount);
	else
		printf("#Will use device %d : %s, compute capability %d.%d (you should compile the program for this compute capability to be more efficient), %d MPs.\n",device,deviceProp.name,deviceProp.major,deviceProp.minor,deviceProp.multiProcessorCount);

	if (number_of_curves==0)
		number_of_curves=deviceProp.multiProcessorCount*CURVES_BY_BLOCK;

	gmp_printf ("#gpu_ecm launched with :\nN=%Zd\nB1=%u\ncurves=%u\nfirstsigma=%Zd\n",N,B1,number_of_curves,sigma);
	if (seed!=0)
		printf("#used seed %lu to generate sigma\n",seed);

	h_xarray=(biguint_t *) malloc(number_of_curves*sizeof(biguint_t));
	h_zarray=(biguint_t *) malloc(number_of_curves*sizeof(biguint_t));
	h_x2array=(biguint_t *) malloc(number_of_curves*sizeof(biguint_t));
	h_z2array=(biguint_t *) malloc(number_of_curves*sizeof(biguint_t));

	//Some precomputation
	//Compute N^-1 mod 2^(32*SIZE_NUMBER)
	mpz_invert(mpz_invmod,N,mpz_max);
	mpz_sub(mpz_invmod,mpz_max,mpz_invmod);
	//Compute  2^-(32*SIZE_NUMBER) mod N
	mpz_invert(mpz_Rinv,mpz_max,N);
	
	mpz_to_biguint(h_N,N);	
	mpz_to_biguint(h_invmod,mpz_invmod);	
	mpz_sub_ui(N2,N,1);
	mpz_divexact_ui(N2,N2,2);
	mpz_to_biguint(h_N2,N2);	
	
	mpz_set_ui(mpz_d,TWO32);
	mpz_invert(mpz_d,mpz_d,N);
	
	//Compute the Montgomery representation of x0, z0 
	mpz_mul_2exp(xp,xp,32*SIZE_NUMBER);
	mpz_mod(xp,xp,N);
	
	mpz_mul_2exp(zp,zp,32*SIZE_NUMBER);
	mpz_mod(zp,zp,N);
	
	for(i=0;i<number_of_curves;i++)
	{
		mpz_set_ui(x2p,9);

		mpz_mul(z2p,mpz_d,sigma);
		mpz_mod(z2p,z2p,N);
		mpz_mul_ui(z2p,z2p,64);
		mpz_add_ui(z2p,z2p,8);
		mpz_mod(z2p,z2p,N);

		mpz_to_biguint(h_xarray[i],xp);	

		mpz_to_biguint(h_zarray[i],zp);	

		mpz_mul_2exp(x2p,x2p,32*SIZE_NUMBER);
		mpz_mod(x2p,x2p,N);
		mpz_to_biguint(h_x2array[i],x2p);	
		mpz_mul_2exp(z2p,z2p,32*SIZE_NUMBER);
		mpz_mod(z2p,z2p,N);
		mpz_to_biguint(h_z2array[i],z2p);	

		mpz_add_ui(sigma,sigma,1);
	}	

	mpz_sub_ui(sigma,sigma,number_of_curves);


	printf("\n#Begin GPU computation...\n");
	time_gpu=cuda_Main(h_N,h_N2,h_invmod,h_xarray,h_zarray,h_x2array,h_z2array,B1,17,number_of_curves);
	printf("#All kernels finished, analysing results...\n");

	for(i=0;i<number_of_curves;i++)
	{
	//printf("Test\n");
		biguint_to_mpz(xfin,h_xarray[i]);	
	//printf("Test\n");
		mpz_mul(xfin,xfin,mpz_Rinv);
		mpz_mod(xfin,xfin,N);

		biguint_to_mpz(zfin,h_zarray[i]);	
		mpz_mul(zfin,zfin,mpz_Rinv);
		mpz_mod(zfin,zfin,N);

#ifndef TEST
		if (i==0 || i==number_of_curves-1)
		{
			//biguint_print(h_zarray[i]);
			printf("\n");
			gmp_printf("#Looking for factors for the curves with sigma=%Zd\n",sigma);
			/*printf("x=");
			biguint_print(h_xarray[i]);
			printf("\nz=");
			biguint_print(h_zarray[i]);
			printf("\nd=");
			biguint_print(h_darray[i]);
			printf("\n");*/
			nbfactor+=findfactor(N,xfin,zfin);
		}

		mpz_add_ui(sigma,sigma,1);
#endif

#ifdef TEST
		if (i==0)
		{
			printf ("res=");
			biguint_print(h_xarray[i]);
			printf ("\ncy=");
			//printf ("+(");
			biguint_print(h_zarray[i]);
			printf ("\n");
			//printf (")*2^1024\n");
			//mpz_mul(mpztemp,mpztemp,mpztemp);
			//gmp_printf("b=%Zd\n",mpztemp);
			//printf("print A;print B;diff=(A+B)*(B-A)-res; print diff.digits(2^32); print 2^(-1024)*(A+B)*(B-A) %% N - res\n");
		}
#endif
	}

	if (nbfactor==0)
		printf("#Results : No factor found\n");
	else if (nbfactor==1)
		printf("#Results : 1 factor found\n");
	else
		printf("#Results : %u curves find a factor (not necessarily different)\n",nbfactor);
	
	printf("\n#Temps gpu : %.3f init&copy=%.3f computation=%.3f\n",(double)(time_gpu.init+time_gpu.computation)/CLOCKS_PER_SEC,(double)(time_gpu.init)/CLOCKS_PER_SEC,(double)(time_gpu.computation)/CLOCKS_PER_SEC);

	mpz_clear (sigma);
  mpz_clear (N);
  mpz_clear (N2);
  mpz_clear (mpz_max);
  mpz_clear (mpz_invmod);
  mpz_clear (mpz_B1);
  mpz_clear (mpz_d);
  mpz_clear (xp);
  mpz_clear (zp);
  mpz_clear (x2p);
  mpz_clear (z2p);
  mpz_clear (xfin);
  mpz_clear (zfin);
  mpz_clear (mpz_Rinv);
  
	
	free((void *)h_xarray);
	free((void *)h_zarray);
	free(h_x2array);
	free(h_z2array);
	
  return 0;
}

