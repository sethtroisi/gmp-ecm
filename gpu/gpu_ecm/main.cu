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
	//unsigned int nbfactor=0;
	int ret;
	unsigned int firstinvd;

	FILE *savefile = NULL;
	char *savefilename = NULL;

	clock2_t time_gpu;
	int device;


	biguint_t *h_xarray;
	biguint_t *h_zarray;
	biguint_t h_invmod;
	biguint_t h_N;
	biguint_t *h_x2array;
	biguint_t *h_z2array;
 
  mpz_t N;
	mpz_t mpz_max; //2^(32*SIZE_NUMBER)
	mpz_t mpz_invmod; //N^-1 mod (2^(32*SIZE_NUMBER))
	mpz_t mpz_Rinv; // 2^-(32*SIZE_NUMBER) mod N
  mpz_t mpz_d;
  mpz_t xp;
  mpz_t zp;
	mpz_t x2p;
	mpz_t z2p;

  mpz_init (N);
  mpz_init (mpz_max);
  //mpz_init (mpz_B1);
  mpz_init (mpz_d);
  mpz_init_set_ui (xp,2);
  mpz_init_set_ui (zp,1);
	mpz_init_set_ui (x2p,9);
	mpz_init (z2p);
  mpz_init (mpz_invmod);
  mpz_init (mpz_Rinv);

	//default values
	number_of_curves=0;
	device=-1;
	firstinvd=0;
	
	while (argc > 1 && argv[1][0]=='-')
	{
		if (strcmp(argv[1], "-h")== 0 || strcmp (argv[1], "--help") ==0)
		{
			usage();
			exit(EXIT_SUCCESS);
		}
		else if ((argc > 2) && (strcmp(argv[1], "-s") == 0))
		{
			sscanf(argv[2], "%u", &firstinvd);
			argc-=2;
			argv+=2;
		}
		else if ((argc > 2) && (strcmp(argv[1],"-n") == 0))
		{
			sscanf(argv[2], "%u", &number_of_curves);
			argc-=2;
			argv+=2;
		}
		else if ((argc > 2) && (strcmp(argv[1],"-d") == 0))
		{
			sscanf(argv[2], "%d", &device);
			argc-=2;
			argv+=2;
		}
		else if ((argc > 2) && (strcmp(argv[1],"-save") == 0))
		{
			savefilename=argv[2];
			argc-=2;
			argv+=2;
		}
		else
		{
			fprintf(stderr,"Unknow option: %s\n",argv[1]);
			exit(EXIT_FAILURE);
		}
	}
	
  if (argc != 2)
  {  
		fprintf(stderr,"Invalid arguments. See gpu_ecm --help.\n");
		exit(EXIT_FAILURE);
	}
	
  //mpz_set_str (N, argv[1], 10); // in base 10 
	sscanf(argv[1], "%u", &B1);

	//argc-=2;
	//argv+=2;

	cudaDeviceProp deviceProp;
				
	fprintf(stdout,"#Compiled for a NVIDIA GPU with compute capability %d.%d.\n",MAJOR,MINOR);
	if (device!=-1)
	{
		fprintf(stdout,"#Device %d is required.\n",device);
		cudaError_t err= cudaSetDevice(device);
		if (err != cudaSuccess)
		{
			fprintf(stdout,"#Error : Could not use device %d\n",device);
			exit(EXIT_FAILURE);
		}
	}
	
	cudaGetDevice(&device);
	cudaGetDeviceProperties(&deviceProp,device);
	if (deviceProp.major < MAJOR || (deviceProp.major== MAJOR && deviceProp.minor< MINOR))
	{
		fprintf(stdout,"#Error : Device %d have a compute capability of %d.%d (required %d.%d).\n",device,deviceProp.major,deviceProp.minor,MAJOR,MINOR);
		exit(EXIT_FAILURE);
	}
	else if (deviceProp.major==MAJOR && deviceProp.minor==MINOR)
		fprintf(stdout,"#Will use device %d : %s, compute capability %d.%d, %d MPs.\n",device,deviceProp.name,deviceProp.major,deviceProp.minor,deviceProp.multiProcessorCount);
	else
		fprintf(stdout,"#Will use device %d : %s, compute capability %d.%d (you should compile the program for this compute capability to be more efficient), %d MPs.\n",device,deviceProp.name,deviceProp.major,deviceProp.minor,deviceProp.multiProcessorCount);

	
	//number_of_curves should be a multiple of CURVES_BY_BLOCK
	number_of_curves=(number_of_curves/CURVES_BY_BLOCK)*CURVES_BY_BLOCK;
	if (number_of_curves==0)
		number_of_curves=deviceProp.multiProcessorCount*CURVES_BY_BLOCK;

//TODO faire une boucle pour boucler sur toutes les lignes de stdin
//while fscanf ==1 ??
	gmp_fscanf(stdin,"%Zd",N);
  
	if (mpz_divisible_ui_p (N, 2))
  {
		//TODO ne pas faire erreur mais sortir normale avec 2 comme facteur
		//faire une fonction pour uniformiser la sortie avec les autres
  	fprintf (stderr, "Error, N should be odd\n");
   	exit (EXIT_FAILURE);
  }

	mpz_ui_pow_ui(mpz_max,2,32*SIZE_NUMBER-1);	
	if (mpz_cmp(N,mpz_max) >=0)
  {
  	fprintf (stderr, "Error, N should be stricly lower than 2^%d\n",32*SIZE_NUMBER-1);
   	exit (EXIT_FAILURE);
  }
	
	if (firstinvd==0)
	{
		srand(time(NULL));
		//generer firstinvd aléatoirement dans les bonnes bornes[2,2^32-1-nbrecurves]
		firstinvd=rand() % TWO31 + 2;
	}

	gmp_fprintf (stderr,"#N=%Zd\n#B1=%u\n#firstinvd=%u\n%u ",
															N,B1,firstinvd,number_of_curves);

	h_xarray=(biguint_t *) malloc(number_of_curves*sizeof(biguint_t));
	h_zarray=(biguint_t *) malloc(number_of_curves*sizeof(biguint_t));
	h_x2array=(biguint_t *) malloc(number_of_curves*sizeof(biguint_t));
	h_z2array=(biguint_t *) malloc(number_of_curves*sizeof(biguint_t));

	//Some precomputation
	//Compute N^-1 mod 2^(32*SIZE_NUMBER)
	mpz_ui_pow_ui(mpz_max,2,32*SIZE_NUMBER);	
	mpz_invert(mpz_invmod,N,mpz_max);
	mpz_sub(mpz_invmod,mpz_max,mpz_invmod);
	//Compute  2^-(32*SIZE_NUMBER) mod N
	mpz_invert(mpz_Rinv,mpz_max,N);
	
	mpz_to_biguint(h_N,N);	
	mpz_to_biguint(h_invmod,mpz_invmod);	
	
	mpz_set_ui(mpz_d,TWO32);
	mpz_invert(mpz_d,mpz_d,N);
	
	//Compute the Montgomery representation
	mpz_mul_2exp(xp,xp,32*SIZE_NUMBER);
	mpz_mod(xp,xp,N);
	
	mpz_mul_2exp(zp,zp,32*SIZE_NUMBER);
	mpz_mod(zp,zp,N);
	
	mpz_mul_2exp(x2p,x2p,32*SIZE_NUMBER);
	mpz_mod(x2p,x2p,N);

	#ifdef TEST
	mpz_set_str(xp,"143805795474180477715551670989706028047729675645425252587319477048578939317229826117558190107889475963854697259734869847815638012718519302222375492556416739307266882805602763960478644417547693257368452083973266767326241158429824355507033182709885084613778430054616797020729738481697557927534650235438",10);
	mpz_set_ui(zp,1);
	mpz_mul(x2p,xp,zp);
	mpz_mod(x2p,x2p,N);
	gmp_printf("\n%Zd\n%Zd\n",xp,zp);

	mpz_mul_2exp(xp,xp,32*SIZE_NUMBER);
	mpz_mod(xp,xp,N);
	
	mpz_mul_2exp(zp,zp,32*SIZE_NUMBER);
	mpz_mod(zp,zp,N);
	#endif

	for(i=0;i<number_of_curves;i++)
	{
		mpz_mul_ui(z2p,mpz_d,firstinvd);
		mpz_mod(z2p,z2p,N);
		mpz_mul_ui(z2p,z2p,64);
		mpz_add_ui(z2p,z2p,8);
		mpz_mod(z2p,z2p,N);

		if (i==0 || i==number_of_curves-1)
			gmp_fprintf(stdout,"\n8+64*d=%Zd",z2p);

		mpz_mul_2exp(z2p,z2p,32*SIZE_NUMBER);
		mpz_mod(z2p,z2p,N);

		mpz_to_biguint(h_xarray[i],xp);	
		mpz_to_biguint(h_zarray[i],zp);	
		mpz_to_biguint(h_x2array[i],x2p);	
		mpz_to_biguint(h_z2array[i],z2p);	

		firstinvd++;
	}	

	firstinvd-=number_of_curves;

	fprintf(stdout,"\n#Begin GPU computation...\n");
	time_gpu=cuda_Main(h_N,h_invmod,h_xarray,h_zarray,h_x2array,h_z2array,B1,firstinvd,number_of_curves);
	fprintf(stdout,"#All kernels finished, analysing results...\n");

	if (savefilename!=NULL)
	{
		savefile=fopen(savefilename,"a");
		if (savefile == NULL)
		{
			fprintf(stderr,"Could not open file %s for writing\n", "tempname");
			exit(EXIT_FAILURE);
		}
	}

	for(i=0;i<number_of_curves;i++)
	{
		biguint_to_mpz(xp,h_xarray[i]);	
		mpz_mul(xp,xp,mpz_Rinv);
		mpz_mod(xp,xp,N);

		biguint_to_mpz(zp,h_zarray[i]);	
		mpz_mul(zp,zp,mpz_Rinv);
		mpz_mod(zp,zp,N);

#ifndef TEST
		if (i==0 || i==number_of_curves-1)
		{
			printf("\n");
			fprintf(stdout,
			 "#Looking for factors for the curves with (d*2^32) mod N = %u\n",
			 firstinvd);
			gmp_fprintf(stdout,"  xfin=%Zd\n  zfin=%Zd\n",xp,zp);
		}
		
		ret=findfactor(N,xp,zp);
		if (ret==ECM_NO_FACTOR_FOUND && savefile != NULL)
			write_resumefile_line (savefile, N, B1, xp, firstinvd, mpz_d);
		//else if (ret==ECM_FACTOR_FOUND)
			//print le A pour GMP-ECM
				
		if (i==0 || i==number_of_curves-1)
		{
			//temporaire juste pour matcher les anciennes sorties
			gmp_fprintf(stdout,"  xunif=%Zd\n",xp);
			fprintf(stdout,"  #No factors found. You shoud try with a bigger B1.\n");
		}

		firstinvd++;
#endif

#ifdef TEST
		if (i==0)
		{
			printf ("xp=");
			biguint_print(h_xarray[i]);
			printf ("\nzp=");
			//printf ("+(");
			biguint_print(h_zarray[i]);
			//printf ("\n");
			gmp_printf("\n%Zd\n%Zd\n",xp,zp);
			mpz_sub(xp,xp,x2p);
			gmp_printf("\n%Zd\n",xp);
			//printf (")*2^1024\n");
			//mpz_mul(mpztemp,mpztemp,mpztemp);
			//gmp_printf("b=%Zd\n",mpztemp);
			//printf("print A;print B;diff=(A+B)*(B-A)-res; print diff.digits(2^32); print 2^(-1024)*(A+B)*(B-A) %% N - res\n");
		}
#endif
	}

//temporaire juste pour matcher les anciennes sorties
	fprintf(stdout,"#Results : No factor found\n");
/*
	if (nbfactor==0)
		fprintf(stdout,"#Results : No factor found\n");
	else if (nbfactor==1)
		fprintf(stdout,"#Results : 1 factor found\n");
	else
		fprintf(stdout,"#Results : %u curves find a factor (not necessarily different)\n",nbfactor);
	*/

	fprintf(stdout,"\n#Temps gpu : %.3f init&copy=%.3f computation=%.3f\n",(double)(time_gpu.init+time_gpu.computation)/CLOCKS_PER_SEC,(double)(time_gpu.init)/CLOCKS_PER_SEC,(double)(time_gpu.computation)/CLOCKS_PER_SEC);

  mpz_clear (N);
  mpz_clear (mpz_max);
  mpz_clear (mpz_invmod);
  mpz_clear (mpz_d);
  mpz_clear (xp);
  mpz_clear (zp);
  mpz_clear (x2p);
  mpz_clear (z2p);
  mpz_clear (mpz_Rinv);
  
	
	free((void *)h_xarray);
	free((void *)h_zarray);
	free(h_x2array);
	free(h_z2array);
	
	if (savefile != NULL)
		fclose(savefile);

  return EXIT_SUCCESS;
}

