#ifdef _MSC_VER
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

#include "def.h"
#include "utils.h"

#define TWO32 4294967296 // 2^32 

 
extern int select_GPU (int, int, FILE*);
extern void cuda_Main (biguint_t, biguint_t, biguint_t*, biguint_t*, biguint_t*,
                       biguint_t*, mpz_t, unsigned int, unsigned int, FILE*, 
                       FILE*);

int main (int argc, char * argv[]) 
{
  int main_ret=EXIT_SUCCESS;

  long begincputime = 0, endcputime = 0, begingputime = 0, endgputime = 0; 

  unsigned int i;
  unsigned int B1;
  unsigned int number_of_curves = 0;

  int ret;
  unsigned int firstinvd;
  unsigned int firstinvd_arg = 0;

  FILE *savefile = NULL;
  char *savefilename = NULL;

  int device = -1;

  int verbose=0;
  FILE *OUTPUT_VERBOSE = NULL;
  FILE *OUTPUT_VVERBOSE = NULL;

  biguint_t *h_xarray;
  biguint_t *h_zarray;
  biguint_t h_invmod;
  biguint_t h_N;
  biguint_t *h_x2array;
  biguint_t *h_z2array;
 
  mpz_t N;
  mpz_t mpz_B; //2^(MAX_BITS)
  mpz_t mpz_invmod; //N^-1 mod (2^(MAX_BITS))
  mpz_t mpz_Rinv; // 2^-(MAX_BITS) mod N
  mpz_t mpz_d;
  mpz_t xp;
  mpz_t zp;
  mpz_t x2p;
  mpz_t z2p;
  mpz_t s;

#ifdef _MSC_VER
  if(!SetPriorityClass(GetCurrentProcess(), IDLE_PRIORITY_CLASS))
	  printf("Could not set process priority\n");
#endif

  /* check ecm is linked with a compatible library */
  if (mp_bits_per_limb != GMP_NUMB_BITS)
  {
    fprintf (stderr, "Error, mp_bits_per_limb and GMP_NUMB_BITS differ\n");
    fprintf (stderr, "Please check your LD_LIBRARY_PATH variable\n");
    exit (EXIT_FAILURE);
  }
  if (GMP_NUMB_BITS != 32 && GMP_NUMB_BITS != 64) 
  {
    fprintf (stderr, "Error, GMP_NUMB_BITS should be either 32 or 64.\n");
    exit (EXIT_FAILURE);
  }

  
  /***********************/
  /* Looking for options */
  /***********************/

  while (argc > 1 && argv[1][0]=='-')
  {
    if (strcmp(argv[1], "-h")== 0 || strcmp (argv[1], "--help") ==0)
    {
      usage();
      exit(EXIT_SUCCESS);
    }
    else if (strcmp (argv[1], "-v") == 0)
    {
      verbose=1;
      argc-=1;
      argv+=1;
    }
    else if (strcmp (argv[1], "-vv") == 0)
    {
      verbose=2;
      argc-=1;
      argv+=1;
    }
    else if ((argc > 2) && (strcmp(argv[1], "-s") == 0))
    {
      sscanf(argv[2], "%u", &firstinvd_arg);
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

  if (verbose < 2)
#ifdef _MSC_VER
    OUTPUT_VVERBOSE=fopen("NUL:","a");
#else
    OUTPUT_VVERBOSE=fopen("/dev/null","a");
#endif
  else
    OUTPUT_VVERBOSE=stdout;

  if (verbose == 1)
    OUTPUT_VERBOSE=stdout;
  else
    OUTPUT_VERBOSE = OUTPUT_VVERBOSE;

  /***********************************************/
  /* Select the GPU and analyse number_of_curves */
  /***********************************************/

  number_of_curves = select_GPU(device, number_of_curves, OUTPUT_VERBOSE);

  /*****************************/
  /* Initialize some variables */
  /*****************************/

  //If needed, open savefile
  if (savefilename!=NULL)
  {
    if ( strcmp (savefilename, "-") == 0 )
      savefile = stdout;
    else
    {
      savefile=fopen(savefilename, "a");
      if (savefile == NULL)
      {
        fprintf(stderr,"Error, could not open file: %s\n", savefilename);
        exit(EXIT_FAILURE);
      }
    }
  }

  if (firstinvd_arg !=0 && 
                (firstinvd_arg < 2 || firstinvd_arg > TWO32-number_of_curves))
  {
    fprintf(stderr,"Error, -s should be bigger than 2 and smaller than %lu.\n",
                                                       TWO32-number_of_curves);
    exit(EXIT_FAILURE);
  }

  mpz_init (N);
  mpz_init (mpz_B);
  mpz_init (mpz_d);
  mpz_init (xp);
  mpz_init (zp);
  mpz_init (x2p);
  mpz_init (z2p);
  mpz_init (mpz_invmod);
  mpz_init (mpz_Rinv);
  mpz_init (s);


  h_xarray=(biguint_t *) malloc(number_of_curves*sizeof(biguint_t));
  h_zarray=(biguint_t *) malloc(number_of_curves*sizeof(biguint_t));
  h_x2array=(biguint_t *) malloc(number_of_curves*sizeof(biguint_t));
  h_z2array=(biguint_t *) malloc(number_of_curves*sizeof(biguint_t));
  
  /****************************/
  /*Some common precomputation*/
  /****************************/

  mpz_ui_pow_ui(mpz_B, 2, MAX_BITS);  
  
  //precompute s from B1
  begincputime = cputime();
  compute_s(s, B1);
  endcputime = cputime();
  gmp_fprintf(OUTPUT_VERBOSE,"#s has %lu bits\n", 
                  mpz_sizeinbase(s,2));
  gmp_fprintf(stdout,"Precomputation of s took %.3fs\n", 
                  (double) (endcputime-begincputime)/1000);
  
  
  /***********************************************/
  /*Computation for each input number in the file*/
  /***********************************************/

  while (!feof(stdin) && read_number(N, stdin)==1)
  {
    begincputime=cputime();
    gmp_fprintf (stdout, "Input number is %Zd (%u digits)\n", 
                                          N, mpz_sizeinbase(N, 10));

    if (mpz_sizeinbase(N, 2) > MAX_BITS-1)
    {
      fprintf (stderr, "Error, input number should be stricly lower than 2^%d\n",
                                                          MAX_BITS-1);
      main_ret= EXIT_FAILURE;
      goto free_memory_and_exit;
    }

    if (mpz_cmp_ui (N, 0) <= 0)
    {
      fprintf (stderr, "Error, input number should be positive\n");
      main_ret= EXIT_FAILURE;
      goto free_memory_and_exit;
    }
    else if (mpz_cmp_ui (N, 1) == 0)
    {
      fprintf (stdout, "********** Factor found in step 1: 1\n");
      goto end_main_loop;
    }
    else if (mpz_divisible_ui_p (N, 2))
    {
      fprintf (stdout, "********** Factor found in step 1: 2\n");
      goto end_main_loop;
      //FIXME cofactor
    }
  
    
    if (firstinvd_arg == 0)
      firstinvd = (get_random_ui() % (TWO32-2-number_of_curves)) + 2;
    else
      firstinvd = firstinvd_arg;
 
    gmp_fprintf (stdout, "Using B1=%u, firstinvd=%u, with %u curves\n", 
                                          B1, firstinvd, number_of_curves);
  
    //Some precomputation
        
    //Compute N^-1 mod 2^(MAX_BITS)
    mpz_invert(mpz_invmod, N, mpz_B);
    mpz_sub(mpz_invmod, mpz_B, mpz_invmod);
    //Compute  2^-(MAX_BITS) mod N
    mpz_invert(mpz_Rinv, mpz_B, N);
    
    mpz_to_biguint(h_N, N); 
    mpz_to_biguint(h_invmod, mpz_invmod); 
#ifdef _MSC_VER
   mpz_set_ux(mpz_d, TWO32);
#else
    mpz_set_ui(mpz_d, TWO32);
#endif
    mpz_invert(mpz_d, mpz_d, N);
    
    //Compute the Montgomery representation
    mpz_set_ui (xp, 2);
    mpz_mul_2exp(xp, xp, MAX_BITS);
    mpz_mod(xp, xp, N);
    
    mpz_set_ui (zp, 1);
    mpz_mul_2exp(zp, zp, MAX_BITS);
    mpz_mod(zp, zp, N);
    
    mpz_set_ui (x2p, 9);
    mpz_mul_2exp(x2p, x2p, MAX_BITS);
    mpz_mod(x2p, x2p, N);
  
    for (i=0; i<number_of_curves; i++)
    {
      mpz_mul_ui(z2p,mpz_d,firstinvd);
      mpz_mod(z2p,z2p,N);
      mpz_mul_ui(z2p,z2p,64);
      mpz_add_ui(z2p,z2p,8);
      mpz_mod(z2p,z2p,N);
  
      if (i==0 || i==number_of_curves-1)
        gmp_fprintf(OUTPUT_VVERBOSE,"8+64*d=%Zd\n",z2p);
  
      mpz_mul_2exp(z2p, z2p, MAX_BITS);
      mpz_mod(z2p, z2p, N);
  
      mpz_to_biguint(h_xarray[i],xp); 
      mpz_to_biguint(h_zarray[i],zp); 
      mpz_to_biguint(h_x2array[i],x2p); 
      mpz_to_biguint(h_z2array[i],z2p); 
  
      firstinvd++;
    } 
  
    firstinvd-=number_of_curves;
  
    begingputime=cputime();
    fprintf(OUTPUT_VERBOSE,"#Begin GPU computation...\n");
    cuda_Main( h_N, h_invmod, h_xarray, h_zarray, h_x2array, 
                          h_z2array, s, firstinvd, number_of_curves,
                          OUTPUT_VERBOSE, OUTPUT_VVERBOSE);
    endgputime=cputime();
  
    for(i=0;i<number_of_curves;i++)
    {
      biguint_to_mpz(xp,h_xarray[i]); 
      mpz_mul(xp,xp,mpz_Rinv);
      mpz_mod(xp,xp,N);
  
      biguint_to_mpz(zp,h_zarray[i]); 
      mpz_mul(zp,zp,mpz_Rinv);
      mpz_mod(zp,zp,N);
  
      if (i==0 || i==number_of_curves-1)
      {
        fprintf(OUTPUT_VVERBOSE,"\n");
        fprintf(OUTPUT_VVERBOSE,
         "#Looking for factors for the curves with (d*2^32) mod N = %u\n",
         firstinvd);
        gmp_fprintf(OUTPUT_VVERBOSE,"  xfin=%Zd\n  zfin=%Zd\n",xp,zp);
      }
      
      ret=findfactor(N,xp,zp);
      if (ret==ECM_NO_FACTOR_FOUND && savefile != NULL)
        write_resumefile_line (savefile, N, B1, xp, firstinvd, mpz_d);
      else if (ret==ECM_FACTOR_FOUND)
        fprintf(stdout,"Factor found with (d*2^32) mod N = %u\n",firstinvd);
            //Maybe print A for GMP-ECM
          
      if (i==0 || i==number_of_curves-1)
      {
        gmp_fprintf(OUTPUT_VVERBOSE,"  xunif=%Zd\n",xp);
      }
  
      firstinvd++;
    }
  
end_main_loop:
    endcputime=cputime();
    fprintf(stdout, "gpu_ecm took : %.3fs (%.3f+%.3f+%.3f)\n",
        (double) (endcputime-begincputime)/1000, 
        (double) (begingputime-begincputime)/1000, 
        (double) (endgputime-begingputime)/1000, 
        (double) (endcputime-endgputime)/1000);
    fprintf(stdout, "Throughput : %.3f\n", 
        1000 * (double)(number_of_curves)/(double)(endcputime-begincputime));
    }

free_memory_and_exit:
  mpz_clear (N);
  mpz_clear (mpz_B);
  mpz_clear (mpz_invmod);
  mpz_clear (mpz_d);
  mpz_clear (xp);
  mpz_clear (zp);
  mpz_clear (x2p);
  mpz_clear (z2p);
  mpz_clear (mpz_Rinv);
  mpz_clear (s);

  
  free((void *) h_xarray);
  free((void *) h_zarray);
  free((void *) h_x2array);
  free((void *) h_z2array);
  
  if (savefile != NULL)
    fclose(savefile);

  return main_ret;
}
