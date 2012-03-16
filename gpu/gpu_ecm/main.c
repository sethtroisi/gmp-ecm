#ifdef _MSC_VER
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

#include "def.h"
#include "utils.h"

#define TWO32 4294967296 // 2^32 

 
extern int select_and_init_GPU (int, int, FILE*);
extern void cuda_Main (biguint_t, biguint_t, biguint_t, digit_t, biguint_t*, 
                       biguint_t*, biguint_t*, biguint_t*, mpz_t, unsigned int, 
                       unsigned int, FILE*, FILE*);

int main (int argc, char * argv[]) 
{
  int main_ret=EXIT_SUCCESS;

  long begincputime = 0, endcputime = 0, begingputime = 0, endgputime = 0; 

  unsigned int i;
  unsigned int B1;
  unsigned int number_of_curves = 0;

  mpcandi_t n;

  int ret;
  unsigned int firstinvd;
  unsigned int firstinvd_arg = 0;
  unsigned int invd;

  char *savefilename = NULL;

  int device = -1;

  int verbose = 0;
  FILE *OUTPUT_STD_VERBOSE = NULL;
  FILE *OUTPUT_STD_VVERBOSE = NULL;

  mpz_t N3; /* N3 = 3*N */
  mpz_t w; /* w = 2^(SIZE_DIGIT) */
  mpz_t invN; /* invN = -N^-1 mod w */
  mpz_t invB; /* invB = 2^(-MAX_BITS) mod N ; B is w^NB_DIGITS */
  mpz_t invw; /* w^(-1) mod N */
  mpz_t M; /* (invN*N+1)/w */
  mpz_t xp;
  mpz_t zp;
  mpz_t x2p;
  mpz_t z2p;
  mpz_t s;

  /* The same variables but for the GPU */
  biguint_t *h_xarray;
  biguint_t *h_zarray;
  digit_t h_invN;
  biguint_t h_N;
  biguint_t h_3N;
  biguint_t h_M;
  biguint_t *h_x2array;
  biguint_t *h_z2array;
 

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
  
  sscanf(argv[1], "%u", &B1);

  if (verbose < 2)
#ifdef _MSC_VER
    OUTPUT_STD_VVERBOSE=fopen("NUL:","a");
#else
    OUTPUT_STD_VVERBOSE=fopen("/dev/null","a");
#endif
  else
    OUTPUT_STD_VVERBOSE=stdout;

  if (verbose == 1)
    OUTPUT_STD_VERBOSE=stdout;
  else
    OUTPUT_STD_VERBOSE = OUTPUT_STD_VVERBOSE;

  /***********************************************/
  /* Select the GPU and analyse number_of_curves */
  /***********************************************/

  begincputime = cputime ();
  number_of_curves = select_and_init_GPU (device, number_of_curves, 
                                                                OUTPUT_STD_VERBOSE);
  fprintf(stdout, "Selecting and initialization of the device s took %.3fs\n", 
                                       (double) (cputime ()-begincputime)/1000);
  /* TRICKS: If initialization of the device is too long (few seconds), */
  /* try running 'nvidia-smi -q -l' on the background .                 */

  /*****************************/
  /* Initialize some variables */
  /*****************************/

  if (firstinvd_arg !=0 && 
                (firstinvd_arg < 2 || firstinvd_arg > TWO32-number_of_curves))
  {
    fprintf(stderr,"Error, -s should be bigger than 2 and smaller than %lu.\n",
                                                       TWO32-number_of_curves);
    exit(EXIT_FAILURE);
  }

  mpcandi_t_init (&n);

  mpz_init (N3);
  mpz_init (invw);
  mpz_init (w);
  mpz_init (M);
  mpz_init (xp);
  mpz_init (zp);
  mpz_init (x2p);
  mpz_init (z2p);
  mpz_init (invN);
  mpz_init (invB);
  mpz_init (s);


  h_xarray=(biguint_t *) malloc(number_of_curves*sizeof(biguint_t));
  h_zarray=(biguint_t *) malloc(number_of_curves*sizeof(biguint_t));
  h_x2array=(biguint_t *) malloc(number_of_curves*sizeof(biguint_t));
  h_z2array=(biguint_t *) malloc(number_of_curves*sizeof(biguint_t));
  
  /****************************/
  /*Some shared precomputation*/
  /****************************/

  /* precompute s from B1 */
  begincputime = cputime ();
  compute_s (s, B1);

  fprintf(OUTPUT_STD_VERBOSE, "#s has %lu bits\n", mpz_sizeinbase (s, 2));
  fprintf(stdout, "Precomputation of s took %.3fs\n", 
                                       (double) (cputime ()-begincputime)/1000);
  
  
  /***********************************************/
  /*Computation for each input number in the file*/
  /***********************************************/

  while (read_number(&n, stdin, 0) == 1)
  {
    /* The integer N we try to factor is in n.n */
    begincputime=cputime();

    if (n.nexprlen == 0)
      gmp_fprintf (stdout, "Input number is %Zd (%u digits)\n", n.n, n.ndigits);
    else
      fprintf (stdout, "Input number is %s (%u digits)\n", n.cpExpr, n.ndigits);

    /* Check that N is not too big */
    if (mpz_sizeinbase(n.n, 2) > MAX_BITS-6)
    {
      fprintf (stderr, "Error, input number should be stricly lower than 2^%d\n",
                                                          MAX_BITS-6);
      main_ret= EXIT_FAILURE;
      goto free_memory_and_exit;
    }

    /* Check that N is positive and handle special cases (N=1 and N is even) */
    if (mpz_cmp_ui (n.n, 0) <= 0)
    {
      fprintf (stderr, "Error, input number should be positive\n");
      main_ret= EXIT_FAILURE;
      goto free_memory_and_exit;
    }
    else if (mpz_cmp_ui (n.n, 1) == 0)
    {
      fprintf (stdout, "********** Factor found in step 1: 1\n");
      goto end_main_loop;
    }
    else if (mpz_divisible_ui_p (n.n, 2))
    {
      mpz_t two;
      mpz_init_set_ui (two, 2);
      print_factor_cofactor (n, two);
      mpz_clear (two);
      goto end_main_loop;
    }
  
    if (firstinvd_arg == 0)
      firstinvd = (get_random_ui() % (TWO32-2-number_of_curves)) + 2;
    else
      firstinvd = firstinvd_arg;
 
    gmp_fprintf (stdout, "Using B1=%u, firstinvd=%u, with %u curves\n", 
                                          B1, firstinvd, number_of_curves);
  
    /*Some precomputation specific to each N */
    mpz_mul_ui (N3, n.n, 3); /* Compute N3 = 3*N */

    mpz_ui_pow_ui (w, 2, SIZE_DIGIT); /* Compute w = 2^SIZE_DIGIT */
    
    mpz_invert (invN, n.n, w);
    mpz_sub (invN, w, invN); /* Compute invN = -N^-1 mod w */
   
    mpz_mul (M, invN, n.n);
    mpz_add_ui (M, M, 1);
    mpz_divexact (M, M, w); /* Compute M = (invN*N+1)/w */

    mpz_to_biguint (h_N, n.n); 
    mpz_to_biguint (h_3N, N3); 
    mpz_to_biguint (h_M, M); 
    h_invN = mpz_get_ui (invN); 
   
    mpz_ui_pow_ui (invB, 2, MAX_BITS); 
    mpz_invert (invB, invB, n.n); /* Compute invB = 2^(-MAX_BITS) mod N */

    mpz_invert (invw, w, n.n); /* Compute inw = 2^-SIZE_DIGIT % N */
    
    /* xp zp x2p are independent of N and the curve */
    mpz_set_ui (xp, 2);
    mpz_set_ui (zp, 1);
    mpz_set_ui (x2p, 9);

    /* Compute their Montgomery representation */
    to_mont_repr (xp, n);
    to_mont_repr (zp, n);
    to_mont_repr (x2p, n);
 
    /* for each curve, compute z2p and put xp, zp, x2p, z2p in the h_*array  */
    for (invd = firstinvd; invd < firstinvd+number_of_curves; invd++)
    {
      i = invd - firstinvd;

      mpz_mul_ui (z2p, invw, invd);
      mpz_mod (z2p, z2p, n.n);
      mpz_mul_2exp (z2p, z2p, 6);
      mpz_add_ui (z2p, z2p, 8);
      mpz_mod (z2p, z2p, n.n); /* z2p = 8+64*d */
  
      if (i == 0 || i == number_of_curves-1)
        gmp_fprintf(OUTPUT_STD_VVERBOSE,"8+64*d=%Zd\n",z2p);
  
      to_mont_repr (z2p, n);

      mpz_to_biguint(h_xarray[i], xp); 
      mpz_to_biguint(h_zarray[i], zp); 
      mpz_to_biguint(h_x2array[i], x2p); 
      mpz_to_biguint(h_z2array[i], z2p); 
    } 
 
    /* Call the wrapper function that call the GPU */
    begingputime=cputime();
    fprintf(OUTPUT_STD_VERBOSE,"#Begin GPU computation...\n");
    cuda_Main( h_N, h_3N, h_M, h_invN, h_xarray, h_zarray, h_x2array, 
                          h_z2array, s, firstinvd, number_of_curves,
                          OUTPUT_STD_VERBOSE, OUTPUT_STD_VVERBOSE);
    endgputime=cputime();
 
    /* Analyse results */
    for (invd = firstinvd; invd < firstinvd+number_of_curves; invd++)
    {
      i = invd - firstinvd;

      biguint_to_mpz(xp, h_xarray[i]); 
      biguint_to_mpz(zp, h_zarray[i]); 
      
      from_mont_repr (xp, n, invB);
      from_mont_repr (zp, n, invB);
  
      if (i==0 || i==number_of_curves-1)
      {
        fprintf(OUTPUT_STD_VVERBOSE,"\n");
        fprintf(OUTPUT_STD_VVERBOSE,
         "#Looking for factors for the curves with (d*2^32) mod N = %u\n", invd);
        gmp_fprintf(OUTPUT_STD_VVERBOSE,"  xfin=%Zd\n  zfin=%Zd\n",xp,zp);
      }
      
      ret = findfactor (n, xp, zp);

      if (ret==ECM_NO_FACTOR_FOUND && savefilename != NULL)
        write_resumefile_wrapper (savefilename, &n, B1, xp, invd, invw);
      else if (ret==ECM_FACTOR_FOUND)
        //Maybe print A for GMP-ECM
        fprintf(OUTPUT_STD_VERBOSE, "Factor found with (d*2^32) mod N = %u\n", invd);
          
      if (i==0 || i==number_of_curves-1)
      {
        gmp_fprintf(OUTPUT_STD_VVERBOSE,"  xunif=%Zd\n",xp);
      }
    }
  
end_main_loop:
    endcputime=cputime();
    if (endgputime == 0)
      endgputime = endcputime;
    if (begingputime == 0)
      begingputime = endcputime;

    fprintf(stdout, "gpu_ecm took : %.3fs (%.3f+%.3f+%.3f)\n",
        (double) (endcputime-begincputime)/1000, 
        (double) (begingputime-begincputime)/1000, 
        (double) (endgputime-begingputime)/1000, 
        (double) (endcputime-endgputime)/1000);
    fprintf(stdout, "Throughput : %.3f\n", 
        1000 * (double)(number_of_curves)/(double)(endcputime-begincputime));
    }

free_memory_and_exit:
  mpcandi_t_free(&n);

  mpz_clear (N3);
  mpz_clear (invN);
  mpz_clear (invw);
  mpz_clear (w);
  mpz_clear (M);
  mpz_clear (xp);
  mpz_clear (zp);
  mpz_clear (x2p);
  mpz_clear (z2p);
  mpz_clear (invB);
  mpz_clear (s);

  
  free((void *) h_xarray);
  free((void *) h_zarray);
  free((void *) h_x2array);
  free((void *) h_z2array);
  
  return main_ret;
}

