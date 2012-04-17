#include "ecm-gpu.h"

#ifdef WITH_GPU
#include "cuda.h"

#define TWO32 4294967296 /* 2^32 */ 

extern void select_and_init_GPU (int, unsigned int*, int);
extern float cuda_Main (biguint_t, biguint_t, biguint_t, digit_t, biguint_t*, 
                        biguint_t*, biguint_t*, biguint_t*, mpz_t, unsigned int, 
                        unsigned int, FILE*, FILE*);

void print_factor_cofactor (mpz_t N, mpz_t factor)
{
    gmp_fprintf(stdout,"********** Factor found in step 1: %Zd\n", factor);
    if (mpz_cmp (factor, N) ==0 )
      fprintf(stdout, "Found input number N\n");
    else
    {
      mpz_t cofactor;
      mpz_init (cofactor);
    
      mpz_divexact(cofactor, N, factor);
      
      gmp_fprintf(stdout,"Found %s factor of %u digits: %Zd\n",
              mpz_probab_prime_p (factor, 5) ? "probable prime" : "composite", 
              mpz_sizeinbase (factor, 10), factor);
     
      gmp_fprintf(stdout,"%s cofactor %Zd has %u digits\n",
            mpz_probab_prime_p (cofactor, 5) ? "Probable prime" : "Composite", 
            cofactor, mpz_sizeinbase (cofactor, 10));

      mpz_clear(cofactor);
    }
}

unsigned int findfactor(mpz_t N, mpz_t xfin, mpz_t zfin)
{
  mpz_t gcd;
  mpz_init (gcd);

  mpz_gcd(gcd, zfin, N);
  
  if (mpz_cmp_ui (gcd, 1)==0)
  {
    mpz_invert (zfin, zfin, N);
    mpz_mul (xfin, xfin, zfin);
    mpz_mod (xfin, xfin, N);
      
    mpz_clear(gcd);
    return ECM_NO_FACTOR_FOUND;
  }
  else //gcd !=1 (and gcd>0 because N>0) so we found a factor
  {
    print_factor_cofactor (N, gcd);
    mpz_clear(gcd);
    return ECM_FACTOR_FOUND_STEP1;
  }
}

void to_mont_repr (mpz_t x, mpz_t n)
{
  mpz_mul_2exp (x, x, ECM_GPU_MAX_BITS);
  mpz_mod(x, x, n);
}

void from_mont_repr (mpz_t x, mpz_t n, mpz_t invB)
{
  mpz_mul(x, x, invB);
  mpz_mod(x, x, n);
}

void mpz_to_biguint (biguint_t a, mpz_t b)
{
  int i;

  for (i=0;i<NB_DIGITS;i++)
  {
#if GMP_NUMB_BITS == 32
    a[i]=mpz_getlimbn(b, i);  
#else // GMP_NUMB_BITS == 64
    if (i%2 == 0)
      a[i]=(mpz_getlimbn(b, i/2) & 0x00000000ffffffff);
    else
      a[i]=(mpz_getlimbn(b, i/2) >> 32);  
#endif
  }
}

void biguint_to_mpz (mpz_t a, biguint_t b)
{
  int i;
  
  mpz_set_ui(a, 0);

  for (i=NB_DIGITS-1;i>=0;i--)
  {
    mpz_mul_2exp(a, a, 32);
	  mpz_add_ui(a , a, b[i]);
  }
}

int gpu_ecm_stage1 (mpz_t N, mpz_t s, unsigned int number_of_curves, 
                    unsigned int firstsigma, float *gputime)
{
  int fct_ret = ECM_NO_FACTOR_FOUND;
  int ret;

  unsigned int sigma;
  unsigned int i;

  mpz_t N3; /* N3 = 3*N */
  mpz_t w; /* w = 2^(SIZE_DIGIT) */
  mpz_t invN; /* invN = -N^-1 mod w */
  mpz_t invB; /* invB = 2^(-MAX_BITS) mod N ; B is w^NB_DIGITS */
  mpz_t invw; /* w^(-1) mod N */
  mpz_t M; /* (invN*N+1)/w */
  mpz_t xp, zp, x2p, z2p;

  /* The same variables but for the GPU */
  biguint_t *h_xarray, *h_zarray, *h_x2array, *h_z2array;
  digit_t h_invN;
  biguint_t h_N, h_3N, h_M;

  /*****************************/
  /* Initialize some variables */
  /*****************************/
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

  h_xarray= (biguint_t *) malloc (number_of_curves * sizeof (biguint_t));
  h_zarray= (biguint_t *) malloc (number_of_curves * sizeof (biguint_t));
  h_x2array= (biguint_t *) malloc (number_of_curves * sizeof (biguint_t));
  h_z2array= (biguint_t *) malloc (number_of_curves * sizeof (biguint_t));

  /*Some computation depending on N */
  mpz_mul_ui (N3, N, 3); /* Compute N3 = 3*N */
  mpz_ui_pow_ui (w, 2, SIZE_DIGIT); /* Compute w = 2^SIZE_DIGIT */
    
  mpz_invert (invN, N, w);
  mpz_sub (invN, w, invN); /* Compute invN = -N^-1 mod w */
   
  mpz_mul (M, invN, N);
  mpz_add_ui (M, M, 1);
  mpz_divexact (M, M, w); /* Compute M = (invN*N+1)/w */

  mpz_to_biguint (h_N, N); 
  mpz_to_biguint (h_3N, N3); 
  mpz_to_biguint (h_M, M); 
  h_invN = mpz_get_ui (invN); 
   
  mpz_ui_pow_ui (invB, 2, ECM_GPU_MAX_BITS); 
  mpz_invert (invB, invB, N); /* Compute invB = 2^(-MAX_BITS) mod N */

  mpz_invert (invw, w, N); /* Compute inw = 2^-SIZE_DIGIT % N */

  /* xp zp x2p are independent of N and the curve */
  mpz_set_ui (xp, 2);
  mpz_set_ui (zp, 1);
  mpz_set_ui (x2p, 9);

  /* Compute their Montgomery representation */
  to_mont_repr (xp, N);
  to_mont_repr (zp, N);
  to_mont_repr (x2p, N);
  
  /* for each curve, compute z2p and put xp, zp, x2p, z2p in the h_*array  */
  for (sigma = firstsigma; sigma < firstsigma+number_of_curves; sigma++)
  {
    i = sigma - firstsigma;

    mpz_mul_ui (z2p, invw, sigma);
    mpz_mod (z2p, z2p, N);
    mpz_mul_2exp (z2p, z2p, 6);
    mpz_add_ui (z2p, z2p, 8);
    mpz_mod (z2p, z2p, N); /* z2p = 8+64*d */

    to_mont_repr (z2p, N);

    mpz_to_biguint(h_xarray[i], xp); 
    mpz_to_biguint(h_zarray[i], zp); 
    mpz_to_biguint(h_x2array[i], x2p); 
    mpz_to_biguint(h_z2array[i], z2p); 
  } 
 
  /* Call the wrapper function that call the GPU */
  *gputime=cuda_Main (h_N, h_3N, h_M, h_invN, h_xarray, h_zarray, h_x2array, 
                     h_z2array, s, firstsigma, number_of_curves,
                     stdout, stdout);

  /* Analyse results */
  for (sigma = firstsigma; sigma < firstsigma+number_of_curves; sigma++)
  {
    i = sigma - firstsigma;

    biguint_to_mpz(xp, h_xarray[i]); 
    biguint_to_mpz(zp, h_zarray[i]); 
    
    from_mont_repr (xp, N, invB);
    from_mont_repr (zp, N, invB);
  
    ret = findfactor (N, xp, zp);

    if (ret != ECM_NO_FACTOR_FOUND)
      fct_ret = ret;
    //if (ret==ECM_NO_FACTOR_FOUND && savefilename != NULL)
    //  write_resumefile_wrapper (savefilename, &n, B1, xp, invd, invw);
    //else if (ret==ECM_FACTOR_FOUND)
      /* Maybe print A (or NU now) for GMP-ECM */
    //  fprintf(OUTPUT_STD_VERBOSE, "Factor found with (d*2^32) mod N = %u\n", invd);
          
    }
  
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
  
  free ((void *) h_xarray);
  free ((void *) h_zarray);
  free ((void *) h_x2array);
  free ((void *) h_z2array);

  return fct_ret;
}
#endif

/*
int gpu_ecm (ATTRIBUTE_UNUSED mpz_t f, ATTRIBUTE_UNUSED mpz_t N, 
             ATTRIBUTE_UNUSED mpz_t s, ATTRIBUTE_UNUSED double B1,
             ATTRIBUTE_UNUSED int device, ATTRIBUTE_UNUSED int *device_init, 
             ATTRIBUTE_UNUSED unsigned int *nb_curves, 
             ATTRIBUTE_UNUSED unsigned int firstinvd)
*/
int
gpu_ecm (mpz_t f, int param, mpz_t firstsigma, mpz_t n, mpz_t go, 
         double *B1done, double B1, mpz_t B2min_parm, mpz_t B2_parm, 
         double B2scale, unsigned long k, const int S, int verbose, int repr,
         int nobase2step2, int use_ntt, int sigma_is_A, FILE *os, FILE* es, 
         char *chkfilename, char *TreeFilename, double maxmem, 
         double stage1time, gmp_randstate_t rng, int (*stop_asap)(void), 
         mpz_t batch_s, double *batch_last_B1_used, int device, 
         int *device_init, unsigned int *nb_curves)

#ifndef WITH_GPU
{
  fprintf(stderr, "This version of libecm does not contain the GPU code.\n"
                  "You should recompile it with ./configure --enable-gpu or\n"
                  "link a version of libecm which contain the GPU code.\n");
  exit(EXIT_FAILURE);
}
#else
{
  int main_ret = ECM_NO_FACTOR_FOUND;
  long st;
  float gputime = 0.0;

  ASSERT((-1 <= sigma_is_A) && (sigma_is_A <= 1));
  ASSERT((GMP_NUMB_BITS == 32) || (GMP_NUMB_BITS == 64));

  set_verbose (verbose);
  ECM_STDOUT = (os == NULL) ? stdout : os;
  ECM_STDERR = (es == NULL) ? stdout : es;


  /* Check that N is not too big */
  if (mpz_sizeinbase (n, 2) > ECM_GPU_MAX_BITS-6)
    {
      outputf (OUTPUT_ERROR, "GPU: Error, input number should be stricly lower"
                             " than 2^%d\n", ECM_GPU_MAX_BITS-6);
      return ECM_ERROR;
    }

  /* Only param = ECM_PARAM_BATCH_SMALL_D is accepted on GPU */
  if (param == ECM_PARAM_DEFAULT)
      param = ECM_PARAM_BATCH_SMALL_D;
    
  if (param != ECM_PARAM_BATCH_SMALL_D)
    {
      outputf (OUTPUT_ERROR, "GPU: Error, only param = ECM_PARAM_BATCH_SMALL_D "
                             "is accepted on GPU.\n");
      return ECM_ERROR;
    }

  /* check that repr == ECM_MOD_DEFAULT */
  if (repr != ECM_MOD_DEFAULT)
    {
      outputf (OUTPUT_ERROR, "GPU: Error, only repr = ECM_MOD_DEFAULT "
                             "is accepted on GPU.\n");
      return ECM_ERROR;
    }
 
  /* Cannot do resume on GPU */
  if (!ECM_IS_DEFAULT_B1_DONE(*B1done) && *B1done < B1)
    {
      outputf (OUTPUT_ERROR, "GPU: Error, cannot resume on GPU.\n");
      return ECM_ERROR;
    }

  /* Compute s */
  if (B1 != *batch_last_B1_used || mpz_cmp_ui (batch_s, 1) <= 0)
    {
      *batch_last_B1_used = B1;

      st = cputime ();
      /* construct the batch exponent */
      compute_s (batch_s, B1);
      outputf (OUTPUT_VERBOSE, "GPU: computing prime product of %zu bits took " 
                   "%ldms\n", mpz_sizeinbase (batch_s, 2), cputime () - st);
    }

  /* Initialiaze the GPU if necessary */
  if (!*device_init)
    {
      st = cputime ();
      if (test_verbose (OUTPUT_VERBOSE))
          select_and_init_GPU (device, nb_curves, 1);
      else
          select_and_init_GPU (device, nb_curves, 0);

      outputf (OUTPUT_VERBOSE, "GPU: Selection and initialization of the device "
                               "took %ldms\n", elltime (st, cputime ()));
      /* TRICKS: If initialization of the device is too long (few seconds), */
      /* try running 'nvidia-smi -q -l' on the background .                 */
      *device_init = 1;
    }

  /* */
  if (sigma_is_A == -1)
    {
      /* Cant do stage 2 on gpu*/
      outputf (OUTPUT_ERROR, "GPU: Error, cannot do stage 2 on GPU.\n");
      return ECM_ERROR;
    }
  else if (sigma_is_A == 1)
    {
      /*compute sigma from A*/
      outputf (OUTPUT_ERROR, "GPU: Not yet implemented.\n");
      return ECM_ERROR;
    }

  if (sigma_is_A == 0 && mpz_sgn(firstsigma) == 0)
    {
      /*generate random one*/
      mpz_set_ui (firstsigma, (get_random_ui() % (TWO32-2-*nb_curves)) + 2 );    
    }
  else /* sigma should be in [2, 2^32-nb_curves] */
    {
      if (mpz_cmp_ui (firstsigma, 2) < 0 || 
          mpz_cmp_ui (firstsigma, TWO32-*nb_curves) > 0)
        {
          outputf (OUTPUT_ERROR, "GPU: Error, sigma should be bigger than 2 "
                                 "and smaller than %lu.\n", TWO32-*nb_curves);
          return ECM_ERROR;
        }
    }

  /* TODO: before beginning stage1:
            go should be NULL (or 1)
            print B1, firstinv, nb_curves (modify a little 
              print_B1_B2_poly) 
  print_B1_B2_poly (OUTPUT_NORMAL, ECM_ECM, B1, *B1done, 
		  mpz_t B2min_param, mpz_t B2min, mpz_t B2, int S, sigma,
		  sigma_is_A, go, param)
  In the meantime, we use this temporary printf:
  */
  gmp_fprintf (stdout, "Using B1=%1.0f, sigma=%Zd, with %u curves\n", 
                                          B1, firstsigma, *nb_curves);
  

  gpu_ecm_stage1 (n, batch_s, *nb_curves, mpz_get_ui(firstsigma), &gputime);

  fprintf (stdout, "time GPU: %.3fs\n", (gputime/1000));
  fprintf (stdout, "Throughput: %.3f\n", 1000 * (*nb_curves)/gputime);

  /* TODO: print CPU time for stage1 */
  /* TODO: Do stage2 (compute the parameters first (like in ecm.c) */
  /* TODO: Save in a file if requested */


  end_gpu_ecm:

  return main_ret;
}
#endif



