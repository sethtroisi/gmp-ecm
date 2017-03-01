#include "ecm-gpu.h"

#ifdef WITH_GPU

#define TWO32 4294967296 /* 2^32 */ 

extern int select_and_init_GPU (int, unsigned int*, int);
extern float cuda_Main (biguint_t, biguint_t, biguint_t, digit_t, biguint_t*, 
                        biguint_t*, biguint_t*, biguint_t*, mpz_t, unsigned int, 
                        unsigned int, int);

int findfactor (mpz_t factor, mpz_t N, mpz_t xfin, mpz_t zfin)
{
  int youpi;
  mpz_t gcd;
  mpz_init (gcd);

  mpz_gcd (gcd, zfin, N);
  
  if (mpz_cmp_ui (gcd, 1) == 0)
  {
    mpz_invert (zfin, zfin, N);
    mpz_mul (xfin, xfin, zfin);
    mpz_mod (xfin, xfin, N);
      
    mpz_set (factor, xfin);
    youpi = ECM_NO_FACTOR_FOUND;
  }
  else //gcd !=1 (and gcd>0 because N>0) so we found a factor
  {
      mpz_set (factor, gcd);
      youpi = ECM_FACTOR_FOUND_STEP1;
    }

  mpz_clear (gcd);
  return youpi;
}

void to_mont_repr (mpz_t x, mpz_t n)
{
  mpz_mul_2exp (x, x, ECM_GPU_MAX_BITS);
  mpz_mod (x, x, n);
}

void from_mont_repr (mpz_t x, mpz_t n, mpz_t invB)
{
  mpz_mul (x, x, invB);
  mpz_mod (x, x, n);
}

void mpz_to_biguint (biguint_t a, mpz_t b)
{
  int i;

  for (i=0;i<ECM_GPU_NB_DIGITS;i++)
  {
#if GMP_NUMB_BITS == 32
    a[i]=mpz_getlimbn (b, i);  
#else // GMP_NUMB_BITS == 64
    if (i%2 == 0)
      a[i]=(mpz_getlimbn (b, i/2) & 0x00000000ffffffff);
    else
      a[i]=(mpz_getlimbn (b, i/2) >> 32);  
#endif
  }
}

void biguint_to_mpz (mpz_t a, biguint_t b)
{
  int i;
  
  mpz_set_ui (a, 0);

  for (i=ECM_GPU_NB_DIGITS-1;i>=0;i--)
  {
    mpz_mul_2exp (a, a, 32);
	  mpz_add_ui (a , a, b[i]);
  }
}

int gpu_ecm_stage1 (mpz_t *factors, int *array_stage_found, mpz_t N, mpz_t s, 
                    unsigned int number_of_curves, unsigned int firstsigma, 
                    float *gputime, int verbose)
{
  int youpi = ECM_NO_FACTOR_FOUND;

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
  mpz_ui_pow_ui (w, 2, ECM_GPU_SIZE_DIGIT); /* Compute w = 2^SIZE_DIGIT */
    
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
  for (i = 0; i < number_of_curves; i++)
  {
    sigma = firstsigma + i;

    mpz_mul_ui (z2p, invw, sigma);
    mpz_mod (z2p, z2p, N);
    mpz_mul_2exp (z2p, z2p, 6);
    mpz_add_ui (z2p, z2p, 8);
    mpz_mod (z2p, z2p, N); /* z2p = 8+64*d */

    to_mont_repr (z2p, N);

    mpz_to_biguint (h_xarray[i], xp); 
    mpz_to_biguint (h_zarray[i], zp); 
    mpz_to_biguint (h_x2array[i], x2p); 
    mpz_to_biguint (h_z2array[i], z2p); 
  } 
 
  /* Call the wrapper function that call the GPU */
  *gputime=cuda_Main (h_N, h_3N, h_M, h_invN, h_xarray, h_zarray, h_x2array, 
                     h_z2array, s, firstsigma, number_of_curves, verbose);

  /* Analyse results */
  for (i = 0; i < number_of_curves; i++)
  {
    sigma = firstsigma + i;

    biguint_to_mpz (xp, h_xarray[i]); 
    biguint_to_mpz (zp, h_zarray[i]); 
    
    from_mont_repr (xp, N, invB);
    from_mont_repr (zp, N, invB);
  
    array_stage_found[i] = findfactor (factors[i], N, xp, zp);

    if (array_stage_found[i] != ECM_NO_FACTOR_FOUND)
      {
        youpi = array_stage_found[i];
        outputf (OUTPUT_NORMAL, "GPU: factor %Zd found in Step 1 with"
                " curve %u (-sigma 3:%u)\n", factors[i], i, sigma);
      }
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

  return youpi;
}

static void
A_from_sigma (mpz_t A, unsigned int sigma, mpz_t n)
{
  mpz_t tmp;
  int i;
  mpz_init_set_ui (tmp, sigma);
  /* Compute d = sigma/2^ECM_GPU_SIZE_DIGIT */
  for (i = 0; i < ECM_GPU_SIZE_DIGIT; i++)
    {
      if (mpz_tstbit (tmp, 0) == 1)
      mpz_add (tmp, tmp, n);
      mpz_div_2exp (tmp, tmp, 1);
    }
  mpz_mul_2exp (tmp, tmp, 2);           /* 4d */
  mpz_sub_ui (tmp, tmp, 2);             /* 4d-2 */
      
  mpz_set (A, tmp);

  mpz_clear (tmp);
}

int
gpu_ecm (mpz_t f, mpz_t x, int param, mpz_t firstsigma, mpz_t n, mpz_t go,
         double *B1done, double B1, mpz_t B2min_parm, mpz_t B2_parm, 
         unsigned long k, const int S, int verbose, int repr,
         int nobase2step2, int use_ntt, int sigma_is_A, FILE *os, FILE* es, 
         char *chkfilename ATTRIBUTE_UNUSED, char *TreeFilename, double maxmem,
         int (*stop_asap)(void), mpz_t batch_s, double *batch_last_B1_used, 
         int device, int *device_init, unsigned int *nb_curves)
{
  unsigned int i;
  int youpi = ECM_NO_FACTOR_FOUND;
  int factor_found = ECM_NO_FACTOR_FOUND;
  long st, st2;
  long tottime; /* at the end, total time in ms */
  unsigned int firstsigma_ui;
  float gputime = 0.0;
  mpz_t tmp_A;
  mpz_t *factors = NULL; /* Contains either a factor of n either end-of-stage-1
                         residue (depending of the value of array_stage_found */
  int *array_stage_found = NULL;
  /* Only for stage 2 */
  int base2 = 0;  /* If n is of form 2^n[+-]1, set base to [+-]n */
  int Fermat = 0; /* If base2 > 0 is a power of 2, set Fermat to base2 */
  int po2 = 0;    /* Whether we should use power-of-2 poly degree */
  /* Use only in stage 2 */
  mpmod_t modulus;
  curve P;
  mpz_t B2min, B2; /* Local B2, B2min to avoid changing caller's values */
  unsigned long dF;
  root_params_t root_params;

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

  /* Only param = ECM_PARAM_BATCH_32BITS_D is accepted on GPU */
  if (param == ECM_PARAM_DEFAULT)
      param = ECM_PARAM_BATCH_32BITS_D;
    
  if (param != ECM_PARAM_BATCH_32BITS_D)
    {
      outputf (OUTPUT_ERROR, "GPU: Error, only param = ECM_PARAM_BATCH_32BITS_D "
                             "is accepted on GPU.\n");
      return ECM_ERROR;
    }

  /* check that repr == ECM_MOD_DEFAULT or ECM_MOD_BASE2 (only for stage 2) */
  if (repr != ECM_MOD_DEFAULT && repr != ECM_MOD_BASE2)
      outputf (OUTPUT_ERROR, "GPU: Warning, the value of repr will be ignored "
      "for step 1 on GPU.\n");

  /* It is only for stage 2, it is not taken into account for GPU code */
  if (mpmod_init (modulus, n, repr) != 0)
    return ECM_ERROR;

  /* See what kind of number we have as that may influence optimal parameter 
     selection. Test for base 2 number. Note: this was already done by
     mpmod_init. */

  if (modulus->repr == ECM_MOD_BASE2)
    base2 = modulus->bits;

  /* For a Fermat number (base2 a positive power of 2) */
  for (Fermat = base2; Fermat > 0 && (Fermat & 1) == 0; Fermat >>= 1);
  if (Fermat == 1) 
    {
      Fermat = base2;
      po2 = 1;
    }
  else
      Fermat = 0;
 
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
      compute_s (batch_s, B1, NULL);
      outputf (OUTPUT_VERBOSE, "Computing batch product (of %" PRIu64
                               " bits) of primes up to B1=%1.0f took %ldms\n",
                               mpz_sizeinbase (batch_s, 2), B1, cputime () - st);
    }

  /* Set parameters for stage 2 */
  mpres_init (P.x, modulus);
  mpres_init (P.y, modulus);
  mpres_init (P.A, modulus);
  mpz_init (tmp_A);
  mpz_init (B2);
  mpz_init (B2min);

  youpi = set_stage_2_params (B2, B2_parm, B2min, B2min_parm, &root_params,
                              B1, &k, S, use_ntt, &po2, &dF,
                              TreeFilename, maxmem, Fermat, modulus);
  if (youpi == ECM_ERROR)
      goto end_gpu_ecm;
  

  /* Initialize the GPU if necessary */
  if (!*device_init)
    {
      st = cputime ();
      youpi = select_and_init_GPU (device, nb_curves,
                                   test_verbose (OUTPUT_VERBOSE));

      if (youpi != 0)
        {
          youpi = ECM_ERROR;
          goto end_gpu_ecm2;
        }

      outputf (OUTPUT_VERBOSE, "GPU: Selection and initialization of the device "
                               "took %ldms\n", elltime (st, cputime ()));
      /* TRICKS: If initialization of the device is too long (few seconds), */
      /* try running 'nvidia-smi -q -l' on the background .                 */
      *device_init = 1;
    }
  
  /* Init arrays */
  factors = (mpz_t *) malloc (*nb_curves * sizeof (mpz_t));
  ASSERT_ALWAYS (factors != NULL);

  array_stage_found = (int *) malloc (*nb_curves * sizeof (int));
  ASSERT_ALWAYS (array_stage_found != NULL);

  for (i = 0; i < *nb_curves; i++)
    {
      mpz_init (factors[i]);
      array_stage_found[i] = ECM_NO_FACTOR_FOUND;
    }


  /* Current code works only for sigma_is_A = 0 */
  if (sigma_is_A != 0)
    {
      outputf (OUTPUT_ERROR, "GPU: Not yet implemented.\n");
      youpi= ECM_ERROR;
      goto end_gpu_ecm;
    }

  ASSERT (sigma_is_A == 0);
  if (mpz_sgn (firstsigma) == 0)
    {
      /* generate random value in [2, 2^32 - nb_curves - 1] */
      mpz_set_ui (firstsigma, (get_random_ul () %
                               (TWO32 - 2 - *nb_curves)) + 2);
    }
  else /* sigma should be in [2, 2^32-nb_curves] */
    {
      if (mpz_cmp_ui (firstsigma, 2) < 0 || 
          mpz_cmp_ui (firstsigma, TWO32 - *nb_curves) >= 0)
        {
          outputf (OUTPUT_ERROR, "GPU: Error, sigma should be in [2,%lu]\n",
                                 TWO32 - *nb_curves - 1);
          youpi= ECM_ERROR;
          goto end_gpu_ecm;
        }
    }
  firstsigma_ui = mpz_get_ui (firstsigma);

  print_B1_B2_poly (OUTPUT_NORMAL, ECM_ECM, B1, *B1done,  B2min_parm, B2min,
                    B2, S, firstsigma, sigma_is_A, ECM_EC_TYPE_MONTGOMERY,
                    go, param, *nb_curves);
  outputf (OUTPUT_VERBOSE, "dF=%lu, k=%lu, d=%lu, d2=%lu, i0=%Zd\n", 
           dF, k, root_params.d1, root_params.d2, root_params.i0);

  if (go != NULL && mpz_cmp_ui (go, 1) > 0)
    {
      outputf (OUTPUT_ERROR, "GPU: Error, option -go is not allowed\n");
      youpi= ECM_ERROR;
      goto end_gpu_ecm;
    }

  if (test_verbose (OUTPUT_VERBOSE))
    {
      if (mpz_cmp_d (B2min, B1) != 0)
        {
          outputf (OUTPUT_VERBOSE, 
            "Can't compute success probabilities for B1 <> B2min\n");
        }
      else
        {
          rhoinit (256, 10);
          print_expcurves (B1, B2, dF, k, root_params.S, param);
        }
    }
  
  st = cputime ();
  youpi = gpu_ecm_stage1 (factors, array_stage_found, n, batch_s, *nb_curves, 
                          firstsigma_ui, &gputime, verbose);

  outputf (OUTPUT_NORMAL, "Computing %u Step 1 took %ldms of CPU time / "
                          "%.0fms of GPU time\n", *nb_curves, 
                                           elltime (st, cputime ()), gputime);
  outputf (OUTPUT_VERBOSE, "Throughput: %.3f curves per second ", 
                                                 1000 * (*nb_curves)/gputime);
  outputf (OUTPUT_VERBOSE, "(on average %.2fms per Step 1)\n", 
                                                        gputime/(*nb_curves));
  tottime = (long) gputime;

  *B1done=B1;

  /* Save stage 1 residues */
  mpz_set_ui (x, 0);
  for (i = 0; i < *nb_curves; i++)
    {
      mpz_mul (x, x, n);
      mpz_add (x, x, factors[i]);
    }

  /* was a factor found in stage 1 ? */
  if (youpi != ECM_NO_FACTOR_FOUND)
      goto end_gpu_ecm_rhotable;

  /* If using 2^k +/-1 modulus and 'nobase2step2' flag is set,
     set default (-nobase2) modular method and remap P.x, P.y, and P.A */
  if (modulus->repr == ECM_MOD_BASE2 && nobase2step2)
    {
      mpmod_clear (modulus);

      repr = ECM_MOD_NOBASE2;
      if (mpmod_init (modulus, n, repr) != 0) /* reset modulus for nobase2 */
        {
          youpi = ECM_ERROR;
          goto end_gpu_ecm_rhotable;
        }
    }

  if (mpz_cmp (B2, B2min) < 0)
      goto end_gpu_ecm_rhotable;

  st2 = cputime ();
  
  P.disc = 0; /* For stage2 this needs to be 0, in order not to use CM stuff */

  for (i = 0; i < *nb_curves; i++)
    {
      if (test_verbose (OUTPUT_RESVERBOSE)) 
        outputf (OUTPUT_RESVERBOSE, "x=%Zd\n", factors[i]);

      if (stop_asap != NULL && (*stop_asap) ())
          goto end_gpu_ecm_rhotable;
    
      mpres_set_z (P.x, factors[i], modulus);
      mpres_set_ui (P.y, 1, modulus);
      A_from_sigma (tmp_A, i+firstsigma_ui, modulus->orig_modulus);
      mpres_set_z (P.A, tmp_A, modulus);
  
      /* compute stage 2 */
      youpi = montgomery_to_weierstrass (factors[i], P.x, P.y, P.A, modulus);
      if (youpi != ECM_NO_FACTOR_FOUND)
        goto next_curve;

      if (test_verbose (OUTPUT_RESVERBOSE) && youpi == ECM_NO_FACTOR_FOUND
          && mpz_cmp (B2, B2min) >= 0)
        {
          mpz_t t;

          mpz_init (t);
          mpres_get_z (t, P.x, modulus);
          outputf (OUTPUT_RESVERBOSE, "After switch to Weierstrass form, "
                                      "P=(%Zd", t);
          mpres_get_z (t, P.y, modulus);
          outputf (OUTPUT_RESVERBOSE, ", %Zd)\n", t);
          mpres_get_z (t, P.A, modulus);
          outputf (OUTPUT_RESVERBOSE, "on curve Y^2 = X^3 + %Zd * X + b\n", 
                       t);
          mpz_clear (t);
        }
  
      /* It is a hack to avoid very verbose Step 2 
        (without it, stage2() prints a least a line by curves) */
      if (!test_verbose (OUTPUT_VERBOSE)) 
        set_verbose (0);
      youpi = stage2 (factors[i], &P, modulus, dF, k, &root_params, use_ntt, 
                      TreeFilename, stop_asap);
      set_verbose (verbose);

    next_curve:
      if (youpi != ECM_NO_FACTOR_FOUND)
        {
          array_stage_found[i] = youpi;
          outputf (OUTPUT_NORMAL, "GPU: factor %Zd found in Step 2 with"
                " curve %u (-sigma 3:%u)\n", factors[i], i, i+firstsigma_ui);
          /* factor_found corresponds to the first factor found */
          if (factor_found == ECM_NO_FACTOR_FOUND)
            factor_found = youpi;
        }
    }

  /* If a factor was found in Step 2, make sure we set
   * our return value "youpi" appropriately
   */
  youpi = factor_found;

  st2 = elltime (st2, cputime ());
  outputf (OUTPUT_NORMAL, "Computing %u Step 2 on CPU took %ldms\n", 
                                                              *nb_curves, st2);
  outputf (OUTPUT_VERBOSE, "Throughput: %.3f Step 2 per second ", 
                                  1000 * ((double)(*nb_curves))/((double)st2));
  outputf (OUTPUT_VERBOSE, "(on average %0.2fms per Step 2)\n", 
                                         ((double) st2)/((double) *nb_curves));
  tottime += st2;


end_gpu_ecm_rhotable:
  if (test_verbose (OUTPUT_VERBOSE))
    {
      if (mpz_cmp_d (B2min, B1) == 0)
        {
          if (youpi == ECM_NO_FACTOR_FOUND && 
              (stop_asap == NULL || !(*stop_asap)()))
              print_exptime (B1, B2, dF, k, root_params.S, 
                             (long) (tottime / *nb_curves), param);
          rhoinit (1, 0); /* Free memory of rhotable */
        }
    }

  /* If f0, ,fk are the factors found (in stage 1 or 2) 
   * f = f0 + f1*n + .. + fk*n^k
   * The purpose of this construction is to be able to return more than one
   * factor if needed without breaking the lib interface (as gcd(f,n)=gcd(f0,n).
   */
  mpz_set_ui (f, 0);
  for (i = 0; i < *nb_curves; i++)
  {
    if (array_stage_found[i] != ECM_NO_FACTOR_FOUND)
      {
        mpz_mul (f, f, n);
        mpz_add (f, f, factors[i]);
      }
  }

end_gpu_ecm:
  mpz_clear (root_params.i0);
  mpz_clear (B2);
  mpz_clear (B2min);

  for (i = 0; i < *nb_curves; i++)
      mpz_clear (factors[i]);

  free (array_stage_found);
  free (factors);

end_gpu_ecm2:
  mpz_clear (tmp_A);
  mpres_clear (P.A, modulus);
  mpres_clear (P.y, modulus);
  mpres_clear (P.x, modulus);
  mpmod_clear (modulus);

  return youpi;
}
#endif



