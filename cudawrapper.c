#include "ecm-gpu.h"

#ifdef WITH_GPU

#include "cudacommon.h"

#include "cgbn_stage1.h"


#define TWO32 4294967296 /* 2^32 */

/* Try to reduce all composite factors to primes.
 * This can be hard if factors overlap e.g. (a*b, a*c*d, b*c)
 */
void reducefactors (mpz_t *factors, int *array_found, unsigned int nb_curves)
{
  unsigned int i, j;
  unsigned int found;
  unsigned int updates;
  mpz_t gcd;
  mpz_init (gcd);

  found = 0;
  mpz_t *reduced = (mpz_t *) malloc (nb_curves * sizeof (mpz_t));
  ASSERT_ALWAYS (reduced != NULL);

  /* Add all unique factors to reduced */
  for (i = 0; i < nb_curves; i++)
    {
      if (array_found[i] == ECM_NO_FACTOR_FOUND)
        continue;

      /* Scan for match */
      updates = 0;
      for (j = 0; j < found; j++) {
          if (mpz_cmp (factors[i], reduced[j]) == 0) {
              updates = 1;
              break;
          }
      }
      if (!updates)
          mpz_init_set (reduced[found++], factors[i]);
    }

  do {
    outputf (OUTPUT_DEVVERBOSE, "GPU: Reducing %d factors\n", found);
    updates = 0;

    /* remove any trivial factor */
    for (i = 0; i < found; i++)
      {
        while (mpz_cmp_ui (reduced[i], 1) == 0) {
          found--;
          mpz_swap (reduced[i], reduced[found]);
          mpz_clear (reduced[found]);
          if (i == found)
              break;
        }
      }

    for (i = 0; i < found; i++)
      {
        /* Try to reduce an existing factor */
        for (j = i+1; j < found; j++)
          {
            /* if i == j remove reduced[j] */
            if (mpz_cmp (reduced[i], reduced[j]) == 0)
              {
                  updates += 1;
                  found--;
                  mpz_swap (reduced[j], reduced[found]);
                  mpz_clear (reduced[found]);
                  if (j == found)
                      break;
              }

            mpz_gcd (gcd, reduced[i], reduced[j]);
            if (mpz_cmp_ui (gcd, 1) > 0)
              {
                /* gcd(2*3, 2*3*5) remove 2*3 from F2 leaving 2*3 and 5 */
                if (mpz_cmp (gcd, reduced[i]) == 0)
                  {
                    updates += 1;
                    assert( mpz_divisible_p (reduced[j], gcd) );
                    mpz_divexact (reduced[j], reduced[j], gcd);
                  }
                /* gcd(2*3*5, 2*3) == 2*3 from F1 leaving 5 and 2*3 */
                else if (mpz_cmp (gcd, reduced[j]) == 0)
                  {
                    updates += 1;
                    assert( mpz_divisible_p (reduced[i], gcd) );
                    mpz_divexact (reduced[i], reduced[i], gcd);
                  }

                /* hard case gcd(2*3, 3*5) = 3, remove 3 from both, add 3 as new factor */
                else if (found < nb_curves)
                  {
                    updates += 1;
                    mpz_divexact (reduced[j], reduced[j], gcd);
                    mpz_divexact (reduced[i], reduced[i], gcd);

                    mpz_init (reduced[found]);
                    mpz_set (reduced[found], gcd);
                    found++;
                  }
              }
            if (mpz_cmp_ui (reduced[i], 1) == 0)
                break;
          }
      }
  } while (updates > 0);

  /* bubble_sort, fast enough because found < num_curves */
  do {
    updates = 0;
    for (j = 1; j < found; j++)
      {
        if (mpz_cmp(reduced[j-1], reduced[j]) > 0)
          {
            updates += 1;
            mpz_swap(reduced[j-1], reduced[j]);
          }
      }
  } while (updates > 0);

  outputf (OUTPUT_DEVVERBOSE, "GPU: Reduced to %d factors\n", found);
  /* write out reduced[i], update array_found */
  for (i = 0; i < found; i++)
    {
      mpz_swap(factors[i], reduced[i]);
      mpz_clear(reduced[i]);
      array_found[i] = ECM_FACTOR_FOUND_STEP1;
      outputf (OUTPUT_DEVVERBOSE, "GPU: Reduced factor %d: %Zd\n", i+1, factors[i]);
    }

  for (i = found; i < nb_curves; i++)
    array_found[i] = ECM_NO_FACTOR_FOUND;

  mpz_clear (gcd);
  free(reduced);
}


static void
A_from_sigma (mpz_t A, unsigned int sigma, mpz_t n)
{
  mpz_t tmp;
  int i;
  mpz_init_set_ui (tmp, sigma);
  /* Compute d = sigma/2^32 */
  for (i = 0; i < 32; i++)
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
                         residue (depending of the value of array_found */
  int *array_found = NULL;
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

  /* Set global VERBOSE to avoid the need to explicitly passing verbose */
  set_verbose (verbose);
  ECM_STDOUT = (os == NULL) ? stdout : os;
  ECM_STDERR = (es == NULL) ? stdout : es;


  /* Check that N is not too big */
  size_t max_bits = ECM_GPU_CGBN_MAX_BITS - 6;
  if (mpz_sizeinbase (n, 2) > max_bits)
    {
      outputf (OUTPUT_ERROR, "GPU: Error, input number should be stricly lower"
                             " than 2^%d\n", max_bits);
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

  /* Initialize the GPU if necessary and determine nb_curves */
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

  array_found = (int *) malloc (*nb_curves * sizeof (int));
  ASSERT_ALWAYS (array_found != NULL);

  for (i = 0; i < *nb_curves; i++)
    {
      mpz_init (factors[i]);
      array_found[i] = ECM_NO_FACTOR_FOUND;
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

  st = cputime ();

  youpi = cgbn_ecm_stage1 (factors, array_found, n, batch_s, *nb_curves,
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

  /* GMP documentation says mpz_sizeinbase(op, 2) is always the exact value. */
  size_t n_bits = mpz_sizeinbase(n, 2);

  /* Save stage 1 residues as x = x0 + x1 * 2^bits + ... + xk * 2^(bits*k) */
  mpz_set_ui (x, 0);
  for (i = 0; i < *nb_curves; i++)
    {
      mpz_mul_2exp (x, x, n_bits);
      mpz_add (x, x, factors[*nb_curves - 1 - i]);
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
      /* hack to reduce verbose Step 2 */
      if (verbose > 0)
        set_verbose (verbose-1);

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
 
      youpi = stage2 (factors[i], &P, modulus, dF, k, &root_params, use_ntt, 
                      TreeFilename, i+1, stop_asap);
      
    next_curve:
      set_verbose (verbose);

      if (youpi != ECM_NO_FACTOR_FOUND)
        {
          array_found[i] = youpi;
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

  reducefactors(factors, array_found, *nb_curves);

  /* If f0, ,fk are the factors found (in stage 1 or 2) 
   * f = f0 + f1*n + .. + fk*n^k
   * The purpose of this construction is to be able to return more than one
   * factor if needed without breaking the lib interface (as gcd(f,n)=gcd(f0,n).
   */
  mpz_set_ui (f, 0);
  for (i = 0; i < *nb_curves; i++)
  {
    /* invert order of factors so they are processed in same order found */
    if (array_found[*nb_curves-1-i] != ECM_NO_FACTOR_FOUND)
      {
        mpz_mul (f, f, n);
        mpz_add (f, f, factors[*nb_curves-1-i]);
      }
  }

end_gpu_ecm:
  mpz_clear (root_params.i0);
  mpz_clear (B2);
  mpz_clear (B2min);

  for (i = 0; i < *nb_curves; i++)
      mpz_clear (factors[i]);

  free (array_found);
  free (factors);

end_gpu_ecm2:
  mpz_clear (tmp_A);
  mpres_clear (P.A, modulus);
  mpres_clear (P.y, modulus);
  mpres_clear (P.x, modulus);
  mpmod_clear (modulus);

  return youpi;
}

#endif /* HAVE_GPU */
