#include "ecm-gpu.h"

#ifdef WITH_GPU

#include "cudacommon.h"

#include "cgbn_stage1.h"
#include "ecm-ecm.h"


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
gpu_ecm (mpz_t f, const ecm_params params, ecm_params mutable_params, mpz_t n, double B1)
{
  unsigned int i;
  int youpi = ECM_NO_FACTOR_FOUND;
  int factor_found = ECM_NO_FACTOR_FOUND;
  long st, st2;
  long tottime; /* at the end, total time in ms */
  unsigned int firstsigma_ui;
  float gputime = 0.0;
  mpz_t tmp_A;
  mpz_t *factors = NULL; /* Contains either a factor(s) of n or end-of-stage-1
                            residue (depending of the value of array_found). */
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
  unsigned int nb_curves = 0; /* Local copy of number of curves */

  /* This helps keep track of when params is changed in this function. */
  assert((void*) params == (void*) mutable_params); /* params != mutable params */

  ASSERT((-1 <= params->sigma_is_A) && (params->sigma_is_A <= 1));
  ASSERT((GMP_NUMB_BITS == 32) || (GMP_NUMB_BITS == 64));

  /* Set global VERBOSE to avoid the need to explicitly passing verbose */
  set_verbose (params->verbose);
  ECM_STDOUT = (params->os == NULL) ? stdout : params->os;
  ECM_STDERR = (params->es == NULL) ? stdout : params->es;

  /* Check that N is not too big */
  size_t max_bits = ECM_GPU_CGBN_MAX_BITS - 6;
  if (mpz_sizeinbase (n, 2) > max_bits)
    {
      outputf (OUTPUT_ERROR, "GPU: Error, input number should be stricly lower"
                             " than 2^%d\n", max_bits);
      return ECM_ERROR;
    }

  /* Only param = ECM_PARAM_BATCH_32BITS_D is accepted on GPU */
  if (params->param == ECM_PARAM_DEFAULT)
      mutable_params->param = ECM_PARAM_BATCH_32BITS_D;

  if (params->param != ECM_PARAM_BATCH_32BITS_D)
    {
      outputf (OUTPUT_ERROR, "GPU: Error, only param = ECM_PARAM_BATCH_32BITS_D "
                             "is accepted on GPU.\n");
      return ECM_ERROR;
    }

  if (params->go != NULL && mpz_cmp_ui (params->go, 1) > 0)
    {
      outputf (OUTPUT_ERROR, "GPU: Error, option -go is not allowed\n");
      return ECM_ERROR;
    }

  /* Cannot do resume on GPU */
  if (!ECM_IS_DEFAULT_B1_DONE(params->B1done) && params->B1done < B1)
    {
      outputf (OUTPUT_ERROR, "GPU: Error, cannot resume on GPU.\n");
      return ECM_ERROR;
    }

  /* Current code works only for sigma_is_A = 0 */
  if (params->sigma_is_A != 0)
    {
      outputf (OUTPUT_ERROR, "GPU: ERROR, sigma_is_A not allowed.\n");
      return ECM_ERROR;
    }

  /* check that repr == ECM_MOD_DEFAULT or ECM_MOD_BASE2 (only for stage 2) */
  if (params->repr != ECM_MOD_DEFAULT && params->repr != ECM_MOD_BASE2)
      outputf (OUTPUT_ERROR, "GPU: Warning, the value of repr will be ignored "
      "for step 1 on GPU.\n");

  /* It is only for stage 2, it is not taken into account for GPU code */
  if (mpmod_init (modulus, n, params->repr) != 0)
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

  /* Set parameters for stage 2 */
  mpres_init (P.x, modulus);
  mpres_init (P.y, modulus);
  mpres_init (P.A, modulus);
  mpz_init (tmp_A);
  mpz_init (B2);
  mpz_init (B2min);

  youpi = set_stage_2_params (
          B2, params->B2,
          B2min, params->B2min,
          &root_params,
          B1, &mutable_params->k, params->S, params->use_ntt, &po2, &dF,
                              params->TreeFilename, params->maxmem, Fermat, modulus);
  if (youpi == ECM_ERROR)
      goto end_gpu_ecm;

  /* Initialize the GPU if necessary and determine nb_curves */
  if (!params->gpu_device_init)
    {
      st = cputime ();
      youpi = select_and_init_GPU (
          params->gpu_device,
          &mutable_params->gpu_number_of_curves,
          test_verbose (OUTPUT_VERBOSE));

      if (youpi != 0)
        {
          youpi = ECM_ERROR;
          goto end_gpu_ecm;
        }

      outputf (OUTPUT_VERBOSE, "GPU: Selection and initialization of the device "
                               "took %ldms\n", elltime (st, cputime ()));
      /* TRICKS: If initialization of the device is too long (few seconds), */
      /* try running 'nvidia-smi -q -l' on the background .                 */
      mutable_params->gpu_device_init = 1;
    }

  /* Number of curves is only set after select_and_init_GPU */
  nb_curves = params->gpu_number_of_curves;

  ASSERT (params->sigma_is_A == 0);
  if (mpz_sgn (params->sigma) == 0)
    {
      /* generate random value in [2, 2^32 - nb_curves - 1] */
      mpz_set_ui (mutable_params->sigma,
                  (get_random_ul () % (TWO32 - 2 - nb_curves)) + 2);
    }
  else /* sigma should be in [2, 2^32-nb_curves] */
    {
      if (mpz_cmp_ui (params->sigma, 2) < 0 ||
          mpz_cmp_ui (params->sigma, TWO32 - nb_curves) >= 0)
        {
          outputf (OUTPUT_ERROR, "GPU: Error, sigma should be in [2,%lu]\n",
                                 TWO32 - nb_curves - 1);
          youpi = ECM_ERROR;
          goto end_gpu_ecm;
        }
    }
  firstsigma_ui = mpz_get_ui (params->sigma);

  print_B1_B2_poly (OUTPUT_NORMAL, ECM_ECM, B1, params->B1done,  params->B2min, B2min,
                    B2, params->S, params->sigma, params->sigma_is_A, ECM_EC_TYPE_MONTGOMERY,
                    params->go, params->param, nb_curves);
  outputf (OUTPUT_VERBOSE, "dF=%lu, k=%lu, d=%lu, d2=%lu, i0=%Zd\n",
           dF, params->k, root_params.d1, root_params.d2, root_params.i0);

  if ((mpz_cmp (B2, B2min) < 0) && modulus->repr == ECM_MOD_BASE2 && params->nobase2step2)
    {
      outputf (OUTPUT_ERROR, "GPU: ERROR, -nobase2s2 not supported on GPU.\n");
      youpi = ECM_ERROR;
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
          print_expcurves (B1, B2, dF, params->k, root_params.S, params->param);
        }
    }

  /* Init arrays */
  factors = (mpz_t *) malloc (nb_curves * sizeof (mpz_t));
  ASSERT_ALWAYS (factors != NULL);

  array_found = (int *) malloc (nb_curves * sizeof (int));
  ASSERT_ALWAYS (array_found != NULL);

  for (i = 0; i < nb_curves; i++)
    {
      mpz_init (factors[i]);
      array_found[i] = ECM_NO_FACTOR_FOUND;
    }


  /* Compute s */
  if (B1 != params->batch_last_B1_used || mpz_cmp_ui (params->batch_s, 1) <= 0)
    {
      mutable_params->batch_last_B1_used = B1;

      st = cputime ();
      /* construct the batch exponent */
      compute_s (mutable_params->batch_s, B1, /* B1done */ 0, NULL);
      outputf (OUTPUT_VERBOSE, "Computing batch product (of %" PRIu64
                               " bits) of primes up to B1=%1.0f took %ldms\n",
                               mpz_sizeinbase (params->batch_s, 2), B1, cputime () - st);
    }


  st = cputime ();

  youpi = cgbn_ecm_stage1 (factors, array_found, n, params->batch_s, nb_curves,
                           firstsigma_ui, &gputime, params->verbose);

  outputf (OUTPUT_NORMAL, "Computing %u Step 1 took %ldms of CPU time / "
                          "%.0fms of GPU time\n", nb_curves,
                          elltime (st, cputime ()), gputime);
  outputf (OUTPUT_VERBOSE, "Throughput: %.3f curves per second ",
                           1000 * nb_curves/gputime);
  outputf (OUTPUT_VERBOSE, "(on average %.2fms per Step 1)\n",
                           gputime/nb_curves);
  tottime = (long) gputime;

  mutable_params->B1done = B1;

  /* GMP documentation says mpz_sizeinbase(op, 2) is always the exact value. */
  size_t n_bits = mpz_sizeinbase(n, 2);

  /* Save stage 1 residues as x = x0 + x1 * 2^bits + ... + xk * 2^(bits*k) */
  mpz_set_ui (mutable_params->x, 0);
  /* Equivalent to using mpz_mul_2exp and mpz_add while avoiding O(n*k) limp copies */
  mpz_realloc2(mutable_params->x, nb_curves * n_bits);
  for (i = 0; i < nb_curves; i++)
    for (size_t j = 0; j < n_bits; j++)
      if (mpz_tstbit (factors[i], j))
        mpz_setbit(mutable_params->x, j + n_bits * i);

  /* was a factor found in stage 1 ? */
  if (youpi != ECM_NO_FACTOR_FOUND)
      goto end_gpu_ecm_factors;

  if (mpz_cmp (B2, B2min) < 0)
      goto end_gpu_ecm_factors;

  st2 = cputime ();

  P.disc = 0; /* For stage2 this needs to be 0, in order not to use CM stuff */

  for (i = 0; i < nb_curves; i++)
    {
      /* hack to reduce verbose Step 2 */
      if (params->verbose > 0)
        set_verbose (params->verbose-1);

      if (test_verbose (OUTPUT_RESVERBOSE))
        outputf (OUTPUT_RESVERBOSE, "x=%Zd\n", factors[i]);

      if (params->stop_asap != NULL && params->stop_asap())
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

      youpi = stage2 (factors[i], &P, modulus, dF, params->k, &root_params,
                      params->use_ntt, params->TreeFilename, i+1, params->stop_asap);

    next_curve:
      set_verbose (params->verbose);

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
                          nb_curves, st2);
  outputf (OUTPUT_VERBOSE, "Throughput: %.3f Step 2 per second ",
                           1000.0 * nb_curves /st2);
  outputf (OUTPUT_VERBOSE, "(on average %0.2fms per Step 2)\n",
                           ((double) st2)/((double) nb_curves));
  tottime += st2;

end_gpu_ecm_factors:

  reducefactors(factors, array_found, nb_curves);

  /* If f0, ,fk are the factors found (in stage 1 or 2)
   * f = f0 + f1*n + .. + fk*n^k
   * The purpose of this construction is to be able to return more than one
   * factor if needed without breaking the lib interface (as gcd(f,n)=gcd(f0,n).
   */
  mpz_set_ui (f, 0);
  for (i = 0; i < nb_curves; i++)
  {
    /* invert order of factors so they are processed in same order found */
    if (array_found[nb_curves-1-i] != ECM_NO_FACTOR_FOUND)
      {
        mpz_mul (f, f, n);
        mpz_add (f, f, factors[nb_curves-1-i]);
      }
  }

  for (i = 0; i < nb_curves; i++)
      mpz_clear (factors[i]);

  free (array_found);
  free (factors);

end_gpu_ecm_rhotable:
  if (test_verbose (OUTPUT_VERBOSE))
    {
      if (mpz_cmp_d (B2min, B1) == 0)
        {
          if (youpi == ECM_NO_FACTOR_FOUND &&
              (params->stop_asap == NULL || !params->stop_asap()))
              print_exptime (B1, B2, dF, params->k, root_params.S,
                             (long) (tottime / nb_curves), params->param);
          rhoinit (1, 0); /* Free memory of rhotable */
        }
    }


end_gpu_ecm:
  mpz_clear (root_params.i0);
  mpres_clear (P.A, modulus);
  mpres_clear (P.y, modulus);
  mpres_clear (P.x, modulus);
  mpmod_clear (modulus);
  mpz_clear (tmp_A);
  mpz_clear (B2);
  mpz_clear (B2min);

  return youpi;
}


// TODO trying to get two arrays out. residuals (in x) and n's in???
/* Input: infile file containing a series of unique numbers
          B1 is the stage 1 bound
          params is ecm_parameters
   Output: n is array of candidates numbers,
           f is array of factors found,
           x is array of stage-1 residual,
           mutable_params->gpu_number_of_curves determines length of array
   Return value: non-zero iff a factor is found (1 for stage 1, 2 for stage 2)
*/
int
gpu_pm1 (FILE * infile, char* savefilename, mpcandi_t **n, mpz_t **f, mpz_t **x,
         const ecm_params params, ecm_params mutable_params, double B1)
{
  unsigned int i;
  int youpi = ECM_NO_FACTOR_FOUND;
  long st;
  float gputime = 0.0;

  /* Local pointers to params->gpu_return1, params->gpu_return2, params->gpu_return3. */
  mpcandi_t *inputs = NULL;
  mpz_t *numbers = NULL;
  mpz_t *factors = NULL; 
  mpz_t *residuals = NULL; 

  unsigned int nb_curves = 0; /* Local copy of number of curves */

  /* This helps keep track of when params is changed in this function. */
  assert((void*) params == (void*) mutable_params); /* params != mutable params */

  ASSERT((GMP_NUMB_BITS == 32) || (GMP_NUMB_BITS == 64));

  /* Set global VERBOSE to avoid the need to explicitly passing verbose */
  set_verbose (params->verbose);
  ECM_STDOUT = (params->os == NULL) ? stdout : params->os;
  ECM_STDERR = (params->es == NULL) ? stdout : params->es;

  /* Check for things GPU P-1 doesn't currently support */
  if (params->go != NULL && mpz_cmp_ui (params->go, 1) > 0)
    {
      outputf (OUTPUT_ERROR, "GPU: Error, option -go is not allowed\n");
      return ECM_ERROR;
    }

  if (mpz_cmp_ui (params->B2, 0) != 0)
    {
      outputf (OUTPUT_ERROR, "GPU P-1: Only supports stage 1, requires B2=0\n");
      return ECM_ERROR;
    }

  /* check that repr == ECM_MOD_DEFAULT or ECM_MOD_BASE2 (only for stage 2) */
  if (params->repr != ECM_MOD_DEFAULT && params->repr != ECM_MOD_BASE2)
      outputf (OUTPUT_ERROR, "GPU: Warning, the value of repr will be ignored "
      "for step 1 on GPU.\n");

  /* Cannot do resume on GPU */
  if (!ECM_IS_DEFAULT_B1_DONE(params->B1done) && params->B1done < B1)
    {
      outputf (OUTPUT_ERROR, "GPU: Error, cannot resume on GPU.\n");
      return ECM_ERROR;
    }

  /* Cannot do resume on GPU */
  if (infile == NULL || feof(infile))
    {
      printf("Hi %d\n", feof(infile));
      outputf (OUTPUT_ERROR, "GPU: Must pass -inp file for GPU P-1.\n");
      return ECM_ERROR;
    }

  /* Initialize the GPU if necessary and determine nb_curves */
  if (!params->gpu_device_init)
    {
      st = cputime ();
      youpi = select_and_init_GPU (
              params->gpu_device,
              &mutable_params->gpu_number_of_curves,
              test_verbose (OUTPUT_VERBOSE));
      if (youpi != 0)
        return ECM_ERROR;

      outputf (OUTPUT_VERBOSE, "GPU: Selection and initialization of the device "
                               "took %ldms\n", elltime (st, cputime ()));
      /* TRICKS: If initialization of the device is too long (few seconds), */
      /* try running 'nvidia-smi -q -l' on the background .                 */
      mutable_params->gpu_device_init = 1;
    }
  // Set local copy of number of curves
  nb_curves = params->gpu_number_of_curves;

  if (mpz_cmp_ui (params->x, 0) == 0)
    {
      mpz_t temp;
      mpz_init_set_ui(temp, 0xFFFFFFFF);
      // TODO figure out random include
      mpz_set_ui(mutable_params->x, 3);
      //__ecm_pm1_random_seed (mutable_params->x, temp, mutable_params->rng);
      mpz_clear(temp);
    }
  outputf (OUTPUT_VERBOSE, "GPU P-1: Using x0=%Zd\n", mutable_params->x);

  /* Init arrays */
  inputs =  (mpcandi_t *) malloc(nb_curves * sizeof (mpcandi_t));
  ASSERT_ALWAYS (inputs != NULL);
  numbers = (mpz_t *) malloc (nb_curves * sizeof (mpz_t));
  ASSERT_ALWAYS (numbers != NULL);
  factors = (mpz_t *) malloc (nb_curves * sizeof (mpz_t));
  ASSERT_ALWAYS (factors != NULL);
  residuals = (mpz_t *) malloc (nb_curves * sizeof (mpz_t));
  ASSERT_ALWAYS (residuals != NULL);

  *n = inputs;
  *f = factors;
  *x = residuals;

  /* TODO figure out interface for this */
  assert (infile != NULL);
  outputf (OUTPUT_VERBOSE, "GPU P-1: Loading numbers from 'TODO'\n");

  for (i = 0; i < nb_curves; i++)
    {
      mpcandi_t_init(&inputs[i]);
      if (read_number(&inputs[i], infile, 1))
        {
          mpz_init_set (numbers[i], inputs[i].n);
          mpz_init (factors[i]);
          mpz_init (residuals[i]);
           
          // Validity checks, > 1, odd
          if (mpz_cmp_ui(numbers[i], 1) <= 0)
            {
              fprintf(ECM_STDERR, "Error, n= should be great than 1.\n");
              gmp_fprintf(ECM_STDERR, "n[%u]=%Zd\n", i, numbers[i]);
              // GOTO FREE
              return ECM_ERROR;
            }
          /* TODO figure out how to divide out 2's here and add to factors[i] */
          if (mpz_tstbit(numbers[i], 0) == 0)
            {
              fprintf(ECM_STDERR, "Error, numbers should all be odd.\n");
              gmp_fprintf(ECM_STDERR, "n[%u]=%Zd\n", i, numbers[i]);
              // GOTO FREE
              return ECM_ERROR;
            }
        }
      else
        {
            mpcandi_t_free(&inputs[i]);

            outputf (OUTPUT_VERBOSE,
                     "GPU P-1: End of input truncating to %i curves\n", i);
            // Reduce to running i curves
            nb_curves = i;
            mutable_params->gpu_number_of_curves = i;
            break;
        }
    }

  /* Compute s */
  // TODO consider bload
  ASSERT_ALWAYS (B1 != params->batch_last_B1_used);
  ASSERT_ALWAYS (mpz_cmp_ui (params->batch_s, 1) <= 0);


  st = cputime ();
  /* construct the batch exponent */
  compute_s (mutable_params->batch_s, B1, NULL);
  outputf (OUTPUT_VERBOSE, "Computing batch product (of %" PRIu64
                           " bits) of primes up to B1=%1.0f took %ldms\n",
                           mpz_sizeinbase (params->batch_s, 2), B1, cputime () - st);
  mutable_params->batch_last_B1_used = B1;

  st = cputime ();

  youpi = cgbn_pm1_stage1 (
      numbers, factors, residuals,
      params->batch_s, params->x, nb_curves, &gputime, params->verbose);

  outputf (OUTPUT_NORMAL, "Computing %u P-1 Step 1 took %ldms of CPU time / "
                          "%.0fms of GPU time\n",
                          nb_curves, elltime (st, cputime ()), gputime);
  outputf (OUTPUT_VERBOSE, "Throughput: %.3f numbers per second ",
                                                 1000 * nb_curves/gputime);
  outputf (OUTPUT_VERBOSE, "(on average %.2fms per Step 1)\n",
                                                        gputime/nb_curves);

  /* These have to be saved and restored on each output. */
  mutable_params->B1done = B1;

  /* Variables needed for write_resumefile */
  mpz_t x0;
  char comment[2] = "";
  mpz_init_set(x0, params->x);

  int factors_found = 0;
  for (i = 0; i < nb_curves; i++)
      factors_found += mpz_cmp_ui(factors[i], 1) > 0;

  if (factors_found)
      fprintf(ECM_STDOUT, "\n\n");

  /* Copy out result from saved i'th P-1 results. */
  for (i = 0; i < nb_curves; i++)
    {
      outputf (OUTPUT_TRACE, "%d -> %Zd -> %Zd | %Zd\n", i, numbers[i], factors[i], residuals[i]);
      ASSERT_ALWAYS (mpz_cmp (numbers[i], inputs[i].n) == 0);

      if (savefilename != NULL)
        {
          if (mpz_cmp_ui (factors[i], 1) > 0) {
              ASSERT_ALWAYS (mpz_divisible_p (inputs[i].n, factors[i]));

              // TODO this or possible just mpcandi_t_addfactor and some print statements

              unsigned int fake_cnt = 1;
              mpz_t fake_lastfac;
              mpz_init(fake_lastfac);
              // TODO not sure what these should be set to
              int resume_wasPrp = 0;

              process_newfactor (
                      factors[i], ECM_FACTOR_FOUND_STEP1, &inputs[i], ECM_PM1,
                      /* returncode */ 0, /* gpu (TODO change var to is_ecm_gpu) */ 0,
                      &fake_cnt, &resume_wasPrp, fake_lastfac, /* resumefile */ NULL,
                      params->verbose, /* deep */ 0);

              mpz_clear (fake_lastfac);
          }


          /* write_resumefile expects residual in params->x */
          mpz_set(mutable_params->x, residuals[i]);

          write_resumefile (savefilename, ECM_PM1, mutable_params, &inputs[i],
                  /* orig_n */ numbers[i], /* orig_x0 */ x0, /* orig_y0 */ NULL, comment);
        }
    }

  mpz_set(mutable_params->x, x0);
  mpz_clear(x0);

  if (factors_found)
    {
      fprintf(ECM_STDOUT, "\n\n");
      return ECM_FACTOR_FOUND_STEP1;
    }

  return ECM_NO_FACTOR_FOUND;
}

#endif /* HAVE_GPU */
