/* gpu_pm1.h: header for GPU (CGBN) P-1stage 1.

Copyright 2021-2026 Seth Troisi

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
more details.

You should have received a copy of the GNU General Public License
along with this program; see the file COPYING.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#ifndef _GPU_PM1_CU
#define _GPU_PM1_CU 1

#ifndef __CUDACC__
#error "This file should only be compiled with nvcc"
#endif

#include "gpu_pm1.h"

#include <cassert>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

// GMP import must proceed cgbn.h
#include <gmp.h>
#include <cgbn.h>
#include <cuda.h>

#include "cudacommon.h"

#include "ecm.h"
#include "ecm-gpu.h"
#include "ecm-ecm.h"
#include "batched_s.h"

#define BATCH_MS_GOAL 100

// See cgbn_error_t enum (cgbn.h:39)
#define cgbn_normalized_error ((cgbn_error_t) 14)
#define cgbn_positive_overflow ((cgbn_error_t) 15)
#define cgbn_negative_overflow ((cgbn_error_t) 16)

// Seems to adds very small overhead (1-10%)
#define VERIFY_NORMALIZED 0
// Adds even less overhead (<1%)
#define CHECK_ERROR 1

// Tested with check_gpuecm.sage
#define CARRY_BITS 6

// Can dramatically change compile time
#if 1
    #define FORCE_INLINE __forceinline__
#else
    #define FORCE_INLINE
#endif

// support routine copied from  "CGBN/samples/utility/support.h"
static
void cgbn_check(cgbn_error_report_t *report, const char *file=NULL, int32_t line=0) {
  // check for cgbn errors

  if(cgbn_error_report_check(report)) {
    fprintf (stderr, "\n");
    fprintf (stderr, "CGBN error occurred: %s\n", cgbn_error_string(report));

    if(report->_instance!=0xFFFFFFFF) {
      fprintf (stderr, "Error reported by instance %d", report->_instance);
      if(report->_blockIdx.x!=0xFFFFFFFF)
        fprintf (stderr, ", blockIdx=(%d, %d, %d)", report->_blockIdx.x, report->_blockIdx.y, report->_blockIdx.z);
      if(report->_threadIdx.x!=0xFFFFFFFF)
        fprintf (stderr, ", threadIdx=(%d, %d, %d)", report->_threadIdx.x, report->_threadIdx.y, report->_threadIdx.z);
      fprintf (stderr, "\n");
    }
    else {
      fprintf (stderr, "Error reported by blockIdx=(%d %d %d)", report->_blockIdx.x, report->_blockIdx.y, report->_blockIdx.z);
      fprintf (stderr, "threadIdx=(%d %d %d)\n", report->_threadIdx.x, report->_threadIdx.y, report->_threadIdx.z);
    }
    if(file!=NULL)
      fprintf (stderr, "file %s, line %d\n", file, line);
    exit(1);
  }
}

#define CGBN_CHECK(report) cgbn_check(report, __FILE__, __LINE__)

static void
to_mpz(mpz_t r, const uint32_t *x, uint32_t count) {
  mpz_import (r, count, -1, sizeof(uint32_t), 0, 0, x);
}

static void
from_mpz(const mpz_t s, uint32_t *x, uint32_t count) {
  size_t words;

  if(mpz_sizeinbase (s, 2) > count * 32) {
    fprintf (stderr, "from_mpz failed -- result does not fit\n");
    exit(EXIT_FAILURE);
  }

  /* Zero all words just to be safe */
  for (words = 0; words < count; words++)
      x[words] = 0;
  mpz_export (x, &words, -1, sizeof(uint32_t), 0, 0, s);
  assert (words <= count);
}



// ---------------------------------------------------------------- //

// The CGBN context uses the following three parameters:
//   TBP             - threads per block (zero means to use the blockDim.x)
//   MAX_ROTATION    - must be small power of 2, empirically, 4 works well
//   CONSTANT_TIME   - require constant time algorithms (currently, constant time algorithms are not available)

// Locally it will also be helpful to have several parameters:
//   TPI             - threads per instance
//   BITS            - number of bits per instance

/* Doesn't seem to have any noticeable impact on performance. */
/* NOTE: >= 512 may not be supported for > 2048 bit kernels */
const uint32_t TPB_DEFAULT = 256;

template<uint32_t tpi, uint32_t bits>
class cgbn_params_t {
  public:
  // parameters used by the CGBN context
  static const uint32_t TPB=TPB_DEFAULT;           // Reasonable default
  static const uint32_t MAX_ROTATION=4;            // good default value
  static const uint32_t SHM_LIMIT=0;               // no shared mem available
#ifndef _MSC_VER
  static const bool     CONSTANT_TIME=false;       // not implemented
#endif

  // parameters used locally in the application
  static const uint32_t TPI=tpi;                   // threads per instance
  static const uint32_t BITS=bits;                 // instance size
  static const uint32_t WINDOW_BITS=7;             // For P-1, 2^N values are pre-computed
};


static
int verify_size_of_n(const mpz_t N, size_t max_bits) {
  size_t n_log2 = mpz_sizeinbase(N, 2);

  /* Using check_gpuecm.sage it looks like 4 bits would suffice. */
  size_t max_usable_bits = max_bits - CARRY_BITS;

  if (n_log2 <= max_usable_bits)
    return ECM_NO_FACTOR_FOUND;

  outputf (OUTPUT_ERROR,
          "GPU P-1: N(%d bits) + carry(%d bits) > BITS(%d)\n",
          n_log2, CARRY_BITS, max_bits);
  outputf (OUTPUT_ERROR,
          "GPU P-1: Error, input number should be stricly lower than 2^%d\n",
          max_usable_bits);
  return ECM_ERROR;
}


static void
set_s_bits_copy_to_gpu(
    const mpz_t s, const uint64_t buffer_bytes,
    uint32_t* const s_bits, uint32_t *gpu_s_bits)
{
  memset(s_bits, 0, buffer_bytes);

   uint64_t countp;
  mpz_export (s_bits, &countp, -1, sizeof(uint32_t), /* endian */ 0, 0, s);
  assert (4 * countp < buffer_bytes);

  CUDA_CHECK(cudaMemcpy(gpu_s_bits, s_bits, buffer_bytes, cudaMemcpyHostToDevice));
}


static
int print_nth_batch(int n)
{
  return ((n < 2) ||
          (n < 20 && n % 10 == 0) ||
          (n < 500 && n % 100 == 0) ||
#if BATCH_MS_GOAL >= 1000
          (n % 1000 == 0));
#else
          (n < 5000 && n % 1000 == 0) ||
          (n % 10000 == 0));
#endif
}

///////////////////////////////////////////////////////////////////////////////
//  GPU P-1 code starts here, I tried to put this in cgbn_pm1.cu but ran     //
//  into issues with duplicate symbols from CGBN.h (monitor and error_report)//
///////////////////////////////////////////////////////////////////////////////

template<class params>
class power_mod_t {
  public:

  typedef cgbn_context_t<params::TPI, params>   context_t;
  typedef cgbn_env_t<context_t, params::BITS>   env_t;
  typedef typename env_t::cgbn_t                bn_t;
  typedef typename env_t::cgbn_local_t          bn_local_t;
  typedef cgbn_mem_t<params::BITS>              mem_t;

  context_t _context;
  env_t     _env;
  int32_t   _instance; // which curve instance is this

  // Constructor
  __device__ FORCE_INLINE power_mod_t(cgbn_monitor_t monitor, cgbn_error_report_t *report, int32_t instance) :
      _context(monitor, report, (uint32_t)instance), _env(_context), _instance(instance) {}

  // Verify 0 <= r < modulus
  __device__ FORCE_INLINE void assert_normalized(bn_t &r, const bn_t &modulus) {
    //if (VERIFY_NORMALIZED && _context.check_errors())
    if (VERIFY_NORMALIZED && CHECK_ERROR) {

        // Negative overflow
        if (cgbn_extract_bits_ui32(_env, r, params::BITS-1, 1)) {
            _context.report_error(cgbn_negative_overflow);
        }
        // Positive overflow
        if (cgbn_compare(_env, r, modulus) >= 0) {
            _context.report_error(cgbn_positive_overflow);
        }
    }
  }

  // One step of normalization, sometimes needed after mont_mul
  __device__ FORCE_INLINE void normalize_addition(bn_t &r, const bn_t &modulus) {
      if (cgbn_compare(_env, r, modulus) >= 0) {
          cgbn_sub(_env, r, r, modulus);
      }
  }
};

// kernel implementation using cgbn
template<class params>
__global__ void kernel_pm1_partial(
        cgbn_error_report_t *report,
        uint64_t s_bits,
        uint64_t s_bits_start,
        uint64_t s_bits_interval,
        uint32_t* gpu_s_bits,
        uint32_t *data,
        uint32_t count
        ) {
  // decode an instance_i number from the blockIdx and threadIdx
  int32_t instance_i = (blockIdx.x*blockDim.x + threadIdx.x)/params::TPI;
  if (instance_i >= count)
      return;

  /* Cast uint32_t array to mem_t */
  typename power_mod_t<params>::mem_t *data_cast = (typename power_mod_t<params>::mem_t*) data;
  cgbn_monitor_t monitor = CHECK_ERROR ? cgbn_report_monitor : cgbn_no_checks;
  power_mod_t<params> exp(monitor, report, instance_i);

  typename power_mod_t<params>::bn_t  modulus, const_base, partial_base, partial_result, temp;


  uint32_t np0;

  { // Setup
      /* const_base is used in WINDOW algorithm, partial_base is used in double and square. */
      cgbn_load(exp._env, modulus,        &data_cast[4*instance_i+0]);
      cgbn_load(exp._env, const_base,     &data_cast[4*instance_i+1]);
      cgbn_load(exp._env, partial_base,   &data_cast[4*instance_i+2]);
      cgbn_load(exp._env, partial_result, &data_cast[4*instance_i+3]);

      // Ignore any n=0 that might have been filtered.
      if (cgbn_equals_ui32(exp._env, modulus, 0))
          return;

      /* Convert to mont, has a miniscule bit of overhead with batching. */
      np0 = cgbn_bn2mont(exp._env, const_base, const_base, modulus);
      cgbn_bn2mont(exp._env, partial_base, partial_base, modulus);
      cgbn_bn2mont(exp._env, partial_result, partial_result, modulus);

      exp.assert_normalized(const_base, modulus);
      exp.assert_normalized(partial_base, modulus);
      exp.assert_normalized(partial_result, modulus);
  }

  if (1) {
      /* WINDOW algorithm, MSB to LSB and mixin bits from s in groups of WINDOW_BITS. */
      static const uint32_t WINDOW_BITS = params::WINDOW_BITS;
      typename power_mod_t<params>::bn_local_t  window[1 << WINDOW_BITS];

      cgbn_set_ui32(exp._env, temp, 1);
      cgbn_store(exp._env, window+0, temp);
      cgbn_set(exp._env, temp, const_base);
      cgbn_store(exp._env, window+1, temp);
      #pragma nounroll
      for(uint16_t index = 2; index < (1 << WINDOW_BITS); index++) {
        cgbn_mont_mul(exp._env, temp, temp, const_base, modulus, np0);
        exp.normalize_addition(temp, modulus);
        cgbn_store(exp._env, window+index, temp);
      }

      /* Partial base is modified so that process_result isn't suspicious. */
      cgbn_set_ui32(exp._env, partial_base, 0);

      uint32_t window_bits = 0;
      for (uint64_t b = s_bits_start; b < s_bits_start + s_bits_interval; b++) {
        cgbn_mont_sqr(exp._env, partial_result, partial_result, modulus, np0);
        /* From https://github.com/NVlabs/CGBN/issues/15
         * it's not clear if the redundant representation has consequences
         * when chained */
        exp.normalize_addition(partial_result, modulus);
        exp.assert_normalized(partial_result, modulus);

        /* Process bits from MSB to LSB
         * b counts from 0 to s_num_bits */
        uint64_t nth = s_bits - 1 - b;
        int bit = (gpu_s_bits[nth/32] >> (nth&31)) & 1;

        window_bits = (window_bits << 1) | bit;
        // Can do simplier test that ((window_bits << 1) + 1) <= (1 << WINDOW_BITS)
        if (nth % WINDOW_BITS == 0) {
            // Can skip the mult by 1 case.
            if (window_bits > 0) {
                cgbn_load(exp._env, temp, window + window_bits);
                cgbn_mont_mul(exp._env, partial_result, partial_result, temp, modulus, np0);
                exp.normalize_addition(partial_result, modulus);
                exp.assert_normalized(partial_result, modulus);
            }
            window_bits = 0;
        }
      }

      // If any left over window_bits handle them here.
      if (window_bits) {
        cgbn_load(exp._env, temp, window + window_bits);
        cgbn_mont_mul(exp._env, partial_result, partial_result, temp, modulus, np0);
        exp.normalize_addition(partial_result, modulus);
        exp.assert_normalized(partial_result, modulus);
      }
  } else {
      for (uint64_t b = s_bits_start; b < s_bits_start + s_bits_interval; b++) {
        /* Process bits from LSB to MSB
         * b counts from 0 to s_num_bits */
        uint64_t nth = b;
        int bit = (gpu_s_bits[nth/32] >> (nth&31)) & 1;
        if (bit) {
            cgbn_mont_mul(exp._env, partial_result, partial_result, partial_base, modulus, np0);
            exp.normalize_addition(partial_result, modulus);
            exp.assert_normalized(partial_result, modulus);
        }
        cgbn_mont_sqr(exp._env, partial_base, partial_base, modulus, np0);
        exp.normalize_addition(partial_base, modulus);
        exp.assert_normalized(partial_base, modulus);
      }
  }

  { // Final output
    // Convert everything back to bn
    cgbn_mont2bn(exp._env, partial_base, partial_base, modulus, np0);
    cgbn_mont2bn(exp._env, partial_result, partial_result, modulus, np0);
    exp.assert_normalized(partial_base, modulus);
    exp.assert_normalized(partial_result, modulus);
    cgbn_store(exp._env, &data_cast[4*instance_i+2], partial_base);
    cgbn_store(exp._env, &data_cast[4*instance_i+3], partial_result);
  }
}

static void
set_gpu_pm1_data(
        const mpz_t *x0, const mpz_t *n, uint32_t *datum,
        uint32_t instances, uint32_t BITS, size_t data_size) {
  /**
   * Store 4 numbers per curve:
   * N, Base, partial Base, partial Result
   */
  const size_t limbs_per = BITS/32;
  mpz_t x;
  mpz_init(x);

  for(int index = 0; index < instances; index++, datum += 4 * limbs_per)
    {
      // Some numbers may be zero at end of batch
      if (mpz_sgn(n[index]) == 0) {
        outputf (OUTPUT_VERBOSE, "GPU P-1: skipping %i n=0\n", index);
        // Set base and result to 0 also
        mpz_set_ui (x, 0);
        from_mpz(x, datum + 0 * limbs_per, BITS/32);
        from_mpz(x, datum + 1 * limbs_per, BITS/32);
        from_mpz(x, datum + 2 * limbs_per, BITS/32);
        from_mpz(x, datum + 3 * limbs_per, BITS/32);
        continue;
      }

      /* Modulo (n) */
      from_mpz(n[index], datum + 0 * limbs_per, BITS/32);
      from_mpz(x0[index], datum + 1 * limbs_per, BITS/32);
      from_mpz(x0[index], datum + 2 * limbs_per, BITS/32);
      /* Result set to 1 as sentinel */
      mpz_set_ui (x, 1);
      from_mpz(x, datum + 3 * limbs_per, BITS/32);

      outputf (OUTPUT_TRACE, "GPU P-1: %5d | N=%Zd x0=%Zd\n",
              index, n[index], x0[index]);
    }
  mpz_clear(x);
}


static
int find_pm1_factor(mpz_t factor, const mpz_t N, const mpz_t a) {
    /* Check if factor found <=> gcd(x0^g, N) > 1 */
    mpz_sub_ui (factor, a, 1);
    mpz_gcd (factor, factor, N);

    if (mpz_cmp_ui (factor, 1) > 0) {
        return ECM_FACTOR_FOUND_STEP1;
    }

    mpz_set_ui (factor, 1);
    return ECM_NO_FACTOR_FOUND;
}


/**
 * Read GPU data back.
 * res[i] = x0[i] ^ s mod numbers[i]
 * return 0 on success, ECM_ERROR on failure.
 */
static
int process_pm1_results(
        const mpz_t *x0,
        const mpz_t *numbers,
        mpz_t *factors,
        mpz_t *res,
        const uint32_t *data, uint32_t cgbn_bits,
        int instances) {
  mpz_t modulo, base, result, tmp;
  mpz_init(modulo);
  mpz_init(base);
  mpz_init(result);
  mpz_init(tmp);

  const uint32_t limbs_per = cgbn_bits / 32;
  int errors = 0;
  for(size_t i = 0; i < instances; i++) {
    const uint32_t *datum = data + (4 * i * limbs_per);;

    // Make sure we were testing the right number.
    to_mpz(modulo, datum + 0 * limbs_per, limbs_per);
    assert(mpz_cmp(numbers[i], modulo) == 0);

    to_mpz(base, datum + 2 * limbs_per, limbs_per);
    to_mpz(result, datum + 3 * limbs_per, limbs_per);

    if (test_verbose (OUTPUT_TRACE) && i <= 10) {
      outputf (OUTPUT_TRACE, "index: %lu, modulo: %Zd\n", i, modulo);
      outputf (OUTPUT_TRACE, "index: %lu, base: %Zd\n", i, base);
      outputf (OUTPUT_TRACE, "index: %lu, result: %Zd\n", i, result);
    }

    if (mpz_cmp_ui(modulo, 0) == 0)
      continue;

    /* Very suspicious for base, result to match initial value */
    if (mpz_cmp (base, x0[i]) == 0 && mpz_cmp_ui(result, 1) == 0) {
      errors += 1;
      if (errors < 10 || errors % 100 == 1)
        outputf (OUTPUT_ERROR, "GPU: curve %d n=%Zd\n", i, modulo);
        outputf (OUTPUT_ERROR, "GPU: curve %d didn't compute or found input?\n", i);
    }

    mpz_set(res[i], result);
    mpz_set(tmp, factors[i]);
    int found = find_pm1_factor(factors[i], modulo, result);
    // Don't print already found factors again.
    if (found != ECM_NO_FACTOR_FOUND && mpz_cmp(factors[i], tmp) != 0) {
      outputf (OUTPUT_NORMAL, "GPU P-1: factor %Zd found in Step 1 with curve %ld\n",
          factors[i], i);
      outputf (OUTPUT_VERBOSE, "Input number is %Zd\n", modulo);
      outputf (OUTPUT_VERBOSE, "********** Factor found in step 1: %Zd\n", factors[i]);
    }
  }

  mpz_clear(modulo);
  mpz_clear(base);
  mpz_clear(result);
  mpz_clear(tmp);

#ifdef IS_DEV_BUILD
  if (errors)
        outputf (OUTPUT_ERROR, "Had %d errors. Try `make clean; make` or reducing TPB_DEFAULT\n",
            errors);
#endif

  if (errors > 2)
      return ECM_ERROR;

  return 0;
}

/* Similiar to writechkfile for batch for results */
void
writechkfile_batch (
     char *chkfilename, const ecm_params params,
     uint32_t instances, const mpz_t *n, const mpz_t *x, const mpz_t *orig_x0)
{
  FILE *chkfile;
  char comment[] = "GPU P-1 checkpoint";

  mpcandi_t tmp_n;
  mpcandi_t_init(&tmp_n);

  outputf (OUTPUT_DEVVERBOSE, "Writing GPU checkpoint to %s at B1 = %" PRIu64 "\n",
           chkfilename, params->B1done);

  chkfile = fopen (chkfilename, "w");
  ASSERT_ALWAYS(chkfile != NULL);

  for (size_t i = 0; i < instances; i++)
    {
      mpz_set(tmp_n.n, n[i]);
      write_resumefile_line(
              chkfile, ECM_PM1, params->B1done,
              /*sigma=*/ NULL, /*sigma_is_A=*/ 0, /*Etype=*/ 0, /*param=*/ 0,
              /*x=*/ x[i], /*y=*/ NULL, /*n=*/ &tmp_n,
              /*x0=*/ orig_x0[i], /*y0=*/ NULL,
              comment);
    }

  fflush (chkfile);
  fclose (chkfile);
  mpcandi_t_free(&tmp_n);
}


/* Performs P-1 from B1done-B1 on many numbers at the same time.
 * Params:
 *   n: numbers
 *   x0: x0 for computation
 *   orig_x0: orig_x0 for resume files.
 *   factors: any factor found during stage 1
 *   res: at return contains x = x0 ^ S, S = lcm(prod(B1)
 * Returns ECM_ERROR if stop_asap or problem during execution
 */
int cgbn_pm1_stage1(
     const mpz_t *n, const mpz_t *x0, const mpz_t *orig_x0,
     mpz_t *factors, mpz_t *res,
     const uint64_t B1, uint32_t instances,
     float *gputime, const ecm_params params, ecm_params mutable_params)
{
  int verbose = params->verbose;

  /* Choose batch so that cascade isn't huge, but also minizing number
     of Host->GPU memcpys */
  uint64_t B1_incr = MIN(B1, 2000000);

  /* Batches are 1.44 bits * B1_incr, allocate 20% extra. */
  const uint64_t s_buffer_bytes = ((B1_incr * 216 / 100 + 127) / 128) * 128;
  uint32_t* const s_bits = (uint32_t*) malloc (s_buffer_bytes);
  assert( s_bits != NULL );
  uint32_t *gpu_s_bits;
  CUDA_CHECK(cudaMalloc((void **)&gpu_s_bits, s_buffer_bytes));
  assert( gpu_s_bits != NULL );

  outputf (OUTPUT_DEVVERBOSE, "GPU P-1: B1 increment = %'lu, allocated %d MB\n",
          B1_incr, s_buffer_bytes >> 20);

  cudaEvent_t global_start, batch_start, stop;
  CUDA_CHECK(cudaEventCreate (&global_start));
  CUDA_CHECK(cudaEventCreate (&batch_start));
  CUDA_CHECK(cudaEventCreate (&stop));
  CUDA_CHECK(cudaEventRecord (global_start));


  cgbn_error_report_t *report;
  // create a cgbn_error_report for CGBN to report back errors
  CUDA_CHECK(cgbn_error_report_alloc(&report));

  // Find Largest N
  size_t n_log2 = 0;
  size_t largest_i = 0;
  for (size_t i = 1; i < instances; i++)
      if (mpz_cmp(n[i], n[largest_i]) > 0)
          largest_i = i;

  n_log2 = mpz_sizeinbase(n[largest_i], 2);
  assert(n_log2 > 1);
  outputf (OUTPUT_VERBOSE, "GPU P-1: Largest number line %lu, %lu bits\n",
           largest_i+1, n_log2);

  uint32_t  BITS = 0;        // kernel bits
  int32_t   TPB=TPB_DEFAULT; // Always the same default
  int32_t   TPI;
  int32_t   IPB;             // IPB = TPB / TPI, instances per block
  size_t    BLOCK_COUNT;     // How many blocks to cover all instances

  /**
   * Smaller TPI is faster, Larger TPI is needed for large inputs.
   * N > 512 TPI=8 | N > 2048 TPI=16 | N > 8192 TPI=32
   *
   * Larger takes longer to compile (and increases binary size)
   * No GPU, No CGBN | ecm 3.4M, 2 seconds to compile
   * GPU, No CGBN    | ecm 3.5M, 3 seconds
   * (8, 1024)       | ecm 3.8M, 12 seconds
   * (16,8192)       | ecm 4.2M, 1 minute
   * (32,16384)      | ecm 4.2M, 1 minute
   * (32,32768)      | ecm 5.2M, 4.7 minutes
   */
  /* NOTE: Custom kernel changes here
   * For "Compiling custom kernel for %d bits should be XX% faster"
   * Change the 512 in cgbn_params_t<4, 512> cgbn_params_small;
   * to the suggested value (a multiple of 32 >= bits + 6).
   * You may need to change the 4 to an 8 (or 16) if bits >512, >2048
   */
  /** TODO: try with const vector for BITs/TPI, see if compiler is happy */
  std::vector<uint32_t> available_kernels;

  /* These are P-1 kernels, don't mistake them for ECM kernels */
  typedef cgbn_params_t<8, 768>   cgbn_params_small;
  typedef cgbn_params_t<8, 1024>  cgbn_params_medium;
  available_kernels.push_back((uint32_t)cgbn_params_small::BITS);
  available_kernels.push_back((uint32_t)cgbn_params_medium::BITS);

#ifndef IS_DEV_BUILD
  /**
   * TPI and BITS have to be set at compile time. Adding multiple cgbn_params
   * (and their associated kernels) allows for better dynamic selection based
   * on the size of N (e.g. N < 1024, N < 2048, N < 4096) but increase compile
   * time and binary size. A few reasonable sizes are included and a verbose
   * warning is printed when a particular N might benefit from a custom sized
   * kernel.
   */
  typedef cgbn_params_t<8, 1536>  cgbn_params_1536;
  typedef cgbn_params_t<8, 2048>  cgbn_params_2048;
  typedef cgbn_params_t<16, 3072> cgbn_params_3072;
  typedef cgbn_params_t<16, 4096> cgbn_params_4096;
  available_kernels.push_back((uint32_t)cgbn_params_1536::BITS);
  available_kernels.push_back((uint32_t)cgbn_params_2048::BITS);
  available_kernels.push_back((uint32_t)cgbn_params_3072::BITS);
  available_kernels.push_back((uint32_t)cgbn_params_4096::BITS);
#endif

  for (int k_i = 0; k_i < available_kernels.size(); k_i++) {
    uint32_t kernel_bits = available_kernels[k_i];
    if (kernel_bits >= n_log2 + CARRY_BITS) {
      BITS = kernel_bits;
      assert( BITS % 32 == 0 );

      /* Print some debug info about kernel. */
      /* TODO: return kernelAttr and validate maxThreadsPerBlock. */
      if (BITS == cgbn_params_small::BITS) {
        TPI = cgbn_params_small::TPI;
        kernel_info((const void*)kernel_pm1_partial<cgbn_params_small>, verbose);
      } else if (BITS == cgbn_params_medium::BITS) {
        TPI = cgbn_params_medium::TPI;
        kernel_info((const void*)kernel_pm1_partial<cgbn_params_medium>, verbose);
#ifndef IS_DEV_BUILD
      } else if (BITS == cgbn_params_1536::BITS) {
        TPI = cgbn_params_1536::TPI;
        kernel_info((const void*)kernel_pm1_partial<cgbn_params_1536>, verbose);
      } else if (BITS == cgbn_params_2048::BITS) {
        TPI = cgbn_params_2048::TPI;
        kernel_info((const void*)kernel_pm1_partial<cgbn_params_2048>, verbose);
      } else if (BITS == cgbn_params_3072::BITS) {
        TPI = cgbn_params_3072::TPI;
        kernel_info((const void*)kernel_pm1_partial<cgbn_params_3072>, verbose);
      } else if (BITS == cgbn_params_4096::BITS) {
        TPI = cgbn_params_4096::TPI;
        kernel_info((const void*)kernel_pm1_partial<cgbn_params_4096>, verbose);
#endif
      } else {
        /* lowercase k to help differentiate this error from one below */
        outputf (OUTPUT_ERROR, "CGBN kernel not found for %d bits\n", BITS);
        return ECM_ERROR;
      }

      IPB = TPB / TPI;
      BLOCK_COUNT = (instances + IPB - 1) / IPB;

      break;
    }
  }

  if (BITS == 0) {
    outputf (OUTPUT_ERROR, "No available CGBN Kernel large enough to process N(%d bits)\n", n_log2);
    return ECM_ERROR;
  }

  /* Alert that recompiling with a smaller kernel would likely improve speed */
  {
    size_t optimized_bits = ((n_log2 + CARRY_BITS + 127)/128) * 128;
    /* Assume speed is roughly O(N) but slightly slower for not being a power of two */
    float pct_faster = 90 * BITS / optimized_bits;

    if (pct_faster > 110) {
      outputf (OUTPUT_VERBOSE, "Compiling custom kernel for %d bits should be ~%.0f%% faster see README.gpu\n",
              optimized_bits, pct_faster);
    }
  }

  int youpi = verify_size_of_n(n[largest_i], BITS);
  if (youpi != ECM_NO_FACTOR_FOUND) {
    return youpi;
  }

  /* Consistency check that struct cgbn_mem_t is byte aligned without extra fields. */
  assert( sizeof(power_mod_t<cgbn_params_small>::mem_t) == cgbn_params_small::BITS/8 );
  assert( sizeof(power_mod_t<cgbn_params_medium>::mem_t) == cgbn_params_medium::BITS/8 );

  outputf (OUTPUT_VERBOSE,
          "CGBN<%d, %d> running kernel<%d block x %d threads> input number is %d bits\n",
          BITS, TPI, BLOCK_COUNT, TPB, n_log2);

  /* Set up CPU and GPU memory for data */
  const size_t limbs_per = BITS/32;
  const size_t data_size = 4 * instances * limbs_per * sizeof(uint32_t);
  uint32_t* const data = (uint32_t*) malloc(data_size);
  uint32_t *gpu_data;
  CUDA_CHECK(cudaMalloc((void **)&gpu_data, data_size));

  outputf (OUTPUT_DEVVERBOSE,
          "GPU P-1: data is %'lu bytes, S is %'lu bytes\n",
          data_size, s_buffer_bytes);

  /* For eta */
  uint64_t B1start = params->B1done;
  /* For computing S in batches. */
  uint64_t B1_partial = (params->B1done) < 2.0 ? 0 : params->B1done;
  mpz_t s;
  mpz_init(s);
  batched_info_t batched;
  batched_info_init(batched);
  advance_to(batched, params->B1done);

  /* Start with small batches and increase till timing is BATCH_MS_GOAL */
  uint64_t batch_size = 2000;
  int batches_complete = 0;
  /* gputime and batch_time are measured in ms */
  float batch_time = 0;
  long last_chkpnt_time = time(0);
  const long GPU_CHKPNT_PERIOD = 180;

  while (B1_partial < B1)
    {
      /* This isn't working the second time around */
      set_gpu_pm1_data(
        batches_complete == 0 ? x0 : res,
        n, data, instances, BITS, data_size);

      CUDA_CHECK(cudaMemcpy(gpu_data, data, data_size, cudaMemcpyHostToDevice));

      /* Set up S */
      B1_partial = std::min<uint64_t>(B1_partial + B1_incr, B1);
      get_batch(batched, s, B1_partial);
      uint64_t s_num_bits = mpz_sizeinbase(s, 2);
      set_s_bits_copy_to_gpu(s, s_buffer_bytes, s_bits, gpu_s_bits);
      outputf (OUTPUT_TRACE, "GPU P-1: Processing B1=%lu-%lu\n", params->B1done+1, B1_partial);

      /* Synchronize, waiting for memCpy to GPU to finish */
      CUDA_CHECK(cudaDeviceSynchronize());

      uint64_t s_partial = 0;
      while (s_partial < s_num_bits)
        {
          if (params->stop_asap != NULL && params->stop_asap())
            {
              outputf (OUTPUT_ERROR, "GPU P-1: Exiting from stop_asap()\n");
              if (batches_complete < 1000)
                youpi = ECM_ERROR;
              goto cleanup;
            }

          /* decrease batch_size for final batch if needed */
          uint64_t remaining = s_num_bits - s_partial;
          if (remaining < 2 * batch_size) {
              batch_size = remaining;
          }

          /* print with decreasing frequency, 3 + 3/1s + 3/10s + 5/100s + every 1000s */
          if (print_nth_batch (batches_complete)) {
            double batch_partial = (double) s_partial / s_num_bits * (B1_partial - params->B1done);
            double fraction = ((double) params->B1done + batch_partial - B1start) / (B1 - B1start);
            outputf (OUTPUT_VERBOSE,
                "GPU P-1: B1=%.0f-%lu, %lu/%lu batch: %lu (%.1f%%)",
                params->B1done+1, B1_partial,
                s_partial, s_num_bits, batch_size,
                100.0 * fraction);

            if (*gputime < 1000 || batches_complete < 2) {
              outputf (OUTPUT_VERBOSE, "\n");
            } else {
              /* TODO need to account for any B1done that was passed in */
              float estimated_total = (*gputime) / fraction;
              float eta = estimated_total - (*gputime);
              outputf (OUTPUT_VERBOSE,
                      ", ETA %.f + %.f = %.f seconds (~%.f ms/instance)\n",
                      eta / 1000, *gputime / 1000, estimated_total / 1000,
                      estimated_total / instances);
            }
          }

          CUDA_CHECK(cudaEventRecord (batch_start));

          if (BITS == cgbn_params_small::BITS) {
            kernel_pm1_partial<cgbn_params_small><<<BLOCK_COUNT, TPB>>>(
                report, s_num_bits, s_partial, batch_size, gpu_s_bits, gpu_data, instances);
          } else if (BITS == cgbn_params_medium::BITS) {
            kernel_pm1_partial<cgbn_params_medium><<<BLOCK_COUNT, TPB>>>(
                report, s_num_bits, s_partial, batch_size, gpu_s_bits, gpu_data, instances);
#ifndef IS_DEV_BUILD
          } else if (BITS == cgbn_params_1536::BITS) {
            kernel_pm1_partial<cgbn_params_1536><<<BLOCK_COUNT, TPB>>>(
                report, s_num_bits, s_partial, batch_size, gpu_s_bits, gpu_data, instances);
          } else if (BITS == cgbn_params_2048::BITS) {
            kernel_pm1_partial<cgbn_params_2048><<<BLOCK_COUNT, TPB>>>(
                report, s_num_bits, s_partial, batch_size, gpu_s_bits, gpu_data, instances);
          } else if (BITS == cgbn_params_3072::BITS) {
            kernel_pm1_partial<cgbn_params_3072><<<BLOCK_COUNT, TPB>>>(
                report, s_num_bits, s_partial, batch_size, gpu_s_bits, gpu_data, instances);
          } else if (BITS == cgbn_params_4096::BITS) {
            kernel_pm1_partial<cgbn_params_4096><<<BLOCK_COUNT, TPB>>>(
                report, s_num_bits, s_partial, batch_size, gpu_s_bits, gpu_data, instances);
#endif
          } else {
            outputf (OUTPUT_ERROR, "CGBN P-1 Kernel not found for %d bits\n", BITS);
            return ECM_ERROR;
          }

          s_partial += batch_size;
          batches_complete++;

          /* error report uses managed memory, sync the device and check for cgbn errors */
          CUDA_CHECK(cudaDeviceSynchronize());
          if (report->_error)
            outputf (OUTPUT_ERROR, "\n\nerror: %d\n", report->_error);
          CGBN_CHECK(report);

          CUDA_CHECK(cudaEventRecord (stop));
          CUDA_CHECK(cudaEventSynchronize (stop));
          cudaEventElapsedTime (&batch_time, batch_start, stop);
          cudaEventElapsedTime (gputime, global_start, stop);

          /* Adjust batch_size to aim for BATCH_MS_GOAL. */
          if (10 * batch_time < 8 * BATCH_MS_GOAL ) {
            batch_size = 11*batch_size/10 + 1;
          } else if (8 * batch_time > 10 * BATCH_MS_GOAL) {
            batch_size = max(100ul, 9*batch_size / 10);
          }
        }

      /* After S is finished copy data back, possible take checkpoint */
      CUDA_CHECK(cudaMemcpy(data, gpu_data, data_size, cudaMemcpyDeviceToHost));

      /* Get x (res) back from data */
      if (process_pm1_results(x0, n, factors, res, data, BITS, instances) == ECM_ERROR)
          goto cleanup;

      /* Update B1done */
      mutable_params->B1done = B1_partial;

      if (params->chkfilename != NULL &&
          (time(0) - last_chkpnt_time) > GPU_CHKPNT_PERIOD)
        {
          writechkfile_batch (params->chkfilename, params, instances, n, res, orig_x0);
          last_chkpnt_time = time(0);
        }
    }

cleanup:

  // clean up
  CUDA_CHECK(cudaFree(gpu_s_bits));
  CUDA_CHECK(cudaFree(gpu_data));
  CUDA_CHECK(cgbn_error_report_free(report));
  CUDA_CHECK(cudaEventDestroy (global_start));
  CUDA_CHECK(cudaEventDestroy (batch_start));
  CUDA_CHECK(cudaEventDestroy (stop));

  free(s_bits);
  free(data);

  return youpi;
}



#ifdef __CUDA_ARCH__
  #if __CUDA_ARCH__ < 350
    #error "Unsupported architecture"
  #endif
#endif

#endif  /* _GPU_PM1_CU */
