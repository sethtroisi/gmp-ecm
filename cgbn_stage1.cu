/* cgbn_stage1.h: header for CGBN (GPU) based ecm stage 1.

Copyright 2021 Seth Troisi

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

#ifndef _CGBN_STAGE1_CU
#define _CGBN_STAGE1_CU 1

#include "cgbn_stage1.h"

#include <cassert>
#include <stdio.h>
#include <stdlib.h>

// GMP import must proceed cgbn.h
#include <gmp.h>
#include <cgbn.h>
#include <cuda.h>

#include "ecm.h"
#include "ecm-gpu.h"

void cuda_check(cudaError_t status, const char *action=NULL, const char *file=NULL, int32_t line=0) {
  // check for cuda errors
  if (status!=cudaSuccess) {
    fprintf (stderr, "CUDA error occurred: %s\n", cudaGetErrorString(status));
    if (action!=NULL)
      fprintf (stderr, "While running %s   (file %s, line %d)\n", action, file, line);
    exit(1);
  }
}

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

// Unify this with cudakernel.cu
#define CUDA_CHECK(action) cuda_check(action, #action, __FILE__, __LINE__)
#define CGBN_CHECK(report) cgbn_check(report, __FILE__, __LINE__)

void to_mpz(mpz_t r, const uint32_t *x, uint32_t count) {
  mpz_import (r, count, -1, sizeof(uint32_t), 0, 0, x);
}

void from_mpz(const mpz_t s, uint32_t *x, uint32_t count) {
  size_t words;

  if(mpz_sizeinbase (s, 2) > count * 32) {
    fprintf (stderr, "from_mpz failed -- result does not fit\n");
    exit(1);
  }

  mpz_export (x, &words, -1, sizeof(uint32_t), 0, 0, s);
  while(words<count)
    x[words++]=0;
}



// ---------------------------------------------------------------- //

// See cgbn_error_t enum (cgbn.h:39)
#define cgbn_normalized_error ((cgbn_error_t) 14)
#define cgbn_positive_overflow ((cgbn_error_t) 15)
#define cgbn_negative_overflow ((cgbn_error_t) 16)

// Seems to adds very small overhead (1-10%)
#define VERIFY_NORMALIZED 1
// Adds even less overhead (<1%)
#define CHECK_ERROR 1

// Can dramatically change compile time
#if 1
    #define FORCE_INLINE __forceinline__
#else
    #define FORCE_INLINE
#endif


// The CGBN context uses the following three parameters:
//   TBP             - threads per block (zero means to use the blockDim.x)
//   MAX_ROTATION    - must be small power of 2, imperically, 4 works well
//   SHM_LIMIT       - number of bytes of dynamic shared memory available to the kernel
//   CONSTANT_TIME   - require constant time algorithms (currently, constant time algorithms are not available)

// Locally it will also be helpful to have several parameters:
//   TPI             - threads per instance
//   BITS            - number of bits per instance

template<uint32_t tpi, uint32_t bits>
class cgbn_params_t {
  public:
  // parameters used by the CGBN context
  static const uint32_t TPB=0;                     // get TPB from blockDim.x
  static const uint32_t MAX_ROTATION=4;            // good default value
  static const uint32_t SHM_LIMIT=0;               // no shared mem available
  static const bool     CONSTANT_TIME=false;       // not implemented

  // parameters used locally in the application
  static const uint32_t TPI=tpi;                   // threads per instance
  static const uint32_t BITS=bits;                 // instance size
};


template<class params>
class curve_t {
  public:

  // define the instance structure
  typedef struct {
    cgbn_mem_t<params::BITS> aX;
    cgbn_mem_t<params::BITS> aY;
    cgbn_mem_t<params::BITS> bX;
    cgbn_mem_t<params::BITS> bY;
    cgbn_mem_t<params::BITS> modulus;
    uint32_t d;
  } instance_t;


  typedef cgbn_context_t<params::TPI, params>   context_t;
  typedef cgbn_env_t<context_t, params::BITS>   env_t;
  typedef typename env_t::cgbn_t                bn_t;

  context_t _context;
  env_t     _env;
  int32_t   _instance; // which curve instance is this

  // Constructor
  __device__ FORCE_INLINE curve_t(cgbn_monitor_t monitor, cgbn_error_report_t *report, int32_t instance) :
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

  // Normalize after addition
  __device__ FORCE_INLINE void normalize_addition(bn_t &r, const bn_t &modulus) {
      if (cgbn_compare(_env, r, modulus) >= 0) {
          cgbn_sub(_env, r, r, modulus);
      }
  }

  // Normalize after subtraction (handled instead by checking carry)
  /*
  __device__ FORCE_INLINE void normalize_subtraction(bn_t &r, const bn_t &modulus) {
      if (cgbn_extract_bits_ui32(_env, r, params::BITS-1, 1)) {
          cgbn_add(_env, r, r, modulus);
      }
  }
  */

  /**
   * Calculate (r * m) / 2^32 mod modulus
   *
   * This removes a factor of 2^32 which is not present in m.
   * Otherwise m (really d) needs to be passed as a bigint not a uint32
   */
  __device__ FORCE_INLINE void special_mult_ui32(bn_t &r, uint32_t m, const bn_t &modulus, uint32_t np0) {
    //uint32_t thread_i = (blockIdx.x*blockDim.x + threadIdx.x)%params::TPI;
    bn_t temp;

    uint32_t carry_t1 = cgbn_mul_ui32(_env, r, r, m);
    uint32_t t1_0 = cgbn_extract_bits_ui32(_env, r, 0, 32);
    uint32_t q = t1_0 * np0;
    uint32_t carry_t2 = cgbn_mul_ui32(_env, temp, modulus, q);

    cgbn_shift_right(_env, r, r, 32);
    cgbn_shift_right(_env, temp, temp, 32);
    // Add back overflow carry
    cgbn_insert_bits_ui32(_env, r, r, params::BITS-32, 32, carry_t1);
    cgbn_insert_bits_ui32(_env, temp, temp, params::BITS-32, 32, carry_t2);

    // This needs to be measured at block containing top bit of modulus
    int32_t carry_q = cgbn_add(_env, r, r, temp);
    carry_q += cgbn_add_ui32(_env, r, r, t1_0 != 0); // add 1
    while (carry_q != 0) {
        carry_q -= cgbn_sub(_env, r, r, modulus);
    }

    // 0 <= r, temp < modulus => r + temp + 1 < 2*modulus
    if (cgbn_compare(_env, r, modulus) >= 0) {
        cgbn_sub(_env, r, r, modulus);
    }
  }


  __device__ FORCE_INLINE void double_add_v2(
          bn_t &q, bn_t &u,
          bn_t &w, bn_t &v,
          uint32_t d,
          uint32_t bit_number,
          const bn_t &modulus,
          const uint32_t np0) {
    // q = xA = aX
    // u = zA = aY
    // w = xB = bX
    // v = zB = bY

    /* Doesn't seem to be a large cost to using many extra variables */
    bn_t t, CB, DA, AA, BB, K, dK;

    /* Can maybe use one more bit if cgbn_add subtracts when carry happens */

    cgbn_add(_env, t, v, w); // t = (bY + bX)
    normalize_addition(t, modulus);
    if (cgbn_sub(_env, v, v, w)) // v = (bY - bX)
        cgbn_add(_env, v, v, modulus);


    cgbn_add(_env, w, u, q); // w = (aY + aX)
    normalize_addition(w, modulus);
    if (cgbn_sub(_env, u, u, q)) // u = (aY - aX)
        cgbn_add(_env, u, u, modulus);
    if (VERIFY_NORMALIZED) {
        assert_normalized(t, modulus);
        assert_normalized(v, modulus);
        assert_normalized(w, modulus);
        assert_normalized(u, modulus);
    }

    cgbn_mont_mul(_env, CB, t, u, modulus, np0); // C*B
        normalize_addition(CB, modulus);
    cgbn_mont_mul(_env, DA, v, w, modulus, np0); // D*A
        normalize_addition(DA, modulus);

    cgbn_mont_sqr(_env, AA, w, modulus, np0);    // AA
    cgbn_mont_sqr(_env, BB, u, modulus, np0);    // BB
    normalize_addition(AA, modulus);
    normalize_addition(BB, modulus);
    if (VERIFY_NORMALIZED) {
        assert_normalized(CB, modulus);
        assert_normalized(DA, modulus);
        assert_normalized(AA, modulus);
        assert_normalized(BB, modulus);
    }

    // q = aX is finalized
    cgbn_mont_mul(_env, q, AA, BB, modulus, np0); // AA*BB
        normalize_addition(q, modulus);
        assert_normalized(q, modulus);

    if (cgbn_sub(_env, K, AA, BB)) // K = AA-BB
        cgbn_add(_env, K, K, modulus);

    // By definition of d = (sigma / 2^32) % MODN
    // K = k*R
    // dK = d*k*R = (K * R * sigma) >> 32
    cgbn_set(_env, dK, K);
    special_mult_ui32(dK, d, modulus, np0); // dK = K*d
        assert_normalized(dK, modulus);

    cgbn_add(_env, u, BB, dK); // BB + dK
    normalize_addition(u, modulus);
    if (VERIFY_NORMALIZED) {
        assert_normalized(K, modulus);
        assert_normalized(dK, modulus);
        assert_normalized(u, modulus);
    }

    // u = aY is finalized
    cgbn_mont_mul(_env, u, K, u, modulus, np0); // K(BB+dK)
        normalize_addition(u, modulus);
        assert_normalized(u, modulus);

    cgbn_add(_env, w, DA, CB); // DA + CB
    normalize_addition(w, modulus);
    if (cgbn_sub(_env, v, DA, CB)) // DA - CB
        cgbn_add(_env, v, v, modulus);
    if (VERIFY_NORMALIZED) {
        assert_normalized(w, modulus);
        assert_normalized(v, modulus);
    }

    // w = bX is finalized
    cgbn_mont_sqr(_env, w, w, modulus, np0); // (DA+CB)^2 mod N
        normalize_addition(w, modulus);
        assert_normalized(w, modulus);

    cgbn_mont_sqr(_env, v, v, modulus, np0); // (DA-CB)^2 mod N
        normalize_addition(v, modulus);
        assert_normalized(v, modulus);

    // v = bY is finalized
    cgbn_shift_left(_env, v, v, 1); // double
    normalize_addition(v, modulus);
        assert_normalized(v, modulus);
  }

  __host__ static instance_t *generate_instances(const mpz_t N, const mpz_t s, ecm_params_t &ecm_params) {
    // P1_x, P1_y = (2,1)
    // 2P_x, 2P_y = (9, 64 * d + 8)

    assert(ecm_params.sigma > 0);

    ecm_params.n_log2 = mpz_sizeinbase(N, 2);
    if (ecm_params.n_log2 + 2 >= params::BITS) {
        printf("N(%d bits)+carry >= BITS %d\n", ecm_params.n_log2, params::BITS);
        exit(1);
    }

    if (ecm_params.n_log2 + 6 > params::BITS) {
        printf("N(%d bits)+carry ~ BITS %d\n", ecm_params.n_log2, params::BITS);
        printf("being caution, feel free to disable this check\n");
        printf("if you do disable, probably should verify a result against gmp-ecm\n");
        exit(1);
    }

    uint32_t num_bits = mpz_sizeinbase(s, 2) - 1;

    assert( 1 <= num_bits <= 1e8 );
    ecm_params.num_bits = num_bits;
    // Use int* so that size can be stored in first element, could pass around extra size.
    ecm_params.s_bits = (char*) malloc(sizeof(char) * num_bits);
    for (int i = 0; i < num_bits; i++) {
        ecm_params.s_bits[i] = mpz_tstbit (s, num_bits - 1 - i);
    }

    instance_t *instances=(instance_t *)malloc(sizeof(instance_t)*ecm_params.curves);

    mpz_t x;
    mpz_init(x);
    for(int index=0;index<ecm_params.curves;index++) {
        instance_t &instance = instances[index];

        // d = (sigma / 2^32) mod N BUT 2^32 handled by special_mul_ui32
        instance.d = ecm_params.sigma + index;

        // mod
        from_mpz(N, instance.modulus._limbs, params::BITS/32);

        // P1 (X, Y)
        mpz_set_ui(x, 2);
        from_mpz(x, instance.aX._limbs, params::BITS/32);
        mpz_set_ui(x, 1);
        from_mpz(x, instance.aY._limbs, params::BITS/32);

        // 2P = P2 (X, Y)
        // P2_y = 64 * d + 8
        mpz_set_ui(x, 9);
        from_mpz(x, instance.bX._limbs, params::BITS/32);

        // d = sigma * mod_inverse(2 ** 32, N)
        mpz_ui_pow_ui(x, 2, 32);
        mpz_invert(x, x, N);
        mpz_mul_ui(x, x, instance.d);
        // P2_y = 64 * d - 2;
        mpz_mul_ui(x, x, 64);
        mpz_add_ui(x, x, 8);
        mpz_mod(x, x, N);

        outputf (OUTPUT_TRACE, "sigma %d => P2_y: %Zd\n", instance.d, x);
        from_mpz(x, instance.bY._limbs, params::BITS/32);
    }

    mpz_clear(x);
    return instances;
  }
};

// kernel implementation using cgbn
template<class params>
__global__ void kernel_double_add(
        cgbn_error_report_t *report,
        uint32_t num_bits,
        char* gpu_s_bits,
        typename curve_t<params>::instance_t *instances,
        uint32_t count) {

  // decode an instance_i number from the blockIdx and threadIdx
  int32_t instance_i=(blockIdx.x*blockDim.x + threadIdx.x)/params::TPI;
  int32_t instance_j=(blockIdx.x*blockDim.x + threadIdx.x)%params::TPI;
  if (instance_j == -123) return;   // avoid unused warning
  if(instance_i >= count)
    return;

  cgbn_monitor_t monitor = CHECK_ERROR ? cgbn_report_monitor : cgbn_no_checks;

  curve_t<params> curve(monitor, report, instance_i);
  typename curve_t<params>::bn_t  aX, aY, bX, bY, modulus;
  typename curve_t<params>::instance_t &instance = instances[instance_i];

  uint32_t np0;
  { // Setup
      cgbn_load(curve._env, aX, &(instance.aX));
      cgbn_load(curve._env, aY, &(instance.aY));
      cgbn_load(curve._env, bX, &(instance.bX));
      cgbn_load(curve._env, bY, &(instance.bY));
      cgbn_load(curve._env, modulus, &(instance.modulus));

      // Convert everything to mont
      np0 = cgbn_bn2mont(curve._env, aX, aX, modulus);
      cgbn_bn2mont(curve._env, aY, aY, modulus);
      cgbn_bn2mont(curve._env, bX, bX, modulus);
      cgbn_bn2mont(curve._env, bY, bY, modulus);

      {
        curve.assert_normalized(aX, modulus);
        curve.assert_normalized(aY, modulus);
        curve.assert_normalized(bX, modulus);
        curve.assert_normalized(bY, modulus);
      }
  }

  for (int b = num_bits; b > 0; b--) {
    /**
     * TODO generates a lot of duplicate inlined code, not sure how to improve
     * Tried with swappings pointers (with a single call to double_add_v2)
     */
    if (gpu_s_bits[num_bits - b] == 0) {
        curve.double_add_v2(aX, aY, bX, bY, instance.d, b, modulus, np0);
    } else {
        curve.double_add_v2(bX, bY, aX, aY, instance.d, b, modulus, np0);
    }
  }

  { // Final output
    // Convert everything back to bn
    cgbn_mont2bn(curve._env, aX, aX, modulus, np0);
    cgbn_mont2bn(curve._env, aY, aY, modulus, np0);
    cgbn_mont2bn(curve._env, bX, bX, modulus, np0);
    cgbn_mont2bn(curve._env, bY, bY, modulus, np0);
    {
      curve.assert_normalized(aX, modulus);
      curve.assert_normalized(aY, modulus);
      curve.assert_normalized(bX, modulus);
      curve.assert_normalized(bY, modulus);
    }
    cgbn_store(curve._env, &(instance.aX), aX);
    cgbn_store(curve._env, &(instance.aY), aY);
    cgbn_store(curve._env, &(instance.bX), bX);
    cgbn_store(curve._env, &(instance.bY), bY);
  }
}


int findfactor(mpz_t factor, const mpz_t N, const mpz_t x_final, const mpz_t y_final) {
    // XXX: combine / refactor logic with cudawrapper.c findfactor

    mpz_t temp;
    mpz_init(temp);

    // Check if factor found

    bool inverted = mpz_invert(temp, y_final, N);    // aY ^ (N-2) % N
    if (inverted) {
        mpz_mul(temp, x_final, temp);         // aX * aY^-1
        mpz_mod(factor, temp, N);             // "Residual"
        mpz_clear(temp);
        return ECM_NO_FACTOR_FOUND;
    }
    mpz_clear(temp);

    mpz_gcd(factor, y_final, N);
    return ECM_FACTOR_FOUND_STEP1;
}

int run_cgbn(mpz_t *factors, int *array_stage_found,
             const mpz_t N, const mpz_t s, float *gputime,
             ecm_params_t *ecm_params) {
  // TODO figure out how to configure this.
  /**
   * TPI and BITS have to be set at compile time.
   * Adding multiple cgbn_params (and kernals) allows for dynamic selection
   * based on the size of N (e.g. N < 1024, N < 2048, N < 4096). If we go that
   * route it would be helpful to do that only during release builds.
   */

  /**
   * Smaller TPI (e.g. 8) is faster (TPI=4 seems worse than TPI=8).
   * Larger TPI (e.g. 32) for running a single curve (or large N).
   * TPI=8  is required for N > 512
   * TPI=16 is required for N > 2048
   * TPI=32 is required for N > 8192
   */
  typedef cgbn_params_t<8, 1024> cgbn_params;

  typedef typename curve_t<cgbn_params>::instance_t instance_t;

  char                *gpu_s_bits;
  instance_t          *instances, *gpu_instances;
  size_t               instance_size = sizeof(instance_t) * ecm_params->curves;
  cgbn_error_report_t *report;
  int32_t              TPB=(cgbn_params::TPB==0) ? 128 : cgbn_params::TPB;    // default threads per block to 128
  int32_t              TPI=cgbn_params::TPI,
                       IPB=TPB/TPI;                                 // IPB is instances per block

  // Keeps CPU from busy waiting during GPU execution.
  CUDA_CHECK(cudaSetDeviceFlags (cudaDeviceScheduleBlockingSync));
  cudaEvent_t start, stop;
  CUDA_CHECK(cudaEventCreate (&start));
  CUDA_CHECK(cudaEventCreate (&stop));
  CUDA_CHECK(cudaEventRecord (start, 0));


  // TODO can IPB / TPB be reduced when curves = 1 for faster execution?
  size_t gpu_block_count = (ecm_params->curves+IPB-1)/IPB;

  instances = curve_t<cgbn_params>::generate_instances(N, s, *ecm_params);
  assert(ecm_params->s_bits != NULL);

  CUDA_CHECK(cudaSetDevice(0));
  // Copy s_bits
  CUDA_CHECK(cudaMalloc((void **)&gpu_s_bits, sizeof(char) * ecm_params->num_bits));
  CUDA_CHECK(cudaMemcpy(gpu_s_bits, ecm_params->s_bits, sizeof(char) * (ecm_params->num_bits), cudaMemcpyHostToDevice));
  // Copy instances
  CUDA_CHECK(cudaMalloc((void **)&gpu_instances, instance_size));
  CUDA_CHECK(cudaMemcpy(gpu_instances, instances, instance_size, cudaMemcpyHostToDevice));

  // create a cgbn_error_report for CGBN to report back errors
  CUDA_CHECK(cgbn_error_report_alloc(&report));

  outputf (OUTPUT_VERBOSE, "Running GPU kernel<%ld,%d> ...\n", gpu_block_count, TPB);
  kernel_double_add<cgbn_params><<<gpu_block_count, TPB>>>(
    report, ecm_params->num_bits, gpu_s_bits, gpu_instances, ecm_params->curves);

  // error report uses managed memory, so we sync the device (or stream) and check for cgbn errors
  CUDA_CHECK(cudaDeviceSynchronize());
  if (report->_error) {
      outputf (OUTPUT_ERROR, "\n\nerror: %d\n", report->_error);
  }
  CGBN_CHECK(report);

  // Copy the instances back from gpuMemory
  outputf (OUTPUT_VERBOSE, "Copying results back to CPU ...\n");
  CUDA_CHECK(cudaMemcpy(instances, gpu_instances, instance_size, cudaMemcpyDeviceToHost));

  CUDA_CHECK(cudaFree(gpu_s_bits));
  CUDA_CHECK(cudaFree(gpu_instances));
  CUDA_CHECK(cgbn_error_report_free(report));

  /* gputime is measured in ms */
  CUDA_CHECK(cudaEventRecord (stop, 0));
  CUDA_CHECK(cudaEventSynchronize (stop));

  cudaEventElapsedTime (gputime, start, stop);
  CUDA_CHECK(cudaEventDestroy (start));
  CUDA_CHECK(cudaEventDestroy (stop));

  int youpi = ECM_NO_FACTOR_FOUND;
  { // Process Results
    mpz_t x_final, y_final, modulo;
    mpz_init(modulo);
    mpz_init(x_final);
    mpz_init(y_final);

    for(size_t i = 0; i < ecm_params->curves; i++) {
      instance_t &instance = instances[i];

      if (test_verbose (OUTPUT_TRACE) && i == 0) {
          to_mpz(x_final, instance.aX._limbs, cgbn_params::BITS/32);
          to_mpz(y_final, instance.aY._limbs, cgbn_params::BITS/32);
          outputf (OUTPUT_TRACE, "index: 0 pA: (%Zd, %Zd)\n", x_final, y_final);

          to_mpz(x_final, instance.bX._limbs, cgbn_params::BITS/32);
          to_mpz(y_final, instance.bY._limbs, cgbn_params::BITS/32);
          outputf (OUTPUT_TRACE, "index: 0 pB: (%Zd, %Zd)\n", x_final, y_final);
      }

      // Make sure we were testing the right number.
      to_mpz(modulo, instance.modulus._limbs, cgbn_params::BITS/32);
      assert(mpz_cmp(modulo, N) == 0);

      to_mpz(x_final, instance.aX._limbs, cgbn_params::BITS/32);
      to_mpz(y_final, instance.aY._limbs, cgbn_params::BITS/32);

      array_stage_found[i] = findfactor(factors[i], N, x_final, y_final);
      if (array_stage_found[i] != ECM_NO_FACTOR_FOUND) {
          youpi = array_stage_found[i];
          outputf (OUTPUT_NORMAL, "GPU: factor %Zd found in Step 1 with curve %ld (-sigma %d:%d)\n",
                  factors[i], i, ECM_PARAM_BATCH_32BITS_D, instance.d);
      }
    }

    mpz_init(modulo);
    mpz_clear(x_final);
    mpz_clear(y_final);

    // clean up
    free(ecm_params->s_bits);
    free(instances);
  }

  return youpi;
}

#endif  /* _CGBN_STAGE1_CU */
