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

#ifndef __CUDACC__
#error "This file should only be compiled with nvcc"
#endif

#include "cgbn_stage1.h"

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

static
void to_mpz(mpz_t r, const uint32_t *x, uint32_t count) {
  mpz_import (r, count, -1, sizeof(uint32_t), 0, 0, x);
}

static
void from_mpz(const mpz_t s, uint32_t *x, uint32_t count) {
  size_t words;

  if(mpz_sizeinbase (s, 2) > count * 32) {
    fprintf (stderr, "from_mpz failed -- result does not fit\n");
    exit(EXIT_FAILURE);
  }

  mpz_export (x, &words, -1, sizeof(uint32_t), 0, 0, s);
  while(words<count)
    x[words++]=0;
}



// ---------------------------------------------------------------- //

// The CGBN context uses the following three parameters:
//   TBP             - threads per block (zero means to use the blockDim.x)
//   MAX_ROTATION    - must be small power of 2, imperically, 4 works well
//   CONSTANT_TIME   - require constant time algorithms (currently, constant time algorithms are not available)

// Locally it will also be helpful to have several parameters:
//   TPI             - threads per instance
//   BITS            - number of bits per instance

/* TODO test how this changes gpu_throughput_test */
/* NOTE: >= 512 may not be supported for > 2048 bit kernels */
const uint32_t TPB_DEFAULT = 256;

template<uint32_t tpi, uint32_t bits>
class cgbn_params_t {
  public:
  // parameters used by the CGBN context
  static const uint32_t TPB=TPB_DEFAULT;           // Reasonable default
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

  typedef cgbn_context_t<params::TPI, params>   context_t;
  typedef cgbn_env_t<context_t, params::BITS>   env_t;
  typedef typename env_t::cgbn_t                bn_t;
  typedef cgbn_mem_t<params::BITS>              mem_t;

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

    /* Roughly 40% of time is spent in these two calls */
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
};

static
uint32_t* set_p_2p(const mpz_t N,
                   uint32_t curves, uint32_t sigma,
                   uint32_t BITS, size_t *data_size) {
  // P1_x, P1_y = (2,1)
  // 2P_x, 2P_y = (9, 64 * d + 8)

  /** Keeps a copy of N (AKA modulo) per curve */
  const size_t limbs_per = BITS/32;
  *data_size = 5 * curves * limbs_per * sizeof(uint32_t);
  uint32_t *data = (uint32_t*) malloc(*data_size);
  uint32_t *datum = data;

  mpz_t x;
  mpz_init(x);
  for(int index = 0; index < curves; index++) {
      // d = (sigma / 2^32) mod N BUT 2^32 handled by special_mul_ui32
      uint32_t d = sigma + index;

      // mod
      from_mpz(N, datum + 0 * limbs_per, BITS/32);

      // P1 (X, Y)
      mpz_set_ui(x, 2);
      from_mpz(x, datum + 1 * limbs_per, BITS/32);
      mpz_set_ui(x, 1);
      from_mpz(x, datum + 2 * limbs_per, BITS/32);

      // 2P = P2 (X, Y)
      // P2_y = 64 * d + 8
      mpz_set_ui(x, 9);
      from_mpz(x, datum + 3 * limbs_per, BITS/32);

      // d = sigma * mod_inverse(2 ** 32, N)
      mpz_ui_pow_ui(x, 2, 32);
      mpz_invert(x, x, N);
      mpz_mul_ui(x, x, d);
      // P2_y = 64 * d - 2;
      mpz_mul_ui(x, x, 64);
      mpz_add_ui(x, x, 8);
      mpz_mod(x, x, N);

      outputf (OUTPUT_TRACE, "sigma %d => P2_y: %Zd\n", d, x);
      from_mpz(x, datum + 4 * limbs_per, BITS/32);
      datum += 5 * limbs_per;
  }
  mpz_clear(x);
  return data;
}


// kernel implementation using cgbn
template<class params>
__global__ void kernel_double_add(
        cgbn_error_report_t *report,
        uint32_t s_bits,
        uint32_t s_bits_start,
        uint32_t s_bits_interval,
        char* gpu_s_bits,
        uint32_t *data,
        uint32_t count,
        uint32_t sigma_0,
        uint32_t np0
        ) {
  // decode an instance_i number from the blockIdx and threadIdx
  int32_t instance_i = (blockIdx.x*blockDim.x + threadIdx.x)/params::TPI;
  if(instance_i >= count)
    return;

  /* Cast uint32_t array to mem_t */
  typename curve_t<params>::mem_t *data_cast = (typename curve_t<params>::mem_t*) data;

  cgbn_monitor_t monitor = CHECK_ERROR ? cgbn_report_monitor : cgbn_no_checks;

  curve_t<params> curve(monitor, report, instance_i);
  typename curve_t<params>::bn_t  aX, aY, bX, bY, modulus;

  { // Setup
      cgbn_load(curve._env, modulus, &data_cast[5*instance_i+0]);
      cgbn_load(curve._env, aX, &data_cast[5*instance_i+1]);
      cgbn_load(curve._env, aY, &data_cast[5*instance_i+2]);
      cgbn_load(curve._env, bX, &data_cast[5*instance_i+3]);
      cgbn_load(curve._env, bY, &data_cast[5*instance_i+4]);

      /* Convert points to mont, has a miniscule bit of overhead with batching. */
      uint32_t np0_test = cgbn_bn2mont(curve._env, aX, aX, modulus);
      assert(np0 == np0_test);

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

  uint32_t d = sigma_0 + instance_i;
  int swapped = 0;
  for (int b = 0; b < s_bits_interval; b++) {
    if (gpu_s_bits[s_bits_start + b] != swapped) {
        swapped = !swapped;
        cgbn_swap(curve._env, aX, bX);
        cgbn_swap(curve._env, aY, bY);
    }
    curve.double_add_v2(aX, aY, bX, bY, d, modulus, np0);
  }

  if (swapped) {
    cgbn_swap(curve._env, aX, bX);
    cgbn_swap(curve._env, aY, bY);
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
    cgbn_store(curve._env, &data_cast[5*instance_i+1], aX);
    cgbn_store(curve._env, &data_cast[5*instance_i+2], aY);
    cgbn_store(curve._env, &data_cast[5*instance_i+3], bX);
    cgbn_store(curve._env, &data_cast[5*instance_i+4], bY);
  }
}

static
int findfactor(mpz_t factor, const mpz_t N, const mpz_t x_final, const mpz_t y_final) {
    // XXX: combine / refactor logic with cudawrapper.c findfactor

    mpz_t temp;
    mpz_init(temp);

    /* Check if factor found */
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


static
int verify_size_of_n(const mpz_t N, size_t max_bits) {
  size_t n_log2 = mpz_sizeinbase(N, 2);

  /* Using check_gpuecm.sage it looks like 4 bits would suffice. */
  size_t max_usable_bits = max_bits - CARRY_BITS;

  if (n_log2 <= max_usable_bits)
    return ECM_NO_FACTOR_FOUND;

  outputf (OUTPUT_ERROR, "GPU: N(%d bits) + carry(%d bits) > BITS(%d)\n",
      n_log2, CARRY_BITS, max_bits);
  outputf (OUTPUT_ERROR, "GPU: Error, input number should be stricly lower than 2^%d\n",
      max_usable_bits);
  return ECM_ERROR;
}


static
uint32_t find_np0(const mpz_t N) {
  uint32_t np0;
  mpz_t temp;
  mpz_init(temp);
  mpz_ui_pow_ui(temp, 2, 32);
  assert(mpz_invert(temp, N, temp));
  np0 = -mpz_get_ui(temp);
  mpz_clear(temp);
  return np0;
}


static
char* allocate_and_set_s_bits(const mpz_t s, int *nbits) {
  uint32_t num_bits = *nbits = mpz_sizeinbase(s, 2) - 1;
  /* XXX: Remove arbitrary limit on S (based on B1) */
  assert( 1 <= num_bits && num_bits <= 100000000 );

  char *s_bits = (char*) malloc(sizeof(char) * num_bits);
  for (int i = 0; i < num_bits; i++) {
      s_bits[i] = mpz_tstbit (s, num_bits - 1 - i);
  }
  return s_bits;
}

static
int process_results(mpz_t *factors, int *array_found,
                    const mpz_t N,
                    const uint32_t *data, uint32_t cgbn_bits,
                    int curves, uint32_t sigma) {
  mpz_t x_final, y_final, modulo;
  mpz_init(modulo);
  mpz_init(x_final);
  mpz_init(y_final);

  const uint32_t limbs_per = cgbn_bits / 32;

  int youpi = ECM_NO_FACTOR_FOUND;
  int errors = 0;
  for(size_t i = 0; i < curves; i++) {
    const uint32_t *datum = data + (5 * i * limbs_per);;

    if (test_verbose (OUTPUT_TRACE) && i == 0) {
      to_mpz(modulo, datum + 0 * limbs_per, limbs_per);
      outputf (OUTPUT_TRACE, "index: 0 modulo: %Zd\n", modulo);

      to_mpz(x_final, datum + 1 * limbs_per, limbs_per);
      to_mpz(y_final, datum + 2 * limbs_per, limbs_per);
      outputf (OUTPUT_TRACE, "index: 0 pA: (%Zd, %Zd)\n", x_final, y_final);

      to_mpz(x_final, datum + 3 * limbs_per, limbs_per);
      to_mpz(y_final, datum + 4 * limbs_per, limbs_per);
      outputf (OUTPUT_TRACE, "index: 0 pB: (%Zd, %Zd)\n", x_final, y_final);
    }

    // Make sure we were testing the right number.
    to_mpz(modulo, datum + 0 * limbs_per, limbs_per);
    assert(mpz_cmp(modulo, N) == 0);

    to_mpz(x_final, datum + 1 * limbs_per, limbs_per);
    to_mpz(y_final, datum + 2 * limbs_per, limbs_per);

    /* Very suspicious for (x_final, y_final) to match (x_0, y_0) == (2, 1)
     * Can happen when
     * 1. block calculation performed incorrectly (and some blocks not run)
     * 2. Kernel didn't run because not enough register
     * 3. nvcc links old version of kernel when something changed
     */
    if (mpz_cmp_ui (x_final, 2) == 0 && mpz_cmp_ui (y_final, 1) == 0) {
      errors += 1;
      if (errors < 10 || errors % 100 == 1)
        outputf (OUTPUT_ERROR, "GPU: curve %d didn't compute?\n", i);
    }

    array_found[i] = findfactor(factors[i], N, x_final, y_final);
    if (array_found[i] != ECM_NO_FACTOR_FOUND) {
      youpi = array_found[i];
      outputf (OUTPUT_NORMAL, "GPU: factor %Zd found in Step 1 with curve %ld (-sigma %d:%d)\n",
          factors[i], i, ECM_PARAM_BATCH_32BITS_D, sigma + i);
    }
  }

  mpz_init(modulo);
  mpz_clear(x_final);
  mpz_clear(y_final);

#ifdef IS_DEV_BUILD
  if (errors)
        outputf (OUTPUT_ERROR, "Had %d errors. Try `make clean; make` or reducing TPB_DEFAULT\n",
            errors);
#endif

  if (errors > 2)
      return ECM_ERROR;

  return youpi;
}

int cgbn_ecm_stage1(mpz_t *factors, int *array_found,
             const mpz_t N, const mpz_t s,
             uint32_t curves, uint32_t sigma,
             float *gputime, int verbose)
{
  assert( sigma > 0 );
  assert( ((uint64_t) sigma + curves) <= 0xFFFFFFFF ); // no overflow

  int s_num_bits;
  char *s_bits = allocate_and_set_s_bits(s, &s_num_bits);
  assert( 1 <= s_num_bits && s_num_bits <= 100000000 );
  assert( s_bits != NULL );

  cudaEvent_t global_start, batch_start, stop;
  CUDA_CHECK(cudaEventCreate (&global_start));
  CUDA_CHECK(cudaEventCreate (&batch_start));
  CUDA_CHECK(cudaEventCreate (&stop));
  CUDA_CHECK(cudaEventRecord (global_start));

  // Copy s_bits
  char     *gpu_s_bits;
  CUDA_CHECK(cudaMalloc((void **)&gpu_s_bits, sizeof(char) * s_num_bits));
  CUDA_CHECK(cudaMemcpy(gpu_s_bits, s_bits, sizeof(char) * s_num_bits, cudaMemcpyHostToDevice));

  cgbn_error_report_t *report;
  // create a cgbn_error_report for CGBN to report back errors
  CUDA_CHECK(cgbn_error_report_alloc(&report));

  size_t    data_size;
  uint32_t *data, *gpu_data;

  uint32_t  BITS = 0;        // kernel bits
  int32_t   TPB=TPB_DEFAULT; // Always the same default
  int32_t   TPI;
  int32_t   IPB;             // IPB = TPB / TPI, instances per block
  size_t    BLOCK_COUNT;     // How many blocks to cover all curves

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
   * For "Compling custom kernel for %d bits should be XX% faster"
   * Change the 512 in cgbn_params_t<4, 512> cgbn_params_512;
   * to the suggested value (a multiple of 32 >= bits + 6).
   * You may need to change the 4 to an 8 (or 16) if bits >512, >2048
   */
  /** TODO: try with const vector for BITs/TPI, see if compiler is happy */
  std::vector<uint32_t> available_kernels;

  typedef cgbn_params_t<4, 512>   cgbn_params_512;
  typedef cgbn_params_t<8, 1024>  cgbn_params_1024;
  available_kernels.push_back((uint32_t)cgbn_params_512::BITS);
  available_kernels.push_back((uint32_t)cgbn_params_1024::BITS);

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

  size_t n_log2 = mpz_sizeinbase(N, 2);
  for (int k_i = 0; k_i < available_kernels.size(); k_i++) {
    uint32_t kernel_bits = available_kernels[k_i];
    if (kernel_bits >= n_log2 + CARRY_BITS) {
      BITS = kernel_bits;
      assert( BITS % 32 == 0 );

      /* Print some debug info about kernel. */
      /* TODO: return kernelAttr and validate maxThreadsPerBlock. */
      if (BITS == cgbn_params_512::BITS) {
        TPI = cgbn_params_512::TPI;
        kernel_info((const void*)kernel_double_add<cgbn_params_512>, verbose);
      } else if (BITS == cgbn_params_1024::BITS) {
        TPI = cgbn_params_1024::TPI;
        kernel_info((const void*)kernel_double_add<cgbn_params_1024>, verbose);
#ifndef IS_DEV_BUILD
      } else if (BITS == cgbn_params_1536::BITS) {
        TPI = cgbn_params_1536::TPI;
        kernel_info((const void*)kernel_double_add<cgbn_params_1536>, verbose);
      } else if (BITS == cgbn_params_2048::BITS) {
        TPI = cgbn_params_2048::TPI;
        kernel_info((const void*)kernel_double_add<cgbn_params_2048>, verbose);
      } else if (BITS == cgbn_params_3072::BITS) {
        TPI = cgbn_params_3072::TPI;
        kernel_info((const void*)kernel_double_add<cgbn_params_3072>, verbose);
      } else if (BITS == cgbn_params_4096::BITS) {
        TPI = cgbn_params_4096::TPI;
        kernel_info((const void*)kernel_double_add<cgbn_params_4096>, verbose);
#endif
      } else {
        /* lowercase k to help differentiate this error from one below */
        outputf (OUTPUT_ERROR, "CGBN kernel not found for %d bits\n", BITS);
        return ECM_ERROR;
      }

      IPB = TPB / TPI;
      BLOCK_COUNT = (curves + IPB - 1) / IPB;

      break;
    }
  }

  if (BITS == 0) {
    outputf (OUTPUT_ERROR, "No available CGBN Kernel large enough to process N(%d bits)\n", n_log2);
    return ECM_ERROR;
  }

  /* Alert that recompiling with a smaller kernel would likely improve speed */
  {
    size_t optimized_bits = ((n_log2 + 5)/128 + 1) * 128;
    if (optimized_bits < BITS && 0.8 * BITS > n_log2 ) {
      /* Assume speed is roughly O(N) but slightly slower for not being a power of two */
      float pct_faster = 90 * BITS / optimized_bits;
      assert(pct_faster > 100);
      outputf (OUTPUT_VERBOSE, "Compiling custom kernel for %d bits should be ~%.0f%% faster\n",
              optimized_bits, pct_faster);
    }
  }

  int youpi = verify_size_of_n(N, BITS);
  if (youpi != ECM_NO_FACTOR_FOUND) {
    return youpi;
  }

  /* Consistency check that struct cgbn_mem_t is byte aligned without extra fields. */
  assert( sizeof(curve_t<cgbn_params_512>::mem_t) == cgbn_params_512::BITS/8 );
  assert( sizeof(curve_t<cgbn_params_1024>::mem_t) == cgbn_params_1024::BITS/8 );
  data = set_p_2p(N, curves, sigma, BITS, &data_size);

  /* np0 is -(N^-1 mod 2**32), used for montgomery representation */
  uint32_t np0 = find_np0(N);

  // Copy data
  outputf (OUTPUT_VERBOSE, "Copying %d bits of data to GPU\n", data_size);
  CUDA_CHECK(cudaMalloc((void **)&gpu_data, data_size));
  CUDA_CHECK(cudaMemcpy(gpu_data, data, data_size, cudaMemcpyHostToDevice));

  outputf (OUTPUT_VERBOSE,
          "CGBN<%d, %d> running kernel<%d block x %d threads> input number is %d bits\n",
          BITS, TPI, BLOCK_COUNT, TPB, n_log2);

  /* Break s_num_bits into smaller chunks */
  int s_partial = 0;
  int batches_complete = 0;
  /* Start with small batches and increase till timing is ~100ms */
  int batch_size = 100;
  /* gputime/batch_time are measured in ms */
  float batch_time = 0;

  while (s_partial < s_num_bits) {
    /* decrease batch_size for final batch if needed */
    batch_size = std::min(s_num_bits - s_partial, batch_size);

    /* print ETA with lessing frequently, 5 early + 5 per 10s + 5 per 100s + every 1000s */
    if ((batches_complete < 3) ||
        (batches_complete < 30 && batches_complete % 10 == 0) ||
        (batches_complete < 500 && batches_complete % 100 == 0) ||
        (batches_complete < 5000 && batches_complete % 1000 == 0) ||
        (batches_complete % 10000 == 0)) {
      outputf (OUTPUT_VERBOSE, "Computing %d bits/call, %d/%d (%.1f%%)",
          batch_size, s_partial, s_num_bits, 100.0 * s_partial / s_num_bits);
      if (batches_complete < 2 || *gputime < 1000) {
        outputf (OUTPUT_VERBOSE, "\n");
      } else {
        float estimated_total = (*gputime) * ((float) s_num_bits) / s_partial;
        float eta = estimated_total - (*gputime);
        outputf (OUTPUT_VERBOSE, ", ETA %.f + %.f = %.f seconds (~%.f ms/curves)\n",
                eta / 1000, *gputime / 1000, estimated_total / 1000,
                estimated_total / curves);
      }
    }

    CUDA_CHECK(cudaEventRecord (batch_start));

    if (BITS == cgbn_params_512::BITS) {
      kernel_double_add<cgbn_params_512><<<BLOCK_COUNT, TPB>>>(
          report, s_num_bits, s_partial, batch_size, gpu_s_bits, gpu_data, curves, sigma, np0);
    } else if (BITS == cgbn_params_1024::BITS) {
      kernel_double_add<cgbn_params_1024><<<BLOCK_COUNT, TPB>>>(
          report, s_num_bits, s_partial, batch_size, gpu_s_bits, gpu_data, curves, sigma, np0);
#ifndef IS_DEV_BUILD
    } else if (BITS == cgbn_params_1536::BITS) {
      kernel_double_add<cgbn_params_1536><<<BLOCK_COUNT, TPB>>>(
          report, s_num_bits, s_partial, batch_size, gpu_s_bits, gpu_data, curves, sigma, np0);
    } else if (BITS == cgbn_params_2048::BITS) {
      kernel_double_add<cgbn_params_2048><<<BLOCK_COUNT, TPB>>>(
          report, s_num_bits, s_partial, batch_size, gpu_s_bits, gpu_data, curves, sigma, np0);
    } else if (BITS == cgbn_params_3072::BITS) {
      kernel_double_add<cgbn_params_3072><<<BLOCK_COUNT, TPB>>>(
          report, s_num_bits, s_partial, batch_size, gpu_s_bits, gpu_data, curves, sigma, np0);
    } else if (BITS == cgbn_params_4096::BITS) {
      kernel_double_add<cgbn_params_4096><<<BLOCK_COUNT, TPB>>>(
          report, s_num_bits, s_partial, batch_size, gpu_s_bits, gpu_data, curves, sigma, np0);
#endif
    } else {
      outputf (OUTPUT_ERROR, "CGBN Kernel not found for %d bits\n", BITS);
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
    /* Adjust batch_size to aim for 100ms */
    if (batch_time < 80) {
      batch_size = 11*batch_size/10;
    } else if (batch_time > 120) {
      batch_size = max(100, 9*batch_size / 10);
    }
  }

  // Copy data back from GPU memory
  outputf (OUTPUT_VERBOSE, "Copying results back to CPU ...\n");
  CUDA_CHECK(cudaMemcpy(data, gpu_data, data_size, cudaMemcpyDeviceToHost));

  cudaEventElapsedTime (gputime, global_start, stop);

  youpi = process_results(factors, array_found, N, data, BITS, curves, sigma);

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
  #if __CUDA_ARCH__ < 300
    #error "Unsupported architecture"
  #endif
#endif

#endif  /* _CGBN_STAGE1_CU */
