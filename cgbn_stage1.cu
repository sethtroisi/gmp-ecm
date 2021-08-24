/***

Copyright (c) 2018-2019, NVIDIA CORPORATION.  All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
IN THE SOFTWARE.

***/

#include "cgbn_stage1.h"

#include <cassert>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

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

// support routine copied from  "CGBN/samples/utility/support.h"
void cgbn_check(cgbn_error_report_t *report, const char *file=NULL, int32_t line=0) {
  // check for cgbn errors

  if(cgbn_error_report_check(report)) {
    fprintf (stderr, "\n");
    fprintf (stderr, "CGBN error occurred: %s\n", cgbn_error_string(report));

    if(report->_instance!=0xFFFFFFFF) {
      fprintf (stderr, "Error reported by instance %d", report->_instance);
      if(report->_blockIdx.x!=0xFFFFFFFF || report->_threadIdx.x!=0xFFFFFFFF)
        fprintf (stderr, ", ");
      if(report->_blockIdx.x!=0xFFFFFFFF)
      fprintf (stderr, "blockIdx=(%d, %d, %d) ", report->_blockIdx.x, report->_blockIdx.y, report->_blockIdx.z);
      if(report->_threadIdx.x!=0xFFFFFFFF)
        fprintf (stderr, "threadIdx=(%d, %d, %d)", report->_threadIdx.x, report->_threadIdx.y, report->_threadIdx.z);
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

static
void to_mpz(mpz_t r, const uint32_t *x, uint32_t count) {
  mpz_import(r, count, -1, sizeof(uint32_t), 0, 0, x);
}

static
void from_mpz(const mpz_t s, uint32_t *x, uint32_t count) {
  size_t words;

  if(mpz_sizeinbase(s, 2)>count*32) {
    fprintf (stderr, "from_mpz failed -- result does not fit\n");
    exit(1);
  }

  mpz_export(x, &words, -1, sizeof(uint32_t), 0, 0, s);
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
//   CONSTANT_TIME   - require constant time algorithms (currently, constant time algorithms are not available)

// Locally it will also be helpful to have several parameters:
//   TPI             - threads per instance
//   BITS            - number of bits per instance

const uint32_t TPB_DEFAULT = 128;

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
};


static
uint32_t* set_p_2p(const mpz_t N, const mpz_t s,
                   uint32_t curves, uint32_t sigma,
                   uint32_t BITS, size_t *data_size) {
  // P1_x, P1_y = (2,1)
  // 2P_x, 2P_y = (9, 64 * d + 8)

  /** Consider curve copies of N (AKA modulo) */
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
        char* gpu_s_bits,
        uint32_t *data,
        uint32_t count,
        uint32_t sigma_0) {
  // decode an instance_i number from the blockIdx and threadIdx
  int32_t instance_i = (blockIdx.x*blockDim.x + threadIdx.x)/params::TPI;
  if(instance_i >= count)
    return;

  /* Cast uint32_t array to mem_t */
  typename curve_t<params>::mem_t *data_cast = (typename curve_t<params>::mem_t*) data;

  cgbn_monitor_t monitor = CHECK_ERROR ? cgbn_report_monitor : cgbn_no_checks;

  curve_t<params> curve(monitor, report, instance_i);
  typename curve_t<params>::bn_t  aX, aY, bX, bY, modulus;
 
  uint32_t np0;
  { // Setup
      cgbn_load(curve._env, modulus, &data_cast[5*instance_i+0]);
      cgbn_load(curve._env, aX, &data_cast[5*instance_i+1]);
      cgbn_load(curve._env, aY, &data_cast[5*instance_i+2]);
      cgbn_load(curve._env, bX, &data_cast[5*instance_i+3]);
      cgbn_load(curve._env, bY, &data_cast[5*instance_i+4]);

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

  uint32_t d = sigma_0 + instance_i;

  for (int b = s_bits; b > 0; b--) {
    /**
     * TODO generates a lot of duplicate inlined code, not sure how to improve
     * Tried with swappings pointers (with a single call to double_add_v2)
     */
    if (gpu_s_bits[s_bits - b] == 0) {
        curve.double_add_v2(aX, aY, bX, bY, d, b, modulus, np0);
    } else {
        curve.double_add_v2(bX, bY, aX, aY, d, b, modulus, np0);
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


static
int verify_size_of_n(const mpz_t N, size_t max_bits) {
  size_t n_log2 = mpz_sizeinbase(N, 2);

  // using check_gpuecm.sage it looks like 4 bits would suffice
  size_t max_usable_bits = max_bits - 6;

  if (n_log2 <= max_usable_bits)
    return ECM_NO_FACTOR_FOUND;

  outputf (OUTPUT_ERROR, "GPU: N(%d bits) + carry(%d bits) > BITS(%d)\n",
      n_log2, max_bits - max_usable_bits, max_bits);
  outputf (OUTPUT_ERROR, "GPU: Error, input number should be stricly lower than 2^%d\n",
      max_usable_bits);
  return ECM_ERROR;
}


static
char* allocate_and_set_s_bits(const mpz_t s, int *sbits) {
  uint32_t num_bits = *sbits = mpz_sizeinbase(s, 2) - 1;
  assert( 1 <= num_bits <= 100'000'000 );

  // Use int* so that size can be stored in first element, could pass around extra size.
  char *s_bits = (char*) malloc(sizeof(char) * num_bits);
  for (int i = 0; i < num_bits; i++) {
      s_bits[i] = mpz_tstbit (s, num_bits - 1 - i);
  }
  return s_bits;
}

static
int process_results(mpz_t *factors, int *array_stage_found,
                    const mpz_t N,
                    const uint32_t *data, uint32_t cgbn_bits,
                    int curves, uint32_t sigma) {
  mpz_t x_final, y_final, modulo;
  mpz_init(modulo);
  mpz_init(x_final);
  mpz_init(y_final);

  const uint32_t limbs_per = cgbn_bits / 32;

  int youpi = ECM_NO_FACTOR_FOUND;
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

    array_stage_found[i] = findfactor(factors[i], N, x_final, y_final);
    if (array_stage_found[i] != ECM_NO_FACTOR_FOUND) {
      youpi = array_stage_found[i];
      outputf (OUTPUT_NORMAL, "GPU: factor %Zd found in Step 1 with curve %ld (-sigma %d:%d)\n",
          factors[i], i, ECM_PARAM_BATCH_32BITS_D, sigma + i);
    }
  }
  
  mpz_init(modulo);
  mpz_clear(x_final);
  mpz_clear(y_final);

  return youpi;
}

int run_cgbn(mpz_t *factors, int *array_stage_found,
             const mpz_t N, const mpz_t s, float *gputime,
             ecm_params_t *ecm_params) {

  const size_t MAX_BITS = 1024;

  size_t curves = ecm_params->curves;  
  assert( ecm_params->sigma > 0 );
  assert( ((uint64_t) ecm_params->sigma + curves) <= 0xFFFFFFFF ); // no overflow

  /* Validate N's size */
  int youpi = verify_size_of_n(N, MAX_BITS);
  if (youpi != ECM_NO_FACTOR_FOUND) {
    return youpi;
  }

  int s_num_bits;
  char *s_bits = allocate_and_set_s_bits(s, &s_num_bits);
  assert( 1 <= s_num_bits <= 100'000'000 );
  assert( s_bits != NULL );

  // Keeps CPU from busy waiting during GPU execution.
  CUDA_CHECK(cudaSetDeviceFlags (cudaDeviceScheduleBlockingSync));
  cudaEvent_t start, stop;
  CUDA_CHECK(cudaEventCreate (&start));
  CUDA_CHECK(cudaEventCreate (&stop));
  CUDA_CHECK(cudaEventRecord (start, 0));

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
   * Smaller TPI (e.g. 8) is faster (TPI=4 seems worse than TPI=8).
   * Larger TPI (e.g. 32) for running a single curve (or large N).
   * TPI=8  is required for N > 512
   * TPI=16 is required for N > 2048
   * TPI=32 is required for N > 8192
   */
  /**
   * Larger takes longer to compile:
   * (32,32768) takes ~10 minutes (ecm was 6.0M)
   * (32,16384) takes ~2  minutes (ecm was 4.3M)
   * (16,8192) takes ~20 seconds (ecm was 4.3M)
   * (8, 1024) takes ~10 seconds (ecm was 3.7M)
   * GPU, No CGBN                (ecm was 3.5M)
   * No GPU, No CGBN             (ecm was 3.4M)
   */
  /**
   * TPI and BITS have to be set at compile time.
   * Adding multiple cgbn_params (and kernals) allows for dynamic selection
   * based on the size of N (e.g. N < 1024, N < 2048, N < 4096). If we go that
   * route it would be helpful to do that only during release builds.
   */

  typedef cgbn_params_t<8, 1024>  cgbn_params_8_1024;
  typedef cgbn_params_t<8, 512>   cgbn_params_8_512;
#ifdef IS_DEV_BUILD
  const std::vector<uint32_t> available_kernels = { 1024 };
#else
  typedef cgbn_params_t<8, 1024>  cgbn_params_8_1024;
  typedef cgbn_params_t<8, 1536>  cgbn_params_8_1536;
  typedef cgbn_params_t<8, 2048>  cgbn_params_8_2048;
  typedef cgbn_params_t<16, 3072> cgbn_params_16_3072;
  const std::vector<uint32_t> available_kernels = { 512, 1024, 1536, 2048, 3072 };
#endif /* IS_DEV_BUILD */

  for (uint32_t kernel_bits : available_kernels) {
    if (kernel_bits + 6 >=  mpz_sizeinbase(N, 2)) {
      BITS = kernel_bits;
      assert( BITS % 32 == 0 );
      TPI = (BITS <= 2048) ? 8 : ((BITS <= 8192) ? 16 : 32);
      IPB = TPB / TPI;
      BLOCK_COUNT = (curves + IPB - 1) / IPB;
      break;
    }
  }
  if (BITS == 0) {
    outputf (OUTPUT_ERROR, "No available CGBN Kernel large enough to process N(%d bits)\n",
        mpz_sizeinbase(N, 2));
    return ECM_ERROR;
  }

  /** Relies on mem_t (AKA struct cgbn_mem_t) being byte aligned without extra fields. */
  assert( sizeof(curve_t<cgbn_params_8_512>::mem_t) == 512/8 );
  assert( sizeof(curve_t<cgbn_params_8_1024>::mem_t) == 1024/8 );
  data = set_p_2p(N, s, curves, ecm_params->sigma, BITS, &data_size);

  // Copy data
  outputf (OUTPUT_VERBOSE, "Copying %d bits of data to GPU\n", data_size);
  CUDA_CHECK(cudaMalloc((void **)&gpu_data, data_size));
  CUDA_CHECK(cudaMemcpy(gpu_data, data, data_size, cudaMemcpyHostToDevice));

  outputf (OUTPUT_VERBOSE, "Running CGBN<%d,%d> kernel<%ld,%d>...\n", BITS, TPI, BLOCK_COUNT, TPB);
  if (BITS == 512) {
    kernel_double_add<cgbn_params_8_512><<<BLOCK_COUNT, TPB>>>(report, s_num_bits, gpu_s_bits, gpu_data, curves, ecm_params->sigma);
  } else if (BITS == 1024) {
    kernel_double_add<cgbn_params_8_1024><<<BLOCK_COUNT, TPB>>>(report, s_num_bits, gpu_s_bits, gpu_data, curves, ecm_params->sigma);
#ifndef IS_DEV_BUILD
  } else if (BITS == 1536) {
    kernel_double_add<cgbn_params_8_1536><<<BLOCK_COUNT, TPB>>>(report, s_num_bits, gpu_s_bits, gpu_data, curves, ecm_params->sigma);
  } else if (BITS == 2048) {
    kernel_double_add<cgbn_params_8_2048><<<BLOCK_COUNT, TPB>>>(report, s_num_bits, gpu_s_bits, gpu_data, curves, ecm_params->sigma);
  } else if (BITS == 3072) {
    kernel_double_add<cgbn_params_16_3072><<<BLOCK_COUNT, TPB>>>(report, s_num_bits, gpu_s_bits, gpu_data, curves, ecm_params->sigma);
#endif
  } else {
    outputf (OUTPUT_ERROR, "CGBN Kernel not found for %d bits\n", BITS);
    return ECM_ERROR;
  }

  /* error report uses managed memory, sync the device and check for cgbn errors */
  CUDA_CHECK(cudaDeviceSynchronize());
  if (report->_error)
    outputf (OUTPUT_ERROR, "\n\nerror: %d\n", report->_error);
  CGBN_CHECK(report);

  /* gputime is measured in ms */
  CUDA_CHECK(cudaEventRecord (stop, 0));
  CUDA_CHECK(cudaEventSynchronize (stop));

  // Copy data back from GPU memory
  outputf (OUTPUT_VERBOSE, "Copying results back to CPU ...\n");
  CUDA_CHECK(cudaMemcpy(data, gpu_data, data_size, cudaMemcpyDeviceToHost));

  cudaEventElapsedTime (gputime, start, stop);

  youpi = process_results(
      factors, array_stage_found, N,
      data, BITS,
      curves, ecm_params->sigma);

  // clean up
  CUDA_CHECK(cudaFree(gpu_s_bits));
  CUDA_CHECK(cudaFree(gpu_data));
  CUDA_CHECK(cgbn_error_report_free(report));
  CUDA_CHECK(cudaEventDestroy (start));
  CUDA_CHECK(cudaEventDestroy (stop));

  free(s_bits);
  free(data);

  return youpi;
}
