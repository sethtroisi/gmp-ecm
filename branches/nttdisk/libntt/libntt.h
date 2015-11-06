#ifndef _LIBNTT_H
#define _LIBNTT_H

#include "basicdefs.h"
#include "gmp-xface.h"

#ifdef __cplusplus
extern "C" {
#endif

/* helper data for primality tests */

extern const uint8_t sprp32_lookup[]; 
extern const uint16_t sprp64_lookup[]; 

/* wrapper for aligned memory */

void * sp_aligned_malloc (size_t len);
void sp_aligned_free (void *newptr);

static INLINE uint64_t
read_clock(void) 
{
#if defined(_MSC_VER)
	LARGE_INTEGER ret;
	QueryPerformanceCounter(&ret);
	return ret.QuadPart;
#else
	uint32_t lo, hi;
	asm volatile("rdtsc":"=d"(hi),"=a"(lo));
	return (uint64_t)hi << 32 | lo;
#endif
}

/* an NTT is built up of one or more passes through
   the input data */

typedef enum
{
  PASS_TYPE_DIRECT,
  PASS_TYPE_PFA,
  PASS_TYPE_TWIDDLE
} pass_type_t;

/* guides for constructing transforms */

typedef struct
{
  uint32_t codelet_size;
  uint32_t group_type;
  pass_type_t pass_type;
} nttplan_t;

#define MAX_PLANS 10

typedef struct
{
  uint32_t num_plans;
  nttplan_t plans[MAX_PLANS];
} nttplangroup_t;

/* top-level interface for a group of transforms that share 
   the same compiler optimizations and prime size */

typedef uint32_t (*mpzspm_get_num_groups_t)(void);

typedef void * (*mpzspm_init_t)(uint32_t max_ntt_size, mpz_t modulus,
				mpz_t P, mpz_t S, uint32_t interleaved, uint32_t *done);
                                
typedef void (*mpzspm_clear_t)(void * mpzspm);

typedef void (*mpzspm_ntt_init_t)(void * mpzspm, uint32_t ntt_size,
    				uint32_t max_ntt_size, nttplangroup_t *p);

typedef void (*mpzspm_ntt_run_t)(void * mpzspm, mpz_t * x, uint32_t ntt_size);

typedef void (*mpzspm_ntt_reset_t)(void * mpzspm);

typedef void (*mpzspm_random_t)(void * mpzspm, uint32_t ntt_size);

typedef void (*mpzspm_test_t)(void * mpzspm, uint32_t ntt_size,
    				uint32_t max_ntt_size);

typedef struct
  {
    const uint32_t sp_bits;
    mpzspm_get_num_groups_t mpzspm_get_num_groups;
    mpzspm_init_t mpzspm_init;
    mpzspm_clear_t mpzspm_clear;
    mpzspm_ntt_init_t mpzspm_ntt_init;
    mpzspm_ntt_run_t mpzspm_ntt_run;
    mpzspm_ntt_reset_t mpzspm_ntt_reset;
    mpzspm_random_t mpzspm_random;
    mpzspm_test_t mpzspm_test;
  } __nttinit_struct;

typedef __nttinit_struct * nttinit_t;

/* all such groups we know of */

#define DECLARE(spbits, wbits) \
extern const __nttinit_struct nttinit_sp##spbits##w##wbits;

#if GMP_LIMB_BITS == 32
DECLARE(30,32)
DECLARE(31,32)
DECLARE(50,32)
DECLARE(62,32)
#else
DECLARE(30,64)
DECLARE(31,64)
DECLARE(50,64)
DECLARE(62,64)
#endif

  /* all the data for convolutions of specified size with a
     single modulus */

typedef struct
{
  mpz_t modulus;
  uint32_t max_ntt_size;
  uint32_t ntt_size;
  uint32_t mpzspm_num;
  nttinit_t *nttinit;
  void **mpzspm;
} __nttwork_struct;

typedef __nttwork_struct * nttwork_t;

nttwork_t nttwork_init(uint32_t max_ntt_size, mpz_t modulus, 
                        uint32_t interleaved,
                        uint32_t *sp_sizes, uint32_t num_sp_sizes);

void nttwork_clear(nttwork_t nttwork);

void nttwork_ntt_init(nttwork_t nttwork, nttplangroup_t *plans);

void nttwork_ntt_run(nttwork_t nttwork, mpz_t * x);

void nttwork_ntt_reset(nttwork_t nttwork);

void nttwork_random(nttwork_t nttwork);

double nttwork_ntt_test(nttwork_t nttwork, uint32_t verify);

#ifdef __cplusplus
}
#endif

#endif /* _LIBNTT_H */
