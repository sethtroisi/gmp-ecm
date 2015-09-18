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

/* top-level interface for a group of transforms that share 
   the same compiler optimizations and prime size */

typedef void * (*mpzspm_init_t)(uint32_t max_len, mpz_t modulus,
				mpz_t P, mpz_t S, uint32_t *done);
                                
typedef void (*mpzspm_clear_t)(void * mpzspm);
                                
typedef struct
  {
    uint32_t sp_bits;
    mpzspm_init_t mpzspm_init;
    mpzspm_clear_t mpzspm_clear;
  } __nttinit_struct;

typedef __nttinit_struct * nttinit_t;

/* all such groups we know of */

#define DECLARE(spbits, wbits) \
extern const __nttinit_struct nttinit_sp##spbits##w##wbits; \
int test_main_sp##spbits##w##wbits(int argc, char **argv);

#if GMP_LIMB_BITS == 32
DECLARE(30,32)
DECLARE(31,32)
DECLARE(50,32)
DECLARE(62,32)
#else
DECLARE(30,64)
DECLARE(50,64)
DECLARE(62,64)
#endif

  /* all the data for convolutions of specified size with a
     single modulus */

typedef struct
{
  mpz_t modulus;
  uint32_t max_len;
  uint32_t mpzspm_num;
  nttinit_t *nttinit;
  void **mpzspm;
} __nttwork_struct;

typedef __nttwork_struct * nttwork_t;

nttwork_t nttwork_init(uint32_t max_len, mpz_t modulus,
                        uint32_t *sp_sizes, uint32_t num_sp_sizes);

void nttwork_clear(nttwork_t nttwork);

#ifdef __cplusplus
}
#endif

#endif /* _LIBNTT_H */
