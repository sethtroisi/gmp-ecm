#ifndef _NTT_IMPL_H
#define _NTT_IMPL_H

#include "sp.h"

#ifdef __cplusplus
extern "C" {
#endif

/* low-level stuff for executing one set of transforms, where
   the transform size, prime size, and arithmetic are all specified */

typedef void (*nttdata_init_t)(spv_t out, 
				sp_t p, sp_t d,
				sp_t primroot, sp_t order,
				sp_t perm);

/* packed transforms */

typedef void (*ntt_run_t)(
    			spv_t in, spv_size_t stride, spv_size_t dist,
    			spv_size_t num_transforms, sp_t p, 
			spv_t ntt_const);

typedef void (*ntt_pfa_run_t)(
    			spv_t in, spv_size_t stride, spv_size_t dist,
    			spv_size_t num_transforms, sp_t p, 
			spv_size_t cofactor, spv_t ntt_const);

typedef void (*ntt_twiddle_run_t)(
    			spv_t in, spv_size_t stride, spv_size_t dist,
    			spv_size_t num_transforms, sp_t p, 
			spv_t ntt_const, spv_t w);

/* interleaved transforms */

typedef void (*ntt_run_interleaved_t)(
    			spv_t in, spv_size_t stride, spv_size_t dist,
    			spv_size_t num_transforms, spv_t p, 
			spv_size_t vsize, spv_t ntt_const);

typedef void (*ntt_pfa_run_interleaved_t)(
    			spv_t in, spv_size_t stride, spv_size_t dist,
    			spv_size_t num_transforms, spv_t p, 
			spv_size_t cofactor, spv_size_t vsize, 
			spv_t ntt_const);

typedef void (*ntt_twiddle_run_interleaved_t)(
    			spv_t in, spv_size_t stride, spv_size_t dist,
    			spv_size_t num_transforms, spv_t p, 
			spv_size_t vsize, spv_t ntt_const,
			spv_t w);

typedef struct
{
  uint32_t size;
  uint32_t num_ntt_const;
  const uint8_t *fixed_ntt_const;
  nttdata_init_t nttdata_init;

  ntt_run_t ntt_run;
  ntt_pfa_run_t ntt_pfa_run;
  ntt_twiddle_run_t ntt_twiddle_run;

  ntt_run_interleaved_t ntt_run_interleaved;
  ntt_pfa_run_interleaved_t ntt_pfa_run_interleaved;
  ntt_twiddle_run_interleaved_t ntt_twiddle_run_interleaved;
} nttconfig_t;

/* a group of transforms that share the same compiler optimizations,
   prime size and instruction set */

typedef const nttconfig_t ** (*get_transform_list_t)(void);

typedef struct
{
  const char * name;
  uint32_t num_transforms;
  uint32_t vsize;
  get_transform_list_t get_transform_list;
} nttgroup_t;

extern const nttgroup_t * X(ntt_master_group_list)[];
uint32_t X(get_master_group_list_size)();

/* all the groups available */

extern const nttgroup_t X(ntt_group);
extern const nttgroup_t MANGLE_SSE2(X(ntt_group_simd));
extern const nttgroup_t MANGLE_SSE42(X(ntt_group_simd));
extern const nttgroup_t MANGLE_AVX(X(ntt_group_simd));
extern const nttgroup_t MANGLE_AVX2(X(ntt_group_simd));
extern const nttgroup_t MANGLE_FMA(X(ntt_group_simd));


/* if the modulus is 2 bits or more smaller than the machine
   word size, the core NTT routines use a redundant representation 
   of the transform elements. Modular multiplies do not do their
   final modular reduction, and inputs to the multiplies are not 
   reduced mod p either. The transform results are in [0, 2p) */

#if SP_NUMB_BITS == SP_TYPE_BITS - 2
#define HAVE_PARTIAL_MOD
#endif

/* an NTT is built up of one or more passes through
   the input data */

typedef struct
{
  pass_type_t pass_type;
  const nttconfig_t *codelet;
  spv_t codelet_const;
  spv_size_t const_size;
  spv_size_t num_transforms;
  spv_size_t vsize;

  union
  {
    struct
    {
    } direct;

    struct
    {
      spv_size_t ntt_size;
      spv_size_t cofactor;
    } pfa;

    struct
    {
      spv_size_t stride;
      spv_size_t twiddle_size;
      spv_t w;
    } twiddle;

  } d;

} nttpass_t;

/* central repository for all precomputed constants
   for computing a group of parallel NTTs */
typedef struct
{
  uint32_t num_passes;
  nttpass_t passes[MAX_PLANS];
} nttdata_t;

/* SPM */

/* small prime modulus - this contains some precomputed constants to
 * calculate modulo a sp */
typedef struct
{
  sp_t sp;		/* value of the sp */
  sp_t mul_c;		/* constant used for reduction mod sp */
  sp_t primroot;
  sp_t inv_primroot;
} __spm_struct;

typedef __spm_struct * spm_t;

spm_t X(spm_init)(uint32_t, sp_t);
void X(spm_clear)(spm_t);

/* MPZSPM */

typedef struct
{
  uint32_t max_vsize;
  uint32_t sp_num;
  uint32_t interleaved;
  spm_t * spm;

  spv_t work;
  spv_t sp;
  nttdata_t nttdata;
  nttdata_t inttdata;
} __mpzspm_struct;

typedef __mpzspm_struct * mpzspm_t;


/* internal routines everybody needs */

void
X(nttdata_init_generic)(const nttconfig_t *c,
		      spv_t out, sp_t p, sp_t d, 
		      sp_t primroot, sp_t order,
		      sp_t perm);

sp_t X(sp_ntt_reciprocal)(sp_t w, sp_t p);

/* routines for nttdata_t */

uint32_t X(ntt_build_passes)(nttplangroup_t *p, mpzspm_t mpzspm, 
    				uint32_t ntt_size, uint32_t max_ntt_size);

void X(ntt_reset)(nttdata_t *data);
void X(ntt_run)(void * m, mpz_t * x, uint32_t ntt_size);

#ifdef __cplusplus
}
#endif

#endif /* _NTT_IMPL_H */
