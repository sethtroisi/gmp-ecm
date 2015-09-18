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

typedef void (*ntt_run_t)(spv_t in, spv_size_t istride, spv_size_t idist,
    			spv_t out, spv_size_t ostride, spv_size_t odist,
    			spv_size_t num_transforms, sp_t p, spv_t ntt_const);

typedef void (*ntt_pfa_run_t)(spv_t x, spv_size_t cofactor, 
			  sp_t p, spv_t ntt_const);

typedef void (*ntt_twiddle_run_t)(
    			spv_t in, spv_size_t istride, spv_size_t idist,
    			spv_t out, spv_size_t ostride, spv_size_t odist,
    			spv_t w, spv_size_t num_transforms, sp_t p, spv_t ntt_const);

typedef struct
{
  uint32_t size;
  uint32_t num_ntt_const;
  const uint8_t *fixed_ntt_const;
  nttdata_init_t nttdata_init;
  ntt_run_t ntt_run;
  ntt_pfa_run_t ntt_pfa_run;
  ntt_twiddle_run_t ntt_twiddle_run;
} nttconfig_t;

/* a group of transforms that share the same compiler optimizations,
   prime size and instruction set */

typedef spv_t (*alloc_twiddle_t)(sp_t primroot, sp_t order, sp_t p, sp_t d,
				spv_size_t rows, spv_size_t cols);

typedef void (*free_twiddle_t)(spv_t twiddle);

typedef const nttconfig_t ** (*get_transform_list_t)(void);

typedef struct
{
  const char * name;
  uint32_t num_transforms;
  get_transform_list_t get_transform_list;
  alloc_twiddle_t alloc_twiddle;
  free_twiddle_t free_twiddle;
} nttgroup_t;

extern const nttgroup_t * X(ntt_master_group_list)[];
extern const uint32_t X(ntt_master_group_list_size);

/* all the groups available */

extern const nttgroup_t X(ntt_group);
extern const nttgroup_t MANGLE_SSE2(X(ntt_group_simd));
extern const nttgroup_t MANGLE_SSE42(X(ntt_group_simd));
extern const nttgroup_t MANGLE_AVX(X(ntt_group_simd));


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

typedef enum
{
  PASS_TYPE_DIRECT,
  PASS_TYPE_PFA,
  PASS_TYPE_TWIDDLE
} pass_type_t;

#define MAX_PASSES 10

typedef struct
{
  pass_type_t pass_type;
  const nttgroup_t *group;
  const nttconfig_t *codelet;
  spv_t codelet_const;

  union
  {
    struct
    {
      spv_size_t num_transforms;
    } direct;

    struct
    {
      spv_size_t cofactor;
    } pfa;

    struct
    {
      spv_size_t num_transforms;
      spv_size_t stride;
      spv_t w;
    } twiddle;

  } d;

} nttpass_t;

/* central repository for all NTT data that shares a
   modulus and primitive root */
typedef struct
{
  uint32_t num_passes;
  nttpass_t passes[MAX_PASSES];
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

  nttdata_t ntt_data;
} __spm_struct;

typedef __spm_struct * spm_t;

spm_t X(spm_init)(uint32_t, sp_t);
void X(spm_clear)(spm_t);

/* MPZSPM */

typedef struct
{
  uint32_t sp_num;
  spm_t * spm;
} __mpzspm_struct;

typedef __mpzspm_struct * mpzspm_t;

/* guides for constructing transforms */
typedef struct
{
  uint32_t codelet_size;
  uint32_t group_type;
  pass_type_t pass_type;
} nttplan_t;

/* internal routines everybody needs */

void
X(nttdata_init_generic)(const nttconfig_t *c,
		      spv_t out, sp_t p, sp_t d, 
		      sp_t primroot, sp_t order,
		      sp_t perm);

sp_t X(sp_ntt_reciprocal)(sp_t w, sp_t p);

/* external interface */

uint32_t X(ntt_build_passes)(nttdata_t *data,
    		nttplan_t *plans, uint32_t num_plans,
		sp_t size, sp_t p, sp_t primroot, sp_t order, sp_t d);

void * X(ntt_init)(sp_t size, sp_t primroot, sp_t p, sp_t d);
void X(ntt_free)(void *data);
void X(ntt_reset)(void *data);
void X(ntt_run)(spv_t x, sp_t p, void *data);

#ifdef __cplusplus
}
#endif

#endif /* _NTT_IMPL_H */
