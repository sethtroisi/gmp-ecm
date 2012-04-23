#ifndef _NTT_IMPL_H
#define _NTT_IMPL_H

#include "sp.h"

typedef uint32_t (*get_num_ntt_const_t)(void);

typedef void (*nttdata_init_t)(spv_t out, 
				sp_t p, sp_t d,
				sp_t primroot, sp_t order);

typedef void (*ntt_run_t)(spv_t x, spv_size_t stride,
			  sp_t p, sp_t d, spv_t ntt_const);

typedef void (*ntt_pfa_run_t)(spv_t x, spv_size_t stride,
			  spv_size_t cofactor, 
			  sp_t p, sp_t d, spv_t ntt_const);


/* a copy of sp_add, but operating on array offsets */

static inline spv_size_t sp_array_inc(spv_size_t a, spv_size_t b, spv_size_t m) 
{
#if (defined(__GNUC__) || defined(__ICL)) && \
    (defined(__x86_64__) || defined(__i386__))

  spv_size_t t = a - m, tr = a + b;

  __asm__ (
    "add %2, %1    # sp_array_inc: t += b\n\t"
    "cmovc %1, %0  # sp_array_inc: if (cy) tr = t \n\t"
    : "+r" (tr), "+&r" (t)
    : "g" (b)
    : "cc"
  );

  return tr;

#elif defined(_MSC_VER) && !defined(_WIN64)

  __asm
    {
        mov     eax, a
        add     eax, b
        mov     edx, eax
        sub     edx, m
        cmovnc  eax, edx
    }

#else

  spv_size_t t = a + b;
  if (t >= m)
    t -= m;
  return t;

#endif
}



typedef struct
{
  uint32_t size;
  get_num_ntt_const_t get_num_ntt_const;
  nttdata_init_t nttdata_init;
  ntt_run_t ntt_run;
  ntt_pfa_run_t ntt_pfa_run;
} nttconfig_t;

extern const nttconfig_t ntt2_config;
extern const nttconfig_t ntt3_config;
extern const nttconfig_t ntt4_config;
extern const nttconfig_t ntt5_config;
extern const nttconfig_t ntt6_config;
extern const nttconfig_t ntt7_config;
extern const nttconfig_t ntt8_config;
extern const nttconfig_t ntt9_config;
extern const nttconfig_t ntt15_config;

#endif /* _NTT_IMPL_H */
