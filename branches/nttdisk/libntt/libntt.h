#ifndef _LIBNTT_H
#define _LIBNTT_H

#include "basicdefs.h"
#include "gmp-xface.h"

void * sp_aligned_malloc (size_t len);
void sp_aligned_free (void *newptr);

/* all the fat binary routines */

#define DECLARE(spbits, wbits) \
void * mpzspm_init_sp##spbits##w##wbits(uint32_t max_len_in, mpz_t modulus); \
void mpzspm_clear_sp##spbits##w##wbits(void *mpzspm); \
int test_main_sp##spbits##w##wbits(int argc, char **argv); \

#if GMP_LIMB_BITS == 32
DECLARE(30,32)
DECLARE(31,32)
DECLARE(62,32)
#else
DECLARE(30,64)
DECLARE(62,64)
#endif
DECLARE(50,64)

#endif /* _LIBNTT_H */
