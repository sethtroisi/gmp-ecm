#ifndef GETPRIME_R_H_
#define GETPRIME_R_H_

#include "ecm_int.h"

#ifdef HAVE_LIBPRIMESIEVE
#include <primesieve.h>
typedef primesieve_iterator prime_info_t[1];
#else
struct prime_info_s {
  ecm_uint offset;  /* offset for current primes */
  ecm_int current;          /* index of previous prime */
  ecm_uint *primes;  /* small primes up to sqrt(p) */
  ecm_uint nprimes; /* length of primes[] */
  unsigned char *sieve;  /* sieving table */
  ecm_int len;              /* length of sieving table */
  ecm_uint *moduli;  /* offset for small primes */
};
typedef struct prime_info_s prime_info_t[1];
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* The getprime_mt function returns successive odd primes, starting with 3. */
void prime_info_init (prime_info_t);
void prime_info_clear (prime_info_t);
ecm_uint getprime_mt (prime_info_t);
/* Returns the next prime >= n. */
ecm_uint getprime_jump_and_next_mt (prime_info_t, ecm_uint n);

#ifdef __cplusplus
}
#endif

#endif	/* GETPRIME_R_H_ */
