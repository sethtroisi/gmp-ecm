#ifndef GETPRIME_R_H_
#define GETPRIME_R_H_

#include "ecm_int.h"

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

#ifdef __cplusplus
extern "C" {
#endif

/* The getprime_mt function returns successive odd primes, starting with 3. */
void prime_info_init (prime_info_t);
void prime_info_clear (prime_info_t);
ecm_uint getprime_mt (prime_info_t);

#ifdef __cplusplus
}
#endif

#endif	/* GETPRIME_R_H_ */
