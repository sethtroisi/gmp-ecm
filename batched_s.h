#ifndef BATCHED_S_H_
#define BATCHED_S_H_

#include "basicdefs.h"
#include "ecm_int.h"
#include "getprime_r.h"
#include <gmp.h>

/* Keep this many iterators */
#define BATCHED_ITERATORS 9
/* Need to keep small primes up to LIMIT_B1 ^ (1/(ITERATORS+1))
   (2^64) ** (1/10) = 84, last prime is 83
 */
#define SMALL_PRIMES 23

/* batched_info_s is made up of several prime iterators
   with each iterator tracking the current prime */
struct prime_iterator_s {
  uint64_t current_prime;
  prime_info_t p_i;
};
typedef struct prime_iterator_s prime_iterator_t[1];

typedef struct {
  unsigned int size;
  mpz_t *val;
} batched_mul_casc;


struct batched_info_s {
  uint32_t B1;
  /* q is p^current_power - 1
     -1 is important so that UINT_MAX can be a sentinel */
  uint64_t small_q[SMALL_PRIMES];
  struct prime_iterator_s iterator[BATCHED_ITERATORS];
  batched_mul_casc *cascade;
};
typedef struct batched_info_s batched_info_t[1];



#ifdef __cplusplus
extern "C" {
#endif

void batched_info_init (batched_info_t);
void batched_info_clear (batched_info_t);
void advance_to (batched_info_t, uint64_t);
void get_batch (batched_info_t, mpz_t, uint64_t);

#ifdef __cplusplus
}
#endif

#endif	/* BATCHED_S_H_ */
