#ifndef BATCHED_S_H_
#define BATCHED_S_H_

#include "ecm_int.h"
#include "getprime_r.h"

struct batched_info_s {
  ecm_uint batch; /* index of batch */
  prime_info_s;   /* current prime iterator */
  ecm_uint prime; /* next prime to process */
  double B1;      /* B1 */
  ecm_uint B1_2;  /* B1^(1/2) */
  ecm_uint B1_3;  /* B1^(1/3) */
};
typedef struct batched_info_s prime_info_t[1];

#ifdef __cplusplus
extern "C" {
#endif

void batched_info_init (batched_info_t);
void batched_info_clear (batched_info_t);

/* Get next batch of primes in 
int  get_batch (prime_info_t);

#ifdef __cplusplus
}
#endif

#endif	/* BATCHED_S_H_ */
