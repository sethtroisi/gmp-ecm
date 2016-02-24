#ifndef GETPRIME_R_H_
#define GETPRIME_R_H_

struct prime_info_s {
  unsigned long offset;  /* offset for current primes */
  long current;          /* index of previous prime */
  unsigned int *primes;  /* small primes up to sqrt(p) */
  unsigned long nprimes; /* length of primes[] */
  unsigned char *sieve;  /* sieving table */
  long len;              /* length of sieving table */
  unsigned int *moduli;  /* offset for small primes */
};
typedef struct prime_info_s prime_info_t[1];

#ifdef __cplusplus
extern "C" {
#endif

/* The getprime function returns successive odd primes, starting with 3. */
void prime_info_init (prime_info_t);
void prime_info_clear (prime_info_t);
unsigned long getprime_mt (prime_info_t);
unsigned long getprime_unsafe (unsigned long);


#ifdef __cplusplus
}
#endif

#endif	/* GETPRIME_R_H_ */
