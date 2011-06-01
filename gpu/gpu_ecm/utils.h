#include <gmp.h>

#define NV 10
#define ADD 6.0 // number of multiplications in an addition 
#define DUP 5.0 // number of multiplications in a duplicate 

unsigned int findfactor(mpz_t N, mpz_t xfin, mpz_t zfin);
void biguint_print (biguint_t a);
void bigint_print (dbigint_t a);
void mpz_to_biguint (biguint_t a, mpz_t b);
void biguint_to_mpz (mpz_t a, biguint_t b);
unsigned long getprime (unsigned long pp);


