#include <gmp.h>

//#define NV 10
//#define ADD 6.0 // number of multiplications in an addition 
//#define DUP 5.0 // number of multiplications in a duplicate 

void usage (void);
long cputime();
unsigned int findfactor(mpz_t N, mpz_t xfin, mpz_t zfin);
void biguint_print (biguint_t a);
void bigint_print (dbigint_t a);
void mpz_to_biguint (biguint_t a, mpz_t b);
void biguint_to_mpz (mpz_t a, biguint_t b);
unsigned long getprime (unsigned long pp);
void compute_s (mpz_t s, unsigned int B1);
void write_resumefile_line (FILE *savefile, mpz_t N, unsigned int B1, mpz_t xp, 
                  unsigned int firstinvd, mpz_t mpz_d);
int read_number (mpz_t n, FILE *fd);
