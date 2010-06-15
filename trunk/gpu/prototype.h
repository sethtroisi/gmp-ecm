# ifndef PROTOTYPE_H
# define PROTOTYPE_H

void calculParam(mpz_t sigma,mpz_t N,mpz_t d,mpz_t x0,mpz_t z0,mpz_t u,mpz_t v);
void duplicate(mpz_t xp,mpz_t zp,mpz_t d,mpz_t N,mpz_t x2p,mpz_t z2p);
void addition(mpz_t xP,mpz_t zP,mpz_t xQ,mpz_t zQ,mpz_t x_PminusQ,mpz_t z_PminusQ,mpz_t N,mpz_t x_PplusQ,mpz_t z_PplusQ);
void mongomerysChains(mpz_t PI);
void multiplication(unsigned long PI,mpz_t xQ,mpz_t zQ,mpz_t d,mpz_t N,mpz_t sigma,mpz_t x_PIQ,mpz_t z_PIQ);
unsigned long int eratosthene(unsigned long min,mpz_t max);
void stageOne(mpz_t B1,mpz_t sigma,mpz_t N,mpz_t Q);

# endif
