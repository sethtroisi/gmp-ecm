typedef double* dpn_t;

unsigned long dpn_size (mpz_t N);
void conversion64to52 (dpn_t b, mpz_t N, unsigned long n);
void conversion52to64 (dpn_t b, unsigned long n, mpz_t M);
void dpn_set (dpn_t a, dpn_t b, unsigned long n);
double dpn_add (dpn_t a, dpn_t b, dpn_t c, unsigned long n);
double dpn_sub (dpn_t a, dpn_t b, dpn_t c, unsigned long n);
void dpn_mul (dpn_t a, dpn_t b, dpn_t c, unsigned long n);
void dpn_mod (dpn_t a, dpn_t mod, dpn_t mu, unsigned long n);
void dpn_add_mod (dpn_t amod, dpn_t b, dpn_t c, dpn_t mod, unsigned long n);
void dpn_sub_mod (dpn_t smod, dpn_t b, dpn_t c, dpn_t mod, unsigned long n);
void mpz_to_montgomery (mpz_t a, mpz_t N, unsigned long b);
void mpz_from_montgomery (mpz_t a, mpz_t N, unsigned long b);
void dpn_print (dpn_t a, unsigned long n);
int dpn_check (dpn_t a, dpn_t N, unsigned long n);
double dpn_add_1 (dpn_t a, dpn_t b, unsigned long n, double cy);
void dpn_zero (dpn_t a, unsigned long n);
