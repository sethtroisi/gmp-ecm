/* See ecm-ecm.h for explainations */

#include "../../ecm-impl.h"

typedef struct
{
  char *cpExpr;	
  mpz_t n;	
  unsigned ndigits;
  unsigned nexprlen;
  int isPrp;	
} mpcandi_t;

void usage (void);
long cputime();
void to_mont_repr (mpz_t, mpcandi_t);
void from_mont_repr (mpz_t, mpcandi_t, mpz_t);
unsigned int findfactor (mpcandi_t n, mpz_t xfin, mpz_t zfin);
void print_factor_cofactor (mpcandi_t n, mpz_t factor);
void biguint_print (biguint_t a);
void bigint_print (dbigint_t a);
void mpz_to_biguint (biguint_t a, mpz_t b);
void biguint_to_mpz (mpz_t a, biguint_t b);
void write_resumefile_wrapper (char *, mpcandi_t *, unsigned int , mpz_t,
                               unsigned int, mpz_t);

/* Require resume.c from trunk/ directory */
int write_resumefile_line (char *, int, double, mpz_t, mpz_t, mpz_t, 
                           mpcandi_t *, mpz_t, const char *);

/* Require mpcandi.c from trunk/ directory */
extern void mpcandi_t_init (mpcandi_t *n);  
extern void mpcandi_t_free (mpcandi_t *n);

/* Require auxi.c and eval.c from trunk/ directory */
extern int read_number (mpcandi_t *, FILE *, int);

/* compute_s (mpz_t s, unsigned int B1); */
/* unsigned int get_random_ui (void);    */
/*      are defined in ecm-impl.h        */
