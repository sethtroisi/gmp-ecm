/* See ecm-ecm.h for explainations */
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
unsigned int findfactor(mpz_t N, mpz_t xfin, mpz_t zfin);
void biguint_print (biguint_t a);
void bigint_print (dbigint_t a);
void mpz_to_biguint (biguint_t a, mpz_t b);
void biguint_to_mpz (mpz_t a, biguint_t b);

void write_resumefile2_line (FILE *savefile, mpz_t N, unsigned int B1, 
                                mpz_t xp, unsigned int firstinvd, mpz_t mpz_d);

/* Require mpcandi.c from trunk/ directory */
extern void mpcandi_t_init (mpcandi_t *n);  
extern void mpcandi_t_free (mpcandi_t *n);

/* Require auxi.c and eval.c from trunk/ directory */
extern int read_number (mpcandi_t *, FILE *, int);

/* Require batch.c and getprime.c from trunk/ directory*/ 
extern void compute_s (mpz_t s, unsigned int B1);

/* Require random.c from trunk/ directory*/ 
extern unsigned int get_random_ui (void);
