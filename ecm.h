typedef mpz_t* listz_t;

typedef struct
{
  int alloc;
  int degree;
  listz_t coeff;
} __polyz_struct;
typedef __polyz_struct polyz_t[1];

double   getprime       (double);
unsigned int nb_digits  (mpz_t);
int          stage1     (mpz_t, mpz_t, double);
int          pm1        (mpz_t, mpz_t, double, double, unsigned int, 
                         unsigned int, int);
int          ecm        (mpz_t, mpz_t, double, double, unsigned int, 
                         unsigned int, int);
unsigned int phi        (unsigned int);
unsigned int bestD      (double);
double       block_size (unsigned int);
unsigned int gcd        (unsigned int, unsigned int);
int          cputime    (void);

/* stage2.c */
int          rootsF     (listz_t, unsigned int, mpz_t, mpz_t, mpz_t *, 
                         unsigned int , mpz_t, int);
void         rootsG     (listz_t, unsigned int, listz_t, listz_t, 
                         listz_t, unsigned int, mpz_t, int);
int          stage2     (mpz_t, mpz_t, double, unsigned int, unsigned int, int);

/* listz.c */
int          list_mul_mem (unsigned int);
listz_t      init_list  (unsigned int);
void         clear_list (listz_t, unsigned int);
void         print_list (listz_t, unsigned int);
void         list_set   (listz_t, listz_t, unsigned int);
void         list_mod   (listz_t, listz_t, unsigned int, mpz_t);
void         list_add   (listz_t, listz_t, listz_t, unsigned int);
void         list_sub   (listz_t, listz_t, listz_t, unsigned int);
int          list_zerop (listz_t, unsigned int);
void         karatsuba  (listz_t, listz_t, listz_t, unsigned int, listz_t);
int          toomcook3  (listz_t, listz_t, listz_t, unsigned int, listz_t);
void         list_mul   (listz_t, listz_t, unsigned int, listz_t, unsigned int,
			 listz_t);
void        list_mulmod (listz_t, listz_t, listz_t, unsigned int, listz_t,
			 mpz_t);
void         buildG     (listz_t, unsigned int, listz_t, int, mpz_t, char);
void  RecursiveDivision (listz_t, listz_t, listz_t, unsigned int, listz_t,
			 mpz_t);
void         Div3by2    (listz_t, listz_t, listz_t, unsigned int, listz_t,
			 mpz_t);
int          list_mod1  (mpz_t, listz_t, listz_t, unsigned int, mpz_t, mpz_t*);
void      poly_submul2 (listz_t, listz_t, listz_t, unsigned int, mpz_t, mpz_t);
int          hgcd_naive (mpz_t, listz_t, listz_t, listz_t, listz_t,
                         listz_t, listz_t, unsigned int, mpz_t, listz_t, int);
int          hgcd       (mpz_t, listz_t, listz_t, listz_t, listz_t,
                         listz_t, listz_t, unsigned int, mpz_t, listz_t, int);
int          list_gcd   (mpz_t, listz_t, listz_t, unsigned int, mpz_t, listz_t);
int          list_invert(listz_t, listz_t, unsigned int, mpz_t, mpz_t);

/* polyz.c */
void init_poly      (polyz_t, int);
void init_poly_list (polyz_t, int, listz_t);
void clear_poly     (polyz_t);
int  poly_zerop     (polyz_t);
void poly_set_ui    (polyz_t, unsigned long int);
int  poly_gcd       (mpz_t, polyz_t, polyz_t, mpz_t, listz_t);
int  poly_mod1      (mpz_t, polyz_t, polyz_t, mpz_t, listz_t);

/* ntl.c */
int NTL_major_version (void);
int NTL_minor_version (void);
int ntl_poly_gcd   (mpz_t, polyz_t, polyz_t, mpz_t);
void NTL_init (void);
void NTL_clear (void);
void NTL_get_factor (mpz_t);
