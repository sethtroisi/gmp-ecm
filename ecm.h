/* maximal stage 1 bound = 2^53 + 4, the next prime being 2^53 + 5 */
#define MAX_B1 9007199254740996.0

/* different methods implemented */
#define EC_METHOD 0
#define PM1_METHOD 1
#define PP1_METHOD 2

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
int      pm1_stage1     (mpz_t, mpz_t, double);
int          pm1        (mpz_t, mpz_t, double, double, unsigned int, 
                         unsigned int, int);
int          ecm        (mpz_t, mpz_t, mpz_t, double, double, unsigned int, 
                         unsigned int, int);
unsigned int phi        (unsigned int);
unsigned int bestD      (double);
double       block_size (unsigned int);
unsigned int gcd        (unsigned int, unsigned int);
int          cputime    (void);

/* pp1.c */
void         pp1_mul    (mpz_t, mpz_t, mpz_t, mpz_t, mpz_t, mpz_t);
void         pp1_mul_ui (mpz_t, mpz_t, mpz_t, mpz_t, mp_limb_t, mpz_t);
int          pp1_stage1 (mpz_t, mpz_t, double);
int          pp1        (mpz_t, mpz_t, double, double, unsigned int, 
                         unsigned int, int);

/* stage2.c */
int          rootsF     (listz_t, unsigned int, mpz_t, mpz_t, mpz_t *, 
                         unsigned int , mpz_t, int, int);
void         rootsG     (listz_t, unsigned int, listz_t, listz_t, 
                         listz_t, unsigned int, mpz_t, int);
int          stage2     (mpz_t, mpz_t, double, unsigned int, unsigned int, int,
                         int, int);

/* listz.c */
int          list_mul_mem (unsigned int);
listz_t      init_list  (unsigned int);
void         clear_list (listz_t, unsigned int);
void         print_list (listz_t, unsigned int);
void         list_set   (listz_t, listz_t, unsigned int);
void         list_neg   (listz_t, listz_t, unsigned int);
void         list_mod   (listz_t, listz_t, unsigned int, mpz_t);
void         list_add   (listz_t, listz_t, listz_t, unsigned int);
void         list_sub   (listz_t, listz_t, listz_t, unsigned int);
void         list_zero  (listz_t, unsigned int);
int          list_zerop (listz_t, unsigned int);
int          karatsuba  (listz_t, listz_t, listz_t, unsigned int, listz_t);
int          list_mul   (listz_t, listz_t, unsigned int, listz_t, unsigned int,
			 listz_t);
int         list_mulmod (listz_t, listz_t, listz_t, unsigned int, listz_t,
			 mpz_t);
int         list_mulmod2(listz_t, listz_t, listz_t, listz_t, unsigned int,
                         listz_t, mpz_t);
int       PolyFromRoots (listz_t, unsigned int, listz_t, int, mpz_t, char);
int          PolyInvert (listz_t, listz_t, unsigned int, listz_t, mpz_t);
int   RecursiveDivision (listz_t, listz_t, listz_t, unsigned int, listz_t,
			 mpz_t);
int   PrerevertDivision (listz_t, listz_t, listz_t, unsigned int, listz_t,
                         mpz_t);
void         Div3by2    (listz_t, listz_t, listz_t, unsigned int, listz_t,
			 mpz_t);
int          list_mod1  (mpz_t, listz_t, listz_t, unsigned int, mpz_t, mpz_t*);
void      poly_submul2 (listz_t, listz_t, listz_t, unsigned int, mpz_t, mpz_t);
int          list_invert (listz_t, listz_t, unsigned int, mpz_t, mpz_t);

/* toomcook.c */
void     mpz_divby3_1op (mpz_t RS);
int           toomcook3 (listz_t, listz_t, listz_t, unsigned int, listz_t);
int           toomcook4 (listz_t, listz_t, listz_t, unsigned int, listz_t);

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

/* a <- b * c where a and b are mpz, c is a double, and t an auxiliary mpz */
#if (BITS_PER_MP_LIMB >= 53)
#define mpz_mul_d(a, b, c, t) \
   mpz_mul_ui (a, b, (unsigned long int) c);
#else
#if (BITS_PER_MP_LIMB >= 32)
#define mpz_mul_d(a, b, c, t) \
   if (c < 4294967296.0) \
      mpz_mul_ui (a, b, (unsigned long int) c); \
   else { \
   mpz_set_d (t, c); \
   mpz_mul (a, b, t); }
#else
#define mpz_mul_d(a, b, c, t) \
   mpz_set_d (t, c); \
   mpz_mul (a, b, t);
#endif
#endif

