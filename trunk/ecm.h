/* ecm.h - header file for gmp-ecm
 
  Copyright 2001, 2002, 2003 Paul Zimmermann and Alexander Kruppa.
 
  This program is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 2 of the License, or (at your
  option) any later version.
 
  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
  more details.
 
  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.
*/

#define ECM_VERSION "5.1.2-beta"

/* To use new polyeval_tellegen */
#define POLYEVALTELLEGEN

/* define top-level multiplication */
#define KARA 2
#define TOOM3 3
#define TOOM4 4
/* compile with -DMULT=2 to override default */
#ifndef MULT
#define MULT TOOM4
#endif

#ifdef POLYEVALTELLEGEN
#define POLYEVAL
#define USE_SHORT_PRODUCT
#endif

#ifndef POLYGCD
#define POLYEVAL
#endif

#include <stdio.h>
#include <gmp.h>

/* thresholds */
#ifndef WANT_GMP_IMPL
#ifndef MUL_KARATSUBA_THRESHOLD
#define MUL_KARATSUBA_THRESHOLD 32
#endif
#endif

#ifndef DIV_DC_THRESHOLD
#define DIV_DC_THRESHOLD    (3 * MUL_KARATSUBA_THRESHOLD)
#endif

#define MPZMOD_THRESHOLD_DEFAULT (3 * DIV_DC_THRESHOLD / 2)
#define REDC_THRESHOLD_DEFAULT   (2 * DIV_DC_THRESHOLD)

/* default number of probable prime tests */
#define PROBAB_PRIME_TESTS 1

#define ABS(x) ((x) >= 0 ? (x) : -(x))

#define mpz_mulmod(a,b,c,n) \
        { mpz_mul (a, b, c); \
        mpz_mod (a, a, n); }

/* maximal stage 1 bound = 2^53 + 4, the next prime being 2^53 + 5 */
#define MAX_B1 9007199254740996.0

/* The checksum for savefile is the product of all mandatory files, modulo
   the greatest prime below 2^32 */
#define CHKSUMMOD 4294967291U

/* different methods implemented */
#define EC_METHOD 0
#define PM1_METHOD 1
#define PP1_METHOD 2

#define MOD_PLAIN 0
#define MOD_BASE2 1
#define MOD_MODMULN 2
#define MOD_REDC 3

typedef mpz_t mpres_t;

typedef mpz_t* listz_t;

typedef struct
{
  mpres_t x;
  mpres_t y;
} __point_struct;
typedef __point_struct point;

typedef struct
{
  mpres_t x;
  mpres_t y;
  mpres_t A;
} __curve_struct;
typedef __curve_struct curve;

typedef struct
{
  int alloc;
  int degree;
  listz_t coeff;
} __polyz_struct;
typedef __polyz_struct polyz_t[1];

typedef struct 
{
  int repr;           /* 0: plain modulus, possibly normalized
                         1: base 2 number
                         2: MODMULN
                         3: REDC representation */
  int bits;           /* in case of a base 2 number, 2^k[+-]1, bits = [+-]k
                         in case of MODMULN or REDC representation, nr. of 
                         bits b so that 2^b > orig_modulus and 
                         mp_bits_per_limb | b */
  mp_limb_t Nprim;    /* For MODMULN */
  mpz_t orig_modulus; /* The original modulus */
  mpz_t aux_modulus;  /* The auxiliary modulus value (i.e. normalized 
                         modulus, or -1/N (mod 2^bits) for REDC */
  mpz_t multiple;     /* The smallest multiple of N that is larger or
			 equal to 2^bits for REDC/MODMULN */
  mpz_t R2, R3;       /* For MODMULN and REDC, R^2 and R^3 (mod orig_modulus), 
                         where R = 2^bits. */
  mpz_t temp1, temp2; /* Temp values used during multiplication etc. */
} __mpmod_struct;
typedef __mpmod_struct mpmod_t[1];

/* Structure for candidate usage.  This is much more powerful than using a
   simple mpz_t to hold the candidate.  This structure also houses the 
   expression (in raw form), and will modify the expression as factors 
   are found (if in looping modes).  Also, since we are warehousing all
   of the data associated with the candidate, we also store whether the
   candidate is PRP here (so testing will cease), along with the length
   of the candidate.  As each factor is found, the candidate will also
   have the factor removed from it */
typedef struct
{
#if defined (CANDI_DEBUG)
  unsigned long magic;	/* used for debugging purposes while writing this code */
#endif
  char *cpExpr;		/* if non-NULL, then this is a "simpler" expression than the 
			   decimal output of n */
  mpz_t n;		/* the cofactor canidate currently being used to find factors from */
  unsigned ndigits;	/* the number of digits (decimal) in n */
  unsigned nexprlen;	/* strlen of expression, 0 if there is NO expression */
  int isPrp;		/* usually 0, but turns 1 if factor found, and the cofactor is PRP, 
			   OR if the original candidate was PRP and the user asked to prp check */
} mpcandi_t;

#if defined (__cplusplus)
extern "C" {
#endif  

double   getprime       (double);

/* auxi.c */
unsigned int nb_digits  (const mpz_t);
unsigned int gcd        (unsigned int, unsigned int);
void         mpz_sub_si (mpz_t, mpz_t, int);
unsigned int ceil_log2  (unsigned int);
int          cputime    ();
unsigned int get_random_ui (void);

/* pm1.c */
void    pm1_random_seed  (mpz_t, mpz_t, gmp_randstate_t);
int          pm1         (mpz_t, mpz_t, mpz_t, double, double, double, double, double, 
                          unsigned int, int, int, int, mpz_t);
int     pm1_rootsF       (mpz_t, listz_t, unsigned int, mpres_t *, listz_t,
                          int, mpmod_t, int, unsigned long *);
mpres_t *pm1_rootsG_init (mpres_t *, double, unsigned int, int, mpmod_t);
void    pm1_rootsG_clear (mpres_t *, int, mpmod_t);
int     pm1_rootsG       (mpz_t, listz_t, unsigned int, mpres_t *, listz_t, 
                          int, mpmod_t, int, unsigned long *);



/* ecm.c */
int          ecm        (mpz_t, mpz_t, mpz_t, mpz_t, double, double, double, double,
                         double, unsigned int, int, int, int, int);
int          cputime    (void);

/* bestd.c */
unsigned long muls_toom3 (unsigned int);
unsigned long muls_gen (unsigned int);
unsigned long muls_gen_short (unsigned int);
unsigned long   phi (unsigned long);
double   block_size (unsigned long);
unsigned long bestD (double, unsigned int, unsigned int *, unsigned int,
                     unsigned long *);

/* trial.c */
int trial_factor (mpcandi_t *n, double maxfact, int deep);

/* ecm2.c */
int     ecm_rootsF       (mpz_t, listz_t, unsigned int, curve *,
                          int, mpmod_t, int, unsigned long *);
point * ecm_rootsG_init  (mpz_t, curve *, double, unsigned int, 
                          int, mpmod_t, int);
void    ecm_rootsG_clear (point *, int, mpmod_t);
int     ecm_rootsG       (mpz_t, listz_t, unsigned int, point *,
                          int, mpmod_t, int, unsigned long *);

/* pp1.c */
int          pp1_mul     (mpres_t, mpres_t, mpz_t, mpmod_t, mpres_t, mpres_t);
int          pp1_mul_prac(mpres_t, unsigned long, mpmod_t, mpres_t, mpres_t,
                          mpres_t, mpres_t, mpres_t);
void    pp1_random_seed  (mpz_t, mpz_t, gmp_randstate_t);
int          pp1         (mpz_t, mpz_t, mpz_t, double, double, double, double, double, 
                          unsigned int, unsigned int, int, int);
int   pp1_rootsF         (listz_t, unsigned int, mpres_t *, listz_t,
                          mpmod_t, int, unsigned long *);
mpres_t *pp1_rootsG_init (mpres_t *, double, unsigned int, mpmod_t);
void  pp1_rootsG_clear   (mpres_t *, mpmod_t);
int   pp1_rootsG         (listz_t, unsigned int, mpres_t *, mpmod_t,
                          unsigned long *);

/* stage2.c */
int          stage2     (mpz_t, void *, mpmod_t, double, double, unsigned int,
                         int, int, int);
void  fin_diff_coeff    (listz_t, double, unsigned int, unsigned int, int);

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
void         list_mul_z (listz_t, listz_t, mpz_t, unsigned int, mpz_t);
int          list_gcd   (mpz_t, listz_t, unsigned int, mpz_t);
void         list_zero  (listz_t, unsigned int);
int          list_zerop (listz_t, unsigned int);
int       list_mul_high (listz_t, listz_t, listz_t, unsigned int, listz_t);
int       toomcook4_low (listz_t, listz_t, listz_t, unsigned int, listz_t);
int      toomcook4_high (listz_t, listz_t, listz_t, unsigned int, listz_t);
int          karatsuba  (listz_t, listz_t, listz_t, unsigned int, listz_t);
int          list_mul   (listz_t, listz_t, unsigned int, int, listz_t,
                         unsigned int, int, listz_t);
int         list_mulmod (listz_t, listz_t, listz_t, unsigned int, listz_t,
			 mpz_t);
int         list_mulmod2(listz_t, listz_t, listz_t, listz_t, unsigned int,
                         listz_t, mpz_t);
int       PolyFromRoots (listz_t, listz_t, unsigned int, listz_t, int, mpz_t,
                         char, listz_t*, unsigned int);
int          PolyInvert (listz_t, listz_t, unsigned int, listz_t, mpz_t);
int   RecursiveDivision (listz_t, listz_t, listz_t, unsigned int,
                         listz_t, mpz_t, int);
int   PrerevertDivision (listz_t, listz_t, listz_t, unsigned int, listz_t,
                         mpz_t);
int          list_mod1  (mpz_t, listz_t, listz_t, unsigned int, mpz_t, mpz_t*);
void      poly_submul2 (listz_t, listz_t, listz_t, unsigned int, mpz_t, mpz_t);
int          list_invert (listz_t, listz_t, unsigned int, mpz_t, mpz_t);

/* polyeval.c */
unsigned int polyeval (listz_t, unsigned int, listz_t*, listz_t, mpz_t, int,
               unsigned int);
unsigned int
polyeval_tellegen (listz_t b, unsigned int k, listz_t *Tree, listz_t T,
                   unsigned int sizeT, listz_t invF, mpz_t n, unsigned int sh);
unsigned int TUpTree (listz_t b, listz_t *Tree, unsigned int k,
              listz_t tmp, unsigned int sh, mpz_t n);
unsigned int TUpTree_space (unsigned int k);
unsigned int muls_tuptree (unsigned int k);
unsigned int muls_polyeval_tellegen (unsigned int k);

/* toomcook.c */
void     mpz_divby3_1op (mpz_t);
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

/* mpmod.c */
int isbase2 (mpz_t, double);
void mpmod_init (mpmod_t, mpz_t, int, int);
void mpmod_init_MPZ (mpmod_t, mpz_t);
void mpmod_init_BASE2 (mpmod_t, int, mpz_t);
void mpmod_init_MODMULN (mpmod_t, mpz_t);
void mpmod_init_REDC (mpmod_t, mpz_t);
void mpmod_clear (mpmod_t);
void mpres_pow(mpres_t, mpres_t, mpres_t, mpmod_t);
void mpres_ui_pow (mpres_t, unsigned int, mpres_t, mpmod_t);
void mpres_mul(mpres_t, mpres_t, mpres_t, mpmod_t);
void mpres_div_2exp(mpres_t, mpres_t, unsigned int, mpmod_t);
void mpres_add_ui (mpres_t, mpres_t, unsigned int, mpmod_t);
void mpres_add (mpres_t, mpres_t, mpres_t, mpmod_t);
void mpres_sub_ui (mpres_t, mpres_t, unsigned int, mpmod_t);
void mpres_sub (mpres_t, mpres_t, mpres_t, mpmod_t);
void mpres_set_z (mpres_t, mpz_t, mpmod_t);
void mpres_get_z (mpz_t, mpres_t, mpmod_t);
void mpres_set_ui (mpres_t, unsigned int, mpmod_t);
void mpres_init (mpres_t, mpmod_t);
void mpres_realloc (mpres_t, mpmod_t);
void mpres_mul_ui (mpres_t, mpres_t, unsigned int, mpmod_t);
void mpres_neg (mpres_t, mpres_t, mpmod_t);
int  mpres_invert (mpres_t, mpres_t, mpmod_t);
void mpres_gcd (mpz_t, mpres_t, mpmod_t);
void mpres_out_str (FILE *, unsigned int, mpres_t, mpmod_t);
#define mpres_clear(a,n) mpz_clear (a)
#define mpres_set(a,b,n) mpz_set (a, b)
#define mpres_swap(a,b,n) mpz_swap (a, b)
#define mpres_is_zero(n) (mpz_sgn (n) == 0)

/* mul_lo.c */
void mpn_mul_lo_n (mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);

/* resume.c */
int  facceptstr (FILE *, char *);
int  freadstrn (FILE *, char *, char, unsigned int);
int  read_resumefile_line (int *, mpz_t, mpcandi_t *, mpz_t, mpz_t, mpz_t, double *,
                           char *, char *, char *, char *, FILE *);
void write_resumefile_line (FILE *, int, double, mpz_t, mpz_t, mpz_t, mpcandi_t *, 
                            mpz_t, const char *);
void write_temp_resumefile (int method, double B1, mpz_t sigma, mpz_t A, mpz_t x, mpz_t n, mpz_t orig_X0, int);
void kill_temp_resume_file (void);

/* main.c */
int read_number (mpcandi_t *n, FILE *, int primetest);
void showscreenticks(int stage, int percentage);   /* for outputting (or not outputting) the 1:98  or stage:percentage_complete  -Q will turn this off */
void showscreenticks_change_stage(int stage);      /* puts up 1:00 or 2:00 and resets the "timeout" clock.  */

/* eval.c */
int eval (mpcandi_t *n, FILE *fd, int bPrp);
int eval_str (mpcandi_t *n, char *cp, int primetest, char **EndChar); /* EndChar can be NULL */
void init_expr (void);
void free_expr (void);

/* candi.c */
void mpcandi_t_init (mpcandi_t *n);  /* damn, a C++ class sure would have been nice :(  */
void mpcandi_t_free (mpcandi_t *n);
int  mpcandi_t_copy (mpcandi_t *to, mpcandi_t *from);
int  mpcandi_t_add_candidate (mpcandi_t *n, mpz_t c, const char *cpExpr, int bPrp);
int  mpcandi_t_addfoundfactor (mpcandi_t *n, mpz_t f, int displaywarning);
int  mpcandi_t_addfoundfactor_d (mpcandi_t *n, double f);

/* smartprp.c */
int smart_probab_prime_p(mpz_t const n, int c);

/* Options for using an external PRPer rather than internal GMP */
extern char *externalprp; /* call out to external program  */  
extern char *externallog; /* where to find output */
extern char *externalinputprpfile; /* Used to place the n value (a temp file). Is deleted after system */
extern char *externalisprp; /* what to match output against */
extern char *externaliscomposite; /* ditto */
extern int externalprplen; /* length where we go external */
extern int externalprpval; /* exit value meaning it's prp, -1 is impossible */

/* b1_ainc.c */
double calc_B1_AutoIncrement(double cur_B1, double incB1val, int calcInc);

/* median.c */
unsigned int
TToomCookMul (listz_t, unsigned int, listz_t, unsigned int, listz_t, 
         unsigned int, listz_t);
unsigned int
TMulGen (listz_t, unsigned int, listz_t, unsigned int, listz_t, 
         unsigned int, listz_t);
unsigned int
TKarMul (listz_t, unsigned int, listz_t, unsigned int, listz_t, 
         unsigned int, listz_t);
void list_add_wrapper (listz_t, listz_t, listz_t, unsigned int,
                       unsigned int);
void list_sub_wrapper (listz_t, listz_t, listz_t, unsigned int,
                       unsigned int);
unsigned int
TKarMul_space (unsigned int n, unsigned int m, unsigned int l);
unsigned int
TMulGen_space (unsigned int n, unsigned int m, unsigned int l);
unsigned int
TToomCookMul_space (unsigned int n, unsigned int m, unsigned int l);
unsigned int muls_tkara (unsigned int n);
unsigned int muls_tgen (unsigned int n);
void list_sub_safe (listz_t ret, listz_t a, listz_t b,
                        unsigned int sizea, unsigned int sizeb,
                        unsigned int needed);
void list_add_safe (listz_t ret, listz_t a, listz_t b,
                        unsigned int sizea, unsigned int sizeb,
                        unsigned int needed);

/* memory.c */
#ifdef MEMORY_DEBUG
void __gmp_default_free (void *, size_t);
void *__gmp_default_allocate (size_t);
void *__gmp_default_reallocate (void *, size_t, size_t);
void tests_memory_start (void);
void tests_memory_end   (void);
void tests_memory_reset (void);
void tests_free (void *, size_t);
void tests_memory_status (void);
#define FREE(ptr,size) tests_free(ptr,size)
#else
#define FREE(ptr,size) free(ptr)
#endif

#if defined (__cplusplus)
}
#endif  

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
