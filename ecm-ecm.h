/* ecm-ecm.h - private header file for GMP-ECM.

  Copyright 2001, 2002, 2003, 2004, 2005 Paul Zimmermann and Alexander Kruppa.

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

#if WANT_ASSERT
#include <assert.h>
#define ASSERT(expr)   assert (expr)
#else
#define ASSERT(expr)   do {} while (0)
#endif

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
  mpz_t n;		/* the cofactor candidate currently being used to find factors from */
  unsigned ndigits;	/* the number of digits (decimal) in n */
  unsigned nexprlen;	/* strlen of expression, 0 if there is NO expression */
  int isPrp;		/* usually 0, but turns 1 if factor found, and the cofactor is PRP, 
			   OR if the original candidate was PRP and the user asked to prp check */
} mpcandi_t;

typedef struct
{
  int  Valid;           /* Is ONLY set to 1 if there is a proper -go <integer> switch.  Otherwise is 0
                           and if 0, then PM1, PP1 and ECM all ignore it */
  char *cpOrigExpr;	/* if non-NULL, then this is a "simpler" expression than the 
			   decimal output of n */
  mpcandi_t Candi;      /* The value after expression checked */
  int containsN;        /* 0 for simple number or expression.  1 if the expression "contains" N as
                           that expression will have to be built for each candidate */
} mpgocandi_t;

/* auxi.c */
unsigned int nb_digits  (const mpz_t);
unsigned int get_random_ui (void);

#define OUTPUT_ALWAYS 0
#define OUTPUT_NORMAL 1
#define OUTPUT_VERBOSE 2
#define OUTPUT_DEVVERBOSE 3
#define OUTPUT_TRACE 4
#define OUTPUT_ERROR -1

/* auxlib.c */
int  test_verbose (int);
void set_verbose (int);
int  inc_verbose ();

/* different methods implemented */
#define ECM_ECM 0
#define ECM_PM1 1
#define ECM_PP1 2

/* getprime2.c */
double getprime (double);
#define WANT_FREE_PRIME_TABLE(p) (p < 0.0)
#define FREE_PRIME_TABLE -1.0

/* b1_ainc.c */
double calc_B1_AutoIncrement(double cur_B1, double incB1val, int calcInc);

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
#endif

/* smartprp.c */
int smart_probab_prime_p (mpz_t const n, int c);

/* trial.c */
int trial_factor (mpcandi_t *n, double maxfact, int deep);

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
void usage (void);

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
/* candi.c   Group Order candidate functions.  */
void mpgocandi_t_init(mpgocandi_t *go);
void mpgocandi_t_free(mpgocandi_t *go);
int  mpgocandi_fixup_with_N(mpgocandi_t *go, mpcandi_t *n);

/* random2.c */
void pp1_random_seed  (mpz_t, mpz_t, gmp_randstate_t);
void pm1_random_seed  (mpz_t, mpz_t, gmp_randstate_t);

/* default number of probable prime tests */
#define PROBAB_PRIME_TESTS 1

/* Options for using an external PRPer rather than internal GMP */
extern char *externalprp; /* call out to external program  */  
extern char *externallog; /* where to find output */
extern char *externalinputprpfile; /* Used to place the n value (a temp file). Is deleted after system */
extern char *externalisprp; /* what to match output against */
extern char *externaliscomposite; /* ditto */
extern int externalprplen; /* length where we go external */
extern int externalprpval; /* exit value meaning it's prp, -1 is impossible */

/* maximal stage 1 bound = 2^53 + 4, the next prime being 2^53 + 5 */
#define MAX_B1 9007199254740996.0

/* The checksum for savefile is the product of all mandatory fields, modulo
   the greatest prime below 2^32 */
#define CHKSUMMOD 4294967291U

#ifdef MEMORY_DEBUG
#define FREE(ptr,size) tests_free(ptr,size)
#else
#define FREE(ptr,size) free(ptr)
#endif

#define ABS(x) ((x) >= 0 ? (x) : -(x))


