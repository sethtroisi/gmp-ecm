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
  Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
  02111-1307, USA.
*/

#ifndef _ECM_ECM_H
#define _ECM_ECM_H 1

#include "config.h"

#include <assert.h>
#define ASSERT_ALWAYS(expr)   assert (expr)
#ifdef WANT_ASSERT
#define ASSERT(expr)   assert (expr)
#else
#define ASSERT(expr)   do {} while (0)
#endif

#include "ecm.h"

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

#define MAX_NUMBER_PRINT_LEN 1000

#define NTT_SIZE_THRESHOLD 30

/* auxlib.c */
int  test_verbose (int);
void set_verbose (int);
int  inc_verbose ();

/* Return codes */
/* Bit coded values: 1: error (for example out of memory)
   2: proper factor found, 4: factor is prime, 
   8: cofactor is prime or 1 */
#define ECM_EXIT_ERROR 1
#define ECM_COMP_FAC_COMP_COFAC 2
#define ECM_PRIME_FAC_COMP_COFAC (2+4)
#define ECM_INPUT_NUMBER_FOUND 8
#define ECM_COMP_FAC_PRIME_COFAC (2+8)
#define ECM_PRIME_FAC_PRIME_COFAC (2+4+8)

/* getprime.c */
double getprime ();
void getprime_clear ();
void getprime_seek (double);
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
void tests_memory_set_location (char *, unsigned int);
#endif

/* trial.c */
int trial_factor (mpcandi_t *n, double maxfact, int deep);

/* resume.c */
int  read_resumefile_line (int *, mpz_t, mpcandi_t *, mpz_t, mpz_t, mpz_t, double *,
                           char *, char *, char *, char *, FILE *);
int write_resumefile_line (char *, int, double, mpz_t, mpz_t, mpz_t, mpcandi_t *, 
                           mpz_t, const char *);

/* main.c */
int read_number (mpcandi_t *n, FILE *, int primetest);
int probab_prime_p (mpz_t, int);
int kbnc_z (double *k, unsigned long *b, unsigned long *n, signed long *c,
            mpz_t z);
int kbnc_str (double *k, unsigned long *b, unsigned long *n, signed long *c,
              char *z, mpz_t num);

/* batch.c */
void compute_s (mpz_t, unsigned long);

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

/* maximal stage 1 bound = 2^53 - 1, the next prime being 2^53 + 5 */
#define MAX_B1 9007199254740991.0

/* The checksum for savefile is the product of all mandatory fields, modulo
   the greatest prime below 2^32 */
#define CHKSUMMOD 4294967291U

#ifdef MEMORY_DEBUG
#define FREE(ptr,size) tests_free(ptr,size)
#define MEMORY_TAG tests_memory_set_location(__FILE__,__LINE__)
#define MEMORY_UNTAG tests_memory_set_location("",0)
#define MPZ_INIT(x)    {MEMORY_TAG;mpz_init(x);MEMORY_UNTAG;}
#define MPZ_INIT2(x,n) {MEMORY_TAG;mpz_init2(x,n);MEMORY_UNTAG;}
#else
#define FREE(ptr,size) free(ptr)
#define MEMORY_TAG do{}while(0)
#define MEMORY_UNTAG do{}while(0)
#define MPZ_INIT(x) mpz_init(x)
#define MPZ_INIT2(x,n) mpz_init2(x,n)
#endif

#define ABS(x) ((x) >= 0 ? (x) : -(x))

/* could go in auxi.c as a function */
#ifdef HAVE_SETPRIORITY
# include <sys/time.h>
# ifdef HAVE_SYS_RESOURCE_H
#  include <sys/resource.h>
# endif
# define NICE10 setpriority (PRIO_PROCESS, 0, 10)
# define NICE20 setpriority (PRIO_PROCESS, 0, 20)

#elif defined(HAVE_NICE)
# ifdef HAVE_UNISTD_H
#  include <unistd.h>
# endif
# define NICE10 nice (10)
# define NICE20 nice (20)

#elif defined(HAVE_WINDOWS_H)
# include <windows.h>
# define NICE10 do { \
   SetPriorityClass (GetCurrentProcess (), BELOW_NORMAL_PRIORITY_CLASS); \
   SetThreadPriority (GetCurrentThread (), THREAD_PRIORITY_BELOW_NORMAL); \
   } while (0)
# define NICE20 do { \
   SetPriorityClass (GetCurrentProcess (), IDLE_PRIORITY_CLASS); \
   SetThreadPriority (GetCurrentThread (), THREAD_PRIORITY_IDLE); \
   } while (0)

#else
# warning "Can't find a way to change priority"
# define NICE10 do {} while (0)
# define NICE20 do {} while (0)
#endif


#endif /* _ECM_ECM_H */
