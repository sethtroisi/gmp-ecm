/* ecm-impl.h - header file for libecm
 
Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011,
2012 Paul Zimmermann, Alexander Kruppa and Cyril Bouvier.
 
This file is part of the ECM Library.

The ECM Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The ECM Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the ECM Library; see the file COPYING.LIB.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#ifndef _ECM_IMPL_H
#define _ECM_IMPL_H 1

#include "config.h"
#include "basicdefs.h"
#include "ecm.h"
#include "sp.h"

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h> /* needed for size_t */
#endif

#if HAVE_STDINT_H
#include <stdint.h>
/* needed for int64_t and uint64_t */
/* or configure will define these for us if possible */
#endif

/* We do not use torsion.[ch] so far since they are not tested enough. */
#define HAVE_TORSION
/* We do not use addlaws.[ch] so far since they are not tested enough. */
#define HAVE_ADDLAWS

#include "ecm_int.h"

#ifndef TUNE
#include "ecm-params.h"
#else
extern size_t MPZMOD_THRESHOLD;
extern size_t REDC_THRESHOLD;
#define TUNE_MULREDC_TABLE {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
#define TUNE_SQRREDC_TABLE {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
#define LIST_MUL_TABLE {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
#endif
extern size_t mpn_mul_lo_threshold[];

#define TUNE_LIST_MUL_N_MAX_SIZE 32

#include <stdio.h> /* needed for "FILE *" */
#include <limits.h>

#if  defined (__STDC__)                                 \
  || defined (__cplusplus)                              \
  || defined (_AIX)                                     \
  || defined (__DECC)                                   \
  || (defined (__mips) && defined (_SYSTYPE_SVR4))      \
  || defined (_MSC_VER)                                 \
  || defined (_WIN32)
#define __ECM_HAVE_TOKEN_PASTE  1
#else
#define __ECM_HAVE_TOKEN_PASTE  0
#endif

#ifndef __ECM
#if __ECM_HAVE_TOKEN_PASTE
#define __ECM(x) __ecm_##x
#else
#define __ECM(x) __ecm_/**/x
#endif
#endif

#define ECM_STDOUT __ecm_stdout
#define ECM_STDERR __ecm_stderr
extern FILE *ECM_STDOUT, *ECM_STDERR;

/* #define TIMING_CRT */

/* default B2 choice: pow (B1 * METHOD_COST / 6.0, DEFAULT_B2_EXPONENT) */
#define DEFAULT_B2_EXPONENT 1.43
#define PM1_COST 1.0 / 6.0
#define PP1_COST 2.0 / 6.0
#define ECM_COST 11.0 / 6.0
/* For new P-/+1 stage 2: */
#define PM1FS2_DEFAULT_B2_EXPONENT 1.7
#define PM1FS2_COST 1.0 / 4.0
#define PP1FS2_COST 1.0 / 4.0

/* define top-level multiplication */
#define KARA 2
#define TOOM3 3
#define TOOM4 4
#define KS 5
#define NTT 6

/* maximal limb size of assembly mulredc */
#define MULREDC_ASSEMBLY_MAX 20

#include "sp.h"

#include <assert.h>
#define ASSERT_ALWAYS(expr)   assert (expr)
#ifdef WANT_ASSERT
#define ASSERT(expr)   assert (expr)
#else
#define ASSERT(expr)   do {} while (0)
#endif

/* thresholds */
#define MPN_MUL_LO_THRESHOLD 32

/* base2mod is used when size(2^n+/-1) <= BASE2_THRESHOLD * size(cofactor) */
#define BASE2_THRESHOLD 1.4

/* default number of probable prime tests */
#define PROBAB_PRIME_TESTS 1

/* threshold for median product */
#define KS_TMUL_THRESHOLD 8e5

#define ABS(x) ((x) >= 0 ? (x) : -(x))

/* getprime */
#define WANT_FREE_PRIME_TABLE(p) (p < 0.0)
#define FREE_PRIME_TABLE -1.0

/* 2^n+-1 with n < MOD_MINBASE2 cannot use base-2 reduction */
#define MOD_MINBASE2 16

/* Various logging levels */
/* OUTPUT_ALWAYS means print always, regardless of verbose value */
#define OUTPUT_ALWAYS 0
/* OUTPUT_NORMAL means print during normal program execution */
#define OUTPUT_NORMAL 1
/* OUTPUT_VERBOSE means print if the user requested more verbosity */
#define OUTPUT_VERBOSE 2
/* OUTPUT_RESVERBOSE is for printing residues (after stage 1 etc) */
#define OUTPUT_RESVERBOSE 3
/* OUTPUT_DEVVERBOSE is for printing internal parameters (for developers) */
#define OUTPUT_DEVVERBOSE 4
/* OUTPUT_TRACE is for printing trace data, produces lots of output */
#define OUTPUT_TRACE 5
/* OUTPUT_ERROR is for printing error messages */
#define OUTPUT_ERROR -1

/* Interval length for writing checkpoints in stage 1, in milliseconds */
#define CHKPNT_PERIOD 600000

/* Does the parametrization imply batch mode ? */
#define IS_BATCH_MODE(p) ( p == ECM_PARAM_BATCH_SQUARE || \
                           p == ECM_PARAM_BATCH_2 || \
                           p == ECM_PARAM_BATCH_32BITS_D )

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
  /* for CM curves */
  int disc;
  mpres_t sq[10];
} __curve_struct;
typedef __curve_struct curve;

typedef struct
{
  unsigned long d1;
  unsigned long d2;
  mpz_t i0;
  int S;
} __root_params_t;
typedef __root_params_t root_params_t;

typedef struct
{
  unsigned long P, s_1, s_2, l;
  mpz_t m_1;
  const char *file_stem;
} __faststage2_param_t;
typedef __faststage2_param_t faststage2_param_t;

typedef struct
{
  unsigned int size_fd; /* How many entries .fd has, always nr * (S+1) */
  unsigned int nr;     /* How many separate progressions there are */
  unsigned int next;   /* From which progression to take the next root */
  unsigned int S;      /* Degree of the polynomials */
  unsigned int dsieve; /* Values not coprime to dsieve are skipped */
  unsigned int rsieve; /* Which residue mod dsieve current .next belongs to */
  int dickson_a;       /* Parameter for Dickson polynomials */
} progression_params_t;

typedef struct
{
  progression_params_t params;
  point *fd;
  unsigned int size_T; /* How many entries T has */
  mpres_t *T;          /* For temp values. FIXME: should go! */
  curve *X;            /* The curve the points are on */
} ecm_roots_state_t;


typedef struct
{
  progression_params_t params;
  mpres_t *fd;
  int invtrick;
} pm1_roots_state_t;

typedef struct
{
  progression_params_t params;
  point *fd;           /* for S != 1 */
  mpres_t tmp[4];      /* for S=1 */
} pp1_roots_state_t;

typedef struct
{
  int alloc;
  int degree;
  listz_t coeff;
} __polyz_struct;
typedef __polyz_struct polyz_t[1];

typedef struct 
{
  int repr;           /* ECM_MOD_MPZ: plain modulus, possibly normalized
                         ECM_MOD_BASE2: base 2 number
                         ECM_MOD_MODMULN: MODMULN
                         ECM_MOD_REDC: REDC representation */
  int bits;           /* in case of a base 2 number, 2^k[+-]1, bits = [+-]k
                         in case of MODMULN or REDC representation, nr. of 
                         bits b so that 2^b > orig_modulus and 
                         GMP_NUMB_BITS | b */
  int Fermat;         /* If repr = 1 (base 2 number): If modulus is 2^(2^m)+1, 
                         i.e. bits = 2^m, then Fermat = 2^m, 0 otherwise.
                         If repr != 1, undefined */
  mp_limb_t *Nprim;   /* For MODMULN */
  mpz_t orig_modulus; /* The original modulus N */
  mpz_t aux_modulus;  /* Used only for MPZ and REDC:
			 - the auxiliary modulus value (i.e. normalized 
                           modulus, or -1/N (mod 2^bits) for REDC,
                         - B^(n + ceil(n/2)) mod N for MPZ,
  			   where B = 2^GMP_NUMB_BITS */
  mpz_t multiple;     /* The smallest multiple of N that is larger or
			 equal to 2^bits for REDC/MODMULN */
  mpz_t R2, R3;       /* For MODMULN and REDC, R^2 and R^3 (mod orig_modulus), 
                         where R = 2^bits. */
  mpz_t temp1, temp2; /* Temp values used during multiplication etc. */
} __mpmod_struct;
typedef __mpmod_struct mpmod_t[1];

#if defined (__cplusplus)
extern "C" {
#endif  

/* getprime.c */
#define getprime __ECM(getprime)
double   getprime       ();
#define getprime_clear __ECM(getprime_clear)
void     getprime_clear ();
#define getprime_seek __ECM(getprime_seek)
void getprime_seek (double);

/* pm1fs2.c */
#define pm1fs2_memory_use __ECM(pm1fs2_ntt_memory_use)
size_t  pm1fs2_memory_use (const unsigned long, const mpz_t, const int);
#define pm1fs2_maxlen __ECM(pm1fs2_maxlen)
unsigned long pm1fs2_maxlen (const size_t, const mpz_t, const int);
#define pp1fs2_memory_use __ECM(pp1fs2_ntt_memory_use)
size_t  pp1fs2_memory_use (const unsigned long, const mpz_t, const int, 
                           const int);
#define pp1fs2_maxlen __ECM(pp1fs2_maxlen)
unsigned long pp1fs2_maxlen (const size_t, const mpz_t, const int, const int);
#define choose_P __ECM(choose_P)
long    choose_P (const mpz_t, const mpz_t, const unsigned long,
                  const unsigned long, faststage2_param_t *, mpz_t, mpz_t,
                  const int, const int);
#define pm1fs2 __ECM(pm1fs2)
int	pm1fs2 (mpz_t, const mpres_t, mpmod_t, const faststage2_param_t *);
#define pm1fs2_ntt __ECM(pm1fs2_ntt)
int	pm1fs2_ntt (mpz_t, const mpres_t, mpmod_t, const faststage2_param_t *);
#define pp1fs2 __ECM(pp1fs2)
int     pp1fs2 (mpz_t, const mpres_t, mpmod_t, const faststage2_param_t *);
#define pp1fs2_ntt __ECM(pp1fs2_ntt)
int     pp1fs2_ntt (mpz_t, const mpres_t, mpmod_t, const faststage2_param_t *,
                    const int);

/* bestd.c */
#define bestD __ECM(bestD)
int     bestD (root_params_t *, unsigned long *, unsigned long *, mpz_t, 
               mpz_t, int, int, double, int, mpmod_t);

/* ecm.c */
#define choose_S __ECM(choose_S)
int  choose_S (mpz_t);
#define add3 __ECM(add3)
void add3 (mpres_t, mpres_t, mpres_t, mpres_t, mpres_t, mpres_t, mpres_t, 
           mpres_t, mpmod_t, mpres_t, mpres_t, mpres_t);
#define duplicate __ECM(duplicate)
void duplicate (mpres_t, mpres_t, mpres_t, mpres_t, mpmod_t, mpres_t, mpres_t,
                mpres_t, mpres_t);

#define ecm_mul __ECM(ecm_mul)
void ecm_mul (mpres_t, mpres_t, mpz_t, mpmod_t, mpres_t);
#define print_B1_B2_poly __ECM(print_B1_B2_poly)
void print_B1_B2_poly (int, int, double, double, mpz_t, mpz_t, mpz_t, int S,  
                       mpz_t, int, int, mpz_t, int, unsigned int);
#define set_stage_2_params __ECM(set_stage_2_params)
int set_stage_2_params (mpz_t, mpz_t, mpz_t, mpz_t, root_params_t *,
                        double, unsigned long *, const int, int, int *,
                        unsigned long *, char *, double, int, mpmod_t);
#define print_expcurves __ECM(print_expcurves)
void print_expcurves (double, const mpz_t, unsigned long, unsigned long, int, 
                      int);
#define print_expwork __ECM(print_expwork)
void print_expwork (double, const mpz_t, unsigned long, unsigned long, int, 
                    int, unsigned int);
#define print_exptime __ECM(print_exptime)
void print_exptime (double, const mpz_t, unsigned long, unsigned long, int, 
                    double, int);
#define montgomery_to_weierstrass __ECM(montgomery_to_weierstrass)
int montgomery_to_weierstrass (mpz_t, mpres_t, mpres_t, mpres_t, mpmod_t);

/* ecm2.c */
#define ecm_rootsF __ECM(ecm_rootsF)
int     ecm_rootsF       (mpz_t, listz_t, root_params_t *, unsigned long, 
                          curve *, mpmod_t);
#define ecm_rootsG_init __ECM(ecm_rootsG_init)
ecm_roots_state_t* ecm_rootsG_init (mpz_t, curve *, root_params_t *, 
                                    unsigned long, unsigned long, mpmod_t);
#define ecm_rootsG __ECM(ecm_rootsG)
int     ecm_rootsG       (mpz_t, listz_t, unsigned long, ecm_roots_state_t *, 
                          mpmod_t);
#define ecm_rootsG_clear __ECM(ecm_rootsG_clear)
void    ecm_rootsG_clear (ecm_roots_state_t *, mpmod_t);

/* lucas.c */
#define pp1_mul_prac __ECM(pp1_mul_prac)
void  pp1_mul_prac     (mpres_t, ecm_uint, mpmod_t, mpres_t, mpres_t,
                        mpres_t, mpres_t, mpres_t);

/* stage2.c */
#define stage2 __ECM(stage2)
int     stage2 (mpz_t, void *, mpmod_t, unsigned long, unsigned long,
                root_params_t *, int, char *, unsigned int, int (*)(void));
#define init_progression_coeffs __ECM(init_progression_coeffs)
listz_t init_progression_coeffs (mpz_t, const unsigned long, const unsigned long, 
				 const unsigned int, const unsigned int, 
				 const unsigned int, const int);
#define init_roots_params __ECM(init_roots_params)
void init_roots_params  (progression_params_t *, const int, 
			 const unsigned long, const unsigned long, 
			 const double);
#define memory_use __ECM(memory_use)
double memory_use (unsigned long, unsigned int, unsigned int, mpmod_t);

/* listz.c */
#define list_mul_mem __ECM(list_mul_mem)
int          list_mul_mem (unsigned int);
#define init_list __ECM(init_list)
listz_t      init_list  (unsigned int);
#define init_list2 __ECM(init_list2)
listz_t      init_list2  (unsigned int, unsigned int);
#define clear_list __ECM(clear_list)
void         clear_list (listz_t, unsigned int);
#define list_inp_raw __ECM(list_inp_raw)
int          list_inp_raw (listz_t, FILE *, unsigned int);
#define list_out_raw __ECM(list_out_raw)
int          list_out_raw (FILE *, listz_t, unsigned int);
#define print_list __ECM(print_list)
void         print_list (listz_t, unsigned int);
#define list_set __ECM(list_set)
void         list_set   (listz_t, listz_t, unsigned int);
#define list_revert __ECM(list_revert)
void         list_revert (listz_t, unsigned int);
#define list_swap __ECM(list_swap)
void         list_swap  (listz_t, listz_t, unsigned int);
#define list_neg __ECM(list_neg)
void         list_neg (listz_t, listz_t, unsigned int, mpz_t);
#define list_mod __ECM(list_mod)
void         list_mod   (listz_t, listz_t, unsigned int, mpz_t);
#define list_add __ECM(list_add)
void         list_add   (listz_t, listz_t, listz_t, unsigned int);
#define list_sub __ECM(list_sub)
void         list_sub   (listz_t, listz_t, listz_t, unsigned int);
#define list_mul_z __ECM(list_mul_z)
void         list_mul_z (listz_t, listz_t, mpz_t, unsigned int, mpz_t);
#define list_mulup __ECM(list_mulup)
void          list_mulup (listz_t, unsigned int, mpz_t, mpz_t);
#define list_zero __ECM(list_zero)
void         list_zero  (listz_t, unsigned int);
#define list_mul __ECM(list_mul)
void         list_mul (listz_t, listz_t, unsigned int, listz_t,
    unsigned int, int, listz_t);
#define list_mul_high __ECM(list_mul_high)
void      list_mul_high (listz_t, listz_t, listz_t, unsigned int);
#define list_mulmod __ECM(list_mulmod)
void        list_mulmod (listz_t, listz_t, listz_t, listz_t, unsigned int,
                         listz_t, mpz_t);
#define PolyFromRoots __ECM(PolyFromRoots)
void      PolyFromRoots (listz_t, listz_t, unsigned int, listz_t, mpz_t);
#define PolyFromRoots_Tree __ECM(PolyFromRoots_Tree)
int       PolyFromRoots_Tree (listz_t, listz_t, unsigned int, listz_t, int, 
                         mpz_t, listz_t*, FILE*, unsigned int);

#define ntt_PolyFromRoots __ECM(ntt_PolyFromRoots)
void	  ntt_PolyFromRoots (mpzv_t, mpzv_t, spv_size_t, mpzv_t, mpzspm_t);
#define ntt_PolyFromRoots_Tree __ECM(ntt_PolyFromRoots_Tree)
int       ntt_PolyFromRoots_Tree (mpzv_t, mpzv_t, spv_size_t, mpzv_t,
                         int, mpzspm_t, mpzv_t *, FILE *);
#define ntt_polyevalT __ECM(ntt_polyevalT)
int  ntt_polyevalT (mpzv_t, spv_size_t, mpzv_t *, mpzv_t, mpzspv_t,
		mpzspm_t, char *);
#define ntt_mul __ECM(ntt_mul)
void  ntt_mul (mpzv_t, mpzv_t, mpzv_t, spv_size_t, mpzv_t, int, mpzspm_t);
#define ntt_PrerevertDivision __ECM(ntt_PrerevertDivision)
void  ntt_PrerevertDivision (mpzv_t, mpzv_t, mpzv_t, mpzspv_t, mpzspv_t,
		spv_size_t, mpzv_t, mpzspm_t);
#define ntt_PolyInvert __ECM(ntt_PolyInvert)
void	     ntt_PolyInvert (mpzv_t, mpzv_t, spv_size_t, mpzv_t, mpzspm_t);

#define PrerevertDivision __ECM(PrerevertDivision)
int   PrerevertDivision (listz_t, listz_t, listz_t, unsigned int, listz_t,
			 mpz_t);
#define PolyInvert __ECM(PolyInvert)
void         PolyInvert (listz_t, listz_t, unsigned int, listz_t, mpz_t);

#define RecursiveDivision __ECM(RecursiveDivision)
void  RecursiveDivision (listz_t, listz_t, listz_t, unsigned int,
                         listz_t, mpz_t, int);

/* polyeval.c */
#define polyeval __ECM(polyeval)
void polyeval (listz_t, unsigned int, listz_t*, listz_t, mpz_t, unsigned int);
#define polyeval_tellegen __ECM(polyeval_tellegen)
int polyeval_tellegen (listz_t, unsigned int, listz_t*, listz_t,
		       unsigned int, listz_t, mpz_t, char *);
#define TUpTree __ECM(TUpTree)
void TUpTree (listz_t, listz_t *, unsigned int, listz_t, int, unsigned int,
		mpz_t, FILE *);

/* ks-multiply.c */
#define list_mul_n_basecase __ECM(list_mul_n_basecase)
void list_mul_n_basecase (listz_t, listz_t, listz_t, unsigned int);
#define list_mul_n_karatsuba __ECM(list_mul_n_karatsuba)
void list_mul_n_karatsuba (listz_t, listz_t, listz_t, unsigned int);
#define list_mul_n_KS1 __ECM(list_mul_n_KS1)
void list_mul_n_KS1 (listz_t, listz_t, listz_t, unsigned int);
#define list_mul_n_KS2 __ECM(list_mul_n_KS2)
void list_mul_n_KS2 (listz_t, listz_t, listz_t, unsigned int);
#define list_mult_n __ECM(list_mult_n)
void list_mult_n (listz_t, listz_t, listz_t, unsigned int);
#define TMulKS __ECM(TMulKS)
int TMulKS     (listz_t, unsigned int, listz_t, unsigned int, listz_t,
                unsigned int, mpz_t, int);
#define ks_wrapmul_m __ECM(ks_wrapmul_m)
unsigned int ks_wrapmul_m (unsigned int, unsigned int, mpz_t);
#define ks_wrapmul __ECM(ks_wrapmul)
unsigned int ks_wrapmul (listz_t, unsigned int, listz_t, unsigned int,
                         listz_t, unsigned int, mpz_t);

/* mpmod.c */
/* Define MPRESN_NO_ADJUSTMENT if mpresn_add, mpresn_sub and mpresn_addsub
   should perform no adjustment step. This yields constraints on N. */
/* #define MPRESN_NO_ADJUSTMENT */
#define isbase2 __ECM(isbase2)
int isbase2 (const mpz_t, const double);
#define mpmod_init __ECM(mpmod_init)
int mpmod_init (mpmod_t, const mpz_t, int);
#define mpmod_init_MPZ __ECM(mpmod_init_MPZ)
void mpmod_init_MPZ (mpmod_t, const mpz_t);
#define mpmod_init_BASE2 __ECM(mpmod_init_BASE2)
int mpmod_init_BASE2 (mpmod_t, const int, const mpz_t);
#define mpmod_init_MODMULN __ECM(mpmod_init_MODMULN)
void mpmod_init_MODMULN (mpmod_t, const mpz_t);
#define mpmod_init_REDC __ECM(mpmod_init_REDC)
void mpmod_init_REDC (mpmod_t, const mpz_t);
#define mpmod_clear __ECM(mpmod_clear)
void mpmod_clear (mpmod_t);
#define mpmod_init_set __ECM(mpmod_init_set)
void mpmod_init_set (mpmod_t, const mpmod_t);
#define mpmod_pausegw __ECM(mpmod_pausegw)
void mpmod_pausegw (const mpmod_t modulus);
#define mpmod_contgw __ECM(mpmod_contgw)
void mpmod_contgw (const mpmod_t modulus);
#define mpres_equal __ECM(mpres_equal)
int mpres_equal (const mpres_t, const mpres_t, mpmod_t);
#define mpres_pow __ECM(mpres_pow)
void mpres_pow (mpres_t, const mpres_t, const mpz_t, mpmod_t);
#define mpres_ui_pow __ECM(mpres_ui_pow)
void mpres_ui_pow (mpres_t, const unsigned long, const mpres_t, mpmod_t);
#define mpres_mul __ECM(mpres_mul)
void mpres_mul (mpres_t, const mpres_t, const mpres_t, mpmod_t) ATTRIBUTE_HOT;
#define mpres_sqr __ECM(mpres_sqr)
void mpres_sqr (mpres_t, const mpres_t, mpmod_t) ATTRIBUTE_HOT;
#define mpres_mul_z_to_z __ECM(mpres_mul_z_to_z)
void mpres_mul_z_to_z (mpz_t, const mpres_t, const mpz_t, mpmod_t);
#define mpres_set_z_for_gcd __ECM(mpres_set_z_for_gcd)
void mpres_set_z_for_gcd (mpres_t, const mpz_t, mpmod_t);
#define mpres_set_z_for_gcd_fix __ECM(mpres_set_z_for_gcd_fix)
void mpres_set_z_for_gcd_fix (mpres_t, const mpres_t, const mpz_t, mpmod_t);
#define mpres_div_2exp __ECM(mpres_div_2exp)
void mpres_div_2exp (mpres_t, const mpres_t, const unsigned int, mpmod_t);
#define mpres_add_ui __ECM(mpres_add_ui)
void mpres_add_ui (mpres_t, const mpres_t, const unsigned long, mpmod_t);
#define mpres_add __ECM(mpres_add)
void mpres_add (mpres_t, const mpres_t, const mpres_t, mpmod_t) ATTRIBUTE_HOT;
#define mpres_sub_ui __ECM(mpres_sub_ui)
void mpres_sub_ui (mpres_t, const mpres_t, const unsigned long, mpmod_t);
#define mpres_ui_sub __ECM(mpres_ui_sub)
void mpres_ui_sub (mpres_t, const unsigned long, const mpres_t, mpmod_t);
#define mpres_sub __ECM(mpres_sub)
void mpres_sub (mpres_t, const mpres_t, const mpres_t, mpmod_t) ATTRIBUTE_HOT;
#define mpres_set_z __ECM(mpres_set_z)
void mpres_set_z (mpres_t, const mpz_t, mpmod_t);
#define mpres_get_z __ECM(mpres_get_z)
void mpres_get_z (mpz_t, const mpres_t, mpmod_t);
#define mpres_set_ui __ECM(mpres_set_ui)
void mpres_set_ui (mpres_t, const unsigned long, mpmod_t);
#define mpres_set_si __ECM(mpres_set_si)
void mpres_set_si (mpres_t, const long, mpmod_t);
#define mpres_init __ECM(mpres_init)
void mpres_init (mpres_t, const mpmod_t);
#define mpres_clear __ECM(mpres_clear)
void mpres_clear (mpres_t, const mpmod_t);
#define mpres_realloc __ECM(mpres_realloc)
void mpres_realloc (mpres_t, const mpmod_t);
#define mpres_mul_ui __ECM(mpres_mul_ui)
void mpres_mul_ui (mpres_t, const mpres_t, const unsigned long, mpmod_t);
#define mpres_mul_2exp __ECM(mpres_mul_2exp)
void mpres_mul_2exp (mpres_t, const mpres_t, const unsigned long, mpmod_t);
#define mpres_muldivbysomething_si __ECM(mpres_muldivbysomething_si)
void mpres_muldivbysomething_si (mpres_t, const mpres_t, const long, mpmod_t);
#define mpres_neg __ECM(mpres_neg)
void mpres_neg (mpres_t, const mpres_t, mpmod_t);
#define mpres_invert __ECM(mpres_invert)
int  mpres_invert (mpres_t, const mpres_t, mpmod_t);
#define mpres_gcd __ECM(mpres_gcd)
void mpres_gcd (mpz_t, const mpres_t, const mpmod_t);
#define mpres_out_str __ECM(mpres_out_str)
void mpres_out_str (FILE *, const unsigned int, const mpres_t, mpmod_t);
#define mpres_is_zero __ECM(mpres_is_zero)
int  mpres_is_zero (const mpres_t, mpmod_t);
#define mpres_set(a,b,n) mpz_set (a, b)
#define mpres_swap(a,b,n) mpz_swap (a, b)
#define mpresn_mul __ECM(mpresn_mul)
void mpresn_mul (mpres_t, const mpres_t, const mpres_t, mpmod_t);
#define mpresn_addsub __ECM(mpresn_addsub)
void mpresn_addsub (mpres_t, mpres_t, const mpres_t, const mpres_t, mpmod_t);
#define mpresn_pad __ECM(mpresn_pad)
void mpresn_pad (mpres_t R, mpmod_t N);
#define mpresn_unpad __ECM(mpresn_unpad)
void mpresn_unpad (mpres_t R);
#define mpresn_sqr __ECM(mpresn_sqr)
void mpresn_sqr (mpres_t, const mpres_t, mpmod_t);
#define mpresn_add __ECM(mpresn_add)
void mpresn_add (mpres_t, const mpres_t, const mpres_t, mpmod_t);
#define mpresn_sub __ECM(mpresn_sub)
void mpresn_sub (mpres_t, const mpres_t, const mpres_t, mpmod_t);
#define mpresn_mul_1 __ECM(mpresn_mul_ui)
void mpresn_mul_1 (mpres_t, const mpres_t, const mp_limb_t, mpmod_t);

/* mul_lo.c */
#define ecm_mul_lo_n __ECM(ecm_mul_lo_n)
void ecm_mul_lo_n (mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);
#define ecm_mul_lo_basecase __ECM(ecm_mul_lo_basecase)
void ecm_mul_lo_basecase (mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);
	
/* median.c */
#define TMulGen __ECM(TMulGen)
int
TMulGen (listz_t, unsigned int, listz_t, unsigned int, listz_t, 
         unsigned int, listz_t, mpz_t);
#define TMulGen_space __ECM(TMulGen_space)
unsigned int TMulGen_space (unsigned int, unsigned int, unsigned int);

/* schoen_strass.c */
#define DEFAULT 0
#define MONIC 1
#define NOPAD 2
#define F_mul __ECM(F_mul)
unsigned int F_mul (mpz_t *, mpz_t *, mpz_t *, unsigned int, int,
                    unsigned int, mpz_t *);
#define F_mul_trans __ECM(F_mul_trans)
unsigned int F_mul_trans (mpz_t *, mpz_t *, mpz_t *, unsigned int,
                          unsigned int, unsigned int, mpz_t *);
#define F_clear __ECM(F_clear)
void F_clear ();

/* rho.c */
#define rhoinit __ECM(rhoinit)
void   rhoinit (int, int);
#define ecmprob __ECM(ecmprob)
double ecmprob (double, double, double, double, int);
double pm1prob (double, double, double, double, int, const mpz_t);

/* pm1.c */
void print_prob (double, const mpz_t, unsigned long, unsigned long, int,
                 const mpz_t);

/* auxlib.c */
#define mpz_add_si __ECM(mpz_add_si)
void         mpz_add_si (mpz_t, mpz_t, long);
#define mpz_sub_si __ECM(mpz_sub_si)
void         mpz_sub_si (mpz_t, mpz_t, long);
#define mpz_divby3_1op __ECM(mpz_divby3_1op)
void         mpz_divby3_1op (mpz_t);
#define double_to_size __ECM(double_to_size)
size_t   double_to_size (double d);
#define cputime __ECM(cputime)
long         cputime    (void);
#define realtime __ECM(realtime)
long         realtime    (void);
#define elltime __ECM(elltime)
long         elltime    (long, long);
#define test_verbose __ECM(test_verbose)
int          test_verbose (int);
#define set_verbose __ECM(set_verbose)
void         set_verbose (int);
#define outputf __ECM(outputf)
int          outputf (int, const char *, ...);
#define writechkfile __ECM(writechkfile)
void writechkfile (char *, int, double, mpmod_t, mpres_t, mpres_t, mpres_t, mpres_t);
#define aux_fseek64 __ECM(aux_fseek64)
int aux_fseek64(FILE *, const int64_t, const int);
#define ecm_tstbit __ECM(ecm_tstbit)
int ecm_tstbit (mpz_srcptr, ecm_uint);

/* Due to GMP (6.x and prior) using long as input to mpz_tstbit, factors would be missed
   on computers with 32-bit longs in batch mode when using B1 > 2977044736UL.
   So, we need to use our own function when long is not 64-bits wide */
#if ULONG_MAX == 0xffffffffUL
#undef mpz_tstbit
#define mpz_tstbit ecm_tstbit
#endif

/* auxarith.c */
#define gcd __ECM(gcd)
unsigned long gcd (unsigned long, unsigned long);
#define eulerphi __ECM(eulerphi)
unsigned long eulerphi (unsigned long);
#define ceil_log2 __ECM(ceil_log2)
unsigned int  ceil_log2  (unsigned long);
#define find_factor __ECM(find_factor)
unsigned long find_factor (const unsigned long);

/* random.c */
#define init_randstate __ECM(init_randstate)
void init_randstate (gmp_randstate_t);
#define pp1_random_seed __ECM(pp1_random_seed)
void pp1_random_seed  (mpz_t, mpz_t, gmp_randstate_t);
#define pm1_random_seed __ECM(pm1_random_seed)
void pm1_random_seed  (mpz_t, mpz_t, gmp_randstate_t);
#define get_random_ul   __ECM(get_random_ul)
unsigned long get_random_ul (void);

/* Fgw.c */
#ifdef HAVE_GWNUM
int  gw_ecm_stage1 (mpz_t, curve *, mpmod_t, double, double *, mpz_t,
                    double, unsigned long, unsigned long, signed long);
#endif

/* batch.c */
#define compute_s  __ECM(compute_s )
void compute_s (mpz_t, ecm_uint, int *);
#define ecm_stage1_batch  __ECM(ecm_stage1_batch)
int ecm_stage1_batch (mpz_t, mpres_t, mpres_t, mpmod_t, double, double *, 
                                                                int,  mpz_t);

/* parametrizations.c */
#define get_curve_from_random_parameter __ECM(get_curve_from_random_parameter)
int get_curve_from_random_parameter (mpz_t, mpres_t, mpres_t, mpz_t, int, 
                                      mpmod_t, gmp_randstate_t);
#define get_curve_from_param0 __ECM(get_curve_from_param0)
int get_curve_from_param0 (mpz_t, mpres_t, mpres_t, mpz_t, mpmod_t);
#define get_curve_from_param1 __ECM(get_curve_from_param1)
int get_curve_from_param1 (mpres_t, mpres_t, mpz_t, mpmod_t);
#define get_curve_from_param2 __ECM(get_curve_from_param2)
int get_curve_from_param2 (mpz_t, mpres_t, mpres_t, mpz_t, mpmod_t);
#define get_curve_from_param3 __ECM(get_curve_from_param3)
int get_curve_from_param3 (mpres_t, mpres_t, mpz_t, mpmod_t);
#define get_default_param __ECM(get_default_param)
int get_default_param (int, double, int);

/* sets_long.c */
/* A set of long ints */
typedef struct {
  unsigned long card;
  long elem[1];
} set_long_t;

/* A set of sets of long ints */
typedef struct {
  unsigned long nr;
  set_long_t sets[1];
} sets_long_t;

#define quicksort_long __ECM(quicksort_long)
void          quicksort_long (long *, unsigned long);
#define sets_print __ECM(sets_print)
void          sets_print (const int, sets_long_t *);
#define sets_max __ECM(sets_max)
void          sets_max (mpz_t, const unsigned long);
#define sets_sumset __ECM(sets_sumset)
void          sets_sumset (set_long_t *, const sets_long_t *);
#define sets_sumset_minmax __ECM(sets_sumset_minmax)
void          sets_sumset_minmax (mpz_t, const sets_long_t *, const int);
#define sets_extract __ECM(sets_extract)
void          sets_extract (sets_long_t *, size_t *, sets_long_t *, 
                            const unsigned long);
#define sets_get_factored_sorted __ECM(sets_get_factored_sorted)
sets_long_t *  sets_get_factored_sorted (const unsigned long);

/* Return the size in bytes of a set of cardinality c */
#define set_sizeof __ECM(set_sizeof)
ATTRIBUTE_UNUSED
static size_t 
set_sizeof (const unsigned long c)
{
  return sizeof (long) + (size_t) c * sizeof (unsigned long);
}


/* Return pointer to the next set in "*sets" */
ATTRIBUTE_UNUSED
static set_long_t *
sets_nextset (const set_long_t *sets)
{
  return (set_long_t *) ((char *)sets + sizeof(unsigned long) + 
                         sets->card * sizeof(long));
}


#if defined (__cplusplus)
}
#endif

/* a <- b * c where a and b are mpz, c is a double, and t an auxiliary mpz */
/* Not sure how the preprocessor handles shifts by more than the integer 
   width on 32 bit machines, so do the shift by 53 in two pieces */
#if (((ULONG_MAX >> 27) >> 26) >= 1)
#define mpz_mul_d(a, b, c, t) \
   mpz_mul_ui (a, b, (unsigned long int) c);
#else
#define mpz_mul_d(a, b, c, t) \
   if (c < (double) ULONG_MAX) \
      mpz_mul_ui (a, b, (unsigned long int) c); \
   else { \
   mpz_set_d (t, c); \
   mpz_mul (a, b, t); }
#endif

#endif /* _ECM_IMPL_H */
