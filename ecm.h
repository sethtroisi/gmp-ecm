/* ecm.h - public interface for libecm.
 
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

typedef struct
{
  int method;     /* factorization method, default is ecm */
  mpz_t x;        /* starting point (if non zero) */
  mpz_t sigma;    /* contains sigma or A (ecm only) */
  int sigma_is_A; /* if non-zero, 'sigma' contains A */
  mpz_t go;       /* initial group order to preload (if NULL: do nothing) */
  double B1done;  /* step 1 was already done up to B1done */
  mpz_t B2min;    /* lower bound for stage 2 (default is B1) */
  mpz_t B2;       /* step 2 bound (chosen automatically if < 0.0) */
  unsigned long k;/* number of blocks in stage 2 */
  int S;          /* degree of the Brent-Suyama's extension for stage 2 */
  int repr;       /* representation for modular arithmetic: ECM_MOD_MPZ=mpz,         
		     ECM_MOD_MODMULN=modmuln (Montgomery's quadratic multiplication),
		     ECM_MOD_REDC=redc (Montgomery's subquadratic multiplication),
		     ECM_MOD_GWNUM=Woltman's gwnum routines (tbd),
		     > 16 : special base-2 representation        
		     MOD_DEFAULT: automatic choice */
  int verbose;    /* verbosity level: 0 no output, 1 normal output,   
		     2 diagnostic output */
  FILE *os;       /* output stream (for verbose messages) */
  FILE *es;       /* error  stream (for error   messages) */
  char *TreeFilename; /* Base filename for storing product tree of F */
  double maxmem;  /* Maximal amount of memory to use in stage 2, in bytes.
                     0. means no limit (optimise only for speed) */
} __ecm_param_struct;
typedef __ecm_param_struct ecm_params[1];

#define ECM_MOD_NOBASE2 -1
#define ECM_MOD_DEFAULT 0
#define ECM_MOD_MPZ 1
#define ECM_MOD_BASE2 2
#define ECM_MOD_MODMULN 3
#define ECM_MOD_REDC 4

int ecm_factor (mpz_t, mpz_t, double, ecm_params);
void ecm_init (ecm_params);
void ecm_clear (ecm_params);

/* the following interface is not supported */
int ecm (mpz_t, mpz_t, mpz_t, mpz_t, mpz_t, double, double, mpz_t, mpz_t,
         double, unsigned long, const int, int, int, int, FILE*, FILE*, 
         char*, double);
int pp1 (mpz_t, mpz_t, mpz_t, mpz_t, double, double, mpz_t, mpz_t, 
         double, unsigned long, const int, int, int, FILE*, FILE*, char*,
         double);
int pm1 (mpz_t, mpz_t, mpz_t, mpz_t, double, double, mpz_t, 
          mpz_t, double, unsigned long, const int, int, int, FILE*, FILE*, 
          char*, double);

/* different methods implemented */
#define ECM_ECM 0
#define ECM_PM1 1
#define ECM_PP1 2

/* return value of ecm, pm1, pp1 */
#define ECM_FACTOR_FOUND_STEP1 1 /* should be positive */
#define ECM_FACTOR_FOUND_STEP2 2 /* should be positive */
#define ECM_NO_FACTOR_FOUND 0 /* should be zero */
#define ECM_ERROR -1 /* should be non-zero */
#define ECM_FACTOR_FOUND_P(x) ((x) > 0)
#define ECM_ERROR_P(x)        ((x) < 0)

#define ECM_DEFAULT_B1_DONE 1.0
#define ECM_IS_DEFAULT_B1_DONE(x) (x <= 1.0)

/* stage 2 bound */
#define ECM_DEFAULT_B2 -1
#define ECM_IS_DEFAULT_B2(x) (mpz_sgn (x) < 0)

#define ECM_DEFAULT_K 0 /* default number of blocks in stage 2. 0 = automatic
                           choice */
#define ECM_DEFAULT_S 0 /* polynomial is chosen automatically */

/* Apple uses '\r' for newlines */
#define IS_NEWLINE(c) (((c) == '\n') || ((c) == '\r'))

