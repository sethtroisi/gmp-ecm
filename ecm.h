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
  int method;  /* factorization method, default is ecm */
  mpz_t x;     /* starting point (if non zero) */
  mpz_t sigma; /* contains sigma or A (ecm only) */
  int sigma_is_A; /* if non-zero, 'sigma' contains A */
  mpz_t go;    /* initial group order to preload (if NULL: do nothing) */
  double B1done; /* step 1 was already done up to B1done */
  double B2min;  /* lower bound for stage 2 (default is B1) */
  double B2;     /* step 2 bound (chosen automatically if < 0.0) */
  unsigned int k; /* number of blocks in stage 2 */
  int S;          /* degree of the Brent-Suyama's extension for stage 2 */
  int repr;       /* representation for modular arithmetic: 1=mpz,         
		     2=modmuln (Montgomery's quadratic multiplication),
		     3=redc (Montgomery's subquadratic multiplication),
		     > 16 : special base-2 representation        
		     otherwise: automatic choice */
  int verbose;    /* verbosity level: 0 no output, 1 normal output,   
		     2 diagnostic output */
  FILE *os;       /* output stream (for verbose messages) */
  FILE *es;       /* error  stream (for error   messages) */
} __ecm_param_struct;
typedef __ecm_param_struct ecm_params[1];

int ecm_factor (mpz_t, mpz_t, double, ecm_params);
void ecm_init (ecm_params);
void ecm_clear (ecm_params);

/* the following interface is not supported */
int ecm (mpz_t, mpz_t, mpz_t, mpz_t, mpz_t, double, double, double, double,
         double, unsigned int, int, int, int, int, FILE*, FILE*);
int pp1 (mpz_t, mpz_t, mpz_t, mpz_t, double, double, double, 
          double, double, unsigned int, unsigned int, int, int, FILE*, FILE*);
int pm1 (mpz_t, mpz_t, mpz_t, mpz_t, double, double, double, 
          double, double, unsigned int, int, int, int, FILE*, FILE*);

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
#define ECM_DEFAULT_B2 -1.0
#define ECM_IS_DEFAULT_B2(x) ((x) < 0.0)

#define ECM_DEFAULT_K 2 /* default number of blocks in stage 2 */
#define ECM_DEFAULT_S 0 /* polynomial is chosen automatically */
#define ECM_DEFAULT_REPR 0 /* automatic choice */
