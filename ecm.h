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

/* Input: x is the starting point or zero
          sigma is sigma value (if x is set to zero) or 
            A parameter (if x is non-zero) of curve
          n is the number to factor
          B1, B2 are the stage 1/stage 2 bounds, respectively
          k is the number of blocks to do in stage 2
          S is the degree of the Suyama-Brent extension for stage 2
          verbose is verbosity level: 0 no output, 1 normal output,
            2 diagnostic output.
          sigma_is_a: If true, the sigma parameter contains the curve's A value
   Output: f is the factor found.
   Return value: non-zero iff a factor was found.
*/
int /* return value is non-zero iff a factor was found */
ecm (mpz_t, /* f: factor found if any */
     mpz_t, /* x: starting point (if non zero) */
     mpz_t, /* sigma: contains sigma or A */
     mpz_t, /* n: number to factor */
     mpz_t, /* go: initial group order to preload (if NULL: do nothing) */
     double, /* B1done: step 1 was already done up to B1done */
     double, /* B1: step 1 bound */
     double, /* B2min: lower bound for stage 2 (default is B1) */
     double, /* B2: step 2 bound (chosen automatically if < 0.0) */
     double, /* B2scale: B2 scale factor */
     unsigned int, /* k: number of blocks in stage 2 */
     int,          /* S: degree of the Brent-Suyama's extension for stage 2 */
     int,          /* verbose: verbosity level: 0 no output, 1 normal output,
                      2 diagnostic output */
     int,          /* representation for modular arithmetic: 1=mpz,
                      2=modmuln (Montgomery's quadratic multiplication),
                      3=redc (Montgomery's subquadratic multiplication),
                      > 16 : special base-2 representation
                      otherwise: automatic choice */
     int,          /* sigma_is_A: if non-zero, 'sigma' contains A */
     FILE*,        /* standard output (for verbose messages) */
     FILE*);       /* standard error  (for error   messages) */

int pp1 (mpz_t, mpz_t, mpz_t, mpz_t, double, double, double, 
          double, double, unsigned int, unsigned int, int, int, FILE*, FILE*);
int pm1 (mpz_t, mpz_t, mpz_t, mpz_t, double, double, double, 
          double, double, unsigned int, int, int, int, FILE*, FILE*);

/* different methods implemented */
#define ECM_ECM 0
#define ECM_PM1 1
#define ECM_PP1 2

/* return value of ecm, pm1, pp1 */
#define ECM_FACTOR_FOUND 1 /* should be non-zero */
#define ECM_NO_FACTOR_FOUND 0 /* should be zero */
#define ECM_ERROR -1 /* should be non-zero */

#define DEFAULT_B1_DONE 1.0
#define IS_DEFAULT_B1_DONE(x) (x <= 1.0)

#define ECM_DEFAULT_K 2 /* default number of blocks in stage 2 */
