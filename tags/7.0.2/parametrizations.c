/* parametrizations.c - functions to compute coefficients of the curve from
parameter and to choose random parameter.
 
Copyright 2012 Cyril Bouvier.
 
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

#include "ecm-gmp.h"
#include "ecm-impl.h"

#if 0
/* this function is useful in debug mode to print residues */
static void
mpres_print (mpres_t x, char* name, mpmod_t n)
{
  mp_size_t m, xn;
  mpres_t t;
  mpres_init(t, n);
  mpz_set_ui(t, 1);
  mpres_mul (t, x, t, n);

  xn = SIZ(t);
  m = ABSIZ(t);
  MPN_NORMALIZE(PTR(t), m);
  SIZ(t) = xn >= 0 ? m : -m;
  gmp_printf ("%s=%Zd\n", name, t);
  SIZ(t) = xn;
  mpres_clear (t, n);
}
#endif

static void 
dbl_param (mpres_t x, mpres_t y, mpres_t z, mpres_t t, mpres_t u, mpres_t v,
                                                                  mpmod_t n)
{
  mpres_mul (z, y, z, n); /* Y1*Z1  */
  mpres_mul_ui (z, z, 2, n); /* Z3 = 2*Y1*Z1  */

  mpres_sqr (u, x, n); /* A = X1*X1  */
  mpres_sqr (t, y, n); /* B = Y1*Y1  */
  mpres_sqr (y, t, n); /* C = B^2  */
  mpres_add (v, x, t, n); /* X1+B  */
  mpres_sqr (v, v, n); /* (X1+B)^2  */
  mpres_sub (v, v, u, n); /* (X1+B)^2-A  */
  mpres_sub (v, v, y, n); /* (X1+B)^2-A-C  */
  mpres_mul_ui (v, v, 2, n); /* D = 2*((X1+B)^2-A-C)  */
  mpres_mul_ui (u, u, 3, n); /* E = 3*A  */
  mpres_sqr (t, u, n); /* F = E^2  */

  mpres_mul_ui (x, v, 2, n); /* 2*D  */
  mpres_sub (x, t, x, n); /* X3 = F-2*D  */

  mpres_sub (v, v, x, n); /* D-X3  */
  mpres_mul_ui (y, y, 8, n); /* 8*C  */
  mpres_mul (t, u, v, n); /* E*(D-X3)  */
  mpres_sub (y, t, y, n); /* Y3 = E*(D-X3)-8*C */
}

/*Add sgn*P=(-3:sgn*3:1) to Q=(x:y:z) */
static void 
add_param (mpres_t x, mpres_t y, mpres_t z, int sgn, mpres_t t, mpres_t u, 
                                          mpres_t v, mpres_t w, mpmod_t n)
{
  mpres_sqr (t, z, n); /* Z1Z1 = Z1^2   */
  mpres_mul_ui (u, t, 3, n); 
  mpres_neg (u, u, n); /* U2 = X2*Z1Z1 with X2=-3 */
  mpres_mul (v, z, t, n); /* Z1*Z1Z1  */
  mpres_mul_ui (v, v, 3, n); /* S2 = Y2*Z1*Z1Z1 with Y2=sgn*3  */
  if (sgn == -1) 
    mpres_neg (v, v, n); /* S2 = Y2*Z1*Z1Z1 with Y2=sgn*3 */
  mpres_sub (u, u, x, n); /* H = U2-X1  */
  mpres_sqr (w, u, n); /* HH = H^2  */

  mpres_add (z, z, u, n); /* Z1+H  */
  mpres_sqr (z, z, n); /* (Z1+H)^2  */
  mpres_sub (z, z, t, n); /* (Z1+H)^2-Z1Z1   */
  mpres_sub (z, z, w, n); /* Z3 = (Z1+H)^2-Z1Z1-HH  */

  mpres_mul_ui (t, w, 4, n); /* I = 4*HH  */
  mpres_mul (u, u, t, n); /* J = H*I  */
  mpres_sub (v, v, y, n); /* S2-Y1  */
  mpres_mul_ui (v, v, 2, n); /* r = 2*(S2-Y1) */
  mpres_mul (t, x, t, n); /* V = X1*I */
  mpres_sqr (x, v, n); /* r^2 */
  mpres_mul_ui (w, t, 2, n); /* 2*V  */
  mpres_sub (x, x, u, n); /* r^2-J  */
  mpres_sub (x, x, w, n); /* X3 = r^2-J-2*V  */

  mpres_sub (w, t, x, n); /* V-X3 */
  mpres_mul (y, y, u, n); /* Y1*J */
  mpres_mul_ui (y, y, 2, n); /* 2*Y1*J   */
  mpres_mul (w, v, w, n); /* r*(V-X3)  */
  mpres_sub (y, w, y, n); /* Y3=r*(V-X3)-2*Y1*J  */
}

/* computes s*(x:y:z) mod n, where t, u, v, w are temporary variables */
static void
addchain_param (mpres_t x, mpres_t y, mpres_t z, mpz_t s, mpres_t t,
                mpres_t u, mpres_t v, mpres_t w, mpmod_t n)
{
  if (mpz_cmp_ui (s, 1) == 0)
    {
      mpres_set_si (x, -3, n);
      mpres_set_ui (y, 3, n);
      mpres_set_ui (z, 1, n);
    }
  else if (mpz_cmp_ui (s, 3) == 0)
    {
      mpz_sub_ui (s, s, 1);
      addchain_param (x, y, z, s, t, u, v, w, n);
      add_param (x, y, z, +1, t, u, v, w, n);
    }
  else if (mpz_divisible_2exp_p (s, 1))
    {
      mpz_tdiv_q_2exp (s, s, 1);
      addchain_param (x, y, z, s, t, u, v, w, n);
      dbl_param (x, y, z, t, u, v, n);
    }
  else if (mpz_congruent_ui_p (s, 1, 4))
    {
      mpz_sub_ui (s, s, 1);
      addchain_param (x, y, z, s, t, u, v, w, n);
      add_param (x, y, z, +1, t, u, v, w, n);
    }
  else /* (s % 4 == 3) and s != 3 */
    {
      mpz_add_ui (s, s, 1);
      addchain_param (x, y, z, s, t, u, v, w, n);
      add_param (x, y, z, -1, t, u, v, w, n);
    }
}

/*
  get_curve_from_param* functions compute curve coefficient A and a starting
  point (x::1) from a given sigma value 
   
  If a factor of n was found during the process, returns 
  ECM_FACTOR_FOUND_STEP1 (and factor in f). 
  If a invalid value of sigma is given, returns ECM_ERROR,
  Returns ECM_NO_FACTOR_FOUND otherwise. 
*/



/* Parametrization ECM_PARAM_SUYAMA */
/* (sigma mod N) should be different from 0, 1, 3, 5, 5/3, -1, -3, -5, -5/3 */
int
get_curve_from_param0 (mpz_t f, mpres_t A, mpres_t x, mpz_t sigma, mpmod_t n)
{
  mpres_t t, u, v, b, z;
  mpz_t tmp;
  
  mpres_init (t, n);
  mpres_init (u, n);
  mpres_init (v, n);
  mpres_init (b, n);
  mpres_init (z, n);
  mpz_init (tmp);

  mpz_mod (tmp, sigma, n->orig_modulus);
  /* TODO add -5 -3 -1 and +/- 5/3 */
  if (mpz_cmp_ui (tmp, 5) == 0 || mpz_cmp_ui (tmp, 3) == 0 || 
      mpz_cmp_ui (tmp, 1) == 0 || mpz_sgn (tmp) == 0)
    return ECM_ERROR;

  mpres_set_z  (u, sigma, n);
  mpres_mul_ui (v, u, 4, n);   /* v = (4*sigma) mod n */
  mpres_sqr (t, u, n);
  mpres_sub_ui (u, t, 5, n);       /* u = (sigma^2-5) mod n */
  mpres_sqr (t, u, n);
  mpres_mul (x, t, u, n);          /* x = (u^3) mod n */
  mpres_sqr (t, v, n);
  mpres_mul (z, t, v, n);          /* z = (v^3) mod n */
  mpres_mul (t, x, v, n);
  mpres_mul_ui (b, t, 4, n);       /* b = (4*x*v) mod n */
  mpres_mul_ui (t, u, 3, n);
  mpres_sub (u, v, u, n);          /* u' = v-u */
  mpres_add (v, t, v, n);          /* v' = (3*u+v) mod n */
  mpres_sqr (t, u, n);
  mpres_mul (u, t, u, n);          /* u'' = ((v-u)^3) mod n */
  mpres_mul (A, u, v, n);          /* a = (u'' * v') mod n = 
                                      ((v-u)^3 * (3*u+v)) mod n */
  
  /* Normalize b and z to 1 */
  mpres_mul (v, b, z, n);
  if (!mpres_invert (u, v, n)) /* u = (b*z)^(-1) (mod n) */
    {
      mpres_gcd (f, v, n);
      mpres_clear (t, n);
      mpres_clear (u, n);
      mpres_clear (v, n);
      mpres_clear (b, n);
      mpres_clear (z, n);
      mpz_clear (tmp);
      if (mpz_cmp (f, n->orig_modulus) == 0)
          return ECM_ERROR;
      else
          return ECM_FACTOR_FOUND_STEP1;
    }
  
  mpres_mul (v, u, b, n);   /* v = z^(-1) (mod n) */
  mpres_mul (x, x, v, n);   /* x = x * z^(-1) */
  
  mpres_mul (v, u, z, n);   /* v = b^(-1) (mod n) */
  mpres_mul (t, A, v, n);
  mpres_sub_ui (A, t, 2, n);
  
  mpres_clear (t, n);
  mpres_clear (u, n);
  mpres_clear (v, n);
  mpres_clear (b, n);
  mpres_clear (z, n);
  mpz_clear (tmp);

  return ECM_NO_FACTOR_FOUND;
}

/* Parametrization ECM_PARAM_BATCH_SQUARE */
/* Only work for 64-bit machines */
/* d = (sigma^2/2^64 mod N) should be different from 0, 1, -1/8 */
int  
get_curve_from_param1 (mpres_t A, mpres_t x0, mpz_t sigma, mpmod_t n)
{
  int i;
  mpz_t tmp;
  mpz_init (tmp);

  ASSERT (GMP_NUMB_BITS == 64);
      
  mpz_mul (tmp, sigma, sigma); /* tmp = sigma^2*/
  
  /* A=4*d-2 with d = sigma^2/2^GMP_NUMB_BITS*/
  /* Compute d = sigma^2/2^GMP_NUMB_BITS */
  for (i = 0; i < GMP_NUMB_BITS; i++)
    {
      if (mpz_tstbit (tmp, 0) == 1)
          mpz_add (tmp, tmp, n->orig_modulus);
      mpz_div_2exp (tmp, tmp, 1);
    }

  mpz_mod (tmp, tmp, n->orig_modulus);
  /* TODO add d!=-1/8*/
  if (mpz_sgn (tmp) == 0 || mpz_cmp_ui (tmp, 1) == 0)
      return ECM_ERROR;

  mpz_mul_2exp (tmp, tmp, 2);           /* 4d */
  mpz_sub_ui (tmp, tmp, 2);             /* 4d-2 */
      
  mpres_set_z (A, tmp, n);
  mpres_set_ui (x0, 2, n);

  mpz_clear(tmp);
  return ECM_NO_FACTOR_FOUND;
}

/* Parametrization ECM_PARAM_BATCH_2 */
/* 2 <= sigma */
/* Compute (x:y:z) = sigma*(-3:3:1) on the elliptic curve y^2 = x^3 + 36
  using Jacobian coordinates; formulae were found at
      https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html
  Then we let x3 = (3*x+y+6)/(2*(y-3)), A = -(3*x3^4+6*x3^2-1)/(4*x3^3) and
  x0 = 2. The Sage code below gives the factored group order:

  def FindGroupOrderParam2(p,sigma):
   K = GF(p)
   E = EllipticCurve(K,[0,36])
   P = sigma*E(-3,3)
   x,y = P.xy()
   x3 = (3*x+y+6)/(2*(y-3))
   A = -(3*x3^4+6*x3^2-1)/(4*x3^3)
   d = K((A+2)/4)
   a = K(4*d-2)
   b = K(16*d+2)
   E = EllipticCurve(K,[0,a/b,0,1/b^2,0])
   return factor(E.cardinality())
*/
int 
get_curve_from_param2 (mpz_t f, mpres_t A, mpres_t x0, mpz_t sigma, mpmod_t n)
{
  mpres_t t, u, v, w, x, y, z;
  mpz_t k;
  int ret = ECM_NO_FACTOR_FOUND;

  mpres_init (t, n);
  mpres_init (u, n);
  mpres_init (v, n);
  mpres_init (w, n);
  mpres_init (x, n);
  mpres_init (y, n);
  mpres_init (z, n);
  mpz_init (k);

  mpz_set (k, sigma);

  if (mpz_cmp_ui (k, 2) < 0)
    {
      ret = ECM_ERROR;
      goto clear_and_exit;
    }

  addchain_param (x, y, z, k, t, u, v, w, n);

  /* Now (x:y:z) = k*P */

  if (!mpres_invert (u, z, n)) 
    {
      mpres_gcd (f, z, n);
      ret = ECM_FACTOR_FOUND_STEP1;
      goto clear_and_exit;
    }

  mpres_sqr (v, u, n);
  mpres_mul (u, v, u, n);
  mpres_mul (x, x, v, n);  
  mpres_mul (y, y, u, n);  

  mpres_sub_ui (t, y, 3, n);
  mpres_mul_ui (t, t, 2, n);

  if (!mpres_invert (u, t, n)) 
    {
      mpres_gcd (f, t, n);
      ret = ECM_FACTOR_FOUND_STEP1;
      goto clear_and_exit;
    }
  
  mpres_mul_ui (w, x, 3, n);
  mpres_add (w, w, y, n);
  mpres_add_ui (w, w, 6, n);
  mpres_mul (x, w, u, n);   /* Now x contains x_3 */  

  /* A=-(3*x3^4+6*x3^2-1)/(4*x3^3) */
  mpres_sqr (u, x, n);
  mpres_mul (v, u, x, n);
  mpres_sqr (w, u, n);

  mpres_mul_ui (u, u, 6, n);
  mpres_neg (u, u, n);
  mpres_mul_ui (v, v, 4, n);
  mpres_mul_ui (w, w, 3, n);
  mpres_neg (w, w, n);

  if (!mpres_invert (t, v, n)) 
    {
      mpres_gcd (f, v, n);
      ret = ECM_FACTOR_FOUND_STEP1;
      goto clear_and_exit;
    }

  mpres_add (w, w, u, n);
  mpres_add_ui (w, w, 1, n);
  mpres_mul (A, w, t, n);
  mpz_mod (A, A, n->orig_modulus); 

  mpres_set_ui (x0, 2, n);

 clear_and_exit:
  mpres_clear (t, n);
  mpres_clear (u, n);
  mpres_clear (v, n);
  mpres_clear (w, n);
  mpres_clear (x, n);
  mpres_clear (y, n);
  mpres_clear (z, n);
  mpz_clear (k);

  return ret;
}

/* Parametrization ECM_PARAM_BATCH_32BITS_D */
/* d = (sigma/2^32 mod N) should be different from 0, 1, -1/8 */
int  
get_curve_from_param3 (mpres_t A, mpres_t x0, mpz_t sigma, mpmod_t n)
{
  int i;
  mpz_t tmp;
  mpz_t two32;
  mpz_init (two32);
  mpz_ui_pow_ui (two32, 2, 32);
  mpz_init (tmp);

  /* sigma < 2^32 (it was generated for 32-bit machines) */
  /* To use it on a 64-bits machines one should multiplied it by 2^32 */
  if (GMP_NUMB_BITS == 64)
      mpz_mul (tmp, sigma, two32);
  else  
      mpz_set (tmp, sigma);
  
  /* A=4*d-2 with d = sigma/2^GMP_NUMB_BITS*/
  /* Compute d = sigma/2^GMP_NUMB_BITS */
  for (i = 0; i < GMP_NUMB_BITS; i++)
    {
      if (mpz_tstbit (tmp, 0) == 1)
      mpz_add (tmp, tmp, n->orig_modulus);
      mpz_div_2exp (tmp, tmp, 1);
    }

  mpz_mod (tmp, tmp, n->orig_modulus);
  /* TODO add d!=-1/8*/
  if (mpz_sgn (tmp) == 0 || mpz_cmp_ui (tmp, 1) == 0)
      return ECM_ERROR;

  mpz_mul_2exp (tmp, tmp, 2);           /* 4d */
  mpz_sub_ui (tmp, tmp, 2);             /* 4d-2 */
      
  mpres_set_z (A, tmp, n);
  mpres_set_ui (x0, 2, n);

  mpz_clear(tmp);
  mpz_clear (two32);
  return ECM_NO_FACTOR_FOUND;
}

int
get_curve_from_random_parameter (mpz_t f, mpres_t A, mpres_t x, mpz_t sigma, 
                                 int param, mpmod_t modulus, gmp_randstate_t rng)
{
  int ret;

  /* initialize the random number generator if not already done */
  init_randstate (rng);
  do
    {
      if (param == ECM_PARAM_SUYAMA)
        {
          mpz_urandomb (sigma, rng, 64);
          ret = get_curve_from_param0 (f, A, x, sigma, modulus);
        }
      else if (param == ECM_PARAM_BATCH_SQUARE)
        {
          mpz_urandomb (sigma, rng, 32);
          ret = get_curve_from_param1 (A, x, sigma, modulus);
        }
      else if (param == ECM_PARAM_BATCH_2)
        {
          mpz_urandomb (sigma, rng, 64);
          ret = get_curve_from_param2 (f, A, x, sigma, modulus);
        }
      else if (param == ECM_PARAM_BATCH_32BITS_D)
        {
          mpz_urandomb (sigma, rng, 32);
          ret = get_curve_from_param3 (A, x, sigma, modulus);
        }
      else
        return ECM_ERROR;
    } while (ret == ECM_ERROR);

  return ret;
}

int 
get_default_param (int sigma_is_A, double B1done, int repr)
{
  /* if B1done is not the default value, use ECM_PARAM_SUYAMA, since
     ECM_PARAM_BATCH* requires B1done is the default */
  if (!ECM_IS_DEFAULT_B1_DONE(B1done))
      return ECM_PARAM_SUYAMA;
  
  if (sigma_is_A == 1 || sigma_is_A == -1)
    {
      /* For now we keep the default values in order not to compute the 
         expected number of curves. But it will do stage 1 as ECM_PARAM_SUYAMA */
      return ECM_PARAM_DEFAULT; 
    }

  /* ECM_PARAM_BATCH* requires ECM_MOD_MODMULN */
  if (repr != ECM_MOD_MODMULN)
    return ECM_PARAM_SUYAMA;

  if (GMP_NUMB_BITS == 64)
    return ECM_PARAM_BATCH_SQUARE;
  else
    return ECM_PARAM_BATCH_32BITS_D;
}
