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

static void
addchain_param (mpres_t x, mpres_t y, mpres_t z, unsigned int s, mpres_t t,
                                    mpres_t u, mpres_t v, mpres_t w, mpmod_t n)
{
  if (s == 1)
    {
      mpres_set_si (x, -3, n);
      mpres_set_ui (y, 3, n);
      mpres_set_ui (z, 1, n);
    }
  else if (s == 3)
    {
      addchain_param(x, y, z, s-1, t, u, v, w, n);
      add_param (x, y, z, +1, t, u, v, w, n);
    }
  else if (s % 2 == 0)
    {
      addchain_param(x, y, z, s/2, t, u, v, w, n);
      dbl_param (x, y, z, t, u, v, n);
    }
  else if (s % 4 == 1)
    {
      addchain_param(x, y, z, s-1, t, u, v, w, n);
      add_param (x, y, z, +1, t, u, v, w, n);
    }
  else /* (s % 4 == 3) and s != 3 */
    {
      addchain_param(x, y, z, s+1, t, u, v, w, n);
      add_param (x, y, z, -1, t, u, v, w, n);
    }
}

/*Parametrization for BATCHMODE 2: generate curves with a point of order 3 and
  starting point (2:1) 
  Compute k*P on y^2=x^3+36 with P=(-3,3); need k>1
  x3 = (3*x+y+6)/(2*(y-3)) and A=-(3*x3^4+6*x3^2-1)/(4*x3^3)*/
int 
get_curve_from_ell_parametrization (mpz_t f, mpres_t A, mpz_t k, mpmod_t n)
{
  mpres_t t, u, v, w, x, y, z;
  unsigned int s;

  MEMORY_TAG;
  mpres_init (t, n);
  MEMORY_TAG;
  mpres_init (u, n);
  MEMORY_TAG;
  mpres_init (v, n);
  MEMORY_TAG;
  mpres_init (w, n);
  MEMORY_TAG;
  mpres_init (x, n);
  MEMORY_TAG;
  mpres_init (y, n);
  MEMORY_TAG;
  mpres_init (z, n);
  MEMORY_UNTAG;

  s = mpz_get_ui (k);

  addchain_param (x, y, z, s, t, u, v, w, n); 

  /* Now (x:y:z) = k*P */

  if (!mpres_invert(u, z, n)) 
    {
      mpres_gcd (f, z, n);
      mpres_clear (t, n);
      mpres_clear (u, n);
      mpres_clear (v, n);
      mpres_clear (w, n);
      mpres_clear (x, n);
      mpres_clear (y, n);
      mpres_clear (z, n);
      return ECM_FACTOR_FOUND_STEP1;
    }

  mpres_sqr (v, u, n);
  mpres_mul (u, v, u, n);
  mpres_mul (x, x, v, n);  
  mpres_mul (y, y, u, n);  

  mpres_sub_ui (t, y, 3, n);
  mpres_mul_ui (t, t, 2, n);

  if (!mpres_invert(u, t, n)) 
    {
      mpres_gcd (f, t, n);
      mpres_clear (t, n);
      mpres_clear (u, n);
      mpres_clear (v, n);
      mpres_clear (w, n);
      mpres_clear (x, n);
      mpres_clear (y, n);
      mpres_clear (z, n);
      return ECM_FACTOR_FOUND_STEP1;
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

  if (!mpres_invert(t, v, n)) 
    {
      mpres_gcd (f, v, n);
      mpres_clear (t, n);
      mpres_clear (u, n);
      mpres_clear (v, n);
      mpres_clear (w, n);
      mpres_clear (x, n);
      mpres_clear (y, n);
      mpres_clear (z, n);
      return ECM_FACTOR_FOUND_STEP1;
    }

  mpres_add (w, w, u, n);
  mpres_add_ui (w, w, 1, n);
  mpres_mul (A, w, t, n);
  mpz_mod (A, A, n->orig_modulus); 

  mpres_clear (t, n);
  mpres_clear (u, n);
  mpres_clear (v, n);
  mpres_clear (w, n);
  mpres_clear (x, n);
  mpres_clear (y, n);
  mpres_clear (z, n);

  return ECM_NO_FACTOR_FOUND;
}
