/* torsions.c - ECM with special torsion curves
   Author: F. Morain
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gmp.h> /* GMP header file */

#include "ecm.h" /* ecm header file */
#include "ecm-impl.h"
#include "ecm-ecm.h"
#include "mpmod.h"

#include "addlaws.h"
#include "torsions.h"

#define DEBUG_TORSION 0

/* We use three variants of Weierstrass parametrization:
   CW (complete): y^2+a1*x*y+a3*y=x^3+a2*x^2+a4*x+a6
   MW (medium)  : y^2=x^3+a2*x^2+a4*x+a6
   SW (short)   : y^2=x^3+a4*x+a6

   A Kubert curve is the special case Y^2+(1-c)*X*Y-b*Y = X^3-b*X^2

   Generally, we build a curve under the SW form, with affine law, meaning
   that constructed points will be [x, y, 1].
 */

/********** utilities **********/

void
mod_div_2(mpz_t x, mpz_t n)
{
    if(mpz_tstbit(x, 0)){
	/* x is odd, x/2 = (x+N)/2 */
	mpz_add(x, x, n);
	mpz_div_2exp(x, x, 1);
    }
    else
	/* x is even, x/2 is easy */
	mpz_div_2exp(x, x, 1);
}

/* r <- q mod N. 
   Return value: 1 if den invertible, 0 if factor found; in this case
   gcd(den(q), N) is put in r.
 */
int
mod_from_rat2(mpz_t r, mpz_t num, mpz_t den, mpz_t N)
{
    int ret = 1;
 
    if(mpz_invert(r, den, N) == 0){
	mpz_gcd(r, den, N);
	ret = 0;
    }
    else{
	mpz_mul(r, r, num);
	mpz_mod(r, r, N);
    }
    return ret;
}

int
mod_from_rat_str(mpz_t r, char *str, mpz_t N)
{
    mpq_t q;
    int ret;

    mpq_init(q);
    mpq_set_str(q, str, 10);
    ret = mod_from_rat2(r, mpq_numref(q), mpq_denref (q), N);
    mpq_clear(q);
    return ret;
}

/* From a curve in Kubert form Y^2+(1-c)*X*Y-b*Y = X^3-b*X^2
   to a Weierstrass form y^2 = X^3 + a2 * X^2 + a4 * X + a6
   where y = Y+((1-c)*X-b)/2
   WE:=[0,(1/4*c^2+1/4-1/2*c-b),0,(1/2*c*b-1/2*b),1/4*b^2]);
   We compute:
   a2 = 1/4*c^2+1/4-1/2*c-b = ((c-1)/2)^2-b
   a4 = 1/2*c*b-1/2*b = b*(c-1)/2
   a6 = (b/2)^2
   TODO: rewrite this with MediumW, etc.
*/
void
KW2W246(mpz_t a2, mpz_t a4, mpz_t a6, mpz_t b, mpz_t c, mpz_t n, int compute_a6)
{
    /** a4 <- (c-1)/2 **/
    mpz_sub_si(a4, c, 1);
    mod_div_2(a4, n);
    /** a2 <- a4^2-b **/
    mpz_mul(a2, a4, a4);
    mpz_sub(a2, a2, b);
    mpz_mod(a2, a2, n);
    /** a4 <- a4*b **/
    mpz_mul(a4, a4, b);
    mpz_mod(a4, a4, n);
    if(compute_a6 != 0){
	mpz_set(a6, b);
	mod_div_2(a6, n);
	mpz_mul(a6, a6, a6);
	mpz_mod(a6, a6, n);
    }
#if DEBUG_TORSION >= 2
    gmp_printf("N:=%Zd;\n", n);
    gmp_printf("b:=%Zd;\n", b);
    gmp_printf("c:=%Zd;\n", c);
    gmp_printf("a2:=%Zd;\n", a2);
    gmp_printf("a4:=%Zd;\n", a4);
    printf("a6:=RatMod(b^2/4, N);\n");
    if(compute_a6 != 0)
	gmp_printf("a6:=%Zd;\n", a6);
#endif
}

static int
check_weierstrass(mpz_t A, mpz_t B, mpz_t X, mpz_t Y, mpz_t tmp1, mpz_t tmp2,
		  mpz_t n)
{
    mpz_mul(tmp1, Y, Y);
    mpz_mul(tmp2, X, X);
    mpz_add(tmp2, tmp2, A);
    mpz_mul(tmp2, tmp2, X);
    mpz_add(tmp2, tmp2, B);
    mpz_sub(tmp1, tmp1, tmp2);
    mpz_mod(tmp1, tmp1, n);
    return mpz_sgn(tmp1) == 0;
}

/* Weierstrass (a2, a4, a6) to (A, B)
   A = (a4-1/3*a2^2)
   B = -1/3*a4*a2 + 2/27*a2^3 + a6
     = -1/3*a2*(a4-2/9*a2^2) + a6
   X = x+a2/3
   Y = y
   INPUT: if x0 == NULL, we have no point to translate
          if B == NULL, we do not need and we do not compute B
   REM: we assume gcd(n, 3) = 1.
*/
void
MediumWeierstrassToShortWeierstrass(mpz_t A, mpz_t B, mpz_t X, mpz_t Y,
				    mpz_t a2, mpz_t a4, mpz_t a6, 
				    mpz_t x0, mpz_t y0, mpz_t n)
{
    mpz_t tmp1, tmp2, tmp3;

    mpz_init(tmp1);
    mpz_init(tmp2);
    /* tmp2 <- a2/3 */
    mpz_init_set_si(tmp3, 3);
    mod_from_rat2(tmp2, a2, tmp3, n);
    if(X != NULL && x0 != NULL){
	/* wx0 = x0 + a2/3 */
	mpz_add(X, tmp2, x0);
	mpz_mod(X, X, n);
    }
    if(Y != NULL && y0 != NULL){
	mpz_set(Y, y0);
	mpz_mod(Y, Y, n);
    }
    /* A = a4-1/3*a2^2 = a4 - a2 * (a2/3) */
    /** tmp1 <- tmp2*a2 = a2^2/3 */
    mpz_mul(tmp1, a2, tmp2);
    mpz_mod(tmp1, tmp1, n);
    mpz_sub(A, a4, tmp1);
    mpz_mod(A, A, n);
    if(B != NULL){
	/* B = -1/3*a2*(a4-2/9*a2^2) + a6 */
	/** B <- 2/9*a2^2 = 2 * (a2^2/3) / 3 **/
	mod_from_rat2(B, tmp1, tmp3, n);
	mpz_mul_si(B, B, 2);
	mpz_sub(B, a4, B);
	mpz_mul(B, B, tmp2);
	mpz_sub(B, a6, B);
	mpz_mod(B, B, n);
    }
#if DEBUG_TORSION >= 2
    gmp_printf("N:=%Zd;\n", n);
    gmp_printf("a2:=%Zd; a4:=%Zd; a6:=%Zd;\n", a2, a4, a6);
    gmp_printf("A:=%Zd; B:=%Zd;\n", A, B);
    if(X != NULL && x0 != NULL){
	gmp_printf("x:=%Zd;\n", x0);
	gmp_printf("X:=%Zd;\n", X);
    }
    if(Y != NULL && y0 != NULL){
	gmp_printf("y:=%Zd;\n", Y);
	printf("(y^2-x^3-a2*x^2-a4*x-a6) mod N;\n");
	printf("(y^2-X^3-A*X-B) mod N;\n");
    }
#endif
    mpz_clear(tmp1);
    mpz_clear(tmp2);
    mpz_clear(tmp3);
}

/* 
   The original Kubert curve E(b, c) is y^2+(1-c)*x*y-b*y = x^3-b*x^2
   The medium Weierstrass form is y^2=x^3+a2*x^2+a4*x+a6 with point (x0, y0);
   we convert this to short Weierstrass form:
   E: Y^2 = X^3 + A * X + B
   and point P=(X, Y) on E.
*/
void
kubert_to_weierstrass(mpz_t A, mpz_t B, mpz_t X, mpz_t Y, 
		      mpz_t b, mpz_t c, mpz_t x0, mpz_t y0, mpz_t n)
{
    mpz_t a2, a4, a6;

    mpz_init(a2);
    mpz_init(a4);
    mpz_init(a6);
    KW2W246(a2, a4, a6, b, c, n, 1);
    /* second conversion */
    MediumWeierstrassToShortWeierstrass(A, B, X, Y, a2, a4, a6, x0, y0, n);
#if DEBUG_TORSION >= 2
    gmp_printf("a2:=%Zd; a4:=%Zd; a6:=%Zd; A:=%Zd; B:=%Zd;\n", a2, a4, a6,A,B);
    gmp_printf("X:=%Zd; Y:=%Zd;\n", X, Y);
#endif
    mpz_clear(a2);
    mpz_clear(a4);
    mpz_clear(a6);
}

static int
forbidden(char *torsion, int u){
    if(strcmp(torsion, "Z10") == 0)
	return u == 1 || u == 2;
    else if(strcmp(torsion, "Z3xZ3") == 0)
	return u == 2;
    return 0;
}

/* Kubert: put b = c. 
   SIDE EFFECT: tE[0..nE[ and tP[0..nE[ receive a curve of torsion Z5
                and a point on it using parameters [smin..smax[.
   OUTPUT: ECM_NO_FACTOR_FOUND or ECM_FACTOR_FOUND_STEP1 if a factor is found.
*/
int
build_curves_with_torsion_Z5(mpz_t f, mpmod_t n,
			     ell_curve_t *tE, ell_point_t *tP,
			     int smin, int smax, int nE)
{
    mpz_t A, B, X, Y;
    int s, ret = ECM_NO_FACTOR_FOUND, nc = 0;
    mpz_t x0, y0, c, tmp;

    mpz_init(A);
    mpz_init(B);
    mpz_init(X);
    mpz_init(Y);
    mpz_init(x0);
    mpz_init(y0);
    mpz_init(c);
    mpz_init(tmp);
    for(s = smin; s < smax; s++){
	mpz_set_si(x0, s);
	/* c:=1/2*x0*(4*x0+1)/(3*x0+1); */
	/* y0 <- 2*(3*x0+1) */
	mpz_mul_si(y0, x0, 3);
	mpz_add_si(y0, y0, 1);
	/* tmp <- x0*(4*x0+1) */
	mpz_add(tmp, y0, x0);
	mpz_mul(tmp, tmp, x0);
	mpz_add(y0, y0, y0);
	if(mod_from_rat2(c, tmp, y0, n->orig_modulus) == 0){
	    printf("factor found during Z5_init\n");
	    mpz_gcd(f, c, n->orig_modulus);
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
	/* y0:=x0*(x0+1)*(4*x0+1)/4/(3*x0+1) = (x0+1)*c/2 */
	mpz_add_si(y0, x0, 1);
	mpz_mul(y0, y0, c);
	mpz_mod(y0, y0, n->orig_modulus);
	mod_div_2(y0, n->orig_modulus);
#if DEBUG_TORSION >= 2
	gmp_printf("x0:=%Zd;\nc:=%Zd;\ny0:=%Zd;\n", x0, c, y0);
	printf("cr:=1/2*x0*(4*x0+1)/(3*x0+1);\n");
#endif
	/* P:=WE![x0, y0, 1]; */
	kubert_to_weierstrass(A, B, X, Y, c, c, x0, y0, n->orig_modulus);
	if(check_weierstrass(A, B, X, Y, tmp, x0, n->orig_modulus) == 0){
	    printf("#!# check_weierstrass false\n");
	    ret = ECM_ERROR;
	    break;
	}
	ell_curve_init(tE[nc], ECM_EC_TYPE_WEIERSTRASS, ECM_LAW_AFFINE,n);
	mpz_set(tE[nc]->a4, A);
	mpz_set(tE[nc]->a6, B);
	ell_point_init(tP[nc], tE[nc], n);
	mpz_set(tP[nc]->x, X);
	mpz_set(tP[nc]->y, Y);
	nc++;
	if(nc >= nE)
	    break;
    }
    mpz_clear(A);
    mpz_clear(B);
    mpz_clear(X);
    mpz_clear(Y);
    mpz_clear(x0);
    mpz_clear(y0);
    mpz_clear(c);
    mpz_clear(tmp);
    return ret;
}

/* 
     E_aux: T^2 = S^3 + A * S + B
   => quartic QC: Y^2 = X^4 - 6 * A2 * X^2 + 4 * A1 * X + A0, with
     X = (T-A1/2)/(S-A2), Y = -X^2 + 2 * S + A2.
   => quartic y^2 = f(x) = a4*x^4+...+a0, where
     x = x0+y0/(X-cte), where cte = f'(x0)/4/y0
     y = Y/y0*(x-x0)^2 = Y*y0/(X-cte)^2
   INPUT: (s, t) is a point on E_aux; (x0, y0) a point on QC.
   SIDE EFFECT: x, y contain a point on the elliptic curve.
   OUTPUT: 1 if no pb occurred,
           0 if a factor was found and put in f
 */
int
cubic_to_quartic(mpz_t f, mpz_t n, mpz_t x, mpz_t y,
		 mpz_t s, mpz_t t, mpz_t A2, mpz_t A1div2,
		 mpz_t x0, mpz_t y0, mpz_t cte)
{
    mpz_t X, Y;
    int ret = 1;

    mpz_init(X);
    mpz_init(Y);
    /* X <- (t-A1/2)/(s-A2) */
    mpz_sub(x, t, A1div2);
    mpz_sub(y, s, A2);
    if(mod_from_rat2(X, x, y, n) == 0){
	mpz_set(f, X);
	ret = 0;
    }
    else{
	/* Y <- -X^2 + 2 * s + A2 */
	mpz_mul(Y, X, X);
	mpz_sub(Y, A2, Y);
	mpz_add(Y, Y, s);
	mpz_add(Y, Y, s);
	mpz_mod(Y, Y, n);
	/* x <- x0+y0/(X-cte) */
	mpz_sub(X, X, cte);
	mpz_mod(X, X, n);
	if(mpz_invert(f, X, n) == 0){
	    mpz_gcd(f, X, n);
	    ret = 0;
	}
	else{
	    /* x <- y0/(X-cte) */
	    mpz_mul(x, f, y0);
	    mpz_mod(x, x, n);
	    /* y <- x/(X-cte) = y0/(X-cte)^2 */
	    mpz_mul(y, x, f);
	    mpz_mod(y, y, n);
	    mpz_mul(y, y, Y);
	    mpz_mod(y, y, n);
	    mpz_add(x, x, x0);
	    mpz_mod(x, x, n);
	}
    }
    mpz_clear(X);
    mpz_clear(Y);
    return ret;
}

int
build_curves_with_torsion_aux(ell_curve_t Eaux, ell_point_t Paux,
			      mpz_t A2, mpz_t A1div2, mpz_t x0, mpz_t y0,
			      mpz_t cte,
			      char *sa4, char *sa6, char *sPx, char *sPy,
			      char *sA2, char *sA1div2, char *sx0, char *sy0,
			      char *scte, mpmod_t n, mpres_t tmp)
{
    mpz_t f;

    mpz_init(f);
    mod_from_rat_str(f, sa4, n->orig_modulus);
    mpres_set_z(tmp, f, n);
    ell_curve_init_set(Eaux, ECM_EC_TYPE_WEIERSTRASS, ECM_LAW_AFFINE, tmp, n);
    mod_from_rat_str(f, sa6, n->orig_modulus);
    mpres_set_z(Eaux->a6, f, n);
    ell_point_init(Paux, Eaux, n);
    mod_from_rat_str(f, sPx, n->orig_modulus);
    mpres_set_z(Paux->x, f, n);
    mod_from_rat_str(f, sPy, n->orig_modulus);
    mpres_set_z(Paux->y, f, n);
#if DEBUG_TORSION >= 2
    printf("Paux:=");
    pt_print(Eaux, Paux, n);
    printf(";\n");
#endif
    mod_from_rat_str(A2, sA2, n->orig_modulus);
    mod_from_rat_str(A1div2, sA1div2, n->orig_modulus);
    mod_from_rat_str(x0, sx0, n->orig_modulus);
    mod_from_rat_str(y0, sy0, n->orig_modulus);
    mod_from_rat_str(cte, scte, n->orig_modulus);
    mpz_clear(f);
    return 1;
}

/* 
   SIDE EFFECT: tE[0..nE[ and tP[0..nE[ receive a curve of torsion Z7
                and a point on it using parameters [umin..umax[.
   OUTPUT: ECM_NO_FACTOR_FOUND or ECM_FACTOR_FOUND_STEP1 if a factor is found.
   tE[i], tP[i] are built in raw modular form, not Montgomery form. 
   REM: we assume gcd(n, 6).
*/
int
build_curves_with_torsion_Z7(mpz_t fac, mpmod_t n, 
			     ell_curve_t *tE, ell_point_t *tP,
			     int umin, int umax, int nE)
{
    int u, ret = ECM_NO_FACTOR_FOUND, nc = 0;
    mpz_t A2, A1div2, x0, y0, cte, d, c, b, kx0, ky0, A, B, X, Y;
    mpres_t tmp;
    ell_curve_t E;
    ell_point_t P, Q;

    mpz_init(A2);
    mpz_init(A1div2);
    mpz_init(cte);
    mpz_init(x0);
    mpz_init(y0);
    mpz_init(A);
    mpz_init(B);
    mpz_init(X);
    mpz_init(Y);
    /* Eaux = "1295/48", "-1079/864" */
    /* Paux = "2185/12", "-2458" */
    /* Y^2 = X^4-1/2*X^2-8*X-1727/16 */
    mpres_init(tmp, n);
    build_curves_with_torsion_aux(E, P, A2, A1div2, x0, y0, cte,
				  "1295/48", "-1079/864",
				  "2185/12", "-2458",
				  "1/12", "-1",
				  "-1", "8", "-7/2",
				  n, tmp);
    mpz_init(d);
    mpz_init(c);
    mpz_init(b);
    mpz_init(kx0);
    mpz_init(ky0);
    ell_point_init(Q, E, n);
    mpz_set_si(d, umin-1);
    if(ell_point_mul(fac, Q, d, P, E, n) == 0){
	printf("found factor during init of Q in Z7\n");
	ret = ECM_FACTOR_FOUND_STEP1;
    }
    for(u = umin; (ret != ECM_FACTOR_FOUND_STEP1) && u < umax; u++){
	/* update Q */
	if(ell_point_add(fac, Q, P, Q, E, n) == 0){
	    printf("found factor during update of Q in Z7\n");
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
	if(ell_point_is_on_curve(Q, E, n) == 0){
	    printf("#!# Q=[%d]P is not on E\n", u);
	    //	    ell_point_print(Q, E, n); printf("\n");
	    ret = ECM_ERROR;
	    break;
	}
	/* come back to plain (not Montgomery) residues */
	mpres_get_z(b, Q->x, n);
	mpres_get_z(c, Q->y, n);
#if DEBUG_TORSION >= 2
	gmp_printf("b:=%Zd; c:=%Zd;\n", b, c);
	printf("(s, t)[%d]:=", u);
	pt_print(E, Q, n);
	printf(";\n");
#endif
	if(cubic_to_quartic(fac, n->orig_modulus, d, ky0, b, c, 
			    A2, A1div2, x0, y0, cte) == 0){
	    printf("found factor in Z7 (cubic_to_quartic)\n");
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
	/* (d, ky0) is a point on y^2 = x^4-18*x^3+13*x^2-28*x+4 */
	/* d:=x; */
	/* x0:=-2*d; */
	mpz_mul_si(kx0, d, -2);
	mpz_mod(kx0, kx0, n->orig_modulus);
	/* y0:=d*y/2; */
	mpz_mul(ky0, ky0, d);
	mpz_mod(ky0, ky0, n->orig_modulus);
	mod_div_2(ky0, n->orig_modulus);
	/* c:=d^2-d; */
	mpz_mul(c, d, d);
	mpz_sub(c, c, d);
	mpz_mod(c, c, n->orig_modulus);
	/* b:=c*d; */
	mpz_mul(b, c, d);
	mpz_mod(b, b, n->orig_modulus);
	/* to short Weierstrass form */
	kubert_to_weierstrass(A, B, X, Y, b, c, kx0, ky0, n->orig_modulus);
	if(check_weierstrass(A, B, X, Y, tmp, x0, n->orig_modulus) == 0){
	    ret = ECM_ERROR;
            break;
	}
	ell_curve_init(tE[nc], ECM_EC_TYPE_WEIERSTRASS, ECM_LAW_AFFINE,n);
        mpz_set(tE[nc]->a4, A);
        mpz_set(tE[nc]->a6, B);
	ell_point_init(tP[nc], tE[nc], n);
	mpz_set(tP[nc]->x, X);
	mpz_set(tP[nc]->y, Y);
#if DEBUG_TORSION >= 2
	gmp_printf("E[%d]:=[%Zd];\n", nc, tE[nc]->a4);
	gmp_printf("P[%d]:=[%Zd, %Zd, %Zd];\n", 
		   nc, tP[nc]->x, tP[nc]->y, tP[nc]->z);
#endif
	nc++;
	if(nc >= nE)
	    break;
    }
    mpz_clear(A2);
    mpz_clear(A1div2);
    mpz_clear(x0);
    mpz_clear(y0);
    mpz_clear(cte);
    mpz_clear(A);
    mpz_clear(B);
    mpz_clear(X);
    mpz_clear(Y);
    ell_point_clear(P, E, n);
    ell_point_clear(Q, E, n);
    ell_curve_clear(E, n);
    mpz_clear(d);
    mpz_clear(c);
    mpz_clear(b);
    mpz_clear(kx0);
    mpz_clear(ky0);
    mpres_clear(tmp, n);
    return ret;
}

/* 
   SIDE EFFECT: tE[0..nE[ and tP[0..nE[ receive a curve of torsion Z9
                and a point on it using parameters [umin..umax[.
   OUTPUT: ECM_NO_FACTOR_FOUND or ECM_FACTOR_FOUND_STEP1 if a factor is found.
   tE[i], tP[i] are built in raw modular form, not Montgomery form. 
   REM: we assume gcd(n, 6).
*/
int
build_curves_with_torsion_Z9(mpz_t fac, mpmod_t n, ell_curve_t *tE, 
			     ell_point_t *tP, int umin, int umax, int nE)
{
    int u, ret = ECM_NO_FACTOR_FOUND, nc = 0;
    mpz_t A2, A1div2, x0, y0, cte, d, c, b, kx0, ky0, A, B, X, Y, f;
    mpres_t tmp;
    ell_curve_t E;
    ell_point_t P, Q;

    mpz_init(A2);
    mpz_init(A1div2);
    mpz_init(cte);
    mpz_init(x0);
    mpz_init(y0);
    mpz_init(A);
    mpz_init(B);
    mpz_init(X);
    mpz_init(Y);
    /* Eaux = [-9, 9] */
    /* Paux = [1, 1, 1] */
    /* Y^2 = X^4-24*X-36 */
    mpres_init(tmp, n);
    build_curves_with_torsion_aux(E, P, A2, A1div2, x0, y0, cte,
				  "-9", "9", "1", "1", 
				  "0", "3", 
				  "2", "3", "0",
				  n, tmp);
    mpz_init(f);
    mpz_init(d);
    mpz_init(c);
    mpz_init(b);
    mpz_init(kx0);
    mpz_init(ky0);
    ell_point_init(Q, E, n);
    mpz_set_si(d, umin-1);
    if(ell_point_mul(fac, Q, d, P, E, n) == 0){
	printf("found factor during init of Q in Z9\n");
	ret = ECM_FACTOR_FOUND_STEP1;
    }
    for(u = umin; (ret != ECM_FACTOR_FOUND_STEP1) && u < umax; u++){
	/* update Q */
	if(ell_point_add(fac, Q, P, Q, E, n) == 0){
	    printf("found factor during update of Q in Z9\n");
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
#if DEBUG_TORSION >= 2
	printf("(s, t)[%d]:=", u);
	pt_print(E, Q, n);
	printf(";\n");
#endif
	if(ell_point_is_on_curve(Q, E, n) == 0){
	    printf("#!# Q=[%d]P is not on E\n", u);
	    ret = ECM_ERROR;
	    break;
	}
	mpres_get_z(b, Q->x, n);
	mpres_get_z(c, Q->y, n);
	if(cubic_to_quartic(fac, n->orig_modulus, f, ky0, b, c, 
			    A2, A1div2, x0, y0, cte) == 0){
	    printf("found factor in Z9 (cubic_2_quartic)\n");
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
	/* f:=x; */
	/* d:=f*(f-1)+1; */
	mpz_sub_si(d, f, 1);
	mpz_mul(d, d, f);
	mpz_add_si(d, d, 1);
	mpz_mod(d, d, n->orig_modulus);
	/* c:=f*(d-1); */
	mpz_sub_si(c, d, 1);
	mpz_mul(c, c, f);
	mpz_mod(c, c, n->orig_modulus);
	/* kx0:=(2*f-1)*f^2; */
	/** b <- f^2 **/
	mpz_mul(b, f, f);
	mpz_mod(b, b, n->orig_modulus);
	mpz_mul_si(kx0, f, 2);
	mpz_sub_si(kx0, kx0, 1);
	mpz_mul(kx0, kx0, b);
	mpz_mod(kx0, kx0, n->orig_modulus);
	/* ky0:=y*f^4/2; */
	/** b <- b^2 = f^4 **/
	mpz_mul(b, b, b);
	mpz_mul(ky0, ky0, b);
	mpz_mod(ky0, ky0, n->orig_modulus);
	mod_div_2(ky0, n->orig_modulus);
	/* b:=c*d; */
	mpz_mul(b, c, d);
	mpz_mod(b, b, n->orig_modulus);
#if DEBUG_TORSION >= 2
	gmp_printf("f=%Zd d=%Zd c=%Zd b=%Zd\n", f, d, c, b);
	gmp_printf("kx0=%Zd ky0=%Zd\n", kx0, ky0);
#endif
	/* to short Weierstrass form */
	kubert_to_weierstrass(A, B, X, Y, b, c, kx0, ky0, n->orig_modulus);
	if(check_weierstrass(A, B, X, Y, tmp, x0, n->orig_modulus) == 0){
            ret = ECM_ERROR;
            break;
        }
	ell_curve_init(tE[nc], ECM_EC_TYPE_WEIERSTRASS, ECM_LAW_AFFINE, n);
	mpz_set(tE[nc]->a4, A);
        mpz_set(tE[nc]->a6, B);
	ell_point_init(tP[nc], tE[nc], n);
        mpz_set(tP[nc]->x, X);
        mpz_set(tP[nc]->y, Y);
	nc++;
	if(nc >= nE)
	    break;
    }
    mpz_clear(A);
    mpz_clear(B);
    mpz_clear(X);
    mpz_clear(Y);
    mpz_clear(A2);
    mpz_clear(A1div2);
    mpz_clear(x0);
    mpz_clear(y0);
    mpz_clear(cte);
    ell_point_clear(P, E, n);
    ell_point_clear(Q, E, n);
    mpz_clear(f);
    mpz_clear(d);
    mpz_clear(c);
    mpz_clear(b);
    mpz_clear(kx0);
    mpz_clear(ky0);
    mpres_clear(tmp, n);
    return ret;
}

int
build_curves_with_torsion_Z10(mpz_t fac, mpmod_t n, ell_curve_t *tE, 
			      ell_point_t *tP, int umin, int umax, int nE)
{
    int u, ret = ECM_NO_FACTOR_FOUND, nc = 0;
    mpz_t A2, A1div2, x0, y0, cte, d, c, b, kx0, ky0, A, B, X, Y;
    mpz_t f;
    mpres_t tmp;
    ell_curve_t E;
    ell_point_t P, Q;

    mpz_init(A2);
    mpz_init(A1div2);
    mpz_init(cte);
    mpz_init(x0);
    mpz_init(y0);
    mpz_init(A);
    mpz_init(B);
    mpz_init(X);
    mpz_init(Y);
    /* Eaux = [2/3, -53/108] */
    /* Paux = [2/3, 1/2, 1] */
    /* Y^2 = X^4-4*X^2-4*X-4 */
    mpres_init(tmp, n);
    build_curves_with_torsion_aux(E, P, A2, A1div2, x0, y0, cte,
				  "2/3", "-53/108", "2/3", "1/2",
				  "2/3", "-1/2", 
				  "0", "1", "-2",
				  n, tmp);
    mpz_init(f);
    mpz_init(d);
    mpz_init(c);
    mpz_init(b);
    mpz_init(kx0);
    mpz_init(ky0);
    ell_point_init(Q, E, n);
    for(u = umin; u < umax; u++){
	if(forbidden("Z10", u))
	    continue;
	/* update Qaux */
	mpz_set_si(d, u);
	if(ell_point_mul(fac, Q, d, P, E, n) == 0){
	    printf("found factor during update of Q in Z10\n");
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
#if DEBUG_TORSION >= 2
	printf("(s, t)[%d]:=", u);
	pt_print(E, Q, n);
	printf(";\n");
#endif
	if(ell_point_is_on_curve(Q, E, n) == 0){
	    printf("#!# Q=[%d]P is not on E\n", u);
	    ret = ECM_ERROR;
	    break;
	}
	mpres_get_z(b, Q->x, n);
	mpres_get_z(c, Q->y, n);
	if(cubic_to_quartic(fac, n->orig_modulus, f, ky0, b, c, 
			    A2, A1div2, x0, y0, cte) == 0){
	    printf("found factor in Z10 (cubic_2_quartic)\n");
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
	/* f:=x; */
	/* d:=f^2/(f-(f-1)^2) = -f^2/(f^2-3*f+1) */
	/** b <- -f^2 **/
	mpz_mul(b, f, f);
	mpz_neg(b, b);
	mpz_mod(b, b, n->orig_modulus);
	/* c = f^2-3*f+1 = f*(f-3)+1 */
	mpz_sub_si(c, f, 3);
	mpz_mul(c, c, f);
	mpz_add_si(c, c, 1);
	mpz_mod(c, c, n->orig_modulus);
	if(mod_from_rat2(d, b, c, n->orig_modulus) == 0){
	    printf("inverse found in Z10 (d)\n");
	    mpz_set(fac, d);
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
	/* c:=f*(d-1); */
	mpz_sub_si(c, d, 1);
	mpz_mul(c, c, f);
	mpz_mod(c, c, n->orig_modulus);
	/* ky0:=y*f^4/(f^2-3*f+1)^2/2; = num/den */
	/* it seems that ky0 = y*d^2/2 */
	mpz_mul(b, ky0, d);
	mpz_mul(b, b, d);
	mpz_mod(b, b, n->orig_modulus);
	mpz_set_si(fac, 2);
	mod_from_rat2(ky0, b, fac, n->orig_modulus);
	/* kx0:=-f*d; */
	mpz_mul(kx0, f, d);
	mpz_neg(kx0, kx0);
	mpz_mod(kx0, kx0, n->orig_modulus);
	/* b:=c*d; */
	mpz_mul(b, c, d);
	mpz_mod(b, b, n->orig_modulus);
#if DEBUG_TORSION >= 2
	gmp_printf("f:=%Zd; d:=%Zd; c:=%Zd; b:=%Zd;\n", f, d, c, b);
	gmp_printf("kx0:=%Zd; ky0:=%Zd;\n", kx0, ky0);
#endif
	/* to short Weierstrass form */
	kubert_to_weierstrass(A, B, X, Y, b, c, kx0, ky0, n->orig_modulus);
	if(check_weierstrass(A, B, X, Y, tmp, x0, n->orig_modulus) == 0){
            ret = ECM_ERROR;
            break;
        }
	ell_curve_init(tE[nc], ECM_EC_TYPE_WEIERSTRASS, ECM_LAW_AFFINE, n);
	mpz_set(tE[nc]->a4, A);
        mpz_set(tE[nc]->a6, B);
	ell_point_init(tP[nc], tE[nc], n);
        mpz_set(tP[nc]->x, X);
        mpz_set(tP[nc]->y, Y);
	nc++;
	if(nc >= nE)
	    break;
    }
#if DEBUG_TORSION >= 2
    if(ret != ECM_ERROR && nc > 0){
	printf("Curves built\n");
	pt_many_print(tE, tP, nE, n);
    }
#endif
    mpz_clear(A);
    mpz_clear(B);
    mpz_clear(X);
    mpz_clear(Y);
    mpz_clear(A2);
    mpz_clear(A1div2);
    mpz_clear(x0);
    mpz_clear(y0);
    mpz_clear(cte);
    ell_point_clear(P, E, n);
    ell_point_clear(Q, E, n);
    ell_curve_clear(E, n);
    mpres_clear(tmp, n);
    mpz_clear(d);
    mpz_clear(c);
    mpz_clear(b);
    mpz_clear(kx0);
    mpz_clear(ky0);
    mpz_clear(f);
    return ret;
}

/* Warning: b and a have the Montgomery meaning in this function. 
   All tE[i] will be in Montgomery form: B*Y^2 = X^3 + A * X^2 + X.
*/
int
build_curves_with_torsion_Z2xZ8(mpz_t fac, mpmod_t n, 
				ell_curve_t *tE, ell_point_t *tP,
				int umin, int umax, int nE)
{
    int u, nc = 0, ret = ECM_NO_FACTOR_FOUND;
    mpz_t tmp, a, b, alpha, beta, c, d, kx0, ky0, wx0, mb;
    mpres_t tmp2;
    ell_curve_t E;
    ell_point_t P, Q;

    mpz_init(alpha);
    mpz_init(beta);
    mpz_init(tmp);
    mpz_init(a);
    mpz_init(b);
    mpz_init(c);
    mpz_init(d);
    mpz_init(kx0);
    mpz_init(ky0);
    mpz_init(wx0);
    mpz_init(mb);

    /* Eaux = [-8, -32] */
    /* Paux = [12, 40, 1] */
    mpres_init(tmp2, n);
    mpz_set_str(fac, "-8", 10); 
    mpres_set_z(tmp2, fac, n);
    ell_curve_init_set(E, ECM_EC_TYPE_WEIERSTRASS, ECM_LAW_AFFINE, tmp2, n);
    ell_point_init(P, E, n);
    mpz_set_str(fac, "12", 10); 
    mpres_set_z(P->x, fac, n);
    mpz_set_str(fac, "40", 10);
    mpres_set_z(P->y, fac, n);
    mpz_set_ui(P->z, 1);

    ell_point_init(Q, E, n);
    mpz_set_si(d, umin-1);
    if(ell_point_mul(fac, Q, d, P, E, n) == 0){
	printf("found factor during init of Q in Z2xZ8\n");
	ret = ECM_FACTOR_FOUND_STEP1;
    }
    for(u = umin; (ret != ECM_FACTOR_FOUND_STEP1) && u < umax; u++){
	/* update Q */
	if(ell_point_add(fac, Q, P, Q, E, n) == 0){
	    printf("found factor during update of Q in Z2xZ8\n");
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
#if DEBUG_TORSION >= 2
	printf("(s, t)[%d]:=", u);
	pt_print(E, Q, n);
	printf(";\n");
#endif
	/* beta <- (y+25)/(x-9) */
	mpres_get_z(a, Q->x, n);
	mpres_get_z(b, Q->y, n);
	mpz_mod(wx0, a, n->orig_modulus);
	mpz_sub_si(a, a, 9);
	mpz_mod(a, a, n->orig_modulus);
	mpz_add_si(b, b, 25);
	mpz_mod(b, b, n->orig_modulus);
	if(mod_from_rat2(beta, b, a, n->orig_modulus) == 0){
            printf("found factor in Z2xZ8 (beta)\n");
	    mpz_set(fac, beta);
	    ret = ECM_FACTOR_FOUND_STEP1;
            break;
	}
	mpz_add_si(tmp, beta, 1);
	mpz_mod(tmp, tmp, n->orig_modulus);
	/* alpha <- 1/(beta+1) */
	if(mpz_invert(alpha, tmp, n->orig_modulus) == 0){
            printf("found factor in Z2xZ8 (alpha)\n");
	    mpz_gcd(fac, tmp, n->orig_modulus);
	    ret = ECM_FACTOR_FOUND_STEP1;
            break;
	}
	/** d <- 8*alpha^2-1; 
	    d = -(beta^2+2*beta-7)/(beta+1)^2
	 **/
	mpz_mul(d, alpha, alpha);
	mpz_mul_si(d, d, 8);
	mpz_sub_si(d, d, 1);
	mpz_mod(d, d, n->orig_modulus);
	/* d:=2*alpha*(4*alpha+1)/d; */
	mpz_mul_si(c, alpha, 4);
	mpz_add_si(c, c, 1);
	mpz_mul(c, c, alpha);
	mpz_mul_si(c, c, 2);
	mpz_mod(c, c, n->orig_modulus);
	if(mod_from_rat2(fac, c, d, n->orig_modulus) == 0){
	    // the only possibility is d = 0 mod p or 8*alpha^2-1 = 0 mod  p
            printf("found factor in Z2xZ8 (d)\n");
	    ret = ECM_FACTOR_FOUND_STEP1;
            break;
	}
	mpz_set(d, fac);
	/* c:=(2*d-1)*(d-1)/d;*/
	mpz_sub_si(fac, d, 1);
	/** kx0 <- 2*d-1 **/
	mpz_mul_si(kx0, d, 2);
	mpz_sub_si(kx0, kx0, 1);
	mpz_mul(fac, fac, kx0);
	mpz_mod(fac, fac, n->orig_modulus);
        if(mod_from_rat2(c, fac, d, n->orig_modulus) == 0){
	    // this is possible only if d = 0 mod p or 
	    // 2*alpha*(4*alpha+1) = 0 mod p
            printf("found factor in Z2xZ8 (d2)\n");
	    mpz_set(fac, c);
            ret = ECM_FACTOR_FOUND_STEP1;
            break;
        }
	/* b = c*d */
	mpz_mul(b, c, d);
	mpz_mod(b, b, n->orig_modulus);
	/* kx0:=-(2*d-1)/4;*/
	mod_div_2(kx0, n->orig_modulus);
	mod_div_2(kx0, n->orig_modulus);
	mpz_mul_si(kx0, kx0, -1);
	mpz_mod(kx0, kx0, n->orig_modulus);
	/* ky0:=(c/8)*(-beta^2+2*uP[1]+9); */
	mpz_mul(fac, beta, beta);
	mpz_set(a, wx0);
	mpz_sub(fac, a, fac);
	mpz_add(fac, fac, a);
	mpz_add_si(fac, fac, 9);
	mpz_mul(fac, fac, c);
	mpz_mod(fac, fac, n->orig_modulus);
	mod_div_2(fac, n->orig_modulus);
	mod_div_2(fac, n->orig_modulus);
	mod_div_2(fac, n->orig_modulus);
	/* ky0:=ky0/(beta^2+2*beta-7); */
	mpz_add_si(tmp, beta, 2);
	mpz_mul(tmp, tmp, beta);
	mpz_sub_si(tmp, tmp, 7);
	mpz_mod(tmp, tmp, n->orig_modulus);
	/* as proven above, we cannot have tmp non invertible at that point */
	mod_from_rat2(ky0, fac, tmp, n->orig_modulus);
	KW2W246(fac, a, NULL, b, c, n->orig_modulus, 0);
#if DEBUG_TORSION >= 2
	gmp_printf("kwx0:=%Zd;\n", kx0);
	gmp_printf("kwy0:=%Zd;\n", ky0);
	printf("(kwy0^2-(kwx0^3+a2*kwx0^2+a4*kwx0+a6)) mod N;\n");
#endif
	/* wx0:=kx0+a2/3; */
        mpz_set_si(tmp, 3);
	mod_from_rat2(wx0, fac, tmp, n->orig_modulus);
	mpz_add(wx0, wx0, kx0);
	mpz_mod(wx0, wx0, n->orig_modulus);
	/* mb:=-1/(d-1)^2; */
	mpz_sub_si(tmp, d, 1);
	mpz_mul(tmp, tmp, tmp);
	mpz_mod(tmp, tmp, n->orig_modulus);
	mpz_neg(tmp, tmp);
	if(mpz_invert(mb, tmp, n->orig_modulus) == 0){
	    printf("found factor in Z2xZ8 (mb)\n");
	    mpz_gcd(fac, tmp, n->orig_modulus);
            ret = ECM_FACTOR_FOUND_STEP1;
            break;
	}
	/* ma:=-1/4*(8*d^4-16*d^3+16*d^2-8*d+1)/(d-1)^2/d^2; 
	     :=mb*(8*d^4-16*d^3+16*d^2-8*d+1)/(4*d^2)
	 */
	mpz_set_si(fac, 8);         /* num */
	mpz_mul(fac, fac, d); mpz_add_si(fac, fac, -16);
	mpz_mul(fac, fac, d); mpz_add_si(fac, fac, 16);
	mpz_mul(fac, fac, d); mpz_add_si(fac, fac, -8);
	mpz_mul(fac, fac, d); mpz_add_si(fac, fac, 1);
#if 0
	mpz_sub_si(tmp, d, 1);    /* den */
	mpz_mul(tmp, tmp, d);
	mpz_mul(tmp, tmp, tmp);
	mpz_mul_si(tmp, tmp, -4);
	mpz_mod(tmp, tmp, n->orig_modulus);
#else
	mpz_mul(fac, fac, mb);
	/* one day, we could save 1/d computation again */
	mpz_mul(tmp, d, d);
	mpz_mul_si(tmp, tmp, 4);
#endif
	/* to Montgomery form */
	ell_curve_init(tE[nc], ECM_EC_TYPE_MONTGOMERY, ECM_LAW_HOMOGENEOUS,n);
	ell_point_init(tP[nc], tE[nc], n);
	/* this cannot yield a factor, since d is invertible at that point */
	mod_from_rat2(tE[nc]->a2, fac, tmp, n->orig_modulus);
	/* not really needed, but useful for debug */
	mpz_set_ui(tE[nc]->a4, 1);
	mpz_set_ui(tE[nc]->a6, 0);
	/* mx:=mb*wx0-ma/3; */
	mpz_mul(fac, mb, wx0);
        mpz_set_si(tmp, 3);
        mod_from_rat2(tP[nc]->x, tE[nc]->a2, tmp, n->orig_modulus);
	mpz_sub(tP[nc]->x, fac, tP[nc]->x);
	mpz_mod(tP[nc]->x, tP[nc]->x, n->orig_modulus);
	/* my:=mb*ky0; */
#if DEBUG_TORSION >= 2
	gmp_printf("N:=%Zd;\n", n->orig_modulus);
	gmp_printf("ma:=%Zd;\n", tE[nc]->a2);
	gmp_printf("mb:=%Zd;\n", mb);
	gmp_printf("kx0:=%Zd;\n", kx0);
	gmp_printf("ky0:=%Zd;\n", ky0);
	gmp_printf("mx0:=%Zd;\n", tP[nc]->x);
	mpz_mul(tmp, mb, ky0);
	mpz_mod(tmp, tmp, n->orig_modulus);
	gmp_printf("my0:=%Zd;\n", tmp);
	printf("chk:=(mb*my0^2-mx0^3-ma*mx0^2-mx0) mod N;\n");
#endif
	nc++;
	if(nc >= nE)
	    break;
    }
#if DEBUG_TORSION >= 2
    printf("Curves built\n");
    pt_many_print(tE, tP, nE, n);
#endif
    ell_point_clear(P, E, n);
    ell_point_clear(Q, E, n);
    ell_curve_clear(E, n);
    mpz_clear(mb);
    mpz_clear(tmp);
    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(c);
    mpz_clear(d);
    mpz_clear(alpha);
    mpz_clear(beta);
    mpz_init(kx0);
    mpz_init(ky0);
    mpz_init(wx0);
    mpres_clear(tmp2, n);
    return ret;
}

/* Z3xZ3 over Q(sqrt(-3)). Interesting if we know that p | N is s.t.
   p = 1 mod 3.
   Source: Dujella and Najman, arxiv:1201.0266v1 
   A more simpler and more efficient stuff, using Hessian form. */
int
build_curves_with_torsion_Z3xZ3(mpz_t f, mpmod_t n, 
				ell_curve_t *tE, ell_point_t *tP,
				int umin, int umax, int nE)
{
    int u, nc = 0, ret = ECM_NO_FACTOR_FOUND;
    mpz_t u0, v0, D, num, den;

    mpz_init(u0);
    mpz_init(num);
    mpz_init(den);
    mpz_init(D);
    mpz_init_set_si(v0, umin-1); /* to prevent u0 = v0 */
    for(u = umin; u < umax; u++){
	if(forbidden("Z3xZ3", u))
	    continue;
	mpz_set_si(u0, u);
	/* D:=RatMod((u0^3+v0^3+1)/(3*u0*v0), N); */
	mpz_mul(num, u0, u0);
	mpz_mul(num, num, u0);
	mpz_mul(den, v0, v0);
	mpz_mul(den, den, v0);
	mpz_add(num, num, den);
	mpz_add_si(num, num, 1);
	mpz_mod(num, num, n->orig_modulus);
	if(mpz_sgn(num) == 0)
	    continue;
	
	mpz_mul(den, u0, v0);
	mpz_mul_si(den, den, 3);
	mpz_mod(den, den, n->orig_modulus);

	if(mod_from_rat2(D, num, den, n->orig_modulus) == 0){
	    printf("found factor in Z3xZ3 (D)\n");
            mpz_set(f, D);
            ret = ECM_FACTOR_FOUND_STEP1;
            break;
	}
	mpz_mul(num, D, D);
	mpz_mul(num, num, D);
	mpz_mod(num, num, n->orig_modulus);
	if(mpz_cmp_ui(num, 1) == 0){
	    printf("D^3=1 => singluar curve\n");
	    ret = ECM_ERROR;
	    break;
	}
	ell_curve_init_set(tE[nc],ECM_EC_TYPE_HESSIAN,ECM_LAW_HOMOGENEOUS,D,n);
	ell_point_init(tP[nc], tE[nc], n);
	mpz_set(tP[nc]->x, u0);
	mpz_set(tP[nc]->y, v0);
	mpz_set_ui(tP[nc]->z, 1);
	nc++;
	if(nc >= nE)
	    break;
    }
    mpz_clear(u0);
    mpz_clear(v0);
    mpz_clear(D);
    mpz_clear(num);
    mpz_clear(den);
    return ret;
}

/* For a small price, add a 2-torsion point, also over Q(sqrt(-3)). */
int
build_curves_with_torsion_Z3xZ6(mpz_t f, mpmod_t n, 
				ell_curve_t *tE, ell_point_t *tP,
				int umin, int umax, int nE)
{
    int u, nc = 0, ret = ECM_NO_FACTOR_FOUND;
    ell_curve_t E;
    ell_point_t P, Q;
    mpres_t tmp, num, den, tk, sk;
    mpz_t t;

    mpz_init(t);
    mpz_init(num);
    mpz_init(den);
    mpz_init(tk);
    mpz_init(sk);
    /* Eaux:=EllipticCurve([0, -4]); */
    /* Paux:=Eaux![2, 2, 1]; */
    mpres_init(tmp, n);
    mpres_set_ui(tmp, 0, n);
    ell_curve_init_set(E, ECM_EC_TYPE_WEIERSTRASS, ECM_LAW_AFFINE, tmp, n);
    ell_point_init(P, E, n);
    mpz_set_str(f, "2", 10);
    mpres_set_z(P->x, f, n);
    mpz_set_str(f, "2", 10);
    mpres_set_z(P->y, f, n);
    mpz_set_ui(P->z, 1);

    ell_point_init(Q, E, n);
    for(u = umin; u < umax; u++){
	/* update Qaux */
	mpz_set_si(f, u);
	if(ell_point_mul(f, Q, f, P, E, n) == 0){
	    printf("found factor in Z3xZ6 (update of Q)\n");
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
#if DEBUG_TORSION >= 2
	printf("(s, t)[%d]:=", u);
	pt_print(E, Q, n);
	printf(";\n");
#endif
	mpres_get_z(tk, Q->x, n);
	mpres_get_z(sk, Q->y, n);
#if 0 /* useless in affine form? */
	mpres_get_z(t, Q->z, n);
	if(mpz_invert(f, t, n->orig_modulus) == 0){
	    printf("found factor in Z3xZ6 (normalization)\n");
	    mpz_gcd(f, t, n->orig_modulus);
	    break;
	}
	mpz_mul(tk, tk, f);
	mpz_mod(tk, tk, n->orig_modulus);
	mpz_mul(sk, sk, f);
	mpz_mod(sk, sk, n->orig_modulus);
#endif
	/* t:=RatMod(-tk/2, N); */
	mpz_mul_si(t, tk, -1);
	mod_div_2(t, n->orig_modulus);
	/* D:=RatMod((2*t^3+1)/3/t^2, N); */
	mpz_mul(den, t, t);
	mpz_mod(den, den, n->orig_modulus);
	mpz_mul(num, den, t);
	mpz_mul_si(num, num, 2);
	mpz_add_si(num, num, 1);
	mpz_mod(num, num, n->orig_modulus);
	mpz_mul_si(den, den, 3);
	mpz_mod(den, den, n->orig_modulus);
	ell_curve_init(tE[nc], ECM_EC_TYPE_HESSIAN, ECM_LAW_HOMOGENEOUS, n);
	ell_point_init(tP[nc], tE[nc], n);
	if(mod_from_rat2(tE[nc]->a4, num, den, n->orig_modulus) == 0){
	    /* only if t = 0, which seems hard */
            printf("found factor in Z3xZ6 (D)\n");
            mpz_set(f, tE[nc]->a4);
            ret = ECM_FACTOR_FOUND_STEP1;
            break;
        }
#if DEBUG_TORSION >= 1
	gmp_printf("D%d:=%Zd;\n", nc, tE[nc]->a4);
#endif
	/* u0:=RatMod(sk/tk, N); 
	   if tk was not invertible, it would have been caught before
	 */
	mod_from_rat2(tP[nc]->x, sk, tk, n->orig_modulus);
	/* v0:=-1; */
	mpz_sub_si(tP[nc]->y, n->orig_modulus, 1);
        mpz_set_ui(tP[nc]->z, 1);
	nc++;
	if(nc >= nE)
	    break;
    }
    mpz_clear(t);
    mpz_clear(num);
    mpz_clear(den);
    mpz_clear(sk);
    mpz_clear(tk);
    mpres_clear(tmp, n);
    return ret;
}

/* JKL: K = Q(sqrt(-3), sqrt(8*t^3+1), t in Q, t != 0, 1, -1/2;
   mu = (2*t^3+1)/(3*t^2) => parameter for Hessian form.
   Tors(E) = Z/6xZ/6.
   A "specified" point is (0:-1:1), but does it have infinite order?
   Also: twisted Hessian is a*X^3+Y^3+Z^3=d*X*Y*Z, d/a=3*mu.
   See JKL-ECM in ANTS-XII.
 */

/* Original source is Brier + Clavier.
   We can build curves in Montgomery form directly... 
   Useful if one knows that all p | n are 1 mod 4 (Cunningham, etc.).
*/
int
build_curves_with_torsion_Z4xZ4(mpz_t f, mpmod_t n, ell_curve_t *tE,
				ell_point_t *tP,
				int smin, int smax, int nE)
{
    mpz_t tau, lambda, nu2, tmp, b, x0;
    int nu, nc = 0, ret = ECM_NO_FACTOR_FOUND;

    mpz_init(tau);
    mpz_init(lambda);
    mpz_init(nu2);
    mpz_init(tmp);
    mpz_init(b);
    mpz_init(x0);
    for(nu = smin; nu < smax; nu++){
	mpz_set_si(nu2, nu*nu);
	/* tau:=(nu^2+3)/2/nu; */
	mpz_add_si(lambda, nu2, 3);
	mpz_set_si(tmp, 2*nu);
	if(mod_from_rat2(tau, lambda, tmp, n->orig_modulus) == 0){
            printf("Factor found during init of Z4xZ4 (tau)\n");
	    mpz_set(f, tau);
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
	/* lambda:=8*nu^3; */
	mpz_mul_si(lambda, nu2, 8*nu);
	mpz_mod(lambda, lambda, n->orig_modulus);
	/* A:=-27*lambda^4*(tau^8+14*tau^4+1); */
	/* B:=54*lambda^6*(tau^12-33*tau^8-33*tau^4+1); */
	/* x0:=3*(3*nu^12+34*nu^10+117*nu^8+316*nu^6+1053*nu^4+2754*nu^2+2187); */
	/* y0:=27*(nu^2-3)*(nu^2+1)*(nu^2+9)*(nu^6+5*nu^4+15*nu^2+27)^2; */
	/* P = (x0, y0) is a point on Y^2 = X^3+A*X+B */

	/* Montgomery form: there are several possible mb */
	/* mb:=1/(9*lambda^2*(tau^4-1));
	   lambda is invertible iff nu is;
	   tau^4-1 = (tau-1)(tau+1)(tau^2+1)
	*/
	mpz_powm_ui(x0, tau, 4, n->orig_modulus);
	mpz_sub_si(x0, x0, 1);
	mpz_mod(x0, x0, n->orig_modulus);
	mpz_mul(tmp, x0, lambda);
	mpz_mul(tmp, tmp, lambda);
	mpz_mul_si(tmp, tmp, 9);
	if(mpz_invert(b, tmp, n->orig_modulus) == 0){
	    printf("Factor found during init of Z4xZ4 (mb)\n");
	    mpz_gcd(f, tmp, n->orig_modulus);
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
	/* ma:=-2*(tau^4+1)/(tau^4-1); at this point: invertible! */
	mpz_add_si(tmp, x0, 2);
	mpz_mul_si(tmp, tmp, -2);
	mpz_mod(tmp, tmp, n->orig_modulus);
        /* to Montgomery form */
        ell_curve_init(tE[nc], ECM_EC_TYPE_MONTGOMERY, ECM_LAW_HOMOGENEOUS, n);
        ell_point_init(tP[nc], tE[nc], n);
	mod_from_rat2(tE[nc]->a4, tmp, x0, n->orig_modulus);
	/* now compute real x0 */
	/* x0:=3*(3*nu^12+34*nu^10+117*nu^8+316*nu^6+1053*nu^4+2754*nu^2+2187); */
	mpz_set_si(x0, 3);
	mpz_mul(x0, x0, nu2); mpz_add_si(x0, x0, 34);
	mpz_mul(x0, x0, nu2); mpz_add_si(x0, x0, 117);
	mpz_mul(x0, x0, nu2); mpz_add_si(x0, x0, 316);
	mpz_mul(x0, x0, nu2); mpz_add_si(x0, x0, 1053);
	mpz_mul(x0, x0, nu2); mpz_add_si(x0, x0, 2754);
	mpz_mul(x0, x0, nu2); mpz_add_si(x0, x0, 2187);
	mpz_mul_si(x0, x0, 3);
	mpz_mod(x0, x0, n->orig_modulus);
#if DEBUG_TORSION >= 2
	gmp_printf("N:=%Zd;\n", n);
	printf("nu:=%d;\n", nu);
	gmp_printf("tau:=%Zd;\n", tau);
	gmp_printf("lambda:=%Zd;\n", lambda);
	gmp_printf("a:=%Zd;\n", tE[nc]->a4);
	gmp_printf("x0:=%Zd;\n", x0);
#endif
	/* x:=b*x0-a/3; not needed: y:=b*y0 */
	mpz_set_si(tmp, 3);
	mod_from_rat2(tP[nc]->x, tE[nc]->a4, tmp, n->orig_modulus);
	mpz_mul(b, b, x0);
	mpz_mod(b, b, n->orig_modulus);
	mpz_sub(tP[nc]->x, b, tP[nc]->x);
	mpz_mod(tP[nc]->x, tP[nc]->x, n->orig_modulus);
	nc++;
	if(nc >= nE)
	    break;
    }
    mpz_clear(tau);
    mpz_clear(lambda);
    mpz_clear(nu2);
    mpz_clear(tmp);
    mpz_clear(b);
    mpz_clear(x0);
    if(ret != ECM_FACTOR_FOUND_STEP1 && nc < nE){
	printf("Not enough curves generated\n");
	return ECM_ERROR;
    }
    return ret;
}

/* Assuming we can generate curves with given torsion using parameter s
   in interval [smin..smax[.
*/
int
build_curves_with_torsion(mpz_t f, mpmod_t n, ell_curve_t *tE, ell_point_t *tP,
			  char *torsion, int smin, int smax, int nE)
{
    int ret = 0;

    /* over Q: see Atkin-Morain, Math. Comp., 1993 */
    if(strcmp(torsion, "Z5") == 0)
	return build_curves_with_torsion_Z5(f, n, tE, tP, smin, smax, nE);
    else if(strcmp(torsion, "Z7") == 0)
	return build_curves_with_torsion_Z7(f, n, tE, tP, smin, smax, nE);
    else if(strcmp(torsion, "Z9") == 0)
	return build_curves_with_torsion_Z9(f, n, tE, tP, smin, smax, nE);
    else if(strcmp(torsion, "Z10") == 0)
	return build_curves_with_torsion_Z10(f, n, tE, tP, smin, smax, nE);
    else if(strcmp(torsion, "Z2xZ8") == 0)
	return build_curves_with_torsion_Z2xZ8(f, n, tE, tP, smin, smax, nE);
    /* no longer over Q */
    /** interesting when p = 1 mod 3 **/
    else if(strcmp(torsion, "Z3xZ3") == 0) /* over Q(sqrt(-3)) */
	return build_curves_with_torsion_Z3xZ3(f, n, tE, tP, smin, smax, nE);
    else if(strcmp(torsion, "Z3xZ6") == 0) /* over Q(sqrt(-3)) */
	return build_curves_with_torsion_Z3xZ6(f, n, tE, tP, smin, smax, nE);
    /** interesting when p = 1 mod 4 **/
    else if(strcmp(torsion, "Z4xZ4") == 0) /* over Q(sqrt(-1)) */
	return build_curves_with_torsion_Z4xZ4(f, n, tE, tP, smin, smax, nE);
    else{
	printf("Unknown torsion group: %s\n", torsion);
	ret = ECM_ERROR;
    }
    return ret;
}

/* E is a curve with given torsion and (x, y) a point on E mod n.
   OUTPUT: ECM_NO_FACTOR_FOUND if everything went ok
           ECM_FACTOR_FOUND_STEP1 in case a factor was found when building E.
   REM: E is defined over Z, not in mpres_t.
 */
int
build_curves_with_torsion2(mpz_t f, mpz_t n, ell_curve_t E, 
			   mpz_t x, mpz_t y, char *torsion, 
			   mpz_t sigma)
{
    ell_curve_t tE[1];
    ell_point_t tP[1];
    mpmod_t modulus;
    int ret, smin, smax;

    smin = (int)mpz_get_si(sigma);
    smax = smin+10;
    mpmod_init(modulus, n, ECM_MOD_DEFAULT);
    ret = build_curves_with_torsion(f, modulus, tE, tP, torsion, smin, smax,1);
    if(ret == ECM_NO_FACTOR_FOUND){
	E->type = tE[0]->type;
	E->law = tE[0]->law;
	mpz_set(E->a2, tE[0]->a2);
	mpz_set(E->a4, tE[0]->a4);
	mpz_set(E->a6, tE[0]->a6);
	mpz_set(x, tP[0]->x);
	mpz_set(y, tP[0]->y);
	ell_point_clear(tP[0], tE[0], modulus);
	ell_curve_clear(tE[0], modulus);
    }
    mpmod_clear(modulus);
    return ret;
}
