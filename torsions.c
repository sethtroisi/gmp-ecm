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
mod_from_rat(mpz_t r, mpq_t q, mpz_t N)
{
    return mod_from_rat2(r, mpq_numref(q), mpq_denref (q), N);
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

/* From a curve in long form y^2+a1*x*y+a3*y = x^3+a2*x^2+a4*x+a6
   to a medium Weiestrass form Y^2 = X^3 + A2 * X^2 + A4 * X + A6
   where X = x, Y = y+(a1*X+a3)/2
   A2 = 1/4*a1^2+a2
   A4 = 1/2*a1*a3+a4
   A6 = 1/4*a3^2+a6
*/
void
LongWeierstrassToMediumWeierstrass(mpz_t A2, mpz_t A4, mpz_t A6, 
				   mpz_t X, mpz_t Y,
				   mpz_t a1, mpz_t a3, mpz_t a2,
				   mpz_t a4, mpz_t a6, mpz_t x, mpz_t y,
				   mpz_t n)
{
    /** A4 <- a1/2 **/
    mpz_set(A4, a1);
    mod_div_2(A4, n);
    /** A2 <- A4^2+a2 **/
    mpz_mul(A2, A4, A4);
    mpz_add(A2, A2, a2);
    mpz_mod(A2, A2, n);
    /** finish A4 **/
    mpz_mul(A4, A4, a3);
    mpz_add(A4, A4, a4);
    mpz_mod(A4, A4, n);
    /** A6 **/
    mpz_mul(A6, a3, a3);
    mpz_mod(A6, A6, n);
    mod_div_2(A6, n);
    mod_div_2(A6, n);
    mpz_add(A6, A6, a6);
    mpz_mod(A6, A6, n);
    if(X != NULL && x != NULL)
	mpz_set(X, x);
    if(Y != NULL && y != NULL){
	mpz_mul(Y, a1, x);
	mpz_add(Y, Y, a3);
	mpz_mod(Y, Y, n);
	mod_div_2(Y, n);
	mpz_add(Y, Y, y);
	mpz_mod(Y, Y, n);
    }
#if DEBUG_TORSION >= 2
    gmp_printf("N:=%Zd;\n", n);
    gmp_printf("a1:=%Zd;a3:=%Zd;a2:=%Zd;a4:=%Zd;a6:=%Zd;\n",a1,a3,a2,a4,a6);
    gmp_printf("x:=%Zd; y:=%Zd;\n", x, y);
    printf("(y^2+a1*x*y+a3*y-(x^3+a2*x^2+a4*x+a6)) mod N;\n");
    gmp_printf("A2:=%Zd; A4:=%Zd; A6:=%Zd;\n", A2, A4, A6);
    gmp_printf("X:=%Zd; Y:=%Zd;\n", X, Y);
    printf("(Y^2-(X^3+A2*X^2+A4*X+A6)) mod N;\n");
#endif
}

/* Weierstrass (a2, a4, a6) to (A, B)
   A = (a4-1/3*a2^2)
   B = -1/3*a4*a2 + 2/27*a2^3 + a6
     = -1/3*a2*(a4-2/9*a2^2) + a6
   X = x+a2/3
   Y = y
   INPUT: if x0 == NULL, we have no point to translate
          if B == NULL, we do not need and we do not compute B
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

/* Eaux, Paux <- shortW */
void
LongWeierstrassToShortWeierstrass(mpz_t A, mpz_t B, mpz_t XX, mpz_t YY,
				  mpz_t a1, mpz_t a3, mpz_t a2,
				  mpz_t a4, mpz_t a6, mpz_t x, mpz_t y,
				  mpz_t n)
{
    mpz_t A2, A4, A6, X, Y;

    mpz_init(A2);
    mpz_init(A4);
    mpz_init(A6);
    mpz_init(X);
    mpz_init(Y);
    LongWeierstrassToMediumWeierstrass(A2,A4,A6,X,Y,a1,a3,a2,a4,a6,x,y,n);
    MediumWeierstrassToShortWeierstrass(A, B, XX, YY, A2, A4, A6, X, Y, n);
    mpz_clear(A2);
    mpz_clear(A4);
    mpz_clear(A6);
    mpz_clear(X);
    mpz_clear(Y);
#if DEBUG_TORSION >= 2
    gmp_printf("A:=%Zd; B:=%Zd;\n", A, B);
    gmp_printf("Px:=%Zd; Py:=%Zd;\n", XX, YY);
    printf("(Py^2-Px^3-A*Px-B) mod N;\n");
#endif
}

/* From a curve in Kubert form Y^2+(1-c)*X*Y-b*Y = X^3-b*X^2
   to a Weiestrass form y^2 = X^3 + a2 * X^2 + a4 * X + a6
   where y = (Y+(1-c)*X-b)/2
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

/* 
   Sends Kubert curve E(b, c): y^2+(1-c)*x*y-b*y = x^3-b*x^2
   with point (x0, y0) to short Weierstrass form:
   Y^2 = X^3 + A * X + B
*/
void
kubert_to_weierstrass(ell_curve_t E, mpz_t B, ell_point_t P, mpz_t b, mpz_t c, 
		      mpz_t x0, mpz_t y0, mpz_t n)
{
    mpz_t a2, a4, a6;
    int compute_a6 = (B != NULL);

    mpz_init(a2);
    mpz_init(a4);
    if(compute_a6)
	mpz_init(a6);
    /* ((y+(1-c)*x-b)/2)^2=x^3+a2*x^2+a4*x+a6 */
    KW2W246(a2, a4, a6, b, c, n, compute_a6);
    /* second conversion */
    MediumWeierstrassToShortWeierstrass(E->a4,B,P->x,P->y,a2,a4,a6,x0,y0,n);
    mpz_clear(a2);
    mpz_clear(a4);
    if(compute_a6)
	mpz_clear(a6);
}

/* Kubert: put b = c. */
int
build_curves_with_torsion_Z5(mpz_t f, mpmod_t n, 
			     ell_curve_t *tE, ell_point_t *tP,
			     int smin, int smax, int nE)
{
    int s, ret = ECM_NO_FACTOR_FOUND, nc = 0;
    mpz_t x0, y0, c, tmp;

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
	    /* factor found! */
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
	/* convert to short Weierstrass form */
	ell_curve_init(tE[nc], ECM_EC_TYPE_WEIERSTRASS, ECM_LAW_HOMOGENEOUS,n);
	ell_point_init(tP[nc], tE[nc], n);
	kubert_to_weierstrass(tE[nc], NULL, tP[nc], c, c, x0, y0, n->orig_modulus);
	nc++;
	if(nc >= nE)
	    break;
    }
    mpz_clear(x0);
    mpz_clear(y0);
    mpz_clear(c);
    mpz_clear(tmp);
    return ret;
}

/* INPUT: 
   T^2 = S^3 + A * S + B
   => quartic Y^2 = X^4 - 6 * A2 * X^2 + 4 * A1 * X + A0, with
   X = (T-A1/2)/(S-A2), Y = -X^2 + 2 * S + A2.
   => quartic y^2 = f(x) = a4*x^4+...+a0, where
   x = x0+y0/(X-cte), where cte = f'(x0)/4/y0
   y = Y/y0*(x-x0)^2 = Y*y0/(X-cte)^2
   OUTPUT: x, y
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
    /* 1st move */
    mpz_sub(x, t, A1div2);
    mpz_sub(y, s, A2);
    if(mod_from_rat2(X, x, y, n) == 0){
	mpz_set(f, X);
	ret = 0;
    }
    else{
	mpz_mul(Y, X, X);
	mpz_sub(Y, A2, Y);
	mpz_add(Y, Y, s);
	mpz_add(Y, Y, s);
	mpz_mod(Y, Y, n);
	/* 2nd move */
	mpz_sub(X, X, cte);
	mpz_mod(X, X, n);
	if(mpz_invert(f, X, n) == 0){
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

/* tE[i], tP[i] are built in raw modular form, not Montgomery form. */
int
build_curves_with_torsion_Z7(mpz_t f, mpmod_t n, 
			     ell_curve_t *tE, ell_point_t *tP,
			     int umin, int umax, int nE)
{
    int u, ret = ECM_NO_FACTOR_FOUND, nc = 0;
    mpz_t A2, A1div2, x0, y0, cte, d, c, b, kx0, ky0;
    mpres_t tmp;
    ell_curve_t E;
    ell_point_t P, Q;

    /* Eaux = "1295/48", "-1079/864" */
    /* Paux = "2185/12", "-2458" */
    /* E and P consist of Montgomery residues */
    mpres_init(tmp, n);
    mod_from_rat_str(f, "1295/48", n->orig_modulus);
    mpres_set_z(tmp, f, n);
    ell_curve_init_set(E, ECM_EC_TYPE_WEIERSTRASS, ECM_LAW_AFFINE, tmp, n);
    ell_point_init(P, E, n);
    mod_from_rat_str(f, "2185/12", n->orig_modulus);
    mpres_set_z(P->x, f, n);
    mpz_set_str(f, "-2458", 10);
    mpres_set_z(P->y, f, n);
    mpz_set_ui(P->z, 1);
#if DEBUG_TORSION >= 2
    printf("Paux:=");
    pt_print(E, P, n);
    printf(";\n");
#endif
    mpz_init(A2);
    mod_from_rat_str(A2, "1/12", n->orig_modulus);
    mpz_init_set_str(A1div2, "-1", 10);
    mpz_mod(A1div2, A1div2, n->orig_modulus);
    mpz_init_set_str(x0, "-1", 10);
    mpz_mod(x0, x0, n->orig_modulus);
    mpz_init_set_str(y0, "8", 10);
    mpz_mod(y0, y0, n->orig_modulus);
    mpz_init(cte);
    mod_from_rat_str(cte, "-7/2", n->orig_modulus);

    mpz_init(d);
    mpz_init(c);
    mpz_init(b);
    mpz_init(kx0);
    mpz_init(ky0);
    ell_point_init(Q, E, n);
    for(u = umin; u < umax; u++){
	/* update Qaux */
	mpz_set_ui(d, u);
	/* TODO: replace with ell_point_add, one of these days */
	if(ell_point_mul(Q, d, P, E, n) == 0){
	    printf("found factor during update of Q\n");
	    mpz_set(f, Q->x);
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
#if DEBUG_TORSION >= 2
	printf("(s, t)[%d]:=", u);
	pt_print(Q, n);
	printf(";\n");
#endif
	/* come back to plain (not Montgomery) residues */
	mpres_get_z(b, Q->x, n);
	mpres_get_z(c, Q->y, n);
#if 0
	mpres_get_z(d, Q->z, n);
	if(mpz_invert(f, d, n->orig_modulus) == 0){
	    printf("found factor in Z7 (normalization)\n");
	    mpz_gcd(f, d, n->orig_modulus);
	    break;
	}
	mpz_mul(b, b, f);
	mpz_mod(b, b, n->orig_modulus);
	mpz_mul(c, c, f);
	mpz_mod(c, c, n->orig_modulus);
#endif
	if(cubic_to_quartic(f, n->orig_modulus, d, ky0, b, c, 
			    A2, A1div2, x0, y0, cte) == 0){
	    printf("found factor in Z7 (cubic_2_quartic)\n");
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
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
	ell_curve_init(tE[nc], ECM_EC_TYPE_WEIERSTRASS, ECM_LAW_HOMOGENEOUS,n);
	ell_point_init(tP[nc], tE[nc], n);
	kubert_to_weierstrass(tE[nc], NULL, tP[nc], b, c, kx0, ky0, n->orig_modulus);
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

int
build_curves_with_torsion_Z9(mpz_t fac, mpmod_t n, ell_curve_t *tE, 
			     ell_point_t *tP, int umin, int umax, int nE)
{
    int u, ret = ECM_NO_FACTOR_FOUND, nc = 0;
    mpz_t A2, A1div2, x0, y0, cte, d, c, b, kx0, ky0;
    mpz_t f;
    mpres_t tmp;
    ell_curve_t E;
    ell_point_t P, Q;

    /* Eaux = [-9, 9] */
    /* Paux = [1, 1, 1] */
    mpres_init(tmp, n);
    mpz_init_set_str(f, "-9", 10);
    mpres_set_z(tmp, f, n);
    ell_curve_init_set(E, ECM_EC_TYPE_WEIERSTRASS, ECM_LAW_HOMOGENEOUS,tmp, n);
    ell_point_init(P, E, n);
    mpz_set_str(f, "1", 10); 
    mpres_set_z(P->x, f, n);
    mpz_set_str(f, "1", 10);
    mpres_set_z(P->y, f, n);
    mpres_init(P->z, n);
    mpz_set_ui(P->z, 1);
#if DEBUG_TORSION >= 2
    printf("Paux:=");
    pt_print(E, P, n);
    printf(";\n");
#endif
    mpz_init_set_ui(A2, 0);
    mpz_init_set_str(A1div2, "3", 10);
    mpz_mod(A1div2, A1div2, n->orig_modulus);
    mpz_init_set_str(x0, "2", 10);
    mpz_mod(x0, x0, n->orig_modulus);
    mpz_init_set_str(y0, "3", 10);
    mpz_mod(y0, y0, n->orig_modulus);
    mpz_init_set_ui(cte, 0);

    mpz_init(d);
    mpz_init(c);
    mpz_init(b);
    mpz_init(kx0);
    mpz_init(ky0);
    ell_point_init(Q, E, n);
    for(u = umin; u < umax; u++){
	/* update Qaux */
	mpz_set_ui(d, u);
        if(ell_point_mul(Q, d, P, E, n) == 0){
	    printf("found factor during update of Q\n");
	    mpz_set(fac, Q->x);
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
#if DEBUG_TORSION >= 2
	printf("(s, t)[%d]:=", u);
	pt_print(Q, n);
	printf(";\n");
#endif
	mpres_get_z(b, Q->x, n);
	mpres_get_z(c, Q->y, n);
#if 0
	mpres_get_z(d, Q->z, n);
	if(mpz_invert(fac, d, n->orig_modulus) == 0){
	    printf("found factor in Z9 (normalization)\n");
	    mpz_gcd(fac, d, n->orig_modulus);
	    break;
	}
	mpz_mul(b, b, fac);
	mpz_mod(b, b, n->orig_modulus);
	mpz_mul(c, c, fac);
	mpz_mod(c, c, n->orig_modulus);
#endif
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
	/* to short Weierstrass form */
	ell_curve_init(tE[nc], ECM_EC_TYPE_WEIERSTRASS, ECM_LAW_HOMOGENEOUS,n);
	ell_point_init(tP[nc], tE[nc], n);
	kubert_to_weierstrass(tE[nc], NULL, tP[nc], b, c, kx0, ky0, n->orig_modulus);
	nc++;
	if(nc >= nE)
	    break;
    }
#if DEBUG_TORSION >= 2
    printf("Curves built\n");
    pt_many_print(tE, tP, nE, n);
#endif
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
    mpz_t A2, A1div2, x0, y0, cte, d, c, b, kx0, ky0;
    mpz_t f;
    mpres_t tmp;
    ell_curve_t E;
    ell_point_t P, Q;

    mpz_init(f);
    /* Eaux = [2/3, -53/108] */
    /* Paux = [2/3, 1/2, 1] */
    mpres_init(tmp, n);
    mod_from_rat_str(f, "2/3", n->orig_modulus);
    mpres_set_z(tmp, f, n);
    ell_curve_init_set(E, ECM_EC_TYPE_WEIERSTRASS, ECM_LAW_HOMOGENEOUS, tmp,n);
    ell_point_init(P, E, n);
    mod_from_rat_str(f, "2/3", n->orig_modulus);
    mpres_set_z(P->x, f, n);
    mod_from_rat_str(f, "1/2", n->orig_modulus);
    mpres_set_z(P->y, f, n);
    mpz_set_ui(P->z, 1);
#if DEBUG_TORSION >= 2
    printf("P:=");
    pt_print(E, P, n);
    printf(";\n");
#endif

    mpz_init(A2);
    mod_from_rat_str(A2, "2/3", n->orig_modulus);
    mpz_init(A1div2);
    mod_from_rat_str(A1div2, "-1/2", n->orig_modulus);
    mpz_mod(A1div2, A1div2, n->orig_modulus);
    mpz_init_set_si(x0, 0);
    mpz_init_set_si(y0, 1);
    mpz_init_set_si(cte, -2);
    mpz_mod(cte, cte, n->orig_modulus);

    mpz_init(d);
    mpz_init(c);
    mpz_init(b);
    mpz_init(kx0);
    mpz_init(ky0);
    ell_point_init(Q, E, n);
    for(u = umin; u < umax; u++){
	/* update Qaux */
	mpz_set_ui(d, u);
	if(ell_point_mul(Q, d, P, E, n) == 0){
	    printf("found factor in Z10 (update of Q)\n");
	    mpz_set(fac, Q->x);
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
#if DEBUG_TORSION >= 2
	printf("(s, t)[%d]:=", u);
	pt_print(Q, n);
	printf(";\n");
#endif
	mpres_get_z(b, Q->x, n);
	mpres_get_z(c, Q->y, n);
#if 0
	mpres_get_z(d, Q->z, n);
	if(mpz_invert(fac, d, n->orig_modulus) == 0){
	    printf("found factor in Z10 (normalization)\n");
	    mpz_gcd(fac, d, n->orig_modulus);
	    break;
	}
	mpz_mul(b, b, fac);
	mpz_mod(b, b, n->orig_modulus);
	mpz_mul(c, c, fac);
	mpz_mod(c, c, n->orig_modulus);
#endif
	if(cubic_to_quartic(fac, n->orig_modulus, f, ky0, b, c, 
			    A2, A1div2, x0, y0, cte) == 0){
	    printf("found factor in Z10 (cubic_2_quartic)\n");
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
	/* f:=x; */
	/* d:=f^2/(f-(f-1)^2); */
	/** b <- f^2 **/
	mpz_mul(b, f, f);
	mpz_mod(b, b, n->orig_modulus);
	mpz_sub_si(c, f, 1);
	mpz_mul(c, c, c);
	mpz_sub(c, f, c);
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
	/** b <- b^2 = f^4 **/
	mpz_mul(b, b, b);
	mpz_mod(b, b, n->orig_modulus);
	mpz_mul(kx0, ky0, b);    /* num */
	mpz_sub_si(fac, f, 3);
	mpz_mul(fac, fac, f);
	mpz_add_si(fac, fac, 1);
	mpz_mul(fac, fac, fac);
	mpz_mul_si(fac, fac, 2); 
	mpz_mod(fac, fac, n->orig_modulus);    /* den */
	if(mod_from_rat2(ky0, kx0, fac, n->orig_modulus) == 0){
            printf("inverse found in Z10 (ky0)\n");
            ret = ECM_FACTOR_FOUND_STEP1;
            break;
        }
	/* kx0:=-f*d; */
	mpz_mul(kx0, f, d);
	mpz_mul_si(kx0, kx0, -1); /* humf */
	mpz_mod(kx0, kx0, n->orig_modulus);
	/* b:=c*d; */
	mpz_mul(b, c, d);
	mpz_mod(b, b, n->orig_modulus);
	/* to short Weierstrass form */
	ell_curve_init(tE[nc], ECM_EC_TYPE_WEIERSTRASS, ECM_LAW_HOMOGENEOUS,n);
	ell_point_init(tP[nc], tE[nc], n);
	kubert_to_weierstrass(tE[nc], NULL, tP[nc], b, c, kx0, ky0, n->orig_modulus);
	nc++;
	if(nc >= nE)
	    break;
    }
#if DEBUG_TORSION >= 2
    printf("Curves built\n");
    pt_many_print(tE, tP, nE, n);
#endif
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

/* Warning: b and a have the Montgomery meaning in this function. */
int
build_curves_with_torsion_Z2xZ8(mpz_t f, mpmod_t n, 
				ell_curve_t *tE, ell_point_t *tP,
				int umin, int umax, int nE)
{
    int u, nc = 0, ret = ECM_NO_FACTOR_FOUND;
    mpz_t tmp, a, b, alpha, beta, c, d, kx0, ky0, wx0;
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

    /* Eaux = [-8, -32] */
    /* Paux = [12, 40, 1] */
    mpres_init(tmp2, n);
    mpz_set_str(f, "-8", 10); 
    mpres_set_z(tmp2, f, n);
    ell_curve_init_set(E, ECM_EC_TYPE_WEIERSTRASS, ECM_LAW_AFFINE, tmp2, n);
    ell_point_init(P, E, n);
    mpz_set_str(f, "12", 10); 
    mpres_set_z(P->x, f, n);
    mpz_set_str(f, "40", 10);
    mpres_set_z(P->y, f, n);
    mpz_set_ui(P->z, 1);

    ell_point_init(Q, E, n);
    for(u = umin; u < umax; u++){
	/* update Qaux */
	mpz_set_ui(d, u);
	if(ell_point_mul(Q, d, P, E, n) == 0){
	    printf("found factor in Z2xZ8 (update of Q)\n");
	    mpz_set(f, Q->x);
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
#if DEBUG_TORSION >= 2
	printf("(s, t)[%d]:=", u);
	pt_print(E, Q, n);
	printf(";\n");
#endif
	mpres_get_z(a, Q->x, n);
	mpres_get_z(b, Q->y, n);
#if 0 /* useless in affine form? */
	mpres_get_z(d, Q->z, n);
	if(mpz_invert(f, d, n->orig_modulus) == 0){
	    printf("found factor in Z2xZ8 (normalization)\n");
	    mpz_gcd(f, d, n->orig_modulus);
	    break;
	}
	mpz_mul(a, a, f);
	mpz_mul(b, b, f);
#endif
	mpz_mod(wx0, a, n->orig_modulus);
	mpz_sub_si(a, a, 9);
	mpz_mod(a, a, n->orig_modulus);
	mpz_add_si(b, b, 25);
	mpz_mod(b, b, n->orig_modulus);
	if(mod_from_rat2(beta, b, a, n->orig_modulus) == 0){
            printf("found factor in Z2xZ8 (beta)\n");
	    mpz_set(f, beta);
	    ret = ECM_FACTOR_FOUND_STEP1;
            break;
	}
	mpz_add_si(tmp, beta, 1);
	mpz_mod(tmp, tmp, n->orig_modulus);
	if(mpz_invert(alpha, tmp, n->orig_modulus) == 0){
            printf("found factor in Z2xZ8 (alpha)\n");
	    mpz_gcd(f, tmp, n->orig_modulus);
	    ret = ECM_FACTOR_FOUND_STEP1;
            break;
	}
	/** d <- 8*alpha^2-1; **/
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
	if(mod_from_rat2(f, c, d, n->orig_modulus) == 0){
            printf("found factor in Z2xZ8 (d)\n");
	    ret = ECM_FACTOR_FOUND_STEP1;
            break;
	}
	mpz_set(d, f);
	/* c:=(2*d-1)*(d-1)/d;*/
	mpz_sub_si(f, d, 1);
	/** kx0 <- 2*d-1 **/
	mpz_mul_si(kx0, d, 2);
	mpz_sub_si(kx0, kx0, 1);
	mpz_mul(f, f, kx0);
	mpz_mod(f, f, n->orig_modulus);
        if(mod_from_rat2(c, f, d, n->orig_modulus) == 0){
            printf("found factor in Z2xZ8 (d)\n");
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
	mpz_mul(f, beta, beta);
	mpz_set(a, wx0);
	mpz_sub(f, a, f);
	mpz_add(f, f, a);
	mpz_add_si(f, f, 9);
	mpz_mul(f, f, c);
	mpz_mod(f, f, n->orig_modulus);
	mod_div_2(f, n->orig_modulus);
	mod_div_2(f, n->orig_modulus);
	mod_div_2(f, n->orig_modulus);
	/* ky0:=ky0/(beta^2+2*beta-7); */
	mpz_add_si(tmp, beta, 2);
	mpz_mul(tmp, tmp, beta);
	mpz_sub_si(tmp, tmp, 7);
	mpz_mod(tmp, tmp, n->orig_modulus);
	if(mod_from_rat2(ky0, f, tmp, n->orig_modulus) == 0){
            printf("found factor in Z2xZ8 (ky0)\n");
	    mpz_set(f, ky0);
            ret = ECM_FACTOR_FOUND_STEP1;
            break;
        }
	KW2W246(f, a, NULL, b, c, n->orig_modulus, 0);
#if DEBUG_TORSION >= 2
	gmp_printf("kwx0:=%Zd;\n", kx0);
	gmp_printf("kwy0:=%Zd;\n", ky0);
	printf("(kwy0^2-(kwx0^3+a2*kwx0^2+a4*kwx0+a6)) mod N;\n");
#endif
	/* wx0:=kx0+a2/3; */
        mpz_set_si(tmp, 3);
	mod_from_rat2(wx0, f, tmp, n->orig_modulus);
	mpz_add(wx0, wx0, kx0);
	mpz_mod(wx0, wx0, n->orig_modulus);
	/* ma:=-1/4*(8*d^4-16*d^3+16*d^2-8*d+1)/(d-1)^2/d^2; */
	mpz_sub_si(tmp, d, 1);    /* num */
	mpz_mul(tmp, tmp, d);
	mpz_mul(tmp, tmp, tmp);
	mpz_mul_si(tmp, tmp, -4);
	mpz_mod(tmp, tmp, n->orig_modulus);
	mpz_set_si(f, 8);         /* den */
	mpz_mul(f, f, d); mpz_add_si(f, f, -16);
	mpz_mul(f, f, d); mpz_add_si(f, f, 16);
	mpz_mul(f, f, d); mpz_add_si(f, f, -8);
	mpz_mul(f, f, d); mpz_add_si(f, f, 1);

	/* to Montgomery form */
	ell_curve_init(tE[nc], ECM_EC_TYPE_MONTGOMERY, ECM_LAW_HOMOGENEOUS,n);
	ell_point_init(tP[nc], tE[nc], n);
	if(mod_from_rat2(tE[nc]->a4, f, tmp, n->orig_modulus) == 0){
            printf("found factor in Z2xZ8 (ma)\n");
	    mpz_set(f, tE[nc]->a4);
            ret = ECM_FACTOR_FOUND_STEP1;
            break;
        }
	/* mb:=-1/(d-1)^2; */
	mpz_sub_si(tmp, d, 1);
	mpz_mul(tmp, tmp, tmp);
	mpz_mod(tmp, tmp, n->orig_modulus);
	if(mpz_invert(f, tmp, n->orig_modulus) == 0){
	    printf("found factor in Z2xZ8 (mb)\n");
	    mpz_gcd(f, tmp, n->orig_modulus);
            ret = ECM_FACTOR_FOUND_STEP1;
            break;
	}
	mpz_set_si(tmp, 0);
	mpz_sub(tmp, tmp, f);
	mpz_mod(tmp, tmp, n->orig_modulus);
	/* mx:=mb*wx0-ma/3; */
	mpz_mul(f, tmp, wx0);
        mpz_set_si(tmp, 3);
        mod_from_rat2(tP[nc]->x, tE[nc]->a4, tmp, n->orig_modulus);
	mpz_sub(tP[nc]->x, f, tP[nc]->x);
	mpz_mod(tP[nc]->x, tP[nc]->x, n->orig_modulus);
	/* my:=mb*ky0; */
#if DEBUG_TORSION >= 2
	gmp_printf("N:=%Zd;\n", n->orig_modulus);
	gmp_printf("ma:=%Zd;\n", tE[nc]->a4);
	gmp_printf("kx0:=%Zd;\n", kx0);
	gmp_printf("ky0:=%Zd;\n", ky0);
	gmp_printf("mx0:=%Zd;\n", tP[nc]->x);
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

/* Source: Dujella and Najman, arxiv:1201.0266v1 */
int
build_curves_with_torsion_Z3xZ3_DuNa(mpmod_t n, ell_curve_t *tE, ell_point_t *tP,
				     int smin, int smax, int nE)
{
    mpz_t a2, a4, x0, y0;
    int T, nc = 0;

    mpz_init(x0);
    mpz_init(y0);
    mpz_init(a2);
    mpz_init(a4);
    for(T = smin; T < smax; T++){
	/* x0 <- T^6 */
	mpz_ui_pow_ui(x0, T, 6);
	mpz_mod(x0, x0, n->orig_modulus);
	/* a2:=T^6+108;*/
	mpz_add_ui(a2, x0, 108);
	mpz_mod(a2, a2, n->orig_modulus);
	/* a4:=144*T^6+3888; */
	mpz_mul_ui(a4, x0, 144);
	mpz_add_ui(a4, a4, 3888);
	mpz_mod(a4, a4, n->orig_modulus);
#if 0
	{
	    mpz_t a6;
	    /* a6:=64*T^12+3456*T^6+46656; */
	    mpz_init(a6);
	    mpz_mul_ui(a6, x0, 64);
	    mpz_add_ui(a6, a6, 3456);
	    mpz_mul(a6, a6, x0);
	    mpz_add_ui(a6, a6, 46656);
	    mpz_mod(a6, a6, n);
	    mpz_clear(a6);
	}
#endif
	/* P:=E![0, 8*T^6+216, 1] has infinite order.*/
	/* convert to short Weierstrass form */
	mpz_mul_ui(y0, x0, 8);
	mpz_add_ui(y0, y0, 216);
	mpz_mod(y0, y0, n->orig_modulus);
	mpz_set_ui(x0, 0);
	ell_curve_init(tE[nc], ECM_EC_TYPE_WEIERSTRASS, ECM_LAW_HOMOGENEOUS,n);
	ell_point_init(tP[nc], tE[nc], n);
	MediumWeierstrassToShortWeierstrass(tE[nc]->a4, NULL,
					    tP[nc]->x, tP[nc]->y,
					    a2, a4, NULL, x0, y0,
					    n->orig_modulus);
	nc++;
	if(nc >= nE)
	    break;
    }
    mpz_clear(x0);
    mpz_clear(y0);
    mpz_clear(a2);
    mpz_clear(a4);
    return ECM_NO_FACTOR_FOUND;
}

/* A more simpler and more efficient stuff, using Hessian form. */
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
	mpz_set_si(u0, u);
	/* D:=RatMod((u0^3+v0^3+1)/(3*u0*v0), N); */
	mpz_mul(num, u0, u0);
	mpz_mul(num, num, u0);
	mpz_mul(den, v0, v0);
	mpz_mul(den, den, v0);
	mpz_add(num, num, den);
	mpz_add_si(num, num, 1);
	mpz_mod(num, num, n->orig_modulus);
	
	mpz_mul(den, u0, v0);
	mpz_mul_si(den, den, 3);
	mpz_mod(den, den, n->orig_modulus);

	if(mod_from_rat2(D, num, den, n->orig_modulus) == 0){
	    printf("found factor in Z3xZ3 (D)\n");
            mpz_set(f, D);
            ret = ECM_FACTOR_FOUND_STEP1;
            break;
	}
	ell_curve_init_set(tE[nc],ECM_EC_TYPE_HESSIAN,ECM_LAW_HOMOGENEOUS,D,n);
	ell_point_init(tP[nc], tE[nc], n);
	mpz_set(tE[nc]->a4, D);
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

/* For a small price, add a 2-torsion point. */
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
	mpz_set_ui(f, u);
	if(ell_point_mul(Q, f, P, E, n) == 0){
	    printf("found factor in Z3xZ6 (update of Q)\n");
	    mpz_set(f, Q->x);
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
            printf("found factor in Z3xZ6 (D)\n");
            mpz_set(f, tE[nc]->a4);
            ret = ECM_FACTOR_FOUND_STEP1;
            break;
        }
#if DEBUG_TORSION >= 1
	gmp_printf("D%d:=%Zd;\n", nc, tE[nc]->a4);
#endif
	/* u0:=RatMod(sk/tk, N); */
	if(mod_from_rat2(tP[nc]->x, sk, tk, n->orig_modulus) == 0){
            printf("found factor in Z3xZ6 (u0)\n");
            mpz_set(f, tP[nc]->x);
            ret = ECM_FACTOR_FOUND_STEP1;
            break;
        }
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
    printf("# s:");
    for(nu = smin; nu < smax; nu++){
	mpz_set_ui(nu2, nu*nu);
	/* tau:=(nu^2+3)/2/nu; */
	mpz_add_si(lambda, nu2, 3);
	mpz_set_si(tmp, 2*nu);
	if(mod_from_rat2(tau, lambda, tmp, n->orig_modulus) == 0){
            printf("\nFactor found durint init of Z4xZ4 (tau)\n");
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

	/* Montgomery form: there are several b possible */
	/* b:=1/9/lambda^2/(tau^4-1); */
	mpz_powm_ui(x0, tau, 4, n->orig_modulus);
	mpz_sub_si(x0, x0, 1);
	mpz_mod(x0, x0, n->orig_modulus);
	mpz_mul(tmp, x0, lambda);
	mpz_mul(tmp, tmp, lambda);
	mpz_mul_si(tmp, tmp, 9);
	if(mpz_invert(b, tmp, n->orig_modulus) == 0){
	    printf("\nFactor found durint init of Z4xZ4 (b)\n");
	    mpz_gcd(f, tmp, n->orig_modulus);
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
	/* a:=-2*(tau^4+1)/(tau^4-1); */
	mpz_add_si(tmp, x0, 2);
	mpz_mul_si(tmp, tmp, -2);
	mpz_mod(tmp, tmp, n->orig_modulus);
        /* to Montgomery form */
        ell_curve_init(tE[nc], ECM_EC_TYPE_MONTGOMERY, ECM_LAW_HOMOGENEOUS, n);
        ell_point_init(tP[nc], tE[nc], n);
	if(mod_from_rat2(tE[nc]->a4, tmp, x0, n->orig_modulus) == 0){
	    printf("\nFactor found durint init of Z4xZ4 (a)\n");
	    mpz_set(f, tE[nc]->a4);
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
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
	printf(" %d", nu);
	nc++;
	if(nc >= nE)
	    break;
    }
    printf("\n");
    mpz_clear(tau);
    mpz_clear(lambda);
    mpz_clear(nu2);
    mpz_clear(tmp);
    mpz_clear(b);
    mpz_clear(x0);
    if(nc < nE){
	printf("Not enough curves generated\n");
	return ECM_ERROR;
    }
    return ret;
}

/* coeffs = {deg, c_deg, ..., c_0} */
void
mpz_eval_poly(mpz_t y, long* coeffs, mpz_t x, mpz_t n)
{
    int deg, i;

    deg = coeffs[0];
    /* Horner's rule */
    mpz_set_si(y, coeffs[1]);
    for(i = 2; i <= deg+1; i++){
	mpz_mul(y, y, x);
	mpz_add_si(y, y, coeffs[i]);
	mpz_mod(y, y, n);
    }
}

void
ec_force_point(ell_curve_t E, ell_point_t P, mpz_t B, long *x0, mpz_t n)
{
    mpz_t lambda;

    mpz_init(lambda);
    while(1){
	*x0 += 1;
	/* lambda <- x0^3+A*x0+B = (x0^2+A)*x0+B */
	mpz_set_si(lambda, *x0);
	mpz_mul_si(lambda, lambda, *x0);
	mpz_add(lambda, lambda, E->a4);
	mpz_mul_si(lambda, lambda, *x0);
	mpz_add(lambda, lambda, B);
	mpz_mod(lambda, lambda, n);
#if 0
	gmp_printf("lambda:=%Zd;\n", lambda);
	printf("# jac=%d\n", mpz_jacobi(lambda, n));
#endif
	/* P */
	mpz_mul_si(P->x, lambda, *x0);
	mpz_mod(P->x, P->x, n);
	mpz_mul(P->y, lambda, lambda);
	mpz_mod(P->y, P->y, n);
        mpz_set_si(P->z, 1);
	/* modify E */
	mpz_mul(E->a4, E->a4, P->y);
	mpz_mod(E->a4, E->a4, n);
	/* not really needed */
	mpz_mul(B, B, P->y);
	mpz_mul(B, B, lambda);
	mpz_mod(B, B, n);
	break;
    }
#if DEBUG_TORSION >= 2
    gmp_printf("newA:=%Zd;\n", E->a4);
    gmp_printf("newB:=%Zd;\n", B);
#endif
    mpz_clear(lambda);
}

/* Source: Brier+Clavier or Kohel or Silverberg or Klein. */
int
build_curves_with_torsion_Z5xZ5(mpmod_t n, ell_curve_t *tE,
				ell_point_t *tP,
				int smin, int smax, int nE)
{
    int u, nc = 0, ret = ECM_NO_FACTOR_FOUND;
    long polyA[] = {4, 1, -228, 494, 228, 1}; /* deg, c_deg, ..., c_0 */
    long polyB[] = {6, 1, 522, -10005, 0, -10005, -522, 1};
    long x0;
    mpz_t t, num, den, B;

    mpz_init(t);
    mpz_init(num);
    mpz_init(den);
    mpz_init(B);
    printf("# s:");
    for(u = smin; u < smax; u++){
	printf(" %d", u);
	mpz_set_ui(t, u);
	mpz_powm_ui(t, t, 5, n->orig_modulus);
	ell_curve_init(tE[nc], ECM_EC_TYPE_WEIERSTRASS, ECM_LAW_HOMOGENEOUS,n);
	/* A = -(t^4-228*t^3+494*t^2+228*t+1)/48; */
	mpz_eval_poly(num, polyA, t, n->orig_modulus);
	mpz_sub_si(den, n->orig_modulus, 48);
	mod_from_rat2(tE[nc]->a4, num, den, n->orig_modulus);
	/* B = (t^6+522*t^5-10005*t^4-10005*t^2-522*t+1)/864 */
	mpz_eval_poly(num, polyB, t, n->orig_modulus);
	mpz_set_si(den, 864);
	mod_from_rat2(B, num, den, n->orig_modulus);
	ell_point_init(tP[nc], tE[nc], n);
	x0 = 0;
	ec_force_point(tE[nc], tP[nc], B, &x0, n->orig_modulus);
	nc++;
	if(nc >= nE)
	    break;
    }
    printf("\n");
    mpz_clear(t);
    mpz_clear(num);
    mpz_clear(den);
    mpz_clear(B);
    return ret;
}

int
build_curves_with_torsion_Z2xZ10(mpz_t f, mpmod_t n, ell_curve_t *tE,
				 ell_point_t *tP,
				 int smin, int smax, int nE)
{
    int u, nc = 0, ret = ECM_NO_FACTOR_FOUND;
    long poly1[] = {3, 2, -3, 1, 0}; /* deg, c_deg, ..., c_0 */
    long poly2[] = {2, 1, -3, 1}; /* deg, c_deg, ..., c_0 */
    ATTRIBUTE_UNUSED long polydisc[] = {3, 8, -8, 0, 1};
    long x0;
    mpz_t t, num, den, b, c, B;

    mpz_init(t);
    mpz_init(num);
    mpz_init(den);
    mpz_init(b);
    mpz_init(c);
    mpz_init(B);
#if DEBUG_TORSION >= 2
    printf("info:=[];\n");
#endif
    for(u = smin; u < smax; u++){
#if DEBUG_TORSION >= 2
	printf("info[%d]:=\"Z2xZ10:%d\";\n", nc+1, u);
#endif
	mpz_set_ui(t, u);
	/* num = t*(2*t^2-3*t+1) */
	mpz_eval_poly(num, poly1, t, n->orig_modulus);
	/* den = t^2-3*t+1 */
	mpz_eval_poly(den, poly2, t, n->orig_modulus);
	/* b:=RatMod(t^3*(2*t^2-3*t+1)/(t^2-3*t+1)^2, N); */
	/* c:=RatMod(-t*(2*t^2-3*t+1)/(t^2-3*t+1), N); */
	if(mpz_invert(f, den, n->orig_modulus) == 0){
	    printf("# factor found in Z2xZ10 (den)\n");
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
	mpz_mul(c, num, f);
	mpz_mod(c, c, n->orig_modulus);
	mpz_mul(b, c, f);
	mpz_mul(b, b, t);
	mpz_mul(b, b, t);
	mpz_mod(b, b, n->orig_modulus);
	mpz_mul_si(c, c, -1);
	mpz_mod(c, c, n->orig_modulus);
#if DEBUG_TORSION >= 2
	gmp_printf("t:=%Zd;\n", t);
	gmp_printf("b:=%Zd;\n", b);
	gmp_printf("c:=%Zd;\n", c);
#endif	
	ell_curve_init(tE[nc], ECM_EC_TYPE_WEIERSTRASS, ECM_LAW_HOMOGENEOUS,n);
	ell_point_init(tP[nc], tE[nc], n);
	/* transform with no input point */
	kubert_to_weierstrass(tE[nc],B,NULL,b,c,NULL,NULL,n->orig_modulus);
#if DEBUG_TORSION >= 2
        gmp_printf("A:=%Zd;\n", tE[nc]->a4);
        gmp_printf("B:=%Zd;\n", B);
	mpz_eval_poly(num, polydisc, t, n->orig_modulus);
	gmp_printf("# disc:=%Zd;\n", num);
#endif
	x0 = 0;
	ec_force_point(tE[nc], tP[nc], B, &x0, n->orig_modulus);
	nc++;
	if(nc >= nE)
	    break;
    }
    mpz_clear(t);
    mpz_clear(num);
    mpz_clear(den);
    mpz_clear(b);
    mpz_clear(c);
    mpz_clear(B);
    return ret;
}

int
build_curves_with_torsion_Z2xZ12(mpz_t f, mpmod_t n, ell_curve_t *tE,
				 ell_point_t *tP,
				 int smin, int smax, int nE)
{
    int u, nc = 0, ret = ECM_NO_FACTOR_FOUND;
    long x0;
    mpz_t t, num, den, b, c, B;

    mpz_init(t);
    mpz_init(num);
    mpz_init(den);
    mpz_init(b);
    mpz_init(c);
    mpz_init(B);
#if DEBUG_TORSION >= 2
    printf("info:=[];\n");
#endif
    for(u = smin; u < smax; u++){
#if DEBUG_TORSION >= 2
	printf("info[%d]:=\"Z2xZ12:%d\";\n", nc+1, u);
#endif
	mpz_set_ui(t, u);
	/* num = 1-t^2 */
	/* den = t^2*(t^2+3) */
	/* c:=RatMod(num/den, N); */
	mpz_mul(num, t, t);
	mpz_mod(num, num, n->orig_modulus);
	mpz_add_si(den, num, 3);
	mpz_mul(den, den, num);
	mpz_mod(den, den, n->orig_modulus);
	mpz_neg(num, num);
	mpz_add_si(num, num, 1);
	mpz_mod(num, num, n->orig_modulus);
	if(mod_from_rat2(c, num, den, n->orig_modulus) == 0){
	    printf("# factor found in Z2xZ12 (den)\n");
	    mpz_set(f, c);
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
	mpz_mul(b, c, c);
	mpz_add(b, b, c);
	mpz_mod(b, b, n->orig_modulus);
#if DEBUG_TORSION >= 2
	gmp_printf("t:=%Zd;\n", t);
	gmp_printf("b:=%Zd;\n", b);
	gmp_printf("c:=%Zd;\n", c);
#endif	
	ell_curve_init(tE[nc], ECM_EC_TYPE_WEIERSTRASS, ECM_LAW_HOMOGENEOUS,n);
	ell_point_init(tP[nc], tE[nc], n);
	/* transform with no input point */
	kubert_to_weierstrass(tE[nc],B,NULL,b,c,NULL,NULL,n->orig_modulus);
#if DEBUG_TORSION >= 2
        gmp_printf("A:=%Zd;\n", tE[nc]->a4);
        gmp_printf("B:=%Zd;\n", B);
#endif
	x0 = 0;
	ec_force_point(tE[nc], tP[nc], B, &x0, n->orig_modulus);
	nc++;
	if(nc >= nE)
	    break;
    }
    mpz_clear(t);
    mpz_clear(num);
    mpz_clear(den);
    mpz_clear(b);
    mpz_clear(c);
    mpz_clear(B);
    return ret;
}

/***** From Rabarison10 *****/

/* Operate over Q(sqrt(d)); one must have sqroots[0] = sqrt(d) mod n. */
int
build_curves_with_torsion_Z11(mpz_t f, mpmod_t n, 
			      ell_curve_t *tE, ell_point_t *tP,
			      int smin, int smax, int nE,
			      int d, mpz_t *sqroots)
{
    int ret = ECM_NO_FACTOR_FOUND, u, nc = 0, i;
    ell_curve_t Eaux;
    ell_point_t Paux, kPaux;
    long x0;
    mpz_t aux1, aux2, aux3, aux4, aux6, xaux, yaux, s, t, b, c, tmp, B;

    /* Eaux = [a1, a3, a2, a4, a6] */
    mpz_init_set_si(aux1, 0);
    mpz_init_set_si(aux3, -1);
    mpz_init_set_si(aux2, -1);
    mpz_init_set_si(aux4, 0);
    mpz_init_set_si(aux6, 0);
    /* use points on Eaux */
    if(d == 2){
	/* Paux = [1/2, (2-Z)/4, 1]; */
	mpz_init(xaux);
	mod_from_rat_str(xaux, "1/2", n->orig_modulus);
	mpz_init_set_si(yaux, 2);
	mpz_sub(yaux, yaux, sqroots[0]);
	mpz_mod(yaux, yaux, n->orig_modulus);
	mod_div_2(yaux, n->orig_modulus);
	mod_div_2(yaux, n->orig_modulus);
    }
    mpz_init(t);
    mpz_init(s);
    mpz_init(b);
    mpz_init(c);
    mpz_init(tmp);
    mpz_init(B);
    /* make Eaux in short Weierstrass form */
    ell_curve_init(Eaux, ECM_EC_TYPE_WEIERSTRASS, ECM_LAW_AFFINE, n);
    ell_point_init(Paux, Eaux, n);
    /* make it mpres */
    mpres_set_z(Eaux->a1, aux1, n);
    mpres_set_z(Eaux->a3, aux3, n);
    mpres_set_z(Eaux->a2, aux2, n);
    mpres_set_z(Eaux->a4, aux4, n);
    mpres_set_z(Eaux->a6, aux6, n);
    mpres_set_z(Paux->x, xaux, n);
    mpres_set_z(Paux->y, yaux, n);
    mpz_set_ui(Paux->z, 1);
    mpz_clear(aux1);
    mpz_clear(aux3);
    mpz_clear(aux2);
    mpz_clear(aux4);
    mpz_clear(aux6);
    mpz_clear(xaux);
    mpz_clear(yaux);
    ell_point_init(kPaux, Eaux, n);
    for(u = smin; u < smax; u++){
#if DEBUG_TORSION >= 2
	printf("info[%d]:=\"Z11:%d\";\n", nc+1, u);
#endif
	mpz_set_si(t, u);
	if(ell_point_mul(kPaux, t, Paux, Eaux, n) == 0){
	    printf("# found factor during update(Z11)\n");
	    mpz_set(f, kPaux->x);
            ret = ECM_FACTOR_FOUND_STEP1;
            break;
	}
	mpres_get_z(t, kPaux->x, n);
	mpres_get_z(s, kPaux->y, n);
	/* b:=s*(s-1)*(s-t)/t mod N;
	   a:=(s*t+t-s^2)/t mod N = (t-s*(s-t))/t; */
	if(mpz_invert(tmp, t, n->orig_modulus) == 0){
	    printf("# found factor during inversion(Z11)\n");
	    mpz_gcd(f, t, n->orig_modulus);
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
	/** b <- s*(s-t) **/
	mpz_sub(b, s, t);
	mpz_mul(b, b, s);
	/** c <- t-s*(s-t) **/
	mpz_sub(c, t, b);
	mpz_mul(c, c, tmp);
	mpz_mod(c, c, n->orig_modulus);
	/** b <- b/t*(s-1) **/
	mpz_mul(b, b, tmp);
	mpz_sub_si(tmp, s, 1);
	mpz_mul(b, b, tmp);
	mpz_mod(b, b, n->orig_modulus);
	/* c:=1-a mod N; */
	mpz_sub_si(c, c, 1);
	mpz_neg(c, c);
	mpz_mod(c, c, n->orig_modulus);
	/* E:=QEcFromKubertForm(-b, c) mod N; */
	mpz_neg(b, b);
	mpz_mod(b, b, n->orig_modulus);
	/* forcing this number of twist candidates */
	x0 = 0;
	for(i = 0; i < 3; i++){
	    /* lazy, isn't it? */
	    ell_curve_init(tE[nc], ECM_EC_TYPE_WEIERSTRASS, 
			   ECM_LAW_HOMOGENEOUS, n);
	    ell_point_init(tP[nc], tE[nc], n);
	    kubert_to_weierstrass(tE[nc], B, NULL, b, c, NULL, NULL,
				  n->orig_modulus);
	    ec_force_point(tE[nc], tP[nc], B, &x0, n->orig_modulus);
	    nc++;
	    if(nc >= nE)
		break;
	}
	if(nc >= nE)
	    break;
    }
    mpz_clear(t);
    mpz_clear(s);
    mpz_clear(b);
    mpz_clear(c);
    mpz_clear(tmp);
    mpz_clear(B);
    ell_point_clear(Paux, Eaux, n);
    ell_point_clear(kPaux, Eaux, n);
    ell_curve_clear(Eaux, n);
    return ret;
}
