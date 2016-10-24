/* addlaws.c - various addition laws for ECM
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

#if DEBUG_ADD_LAWS >= 1
void
print_mpz_from_mpres(mpres_t x, mpmod_t n)
{
    mpz_t tmp;

    mpz_init(tmp);
    mpres_get_z(tmp, x, n);
    gmp_printf("%Zd", tmp);
    mpz_clear(tmp);
}
#endif

/******************** Weierstrass section ********************/

void
pt_w_set_to_zero(ell_point_t P, mpmod_t n)
{
    mpres_set_ui(P->x, 0, n);
    mpres_set_ui(P->y, 1, n);
    mpres_set_ui(P->z, 0, n);
}

int
pt_w_is_zero(mpres_t z, mpmod_t n)
{
    return mpres_is_zero(z, n);
}

void
pt_w_set(mpres_t x0, mpres_t y0, mpres_t z0,
	    mpres_t x, mpres_t y, mpres_t z,
	    ATTRIBUTE_UNUSED mpmod_t n)
{
  mpres_set(x0, x, n);
  mpres_set(y0, y, n);
  mpres_set(z0, z, n);
}

#if DEBUG_ADD_LAWS >= 1
void
pt_w_print(mpres_t x, mpres_t y, mpres_t z, ell_curve_t E, mpmod_t n)
{
    printf("[");
    print_mpz_from_mpres(x, n);
    printf(", ");
    print_mpz_from_mpres(y, n);
    printf(", ");
    if(E->type == ECM_EC_TYPE_WEIERSTRASS && E->law == ECM_LAW_AFFINE)
	gmp_printf("%Zd", z);
    else
	print_mpz_from_mpres(z, n);
    printf("]");
}
#endif

/* [x0, y0, z0] <- [x1, y1, z1] + [x2, y2, z2] using lambda=num/den
   with buffer inv.

   (lambda*x+mu)^2+a1*x*(lambda*x+mu)+a3*(lambda*x+mu)=x^3+a2*x^2+...
   x^3+(a2-lambda^2-a1*lambda)*x^2+... = 0
   x1+x2+x3 = lambda^2+a1*lambda-a2.
   y3 = lambda*(x1-x3)-y1-a1*x3-a3
 */
static int
pt_w_common_aff(mpz_t f, mpres_t x0, mpres_t y0, mpres_t z0,
		mpres_t x1, mpres_t y1,
		mpres_t x2, mpres_t a1, mpres_t a3, mpres_t a2,
		mpmod_t n, mpres_t num, mpres_t den, mpres_t lambda)
{
    if(mpres_invert(lambda, den, n) == 0){
	mpres_gcd(f, den, n);
	return 0;
    }
    /** lambda = num/den **/
    mpres_mul(lambda, lambda, num, n);
    /** num <- (lambda+a1)*lambda **/
    mpres_add(num, lambda, a1, n);
    mpres_mul(num, num, lambda, n);
    mpres_sub(num, num, a2, n);
    /** x0 = den <- num-x1-x2 **/
    mpres_sub(den, num, x1, n);
    mpres_sub(den, den, x2, n);
    /** y0 = num <- lambda*(x1-x0)-(y1+a1*x0+a3) **/
    mpres_sub(num, x1, den, n);
    mpres_mul(num, num, lambda, n);
    mpres_sub(y0, num, y1, n);
    mpres_sub(y0, y0, a3, n);
    mpres_mul(x0, a1, den, n);
    mpres_sub(y0, y0, x0, n);
    /** finish **/
    mpres_set(x0, den, n);
    mpz_set_ui(z0, 1); /* just in case */
    return 1;
}

/* [x3, y3, z3] <- [2] * [x1, y1, z1] */
int
pt_w_duplicate(mpz_t f, mpres_t x3, mpres_t y3, mpres_t z3,
	       mpres_t x1, mpres_t y1, mpres_t z1,
	       mpmod_t n, ell_curve_t E)
{
    if(pt_w_is_zero(z1, n) == 1){
      pt_w_set(x3, y3, z3, x1, y1, z1, n);
      return 1;
    }
    if(E->type == ECM_EC_TYPE_WEIERSTRASS && E->law == ECM_LAW_AFFINE){
	/* buf[1] <- 2*y1+a1*x1+a3 */
	mpres_mul(E->buf[1], E->a1, x1, n);
	mpres_add(E->buf[1], E->buf[1], E->a3, n);
	mpres_add(E->buf[1], E->buf[1], y1, n);
	mpres_add(E->buf[1], E->buf[1], y1, n);
	if(mpres_is_zero(E->buf[1], n)){
	    /* buf1 = 0 <=> P is a [2]-torsion point */
	    mpres_set_ui(x3, 0, n);
	    mpres_set_ui(y3, 1, n);
	    mpres_set_ui(z3, 0, n);
	    return 1;
	}
	/* buf[0] <- 3*x^2+2*a2*x+a4-a1*y = (3*x+2*a2)*x+a4-a1*y */
	mpres_mul_ui(E->buf[0], x1, 3, n);
	mpres_add(E->buf[0], E->buf[0], E->a2, n);
	mpres_add(E->buf[0], E->buf[0], E->a2, n);
	mpres_mul(E->buf[0], E->buf[0], x1, n);
	mpres_add(E->buf[0], E->buf[0], E->a4, n);
	mpres_mul(E->buf[2], E->a1, y1, n);
	mpres_sub(E->buf[0], E->buf[0], E->buf[2], n);
	return pt_w_common_aff(f, x3, y3, z3, x1, y1, x1, 
			       E->a1, E->a3, E->a2, n, 
			       E->buf[0], E->buf[1], E->buf[2]);
    }
    else if(E->type == ECM_EC_TYPE_WEIERSTRASS 
	    && E->law == ECM_LAW_HOMOGENEOUS){
	/* source is dbl-2007-bl: 5M + 6S + 1*a + 7add + 3*2 + 1*3 */
	/* mapping: h = buf[0], w = buf[1], s = buf[2], RR = buf[3], B = buf[4];*/
	/*  h:=X1^2 mod p;	    # S*/
	mpres_sqr(E->buf[0], x1, n);
	/*	w:=Z1^2 mod p;*/
	mpres_sqr(E->buf[1], z1, n);
	/*	w:=a*w mod p;*/
	mpres_mul(E->buf[1], E->buf[1], E->a4, n);
	/*	s:=3*h mod p;	    # *3*/
	mpres_mul_ui(E->buf[2], E->buf[0], 3, n);
	/*	w:=w+s mod p;*/
	mpres_add(E->buf[1], E->buf[1], E->buf[2], n);
	/*	s:=Y1*Z1 mod p;*/
	mpres_mul(E->buf[2], y1, z1, n);
	/*	s:=2*s mod p;*/
	mpres_mul_ui(E->buf[2], E->buf[2], 2, n);
	/*	Z3:=s^2 mod p;*/
	mpres_sqr(z3, E->buf[2], n);
	/*	Z3:=s*Z3 mod p;*/
	mpres_mul(z3, z3, E->buf[2], n);
	/*	RR:=Y1*s mod p;	    # M*/
	mpres_mul(E->buf[3], y1, E->buf[2], n);
	/*	B:=X1+RR mod p;    # add*/
	mpres_add(E->buf[4], x1, E->buf[3], n);
	/*	B:=B^2 mod p;*/
	mpres_sqr(E->buf[4], E->buf[4], n);
	/*	RR:=RR^2 mod p;	    # S*/
	mpres_sqr(E->buf[3], E->buf[3], n);
	/*	B:=B-h mod p;*/
	mpres_sub(E->buf[4], E->buf[4], E->buf[0], n);
	/*	B:=B-RR mod p;*/
	mpres_sub(E->buf[4], E->buf[4], E->buf[3], n);
	/*	h:=w^2 mod p;*/
	mpres_sqr(E->buf[0], E->buf[1], n);
	/*	X3:=2*B mod p;*/
	mpres_mul_ui(x3, E->buf[4], 2, n);
	/*	h:=h-X3 mod p;*/
	mpres_sub(E->buf[0], E->buf[0], x3, n);
	/*	X3:=h*s mod p;	    # M*/
	mpres_mul(x3, E->buf[0], E->buf[2], n);
	/*	s:=B-h mod p;*/
	mpres_sub(E->buf[2], E->buf[4], E->buf[0], n);
	/*	s:=w*s mod p;*/
	mpres_mul(E->buf[2], E->buf[2], E->buf[1], n);
	/*	Y3:=2*RR mod p;*/
	mpres_mul_ui(y3, E->buf[3], 2, n);
	/*	Y3:=s-Y3 mod p;*/
	mpres_sub(y3, E->buf[2], y3, n);
	return 1;
    }
    return 0;
}

/* [x3, y3, z3] <- [x1, y1, z1] + [x2, y2, z2]; P3 can be either P1 or P2. */
int
pt_w_add(mpz_t f, mpres_t x3, mpres_t y3, mpres_t z3,
	 mpres_t x1, mpres_t y1, mpres_t z1,
	 mpres_t x2, mpres_t y2, mpres_t z2,
	 mpmod_t n, ell_curve_t E)
{
    if(pt_w_is_zero(z1, n)){
	pt_w_set(x3, y3, z3, x2, y2, z2, n);
	return 1;
    }
    else if(pt_w_is_zero(z2, n)){
	pt_w_set(x3, y3, z3, x1, y1, z1, n);
	return 1;
    }
    if(E->type == ECM_EC_TYPE_WEIERSTRASS && E->law == ECM_LAW_AFFINE)
	if(mpres_equal(x1, x2, n) && mpres_equal(y1, y2, n))
	    return pt_w_duplicate(f, x3, y3, z3, x1, y1, z1, n, E);
	else{
	    mpres_sub(E->buf[0], y1, y2, n);
	    mpres_sub(E->buf[1], x1, x2, n);
	    return pt_w_common_aff(f, x3, y3, z3, x1, y1, x2, 
				   E->a1, E->a3, E->a2,
				   n, E->buf[0], E->buf[1], E->buf[2]);
	}
    else if(E->type == ECM_EC_TYPE_WEIERSTRASS 
	    && E->law == ECM_LAW_HOMOGENEOUS){
	/* Cohen-Miyaji-Ono: 12M+2S+6add+1*2 */
	/* mapping: y1z2 = buf, AA = buf+1, u = buf+2, v = buf+3, R = buf+4, */
	/* vvv = buf+5; */
#if DEBUG_ADD_LAWS >= 2
	printf("y1="); print_mpz_from_mpres(y1, n); printf("\n");
	printf("y2="); print_mpz_from_mpres(y2, n); printf("\n");
	printf("z1="); print_mpz_from_mpres(z1, n); printf("\n");
	printf("z2="); print_mpz_from_mpres(z2, n); printf("\n");
#endif
	/*  Y1Z2:=Y1*Z2 mod p;	# M*/
	mpres_mul(E->buf[0], y1, z2, n);
	/*	A:=X1*Z2 mod p;	# M*/
	mpres_mul(E->buf[1], x1, z2, n);
	/*	u:=Y2*Z1 mod p;*/
	mpres_mul(E->buf[2], y2, z1, n);
	/*	u:=u-Y1Z2 mod p;*/
	mpres_sub(E->buf[2], E->buf[2], E->buf[0], n);
	/*	v:=X2*Z1 mod p;*/
	mpres_mul(E->buf[3], x2, z1, n);
	/*	v:=v-A mod p;*/
	mpres_sub(E->buf[3], E->buf[3], E->buf[1], n);
	if(mpz_sgn(E->buf[2]) == 0 && mpz_sgn(E->buf[3]) == 0){
	    /* u = 0 <=> Y2*Z1 = Y1*Z2 <=> Y2/Z2 = Y1/Z1*/
	    /* v = 0 <=> X2*Z1 = X1*Z2 <=> X2/Z2 = X1/Z1*/
	    return pt_w_duplicate(f, x3, y3, z3, x1, y1, z1, n, E);
	}
	/*	Z3:=Z1*Z2 mod p;	# M*/
	mpres_mul(z3, z1, z2, n);
	/*	X3:=u^2 mod p;*/
	mpres_sqr(x3, E->buf[2], n);
	/*	X3:=X3*Z3 mod p;*/
	mpres_mul(x3, x3, z3, n);
	/*	R:=v^2 mod p;*/
	mpres_sqr(E->buf[4], E->buf[3], n);
	/*	vvv:=v*R mod p;*/
	mpres_mul(E->buf[5], E->buf[3], E->buf[4], n);
	/*	R:=R*A mod p;*/
	mpres_mul(E->buf[4], E->buf[4], E->buf[1], n);
	/*	Y3:=2*R mod p; 		# *2*/
	mpres_mul_ui(y3, E->buf[4], 2, n);
	/*	A:=X3-vvv mod p;*/
	mpres_sub(E->buf[1], x3, E->buf[5], n);
	/*	A:=A-Y3 mod p;*/
	mpres_sub(E->buf[1], E->buf[1], y3, n);
	/*	X3:=v*A mod p;		# M*/
	mpres_mul(x3, E->buf[3], E->buf[1], n);
	/*	Y3:=R-A mod p;*/
	mpres_sub(y3, E->buf[4], E->buf[1], n);
	/*	Y3:=u*Y3 mod p;*/
	mpres_mul(y3, y3, E->buf[2], n);
	/*	A:=vvv*Y1Z2 mod p;*/
	mpres_mul(E->buf[1], E->buf[5], E->buf[0], n);
	/*	Y3:=Y3-A mod p;*/
	mpres_sub(y3, y3, E->buf[1], n);
	/*	Z3:=vvv*Z3 mod p;	# M*/
	mpres_mul(z3, z3, E->buf[5], n);
	return 1;
    }
    return 0;
}

#if USE_ADD_SUB_CHAINS > 0
/* [x3, y3, z3] <- [x1, y1, z1] - [x2, y2, z2]; P3 != P1, P3 != P2. 
   -P2 ~ -(x2/z2, y2/z2, 1) = (x2/z2, -y2/z2-a1*x/z2-a3, 1) 
                            ~ (x2, -y2-a1*x2-a3*z2, z2).
*/
int
pt_w_sub(mpz_t f, mpres_t x3, mpres_t y3, mpres_t z3,
	 mpres_t x1, mpres_t y1, mpres_t z1,
	 mpres_t x2, mpres_t y2, mpres_t z2,
	 mpmod_t n, ell_curve_t E)
{
    int res = 1;

    if(E->law == ECM_LAW_HOMOGENEOUS){
	/* FIXME: does not work for complete equation! */
	mpres_neg(y2, y2, n);
	res = pt_w_add(f, x3, y3, z3, x1, y1, z1, x2, y2, z2, n, E);
	mpres_neg(y2, y2, n);
    }
    else if(E->law == ECM_LAW_AFFINE){
	/* buf[3] not used in law, so use it */
	mpres_mul(E->buf[3], E->a1, x2, n);
	mpres_add(E->buf[3], E->buf[3], E->a3, n);
	mpres_add(E->buf[3], E->buf[3], y2, n);
	mpres_neg(E->buf[3], E->buf[3], n);
	res = pt_w_add(f, x3, y3, z3, x1, y1, z1, x2, E->buf[3], z2, n, E);
    }
    return res;
}
#endif

/******************** projective Hessian form ********************/

/* U^3+V^3+W^3 = 3*D*U*V*W, D^3 <> 1.
   O_H = [1:-1:0]
   -[u:v:w] = [v:u:w]
   Warning: there can exist two other points at infinity, namely
   [1:-omega:0] and [1:-omega^2:0] where omega^3 = 1.
*/
int
hessian_is_zero(ell_point_t P, ATTRIBUTE_UNUSED ell_curve_t E, mpmod_t n)
{
    mpres_t tmp;
    int ret;

    if(mpz_sgn(P->z) != 0)
	return 0;
    mpres_init(tmp, n);
    mpres_add(tmp, P->x, P->y, n);
    ret = mpz_sgn(tmp) == 0;
#if 0
    if(ret)
	gmp_printf("found a third root of unity? %Zd/%Zd\n", P->x, P->y);
#endif
    mpres_clear(tmp, n);
    return ret;
}

void
hessian_set_to_zero(ell_point_t P, ATTRIBUTE_UNUSED ell_curve_t E, mpmod_t n)
{
    mpres_set_si(P->x,  1, n);
    mpres_set_si(P->y, -1, n);
    mpres_set_si(P->z,  0, n);
}

#if DEBUG_ADD_LAWS >= 1
void
hessian_print(ell_point_t P, ell_curve_t E, mpmod_t n)
{
    pt_w_print(P->x, P->y, P->z, E, n);
}
#endif

#if USE_ADD_SUB_CHAINS > 0
/* -[u:v:w] = [v:u:w] */
void
hessian_negate(ell_point_t P, ATTRIBUTE_UNUSED ell_curve_t E, ATTRIBUTE_UNUSED mpmod_t n)
{
    mpz_swap(P->x, P->y); /* humf */
}
#endif

/* TODO: decrease the number of buffers? */
int
hessian_duplicate(ell_point_t R, ell_point_t P, 
		  ATTRIBUTE_UNUSED ell_curve_t E, mpmod_t n)
{
    /* A = buf[0], ..., G = buf[6], H = buf[7], J = buf[8] */
    /* A:=P[1]^2 mod N; */
    mpres_mul(E->buf[0], P->x, P->x, n);
    /* B:=P[2]^2 mod N; */
    mpres_mul(E->buf[1], P->y, P->y, n);
    /* C:=P[3]^2 mod N; */
    mpres_mul(E->buf[2], P->z, P->z, n);
    /* D:=(A+B) mod N; */
    mpres_add(E->buf[3], E->buf[0], E->buf[1], n);
    /* E:=(A+C) mod N; */
    mpres_add(E->buf[4], E->buf[0], E->buf[2], n);
    /* F:=(B+C) mod N; */
    mpres_add(E->buf[5], E->buf[1], E->buf[2], n);
    /* G:=((P[1]+P[2])^2-D) mod N; */
    mpres_add(E->buf[6], P->x, P->y, n);
    mpres_mul(E->buf[6], E->buf[6], E->buf[6], n);
    mpres_sub(E->buf[6], E->buf[6], E->buf[3], n);
    /* H:=((P[1]+P[3])^2-E) mod N; */
    mpres_add(E->buf[7], P->x, P->z, n);
    mpres_mul(E->buf[7], E->buf[7], E->buf[7], n);
    mpres_sub(E->buf[7], E->buf[7], E->buf[4], n);
    /* J:=((P[2]+P[3])^2-F) mod N; */
    mpres_add(E->buf[8], P->y, P->z, n);
    mpres_mul(E->buf[8], E->buf[8], E->buf[8], n);
    mpres_sub(E->buf[8], E->buf[8], E->buf[5], n);
    /* R->x = ((J-G)*(H+2*E)) mod N */
    mpres_sub(E->buf[0], E->buf[8], E->buf[6], n);
    mpres_add(E->buf[1], E->buf[7], E->buf[4], n);
    mpres_add(E->buf[1], E->buf[1], E->buf[4], n);
    mpres_mul(R->x, E->buf[0], E->buf[1], n);
    /* R->y = ((G-H)*(J+2*F)) mod N */
    mpres_sub(E->buf[0], E->buf[6], E->buf[7], n);
    mpres_add(E->buf[1], E->buf[8], E->buf[5], n);
    mpres_add(E->buf[1], E->buf[1], E->buf[5], n);
    mpres_mul(R->y, E->buf[0], E->buf[1], n);
    /* R->z = ((H-J)*(G+2*D)) mod N */
    mpres_sub(E->buf[0], E->buf[7], E->buf[8], n);
    mpres_add(E->buf[1], E->buf[6], E->buf[3], n);
    mpres_add(E->buf[1], E->buf[1], E->buf[3], n);
    mpres_mul(R->z, E->buf[0], E->buf[1], n);
    return 1;
}

/* TODO: reduce the number of buffers? */
int
hessian_plus(ell_point_t R, ell_point_t P, ell_point_t Q, 
	     ATTRIBUTE_UNUSED ell_curve_t E, mpmod_t n)
{
    /* P = [T1,T2,T3], Q = [T4,T5,T6] */
    /* P = Q <=> T1/T3=T4/T6 and T2/T3=T5/T6 
             <=> T1*T6=T3*T4 and T2*T6=T3*T5
     */
    /* T1 = buf[0], ..., T7 = buf[6] */
    /* T7:=(T1*T6) mod N; */
    mpres_mul(E->buf[6], P->x, Q->z, n);
    /* T1:=(T1*T5) mod N; */
    mpres_mul(E->buf[0], P->x, Q->y, n);
    /* T5:=(T3*T5) mod N; */
    mpres_mul(E->buf[4], P->z, Q->y, n);
    /* T3:=(T3*T4) mod N; */
    mpres_mul(E->buf[2], P->z, Q->x, n);
    /* T4:=(T2*T4) mod N; */
    mpres_mul(E->buf[3], P->y, Q->x, n);
    /* T2:=(T2*T6) mod N; */
    mpres_mul(E->buf[1], P->y, Q->z, n);

    if(mpres_equal(E->buf[6], E->buf[2], n)
       && mpres_equal(E->buf[4], E->buf[1], n))
	/* as a matter of that, P = Q and we need duplicate */
	return hessian_duplicate(R, P, E, n);

    /* T6:=(T2*T7) mod N; */
    mpres_mul(E->buf[5], E->buf[1], E->buf[6], n);
    /* T2:=(T2*T4) mod N; */
    mpres_mul(E->buf[1], E->buf[1], E->buf[3], n);
    /* T4:=(T3*T4) mod N; */
    mpres_mul(E->buf[3], E->buf[2], E->buf[3], n);
    /* T3:=(T3*T5) mod N; */
    mpres_mul(E->buf[2], E->buf[2], E->buf[4], n);
    /* T5:=(T1*T5) mod N; */
    mpres_mul(E->buf[4], E->buf[0], E->buf[4], n);
    /* T1:=(T1*T7) mod N; */
    mpres_mul(E->buf[0], E->buf[0], E->buf[6], n);
    /* T1:=(T1-T4) mod N; */
    mpres_sub(R->y, E->buf[0], E->buf[3], n);
    /* T2:=(T2-T5) mod N; */
    mpres_sub(R->x, E->buf[1], E->buf[4], n);
    /* T3:=(T3-T6) mod N; */
    mpres_sub(R->z, E->buf[2], E->buf[5], n);
    /* return [T2, T1, T3]; */
    return 1;
}

int
hessian_add(ell_point_t R, ell_point_t P, ell_point_t Q, ell_curve_t E, mpmod_t n)
{
    if(hessian_is_zero(P, E, n)){
	ell_point_set(R, Q, E, n);
	return 1;
    }
    else if(hessian_is_zero(Q, E, n)){
	ell_point_set(R, P, E, n);
	return 1;
    }
    else
	return hessian_plus(R, P, Q, E, n);
}

#if USE_ADD_SUB_CHAINS > 0
int
hessian_sub(ell_point_t R, ell_point_t P, ell_point_t Q, ell_curve_t E, mpmod_t n)
{
    int ret;

    hessian_negate(Q, E, n);
    ret = hessian_add(R, P, Q, E, n);
    hessian_negate(Q, E, n);
    return ret;
}
#endif

/* switch from X^3+Y^3+1=3*D*X*Y to Y^2=X^3+A*X+B
   A:=-27*D*(D^3+8);
   B:=54*(D^6-20*D^3-8);
   xi:=12*(D^3-1)/(D*u+v+1);
   x:=-9*D^2+xi*u;
   y:=3*xi*(v-1);
   OUTPUT: If a factor is found during the inversion, it is put in f and
   ECM_FACTOR_FOUND_STEP1 is returned. Otherwise, ECM_NO_FACTOR_FOUND is
   returned.
   SIDE-EFFECT: (x, y, D) <- (x_on_W, y_on_W, A_of_W)
 */
int
hessian_to_weierstrass(mpz_t f, mpres_t x, mpres_t y, mpres_t D, mpmod_t n)
{
    mpres_t D3, A, xi, tmp1, tmp2;
    int ret = ECM_NO_FACTOR_FOUND;

#if DEBUG_ADD_LAWS >= 1
    printf("P:=[");
    print_mpz_from_mpres(x, n);
    printf(", ");
    print_mpz_from_mpres(y, n);
    printf(", 1];\n");
    printf("D:=");
    print_mpz_from_mpres(D, n);
    printf(";\n");
#endif
    /* D3 <- D^3 */
    mpres_init(D3, n);
    mpres_mul(D3, D, D, n);
    mpres_mul(D3, D3, D, n);
    /* finish A */
    mpres_init(A, n);
    mpres_add_ui(A, D3, 8, n);
    mpres_mul(A, A, D, n);
    mpres_mul_ui(A, A, 27, n);
    mpres_neg(A, A, n);
    /* compute xi */
    mpres_init(xi, n);
    mpres_init(tmp1, n);
    mpres_mul(tmp1, D, x, n);
    mpres_add(tmp1, tmp1, y, n);
    mpres_add_ui(tmp1, tmp1, 1, n);
    mpres_init(tmp2, n);
    mpres_sub_ui(tmp2, D3, 1, n);
    mpres_mul_ui(tmp2, tmp2, 12, n);
    if(mpres_invert(xi, tmp1, n) == 0){
	mpres_gcd(f, tmp1, n);
	ret = ECM_FACTOR_FOUND_STEP1;
    }
    else{
	mpres_mul(xi, xi, tmp2, n);
	/* compute x */
	mpres_mul(tmp1, D, D, n);
	mpres_mul_ui(tmp1, tmp1, 9, n);
	mpres_mul(tmp2, xi, x, n);
	mpres_sub(x, tmp2, tmp1, n);
	/* compute y */
	mpres_sub_ui(tmp1, y, 1, n);
	mpres_mul(tmp1, tmp1, xi, n);
	mpres_mul_ui(y, tmp1, 3, n);
	mpres_set(D, A, n);
#if DEBUG_ADD_LAWS >= 1
	printf("WP:=[");
	print_mpz_from_mpres(x, n);
	printf(", ");
	print_mpz_from_mpres(y, n);
	printf(", 1];\n");
	printf("WA:=");
	print_mpz_from_mpres(D, n);
	printf(";\nWB:=(WP[2]^2-WP[1]^3-WA*WP[1]) mod N;WE:=[WA, WB];\n");
#endif
    }
    mpres_clear(A, n);
    mpres_clear(D3, n);
    mpres_clear(xi, n);
    mpres_clear(tmp1, n);
    mpres_clear(tmp2, n);
    return ret;
}

int
mult_by_3(mpz_t f, mpres_t x, mpres_t y, mpres_t A, mpmod_t n)
{
    ell_curve_t E;
    ell_point_t P, Q;
    int ret = ECM_NO_FACTOR_FOUND;
    mpz_t e;

    ell_curve_init_set(E, ECM_EC_TYPE_WEIERSTRASS, ECM_LAW_AFFINE, A, n);
    ell_point_init(P, E, n);
    mpres_set(P->x, x, n);
    mpres_set(P->y, y, n);
    mpres_set_ui(P->z, 1, n);
    ell_point_init(Q, E, n);
    mpz_init_set_ui(e, 3);
    if(ell_point_mul(f, Q, e, P, E, n) != 0){
	mpres_set(x, Q->x, n);
	mpres_set(y, Q->y, n);
    }
    mpz_clear(e);
    ell_point_clear(Q, E, n);
    ell_point_clear(P, E, n);
    ell_curve_clear(E, n);
    return ret;
}

/******************** projective twisted Hessian form ********************/

/* a*U^3+V^3+W^3 = d*U*V*W
   O_E = [0:-1:1]
   -[U:V:W]=[U:W:V]
*/
int
twisted_hessian_is_zero(ell_point_t P, ATTRIBUTE_UNUSED ell_curve_t E, mpmod_t n)
{
    mpres_t tmp;
    int ret;

    if(mpz_sgn(P->x) != 0)
	return 0;
    mpres_init(tmp, n);
    mpres_add(tmp, P->y, P->z, n);
    ret = mpz_sgn(tmp) == 0;
#if 0
    if(ret)
	gmp_printf("found a third root of unity? %Zd/%Zd\n", P->x, P->y);
#endif
    mpres_clear(tmp, n);
    return ret;
}

void
twisted_hessian_set_to_zero(ell_point_t P, ATTRIBUTE_UNUSED ell_curve_t E, mpmod_t n)
{
    mpres_set_si(P->x,  0, n);
    mpres_set_si(P->y, -1, n);
    mpres_set_si(P->z,  1, n);
}

#if DEBUG_ADD_LAWS >= 1
void
twisted_hessian_print(ell_point_t P, ell_curve_t E, mpmod_t n)
{
    pt_w_print(P->x, P->y, P->z, E, n);
}
#endif

#if USE_ADD_SUB_CHAINS > 0
/* -[u:v:w] = [u:w:v] */
void
twisted_hessian_negate(ell_point_t P, ATTRIBUTE_UNUSED ell_curve_t E, ATTRIBUTE_UNUSED mpmod_t n)
{
    mpz_swap(P->y, P->z); /* humf */
}
#endif

/* TODO: decrease the number of buffers? */
/* 6M+2S+1M_d: better when d is small */
int
twisted_hessian_duplicate(ell_point_t R, ell_point_t P, 
		  ATTRIBUTE_UNUSED ell_curve_t E, mpmod_t n)
{
    /* R = buf[0], ..., W = buf[5], C = buf[6], D = buf[7], E = buf[8] */
    /*  R:=Y1+Z1;*/
    mpres_add(E->buf[0], P->y, P->z, n);
    /*	S:=Y1-Z1;*/
    mpres_sub(E->buf[1], P->y, P->z, n);
    /*	T:=R^2 mod N;*/
    mpres_sqr(E->buf[2], E->buf[0], n);
    /*	U:=S^2 mod N;*/
    mpres_sqr(E->buf[3], E->buf[1], n);
    /*	V:=T+3*U;*/
    mpres_add(E->buf[4], E->buf[2], E->buf[3], n);
    mpres_add(E->buf[4], E->buf[4], E->buf[3], n);
    mpres_add(E->buf[4], E->buf[4], E->buf[3], n);
    /*	W:=3*T+U;*/
    mpres_add(E->buf[5], E->buf[3], E->buf[2], n);
    mpres_add(E->buf[5], E->buf[5], E->buf[2], n);
    mpres_add(E->buf[5], E->buf[5], E->buf[2], n);
    /*	C:=(R*V) mod N;*/
    mpres_mul(E->buf[6], E->buf[0], E->buf[4], n);
    /*	D:=(S*W) mod N;*/
    mpres_mul(E->buf[7], E->buf[1], E->buf[5], n);
    /*	E:=(3*C-E0[2]*X1*(W-V)) mod N;*/
    mpres_sub(E->buf[8], E->buf[5], E->buf[4], n);
    mpres_mul(E->buf[8], E->buf[8], P->x, n);
    mpres_mul(E->buf[8], E->buf[8], E->a6, n);
    mpres_sub(E->buf[8], E->buf[6], E->buf[8], n);
    mpres_add(E->buf[8], E->buf[8], E->buf[6], n);
    mpres_add(E->buf[8], E->buf[8], E->buf[6], n);
    /*	X3:=(-2*X1*D) mod N;*/
    mpres_mul(R->x, P->x, E->buf[7], n);
    mpres_add(R->x, R->x, R->x, n);
    mpres_neg(R->x, R->x, n);
    /*	Y3:=((D+E)*Z1) mod N;*/
    mpres_add(E->buf[0], E->buf[7], E->buf[8], n);
    mpres_mul(E->buf[1], E->buf[0], P->z, n);
    /*	Z3:=((D-E)*Y1) mod N;*/
    mpres_sub(E->buf[0], E->buf[7], E->buf[8], n);
    mpres_mul(R->z, E->buf[0], P->y, n);
    mpres_set(R->y, E->buf[1], n);
    return 1;
}

/* TODO: reduce the number of buffers? */
int
twisted_hessian_plus(ell_point_t R, ell_point_t P, ell_point_t Q, 
	     ATTRIBUTE_UNUSED ell_curve_t E, mpmod_t n)
{
    /* A = buf[0], ... F = buf[5], G = [6], H = [7], J = [8] */
    // A:=X1*Z2 mod N;
    mpres_mul(E->buf[0], P->x, Q->z, n);
    // B:=Z1*Z2 mod N;
    mpres_mul(E->buf[1], P->z, Q->z, n);
    // C:=Y1*X2 mod N;
    mpres_mul(E->buf[2], P->y, Q->x, n);
    // D:=Y1*Y2 mod N;
    mpres_mul(E->buf[3], P->y, Q->y, n);
    // E:=Z1*Y2 mod N;
    mpres_mul(E->buf[4], P->z, Q->y, n);
    // F:=E0[1]*X1*X2 mod N;
    mpres_mul(E->buf[5], P->x, Q->x, n);
    mpres_mul(E->buf[5], E->buf[5], E->a4, n);
    // Hisil
    // G := (D+B)*(A-C) mod N;
    mpres_add(E->buf[9], E->buf[3], E->buf[1], n);
    mpres_sub(E->buf[6], E->buf[0], E->buf[2], n);
    mpres_mul(E->buf[6], E->buf[6], E->buf[9], n);
    // H := (D-B)*(A+C) mod N;
    mpres_sub(E->buf[9], E->buf[3], E->buf[1], n);
    mpres_add(E->buf[7], E->buf[0], E->buf[2], n);
    mpres_mul(E->buf[7], E->buf[7], E->buf[9], n);
    // J := (D+F)*(A-E) mod N;
    mpres_add(E->buf[9], E->buf[3], E->buf[5], n);
    mpres_sub(E->buf[8], E->buf[0], E->buf[4], n);
    mpres_mul(E->buf[8], E->buf[8], E->buf[9], n);
    // K := (D-F)*(A+E) mod N;
    // this is the last use of A, so that K -> buf[0]
    mpres_sub(E->buf[9], E->buf[3], E->buf[5], n);
    mpres_add(E->buf[0], E->buf[0], E->buf[4], n);
    mpres_mul(E->buf[0], E->buf[0], E->buf[9], n);
    // X3 := G-H
    mpres_sub(R->x, E->buf[6], E->buf[7], n);
    // Y3 := K-J
    mpres_sub(R->y, E->buf[0], E->buf[8], n);
    // Z3 := (J+K-G-H-2*(B-F)*(C+E)) mod N;
    mpres_sub(E->buf[9], E->buf[1], E->buf[5], n);
    mpres_add(R->z, E->buf[2], E->buf[4], n);
    mpres_mul(R->z, R->z, E->buf[9], n);
    mpres_add(R->z, R->z, R->z, n);
    mpres_add(R->z, R->z, E->buf[7], n);
    mpres_add(R->z, R->z, E->buf[6], n);
    mpres_sub(R->z, E->buf[0], R->z, n);
    mpres_add(R->z, R->z, E->buf[8], n);
    if(mpz_sgn(R->x) == 0 && mpz_sgn(R->y) == 0 && mpz_sgn(R->z) == 0){
	// iff (X2:Y2:Z2)=(Z1:gamma^2*X1:gamma*Y1), gamma^3 = a
	fprintf(stderr, "GASP: X3, Y3 and Z3 are 0\n");
	exit(-1);
#if 0
	    // TODO: rewrite with above quantities!
	    X3:=(X1^2*Y2*Z2-X2^2*Y1*Z1) mod N;
	    // A*X1*Y2-C*X2*Z1 = A*U-C*V
	    Y3:=(Z1^2*X2*Y2-Z2^2*X1*Y1) mod N;
	    // E*Z1*X2-A*Z2*Y1 = E*V-A*W
	    Z3:=(Y1^2*X2*Z2-Y2^2*X1*Z1) mod N;
	    // C*Y1*Z2-E*Y2*X1 = C*W-E*U

	    // X3 =     Y1*(a*X1^3-Z1^3)
	    // Y3 = g^2*X1*(Z1^3-Y1^3)
	    // Z3 =   g*Z1*(Y1^3-Z1^3)

	    // can be made faster with a = aa^3, since then g = aa and we
	    // can share many things

#endif
    }
    return 1;
}

int
twisted_hessian_add(ell_point_t R, ell_point_t P, ell_point_t Q, ell_curve_t E, mpmod_t n)
{
    if(twisted_hessian_is_zero(P, E, n)){
	ell_point_set(R, Q, E, n);
	return 1;
    }
    else if(twisted_hessian_is_zero(Q, E, n)){
	ell_point_set(R, P, E, n);
	return 1;
    }
    else
	return twisted_hessian_plus(R, P, Q, E, n);
}

#if USE_ADD_SUB_CHAINS > 0
int
twisted_hessian_sub(ell_point_t R, ell_point_t P, ell_point_t Q, ell_curve_t E, mpmod_t n)
{
    int ret;

    twisted_hessian_negate(Q, E, n);
    ret = twisted_hessian_add(R, P, Q, E, n);
    twisted_hessian_negate(Q, E, n);
    return ret;
}
#endif

/* INPUT: a*x^3+y^3+1 = d*x*y
   OUTPUT: Y^2 = X^3+A*X+B
   If a=c^3, then curve isom to Hessian (c*x)^3+y^3+1=3*(d/(3*c))*(c*x)*y
   SIDE EFFECT: (x, y, c) <- (x_on_W, y_on_W, A_of_W)
 */
int
twisted_hessian_to_weierstrass(mpz_t f, mpres_t x, mpres_t y, mpres_t c, mpres_t d, mpmod_t n)
{
    int ret = ECM_NO_FACTOR_FOUND;
    mpres_t tmp;

#if DEBUG_ADD_LAWS >= 2
    printf("x_tH="); print_mpz_from_mpres(x, n); printf("\n");
    printf("y_tH="); print_mpz_from_mpres(y, n); printf("\n");
    printf("c_tH="); print_mpz_from_mpres(c, n); printf("\n");
    printf("d_tH="); print_mpz_from_mpres(d, n); printf("\n");
#endif
    mpres_init(tmp, n);
    mpres_mul_ui(tmp, c, 3, n);
    if(mpres_invert(tmp, tmp, n) == 0){
        mpres_gcd(f, tmp, n);
        ret = ECM_FACTOR_FOUND_STEP1;
    }
    else{
	mpres_mul(x, x, c, n);
	mpres_mul(c, tmp, d, n);
	/* from x^3+y^3+1=3*c*x*y to Weierstrass stuff */
	ret = hessian_to_weierstrass(f, x, y, c, n);
#if DEBUG_ADD_LAWS >= 2
	printf("A_W="); print_mpz_from_mpres(c, n); printf("\n");
	printf("x_W="); print_mpz_from_mpres(x, n); printf("\n");
	printf("y_W="); print_mpz_from_mpres(y, n); printf("\n");
#endif
    }
    mpres_clear(tmp, n);
    return ret;
}

/******************** generic ec's ********************/

void
ell_point_init(ell_point_t P, ell_curve_t E, mpmod_t n)
{
    mpres_init(P->x, n);
    mpres_init(P->y, n);
    mpres_init(P->z, n);
    if(E->type == ECM_EC_TYPE_WEIERSTRASS){
      if(E->law == ECM_LAW_AFFINE)
	mpz_set_ui(P->z, 1); /* humf */
      else if(E->law == ECM_LAW_HOMOGENEOUS)
	mpres_set_ui(P->z, 1, n);
    }
    else if(E->type == ECM_EC_TYPE_HESSIAN 
	    || E->type == ECM_EC_TYPE_TWISTED_HESSIAN)
	mpres_set_ui(P->z, 1, n);
}

/* TODO: change this according to E->type */
void
ell_point_clear(ell_point_t P, ATTRIBUTE_UNUSED ell_curve_t E, mpmod_t n)
{
    mpres_clear(P->x, n);
    mpres_clear(P->y, n);
    mpres_clear(P->z, n);
}

#if DEBUG_ADD_LAWS >= 1
void
ell_point_print(ell_point_t P, ell_curve_t E, mpmod_t n)
{
    if(E->type == ECM_EC_TYPE_WEIERSTRASS)
	pt_w_print(P->x, P->y, P->z, E, n);
    else if(E->type == ECM_EC_TYPE_HESSIAN)
	hessian_print(P, E, n);
    else if(E->type == ECM_EC_TYPE_TWISTED_HESSIAN)
	twisted_hessian_print(P, E, n);
}
#endif

/* TODO: should depend on E->type... */
void
ell_point_set(ell_point_t Q, ell_point_t P,
	     ATTRIBUTE_UNUSED ell_curve_t E, ATTRIBUTE_UNUSED mpmod_t n)
{
    mpres_set(Q->x, P->x, n);
    mpres_set(Q->y, P->y, n);
    mpres_set(Q->z, P->z, n);
}

void
ell_curve_init(ell_curve_t E, int etype, int law, mpmod_t n)
{
    int i;

    E->type = etype;
    E->law = law;
    mpres_init(E->a1, n);
    mpres_init(E->a3, n);
    mpres_init(E->a2, n);
    mpres_init(E->a4, n);
    mpres_init(E->a6, n);
    mpres_set_ui(E->a1, 0, n);
    mpres_set_ui(E->a3, 0, n);
    mpres_set_ui(E->a2, 0, n);
    mpres_set_ui(E->a4, 0, n);
    mpres_set_ui(E->a6, 0, n);
    for(i = 0; i < EC_W_NBUFS; i++)
	mpres_init (E->buf[i], n);
}

void
ell_curve_init_set(ell_curve_t E, int etype, int law, mpres_t A, mpmod_t n)
{
    ell_curve_init(E, etype, law, n);
    mpres_set(E->a4, A, n);
}

void
ell_curve_set_z(ell_curve_t E, ell_curve_t zE, mpmod_t n)
{
    ell_curve_init(E, zE->type, zE->law, n);
    mpres_set_z(E->a1, zE->a1, n);
    mpres_set_z(E->a3, zE->a3, n);
    mpres_set_z(E->a2, zE->a2, n);
    mpres_set_z(E->a4, zE->a4, n);
    mpres_set_z(E->a6, zE->a6, n);
#if 0
    E->disc = zE->disc;
    if(E->disc != 0){
	mpres_init(E->sq[0], n);
	mpres_set_z(E->sq[0], zE->sq[0], n);
    }
#endif
}

void
ell_curve_clear(ell_curve_t E, mpmod_t n)
{
    int i;

    mpres_clear(E->a4, n);
    for(i = 0; i < EC_W_NBUFS; i++)
	mpres_clear (E->buf[i], n);
    /* TODO: case of sq */
}

#if DEBUG_ADD_LAWS >= 1
void
ell_curve_print(ell_curve_t E, mpmod_t n)
{
    if(E->type == ECM_EC_TYPE_WEIERSTRASS){
	printf("["); print_mpz_from_mpres(E->a1, n);
	printf(", "); print_mpz_from_mpres(E->a3, n);
	printf(", "); print_mpz_from_mpres(E->a2, n);
	printf(", "); print_mpz_from_mpres(E->a4, n);
	printf(", "); print_mpz_from_mpres(E->a6, n); printf("];\n");
    }
    else if(E->type == ECM_EC_TYPE_HESSIAN){
	printf("D:="); print_mpz_from_mpres(E->a4, n); printf(";\n");
	printf("E:=[D];\n");
    }
    else if(E->type == ECM_EC_TYPE_TWISTED_HESSIAN){
	printf("a:="); print_mpz_from_mpres(E->a4, n); printf(";\n");
	printf("d:="); print_mpz_from_mpres(E->a6, n); printf(";\n");
	printf("E:=[a, d];\n");
    }
}
#endif

/* OUTPUT: 1 if P = O_E, 0 otherwise. */
int
ell_point_is_zero(ell_point_t P, ell_curve_t E, mpmod_t n)
{
    if(E->type == ECM_EC_TYPE_WEIERSTRASS)
	return pt_w_is_zero(P->z, n);
    else if(E->type == ECM_EC_TYPE_HESSIAN)
	return hessian_is_zero(P, E, n);
    else if(E->type == ECM_EC_TYPE_TWISTED_HESSIAN)
	return twisted_hessian_is_zero(P, E, n);
    return ECM_ERROR;
}

void
ell_point_set_to_zero(ell_point_t P, ell_curve_t E, mpmod_t n)
{
    if(E->type == ECM_EC_TYPE_WEIERSTRASS)
	pt_w_set_to_zero(P, n);
    else if(E->type == ECM_EC_TYPE_HESSIAN)
	hessian_set_to_zero(P, E, n);
    else if(E->type == ECM_EC_TYPE_TWISTED_HESSIAN)
	twisted_hessian_set_to_zero(P, E, n);
}

int
ell_point_is_on_curve(ell_point_t P, ell_curve_t E, mpmod_t n)
{
    int ok = 1;

    if(ell_point_is_zero(P, E, n))
	return 1;
    if(E->type == ECM_EC_TYPE_WEIERSTRASS){
	mpres_t tmp1, tmp2;

	mpres_init(tmp1, n);
	mpres_init(tmp2, n);
	if(E->law == ECM_LAW_AFFINE){
	    /* y^2+a1*x*y+a3*y = x^3+a2*x^2+a4*x+a6? */
	    mpres_mul(tmp1, E->a1, P->x, n);
	    mpres_add(tmp1, tmp1, P->y, n);
	    mpres_add(tmp1, tmp1, E->a3, n);
	    mpres_mul(tmp1, tmp1, P->y, n);
	    
	    mpres_add(tmp2, E->a2, P->x, n);
	    mpres_mul(tmp2, tmp2, P->x, n);
	    mpres_add(tmp2, tmp2, E->a4, n);
	    mpres_mul(tmp2, tmp2, P->x, n);
	    mpres_add(tmp2, tmp2, E->a6, n);
	}
#if 0 // useless for the time being
	else{
	    /* y^2*z+a1*x*y*z+a3*y*z^2 = x^3+a2*x^2*z+a4*x*z^2+a6*z^3? */
	    /* y*z*(y+a1*x+a3*z) = ((x+a2*z)*x+a4*z^2)*x+a6*z^3? */
	    mpres_t tmp3;

	    mpres_mul(tmp1, E->a1, P->x, n);  /* a1*x */
	    mpres_add(tmp1, tmp1, P->y, n);   /* a1*x+y */
	    mpres_mul(tmp2, E->a3, P->z, n);  /* a3*z */
	    mpres_add(tmp1, tmp1, tmp2, n);   /* y+a1*x+a3*z */
	    mpres_mul(tmp1, tmp1, P->y, n);   /* y*(...) */
	    mpres_mul(tmp1, tmp1, P->z, n);   /* lhs */

	    mpres_init(tmp3, n);
	    mpres_mul(tmp2, E->a2, P->z, n);  /* a2*z */
	    mpres_add(tmp2, tmp2, P->x, n);   /* x+a2*z */
	    mpres_mul(tmp2, tmp2, P->x, n);   /* (x+a2*z)*x */
	    mpres_mul(tmp3, E->a4, P->z, n);  /* a4*z */
	    mpres_mul(tmp3, tmp3, P->z, n);   /* a4*z^2 */
	    mpres_add(tmp2, tmp2, tmp3, n);   /* (x+a2*z)*x+a4*z^2 */
	    mpres_mul(tmp2, tmp2, P->x, n);   /* (...)*x */
	    mpres_mul(tmp3, P->z, P->z, n);   /* z^2 */
	    mpres_mul(tmp3, tmp3, P->z, n);   /* z^3 */
	    mpres_mul(tmp3, tmp3, E->a6, n);  /* a6*z^3 */
	    mpres_add(tmp2, tmp2, tmp3, n);   /* rhs */
	    mpres_clear(tmp3, n);
	}
#endif
	ok = mpres_equal(tmp1, tmp2, n);

	mpres_clear(tmp1, n);
	mpres_clear(tmp2, n);
    }
    else if(E->type == ECM_EC_TYPE_HESSIAN){
	/* TODO */
    }
    else if(E->type == ECM_EC_TYPE_TWISTED_HESSIAN){
	/* TODO */
    }
    return ok;
}

#if DEBUG_ADD_LAWS >= 1
static void
ell_point_check(ell_point_t P, ell_curve_t E, mpmod_t n)
{
    if(ell_point_is_on_curve(P, E, n) == 0){
	printf("Point not on curve\n");
	printf("E:=");
	ell_curve_print(E, n);
	printf("P:=");
	pt_print(E, P, n);
	printf("\n");
	exit(-1);
    }
}
#endif

#if DEBUG_ADD_LAWS >= 1
int
ell_point_equal(ell_point_t P, ell_point_t Q, ell_curve_t E, mpmod_t n)
{
    int ret = 1;

    if(E->type == ECM_EC_TYPE_WEIERSTRASS){
	if(E->law == ECM_LAW_AFFINE)
	    return mpres_equal(P->x, Q->x, n) 
		   && mpres_equal(P->y, Q->y, n)
		   && mpres_equal(P->z, Q->z, n);
	else if(E->law == ECM_LAW_HOMOGENEOUS){
	    mpres_t tmp1, tmp2;

	    mpres_init(tmp1, n);
	    mpres_init(tmp2, n);
	    mpres_mul(tmp1, P->x, Q->z, n);
	    mpres_mul(tmp2, P->z, Q->x, n);
	    if(mpres_equal(tmp1, tmp2, n) == 0){
		printf("Px/Pz != Qx/Qz\n");
		ret = 0;
		exit(-1);
	    }
	    else{
		mpres_mul(tmp1, P->y, Q->z, n);
		mpres_mul(tmp2, P->z, Q->y, n);
		if(mpres_equal(tmp1, tmp2, n) == 0){
		    printf("Py/Pz != Qy/Qz\n");
		    ret = 0;
		    exit(-1);
		}
	    }
	    mpres_clear(tmp1, n);
	    mpres_clear(tmp2, n);
	}
    }
    return ret;
}
#endif

/* OUTPUT: 1 if everything ok, 0 otherwise */
int
ell_point_add(mpz_t f, ell_point_t R, ell_point_t P, ell_point_t Q, ell_curve_t E, mpmod_t n)
{
    if(E->type == ECM_EC_TYPE_WEIERSTRASS)
	return pt_w_add(f, R->x, R->y, R->z, P->x, P->y, P->z, 
			Q->x, Q->y, Q->z, n, E);
    else if(E->type == ECM_EC_TYPE_HESSIAN)
	return hessian_add(R, P, Q, E, n);
    else if(E->type == ECM_EC_TYPE_TWISTED_HESSIAN)
	return twisted_hessian_add(R, P, Q, E, n);
    else
	return ECM_ERROR;
}

#if USE_ADD_SUB_CHAINS > 0
/* R <- P-Q */
int
ell_point_sub(mpz_t f, ell_point_t R, ell_point_t P, ell_point_t Q, ell_curve_t E, mpmod_t n)
{
    if(E->type == ECM_EC_TYPE_WEIERSTRASS)
	return pt_w_sub(f, R->x, R->y, R->z, P->x, P->y, P->z,
			Q->x, Q->y, Q->z, n, E);
    else if(E->type == ECM_EC_TYPE_HESSIAN)
	return hessian_sub(R, P, Q, E, n);
    else if(E->type == ECM_EC_TYPE_TWISTED_HESSIAN)
	return twisted_hessian_sub(R, P, Q, E, n);
    else
	return ECM_ERROR;
}
#endif

int
ell_point_duplicate(mpz_t f, ell_point_t R, ell_point_t P, ell_curve_t E, mpmod_t n)
{
#if DEBUG_ADD_LAWS >= 2
    printf("E:=");
    ell_curve_print(E, n);
#endif
    if(E->type == ECM_EC_TYPE_WEIERSTRASS)
	return pt_w_duplicate(f, R->x, R->y, R->z, P->x, P->y, P->z, n, E);
    else if(E->type == ECM_EC_TYPE_HESSIAN)
	return hessian_duplicate(R, P, E, n);
    else if(E->type == ECM_EC_TYPE_TWISTED_HESSIAN)
	return twisted_hessian_duplicate(R, P, E, n);
    else
	return ECM_ERROR;
}

void
ell_point_negate(ell_point_t P, ell_curve_t E, mpmod_t n)
{
#if DEBUG_ADD_LAWS >= 2
    printf("P:="); ell_point_print(P, E, n); printf(";\n");
#endif
    if(ell_point_is_zero(P, E, n) == 0){
	if(E->type == ECM_EC_TYPE_WEIERSTRASS){
	    if(E->law == ECM_LAW_HOMOGENEOUS){
		/* FIXME: does not work for complete equation! */
		mpres_neg(P->y, P->y, n);
	    }
	    else if(E->law == ECM_LAW_AFFINE){
		/* (-P).y = -P.y-a1*P.x-a3 */
		if(mpz_sgn(E->a1) != 0
		   || mpz_sgn(E->a3) != 0
		   || mpz_sgn(E->a2) != 0){ /* FIXME */
		    printf("GROUMF\n");
		    exit(-1);
		}
		mpres_neg(P->y, P->y, n);
	    }
	}
#if USE_ADD_SUB_CHAINS > 0
	else if(E->type == ECM_EC_TYPE_HESSIAN)
	    hessian_negate(P, E, n);
	else if(E->type == ECM_EC_TYPE_TWISTED_HESSIAN)
	    twisted_hessian_negate(P, E, n);
#endif
    }
#if DEBUG_ADD_LAWS >= 2
    printf("neg(P):="); ell_point_print(P, E, n); printf(";\n");
#endif
}

/* Q <- [e]*P
   Return value: 0 if a factor is found, and the factor is in Q->x,
                 1 otherwise.
*/
int
ell_point_mul_plain (mpz_t f, ell_point_t Q, mpz_t e, ell_point_t P, ell_curve_t E, mpmod_t n)
{
  size_t l;
  int negated = 0, status = 1;
  ell_point_t P0;

  if(ell_point_is_zero(P, E, n) != 0){
      ell_point_set(Q, P, E, n);
      return 1;
  }

  if (mpz_sgn (e) == 0)
    {
      ell_point_set_to_zero(Q, E, n);
      return 1;
    }

  if (mpz_sgn (e) < 0)
    {
      negated = 1;
      mpz_neg (e, e);
      ell_point_negate(P, E, n); /* since the point is non-zero */
    }

  if (mpz_cmp_ui (e, 1) == 0){
      ell_point_set(Q, P, E, n);
      goto ell_point_mul_plain_end;
  }

  l = mpz_sizeinbase (e, 2) - 1; /* l >= 1 */

  ell_point_init(P0, E, n);
  ell_point_set(P0, P, E, n);

#if DEBUG_ADD_LAWS >= 2
  printf("P:="); ell_point_print(P, E, n); printf(";\n");
#endif
  while (l-- > 0)
    {
#if DEBUG_ADD_LAWS >= 2
	printf("P0:="); ell_point_print(P0, E, n); printf(";\n");
#endif
	if(ell_point_duplicate (f, P0, P0, E, n) == 0)
	  {
	    status = 0;
	    break;
	  }
#if DEBUG_ADD_LAWS >= 2
	printf("Rdup:="); ell_point_print(P0, E, n); printf(";\n");
	printf("dup:=ProjEcmDouble(P0, E, N); ProjEcmEqual(dup, Rdup, N);\n");
#endif
	if (mpz_tstbit (e, l))
	  {
	      if(ell_point_add (f, P0, P0, P, E, n) == 0)
	      {
		status = 0;
		break;
	      }
#if DEBUG_ADD_LAWS >= 2
	      printf("Radd:="); ell_point_print(P0, E, n); printf(";\n");
	      printf("Padd:=ProjEcmAdd(P, Rdup, E, N); ProjEcmEqual(Padd, Radd, N);\n");
#endif
	  }
    }

  ell_point_set(Q, P0, E, n);
  ell_point_clear(P0, E, n);
ell_point_mul_plain_end:

  /* Undo negation to avoid changing the caller's e value */
  if (negated){
    mpz_neg (e, e);
    ell_point_negate(P, E, n);
  }
  return status;
}

int
ell_point_mul(mpz_t f, ell_point_t Q, mpz_t e, ell_point_t P, ell_curve_t E, mpmod_t n)
{
#if 1 /* keeping it simple */
    return ell_point_mul_plain(f, Q, e, P, E, n);
#else
    return ell_point_mul_add_sub(f, Q, e, P, E, n);
#endif
}

