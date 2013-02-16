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

#define DEBUG_ADD_LAWS 0

/********** group law on points **********/

int
pt_is_zero(ell_point_t P, ATTRIBUTE_UNUSED mpmod_t n)
{
    return mpz_sgn(P->z) == 0;
}

void
pt_set_to_zero(ell_point_t P, mpmod_t n)
{
    mpz_set_ui(P->x, 0);
    mpres_set_ui(P->y, 1, n);
    mpz_set_ui(P->z, 0);
}

void
pt_assign(ell_point_t Q, ell_point_t P, ATTRIBUTE_UNUSED mpmod_t n)
{
    mpres_set(Q->x, P->x, n);
    mpres_set(Q->y, P->y, n);
    mpres_set(Q->z, P->z, n);
}

void
pt_neg(ell_point_t P, mpmod_t n)
{
    if(pt_is_zero(P, n) == 0)
	mpres_neg(P->y, P->y, n);
}

void
pt_many_set_to_zero(ell_point_t *tP, int nE, mpmod_t n)
{
    int i;

    for(i = 0; i < nE; i++)
	pt_set_to_zero(tP[i], n);
}

void
pt_many_neg(ell_point_t *tP, int nE, mpmod_t n)
{
    int i;

    for(i = 0; i < nE; i++)
	pt_neg(tP[i], n);
}

void
pt_many_assign(ell_point_t *tQ, ell_point_t *tP, int nE, mpmod_t n)
{
    int i;

    for(i = 0; i < nE; i++)
	pt_assign(tQ[i], tP[i], n);
}

void
print_mpz_from_mpres(mpres_t x, mpmod_t n)
{
    mpz_t tmp;

    mpz_init(tmp);
    mpres_get_z(tmp, x, n);
    gmp_printf("%Zd", tmp);
    mpz_clear(tmp);
}

void
pt_print(ell_point_t P, mpmod_t n)
{
    printf("[");
    print_mpz_from_mpres(P->x, n);
    printf(", ");
    print_mpz_from_mpres(P->y, n);
    printf(", ");
    print_mpz_from_mpres(P->z, n);
    printf("]");
}

void
pt_many_print(ell_curve_t *tE, ell_point_t *tP, int nE, mpmod_t n)
{
    int i;

    for(i = 0; i < nE; i++){
	printf("%d: ", i);
	pt_print(tP[i], n);
	printf(" on E.A=");
	print_mpz_from_mpres(tE[i]->A, n);
	printf("\n");
    }
}

/* Computes inv[i] = 1/x[i] using only one inversion, a la Montgomery.
   If takeit[i] != 1, do not compute 1/x[i] (it is probably 0, or irrelevant).
   We should have inv != x.
   x[nx] is a buffer.
   When a factor is found, the i s.t. x[i] is not invertible are looked for
   and the corresponding values of takeit put to 2.
*/
int
compute_all_inverses(mpres_t *inv, mpres_t *x, int nx, mpmod_t n, char *takeit)
{
    int i;

#if 0
    /* plain version, to debug the architecture */
    for(i = 0; i < nx; i++){
	if(takeit[i] != 1)
	    continue;
	if(!mpres_invert(inv[i], x[i], n)){
	    mpres_gcd(inv[0], x[i], n);
#if DEBUG_ADD_LAWS >= 1
	    printf("Factor[%d]: ", i);
            mpz_out_str (stdout, 10, inv[0]);
            printf ("\n");
#endif
	    return 0;
	}
    }
#else
    /* Montgomery's trick */
    for(i = 0; i < nx; i++){
	if(takeit[i] != 1){
	    if(i == 0)
		mpres_set_ui(inv[i], 1, n);
	    else
		mpres_set(inv[i], inv[i-1], n);
	}
	else{
	    if(i == 0)
		mpres_set(inv[i], x[i], n);
	    else
		mpres_mul(inv[i], inv[i-1], x[i], n);
	}
    }
    /* invert */
    if(!mpres_invert(x[nx], inv[nx-1], n)){
	mpres_gcd(inv[0], inv[nx-1], n);
#if DEBUG_ADD_LAWS >= 1
	printf("Factor[%d]: ", i);
	mpz_out_str (stdout, 10, inv[0]);
	printf ("\n");
#endif
	/* identifying the x[i]'s */
	for(i = 0; i < nx; i++){
	    mpres_gcd(x[nx], x[i], n);
	    if(mpz_cmp_ui(x[nx], 1) != 0){
#if DEBUG_ADD_LAWS >= 0
		printf("# x[%d] not invertible: ", i);
		mpz_out_str (stdout, 10, x[nx]);
		printf ("\n");
#endif
		/* ONE DAY: if x[nx] != inv[0], we have another factor! */
		takeit[i] = 2;
	    }
	}
	return 0;
    }
    /* get inverses back */
    /* say inv = 1/(x1*x2*x3) */
    for(i = nx-1; i > 0; i--)
	if(takeit[i] == 1){
	    mpres_mul(inv[i], x[nx], inv[i-1], n); /* 1/x3 = inv * (x1*x2) */
	    mpres_mul(x[nx], x[nx], x[i], n); /* inv = 1/(x1*x2) */
	}
    mpres_set(inv[0], x[nx], n);
#endif
#if DEBUG_ADD_LAWS >= 1
    /*    printf("# checking inverses\n"); */
    mpres_t tmp;
    mpres_init(tmp, n);
    for(i = 0; i < nx; i++){
	mpres_mul(tmp, inv[i], x[i], n);
	mpres_get_z(tmp, tmp, n);
	if(mpz_cmp_ui(tmp, 1) != 0)
	    printf("ERROR in compute_all_inverses[%d]\n", i);
    }
    mpres_clear(tmp, n);
#endif
    return 1;
}

/* NOTE: we can have tR = tP or tQ.
   In case a factor is found, it is put in num[nE].
 */
int
pt_many_common(ell_point_t *tR, ell_point_t *tP, ell_point_t *tQ, int nE, 
	       mpmod_t n, 
	       mpres_t *num, mpres_t *den, mpres_t *inv, char *takeit)
{
    int i;

    if(compute_all_inverses(inv, den, nE, n, takeit) == 0){
	mpz_set(num[nE], inv[0]);
	return 0;
    }
    for(i = 0; i < nE; i++){
	if(takeit[i] != 1)
	    continue;
	/* l:=(inv[i]*num[i]) mod N; */
	mpres_mul(num[i], num[i], inv[i], n);
	/* x:=(l^2-P[1]-Q[1]) mod N; */
	mpres_sqr(den[i], num[i], n);
	mpres_sub(den[i], den[i], tP[i]->x, n);
	mpres_sub(den[i], den[i], tQ[i]->x, n);
	/* tR[i]:=[x, (l*(P[1]-x)-P[2]) mod N, 1]; */
	mpres_sub(tR[i]->x, tP[i]->x, den[i], n);
	mpres_mul(tR[i]->x, tR[i]->x, num[i], n);
	mpres_sub(tR[i]->y, tR[i]->x, tP[i]->y, n);
	mpres_set(tR[i]->x, den[i], n);
    }
    return 1;
}

/*   In case a factor is found, it is put in num[nE]. */
int
pt_many_duplicate(ell_point_t *tQ, ell_point_t *tP, ell_curve_t *tE, int nE, 
		  mpmod_t n, 
		  mpres_t *num, mpres_t *den, mpres_t *inv, char *ok)
{
    char *takeit = (char *)malloc(nE * sizeof(char));
    int i, res;

    memcpy(takeit, ok, nE);
    for(i = 0; i < nE; i++){
	if(ok[i] == 0)
	    continue; /* takeit[i] = 0 */
	if(pt_is_zero(tP[i], n)){
	    takeit[i] = 0;
	    pt_set_to_zero(tQ[i], n);
	}
	else if(mpz_sgn(tP[i]->y) == 0){
	    /* 2 * P[i] = O_E */
	    takeit[i] = 0;
	    pt_set_to_zero(tP[i], n);
	    printf("# [2] * P[%d] = O_E\n", i);
	}
	else{
	    mpres_sqr(num[i], tP[i]->x, n);
	    mpres_mul_ui(num[i], num[i], 3, n);
	    mpres_add(num[i], num[i], tE[i]->A, n);
	    mpres_mul_ui(den[i], tP[i]->y, 2, n);
	}
    }
    res = pt_many_common(tQ, tP, tP, nE, n, num, den, inv, takeit);
    /* TODO: case takeit[i] == 2 */
    free(takeit);
    return res;
}

/* R[i] <- P[i] + Q[i], or a factor is found which is put in num[nE]. */
int
pt_many_add(ell_point_t *tR, ell_point_t *tP, ell_point_t *tQ, ell_curve_t *tE, 
	    int nE, mpmod_t n, 
	    mpres_t *num, mpres_t *den, mpres_t *inv, char *ok)
{
    char *takeit = (char *)malloc(nE * sizeof(char));
    int i, res;

    memcpy(takeit, ok, nE);
#if DEBUG_ADD_LAWS >= 2
    printf("In pt_many_add, adding\n");
    pt_many_print(tP, nE, n);
    printf("and\n");
    pt_many_print(tQ, nE, n);
#endif
    for(i = 0; i < nE; i++){
	if(ok[i] == 0)
	    continue; /* takeit[i] = 0 */
	if(pt_is_zero(tP[i], n)){
#if DEBUG_ADD_LAWS >= 2
	    printf("# tEP[%d] = O_{E[%d]}\n", i, i);
#endif
	    takeit[i] = 0;
	    pt_assign(tR[i], tQ[i], n);
	}
	else if(pt_is_zero(tQ[i], n)){
#if DEBUG_ADD_LAWS >= 2
	    printf("# tEQ[%d] = O_{E[%d]}\n", i, i);
#endif
	    takeit[i] = 0;
	    pt_assign(tR[i], tP[i], n);
	}
	else if(pt_is_equal(tP[i], tQ[i])){
	    /* we should double */
	    if(mpz_sgn(tP[i]->y) == 0){
#if DEBUG_ADD_LAWS >= 2
		printf("# 2 * P[%d] = O_{E[%d]}\n", i, i);
#endif
		takeit[i] = 0;
		pt_set_to_zero(tR[i], n);
	    }
	    else{
		/* ordinary doubling */
		mpres_sqr(num[i], tP[i]->x, n);
		mpres_mul_ui(num[i], num[i], 3, n);
		mpres_add(num[i], num[i], tE[i]->A, n);
		mpres_mul_ui(den[i], tP[i]->y, 2, n);
	    }
	}
	else if(mpz_cmp(tQ[i]->x, tP[i]->x) == 0){
	    mpres_add(num[i], tQ[i]->x, tP[i]->x, n);
	    if(mpz_sgn(num[i]) == 0){
		takeit[i] = 0;
		pt_set_to_zero(tR[i], n);
	    }
	}
	else{
	    mpres_sub(num[i], tQ[i]->y, tP[i]->y, n);
	    mpres_sub(den[i], tQ[i]->x, tP[i]->x, n);
	}
    }
    res = pt_many_common(tR, tP, tQ, nE, n, num, den, inv, takeit);
    /* TODO: case takeit[i] == 2 */
    free(takeit);
    return res;
}

/* tER != tEP */
static int
pt_many_sub(ell_point_t *tR, ell_point_t *tQ, ell_point_t *tP, ell_curve_t *tE,
	    int nE, mpmod_t n, 
	    mpres_t *num, mpres_t *den, mpres_t *inv, char *ok)
{
    int i, res;

    for(i = 0; i < nE; i++)
	if(ok[i] == 1)
	    pt_neg(tP[i], n);
    res = pt_many_add(tR, tQ, tP, tE, nE, n, num, den, inv, ok);
    for(i = 0; i < nE; i++)
	if(ok[i] == 1)
	    pt_neg(tP[i], n);
    return res;
}

/* Ordinary binary left-right addition */
static int
pt_many_mul_plain(ell_point_t *tQ, ell_point_t *tP, ell_curve_t *tE,
		  int nE, mpz_t e, mpmod_t n, 
		  mpres_t *num, mpres_t *den, mpres_t *inv, char *ok)
{
  size_t l = mpz_sizeinbase (e, 2) - 1; /* l >= 1 */
  int status = 1;

  pt_many_assign(tQ, tP, nE, n);
  while (l-- > 0)
    {
	if(pt_many_duplicate (tQ, tQ, tE, nE, n, num, den, inv, ok) == 0)
	  {
	    status = 0;
	    break;
	  }
#if DEBUG_ADD_LAWS >= 2
	printf("Rdup:="); pt_many_print(tQ, nE, n); printf(";\n");
#endif
	if (mpz_tstbit (e, l))
	  {
	      if(pt_many_add (tQ, tP, tQ, tE, nE, n, num, den, inv, ok) == 0)
	      {
		status = 0;
		break;
	      }
#if DEBUG_ADD_LAWS >= 2
	      printf("Radd:="); pt_many_print(tQ, nE, n); printf(";\n");
#endif
	  }
    }
  return status;
}

/* Ordinary binary left-right addition; see Solinas00. Morally, we use
 w = 2. */
static int
pt_many_mul_add_sub_si(ell_point_t *tQ, ell_point_t *tP, ell_curve_t *tE, int nE,
		       long c, mpmod_t n, 
		       mpres_t *num, mpres_t *den, mpres_t *inv, char *ok)
{
    long u, S[64];
    int j, iS = 0, status = 1;
    ATTRIBUTE_UNUSED int w = 2;

    /* build NAF_w(c) */
    while(c > 0){
	if((c & 1) == 1){
	    /* c is odd */
	    u = c & (long)3;
	    if(u == 3)
		u = -1;
	}
	else
	    u = 0;
	S[iS++] = u;
	c >>= 1;
    }
    /* use it */
    pt_many_set_to_zero(tQ, nE, n);
    for(j = iS-1; j >= 0; j--){
	if(pt_many_duplicate(tQ, tQ, tE, nE, n, num, den, inv, ok) == 0){
	    status = 0;
	    break;
	}
#if DEBUG_ADD_LAWS >= 2
	printf("Rdup:="); pt_many_print(tQ, nE, n); printf(";\n");
#endif
	if(S[j] == 1){
	    if(pt_many_add(tQ, tQ, tP, tE, nE, n, num, den, inv, ok) == 0){
		status = 0;
		break;
	    }
#if DEBUG_ADD_LAWS >= 2
	    printf("Radd:="); pt_many_print(tQ, nE, n); printf(";\n");
#endif
	}
	else if(S[j] == -1){
	    if(pt_many_sub(tQ, tQ, tP, tE, nE, n, num, den, inv, ok) == 0){
		status = 0;
		break;
	    }
#if DEBUG_ADD_LAWS >= 2
	    printf("Rsub:="); pt_many_print(tQ, nE, n); printf(";\n");
#endif
	}
    }
  return status;
}

/* tEQ[i] <- e * tEP[i]; we must have tEQ != tEP */
/* If a factor is found, it is put back in num[nE]. */
int
pt_many_mul(ell_point_t *tQ, ell_point_t *tP, ell_curve_t *tE, int nE,
	    mpz_t e, mpmod_t n, 
	    mpres_t *num, mpres_t *den, mpres_t *inv, char *ok)
{
  size_t l;
  int negated = 0, status = 1;

  if (mpz_sgn (e) == 0)
    {
	pt_many_set_to_zero(tQ, nE, n);
	return 1;
    }

  /* The negative of a point (x:y:z) is (x:-y:z) */
  if (mpz_sgn (e) < 0)
    {
      negated = 1;
      mpz_neg (e, e);
      pt_many_neg(tP, nE, n);
    }

  if (mpz_cmp_ui (e, 1) == 0)
    goto pt_many_mul_end;

  l = mpz_sizeinbase (e, 2) - 1; /* l >= 1 */
  if(l < 32)
      status = pt_many_mul_add_sub_si(tQ, tP, tE, nE, mpz_get_si(e), n,
				      num, den, inv, ok);
  else
      status = pt_many_mul_plain(tQ, tP, tE, nE, e, n, num, den, inv, ok);


pt_many_mul_end:

  /* Undo negation to avoid changing the caller's e value */
  if (negated){
    mpz_neg (e, e);
    pt_many_neg(tP, nE, n);
  }
  return status;
}

/******************** Weierstrass section ********************/

void
pt_w_set_to_zero(ell_point_t P, mpmod_t n)
{
    mpres_set_ui (P->x, 0, n);
    mpres_set_ui (P->y, 1, n);
    mpres_set_ui (P->z, 0, n);
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
  mpres_set (x0, x, n);
  mpres_set (y0, y, n);
  mpres_set (z0, z, n);
}

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

int
pt_w_common_aff(mpres_t x0, mpres_t y0, mpres_t z0,
		mpres_t x1, mpres_t y1,
		mpres_t x2,
		mpmod_t n, mpres_t num, mpres_t den, mpres_t inv)
{
    if(mpres_invert(inv, den, n) == 0){
	mpres_gcd(x0, den, n);
	return 0;
    }
    mpres_mul(inv, inv, num, n);
    mpres_mul(num, inv, inv, n);
    mpres_sub(den, num, x1, n);
    mpres_sub(den, den, x2, n);
    mpres_sub(num, x1, den, n);
    mpres_mul(num, num, inv, n);
    mpres_sub(y0, num, y1, n);
    mpres_set(x0, den, n);
    mpz_set_ui(z0, 1); /* just in case */
    return 1;
}

/* [x0, y0, z0] <- [2] * [x, y, z] */
int
pt_w_duplicate(mpres_t x3, mpres_t y3, mpres_t z3,
	       mpres_t x1, mpres_t y1, mpres_t z1,
	       mpmod_t n, ell_curve_t E)
{
    if(pt_w_is_zero(z1, n))
      {
	pt_w_set(x3, y3, z3, x1, y1, z1, n);
	return 1;
      }
    if(E->type == ECM_EC_TYPE_WEIERSTRASS && E->law == ECM_LAW_AFFINE){
	/* buf[1] <- 2*y1 */
	mpres_add(E->buf[1], y1, y1, n);
	if(mpres_is_zero(y1, n)){
	    /* y1 = 0 <=> P is a [2]-torsion point */
	    pt_w_set(x3, y3, z3, x1, y1, z1, n);
	    return 1;
	}
	/* buf[0] <- 3*x^2+A */
	mpres_mul_ui(E->buf[0], x1, 3, n);
	mpres_mul(E->buf[0], E->buf[0], x1, n);
	mpres_add(E->buf[0], E->buf[0], E->A, n);
	return pt_w_common_aff(x3, y3, z3, x1, y1, x1, n, 
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
	mpres_mul(E->buf[1], E->buf[1], E->A, n);
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
pt_w_add(mpres_t x3, mpres_t y3, mpres_t z3,
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
	if(mpz_cmp(x1, x2) == 0 && mpz_cmp(y1, y2) == 0)
	    return pt_w_duplicate(x3, y3, z3, x1, y1, z1, n, E);
	else{
	    mpres_sub(E->buf[0], y1, y2, n);
	    mpres_sub(E->buf[1], x1, x2, n);
	    return pt_w_common_aff(x3, y3, z3, x1, y1, x2, n, 
				   E->buf[0], E->buf[1], E->buf[2]);
	}
    else if(E->type == ECM_EC_TYPE_WEIERSTRASS 
	    && E->law == ECM_LAW_HOMOGENEOUS){
	/* Cohen-Miyaji-Ono: 12M+2S+6add+1*2 */
	/* mapping: y1z2 = buf, AA = buf+1, u = buf+2, v = buf+3, R = buf+4, */
	/* vvv = buf+5; */
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
	    return pt_w_duplicate(x3, y3, z3, x1, y1, z1, n, E);
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

/* [x3, y3, z3] <- [x1, y1, z1] - [x2, y2, z2]; P3 != P1, P3 != P1. */
int
pt_w_sub(mpres_t x3, mpres_t y3, mpres_t z3,
	 mpres_t x1, mpres_t y1, mpres_t z1,
	 mpres_t x2, mpres_t y2, mpres_t z2,
	 mpmod_t n, ell_curve_t E)
{
    int res;

    mpres_neg(y2, y2, n);
    res = pt_w_add(x3, y3, z3, x1, y1, z1, x2, y2, z2, n, E);
    mpres_neg(y2, y2, n);
    return res;
}

/* Morain/Olivos */
int
MO_automaton(char *T, size_t *IT, mpz_t e, size_t le)
{
    size_t ie, iT = 0;
    int state = 0, res = 1, bz;

    for(ie = 0; ie < le; ie++){
	bz = mpz_tstbit(e, ie) == 0;
	switch(state){
	case 0:
	    if(bz)
		T[iT++] = 0;
	    else
		state = 1;
	    break;
	case 1:
	    if(bz){
		T[iT++] = 1;
		T[iT++] = 0;
		state = 0;
	    }
	    else{
		T[iT++] = -1;
                T[iT++] = 0;
                state = 11;
	    }
	    break;
	case 11:
	    if(bz)
		state = 110;
	    else{
		T[iT++] = 0;
		state = 11;
	    }
	    break;
	case 110:
	    if(bz){
		T[iT++] = 1;
                T[iT++] = 0;
                state = 0;
	    }
	    else{
		T[iT++] = -1;
                T[iT++] = 0;
                state = 11;
	    }
	}
    }
    if(state == 1 || state == 11)
	T[iT++] = 1;
    else if(state == 110)
	res = 0;
    *IT = iT;
    return res;
}

/* Do we have eval(T) == e? */
int
MO_check(char *T, size_t iT, mpz_t e)
{
    mpz_t tmp;
    int i, ok;

    mpz_init_set_ui(tmp, 0);
    for(i = ((int)iT)-1; i >= 0; i--){
	mpz_mul_2exp(tmp, tmp, 1);
	mpz_add_si(tmp, tmp, (int)T[i]);
    }
    gmp_printf("e:=%Zd;\n", e);
    gmp_printf("t:=%Zd;\n", tmp);
    ok = mpz_cmp(tmp, e) == 0;
    mpz_clear(tmp);
    return ok;
}

/* Do we have eval(T) == e? */
int
Split_check(short *S, int iS, mpz_t e)
{
    mpz_t tmp;
    int i, ok;

    mpz_init_set_ui(tmp, 0);
    for(i = 0; i < iS; i += 2){
	mpz_add_si(tmp, tmp, (int)S[i+1]);
	mpz_mul_2exp(tmp, tmp, (int)S[i]);
    }
    gmp_printf("e:=%Zd;\n", e);
    gmp_printf("t:=%Zd;\n", tmp);
    ok = mpz_cmp(tmp, e) == 0;
    mpz_clear(tmp);
    return ok;
}

/* Adapted from Koyama and Tsuroka, CRYPTO'92, using less space. 
   S is filled in left-to-right from T.
*/
int
Split(short *S, int Slen, char *T, size_t iT, int w)
{
    int i = (int)iT-1, iS = 0, gap, j, k, lW;
    short W;

    while(i >= w-1){
	/* next scan: T[i-w+1..i] */
	/* exclude right zeros */
	gap = 0;
	for(j = i-w+1; j <= i; j++){
	    if(T[j] != 0)
		break;
	    gap++;
	}
	lW = i-j+1;
	W = 0;
	/* at this point, T[j] <> 0 */
	for(k = j; k <= i; k++)
	    W = W+T[k]*(((short)1<<(k-j)));
	i = i-w;
	/* exclude left zeros and update power of 2 */
	while((i >= 0) && (T[i] == 0)){
	    i--;
	    gap++;
	}
	S[iS] = gap;
	S[iS+1] = W;
	if(iS >= 2)
	    S[iS-2] += lW;
	iS += 2;
	if(iS > Slen)
	    return -1;
    }
    /* at this point, we have to examine T[0..i] */
    if(i >= 0){
	/* exclude right zeros */
	gap = 0;
	for(j = 0; j <= i; j++){
	    if(T[j] != 0)
		break;
	    gap++;
	}
	lW = i-j+1;
	W = 0;
	/* at this point, T[j] <> 0 */
	for(k = j; k <= i; k++)
	    W = W+T[k]*(((short)1) << (k-j));
	S[iS] = gap;
	S[iS+1] = W;
	if(iS >= 2)
	    S[iS-2] += lW;
	iS += 2;
	if(iS > Slen)
	    return -1;
    }
    return iS;
}

/*
   OUTPUT: iS such that S[0..iS[ was filled
           -1 if Slen is too small
   The Solinas version is too slow for big entries, since it requires too
   many shifts.
   At the end of the process, we will have written
   e = 2^t0 * (2*d0+1 + 2^t1 *(2*d1+1 + 2^t2 * (2*d2+1+... + 2^ts*(2*ds+1) )
   where ti >= w and -2^(w-1)+1 <= 2*di+1 < 2^(w-1)+1.
   S will contain: [[ts, 2*ds+1], ..., [t1, 2*d1+1], [t0, 2*d0+1]].
*/
int
build_MO_chain(short *S, int Slen, mpz_t e, int w)
{
    /* first use automata */
    size_t le = mpz_sizeinbase(e, 2), iT = 0;
#if DEBUG_ADD_LAWS >= 2
    long tp = cputime();
    int i;
#endif
    char *T = (char *)malloc((2*le) * sizeof(char)); /* humf */
    int iS;

    MO_automaton(T, &iT, e, le);
#if DEBUG_ADD_LAWS >= 2
    /* check value of T */
    gmp_printf("# signed digits(%Zd):", e);
    for(i = 0; i < (int)iT; i++)
	printf(" %d", T[i]);
    printf("\n");
    if(MO_check(T, iT, e) == 0)
	printf("#!# Error in MO\n");
    else
	printf("# good check in MO\n");
    printf("# le = %ld, iT = %ld, time = %ldms\n",le,iT,elltime(tp,cputime()));
    tp = cputime();
#endif
    /* compact T to fill in S */
    iS = Split(S, Slen, T, iT, w);
#if DEBUG_ADD_LAWS >= 2
    printf("# time = %ldms\n", elltime(tp, cputime()));
    printf("S =");
    for(i = 0; i < iS; i++)
	printf(" %d", S[i]);
    printf("\n");
    if(Split_check(S, iS, e) == 0)
        printf("#!# Error in Split\n");
    else
        printf("# good check in Split\n");

#endif
    free(T);
    return iS;
}

int
build_add_sub_chain(short *S, int Slen, mpz_t e, int w)
{
    return build_MO_chain(S, Slen, e, w);
}

/* Checks that x1/z1 = x2/z2 and y1/z1 = y2/z2.
   OUTPUT: 1 if equals, 0 otherwise.
 */
int
pt_w_cmp(mpres_t x1, mpres_t y1, mpres_t z1,
	 mpres_t x2, mpres_t y2, mpres_t z2,
	 mpmod_t n)
{
    if(pt_w_is_zero(z1, n))
	return pt_w_is_zero(z2, n);
    else if(pt_w_is_zero(z2, n))
	return pt_w_is_zero(z1, n);
    else{
	mpres_t tmp1, tmp2;
	int cmp = 1;

	mpres_init(tmp1, n);
	mpres_init(tmp2, n);
	mpres_mul(tmp1, x1, z2, n);
	mpres_mul(tmp2, x2, z1, n);
	mpres_sub(tmp1, tmp1, tmp2, n);
	if(mpres_is_zero(tmp1, n) == 0){
	    printf("x1/z1 != x2/z2\n");
	    cmp = 0;
	}
	else{
	    mpres_mul(tmp1, y1, z2, n);
	    mpres_mul(tmp2, y2, z1, n);
	    mpres_sub(tmp1,tmp1, tmp2, n);
	    cmp = mpres_is_zero(tmp1, n);
	    if(cmp == 0)
		printf("y1/z1 != y2/z2\n");
	}
	mpres_clear(tmp1, n);
	mpres_clear(tmp2, n);
	return cmp;
    }
}

/******************** Hessian form ********************/

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

void
hessian_print(ell_point_t P, ell_curve_t E, mpmod_t n)
{
    pt_w_print(P->x, P->y, P->z, E, n);
}

/* -[u:v:w] = [v:u:w] */
void
hessian_negate(ell_point_t P, ATTRIBUTE_UNUSED ell_curve_t E, ATTRIBUTE_UNUSED mpmod_t n)
{
    mpz_swap(P->x, P->y); /* humf */
}

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

    if(mpz_cmp(E->buf[6], E->buf[2]) == 0 
       && mpz_cmp(E->buf[4], E->buf[1]) == 0)
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

int
hessian_sub(ell_point_t R, ell_point_t P, ell_point_t Q, ell_curve_t E, mpmod_t n)
{
    int ret;

    hessian_negate(Q, E, n);
    ret = hessian_add(R, P, Q, E, n);
    hessian_negate(Q, E, n);
    return ret;
}

/* switch from X^3+Y^3+1=3*D*X*Y to Y^2=X^3+A*X+B
   A:=-27*D*(D^3+8);
   B:=54*(D^6-20*D^3-8);
   xi:=12*(D^3-1)/(D*u+v+1);
   x:=-9*D^2+xi*u;
   y:=3*xi*(v-1);
   If a factor is found during the inversion, it is put in f and
   ECM_FACTOR_FOUND_STEP1 is returned. Otherwise, ECM_NO_FACTOR_FOUND is
   returned.
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
    printf("A:=");
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
    if(ell_point_mul(Q, e, P, E, n) == 0)
	mpz_set(f, Q->x);
    else{
	mpres_set(x, Q->x, n);
	mpres_set(y, Q->y, n);
    }
    mpz_clear(e);
    ell_point_clear(Q, E, n);
    ell_point_clear(P, E, n);
    ell_curve_clear(E, n);
    return ret;
}

/******************** generic ec's ********************/

void
ell_point_init(ell_point_t P, ell_curve_t E, mpmod_t n)
{
    mpres_init(P->x, n);
    mpres_init(P->y, n);
    mpres_init(P->z, n);
    if(E->type == ECM_EC_TYPE_WEIERSTRASS && E->law == ECM_LAW_AFFINE)
	mpz_set_ui (P->z, 1);
    else if(E->type == ECM_EC_TYPE_WEIERSTRASS 
	    && E->law == ECM_LAW_HOMOGENEOUS)
	mpres_set_ui(P->z, 1, n);
    else if(E->type == ECM_EC_TYPE_HESSIAN)
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

void
ell_point_print(ell_point_t P, ell_curve_t E, mpmod_t n)
{
    if(E->type == ECM_EC_TYPE_WEIERSTRASS)
	pt_w_print(P->x, P->y, P->z, E, n);
    else if(E->type == ECM_EC_TYPE_HESSIAN)
	hessian_print(P, E, n);
}

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
    mpres_init(E->A, n);
    for(i = 0; i < EC_W_NBUFS; i++)
	mpres_init (E->buf[i], n);
}

void
ell_curve_init_set(ell_curve_t E, int etype, int law, mpres_t A, mpmod_t n)
{
    ell_curve_init(E, etype, law, n);
    mpres_set(E->A, A, n);
}

void
ell_curve_set_z(ell_curve_t E, ell_curve_t zE, mpmod_t n)
{
    ell_curve_init(E, zE->type, zE->law, n);
    mpres_set_z(E->A, zE->A, n);
    E->disc = zE->disc;
    if(E->disc != 0){
	mpres_init(E->sq[0], n);
	mpres_set_z(E->sq[0], zE->sq[0], n);
    }
}

void
ell_curve_clear(ell_curve_t E, mpmod_t n)
{
    int i;

    mpres_clear(E->A, n);
    for(i = 0; i < EC_W_NBUFS; i++)
	mpres_clear (E->buf[i], n);
    /* TODO: case of sq */
}

void
ell_curve_print(ell_curve_t E, mpmod_t n)
{
    if(E->type == ECM_EC_TYPE_WEIERSTRASS){
	printf("A:="); print_mpz_from_mpres(E->A, n); printf(";\n");
	printf("E:=[A, y0^2-x0^3-A*x0];\n");
    }
    else if(E->type == ECM_EC_TYPE_HESSIAN){
	printf("D:="); print_mpz_from_mpres(E->A, n); printf(";\n");
	printf("E:=[D];\n");
    }
}

int
ell_point_is_zero(ell_point_t P, ell_curve_t E, mpmod_t n)
{
    if(E->type == ECM_EC_TYPE_WEIERSTRASS)
	return pt_w_is_zero(P->z, n);
    else if(E->type == ECM_EC_TYPE_HESSIAN)
	return hessian_is_zero(P, E, n);
    return 0;
}

void
ell_point_set_to_zero(ell_point_t P, ell_curve_t E, mpmod_t n)
{
    if(E->type == ECM_EC_TYPE_WEIERSTRASS)
	pt_w_set_to_zero(P, n);
    else if(E->type == ECM_EC_TYPE_HESSIAN)
	hessian_set_to_zero(P, E, n);
}

int
ell_point_add(ell_point_t R, ell_point_t P, ell_point_t Q, ell_curve_t E, mpmod_t n)
{
    if(E->type == ECM_EC_TYPE_WEIERSTRASS)
	return pt_w_add(R->x, R->y, R->z, P->x, P->y, P->z, Q->x, Q->y, Q->z,
			n, E);
    else if(E->type == ECM_EC_TYPE_HESSIAN)
	return hessian_add(R, P, Q, E, n);
    else
	return -1;
}

/* R <- P-Q */
int
ell_point_sub(ell_point_t R, ell_point_t P, ell_point_t Q, ell_curve_t E, mpmod_t n)
{
    if(E->type == ECM_EC_TYPE_WEIERSTRASS)
	return pt_w_sub(R->x, R->y, R->z, P->x, P->y, P->z, Q->x, Q->y, Q->z,
			n, E);
    else if(E->type == ECM_EC_TYPE_HESSIAN)
	return hessian_sub(R, P, Q, E, n);
    else
	return -1;
}

int
ell_point_duplicate(ell_point_t R, ell_point_t P, ell_curve_t E, mpmod_t n)
{
    if(E->type == ECM_EC_TYPE_WEIERSTRASS)
	return pt_w_duplicate(R->x, R->y, R->z, P->x, P->y, P->z, n, E);
    else if(E->type == ECM_EC_TYPE_HESSIAN)
	return hessian_duplicate(R, P, E, n);
    else
	return -1;
}

void
ell_point_negate(ell_point_t P, ell_curve_t E, mpmod_t n)
{
    if(ell_point_is_zero(P, E, n) != 0){
	if(E->type == ECM_EC_TYPE_WEIERSTRASS)
	    mpres_neg(P->y, P->y, n);
	else if(E->type == ECM_EC_TYPE_HESSIAN)
	    hessian_negate(P, E, n);
    }
}

/* Q <- [e]*P
   Return value: 0 if a factor is found, and the factor is in Q->x,
                 1 otherwise.
*/
int
ell_point_mul_plain (ell_point_t Q, mpz_t e, ell_point_t P, ell_curve_t E, mpmod_t n)
{
  size_t l;
  int negated = 0, status = 1;
  ell_point_t P0;

  if(ell_point_is_zero(P, E, n)){
      ell_point_set(Q, P, E, n);
      return 1;
  }

  if (mpz_sgn (e) == 0)
    {
      ell_point_set_to_zero(Q, E, n);
      return 1;
    }

  /* The negative of a point (x:y:z) is (x:-y:z) */
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
	if(ell_point_duplicate (P0, P0, E, n) == 0)
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
	      if(ell_point_add (P0, P0, P, E, n) == 0)
	      {
		status = 0;
		break;
	      }
#if DEBUG_ADD_LAWS >= 2
	      printf("Radd:="); ell_point_print(P0, E, n);printf(";\n");
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

#define EC_ADD_SUB_WMAX 10
#define EC_ADD_SUB_2_WMAX (1 << EC_ADD_SUB_WMAX)

/* TODO: do better */
int
get_add_sub_w(mpz_t e)
{
    size_t l = mpz_sizeinbase(e, 2);

    if(l <= 16)
	return 2;
    else if(l <= 32)
	return 3;
    else if(l <= 128)
	return 4;
    else if(l <= 1024)
	return 5;
    else if(l <= 10240)
	return 6;
    else if(l <= 102400)
	return 7;
    else if(l <= 1024000)
	return 8;
    else
	return 9;
}

/* pack everybody */
void
add_sub_pack(mpz_t s, int w, short *S, int iS)
{
    int nsh, cte = sizeof(mp_limb_t)/sizeof(short);
    short *tmp;

    nsh = iS / cte;
    if(iS % cte != 0)
	nsh++;
    nsh *= cte;
    nsh += 4;
    /* coding */
    tmp = (short *)malloc(nsh * sizeof(short));
    tmp[0] = w;
    tmp[1] = iS / (1 << 16);
    tmp[2] = iS % (1 << 16);
    memcpy(tmp+4, S, iS * sizeof(short));
    s->_mp_d = (mp_limb_t *)tmp; /* humf */
}

void add_sub_unpack(int *w, short **S, int *iS, mpz_t s)
{
    short *T;

    T = (short *)s->_mp_d; /* humf */
    *w = (int)T[0];
    *iS = (((int)T[1]) << 16) + (int)T[2];
    *S = T+4;
}

/* INPUT: S = [[ts, 2*ds+1], ..., [t1, 2*d1+1], [t0, 2*d0+1]] for
   e = 2^t0 * (2*d0+1 + 2^t1 *(2*d1+1 + 2^t2 * (2*d2+1+... + 2^ts*(2*ds+1) ),
   with -2^(w-1)+1 <= 2*di+1 < 2^{w-1}.
*/
int
ell_point_mul_add_sub_with_S(ell_point_t Q, ell_point_t P, ell_curve_t E,
			    mpmod_t n, int w, short *S, int iS)
{
    ell_point_t P0;
    ell_point_t iP[EC_ADD_SUB_2_WMAX];
    int status = 1, i, j, k;

    /* iP[i] <- (2*i+1) * P */
    k = (1 << (w-1)) - 1;
    for(i = 0; i <= k; i++)
	ell_point_init(iP[i], E, n);
    ell_point_set(iP[0], P, E, n);
    if(k > 0){
	/* P[k] <- [2]*P */
	if(ell_point_duplicate(iP[k], P, E, n) == 0){
	    mpres_set(P0->x, iP[k]->x, n);
            status = 0;
	    goto ell_point_mul_add_sub_end;
        }
	for(i = 1; i <= k; i++){
	    if(ell_point_add(iP[i], iP[i-1], iP[k], E, n) == 0){
		mpres_set(P0->x, iP[i]->x, n);
                status = 0;
		goto ell_point_mul_add_sub_end;
	    }
	}
	/* at this point, P[i] = (2*i+1) P */
    }
  
    ell_point_init(P0, E, n);
    ell_point_set_to_zero(P0, E, n);

#if DEBUG_ADD_LAWS >= 2
    printf("P:="); ell_point_print(P, E, n); printf(";\n");
#endif
    /* S = [[ts, 2*ds+1], ... */
    for(j = 0; j < iS; j += 2){
#if DEBUG_ADD_LAWS >= 2
	printf("P0:="); ell_point_print(P0, E, n); printf(";\n");
#endif
	i = abs(S[j+1]) >> 1; /* (abs(S[j+1])-1)/2, S[j+1] is always odd */
	assert(i <= k);
	if(S[j+1] > 0){
	    if(ell_point_add(P0, P0, iP[i], E, n) == 0){
		status = 0;
		break;
	    }
#if DEBUG_ADD_LAWS >= 2
	    printf("iP%d:=", i); ell_point_print(iP[i], E, n); printf(";\n");
	    printf("Radd:="); ell_point_print(P0, E, n); printf(";\n");
	    printf("Q:=ProjEcmAdd(P0, iP%d, E, N); ProjEcmEqual(Q, Radd, N);\n", i);
#endif
	}
	else{
	    /* add(-P) = sub(P) */
	    if(ell_point_sub(P0, P0, iP[i], E, n) == 0){
		status = 0;
		break;
	    }
#if DEBUG_ADD_LAWS >= 2
	    printf("Rsub:="); ell_point_print(P0, E, n); printf(";\n");
#endif
	}
	/* now multiply */
	for(i = 0; i < S[j]; i++){
	    if(ell_point_duplicate(P0, P0, E, n) == 0){
		status = 0;
		break;
	    }
#if DEBUG_ADD_LAWS >= 2
	    printf("Rdup:="); ell_point_print(P0, E, n); printf(";\n");
#endif
	}
	if(status == 0)
	    break;
    }
 ell_point_mul_add_sub_end:
    ell_point_set(Q, P0, E, n);
    ell_point_clear(P0, E, n);
    for(i = 0; i <= k; i++)
	ell_point_clear(iP[i], E, n);
    return status;
}

/* multiply P=(x:y:z) by e and puts the result in Q.
   Return value: 0 if a factor is found, and the factor is in Q->x,
                 1 otherwise.
   See Solinas 2000 for the most plug-and-play presentation.		 
*/
int
ell_point_mul_add_sub (ell_point_t Q, mpz_t e, ell_point_t P,
		      ell_curve_t E, mpmod_t n)
{
    int negated = 0, status = 1, iS = 0, w, Slen;
#if DEBUG_ADD_LAWS >= 2
    int j;
#endif
    short *S;
    
    if(ell_point_is_zero(P, E, n)){
	ell_point_set(Q, P, E, n);
	return 1;
    }
    
    if(mpz_sgn(e) == 0){
	ell_point_set_to_zero(Q, E, n);
	return 1;
    }
    
    if(mpz_sgn (e) < 0){
	negated = 1;
	mpz_neg(e, e);
	ell_point_negate(P, E, n);
    }
    
    if (mpz_cmp_ui(e, 1) == 0){
	ell_point_set(Q, P, E, n);
	return 1;
    }

    w = get_add_sub_w(e);

    Slen = 2 * mpz_sizeinbase(e, 2);
    /*    printf("# Slen=%d\n", Slen); */
    S = (short *)malloc(Slen * sizeof(short));
    iS = build_add_sub_chain(S, Slen, e, w);
    if(iS == -1){
	printf("build_NAF: Slen=%d too small\n", Slen);
	return -1;
    }
#if DEBUG_ADD_LAWS >= 2
    gmp_printf("addsub[%Zd=>%d]:", e, iS);
    for(j = iS-1; j >= 0; j--)
	printf(" %d", S[j]);
    printf("\n");
    printf("P:="); ell_point_print(P, E, n); printf(";\n");
#endif
    status = ell_point_mul_add_sub_with_S(Q, P, E, n, w, S, iS);
    free(S);
#if DEBUG_ADD_LAWS >= 2
    if(status == 0){
	printf("Not checking, since a factor was found!\n");
    }
    else{
	ell_point_t PP;
	mpz_t f;
	int res;

	mpz_init(f);
	ell_point_init(PP, E, n);
	res = ell_point_mul_plain(PP, e, P, E, n);
	if(res == 0){
	    printf("Factor found during ell_point_mul_plain...!\n");
	}
	else if(pt_w_cmp(Q->x, Q->y, Q->z, PP->x, PP->y, PP->z, n) != 1){
	    mpz_gcd(f, PP->z, n->orig_modulus);
	    if(mpz_cmp_ui(f, 1) != 0){
		gmp_printf("non trivial gcd from plain: %Zd\n", f);
		mpz_gcd(f, Q->z, n->orig_modulus);
		gmp_printf("gcd from addsub: %Zd\n", f);
	    }
	    else{
		printf("PB\n");
		gmp_printf("N:=%Zd;\n", n->orig_modulus);
		gmp_printf("e:=%Zd;\n", e);
		printf("P:="); ell_point_print(P, E, n); printf(";\n");
		printf("x0:=P[1]/P[3] mod N; y0:=P[2]/P[3] mod N;\n");
		ell_curve_print(E, n); printf("E:=E mod N;\n");
		printf("addsub:="); ell_point_print(Q, E, n); printf(";\n");
		printf("plain:="); ell_point_print(PP, E, n); printf(";\n");
		exit(-1);
	    }
	}
	ell_point_clear(PP, E, n);
	mpz_clear(f);
    }
#endif
    /* Undo negation to avoid changing the caller's e value */
    if (negated){
	ell_point_negate(P, E, n);
	mpz_neg(e, e);
    }
    return status;
}

int
ell_point_mul(ell_point_t Q, mpz_t e, ell_point_t P, ell_curve_t E, mpmod_t n)
{
#if 0 /* keeping it simple */
    return ell_point_mul_plain(Q, e, P, E, n);
#else
    return ell_point_mul_add_sub(Q, e, P, E, n);
#endif
}

int *
compute_forbidden_res(int disc)
{
    int *t = NULL;

    if(disc == 0)
	return NULL;
    if(disc == -3){
	/* we do not want p = 2 mod 3 */
	t = (int *)malloc(3 * sizeof(int));
	t[0] = 3;
	t[1] = 2;
	t[2] = -1;
    }
    else if(disc == -4){
	/* we do not want p = 3 mod 4 */
	t = (int *)malloc(3 * sizeof(int));
	t[0] = 4;
	t[1] = 3;
	t[2] = -1;
    }
    else if(disc == -8){
	/* (-2/p) = -1 <=> p = 5 or 7 mod 8 */
	t = (int *)malloc(4 * sizeof(int));
	t[0] = 8;
	t[1] = 5;
	t[2] = 7;
	t[3] = -1;
    }
    else if(disc == -7 || disc == -11){
	/* (-d/p) = (p/d) when d is 3 mod 4 */
	int x, i, d = -disc;

	/* initialize */
	t = (int *)malloc(d * sizeof(int));
	memset(t, 0, d * sizeof(int));
	/* crude, but sufficient */
	for(x = 0; x < d; x++)
	    t[(x*x)%d] = 1;
	/* x = 0 is always ok */
	t[0] = d;
	i = 1;
	for(x = 1; x < d; x++)
	    if(t[x] == 0)
		t[i++] = x;
	t[i++] = -1;
#if 0
	for(x = 0; x < i; x++)
	    printf(" %d", t[x]);
	printf("\n");
#endif
    }
    return t;
}

/* We can probably hack so that s contains the coding of a NAF, containing
   w, iS, S.
*/
int
compute_s_4_add_sub(mpz_t s, unsigned long B1, int disc)
{
    mpz_t t;
    long tp;
    short *S;
    int iS, Slen, w, *forbiddenres = compute_forbidden_res(disc);

    mpz_init(t);
    tp = cputime();
    compute_s(t, B1, forbiddenres);
    free(forbiddenres);
    printf("# computing prod(p^e <= %lu): %ldms\n", B1, elltime(tp,cputime()));
#if 0 /* keeping it simple for the time being */
    mpz_set(s, t);
#else
    tp = cputime();
    w = get_add_sub_w(t);
    /* Slen = 2 * log_{2^w}(t) = 2*log_2(t)/w = 2 * 64 * size(t)/w */
    Slen = (2 * GMP_NUMB_BITS * mpz_size(t)) / w;
    S = (short *)malloc(Slen * sizeof(short));
    iS = build_add_sub_chain(S, Slen, t, w);
    printf("# NAF has %d terms (w=%d, Slen=%d): %ldms\n", iS, w, Slen,
	   elltime(tp,cputime()));
    if(iS == -1){
	printf("build_NAF: Slen=%d too small\n", Slen);
	return 0;
    }
    add_sub_pack(s, w, S, iS);
    free(S);
#endif
    mpz_clear(t);
    return 1;
}

