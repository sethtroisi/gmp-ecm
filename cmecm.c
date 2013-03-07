/* cmecm.c - ECM with CM curves for specific numbers
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
#include "cmecm.h"

#define DEBUG_CMECM 0

static int
adjust_CM(mpz_t f, ell_curve_t E, ell_point_t P, mpz_t N, mpz_t j)
{
    mpz_t k;
    long x0;
    int ret = ECM_NO_FACTOR_FOUND;

    /* k = j/(1728-j) */
    mpz_init_set_si(k, 1728);
    mpz_sub(k, k, j);
    mpz_mod(k, k, N);
    if(mpz_invert(f, k, N) == 0){
	printf("# factor found (adjust_CM)\n");
	mpz_gcd(f, k, N);
	ret = ECM_FACTOR_FOUND_STEP1;
    }
    else{
	mpz_mul(k, f, j);
	mpz_mod(k, k, N);
	/* a = 3*k, b = 2*k */
	mpz_set(E->a4, k);
	mpz_add(j, k, k);
	mpz_mod(j, j, N);
	mpz_add(E->a4, E->a4, j);
	mpz_mod(E->a4, E->a4, N);
	mpz_set(E->a6, j);
	x0 = 0;
	ec_force_point(E, P, &x0, N);
    }
    return ret;
}

/* Curves are built mod N, not in nres */
int
build_curves_with_CM(mpz_t f, int *nE, ell_curve_t *tE, ell_point_t *tP, 
		     int disc, mpmod_t n, mpz_t *sqroots)
{
    mpz_t j, tmp;
    long x0;
    int i, ret = ECM_NO_FACTOR_FOUND, imax, found = 0;

    ell_curve_init(tE[0], ECM_EC_TYPE_WEIERSTRASS, ECM_LAW_AFFINE, n);
    ell_point_init(tP[0], tE[0], n);
    tE[0]->disc = disc;
    mpz_init_set_ui(tE[0]->sq[0], 1);
    *nE = 1;
    if(disc == -3){
	/* D = -3 => E: Y^2 = X^3 + 8 has rank 1, generator (2 : 4 : 1)
	   f4 = P2*P4, disc(P2) = 3*2^4
	*/
	found = 1;
	imax = (sqroots == NULL ? 1 : 6);

	for(i = 0; i < imax; i++){
	    if(i == 0){
		tE[i]->type = ECM_EC_TYPE_WEIERSTRASS;
		tE[i]->law = ECM_LAW_HOMOGENEOUS;
	    }
	    else{
		ell_curve_init(tE[i], ECM_EC_TYPE_WEIERSTRASS, ECM_LAW_HOMOGENEOUS, n);
		ell_point_init(tP[i], tE[i], n);
	    }
	    mpz_set_ui(tE[i]->a4, 0);
	    tE[i]->disc = -3;
	}
	mpz_set_si(tP[0]->x, 2);
	mpz_set_si(tP[0]->y, 4);
	mpz_set_si(tP[0]->z, 1);

	if(sqroots != NULL){
	    /* TODO: use complex twists? */
	    mpz_t zeta6;
	    /* compute zeta6 from sqrt(-3) */
	    mpz_init(zeta6);
	    mpz_init(tmp);
	    mpz_add_si(zeta6, sqroots[0], 1);
	    mod_div_2(zeta6, n->orig_modulus);
	    mpz_powm_ui(tmp, zeta6, 3, n->orig_modulus);
	    mpz_add_si(tmp, tmp, 1);
	    mpz_mod(tmp, tmp, n->orig_modulus);
	    if(mpz_sgn(tmp) != 0){
		gmp_printf("# zeta6^3+1=%Zd\n", tmp);
		exit(-1);
	    }
	    mpz_set_ui(tmp, 8);
	    for(i = 1; i < imax; i++){
		mpz_mul(tmp, tmp, zeta6);
		mpz_mod(tmp, tmp, n->orig_modulus);
		x0 = 0;
		/* works since tE[i]->a4 is always 0... */
		mpz_set(tE[i]->a6, tmp);
		ec_force_point(tE[i], tP[i], &x0, n->orig_modulus);
	    }
	    *nE = 6;
	    mpz_clear(zeta6);
	    mpz_clear(tmp);
            for(i = 0; i < imax; i++)
		mpz_init_set(tE[i]->sq[0], sqroots[0]);
	}
    }
    else if(disc == -4){
	found = 1;
	mpz_init(tmp);
	if(sqroots != NULL){
	    /* sqroots[0] = sqrt(-1) */
	    mpz_mul(tmp, sqroots[0], sqroots[0]);
	    mpz_add_si(tmp, tmp, 1);
	    mpz_mod(tmp, tmp, n->orig_modulus);
	    if(mpz_sgn(tmp) != 0){
		gmp_printf("# zeta4^2+1=%Zd\n", tmp);
		exit(-1);
	    }
	}
#if 0
	/* Y^2 = X^3 + 9 * X has rank 1 and a 4-torsion point */
	/* a generator is (4 : 10 : 1) */
	imax = (sqroots == NULL ? 1 : 4);

	for(i = 0; i < imax; i++){
	    if(i == 0){
		tE[i]->type = ECM_EC_TYPE_WEIERSTRASS;
		tE[i]->law = ECM_LAW_HOMOGENEOUS;
	    }
	    else{
		ell_curve_init(tE[i], ECM_EC_TYPE_WEIERSTRASS, ECM_LAW_HOMOGENEOUS, n);
		ell_point_init(tP[i], tE[i], n);
	    }
	    tE[i]->disc = -4;
	}
        mpz_set_ui(tE[0]->a4, 9);
        mpz_set_si(tP[0]->x, 4);
        mpz_set_si(tP[0]->y, 10);
        mpz_set_si(tP[0]->z, 1);
	if(sqroots != NULL){
	    /* sqroots[0] = sqrt(-1) */
	    for(i = 1; i < imax; i++){
		mpz_mul(tmp, tE[i-1]->a4, sqroots[0]);
		mpz_mod(tE[i]->a4, tmp, n->orig_modulus);
		x0 = 1; /* x0 = 0 is bad, since this is a 2-torsion point */
		mpz_set_si(tE[i]->a6, 0);
		ec_force_point(tE[i], tP[i], &x0, n->orig_modulus);
	    }
	    *nE = 4;
	    for(i = 0; i < imax; i++)
		mpz_init_set(tE[i]->sq[0], sqroots[0]);
	}
#else 
	/* Montgomery form curves, we do not need zeta4
	   Y^2 = X^3+k^2*X -> (1/k)*y^2 = x^3+x, x=X/k, y=Y/k
	   k=3, gen=(4, 10) => 1/3*y^2 = x^3 + x, gen = (4/3, 10/3)
	   {k, num, den} on E.W; divide by k to obtain E.M.
	   Remark this imposes 4 | #E, so not all classes are possible...!
	*/
	{
	long data4[][3] = {{3, 4, 1}, {7, 16, 9}, {10, 5, 1}, {11, 4900, 9},
#if 0
			   {14, 2, 1}, {15, 36, 1}, {17, 144, 1},
			   {19, 722500, 77841}, {23, 7056, 121},
			   {26, 13, 9}, {31, 1787598400, 32798529},
			   {35, 225, 4}, {39, 2025, 4},
#endif
			   {0, 0, 0}};
	for(i = 0; data4[i][0] != 0; i++){
	    ell_curve_init(tE[i], ECM_EC_TYPE_MONTGOMERY, ECM_LAW_HOMOGENEOUS, n);
	    ell_point_init(tP[i], tE[i], n);
	    tE[i]->disc = -4;
	    mpz_init_set_ui(tE[i]->sq[0], 1);
	    /* compute abscissa of generator in Montgomery form */
	    mpz_set_ui(tE[i]->a4, data4[i][1]);
	    mpz_set_ui(f, data4[i][2]);
	    /* do not forget to divide by k */
	    mpz_mul_ui(f, f, data4[i][0]);
	    if(mod_from_rat2(tP[i]->x, tE[i]->a4, f, n->orig_modulus) == 0){
		printf("# factor found during Montgomery preparation\n");
		mpz_set(f, tP[i]->x);
	    }
	    mpz_set_ui(tE[i]->a4, 0);
	}
	}
	printf("# using %d curves in Montgomery form for disc=-4\n", i);
	*nE = i;
#endif
	mpz_clear(tmp);
    }
    if(found == 0){
	/* {D, A, B, x0, y0} */
	/* D = -8: E_c: Y^2 = X^3+4*c*X^2+2*c^2*X => Montgomery when 2 = z^2 
	   c = 1 => rank = 1, generator is (-1 : -1 : 1)
	   alt.: [-2*3*5*cc^2, 2^3*7*cc^3], here with cc=12
	*/
	long h1_data[][5] = {{-7, -2222640, -1568294784, 2052, 50112},
			     {-8, -4320, 96768, 12, -216},
			     {-11, -1056, 13552, 33, 121},
			     {0, 0, 0, 0, 0}};
	int i1;

	for(i1 = 0; h1_data[i1][0] != 0; i1++){
	    if(h1_data[i1][0] == disc){
		tE[0]->type = ECM_EC_TYPE_WEIERSTRASS;
		tE[0]->law = ECM_LAW_HOMOGENEOUS;
		mpz_set_si(tE[0]->a4, h1_data[i1][1]);
		mpz_set_si(tP[0]->x, h1_data[i1][3]);
		mpz_set_si(tP[0]->y, h1_data[i1][4]);
		mpz_set_si(tP[0]->z, 1);
		found = 1;
		break;
	    }
	}
    }
    if(found == 0){
	/* class number 2 */
	/* {D, q1, den, a0, a1} s.t. D = q1*(-q2), j = (a0+a1*sqrt(q1))/den */
	long h2_data[][5] = {{-15, 5, 2, -191025, -85995},
			     {0, 0, 0, 0, 0}};
	int i2;

	if(sqroots == NULL){
	    printf("# missing sqroot for disc=%d\n", disc);
	    return ret;
	}
	for(i2 = 0; h2_data[i2][0] != 0; i2++){
	    if(h2_data[i2][0] == disc){
		/* it must be that sqroots[0] contains sqrt(q1) mod N */
		/* j = (a0+a1*sqrt(q1))/den */
		tE[0]->type = ECM_EC_TYPE_WEIERSTRASS;
		tE[0]->type = ECM_LAW_HOMOGENEOUS;
		mpz_set(tE[0]->sq[0], sqroots[0]);
		mpz_init_set_si(j, h2_data[i2][4]);
		mpz_mul(j, j, sqroots[0]);
		mpz_add_si(j, j, h2_data[i2][3]);
		mpz_mod(j, j, n->orig_modulus);
		if(h2_data[i2][2] == 2)
		    mod_div_2(j, n->orig_modulus);
		ret = adjust_CM(f, tE[0], tP[0], n->orig_modulus, j);
		mpz_clear(j);
		found = 1;
		break;
	    }
	}
    }
    if(found == 0){
	printf("Unknown discriminant: %d\n", disc);
	ret = ECM_ERROR;
    }
    return ret;
}

/******************** Stage 2 ********************/

/* General idea: we want to match [U]P = [V][omega]P, where 
   4*B1 < U^2+D*V^2 <= 4*B2.
   F could have all ([u]P)_x as roots;
   G could have all ([v][omega]P)_x.
   We do not need to test Gcd(F, F').
*/

/********** NAIVE version **********/

/* Find if some x_u = xt_v where [u]P=(x_u, ...), [v]Q=(xt_v, ...). 
   A priori, #v < #u, so we precompute all v.
   One looks for u^2+d*v^2 with umin <= u <= umax, vmin <= v <= vmax,
   where umax = sqrt(cof*B2-d*vmin^2), vmax = sqrt((cof*B2-umin^2)/d).
*/
static
int LoopCM(mpz_t f, ell_curve_t E, ell_point_t P, ell_point_t Q,
	   mpmod_t modulus, mpz_t B2, int cof, int d, mpres_t omega,
	   unsigned long umin, unsigned long du,
	   unsigned long vmin, unsigned long dv)
{
    int ret = ECM_NO_FACTOR_FOUND, compute_uP = 1;
    ell_point_t uP, duP, vQ, dvQ;
    mpz_ptr xtv, ztv; /* suboptimal? */
    mpz_t tmp, xu, zu, tmp2, omega5;
    long tp = cputime();
    unsigned long umax, vmax, u, v, iv, vsup, nv;

    mpz_init(tmp);
    mpz_init(tmp2);
    /* u <= sqrt(cof*B2-d*vmin^2) */
    mpz_set(tmp, B2);
    mpz_mul_si(tmp, tmp, cof);
    mpz_set_ui(tmp2, vmin);
    mpz_mul_ui(tmp2, tmp2, vmin);
    mpz_mul_ui(tmp2, tmp2, d);
    mpz_sub(tmp, tmp, tmp2);
    mpz_sqrt(tmp, tmp);
    gmp_printf("# umax=%Zd\n", tmp);
    if(mpz_fits_ulong_p(tmp) == 0)
	gmp_printf("#!# Gasp: umax=%Zd too large\n", tmp);
    umax = mpz_get_ui(tmp);
    /* v <= sqrt((cof*B2-umin^2)/d) */
    mpz_set(tmp, B2);
    mpz_mul_si(tmp, tmp, cof);
    mpz_set_ui(tmp2, umin);
    mpz_mul_ui(tmp2, tmp2, umin);
    mpz_sub(tmp, tmp, tmp2);
    mpz_fdiv_q_ui(tmp, tmp, d);
    mpz_sqrt(tmp, tmp);
    gmp_printf("# vmax=%Zd\n", tmp);
    if(mpz_fits_ulong_p(tmp) == 0)
	gmp_printf("#!# Gasp: vmax=%Zd too large\n", tmp);
    vmax = mpz_get_ui(tmp);

    /* precompute all v's */
    nv = 1 + ((vmax-vmin)/dv);
    xtv = (mpz_ptr)malloc(nv * sizeof(__mpz_struct));
    ztv = (mpz_ptr)malloc(nv * sizeof(__mpz_struct));
    ell_point_init(vQ, E, modulus);
    mpz_set_ui(tmp, vmin);
    ell_point_mul(vQ, tmp, Q, E, modulus);
    ell_point_init(dvQ, E, modulus);
    mpz_set_ui(tmp, dv);
    ell_point_mul(dvQ, tmp, Q, E, modulus);
    for(v = vmin, iv = 0; v <= vmax; v += dv, iv++){
	mpz_init(xtv+iv);
	mpres_get_z(xtv+iv, vQ->x, modulus);
	mpz_init(ztv+iv);
	mpres_get_z(ztv+iv, vQ->z, modulus);
	ell_point_add(vQ, vQ, dvQ, E, modulus);
    }
    ell_point_clear(vQ, E, modulus);
    ell_point_clear(dvQ, E, modulus);
    printf("# computing %lu (// %lu) additions: %ldms\n", iv, nv,
	   elltime(tp, cputime()));

    /* compute all uP */
    tp = cputime();
    compute_uP = !(d == 3 && du == 1 && dv == 1);
    ell_point_init(uP, E, modulus);
    mpz_set_ui(tmp, umin);
    ell_point_mul(uP, tmp, P, E, modulus);
    ell_point_init(duP, E, modulus);
    mpz_set_ui(tmp, du);
    ell_point_mul(duP, tmp, P, E, modulus);
    if(compute_uP == 0){
	/* omega5 = omega^5 = 1/omega */
	mpz_init(omega5);
	mpres_get_z(omega5, omega, modulus);
	mpz_powm_ui(omega5, omega5, 5, modulus->orig_modulus);
    }

    mpz_init(xu);
    mpz_init(zu);
    mpz_set_ui(f, 1);
    for(u = umin; u <= umax; u += du){
	if(compute_uP || u > vmax){
	    mpres_get_z(xu, uP->x, modulus);
	    mpres_get_z(zu, uP->z, modulus);
	}
	else{
	    /* [omega5]Q = [omega^6]P = P */
	    mpz_set(xu, xtv+(u-1));
	    mpz_mul(xu, xu, omega5);
	    mpz_mod(xu, xu, modulus->orig_modulus);
	    mpz_set(zu, ztv+(u-1));
	}
	/* d*v^2 <= B2-u^2 */
	mpz_set_ui(tmp, u);
	mpz_mul_ui(tmp, tmp, u);
	mpz_sub(tmp, B2, tmp);
	if(d != 1)
	    mpz_fdiv_q_ui(tmp, tmp, d);
	mpz_sqrt(tmp, tmp);
	if(mpz_fits_ulong_p(tmp) == 0)
	    gmp_printf("#!# Gasp: vmax=%Zd too large\n", tmp);
	vsup = mpz_get_ui(tmp);
	/*	printf("# u = %lu => vsup = %lu\n", u, vsup);*/
	for(v = vmin, iv = 0; v <= vsup; v += dv, iv++){
	    /* xu/zu = xtv/ztv <=> xu*ztv-xtv*zu = 0 */
	    mpz_mul(tmp, xu, ztv+iv);
	    mpz_mul(tmp2, zu, xtv+iv);
	    mpz_sub(tmp, tmp, tmp2);
	    if(mpz_sgn(tmp) == 0){
		printf("[%ld]P = [%ld]Q?\n", u, v);
		ell_curve_print(E, modulus);
		exit(-1);
	    }
	    mpz_mul(f, f, tmp);
	    mpz_mod(f, f, modulus->orig_modulus);
	}
	mpz_gcd(f, f, modulus->orig_modulus);
	if(mpz_cmp_ui(f, 1) != 0){
	    gmp_printf("# Youpi: u=%ld, f=%Zd\n", u, f);
	    ret = ECM_FACTOR_FOUND_STEP2;
	    break;
	}
	if(compute_uP || u > vmax)
	    ell_point_add(uP, uP, duP, E, modulus);
	if(compute_uP == 0 && u == (vmax+1)){
	    /* compute next point */
	    mpz_set_ui(tmp, vmax+1);
	    ell_point_mul(uP, tmp, P, E, modulus);
	}
    }
    mpz_clear(xu);
    mpz_clear(zu);
    mpz_clear(tmp);
    mpz_clear(tmp2);
    ell_point_clear(uP, E, modulus);
    ell_point_clear(duP, E, modulus);
    if(compute_uP == 0)
	mpz_clear(omega5);
    for(v = vmin, iv = 0; v <= vmax; v += dv, iv++){
	mpz_clear(xtv+iv);
	mpz_clear(ztv+iv);
    }
    free(xtv);
    free(ztv);
    printf("# using all [u]P: %ldms\n", elltime(tp, cputime()));
    return ret;
}

/* Testing CM stage2 in a very naive way, for the time being. */
int stage2_CM(mpz_t f, ell_curve_t E, ell_point_t P, mpmod_t modulus, 
	      unsigned long dF, mpz_t B2, ATTRIBUTE_UNUSED mpz_t B2min)
{
    int ret = ECM_NO_FACTOR_FOUND;
    ell_point_t Q;
    mpz_t tmp;
    unsigned long umin, du, vmin, dv;
    mpres_t omega;

    printf("PE:="); ell_point_print(P, E, modulus); printf(";\n");
    if (dF == 0)
	return ECM_NO_FACTOR_FOUND;

    mpz_init(tmp);
    if(E->disc == -3){
	/* case where pi = U+V*omega, U and V with the same parity */
	/* omega = (1+sqrt(-3))/2 */
	mpres_init(omega, modulus);
	mpres_set(omega, E->sq[0], modulus);
	mpres_add_ui(omega, omega, 1, modulus);
	mpres_div_2exp(omega, omega, 1, modulus);
	/* 1st case: u^2+3*v^2 <= B2, no restriction on u, v */
	/* we can easily compute [v][omega]P as [omega]([v]P) */
	umin = 1;
	du = 1;
	vmin = 1;
	dv = 1;

	/* Q = [omega](P) = [omega*P.x, P.y, 1] */
#if DEBUG_CMECM >= 0
	gmp_printf("N:=%Zd;\n", modulus->orig_modulus);
	printf("zeta3:=");print_mpz_from_mpres(E->sq[0],modulus);printf(";\n");
#endif
	ell_point_init(Q, E, modulus);
	mpres_mul(Q->x, P->x, omega, modulus);
	mpres_set(Q->y, P->y, modulus);
	mpres_set_ui(Q->z, 1, modulus);
	ret = LoopCM(f,E,P,Q,modulus,B2,1,3,omega,umin,du,vmin,dv);

	if(ret == ECM_NO_FACTOR_FOUND){
	    /* 2nd case: u^2+3*v^2 <= 4*B2, u and v odd */
	    /* we can easily compute [v][omega]P as [omega]([v]P) */
	}
	ell_point_clear(Q, E, modulus);
	mpres_clear(omega, modulus);
    }
    else if(E->disc == -4){
	/* case where pi = u+v*zeta4, u odd, v even */
	umin = 1;
	du = 2;
	vmin = 2;
	dv = 2;

	/* Q = [zeta4](P) = [-P.x, zeta4*P.y, 1] */
	ell_point_init(Q, E, modulus);
	mpres_neg(Q->x, P->x, modulus);
	mpres_mul(Q->y, P->y, E->sq[0], modulus);
	mpres_set_ui(Q->z, 1, modulus);
	ret = LoopCM(f,E,P,Q,modulus,B2,1,1,E->sq[0],umin,du,vmin,dv);
	ell_point_clear(Q, E, modulus);
    }
    else{
	printf("# stage2_CM not ready for disc=%d\n", E->disc);
    }
    mpz_clear(tmp);
    return ret;
}

/********** fast version **********/

/* u^2+d*v^2 <= cof*B2 => number of v's is <= number of u's. 
   As a matter of fact, dF should be a power of 2.
*/
void
set_stage2_params_CM(unsigned long *pdF, unsigned long *pk, mpz_t B2, int disc)
{
    mpz_t tmp;
    unsigned long umax, vmax, dF = 0, k = 0;

    mpz_init(tmp);
    if(disc == -3){
	/* U^2+3*V^2 <= 4*B2 */
	/* tmp <- floor(sqrt(4*B2-3)) */
	mpz_mul_2exp(tmp, B2, 2);
	mpz_sub_si(tmp, tmp, 3);
	mpz_sqrt(tmp, tmp);
	umax = mpz_get_ui(tmp);
	/* find smallest power of 2 >= umax */
	for(dF = 1; dF < umax; dF <<= 1);
	k = 1;
    }
    else if(disc == -4){
	/* tmp <- floor(sqrt(B2-1)) */
	mpz_set(tmp, B2);
	mpz_sub_si(tmp, tmp, 1);
	mpz_sqrt(tmp, tmp);
#if CMECM_FAST == 0 /* standard version */
	/* u^2 + v^2 <= B2, u odd, v even */
	vmax = mpz_get_ui(tmp);
	/* largest v must be even */
	if(vmax % 2 == 1)
	    vmax--;
	/* 2 <= v <= vmax => 0 <= (v-2)/2 <= (vmax-2)/2=(vmax/2)-1 */
	vmax >>= 1;
	/* find smallest power of 2 >= vmax */
	for(dF = 1; dF < vmax; dF <<= 1);

	mpz_init_set(tmp, B2);
	mpz_sub_si(tmp, tmp, 4);
	mpz_sqrt(tmp, tmp);
	umax = mpz_get_ui(tmp);
	if(umax % 2 == 0)
	    umax--;
	for(k = 1; k*dF < umax; k++);
#else /* faster */
	/* u^2 + v^2 <= B2, u > v > 0 */
	umax = mpz_get_ui(tmp);
	/* find smallest power of 2 >= umax */
	for(dF = 1; dF < umax; dF <<= 1);
	k = 1;
#endif
    }
    else{
	int D8 = abs(disc) % 8, Dp;

	if(D8 == 0){
	    Dp = abs(disc)/4;
	    /* u^2+D'*v^2 <= B2, D' even, u odd, v any */
	    /* vmax = floor(sqrt((B2-1)/Dp)) */
	    mpz_set(tmp, B2);
	    mpz_sub_si(tmp, tmp, 1);
	    mpz_fdiv_q_ui(tmp, tmp, Dp);
	    mpz_sqrt(tmp, tmp);
	    vmax = mpz_get_ui(tmp);
	    for(dF = 1; dF < vmax; dF <<= 1);
	    /* umax = floor(sqrt(B2-Dp)) */
	    mpz_set(tmp, B2);
            mpz_sub_si(tmp, tmp, Dp);
	    mpz_sqrt(tmp, tmp);
            umax = mpz_get_ui(tmp);
	    if(umax % 2 == 0)
		umax--;
	    for(k = 1; k*dF < umax; k++);
	}
    }
    mpz_clear(tmp);
    *pdF = dF;
    *pk = k;
}

/* F[i] <- ([kmin+i*dk]P)_x for 0 <= i < dF. */
int all_multiples(mpz_t f, listz_t F, unsigned long dF, 
		  ell_curve_t E, ell_point_t P, mpmod_t modulus,
		  unsigned long kmin, unsigned long dk)
{
    int ret = ECM_NO_FACTOR_FOUND;
    ell_point_t kP, dkP;
    unsigned long k, ik;
    mpz_t tmp;

    mpz_init(tmp);
    ell_point_init(kP, E, modulus);
    mpz_set_ui(tmp, kmin);
    if(ell_point_mul(kP, tmp, P, E, modulus) == 0){
	printf("# factor found at strange place 1\n");
	mpz_set(f, kP->x);
	ret = ECM_FACTOR_FOUND_STEP2;
	goto end_of_all_multiples;
    }
    ell_point_init(dkP, E, modulus);
    mpz_set_ui(tmp, dk);
    if(ell_point_mul(dkP, tmp, P, E, modulus) == 0){
	printf("# factor found at strange place 2\n");
	mpz_set(f, dkP->x);
	ret = ECM_FACTOR_FOUND_STEP2;
	goto end_of_all_multiples;
    }
    for(k = kmin, ik = 0; ik < dF; k += dk, ik++){
	/* F[ik] <= kP_x */
#if 0
	if(k == 4570){
	    printf("[4570]*[i]*P="); ell_point_print(kP, E, modulus); printf("\n");
	}
#endif
	mpres_get_z(F[ik], kP->x, modulus);
	if(ell_point_add(kP, kP, dkP, E, modulus) == 0){
	    printf("# factor found when adding dkP for k=%ld\n", k);
	    mpz_set(f, kP->x);
	    ret = ECM_FACTOR_FOUND_STEP2;
	    break;
	}
    }
    ell_point_clear(dkP, E, modulus);
 end_of_all_multiples:
    ell_point_clear(kP, E, modulus);
    mpz_clear(tmp);
    return ret;
}

/* coeffs = {deg, c_deg, ..., c_0} */
void
mpres_eval_poly(mpres_t y, long* coeffs, mpres_t x, mpmod_t modulus)
{
    int deg, i;

    deg = coeffs[0];
    /* Horner's rule */
    mpres_set_si(y, coeffs[1], modulus);
    for(i = 2; i <= deg+1; i++){
	mpres_mul(y, y, x, modulus);
	if(coeffs[i] > 0)
	    mpres_add_ui(y, y, coeffs[i], modulus);
	else
	    mpres_sub_ui(y, y, -coeffs[i], modulus);
    }
}

/* coeffs = {deg, c_deg, ..., c_0} */
void
mpres_eval_diff(mpres_t y, long* coeffs, mpres_t x, mpmod_t modulus)
{
    long *tmp, i, deg;

    deg = coeffs[0];
    tmp = (long *)malloc((deg+1) * sizeof(long));
    tmp[0] = deg-1;
    for(i = 1; i <= deg; i++)
	tmp[i] = coeffs[i] * (deg-i+1);
    mpres_eval_poly(y, tmp, x, modulus);
    free(tmp);
}

int
eval_CM_mult(ell_point_t omegaP, ell_point_t P, long *num, long *den,
	     mpz_t sq[], mpmod_t modulus)
{
    int ret = 1;
    mpres_t nn, dd, inv, tmp;

    mpres_init(nn, modulus);
    mpres_eval_poly(nn, num, P->x, modulus);
    mpres_init(dd, modulus);
    mpres_eval_poly(dd, den, P->x, modulus);
#if DEBUG_CMECM >= 1
    printf("x:="); print_mpz_from_mpres(P->x, modulus); printf(";\n");
    printf("nn:="); print_mpz_from_mpres(nn, modulus); printf(";\n");
    printf("dd:="); print_mpz_from_mpres(dd, modulus); printf(";\n");
#endif
    mpres_init(inv, modulus);
    mpres_init(tmp, modulus);
    if(mpres_invert(inv, dd, modulus) == 0){
	mpres_gcd(omegaP->x, dd, modulus);
	ret = 0;
    }
    else{
	/* omX <- num/den */
	mpres_mul(omegaP->x, nn, inv, modulus);
	/* omY <- dnum/den-num*dden/den^2 */
	mpres_eval_diff(omegaP->y, num, P->x, modulus);
	mpres_eval_diff(tmp, den, P->x, modulus);
	mpres_mul(tmp, tmp, nn, modulus);
	mpres_mul(tmp, tmp, inv, modulus);
	mpres_sub(omegaP->y, omegaP->y, tmp, modulus);
	mpres_mul(omegaP->y, omegaP->y, inv, modulus);
	mpres_invert(inv, sq[0], modulus);
	mpres_mul(omegaP->y, omegaP->y, inv, modulus);
	mpres_mul(omegaP->y, omegaP->y, P->y, modulus);
#if DEBUG_CMECM >= 1
	printf("omX:="); 
	print_mpz_from_mpres(omegaP->x, modulus);
	printf(";\n");
	printf("y:="); print_mpz_from_mpres(P->y, modulus); printf(";\n");
	printf("omY:="); 
	print_mpz_from_mpres(omegaP->y, modulus);
	printf(";\n");
#endif
    }
    mpres_clear(nn, modulus);
    mpres_clear(dd, modulus);
    mpres_clear(inv, modulus);
    mpres_clear(tmp, modulus);
    return ret;
}

int
apply_CM(ell_point_t omegaP, int disc, mpz_t sq[], ell_point_t P, mpmod_t modulus)
{
    int ret = 1;

    if(disc == -3){
	/* Q = [omega](P) = [omega*P.x, P.y, 1] */
    }
    else if(disc == -4){
	/* Q = [zeta4](P) = [-P.x, zeta4*P.y, 1] */
#if DEBUG_CMECM >= 1
	gmp_printf("N:=%Zd;\n", modulus->orig_modulus);
	printf("zeta4:=");print_mpz_from_mpres(sq[0],modulus);printf(";\n");
#endif
	mpres_neg(omegaP->x, P->x, modulus);
	mpres_mul(omegaP->y, P->y, sq[0], modulus);
	mpres_set_ui(omegaP->z, 1, modulus);
    }
    else if(disc == -8){
	long num[] = {2, -1, 4*12, -18*12*12};
	long den[] = {1, 2, -8*12};

	mpres_t tmp;
	mpres_init(tmp, modulus);
	mpres_mul(tmp, sq[0], sq[0], modulus);
	mpres_add_ui(tmp, tmp, 2, modulus);
	printf("sq^2+2="); print_mpz_from_mpres(tmp, modulus); printf(";\n");
	ret = eval_CM_mult(omegaP, P, num, den, sq, modulus);
    }
    return ret;
}

int ecm_rootsF_CM(mpz_t f, listz_t F, unsigned long dF, curve *C, 
		  mpmod_t modulus)
{
    int ret = ECM_NO_FACTOR_FOUND;
    long tp;
    ell_curve_t E;
    ell_point_t P, omegaP;
    unsigned long vmin = 1, dv = 1;

#if 0
    printf("# Entering ecm_rootsF_CM with disc=%d dF=%ld\n", C->disc, dF);
#endif
    ell_curve_init(E, ECM_EC_TYPE_WEIERSTRASS, ECM_LAW_AFFINE, modulus);
    mpres_set(E->a4, C->A, modulus);
    ell_point_init(P, E, modulus);
    mpres_set(P->x, C->x, modulus);
    mpres_set(P->y, C->y, modulus);
    ell_point_init(omegaP, E, modulus);
    tp = cputime();
    if(C->disc == -3 || C->disc == -4){
#if CMECM_FAST == 0
	printf("#!# not working anymore\n");
#else /* we don't need another hero! */
	ret = all_multiples(f, F, dF, E, P, modulus, 1, 1);
#endif
    }
    else{
	apply_CM(omegaP, C->disc, C->sq, P, modulus);
	if(abs(C->disc) % 8 == 0){
	    vmin = 1;
	    dv = 1;
	}
	ret = all_multiples(f, F, dF, E, omegaP, modulus, vmin, dv);
    }
    printf("# computing %ld ec_adds: %ldms\n", dF, elltime(tp, cputime()));
    ell_point_clear(omegaP, E, modulus);
    ell_point_clear(P, E, modulus);
    ell_curve_clear(E, modulus);
    return ret;
}

ecm_roots_state_t *
ecm_rootsG_init_CM (mpz_t f, curve *X, root_params_t *root_params, 
		    ATTRIBUTE_UNUSED unsigned long dF, 
		    ATTRIBUTE_UNUSED unsigned long blocks, mpmod_t modulus)
{
    ecm_roots_state_t *state;
    progression_params_t *params; /* for less typing */
    ell_curve_t E;
    ell_point_t P, duP;
    unsigned long umin = 1, du = 1, k;
    mpz_t tmp;
    
    state = (ecm_roots_state_t *) malloc (sizeof (ecm_roots_state_t));
    if(state == NULL){
	mpz_set_si (f, -1);
	return NULL;
    }
    /* really needed for CM? */
    params = &(state->params);
    /* If S < 0, use degree |S| Dickson poly, otherwise use x^S */
    params->dickson_a = (root_params->S < 0) ? -1 : 0;
    params->S = abs (root_params->S);
    state->X = X;

    state->T = NULL;
    state->size_T = 0;
    state->params.size_fd = 2;
    state->fd = (point *) malloc (state->params.size_fd * sizeof (point));
    if(state->fd == NULL){
	free (state);
	mpz_set_si (f, -1);
	return NULL;
    }
    for(k = 0; k < state->params.size_fd; k++){
	mpres_init (state->fd[k].x, modulus);
	mpres_init (state->fd[k].y, modulus);
    }
    
    if(X->disc == -3 || X->disc == -4){
#if CMECM_FAST == 0
	umin = 1;
	du = 2;
#else
    /*    printf("# Skipping stuff in rootsG_init_CM for disc=%d\n", X->disc);*/
    return state;
#endif
    }
    else if(abs(X->disc) % 8 == 0){
	umin = 1;
	du = 2;
    }
    
    /* conversions */
    ell_curve_init(E, ECM_EC_TYPE_WEIERSTRASS, ECM_LAW_AFFINE, modulus);
    mpres_set(E->a4, X->A, modulus);
    ell_point_init(P, E, modulus);
    ell_point_init(duP, E, modulus);

    mpres_set(P->x, X->x, modulus);
    mpres_set(P->y, X->y, modulus);
    mpres_set_ui(P->z, 1, modulus);

    /* fd[0] <- [umin]*P */
    mpz_init_set_ui(tmp, umin);
    ell_point_mul(duP, tmp, P, E, modulus);
    mpres_set(state->fd[0].x, duP->x, modulus);
    mpres_set(state->fd[0].y, duP->y, modulus);

    /* fd[1] <- [du]*P */
    mpz_set_ui(tmp, du);
    ell_point_mul(duP, tmp, P, E, modulus);
    mpres_set(state->fd[1].x, duP->x, modulus);
    mpres_set(state->fd[1].y, duP->y, modulus);

    ell_point_clear(P, E, modulus);
    ell_point_clear(duP, E, modulus);
    ell_curve_clear(E, modulus);
    mpz_clear(tmp);
    return state;
}

/* Compute "next" polynomial G. For block i, this means computing
   [umin+du*(i*dF)..umin+du*((i+1)*dF)[ * P on E, or adding duP.
   We need to store [umin+du*(i*dF)]*P and duP in 'state'.
*/
int 
ecm_rootsG_CM (mpz_t f, listz_t G, unsigned long dF, ecm_roots_state_t *state, 
	       mpmod_t modulus)
{
    ell_curve_t E;
    ell_point_t uP, duP;
    unsigned long i;
    int ret = ECM_NO_FACTOR_FOUND;

    /* conversions */
    ell_curve_init(E, ECM_EC_TYPE_WEIERSTRASS, ECM_LAW_AFFINE, modulus);
    mpres_set(E->a4, state->X->A, modulus);
    ell_point_init(uP, E, modulus);
    ell_point_init(duP, E, modulus);

    mpres_set(uP->x, state->fd[0].x, modulus);
    mpres_set(uP->y, state->fd[0].y, modulus);
    mpres_set_ui(uP->z, 1, modulus);

    mpres_set(duP->x, state->fd[1].x, modulus);
    mpres_set(duP->y, state->fd[1].y, modulus);
    mpres_set_ui(duP->z, 1, modulus);

    for(i = 0; i < dF; i++){
	mpres_get_z(G[i], uP->x, modulus);
	if(ell_point_add(uP, uP, duP, E, modulus) == 0){
	    printf("# factor found for i=%ld\n", i);
	    mpz_set(f, uP->x);
	    ret = ECM_FACTOR_FOUND_STEP2;
	    break;
	}
    }

    /* put back last point */
    mpres_set(state->fd[0].x, uP->x, modulus);
    mpres_set(state->fd[0].y, uP->y, modulus);

    ell_point_clear(uP, E, modulus);
    ell_point_clear(duP, E, modulus);
    ell_curve_clear(E, modulus);
    return ret;
}

/* 
   Given F(X) = prod(X-x_u), compute G(X) = prod(X-(omega*P)_x), which
   is G(X) = F([omega]*X) made monic.
 */
int
compute_G_from_F(listz_t G, listz_t F, unsigned long dF, curve *X,
		 mpmod_t modulus)
{
    int ret = ECM_NO_FACTOR_FOUND;
    unsigned long j, ex;
    mpz_t tmp, tmp2;

    if(X->disc == -3){
	/* [omega](X, Y) = (omega*X, Y), hence G(X) = F(omega*X)/omega^dF */
	/*printf("# making G(X) = F(omega*X)\n");*/
	mpz_init(tmp);
	mpz_init(tmp2);
	/* get back sqrt(-3) */
	mpres_get_z(tmp, X->sq[0], modulus);
	/* omega = (-1+sqrt(-3))/2 */
	mpz_sub_si(tmp, tmp, 1);
	mod_div_2(tmp, modulus->orig_modulus);
#if 0
	{
	    gmp_printf("omega=%Zd\n", tmp);
	    mpz_powm_ui(tmp2, tmp, 3, modulus->orig_modulus);
	    gmp_printf("omega^3=%Zd\n", tmp2);
	}
#endif
	/* tmp2 <- omega^2 */
	mpz_powm_ui(tmp2, tmp, 2, modulus->orig_modulus);
	/* F(omega*X) = omega^dF*X^dF+...+omega*c1*X+c0 */
	/* => G(X) = X^dF + ... + omega^(1-dF)*c1*X+c0*omega^(-dF) */
	ex = dF % 3;
	if(ex) ex = 3-ex; /* ex = -dF mod 3 */
	/* [X^j]G = omega^(-dFmod3+j)[X^j]F */
	for(j = 0; j < dF; j++){
	    if(ex == 0)
		mpz_set(G[j], F[j]);
	    else if(ex == 1){
		mpz_mul(G[j], F[j], tmp);
		mpz_mod(G[j], G[j], modulus->orig_modulus);
	    }
	    else{ /* ex == 2 */
		mpz_mul(G[j], F[j], tmp2);
		mpz_mod(G[j], G[j], modulus->orig_modulus);
	    }
	    ex++;
	    if(ex == 3) ex = 0;
	}
	mpz_clear(tmp);
	mpz_clear(tmp2);
    }    
    else if(X->disc == -4){
	/* [omega](X, Y) = (-X, zeta4*Y), hence G(X) = F(-X) */
	assert(dF % 2 == 0); /* lazy? */
	/*printf("# making G(X) = F(-X)\n");*/
	for (j = 0; j < dF; j ++){
	    mpz_set(G[j], F[j]);
	    if(j & 1){
		/* [X^j]G = -[X^j]F */
		if(mpz_sgn(G[j]) != 0)
		    mpz_sub(G[j], modulus->orig_modulus, G[j]);
	    }
	}
    }
    return ret;
}
