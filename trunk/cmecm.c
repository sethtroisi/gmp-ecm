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

#include "torsions.h"
#include "cmecm.h"

static int
adjust_CM(mpz_t f, ec_curve_t E, ec_point_t P, mpz_t N, mpz_t j)
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
	mpz_set(E->A, k);
	mpz_add(j, k, k);
	mpz_mod(j, j, N);
	mpz_add(E->A, E->A, j);
	mpz_mod(E->A, E->A, N);
	x0 = 0;
	ec_force_point(E, P, j, &x0, N);
    }
    return ret;
}

/* Curves are built mod N, not in nres */
int
build_curves_with_CM(mpz_t f, int *nE, ec_curve_t *tE, ec_point_t *tP, 
		     int disc, mpmod_t n, mpz_t *sqroots)
{
    mpz_t j, tmp;
    long x0;
    int i, ret = ECM_NO_FACTOR_FOUND, imax;

    ec_curve_init(tE[0], n);
    ec_point_init(tP[0], tE[0], n);
    tE[0]->disc = disc;
    *nE = 1;
    if(disc == -3){
	/* D = -3 => E: Y^2 = X^3 + 8 has rank 1, generator (2 : 4 : 1)
	   f4 = P2*P4, disc(P2) = 3*2^4
	*/
	imax = (sqroots == NULL ? 1 : 6);

	for(i = 0; i < imax; i++){
	    if(i > 0){
		ec_curve_init(tE[i], n);
		ec_point_init(tP[i], tE[i], n);
	    }
	    tE[i]->type = ECM_EC_TYPE_WEIERSTRASS;
	    mpz_set_ui(tE[i]->A, 0);
	    tE[i]->disc = -3;
	}
	mpz_set_si(tP[0]->x, 2);
	mpz_set_si(tP[0]->y, 4);

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
	    gmp_printf("# zeta6^3+1=%Zd\n", tmp);
	    mpz_set_ui(tmp, 8);
	    for(i = 1; i < imax; i++){
		mpz_mul(tmp, tmp, zeta6);
		mpz_mod(tmp, tmp, n->orig_modulus);
		x0 = 0;
		/* works since tE[i]->A is always 0... */
		ec_force_point(tE[i], tP[i], tmp, &x0, n->orig_modulus);
	    }
	    *nE = 6;
	    mpz_clear(zeta6);
	    mpz_clear(tmp);
	    mpz_set(tE[0]->sq[0], sqroots[0]);
	}
    }
    else if(disc == -4){
#if 1
	/* Y^2 = X^3 + 9 * X has rank 1 and a 4-torsion point */
	/* a generator is (4 : 10 : 1) */
	imax = (sqroots == NULL ? 1 : 4);

	for(i = 0; i < imax; i++){
	    if(i > 0){
		ec_curve_init(tE[i], n);
		ec_point_init(tP[i], tE[i], n);
	    }
	    tE[i]->type = ECM_EC_TYPE_WEIERSTRASS;
	    tE[i]->disc = -4;
	}
        mpz_set_ui(tE[0]->A, 9);
        mpz_set_si(tP[0]->x, 4);
        mpz_set_si(tP[0]->y, 10);
	if(sqroots != NULL){
	    /* sqroots[0] = sqrt(-1) */
	    mpz_init(tmp);
	    mpz_mul(tmp, sqroots[0], sqroots[0]);
	    mpz_add_si(tmp, tmp, 1);
	    mpz_mod(tmp, tmp, n->orig_modulus);
	    gmp_printf("# zeta4^2+1=%Zd\n", tmp);
	    for(i = 1; i < imax; i++){
		mpz_mul(tmp, tE[i-1]->A, sqroots[0]);
		mpz_mod(tE[i]->A, tmp, n->orig_modulus);
		x0 = 1; /* x0 = 0 is bad, since this is a 2-torsion point */
		mpz_set_si(tmp, 0);
		ec_force_point(tE[i], tP[i], tmp, &x0, n->orig_modulus);
	    }
	    *nE = 4;
	    mpz_clear(tmp);
	    mpz_set(tE[0]->sq[0], sqroots[0]);
	}
#else /* one day, use this? */
	/* => 1/3*y^2 = x^3 + x, gen = (4/3, 10/3) */
	tE[0]->type = ECM_EC_TYPE_MONTGOMERY;
	mpres_set_ui(tE[0]->A, 0, n);
	mod_from_rat_str(f, "4/3", n->orig_modulus);
	mpres_set_z(tP[0]->x, f, n);
#endif
    }
    else if(disc == -7){
	/* E = y^2 = x^3 - 2222640*x - 1568294784
	           = x^3 - 5*7*(2^2*3^2*7)^2*x - 2*7^2*(2^2*3^2*7)^3
	   P = (2052 : 50112 : 1) */
	tE[0]->type = ECM_EC_TYPE_WEIERSTRASS;
	mpz_set_si(tE[0]->A, -2222640);
        mpz_set_si(tP[0]->x, 2052);
        mpz_set_si(tP[0]->y, 50112);
    }
    else if(disc == -8){
	/* D = -8: E_c: Y^2 = X^3+4*c*X^2+2*c^2*X => Montgomery when 2 = z^2 
	   c = 1 => rank = 1, generator is (-1 : -1 : 1)
	   alt.: [-2*3*5*cc^2, 2^3*7*cc^3], here with cc=12
	*/
	tE[0]->type = ECM_EC_TYPE_WEIERSTRASS;
	mpz_set_si(tE[0]->A, -4320);
        mpz_set_si(tP[0]->x, 12);
        mpz_set_si(tP[0]->y, -216);
    }
    else if(disc == -11){
	/*     E:=EllipticCurve([0, 0, 0, -2^5*3*11, 2^4*7*11^2]);
	       [ (33 : 121 : 1) ] */
	tE[0]->type = ECM_EC_TYPE_WEIERSTRASS;
	mpz_set_si(tE[0]->A, -1056);
	mpz_set_si(tP[0]->x, 33);
        mpz_set_si(tP[0]->y, 121);
    }
    /* class number 2 */
    else if(disc == -15){
	/* it must be that sqroots[0] contains sqrt(5) mod N */
	/* j = -(191025+85995*sqrt(5))/2 */
	tE[0]->type = ECM_EC_TYPE_WEIERSTRASS;
	mpz_init_set_si(j, 85995);
	mpz_mul(j, j, sqroots[0]);
	mpz_add_si(j, j, 191025);
	mpz_neg(j, j);
	mpz_mod(j, j, n->orig_modulus);
	mod_div_2(j, n->orig_modulus);
	ret = adjust_CM(f, tE[0], tP[0], n->orig_modulus, j);
	mpz_clear(j);
    }
    else{
	printf("Unknown discriminant: %d\n", disc);
	ret = ECM_ERROR;
    }
    return ret;
}

