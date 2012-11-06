/* manyecm.c - ECM with many curves in parallel */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gmp.h> /* GMP header file */
#include "ecm.h" /* ecm header file */

#include "ecm-impl.h"
#include "ecm-ecm.h"
#include "mpmod.h"

#define DEBUG_MANY_EC 0
/* #define ONE_CURVE_AT_A_TIME*/ /* PROTECTION */

#define NCURVE_MAX 1000

/********** group law on points **********/

#define pt_is_equal(EP, EQ) (mpz_cmp((EP)->x, (EQ)->x) == 0 \
	                     && mpz_cmp((EP)->y, (EQ)->y) == 0 \
			     && mpz_cmp((EP)->z, (EQ)->z) == 0)

int
pt_is_zero(curve *EP, mpmod_t n)
{
    return mpz_sgn(EP->z) == 0;
}

void
pt_set_to_zero(curve *EP, mpmod_t n)
{
    mpz_set_ui(EP->x, 0);
    mpres_set_ui(EP->y, 1, n);
    mpz_set_ui(EP->z, 0);
}

void
pt_assign(curve *EQ, curve *EP, mpmod_t n)
{
    mpres_set(EQ->x, EP->x, n);
    mpres_set(EQ->y, EP->y, n);
    mpres_set(EQ->z, EP->z, n);
}

void
pt_neg(curve *EP, mpmod_t n)
{
    if(pt_is_zero(EP, n) == 0)
	mpres_neg(EP->y, EP->y, n);
}

void
pt_many_set_to_zero(curve *tEP, int nEP, mpmod_t n)
{
    int i;

    for(i = 0; i < nEP; i++)
	pt_set_to_zero(tEP+i, n);
}

void
pt_many_neg(curve *tEP, int nEP, mpmod_t n)
{
    int i;

    for(i = 0; i < nEP; i++)
	pt_neg(tEP+i, n);
}

void
pt_many_assign(curve *tEQ, curve *tEP, int nEP, mpmod_t n)
{
    int i;

    for(i = 0; i < nEP; i++)
	pt_assign(tEQ+i, tEP+i, n);
}

static void
print_mpz_from_mpres(mpres_t x, mpmod_t n)
{
    mpz_t tmp;

    mpz_init(tmp);
    mpres_get_z(tmp, x, n);
    gmp_printf("%Zd", tmp);
    mpz_clear(tmp);
}

void
pt_print(curve EP, mpmod_t n)
{
    printf("[");
    print_mpz_from_mpres(EP.x, n);
    printf(", ");
    print_mpz_from_mpres(EP.y, n);
    printf(", ");
    print_mpz_from_mpres(EP.z, n);
    printf("]");
}

void
pt_many_print(curve *tEP, int nEP, mpmod_t n)
{
    int i;

    for(i = 0; i < nEP; i++){
	printf("%d: ", i);
	pt_print(tEP[i], n);
	printf(" on E.A=");
	print_mpz_from_mpres(tEP[i].A, n);
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
#if DEBUG_MANY_EC >= 1
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
#if DEBUG_MANY_EC >= 1
	printf("Factor[%d]: ", i);
	mpz_out_str (stdout, 10, inv[0]);
	printf ("\n");
#endif
	/* identifying the x[i]'s */
	for(i = 0; i < nx; i++){
	    mpres_gcd(x[nx], x[i], n);
	    if(mpz_cmp_ui(x[nx], 1) != 0){
#if DEBUG_MANY_EC >= 0
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
#if DEBUG_MANY_EC >= 1
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

int
pt_many_common(curve *tER, curve *tEP, curve *tEQ, int nEP, mpmod_t n, 
	       mpres_t *num, mpres_t *den, mpres_t *inv, char *takeit)
{
    int i;

    if(compute_all_inverses(inv, den, nEP, n, takeit) == 0){
	mpz_set(num[nEP], inv[0]);
	return 0;
    }
    for(i = 0; i < nEP; i++){
	if(takeit[i] != 1)
	    continue;
	/* l:=(inv[i]*num[i]) mod N; */
	mpres_mul(num[i], num[i], inv[i], n);
	/* x:=(l^2-P[1]-Q[1]) mod N; */
	mpres_sqr(den[i], num[i], n);
	mpres_sub(den[i], den[i], tEP[i].x, n);
	mpres_sub(den[i], den[i], tEQ[i].x, n);
	/* tR[i]:=[x, (l*(P[1]-x)-P[2]) mod N, 1]; */
	mpres_sub(tER[i].x, tEP[i].x, den[i], n);
	mpres_mul(tER[i].x, tER[i].x, num[i], n);
	mpres_sub(tER[i].y, tER[i].x, tEP[i].y, n);
	mpres_set(tER[i].x, den[i], n);
    }
    return 1;
}

int
pt_many_duplicate(curve *tEQ, curve *tEP, int nEP, mpmod_t n, 
		  mpres_t *num, mpres_t *den, mpres_t *inv, char *ok)
{
    char *takeit = (char *)malloc(nEP * sizeof(char));
    int i, res;

    memcpy(takeit, ok, nEP);
    for(i = 0; i < nEP; i++){
	if(ok[i] == 0)
	    continue; /* takeit[i] = 0 */
	if(pt_is_zero(tEP+i, n)){
	    takeit[i] = 0;
	    pt_set_to_zero(tEQ+i, n);
	}
	else if(mpz_sgn(tEP[i].y) == 0){
	    /* 2 * P[i] = O_E */
	    takeit[i] = 0;
	    pt_set_to_zero(tEP+i, n);
	    printf("# [2] * P[%d] = O_E\n", i);
	}
	else{
	    mpres_sqr(num[i], tEP[i].x, n);
	    mpres_mul_ui(num[i], num[i], 3, n);
	    mpres_add(num[i], num[i], tEP[i].A, n);
	    mpres_mul_ui(den[i], tEP[i].y, 2, n);
	}
    }
    res = pt_many_common(tEQ, tEP, tEP, nEP, n, num, den, inv, takeit);
    /* TODO: case takeit[i] == 2 */
    free(takeit);
    return res;
}

/* R[i] <- P[i] + Q[i], or a factor is found which is put in tER[0]->x. */
int
pt_many_add(curve *tER, curve *tEP, curve *tEQ, int nEP, mpmod_t n, 
	    mpres_t *num, mpres_t *den, mpres_t *inv, char *ok)
{
    char *takeit = (char *)malloc(nEP * sizeof(char));
    int i, res;

    memcpy(takeit, ok, nEP);
#if DEBUG_MANY_EC >= 2
    printf("In pt_many_add, adding\n");
    pt_many_print(tEP, nEP, n);
    printf("and\n");
    pt_many_print(tEQ, nEP, n);
#endif
    for(i = 0; i < nEP; i++){
	if(ok[i] == 0)
	    continue; /* takeit[i] = 0 */
	if(pt_is_zero(tEP+i, n)){
#if DEBUG_MANY_EC >= 2
	    printf("# tEP[%d] = O_{E[%d]}\n", i, i);
#endif
	    takeit[i] = 0;
	    pt_assign(tER+i, tEQ+i, n);
	}
	else if(pt_is_zero(tEQ+i, n)){
#if DEBUG_MANY_EC >= 2
	    printf("# tEQ[%d] = O_{E[%d]}\n", i, i);
#endif
	    takeit[i] = 0;
	    pt_assign(tER+i, tEP+i, n);
	}
	else if(pt_is_equal(tEP+i, tEQ+i)){
	    /* we should double */
	    if(mpz_sgn(tEP[i].y) == 0){
#if DEBUG_MANY_EC >= 2
		printf("# 2 * P[%d] = O_{E[%d]}\n", i, i);
#endif
		takeit[i] = 0;
		pt_set_to_zero(tEP+i, n);
	    }
	    else{
		/* ordinary doubling */
		mpres_sqr(num[i], tEP[i].x, n);
		mpres_mul_ui(num[i], num[i], 3, n);
		mpres_add(num[i], num[i], tEP[i].A, n);
		mpres_mul_ui(den[i], tEP[i].y, 2, n);
	    }
	}
	else if(mpz_cmp(tEQ[i].x, tEP[i].x) == 0){
	    mpres_add(num[i], tEQ[i].x, tEP[i].x, n);
	    if(mpz_sgn(num[i]) == 0){
		takeit[i] = 0;
		pt_set_to_zero(tER+i, n);
	    }
	}
	else{
	    mpres_sub(num[i], tEQ[i].y, tEP[i].y, n);
	    mpres_sub(den[i], tEQ[i].x, tEP[i].x, n);
	}
    }
    res = pt_many_common(tER, tEP, tEQ, nEP, n, num, den, inv, takeit);
    /* TODO: case takeit[i] == 2 */
    free(takeit);
    return res;
}

/* tEQ[i] <- e * tEP[i]; we must have tEQ != tEP */
/* If a factor is found, it is put back in num[nEP]. */
int
pt_many_mul(curve *tEQ, curve *tEP, int nEP, mpz_t e, mpmod_t n, 
	    mpres_t *num, mpres_t *den, mpres_t *inv, char *ok)
{
  size_t l;
  int negated = 0, status = 1;

  if (mpz_sgn (e) == 0)
    {
	pt_many_set_to_zero(tEQ, nEP, n);
	return 1;
    }

  /* The negative of a point (x:y:z) is (x:-y:z) */
  if (mpz_sgn (e) < 0)
    {
      negated = 1;
      mpz_neg (e, e);
      pt_many_neg(tEP, nEP, n);
    }

  if (mpz_cmp_ui (e, 1) == 0)
    goto pt_many_mul_end;

  l = mpz_sizeinbase (e, 2) - 1; /* l >= 1 */

  pt_many_assign(tEQ, tEP, nEP, n);

  while (l-- > 0)
    {
	if(pt_many_duplicate (tEQ, tEQ, nEP, n, num, den, inv, ok) == 0)
	  {
	    status = 0;
	    break;
	  }
#if DEBUG_MANY_EC >= 2
	printf("Rdup:="); pt_many_print(tEQ, nEP, n); printf(";\n");
#endif
	if (mpz_tstbit (e, l))
	  {
	      if(pt_many_add (tEQ, tEP, tEQ, nEP, n, num, den, inv, ok) == 0)
	      {
		status = 0;
		break;
	      }
#if DEBUG_MANY_EC >= 2
	      printf("Radd:="); pt_many_print(tEQ, nEP, n); printf(";\n");
#endif
	  }
    }

pt_many_mul_end:

  /* Undo negation to avoid changing the caller's e value */
  if (negated){
    mpz_neg (e, e);
    pt_many_neg(tEP, nEP, n);
  }
  return status;
}

/* r <- q mod N. 
   Return value: 1 if den invertible, 0 if factor found; in this case
   gcd(den(q), N) is put in r.
 */
int
mod_from_rat(mpz_t r, mpq_t q, mpz_t N)
{
    int ret = 1;
 
    if(mpz_invert(r, mpq_denref (q), N) == 0){
	mpz_gcd(r, mpq_denref (q), N);
	ret = 0;
    }
    else{
	mpz_mul(r, r, mpq_numref(q));
	mpz_mod(r, r, N);
    }
    return ret;
}

/* fall back on traditional ECM.
 */
int
process_one_curve(mpz_t f, mpz_t N, double B1, ecm_params params, curve EP)
{
    double B2scale = 1.0;
    int ret;

    /* Taken from main.c; no comment */
    /* Here's an ugly hack to pass B2scale to the library somehow.
       It gets piggy-backed onto B1done */
    params->B1done = params->B1done + floor (B2scale * 128.) / 134217728.; 

    mpz_set_si(params->B2, ECM_DEFAULT_B2); /* compute it automatically from B1 */
    mpz_set_si(params->B2min, ECM_DEFAULT_B2); /* will be set to B1 */
#if DEBUG_MANY_EC >= 2
    params->verbose = 2;
#endif
    mpz_set(params->x, EP.x);
    mpz_set(params->y, EP.y);
    mpz_set(params->sigma, EP.A); /* humf */
    ret = ecm_factor(f, N, B1, params);
    return ret;
}

int
conclude_on_factor(mpz_t N, mpz_t f)
{
    mpz_t C;
    int factor_is_prime, cofactor_is_prime, ret;

    if(mpz_cmp(N, f) == 0){
	printf("# found input number, proceeding to next curve\n");
	return ECM_INPUT_NUMBER_FOUND;
    }
    factor_is_prime = mpz_probab_prime_p (f, PROBAB_PRIME_TESTS);
    mpz_init(C);
    mpz_tdiv_q(C, N, f);
    cofactor_is_prime = mpz_probab_prime_p (C, PROBAB_PRIME_TESTS);
    if (factor_is_prime)
	ret = cofactor_is_prime ? ECM_PRIME_FAC_PRIME_COFAC :
	    ECM_PRIME_FAC_COMP_COFAC;
    else
	ret = cofactor_is_prime ? ECM_COMP_FAC_PRIME_COFAC :
	    ECM_COMP_FAC_COMP_COFAC;
    gmp_printf("f=%Zd is ", f);
    if(factor_is_prime)
	printf("PRP\n");
    else
	printf("COMPOSITE\n");
    gmp_printf("C=%Zd is ", C);
    if(cofactor_is_prime)
	printf("PRP\n");
    else
	printf("COMPOSITE\n");
    mpz_clear(C);
    return ret;
}

#ifdef ONE_CURVE_AT_A_TIME
int
one_curve_at_a_time(mpz_t f, char *ok, curve *tEP, int nc, mpz_t N, double B1)
{
    ecm_params params;
    int ret = 0, i;
    mpz_t C;

    ecm_init(params);
    params->sigma_is_A = -1;
    mpz_init (C);
    /* process curves one at a time */
    for(i = 0; i < nc; i++){
	params->B1done = 1.0;
	ret = process_one_curve(f, N, B1, params, tEP[i]);
	if(ret > 0){ /* humf */
	    ok[i] = 0;
	    ret = conclude_on_factor(N, f);
	}
	else if(ret == ECM_ERROR){
	    printf("Error for curve %d\n", i);
	}
    }
    mpz_clear (C);
    ecm_clear(params);
    return ret;
}
#else
/* Using parallelism.
   Copied from classical ecm_stage1.
 */
int
all_curves_at_once(mpz_t f, char *ok, curve *tEP, int nEP, mpmod_t n, 
		   double B1, double *B1done, int (*stop_asap)(void),
		   char *chkfilename)
{
    curve tEQ[NCURVE_MAX], tER[NCURVE_MAX];
    mpz_t num[NCURVE_MAX+1], den[NCURVE_MAX+1], inv[NCURVE_MAX], e;
    double p = 0.0, r, last_chkpnt_p;
    int ret = ECM_NO_FACTOR_FOUND;
    long last_chkpnt_time;
    int i;
    
    mpz_init(e);
    for(i = 0; i < nEP; i++){
	mpres_init(tEQ[i].A, n); mpres_set(tEQ[i].A, tEP[i].A, n);
	mpres_init(tEQ[i].x, n); mpres_set(tEQ[i].x, tEP[i].x, n);
	mpres_init(tEQ[i].y, n); mpres_set(tEQ[i].y, tEP[i].y, n);
	mpres_init(tEQ[i].z, n); mpres_set(tEQ[i].z, tEP[i].z, n);
	
	mpres_init(tER[i].A, n); mpres_set(tER[i].A, tEP[i].A, n);
	mpres_init(tER[i].x, n);
	mpres_init(tER[i].y, n);
	mpres_init(tER[i].z, n);
	
	mpres_init(num[i], n);
	mpres_init(den[i], n);
	mpres_init(inv[i], n);
    }
    mpres_init(num[nEP], n); /* to be used as buffer in compute_all_inverses */
    mpres_init(den[nEP], n); /* to be used as buffer in compute_all_inverses */
    
    last_chkpnt_time = cputime ();
    
#if DEBUG_MANY_EC >= 2
    printf("Initial points:\n");
    pt_many_print(tEP, nEP, n);
#endif
    for (r = 2.0; r <= B1; r *= 2.0)
	if (r > *B1done){
	    if(pt_many_duplicate (tEQ, tEQ, nEP, n, num, den, inv, ok) == 0){
		mpz_set(f, num[nEP]);
		ret = ECM_FACTOR_FOUND_STEP1;
		goto end_of_all;
	    }
#if DEBUG_MANY_EC >= 2
	    printf("P%ld:=", (long)r); pt_many_print(tEQ, nEP, n); printf(";\n");
#endif
	}

    last_chkpnt_p = 3.;
    for (p = getprime (); p <= B1; p = getprime ()){
	for (r = p; r <= B1; r *= p){
#if DEBUG_MANY_EC >= 1
	    printf("## p = %ld\n", (long)p);
#endif
	    if (r > *B1done){
		mpz_set_ui(e, (ecm_uint) p);
		if(pt_many_mul(tER, tEQ, nEP, e, n, num, den, inv, ok) == 0){
		    mpz_set(f, num[nEP]);
		    ret = ECM_FACTOR_FOUND_STEP1;
		    goto end_of_all;
		}
#if DEBUG_MANY_EC >= 2
		pt_many_print(tER, nEP, n);
#endif
		for(i = 0; i < nEP; i++)
		    if(pt_is_zero(tER+i, n))
			ok[i] = 0;
		pt_many_assign(tEQ, tER, nEP, n); /* TODO: use pointers */
	    }
	    if (stop_asap != NULL && (*stop_asap) ()){
		outputf (OUTPUT_NORMAL, "Interrupted at prime %.0f\n", p);
		break;
	    }
	    
#if 0
	    /* WARNING: not activated yet */
	    if (chkfilename != NULL && p > last_chkpnt_p + 10000. && 
		elltime (last_chkpnt_time, cputime ()) > CHKPNT_PERIOD){
		writechkfile (chkfilename, ECM_ECM, MAX(p, *B1done), n, A, x, y, z);
		last_chkpnt_p = p;
		last_chkpnt_time = cputime ();
	    }
#endif
	}
    }
 end_of_all:
    /* If stage 1 finished normally, p is the smallest prime > B1 here.
       In that case, set to B1 */
    if (p > B1)
	p = B1;
    
    if (p > *B1done)
	*B1done = p;
    
#if 0
    if (chkfilename != NULL)
	writechkfile (chkfilename, ECM_ECM, *B1done, n, A, x, y, z);
#endif
    getprime_clear (); /* free the prime tables, and reinitialize */
    
    /* put results back */
    pt_many_assign(tEP, tEQ, nEP, n);
    /* normalize all points */
    for(i = 0; i < nEP; i++)
	if(pt_is_zero(tEP+i, n))
	    pt_set_to_zero(tEP+i, n);
    /* clear temporary variables */
    mpz_clear(e);
    for(i = 0; i < nEP; i++){
	mpres_clear(tEQ[i].A, n);
	mpres_clear(tEQ[i].x, n);
	mpres_clear(tEQ[i].y, n);
	mpres_clear(tEQ[i].z, n);
	mpres_clear(num[i], n);
	mpres_clear(den[i], n);
	mpres_clear(inv[i], n);
    }
    mpres_clear(num[nEP], n);
    mpres_clear(den[nEP], n);
    return ret;
}
#endif

int
read_and_prepare(mpz_t f, mpz_t x, mpq_t q, char *buf, mpz_t n)
{
    mpq_set_str(q, buf, 10);
    if(mod_from_rat(x, q, n) == 0){
	mpz_set(f, x);
	return 0;
    }
    return 1;
}

int
process_many_curves(mpz_t f, mpz_t n, double B1, curve *tEP, int nEP)
{
#ifndef ONE_CURVE_AT_A_TIME
    mpmod_t modulus;
    double B1done;
    curve tEQ[NCURVE_MAX];
#endif
    ecm_params params;
    char *ok = (char *)malloc(nEP * sizeof(char));
    int ret = 0, i;
    long st = cputime ();
    
    memset(ok, 1, nEP);
#ifdef ONE_CURVE_AT_A_TIME
    ret = one_curve_at_a_time(f, ok, tEP, nEP, n, B1);
#else
    mpmod_init(modulus, n, ECM_MOD_DEFAULT);
    for(i = 0; i < nEP; i++){
       mpres_init(tEQ[i].x, modulus); mpres_set_z(tEQ[i].x, tEP[i].x, modulus);
       mpres_init(tEQ[i].y, modulus); mpres_set_z(tEQ[i].y, tEP[i].y, modulus);
       mpres_init(tEQ[i].z, modulus); mpres_set_z(tEQ[i].z, tEP[i].z, modulus);
       mpres_init(tEQ[i].A, modulus); mpres_set_z(tEQ[i].A, tEP[i].A, modulus);
    }
    B1done = 1.0;
    ret = all_curves_at_once(f, ok, tEQ, nEP, modulus, B1, &B1done, NULL, NULL);
    printf("# Step 1 took %ldms\n", elltime (st, cputime ()));

    if(ret != ECM_NO_FACTOR_FOUND){
	ret = conclude_on_factor(n, f);
    }
    else{
	ecm_init(params);
	params->sigma_is_A = -1;
	params->B1done = B1;
	for(i = 0; i < nEP; i++){
	    if(ok[i] == 0)
		continue;
#if DEBUG_MANY_EC >= 1
	    printf("# Entering Step 2 for E[%d]\n", i);
#endif
	    st = cputime ();
	    mpres_get_z(tEP[i].x, tEQ[i].x, modulus);
	    mpres_get_z(tEP[i].y, tEQ[i].y, modulus);
	    mpres_get_z(tEP[i].z, tEQ[i].z, modulus);
	    mpres_get_z(tEP[i].A, tEQ[i].A, modulus);
	    ret = process_one_curve(f, n, B1, params, tEP[i]);
	    printf("# Step 2 for E[%d] took %ldms\n",i,elltime(st, cputime()));
	    if(ret != ECM_NO_FACTOR_FOUND){
		printf("## factor found in Step 2: ");
		mpz_out_str (stdout, 10, f);
		printf ("\n");
		ret = conclude_on_factor(n, f);
		break;
	    }
	}
	ecm_clear(params);
    }
    for(i = 0; i < nEP; i++){
	mpres_clear(tEQ[i].x, modulus);
	mpres_clear(tEQ[i].y, modulus); 
	mpres_clear(tEQ[i].z, modulus); 
	mpres_clear(tEQ[i].A, modulus); 
    }
#endif
    free(ok);
    return ret;
}

int
process_many_curves_from_file(mpz_t f, mpz_t n, double B1, char *fic_EP,
			      int ncurves)
{
    curve tEP[NCURVE_MAX];
    FILE *ifile = fopen(fic_EP, "r");
    int nEP = 0, i, ret = 0;
    char bufA[1024], bufx[1024], bufy[1024];
    mpq_t q;

    mpq_init(q);
    while(fscanf(ifile, "%s %s %s", bufA, bufx, bufy) != EOF){
	mpz_init(tEP[nEP].A);
	if(read_and_prepare(f, tEP[nEP].A, q, bufA, n) == 0)
	    goto process_end;
	mpz_init(tEP[nEP].x);
	if(read_and_prepare(f, tEP[nEP].x, q, bufx, n) == 0)
	    goto process_end;
	mpz_init(tEP[nEP].y);
	if(read_and_prepare(f, tEP[nEP].y, q, bufy, n) == 0)
	    goto process_end;
	mpz_init_set_ui(tEP[nEP].z, 1);
	nEP++;
	if(ncurves != 0 && nEP == ncurves)
	    break;
    }
    ret = process_many_curves(f, n, B1, tEP, nEP);
    for(i = 0; i < nEP; i++){
	mpz_clear(tEP[i].x);
	mpz_clear(tEP[i].y);
	mpz_clear(tEP[i].z);
	mpz_clear(tEP[i].A);
    }
 process_end:
    fclose(ifile);
    mpq_clear(q);
    return ret;
}

/* Assuming we can generate curves with given torsion using parameter s
   in interval [smin..smax].
*/
int
build_curves_with_torsion(mpz_t f, mpz_t n, curve *tEP, char *torsion, 
			  int smin, int smax, int nEP)
{
    curve E;
    int ret = 0;

    /* over Q: see Atkin-Morain, Math. Comp., 1993 */
    if(strcmp(torsion, "Z5") == 0){
    }
    else if(strcmp(torsion, "Z7") == 0){
    }
    else if(strcmp(torsion, "Z9") == 0){
    }
    else if(strcmp(torsion, "Z10") == 0){
    }
    else if(strcmp(torsion, "Z2xZ8") == 0){
    }
    /* no longer over Q */
    else if(strcmp(torsion, "Z3xZ3") == 0){ /* over Q(sqrt(-3)) */
    }
    else{
	printf("Unknown torsion group: %s\n", torsion);
	ret = ECM_ERROR;
    }
    return ret;
}

int
process_curves_with_torsion(mpz_t f, mpz_t n, double B1, char *torsion,
			    int smin, int smax, int nEP)
{
    curve tEP[NCURVE_MAX];
    int ret = 0, i;

    for(i = 0; i < nEP; i++){
	mpz_init(tEP[i].A);
	mpz_init(tEP[i].x);
	mpz_init(tEP[i].y);
	mpz_init(tEP[i].z);
    }
    if(build_curves_with_torsion(f, n, tEP, torsion, smin, smax, nEP) == 0){
	/* a factor found by sheer luck */
	ret = 0;
    }
    else
	ret = process_many_curves(f, n, B1, tEP, nEP);
    for(i = 0; i < nEP; i++){
	mpz_clear(tEP[i].A);
	mpz_clear(tEP[i].x);
	mpz_clear(tEP[i].y);
	mpz_clear(tEP[i].z);
    }
    return ret;
}

static void
usage (char *cmd)
{
    printf("Usage: %s -inp file_N -B1 B1 -curves file_C", cmd);
    printf(" -torsion T -smin smin -smax smax\n");
    printf("  -inp file_N    numbers to be factored, one per line\n");
    printf("                 file_N can be '-', in which case stdin is used\n");
    printf("  -curves file_C curves to be used, format 'A x0 y0' per line\n");
    printf("  -h, --help   Prints this help and exit.\n");
}

int
main (int argc, char *argv[])
{
  mpz_t n, f;
  int res = 0, smin = -1, smax = -1, ncurves = 0;
  double B1 = 0.0;
  char *infilename = NULL, *curvesname = NULL, *torsion = NULL;
  char buf[10000];
  FILE *infile = NULL;

  /* first look for options */
  while ((argc > 1) && (argv[1][0] == '-')){
      if (strcmp (argv[1], "-h") == 0 || strcmp (argv[1], "--help") == 0){
          usage (argv[0]);
          exit (EXIT_SUCCESS);
      }
      else if ((argc > 2) && (strcmp (argv[1], "-B1") == 0)){
	  B1 = atof(argv[2]);
	  argv += 2;
	  argc -= 2;
      }
      else if ((argc > 2) && (strcmp (argv[1], "-inp") == 0)){
	  infilename = argv[2];
	  if(strcmp(infilename, "-") == 0)
	      infile = stdin;
	  else{
	      infile = fopen (infilename, "r");
	      if (!infile){
		  fprintf (stderr, "Can't find input file %s\n", infilename);
		  exit (EXIT_FAILURE);
	      }
	  }
	  argv += 2;
	  argc -= 2;
      }
      else if ((argc > 2) && (strcmp (argv[1], "-curves") == 0)){
	  curvesname = argv[2];
	  argv += 2;
	  argc -= 2;
      }
      /* one may restrict the number of curves used from file */
      else if ((argc > 2) && (strcmp (argv[1], "-ncurves") == 0)){
	  ncurves = atoi(argv[2]);
	  argv += 2;
	  argc -= 2;
      }
      /** torsion related parameters **/
      else if ((argc > 2) && (strcmp (argv[1], "-torsion") == 0)){
	  torsion = argv[2];
	  argv += 2;
	  argc -= 2;
      }
      else if ((argc > 2) && (strcmp (argv[1], "-smin") == 0)){
	  smin = atoi(argv[2]);
	  argv += 2;
	  argc -= 2;
      }
      else if ((argc > 2) && (strcmp (argv[1], "-smax") == 0)){
	  smax = atoi(argv[2]);
	  argv += 2;
	  argc -= 2;
      }
      else{
	  fprintf (stderr, "Unknown option: %s\n", argv[1]);
	  exit (EXIT_FAILURE);
      }
  }
  if(infile == NULL){
      fprintf (stderr, "No input file given\n");
      exit (EXIT_FAILURE);
  }
  if(curvesname == NULL && torsion == NULL){
      fprintf (stderr, "No curve file given\n");
      exit (EXIT_FAILURE);
  }
  if(curvesname != NULL && torsion != NULL){
      fprintf (stderr, "Cannot have -curves and -torsion at the same time.\n");
      exit (EXIT_FAILURE);
  }
  if(torsion != NULL && ncurves == 0){
      fprintf (stderr, "You must provide ncurves != 0 with -torsion.\n");
      exit (EXIT_FAILURE);
  }

  mpz_init (n);
  mpz_init (f); /* for potential factor */
  while(fscanf(infile, "%s", buf) != EOF){
      /* read number */
      if(buf[0] == '#'){
	  /* print till end of line */
	  printf("%s", buf);
	  char c;
	  while((c = getc(infile)) != '\n')
	      printf("%c", c);
	  printf("\n");
	  continue;
      }
      if(mpz_set_str (n, buf, 10)){
	  fprintf (stderr, "Invalid number: %s\n", argv[1]);
	  exit (1);
      }
      if(curvesname != NULL){
	  if(ncurves == 0)
	      printf ("# Using all");
	  else
	      printf("# Using only %d", ncurves);
	  printf(" curves from %s with B1=%1.0f\n", curvesname, B1);
	  res = process_many_curves_from_file(f, n, B1, curvesname, ncurves);
      }
      else if(torsion != NULL){
	  res = process_curves_with_torsion(f, n, B1, torsion, smin, smax, ncurves);
      }
      fflush(stdout);
  }
  if(infile != stdin)
      fclose(infile);

  mpz_clear (f);
  mpz_clear (n);

  return 0;
}
