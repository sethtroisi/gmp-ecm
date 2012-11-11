/* manyecm.c - ECM with many curves in parallel 
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

#define DEBUG_MANY_EC 0

#define NCURVE_MAX 1000

/********** stolen by lazyness **********/

/* returns the number of decimal digits of n */
unsigned int
nb_digits (const mpz_t n)
{
  mpz_t x;
  unsigned int size;

  size = mpz_sizeinbase (n, 10);

  /* the GMP documentation says mpz_sizeinbase returns the exact value,
     or one too big, thus:
     (a) either n < 10^(size-1), and n has size-1 digits
     (b) or n >= size-1, and n has size digits
     Note: mpz_sizeinbase returns 1 for n=0, thus we always have size >= 1.
  */
				    
  mpz_init (x);
  mpz_ui_pow_ui (x, 10, size - 1);
  if (mpz_cmpabs (n, x) < 0)
    size --;
  mpz_clear (x);

  return size;
}

/********** group law on points **********/

#define pt_is_equal(EP, EQ) (mpz_cmp((EP)->x, (EQ)->x) == 0 \
	                     && mpz_cmp((EP)->y, (EQ)->y) == 0 \
			     && mpz_cmp((EP)->z, (EQ)->z) == 0)

int
pt_is_zero(curve *EP, ATTRIBUTE_UNUSED mpmod_t n)
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
pt_assign(curve *EQ, curve *EP, ATTRIBUTE_UNUSED mpmod_t n)
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

/* NOTE: we can have tER = tEP or tEQ */
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
		pt_set_to_zero(tER+i, n);
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

/* fall back on traditional ECM.
   We decide between Weierstrass and Montgomery by inspection of y: it is
   Weierstrass iff y == NULL.
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
    mpz_set(params->x, EP.x);
    mpz_set(params->sigma, EP.A); /* humf */

    if(EP.y->_mp_alloc == 0) /* humf */
	params->sigma_is_A = 1;
    else{
	params->sigma_is_A = -1;
	mpz_set(params->y, EP.y);
    }
    ret = ecm_factor(f, N, B1, params);
    return ret;
}

int
conclude_on_factor(mpz_t N, mpz_t f, int verbose)
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
    if (verbose >= 1)
      {
        printf ("Found %s factor of %2u digits: ", 
		factor_is_prime ? "probable prime" : "composite",
		nb_digits (f));
	mpz_out_str (stdout, 10, f);
	printf ("\n");
	printf ("%s cofactor ",
		cofactor_is_prime ? "Probable prime" : "Composite");
	mpz_out_str (stdout, 10, C);
	printf (" has %u digits\n", nb_digits(C));
      }
    mpz_clear(C);
    return ret;
}

int
one_curve_at_a_time(mpz_t f, char *ok, curve *tEP, int nc, mpz_t N, double B1)
{
    ecm_params params;
    int ret = 0, i;
    mpz_t C;

    ecm_init(params);
    params->verbose = 1;
    mpz_init (C);
    /* process curves one at a time */
    for(i = 0; i < nc; i++){
	params->B1done = 1.0;
	ret = process_one_curve(f, N, B1, params, tEP[i]);
	if(ret > 0){ /* humf */
	    ok[i] = 0;
	    ret = conclude_on_factor(N, f, params->verbose);
	    break;
	}
	else if(ret == ECM_ERROR){
	    printf("Error for curve %d\n", i);
	}
    }
    mpz_clear (C);
    ecm_clear(params);
    return ret;
}

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
	    printf("## p = %ld at %ldms\n", (long)p, cputime());
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

/* f is a (probable) prime factor of n. tEP is in plain mod n form. */
void
dump_curves(curve *tEP, int nEP, mpz_t f)
{
    int i;

    gmp_printf("p:=%Zd; F:=GF(p); P:=[]; A:=[]; B:=[]; E:=[];\n", f);
    for(i = 0; i < nEP; i++){
	gmp_printf("P[%d]:=[%Zd, %Zd, %Zd];\n", i+1, 
		   tEP[i].x, tEP[i].y, tEP[i].z); 
	gmp_printf("A[%d]:=%Zd;\n", i+1, tEP[i].A);
	printf("B[%d]:=P[%d][2]^2-P[%d][1]^3-A[%d]*P[%d][1];\n", 
	       i+1, i+1, i+1, i+1, i+1);
	printf("E[%d]:=EllipticCurve([F!A[%d], F!B[%d]]);\n", i+1, i+1, i+1);
	printf("Factorization(#E[%d]);\n", i+1);
    }
}

int
process_many_curves(mpz_t f, mpz_t n, double B1, curve *tEP, int nEP,
		    int onebyone)
{
    mpmod_t modulus;
    double B1done;
    curve tEQ[NCURVE_MAX];
    ecm_params params;
    char *ok = (char *)malloc(nEP * sizeof(char));
    int ret = 0, i;
    long st = cputime ();
    
    ecm_init(params);
#if DEBUG_MANY_EC >= 2
    params->verbose = 2;
#else
    params->verbose = 1;
#endif
    memset(ok, 1, nEP);
    if(onebyone)
	return one_curve_at_a_time(f, ok, tEP, nEP, n, B1);
    /* take everybody */
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
	ret = conclude_on_factor(n, f, params->verbose);
#if DEBUG_MANY_EC >= 2
	if(ret == ECM_PRIME_FAC_PRIME_COFAC || ret == ECM_PRIME_FAC_COMP_COFAC)
	    /* output Magma lines to check #E's mod f */
	    dump_curves(tEP, nEP, f);
#endif
    }
    else{
	params->sigma_is_A = -1;
	params->B1done = B1;
	for(i = 0; i < nEP; i++){
	    if(ok[i] == 0)
		continue;
#if DEBUG_MANY_EC >= 1
	    printf("# Entering Step 2 for E[%d]\n", i);
#endif
	    mpres_get_z(tEP[i].x, tEQ[i].x, modulus);
	    mpres_get_z(tEP[i].y, tEQ[i].y, modulus);
	    mpres_get_z(tEP[i].z, tEQ[i].z, modulus);
	    mpres_get_z(tEP[i].A, tEQ[i].A, modulus);
	    ret = process_one_curve(f, n, B1, params, tEP[i]);
	    if(ret != ECM_NO_FACTOR_FOUND){
		printf("## factor found in Step 2: ");
		mpz_out_str (stdout, 10, f);
		printf ("\n");
		ret = conclude_on_factor(n, f, params->verbose);
		break;
	    }
	}
    }
    for(i = 0; i < nEP; i++){
	mpres_clear(tEQ[i].x, modulus);
	mpres_clear(tEQ[i].y, modulus); 
	mpres_clear(tEQ[i].z, modulus); 
	mpres_clear(tEQ[i].A, modulus); 
    }
    ecm_clear(params);
    mpmod_clear(modulus);
    free(ok);
    return ret;
}

int
process_many_curves_from_file(mpz_t f, mpz_t n, double B1, char *fic_EP,
			      int ncurves, int onebyone)
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
    ret = process_many_curves(f, n, B1, tEP, nEP, onebyone);
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

/* A = (a4-1/3*a2^2)
   B = -1/3*a4*a2+a6+2/27*a2^3 (unused, really)
   X = x+a2/3
*/
void
W2W(curve *EP, mpz_t a2, mpz_t a4, mpz_t x0, mpz_t y0, mpz_t n)
{
    /* x <- a2/3 */
    mpz_set_si(EP->y, 3);
    mod_from_rat2(EP->x, a2, EP->y, n);
    /* A = a4-1/3*a2^2 = a4 - a2 * (a2/3) */
    /** a2 <- a2^2/3 **/
    mpz_mul(a2, a2, EP->x);
    mpz_mod(a2, a2, n);
    mpz_sub(EP->A, a4, a2);
    mpz_mod(EP->A, EP->A, n);
    /* wx0 = x0 + a2/3 */
    mpz_add(EP->x, EP->x, x0);
    mpz_mod(EP->x, EP->x, n);
    mpz_set(EP->y, y0);
    mpz_mod(EP->y, EP->y, n);
    mpz_set_ui(EP->z, 1);
#if DEBUG_MANY_EC >= 2
    gmp_printf("N:=%Zd;\n", n);
    gmp_printf("A:=%Zd;\n", EP->A);
    gmp_printf("x0:=%Zd;\n", x0);
    gmp_printf("wx0:=%Zd;\n", EP->x);
    gmp_printf("y0:=%Zd;\n", EP->y);
    printf("B:=(y0^2-wx0^3-A*wx0) mod N;\n");
    exit(-1);
#endif
}

/* From a curve in Weierstrass form to a short form 
   WE:=[0,(1/4*c^2+1/4-1/2*c-b),0,(1/2*c*b-1/2*b),1/4*b^2]);
   We compute:
   a2 = 1/4*c^2+1/4-1/2*c-b = ((c-1)/2)^2-b
   a4 = 1/2*c*b-1/2*b = b*(c-1)/2
   a6 = (b/2)^2
*/
void
K2W24(mpz_t a2, mpz_t a4, mpz_t b, mpz_t c, mpz_t n)
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
#if 0 /* not needed */
    {
	mpz_init(a6);
	mpz_set(a6, b);
	mod_div_2(a6, n);
	mpz_mul(a6, a6, a6);
	mpz_mod(a6, a6, n);
	mpz_clear(a6);
    }
#endif
#if DEBUG_MANY_EC >= 2
    gmp_printf("N:=%Zd;\n", n);
    gmp_printf("b:=%Zd;\n", b);
    gmp_printf("c:=%Zd;\n", c);
    gmp_printf("a2:=%Zd;\n", a2);
    gmp_printf("a4:=%Zd;\n", a4);
    printf("a6:=RatMod(b^2/4, N);\n");
#endif
}

/* 
   Send Kubert curve with point (x0, y0) to short Weierstrass form:
   Y^2 = X^3 + A * X + B
*/
void
K2W4(curve *EP, mpz_t b, mpz_t c, mpz_t x0, mpz_t y0, mpz_t n)
{
    mpz_t a2, a4;

    mpz_init(a2);
    mpz_init(a4);
    K2W24(a2, a4, b, c, n);
    /* second conversion */
    W2W(EP, a2, a4, x0, y0, n);
    mpz_clear(a2);
    mpz_clear(a4);
}

/* Kubert: put b = c. */
int
build_curves_with_torsion_Z5(mpz_t f, mpz_t n, curve *tEP,
			     int smin, int smax, int nEP)
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
	if(mod_from_rat2(c, tmp, y0, n) == 0){
	    /* factor found! */
	    mpz_gcd(f, c, n);
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
	/* y0:=x0*(x0+1)*(4*x0+1)/4/(3*x0+1) = (x0+1)*c/2 */
	mpz_add_si(y0, x0, 1);
	mpz_mul(y0, y0, c);
	mpz_mod(y0, y0, n);
	mod_div_2(y0, n);
#if DEBUG_MANY_EC >= 2
	gmp_printf("x0:=%Zd;\nc:=%Zd;\ny0:=%Zd;\n", x0, c, y0);
	printf("cr:=1/2*x0*(4*x0+1)/(3*x0+1);\n");
#endif
	/* P:=WE![x0, y0, 1]; */
	/* convert to short Weierstrass form */
	K2W4(tEP+nc, c, c, x0, y0, n);
	nc++;
	if(nc >= nEP)
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

int
build_curves_with_torsion_Z7(mpz_t f, mpz_t n, curve *tEP,
			     int umin, int umax, int nEP)
{
    int u, ret = ECM_NO_FACTOR_FOUND, nc = 0;
    mpz_t A2, A1div2, x0, y0, cte, num[2], den[2], inv[1], d, c, b, kx0, ky0;
    mpmod_t modulus;
    curve EP[1], EQ[1]; /* blourk */
    char ok[1];

    ok[0] = 1;
    mpmod_init(modulus, n, ECM_MOD_DEFAULT);
    for(u = 0; u < 2; u++){
	mpres_init(num[u], modulus);
	mpres_init(den[u], modulus);
    }
    mpres_init(inv[0], modulus);
    /* Eaux = "1295/48", "-1079/864" */
    /* Paux = "2185/12", "-2458" */
    mpres_init(EP[0].A, modulus);
    mod_from_rat_str(f, "1295/48", n); mpres_set_z(EP[0].A, f, modulus);
    mpres_init(EP[0].x, modulus);
    mod_from_rat_str(f, "2185/12", n); mpres_set_z(EP[0].x, f, modulus);
    mpres_init(EP[0].y, modulus);
    mpz_set_str(f, "-2458", 10); mpres_set_z(EP[0].y, f, modulus);
    mpres_init(EP[0].z, modulus);
    mpres_set_ui(EP[0].z, 1, modulus);
#if DEBUG_MANY_EC >= 2
    printf("P:=");
    pt_print(EP[0], modulus);
    printf(";\n");
#endif

    mpres_init(EQ[0].x, modulus);
    mpres_init(EQ[0].y, modulus);
    mpres_init(EQ[0].z, modulus);
    mpres_init(EQ[0].A, modulus);
    mpres_set(EQ[0].A, EP[0].A, modulus);

    mpz_init(A2);
    mod_from_rat_str(A2, "1/12", n);
    mpz_init_set_str(A1div2, "-1", 10);
    mpz_mod(A1div2, A1div2, n);
    mpz_init_set_str(x0, "-1", 10);
    mpz_mod(x0, x0, n);
    mpz_init_set_str(y0, "8", 10);
    mpz_mod(y0, y0, n);
    mpz_init(cte);
    mod_from_rat_str(cte, "-7/2", n);

    mpz_init(d);
    mpz_init(c);
    mpz_init(b);
    mpz_init(kx0);
    mpz_init(ky0);
    for(u = umin; u < umax; u++){
	/* update Qaux */
	mpz_set_ui(d, u);
	if(pt_many_mul(EQ, EP, 1, d, modulus, num, den, inv, ok) == 0){
	    printf("found factor during update of Q\n");
	    mpz_set(f, num[1]);
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
#if DEBUG_MANY_EC >= 2
	printf("(s, t)[%d]:=", u);
	pt_print(EQ[0], modulus);
	printf(";\n");
#endif
	mpres_get_z(b, EQ[0].x, modulus);
	mpres_get_z(c, EQ[0].y, modulus);
	if(cubic_to_quartic(f, n, d, ky0, b, c, A2, A1div2, x0, y0, cte) == 0){
	    printf("found factor in Z7 (cubic_2_quartic)\n");
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
	/* d:=x; */
	/* x0:=-2*d; */
	mpz_mul_si(kx0, d, -2);
	mpz_mod(kx0, kx0, n);
	/* y0:=d*y/2; */
	mpz_mul(ky0, ky0, d);
	mpz_mod(ky0, ky0, n);
	mod_div_2(ky0, n);
	/* c:=d^2-d; */
	mpz_mul(c, d, d);
	mpz_sub(c, c, d);
	mpz_mod(c, c, n);
	/* b:=c*d; */
	mpz_mul(b, c, d);
	mpz_mod(b, b, n);
	K2W4(tEP+nc, b, c, kx0, ky0, n);
	nc++;
	if(nc >= nEP)
	    break;
	pt_many_assign(EP, EQ, 1, modulus);
    }
#if DEBUG_MANY_EC >= 2
    printf("Curves built\n");
    pt_many_print(tEP, nEP, modulus);
#endif
    mpz_clear(A2);
    mpz_clear(A1div2);
    mpz_clear(x0);
    mpz_clear(y0);
    mpz_clear(cte);
    mpres_clear(EP[0].A, modulus);
    mpres_clear(EP[0].x, modulus);
    mpres_clear(EP[0].y, modulus);
    mpres_clear(EP[0].z, modulus);
    for(u = 0; u < 2; u++){
	mpres_clear(num[u], modulus);
	mpres_clear(den[u], modulus);
    }
    mpres_clear(inv[0], modulus);
    mpres_clear(EP[0].x, modulus);
    mpres_clear(EP[0].y, modulus);
    mpres_clear(EP[0].z, modulus);
    mpres_clear(EP[0].A, modulus);
    mpres_clear(EQ[0].x, modulus);
    mpres_clear(EQ[0].y, modulus);
    mpres_clear(EQ[0].z, modulus);
    mpres_clear(EQ[0].A, modulus);
    mpmod_clear(modulus);
    mpz_clear(d);
    mpz_clear(c);
    mpz_clear(b);
    mpz_clear(kx0);
    mpz_clear(ky0);
    return ret;
}

int
build_curves_with_torsion_Z9(mpz_t fac, mpz_t n, curve *tEP,
			     int umin, int umax, int nEP)
{
    int u, ret = ECM_NO_FACTOR_FOUND, nc = 0;
    mpz_t A2, A1div2, x0, y0, cte, num[2], den[2], inv[1], d, c, b, kx0, ky0;
    mpz_t f;
    mpmod_t modulus;
    curve EP[1], EQ[1]; /* blourk */
    char ok[1];

    ok[0] = 1;
    mpmod_init(modulus, n, ECM_MOD_DEFAULT);
    for(u = 0; u < 2; u++){
	mpres_init(num[u], modulus);
	mpres_init(den[u], modulus);
    }
    mpres_init(inv[0], modulus);
    /* Eaux = [-9, 9] */
    /* Paux = [1, 1, 1] */
    mpz_init(f);
    mpres_init(EP[0].A, modulus);
    mpz_set_str(f, "-9", 10); mpres_set_z(EP[0].A, f, modulus);
    mpres_init(EP[0].x, modulus);
    mpz_set_str(f, "1", 10); mpres_set_z(EP[0].x, f, modulus);
    mpres_init(EP[0].y, modulus);
    mpz_set_str(f, "1", 10); mpres_set_z(EP[0].y, f, modulus);
    mpres_init(EP[0].z, modulus);
    mpres_set_ui(EP[0].z, 1, modulus);
#if DEBUG_MANY_EC >= 2
    printf("P:=");
    pt_print(EP[0], modulus);
    printf(";\n");
#endif

    mpres_init(EQ[0].x, modulus);
    mpres_init(EQ[0].y, modulus);
    mpres_init(EQ[0].z, modulus);
    mpres_init(EQ[0].A, modulus);
    mpres_set(EQ[0].A, EP[0].A, modulus);

    mpz_init_set_ui(A2, 0);
    mpz_init_set_str(A1div2, "3", 10);
    mpz_mod(A1div2, A1div2, n);
    mpz_init_set_str(x0, "2", 10);
    mpz_mod(x0, x0, n);
    mpz_init_set_str(y0, "3", 10);
    mpz_mod(y0, y0, n);
    mpz_init_set_ui(cte, 0);

    mpz_init(d);
    mpz_init(c);
    mpz_init(b);
    mpz_init(kx0);
    mpz_init(ky0);
    for(u = umin; u < umax; u++){
	/* update Qaux */
	mpz_set_ui(d, u);
	if(pt_many_mul(EQ, EP, 1, d, modulus, num, den, inv, ok) == 0){
	    printf("found factor during update of Q\n");
	    mpz_set(fac, num[1]);
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
#if DEBUG_MANY_EC >= 2
	printf("(s, t)[%d]:=", u);
	pt_print(EQ[0], modulus);
	printf(";\n");
#endif
	mpres_get_z(b, EQ[0].x, modulus);
	mpres_get_z(c, EQ[0].y, modulus);
	if(cubic_to_quartic(fac, n, f, ky0, b, c, A2, A1div2, x0, y0, cte) == 0){
	    printf("found factor in Z9 (cubic_2_quartic)\n");
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
	/* f:=x; */
	/* d:=f*(f-1)+1; */
	mpz_sub_si(d, f, 1);
	mpz_mul(d, d, f);
	mpz_add_si(d, d, 1);
	mpz_mod(d, d, n);
	/* c:=f*(d-1); */
	mpz_sub_si(c, d, 1);
	mpz_mul(c, c, f);
	mpz_mod(c, c, n);
	/* kx0:=(2*f-1)*f^2; */
	/** b <- f^2 **/
	mpz_mul(b, f, f);
	mpz_mod(b, b, n);
	mpz_mul_si(kx0, f, 2);
	mpz_sub_si(kx0, kx0, 1);
	mpz_mul(kx0, kx0, b);
	mpz_mod(kx0, kx0, n);
	/* ky0:=y*f^4/2; */
	/** b <- b^2 = f^4 **/
	mpz_mul(b, b, b);
	mpz_mul(ky0, ky0, b);
	mpz_mod(ky0, ky0, n);
	mod_div_2(ky0, n);
	/* b:=c*d; */
	mpz_mul(b, c, d);
	mpz_mod(b, b, n);
	K2W4(tEP+nc, b, c, kx0, ky0, n);
	nc++;
	if(nc >= nEP)
	    break;
	pt_many_assign(EP, EQ, 1, modulus);
    }
#if DEBUG_MANY_EC >= 2
    printf("Curves built\n");
    pt_many_print(tEP, nEP, modulus);
#endif
    mpz_clear(A2);
    mpz_clear(A1div2);
    mpz_clear(x0);
    mpz_clear(y0);
    mpz_clear(cte);
    mpres_clear(EP[0].A, modulus);
    mpres_clear(EP[0].x, modulus);
    mpres_clear(EP[0].y, modulus);
    mpres_clear(EP[0].z, modulus);
    for(u = 0; u < 2; u++){
	mpres_clear(num[u], modulus);
	mpres_clear(den[u], modulus);
    }
    mpres_clear(inv[0], modulus);
    mpres_clear(EP[0].x, modulus);
    mpres_clear(EP[0].y, modulus);
    mpres_clear(EP[0].z, modulus);
    mpres_clear(EP[0].A, modulus);
    mpres_clear(EQ[0].x, modulus);
    mpres_clear(EQ[0].y, modulus);
    mpres_clear(EQ[0].z, modulus);
    mpres_clear(EQ[0].A, modulus);
    mpmod_clear(modulus);
    mpz_clear(f);
    mpz_clear(d);
    mpz_clear(c);
    mpz_clear(b);
    mpz_clear(kx0);
    mpz_clear(ky0);
    return ret;
}

int
build_curves_with_torsion_Z10(mpz_t fac, mpz_t n, curve *tEP,
			      int umin, int umax, int nEP)
{
    int u, ret = ECM_NO_FACTOR_FOUND, nc = 0;
    mpz_t A2, A1div2, x0, y0, cte, num[2], den[2], inv[1], d, c, b, kx0, ky0;
    mpz_t f;
    mpmod_t modulus;
    curve EP[1], EQ[1]; /* blourk */
    char ok[1];

    ok[0] = 1;
    mpz_init(f);
    mpmod_init(modulus, n, ECM_MOD_DEFAULT);
    for(u = 0; u < 2; u++){
	mpres_init(num[u], modulus);
	mpres_init(den[u], modulus);
    }
    mpres_init(inv[0], modulus);
    /* Eaux = [2/3, -53/108] */
    /* Paux = [2/3, 1/2, 1] */
    mpres_init(EP[0].A, modulus);
    mod_from_rat_str(f, "2/3", n); mpres_set_z(EP[0].A, f, modulus);
    mpres_init(EP[0].x, modulus);
    mod_from_rat_str(f, "2/3", n); mpres_set_z(EP[0].x, f, modulus);
    mpres_init(EP[0].y, modulus);
    mod_from_rat_str(f, "1/2", n); mpres_set_z(EP[0].y, f, modulus);
    mpres_init(EP[0].z, modulus);
    mpres_set_ui(EP[0].z, 1, modulus);
#if DEBUG_MANY_EC >= 2
    printf("P:=");
    pt_print(EP[0], modulus);
    printf(";\n");
#endif

    mpres_init(EQ[0].x, modulus);
    mpres_init(EQ[0].y, modulus);
    mpres_init(EQ[0].z, modulus);
    mpres_init(EQ[0].A, modulus);
    mpres_set(EQ[0].A, EP[0].A, modulus);

    mpz_init(A2);
    mod_from_rat_str(A2, "2/3", n);
    mpz_init(A1div2);
    mod_from_rat_str(A1div2, "-1/2", n);
    mpz_mod(A1div2, A1div2, n);
    mpz_init_set_si(x0, 0);
    mpz_init_set_si(y0, 1);
    mpz_init_set_si(cte, -2);
    mpz_mod(cte, cte, n);

    mpz_init(d);
    mpz_init(c);
    mpz_init(b);
    mpz_init(kx0);
    mpz_init(ky0);
    for(u = umin; u < umax; u++){
	/* update Qaux */
	mpz_set_ui(d, u);
	if(pt_many_mul(EQ, EP, 1, d, modulus, num, den, inv, ok) == 0){
	    printf("found factor in Z10 (update of Q)\n");
	    mpz_set(fac, num[1]);
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
#if DEBUG_MANY_EC >= 2
	printf("(s, t)[%d]:=", u);
	pt_print(EQ[0], modulus);
	printf(";\n");
#endif
	mpres_get_z(b, EQ[0].x, modulus);
	mpres_get_z(c, EQ[0].y, modulus);
	if(cubic_to_quartic(fac, n, f, ky0, b, c, A2, A1div2, x0, y0, cte) == 0){
	    printf("found factor in Z10 (cubic_2_quartic)\n");
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
	/* f:=x; */
	/* d:=f^2/(f-(f-1)^2); */
	/** b <- f^2 **/
	mpz_mul(b, f, f);
	mpz_mod(b, b, n);
	mpz_sub_si(c, f, 1);
	mpz_mul(c, c, c);
	mpz_sub(c, f, c);
	mpz_mod(c, c, n);
	if(mod_from_rat2(d, b, c, n) == 0){
	    printf("inverse found in Z10 (d)\n");
	    mpz_set(fac, d);
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
	/* c:=f*(d-1); */
	mpz_sub_si(c, d, 1);
	mpz_mul(c, c, f);
	mpz_mod(c, c, n);
	/* ky0:=y*f^4/(f^2-3*f+1)^2/2; = num/den */
	/** b <- b^2 = f^4 **/
	mpz_mul(b, b, b);
	mpz_mod(b, b, n);
	mpz_mul(kx0, ky0, b);    /* num */
	mpz_sub_si(fac, f, 3);
	mpz_mul(fac, fac, f);
	mpz_add_si(fac, fac, 1);
	mpz_mul(fac, fac, fac);
	mpz_mul_si(fac, fac, 2); 
	mpz_mod(fac, fac, n);    /* den */
	if(mod_from_rat2(ky0, kx0, fac, n) == 0){
            printf("inverse found in Z10 (ky0)\n");
            ret = ECM_FACTOR_FOUND_STEP1;
            break;
        }
	/* kx0:=-f*d; */
	mpz_mul(kx0, f, d);
	mpz_mul_si(kx0, kx0, -1); /* humf */
	mpz_mod(kx0, kx0, n);
	/* b:=c*d; */
	mpz_mul(b, c, d);
	mpz_mod(b, b, n);
	K2W4(tEP+nc, b, c, kx0, ky0, n);
	nc++;
	if(nc >= nEP)
	    break;
	pt_many_assign(EP, EQ, 1, modulus);
    }
#if DEBUG_MANY_EC >= 2
    printf("Curves built\n");
    pt_many_print(tEP, nEP, modulus);
#endif
    mpz_clear(A2);
    mpz_clear(A1div2);
    mpz_clear(x0);
    mpz_clear(y0);
    mpz_clear(cte);
    mpres_clear(EP[0].A, modulus);
    mpres_clear(EP[0].x, modulus);
    mpres_clear(EP[0].y, modulus);
    mpres_clear(EP[0].z, modulus);
    for(u = 0; u < 2; u++){
	mpres_clear(num[u], modulus);
	mpres_clear(den[u], modulus);
    }
    mpres_clear(inv[0], modulus);
    mpres_clear(EP[0].x, modulus);
    mpres_clear(EP[0].y, modulus);
    mpres_clear(EP[0].z, modulus);
    mpres_clear(EP[0].A, modulus);
    mpres_clear(EQ[0].x, modulus);
    mpres_clear(EQ[0].y, modulus);
    mpres_clear(EQ[0].z, modulus);
    mpres_clear(EQ[0].A, modulus);
    mpmod_clear(modulus);
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
build_curves_with_torsion_Z2xZ8(mpz_t f, mpz_t n, curve *tEP,
				int umin, int umax, int nEP)
{
    int u, nc = 0, ret = ECM_NO_FACTOR_FOUND;
    mpz_t num[2], den[2], inv[1], tmp, a, b, alpha, beta, c, d, kx0, ky0, wx0;
    mpmod_t modulus;
    curve EP[1], EQ[1]; /* blourk */
    char ok[1];

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
    ok[0] = 1;
    mpmod_init(modulus, n, ECM_MOD_DEFAULT);
    for(u = 0; u < 2; u++){
	mpres_init(num[u], modulus);
	mpres_init(den[u], modulus);
    }
    mpres_init(inv[0], modulus);
    /* Eaux = [-8, -32] */
    /* Paux = [12, 40, 1] */
    mpres_init(EP[0].A, modulus);
    mpz_set_str(f, "-8", 10); mpres_set_z(EP[0].A, f, modulus);
    mpres_init(EP[0].x, modulus);
    mpz_set_str(f, "12", 10); mpres_set_z(EP[0].x, f, modulus);
    mpres_init(EP[0].y, modulus);
    mpz_set_str(f, "40", 10); mpres_set_z(EP[0].y, f, modulus);
    mpres_init(EP[0].z, modulus);
    mpres_set_ui(EP[0].z, 1, modulus);

    mpres_init(EQ[0].x, modulus);
    mpres_init(EQ[0].y, modulus);
    mpres_init(EQ[0].z, modulus);
    mpres_init(EQ[0].A, modulus);
    mpres_set(EQ[0].A, EP[0].A, modulus);

    for(u = umin; u < umax; u++){
	/* update Qaux */
	mpz_set_ui(d, u);
	if(pt_many_mul(EQ, EP, 1, d, modulus, num, den, inv, ok) == 0){
	    printf("found factor in Z10 (update of Q)\n");
	    mpz_set(f, num[1]);
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
#if DEBUG_MANY_EC >= 2
	printf("(s, t)[%d]:=", u);
	pt_print(EQ[0], modulus);
	printf(";\n");
#endif
	mpres_get_z(a, EQ[0].x, modulus);
	mpz_sub_si(a, a, 9);
	mpres_get_z(b, EQ[0].y, modulus);
	mpz_add_si(b, b, 25);
	if(mod_from_rat2(beta, b, a, n) == 0){
            printf("found factor in Z2xZ8 (beta)\n");
	    mpz_set(f, beta);
	    ret = ECM_FACTOR_FOUND_STEP1;
            break;
	}
	mpz_add_si(tmp, beta, 1);
	if(mpz_invert(alpha, tmp, n) == 0){
            printf("found factor in Z2xZ8 (alpha)\n");
	    mpz_gcd(f, tmp, n);
	    ret = ECM_FACTOR_FOUND_STEP1;
            break;
	}
	/** d <- 8*alpha^2-1; **/
	mpz_mul(d, alpha, alpha);
	mpz_mul_si(d, d, 8);
	mpz_sub_si(d, d, 1);
	mpz_mod(d, d, n);
	/* d:=2*alpha*(4*alpha+1)/d; */
	mpz_mul_si(c, alpha, 4);
	mpz_add_si(c, c, 1);
	mpz_mul(c, c, alpha);
	mpz_mul_si(c, c, 2);
	mpz_mod(c, c, n);
	if(mod_from_rat2(f, c, d, n) == 0){
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
	mpz_mod(f, f, n);
        if(mod_from_rat2(c, f, d, n) == 0){
            printf("found factor in Z2xZ8 (d)\n");
            ret = ECM_FACTOR_FOUND_STEP1;
            break;
        }
	/* b = c*d */
	mpz_mul(b, c, d);
	mpz_mod(b, b, n);
	/* kx0:=-(2*d-1)/4;*/
	mod_div_2(kx0, n);
	mod_div_2(kx0, n);
	mpz_mul_si(kx0, kx0, -1);
	mpz_mod(kx0, kx0, n);
	/* ky0:=(c/8)*(-beta^2+2*uP[1]+9); */
	mpz_mul(f, beta, beta);
	mpres_get_z(a, EQ[0].x, modulus);
	mpz_sub(f, a, f);
	mpz_add(f, f, a);
	mpz_add_si(f, f, 9);
	mpz_mul(f, f, c);
	mpz_mod(f, f, n);
	mod_div_2(f, n);
	mod_div_2(f, n);
	mod_div_2(f, n);
	/* ky0:=ky0/(beta^2+2*beta-7); */
	mpz_add_si(tmp, beta, 2);
	mpz_mul(tmp, tmp, beta);
	mpz_sub_si(tmp, tmp, 7);
	mpz_mod(tmp, tmp, n);
	if(mod_from_rat2(ky0, f, tmp, n) == 0){
            printf("found factor in Z2xZ8 (ky0)\n");
	    mpz_set(f, ky0);
            ret = ECM_FACTOR_FOUND_STEP1;
            break;
        }
	K2W24(f, a, b, c, n);
	/* wx0:=kx0+a2/3; */
        mpz_set_si(tmp, 3);
	mod_from_rat2(wx0, f, tmp, n);
	mpz_add(wx0, wx0, kx0);
	mpz_mod(wx0, wx0, n);
	/* ma:=-1/4*(8*d^4-16*d^3+16*d^2-8*d+1)/(d-1)^2/d^2; */
	mpz_sub_si(tmp, d, 1);    /* num */
	mpz_mul(tmp, tmp, d);
	mpz_mul(tmp, tmp, tmp);
	mpz_mul_si(tmp, tmp, -4);
	mpz_mod(tmp, tmp, n);
	mpz_set_si(f, 8);         /* den */
	mpz_mul(f, f, d); mpz_add_si(f, f, -16);
	mpz_mul(f, f, d); mpz_add_si(f, f, 16);
	mpz_mul(f, f, d); mpz_add_si(f, f, -8);
	mpz_mul(f, f, d); mpz_add_si(f, f, 1);
	if(mod_from_rat2(tEP[nc].A, f, tmp, n) == 0){
            printf("found factor in Z2xZ8 (ma)\n");
	    mpz_set(f, tEP[nc].A);
            ret = ECM_FACTOR_FOUND_STEP1;
            break;
        }
	/* mb:=-1/(d-1)^2; */
	mpz_sub_si(tmp, d, 1);
	mpz_mul(tmp, tmp, tmp);
	mpz_mod(tmp, tmp, n);
	if(mpz_invert(f, tmp, n) == 0){
	    printf("found factor in Z2xZ8 (mb)\n");
	    mpz_gcd(f, tmp, n);
            ret = ECM_FACTOR_FOUND_STEP1;
            break;
	}
	mpz_set_si(tmp, 0);
	mpz_sub(tmp, tmp, f);
	mpz_mod(tmp, tmp, n);
	/* mx:=mb*wx0-ma/3; */
	mpz_mul(f, tmp, wx0);
        mpz_set_si(tmp, 3);
        mod_from_rat2(tEP[nc].x, tEP[nc].A, tmp, n);
	mpz_sub(tEP[nc].x, f, tEP[nc].x);
	mpz_mod(tEP[nc].x, tEP[nc].x, n);
	/* my:=mb*ky0; */
#if DEBUG_MANY_EC >= 2
	gmp_printf("N:=%Zd;\n", n);
	gmp_printf("ma:=%Zd;\n", tEP[nc].A);
	gmp_printf("kx0:=%Zd;\n", kx0);
	gmp_printf("ky0:=%Zd;\n", ky0);
	gmp_printf("mx0:=%Zd;\n", tEP[nc].x);
#endif
	nc++;
	if(nc >= nEP)
	    break;
	pt_many_assign(EP, EQ, 1, modulus);
    }
#if DEBUG_MANY_EC >= 2
    printf("Curves built\n");
    pt_many_print(tEP, nEP, modulus);
#endif
    mpres_clear(EP[0].A, modulus);
    mpres_clear(EP[0].x, modulus);
    mpres_clear(EP[0].y, modulus);
    mpres_clear(EP[0].z, modulus);
    for(u = 0; u < 2; u++){
	mpres_clear(num[u], modulus);
	mpres_clear(den[u], modulus);
    }
    mpres_clear(inv[0], modulus);
    mpres_clear(EP[0].x, modulus);
    mpres_clear(EP[0].y, modulus);
    mpres_clear(EP[0].z, modulus);
    mpres_clear(EP[0].A, modulus);
    mpres_clear(EQ[0].x, modulus);
    mpres_clear(EQ[0].y, modulus);
    mpres_clear(EQ[0].z, modulus);
    mpres_clear(EQ[0].A, modulus);
    mpmod_clear(modulus);
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
    return ret;
}

int
build_curves_with_torsion_Z3xZ3(mpz_t n, curve *tEP,
				int smin, int smax, int nEP)
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
	mpz_mod(x0, x0, n);
	/* a2:=T^6+108;*/
	mpz_add_ui(a2, x0, 108);
	mpz_mod(a2, a2, n);
	/* a4:=144*T^6+3888; */
	mpz_mul_ui(a4, x0, 144);
	mpz_add_ui(a4, a4, 3888);
	mpz_mod(a4, a4, n);
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
	mpz_mod(y0, y0, n);
	mpz_set_ui(x0, 0);
	W2W(tEP+nc, a2, a4, x0, y0, n);
	nc++;
	if(nc >= nEP)
	    break;
    }
    mpz_clear(x0);
    mpz_clear(y0);
    mpz_clear(a2);
    mpz_clear(a4);
    return ECM_NO_FACTOR_FOUND;
}

/* Original source is Brier + Clavier.
   We can build curves in Montgomery form directly... 
   Useful if one knows that all p | n are 1 mod 4 (Cunningham, etc.).
*/
int
build_curves_with_torsion_Z4xZ4(mpz_t f, mpz_t n, curve *tEP,
				int smin, int smax, int nEP)
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
	mpz_set_ui(nu2, nu*nu);
	/* tau:=(nu^2+3)/2/nu; */
	mpz_add_si(lambda, nu2, 3);
	mpz_set_si(tmp, 2*nu);
	mod_from_rat2(tau, lambda, tmp, n);
	/* lambda:=8*nu^3; */
	mpz_mul_si(lambda, nu2, 8*nu);
	mpz_mod(lambda, lambda, n);
	/* A:=-27*lambda^4*(tau^8+14*tau^4+1); */
	/* B:=54*lambda^6*(tau^12-33*tau^8-33*tau^4+1); */
	/* x0:=3*(3*nu^12+34*nu^10+117*nu^8+316*nu^6+1053*nu^4+2754*nu^2+2187); */
	/* y0:=27*(nu^2-3)*(nu^2+1)*(nu^2+9)*(nu^6+5*nu^4+15*nu^2+27)^2; */
	/* P = (x0, y0) is a point on Y^2 = X^3+A*X+B */

	/* Montgomery form: there are several b possible */
	/* b:=1/9/lambda^2/(tau^4-1); */
	mpz_powm_ui(x0, tau, 4, n);
	mpz_sub_si(x0, x0, 1);
	mpz_mod(x0, x0, n);
	mpz_mul(tmp, x0, lambda);
	mpz_mul(tmp, tmp, lambda);
	mpz_mul_si(tmp, tmp, 9);
	if(mpz_invert(b, tmp, n) == 0){
	    printf("Factor found durint init of Z4xZ4\n");
	    mpz_gcd(f, tmp, n);
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
	/* a:=-2*(tau^4+1)/(tau^4-1); */
	mpz_add_si(tmp, x0, 2);
	mpz_mul_si(tmp, tmp, -2);
	mpz_mod(tmp, tmp, n);
	mod_from_rat2(tEP[nc].A, tmp, x0, n);
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
	mpz_mod(x0, x0, n);
#if DEBUG_MANY_EC >= 2
	gmp_printf("N:=%Zd;\n", n);
	printf("nu:=%d;\n", nu);
	gmp_printf("tau:=%Zd;\n", tau);
	gmp_printf("lambda:=%Zd;\n", lambda);
	gmp_printf("a:=%Zd;\n", tEP[nc].A);
	gmp_printf("x0:=%Zd;\n", x0);
#endif
	/* x:=b*x0-a/3; not needed: y:=b*y0 */
	mpz_set_si(tmp, 3);
	mod_from_rat2(tEP[nc].x, tEP[nc].A, tmp, n);
	mpz_mul(b, b, x0);
	mpz_mod(b, b, n);
	mpz_sub(tEP[nc].x, b, tEP[nc].x);
	mpz_mod(tEP[nc].x, tEP[nc].x, n);
	nc++;
	if(nc >= nEP)
	    break;
    }
    mpz_clear(tau);
    mpz_clear(lambda);
    mpz_clear(nu2);
    mpz_clear(tmp);
    mpz_clear(b);
    mpz_clear(x0);
    return ret;
}

/* Assuming we can generate curves with given torsion using parameter s
   in interval [smin..smax].
*/
int
build_curves_with_torsion(mpz_t f, mpz_t n, curve *tEP, char *torsion, 
			  int smin, int smax, int nEP)
{
    int ret = 0;

    /* over Q: see Atkin-Morain, Math. Comp., 1993 */
    if(strcmp(torsion, "Z5") == 0)
	return build_curves_with_torsion_Z5(f, n, tEP, smin, smax, nEP);
    else if(strcmp(torsion, "Z7") == 0)
	return build_curves_with_torsion_Z7(f, n, tEP, smin, smax, nEP);
    else if(strcmp(torsion, "Z9") == 0)
	return build_curves_with_torsion_Z9(f, n, tEP, smin, smax, nEP);
    else if(strcmp(torsion, "Z10") == 0)
	return build_curves_with_torsion_Z10(f, n, tEP, smin, smax, nEP);
    else if(strcmp(torsion, "Z2xZ8") == 0)
	return build_curves_with_torsion_Z2xZ8(f, n, tEP, smin, smax, nEP);
    /* no longer over Q */
    else if(strcmp(torsion, "Z3xZ3") == 0) /* over Q(sqrt(-3)) */
	return build_curves_with_torsion_Z3xZ3(n, tEP, smin, smax, nEP);
    else if(strcmp(torsion, "Z4xZ4") == 0) /* over Q(sqrt(-1)) */
	return build_curves_with_torsion_Z4xZ4(f, n, tEP, smin, smax, nEP);
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
    int ret = 0, i, mtyform;

    mtyform = strcmp(torsion, "Z4xZ4") == 0 || strcmp(torsion, "Z2xZ8") == 0;
    for(i = 0; i < nEP; i++){
	mpz_init(tEP[i].A);
	mpz_init(tEP[i].x);
	if(mtyform == 0){
	    mpz_init(tEP[i].y);
	    mpz_init(tEP[i].z);
	}
    }
    if(build_curves_with_torsion(f, n, tEP, torsion, smin, smax, nEP) 
       == ECM_NO_FACTOR_FOUND){
	if(mtyform)
	    ret = process_many_curves(f, n, B1, tEP, nEP, 1);
	else
	    ret = process_many_curves(f, n, B1, tEP, nEP, 0);
    }
    for(i = 0; i < nEP; i++){
	mpz_clear(tEP[i].A);
	mpz_clear(tEP[i].x);
	if(mtyform == 0){
	    mpz_clear(tEP[i].y);
	    mpz_clear(tEP[i].z);
	}
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
  int res = 0, smin = -1, smax = -1, ncurves = 0, method = ECM_ECM, onebyone;
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
      else if (strcmp (argv[1], "-pm1") == 0)
	{
	  method = ECM_PM1;
	  argv++;
	  argc--;
	}
      else if (strcmp (argv[1], "-pp1") == 0)
	{
	  method = ECM_PP1;
	  argv++;
	  argc--;
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
  setlinebuf(stdout);
  while(fscanf(infile, "%s", buf) != EOF){
      /* read number */
      if(buf[0] == '#'){
	  char c;
	  /* print till end of line */
	  printf("%s", buf);
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
	  onebyone = 0;
	  res = process_many_curves_from_file(f, n, B1, curvesname, ncurves,
					      onebyone);
      }
      else if(torsion != NULL){
	  res = process_curves_with_torsion(f, n, B1, torsion, smin, smax, ncurves);
      }
  }
  if(infile != stdin)
      fclose(infile);

  mpz_clear (f);
  mpz_clear (n);

  return res;
}
