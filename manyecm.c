/* manyecm.c - ECM with many curves with many torsion and/or in parallel 
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

/********** group law on points **********/

#define pt_is_equal(P, Q) (mpz_cmp((P)->x, (Q)->x) == 0 \
	                     && mpz_cmp((P)->y, (Q)->y) == 0 \
			     && mpz_cmp((P)->z, (Q)->z) == 0)

int
pt_is_zero(ec_point_t P, ATTRIBUTE_UNUSED mpmod_t n)
{
    return mpz_sgn(P->z) == 0;
}

void
pt_set_to_zero(ec_point_t P, mpmod_t n)
{
    mpz_set_ui(P->x, 0);
    mpres_set_ui(P->y, 1, n);
    mpz_set_ui(P->z, 0);
}

void
pt_assign(ec_point_t Q, ec_point_t P, ATTRIBUTE_UNUSED mpmod_t n)
{
    mpres_set(Q->x, P->x, n);
    mpres_set(Q->y, P->y, n);
    mpres_set(Q->z, P->z, n);
}

void
pt_neg(ec_point_t P, mpmod_t n)
{
    if(pt_is_zero(P, n) == 0)
	mpres_neg(P->y, P->y, n);
}

void
pt_many_set_to_zero(ec_point_t *tP, int nE, mpmod_t n)
{
    int i;

    for(i = 0; i < nE; i++)
	pt_set_to_zero(tP[i], n);
}

void
pt_many_neg(ec_point_t *tP, int nE, mpmod_t n)
{
    int i;

    for(i = 0; i < nE; i++)
	pt_neg(tP[i], n);
}

void
pt_many_assign(ec_point_t *tQ, ec_point_t *tP, int nE, mpmod_t n)
{
    int i;

    for(i = 0; i < nE; i++)
	pt_assign(tQ[i], tP[i], n);
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
pt_print(ec_point_t P, mpmod_t n)
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
pt_many_print(ec_curve_t *tE, ec_point_t *tP, int nE, mpmod_t n)
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

/* NOTE: we can have tR = tP or tQ.
   In case a factor is found, it is put in num[nE].
 */
int
pt_many_common(ec_point_t *tR, ec_point_t *tP, ec_point_t *tQ, int nE, 
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
pt_many_duplicate(ec_point_t *tQ, ec_point_t *tP, ec_curve_t *tE, int nE, 
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
pt_many_add(ec_point_t *tR, ec_point_t *tP, ec_point_t *tQ, ec_curve_t *tE, 
	    int nE, mpmod_t n, 
	    mpres_t *num, mpres_t *den, mpres_t *inv, char *ok)
{
    char *takeit = (char *)malloc(nE * sizeof(char));
    int i, res;

    memcpy(takeit, ok, nE);
#if DEBUG_MANY_EC >= 2
    printf("In pt_many_add, adding\n");
    pt_many_print(tP, nE, n);
    printf("and\n");
    pt_many_print(tQ, nE, n);
#endif
    for(i = 0; i < nE; i++){
	if(ok[i] == 0)
	    continue; /* takeit[i] = 0 */
	if(pt_is_zero(tP[i], n)){
#if DEBUG_MANY_EC >= 2
	    printf("# tEP[%d] = O_{E[%d]}\n", i, i);
#endif
	    takeit[i] = 0;
	    pt_assign(tR[i], tQ[i], n);
	}
	else if(pt_is_zero(tQ[i], n)){
#if DEBUG_MANY_EC >= 2
	    printf("# tEQ[%d] = O_{E[%d]}\n", i, i);
#endif
	    takeit[i] = 0;
	    pt_assign(tR[i], tP[i], n);
	}
	else if(pt_is_equal(tP[i], tQ[i])){
	    /* we should double */
	    if(mpz_sgn(tP[i]->y) == 0){
#if DEBUG_MANY_EC >= 2
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
pt_many_sub(ec_point_t *tR, ec_point_t *tQ, ec_point_t *tP, ec_curve_t *tE,
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
pt_many_mul_plain(ec_point_t *tQ, ec_point_t *tP, ec_curve_t *tE,
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
#if DEBUG_MANY_EC >= 2
	printf("Rdup:="); pt_many_print(tQ, nE, n); printf(";\n");
#endif
	if (mpz_tstbit (e, l))
	  {
	      if(pt_many_add (tQ, tP, tQ, tE, nE, n, num, den, inv, ok) == 0)
	      {
		status = 0;
		break;
	      }
#if DEBUG_MANY_EC >= 2
	      printf("Radd:="); pt_many_print(tQ, nE, n); printf(";\n");
#endif
	  }
    }
  return status;
}

/* Ordinary binary left-right addition; see Solinas00. Morally, we use
 w = 2. */
static int
pt_many_mul_add_sub_si(ec_point_t *tQ, ec_point_t *tP, ec_curve_t *tE, int nE,
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
#if DEBUG_MANY_EC >= 2
	printf("Rdup:="); pt_many_print(tQ, nE, n); printf(";\n");
#endif
	if(S[j] == 1){
	    if(pt_many_add(tQ, tQ, tP, tE, nE, n, num, den, inv, ok) == 0){
		status = 0;
		break;
	    }
#if DEBUG_MANY_EC >= 2
	    printf("Radd:="); pt_many_print(tQ, nE, n); printf(";\n");
#endif
	}
	else if(S[j] == -1){
	    if(pt_many_sub(tQ, tQ, tP, tE, nE, n, num, den, inv, ok) == 0){
		status = 0;
		break;
	    }
#if DEBUG_MANY_EC >= 2
	    printf("Rsub:="); pt_many_print(tQ, nE, n); printf(";\n");
#endif
	}
    }
  return status;
}

/* tEQ[i] <- e * tEP[i]; we must have tEQ != tEP */
/* If a factor is found, it is put back in num[nE]. */
int
pt_many_mul(ec_point_t *tQ, ec_point_t *tP, ec_curve_t *tE, int nE,
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
   TODO: use chkfile also.
 */
int
process_one_curve(mpz_t f, mpz_t N, double B1,
		  ecm_params params, ec_curve_t E, ec_point_t P)
{
    int ret;
    double B2scale = 1.0;

    /* Taken from main.c; no comment */
    /* Here's an ugly hack to pass B2scale to the library somehow.
       It gets piggy-backed onto B1done */
    params->B1done = params->B1done + floor (B2scale * 128.) / 134217728.; 

#if 0
    /* compute it automatically from B1: no freedom on B2! */
    mpz_set_si(params->B2, ECM_DEFAULT_B2); 
#endif
    /* will be set to B1 */
    mpz_set_si(params->B2min, ECM_DEFAULT_B2);

    mpz_set(params->x, P->x);
    mpz_set(params->sigma, E->A); /* humf */

    if(E->type == ECM_EC_TYPE_MONTGOMERY)
	params->sigma_is_A = 1;
    else{
	params->sigma_is_A = -1;
	params->Etype = E->type;
	mpz_set(params->y, P->y);
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
	printf("# found input number\n");
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

/* f is a (probable) prime factor of n. tP is in plain mod n form. */
void
dump_curves(ec_curve_t *tE, ec_point_t *tP, int nE, mpz_t f)
{
    int i;

    gmp_printf("p:=%Zd; F:=GF(p); P:=[]; A:=[]; B:=[]; E:=[];\n", f);
    printf("CheckE:=procedure(E, info)\n");
    printf("    printf \"#E[%%o]=%%o\\n\", info, Factorization(#E);\n");
    printf("    gen:=Generators(E);\n");
    printf("    printf \"ords=%%o\\n\", ");
    printf("[Factorization(Order(g)) : g in gen];\n");
    printf("end procedure;\n");
    for(i = 0; i < nE; i++){
	if(tE[i]->type == ECM_EC_TYPE_MONTGOMERY){
	    mpmod_t fmod;
	    mpres_t x, y, A;
	    mpz_t tmp;

	    mpz_init(tmp);
	    mpmod_init(fmod, f, ECM_MOD_DEFAULT);
	    mpres_init(x, fmod);
	    mpres_set_z(x, tP[i]->x, fmod);
	    mpres_init(y, fmod);
	    mpres_init(A, fmod);
	    mpres_set_z(A, tE[i]->A, fmod);
	    if(montgomery_to_weierstrass(tmp, x, y, A, fmod) == ECM_FACTOR_FOUND_STEP1){
		printf("GASP in dump!\n");
	    }
	    printf("P[%d]:=[", i+1); print_mpz_from_mpres(x, fmod);
	    printf(", "); print_mpz_from_mpres(y,fmod);
	    printf(", 1];\n");
	    printf("A[%d]:=", i+1);
	    print_mpz_from_mpres(A, fmod);
	    printf(";\n");
	    mpres_clear(x, fmod);
	    mpres_clear(y, fmod);
	    mpres_clear(A, fmod);
	    mpmod_clear(fmod);
	    mpz_clear(tmp);
	}
	else if(tE[i]->type == ECM_EC_TYPE_WEIERSTRASS){
	    gmp_printf("P[%d]:=[%Zd, %Zd, %Zd];\n", i+1, 
		       tP[i]->x, tP[i]->y, tP[i]->z); 
	    gmp_printf("A[%d]:=%Zd;\n", i+1, tE[i]->A);
	}
	else{
	    printf("Case %d NYI in dump_curves\n", tE[i]->type);
	    break;
	}
	printf("B[%d]:=P[%d][2]^2-P[%d][1]^3-A[%d]*P[%d][1];\n", 
	       i+1, i+1, i+1, i+1, i+1);
	printf("E[%d]:=EllipticCurve([F!A[%d], F!B[%d]]);\n", i+1, i+1, i+1);
	printf("CheckE(E[%d], infos[%d]);\n", i+1, i+1);
    }
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

/* TODO: better control of B2 + dichotomy (cf. #B2) */
int
one_curve_at_a_time(mpz_t f, char *ok, ec_curve_t *tE, ec_point_t *tP, int nE,
		    mpz_t N, ecm_params params, double B1, 
		    char *savefilename)
{
    double tmpB1, tmpB2, B2g = 0, B2d = 0, B2 = 0; /* 1e9; #B2 */
    int ret = 0, i, saveit, nhit, nhitmax = 1; /* #B2 */
    mpcandi_t candi;
    char comment[256] = "";
    mpz_t C;

    mpcandi_t_init(&candi);
    mpcandi_t_add_candidate(&candi, N, NULL, 0);
    mpz_init(C);
    /* process curves one at a time */
    for(i = 0; i < nE; i++){
	tmpB1 = B1;
	tmpB2 = B2;
	nhit = 0;
	while(1){
#if DEBUG_MANY_EC >= 2
	    printf("infos:=[\"E%d\"];\n", i);
	    dump_curves(tE+i, tP+i, 1, N);
#endif
	    params->B1done = 1.0;
#if 0 /* #B2 */
	    mpz_set_d(params->B2, tmpB2);
#endif
	    if(nhit > 0){
		tmpB2 = (B2d+B2g)/2;
		printf("# trying new B2[%d]=%f\n", nhit, tmpB2);
	    }
	    ret = process_one_curve(f,N,tmpB1,params,tE[i],tP[i]);
	    if(ret == ECM_NO_FACTOR_FOUND){
		if(nhit == 0)
		    /* no factor found in any step */
		    break;
		else{
		    /* we are in some recursive step */
		    printf("# B1done=%.0f\n", params->B1done);
		    if(params->B1done == tmpB1)
			/* dichotomy for step 2 */
			B2g = tmpB2;
		}
	    }
	    else if(ret == ECM_FACTOR_FOUND_STEP1){
		if(mpz_cmp(f, N) != 0)
		    /* non-trivial factor found */
		    break;
		else{
		    tmpB1 = params->B1done - 1;
		    printf("# trying again with B1=%.0f\n", tmpB1);
		}
	    }
	    else if(ret == ECM_FACTOR_FOUND_STEP2){
		if(mpz_cmp(f, N) != 0)
		    /* non-trivial factor found */
		    break;
		else{
		    if(nhit == 0)
			B2g = 0;
		    B2d = tmpB2;
		}
	    }
	    else
		break;
	    nhit++;
	    if(nhit == nhitmax) /* caution, Lemmy! */
		break;
	}
	saveit = (savefilename != NULL);
	if(ret > 0){ /* humf */
	    ok[i] = 0;
	    ret = conclude_on_factor(N, f, params->verbose);
	    if(ret == ECM_INPUT_NUMBER_FOUND){
		printf("# B1done=%.0f\n", params->B1done);
		printf("# proceeding to next curve\n");
		saveit = 0;
	    }
	    else{
#if DEBUG_MANY_EC >= 2
		if(ret == ECM_PRIME_FAC_PRIME_COFAC 
		   || ret == ECM_PRIME_FAC_COMP_COFAC)
		    /* output Magma lines to check #E's mod f */
		    dump_curves(tE, tP, nE, f);
#endif
		break;
	    }
	}
	else if(ret == ECM_ERROR){
	    printf("Error for curve %d\n", i);
	}
	if(saveit){
	    write_resumefile(savefilename, ECM_ECM, N, params, &candi, 
			     tP[i]->x, tP[i]->y, comment);
	}

    }
    mpz_clear (C);
    mpcandi_t_free(&candi);
    return ret;
}

/* Using parallelism.
   Copied from classical ecm_stage1.
 */
int
all_curves_at_once(mpz_t f, char *ok, ec_curve_t *tE, ec_point_t *tP, int nE,
		   mpmod_t n, double B1, double *B1done, 
		   int (*stop_asap)(void), char *chkfilename)
{
    ec_point_t tQ[NCURVE_MAX], tR[NCURVE_MAX];
    mpz_t num[NCURVE_MAX+1], den[NCURVE_MAX+1], inv[NCURVE_MAX], e;
    double p = 0.0, r, last_chkpnt_p;
    int ret = ECM_NO_FACTOR_FOUND;
    long last_chkpnt_time;
    int i;
    
    mpz_init(e);
    for(i = 0; i < nE; i++){
	mpres_init(tQ[i]->x, n); mpres_set(tQ[i]->x, tP[i]->x, n);
	mpres_init(tQ[i]->y, n); mpres_set(tQ[i]->y, tP[i]->y, n);
	mpres_init(tQ[i]->z, n); mpres_set(tQ[i]->z, tP[i]->z, n);
	
	mpres_init(tR[i]->x, n);
	mpres_init(tR[i]->y, n);
	mpres_init(tR[i]->z, n);
	
	mpres_init(num[i], n);
	mpres_init(den[i], n);
	mpres_init(inv[i], n);
    }
    mpres_init(num[nE], n); /* to be used as buffer in compute_all_inverses */
    mpres_init(den[nE], n); /* to be used as buffer in compute_all_inverses */
    
    last_chkpnt_time = cputime ();
    
#if DEBUG_MANY_EC >= 2
    printf("Initial points:\n");
    pt_many_print(tP, nE, n);
#endif
    for (r = 2.0; r <= B1; r *= 2.0)
	if (r > *B1done){
	    if(pt_many_duplicate (tQ, tQ, tE, nE, n, num, den, inv, ok) == 0){
		mpz_set(f, num[nE]);
		ret = ECM_FACTOR_FOUND_STEP1;
		goto end_of_all;
	    }
#if DEBUG_MANY_EC >= 2
	    printf("P%ld:=", (long)r); pt_many_print(tQ, nE, n); printf(";\n");
#endif
	}

    last_chkpnt_p = 3.;
    for (p = getprime (); p <= B1; p = getprime ()){
	for (r = p; r <= B1; r *= p){
#if DEBUG_MANY_EC >= 2
	    printf("## p = %ld at %ldms\n", (long)p, cputime());
#endif
	    if (r > *B1done){
		mpz_set_ui(e, (ecm_uint) p);
		if(pt_many_mul(tR, tQ, tE, nE, e, n, num, den, inv, ok) == 0){
		    mpz_set(f, num[nE]);
		    ret = ECM_FACTOR_FOUND_STEP1;
		    goto end_of_all;
		}
#if DEBUG_MANY_EC >= 2
		pt_many_print(tR, nE, n);
#endif
		for(i = 0; i < nE; i++)
		    if(pt_is_zero(tR[i], n))
			ok[i] = 0;
		pt_many_assign(tQ, tR, nE, n); /* TODO: use pointers */
	    }
	    if (stop_asap != NULL && (*stop_asap) ()){
		outputf (OUTPUT_NORMAL, "Interrupted at prime %.0f\n", p);
		break;
	    }
	    
	    /* WARNING: not activated yet */
	    if (chkfilename != NULL && p > last_chkpnt_p + 10000. && 
		elltime (last_chkpnt_time, cputime ()) > CHKPNT_PERIOD){
#if 0 /* TODO: make this work for many curves */
		writechkfile (chkfilename, ECM_ECM, MAX(p, *B1done), n, A, x, y, z);
#endif
		last_chkpnt_p = p;
		last_chkpnt_time = cputime ();
	    }
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
    pt_many_assign(tP, tQ, nE, n);
    /* normalize all points */
    for(i = 0; i < nE; i++)
	if(pt_is_zero(tP[i], n))
	    pt_set_to_zero(tP[i], n);
    /* clear temporary variables */
    mpz_clear(e);
    for(i = 0; i < nE; i++){
	mpres_clear(tQ[i]->x, n);
	mpres_clear(tQ[i]->y, n);
	mpres_clear(tQ[i]->z, n);
	mpres_clear(num[i], n);
	mpres_clear(den[i], n);
	mpres_clear(inv[i], n);
    }
    mpres_clear(num[nE], n);
    mpres_clear(den[nE], n);
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

/* 
   OUTPUT: ECM_NO_FACTOR_FOUND
           ECM_INPUT_NUMBER_FOUND
           ECM_PRIME_FAC_PRIME_COFAC
	   ECM_PRIME_FAC_COMP_COFAC
	   ECM_COMP_FAC_COMP_COFAC
	   ECM_COMP_FAC_PRIME_COFAC
*/
int
process_many_curves(mpz_t f, mpmod_t n, double B1, 
		    ec_curve_t *tE, ec_point_t *tP, int nE, 
		    ecm_params params, int onebyone, char *savefilename)
{
    double B1done;
    ec_point_t tQ[NCURVE_MAX];
    char *ok = (char *)malloc(nE * sizeof(char));
    int ret = 0, i;
    long st = cputime ();
    
    memset(ok, 1, nE);
    if(onebyone != 0){
	ret = one_curve_at_a_time(f, ok, tE, tP, nE, n->orig_modulus, params,
				  B1, savefilename);
	free(ok);
	return ret;
    }
    /* take everybody */
    for(i = 0; i < nE; i++){
	ec_point_init(tQ[i], tE[i], n);
	ec_point_set(tQ[i], tP[i], tE[i], n);
    }
    B1done = 1.0;
    ret = all_curves_at_once(f, ok, tE, tQ, nE, n, B1, &B1done, NULL, NULL);
    printf("# Step 1 took %ldms\n", elltime (st, cputime ()));

    if(ret != ECM_NO_FACTOR_FOUND){
	ret = conclude_on_factor(n->orig_modulus, f, params->verbose);
#if DEBUG_MANY_EC >= 2
	if(ret == ECM_PRIME_FAC_PRIME_COFAC || ret == ECM_PRIME_FAC_COMP_COFAC)
	    /* output Magma lines to check properties of E mod f */
	    dump_curves(tE, tP, nE, f);
#endif
    }
    else{
	params->sigma_is_A = -1;
	params->B1done = B1;
	for(i = 0; i < nE; i++){
	    if(ok[i] == 0)
		continue;
#if DEBUG_MANY_EC >= 1
	    printf("# Entering Step 2 for E[%d]\n", i);
#endif
	    mpres_get_z(tP[i]->x, tQ[i]->x, n);
	    mpres_get_z(tP[i]->y, tQ[i]->y, n);
	    mpres_get_z(tP[i]->z, tQ[i]->z, n);
	    ret = process_one_curve(f, n->orig_modulus, B1, params,
				    tE[i], tP[i]);
	    if(ret != ECM_NO_FACTOR_FOUND){
		printf("## factor found in Step 2: ");
		mpz_out_str (stdout, 10, f);
		printf ("\n");
		ret = conclude_on_factor(n->orig_modulus, f, params->verbose);
		break;
	    }
	}
    }
    for(i = 0; i < nE; i++)
	ec_point_clear(tQ[i], tE[i], n);
    free(ok);
    return ret;
}

int
read_curves_from_file(int *nE, ec_curve_t *tE, ec_point_t *tP, 
		      mpz_t *tf, int *nf,
		      mpmod_t n, char *fic_EP, int ncurves)
{
    FILE *ifile = fopen(fic_EP, "r");
    char bufA[1024], bufx[1024], bufy[1024], c, Etype;
    mpq_t q;
    int ret = ECM_NO_FACTOR_FOUND;

    *nE = 0;
    mpq_init(q);
    while(fscanf(ifile, "%s", bufA) != EOF){
	if(bufA[0] == '#'){
	    /* skip line and print it */
	    printf("%s", bufA);
	    while((c = getc(ifile)) != '\n')
		printf("%c", c);
	    printf("\n");
	    continue;
	}
	else
	    Etype = bufA[0];
	ec_curve_init(tE[*nE], n);
	if(Etype == 'W'){
	    if(fscanf(ifile, "%s %s %s", bufA, bufx, bufy) == EOF)
		break;
	    tE[*nE]->type = ECM_EC_TYPE_WEIERSTRASS;
	}
	else if(Etype == 'H'){
	    if(fscanf(ifile, "%s %s %s", bufA, bufx, bufy) == EOF)
		break;
	    tE[*nE]->type = ECM_EC_TYPE_HESSIAN;
	}
	else if(Etype == 'M'){
	    if(fscanf(ifile, "%s %s", bufA, bufx) == EOF)
		break;
	    tE[*nE]->type = ECM_EC_TYPE_MONTGOMERY;
	}
	else{
	    printf("Unknown curve type: %c\n", Etype);
	    return ECM_ERROR;
	}
	mpz_init(tE[*nE]->A);
	if(read_and_prepare(tf[*nf], tE[*nE]->A, q, bufA, n->orig_modulus) == 0){
	    ret = 0;
	    *nf += 1;
	    goto process_end;
	}
	ec_point_init(tP[*nE], tE[*nE], n);
	mpz_init(tP[*nE]->x);
	if(read_and_prepare(tf[*nf], tP[*nE]->x, q, bufx, n->orig_modulus) == 0){
	    ret = 0;
            *nf+= 1;
	    goto process_end;
	}
	mpz_init(tP[*nE]->y);
	if((Etype == 'W') || (Etype == 'H')){
	    if(read_and_prepare(tf[*nf], tP[*nE]->y, q, bufy, n->orig_modulus) == 0){
		ret = 0;
		*nf+= 1;
		goto process_end;
	    }
	}
	mpz_init_set_ui(tP[*nE]->z, 1);
	*nE += 1;
	if(ncurves != 0 && *nE == ncurves)
	    break;
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

/* Weierstrass (a2, a4, a6) to (A, B)
   A = (a4-1/3*a2^2)
   B = -1/3*a4*a2 + 2/27*a2^3 + a6
     = -1/3*a2*(a4-2/9*a2^2) + a6
   X = x+a2/3
   INPUT: if x0 == NULL, we have no point to translate
          if B == NULL, we do not need and we do not compute B
*/
void
W2W(ec_curve_t E, mpz_t B, ec_point_t P, mpz_t a2, mpz_t a4, mpz_t a6, 
    mpz_t x0, mpz_t y0, mpz_t n)
{
    mpz_t tmp1, tmp2, tmp3;

    mpz_init(tmp1);
    mpz_init(tmp2);
    /* tmp2 <- a2/3 */
    mpz_init_set_si(tmp3, 3);
    mod_from_rat2(tmp2, a2, tmp3, n);
    if(x0 != NULL){
	/* wx0 = x0 + a2/3 */
	mpz_add(P->x, tmp2, x0);
	mpz_mod(P->x, P->x, n);
	mpz_set(P->y, y0);
	mpz_mod(P->y, P->y, n);
	mpz_set_ui(P->z, 1);
    }
    /* A = a4-1/3*a2^2 = a4 - a2 * (a2/3) */
    /** tmp1 <- tmp2*a2 = a2^2/3 */
    mpz_mul(tmp1, a2, tmp2);
    mpz_mod(tmp1, tmp1, n);
    mpz_sub(E->A, a4, tmp1);
    mpz_mod(E->A, E->A, n);
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
#if DEBUG_MANY_EC >= 2
    gmp_printf("N:=%Zd;\n", n);
    gmp_printf("A:=%Zd;\n", E->A);
    if(x0 != NULL){
	gmp_printf("x0:=%Zd;\n", x0);
	gmp_printf("wx0:=%Zd;\n", P->x);
	gmp_printf("y0:=%Zd;\n", P->y);
	printf("B:=(y0^2-wx0^3-A*wx0) mod N;\n");
    }
#endif
    mpz_clear(tmp1);
    mpz_clear(tmp2);
    mpz_clear(tmp3);
}

/* From a curve in Weierstrass form to a short form 
   WE:=[0,(1/4*c^2+1/4-1/2*c-b),0,(1/2*c*b-1/2*b),1/4*b^2]);
   We compute:
   a2 = 1/4*c^2+1/4-1/2*c-b = ((c-1)/2)^2-b
   a4 = 1/2*c*b-1/2*b = b*(c-1)/2
   a6 = (b/2)^2
*/
void
K2W246(mpz_t a2, mpz_t a4, mpz_t a6, mpz_t b, mpz_t c, mpz_t n, int compute_a6)
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
#if DEBUG_MANY_EC >= 2
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
   Sends Kubert curve E(b, c) with point (x0, y0) to short Weierstrass form:
   Y^2 = X^3 + A * X + B
*/
void
kubert_to_weierstrass(ec_curve_t E, mpz_t B, ec_point_t P, mpz_t b, mpz_t c, 
		      mpz_t x0, mpz_t y0, mpz_t n)
{
    mpz_t a2, a4, a6;
    int compute_a6 = (B != NULL);

    mpz_init(a2);
    mpz_init(a4);
    if(compute_a6)
	mpz_init(a6);
    K2W246(a2, a4, a6, b, c, n, compute_a6);
    /* second conversion */
    W2W(E, B, P, a2, a4, a6, x0, y0, n);
    mpz_clear(a2);
    mpz_clear(a4);
    if(compute_a6)
	mpz_clear(a6);
}

/* Kubert: put b = c. */
int
build_curves_with_torsion_Z5(mpz_t f, mpmod_t n, 
			     ec_curve_t *tE, ec_point_t *tP,
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
#if DEBUG_MANY_EC >= 2
	gmp_printf("x0:=%Zd;\nc:=%Zd;\ny0:=%Zd;\n", x0, c, y0);
	printf("cr:=1/2*x0*(4*x0+1)/(3*x0+1);\n");
#endif
	/* P:=WE![x0, y0, 1]; */
	/* convert to short Weierstrass form */
	ec_curve_init(tE[nc], n);
	tE[nc]->type = ECM_EC_TYPE_WEIERSTRASS;
	ec_point_init(tP[nc], tE[nc], n);
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
			     ec_curve_t *tE, ec_point_t *tP,
			     int umin, int umax, int nE)
{
    int u, ret = ECM_NO_FACTOR_FOUND, nc = 0;
    mpz_t A2, A1div2, x0, y0, cte, d, c, b, kx0, ky0;
    mpres_t tmp;
    ec_curve_t E;
    ec_point_t P, Q;

    /* Eaux = "1295/48", "-1079/864" */
    /* Paux = "2185/12", "-2458" */
    /* E and P consist of Montgomery residues */
    mpres_init(tmp, n);
    mod_from_rat_str(f, "1295/48", n->orig_modulus);
    mpres_set_z(tmp, f, n);
    ec_curve_init_set(E, tmp, ECM_EC_TYPE_WEIERSTRASS, n);
    ec_point_init(P, E, n);
    mod_from_rat_str(f, "2185/12", n->orig_modulus);
    mpres_set_z(P->x, f, n);
    mpz_set_str(f, "-2458", 10);
    mpres_set_z(P->y, f, n);
    mpres_set_ui(P->z, 1, n);
#if DEBUG_MANY_EC >= 2
    printf("Paux:=");
    pt_print(P, n);
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
    ec_point_init(Q, E, n);
    for(u = umin; u < umax; u++){
	/* update Qaux */
	mpz_set_ui(d, u);
	/* TODO: replace with ec_point_add, one of these days */
	if(ec_point_mul(Q, d, P, E, n) == 0){
	    printf("found factor during update of Q\n");
	    mpz_set(f, Q->x);
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
#if DEBUG_MANY_EC >= 2
	printf("(s, t)[%d]:=", u);
	pt_print(Q, n);
	printf(";\n");
#endif
	/* come back to plain (not Montgomery) residues */
	mpres_get_z(b, Q->x, n);
	mpres_get_z(c, Q->y, n);
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
	ec_curve_init(tE[nc], n);
	tE[nc]->type = ECM_EC_TYPE_WEIERSTRASS;
	ec_point_init(tP[nc], tE[nc], n);
	kubert_to_weierstrass(tE[nc], NULL, tP[nc], b, c, kx0, ky0, n->orig_modulus);
#if DEBUG_MANY_EC >= 2
	gmp_printf("E[%d]:=[%Zd];\n", nc, tE[nc]->A);
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
    ec_point_clear(P, E, n);
    ec_point_clear(Q, E, n);
    ec_curve_clear(E, n);
    mpz_clear(d);
    mpz_clear(c);
    mpz_clear(b);
    mpz_clear(kx0);
    mpz_clear(ky0);
    mpres_clear(tmp, n);
    return ret;
}

int
build_curves_with_torsion_Z9(mpz_t fac, mpmod_t n, ec_curve_t *tE, 
			     ec_point_t *tP, int umin, int umax, int nE)
{
    int u, ret = ECM_NO_FACTOR_FOUND, nc = 0;
    mpz_t A2, A1div2, x0, y0, cte, d, c, b, kx0, ky0;
    mpz_t f;
    mpres_t tmp;
    ec_curve_t E;
    ec_point_t P, Q;

    /* Eaux = [-9, 9] */
    /* Paux = [1, 1, 1] */
    mpres_init(tmp, n);
    mpz_init_set_str(f, "-9", 10);
    mpres_set_z(tmp, f, n);
    ec_curve_init_set(E, tmp, ECM_EC_TYPE_WEIERSTRASS, n);
    ec_point_init(P, E, n);
    mpz_set_str(f, "1", 10); 
    mpres_set_z(P->x, f, n);
    mpz_set_str(f, "1", 10);
    mpres_set_z(P->y, f, n);
    mpres_init(P->z, n);
    mpres_set_ui(P->z, 1, n);
#if DEBUG_MANY_EC >= 2
    printf("Paux:=");
    pt_print(P, n);
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
    ec_point_init(Q, E, n);
    for(u = umin; u < umax; u++){
	/* update Qaux */
	mpz_set_ui(d, u);
        if(ec_point_mul(Q, d, P, E, n) == 0){
	    printf("found factor during update of Q\n");
	    mpz_set(fac, Q->x);
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
#if DEBUG_MANY_EC >= 2
	printf("(s, t)[%d]:=", u);
	pt_print(Q, n);
	printf(";\n");
#endif
	mpres_get_z(b, Q->x, n);
	mpres_get_z(c, Q->y, n);
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
	ec_curve_init(tE[nc], n);
	tE[nc]->type = ECM_EC_TYPE_WEIERSTRASS;
	ec_point_init(tP[nc], tE[nc], n);
	kubert_to_weierstrass(tE[nc], NULL, tP[nc], b, c, kx0, ky0, n->orig_modulus);
	nc++;
	if(nc >= nE)
	    break;
    }
#if DEBUG_MANY_EC >= 2
    printf("Curves built\n");
    pt_many_print(tE, tP, nE, n);
#endif
    mpz_clear(A2);
    mpz_clear(A1div2);
    mpz_clear(x0);
    mpz_clear(y0);
    mpz_clear(cte);
    ec_point_clear(P, E, n);
    ec_point_clear(Q, E, n);
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
build_curves_with_torsion_Z10(mpz_t fac, mpmod_t n, ec_curve_t *tE, 
			      ec_point_t *tP, int umin, int umax, int nE)
{
    int u, ret = ECM_NO_FACTOR_FOUND, nc = 0;
    mpz_t A2, A1div2, x0, y0, cte, d, c, b, kx0, ky0;
    mpz_t f;
    mpres_t tmp;
    ec_curve_t E;
    ec_point_t P, Q;

    mpz_init(f);
    /* Eaux = [2/3, -53/108] */
    /* Paux = [2/3, 1/2, 1] */
    mpres_init(tmp, n);
    mod_from_rat_str(f, "2/3", n->orig_modulus);
    mpres_set_z(tmp, f, n);
    ec_curve_init_set(E, tmp, ECM_EC_TYPE_WEIERSTRASS, n);
    ec_point_init(P, E, n);
    mod_from_rat_str(f, "2/3", n->orig_modulus);
    mpres_set_z(P->x, f, n);
    mod_from_rat_str(f, "1/2", n->orig_modulus);
    mpres_set_z(P->y, f, n);
    mpres_set_ui(P->z, 1, n);
#if DEBUG_MANY_EC >= 2
    printf("P:=");
    pt_print(P, n);
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
    ec_point_init(Q, E, n);
    for(u = umin; u < umax; u++){
	/* update Qaux */
	mpz_set_ui(d, u);
	if(ec_point_mul(Q, d, P, E, n) == 0){
	    printf("found factor in Z10 (update of Q)\n");
	    mpz_set(fac, Q->x);
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
#if DEBUG_MANY_EC >= 2
	printf("(s, t)[%d]:=", u);
	pt_print(Q, n);
	printf(";\n");
#endif
	mpres_get_z(b, Q->x, n);
	mpres_get_z(c, Q->y, n);
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
	ec_curve_init(tE[nc], n);
	tE[nc]->type = ECM_EC_TYPE_WEIERSTRASS;
	ec_point_init(tP[nc], tE[nc], n);
	kubert_to_weierstrass(tE[nc], NULL, tP[nc], b, c, kx0, ky0, n->orig_modulus);
	nc++;
	if(nc >= nE)
	    break;
    }
#if DEBUG_MANY_EC >= 2
    printf("Curves built\n");
    pt_many_print(tE, tP, nE, n);
#endif
    mpz_clear(A2);
    mpz_clear(A1div2);
    mpz_clear(x0);
    mpz_clear(y0);
    mpz_clear(cte);
    ec_point_clear(P, E, n);
    ec_point_clear(Q, E, n);
    ec_curve_clear(E, n);
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
				ec_curve_t *tE, ec_point_t *tP,
				int umin, int umax, int nE)
{
    int u, nc = 0, ret = ECM_NO_FACTOR_FOUND;
    mpz_t tmp, a, b, alpha, beta, c, d, kx0, ky0, wx0;
    mpres_t tmp2;
    ec_curve_t E;
    ec_point_t P, Q;

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
    ec_curve_init_set(E, tmp2, ECM_EC_TYPE_WEIERSTRASS, n);
    ec_point_init(P, E, n);
    mpz_set_str(f, "12", 10); 
    mpres_set_z(P->x, f, n);
    mpz_set_str(f, "40", 10);
    mpres_set_z(P->y, f, n);
    mpres_set_ui(P->z, 1, n);

    ec_point_init(Q, E, n);
    for(u = umin; u < umax; u++){
	/* update Qaux */
	mpz_set_ui(d, u);
	if(ec_point_mul(Q, d, P, E, n) == 0){
	    printf("found factor in Z2xZ8 (update of Q)\n");
	    mpz_set(f, Q->x);
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
#if DEBUG_MANY_EC >= 2
	printf("(s, t)[%d]:=", u);
	pt_print(Q, n);
	printf(";\n");
#endif
	mpres_get_z(a, Q->x, n);
	mpres_get_z(b, Q->y, n);
	mpres_get_z(d, Q->z, n);
	if(mpz_invert(f, d, n->orig_modulus) == 0){
	    printf("found factor in Z2xZ8 (normalization)\n");
	    mpz_gcd(f, d, n->orig_modulus);
	    break;
	}
	mpz_mul(a, a, f);
	mpz_mod(wx0, a, n->orig_modulus);
	mpz_mul(b, b, f);
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
	K2W246(f, a, NULL, b, c, n->orig_modulus, 0);
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
	ec_curve_init(tE[nc], n);
	tE[nc]->type = ECM_EC_TYPE_MONTGOMERY;
	ec_point_init(tP[nc], tE[nc], n);
	if(mod_from_rat2(tE[nc]->A, f, tmp, n->orig_modulus) == 0){
            printf("found factor in Z2xZ8 (ma)\n");
	    mpz_set(f, tE[nc]->A);
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
        mod_from_rat2(tP[nc]->x, tE[nc]->A, tmp, n->orig_modulus);
	mpz_sub(tP[nc]->x, f, tP[nc]->x);
	mpz_mod(tP[nc]->x, tP[nc]->x, n->orig_modulus);
	/* my:=mb*ky0; */
#if DEBUG_MANY_EC >= 2
	gmp_printf("N:=%Zd;\n", n->orig_modulus);
	gmp_printf("ma:=%Zd;\n", tE[nc]->A);
	gmp_printf("kx0:=%Zd;\n", kx0);
	gmp_printf("ky0:=%Zd;\n", ky0);
	gmp_printf("mx0:=%Zd;\n", tP[nc]->x);
#endif
	nc++;
	if(nc >= nE)
	    break;
    }
#if DEBUG_MANY_EC >= 2
    printf("Curves built\n");
    pt_many_print(tE, tP, nE, n);
#endif
    ec_point_clear(P, E, n);
    ec_point_clear(Q, E, n);
    ec_curve_clear(E, n);
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
build_curves_with_torsion_Z3xZ3_DuNa(mpmod_t n, ec_curve_t *tE, ec_point_t *tP,
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
	ec_curve_init(tE[nc], n);
	tE[nc]->type = ECM_EC_TYPE_WEIERSTRASS;
	ec_point_init(tP[nc], tE[nc], n);
	W2W(tE[nc], NULL, tP[nc], a2, a4, NULL, x0, y0, n->orig_modulus);
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
				ec_curve_t *tE, ec_point_t *tP,
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
	ec_curve_init_set(tE[nc], D, ECM_EC_TYPE_HESSIAN, n);
	ec_point_init(tP[nc], tE[nc], n);
	mpz_set(tE[nc]->A, D);
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
				ec_curve_t *tE, ec_point_t *tP,
				int umin, int umax, int nE)
{
    int u, nc = 0, ret = ECM_NO_FACTOR_FOUND;
    ec_curve_t E;
    ec_point_t P, Q;
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
    ec_curve_init_set(E, tmp, ECM_EC_TYPE_WEIERSTRASS, n);
    ec_point_init(P, E, n);
    mpz_set_str(f, "2", 10);
    mpres_set_z(P->x, f, n);
    mpz_set_str(f, "2", 10);
    mpres_set_z(P->y, f, n);
    mpres_set_ui(P->z, 1, n);

    ec_point_init(Q, E, n);
    for(u = umin; u < umax; u++){
	/* update Qaux */
	mpz_set_ui(f, u);
	if(ec_point_mul(Q, f, P, E, n) == 0){
	    printf("found factor in Z3xZ6 (update of Q)\n");
	    mpz_set(f, Q->x);
	    ret = ECM_FACTOR_FOUND_STEP1;
	    break;
	}
#if DEBUG_MANY_EC >= 2
	printf("(s, t)[%d]:=", u);
	pt_print(Q, n);
	printf(";\n");
#endif
	mpres_get_z(tk, Q->x, n);
	mpres_get_z(sk, Q->y, n);
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
	ec_curve_init(tE[nc], n);
	tE[nc]->type = ECM_EC_TYPE_HESSIAN;
	ec_point_init(tP[nc], tE[nc], n);
	if(mod_from_rat2(tE[nc]->A, num, den, n->orig_modulus) == 0){
            printf("found factor in Z3xZ6 (D)\n");
            mpz_set(f, tE[nc]->A);
            ret = ECM_FACTOR_FOUND_STEP1;
            break;
        }
#if DEBUG_MANY_EC >= 1
	gmp_printf("D%d:=%Zd;\n", nc, tE[nc]->A);
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
build_curves_with_torsion_Z4xZ4(mpz_t f, mpmod_t n, ec_curve_t *tE,
				ec_point_t *tP,
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
        ec_curve_init(tE[nc], n);
        tE[nc]->type = ECM_EC_TYPE_MONTGOMERY;
        ec_point_init(tP[nc], tE[nc], n);
	if(mod_from_rat2(tE[nc]->A, tmp, x0, n->orig_modulus) == 0){
	    printf("\nFactor found durint init of Z4xZ4 (a)\n");
	    mpz_set(f, tE[nc]->A);
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
#if DEBUG_MANY_EC >= 2
	gmp_printf("N:=%Zd;\n", n);
	printf("nu:=%d;\n", nu);
	gmp_printf("tau:=%Zd;\n", tau);
	gmp_printf("lambda:=%Zd;\n", lambda);
	gmp_printf("a:=%Zd;\n", tE[nc]->A);
	gmp_printf("x0:=%Zd;\n", x0);
#endif
	/* x:=b*x0-a/3; not needed: y:=b*y0 */
	mpz_set_si(tmp, 3);
	mod_from_rat2(tP[nc]->x, tE[nc]->A, tmp, n->orig_modulus);
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
ec_force_point(ec_curve_t E, ec_point_t P, mpz_t B, long *x0, mpz_t n)
{
    mpz_t lambda;

    mpz_init(lambda);
    while(1){
	*x0 += 1;
	/* lambda <- x0^3+A*x0+B = (x0^2+A)*x0+B */
	mpz_set_si(lambda, *x0);
	mpz_mul_si(lambda, lambda, *x0);
	mpz_add(lambda, lambda, E->A);
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
	/* modify E */
	mpz_mul(E->A, E->A, P->y);
	mpz_mod(E->A, E->A, n);
	/* not really needed */
	mpz_mul(B, B, P->y);
	mpz_mul(B, B, lambda);
	mpz_mod(B, B, n);
	break;
    }
#if DEBUG_MANY_EC >= 2
    gmp_printf("newA:=%Zd;\n", E->A);
    gmp_printf("newB:=%Zd;\n", B);
#endif
    mpz_clear(lambda);
}

/* Source: Brier+Clavier or Kohel or Silverberg or Klein. */
int
build_curves_with_torsion_Z5xZ5(mpmod_t n, ec_curve_t *tE,
				ec_point_t *tP,
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
	ec_curve_init(tE[nc], n);
	tE[nc]->type = ECM_EC_TYPE_WEIERSTRASS;
	/* A = -(t^4-228*t^3+494*t^2+228*t+1)/48; */
	mpz_eval_poly(num, polyA, t, n->orig_modulus);
	mpz_sub_si(den, n->orig_modulus, 48);
	mod_from_rat2(tE[nc]->A, num, den, n->orig_modulus);
	/* B = (t^6+522*t^5-10005*t^4-10005*t^2-522*t+1)/864 */
	mpz_eval_poly(num, polyB, t, n->orig_modulus);
	mpz_set_si(den, 864);
	mod_from_rat2(B, num, den, n->orig_modulus);
	ec_point_init(tP[nc], tE[nc], n);
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
build_curves_with_torsion_Z2xZ10(mpz_t f, mpmod_t n, ec_curve_t *tE,
				 ec_point_t *tP,
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
#if DEBUG_MANY_EC >= 2
    printf("info:=[];\n");
#endif
    for(u = smin; u < smax; u++){
#if DEBUG_MANY_EC >= 2
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
#if DEBUG_MANY_EC >= 2
	gmp_printf("t:=%Zd;\n", t);
	gmp_printf("b:=%Zd;\n", b);
	gmp_printf("c:=%Zd;\n", c);
#endif	
	ec_curve_init(tE[nc], n);
	tE[nc]->type = ECM_EC_TYPE_WEIERSTRASS;
	ec_point_init(tP[nc], tE[nc], n);
	/* transform with no input point */
	kubert_to_weierstrass(tE[nc],B,NULL,b,c,NULL,NULL,n->orig_modulus);
#if DEBUG_MANY_EC >= 2
        gmp_printf("A:=%Zd;\n", tE[nc]->A);
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
build_curves_with_torsion_Z2xZ12(mpz_t f, mpmod_t n, ec_curve_t *tE,
				 ec_point_t *tP,
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
#if DEBUG_MANY_EC >= 2
    printf("info:=[];\n");
#endif
    for(u = smin; u < smax; u++){
#if DEBUG_MANY_EC >= 2
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
#if DEBUG_MANY_EC >= 2
	gmp_printf("t:=%Zd;\n", t);
	gmp_printf("b:=%Zd;\n", b);
	gmp_printf("c:=%Zd;\n", c);
#endif	
	ec_curve_init(tE[nc], n);
	tE[nc]->type = ECM_EC_TYPE_WEIERSTRASS;
	ec_point_init(tP[nc], tE[nc], n);
	/* transform with no input point */
	kubert_to_weierstrass(tE[nc],B,NULL,b,c,NULL,NULL,n->orig_modulus);
#if DEBUG_MANY_EC >= 2
        gmp_printf("A:=%Zd;\n", tE[nc]->A);
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

/* Assuming we can generate curves with given torsion using parameter s
   in interval [smin..smax].
*/
int
build_curves_with_torsion(mpz_t f, mpmod_t n, ec_curve_t *tE, ec_point_t *tP,
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
    else if(strcmp(torsion, "Z3xZ3_DuNa") == 0) /* over Q(sqrt(-3)) */
	return build_curves_with_torsion_Z3xZ3_DuNa(n, tE, tP, smin, smax, nE);
    else if(strcmp(torsion, "Z3xZ3") == 0) /* over Q(sqrt(-3)) */
	return build_curves_with_torsion_Z3xZ3(f, n, tE, tP, smin, smax, nE);
    else if(strcmp(torsion, "Z3xZ6") == 0) /* over Q(sqrt(-3)) */
	return build_curves_with_torsion_Z3xZ6(f, n, tE, tP, smin, smax, nE);
    /** interesting when p = 1 mod 4 **/
    else if(strcmp(torsion, "Z4xZ4") == 0) /* over Q(sqrt(-1)) */
	return build_curves_with_torsion_Z4xZ4(f, n, tE, tP, smin, smax, nE);
    /** interesting when p = 1 mod 5 **/
    else if(strcmp(torsion, "Z5xZ5") == 0) /* over Q(zeta5) */
	return build_curves_with_torsion_Z5xZ5(n, tE, tP, smin, smax, nE);
    /** forcing points: is this really interesting? **/
    else if(strcmp(torsion, "Z2xZ10") == 0)
	return build_curves_with_torsion_Z2xZ10(f, n, tE, tP, smin, smax, nE);
    else if(strcmp(torsion, "Z2xZ12") == 0)
	return build_curves_with_torsion_Z2xZ12(f, n, tE, tP, smin, smax, nE);
    else{
	printf("Unknown torsion group: %s\n", torsion);
	ret = ECM_ERROR;
    }
    return ret;
}

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

/* 
   OUTPUT: ECM_NO_FACTOR_FOUND
           ECM_PRIME_FAC_PRIME_COFAC
	   ECM_PRIME_FAC_COMP_COFAC
	   ECM_COMP_FAC_COMP_COFAC
	   ECM_COMP_FAC_PRIME_COFAC
  One ring to run them all.
*/
int
process_many_curves_loop(mpz_t tf[], int *nf, mpz_t n, double B1, 
			 ecm_params params,
			 char *fic_EP,
			 char *torsion, int smin, int smax, int nE,
			 int disc, mpz_t *sqroots,
			 char *savefilename)
{
    ec_curve_t tE[NCURVE_MAX];
    ec_point_t tP[NCURVE_MAX];
    mpmod_t modulus;
    int ret = 0, i, onebyone;

    onebyone = 1; /* mtyform; */
    while(1){
	/* cheating with the content of tE and tP */
	mpmod_init(modulus, n, ECM_MOD_DEFAULT);
	if(fic_EP != NULL)
	    ret = read_curves_from_file(&nE, tE, tP, tf, nf, modulus, 
					fic_EP, nE);
	else if(torsion != NULL)
	    ret = build_curves_with_torsion(tf[*nf],modulus,tE,tP,
					    torsion,smin,smax,nE);
	else if(disc != 0)
	    ret = build_curves_with_CM(tf[*nf],&nE,tE,tP,disc,modulus,sqroots);
	if(ret == ECM_NO_FACTOR_FOUND)
	    ret = process_many_curves(tf[*nf],modulus,B1,tE,tP,nE,params,
				      onebyone,savefilename);
	else{
	    printf("Quid? %d\n", ret);
	    break;
	}
	/* clear curves */
	for(i = 0; i < nE; i++){
	    ec_point_clear(tP[i], tE[i], modulus);
	    ec_curve_clear(tE[i], modulus);
	}
	mpmod_clear(modulus);
	/* inspect result */
	if(ret == ECM_PRIME_FAC_PRIME_COFAC){
	    *nf += 1;
	    break;
	}
	else if(ret == ECM_PRIME_FAC_COMP_COFAC){
	    printf("# start again with n/f\n");
	    mpz_tdiv_q(n, n, tf[*nf]);
	    *nf += 1;
	}
	else if(ret == ECM_COMP_FAC_PRIME_COFAC){
	    mpz_t C;

	    printf("# start again with f\n");
	    mpz_init(C);
	    mpz_tdiv_q(C, n, tf[*nf]);
	    mpz_set(n, tf[*nf]);
	    mpz_set(tf[*nf], C);
	    mpz_clear(C);
	    *nf += 1;
	}
	else if(ret == ECM_COMP_FAC_COMP_COFAC){
	    mpz_t f;

	    mpz_init_set(f, tf[*nf]);
	    /* update n right now */
	    mpz_tdiv_q(n, n, f);
	    gmp_printf("# recursive call for f=%Zd\n", f);
	    process_many_curves_loop(tf, nf, f, B1, params, fic_EP,
				     torsion, smin, smax, nE,
				     disc, sqroots, savefilename);
	    /* there is always some cofactor to store */
	    mpz_set(tf[*nf], f);
	    *nf += 1;
	    printf("# start again with n/f\n");
	}
	else /* something happened */
	    break;
    }
    return ret;
}

/* Assume b^n = 1 mod N. */
int
odd_square_root_mod_N(mpz_t f, mpz_t *sqroots, int b, int n, int q, mpz_t N)
{
    mpz_t zeta, tmp, tmp2;
    int np = n, e = 0, *tab, x, ret = ECM_NO_FACTOR_FOUND;

    while(np % q == 0){
	e++;
	np /= q;
    }
    printf("# n = %d = %d^%d * %d\n", n, q, e, np);
    mpz_init_set_ui(zeta, b);
    mpz_powm_ui(zeta, zeta, np, N);
    if(mpz_cmp_ui(zeta, 1) == 0){
	printf("# missed: zeta == 1\n");
    }
    else{
	/* look for k s.t. zeta^{q^k} = 1 */
	mpz_init_set(tmp, zeta);
	do{
	    mpz_set(zeta, tmp);
	    mpz_powm_ui(tmp, tmp, q, N);
	} while(mpz_cmp_ui(tmp, 1) != 0);
	/*	gmp_printf("# zeta_%d = %Zd\n", q, zeta);*/
	mpz_sub_si(f, zeta, 1);
	mpz_gcd(f, f, N);
	if(mpz_cmp_ui(f, 1) != 0){
	    printf("# Factor found (gcd(zeta-1, N)): ");
	    mpz_out_str(stdout, 10, f);
	    printf("\n");
	    ret = ECM_FACTOR_FOUND_STEP1;
	    goto end_of_odd_sqrt;
	}
	/* compute eta0 = sum zeta^R */
	tab = (int *)malloc((q+1) * sizeof(int));
	memset(tab, 0, (q+1) * sizeof(int));
	for(x = 1; x < q; x++)
	    tab[(x*x) % q] = 1;
	mpz_set_ui(tmp, 0);
	mpz_init(tmp2);
	for(x = 1; x < q; x++)
	    if(tab[x] == 1){
		mpz_powm_ui(tmp2, zeta, x, N);
		mpz_add(tmp, tmp, tmp2);
	    }
	mpz_add(tmp, tmp, tmp);
	mpz_add_ui(tmp, tmp, 1);
	mpz_mod(tmp, tmp, N);
	mpz_mul(tmp2, tmp, tmp);
	if(q % 4 == 1){
	    gmp_printf("# sqrt(%d) = %Zd\n", q, tmp);
	    mpz_sub_si(tmp2, tmp2, q);
	}
	else{
	    gmp_printf("# sqrt(-%d) = %Zd\n", q, tmp);
	    mpz_add_si(tmp2, tmp2, q);
	}
	mpz_mod(tmp2, tmp2, N);
	if(mpz_sgn(tmp2) == 0)
	    mpz_init_set(sqroots[0], tmp);
	else{
	    gmp_printf("Bad check: %Zd\n", tmp2);
	    gmp_printf("N:=%Zd;\n", N);
	    ret = ECM_ERROR;
	}
	mpz_clear(tmp);
	mpz_clear(tmp2);
	free(tab);
    }
 end_of_odd_sqrt:
    mpz_clear(zeta);
    return ret;
}

/* Consider M = b^n+1 if n > 0, M = b^(-n)-1 otherwise.
   N is supposed to be a *primitive* cofactor of M.
   Then find a special cocktail of CM curves a` la Atkin.
   TODO: solve the B1 choice problem...! Force disc?
 */
int
process_special_blend(mpz_t tf[], int *nf, mpz_t N, int b, int n, double B1, 
		      ecm_params params, char *savefilename)
{
    int sgn = 1, disc1 = 0, q, nn, disc, i, j, ret = ECM_NO_FACTOR_FOUND;
    int tabd[][3] = {{-3, 5, 0},
#if 0
 {-4, 5, 0}, {8, -3, 0}, {5, -7, 0},
		     {-8, 5, 0}, {-3, 17, 0}, {-4, 13, 0}, {8, -11, 0},
		     {-7, 13, 0}, {5, -23, 0}, {-3, 41, 0},
148, 187, 232, 235, 267, 403, 427,
		     /* h = g = 4 */
		  84, 120, 132, 168, 195, 228, 280, 312, 340, 372, 408, 
		  435, 483, 520, 532, 555, 595, 627, 708, 715, 760, 795,
		  1012, 1435,
		     /* h = g = 8 */
		  420, 660, 840, 1092, 1155, 1320, 1380, 1428, 1540, 1848, 
		  1995, 3003, 3315,
		     /* h = g = 16 */
		  5460
#endif
		     {0, 0, 0}
    };
    mpz_t sqroots[10];

    if(n < 0){
	sgn = -1;
	n = -n;
    }
    /* look for discriminants of class number 1, hence rational CM curves */
    if(sgn == -1){
	if(n % 2 == 1){
	    /* b^(2*k+1) = 1 mod N => (b^(k+1))^2 = b mod N */
	}
    }
    else{
	nn = 2 * n;
	if(n % 2 == 0){
	    /* b^(2*k) = -1 mod N => (b^k)^2 = -1 mod N */
	    disc1 = -4;
	    /* set squareroot of -1 */
	    mpz_init_set_si(sqroots[0], b);
	    mpz_powm_ui(sqroots[0], sqroots[0], n>>1, N);
	}
	else{
	    /* b^(2*k+1) = -1 mod N => (b^(k+1))^2 = -b mod N */
	    if(b == 2)
		disc1 = -8;
	    else if(b == 3 || b == 7 || b == 11)
		disc1 = -b;
	    if(disc1 != 0){
		/* set squareroot of -b */
		mpz_init_set_si(sqroots[0], b);
		mpz_powm_ui(sqroots[0], sqroots[0], (n>>1)+1, N);
	    }
	}
    }
#if 0 /* TMP, TMP, time to try the use of twists */
    /* each odd prime factor q of n can lead to sqrt(q*) */
    printf("# Using factors of n=%d\n", n);
    for(i = 0; tabd[i][0] != 0; i++){
	disc = tabd[i][0];
	for(j = 1; tabd[i][j] != 0; j++)
	    disc *= tabd[i][j];
	if(n % abs(disc) == 0){
	    /* for the time being, works only for h = 2 */
	    printf("# I can use tabd[%d] = %d\n", i, disc);
	    for(j = 0; tabd[i][j] != 0; j++){
		q = tabd[i][j];
		if(abs(q) % 2 == 1){
		    ret = odd_square_root_mod_N(tf[*nf],sqroots,b,nn,abs(q),N);
		    if(ret != ECM_NO_FACTOR_FOUND)
			break;
		}
		else{
		    printf("# I cannot use even primes right now!\n");
		}
	    }
	    if(ret == ECM_NO_FACTOR_FOUND){
		printf("# Let us use disc=%d\n", disc);
		ret = process_many_curves_loop(tf, nf, N, B1, params,
					       NULL, NULL, 0, 0, 1,
					       disc, sqroots,
					       savefilename);
	    }
	}
	if(ret != ECM_NO_FACTOR_FOUND)
	    break;
    }
#endif
    if(ret == ECM_NO_FACTOR_FOUND && disc1 != 0){
	printf("# Let us use disc=%d\n", disc1);
	return process_many_curves_loop(tf, nf, N, B1, params,
					NULL, NULL, 0, 0, 1,
					disc1, sqroots,
					savefilename);
    }
    return ret;
}

static void
usage(char *cmd)
{
    printf("Usage: %s -inp file_N -B1 B1 -curves file_C", cmd);
    printf(" -torsion T -smin smin -smax smax\n");
    printf("  -inp file_N    numbers to be factored, one per line\n");
    printf("                 file_N can be '-', in which case stdin is used\n");
    printf("  -curves file_C curves to be used, format '[M|W|H] A x0 y0' per line\n");
    printf("                 M=Montgomery, W=Weierstrass, H=Hessian\n");
    printf("  -disc D        uses CM curves with discriminant D\n");
    printf("  -b b           for numbers b^n+/-1 (activates some special code)\n");
    printf("  -h, --help     Prints this help and exit.\n");
}

#define NFMAX 100

int
main(int argc, char *argv[])
{
    mpz_t N, tf[NFMAX];
    int res = 0, smin = -1, smax = -1, ncurves = 0, method = ECM_ECM;
    int nf = 0, i;
    double B1 = 0.0;
    int disc = 0, b = 0, n = 0;
    char *infilename = NULL, *curvesname = NULL, *torsion = NULL;
    char buf[10000];
    FILE *infile = NULL;
    char *savefilename = NULL;
    ecm_params params;

    /* print args */
    printf("# ARGS: %s", argv[0]);
    for(i = 1; i < argc; i++)
	printf(" %s", argv[i]);
    printf("\n");
    /* look for options */
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
	else if ((argc > 2) && (strcmp (argv[1], "-disc") == 0)){
	    disc = atol(argv[2]);
	    ncurves = 1;
	    argv += 2;
	    argc -= 2;
	}
	else if ((argc > 2) && (strcmp (argv[1], "-b") == 0)){
	    b = atoi(argv[2]);
	    argv += 2;
	    argc -= 2;
	}
	else if ((argc > 2) && (strcmp (argv[1], "-save") == 0)){
	    savefilename = argv[2];
	    argv += 2;
	    argc -= 2;
	}
	else if (strcmp (argv[1], "-pm1") == 0){
	    method = ECM_PM1;
	    argv++;
	    argc--;
	}
	else if (strcmp (argv[1], "-pp1") == 0){
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
    if(curvesname == NULL && torsion == NULL && disc == 0 && b == 0){
	fprintf (stderr, "No curve file given\n");
	exit (EXIT_FAILURE);
    }
    if(curvesname != NULL && torsion != NULL){
	fprintf (stderr, "Cannot have -curves and -torsion at the same time.\n");
	exit (EXIT_FAILURE);
    }
    if((disc != 0) && ((torsion != NULL) || (curvesname != NULL))){
	fprintf (stderr, "Cannot have -disc and -curves or -torsion at the same time.\n");
	exit (EXIT_FAILURE);
    }
    if(torsion != NULL && ncurves == 0){
	fprintf (stderr, "You must provide ncurves != 0 with -torsion.\n");
	exit (EXIT_FAILURE);
    }
    
    if(torsion != NULL)
	printf("GMP-ECM [torsion=%s:%d-%d]\n", torsion, smin, smax);
    if(disc != 0)
	printf("GMP-ECM [CM=%d]\n", disc);
    
    mpz_init (N);
    for(i = 0; i < NFMAX; i++)
	mpz_init(tf[i]); /* for potential factors */
    ecm_init(params);
#if DEBUG_MANY_EC >= 2
    params->verbose = 2;
#else
    params->verbose = 1;
#endif
#if 1
    compute_s_4_add_sub(params->batch_s, (unsigned long)B1, disc);
#endif
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
      if(b != 0){
	  /* line should be: "b n N" */
	  b = atoi(buf);
	  fscanf(infile, "%d", &n);
	  fscanf(infile, "%s", buf);
      }
      if(mpz_set_str (N, buf, 10)){
	  fprintf (stderr, "Invalid number: %s\n", argv[1]);
	  exit (1);
      }
      if(b != 0){
	  nf = 0;
	  res = process_special_blend(tf,&nf,N,b,n,B1,params,savefilename);
      }
      if(method == ECM_ECM){
	  nf = 0;
	  res = process_many_curves_loop(tf, &nf, N, B1, params,
					 curvesname,
					 torsion, smin, smax, ncurves,
					 disc, NULL,
					 savefilename);
      }
    }
#if 0	  
    printf("List of factors:\n");
    for(i = 0; i < nf; i++)
	gmp_printf("%Zd\n", tf[i]);
#endif
    ecm_clear(params);
    if(infile != stdin)
	fclose(infile);
    for(i = 0; i < NFMAX; i++)
	mpz_clear (tf[i]);
    mpz_clear (N);
    
    return res;
}
