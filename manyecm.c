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

#include "addlaws.h"
#include "torsions.h"
#include "cmecm.h"

#define DEBUG_MANY_EC 0

#define NCURVE_MAX 1000

/* fall back on traditional ECM.
   TODO: use chkfile also.
 */
int
process_one_curve(mpz_t f, mpz_t N, double B1, mpz_t B2,
		  ecm_params params, ec_curve_t E, ec_point_t P)
{
    int ret;
    double B2scale = 1.0;

    /* Taken from main.c; no comment */
    /* Here's an ugly hack to pass B2scale to the library somehow.
       It gets piggy-backed onto B1done */
    params->B1done = params->B1done + floor (B2scale * 128.) / 134217728.; 

    /* if B2 = ECM_DEFAULT_B2, compute it automatically from B1: 
       no freedom on B2! */
    mpz_set(params->B2, B2);
    /* will be set to B1 */
    mpz_set_si(params->B2min, ECM_DEFAULT_B2);

    mpz_set(params->x, P->x);
    mpz_set(params->sigma, E->A); /* humf */

    if(E->type == ECM_EC_TYPE_MONTGOMERY)
	params->sigma_is_A = 1;
    else{
	params->sigma_is_A = -1;
	mpz_set(params->y, P->y);
    }
    params->E = E;
    ret = ecm_factor(f, N, B1, params);
    return ret;
}

/* OUTPUT: ECM_PRIME_FAC_PRIME_COFAC if f prp, N/f prp
           ECM_PRIME_FAC_COMP_COFAC if f prp, N/f composite
	   ECM_COMP_FAC_PRIME_COFAC if f composite, N/f prp
           ECM_COMP_FAC_COMP_COFAC if f composite, N/f composite
 */
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
    printf("CheckE:=procedure(E, P, info)\n");
    printf("    printf \"#E[%%o]=%%o\\n\", info, Factorization(#E);\n");
    printf("    gen:=Generators(E);\n");
    printf("    printf \"ords=%%o\\n\", ");
    printf("[Factorization(Order(g)) : g in gen];\n");
    printf("    printf \"ord(P)=%%o\\n\", ");
    printf("Factorization(Order(E!P));\n");
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
	else if(tE[i]->type == ECM_EC_TYPE_WEIERSTRASS_AFF
		|| tE[i]->type == ECM_EC_TYPE_WEIERSTRASS_HOM){
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
	printf("CheckE(E[%d], P[%d], infos[%d]);\n", i+1, i+1, i+1);
    }
}

/* TODO: better control of B2 + dichotomy (cf. #B2) */
int
one_curve_at_a_time(mpz_t f, char *ok, ec_curve_t *tE, ec_point_t *tP, int nE,
		    mpz_t N, ecm_params params, double B1, mpz_t B2,
		    char *savefilename)
{
    double tmpB1, tmpB2, B2g = 0, B2d = 0, dB2 = 0; /* 1e9; #B2 */
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
	tmpB2 = dB2;
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
	    ret = process_one_curve(f,N,tmpB1,B2,params,tE[i],tP[i]);
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
#if DEBUG_MANY_EC >= 0
		if(ret == ECM_PRIME_FAC_PRIME_COFAC 
		   || ret == ECM_PRIME_FAC_COMP_COFAC){
		    /* output Magma lines to check #E's mod f */
		    /*dump_curves(tE, tP, nE, f);*/
		    printf("infos:=[\"E%d\"];\n", i);
		    dump_curves(tE+i, tP+i, 1, f);
		}
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
process_many_curves(mpz_t f, mpmod_t n, double B1, mpz_t B2,
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
				  B1, B2, savefilename);
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
	    ret = process_one_curve(f, n->orig_modulus, B1, B2, params,
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
	ec_curve_init(tE[*nE], ECM_EC_TYPE_WEIERSTRASS_AFF, n);
	if(Etype == 'W'){
	    if(fscanf(ifile, "%s %s %s", bufA, bufx, bufy) == EOF)
		break;
	    tE[*nE]->type = ECM_EC_TYPE_WEIERSTRASS_HOM;
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

/* 
   OUTPUT: ECM_NO_FACTOR_FOUND
           ECM_PRIME_FAC_PRIME_COFAC
	   ECM_PRIME_FAC_COMP_COFAC
	   ECM_COMP_FAC_COMP_COFAC
	   ECM_COMP_FAC_PRIME_COFAC
  One ring to run them all.
*/
int
process_many_curves_loop(mpz_t tf[], int *nf, mpz_t n, double B1, mpz_t B2,
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
	/* cheating with the content of tE and tP that are first defined
	   over Z/nZ without residues
	 */
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
	    ret = process_many_curves(tf[*nf],modulus,B1,B2,tE,tP,nE,params,
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
	    process_many_curves_loop(tf, nf, f, B1, B2, params, fic_EP,
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
	    printf("# Factor found (gcd(zeta_%d-1, N)): ", q);
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

/* b^(2*k) = -1 => (b^k)^2 = -1 mod N */
static
int psb_minus_even(mpz_t sqroots[], int b, int k, mpz_t N)
{
    /* set squareroot of -1 */
    mpz_init_set_si(sqroots[0], b);
    mpz_powm_ui(sqroots[0], sqroots[0], k, N);
    if(k % 2 == 0){
	/* zeta8 = b^(k/2) = (1+zeta4)/sqrt(2) 
	   => sqrt(2) = (1+zeta4)/zeta8 */
	mpz_init_set_si(sqroots[1], b);
	mpz_powm_ui(sqroots[1], sqroots[1], k>>1, N);
	mpz_invert(sqroots[1], sqroots[1], N);
	mpz_add_si(sqroots[0], sqroots[0], 1);
	mpz_mul(sqroots[1], sqroots[1], sqroots[0]);
	mpz_mod(sqroots[1], sqroots[1], N);
	mpz_sub_si(sqroots[0], sqroots[0], 1);
    }
    return -4;
}

/* b^(2*k+1) = -1 => (b^(k+1))^2 = -b mod N */
static
int psb_minus_odd(mpz_t sqroots[], int b, int k, mpz_t N)
{
    int disc1 = 0;

    if(b == 2)
	disc1 = -8;
    else if(b == 3 || b == 7 || b == 11)
	disc1 = -b;
    if(disc1 != 0){
	/* set squareroot of -b */
	mpz_init_set_si(sqroots[0], b);
	mpz_powm_ui(sqroots[0], sqroots[0], k+1, N);
    }
    return disc1;
}

/* OUTPUT: ECM_NO_FACTOR_FOUND or ECM_FACTOR_FOUND_STEP1 in very rare cases! */
static
int prepare_squareroots(int *disc1, mpz_t f, mpz_t sqroots[], int b, int n, int sgn, mpz_t N)
{
    int k, discr = 0, ret = ECM_NO_FACTOR_FOUND;
    mpz_t tmp, tmp2;

    if(sgn == -1){
	/* b^n = 1 mod N */
	if(n % 2 == 0){
	    /* b^(2*k) = 1 => try to find smallest power */
	    k = n >> 1;
	    while(k % 2 == 0)
		k >>= 1;
	    mpz_init_set_si(tmp, b);
	    mpz_powm_ui(tmp, tmp, k, N);
	    mpz_init_set_si(tmp2, 0);
	    while(mpz_cmp_ui(tmp, 1) != 0){
		mpz_set(tmp2, tmp);
		mpz_mul(tmp, tmp, tmp);
		mpz_mod(tmp, tmp, N);
		k <<= 1;
	    }
	    /* at this point, b^k == 1 */
	    gmp_printf("# %d^%d = 1 mod %Zd;\n", b, k, N);
	    if(k % 2 == 1){
		/* b^(2*r+1) = 1 mod N => (b^(r+1))^2 = b mod N */
		printf("# case k=%d odd\n", k);
	    }
	    else{
		/* b^(2*r) = 1 */
		mpz_add_si(tmp2, tmp2, 1);
		if(mpz_cmp(tmp2, N) == 0){
		    /* case b^r = -1 */
		    printf("# %d^%d = -1 mod N;\n", b, k>>1);
		    if(k % 4 == 0)
			/* b^(2*s) = -1 */
			discr = psb_minus_even(sqroots, b, k>>2, N);
		    else
			/* b^(2*s+1) = -1 */
			discr = psb_minus_odd(sqroots, b, k>>2, N);
		}
		else{
		    /* we have a factor, since tmp2^2 = 1, tmp2 != -1, +1 */
		    mpz_sub_si(tmp2, tmp2, 1);
		    mpz_gcd(f, tmp2, N);
		    gmp_printf("Factor!! %Zd\n", f);
		    return ECM_FACTOR_FOUND_STEP1;
		}
	    }
	    mpz_clear(tmp);
	    mpz_clear(tmp2);
	}
	else{
	    /* b^(2*k+1) = 1 mod N => (b^(k+1))^2 = b mod N */
	    k = n>>1;
	    mpz_init_set_si(sqroots[0], b);
	    mpz_powm_ui(sqroots[0], sqroots[0], k+1, N);
	}
    }
    else{
	/* b^n = -1 mod N */
	if(n % 2 == 0)
	    /* b^(2*k) = -1 mod N => (b^k)^2 = -1 mod N */
	    discr = psb_minus_even(sqroots, b, n>>1, N);
	else
	    /* b^(2*k+1) = -1 mod N => (b^(k+1))^2 = -b mod N */
	    discr = psb_minus_odd(sqroots, b, n>>1, N);
    }
    *disc1 = discr;
    return ret;
}

/* Consider M = b^n+1 if n > 0, M = b^(-n)-1 otherwise.
   N is supposed to be a *primitive* cofactor of M.
   Then find a special cocktail of CM curves a` la Atkin.
   To solve the B1 problem, only consider (b, n)'s s.t. disc(b, n) = discref.
 */
int
process_special_blend(mpz_t tf[], int *nf, mpz_t N, int b, int n, int discref,
		      double B1, mpz_t B2, 
		      ecm_params params, char *savefilename)
{
    int sgn = 1, disc1 = 0, q, nn, disc, i, j, k;
    int ret = ECM_NO_FACTOR_FOUND;
    int tabd[][4] = {{-15, -3, 5, 0}, 
#if 0
		     {-20, -4, 5, 0}, {-24, 8, -3, 0},
		     {-35, 5, -7, 0}, {-40, -8, 5, 0}, {-51, -3, 17, 0},
		     {-52, -4, 13, 0}, {-88, 8, -11, 0},
		     {-91, -7, 13, 0}, {-115, 5, -23, 0}, {-123, -3, 41, 0},
		     {-148, 0, 0, 0}, {-187, 0, 0, 0}, 
		     {-232, 0, 0, 0}, {-235, 0, 0, 0},
		     {-267, 0, 0, 0}, {-403, 0, 0, 0}, {-427, 0, 0, 0},
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
    /* try discriminants of class number 1, hence rational CM curves */
    ret = prepare_squareroots(&disc1, tf[0], sqroots, b, n, sgn, N);
    if(ret != ECM_NO_FACTOR_FOUND)
	return ret;
    if(disc1 == discref){
	gmp_printf("# Let us use disc=%d with B1=%1.0f B2=%Zd\n",disc1,B1,B2);
	ret = process_many_curves_loop(tf, nf, N, B1, B2, params,
				       NULL, NULL, 0, 0, 1,
				       disc1, sqroots,
				       savefilename);
	if(ret != ECM_NO_FACTOR_FOUND){
	    if(disc1 != 0)
		mpz_clear(sqroots[0]); /* really? */
	    return ret;
	}
    }
    if(disc1 != 0)
	mpz_clear(sqroots[0]); /* really? */
#if 0 /* not ready yet */
	nn = 2 * n; /* used below */
    /* each odd prime factor q of n can lead to sqrt(q*) */
    for(i = 0; tabd[i][0] != 0; i++){
	disc = tabd[i][0];
	if(disc != discref)
	    continue;
	if(n % abs(disc) == 0){
	    /* for the time being, works only for h = 2 */
	    printf("# I can use tabd[%d] = %d\n", i, disc);
	    for(j = 1; tabd[i][j] != 0; j++){
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
		ret = process_many_curves_loop(tf, nf, N, B1, B2, params,
					       NULL, NULL, 0, 0, 1,
					       disc, sqroots,
					       savefilename);
	    }
	}
	if(ret != ECM_NO_FACTOR_FOUND)
	    break;
    }
#endif
    return ret;
}

static void
usage(char *cmd)
{
    printf("Usage: %s -inp file_N -B1 B1 -B2 B2 -curves file_C", cmd);
    printf(" -torsion T -smin smin -smax smax\n");
    printf("  -inp file_N    numbers to be factored, one per line\n");
    printf("                 file_N can be '-', in which case stdin is used\n");
    printf("  -curves file_C curves to be used, format '[M|W|H] A x0 y0' per line\n");
    printf("                 M=Montgomery, W=Weierstrass, H=Hessian\n");
    printf("  -disc D        uses CM curves with discriminant D\n");
    printf("  -b b           for numbers b^n+/-1 (activates some special code; b=1 for any b in the file)\n");
    printf("  -h, --help     Prints this help and exit.\n");
}

#define NFMAX 100

int
main(int argc, char *argv[])
{
    mpz_t N, tf[NFMAX], B2;
    int res = 0, smin = -1, smax = -1, ncurves = 0, method = ECM_ECM;
    int nf = 0, i, bb = 0;
    double B1 = 0.0, dB2 = 0.0;
    int disc = 0, b = 0, n = 0;
    char *infilename = NULL, *curvesname = NULL, *torsion = NULL;
    char buf[10000], c;
    FILE *infile = NULL;
    char *savefilename = NULL;
    ecm_params params;

    /* print args */
    printf("# ARGS: %s", argv[0]);
    for(i = 1; i < argc; i++)
	printf(" %s", argv[i]);
    printf("\n");
    mpz_init_set_si(B2, ECM_DEFAULT_B2);
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
	else if ((argc > 2) && (strcmp (argv[1], "-B2") == 0)){
	    dB2 = atof(argv[2]);
	    mpz_set_d(B2, dB2);
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
    /*    params->verbose = OUTPUT_DEVVERBOSE;*/
    params->verbose = OUTPUT_NORMAL;
#endif
#if 1
    compute_s_4_add_sub(params->batch_s, (unsigned long)B1, disc);
#endif
    while(fscanf(infile, "%s", buf) != EOF){
      /* read number */
      if(buf[0] == '#'){
	  /* print till end of line */
	  printf("%s", buf);
	  while((c = getc(infile)) != '\n')
	      printf("%c", c);
	  printf("\n");
	  continue;
      }
      if(b != 0){
	  /* line should be: "b n[+/-/L/M] N" */
	  bb = atoi(buf);
	  /* decode */
	  fscanf(infile, "%s", buf);
	  c = buf[strlen(buf)-1];
	  buf[strlen(buf)-1] = '\0';
	  n = atoi(buf);
	  if(c == '-')
	      n = -n;
	  else if(c == 'L' || c == 'M'){
	      if(bb == 5)
		  n = -n;
	  }
	  else if(c != '+'){
	      printf("#!# unknown suffix: %c\n", c);
	      break;
	  }
	  /* read N */
	  fscanf(infile, "%s", buf);
	  printf("# I read: b=%d n=%d c=%c\n", bb, n, c);
	  if((b > 1) && (bb != b))
	      continue;
      }
      if(mpz_set_str (N, buf, 10)){
	  fprintf (stderr, "Invalid number: %s\n", argv[1]);
	  exit (1);
      }
      if(method == ECM_ECM){
	  if(b != 0){
	      nf = 0;
	      res = process_special_blend(tf,&nf,N,bb,n,disc,B1,B2,
					  params,savefilename);
	  }
	  if(res == ECM_NO_FACTOR_FOUND
	     && (curvesname != NULL || torsion != NULL)){
	      nf = 0;
	      res = process_many_curves_loop(tf, &nf, N, B1, B2, params,
					     curvesname,
					     torsion, smin, smax, ncurves,
					     disc, NULL,
					     savefilename);
	  }
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
	mpz_clear(tf[i]);
    mpz_clear(N);
    mpz_clear(B2);    
    return res;
}
