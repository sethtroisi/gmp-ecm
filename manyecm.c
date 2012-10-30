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
#define ONE_CURVE_AT_A_TIME /* PROTECTION */

#define NCURVE_MAX 1000

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

int
compute_all_inverses(mpres_t *inv, mpres_t *x, int nx, mpmod_t n, char *takeit)
{
    int i;

    /* plain version, to debug the architecture */
    for(i = 0; i < nx; i++){
	if(takeit[i] == 0)
	    continue;
	if(!mpres_invert(inv[i], x[i], n)){
	    mpres_gcd(inv[0], x[i], n);
	    return 0;
	}
    }
    return 1;
}

int
pt_many_common(curve *tER, curve *tEP, curve *tEQ, int nEP, mpmod_t n, 
	       mpres_t *num, mpres_t *den, mpres_t *inv, char *takeit)
{
    int i;

    if(compute_all_inverses(inv, den, nEP, n, takeit) == 0){
	mpz_set(tER[0].x, inv[0]);
	return 0;
    }
    for(i = 0; i < nEP; i++){
	if(takeit[i] == 0)
	    continue;
	/* l:=(inv[i]*num[i]) mod N; */
	mpres_mul(num[i], num[i], inv[i], n);
	mpres_sqr(den[i], num[i], n);
	/* x:=(l^2-P[1]-Q[1]) mod N; */
	mpres_sub(den[i], den[i], tEP[i].x, n);
	mpres_sub(den[i], den[i], tEQ[i].x, n);
	/* tR[i]:=[x, (l*(P[1]-x)-P[2]) mod N, 1]; */
	mpres_sub(tER[i].x, tEP[i].x, den[i], n);
	mpres_mul(tER[i].x, tER[i].x, num[i], n);
	mpres_sub(tER[i].y, tER[i].x, tER[i].y, n);
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
    free(takeit);
    return res;
}

int
pt_many_add(curve *tER, curve *tEP, curve *tEQ, int nEP, mpmod_t n, 
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
	    pt_assign(tER+i, tEQ+i, n);
	}
	else if(pt_is_equal(tEP+i, tEQ+i)){
	    /* we should double */
	    if(mpz_sgn(tEP[i].y) == 0){
		/* 2 * P[i] = O_E */
		takeit[i] = 0;
		pt_set_to_zero(tEP+i, n);
		printf("# [2] * P[%d] = O_E\n", i);
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
	mpres_sub(num[i], tEQ[i].y, tEP[i].y, n);
	mpres_sub(den[i], tEQ[i].x, tEP[i].x, n);
    }
    res = pt_many_common(tER, tEP, tEQ, nEP, n, num, den, inv, takeit);
    free(takeit);
    return res;
}

/* tEQ[i] <- e * tEP[i] */
/* If a factor is found, it is put back in tEP[0].x */
int
pt_many_mul(curve *tEQ, curve *tEP, int nEP, mpz_t e, mpmod_t n, 
	    mpres_t *num, mpres_t *den, mpres_t *inv, char *ok)
{
  size_t l;
  int negated = 0, status = 1;

  if (mpz_sgn (e) == 0)
    {
	pt_many_set_to_zero(tEP, nEP, n);
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
  if (negated)
    mpz_neg (e, e);
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

#ifdef ONE_CURVE_AT_A_TIME
int
one_curve_at_a_time(mpz_t f, curve *tEP, int nc, mpz_t n, double B1)
{
    ecm_params params;
    int ret, i;

    ecm_init(params);
    params->sigma_is_A = -1;
    /* process curves one at a time */
    for(i = 0; i < nc; i++){
	params->B1done = 1.0;
	mpz_set (params->x, tEP[i].x);
	mpz_set (params->y, tEP[i].y);
	mpz_set (params->sigma, tEP[i].A); /* humf */
	ret = ecm_factor(f, n, B1, params);
	if(ret > 0){ /* humf */
	    printf("******************** Factor found with E[%d]: ", i);
	    mpz_out_str (stdout, 10, f);
	    printf ("\n");
	}
	else if(ret == ECM_ERROR){
	    printf("Error for curve %d\n", i);
	}
    }
    ecm_clear(params);
    return ret;
}
#else
/* Using parallelism.
   Copied from classical ecm_stage1.
 */
int
all_curves_at_once(mpz_t f, curve *tEP, int nEP, mpmod_t n, 
		   double B1, double *B1done, int (*stop_asap)(void),
		   char *chkfilename)
{
  curve tEQ[NCURVE_MAX];
  mpz_t num[NCURVE_MAX], den[NCURVE_MAX], inv[NCURVE_MAX], e;
  double p = 0.0, r, last_chkpnt_p;
  int ret = ECM_NO_FACTOR_FOUND;
  long last_chkpnt_time;
  int i;
  char *ok = (char *)malloc(nEP * sizeof(char));

  mpz_init(e);
  memset(ok, 1, nEP);
  for(i = 0; i < nEP; i++){
      mpres_init(tEQ[i].A, n); mpres_set(tEQ[i].A, tEP[i].A, n);
      mpres_init(tEQ[i].x, n); mpres_set(tEQ[i].x, tEP[i].x, n);
      mpres_init(tEQ[i].y, n); mpres_set(tEQ[i].y, tEP[i].y, n);
      mpres_init(tEQ[i].z, n); mpres_set(tEQ[i].z, tEP[i].z, n);
      mpres_init(num[i], n);
      mpres_init(den[i], n);
      mpres_init(inv[i], n);
  }
  
  last_chkpnt_time = cputime ();

#if DEBUG_MANY_EC >= 2
  printf("Initial points:\n");
  pt_many_print(tEP, nEP, n);
#endif
  for (r = 2.0; r <= B1; r *= 2.0)
      if (r > *B1done)
	{
	  if(pt_many_duplicate (tEQ, tEQ, nEP, n, num, den, inv, ok) == 0)
	    {
		mpz_set(f, tEQ[0].x);
		ret = ECM_FACTOR_FOUND_STEP1;
		goto end_of_all;
	    }
#if DEBUG_MANY_EC >= 1
	  printf("P%ld:=", (long)r); pt_many_print(tEQ, nEP, n); printf(";\n");
#endif
	}
  
  last_chkpnt_p = 3.;
  for (p = getprime (); p <= B1; p = getprime ())
    {
      for (r = p; r <= B1; r *= p)
	  printf("## p = %ld\n", (long)p);
	if (r > *B1done)
	  {
	    mpz_set_ui(e, (ecm_uint) p);
	    if(pt_many_mul(tEQ, tEP, nEP, e, n, num, den, inv, ok) == 0)
	      {
		mpz_set(f, tEQ[0].x);
		ret = ECM_FACTOR_FOUND_STEP1;
		goto end_of_all;
	      }
	    /* TODO: update ok if needed */
#if DEBUG_MANY_EC >= 2
	    pt_many_print(tEQ, nEP, n);
#endif
	  }

      if (stop_asap != NULL && (*stop_asap) ())
        {
          outputf (OUTPUT_NORMAL, "Interrupted at prime %.0f\n", p);
          break;
        }

#if 0
      /* WARNING: not activated yet */
      if (chkfilename != NULL && p > last_chkpnt_p + 10000. && 
          elltime (last_chkpnt_time, cputime ()) > CHKPNT_PERIOD)
        {
	  writechkfile (chkfilename, ECM_ECM, MAX(p, *B1done), n, A, x, y, z);
          last_chkpnt_p = p;
          last_chkpnt_time = cputime ();
        }
#endif
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

  /* normalize all points */

  mpz_clear(e);
  free(ok);
  for(i = 0; i < nEP; i++){
      mpres_clear(tEQ[i].A, n);
      mpres_clear(tEQ[i].x, n);
      mpres_clear(tEQ[i].y, n);
      mpres_clear(tEQ[i].z, n);
      mpres_clear(num[i], n);
      mpres_clear(den[i], n);
      mpres_clear(inv[i], n);
  }
  return ret;
}
#endif

int
read_and_prepare(mpz_t f, mpres_t x, mpq_t q, char *buf, mpz_t n)
{
    mpq_set_str(q, buf, 10);
    if(mod_from_rat(x, q, n) == 0){
	mpz_set(f, x);
	return 0;
    }
    return 1;
}

int
process_many_curves(mpz_t f, mpz_t n, double B1, char *fic_EP)
{
    curve tEP[NCURVE_MAX];
    FILE *ifile = fopen(fic_EP, "r");
    int nc = 0, i, ret = 0;
    char bufA[1024], bufx[1024], bufy[1024];
    ecm_params params;
    mpq_t q;
#ifndef ONE_CURVE_AT_A_TIME
    mpmod_t modulus;
#endif
    double B1done;

    mpq_init(q);
    while(fscanf(ifile, "%s %s %s", bufA, bufx, bufy) != EOF){
	mpz_init(tEP[nc].A);
	if(read_and_prepare(f, tEP[nc].A, q, bufA, n) == 0)
	    goto process_end;
	mpz_init(tEP[nc].x);
	if(read_and_prepare(f, tEP[nc].x, q, bufx, n) == 0)
	    goto process_end;
	mpz_init(tEP[nc].y);
	if(read_and_prepare(f, tEP[nc].y, q, bufy, n) == 0)
	    goto process_end;
	mpz_init_set_ui(tEP[nc].z, 1);
	nc++;
    }
#ifdef ONE_CURVE_AT_A_TIME
    ret = one_curve_at_a_time(f, tEP, nc, n, B1);
#else
    mpmod_init(modulus, n, ECM_MOD_DEFAULT);
    for(i = 0; i < nc; i++){
	mpres_set_z(tEP[i].x, tEP[i].x, modulus);
	mpres_set_z(tEP[i].y, tEP[i].y, modulus);
	mpres_set_z(tEP[i].z, tEP[i].z, modulus);
	mpres_set_z(tEP[i].A, tEP[i].A, modulus);
    }
    B1done = 1.0;
    ret = all_curves_at_once(f, tEP, nc, modulus, B1, &B1done, NULL, NULL);
#endif
    for(i = 0; i < nc; i++){
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

int
main (int argc, char *argv[])
{
  mpz_t n, f;
  int res;
  double B1;

  if (argc != 4)
    {
      fprintf (stderr, "Usage: ecmfactor <number> <fic_EP> <B1>\n");
      exit (1);
    }

  mpz_init (n);

  /* read number on command line */
  if (mpz_set_str (n, argv[1], 10))
    {
      fprintf (stderr, "Invalid number: %s\n", argv[1]);
      exit (1);
    }

  B1 = atof (argv[3]);

  mpz_init (f); /* for potential factor */

  printf ("Using all curves from %s with B1=%1.0f\n", argv[2], B1);
  res = process_many_curves(f, n, B1, argv[2]);

  if (res > 0)
    {
      printf ("found factor in step %u: ", res);
      mpz_out_str (stdout, 10, f);
      printf ("\n");
    }
  else if (res == ECM_NO_FACTOR_FOUND)
    printf ("found no factor\n");
  else
    printf ("error\n");

  mpz_clear (f);
  mpz_clear (n);

  return 0;
}
