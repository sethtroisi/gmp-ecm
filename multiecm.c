/* multiecm.c - ECM with many curves with many torsion and/or in parallel 
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
#include "getprime_r.h"

#ifdef HAVE_ADDLAWS
#include "addlaws.h"
#endif

#ifdef HAVE_TORSION
#include "torsions.h"
#endif

#define DEBUG_MULTI_EC 0
#define MULTI_USE_ADD_SUB 1

#define NCURVE_MAX 2000

/* Morain/Olivos 
   OUTPUT: 1 if ok, 0 if pb (bad state reached).
   SIDE EFFECT: fills in T[0..IT[
*/
int
MO_automaton(signed char *T, size_t *IT, mpz_t e, size_t le)
{
    size_t ie, iT = 0;
    int state = 0, res = 1, bz;

    for(ie = 0; ie < le; ie++){
	bz = ecm_tstbit(e, ie) == 0;
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
MO_check(signed char *T, size_t iT, mpz_t e)
{
    mpz_t tmp;
    int i, ok;

    printf("# Entering MO_check...\n");
    mpz_init_set_ui(tmp, 0);
    for(i = ((int)iT)-1; i >= 0; i--){
	mpz_mul_2exp(tmp, tmp, 1);
	mpz_add_si(tmp, tmp, (int)T[i]);
    }
#if DEBUG_ADD_LAWS >= 2
    gmp_printf("e:=%Zd;\n", e);
    gmp_printf("t:=%Zd;\n", tmp);
#endif
    ok = mpz_cmp(tmp, e) == 0;
    mpz_clear(tmp);
    return ok;
}

/* Do we have eval(T) == e? */
int
Split_check(short *S, size_t iS, mpz_t e)
{
    mpz_t tmp;
    size_t i;
    int ok;

    mpz_init_set_ui(tmp, 0);
    for(i = 0; i < iS; i += 2){
	mpz_add_si(tmp, tmp, (int)S[i+1]);
	mpz_mul_2exp(tmp, tmp, (int)S[i]);
    }
#if DEBUG_ADD_LAWS >= 2
    gmp_printf("e:=%Zd;\n", e);
    gmp_printf("t:=%Zd;\n", tmp);
#endif
    ok = mpz_cmp(tmp, e) == 0;
    mpz_clear(tmp);
    return ok;
}

/* Adapted from Koyama and Tsuroka, CRYPTO'92, using less space. 
   S is filled in left-to-right from T.
*/
size_t
Split(short *S, size_t Slen, signed char *T, size_t iT, int w0)
{
    size_t i = iT-1, j, k, lW, iS = 0, w = (size_t)w0;
    int gap;
    short W, nomorei = 0;

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
	    W += T[k]*(((short)1<<(k-j)));
#if 0
	i = i-w; // new i >= -1
#endif
	if(i >= w){
	    i = i-w;
	    /* exclude left zeros and update power of 2 */
#if 0
	    while((i >= 0) && (T[i] == 0)){
#endif
	    while(T[i] == 0){
		gap++;
		if(i == 0){
		    nomorei = 1;
		    break;
		}
		i--;
	    }
	}
	else{
	    // case i = w-1
	    nomorei = 1;
	}
	S[iS] = gap;
	S[iS+1] = W;
	if(iS >= 2)
	    S[iS-2] += lW;
	iS += 2;
	if(iS > Slen)
	    return 0;
	if(nomorei)
	    break;
    }
    /* at this point, we have to examine T[0..i] */
    if(nomorei == 0){
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
	    W += T[k]*(((short)1) << (k-j));
	S[iS] = gap;
	S[iS+1] = W;
	if(iS >= 2)
	    S[iS-2] += lW;
	iS += 2;
	if(iS > Slen)
	    return 0;
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
size_t
build_MO_chain(short *S, size_t Slen, mpz_t e, int w)
{
    /* first use automata */
    size_t le = mpz_sizeinbase(e, 2), iT = 0, iS = 0;
#if DEBUG_ADD_LAWS >= 0
    long tp = cputime();
    size_t i;
#endif
    signed char *T = (signed char *)malloc((2*le) * sizeof(signed char)); /* humf */

    MO_automaton(T, &iT, e, le);
    if(iT > (2*le)){
	printf("#!# Error: iT too large %lu\n", iT);
    }
#if DEBUG_ADD_LAWS >= 2
    /* check value of T */
    gmp_printf("# signed digits(%Zd):", e);
    for(i = 0; i < (int)iT; i++)
	printf(" %d", T[i]);
    printf("\n");
#endif
#if DEBUG_ADD_LAWS >= 2
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
#endif
#if DEBUG_ADD_LAWS >= 2
    if(Split_check(S, iS, e) == 0){
        printf("#!# Error in Split\n");
	exit(-1);
    }
    else
        printf("# good check in Split\n");

#endif
    free(T);
    return iS;
}

size_t
build_add_sub_chain(short *S, size_t Slen, mpz_t e, int w)
{
    return build_MO_chain(S, Slen, e, w);
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
	return 7;
    else
	return 7;
}

/* pack everybody */
void
add_sub_pack(mpz_t s, int w, short *S, size_t iS)
{
    size_t nsh, cte = sizeof(mp_limb_t)/sizeof(short);
    size_t cte2 = ((size_t) 1) << 16;
    unsigned short *tmp;

    nsh = iS / cte;
    if(iS % cte != 0)
	nsh++;
    nsh *= cte;
    nsh += 4;
    /* coding */
    tmp = (unsigned short *)malloc(nsh * sizeof(unsigned short));
    tmp[0] = w;
    tmp[1] = iS / cte2;
    tmp[2] = iS % cte2;
    memcpy(tmp+4, S, iS * sizeof(unsigned short));
    s->_mp_d = (mp_limb_t *)tmp; /* humf */
}

void add_sub_unpack(int *w, short **S, size_t *iS, mpz_t s)
{
    unsigned short *T;

    T = (unsigned short *)s->_mp_d; /* humf */
    *w = (int)T[0];
    *iS = (size_t)((((size_t)T[1]) << 16) + (size_t)T[2]);
    *S = (short *)(T+4);
#if DEBUG_ADD_LAWS >= 2
    printf("# iS_unpack = %lu\n", *iS);
#endif
}

/* INPUT: S = [[ts, 2*ds+1], ..., [t1, 2*d1+1], [t0, 2*d0+1]] for
   e = 2^t0 * (2*d0+1 + 2^t1 *(2*d1+1 + 2^t2 * (2*d2+1+... + 2^ts*(2*ds+1) ),
   with -2^(w-1)+1 <= 2*di+1 < 2^{w-1}.
*/
int
ell_point_mul_add_sub_with_S(mpz_t f, ell_point_t Q, ell_point_t P, 
			     ell_curve_t E, mpmod_t n, int w, short *S, int iS)
{
    ell_point_t P0;
    ell_point_t iP[EC_ADD_SUB_2_WMAX];
    size_t j;
    int status = 1, i, k;
#if DEBUG_ADD_LAWS >= 2
    mpz_t ex;
    unsigned long cpt = 0;
    long tp;
    ell_point_t QQ;

    tp = cputime();
    /* to reconstruct the exponent to check mults */
    mpz_init_set_ui(ex, 0);
    ell_point_init(QQ, E, n);
#endif

    /* iP[i] <- (2*i+1) * P */
    k = (1 << (w-1)) - 1;
    for(i = 0; i <= k; i++)
	ell_point_init(iP[i], E, n);
    ell_point_set(iP[0], P, E, n);
#if DEBUG_ADD_LAWS >= 2
    ell_point_check(iP[0], E, n);
    gmp_printf("N:=%Zd;\n", n->orig_modulus);
    printf("E:="); ell_curve_print(E, n); printf(";\n");
    printf("P:="); ell_point_print(P, E, n); printf(";\n");
#endif
    if(k > 0){
	/* P[k] <- [2]*P */
	if(ell_point_duplicate(f, iP[k], P, E, n) == 0){
	    mpres_set(P0->x, iP[k]->x, n);
            status = 0;
	    goto ell_point_mul_add_sub_end;
        }
#if DEBUG_ADD_LAWS >= 2
	ell_point_check(iP[k], E, n);
#endif
	for(i = 1; i <= k; i++){
	    if(ell_point_add(f, iP[i], iP[i-1], iP[k], E, n) == 0){
                status = 0;
		goto ell_point_mul_add_sub_end;
	    }
#if DEBUG_ADD_LAWS >= 2
	    ell_point_check(iP[i], E, n);
#endif
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
	printf("P0:="); ell_point_print(P0, E, n); printf(":\n");
#endif
	i = abs(S[j+1]) >> 1; /* (abs(S[j+1])-1)/2, S[j+1] is always odd */
	assert(i <= k);
	if(S[j+1] > 0){
	    if(ell_point_add(f, P0, P0, iP[i], E, n) == 0){
		status = 0;
		break;
	    }
#if DEBUG_ADD_LAWS >= 2
	    printf("iP%d:=", i); ell_point_print(iP[i], E, n); printf(":\n");
	    printf("Radd:="); ell_point_print(P0, E, n); printf(":\n");
	    printf("Q:=ProjEcmAdd(P0, iP%d, E, N): ", i);
	    printf("printf(\"CHK%lu: %%a\\n\", ProjEcmEqual(Q, Radd, N));\n",
		   cpt++);
#endif
	}
	else{
	    /* add(-P) = sub(P) */
	    if(ell_point_sub(f, P0, P0, iP[i], E, n) == 0){
		status = 0;
		break;
	    }
#if DEBUG_ADD_LAWS >= 2
	    printf("iP%d:=", i); ell_point_print(iP[i], E, n); printf(":\n");
	    printf("Rsub:="); ell_point_print(P0, E, n); printf(":\n");
	    printf("S:=ProjEcmNegate(iP%d, N):\n", i);
	    printf("Q:=ProjEcmAdd(P0, S, E, N): ");
	    printf("printf(\"CHK%lu: %%a\\n\", ProjEcmEqual(Q, Rsub, N));\n",
		   cpt++);
#endif
	}
#if DEBUG_ADD_LAWS >= 2
	ell_point_check(P0, E, n);
#endif
	/* now multiply */
	for(i = 0; i < S[j]; i++){
	    if(ell_point_duplicate(f, P0, P0, E, n) == 0){
		status = 0;
		break;
	    }
#if DEBUG_ADD_LAWS >= 2
	    printf("Rdup:="); ell_point_print(P0, E, n); printf(":\n");
	    printf("Q:=ProjEcmDouble(P0, E, N): ");
	    printf("printf(\"CHK%lu: %%a\\n\", ProjEcmEqual(Q, Rdup, N));\n",
		   cpt++);
#endif
	}
#if DEBUG_ADD_LAWS >= 2
	mpz_add_si(ex, ex, (int)S[j+1]);
	mpz_mul_2exp(ex, ex, (int)S[j]);
#endif
#if DEBUG_ADD_LAWS >= 2
	ell_point_check(P0, E, n);
	ell_point_mul_plain(QQ, ex, P, E, n);
	if(ell_point_equal(QQ, P0, E, n) == 0){
	    gmp_printf("ex:=%Zd;\n", ex);
	    printf("P0:="); ell_point_print(P0, E, n); printf(";\n");
	    printf("QQ:="); ell_point_print(QQ, E, n); printf(";\n");
	    exit(-1);
	}
#endif
	if(status == 0)
	    break;
    }
 ell_point_mul_add_sub_end:
#if DEBUG_ADD_LAWS >= 2
    printf("# time[addsub%d] = %ldms\n", w, elltime(tp, cputime()));
    printf("Checking with [ex]*P with %ld bits\n", (long)mpz_sizeinbase(ex,2));
    tp = cputime();
    ell_point_mul_plain(QQ, ex, P, E, n);
    if(ell_point_equal(QQ, P0, E, n) == 0){
	gmp_printf("ex:=%Zd;\n", ex);
	printf("P0:="); ell_point_print(P0, E, n); printf(";\n");
	printf("QQ:="); ell_point_print(QQ, E, n); printf(";\n");
	exit(-1);
    }
    printf("# time[plain] = %ldms\n", elltime(tp, cputime()));
    mpz_clear(ex);
#endif
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
ell_point_mul_add_sub(mpz_t f, ell_point_t Q, mpz_t e, ell_point_t P,
		      ell_curve_t E, mpmod_t n)
{
    size_t iS = 0, Slen, w;
    int negated = 0, status = 1;
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
#if DEBUG_ADD_LAWS >= 2
    printf("# Slen=%lu\n", Slen);
#endif
    S = (short *)malloc(Slen * sizeof(short));
    iS = build_add_sub_chain(S, Slen, e, w);
    if(iS == 0){
	printf("build_NAF: Slen=%"PRIu64" too small\n", Slen);
	return -1;
    }
#if DEBUG_ADD_LAWS >= 2
    gmp_printf("addsub[%Zd=>%d]:", e, iS);
    for(j = iS-1; j >= 0; j--)
	printf(" %d", S[j]);
    printf("\n");
    printf("P:="); ell_point_print(P, E, n); printf(";\n");
#endif
    status = ell_point_mul_add_sub_with_S(f, Q, P, E, n, w, S, iS);
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
		mpz_gcd(f, Q->z, n->orig_modulus); // FIXME: why?
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
compute_s_4_add_sub(mpz_t s, ecm_uint B1, int disc)
{
    mpz_t t;
    long tp;
    short *S;
    size_t Slen, iS;
    int w, *forbiddenres = compute_forbidden_res(disc);

    mpz_init(t);
    tp = cputime();
    compute_s(t, B1, forbiddenres);
    free(forbiddenres);
    printf("# computing prod(p^e <= %lu): %ldms\n", B1, elltime(tp,cputime()));
#if USE_ADD_SUB_CHAINS == 0 /* keeping it simple for the time being */
    mpz_set(s, t);
#else
    tp = cputime();
    w = get_add_sub_w(t);
    /* Slen = 2 * log_{2^w}(t) = 2*log_2(t)/w = 2 * 64 * size(t)/w */
    Slen = (2 * GMP_NUMB_BITS * mpz_size(t)) / w;
    S = (short *)malloc(Slen * sizeof(short));
    iS = build_add_sub_chain(S, Slen, t, w);
    printf("# NAF has %"PRIu64" terms (w=%d, Slen=%"PRIu64"): %ldms\n", iS, w, Slen,
	   elltime(tp,cputime()));
    if(iS == 0){
	printf("build_NAF: Slen=%"PRIu64" too small\n", Slen);
	return 0;
    }
    add_sub_pack(s, w, S, iS);
    free(S);
#endif
    mpz_clear(t);
    return 1;
}

/* fall back on traditional ECM.
   TODO: use chkfile also.
 */
int
process_one_curve(mpz_t f, mpz_t N, double B1, mpz_t B2,
		  ecm_params params, ell_curve_t E, ell_point_t P)
{
    int ret;

    /* if B2 = ECM_DEFAULT_B2, compute it automatically from B1: 
       no freedom on B2! */
    mpz_set(params->B2, B2);
    /* will be set to B1 */
    mpz_set_si(params->B2min, ECM_DEFAULT_B2);

    mpz_set(params->x, P->x);
    mpz_set(params->sigma, E->a4); /* humf */

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
        printf ("Found %s factor of %u digits: ", 
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

#if DEBUG_MULTI_EC >= 2
/* f is a (probable) prime factor of n. tP is in plain mod n form. */
void
dump_curves(ell_curve_t *tE, ell_point_t *tP, int nE, mpz_t f)
{
    int i;

    printf("CheckE:=procedure(E, D, P, info)\n");
    printf("    K:=QuadraticField(D); OK:=MaximalOrder(K);\n");
    printf("    printf \"#E[%%o]=%%o\\n\", info, Factorization(#E);\n");
    printf("    tw:=Twists(E); Et:=tw[2];\n");
    printf("    printf \"#Et[%%o]=%%o\\n\", info, Factorization(#Et);\n");
    printf("    gen:=Generators(E); printf \"ords=%%o\\n\", ");
    printf("[Factorization(Order(g)):g in gen];\n");
    printf("    lf:=Factorization(Order(E!P)); printf \"ord(P)=%%o\\n\", lf;\n");
    printf("    for i:=1 to #lf do\n");
    printf("        lfi:=Factorization(lf[i][1]*OK);\n");
    printf("        ok,gen:=IsPrincipal(lfi[1][1]); print lf[i], ok, gen;\n");
    printf("    end for;\n");
    printf("end procedure;\n");
    gmp_printf("p:=%Zd; F:=GF(p); P:=[]; A:=[]; B:=[]; E:=[]; D:=[];\n", f);
    for(i = 0; i < nE; i++){
	printf("D[%d]:=%d;\n", i+1, tE[i]->disc);
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
	    mpres_set_z(A, tE[i]->a4, fmod);
	    if(montgomery_to_weierstrass(tmp, x, y, A, fmod) 
	       == ECM_FACTOR_FOUND_STEP1){
		printf("GASP while dumping a Montgomery form curve!\n");
	    }
	    printf("P[%d]:=[", i+1); print_mpz_from_mpres(x, fmod);
	    printf(", "); print_mpz_from_mpres(y,fmod);
	    printf(", 1];\n");
	    printf("A[%d]:=", i+1);
	    print_mpz_from_mpres(A, fmod);
	    printf(";\n");
	    printf("B[%d]:=(P[%d][2]^2-P[%d][1]^3-A*P[%d][1]) mod N;\n", 
		   i+1, i+1, i+1, i+1);
	    mpres_clear(x, fmod);
	    mpres_clear(y, fmod);
	    mpres_clear(A, fmod);
	    mpmod_clear(fmod);
	    mpz_clear(tmp);
	}
	else if(tE[i]->type == ECM_EC_TYPE_WEIERSTRASS){
	    gmp_printf("P[%d]:=[%Zd, %Zd, %Zd];\n", i+1, 
		       tP[i]->x, tP[i]->y, tP[i]->z); 
	    gmp_printf("A[%d]:=%Zd;\n", i+1, tE[i]->a4);
	    gmp_printf("B[%d]:=%Zd;\n", i+1, tE[i]->a6);
	}
	else{
	    printf("Case %d NYI in dump_curves\n", tE[i]->type);
	    break;
	}
	printf("E[%d]:=EllipticCurve([F!A[%d], F!B[%d]]);\n", i+1, i+1, i+1);
	printf("CheckE(E[%d], D[%d], P[%d], infos[%d]);\n",i+1,i+1,i+1,i+1);
    }
}
#endif /* DEBUG_MULTI_EC >= 2 */

/* TODO: better control of B2 + dichotomy (cf. #B2) */
int
one_curve_at_a_time(mpz_t f, char *ok, ell_curve_t *tE, ell_point_t *tP, int nE,
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
#if DEBUG_MULTI_EC >= 2
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
#if DEBUG_MULTI_EC >= 2
		if(ret == ECM_PRIME_FAC_PRIME_COFAC 
		   || ret == ECM_PRIME_FAC_COMP_COFAC){
		    /* output Magma lines to check #E's mod f */
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
#if DEBUG_MULTI_EC >= 2
    printf("# let's debug all curves\n");
    dump_curves(tE, tP, nE, N);
#endif
    mpz_clear (C);
    mpcandi_t_free(&candi);
    return ret;
}

/**********************************************************************
 Using parallelism.
**********************************************************************/

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
pt_print(ell_curve_t E, ell_point_t P, mpmod_t n)
{
    printf("[");
    print_mpz_from_mpres(P->x, n);
    printf(", ");
    print_mpz_from_mpres(P->y, n);
    printf(", ");
    if(E->type == ECM_EC_TYPE_WEIERSTRASS && E->law == ECM_LAW_AFFINE)
	gmp_printf("%Zd", P->z);
    else
	print_mpz_from_mpres(P->z, n);
    printf("]");
}

void
pt_many_print(ell_curve_t *tE, ell_point_t *tP, int nE, mpmod_t n)
{
    int i;

    for(i = 0; i < nE; i++){
	printf("%d: ", i);
	pt_print(tE[i], tP[i], n);
	printf(" on E.A=");
	print_mpz_from_mpres(tE[i]->a4, n);
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
compute_all_inverses(mpz_t f, mpres_t *inv, mpres_t *x, int nx, mpmod_t n, char *takeit)
{
    int i;

#if 0
    /* plain version, to debug the architecture */
    for(i = 0; i < nx; i++){
	if(takeit[i] != 1)
	    continue;
	if(!mpres_invert(inv[i], x[i], n)){
	    mpres_gcd(inv[0], x[i], n); // FIXME!!
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
	mpres_gcd(f, inv[nx-1], n);
#if DEBUG_ADD_LAWS >= 1
	printf("Factor[%d]: ", i);
	mpz_out_str (stdout, 10, f);
	printf ("\n");
#endif
	/* identifying the x[i]'s */
	for(i = 0; i < nx; i++){
	    mpres_gcd(f, x[i], n);
	    if(mpz_cmp_ui(f, 1) != 0){
#if DEBUG_ADD_LAWS >= 0
		printf("# x[%d] not invertible: ", i);
		mpz_out_str (stdout, 10, f);
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
pt_many_common(mpz_t f, ell_point_t *tR, ell_point_t *tP, ell_point_t *tQ, 
	       int nE, mpmod_t n, mpres_t *num, mpres_t *den, mpres_t *inv,
	       char *takeit)
{
    int i;

    if(compute_all_inverses(f, inv, den, nE, n, takeit) == 0){
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
pt_many_duplicate(mpz_t f, ell_point_t *tQ, ell_point_t *tP, ell_curve_t *tE,
		  int nE, mpmod_t n, 
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
	    mpres_add(num[i], num[i], tE[i]->a4, n);
	    mpres_mul_ui(den[i], tP[i]->y, 2, n);
	}
    }
    res = pt_many_common(f, tQ, tP, tP, nE, n, num, den, inv, takeit);
    /* TODO: case takeit[i] == 2 */
    free(takeit);
    return res;
}

/* R[i] <- P[i] + Q[i], or a factor is found which is put in num[nE]. */
int
pt_many_add(mpz_t f, ell_point_t *tR, ell_point_t *tP, ell_point_t *tQ, ell_curve_t *tE, 
	    int nE, mpmod_t n, 
	    mpres_t *num, mpres_t *den, mpres_t *inv, char *ok)
{
    char *takeit = (char *)malloc(nE * sizeof(char));
    int i, res;

    memcpy(takeit, ok, nE);
#if DEBUG_ADD_LAWS >= 2
    printf("In pt_many_add, adding\n");
    pt_many_print(tE, tP, nE, n);
    printf("and\n");
    pt_many_print(tE, tQ, nE, n);
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
		mpres_add(num[i], num[i], tE[i]->a4, n);
		mpres_mul_ui(den[i], tP[i]->y, 2, n);
	    }
	}
	else if(mpres_equal(tQ[i]->x, tP[i]->x, n)){
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
    res = pt_many_common(f, tR, tP, tQ, nE, n, num, den, inv, takeit);
    /* TODO: case takeit[i] == 2 */
    free(takeit);
    return res;
}

/* tER != tEP */
static int
pt_many_sub(mpz_t f, ell_point_t *tR, ell_point_t *tQ, ell_point_t *tP, ell_curve_t *tE,
	    int nE, mpmod_t n, 
	    mpres_t *num, mpres_t *den, mpres_t *inv, char *ok)
{
    int i, res;

    for(i = 0; i < nE; i++)
	if(ok[i] == 1)
	    pt_neg(tP[i], n);
    res = pt_many_add(f, tR, tQ, tP, tE, nE, n, num, den, inv, ok);
    for(i = 0; i < nE; i++)
	if(ok[i] == 1)
	    pt_neg(tP[i], n);
    return res;
}

/* Ordinary binary left-right addition */
static int
pt_many_mul_plain(mpz_t f, ell_point_t *tQ, ell_point_t *tP, ell_curve_t *tE,
		  int nE, mpz_t e, mpmod_t n, 
		  mpres_t *num, mpres_t *den, mpres_t *inv, char *ok)
{
  size_t l = mpz_sizeinbase (e, 2) - 1; /* l >= 1 */
  int status = 1;

  pt_many_assign(tQ, tP, nE, n);
  while (l-- > 0)
    {
	if(pt_many_duplicate (f, tQ, tQ, tE, nE, n, num, den, inv, ok) == 0)
	  {
	    status = 0;
	    break;
	  }
#if DEBUG_ADD_LAWS >= 2
	printf("Rdup:="); pt_many_print(tE, tQ, nE, n); printf(";\n");
#endif
	if (ecm_tstbit (e, l))
	  {
	      if(pt_many_add (f, tQ, tP, tQ, tE, nE, n, num, den, inv, ok) == 0)
	      {
		status = 0;
		break;
	      }
#if DEBUG_ADD_LAWS >= 2
	      printf("Radd:="); pt_many_print(tE, tQ, nE, n); printf(";\n");
#endif
	  }
    }
  return status;
}

/* Ordinary binary left-right addition; see Solinas00. Morally, we use
 w = 2. */
static int
pt_many_mul_add_sub_si(mpz_t f, ell_point_t *tQ, ell_point_t *tP, ell_curve_t *tE, int nE,
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
	if(pt_many_duplicate(f, tQ, tQ, tE, nE, n, num, den, inv, ok) == 0){
	    status = 0;
	    break;
	}
#if DEBUG_ADD_LAWS >= 2
	printf("Rdup:="); pt_many_print(tE, tQ, nE, n); printf(";\n");
#endif
	if(S[j] == 1){
	    if(pt_many_add(f, tQ, tQ, tP, tE, nE, n, num, den, inv, ok) == 0){
		status = 0;
		break;
	    }
#if DEBUG_ADD_LAWS >= 2
	    printf("Radd:="); pt_many_print(tE, tQ, nE, n); printf(";\n");
#endif
	}
	else if(S[j] == -1){
	    if(pt_many_sub(f, tQ, tQ, tP, tE, nE, n, num, den, inv, ok) == 0){
		status = 0;
		break;
	    }
#if DEBUG_ADD_LAWS >= 2
	    printf("Rsub:="); pt_many_print(tE, tQ, nE, n); printf(";\n");
#endif
	}
    }
  return status;
}

/* tEQ[i] <- e * tEP[i]; we must have tEQ != tEP */
/* If a factor is found, it is put back in num[nE]. */
int
pt_many_mul(mpz_t f, ell_point_t *tQ, ell_point_t *tP, ell_curve_t *tE, int nE,
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
      status = pt_many_mul_add_sub_si(f, tQ, tP, tE, nE, mpz_get_si(e), n,
				      num, den, inv, ok);
  else
      status = pt_many_mul_plain(f, tQ, tP, tE, nE, e, n, num, den, inv, ok);


pt_many_mul_end:

  /* Undo negation to avoid changing the caller's e value */
  if (negated){
    mpz_neg (e, e);
    pt_many_neg(tP, nE, n);
  }
  return status;
}

/* Copied from classical ecm_stage1. */
int
all_curves_at_once(mpz_t f, char *ok, ell_curve_t *tE, ell_point_t *tP, int nE,
		   mpmod_t n, double B1, double *B1done, 
		   int (*stop_asap)(void), char *chkfilename)
{
    ell_point_t tQ[NCURVE_MAX], tR[NCURVE_MAX];
    mpz_t num[NCURVE_MAX+1], den[NCURVE_MAX+1], inv[NCURVE_MAX], e;
    double p = 0.0, r, last_chkpnt_p;
    int ret = ECM_NO_FACTOR_FOUND;
    long last_chkpnt_time;
    int i;
    prime_info_t prime_info;
    
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
    
#if DEBUG_MULTI_EC >= 2
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
#if DEBUG_MULTI_EC >= 2
	    printf("P%ld:=", (long)r); pt_many_print(tQ, nE, n); printf(";\n");
#endif
	}

    last_chkpnt_p = 3.;
    prime_info_init (prime_info);
    for (p = getprime_mt (prime_info); p <= B1; p = getprime_mt (prime_info)){
	for (r = p; r <= B1; r *= p){
#if DEBUG_MULTI_EC >= 2
	    printf("## p = %ld at %ldms\n", (long)p, cputime());
#endif
	    if (r > *B1done){
		mpz_set_ui(e, (ecm_uint) p);
		if(pt_many_mul(tR, tQ, tE, nE, e, n, num, den, inv, ok) == 0){
		    mpz_set(f, num[nE]);
		    ret = ECM_FACTOR_FOUND_STEP1;
		    goto end_of_all;
		}
#if DEBUG_MULTI_EC >= 2
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
    prime_info_clear (prime_info); /* free the prime table */
    
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
read_and_prepare(mpz_t f ATTRIBUTE_UNUSED, mpz_t x ATTRIBUTE_UNUSED, mpq_t q,
                 char *buf, mpz_t n ATTRIBUTE_UNUSED)
{
    mpq_set_str(q, buf, 10);
#ifdef HAVE_TORSION
    if(mod_from_rat(x, q, n) == 0){
	mpz_set(f, x);
	return 0;
    }
#endif
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
		    ell_curve_t *tE, ell_point_t *tP, int nE, 
		    ecm_params params, int onebyone, char *savefilename)
{
    double B1done;
    ell_point_t tQ[NCURVE_MAX];
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
	ell_point_init(tQ[i], tE[i], n);
	ell_point_set(tQ[i], tP[i], tE[i], n);
    }
    B1done = 1.0;
    ret = all_curves_at_once(f, ok, tE, tQ, nE, n, B1, &B1done, NULL, NULL);
    printf("# Step 1 took %ldms\n", elltime (st, cputime ()));

    if(ret != ECM_NO_FACTOR_FOUND){
	ret = conclude_on_factor(n->orig_modulus, f, params->verbose);
#if DEBUG_MULTI_EC >= 2
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
#if DEBUG_MULTI_EC >= 1
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
	ell_point_clear(tQ[i], tE[i], n);
    free(ok);
    return ret;
}

int
read_curves_from_file(int *nE, ell_curve_t *tE, ell_point_t *tP, 
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
	ell_curve_init(tE[*nE],ECM_EC_TYPE_WEIERSTRASS,ECM_LAW_HOMOGENEOUS,n);
	if(Etype == 'W'){
	    if(fscanf(ifile, "%s %s %s", bufA, bufx, bufy) == EOF)
		break;
	    tE[*nE]->type = ECM_EC_TYPE_WEIERSTRASS;
	    tE[*nE]->law = ECM_LAW_AFFINE;
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
	mpz_init(tE[*nE]->a4);
	if(read_and_prepare(tf[*nf], tE[*nE]->a4, q, bufA, n->orig_modulus) == 0){
	    ret = 0;
	    *nf += 1;
	    goto process_end;
	}
	ell_point_init(tP[*nE], tE[*nE], n);
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
    ell_curve_t tE[NCURVE_MAX];
    ell_point_t tP[NCURVE_MAX];
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
#ifdef HAVE_TORSION
	else if(torsion != NULL)
	    ret = build_curves_with_torsion(tf[*nf],modulus,tE,tP,
					    torsion,smin,smax,nE,disc,sqroots);
#endif
	else if(disc != 0){
#if 0
	    ret = build_curves_with_CM(tf[*nf],&nE,tE,tP,disc,modulus,sqroots);
#else
	    printf("Sorry, disabled right now!\n");
	    exit(-1);
#endif
	}
	if(ret == ECM_NO_FACTOR_FOUND)
	    ret = process_many_curves(tf[*nf],modulus,B1,B2,tE,tP,nE,params,
				      onebyone,savefilename);
	else{
	    printf("Quid? %d\n", ret);
	    break;
	}
	/* clear curves */
	for(i = 0; i < nE; i++){
	    ell_point_clear(tP[i], tE[i], modulus);
	    ell_curve_clear(tE[i], modulus);
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

/* Assume b^n = 1 mod N.
   status =  0 if the squareroot could not be computed,
            -1 if the check is bad,
	     1 otherwise.
 */
int
odd_square_root_mod_N(mpz_t f, int *status, mpz_t *sqroots,
		      int b, int n, int q, mpz_t N)
{
    mpz_t zeta, tmp, tmp2;
    int np = n, e = 0, *tab, x, ret = ECM_NO_FACTOR_FOUND;

    *status = 1;
    while(np % q == 0){
	e++;
	np /= q;
    }
    /*printf("# n = %d = %d^%d * %d\n", n, q, e, np);*/
    mpz_init_set_ui(zeta, b);
    mpz_powm_ui(zeta, zeta, np, N);
    if(mpz_cmp_ui(zeta, 1) == 0){
	printf("# missed: zeta == 1\n");
	*status = 0;
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
	    *status = -1;
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
static void
psb_minus_even(int *tsq, mpz_t sqroots[], int b, int k, mpz_t N)
{
    int isq = 0;

    /*    printf("# got sqrt(-1)\n");*/
    tsq[isq++] = -1;
    mpz_init_set_si(sqroots[0], b);
    mpz_powm_ui(sqroots[0], sqroots[0], k, N);
    if(k % 2 == 0){
	/* zeta8 = b^(k/2) = (1+zeta4)/sqrt(2) 
	   => sqrt(2) = (1+zeta4)/zeta8 */
	/*	printf("# got sqrt(2)\n");*/
	mpz_init_set_si(sqroots[1], b);
	mpz_powm_ui(sqroots[1], sqroots[1], k>>1, N);
	mpz_invert(sqroots[1], sqroots[1], N);
	mpz_add_si(sqroots[0], sqroots[0], 1);
	mpz_mul(sqroots[1], sqroots[1], sqroots[0]);
	mpz_mod(sqroots[1], sqroots[1], N);
	mpz_sub_si(sqroots[0], sqroots[0], 1);
	tsq[isq++] = 2;
    }
    tsq[isq] = 0;
}

/* b^(2*k+1) = -1 => (b^(k+1))^2 = -b mod N */
static void
psb_minus_odd(int *tsq, mpz_t sqroots[], int b, int k, mpz_t N)
{
    /* printf("# got sqrt(-%d)\n", b);*/
    tsq[0] = -b;
    mpz_init_set_si(sqroots[0], b);
    mpz_powm_ui(sqroots[0], sqroots[0], k+1, N);
}

/* b^(2*k+1) = 1 mod N => (b^(k+1))^2 = b mod N */
static void
psb_plus_odd(int *tsq, mpz_t sqroots[], int b, int k, mpz_t N)
{
    /*    printf("# got sqrt(%d)\n", b);*/
    mpz_init_set_si(sqroots[0], b);
    mpz_powm_ui(sqroots[0], sqroots[0], k+1, N);
    tsq[0] = b;
}

/* N | b^n+c
   OUTPUT: ECM_NO_FACTOR_FOUND or ECM_FACTOR_FOUND_STEP1 in very rare cases! 
*/
static int
prepare_squareroots_from_powers(mpz_t f, int *tsq, mpz_t sqroots[], 
				int b, int n, int c, mpz_t N)
{
    int k, ret = ECM_NO_FACTOR_FOUND;
    mpz_t tmp, tmp2;

    tsq[0] = 0;
    tsq[1] = 0;
    if(c == -1){
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
	    if(k % 2 == 1)
		/* b^(2*r+1) = 1 mod N => (b^(r+1))^2 = b mod N */
		psb_plus_odd(tsq, sqroots, b, k>>1, N);
	    else{
		/* b^(2*r) = 1 */
		mpz_add_si(tmp2, tmp2, 1);
		if(mpz_cmp(tmp2, N) == 0){
		    /* case b^r = -1 */
		    printf("# %d^%d = -1 mod N;\n", b, k>>1);
		    if(k % 4 == 0)
			/* b^(2*s) = -1 */
			psb_minus_even(tsq, sqroots, b, k>>2, N);
		    else
			/* b^(2*s+1) = -1 */
			psb_minus_odd(tsq, sqroots, b, k>>2, N);
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
	else
	    /* b^(2*k+1) = 1 mod N => (b^(k+1))^2 = b mod N */
	    psb_plus_odd(tsq, sqroots, b, n>>1, N);
    }
    else if(c == 1){
	/* b^n = -1 mod N */
	if(n % 2 == 0)
	    /* b^(2*k) = -1 mod N => (b^k)^2 = -1 mod N */
	    psb_minus_even(tsq, sqroots, b, n>>1, N);
	else
	    /* b^(2*k+1) = -1 mod N => (b^(k+1))^2 = -b mod N */
	    psb_minus_odd(tsq, sqroots, b, n>>1, N);
    }
    else{
	/* b^n = -c mod N */
	if(n % 2 == 0){
	    /* (b^k)^2 = -c */
	    tsq[0] = -c; /* FIXME: case c non squarefree? */
	    mpz_init_set_si(sqroots[0], b);
	    mpz_powm_ui(sqroots[0], sqroots[0], n>>1, N);
	}
	else{
	    /* (b^(k+1))^2 = -c*b */
	    tsq[0] = -c*b; /* FIXME: case c non squarefree? */
	    mpz_init_set_si(sqroots[0], b);
	    mpz_powm_ui(sqroots[0], sqroots[0], (n+1) >>1, N);
	}
    }
    return ret;
}

/* N is a cofactor of b^n+c. */
static int
prepare_squareroots(mpz_t f, int *tsq, mpz_t sqroots[], 
		    int b, int n, int c, mpz_t N)
{
    int ret, tabq[] = {3, 5, 7, 11, 13, 19, 0}, q, iq, qs, isq, nn, status;

    tsq[0] = 0;
    ret = prepare_squareroots_from_powers(f,tsq,sqroots,b,n,c,N);
    if(ret != ECM_NO_FACTOR_FOUND)
	return ret;
    printf("# I already found squareroots for:");
    for(isq = 0; tsq[isq] != 0; isq++)
	printf(" %d", tsq[isq]);
    printf("\n");
    /* let's find some squareroots using small odd prime divisors of n */
    if(abs(c) == 1){
	for(iq = 0; tabq[iq] != 0; iq++){
	    q = tabq[iq];
	    qs = (q % 4 == 1 ? q : -q);
	    if(n % q == 0){
		/*	    printf("# I can find sqrt(%d)\n", qs);*/
		/* make sure that b^nn = 1 */
		nn = (c == -1 ? n : 2*n);
		ret = odd_square_root_mod_N(f,&status,sqroots+isq,b,nn,q,N);
		if(ret != ECM_NO_FACTOR_FOUND)
		    break;
		if(status == 1)
		    tsq[isq++] = qs;
	    }
	}
	tsq[isq] = 0;
    }
    return ret;
}

static int
rebuild_squareroot(mpz_t sq2[], int *tsq, mpz_t sqroots[], int *tqs, mpz_t N)
{
    mpz_t tmp;
    int isq, iqs, disc = tqs[0], ret, qs;
    
    mpz_set_ui(sq2[0], 1);
    for(iqs = 1; tqs[iqs] != 0; iqs++){
	for(isq = 0; tsq[isq] != 0; isq++){
	    if(tsq[isq] == -1)
		qs = -4;
	    else if(tsq[isq] == -2)
		qs = -8;
	    else
		qs = tsq[isq];
	    if(qs == tqs[iqs]){
		disc /= qs;
		mpz_mul(sq2[0], sq2[0], sqroots[isq]);
		mpz_mod(sq2[0],sq2[0], N);
		break;
	    }
	}
    }
    if(disc != 1){
	printf("#!# Pb: disc != 1: %d\n", disc);
	return 0;
    }
    /* check */
    disc = tqs[0];
    if(disc % 4 == 0) disc >>= 2;
    mpz_init_set(tmp, sq2[0]);
    mpz_mul(tmp, tmp, sq2[0]);
    mpz_sub_si(tmp, tmp, disc);
    mpz_mod(tmp, tmp, N);
    if(mpz_sgn(tmp) == 0){
	printf("# good check for sqrt(%d)\n", tqs[0]);
	ret = 1;
    }
    else{
	printf("# bad check for sqrt(%d)\n", tqs[0]);
	ret = 0;
    }
    mpz_clear(tmp);
    return ret;
}

/* Consider M = b^n+1 if n > 0, M = b^(-n)-1 otherwise.
   N is supposed to be a *primitive* cofactor of M.
   Then find a special cocktail of CM curves a` la Atkin.
   To solve the B1 problem, only consider (b, n)'s s.t. disc(b, n) = discref.

   When torsion != NULL, this means we are using some curves over
   Q(sqrt(discref)).

 */
int
process_special_blend(mpz_t tf[], int *nf, int *tried, 
		      mpz_t N, int b, int n, int c, double B1, mpz_t B2, 
		      ecm_params params, char *savefilename,
		      int discref,
		      char *torsion, int smin, int smax, int ncurves)
{
    int i;
    int ret = ECM_NO_FACTOR_FOUND;
    int tabd[][4] = {{-3, -3, 0, 0}, {-4, -4, 0, 0}, {-7, -7, 0, 0}, 
		     {-8, -8, 0, 0}, {-11, -11, 0, 0},
		     /* h = g = 2 */
		     {-15, -3, 5, 0}, 
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
		     {0, 0, 0, 0}
    };
    mpz_t sqroots[10], sqd[10];
    int tsq[10], disc;

    tsq[0] = 0;
    ret = prepare_squareroots(tf[0], tsq, sqroots, b, n, c, N);
    if(ret != ECM_NO_FACTOR_FOUND)
	return conclude_on_factor(N, tf[0], 1);
    if(torsion != NULL){
	if(tsq[0] != discref)
	    printf("#W# tsq[0]=%d != discref\n", tsq[0]);
	else{
	    printf("# Using curves with torsion %s and disc=%d",
		   torsion, discref);
	    gmp_printf(" together with B1=%1.0f B2=%Zd\n", B1, B2);
	    *tried = 1;
	    ret = process_many_curves_loop(tf, nf, N, B1, B2, params, NULL,
					   torsion, smin, smax, ncurves,
					   discref, sqroots, savefilename);
	    /* TODO: improve this? */
	}
	return ret;
    }
    mpz_init_set_ui(sqd[0], 1);
    for(i = 0; tabd[i][0] != 0; i++){
	disc = tabd[i][0];
	if(disc != discref || n % abs(disc) != 0)
	    continue;
	/* rebuild sqrt(disc) */
	if(rebuild_squareroot(sqd, tsq, sqroots, tabd[i], N)){
	    printf("# Using CM curves with disc=%d", disc);
	    gmp_printf(" together  with B1=%1.0f B2=%Zd\n", B1, B2);
	    *tried = 1;
	    ret = process_many_curves_loop(tf, nf, N, B1, B2, params,NULL, 
					   NULL, 0, 0, 1,
					   disc, sqd, savefilename);
	    if(ret != ECM_NO_FACTOR_FOUND)
		break;
	}
    }
    mpz_clear(sqd[0]);
    for(i = 0; tsq[i] != 0; i++)
	mpz_clear(sqroots[i]);
    return ret;
}

/* for N | b^n+c */
static char *
best_M_d(int *disc, int b, int n, int c)
{
    *disc = 0;
    if(c == -1)
	/* b^(2*k+1) = 1 mod N => (b^(k+1))^2 = b mod N */
	*disc = b;
    else if(c == 1 && (n % 2 == 1))
	/* b^(2*k+1) = -1 mod N => (b^(k+1))^2 = -b mod N */
	*disc = -b;
    else{ /* FIXME: squarefree part of -c or -b*c */
	if(n % 2 == 0)
	    *disc = -c;
	else
	    *disc = -b*c;
    }
#ifdef HAVE_TORSION
    /* TODO: case of b with a square prime factor */
    if(*disc != 0){
        int i, M = -1, Mi, di;
	for(i = 0; strcmp(XM_data[i][0] , "0") != 0; i++){
	    Mi = atoi(XM_data[i][0]);
	    di = atoi(XM_data[i][1]);
	    if(di == *disc)
		M = Mi;
	}
	if(M != -1){
	    char *tmp = (char *)malloc(4 * sizeof(char));
	    sprintf(tmp, "Z%d", M);
	    return tmp;
	}
    }
#endif
    return NULL;
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
    printf("  -format format where format = \"bn\" or \"plain\" (default)\n");
    printf("  -X1            select best X1(M) for b^n+/-1\n");
    printf("  -h, --help     Prints this help and exit.\n");
}

#define NFMAX 100

int
main(int argc, char *argv[])
{
    mpz_t N, tf[NFMAX], B2;
    int res = 0, smin = -1, smax = -1, ncurves = 0, method = ECM_ECM, tried;
    int nf = 0, i, bb = 0;
    double B1 = 0.0, dB2 = 0.0;
    int disc = 0, b = 0, n = 0, useX1 = 0, c;
    char *infilename = NULL, *curvesname = NULL, *torsion = NULL;
    char buf[10000], ch;
    FILE *infile = NULL;
    char *savefilename = NULL, format[20];
    ecm_params params;
    int ret;

    /* print args */
    sprintf(format, "plain");
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
	    disc = atoi(argv[2]);
	    argv += 2;
	    argc -= 2;
	}
	else if ((argc > 2) && (strcmp (argv[1], "-b") == 0)){
	    b = atoi(argv[2]);
	    argv += 2;
	    argc -= 2;
	}
	else if ((argc > 2) && (strcmp (argv[1], "-format") == 0)){
	    sprintf(format, "%s", argv[2]);
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
	else if (strcmp (argv[1], "-X1") == 0){
	    useX1 = 1;
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
    if(curvesname == NULL && torsion == NULL && disc == 0 && b == 0 
       && useX1 == 0){
	fprintf (stderr, "Not enough parameters\n");
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
    if(ncurves > NCURVE_MAX){
      fprintf(stderr, "Too many curves: %d\n", ncurves);
      exit (EXIT_FAILURE);
    }

    if(torsion != NULL || useX1){
	if(disc == 0)
	    printf("# GMP-ECM [torsion=%s:%d-%d]\n", torsion, smin, smax);
	else
	    printf("# GMP-ECM [torsion=%s:%d-%d;d=%d]\n",
		   torsion, smin, smax, disc);
    }
    else if(disc != 0){
	printf("# GMP-ECM [CM=%d]\n", disc);
	ncurves = 1; /* FIXME */
    }
    else if(useX1)
	printf("# GMP-ECM [X1:%d-%d]\n", smin, smax);
    
    mpz_init (N);
    for(i = 0; i < NFMAX; i++)
	mpz_init(tf[i]); /* for potential factors */
    ecm_init(params);
#if DEBUG_MULTI_EC >= 2
    params->verbose = 2;
#else
    /*    params->verbose = OUTPUT_DEVVERBOSE;*/
    params->verbose = OUTPUT_NORMAL;
#endif
#if MULTI_USE_ADD_SUB
    if(torsion != NULL || useX1)
	compute_s_4_add_sub(params->batch_s, (ecm_uint)B1, 0);
    else
	compute_s_4_add_sub(params->batch_s, (ecm_uint)B1, disc);
#endif
    while(fscanf(infile, "%s", buf) != EOF){
	/* read number */
	if(buf[0] == '#'){
	    /* print till end of line */
	    printf("%s", buf);
	    while((ch = getc(infile)) != '\n')
		printf("%c", ch);
	    printf("\n");
	    continue;
	}
	if(strcmp(format, "bn") == 0 || strcmp(format, "bnc") == 0){
	    /* line should be: "b n[+/-/L/M] N" 
	       or "b n[+/-]c N"
	     */
	    bb = atoi(buf);
	    /* decode */
	    ret = fscanf(infile, "%s", buf);
            ASSERT_ALWAYS(ret != EOF);
	    if(strcmp(format, "bn") == 0){
		/* buf = "n[+/-/L/M]" */
		ch = buf[strlen(buf)-1];
		buf[strlen(buf)-1] = '\0';
		n = atoi(buf);
		c = 1;
		if(ch == '-')
		    c = -1;
		else if(ch == 'L' || ch == 'M'){
		    if(bb == 5)
			c = -1;
		}
		else if(ch != '+'){
		    printf("#!# unknown suffix: %c\n", ch);
		    break;
		}
#if DEBUG_MULTI_EC >= 2
		printf("# I read: b=%d n=%d c=%d\n", bb, n, c);
#endif
	    }
	    else if(strcmp(format, "bnc") == 0){
		/* buf = "n[+/-]c" */
		sscanf(buf, "%d%c%d", &n, &ch, &c);
		if(ch == '-')
		    c = -c;
#if DEBUG_MULTI_EC >= 2
		printf("# I read: b=%d n=%d c=%d\n", bb, n, c);
#endif
	    }
	    /* read N */
	    ret = fscanf(infile, "%s", buf);
            ASSERT_ALWAYS(ret != EOF);
	    if((b > 1) && (bb != b))
		continue;
	    if(useX1){
		/* select best (M, d) */
		torsion = best_M_d(&disc, bb, n, c);
		if(torsion == NULL){
		    printf("# no level found for this number, sorry\n");
		    continue;
		}
	    }
	}
	if(mpz_set_str (N, buf, 10)){
	    fprintf (stderr, "Invalid number: %s\n", argv[1]);
	    exit (1);
	}
	res = ECM_NO_FACTOR_FOUND;
	tried = 0;
	if(method == ECM_ECM){
	    nf = 0;
	    if((strcmp(format, "bn") == 0 || strcmp(format, "bnc") == 0)
	       && disc != 0){
		res = process_special_blend(tf,&nf,&tried,N,bb,n,c,B1,B2,
					    params,savefilename,disc,
					    torsion, smin, smax, ncurves);
	    }
	    if(res == ECM_NO_FACTOR_FOUND && !tried
	       && (curvesname != NULL || torsion != NULL || disc != 0)){
		res = process_many_curves_loop(tf, &nf, N, B1, B2, params,
					       curvesname,
					       torsion, smin, smax, ncurves,
					       disc, NULL,
					       savefilename);
	    }
	}
#if DEBUG_MULTI_EC >= 2
	printf("List of factors:\n");
	for(i = 0; i < nf; i++)
	    gmp_printf("%Zd\n", tf[i]);
#endif
    }
    ecm_clear(params);
    if(infile != stdin)
	fclose(infile);
    for(i = 0; i < NFMAX; i++)
	mpz_clear(tf[i]);
    mpz_clear(N);
    mpz_clear(B2);    
    return res;
}
