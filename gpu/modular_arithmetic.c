#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <sys/cdefs.h>


#define ONE ((mp_limb_t)1)
#define MASK ((ONE<<52)-ONE)
#define MASK64 (~0UL)
#define TWO52 4503599627370496.0 /* 2^52 */
#define TWO104 20282409603651670423947251286016.0 /* 2^104 */
#define TWOm52 1/4503599627370496 /* 2^(-52) */
#define C0 67108864.0                         /* 2^26 */
#define C1 302231454903657293676544.0         /* 2^78 */
#define C2 20282409603651670423947251286016.0 /* 2^104 */


/* computes size of N in base 2^52 */
unsigned long sizeAlloc(mpz_t N) {
	mpz_t tmp; mpz_init(tmp);
	mpz_t exp; mpz_init(exp);
	mpz_t exp_tmp; mpz_init(exp_tmp);
	unsigned long n;

	mpz_ui_pow_ui(exp,2,52);
	n = mpz_size(N); /* starts with (2^52)^sizeN */
	while (mpz_cmp_ui(tmp,0)>=0) {
		mpz_pow_ui(exp_tmp,exp,n);
		mpz_sub(tmp,N,exp_tmp); /* N - (2^52)^i */
		n++;
	}
	mpz_clear(tmp); mpz_clear(exp); mpz_clear(exp_tmp);
	return n-1;
}



/* converts N from base 2^64 to base 2^52 */
void
conversion64to52(double *b,mpz_t N,unsigned long n) {
	unsigned long tabr[14] = {52,40,28,16,4,56,44,32,20,8,60,48,36,24};
	unsigned long tabl[14] = {12,24,36,48,60,8,20,32,44,56,4,16,28,40};
	unsigned long i,j,k;

	for (i=0; i<n; i++)
		b[i] = 0.0;
	
	for (i=mpz_size(N); i<n; i++)
		N->_mp_d[i] = 0;
	
	j = 0; k = 0;
	b[0] = (N->_mp_d[0])&MASK;
	for (i=1; i<n; i++) {
		if ((j%5 == 0) && (j != 0))
			k++;
			
		if (i%15 == 0) {
			j = 0;
			b[i] = ((N->_mp_d[i-k-1])>>12)&MASK;
			b[i+1] = (N->_mp_d[i-k])&MASK;
			i++; k++;
		}
		else {
			b[i] = (((N->_mp_d[i-k-1])>>tabr[j])|((N->_mp_d[i-k])<<tabl[j]))&MASK;
			j++;
		}
	}
}


/* converts N from base 2^52 to mpz_t */
void
conversion52to64(double * b,unsigned long n,mpz_t M) {
	unsigned long i;
	mpz_t exp; mpz_init(exp);
	mpz_t exp_tmp; mpz_init(exp_tmp);
	
	mpz_ui_pow_ui(exp,2,52);
	for (i=0; i<n; i++) {
		mpz_pow_ui(exp_tmp,exp,i);
		mpz_mul_ui(exp_tmp,exp_tmp,b[i]);
		mpz_add(M,M,exp_tmp);
	}

	mpz_clear(exp);
	mpz_clear(exp_tmp);
}


/* converts N from base 2^64 to base 2^52 */
void
conversion52to64bis(unsigned long *a,double *b,unsigned long sizeb,unsigned long n) {
	unsigned long tabr[12] = {12,24,36,48,8,20,32,44,4,16,28,40};
	unsigned long tabl[12] = {40,28,16,56,44,32,20,60,48,36,24,12};
	unsigned long i,j,k,cmp;
	unsigned long *tmp = malloc(sizeb*sizeof(unsigned long));
	/* initialization */
	for (i=0; i<sizeb; i++) {
		tmp[i] = 0.0; a[i] = 0;
	}
	
	/* converts b to double and sets */
	for (i=0; i<sizeb; i++)
		tmp[i] = b[i];
	
	j = 0; k = 0; cmp = 1;
	a[0] = (tmp[0]|((tmp[1]&((1UL<<12)-1UL))<<52));
	for (i=1; i<n; i++) {
		if ((j%4 == 0) && (j != 0))
			k++;

		if (i%13 == 0) {
			j = 0;
			a[i] = (tmp[i+k]|((tmp[i+k+1]&((1UL<<12)-1UL))<<52))&MASK64;
			cmp = 0;
		}
		else {
			if ((cmp == 4) && ((i+1)%13 != 0)) {
				a[i] = ((tmp[i+k]>>tabr[j])|((tmp[i+k+1]&((1UL<<52)-1UL))<<tabr[11-j])|(tmp[i+k+2]<<tabl[j]))&MASK64;
				cmp = 0;
				j++;
			}
			else {
				a[i] = (((tmp[i+k]>>tabr[j])|(tmp[i+k+1]<<tabl[j])))&MASK64;
				j++;
			}
		}
		cmp++;
	}
	
	free(tmp);
}



/* {a, n} <- {b, n} + {c, n}, returns carry out (0 or 1) */
double
dpn_add(double *a,double *b,double *c,unsigned long n) {
	unsigned long i;
	double cy = 0.0;
	for (i = 0; i < n; i++) {
		cy += b[i] + c[i]; /* since 0 <= b[i], c[i] < 2^52, 0 <= cy < 2^53 */
		if (cy >= TWO52) {
			a[i] = cy - TWO52;
			cy = 1.0;
		}
		else {
			a[i] = cy;
			cy = 0.0;
		}
	}
	return cy;
}



/* {a, n} <- {b, n} - {c, n}, returns carry out (0 or 1) */
double
dpn_sub(double *a,double *b,double *c,unsigned long n) {
	unsigned long i;
	double cy = 0.0;
	double t;

	for (i=0; i<n; i++) {
		t = b[i]-c[i]-cy;
		if (t >= 0) {
			a[i] = t;
			cy = 0.0;
		}
		else {
			a[i] = t+TWO52;
			cy = 1.0;
		}
	}

	return cy;
}



/* ({a,n} > {b,n})?, returns 1(true), -1(false) or 0(a=b) */
int
dpn_cmp(double *a,double *b,unsigned long n) {
	unsigned long i;
	i = n;
	while (i--) {
		if (a[i] > b[i])
			return 1;
		if (a[i] < b[i])
			return -1;
		else
			i = i;
	}
	return 0;
}

/* {a, n} <- {b, n} + {c, n} mod {mod, n}, returns carry out (0 or 1) */
double
dpn_add_mod(double *amod,double *b,double *c,double *mod,unsigned long n) {
	double cy;
	
	cy = dpn_add(amod,b,c,n);
	
	while (dpn_cmp(amod,mod,n) >= 0) /* (a >= mod)? */
		cy = dpn_sub(amod,amod,mod,n); /* {a, n} <- {a, n} - {mod, n} */
		
	return cy;	
}


		
/* {a, n} <- {b, n} - {c, n} mod {mod, n}, returns carry out (0 or 1) */
double
dpn_sub_mod(double *smod,double *b,double *c,double *mod,unsigned long n) {
	unsigned long i;
	double cy;
	
	cy = dpn_sub(smod,b,c,n);
	
	if (cy != 0) {
		cy = dpn_add(smod,smod,mod,n);
		if (cy == 0)
			abort();
	}
	return cy;	
}



/* Assume b, c, h are integers with 0 <= b, c < 2^52.
  Return h, l such that h*2^52 + l = b*c */
void
mul (double *h, double *l, double b, double c)
{
	double bh, bl, ch, cl, m;

	bh = C1 + b; /* ulp(b1) = 2^26 */
	bh = bh - C1;
	bl = b - bh;
	if (bl < 0.0)
	{
		bh -= C0;
		bl += C0;
	}

	ch = C1 + c; /* ulp(b1) = 2^26 */
	ch = ch - C1;
	cl = c - ch;
	if (cl < 0.0)
	{
		ch -= C0;
		cl += C0;
	}

	*h = bh * ch;
	m = bh * cl + bl * ch;
	*l = bl * cl;

	/* b*c = h + m + l where h is multiple of 2^52, m is multiple of 2^26,
	0 <= h/2^52 <= (2^26-1)^2, 0 <= m/2^26 <= 2*(2^26-1)^2,
	0 <= l <= (2^26-1)^2 */

	bh = C2 + m;
	bh = bh - C2;
	bl = m - bh;

	/* now ulp(bh) = 2^52, ulp(bl) = 2^26, 0 <= bh/2^52 <= 2^27-4
	0 <= bl/2^26 < 2^26 */

	*h += bh; /* 0 <= h/2^52 < 2^52 */
	*l += bl; /* 0 <= l < 2^52 */
}


/* precomputes b[i]*c[j], return hbc, lbc such that hbc[i][j]*2^52 + lbc[i][j] = b[i]*c[j]
   with 0 < b[i],c[j] < 2^52, 0 < hbc[i][j],lbc[i][j] < 2^52 */
void
coef(double **hbc,double **lbc,double *b,double *c,unsigned long n) {
	unsigned long i,j,k;
	unsigned long nn = n*n;
	double *h = (double*)malloc(nn*sizeof(double));
	double *l = (double*)malloc(nn*sizeof(double));
	/* initialization */
	for (i=0; i<nn; i++) {
		h[i] = 0.0; l[i] = 0.0;
	}

	k = 0;
	for (i=0; i<n; i++) {
		for (j=0; j<n; j++) {
			mul(h+k,l+k,b[j],c[i]);
			if (l[k] < 0) {
				l[k] = l[k]+TWO52;
				h[k] = h[k]-TWO52;
			}
			h[k] = h[k]*TWOm52;
			hbc[j][i] = h[k];
			lbc[j][i] = l[k];
			k++;
		}
	}
	free(h); free(l);
}



/* {a, n} <- {b, n} * {c, n}, returns carry out (0 or 1) */
double
dpn_mul(double *a,double *b,double *c,unsigned long n) {
	unsigned long i,j,k;
	unsigned long nn = 2*n;
	double cy = 0.0;
	
	double ** hbc = malloc(nn*sizeof(double));
	for (i=0; i<nn; i++)
		hbc[i] = malloc(nn*sizeof(double));
	double ** lbc = malloc(nn*sizeof(double));
	for (i=0; i<nn; i++)
		lbc[i] = malloc(nn*sizeof(double));
	/* initialization */
	for (i=0; i<nn; i++) {
		for (j=0; j<nn; j++) {
			hbc[i][j] = 0.0; lbc[i][j] = 0.0;
		}
	}

	/* if b > c precomputes b[i]*c[j], else computes c[i]*b[j] */
	int cmp;
	cmp = dpn_cmp(b,c,n);
	if (cmp == 1)
		coef(hbc,lbc,b,c,n);
	else
		coef(hbc,lbc,c,b,n);

	a[0] = lbc[0][0];
	k = 0;
	/* the coefficients b[i]*c[j] are already calculated,
	   it remains to add the coefficients as in schoolbook multiplication */
	for(i=1; i<nn; i++) {
		a[i] = cy;
		cy = 0.0;
		if (i < n) {
			for (j=0; j<=i; j++) {
				a[i] += lbc[i-j][j]; /* since 0 <= a[i], lbc[i-j][j] < 2^52, 0 <= a+lbc < 2^53 */
				if (TWO52 < a[i]) {
					a[i] = a[i]-TWO52;
					cy++;
				}
			}
			for (j=0; j<=i-1; j++) {
				a[i] += hbc[i-1-j][j]; /* since 0 <= a[i], hbc[i-1-j][j] < 2^52, 0 <= a+hbc < 2^53 */
				if (TWO52 < a[i]) {
					a[i] = a[i]-TWO52;
					cy++;
				}
			}
		}
		else {
			k++;
			for (j=k; j<=n; j++) {
				a[i] += hbc[i-j][j-1];
				if (TWO52 < a[i]) {
					a[i] = a[i]-TWO52;
					cy++;
				}
			}
			for (j=k; j<=n-1; j++) {
				a[i] += lbc[i-j][j];
				if (TWO52 < a[i]) {
					a[i] = a[i]-TWO52;
					cy++;
				}
			}
		}
	}
	for(i=0; i<nn; i++) {
		free(hbc[i]); free(lbc[i]);
	}
	free(hbc); free(lbc);
	
	return cy;
}



/* REDC algorithm computes {a, n} <- {b, n} * {c, n} * (2^52)^(-n),
   returns carry out (0 or 1) */
double
REDC(double *a,double *b,double *c,double *mod,double *mu,unsigned long n) {
	unsigned long i,j,k;
	double cy;
	unsigned long nn = 2*n;
	double *q = (double*)malloc(nn*sizeof(double));
	double *r = (double*)malloc(nn*sizeof(double));
	double *h = (double*)malloc(nn*sizeof(double));
	double *l = (double*)malloc(nn*sizeof(double));
	double *tmp = (double*)malloc(nn*sizeof(double));
	double *tmp1 = (double*)malloc(nn*sizeof(double));
	/* initialization */
	for (i=0; i<nn; i++) {
		q[i] = 0.0; r[i] = 0.0; h[i] = 0.0; l[i] = 0.0; tmp[i] = 0.0; tmp1[i] = 0.0;
	}
	
	cy = dpn_mul(a,b,c,n); /* a = b*c */
		
	for (i=0; i<n; i++) {
		mul(h+i,q+i,mu[0],a[i]);
		if (q[i] < 0)
			q[i] = q[i]+TWO52;
		/* tmp <- q[i]*mod */
		mul(h+0,l+0,q[i],mod[0]);
		if (l[0] < 0) {
			l[0] = l[0]+TWO52;
			h[0] = h[0]-TWO52;
		}
		h[0] = h[0]*TWOm52;
		tmp[0] = l[0];

		for (j=1; j<=n; j++) {
			mul(h+j,l+j,q[i],mod[j]);
			if (l[j] < 0) {
				l[j] = l[j]+TWO52;
				h[j] = h[j]-TWO52;
			}
			h[j] = h[j]*TWOm52;
			tmp[j] = l[j]+h[j-1]+cy;
			cy = 0.0;
			if (TWO52 < tmp[j]) {
				tmp[j] = tmp[j]-TWO52;
				cy++;
			}
		}
		/* q[i]*mod*(2^52)^i */
		k = i;
		for (j=0; j<nn-i; j++) {
			tmp1[k] = tmp[j];
			k++;
		}
		/* a <- a+q[i]*mod*(2^52)^i */
		cy = dpn_add(a,a,tmp1,nn);
		for (j=0; j<nn; j++) {
			tmp[j] = 0.0; tmp1[j] = 0.0;
		}
	}
	/* r <- a*(2^52)^(-n) */
	k = n;
	for (i=0; i<n; i++) {
		r[i] = a[k];
		k++;
	}
	/* r >= (2^52)^n? */
	tmp[n+1]=1; /* tmp=(2^52)^n */
	if (dpn_cmp(r,tmp,nn) == 1)
		cy = dpn_sub(a,r,mod,nn); /* a <- r-mod */
	else {
		for (i=0; i<nn; i++)
			a[i] = r[i]; /* a <- r */
	}
	
	free(q); free(r); free(h); free(l); free(tmp); free(tmp1);
	return cy;
}



/* {a, n} <- {b, n} * {c, n} mod {mod, n}, returns carry out (0 or 1) */
void
dpn_mul_mod(double *a,double *b,double *c,mpz_t mod,double *mu,unsigned long n) {
	double cy;
	unsigned long i;
	unsigned long nn = 2*n;
	double *redc = (double*)malloc(nn*sizeof(double));
	for (i=0; i<nn; i++)
		redc[i] = 0.0;
	double *un = (double*)malloc(n*sizeof(double));
	un[0] = 1.0;
	for (i=1; i<n; i++)
		un[i] = 0.0;
	double *d = (double*)malloc(n*sizeof(double));
	for (i=0; i<n; i++)
		d[i] = 0.0;

	conversion64to52(d,mod,n);
	
	mpz_t lambda; mpz_init(lambda);
	mpz_t exp; mpz_init(exp);
	mpz_t A; mpz_init(A);
		
	cy = dpn_mul(a,b,c,n); /* a = b*c */
	conversion52to64(a,nn,A); /* converts a from base 2^52 to mpz_t (a=A) */
	mpz_ui_pow_ui(exp,2,52);
	mpz_pow_ui(lambda,exp,n);
	mpz_mul(lambda,lambda,A);
	mpz_mod(lambda,lambda,mod); /* lambda*A mod[mod] */
	conversion64to52(redc,lambda,n); /* converts lambda*A mod[mod] to base 2^52 */
	
	for (i=0; i<nn; i++)
		a[i] = 0.0;
	cy = REDC(a,redc,un,d,mu,n); /* computes (lambda*A mod[mod])*(2^52)^(-n) = b*c mod[mod] */
		
	mpz_clear(lambda); mpz_clear(exp); mpz_clear(A);
	free(redc); free(un); free(d);
}




main (int argc, char * argv[]) {
	/* assume 0 < N,N1 < N2 */
	mpz_t N; mpz_init(N); mpz_set_str(N,argv[1],10);
	mpz_t N1; mpz_init(N1); mpz_set_str(N1,argv[2],10);
	mpz_t N2; mpz_init(N2); mpz_set_str(N2,argv[3],10); /* N2 -> modulo */

	unsigned long nb,nc,nd;
	nb = sizeAlloc(N); nc = sizeAlloc(N1); nd = sizeAlloc(N2);
	
	unsigned long n;
	n = (nc>nb)?nc:nb;
	n = (n>nd)?n:nd;
	
	unsigned long i;
	
	double *b = (double*)malloc(n*sizeof(double));
	double *c = (double*)malloc(n*sizeof(double));
	double *d = (double*)malloc(n*sizeof(double));
	for (i=0; i<n; i++) {
		b[i] = 0.0; c[i] = 0.0; d[i] = 0.0;
	}
	conversion64to52(b,N,nb);
	conversion64to52(c,N1,nc);
	conversion64to52(d,N2,nd);
	for (i=0; i<n; i++)
		printf("b[%lu]=%f\n",i,b[i]);
	printf("\n");
	for (i=0; i<n; i++)
		printf("c[%lu]=%f\n",i,c[i]);
	printf("\n");
	for (i=0; i<n; i++)
		printf("d[%lu]=%f\n",i,d[i]);
	printf("\n");

	unsigned long *test = (unsigned long*)malloc(n*sizeof(unsigned long));
	for (i=0; i<n; i++)
		test[i] = 0;
	conversion52to64bis(test,b,n,mpz_size(N));
	for (i=0; i<n; i++)
		printf("test[%lu]=%lu\n",i,test[i]);
	printf("\n");
		
		
	/* add */
	double cy;
	double *a = (double*)malloc(n*sizeof(double));
	for (i=0; i<n; i++)
		a[i] = 0.0;
	cy = dpn_add(a,b,c,n);
	for (i=0; i<n; i++)
		printf("a[%lu]=%f\n",i,a[i]);
	/* verif */
	mpz_t add; mpz_init(add);
	mpz_add(add,N,N1);
	gmp_printf("add=%Zd\n\n",add);
	
	
	/* sub */
	double *s = (double*)malloc(n*sizeof(double));
	for (i=0; i<n; i++)
		s[i] = 0.0;
	cy = dpn_sub(s,b,c,n);
	for (i=0; i<n; i++)
		printf("s[%lu]=%f\n",i,s[i]);
	/* verif */
	mpz_t sub; mpz_init(sub);
	mpz_sub(sub,N,N1);
	gmp_printf("sub=%Zd\n\n",sub);
	
	
	/* add mod */
	double *ma = (double*)malloc(n*sizeof(double));
	for (i=0; i<n; i++)
		ma[i] = 0.0;
	cy = dpn_add_mod(ma,b,c,d,n);
	for (i=0; i<n; i++)
		printf("ma[%lu]=%f\n",i,ma[i]);
	/* verif */
	mpz_t moda; mpz_init(moda);
	mpz_mod(moda,add,N2);
	gmp_printf("addmod=%Zd\n\n",moda);
	
	
	/* sub mod */
	double *ms = (double*)malloc(n*sizeof(double));
	for (i=0; i<n; i++)
		ms[i] = 0.0;
	cy = dpn_sub_mod(ms,b,c,d,n);
	for (i=0; i<n; i++)
		printf("ms[%lu]=%f\n",i,ms[i]);
	/* verif */
	mpz_t mods; mpz_init(mods);
	mpz_mod(mods,sub,N2);
	gmp_printf("submod=%Zd\n\n",mods);
	
	
	/* mul */
	unsigned long maxm = 2*n;
	double *mm = (double*)malloc(maxm*sizeof(double));
	for (i=0; i<maxm; i++)
		mm[i] = 0.0;
	cy = dpn_mul(mm,b,c,n);
	for (i=0; i<maxm; i++)
		printf("mm[%lu]=%f\n",i,mm[i]);
	/* verif */
	mpz_t mult; mpz_init(mult);
	mpz_mul(mult,N,N1);
	gmp_printf("mult=%Zd\n\n",mult);
	
	
	/* mul mod */
	double *mmod = (double*)malloc(maxm*sizeof(double));
	for (i=0; i<maxm; i++)
		mmod[i] = 0.0;
	double *mu = (double*)malloc(sizeof(double));
	double *un = (double*)malloc(n*sizeof(double));
	un[0] = 1.0;
	
	mpz_t lambda; mpz_init(lambda);
	mpz_t tmp; mpz_init(tmp);
	mpz_t exp; mpz_init(exp);
	mpz_ui_pow_ui(exp,2,52);
		/* mu = -N2^(-1) mod[2^52] */
	mpz_mul_ui(tmp,N2,-1);
	mpz_powm_ui(tmp,tmp,-1,exp);
	conversion64to52(mu,tmp,1);
	
	dpn_mul_mod(mmod,b,c,N2,mu,n);
	for (i=0; i<maxm; i++)
		printf("mmod[%lu]=%f\n",i,mmod[i]);
	/* verif */
	mpz_t mulmod; mpz_init(mulmod);
	mpz_mul(mulmod,N,N1);
	mpz_mod(mulmod,mulmod,N2);
	gmp_printf("mulmod=%Zd\n",mulmod);
	
	mpz_clear(N); mpz_clear(N1); mpz_clear(N2);
	mpz_clear(add); mpz_clear(sub); mpz_clear(moda); mpz_clear(mods);
	mpz_clear(mult); mpz_clear(tmp); mpz_clear(exp); mpz_clear(lambda); mpz_clear(mulmod);
	free(b); free(c); free(d); free(test);
	free(a); free(s); free(ma); free(ms);
	free(mm); free(mmod); free(mu); free(un);
	
}

