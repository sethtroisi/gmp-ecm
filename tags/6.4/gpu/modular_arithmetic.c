/* Plain C stage 1 (without GMP for the critical loops).

  Copyright 2010 Julie Feltin and Paul Zimmermann.

  This file is part of the ECM Library.

  The ECM Library is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or (at your
  option) any later version.

  The ECM Library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
  License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with the ECM Library; see the file COPYING.LIB.  If not, write to
  the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
  MA 02110-1301, USA.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <sys/cdefs.h>
#include <assert.h>
#include "modular_arithmetic.h"

#define ONE ((mp_limb_t)1)
#define MASK ((ONE<<52)-ONE)
#define MASK64 (~0UL)
#define TWO52 4503599627370496.0 /* 2^52 */
#define TWO104 20282409603651670423947251286016.0 /* 2^104 */
#define TWOm52 1/4503599627370496 /* 2^(-52) */
#define C0 67108864.0                         /* 2^26 */
#define C1 302231454903657293676544.0         /* 2^78 */
#define C2 20282409603651670423947251286016.0 /* 2^104 */


/* returns number of words of N in base 2^52 */
unsigned long
dpn_size (mpz_t N)
{
  unsigned long n = mpz_sizeinbase (N, 2);

  return (n + 51) / 52;
}

/* converts N from base 2^64 to base 2^52, where n is the (maximal) number of
   words of N in base 2^52, and b[] was allocated with n words */
void
conversion64to52 (dpn_t b, mpz_t N, unsigned long n)
{
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


/* converts {b, n} from base 2^52 to mpz_t M */
void
conversion52to64 (dpn_t b, unsigned long n, mpz_t M)
{
  unsigned long i;

  mpz_set_ui (M, (unsigned long) b[n-1]);
  for (i = n - 1; i > 0;)
    {
      mpz_mul_2exp (M, M, 52);
      mpz_add_ui (M, M, (unsigned long) b[--i]);
    }
}


/* converts N from base 2^64 to base 2^52 */
void
conversion52to64bis (unsigned long *a, dpn_t b, unsigned long sizeb,
                     unsigned long n)
{
	unsigned long tabr[12] = {12,24,36,48,8,20,32,44,4,16,28,40};
	unsigned long tabl[12] = {40,28,16,56,44,32,20,60,48,36,24,12};
	unsigned long i,j,k,cmp;
	unsigned long *tmp = malloc((sizeb+1)*sizeof(unsigned long));
	/* initialization */
	for (i=0; i<sizeb+1; i++)
		tmp[i] = 0.0;
	
	/* converts b to unsigned long */
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
dpn_add (dpn_t a, dpn_t b, dpn_t c, unsigned long n)
{
  unsigned long i;
  double cy = 0.0;

  for (i = 0; i < n; i++)
    {
      cy += b[i] + c[i]; /* since 0 <= b[i], c[i] < 2^52, 0 <= cy < 2^53 */
      if (cy >= TWO52)
        {
          a[i] = cy - TWO52;
          cy = 1.0;
        }
      else
        {
          a[i] = cy;
          cy = 0.0;
        }
    }
  return cy;
}

/* {a, n} <- {b, n} + cy, where 0 <= cy <= 1, return carry out */
double
dpn_add_1 (dpn_t a, dpn_t b, unsigned long n, double cy)
{
  unsigned long i;

  for (i = 0; i < n; i++)
    {
      cy += b[i];
      if (cy == TWO52)
        {
          a[i] = 0.0;
          cy = 1.0;
        }
      else
        {
          a[i] = cy;
          cy = 0.0;
        }
    }
  return cy;
}

/* {a, n} <- {b, n} - {c, n}, returns borrow out (0 or 1) */
double
dpn_sub (dpn_t a, dpn_t b, dpn_t c, unsigned long n)
{
  unsigned long i;
  double cy = 0.0;

  for (i = 0; i < n; i++)
    {
      cy = b[i] - c[i] + cy; /* since 0 <= b[i], c[i] < 2^52 and -1 <= cy <= 0,
                                we have -2^52 <= cy < 2^52 */
                                
      if (cy >= 0)
        {
          a[i] = cy;
          cy = 0.0;
        }
      else
        {
          a[i] = cy + TWO52;
          cy = -1.0;
        }
    }

  return -cy;
}

/* ({a,n} > {b,n})?, returns 1(true), -1(false) or 0(a=b) */
int
dpn_cmp (dpn_t a, dpn_t b, unsigned long n)
{
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

/* {a, n} <- {b, n} + {c, n} mod {mod, n}
   Assumes 0 <= {b, n}, {c, n} < {mod, n} */
void
dpn_add_mod (dpn_t amod, dpn_t b, dpn_t c, dpn_t mod, unsigned long n)
{
  double cy;

  assert (dpn_check (b, mod, n));
  assert (dpn_check (c, mod, n));
  cy = dpn_add (amod, b, c, n);
	
  if (dpn_cmp (amod, mod, n) >= 0) /* (a >= mod)? */
    cy = dpn_sub (amod, amod, mod, n); /* {a, n} <- {a, n} - {mod, n} */

  assert (cy == 0.0);
  assert (dpn_check (amod, mod, n));
}
		
/* {a, n} <- {b, n} - {c, n} mod {mod, n}
   Assumes 0 <= {b, n}, {c, n} < {mod, n} */
void
dpn_sub_mod (dpn_t smod, dpn_t b, dpn_t c, dpn_t mod, unsigned long n)
{
  double cy;

  assert (dpn_check (b, mod, n));
  assert (dpn_check (c, mod, n));
  cy = dpn_sub (smod, b, c, n); /* 0 <= cy <= 1 */
	
  if (cy != 0.0) /* we should subtract 1 at smod[n] */
    {
      cy = dpn_add (smod, smod, mod, n);
      assert (cy != 0.0);
    }
  assert (dpn_check (smod, mod, n));
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
	*l += bl;

        if (*l >= TWO52)
          {
            *h += TWO52;
            *l -= TWO52;
          }
        else if (*l < 0.0)
          {
            *h -= TWO52;
            *l += TWO52;
          }
        assert (0.0 <= *l && *l < TWO52);
}


/* precomputes b[i]*c[j], return hbc, lbc such that hbc[i][j]*2^52 + lbc[i][j] = b[i]*c[j]
   with 0 < b[i],c[j] < 2^52, 0 < hbc[i][j],lbc[i][j] < 2^52 */
void
coef (dpn_t *hbc, dpn_t *lbc, dpn_t b, dpn_t c, unsigned long n)
{
	unsigned long i,j,k;
	unsigned long nn = n*n;
	dpn_t h = (dpn_t)malloc(nn*sizeof(double));
	dpn_t l = (dpn_t)malloc(nn*sizeof(double));
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



/* {a, 2n} <- {b, n} * {c, n} */
void
dpn_mul (dpn_t a, dpn_t b, dpn_t c, unsigned long n)
{
	unsigned long i,j,k;
	unsigned long nn = 2*n;
	double cy = 0.0;
	
	dpn_t *hbc = malloc(nn*sizeof(double));
	dpn_t *lbc = malloc(nn*sizeof(double));
	int cmp;

	for (i=0; i<nn; i++)
		hbc[i] = malloc(nn*sizeof(double));
	for (i=0; i<nn; i++)
		lbc[i] = malloc(nn*sizeof(double));
	/* initialization */
	for (i=0; i<nn; i++) {
		for (j=0; j<nn; j++) {
			hbc[i][j] = 0.0; lbc[i][j] = 0.0;
		}
	}

	/* if b > c precomputes b[i]*c[j], else computes c[i]*b[j] */
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

        assert (cy == 0.0);
}

/* {a, n} = {a, n} + {b, n} * c, returns carry out */
double
dpn_addmul_ui (dpn_t a, dpn_t b, unsigned long n, double c) 
{
  unsigned long i;
  double cy = 0.0, h[1], l[1], t;
  
  for (i = 0; i < n; i++)
    {
      mul (h, l, b[i], c);
      t = l[0] + cy; /* 0 <= cy < 2^53 */
      if (t >= TWO52)
        {
          t -= TWO52;
          cy = 1.0;
        }
      else
        cy = 0.0;
      t += a[i];
      if (t >= TWO52)
        {
          t -= TWO52;
          cy += 1.0;
        }
      a[i] = t;
      cy += h[0];
    }
  return cy;
}

/* dpn_mod performs {a, n} <- {a, 2n} * (2^52)^(-n) mod N,
   where N = {mod, n}, returns carry out (0 or 1) */
void
dpn_mod (dpn_t a, dpn_t mod, dpn_t mu, unsigned long n)
{
  unsigned long i,j,k;
  double cy = 0.0;
  unsigned long nn = 2*n;
  dpn_t q = (dpn_t)malloc(nn*sizeof(double));
  dpn_t r = (dpn_t)malloc(nn*sizeof(double));
  dpn_t h = (dpn_t)malloc(nn*sizeof(double));
  dpn_t l = (dpn_t)malloc(nn*sizeof(double));
  dpn_t tmp = (dpn_t)malloc(nn*sizeof(double));
  dpn_t tmp1 = (dpn_t)malloc(nn*sizeof(double));

  /* initialization */
  for (i=0; i<nn; i++) {
    q[i] = 0.0; r[i] = 0.0; h[i] = 0.0; l[i] = 0.0; tmp[i] = 0.0; tmp1[i] = 0.0;
  }

  for (i = 0; i < n; i++)
    {
      mul (h+i, q+i, mu[0], a[i]);
      /* tmp <- q[i]*mod */
      mul (h+0, l+0, q[i], mod[0]);
      h[0] = h[0] * TWOm52;
    tmp[0] = l[0];
    if (TWO52 < tmp[0]) {
      tmp[0] = tmp[0]-TWO52;
      cy++;
    }
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

  /* {a+n,n} >= mod? */
  if (dpn_cmp (a + n, mod, n) >= 0)
    {
      cy = dpn_sub (a, a + n, mod, n); /* {a,n} <- {a+n,n} - {mod,n} */
      assert (cy == 0.0);
    }
  else
    dpn_set (a, a + n, n);

  assert (dpn_check (a, mod, n));

  free (q);
  free (r);
  free (h);
  free (l);
  free (tmp);
  free (tmp1);
}

/* REDC algorithm computes {a, n} <- {b, n} * {c, n} * (2^52)^(-n),
   returns carry out (0 or 1) */
void
REDC (dpn_t a, dpn_t b, dpn_t c, dpn_t mod, dpn_t mu, unsigned long n)
{
  dpn_mul (a, b, c, n); /* a = b*c */
  dpn_mod (a, mod, mu, n);
}

/* {a, n} <- {b, n} * {c, n} mod {mod, n}, returns carry out (0 or 1) */
void
dpn_mul_mod (dpn_t a, dpn_t b, dpn_t c, mpz_t mod, dpn_t mu, unsigned long n)
{
	unsigned long i;
	unsigned long nn = 2*n;
	dpn_t redc = (dpn_t)malloc(nn*sizeof(double));
	dpn_t un = (dpn_t)malloc(n*sizeof(double));
	dpn_t d = (dpn_t)malloc(n*sizeof(double));
	mpz_t lambda;
	mpz_t exp;
	mpz_t A;

	for (i=0; i<nn; i++)
		redc[i] = 0.0;
	un[0] = 1.0;
	for (i=1; i<n; i++)
		un[i] = 0.0;
	for (i=0; i<n; i++)
		d[i] = 0.0;

	conversion64to52(d,mod,n);

        mpz_init(lambda);
        mpz_init(exp);
        mpz_init(A);
		
	dpn_mul(a,b,c,n); /* a = b*c */

	conversion52to64(a,nn,A); /* converts a from base 2^52 to mpz_t (a=A) */
	mpz_ui_pow_ui(exp,2,52);
	mpz_pow_ui(lambda,exp,n); /* lambda = (2^52)^n */
	mpz_mul(lambda,lambda,A); /* lambda*A */
	mpz_mod(lambda,lambda,mod); /* lambda*A mod[mod] */
	conversion64to52(redc,lambda,n); /* converts lambda*A mod[mod] to base 2^52 */
	
	for (i=0; i<nn; i++)
		a[i] = 0.0;
	REDC(a,redc,un,d,mu,n); /* computes (lambda*A mod[mod])*(2^52)^(-n) = b*c mod[mod] */
	
	mpz_clear(lambda); mpz_clear(exp); mpz_clear(A);
	free(redc); free(un); free(d);
}

/* {a, n} <- {b, n} */
void
dpn_set (dpn_t a, dpn_t b, unsigned long n)
{
  unsigned long i;

  for (i = 0; i < n; i++)
    a[i] = b[i];
}

/* {a, n} <- 0 */
void
dpn_zero (dpn_t a, unsigned long n)
{
  unsigned long i;

  for (i = 0; i < n; i++)
    a[i] = 0.0;
}

void
dpn_print (dpn_t a, unsigned long n)
{
  unsigned long i;

  printf ("%1.0f", a[0]);
  for (i = 1; i < n; i++)
    printf ("+%1.0f*2^%lu", a[i], 52*i);
  printf ("\n");
}

int
dpn_check (dpn_t a, dpn_t N, unsigned long n)
{
  unsigned long i;

  for (i = 0; i < n; i++)
    if (a[i] < 0.0 || a[i] >= TWO52 || a[i] != (double) (unsigned long) a[i])
      return 0;
  if (dpn_cmp (a, N, n) >= 0)
    return 0;
  return 1;
}

/* puts a into b-bit word Montgomery form, i.e., a <- a*2^(b*n) mod N,
   where n is the number of words of N in base 2^b */
void
mpz_to_montgomery (mpz_t a, mpz_t N, unsigned long b)
{
  unsigned long n = mpz_sizeinbase (N, 2);

  n = (n + b - 1) / b; /* ceil(n/b) */
  mpz_mul_2exp (a, a, n * b);
  mpz_mod (a, a, N);
}

/* puts a from b-bit word Montgomery form to classical form,
   i.e., a <- a/2^(b*n) mod N,
   where n is the number of words of N in base 2^b
   Assumes N is odd.
*/
void
mpz_from_montgomery (mpz_t a, mpz_t N, unsigned long b)
{
  unsigned long n = mpz_sizeinbase (N, 2), i;

  n = (n + b - 1) / b; /* ceil(n/b) */
  for (i = 0; i < n * b; i++)
    {
      if (mpz_divisible_ui_p (a, 2) == 0)
        mpz_add (a, a, N);
      /* now a is even */
      mpz_div_2exp (a, a, 1);
    }
}

#ifdef MAIN
int
main (int argc, char * argv[]) {
	/* assume 0 < N,N1 < N2 */
	mpz_t N;
	mpz_t N1;
	mpz_t N2;
	unsigned long nb,nc,nd;
	unsigned long n;
	unsigned long i;
        dpn_t a, b, c, d, s, ma, ms, mm, mmod, mu, un;
        unsigned long *test;
	double cy;
	mpz_t add;
	mpz_t sub;
	mpz_t moda;
	mpz_t mods;
        unsigned long maxm;
	mpz_t mult;
	mpz_t lambda;
	mpz_t tmp;
	mpz_t exp;
	mpz_t mulmod;

        if (argc != 4)
          {
            fprintf (stderr, "Usage: modular_arithmetic N N1 N2\n");
            exit (1);
          }

        mpz_init(N);
        mpz_set_str(N,argv[1],10);
        mpz_init(N1);
        mpz_set_str(N1,argv[2],10);
        mpz_init(N2);
        mpz_set_str(N2,argv[3],10); /* N2 -> modulo */

	nb = dpn_size(N); nc = dpn_size(N1); nd = dpn_size(N2);
	
	n = (nc>nb)?nc:nb;
	n = (n>nd)?n:nd;
	
	b = (dpn_t)malloc(n*sizeof(double));
	c = (dpn_t)malloc(n*sizeof(double));
	d = (dpn_t)malloc(n*sizeof(double));
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

	test = (unsigned long*)malloc(n*sizeof(unsigned long));
	for (i=0; i<n; i++)
		test[i] = 0;
	conversion52to64bis(test,b,n,mpz_size(N));
	for (i=0; i<n; i++)
		printf("test[%lu]=%lu\n",i,test[i]);
	printf("\n");
		
		
	/* add */
	a = (dpn_t)malloc(n*sizeof(double));
	for (i=0; i<n; i++)
		a[i] = 0.0;
	cy = dpn_add(a,b,c,n);
	for (i=0; i<n; i++)
		printf("a[%lu]=%f\n",i,a[i]);
	/* verif */
        mpz_init(add);
	mpz_add(add,N,N1);
	gmp_printf("add=%Zd\n\n",add);
	
	
	/* sub */
	s = (dpn_t)malloc(n*sizeof(double));
	for (i=0; i<n; i++)
		s[i] = 0.0;
	cy = dpn_sub(s,b,c,n);
	for (i=0; i<n; i++)
		printf("s[%lu]=%f\n",i,s[i]);
	/* verif */
        mpz_init(sub);
	mpz_sub(sub,N,N1);
	gmp_printf("sub=%Zd\n\n",sub);
	
	
	/* add mod */
	ma = (dpn_t)malloc(n*sizeof(double));
	for (i=0; i<n; i++)
		ma[i] = 0.0;
	dpn_add_mod (ma, b, c, d, n);
	for (i=0; i<n; i++)
		printf("ma[%lu]=%f\n",i,ma[i]);
	/* verif */
        mpz_init(moda);
	mpz_mod(moda,add,N2);
	gmp_printf("addmod=%Zd\n\n",moda);
	
	
	/* sub mod */
	ms = (dpn_t)malloc(n*sizeof(double));
	for (i=0; i<n; i++)
		ms[i] = 0.0;
	cy = dpn_sub_mod(ms,b,c,d,n);
	for (i=0; i<n; i++)
		printf("ms[%lu]=%f\n",i,ms[i]);
	/* verif */
        mpz_init(mods);
	mpz_mod(mods,sub,N2);
	gmp_printf("submod=%Zd\n\n",mods);
	
	
	/* mul */
	maxm = 2*n;
	mm = (dpn_t)malloc(maxm*sizeof(double));
	for (i=0; i<maxm; i++)
		mm[i] = 0.0;
	cy = dpn_mul(mm,b,c,n);
	for (i=0; i<maxm; i++)
		printf("mm[%lu]=%f\n",i,mm[i]);
	/* verif */
        mpz_init(mult);
	mpz_mul(mult,N,N1);
	gmp_printf("mult=%Zd\n\n",mult);
	
	
	/* mul mod */
	mmod = (dpn_t)malloc(maxm*sizeof(double));
	for (i=0; i<maxm; i++)
		mmod[i] = 0.0;
	mu = (dpn_t)malloc(sizeof(double));
	un = (dpn_t)malloc(n*sizeof(double));
	un[0] = 1.0;
	
        mpz_init(lambda);
        mpz_init(tmp);
        mpz_init(exp);
	mpz_ui_pow_ui(exp,2,52);
		/* mu = -1/N2 mod[2^52] */
	mpz_mul_ui(tmp,N2,-1);
	mpz_powm_ui(tmp,tmp,-1,exp);
	conversion64to52(mu,tmp,1);

	dpn_mul_mod(mmod,b,c,N2,mu,n);
	for (i=0; i<n; i++)
		printf("mmod[%lu]=%f\n",i,mmod[i]);
	/* verif */
        mpz_init(mulmod);
	mpz_mul(mulmod,N,N1);
	mpz_mod(mulmod,mulmod,N2);
	gmp_printf("mulmod=%Zd\n",mulmod);
	
	mpz_clear(N); mpz_clear(N1); mpz_clear(N2);
	mpz_clear(add); mpz_clear(sub); mpz_clear(moda); mpz_clear(mods);
	mpz_clear(mult); mpz_clear(tmp); mpz_clear(exp); mpz_clear(lambda); mpz_clear(mulmod);
	free(b); free(c); free(d); free(test);
	free(a); free(s); free(ma); free(ms);
	free(mm); free(mmod); free(mu); free(un);

        return 0;
}
#endif
