/* Plain C stage 1 (not using GMP for the critical loops).

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

/* Example:
   $ ./stage1 29799904256775982671863388319999573561548825027149399972531599612392671227006866151136667908641695103422986028076864929902803267437351318167549013218980573566942647077444419419003164546362008247462049 17 1000000
a=9122032453422309303975228503303119186594096764565570324571255124462257150200764885675475282720699595928510737885986837078892475054042869967855173419156369123531441860593128494129801744621300240933171
Starting point: x=15103147079049907765711344782576620836201216022308570309709685128202269764959718444195054567986529779670921068575492301131280304481783998012863228763537884832421573625653695166561750401066843178542817
After stage 1, x=17628830287311089492501926991334045672323900629597277022267290291143065118509513624278744509458318241685541031462343435507840406720158512915044990777024060216890526386080098985253934838627077872164414
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <assert.h>
#include "getprime.h"
#include "modular_arithmetic.h"

/* calculParam calculates the values of parameters of Suyama's parametrisation
   input : a random integer sigma and N
   output : parameters of Suyama's parametrisation -> a, x0, z0, u and v */
static void
calculParam (mpz_t sigma, mpz_t N, mpz_t d, mpz_t x0, mpz_t z0, mpz_t u,
             mpz_t v)
{
	mpz_t a; 
	mpz_t tmp; 
	mpz_t bezout1;
	mpz_t bezout2; 
	mpz_t gcd; 

        mpz_init(a);
        mpz_init(tmp);
        mpz_init(bezout1);
        mpz_init(bezout2);
        mpz_init(gcd);

	/* u */
	mpz_pow_ui(u,sigma,2);
	mpz_sub_ui(u,u,5);
	/* v */
	mpz_mul_ui(v,sigma,4);
	/* x0 */
	mpz_powm_ui(x0,u,3,N);
	/* z0 */
	mpz_powm_ui(z0,v,3,N);
	/* a */
	mpz_sub(a,v,u);
	mpz_pow_ui(a,a,3);
	mpz_mul_ui(tmp,u,3);
	mpz_add(tmp,tmp,v);
	mpz_mul(a,a,tmp);
	mpz_mod(a,a,N);
	mpz_pow_ui(tmp,u,3);
	mpz_mul(tmp,tmp,v);
	mpz_mul_ui(tmp,tmp,4);
	mpz_mod(tmp,tmp,N);
	mpz_gcdext(gcd,bezout1,bezout2,tmp,N); /* set gcd to the greatest
                                                  common divisor of tmp and N,
                                                  and in addition set bezout1
                                                  and bezout2 to coefficients
                                                  satisfying tmp*bezout1 +
                                                  N*bezout2 = gcd -> bezout1 =
                                                  (1/tmp)%N */
	mpz_mod(tmp,bezout1,N);
	mpz_mul(a,a,tmp);
	mpz_mod(a,a,N);
	mpz_sub_ui(a,a,2);
	mpz_mod(a,a,N);

	gmp_printf("a=%Zd\n",a);
	
	/* d = (a+2)/4 mod N */
	mpz_add_ui(d,a,2);
        while (mpz_divisible_ui_p (d, 4) == 0)
          mpz_add (d, d, N);
	mpz_divexact_ui(d,d,4);
	mpz_mod(d,d,N);
	
	/* calculation of the starting point x = (x0/z0)%N */
	mpz_gcdext(gcd,bezout1,bezout2,z0,N);
        /* set gcd to the greatest common divisor of z0 and N, and in addition
           set bezout1 and bezout2 to coefficients satisfying
           z0*bezout1 + N*bezout2 = gcd -> bezout1=(1/z0)%N */
	mpz_mod(tmp,bezout1,N);
	mpz_mul(tmp,x0,tmp); /* x0/z0 */
	mpz_mod(tmp,tmp,N);
	gmp_printf("Starting point: x=%Zd\n",tmp);
        
        /* x0 <- x0/z0, z0 <- 1 */
        mpz_set (x0, tmp);
        mpz_set_ui (z0, 1);

	mpz_clear(a);
	mpz_clear(tmp);
	mpz_clear(bezout1);
	mpz_clear(bezout2);
	mpz_clear(gcd);
}

/* duplicates the point P = (xP::zP) on the elliptic curve
   Input : the coordinates (xP::zP) of the point P to duplicate,
   d = (a+2)/4 the Suyama's parameter and N the modulus
   Output : the coordinates (x2P::z2P) of the point 2P
   Assumes the buffers x2P and z2P are allocated to 2n words.
*/
static void
duplicate (dpn_t xP, dpn_t zP, dpn_t d, dpn_t N, dpn_t mu, unsigned long n,
           dpn_t x2P, dpn_t z2P)
{
  dpn_t u, v, t;

  u = malloc (2 * n * sizeof (double));
  v = malloc (2 * n * sizeof (double));
  t = malloc (2 * n * sizeof (double));

  assert (dpn_check (xP, N, n));
  dpn_add_mod (u, xP, zP, N, n);
  dpn_mul (u, u, u, n);
  dpn_mod (u, N, mu, n);
  dpn_sub_mod (v, xP, zP, N, n);
  dpn_mul (v, v, v, n);
  dpn_mod (v, N, mu, n);
  dpn_sub_mod (t, u, v, N, n);
  dpn_mul (t, t, d, n);
  dpn_mod (t, N, mu, n); /* we have to reduce t right now, since if we reduce
                            it after adding v, we would have to shift v */
  dpn_add_mod (t, t, v, N, n);
  /* 0 <= t <= (N-1)^2+(N-1) = N*(N-1) */
  dpn_mul (x2P, u, v, n);
  dpn_mod (x2P, N, mu, n);
  dpn_sub_mod (z2P, u, v, N, n);
  dpn_mul (z2P, z2P, t, n);
  dpn_mod (z2P, N, mu, n);

  assert (dpn_check (x2P, N, n));
  assert (dpn_check (z2P, N, n));

  free (u);
  free (v);
  free (t);
}

/* differential addition P+Q given the coordinates of P, Q, and P-Q
   Input : the coordinates of points P=(xP::zP) and Q=(xQ::zQ),
   x_PminusQ ,z_PminusQ corresponding to P-Q, and N
   Output : the coordinates of P+Q=(x_PplusQ::z_PplusQ)
   Assumes the buffers x_PplusQ and z_PplusQ have 2n words.
*/
static void
addition (dpn_t xP, dpn_t zP, dpn_t xQ, dpn_t zQ, dpn_t x_PminusQ,
          dpn_t z_PminusQ, dpn_t N, dpn_t mu, unsigned long n,
          dpn_t x_PplusQ, dpn_t z_PplusQ)
{
  dpn_t u, v, uv, w, t;

  u = malloc (2 * n * sizeof (double));
  v = malloc (2 * n * sizeof (double));
  uv = malloc (n * sizeof (double));
  w = malloc (2 * n * sizeof (double));
  t = malloc (2 * n * sizeof (double));

  dpn_add_mod (u, xP, zP, N, n);
  dpn_sub_mod (uv, xQ, zQ, N, n);
  dpn_mul (u, u, uv, n);
  dpn_mod (u, N, mu, n);
  dpn_sub_mod (v, xP, zP, N, n);
  dpn_add_mod (uv, xQ, zQ, N, n);
  dpn_mul (v, v, uv, n);
  dpn_mod (v, N, mu, n);
  dpn_add_mod (w, u, v, N, n);
  dpn_mul (w, w, w, n);
  dpn_mod (w, N, mu, n);
  dpn_sub_mod (t, u, v, N, n);
  dpn_mul (t, t, t, n);
  dpn_mod (t, N, mu, n);
  assert (dpn_check (z_PminusQ, N, n));
  assert (dpn_check (w, N, n));
  if (w[0] == 1032067344955804.0)
    {
      printf ("z_PminusQ="); dpn_print (z_PminusQ, n);
      printf ("w="); dpn_print (w, n);
    }
  dpn_mul (x_PplusQ, z_PminusQ, w, n);
  dpn_mod (x_PplusQ, N, mu, n);
  if (w[0] == 1032067344955804.0)
    {
      printf ("x_PplusQ="); dpn_print (x_PplusQ, n);
    }
  assert (dpn_check (x_PplusQ, N, n));
  dpn_mul (z_PplusQ, x_PminusQ, t, n);
  dpn_mod (z_PplusQ, N, mu, n);

  assert (dpn_check (z_PplusQ, N, n));

  free (u);
  free (v);
  free (uv);
  free (w);
  free (t);
}

/* 'multiplication' allows to calculate the multiplication pi.Q on the elliptique curve
   input : integer PI which whom we want to multiplicate the point Q, the coordinates of the point Q = (xQ::zQ),
   Suyama's parameter a and N the number to be factored
   output : the coordinates of the point PI.Q = (x_PIQ::z_PIQ) */
void
multiplication (unsigned long PI, mpz_t xQ, mpz_t zQ, mpz_t d, mpz_t N,
                mpz_t muz, mpz_t x_PIQ, mpz_t z_PIQ)
{
  int i, bits=0;
  unsigned long PIm1,tmp; 
  unsigned long Kplus1;
  dpn_t Nd, x2Qd, z2Qd, x3Qd, z3Qd, x4Qd, z4Qd, x_PIQm1d, z_PIQm1d, x_2KQd,
    z_2KQd, x_2Kp1Qd, z_2Kp1Qd, x_2Kp2Qd, z_2Kp2Qd, xQd, zQd, dd, x_PIQd,
    z_PIQd, mu;
  unsigned long dn;

  dn = dpn_size (N);

  /* inputs */
  Nd = malloc (dn * sizeof(double));
  conversion64to52 (Nd, N, dn);
  xQd = malloc (dn * sizeof(double));
  zQd = malloc (dn * sizeof(double));
  dd = malloc (dn * sizeof(double));
  conversion64to52 (xQd, xQ, dn);
  conversion64to52 (zQd, zQ, dn);
  conversion64to52 (dd, d, dn);
  mu = malloc (dn * sizeof(double));
  conversion64to52 (mu, muz, dn);

  /* outputs */
  x_PIQd = malloc (dn * sizeof(double));
  z_PIQd = malloc (dn * sizeof(double));

  /* auxiliary buffers */
  x2Qd = malloc (2 * dn * sizeof(double));
  z2Qd = malloc (2 * dn * sizeof(double));
  x3Qd = malloc (2 * dn * sizeof(double));
  z3Qd = malloc (2 * dn * sizeof(double));
  x4Qd = malloc (2 * dn * sizeof(double));
  z4Qd = malloc (2 * dn * sizeof(double));
  x_PIQm1d = malloc (dn * sizeof(double));
  z_PIQm1d = malloc (dn * sizeof(double));
  x_2KQd = malloc (2 * dn * sizeof(double));
  z_2KQd = malloc (2 * dn * sizeof(double));
  x_2Kp1Qd = malloc (2 * dn * sizeof(double));
  z_2Kp1Qd = malloc (2 * dn * sizeof(double));
  x_2Kp2Qd = malloc (2 * dn * sizeof(double));
  z_2Kp2Qd = malloc (2 * dn * sizeof(double));

  PIm1 = tmp = PI-1;
	
  /* calculate size of PI-1 in base 2 */
  while (tmp)
    {
      tmp >>= 1;
      bits++;
    }

  dpn_zero (x2Qd, dn);
  dpn_zero (z2Qd, dn);

  assert (dpn_check (xQd, Nd, dn));
  duplicate (xQd, zQd, dd, Nd, mu, dn, x2Qd, z2Qd); /* x2Q and z2Q */
	
  /* if PI=2, returns 2Q=(x2Q::z2Q) */
  if (PI == 2)
    {
      dpn_set (x_PIQd, x2Qd, dn);
      dpn_set (z_PIQd, z2Qd, dn);
    }
  else
    {
      addition (xQd, zQd, x2Qd, z2Qd, xQd, zQd, Nd, mu, dn, x3Qd, z3Qd);
      /* x3Q and z3Q */
      assert (dpn_check (x2Qd, Nd, dn));
      duplicate (x2Qd, z2Qd, dd, Nd, mu, dn, x4Qd, z4Qd); /* x4Q and z4Q */
      tmp = PIm1>>(bits-2); /* read bits until bits-2 (left towards right) */
      tmp = tmp+1; /* calculation of the 2nd value of Montgomery's chain */
      /* test if the value of Montgomery's chain which comes later 2 is 3 or 4 */
      if (tmp == 3)
        {
          /* if the value of Montgomery's chain which comes later 2 is 3,
             set coordinates of 2Q to coordinates of PIQm1 and set coordinates of 3Q to coordinates of PIQ */
          dpn_set (x_PIQm1d, x2Qd, dn);
          dpn_set (z_PIQm1d, z2Qd, dn);
          dpn_set (x_PIQd, x3Qd, dn);
          dpn_set (z_PIQd, z3Qd, dn);
        }
      else
        {
          /* if the value of Montgomery's chain which comes later 2 is 4,
             set coordinates of 3Q to coordinates of PIQm1 and set coordinates of 4Q to coordinates of PIQ */
          dpn_set (x_PIQm1d, x3Qd, dn);
          dpn_set (z_PIQm1d, z3Qd, dn);
          dpn_set (x_PIQd, x4Qd, dn);
          dpn_set (z_PIQd, z4Qd, dn);
        }
      for (i=3; i<bits+1; i++)
        { /* if PI > 4 */
          Kplus1 = PIm1>>(bits-i); /* read i bits of PIm1 left towards right */
          Kplus1 = Kplus1+1;
          if (Kplus1%2 == 0)
            { /* Kplus1 is even */
              addition (x_PIQm1d, z_PIQm1d, x_PIQd, z_PIQd, xQd, zQd, Nd, mu,
                        dn, x_2Kp1Qd, z_2Kp1Qd); /* (2K+1)Q */
              assert (dpn_check (x_PIQd, Nd, dn));
              duplicate (x_PIQd, z_PIQd, dd, Nd, mu, dn, x_2Kp2Qd, z_2Kp2Qd);
              /* (2K+2)Q */
              dpn_set (x_PIQm1d, x_2Kp1Qd, dn);
              dpn_set (z_PIQm1d, z_2Kp1Qd, dn);
              dpn_set (x_PIQd, x_2Kp2Qd, dn);
              dpn_set (z_PIQd, z_2Kp2Qd, dn);
            }
          else
            {
              assert (dpn_check (x_PIQm1d, Nd, dn));
              duplicate (x_PIQm1d, z_PIQm1d, dd, Nd, mu, dn, x_2KQd, z_2KQd);
              /* (2K)Q */
              addition (x_PIQm1d, z_PIQm1d, x_PIQd, z_PIQd, xQd, zQd, Nd, mu,
                        dn, x_2Kp1Qd, z_2Kp1Qd); /* (2K+1)Q */
              dpn_set (x_PIQm1d, x_2KQd, dn);
              dpn_set (z_PIQm1d, z_2KQd, dn);
              dpn_set (x_PIQd, x_2Kp1Qd, dn);
              dpn_set (z_PIQd, z_2Kp1Qd, dn);
            }
        }
      
    }

  conversion52to64 (x_PIQd, dn, x_PIQ);
  conversion52to64 (z_PIQd, dn, z_PIQ);

  free (Nd);
  free (x2Qd);
  free (z2Qd);
  free (x3Qd);
  free (z3Qd);
  free (x4Qd);
  free (z4Qd);
  free (x_PIQm1d);
  free (z_PIQm1d);
  free (x_2KQd);
  free (z_2KQd);
  free (x_2Kp1Qd);
  free (z_2Kp1Qd);
  free (x_2Kp2Qd);
  free (z_2Kp2Qd);
  free (xQd);
  free (zQd);
  free (dd);
  free (x_PIQd);
  free (z_PIQd);
  free (mu);
}

/* 'stageOne' allows to calculate the stage one of the algorithm ECM
   input : a border B1, a random integer sigma and N
   output : the coordinates of the point Q=(xQ::zQ) */
static void
stageOne (mpz_t B1, mpz_t sigma, mpz_t N, mpz_t Q)
{
  mpz_t d;
  mpz_t x0;
  mpz_t z0;
  mpz_t u;
  mpz_t v;
  mpz_t xQ;
  mpz_t zQ;
  mpz_t xq;
  mpz_t zq;
  unsigned long PI, PIj, PIjp1;
  int j, jp1, k, i;
  mpz_t gcd;
  mpz_t bezout1, bezout2, mu;

  mpz_init (mu);
  mpz_init (d);
  mpz_init (x0);
  mpz_init (z0);
  mpz_init (u);
  mpz_init (v);
  calculParam (sigma, N, d, x0, z0, u, v); /* calculate a, x0, z0, u and v */

  /* compute mu = -1/N mod 2^(52*n) */
  mpz_set_ui (mu, 1);
  mpz_mul_2exp (mu, mu, 52 * dpn_size (N));
  mpz_set (v, mu);
  mpz_invert (mu, N, mu);
  mpz_sub (mu, v, mu);

  mpz_init(xQ);
  mpz_init(zQ);
  mpz_init(xq);
  mpz_init(zq);

  mpz_set(xQ,x0); /* initialisation : xQ = x0 */
  mpz_set(zQ,z0); /* initialisation : zQ = z0 */

  /* puts (xQ::zQ) and d into 52-bit word Montgomery form */
  mpz_to_montgomery (xQ, N, 52);
  mpz_to_montgomery (zQ, N, 52);
  mpz_to_montgomery (d, N, 52);

  for (PI = 2; mpz_cmp_ui(B1,PI)>=0; PI=getprime(PI))
    { /* returns successive odd primes, starting with 3 */
      j = 0;
      while (mpz_cmp_ui(B1,j) > 0)
        {
          jp1 = j+1;
          PIj = pow(PI,j);
          PIjp1 = pow(PI,jp1);
          if ((mpz_cmp_ui(B1,PIj)>=0)&&(mpz_cmp_ui(B1,PIjp1)<0))
            break; /* k has been found -> quit 'if' */
          j++;
        }
      k = j;
      i = 0;
      while (i < k)
        {
          /* calculate PI.Q with the function 'multiplication' */
          /* repeat the operation k times */
          multiplication (PI, xQ, zQ, d, N, mu, xq, zq);
          mpz_set(xQ,xq);
          mpz_set(zQ,zq);
          i++;
        }
    }

  /* Change (xQ::zQ) from 52-bit word Montgomery form to classical form.
     FIXME: this might not be needed, since when we divide x/z, the extra
     factor disappears */
  mpz_from_montgomery (xQ, N, 52);
  mpz_from_montgomery (zQ, N, 52);

  getprime (0);  /* { free the memory used by getprime() } */
	
  /* calculate Q = (xQ/zQ)%N, the value of the point at the end of the first step of ECM */
  mpz_init(gcd);
  mpz_init(bezout1);
  mpz_init(bezout2);
	
  mpz_gcdext(gcd,bezout1,bezout2,zQ,N);/* set gcd to the greatest common divisor of zQ and N, and in addition set bezout1 and bezout2 to coefficients satisfying zQ*bezout1 + N*bezout2 = gcd -> bezout1=(1/zQ)%N */
  mpz_mod(zQ,bezout1,N);
  mpz_mul(Q,xQ,zQ); /* xQ/zQ */
  mpz_mod(Q,Q,N);
	
  mpz_clear(d);
  mpz_clear(x0);
  mpz_clear(z0);
  mpz_clear(u);
  mpz_clear(v);
  mpz_clear(xq);
  mpz_clear(zq);
  mpz_clear (mu);
}

int
main (int argc, char*argv[])
{
  if (argc != 4)
    printf ("Error in call function\n./stage1-c N sigma B1\n");
  else
    {
      mpz_t N;
      mpz_t sigma;
      mpz_t B1;
      mpz_t x;

      mpz_init (N);
      mpz_init (sigma);
      mpz_init (B1);
      mpz_init (x);

      mpz_set_str (N, argv[1], 10); /* in base 10 */
      mpz_set_str (sigma, argv[2], 10);
      mpz_set_str (B1, argv[3], 10);

      /* check N is odd */
      if (mpz_divisible_ui_p (N, 2))
        {
          fprintf (stderr, "Error, N should be odd\n");
          exit (1);
        }
		
      stageOne (B1, sigma, N, x);
      gmp_printf ("After stage 1, x=%Zd\n",x);
		
      mpz_clear (sigma);
      mpz_clear (N);
      mpz_clear (B1);
      mpz_clear (x);
    }

  return 0;
}
