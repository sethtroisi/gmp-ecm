/* 

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
#include <math.h>
#include <gmp.h>


/* calculParam calculates the values of parameters of Suyama's parametrisation
   input : a random integer sigma and N
   output : parameters of Suyama's parametrisation -> a, x0, z0, u and v */
void
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



/* 'duplicate' allows to duplicate a point P = (xP::zP) in the multiplication pi.Q on the elliptique curve
   input : the coordonates (xP::zP) of the point P to duplicate, a the Suyama's parameter and N
   output : the coordonates (x2P::z2P) of the point 2P */
void duplicate(mpz_t xP,mpz_t zP,mpz_t d,mpz_t N,mpz_t x2P,mpz_t z2P) {
	mpz_t u; 
	mpz_t v; 
	mpz_t t; 

        mpz_init(u);
        mpz_init(v);
        mpz_init(t);

	/* u */
	mpz_add(u,xP,zP);
	mpz_mul(u,u,u);
	mpz_mod(u,u,N);
	/* v */
	mpz_sub(v,xP,zP);
	mpz_mul(v,v,v);
	mpz_mod(v,v,N);
	/* t */
	mpz_sub(t,u,v);
	mpz_mul(t,t,d);
	mpz_add(t,t,v);
	mpz_mod(t,t,N);
	/* x2p */
	mpz_mul(x2P,u,v);
	mpz_mod(x2P,x2P,N);
	/* z2p */
	mpz_sub(z2P,u,v);
	mpz_mul(z2P,z2P,t);
	mpz_mod(z2P,z2P,N);

	mpz_clear(u);
	mpz_clear(v);
	mpz_clear(t);
}



/* 'addition' allows to add (P,Q->P+Q) in the multiplication pi.Q on the elliptique curve
   input : the coordinates of points P=(xP::zP) and Q=(xQ::zQ),
   x_PminusQ ,z_PminusQ corresponding in the coordonates of the difference P-Q=(x_P-Q::z_P-Q) and N
   output : the coordinates of P+Q=(x_P+Q::z_P+Q) */
void addition(mpz_t xP,mpz_t zP,mpz_t xQ,mpz_t zQ,mpz_t x_PminusQ,mpz_t z_PminusQ,mpz_t N,mpz_t x_PplusQ,mpz_t z_PplusQ) {
	mpz_t u; 
	mpz_t v; 
	mpz_t uv; 
	mpz_t w; 
	mpz_t t; 

        mpz_init(u);
        mpz_init(v);
        mpz_init(uv);
        mpz_init(w);
        mpz_init(t);

	/* u */
	mpz_add(u,xP,zP);
	mpz_sub(uv,xQ,zQ);
	mpz_mul(u,u,uv);
	mpz_mod(u,u,N);
	/* v */
	mpz_sub(v,xP,zP);
	mpz_add(uv,xQ,zQ);
	mpz_mul(v,v,uv);
	mpz_mod(v,v,N);
	/* w */
	mpz_add(w,u,v);
	mpz_mul(w,w,w);
	mpz_mod(w,w,N);
	/* t */
	mpz_sub(t,u,v);
	mpz_mul(t,t,t);
	mpz_mod(t,t,N);
	/* x_PplusQ */
	mpz_mul(x_PplusQ,z_PminusQ,w);
	mpz_mod(x_PplusQ,x_PplusQ,N);
	/* z_PplusQ */
	mpz_mul(z_PplusQ,x_PminusQ,t);
	mpz_mod(z_PplusQ,z_PplusQ,N);

	mpz_clear(u);
	mpz_clear(v);
	mpz_clear(uv);
	mpz_clear(w);
	mpz_clear(t);
}



/* 'multiplication' allows to calculate the multiplication pi.Q on the elliptique curve
   input : integer PI which whom we want to multiplicate the point Q, the coordonates of the point Q = (xQ::zQ),
   Suyama's parameter a and N the number to be factored
   output : the coordonates of the point PI.Q = (x_PIQ::z_PIQ) */
void
multiplication (unsigned long PI, mpz_t xQ, mpz_t zQ, mpz_t d, mpz_t N,
                mpz_t sigma, mpz_t x_PIQ, mpz_t z_PIQ)
{
  int i, bits=0;
  unsigned long PIm1,tmp; 
  unsigned long Kplus1;
  mpz_t x2Q; 
  mpz_t z2Q; 
  mpz_t x3Q; 
  mpz_t z3Q; 
  mpz_t x4Q; 
  mpz_t z4Q; 
  mpz_t x_PIQm1;
  mpz_t z_PIQm1; 
  mpz_t x_2KQ; 
  mpz_t z_2KQ; 
  mpz_t x_2Kp1Q; 
  mpz_t z_2Kp1Q; 
  mpz_t x_2Kp2Q; 
  mpz_t z_2Kp2Q;  

  PIm1 = tmp = PI-1;
	
  mpz_init(x2Q);
  mpz_init(z2Q);
  mpz_init(x3Q);
  mpz_init(z3Q);
  mpz_init(x4Q);
  mpz_init(z4Q);

  mpz_init(x_PIQm1);
  mpz_init(z_PIQm1);

  mpz_init(x_2KQ);
  mpz_init(z_2KQ);
  mpz_init(x_2Kp1Q);
  mpz_init(z_2Kp1Q);
  mpz_init(x_2Kp2Q);
  mpz_init(z_2Kp2Q);
	
  /* calculate size of PI-1 in base 2 */
  while (tmp)
    {
      tmp >>= 1;
      bits++;
    }

  duplicate(xQ,zQ,d,N,x2Q,z2Q); /* x2Q and z2Q */
	
  /* if PI=2, returns 2Q=(x2Q::z2Q) */
  if (PI == 2)
    {
      mpz_set(x_PIQ,x2Q);
      mpz_set(z_PIQ,z2Q);
    }
  else
    {
      addition(xQ,zQ,x2Q,z2Q,xQ,zQ,N,x3Q,z3Q); /* x3Q and z3Q */
      duplicate(x2Q,z2Q,d,N,x4Q,z4Q); /* x4Q and z4Q */
      tmp = PIm1>>(bits-2); /* read bits until bits-2 (left towards right) */
      tmp = tmp+1; /* calculation of the 2nd value of Montgomery's chain */
      /* test if the value of Montgomery's chain which comes later 2 is 3 or 4 */
      if (tmp == 3)
        {
          /* if the value of Montgomery's chain which comes later 2 is 3,
             set coordonates of 2Q to coordonates of PIQm1 and set coordonates of 3Q to coordonates of PIQ */
          mpz_set(x_PIQm1,x2Q);
          mpz_set(z_PIQm1,z2Q);
          mpz_set(x_PIQ,x3Q);
          mpz_set(z_PIQ,z3Q);
        }
      else
        {
          /* if the value of Montgomery's chain which comes later 2 is 4,
             set coordonates of 3Q to coordonates of PIQm1 and set coordonates of 4Q to coordonates of PIQ */
          mpz_set(x_PIQm1,x3Q);
          mpz_set(z_PIQm1,z3Q);
          mpz_set(x_PIQ,x4Q);
          mpz_set(z_PIQ,z4Q);
        }
      for (i=3; i<bits+1; i++)
        { /* if PI > 4 */
          Kplus1 = PIm1>>(bits-i); /* read i bits of PIm1 left towards right */
          Kplus1 = Kplus1+1;
          if (Kplus1%2 == 0) { /* Kplus1 is an even */
            addition(x_PIQm1,z_PIQm1,x_PIQ,z_PIQ,xQ,zQ,N,x_2Kp1Q,z_2Kp1Q); /* (2K+1)Q */
            duplicate(x_PIQ,z_PIQ,d,N,x_2Kp2Q,z_2Kp2Q); /* (2K+2)Q */
            mpz_set(x_PIQm1,x_2Kp1Q);
            mpz_set(z_PIQm1,z_2Kp1Q);
            mpz_set(x_PIQ,x_2Kp2Q);
            mpz_set(z_PIQ,z_2Kp2Q);
          }
          else
            {
              duplicate(x_PIQm1,z_PIQm1,d,N,x_2KQ,z_2KQ); /* (2K)Q */
              addition(x_PIQm1,z_PIQm1,x_PIQ,z_PIQ,xQ,zQ,N,x_2Kp1Q,z_2Kp1Q); /* (2K+1)Q */
              mpz_set(x_PIQm1,x_2KQ);
              mpz_set(z_PIQm1,z_2KQ);
              mpz_set(x_PIQ,x_2Kp1Q);
              mpz_set(z_PIQ,z_2Kp1Q);
            }
        }
      
    }

  mpz_clear(x2Q);
  mpz_clear(z2Q);
  mpz_clear(x3Q);
  mpz_clear(z3Q);
  mpz_clear(x4Q);
  mpz_clear(z4Q);
  mpz_clear(x_PIQm1);
  mpz_clear(z_PIQm1);
  mpz_clear(x_2KQ);
  mpz_clear(z_2KQ);
  mpz_clear(x_2Kp1Q);
  mpz_clear(z_2Kp1Q);
  mpz_clear(x_2Kp2Q);
  mpz_clear(z_2Kp2Q);
}



/* 'stageOne' allows to calculate the stage one of the algorithm ECM
   input : a border B1, a random integer sigma and N
   output : the coordonates of the point Q=(xQ::zQ) */
void
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
  unsigned long PI,PIj,PIjp1;
  int j,jp1,k,i;
  mpz_t gcd;
  mpz_t bezout1;
  mpz_t bezout2;

  mpz_init (d);
  mpz_init (x0);
  mpz_init (z0);
  mpz_init (u);
  mpz_init (v);
  calculParam (sigma, N, d, x0, z0, u, v); /* calculate a, x0, z0, u and v */

  mpz_init(xQ); mpz_set(xQ,x0); /* initialisation : xQ = x0 */
  mpz_init(zQ); mpz_set(zQ,z0); /* initialisation : zQ = z0 */
	
  mpz_init(xq);
  mpz_init(zq);

  for (PI = 2; mpz_cmp_ui(B1,PI)>=0; PI=getprime(PI))
    { /* returns successive odd primes, starting with 3 */
      j = 0;
      while (mpz_cmp_ui(B1,j) > 0)
        {
          jp1 = j+1;
          PIj = pow(PI,j);
          PIjp1 = pow(PI,jp1);
          if ((mpz_cmp_ui(B1,PIj)>=0)&&(mpz_cmp_ui(B1,PIjp1)<0))
            {
              k = j;
              break; /* k has been found -> quit 'if' */
            }
          j++;
        }
      i = 0;
      while (i < k)
        {
          /* calculate PI.Q with the function 'multiplication' */
          /* repeat the operation k times */
          multiplication (PI,xQ,zQ,d,N,sigma,xq,zq);
          mpz_set(xQ,xq);
          mpz_set(zQ,zq);
          i++;
        }
    }
	
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
}

