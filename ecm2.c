/* Elliptic Curve Method implementation: stage 2 routines.

  Copyright (C) 2002 Alexander Kruppa and Paul Zimmermann.

  This program is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 2 of the License, or (at your
  option) any later version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
  more details.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.
*/

#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "ecm.h"

int duplicateW (mpz_t, mpz_t, mpz_t, mpz_t, mpz_t, mpz_t, mpz_t, mpz_t, mpz_t);
int addW (mpz_t, mpz_t, mpz_t, mpz_t, mpz_t, mpz_t, mpz_t, mpz_t, mpz_t,
          mpz_t);
int multiplyW (mpz_t, mpz_t, mpz_t, mpz_t, mpz_t, unsigned long, mpz_t, mpz_t,
               mpz_t, mpz_t);

#define PTR(x) ((x)->_mp_d)
#define getbit(x,i) (PTR(x)[i/mp_bits_per_limb] & ((mp_limb_t)1<<(i%mp_bits_per_limb)))

/*
  (x1:y1) <- 2*(x:y) where (x:y) can be identical to (x1:y1).
  a is the Weierstrass parameter, u and v are auxiliary variables.
*/
int
duplicateW (mpz_t p, mpz_t x1, mpz_t y1, mpz_t x, mpz_t y, mpz_t n, mpz_t a,
            mpz_t u, mpz_t v)
{
  mpz_mul_2exp (u, y, 1);
  mpz_gcdext (p, v, NULL, u, n);
  if (mpz_cmp_ui(p,1) != 0)
    return 1;
  mpz_mul (u, x, x);
  mpz_mul_ui (u, u, 3);
  mpz_add (u, u, a);
  mpz_mod (u, u, n);
  mpz_mulmod (p, u, v, n);
  mpz_mul (u, p, p);
  mpz_mul_2exp (v, x, 1);
  mpz_sub (u, u, v);
  mpz_mod (u, u, n);
  mpz_sub (v, x, u);
  mpz_mul (v, v, p);
  mpz_sub (y1, v, y);
  mpz_mod (y1, y1, n);
  mpz_set (x1, u); /* can we avoid this mpz_set (using mpz_swap perhaps)? */

  return 0;
}

/* Computes (x:y) <- (x1:y1) + (x2:y2).
   Returns non-zero iff a factor is found (then it is stored in p).
   n is the number to factor.
   u, v are auxiliary variables.
*/
int
addW (mpz_t p, mpz_t x, mpz_t y, mpz_t x1, mpz_t y1, mpz_t x2, mpz_t y2,
      mpz_t n, mpz_t u, mpz_t v)
{
  mpz_sub (u, x2, x1);
  mpz_gcdext (p, v, NULL, u, n);
  if (mpz_cmp_ui (p, 1) != 0)
    return 1;
  mpz_sub (p, y2, y1);
  mpz_mulmod (p, v, p, n);
  mpz_mul (u, p, p);
  mpz_sub (u, u, x1);
  mpz_sub (v, u, x2);
  mpz_mod (v, v, n);
  mpz_sub (u, x1, v);
  mpz_mul (u, u, p);
  mpz_sub (y, u, y1);
  mpz_mod (y, y, n);
  mpz_set (x, v); /* can we avoid this copy? */

  return 0;
}

#if 0
/* (x1:y1) <- q*(x:y) where q is a small integer.
   a is the Weierstrass parameter.
   u and v are auxiliary variables.
   Should work when (x1:y1) and (x:y) are the same variables.
   Return non-zero iff a factor was found.
*/
int
multiplyW (mpz_t p, mpz_t x1, mpz_t y1, mpz_t x, mpz_t y, unsigned long q,
           mpz_t n, mpz_t a, mpz_t u, mpz_t v)
{
  unsigned long j, r, restore;
  mpz_t x2, y2;

  restore = (x1 == x);
  if (restore)
    {
      mpz_init (x2);
      mpz_init (y2);
      x1 = x2;
      y1 = y2;
    }
  for (r = q, j = 1; r != 1; r /= 2, j <<= 1);
  j >>= 1; 
  r = duplicateW (p, x1, y1, x, y, n, a, u, v);
  if (r)
    return r;
  if (q & j)
    r = addW (p, x1, y1, x1, y1, x, y, n, u, v);
  if (r)
    return r;
  j >>= 1;
  while (j != 0)
    {
      if (duplicateW (p, x1, y1, x1, y1, n, a, u, v))
        return 1;
      if (q & j)
        if (addW (p, x1, y1, x1, y1, x, y, n, u, v))
          return 1;
    j >>= 1;
  }
  if (restore)
    {
      mpz_set (x, x1);
      mpz_set (y, y1);
      mpz_clear (x2);
      mpz_clear (y2);
    }

  return 0;
}
#endif

/* (x1:y1) <- q*(x:y) where q is a large integer */
int
multiplyW2 (mpz_t p, mpz_t x1, mpz_t y1, mpz_t x, mpz_t y, mpz_t q, mpz_t n,
            mpz_t a, mpz_t u, mpz_t v)
{
  long j;
  int sign_q;

  sign_q = mpz_cmp_ui (q, 0);

  if (sign_q == 0)
    {
      fprintf (stderr, "Error: q=0 in multiplyW2\n");
      exit (1);
    }

  if (sign_q < 0)
    mpz_neg (q, q);

  if (mpz_cmp_ui (q, 1) == 0)
    {
      mpz_set (x1, x);
      mpz_set (y1, y);
      goto exit_multiplyW2;
    }
  j = mpz_sizeinbase (q, 2) - 2;
  if (duplicateW (p, x1, y1, x, y, n, a, u, v))
    return 1;
  if (getbit(q,j) && addW (p, x1, y1, x1, y1, x, y, n, u, v))
    return 1;
  while (--j >= 0)
    {
      if (duplicateW (p, x1, y1, x1, y1, n, a, u, v))
        return 1;
      if (getbit(q,j) && addW (p, x1, y1, x1, y1, x, y, n, u, v))
        return 1;
  }

 exit_multiplyW2:
  if (sign_q < 0)
    {
      mpz_neg (y1, y1);
      mpz_neg (q, q);
    }

  return 0;
}

/* Input: x[0]..x[e], y[0]..y[e].
   Assumes x[e+1]..x[2e] and y[e+1]..y[2e+1] contains e and e+1 more cells.

   Performs the following loop with only one gcdext, using Montgomery's trick:
   for (j=0;j<e;j++) {
       res=addW(p,x[j],y[j],x[j],y[j],x[j+1],y[j+1],n,u[0],v[0]);
       if (res) return(1); }
   return(0);

   Uses one inversion and 6*e multiplications for e>1 (3 muls for e=1)
*/
int
addWn (mpz_t p, point *X, mpz_t n, long e)
{
  long j;
  point *u, *v;

  if (e == 1)
    return addW (p, X[0].x, X[0].y, X[0].x, X[0].y, X[1].x, X[1].y, n, X[2].x,
                 X[2].y);

  u = X + (e + 1);
  v = X + (e + 1);

  mpz_sub (u[e-1].x, X[e].x, X[e-1].x);
  mpz_set (v[e-1].y, u[e-1].x);
  for (j = e - 2; j >= 0; j--)
    {
      mpz_sub (u[j].x, X[j+1].x, X[j].x);
      mpz_mulmod (v[j].y, u[j].x, v[j+1].y, n); /* v[j] = u[j]*u[j+1]*...*u[e-1] */
    }

  mpz_gcdext (p, v[e].y, NULL, v[0].y, n);
  if (mpz_cmp_ui (p, 1) != 0)
    return 1;

  for (j = 0; j < e; j++)
    {
      /* loop invariant: v[e] = 1/(u[j]*u[j+1]*...*u[e-1]) */
      if (j != e - 1)
        {
          mpz_mulmod (v[j+1].y, v[j+1].y, v[e].y, n);
          /* restore v[e] for next loop and make u[j] free */
          mpz_mulmod (v[e].y, v[e].y, u[j].x, n);
        }
      /* now v[j+1] = 1/(x[j+1]-x[j]) mod n */
      mpz_sub (p, X[j+1].y, X[j].y);
      mpz_mulmod (p, v[j+1].y, p, n);
      mpz_mul (u[j].x, p, p);
      mpz_sub (u[j].x, u[j].x, X[j].x);
      mpz_sub (X[j].x, u[j].x, X[j+1].x);
      mpz_mod (X[j].x, X[j].x, n);
      mpz_sub (u[j].x, X[j+1].x, X[j].x);
      mpz_mul (u[j].x, u[j].x, p);
      mpz_sub (X[j].y, u[j].x, X[j+1].y);
      mpz_mod (X[j].y, X[j].y, n);
    }

  return 0;
}

#if 0
/* (x,y) can be identical to (x1,y1) */
int subW(p,x,y,x1,y1,x2,y2,n,u,v) mpz_t p,x,y,x1,y1,x2,y2,n,u,v;
{
  mpz_sub(u,x1,x2);
  mpz_gcdext(p,v,NULL,u,initial_n); if (mpz_cmp_ui(p,1)!=0) return(1);
  mpz_add(p,y1,y2); mpz_mulmod(p,v,p,n);
  mpz_mul(u,p,p); mpz_sub(u,u,x1); mpz_sub(v,u,x2); mod(v,v,n);
  mpz_sub(u,x1,v); mpz_mul(u,u,p); mpz_sub(y,u,y1); mod(y,y,n);
  mpz_set(x,v);
  mul+=3; gcdexts++;
  return(0);
}
#endif
