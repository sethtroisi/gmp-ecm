/* Elliptic Curve Method implementation: stage 2 routines.

  Copyright 2002, 2003 Alexander Kruppa and Paul Zimmermann.

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

int duplicateW (mpz_t, mpres_t, mpres_t, mpres_t, mpres_t, mpmod_t, 
                mpres_t, mpres_t, mpres_t);
int addW (mpz_t, mpres_t, mpres_t, mpres_t, mpres_t, mpres_t, mpres_t, 
          mpmod_t, mpres_t, mpres_t);
int multiplyW2 (mpz_t, mpres_t, mpres_t, mpres_t, mpres_t, mpz_t, mpmod_t, 
                mpres_t, mpres_t, mpres_t);
int addWn (mpz_t, point *, mpmod_t, long);


#define PTR(x) ((x)->_mp_d)
#define getbit(x,i) (PTR(x)[i/mp_bits_per_limb] & ((mp_limb_t)1<<(i%mp_bits_per_limb)))

/*
  (x1:y1) <- 2*(x:y) where (x:y) can be identical to (x1:y1).
  a is the Weierstrass parameter, u and v are auxiliary variables.
*/
int
duplicateW (mpz_t p, mpres_t x1, mpres_t y1, mpres_t x, mpres_t y, mpmod_t n, 
            mpres_t a, mpres_t u, mpres_t v)
{
  mpres_add (u, y, y, n);
  if (!mpres_invert (v, u, n))
    {
      mpres_gcd (p, u, n);
      return 1;
    }
  mpres_mul (u, x, x, n);
  mpres_mul_ui (u, u, 3, n);
  mpres_add (u, u, a, n);
  mpres_mul (p, u, v, n);
  mpres_mul (u, p, p, n);
  mpres_add (v, x, x, n);
  mpres_sub (u, u, v, n);
  mpres_sub (v, x, u, n);
  mpres_mul (v, v, p, n);
  mpres_sub (y1, v, y, n);
  mpres_set (x1, u, n); /* can we avoid this mpz_set (using mpz_swap perhaps)? */

  return 0;
}

/* Computes (x:y) <- (x1:y1) + (x2:y2).
   Returns non-zero iff a factor is found (then it is stored in p).
   n is the number to factor.
   u, v are auxiliary variables.
*/
int
addW (mpz_t p, mpres_t x, mpres_t y, mpres_t x1, mpres_t y1, mpres_t x2, 
      mpres_t y2, mpmod_t n, mpres_t u, mpres_t v)
{
  mpres_sub (u, x2, x1, n);
  if (!mpres_invert (v, u, n))
    {
      mpres_gcd (p, u, n);
      return 1;
    }
  mpres_sub (p, y2, y1, n);
  mpres_mul (p, v, p, n);
  mpres_mul (u, p, p, n);
  mpres_sub (u, u, x1, n);
  mpres_sub (v, u, x2, n);
  mpres_sub (u, x1, v, n);
  mpres_mul (u, u, p, n);
  mpres_sub (y, u, y1, n);
  mpres_set (x, v, n); /* can we avoid this copy? */

  return 0;
}

/* (x1:y1) <- q*(x:y) where q is a large integer */
int
multiplyW2 (mpz_t p, mpres_t x1, mpres_t y1, mpres_t x, mpres_t y, mpz_t q, 
            mpmod_t n, mpres_t a, mpres_t u, mpres_t v)
{
  long j;
  int sign_q;

  sign_q = mpz_sgn (q);

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
  if (getbit (q, j) && addW (p, x1, y1, x1, y1, x, y, n, u, v))
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
addWn (mpz_t p, point *X, mpmod_t n, long e)
{
  long j;
  point *u, *v;

  if (e == 1)
    return addW (p, X[0].x, X[0].y, X[0].x, X[0].y, X[1].x, X[1].y, n, X[2].x,
                 X[2].y);

  u = X + (e + 1);
  v = X + (e + 1);

  mpres_sub (u[e-1].x, X[e].x, X[e-1].x, n);
  mpres_set (v[e-1].y, u[e-1].x, n);
  for (j = e - 2; j >= 0; j--)
    {
      mpres_sub (u[j].x, X[j+1].x, X[j].x, n);
      mpres_mul (v[j].y, u[j].x, v[j+1].y, n); /* v[j] = u[j]*u[j+1]*...*u[e-1] */
    }

  if (!mpres_invert (v[e].y, v[0].y, n))
    {
      mpres_gcd (p, v[0].y, n);
      return 1;
    }

  for (j = 0; j < e; j++)
    {
      /* loop invariant: v[e] = 1/(u[j]*u[j+1]*...*u[e-1]) */
      if (j != e - 1)
        {
          mpres_mul (v[j+1].y, v[j+1].y, v[e].y, n);
          /* restore v[e] for next loop and make u[j] free */
          mpres_mul (v[e].y, v[e].y, u[j].x, n);
        }
      /* now v[j+1] = 1/(x[j+1]-x[j]) mod n */
      mpres_sub (p, X[j+1].y, X[j].y, n);
      mpres_mul (p, v[j+1].y, p, n);
      mpres_mul (u[j].x, p, p, n);
      mpres_sub (u[j].x, u[j].x, X[j].x, n);
      mpres_sub (X[j].x, u[j].x, X[j+1].x, n);
      mpres_sub (u[j].x, X[j+1].x, X[j].x, n);
      mpres_mul (u[j].x, u[j].x, p, n);
      mpres_sub (X[j].y, u[j].x, X[j+1].y, n);
    }

  return 0;
}

/* puts in F[0..dF-1] the successive values of 

   (j^S) * P where P=(s : : 1) is a point on the elliptic curve

   for 0 < j = 1 mod 7 < d, j and d coprime.
   Returns non-zero iff a factor was found (then stored in f).
*/

int
ecm_rootsF (mpz_t f, listz_t F, unsigned int d, curve *s, listz_t t,
        int S, mpmod_t modulus, int verbose)
{
  unsigned int i, j;
  int st, st2;
  point *fd;
  int youpi = 0, dickson_a;
  listz_t coeffs;
  mpres_t u, v;
  
  st = cputime ();

  /* If S < 0, use degree |S| Dickson poly, otherwise use x^S */
  dickson_a = (S < 0) ? -1 : 0;
  S = abs (S);

  mpres_get_z (F[0], s->x, modulus); /* (1*P)=P for ECM */
  i = 1;

  if (d > 7)
    {
      st2 = cputime ();

      mpres_init (u, modulus);
      mpres_init (v, modulus);

      coeffs = init_list (S + 1);
      fin_diff_coeff (coeffs, 7, 6, S, dickson_a);

      fd = (point *) malloc ((2 * S + 2) * sizeof (point));

      for (j = 0; j <= (unsigned) S && youpi == 0; j++)
        {
          mpres_init (fd[j].x, modulus);
          mpres_init (fd[j].y, modulus);
          youpi = multiplyW2 (f, fd[j].x, fd[j].y, s->x, s->y, coeffs[j], 
                              modulus, s->A, u, v);
          if (youpi && verbose >= 2)
            printf ("Found factor while computing fd[%d]\n", j);
        }
      
      clear_list (coeffs, S + 1);

      /* Allocate workspace for addWn */
      for ( ; j < 2 * (unsigned) S + 2; j++)
        {
          mpres_init (fd[j].x, modulus);
          mpres_init (fd[j].y, modulus);
        }

      if (verbose >= 2)
        printf ("Initializing table of differences for F took %dms\n", 
                cputime () - st2);

      j = 7;
      while (j < d && youpi == 0)
        {
          if (gcd (j, d) == 1)
            mpres_get_z (F[i++], fd[0].x, modulus);

          youpi = addWn (f, fd, modulus, S);

          if (youpi && verbose >= 2)
            printf ("Found factor while computing F[%d]\n", i);

          j += 6;
        }

      for (j = 0; j < 2 * (unsigned) S + 2; j++)
        {
          mpres_clear (fd[j].x, modulus);
          mpres_clear (fd[j].y, modulus);
        }
      free (fd);

      mpres_clear (u, modulus);
      mpres_clear (v, modulus);
    }

  if (youpi)
    return 1;
  
  if (verbose >= 2)
    printf ("Computing roots of F took %dms\n", cputime () - st);

  return 0;
}

/* Perform the neccessary initialisation to allow computation of
   
     Dickson_{S, a}(s+n*d) * P , where P is a point on the elliptic curve
   
   for successive n, where Dickson_{S, a} is the degree S Dickson
   polynomial with parameter a. For a == 0, Dickson_{S, a} (x) = x^S.
   
   If a factor is found during the initialisation, NULL is returned and the
     factor in f.
*/


point *
ecm_rootsG_init (mpz_t f, curve *X, unsigned int s, unsigned int d, 
                 int S, mpmod_t modulus, int verbose)
{
  unsigned int k;
  mpres_t u, v;
  listz_t coeffs;
  point *fd;
  int youpi = 0;
  int dickson_a;
  
  /* If S < 0, use degree |S| Dickson poly, otherwise use x^S */
  dickson_a = (S < 0) ? -1 : 0;
  S = abs (S);
  
  coeffs = init_list (S + 1);
  
  fin_diff_coeff (coeffs, s, d, S, dickson_a);
  
  fd = (point *) malloc ((2 * S + 2) * sizeof (point));
  
  mpres_init (u, modulus);
  mpres_init (v, modulus);
  
  for (k = 0; k <= (unsigned) S && youpi == 0; k++)
    {
      mpres_init (fd[k].x, modulus);
      mpres_init (fd[k].y, modulus);
      youpi = multiplyW2 (f, fd[k].x, fd[k].y, X->x, X->y, coeffs[k], modulus, 
                          X->A, u, v);
    }
  
  clear_list (coeffs, S + 1);
  mpres_clear (v, modulus);
  mpres_clear (u, modulus);
  
  if (youpi) 
    {
      unsigned int i;
      
      if (verbose >= 2)
        printf ("Found factor while computing fd[%d]\n", k);

      /* fd[0 .. k-1] have been initialized */

      for (i = 0; i < k; i++)
        {
          mpres_clear (fd[k].x, modulus);
          mpres_clear (fd[k].y, modulus);
        }
      free (fd);
      
      /* Signal that a factor was found */
      return NULL;
    }
    
  /* Allocate workspace for addWn */
  for ( ; k < 2 * (unsigned) S + 2; k++)
    {
      mpres_init (fd[k].x, modulus);
      mpres_init (fd[k].y, modulus);
    }
  
  return fd;
}

void 
ecm_rootsG_clear (point *fd, int S, mpmod_t modulus)
{
  unsigned int k;
  
  S = abs (S);
  for (k = 0; k < 2 * (unsigned) S + 2; k++)
    {
      mpres_clear (fd[k].x, modulus);
      mpres_clear (fd[k].y, modulus);
    }
  
  free (fd);
}

/* Puts in G the successive values of

     Dickson_{S, a}(s+j*k) P
    
     where P is a point on the elliptic curve,
     1<= j <= d, k is the 'd' value from ecm_rootsG_init()
     and s is the 's' value of ecm_rootsG_init() or where a previous
     call to ecm_rootsG has left off.

   Returns non-zero iff a factor was found (then stored in f).
*/

int 
ecm_rootsG (mpz_t f, listz_t G, unsigned int d, point *fd, listz_t t, 
            int S, mpmod_t modulus, int verbose)
{
  unsigned int i;
  int youpi = 0;
  
  S = abs (S);
  for (i = 0; i < d && youpi == 0; i++)
    {
      mpres_get_z (G[i], fd[0].x, modulus);
      youpi = addWn (f, fd, modulus, S);

      if (youpi && verbose >= 2)
        printf ("Found factor while computing G[%d]\n", i);
    }
  
  return youpi;
}
