/* Elliptic Curve Method implementation: stage 2 routines.

  Copyright 2002, 2003, 2004, 2005 Paul Zimmermann and Alexander Kruppa.

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
  the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
  MA 02111-1307, USA.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gmp.h"
#include "ecm.h"
#include "ecm-impl.h"

#if WANT_ASSERT
#include <assert.h>
#define ASSERT(expr)   assert (expr)
#else
#define ASSERT(expr)   do {} while (0)
#endif

#define min(a,b) (((a)<(b))?(a):(b))

/* R_i <- q_i * S, 0 <= i < n, where q_i are large positive integers, S 
   is a point on an elliptic curve. Uses max(bits in q_i) modular 
   inversions (one less if max(q_i) is a power of 2).
   Needs up to n+2 cells in T.
   Returns factor found, not found, error.
*/

static int
multiplyW2n (mpz_t p, point *R, curve *S, mpz_t *q, unsigned int n, 
              mpmod_t modulus, mpres_t u, mpres_t v, mpres_t *T,
              unsigned long *tot_muls, unsigned long *tot_gcds, FILE *es)
{
  unsigned int i, maxbit, k; /* k is the number of values to batch invert */
  unsigned int l, t, muls = 0, gcds = 0;
#ifdef WANT_EXPCOST
  unsigned int hamweight = 0;
#endif
  int youpi = ECM_NO_FACTOR_FOUND;
  mpz_t flag; /* Used as bit field, keeps track of which R[i] contain partial results */
  point s;    /* 2^t * S */

  if (n == 0)
    return ECM_NO_FACTOR_FOUND;
  
  /* Is S the neutral element ? */
  if (mpres_is_zero (S->x, modulus) && mpres_is_zero (S->y, modulus))
    {
      for (i = 0; i < n; i++)
        {
          mpres_set (R[i].x, S->x, modulus);
          mpres_set (R[i].y, S->y, modulus);
        }
      return ECM_NO_FACTOR_FOUND;
    }
  
  mpz_init (flag);
  mpres_init (s.x, modulus);
  mpres_init (s.y, modulus);
  mpres_set (s.x, S->x, modulus);
  mpres_set (s.y, S->y, modulus);

  /* Set maxbit to index of highest set bit among all the q[i] */
  /* Index of highest bit of q is sizeinbase(q, 2) - 1 */
  maxbit = 0;
  for (i = 0; i < n; i++)
    {
      if (mpz_sgn (q[i]) < 0)
        {
          fprintf (es, "multiplyW2n: multiplicand q[%d] < 0, negatives not supported\n", i);
	  youpi = ECM_ERROR;
	  goto clear_w2;
        }
      /* Multiplier == 0? Then set result to neutral element */
      if (mpz_sgn (q[i]) == 0)
        {
           mpres_set_ui (R[i].x, 0, modulus);
           mpres_set_ui (R[i].y, 0, modulus);
        }
#ifdef WANT_EXPCOST
      else
        hamweight += mpz_popcount (q[i]) - 1;
#endif
      if ((t = mpz_sizeinbase (q[i], 2) - 1) > maxbit)
          maxbit = t;
    }

#ifdef WANT_EXPCOST
  fprintf (es, "Expecting %d multiplications and %d extgcds\n", 
          4 * (maxbit) + 6 * hamweight - 3, maxbit + 1);      /* maxbit is floor(log_2(max(q_i))) */
#endif

  for (t = 0; t <= maxbit && !youpi; t++)    /* Examine t-th bit of the q[i] */
    {
      /* See which values need inverting and put them into T[]. Keep number
         of those values in k */
      k = 0;
      
      /* Will we have to double s at the end of this pass? If yes,  
         schedule 2*s.y for inverting */
      if (t < maxbit)
        mpres_add (T[k++], s.y, s.y, modulus);
      
      for (i = 0; i < n && !youpi; i++) 
        if (mpz_tstbit (q[i], t))       /* If q[i] & (1<<t), we'll add s to R[i] */
          if (mpz_tstbit (flag, i))     /* Does R[i] contain a partial result yet ? */
            {                           /* If Yes: need actual point addition so */
              mpres_sub (T[k], s.x, R[i].x, modulus); /* schedule (s.x-R[i].x) for inverting */
              if (k > 0)
                mpres_mul (T[k], T[k], T[k - 1], modulus);
              k++;
            }                           /* If No: we'll simply set R[i] to s later on, nothing tbd here */
      
      /* So there are k values in need of inverting, call them v[m], 0 <= m < k. */      
      /* Here T[m], 0 <= m < k, contains v[0]*...*v[m] */
      
      /* Put inverse of the product of all scheduled values in T[k]*/
      if (k > 0)
        {
          muls += 3 * (k - 1);
          gcds++;
          if (!mpres_invert (T[k], T[k - 1], modulus))
            {
              /* If a factor was found, put factor in p, 
                 flag success and bail out of loop */
              mpres_gcd (p, T[k - 1], modulus);
              youpi = ECM_FACTOR_FOUND;
              break;
            }
        }
      
      /* T[k] now contains 1/(v[0]*...*v[k - 1]), 
         T[m], 0 <= m < k, still contain v[0]*...*v[m] */
      
      l = k - 1;

      for (i = n; i-- > 0; ) /* Go through the R[i] again, backwards */
        if (mpz_tstbit (q[i], t))
          {
            if (mpz_tstbit (flag, i))
              {
                /* T[k] contains 1/(v[0]*...*v[l]) */
                if (l > 0) /* need to separate the values */
                  {
                    /* T[l - 1] has v[0]*...*v[l-1] */
                    mpres_mul (T[l], T[l - 1], T[k], modulus); /* So T[l] now has 1/v[l] == 1/(s.x - R[i].x) */
                    mpres_sub (u, s.x, R[i].x, modulus);
                    mpres_mul (T[k], T[k], u, modulus);        /* T[k] now has 1/(v[0]*...*v[l - 1]) */
                  }
                else
                  {
                    /* T[k] contains 1/v[0] */
                    mpres_set (T[0], T[k], modulus); 
                  }
                
                /* 1/(s.x - R[i].x) is in T[l] */
#ifdef WANT_ASSERT
                mpres_sub (u, s.x, R[i].x, modulus);
                mpres_mul (u, u, T[l], modulus);
                mpres_get_z (p, u, modulus);
                mpz_mod (p, p, modulus->orig_modulus);
                if (mpz_cmp_ui (p, 1) != 0) 
                  gmp_fprintf (es, "Error, (s.x - R[%d].x) * T[%d] == %Zd\n", i, l, T[l - 1]);
#endif
                
                mpres_sub (u, s.y, R[i].y, modulus);   /* U    = y2 - y1 */
                mpres_mul (T[l], T[l], u, modulus);    /* T[l] = (y2-y1)/(x2-x1) = lambda */
                mpres_mul (u, T[l], T[l], modulus);    /* U    = lambda^2 */
                mpres_sub (u, u, R[i].x, modulus);     /* U    = lambda^2 - x1 */
                mpres_sub (R[i].x, u, s.x, modulus);   /* x3   = lambda^2 - x1 - x2 */
                mpres_sub (u, s.x, R[i].x, modulus);   /* U    = x2 - x3 */
                mpres_mul (u, u, T[l], modulus);       /* U    = lambda*(x2 - x3) */
                mpres_sub (R[i].y, u, s.y, modulus);   /* y3   = lambda*(x2 - x3) - y2 */
                muls += 3;
                l--;
              }
            else /* R[i] does not contain a partial result. */
              {
                mpres_set (R[i].x, s.x, modulus);   /* Just set R[i] to s */
                mpres_set (R[i].y, s.y, modulus);
                mpz_setbit (flag, i);               /* and flag it as used */
              }
          }
      
      if (t < maxbit) /* Double s */
        { 
          ASSERT(l==0);
#ifdef WANT_ASSERT
          mpres_add (u, s.y, s.y, modulus);
          mpres_mul (u, u, T[k], modulus);
          mpres_get_z (p, u, modulus);
          mpz_mod (p, p, modulus->orig_modulus);
          if (mpz_cmp_ui (p, 1) != 0)
            gmp_fprintf (es, "Error, at t==%d, 2*s.y / (2*s.y) == %Zd\n", t, p);
#endif          

                                               /* 1/(2*s.y) is in T[k] */
          mpres_mul (u, s.x, s.x, modulus);    /* U = X^2 */
          mpres_mul_ui (u, u, 3, modulus);     /* U = 3*X^2 */
          mpres_add (u, u, S->A, modulus);     /* U = 3*X^2 + A */
          mpres_mul (T[k], T[k], u, modulus);  /* T = (3*X^2 + A) / (2*Y) = lambda */
          mpres_mul (u, T[k], T[k], modulus);  /* U = lambda^2 */
          mpres_sub (u, u, s.x, modulus);      /* U = lambda^2 - X */
          mpres_sub (u, u, s.x, modulus);      /* U = lambda^2 - 2*X = s.x' */
          mpres_sub (v, s.x, u, modulus);      /* V = s.x - s.x' */
          mpres_mul (v, v, T[k], modulus);     /* V = lambda*(s.x - s.x') */
          mpres_sub (s.y, v, s.y, modulus);    /* s.y' = lambda*(s.x - s.x') - s.y */
          mpres_set (s.x, u, modulus);
          muls += 4;
        }
    }

 clear_w2:
  mpres_clear (s.y, modulus);
  mpres_clear (s.x, modulus);
  mpz_clear (flag);

  if (tot_muls != NULL)
    *tot_muls += muls;
  if (tot_gcds != NULL)
    *tot_gcds += gcds;
  
  return youpi;
}


/* Input: Points X[0]..X[(n+1)*m-1]
   T is used for temporary values and needs to have (n-1)*m+1 entries.

   Performs the following loop with only one gcdext, using Montgomery's trick:
   for (i=0;i<m;i++)
     for (j=0;j<n;j++) {
         res=addW(p,x[j+n*i],y[j+n*i],x[j+n*i],y[j+n*i],x[j+1+n*i],y[j+1+n*i],
                  n,u[0],v[0]);
         if (res) return(1); }
   return(0);

   Uses one inversion and 6*n*m-3 multiplications for n*m > 0

   Return factor found or not (no error can occur here).
*/

static int
addWnm (mpz_t p, point *X, curve *S, mpmod_t modulus, unsigned int m, 
        unsigned int n, mpres_t *T, unsigned long *tot_muls, 
        unsigned long *tot_gcds)
{
  unsigned int k, l;
  int i, j;

  if (n == 0 || m == 0)
    return ECM_NO_FACTOR_FOUND;

  k = 0;
  for (i = m - 1; i >= 0; i--)    /* Go through the m different lists */
    for (j = n - 1; j >= 0; j--)  /* Go through each list backwards */
      {                           /* And prepare the values to be inverted */
        point *X1, *X2;
        X1 = X + i * (n + 1) + j;
        X2 = X + i * (n + 1) + j + 1;
        
        /* If either element is the neutral element, nothing tbd here */
        if ((mpres_is_zero (X1->x, modulus) && mpres_is_zero (X1->y, modulus)) ||
            (mpres_is_zero (X2->x, modulus) && mpres_is_zero (X2->y, modulus)))
          continue;
        
        mpres_sub (T[k], X2->x, X1->x, modulus); /* Schedule X2.x - X1.x */

        if (mpres_is_zero (T[k], modulus))  /* If both x-cordinates are identical */
          {
            /* Are the points identical? */
            mpres_sub (T[k], X2->y, X1->y, modulus);
            if (mpres_is_zero (T[k], modulus))
              {
                /* Yes, we need to double. Schedule 2*X[...].y */
                mpres_add (T[k], X1->y, X1->y, modulus); 
              }
            else
              continue; /* No, they are inverses. Nothing tbd here */
          }

        if (k > 0)
          mpres_mul (T[k], T[k], T[k - 1], modulus);
        k++;
      }

  /* v_m = X[i * (n + 1) + j] - X[i * (n + 1) + j + 1], 0 <= j < n,
     and m = i * n + j */
  /* Here T[m] = v_0 * ... * v_m, 0 <= m < k */

  if (k > 0 && !mpres_invert (T[k], T[k - 1], modulus))
    {
      mpres_gcd (p, T[k - 1], modulus);
      (*tot_muls) += m * n - 1;
      (*tot_gcds) ++;
      return ECM_FACTOR_FOUND;
    }

  /* T[k] = 1/(v_0 * ... * v_m), 0 <= m < k */

  l = k - 1;

  for (i = 0; (unsigned) i < m; i++)
    for (j = 0; (unsigned) j < n; j++)
      {
        point *X1, *X2;
        X1 = X + i * (n + 1) + j;
        X2 = X + i * (n + 1) + j + 1;
        
        /* Is X1 the neutral element? */
        if (mpres_is_zero (X1->x, modulus) && mpres_is_zero (X1->y, modulus))
          {
            /* Yes, set X1 to X2 */
            mpres_set (X1->x, X2->x, modulus);
            mpres_set (X1->y, X2->y, modulus);
            continue;
          }
        
        /* Is X2 the neutral element? If so, X1 stays the same */
        if (mpres_is_zero (X2->x, modulus) && mpres_is_zero (X2->y, modulus))
          continue;
        
        /* Are the x-coordinates identical? */
        mpres_sub (T[k + 1], X2->x, X1->x, modulus);
        if (mpres_is_zero (T[k + 1], modulus))
          {
            /* Are the points inverses of each other? */
            mpres_sub (T[k + 1], X2->y, X1->y, modulus);
            if (!mpres_is_zero (T[k + 1], modulus))
              {
                /* Yes. Set X1 to neutral element */
                mpres_set_ui (X1->x, 0, modulus);
                mpres_set_ui (X1->y, 0, modulus);
                continue;
              }
            /* No, we need to double. Restore T[k+1] */
            mpres_sub (T[k + 1], X2->x, X1->x, modulus);
          }

        if (l == 0)
          mpz_set (T[0], T[k]);
        else
          mpres_mul (T[l], T[k], T[l - 1], modulus); 
          /* T_l = 1/(v_0 * ... * v_l) * (v_0 * ... * v_{l-1}) = 1/v_l */


        if (mpres_is_zero (T[k + 1], modulus)) /* Identical points, so double X1 */
          {
            if (l > 0)
              {
                mpres_add (T[k + 1], X1->y, X1->y, modulus); /* T[k+1] = v_{l} */
                mpres_mul (T[k], T[k], T[k + 1], modulus);
                /* T_k = 1/(v_0 * ... * v_l) * v_l = 1/(v_0 * ... * v_{l-1}) */
              }
            
            mpres_mul (T[k + 1], X1->x, X1->x, modulus);
            mpres_mul_ui (T[k + 1], T[k + 1], 3, modulus);
            mpres_add (T[k + 1], T[k + 1], S->A, modulus);
            mpres_mul (T[l], T[k + 1], T[l], modulus); /* T[l] = lambda */
            mpres_mul (T[k + 1], T[l], T[l], modulus);       /* T1   = lambda^2 */
            mpres_sub (T[k + 1], T[k + 1], X1->x, modulus);  /* T1   = lambda^2 - x1 */
            mpres_sub (X1->x, T[k + 1], X2->x, modulus);     /* X1.x = lambda^2 - x1 - x2 = x3 */
            mpres_sub (T[k + 1], X2->x, X1->x, modulus);     /* T1   = x2 - x3 */
            mpres_mul (T[k + 1], T[k + 1], T[l], modulus);   /* T1   = lambda*(x2 - x3) */
            mpres_sub (X1->y, T[k + 1], X2->y, modulus);     /* Y1   = lambda*(x2 - x3) - y2 = y3 */
          }
        else
          {
            if (l > 0)
              {
                mpres_mul (T[k], T[k], T[k + 1], modulus);
                /* T_k = 1/(v_0 * ... * v_l) * v_l = 1/(v_0 * ... * v_{l-1}) */
              }

            mpres_sub (T[k + 1], X2->y, X1->y, modulus);     /* T1   = y2 - y1 */
            mpres_mul (T[l], T[l], T[k + 1], modulus);       /* Tl   = (y2 - y1) / (x2 - x1) = lambda */
            mpres_mul (T[k + 1], T[l], T[l], modulus);       /* T1   = lambda^2 */
            mpres_sub (T[k + 1], T[k + 1], X1->x, modulus);  /* T1   = lambda^2 - x1 */
            mpres_sub (X1->x, T[k + 1], X2->x, modulus);     /* X1.x = lambda^2 - x1 - x2 = x3 */
            mpres_sub (T[k + 1], X2->x, X1->x, modulus);     /* T1   = x2 - x3 */
            mpres_mul (T[k + 1], T[k + 1], T[l], modulus);   /* T1   = lambda*(x2 - x3) */
            mpres_sub (X1->y, T[k + 1], X2->y, modulus);     /* Y1   = lambda*(x2 - x3) - y2 = y3 */
          }
        
        l--;
      }

  if (tot_muls != NULL)
    (*tot_muls) += 6 * m * n - 3;
  if (tot_gcds != NULL)
    (*tot_gcds) ++;

  return ECM_NO_FACTOR_FOUND;
}

/* puts in F[0..dF-1] the successive values of 

   Dickson_{S, a} (j * d2) * s  where s is a point on the elliptic curve

   for j == 1 mod 6, j and d1 coprime.
   Returns non-zero iff a factor was found (then stored in f)
   or an error occurred.
*/

int
ecm_rootsF (mpz_t f, listz_t F, unsigned int d1, unsigned int d2, 
            unsigned int dF, curve *s, int S, mpmod_t modulus, int verbose,
	    FILE *os, FILE *es)
{
  unsigned int i;
  unsigned long muls = 0, gcds = 0;
  int st;
  int youpi = ECM_NO_FACTOR_FOUND;
  listz_t coeffs;
  ecm_roots_state state;
  
  if (dF == 0)
    return ECM_NO_FACTOR_FOUND;

  if (verbose >= 2)
    st = cputime ();

  /* Relative cost of point add during init and computing roots assumed =1 */
  init_roots_state (&state, S, d1, d2, 1.0);

  if (verbose >= 3)
     fprintf (os, "ecm_rootsF: state: nr = %d, dsieve = %d, size_fd = %d, S = %d, "
              "dickson_a = %d\n", state.nr, state.dsieve, state.size_fd, 
              state.S, state.dickson_a);

  /* Init finite differences tables */
  coeffs = init_progression_coeffs (0.0, state.dsieve, d2, 1, 6, state.S, 
                                    state.dickson_a, es);

  if (coeffs == NULL) /* error */
    {
      youpi = ECM_ERROR;
      goto clear;
    }

  /* The highest coefficient is the same for all progressions, so set them
     to one for all but the first progression, later we copy the point */
  for (i = state.S + 1; i < state.size_fd; i += state.S + 1)
    mpz_set_ui (coeffs[i + state.S], 1);

  /* Allocate memory for fd[] and T[] */

  state.fd = (point *) xmalloc (state.size_fd * sizeof (point));
  for (i = 0; i < state.size_fd; i++)
    {
     if (verbose >= 4)
       gmp_fprintf (os, "ecm_rootsF: coeffs[%d] = %Zd\n", i, coeffs[i]);
      mpres_init (state.fd[i].x, modulus);
      mpres_init (state.fd[i].y, modulus);
    }

  state.T = (mpres_t *) xmalloc ((state.size_fd + 4) * sizeof (mpres_t));
  for (i = 0 ; i < state.size_fd + 4; i++)
    mpres_init (state.T[i], modulus);

  /* Multiply fd[] = s * coeffs[] */

  youpi = multiplyW2n (f, state.fd, s, coeffs, state.size_fd, modulus, 
                       state.T[0], state.T[1], state.T + 2, &muls, &gcds, es);
  if (youpi == ECM_FACTOR_FOUND && verbose >= 2)
    fprintf (os, "Found factor while computing coeff[] * X\n");  

  if (youpi == ECM_ERROR)
    goto clear;

  /* Copy the point corresponding to the highest coefficient of the first 
     progression to the other progressions */
  for (i = state.S + 1; i < state.size_fd; i += state.S + 1)
    {
      mpres_set (state.fd[i + state.S].x, state.fd[state.S].x, modulus);
      mpres_set (state.fd[i + state.S].y, state.fd[state.S].y, modulus);
    }

  clear_list (coeffs, state.size_fd);
  coeffs = NULL;

  if (verbose >= 2)
    {
      int st1 = cputime ();
      fprintf (os, "Initializing tables of differences for F took %dms",
	       st1 - st);
      if (verbose > 2)
        fprintf (os, ", %lu muls and %lu extgcds", muls, gcds);
      fprintf (os, "\n");
      st = st1;
      muls = 0;
      gcds = 0;
    }

  /* Now for the actual calculation of the roots. */

  for (i = 0; i < dF && !youpi;)
    {
      /* Is this a rsieve value where we computed Dickson(j * d2) * X? */
      if (gcd (state.rsieve, state.dsieve) == 1) 
        {
          /* Did we use every progression since the last update? */
          if (state.next == state.nr)
            {
              /* Yes, time to update again */
              youpi = addWnm (f, state.fd, s, modulus, state.nr, state.S, 
                              state.T, &muls, &gcds);
	      ASSERT(youpi != ECM_ERROR); /* no error can occur in addWnm */
              state.next = 0;
              if (youpi == ECM_FACTOR_FOUND && verbose >= 2)
                fprintf (os, "Found factor while computing roots of F\n");
            }
          
          /* Is this a j value where we want Dickson(j * d2) * X as a root? */
          if (gcd (state.rsieve, d1) == 1) 
            mpres_get_z (F[i++], state.fd[state.next * (state.S + 1)].x, 
                         modulus);

          state.next ++;
        }
      state.rsieve += 6;
    }

 clear:
  for (i = 0 ; i < state.size_fd + 4; i++)
    mpres_clear (state.T[i], modulus);
  free (state.T);
  
  for (i = 0; i < state.size_fd; i++)
    {
      mpres_clear (state.fd[i].x, modulus);
      mpres_clear (state.fd[i].y, modulus);
    }
  free (state.fd);

  if (youpi)
    return youpi; /* error or factor found */
  
  if (verbose >= 2)
    {
      fprintf (os, "Computing roots of F took %dms", cputime () - st);
      if (verbose > 2)
        fprintf (os, ", %ld muls and %ld extgcds", muls, gcds);
      fprintf (os, "\n");
    }

  return ECM_NO_FACTOR_FOUND;
}

/* Perform the necessary initialization to allow computation of
   
     Dickson_{S, a}(s+n*d) * P , where P is a point on the elliptic curve
   
   for successive n, where Dickson_{S, a} is the degree S Dickson
   polynomial with parameter a. For a == 0, Dickson_{S, a} (x) = x^S.
   
   If a factor is found during the initialisation, NULL is returned and the
   factor in f. If an error occurred, NULL is returned and f is -1.
*/
ecm_roots_state *
ecm_rootsG_init (mpz_t f, curve *X, double s, unsigned int d1, unsigned int d2,
                 unsigned int dF, unsigned int blocks, int S, mpmod_t modulus, 
                 int verbose, FILE *os, FILE *es)
{
  unsigned int k, lenT, phid2;
  unsigned long muls = 0, gcds = 0;
  listz_t coeffs;
  ecm_roots_state *state;
  int youpi = 0;
  int dickson_a;
  unsigned int T_inv;
  double bestnr;
  int st = 0;

  ASSERT (gcd (d1, d2) == 1);

  if (verbose >= 2)
    st = cputime ();
  
  /* If S < 0, use degree |S| Dickson poly, otherwise use x^S */
  dickson_a = (S < 0) ? -1 : 0;
  S = abs (S);

  /* Estimate the cost of a modular inversion (in unit of time per 
     modular multiplication) */
  if (modulus->repr == MOD_BASE2)
    T_inv = 18;
  else
    T_inv = 6;
  
  /* Guesstimate a value for the number of disjoint progressions to use */
  bestnr = -(4. + T_inv) + sqrt(12. * (double) dF * (double) blocks * 
        (T_inv - 3.) * log (2. * d1) / log (2.) - (4. + T_inv) * (4. + T_inv));
  bestnr /= 6. * (double) S * log (2. * d1) / log (2);
  
  if (verbose >= 4)
    fprintf (os, "ecm_rootsG_init: bestnr = %f\n", bestnr);
  
  state = (ecm_roots_state *) xmalloc (sizeof (ecm_roots_state));

  if (bestnr < 1.)
    state->nr = 1;
  else
    state->nr = (unsigned int) (bestnr + .5);

  phid2 = phi (d2);

  /* Round up state->nr to multiple of phi(d2) */
  if (phid2 > 1)
    state->nr = ((state->nr + (phid2 - 1)) / phid2) * phid2;

  state->S = S;
  state->size_fd = state->nr * (state->S + 1);

  if (verbose >= 3)
    fprintf (os, "ecm_rootsG_init: s=%f, d1=%u, d2=%d, dF=%d, blocks=%d, S=%u, T_inv = %d, nr=%d\n", 
	     s, d1, d2, dF, blocks, S, T_inv, state->nr);
  
  state->X = X;
  state->next = 0;
  state->dsieve = 1; /* We only init progressions coprime to d2, so nothing to be skipped */
  state->rsieve = 1;

  coeffs = init_progression_coeffs (s, d2, d1, state->nr / phid2, 1, S, 
                                    dickson_a, es);

  if (coeffs == NULL) /* error */
    {
      free (state);
      mpz_set_si (f, -1);
      return NULL;
    }

  state->fd = (point *) xmalloc (state->size_fd * sizeof (point));
  for (k = 0; k < state->size_fd; k++)
    {
      mpres_init (state->fd[k].x, modulus);
      mpres_init (state->fd[k].y, modulus);
    }
  
  lenT = state->size_fd + 4;
  state->T = (mpres_t *) xmalloc (lenT * sizeof (mpres_t));
  for (k = 0; k < lenT; k++)
    mpres_init (state->T[k], modulus);

  for (k = S + 1; k < state->size_fd; k += S + 1)
     mpz_set_ui (coeffs[k + S], 1);

  if (verbose >= 4)
    for (k = 0; k < state->size_fd; k++)
      gmp_fprintf (os, "ecm_rootsG_init: coeffs[%d] == %Zd\n", k, coeffs[k]);

  youpi = multiplyW2n (f, state->fd, X, coeffs, state->size_fd, modulus, 
                     state->T[0], state->T[1], state->T + 2, &muls, &gcds, es);
  if (youpi == ECM_ERROR)
    mpz_set_si (f, -1); /* fall through */

  for (k = S + 1; k < state->size_fd; k += S + 1)
    {
      mpres_set (state->fd[k + S].x, state->fd[S].x, modulus);
      mpres_set (state->fd[k + S].y, state->fd[S].y, modulus);
    }
  
  clear_list (coeffs, state->size_fd);
  coeffs = NULL;
  
  if (youpi != ECM_NO_FACTOR_FOUND) /* factor found or error */
    {
      if (youpi == ECM_FACTOR_FOUND && verbose >= 2)
        fprintf (os, "Found factor while computing fd[]\n");

      ecm_rootsG_clear (state, S, modulus);
      
      /* Signal that a factor was found, or an error occurred (f=-1) */
      state = NULL;
    }
  else
    {
      if (verbose >= 2)
        {
          st = cputime () - st;
          fprintf (os, "Initializing table of differences for G took %dms", st);
          if (verbose > 2)
            fprintf (os, ", %lu muls and %lu extgcds", muls, gcds);
          fprintf (os, "\n");
        }
    }
  
  return state;
}

void 
ecm_rootsG_clear (ecm_roots_state *state, ATTRIBUTE_UNUSED int S, 
                  ATTRIBUTE_UNUSED mpmod_t modulus)
{
  unsigned int k;
  
  for (k = 0; k < state->size_fd; k++)
    {
      mpres_clear (state->fd[k].x, modulus);
      mpres_clear (state->fd[k].y, modulus);
    }
  free (state->fd);
  
  for (k = 0; k < state->size_fd + 4; k++)
    mpres_clear (state->T[k], modulus);
  free (state->T);
  
  free (state);
}

/* Puts in G the successive values of

     Dickson_{S, a}(s+j*k) P
    
     where P is a point on the elliptic curve,
     0<= j <= dF-1, k is the 'd' value from ecm_rootsG_init()
     and s is the 's' value of ecm_rootsG_init() or where a previous
     call to ecm_rootsG has left off.

   Returns non-zero iff a factor was found (then stored in f).
   Cannot return an error.
*/

int 
ecm_rootsG (mpz_t f, listz_t G, unsigned int dF, ecm_roots_state *state, 
            mpmod_t modulus, int verbose, FILE *os)
{
  unsigned int i;
  unsigned long muls = 0, gcds = 0;
  int youpi = ECM_NO_FACTOR_FOUND, st;
  
  st = cputime ();
  
  for (i = 0; i < dF;)
    {
      /* Did we use every progression since the last update? */
      if (state->next == state->nr)
        {
          /* Yes, time to update again */
          youpi = addWnm (f, state->fd, state->X, modulus, state->nr, 
                          state->S, state->T, &muls, &gcds);
	  ASSERT(youpi != ECM_ERROR); /* no error can occur in addWnm */
          state->next = 0;
          
          if (youpi == ECM_FACTOR_FOUND)
            {
              if (verbose >= 2)
                fprintf (os, "Found factor while computing G[]\n");
              break;
            }
        }
      
      /* Is this a root we should skip? (Take only if gcd == 1) */
      if (gcd (state->rsieve, state->dsieve) == 1)
        mpres_get_z (G[i++], (state->fd + state->next * (state->S + 1))->x, 
                       modulus);

      state->next ++;
      state->rsieve ++;
    }
  
  if (verbose >= 2)
    {
      fprintf (os, "Computing roots of G took %dms", cputime () - st);
      if (verbose > 2)
        fprintf (os, ", %lu muls and %lu extgcds", muls, gcds);
      fprintf (os, "\n");
    }
  
  return youpi;
}
