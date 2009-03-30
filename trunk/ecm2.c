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
  the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
  MA 02110-1301, USA.
*/

#include <stdlib.h>
#include <math.h>
#include "ecm-impl.h"

/* R_i <- q_i * S, 0 <= i < n, where q_i are large integers, S is a point on
   an elliptic curve. Uses max(bits in q_i) modular inversions (one less if 
   max(q_i) is a power of 2). Needs up to n+2 cells in T.
   Returns whether factor was found or not found, factor goes into p. 
   No error can occur.
*/

static int
multiplyW2n (mpz_t p, point *R, curve *S, mpz_t *q, const unsigned int n, 
              mpmod_t modulus, mpres_t u, mpres_t v, mpres_t *T,
              unsigned long *tot_muls, unsigned long *tot_gcds)
{
  unsigned int i, maxbit, k; /* k is the number of values to batch invert */
  unsigned int l, t, muls = 0, gcds = 0;
#ifdef WANT_EXPCOST
  unsigned int hamweight = 0;
#endif
  int youpi = ECM_NO_FACTOR_FOUND;
  mpz_t flag; /* Used as bit field, keeps track of which R[i] contain partial results */
  point s;    /* 2^t * S */
  mpz_t signs; /* Used as bit field, i-th bit is set iff q[i]<0 */
#ifdef WANT_ASSERT
  mpz_t __dummy; /* used for local computations */
#endif

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
  
  MPZ_INIT2 (flag, n);
  MPZ_INIT2 (signs, n);
  mpres_init (s.x, modulus);
  mpres_init (s.y, modulus);
  mpres_set (s.x, S->x, modulus);
  mpres_set (s.y, S->y, modulus);

  /* Set maxbit to index of highest set bit among all the q[i] */
  /* Index of highest bit of q is sizeinbase(q, 2) - 1 */
  maxbit = 0;
  for (i = 0; i < n; i++)
    {
      /* We'll first compute positive multiples and change signs later */
      if (mpz_sgn (q[i]) < 0)
        {
          mpz_setbit (signs, i);;
          mpz_neg (q[i], q[i]);
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
  outputf (OUTPUT_ALWAYS, "Expecting %d multiplications and %d extgcds\n", 
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
              if (p != NULL)
                mpres_gcd (p, T[k - 1], modulus);
              youpi = ECM_FACTOR_FOUND_STEP2;
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
		mpz_init(__dummy);
                mpres_get_z (__dummy, u, modulus);
                mpz_mod (__dummy, __dummy, modulus->orig_modulus);
                if (mpz_cmp_ui (__dummy, 1) != 0) 
                  outputf (OUTPUT_ERROR, "Error, (s.x - R[%d].x) * T[%d] == "
                           "%Zd\n", i, l, __dummy);
		mpz_clear(__dummy);
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
	  mpz_init(__dummy);
          mpres_get_z (__dummy, u, modulus);
          mpz_mod (__dummy, __dummy, modulus->orig_modulus);
          if (mpz_cmp_ui (__dummy, 1) != 0)
            outputf (OUTPUT_ERROR, "Error, at t==%d, 2*s.y / (2*s.y) == %Zd\n", 
                     t, __dummy);
	  mpz_clear(__dummy);
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

  mpres_clear (s.y, modulus);
  mpres_clear (s.x, modulus);
  mpz_clear (flag);

  if (tot_muls != NULL)
    *tot_muls += muls;
  if (tot_gcds != NULL)
    *tot_gcds += gcds;
  
  /* Now take inverse points (negative y-coordinate) where q[i] was < 0 */
  for (i = 0; i < n; i++)
    if (mpz_tstbit (signs, i))
      {
        mpz_neg (R[i].y, R[i].y);
        mpz_neg (q[i], q[i]);
      }

  mpz_clear (signs);
  return youpi;
}


/* Input: Points X[0]..X[(n+1)*m-1]
   T is used for temporary values and needs to have (n-1)*m+1 entries.

   Performs the following loop with only one gcdext, using Montgomery's trick:
   for (i=0;i<m;i++)
     for (j=0;j<n;j++)
         (x[j+(n+1)*i] : y[j+(n+1)*i]) += (x[j+1+(n+1)*i] : y[j+1+(n+1)*i])

   Uses one inversion and 6*n*m-3 multiplications for n*m > 0.
   Processes neutral (zero), identical and negative points correctly.

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
            /* Are the points identical? Compare y coordinates: */
            mpres_sub (T[k], X2->y, X1->y, modulus);
            if (mpres_is_zero (T[k], modulus))
              {
                /* Yes, we need to double. Schedule 2*X[...].y */
                mpres_add (T[k], X1->y, X1->y, modulus); 
              }
            else /* No, they are inverses. Nothing tbd here */
              {
#ifdef WANT_ASSERT
                /* Check that the y coordinates are mutual negatives */
                mpres_add (T[k], X2->y, X1->y, modulus);
                ASSERT (mpres_is_zero (T[k], modulus));
#endif
                continue; 
              }
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
      if (p != NULL)
        mpres_gcd (p, T[k - 1], modulus);
      if (tot_muls != NULL)
        (*tot_muls) += m * n - 1;
      if (tot_gcds != NULL)
        (*tot_gcds) ++;
      return ECM_FACTOR_FOUND_STEP2;
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
ecm_rootsF (mpz_t f, listz_t F, root_params_t *root_params, 
            unsigned long dF, curve *s, mpmod_t modulus)
{
  unsigned long i;
  unsigned long muls = 0, gcds = 0;
  long st;
  int youpi = ECM_NO_FACTOR_FOUND;
  listz_t coeffs;
  ecm_roots_state_t state;
  progression_params_t *params = &state.params; /* for less typing */
  mpz_t t;
  
  if (dF == 0)
    return ECM_NO_FACTOR_FOUND;

  st = cputime ();

  /* Relative cost of point add during init and computing roots assumed =1 */
  init_roots_params (params, root_params->S, root_params->d1, root_params->d2, 
		     1.0);

  outputf (OUTPUT_DEVVERBOSE, "ecm_rootsF: state: nr = %d, dsieve = %d, "
	   "size_fd = %d, S = %d, dickson_a = %d\n", 
	   params->nr, params->dsieve, params->size_fd, params->S, 
	   params->dickson_a);

  /* Init finite differences tables */
  MPZ_INIT (t); /* t = 0 */
  coeffs = init_progression_coeffs (t, params->dsieve, root_params->d2, 
				    1, 6, params->S, params->dickson_a);
  mpz_clear (t);

  if (coeffs == NULL) /* error */
    {
      youpi = ECM_ERROR;
      goto clear;
    }

  /* The highest coefficient is the same for all progressions, so set them
     to one for all but the first progression, later we copy the point.
     FIXME: can we avoid the multiplication of those points in multiplyW2n()
     below?
  */
  for (i = params->S + 1; i < params->size_fd; i += params->S + 1)
    mpz_set_ui (coeffs[i + params->S], 1);

  /* Allocate memory for fd[] and T[] */

  state.fd = (point *) malloc (params->size_fd * sizeof (point));
  if (state.fd == NULL)
    {
      youpi = ECM_ERROR;
      goto exit_ecm_rootsF;
    }
  for (i = 0; i < params->size_fd; i++)
    {
      outputf (OUTPUT_TRACE, "ecm_rootsF: coeffs[%d] = %Zd\n", i, coeffs[i]);
      MEMORY_TAG;
      mpres_init (state.fd[i].x, modulus);
      MEMORY_TAG;
      mpres_init (state.fd[i].y, modulus);
      MEMORY_UNTAG;
    }

  state.T = (mpres_t *) malloc ((params->size_fd + 4) * sizeof (mpres_t));
  if (state.T == NULL)
    {
      youpi = ECM_ERROR;
      goto ecm_rootsF_clearfdi;
    }
  for (i = 0 ; i < params->size_fd + 4; i++)
    {
      MEMORY_TAG;
      mpres_init (state.T[i], modulus);
      MEMORY_UNTAG;
    }

  /* Multiply fd[] = s * coeffs[] */

  youpi = multiplyW2n (f, state.fd, s, coeffs, params->size_fd, modulus,
                       state.T[0], state.T[1], state.T + 2, &muls, &gcds);
  if (youpi == ECM_FACTOR_FOUND_STEP2)
    outputf (OUTPUT_VERBOSE, "Found factor while computing coeff[] * X\n");  

  if (youpi == ECM_ERROR)
    goto clear;

  /* Copy the point corresponding to the highest coefficient of the first 
     progression to the other progressions */
  for (i = params->S + 1; i < params->size_fd; i += params->S + 1)
    {
      mpres_set (state.fd[i + params->S].x, state.fd[params->S].x, modulus);
      mpres_set (state.fd[i + params->S].y, state.fd[params->S].y, modulus);
    }

  clear_list (coeffs, params->size_fd);
  coeffs = NULL;

  if (test_verbose (OUTPUT_VERBOSE))
    {
      unsigned int st1 = cputime ();
      outputf (OUTPUT_VERBOSE,
	       "Initializing tables of differences for F took %ldms",
	       elltime (st, st1));
      outputf (OUTPUT_DEVVERBOSE, ", %lu muls and %lu extgcds", muls, gcds);
      outputf (OUTPUT_VERBOSE, "\n");
      st = st1;
      muls = 0;
      gcds = 0;
    }

  /* Now for the actual calculation of the roots. */

  for (i = 0; i < dF && !youpi;)
    {
      /* Is this a rsieve value where we computed Dickson(j * d2) * X? */
      if (gcd ((unsigned long) params->rsieve, 
	       (unsigned long) params->dsieve) == 1UL) 
        {
          /* Did we use every progression since the last update? */
          if (params->next == params->nr)
            {
              /* Yes, time to update again */
              youpi = addWnm (f, state.fd, s, modulus, params->nr, params->S, 
                              state.T, &muls, &gcds);
	      ASSERT(youpi != ECM_ERROR); /* no error can occur in addWnm */
              params->next = 0;
              if (youpi == ECM_FACTOR_FOUND_STEP2)
                outputf (OUTPUT_VERBOSE,
			 "Found factor while computing roots of F\n");
            }
          
          /* Is this a j value where we want Dickson(j * d2) * X as a root? */
          if (gcd ((unsigned long) params->rsieve, root_params->d1) 
	      == 1UL) 
            mpres_get_z (F[i++], 
			 state.fd[params->next * (params->S + 1)].x, modulus);

          params->next ++;
        }
      params->rsieve += 6;
    }

 clear:
  for (i = 0 ; i < params->size_fd + 4; i++)
    mpres_clear (state.T[i], modulus);
  free (state.T);

 ecm_rootsF_clearfdi:
  for (i = 0; i < params->size_fd; i++)
    {
      mpres_clear (state.fd[i].x, modulus);
      mpres_clear (state.fd[i].y, modulus);
    }
  free (state.fd);

 exit_ecm_rootsF:
  if (youpi)
    return youpi; /* error or factor found */
  
  outputf (OUTPUT_VERBOSE, "Computing roots of F took %ldms",
	   elltime (st, cputime ()));
  outputf (OUTPUT_DEVVERBOSE, ", %ld muls and %ld extgcds", muls, gcds);
  outputf (OUTPUT_VERBOSE, "\n");

  return ECM_NO_FACTOR_FOUND;
}

/* Perform the necessary initialization to allow computation of
   
     Dickson_{S, a}(s+n*d) * P , where P is a point on the elliptic curve
   
   for successive n, where Dickson_{S, a} is the degree S Dickson
   polynomial with parameter a. For a == 0, Dickson_{S, a} (x) = x^S.
   
   If a factor is found during the initialisation, NULL is returned and the
   factor in f. If an error occurred, NULL is returned and f is -1.
*/

ecm_roots_state_t *
ecm_rootsG_init (mpz_t f, curve *X, root_params_t *root_params, 
                 unsigned long dF, unsigned long blocks, mpmod_t modulus)
{
  unsigned int k, phid2;
  unsigned long muls = 0, gcds = 0;
  listz_t coeffs;
  ecm_roots_state_t *state;
  progression_params_t *params; /* for less typing */
  int youpi = 0;
  unsigned int T_inv;
  double bestnr;
  long st = 0;

  ASSERT (gcd (root_params->d1, root_params->d2) == 1UL);

  if (test_verbose (OUTPUT_VERBOSE))
    st = cputime ();
  
  state = (ecm_roots_state_t *) malloc (sizeof (ecm_roots_state_t));
  if (state == NULL)
    {
      mpz_set_si (f, -1);
      return NULL;
    }
  params = &(state->params);

  /* If S < 0, use degree |S| Dickson poly, otherwise use x^S */
  params->dickson_a = (root_params->S < 0) ? -1 : 0;
  params->S = abs (root_params->S);

  /* Estimate the cost of a modular inversion (in unit of time per 
     modular multiplication) */
  if (modulus->repr == ECM_MOD_BASE2)
    T_inv = 18;
  else
    T_inv = 6;
  
  /* Guesstimate a value for the number of disjoint progressions to use */
  bestnr = -(4. + T_inv) + sqrt(12. * (double) dF * (double) blocks * 
        (T_inv - 3.) * log (2. * root_params->d1) / log (2.) - (4. + T_inv) * 
        (4. + T_inv));
  bestnr /= 6. * (double) (params->S) * log (2. * root_params->d1) / log (2.0);
  
  outputf (OUTPUT_TRACE, "ecm_rootsG_init: bestnr = %f\n", bestnr);
  
  if (bestnr < 1.)
    params->nr = 1;
  else
    params->nr = (unsigned int) (bestnr + .5);

  phid2 = eulerphi (root_params->d2);

  /* Round up params->nr to multiple of eulerphi(d2) */
  if (phid2 > 1)
    params->nr = ((params->nr + (phid2 - 1)) / phid2) * phid2;

  params->size_fd = params->nr * (params->S + 1);

  outputf (OUTPUT_DEVVERBOSE, "ecm_rootsG_init: i0=%Zd, d1=%lu, d2=%lu, "
           "dF=%lu, blocks=%lu, S=%u, T_inv = %d, nr=%d\n", root_params->i0, 
           root_params->d1, root_params->d2, dF, blocks, params->S, 
	   T_inv, params->nr);
  
  state->X = X;
  params->next = 0;
  params->dsieve = 1; /* We only init progressions coprime to d2, 
			       so nothing to be skipped */
  params->rsieve = 0;

  coeffs = init_progression_coeffs (root_params->i0, root_params->d2, 
				    root_params->d1, params->nr / phid2, 
				    1, params->S, params->dickson_a);

  if (coeffs == NULL) /* error */
    {
      free (state);
      mpz_set_si (f, -1);
      return NULL;
    }

  state->fd = (point *) malloc (params->size_fd * sizeof (point));
  if (state->fd == NULL)
    {
      clear_list (coeffs, params->size_fd);
      free (state);
      mpz_set_si (f, -1);
      return NULL;
    }
  for (k = 0; k < params->size_fd; k++)
    {
      MEMORY_TAG;
      mpres_init (state->fd[k].x, modulus);
      MEMORY_TAG;
      mpres_init (state->fd[k].y, modulus);
      MEMORY_UNTAG;
    }
  
  state->size_T = params->size_fd + 4;
  state->T = (mpres_t *) malloc (state->size_T * sizeof (mpres_t));
  if (state->T == NULL)
    {
      for (k = 0; k < params->size_fd; k++)
        {
          mpres_clear (state->fd[k].x, modulus);
          mpres_clear (state->fd[k].y, modulus);
        }
      clear_list (coeffs, params->size_fd);
      free (state);
      mpz_set_si (f, -1);
      return NULL;
    }
  for (k = 0; k < state->size_T; k++)
    {
      MEMORY_TAG;
      mpres_init (state->T[k], modulus);
      MEMORY_UNTAG;
    }

  for (k = params->S + 1; k < params->size_fd; k += params->S + 1)
     mpz_set_ui (coeffs[k + params->S], 1);

  if (test_verbose (OUTPUT_TRACE))
    for (k = 0; k < params->size_fd; k++)
      outputf (OUTPUT_TRACE, "ecm_rootsG_init: coeffs[%d] == %Zd\n", 
               k, coeffs[k]);

  youpi = multiplyW2n (f, state->fd, X, coeffs, params->size_fd, modulus, 
                     state->T[0], state->T[1], state->T + 2, &muls, &gcds);
  if (youpi == ECM_ERROR)
    mpz_set_si (f, -1); /* fall through */

  for (k = params->S + 1; k < params->size_fd; k += params->S + 1)
    {
      mpres_set (state->fd[k + params->S].x, state->fd[params->S].x, modulus);
      mpres_set (state->fd[k + params->S].y, state->fd[params->S].y, modulus);
    }
  
  clear_list (coeffs, params->size_fd);
  coeffs = NULL;
  
  if (youpi != ECM_NO_FACTOR_FOUND) /* factor found or error */
    {
      if (youpi == ECM_FACTOR_FOUND_STEP2)
        outputf (OUTPUT_VERBOSE, "Found factor while computing fd[]\n");

      ecm_rootsG_clear (state, modulus);
      
      /* Signal that a factor was found, or an error occurred (f=-1) */
      state = NULL;
    }
  else
    {
      if (test_verbose (OUTPUT_VERBOSE))
        {
          st = elltime (st, cputime ());
          outputf (OUTPUT_VERBOSE,
		   "Initializing table of differences for G took %ldms", st);
	  outputf (OUTPUT_DEVVERBOSE, ", %lu muls and %lu extgcds",
		   muls, gcds);
          outputf (OUTPUT_VERBOSE, "\n");
        }
    }
  
  return state;
}

void 
ecm_rootsG_clear (ecm_roots_state_t *state, ATTRIBUTE_UNUSED mpmod_t modulus)
{
  unsigned int k;
  
  for (k = 0; k < state->params.size_fd; k++)
    {
      mpres_clear (state->fd[k].x, modulus);
      mpres_clear (state->fd[k].y, modulus);
    }
  free (state->fd);
  
  for (k = 0; k < state->size_T; k++)
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
ecm_rootsG (mpz_t f, listz_t G, unsigned long dF, ecm_roots_state_t *state, 
            mpmod_t modulus)
{
  unsigned long i;
  unsigned long muls = 0, gcds = 0;
  int youpi = ECM_NO_FACTOR_FOUND;
  long st;
  point *fd = state->fd; /* to save typing */
  progression_params_t *params = &(state->params); /* for less typing */
  
  st = cputime ();

  outputf (OUTPUT_TRACE, "ecm_rootsG: dF = %lu, state: nr = %u, next = %u, "
           "S = %u, dsieve = %u, rsieve = %u,\n\tdickson_a = %d\n", 
           dF, params->nr, params->next, params->S, params->dsieve, 
	   params->rsieve, params->dickson_a);
  
  for (i = 0; i < dF;)
    {
      /* Did we use every progression since the last update? */
      if (params->next == params->nr)
        {
          /* Yes, time to update again */
	  youpi = addWnm (f, fd, state->X, modulus, params->nr, params->S, 
			  state->T, &muls, &gcds);
	  ASSERT(youpi != ECM_ERROR); /* no error can occur in addWnm */
	  params->next = 0;
          
          if (youpi == ECM_FACTOR_FOUND_STEP2)
            {
	      outputf (OUTPUT_VERBOSE, "Found factor while computing G[]\n");
              break;
            }
        }
      
      /* Is this a root we should skip? (Take only if gcd == 1) */
      if (gcd ((unsigned long) params->rsieve,  
               (unsigned long) params->dsieve) == 1UL)
	{
	  mpres_get_z (G[i++], (fd + params->next * (params->S + 1))->x, 
		       modulus);
	  outputf (OUTPUT_TRACE, 
		   "ecm_rootsG: storing d1*%u*X = %Zd in G[%lu]\n",
		   params->rsieve, G[i - 1], i);
	}
      
      params->next ++;
      params->rsieve ++;
    }
  
  outputf (OUTPUT_VERBOSE, "Computing roots of G took %ldms",
	   elltime (st, cputime ()));
  outputf (OUTPUT_DEVVERBOSE, ", %lu muls and %lu extgcds", muls, gcds);
  outputf (OUTPUT_VERBOSE, "\n");
  
  return youpi;
}


/* Find smallest i >= 0 such that 
   f(j * d2)*X = +-f((i0 + i) * d1)*X over GF(p).
   If "+" holds, return 1, if "-" holds, return -1.
   If the correct i could not be determined (because a non-invertible 
   residue appeared during initialisation) return 0. */
int 
ecm_findmatch (unsigned long *I, const unsigned long j, 
               root_params_t *root_params, const curve *X, mpmod_t n, 
               const mpz_t p)
{
  const int dickson_a = root_params->S < 0 ? -1 : 0;
  const unsigned int S = abs (root_params->S);
  const unsigned int sizeT = S + 3;
  unsigned int k;
  unsigned long i;
  int r, sgn = 0;
  point iX, jX;
  curve Xp; /* The point and curve over GF(p) */
  mpmod_t modulus;
  mpz_t s, t; /* temp vars */
  mpres_t u, v; /* temp vars */
  listz_t coeffs;
  point *fd;
  mpres_t *T;
  
  outputf (OUTPUT_RESVERBOSE, "Looking for i such that "
           "f((i+%Zd)*%lu)*X = f(%lu*%lu)*X\n", root_params->i0, 
           root_params->d1, j, root_params->d2);
  
  mpmod_init (modulus, p, ECM_MOD_DEFAULT);
  mpz_init (s);
  mpz_init (t);
  mpres_init (u, modulus);
  mpres_init (v, modulus);
  mpres_init (Xp.x, modulus);
  mpres_init (Xp.y, modulus);
  mpres_init (Xp.A, modulus);
  mpres_init (iX.x, modulus);
  mpres_init (iX.y, modulus);
  mpres_init (jX.x, modulus);
  mpres_init (jX.y, modulus);
  T = malloc (sizeT * sizeof (mpres_t));
  if (T == NULL)
    goto clear_and_exit;
  for (k = 0; k < sizeT; k++)
    mpres_init (T[k], modulus);
  fd = malloc ((S + 1) * sizeof (point));
  if (fd == NULL)
    goto clear_T_and_exit;
  for (k = 0; k < S + 1; k++)
    {
      mpres_init (fd[k].x, modulus);
      mpres_init (fd[k].y, modulus);
    }

  /* Copy the parameters of the curve over Z/ZN to the curve over GF(p) */
  mpres_get_z (t, X->x, n);
  mpres_set_z (Xp.x, t, modulus);
  mpres_get_z (t, X->y, n);
  mpres_set_z (Xp.y, t, modulus);
  mpres_get_z (t, X->A, n);
  mpres_set_z (Xp.A, t, modulus);
  
  /* We use init_progression_coeffs() to compute f(j * d2) */
  mpz_set_ui (t, j);
  coeffs = init_progression_coeffs (t, 1UL, root_params->d2, 1U, 1U, S, 
                                    dickson_a);
  if (coeffs == NULL)
    goto clear_fd_and_exit;
  
  /* Now compute f(j * d2) X */
  r = multiplyW2n (NULL, &jX, &Xp, coeffs, 1U, modulus, u, v, T, NULL, NULL);
  clear_list (coeffs, S + 1);
  if (r != ECM_NO_FACTOR_FOUND)
    goto clear_fd_and_exit;

  /* We'll keep {f(j * d2) X}_x in s */
  mpres_get_z (s, jX.x, modulus);
  outputf (OUTPUT_DEVVERBOSE, "ecm_findmatch: (f(j * d2) X)_x = %Zd\n", s);

  /* Now compute {f((i0 + i) d1) X}_x one at a time and put them in t, 
     until s == t */
  
  /* Init the progression */
  coeffs = init_progression_coeffs (root_params->i0, 1UL, root_params->d1, 
                                    1U, 1U, S, dickson_a);
  if (coeffs == NULL)
    goto clear_fd_and_exit;
  r = multiplyW2n (NULL, fd, &Xp, coeffs, S + 1, modulus, u, v, T, NULL, NULL);
  clear_list (coeffs, S + 1);
  if (r != ECM_NO_FACTOR_FOUND)
    goto clear_fd_and_exit;
  
  mpres_get_z (t, fd[0].x, modulus);
  for (i = 0; mpz_cmp (s, t) != 0; i++)
    {
      r = addWnm (NULL, fd, &Xp, modulus, 1, S, T, NULL, NULL);
      if (r != ECM_NO_FACTOR_FOUND)
        goto clear_fd_and_exit;
      mpres_get_z (t, fd[0].x, modulus);
    }

  outputf (OUTPUT_DEVVERBOSE, "ecm_findmatch: i - i0 = %lu, "
           "{f(i * d1) X}_x = %Zd\n", i, t);

  /* We'll compute f(i * d1)*X and compare it to f(j * d2)*X to verify 
     correctness of the result, and to determine whether it was
     f(i * d1)-f(j * d2) or f(i * d1)+f(j * d2) that found the factor */
  /* We use init_progression_coeffs() to compute f(i * d1) */
  mpz_add_ui (t, root_params->i0, i);
  coeffs = init_progression_coeffs (t, 1UL, root_params->d1, 1U, 1U, S, 
                                    dickson_a);
  if (coeffs == NULL)
    goto clear_fd_and_exit;
  
  /* Now compute iX = f(i * d1)*X */
  r = multiplyW2n (NULL, &iX, &Xp, coeffs, 1U, modulus, u, v, T, NULL, NULL);
  clear_list (coeffs, S + 1);
  if (r != ECM_NO_FACTOR_FOUND)
    goto clear_fd_and_exit;
  
  mpres_get_z (t, iX.x, modulus);
  if (mpz_cmp (s, t) != 0)
    {
      outputf (OUTPUT_ERROR, "ecm_findmatch: ERROR, (f(i*d1) X)_x != "
               "(f(j*d2) X)_x\n(f(i*d1) X)_x = %Zd\n", t);
      goto clear_fd_and_exit;
    }

  mpres_get_z (s, jX.y, modulus);
  mpres_get_z (t, iX.y, modulus);
  if (mpz_cmp (s, t) == 0)
    {
      *I = i;
      sgn = 1;
    }
  else
    {
      mpz_sub (t, p, t);
      if (mpz_cmp (s, t) == 0)
        {
          *I = i;
          sgn = -1;
        }
      else
        {
          mpz_sub (t, p, t);
          outputf (OUTPUT_ERROR, "ecm_findmatch: ERROR, (f(i*d1) X)_y != "
                   "+-(f(j*d2) X)_y\n");
          outputf (OUTPUT_ERROR, "(f(i*d1) X)_y = %Zd\n", t);
          outputf (OUTPUT_ERROR, "(f(j*d2) X)_y = %Zd\n", s);
        }
    }

clear_fd_and_exit:
  for (k = 0; k < S + 1; k++)
    {
      mpres_clear (fd[k].x, modulus);
      mpres_clear (fd[k].y, modulus);
    }
  free(fd);
clear_T_and_exit:
  for (k = 0; k < sizeT; k++)
    mpres_clear (T[k], modulus);
  free (T);
clear_and_exit:
  mpz_clear (s);
  mpz_clear (t);
  mpres_clear (u, modulus);
  mpres_clear (v, modulus);
  mpres_clear (Xp.x, modulus);
  mpres_clear (Xp.y, modulus);
  mpres_clear (Xp.A, modulus);
  mpres_clear (iX.x, modulus);
  mpres_clear (iX.y, modulus);
  mpres_clear (jX.x, modulus);
  mpres_clear (jX.y, modulus);
  mpmod_clear (modulus);
  
  return sgn;
}
