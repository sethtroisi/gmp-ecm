/* Common stage 2 for ECM, P-1 and P+1 (improved standard continuation).

  Copyright (C) 2001 Paul Zimmermann,
  LORIA/INRIA Lorraine, zimmerma@loria.fr
  See http://www.loria.fr/~zimmerma/records/ecmnet.html

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

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "gmp.h"
#include "ecm.h"
#include "cputime.h"

void dickson_ui(mpz_t r, unsigned int x, unsigned int n, int a);
int fin_diff_init (point **, curve, unsigned int, unsigned int, unsigned int,
                   mpz_t, int, int, mpz_t);
int   fin_diff_next  (mpz_t, point *, unsigned int, mpz_t, int);
void  fin_diff_clear (point *, unsigned int, int);
int          rootsF     (mpz_t, listz_t, unsigned int, curve, listz_t, 
                         unsigned int , mpz_t, int, int, int);
void         rootsG     (mpz_t, listz_t, unsigned int, point *, point *,
                         listz_t, unsigned int, mpz_t, int, int);

#define INVF /* precompute 1/F for divisions by F */

void 
dickson_ui(mpz_t r, unsigned int x, unsigned int n, int a)
{
  unsigned int i, b = 0;
  mpz_t t, u;

  if (n == 0)
    {
      mpz_set_ui (r, 2);
      return;
    }
  
  while (n > 2 && (n & 1) == 0)
    {
      b++;
      n >>= 1;
    }
  
  mpz_set_ui (r, x);
  
  mpz_init(t);
  mpz_init(u);

  if (n > 1)
    {
      mpz_set_ui (r, x);
      mpz_mul_ui (r, r, x);
      mpz_sub_si (r, r, a);
      mpz_sub_si (r, r, a); /* r = dickson(x, 2, a) */
      
      mpz_set_ui (t, x);    /* t = dickson(x, 1, a) */
      
      for (i = 2; i < n; i++)
        {
          mpz_mul_si (u, t, a);
          mpz_set (t, r);     /* t = dickson(x, i, a) */
          mpz_mul_ui (r, r, x);
          mpz_sub (r, r, u);  /* r = dickson(x, i+1, a) */
        }
    }
  
  for ( ; b > 0; b--)
    {
      mpz_mul (t, r, r); /* t = dickson(x, n, a) ^ 2 */
      mpz_ui_pow_ui (u, abs(a), n);
      if (n & 1 && a < 0)
        mpz_neg (u, u);
      mpz_mul_2exp (u, u, 1); /* u = 2 * a^n */
      mpz_sub (r, t, u); /* r = dickson(x, 2*n, a) */
      n <<= 1;
    }
  
  mpz_clear(t);
  mpz_clear(u);
  
}

/* Init table to allow successive computation of:
   X^((s + n*D)^E) mod N for P-1,
   V_{(s + n*D)^E}(X) for P+1,
   ((s + n*D)^E)*X for ECM.
   See Knuth, TAOCP vol.2, 4.6.4 and exercise 7 in 4.6.4.
   For s=7, D=6, E=1: fd[0] = X^7, fd[1] = X^6.
   For s=0, D=d, E=1: fd[0] = X^0, fd[1] = X^d.
   Return non-zero iff a factor was found (then stored in f).
*/
int
fin_diff_init (point **Fd,
               curve X, unsigned int s, unsigned int D, unsigned int E,
               mpz_t N, int method, int use_dickson, mpz_t f)
{
  unsigned int i, k, allocated;
  mpz_t P, Q;      /* for P+1 */
  mpz_t p, x1, u, v; /* for ECM */
  point *fd;
  int youpi = 0;

  allocated = E + 1;

  if (method == PP1_METHOD)
    {
      mpz_init (P);
      mpz_init (Q);
      allocated += 2;
    }
  else if (method == EC_METHOD)
    {
      mpz_init (p);
      mpz_init (x1);
      mpz_init (u);
      mpz_init (v);
      allocated += E + 1; /* auxiliary space for addWn */
    }

  fd = (point *) malloc (allocated * sizeof(point));
  for (i = 0; i < allocated; i++)
    mpz_init (fd[i].x);
  if (method == EC_METHOD)
    for (i = 0; i < allocated; i++)
      mpz_init (fd[i].y);
  
  for (i = 0; i <= E; i++)
    if (use_dickson)
      dickson_ui(fd[i].x, s + i * D, E, -1);
    else
      mpz_ui_pow_ui (fd[i].x, s + i * D, E); /* fd[i] = (s+i*D)^E */
  
  for (k = 1; k <= E; k++)
    for (i = E; i >= k; i--)
      mpz_sub (fd[i].x, fd[i].x, fd[i-1].x);
  
  for (i = 0; i <= E && youpi == 0; i++)
    if (method == PM1_METHOD)
      mpz_powm (fd[i].x, X.x, fd[i].x, N);
    else if (method == PP1_METHOD)
      pp1_mul (fd[i].x, X.x, fd[i].x, N, P, Q);
    else /* ECM */
      {
        /* copy in x1 since multiplyW2 doesn't allow x1 and q to be equal */
        youpi = multiplyW2 (f, x1, fd[i].y, X.x, X.y, fd[i].x, N, X.A, u, v);
        mpz_set (fd[i].x, x1);
      }

  if (method == PP1_METHOD) /* necessarily E=1 */
    {
      /* fd[0] = V_s(x), fd[1] = V_D(x) */
      mpz_set_ui (fd[2].x, (s > D) ? (s - D) : (D - s));
      pp1_mul (fd[2].x, X.x, fd[2].x, N, P, Q); /* V_{s-D}(x) */
      mpz_clear (P);
      mpz_clear (Q);
    }

  if (method == EC_METHOD)
    {
      mpz_clear (p);
      mpz_clear (x1);
      mpz_clear (u);
      mpz_clear (v);
    }

  *Fd = fd;
  return youpi;
}

/* P-1: Computes x^((s + (n+1)*D)^E and stores in fd[0]
   P+1: Computes V_{j+D} from fd[0] = V_j, fd[1] = V_D, and fd[2] = V_{j-D}.
   Return non-zero iff a factor was found (can happen only with ECM), in which
   case the factor is stored in f.
*/
int
fin_diff_next (mpz_t f, point *fd, unsigned int E, mpz_t N, int method)
{
  unsigned int i;

  if (method == PP1_METHOD)
    {
      mpz_swap (fd[0].x, fd[2].x);
      mpz_mul (fd[3].x, fd[2].x, fd[1].x);
      mpz_sub (fd[0].x, fd[3].x, fd[0].x);
      mpz_mod (fd[0].x, fd[0].x, N);
      return 0;
    }

  if (method == EC_METHOD)
    return addWn (f, fd, N, E);
  else /* P-1 */
    for (i = 0; i < E; i++)
      {
        mpz_mul (fd[i].x, fd[i].x, fd[i+1].x);
        mpz_mod (fd[i].x, fd[i].x, N);
      }

  return 0;
}

void 
fin_diff_clear (point *fd, unsigned int E, int method)
{
  unsigned int i, allocated = E + 1;

  if (method == PP1_METHOD)
    allocated += 2;
  else if (method == EC_METHOD)
    allocated += E + 1;
  
  for (i = 0; i < allocated; i++)
    mpz_clear (fd[i].x);
  if (method == EC_METHOD)
    for (i = 0; i < allocated; i++)
      mpz_clear (fd[i].y);

  free (fd);
  return;
}

/* puts in F[0..dF-1] the successive values of 

   s^(j^S) for Pollard P-1
   V_j(P) for Williams P+1 [here S=1]
   (j^S) * P for ECM where P=(s : : 1) is a point on the elliptic curve

   For P+1, we have V_{j+6} = V_j * V_6 - V_{j-6}.

   for 0 < j = 1 mod 7 < d, j and d coprime.
   Returns non-zero iff a factor was found (then stored in f).
   P-1 only: If s.y=0, don't use the x+1/x trick, otherwise s.y=1/(s.x) mod n.
   Requires (dF+1) cells in t.
*/
int
rootsF (mpz_t f, listz_t F, unsigned int d, curve s, listz_t t,
        unsigned int S, mpz_t n, int verbose, int method, int use_dickson)
{
  unsigned int i, j;
  int st, st2;
  point *fd;
  int youpi = 0;
  
  st = cputime ();

  mpz_set (F[0], s.x); /* s^1 for P-1, V_1(P)=P for P+1, (1*P)=P for ECM */
  i = 1;
  if (d > 7)
    {
      st2 = cputime ();
      youpi = fin_diff_init (&fd, s, 7, 6, S, n, method, use_dickson, f);
      /* for P+1, fd[0] = V_7(P), fd[1] = V_6(P) */
      if (verbose >= 2)
        printf ("Initializing table of differences for F took %dms\n", cputime () - st2);
      j = 7;
      while (j < d && youpi == 0)
        {
          if (gcd (j, d) == 1)
            {
              mpz_set (F[i], fd[0].x);
              i++;
            }
          youpi = fin_diff_next (f, fd, S, n, method);
          j += 6;
        }
      fin_diff_clear (fd, S, method);
    }

  if (youpi)
    return 1;

  if (method == PM1_METHOD && mpz_cmp_ui (s.y, 0) != 0)
    {
      if ((d/6)*S > 3*(i-1)) /* Batch inversion is cheaper */
        {
          if (list_invert (t, F, i, t[i], n)) 
            {
              mpz_set (f, t[i]);
              return 1;
            }
     
          for (j = 0; j < i; j++) 
            {
              mpz_add (F[j], F[j], t[j]);
              mpz_mod (F[j], F[j], n);
            }
     
          mpz_set (s.y, t[0]); /* Save s^(-1) in s.y */
    
        }
      else
        { /* fin_diff code is cheaper */

#if 0 /* s.y already contains 1/(s.x) */
          mpz_gcdext (f, s.y, NULL, s.x, n);

          if (mpz_cmp_ui (f, 1) != 0)
            return 1;
#endif

          mpz_add (F[0], F[0], s.y);
          mpz_mod (F[0], F[0], n);

          i = 1;
       
          if (d > 7) 
            {
              curve Y; /* for 1/s */
              mpz_init (Y.x);
              mpz_set (Y.x, s.y);
              youpi = fin_diff_init (&fd, Y, 7, 6, S, n, method, use_dickson, f);
              mpz_clear (Y.x);
              j = 7;
              while (j < d && youpi == 0)
                {
                  if (gcd (j, d) == 1)
                    {
                      mpz_add (F[i], F[i], fd[0].x);
                      mpz_mod (F[i], F[i], n);
                      i++;
                    }
                  youpi = fin_diff_next (f, fd, S, n, method);
                  j += 6;
                }
              fin_diff_clear (fd, S, method);
            }
        }
    }
  
  if (verbose >= 2)
    printf ("Computing roots of F took %dms\n", cputime () - st);

  return youpi;
}


/* puts in G the successive values of t0*s^j, for 1 <= j <= d
   returns in t the value of t0*s^d, in u the value of u0*invs^d.
   If listz_t=NULL, don't use the x+1/x trick.
   Needs d+1 cells in t.
*/
void
rootsG (mpz_t f, listz_t G, unsigned int d, point *fd, point *fdinv, 
        listz_t t, unsigned int S, mpz_t n, int verbose, int method)
{
  unsigned int i;
  int st;

  st = cputime ();

  if (fdinv != NULL && S > 3) /* Montgomery's trick to perform only
                                 one inversion */
    {
      mpz_set (G[0], fd[0].x);
      mpz_set (t[0], fd[0].x);

      for (i = 1; i < d; i++)
        {
          fin_diff_next (f, fd, S, n, method);
          mpz_set (G[i], fd[0].x);
          mpz_mul (t[i], t[i-1], fd[0].x);
          mpz_mod (t[i], t[i], n);
        }

      mpz_gcdext (t[d], t[d-1], NULL, t[d-1], n);

      for (i = d-1; i > 0; i--) 
        {
          mpz_mul (t[d], t[i], t[i-1]);
          mpz_mul (t[i-1], t[i], G[i]);
          mpz_mod (t[i-1], t[i-1], n);
          mpz_add (t[d], t[d], G[i]);
          mpz_mod (G[i], t[d], n);
        }

    }
  else /* fdinv=NULL or S <= 2 */
    {
      for (i = 0; i < d; i++)
        {
          if (fdinv != NULL)
            {
              mpz_add (G[i], fd[0].x, fdinv[0].x);
              mpz_mod (G[i], G[i], n);
              fin_diff_next (f, fdinv, S, n, method);
            }
          else
            mpz_set (G[i], fd[0].x);
          fin_diff_next (f, fd, S, n, method);
        }
    }

  if (verbose >= 2)
    printf ("Computing roots of G took %dms\n", cputime () - st);
}

/* Input:  X is the point at end of stage 1
           n is the number to factor
           B2 is the stage 2 bound
           k is the number of blocks
           S is the exponent for Brent-Suyama's extension
           verbose is the verbose level
           invtrick is non-zero iff one uses x+1/x instead of x.
           method: EC_METHOD, PM1_METHOD or PP1_METHOD
           Cf "Speeding the Pollard and Elliptic Curve Methods
               of Factorization", Peter Montgomery, Math. of Comp., 1987,
               page 257: using x^(i^e)+1/x^(i^e) instead of x^(i^(2e))
               reduces the cost of Brent-Suyama's extension from 2*e
               to e+3 multiplications per value of i.
   Output: f is the factor found
   Return value: non-zero iff a factor was found.
*/
int
stage2 (mpz_t f, curve *X, mpz_t n, double B2, unsigned int k, unsigned int S, 
        int verbose, int method, double B1)
{
  int invtrick = method == PM1_METHOD;
  /*  int use_dickson = !invtrick; */
  int use_dickson = 0;
  double b2;
  unsigned int i, d, dF, sizeT;
  unsigned long muls;
  listz_t F, G, H, T;
  point *fd_x, *fd_invx = NULL;
  polyz_t polyF, polyT;
  mpz_t invx;
  int youpi = 0, st, st0;
#ifdef INVF
  listz_t invF = NULL;
#endif

  if (B2 <= B1)
    return 0;

  st0 = cputime ();

  if (verbose >= 2)
    {
      printf ("starting stage 2 with x=");
      mpz_out_str (stdout, 10, X->x);
      putchar ('\n');
    }

  b2 = ceil(B2 / k); /* b2 = ceil(B2/k): small block size */

  d = bestD (b2);

#if 0
  if (2.0 * (double) d > B1)
    {
      fprintf (stderr, "Error: 2*d > B1\n");
      exit (1);
    }
#endif

  b2 = block_size (d);

  B2 = (double) k * b2;

  dF = phi (d) / 2;

  if (verbose >= 2)
    printf ("B2=%1.0f k=%u b2=%1.0f d=%u dF=%u\n", B2, k, b2, d, dF);

  F = init_list (dF + 1);

  /* if method <> PM1_METHOD, invtrick=0 thus invx is 0 */
  if (invtrick)
    {
      mpz_init (invx);
      mpz_gcdext (f, invx, NULL, X->x, n);
      mpz_set (X->y, invx);
    }

  sizeT = 3 * dF - 1 + list_mul_mem (dF);
#ifdef INVF
  if (dF > 3)
    sizeT += dF - 3;
#endif
  T = init_list (sizeT);
  H = T;

  /* needs dF+1 cells in T */
  if (rootsF (f, F, d, *X, T, S, n, verbose, method, use_dickson))
    {
      mpz_set (f, X->x);
      youpi = 2;
      goto clear_F;
    }

  /* ----------------------------------------------
     |   F    |  invF  |   G   |         T        |
     ----------------------------------------------
     | rootsF |  ???   |  ???  |      ???         |
     ---------------------------------------------- */

  PolyFromRoots (F, dF, T, verbose, n, 'F'); /* needs dF+list_mul_mem(dF/2) cells in T */
  mpz_set_ui (F[dF], 1); /* the leading monic coefficient needs to be stored
                             explicitely for PrerevertDivision and polygcd */

  /* ----------------------------------------------
     |   F    |  invF  |   G   |         T        |
     ----------------------------------------------
     |  F(x)  |  ???   |  ???  |      ???         |
     ---------------------------------------------- */

#ifdef INVF
  /* G*H has degree 2*dF-2, hence we must cancel dF-1 coefficients
     to get degree dF-1 */
  if (dF > 1)
    {
      invF = init_list (dF - 1);
      st = cputime ();
#if 0
      list_zero (T, 2 * dF - 3);
      mpz_set_ui (T[2 * dF - 3], 1); /* T[0..2dF-3] = x^(2dF-3) */
      muls = RecursiveDivision (invF, T, F + 1, dF - 1, T + 2 * dF - 2, n);
#else
      muls = PolyInvert (invF, F + 2, dF - 1, T, n);
#endif
      /* now invF[0..K-2] = Quo(x^(2dF-3), F) */
      if (verbose >= 2)
        printf ("Computing 1/F took %ums and %lumuls\n", cputime() - st, muls);
      
      /* ----------------------------------------------
         |   F    |  invF  |   G   |         T        |
         ----------------------------------------------
         |  F(x)  | 1/F(x) |  ???  |      ???         |
         ---------------------------------------------- */
    }
#endif

  G = init_list (dF);
  st = cputime ();
  if ((youpi = fin_diff_init (&fd_x, *X, 2*d, d, S, n, method, use_dickson, f)))
    goto clear_fd;

  if (verbose >= 2)
    printf ("Initializing table of differences for G took %dms\n", cputime () - st);

  if (invtrick) /* P-1: only x is needed */
    {
      curve Y;
      mpz_init_set (Y.x, invx);
      youpi = fin_diff_init (&fd_invx, Y, 2*d, d, S, n, method, use_dickson, f);
      mpz_clear (Y.x);
      if (youpi)
        goto clear_fd_invx;
    }

  for (i=0; i<k; i++)
    {
      /* needs dF+1 cells in T+dF */
      rootsG (f, G, dF, fd_x, fd_invx, T + dF, S, n, verbose, method);

  /* -----------------------------------------------
     |   F    |  invF  |   G    |         T        |
     -----------------------------------------------
     |  F(x)  | 1/F(x) | rootsG |      ???         |
     ----------------------------------------------- */

      PolyFromRoots (G, dF, T + dF, verbose, n, 'G'); /* needs 2*dF+list_mul_mem(dF/2) cells in T */

  /* -----------------------------------------------
     |   F    |  invF  |   G    |         T        |
     -----------------------------------------------
     |  F(x)  | 1/F(x) |  G(x)  |      ???         |
     ----------------------------------------------- */

      if (i == 0)
        {
          list_sub (H, G, F, dF); /* coefficients 1 of degree cancel,
                                     thus T is of degree < dF */
          /* ------------------------------------------------
             |   F    |  invF  |    G    |         T        |
             ------------------------------------------------
             |  F(x)  | 1/F(x) |  ???    |G(x)-F(x)|  ???   |
             ------------------------------------------------ */
        }
      else
	{
          /* since F and G are monic of same degree, G mod F = G - F */
          list_sub (G, G, F, dF);

          /* ------------------------------------------------
             |   F    |  invF  |    G    |         T        |
             ------------------------------------------------
             |  F(x)  | 1/F(x) |G(x)-F(x)|  H(x)  |         |
             ------------------------------------------------ */

	  st = cputime ();
	  /* previous G mod F is in H, with degree < dF, i.e. dF coefficients:
	     requires 3dF-1+list_mul_mem(dF) cells in T */
	  muls = list_mulmod2 (H, T + dF, G, H, dF, T + 3 * dF - 1, n);
          if (verbose >= 2)
            printf ("Computing G * H took %ums and %lumuls\n", cputime() - st,
                    muls);

          /* ------------------------------------------------
             |   F    |  invF  |    G    |         T        |
             ------------------------------------------------
             |  F(x)  | 1/F(x) |G(x)-F(x)| G * H  |         |
             ------------------------------------------------ */

	  st = cputime ();
#ifdef INVF
          muls = PrerevertDivision (H, F, invF, dF, T + 2 * dF - 1, n);
#else
          mpz_set_ui (T[2*dF-1], 0); /* since RecursiveDivision expects a
                                        dividend of 2*dF coefficients */
	  muls = RecursiveDivision (G, H, F, dF, T + 2 * dF, n);
#endif
          if (verbose >= 2)
            printf ("Reducing G * H mod F took %ums and %lumuls\n",
                    cputime() - st, muls);
	}
    }

  st = cputime ();
  init_poly_list (polyF, dF, F);
  init_poly_list (polyT, dF - 1, T);
  if ((youpi = poly_gcd (f, polyF, polyT, n, T + dF)))
    NTL_get_factor (f);
  if (verbose >= 2)
    printf ("Computing gcd of F and G took %dms\n", cputime() - st);

 clear_fd_invx:
  if (invtrick)
    fin_diff_clear (fd_invx, S, method);
 clear_fd:
  fin_diff_clear (fd_x, S, method);
  clear_list (G, dF);

#ifdef INVF
  if (dF > 1)
    clear_list (invF, dF - 1);
#endif

 clear_F:
  clear_list (T, sizeT);
  if (invtrick)
    mpz_clear (invx);
  clear_list (F, dF + 1);

  if (verbose >= 1)
    printf ("Stage 2 took %dms\n", cputime() - st0);

  return youpi;
}
