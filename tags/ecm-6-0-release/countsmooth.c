/* Enumerate smooth values from a range of numbers 
   (with Brent-Suyama extension)

  Copyright 2003 Alexander Kruppa.

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
#include <string.h>
#include <math.h>
#include <gmp.h>

/* uncomment the following to use primegen,
   cf http://cr.yp.to/primegen/primegen-0.97.tar.gz */
/* #define PRIMEGEN */

#ifdef PRIMEGEN
#include <primegen.h>
#else
#include "ecm.h" /* for getprime() */
#endif

#define mulmod(r,u,v,n) mpz_mul(r,u,v);mpz_mod(r,r,n);

#define DEBUG

void          dicksonmod(mpz_t, mpz_t, unsigned int, mpz_t, mpz_t);
int           is_P_minus_i (mpz_t, unsigned int);
unsigned int  eulerphi (unsigned int);
unsigned int  get_lenF (unsigned int);
unsigned int  gcd (unsigned int, unsigned int);
void          quicksort (mpz_t *, unsigned int);
void          quicksort_with_index (mpz_t *, unsigned int *, unsigned int);
int           issorted (mpz_t *, unsigned int);
int           getparm_ui (int, char **, int, char *, unsigned int *);
int           getparm_d (int, char **, int, char *, double *);
int           getparm_mpz (int, char **, int, char *, mpz_t *);
int           brent_suyama_match (mpz_t, unsigned int, unsigned int, 
                 unsigned int, unsigned int, mpz_t, mpz_t *, mpz_t *,
                 unsigned int *, unsigned int *); 
double        brent_suyama_theo (mpz_t, unsigned int, unsigned int, 
                  unsigned int, unsigned int, mpz_t);
void          print_help (void);

/* Computes Dickson_{n, a}(x), the degree n Dickson polynomial with
   parameter a, evaluated at point x, and returns value modulo p in r */
void dicksonmod(mpz_t r, mpz_t x, unsigned int n, mpz_t a, 
	mpz_t p)
{
  unsigned int i, b = 0;
  mpz_t t, u;

  mpz_init(t);
  mpz_init(u);

  if (n == 0)
    {
      mpz_set_ui (r, 1);
      return;
    }
  
  if (n == 1)
    {
      mpz_set (r, x);
      return;
    }

  /* Now n >= 2 */
  while (n > 2 && (n & 1) == 0)
   {
     b++;
     n >>= 1;
   }
  
  mpz_mul (t, x, x);   /* r = x^2 */
  mpz_sub (t, t, a);
  mpz_sub (t, t, a);   /* r = x^2 - 2*a = Dickson_{2,a}(x) */
  mpz_mod (r, t, p);
  mpz_set (t, x);      /* t = Dickson_{1,a}(x) */
  
  for (i = 2; i < n; i++)
    {
      mulmod (u, t, a, p); /* u = a * Dickson_{i-1,a}(x) */
      mpz_set(t, r);       /* t = Dickson_{i,a}(x) */
      mpz_mul (r, r, x);   /* r = x * Dickson_{i,a}(x) */
      mpz_sub (r, r, u);   /* r = x * Dickson_{i,a}(x) - a * Dickson_{i-1,a}(x) */
      mpz_mod (r, r, p);   /*   = Dickson_{i+1,a}(x) */
    }
  
  for ( ; b > 0; b--)
    {
      mulmod (t, r, r, p); /* t = (Dickson_{n,a}(x))^2 */
      mpz_powm_ui (u, a, n, p);
      mpz_mul_2exp (u, u, 1);
      mpz_sub (r, t, u);   /* r = (Dickson_{n,a}(x))^2 - 2 * a^n */
      mpz_mod (r, r, p);   /*   = Dickson_{2*n,a}(x) */
    }

  mpz_clear (u);
  mpz_clear (t);
}


/* Test if N+1 is a probable prime */

int 
is_P_minus_i (mpz_t N, unsigned int i)
{
  int r;
  
  mpz_add_ui (N, N, i);
  r = mpz_probab_prime_p (N, 2);
  mpz_sub_ui (N, N, i);
  
  return r;
}


/* Euler Phi(n) function, number of residues coprime to n */
unsigned int 
eulerphi (unsigned int n) 
{
  unsigned int i, r=1;
  
  for (i=2; i*i<=n; i++)
    if (n%i == 0) 
      {
        r *= i-1;
        n /= i;
        while (n%i == 0) 
          {
            r *= i;
            n /= i;
          }
      }

  if (n>1) r *= n-1;
  
  return r;
}

unsigned int 
get_lenF (unsigned int D)
{
  return (eulerphi (D) / 2);
}

unsigned int 
gcd (unsigned int a, unsigned int b) 
{
  unsigned int t;
  
  while (b != 0) {
    t = a % b;
    a = b;
    b = t;
  }
  
  return a;
}

#define swap(a,b,t) {(t)=(a);(a)=(b);(b)=(t);}

/* Sorts list of mpz_t values in a of length len. pivot is a temp var */

void 
quicksort_with_index (mpz_t *data, unsigned int *index, unsigned int len) 
{
  unsigned int i, j, t;
  
  if (len <= 1) return;
  
  if (len == 2)
    {
      if (mpz_cmp (data[0], data[1]) > 0)
        {
          mpz_swap (data[0], data[1]);
          swap(index[0],index[1],t);
        }
      return;
    }

  i = 0; j = len;
  while (1) 
    {          /* Top half gets everything greater than pivot */
      while (mpz_cmp(data[++i], data[0]) <= 0 && i < j);
      while (mpz_cmp(data[--j], data[0]) > 0 && i < j);
      if (i >= j) break;
      
      mpz_swap (data[i], data[j]);
      swap(index[i],index[j],t);
    }

  mpz_swap (data[0], data[i-1]); 
  swap(index[0],index[i-1],t);

  quicksort_with_index (data, index, i-1);
  quicksort_with_index (data+i, index+i, len-i);
}

#ifdef NO_INDEX

/* Sorts list of mpz_t values in a of length len. pivot is a temp var */

void 
quicksort (mpz_t *data, unsigned int len) 
{
  unsigned int i, j;
  
  if (len <= 1) return;
  
  if (len == 2)
    {
      if (mpz_cmp (data[0], data[1]) > 0)
        mpz_swap (data[0], data[1]);
      return;
    }

  i = 0; j = len;
  while (1) 
    {          /* Top half gets everything greater than pivot */
      while (mpz_cmp(data[++i], data[0]) <= 0 && i < j);
      while (mpz_cmp(data[--j], data[0]) > 0 && i < j);
      if (i >= j) break;
      
      mpz_swap (data[i], data[j]);
    }

  mpz_swap (data[0], data[i-1]); 

  quicksort (data, i-1);
  quicksort (data+i, len-i);
}

#endif

int 
issorted (mpz_t *a, unsigned int len)
{
  unsigned int i;
  for (i = 0; i < len - 1; i++)
    {
      if (mpz_cmp (a[i], a[i+1]) > 0)
        return 0;
    }
  
  return 1;
}

int 
getparm_ui (int argc, char **argv, int i, char *parm, unsigned int *res) 
{
  if (strcmp(argv[i], parm) == 0) {
    i++;
    if (i >= argc) {
      printf("%s needs parameter\n", parm);
      exit(EXIT_FAILURE);
    }
    *res = atoi (argv[i]);
    return 1;
  }
  
  return 0;
}

int 
getparm_d (int argc, char **argv, int i, char *parm, double *res) 
{
  if (strcmp(argv[i], parm) == 0) {
    i++;
    if (i >= argc) {
      printf("%s needs parameter\n", parm);
      exit(EXIT_FAILURE);
    }
    *res = atof (argv[i]);
    return 1;
  }
  
  return 0;
}

int 
getparm_mpz (int argc, char **argv, int i, char *parm, mpz_t *res) 
{
  if (strcmp(argv[i], parm) == 0) {
    i++;
    if (i >= argc) {
      printf("%s needs parameter\n", parm);
      exit(EXIT_FAILURE);
    }
    
    if (index (argv[i], 'e'))
      {                      /* Hmm, looks like scientific notation */
        double t;
        t = atof (argv[i]);
        mpz_set_d (*res, t);
        mpz_add_ui (*res, *res, 1); /* Avoid getting a value with lots of 2
                                       factors because of mantissa truncation */
      } else
        mpz_set_str (*res, argv[i], 10);

    return 1;
  }
  
  return 0;
}


/* Searches for matches in the lists 
   F={f_S(i) % modulus : 0<i<D, i == 1 (mod 6), gcd(i,D)=1} 
   G={f_S(j*D) % modulus : startG<=j<=endG}. f_S(x) is the S-th power or the 
   degree S Dickson polynomial with parameter a.
   If case of a match, the i and j of the first match are written to luckyF 
   and luckyG, if they are not NULL.
*/

int 
brent_suyama_match (mpz_t modulus, unsigned int S, unsigned int D, 
        unsigned int startG, unsigned int endG, mpz_t a, mpz_t *F_param, 
        mpz_t *G_param, unsigned int *luckyF, unsigned int *luckyG)
{
  unsigned int lenF, lenG, i, j, nrfactors = 0;
  mpz_t *F, *G;
  mpz_t t;
  unsigned int *indexF, *indexG;
  
  if (endG < startG || D <= 1 || D % 6 != 0 || S == 0 || S & 1)
    return 0;
  
  lenF = get_lenF (D);
  lenG = endG - startG + 1;
  mpz_init (t);
  
  /* Allocate memory for F, if not passed by caller */

  if (F_param == NULL)
    {
      F = (mpz_t *) malloc (lenF * sizeof (mpz_t));
      for (i = 0; i < lenF; i++)
        mpz_init (F[i]);
    } else
      F = F_param;
  

  /* Allocate memory for G, if not passed by caller */

  if (G_param == NULL)
    {
      G = (mpz_t *) malloc (lenG * sizeof (mpz_t));
      for (i = 0; i < lenG; i++)
        mpz_init (G[i]);
    } else
      G = G_param;
  

  /* Build lists F and G */
  indexF = (unsigned int *) malloc (lenF * sizeof (unsigned int));
  for (i = 1, j = 0; i < D; i += 6)
    if (gcd (i, D) == 1)
      {
        mpz_set_ui (t, i);
        indexF[j] = i;
        if (mpz_sgn (a) == 0) /* If a == 0, do S-th power */
          mpz_powm_ui (F[j], t, S, modulus);
        else
          dicksonmod (F[j], t, S, a, modulus);
        j++;
      }

  mpz_set_ui (t, startG);
  mpz_mul_ui (t, t, D);
  indexG = (unsigned int *) malloc (lenG * sizeof (unsigned int));
  for (i = 0; i < lenG; i++) 
    {
      indexG[i] = startG + i;
      if (mpz_sgn (a) == 0)
        mpz_powm_ui (G[i], t, S, modulus);
      else
        dicksonmod (G[i], t, S, a, modulus);
      mpz_add_ui (t, t, D);
    }

  /* Sort lists */
  quicksort_with_index (F, indexF, lenF);
  quicksort_with_index (G, indexG, lenG);

#ifdef DEBUG
  /* Check sorted lists */
  if (!issorted (F, lenF) || !issorted (G, lenG))
    {
      fprintf (stderr, "Error, quicksort returns unsorted list\n");
      exit (EXIT_FAILURE);
    }
#endif
  
  /* Look for matches */
  for (i = 0, j = 0; i < lenF && j < lenG; i++) 
    {
      while (j < lenG && mpz_cmp (F[i], G[j]) > 0) 
        j++;
      if (j < lenG && mpz_cmp (F[i], G[j]) == 0 && 
          (startG > 0 || i+j > 0)) /* Eliminate spurious factors in Dickson */
        {
          unsigned int tj = j;
          if (nrfactors == 0 && luckyF != NULL  && luckyG != NULL)
            {
              *luckyF = indexF[i];
              *luckyG = indexG[j];
            }
          while (tj < lenG && mpz_cmp (F[i], G[tj]) == 0) 
            {
              nrfactors++;
              tj++;
            }
        }
  }
  
  if (F_param == NULL)
    {
      for (i = 0; i < lenF; i++)
        mpz_clear (F[i]);
      free (F);
    }

  free (indexF);
  
  if (G_param == NULL)
    {
      for (i = 0; i < lenG; i++)
        mpz_clear (G[i]);
      free (G);
    }
        
  free (indexG);
  
  mpz_clear (t);

  return nrfactors;
}

/* Compute probability that value "modulus" is found by Brent-Suyama extension */
double 
brent_suyama_theo (mpz_t modulus, unsigned int S, unsigned int D, 
        unsigned int startG, unsigned int endG, mpz_t a)
{
  unsigned int mod_S, gcd_S;
  double nr_points, p;
  
  if (endG < startG || D <= 1 || S == 0 || S & 1)
    return 0.;
  
  /* See if the linear factors catch it */
  if (mpz_cmp_d (modulus, (double)D * (double)(startG)) > 0 &&
      mpz_cmp_d (modulus, (double)D * (double)(endG)) < 0)
    return 1.0;
  
  mod_S = mpz_fdiv_ui (modulus, S);
  if (mpz_sgn (a) == 0)
    gcd_S = gcd (mod_S - 1, S) - 2;
  else
    gcd_S = (gcd (mod_S - 1, S) + gcd (mod_S + 1, S)) / 2 - 2;
  
  nr_points = (double) get_lenF (D) * (double) (endG - startG + 1);
  p = mpz_get_d (modulus);
  
  return 1.0 - pow(1.0 - (double)gcd_S / p, nr_points);
}

void 
print_help ()
{
  printf ("countsmooth [-N <N>] [-tests <nrtests>] [-blocks <blocks>] [-B1 <B1>]\n"
          " [-B2 <B2>] [-S <S>] [-D <D>] [-startG <startG>] [-endG <endG>] [-a <a>]\n"
          " [-maxBS <maxBS>] [-v] [-theo] [-pm1] [-ecm]\n\n\n");

  printf (" Determines which values in [N, N+nrtests-1] are smooth according to given\n"
          " parameters.\n");
  printf (" A value is smooth if its factorization contains only primes and prime\n"
          " powers <= B1 (B1-smooth), or if one prime factor is >B1 but <=B2 and all\n"
          " others are B1-smooth, or if the largest prime factor can be expressed as\n"
          " f_{S,a}(k*D) - f_{S,a}(i), where startG<=k<=endG and 0<i<D, gcd(i,D)==1,\n"
          " i==1 (mod 6). f_{S,a}(x) is the S-th Dickson polynomial with parameter a,\n"
          " or the S-th power for a==0.\n\n");
          

  printf (" <N>         Start of range of numbers to invesigate.\n");
  printf (" <nrtests>   Length of range of numbers to investigate.\n");
  printf (" <blocks>    Split range into this many blocks to conserve memory.\n"
          "             By default <blocks> = 1\n");
  printf (" <B1>        B1-smoothness limit.\n");
  printf (" <B2>        B2-smoothness limit.\n");
  printf (" <S>         Degree of Brent-Suyama polynomial in stage 2.\n");
  printf (" <D>         Stride for roots of G in stage 2.\n");
  printf (" <startG>    Starting value for roots of G is <D> * <startG>.\n"
          "             Default is floor(<B1> / <D>)\n");
  printf (" <endG>      Ending value for roots of G is <D> * <endG>.\n"
          "             Default is floor(<B2> / <D>) + 1\n");
  printf (" <a>         Parameter for Dickson polynomial. If <a> == 0, <S>-th powers\n"
          "             are used. <a> = 0 is the default.\n");
  printf (" <maxBS>     Brent-Suyama extension is only tried on cofactors <= <maxBS>\n"
          "             Default is <B2> * 1000\n");
  printf (" -v          Verbose. Print info on each value found smooth.\n");
  printf (" -theo       Don't really compute and match Brent-Suyama values, just\n"
          "             calculate the probability of success and add up expected values.\n");
  printf (" -pm1        Investigate only those values in [N, N+nrtest-1] that are\n"
          "             a prime minus 1\n");
  printf( " -ecm        Adjust for the ECM algorithm (divides N by 12, doubles S)\n");
}

int 
main(int argc, char **argv)
{
  mpz_t N, a, *cofac;
  unsigned int len_cofac;
  unsigned int p, B1=0, i, D, S, startG , endG, Nmod12, blocklen;
  unsigned int nr_tests = 0, nr_blocks = 1, blockstart = 0, nr_primes = 0;
#ifdef PRIMEGEN
  primegen pg[1];
#endif
  double ppow, 
         B2=0., 
         maxBS = 0., /* Try Brent-Suyama only on cofactors <= maxBS */
         nrBSsmooth = 0.; /* Nr (or expected value) of Brent-Suyama successes */
  unsigned int nrB1smooth = 0, nrB2smooth = 0;
  int verbose = 0, theo = 0, ecm_style = 0, pminus1 = 0;
  mpz_t *F, *G;
  
  mpz_init (N);
  mpz_init (a);
  D = endG = S = 0;
  startG = 0xffffffff;
  B1 = 1000000;
  B2 = 100000000;
  nr_tests = 10000;
  mpz_set_str (N, "8333333333333333333333333331", 10);
  
  /* Get parameters */
  
  for (i = 1; i < (unsigned int) argc; i++)
    {
      if (strcmp (argv[i], "-v") == 0)
        verbose = 1;
      else if (strcmp (argv[i], "-theo") == 0)
        theo = 1;
      else if (strcmp (argv[i], "-ecm") == 0)
        ecm_style = 1;
      else if (strcmp (argv[i], "-pm1") == 0)
        pminus1 = 1;
      else if (getparm_ui (argc, argv, i, "-tests", &nr_tests))
        i++;
      else if (getparm_ui (argc, argv, i, "-blocks", &nr_blocks))
        i++;
      else if (getparm_ui (argc, argv, i, "-B1", &B1))
        i++;
      else if (getparm_d (argc, argv, i, "-B2", &B2))
        i++;
      else if (getparm_ui (argc, argv, i, "-D", &D))
        i++;
      else if (getparm_ui (argc, argv, i, "-S", &S))
        i++;
      else if (getparm_ui (argc, argv, i, "-startG", &startG))
        i++;
      else if (getparm_ui (argc, argv, i, "-endG", &endG))
        i++;
      else if (getparm_mpz (argc, argv, i, "-N", &N))
        i++;
      else if (getparm_mpz (argc, argv, i, "-a", &a))
        i++;
      else if (getparm_d (argc, argv, i, "-maxBS", &maxBS))
        i++;
      else
        {
          printf ("Don't know parameter %s\n", argv[i]);
          print_help();
          exit (EXIT_FAILURE);
        }
    }
  
  if (maxBS == 0.)
    maxBS = B2 * 1000.;
  
  if (startG == 0xffffffff && D > 1)
    startG = B1 / D;
  
  if (endG == 0 && D > 1)
    {
      if (B2 > (double)(~1U)*(double)D)
        {
          fprintf (stderr, 
            "Overflow error computing default endG. Please set explicitly\n");
          exit (EXIT_FAILURE);
        }
      endG = (unsigned int) (B2 / (double) D);
    }
  
  if (ecm_style)
    {
      /* ECM is known to have a group order divisible by 12, which 
         effectively makes the part that has yet to be smooth N/12.
         Instead of generating multiples of 12 and sieving those,
         we simply divide N by 12 here. (This leads to rounding down if N
         was not a multiple of 12, which we shouls correct when printing.)
         
         ECM finds a factor if p|f_S(x)-f_S(y) or p|f_S(x)+f_S(y), (because
         a point and it's negative have the same x-coordinate) where p is 
         the missing group order factor, and f_S(x) and f_S(y) are the 
         polynomial values examined by the Brent-Suyma extension in stage 2. 
         If p is prime, this is equivalent to the condition 
         p|(f_S(x)-f_S(y))(f_S(x)+f_S(y)), and for both f_S(x)=x^S and 
         f_S(x) = S-th Dickson polynomials, 
         (f_S(x)-f_S(y))(f_S(x)+f_S(y)) = f_{2S}(x)-f_{2S}(y).
         So effectively, ECM behaves as if S were twice as large, which we
         account for by simply doubling S here. */

      Nmod12 = mpz_fdiv_q_ui (N, N, 12);
      S = S * 2;
    }

  if (verbose)
    {
      printf ("B1=%u, B2=%.0f, %s%u, D=%u, %u<=G<=%u\nN=",
               B1, B2, mpz_sgn(a) ? "Dickson_" : "X^", S, D, startG, endG);
      mpz_out_str (stdout, 10, N);
      printf ("\n");
      fflush (stdout);
    }
  
  len_cofac = (nr_tests + nr_blocks - 1) / nr_blocks;
  cofac = (mpz_t *) malloc (len_cofac * sizeof (mpz_t));
  for (i = 0; i < len_cofac; i++)
    mpz_init (cofac[i]);

  F = (mpz_t *) malloc (get_lenF (D) * sizeof (mpz_t));
  for (i = 0; i < get_lenF (D); i++)
    mpz_init (F[i]);

  G = (mpz_t *) malloc ((endG - startG + 1) * sizeof (mpz_t));
  for (i = 0; i < (endG - startG + 1); i++)
    mpz_init (G[i]);

  while (nr_tests > 0)
    {
      if (nr_tests < len_cofac)
        blocklen = nr_tests;
      else
        blocklen = len_cofac;
      
      /* Init cofac */
      for (i = 0; i < blocklen; i++)
        mpz_add_ui (cofac[i], N, i);

      /* Do the sieving */
#ifdef PRIMEGEN
      primegen_init (pg);
      for (p = primegen_next (pg); p <= B1; p = primegen_next (pg))
#else
      for (p = 2; p <= B1; p = getprime (p))
#endif
        {
          /* Compute first sieve location where p divides */
          i = mpz_fdiv_ui (N, p);
          if (i > 0)
            i = p - i;

          while (i < blocklen)
            {
              if (mpz_fdiv_q_ui (cofac[i], cofac[i], p) != 0)
                {
                  fprintf (stderr, "%d does not divide cofac[%d]\n", p, i);
                  exit (EXIT_FAILURE);
                }
              /* Take out prime powers. Not very efficient. */
              ppow = (double)p * (double)p;
              while (ppow <= B1 && mpz_fdiv_ui (cofac[i], p) == 0)
                {
                  mpz_fdiv_q_ui (cofac[i], cofac[i], p);
                  ppow *= (double) p;
                }
              i += p;
            }
        }

      /* Check which cofactors are small enough and report */
      for (i = 0; i < blocklen; i++)
        {
          if (pminus1) /* Report only primes - 1 if pminus1 is set */
            { 
              /* Test if N+i+1 is prime. If it is, obviously no factors
                 can have been divided out of it, and the cofactor must
                 still be larger than N. */
              if ((i + 1 == blocklen || mpz_cmp (cofac[i + 1], N) > 0) &&
                  is_P_minus_i (N, i+1))
                nr_primes++;
              else
                continue;
            }
          
          if (mpz_cmp_ui (cofac[i], 1) == 0)
            {
              if (verbose)
                printf ("N%s+%d: B1-smooth\n", ecm_style ? "/12" : "", i+blockstart);

              nrB1smooth++;
              continue;
            }

          if (mpz_cmp_d (cofac[i], B2) <= 0)
            {
              if (verbose)
                {
                  printf ("N%s+%d: ", ecm_style ? "/12" : "", i + blockstart);
                  mpz_out_str (stdout, 10, cofac[i]);
                  printf ("\n");
                }

              nrB2smooth++;
              continue;
            }

          if (mpz_cmp_d (cofac[i], maxBS) <= 0)
            {
              if (theo)
                {
                  double pr;
                  pr = brent_suyama_theo (cofac[i], S, D, startG, endG, a);
                  nrBSsmooth += pr;
                  if (verbose && pr >= 0.3)
                    {
                      printf ("N%s+%d: ", ecm_style ? "/12" : "", i + blockstart);
                      mpz_out_str (stdout, 10, cofac[i]);
                      printf (" (Brent-Suyama, pr=%f)\n", pr);
                    }
                } else {
                  unsigned int luckyF, luckyG;
                  if (brent_suyama_match (cofac[i], S, D, startG, endG, a, F, G, &luckyF, &luckyG))
                    {
                      unsigned int j;
                      
                      for (j = 2; j < S; j += 2)
                        if (S % j == 0 && 
                            brent_suyama_match (cofac[i], j, D, startG, endG, a, F, G, &luckyF, &luckyG))
                          break;
                      
                      if (verbose)
                        {
                          printf ("N%s+%d: ", ecm_style ? "/12" : "", i+blockstart);
                          mpz_out_str (stdout, 10, cofac[i]);
                          printf (" (Brent-Suyama, divides ");
                          if (mpz_sgn(a))
                            printf( "Dickson_%u(D*%u)-Dickson_%u(%u))\n", j, luckyG, j, luckyF);
                          else
                            printf( "(D*%u)^%u-%u^%u)\n", luckyG, j, luckyF, j);
                        }

                      nrBSsmooth++;
                      continue;
                    }
                }
            }
        }
      
      nr_tests -= blocklen;
      blockstart += blocklen;
      mpz_add_ui (N, N, blocklen);
#ifndef PRIMEGEN
      p = getprime (FREE_PRIME_TABLE);
#endif
    }
  
  for (i = 0; i < get_lenF (D); i++)
    mpz_clear (F[i]);
  free (F);

  for (i = 0; i < (endG - startG + 1); i++)
    mpz_clear (G[i]);
  free (G);

  printf ("B1-smooth: %d, B2-smooth: %d, found by Brent-Suyama: %f, Total: %d\n", 
    nrB1smooth, nrB2smooth, nrBSsmooth ,nrB1smooth + nrB2smooth + 
    (unsigned int) nrBSsmooth);

  if (pminus1)
    printf ("Number of N that are a prime - 1: %d\n", nr_primes);

  return 0;
}
