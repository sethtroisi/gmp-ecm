/* Dickman's rho function (to compute probability of success of ecm).

  Copyright 2004, 2005, 2007, 2008 Alexander Kruppa.

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

#include "config.h"
#if defined(TESTDRIVE)
#define _ISOC99_SOURCE 1
#endif
#if defined(DEBUG_NUMINTEGRATE) || defined(TESTDRIVE)
# include <stdio.h>
#endif
#include <stdlib.h>
#include <math.h>
#if defined(TESTDRIVE)
#include <string.h>
#include "primegen.h"
#endif
#ifdef HAVE_LIBGSL
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_integration.h>
#endif
#include "ecm-impl.h"

#ifndef ECM_EXTRA_SMOOTHNESS
#define ECM_EXTRA_SMOOTHNESS 23.4
#endif

#define M_PI_SQR   9.869604401089358619 /* Pi^2 */
#define M_PI_SQR_6 1.644934066848226436 /* Pi^2/6 */
/* gsl_math.h defines M_EULER */
#ifndef M_EULER
#define M_EULER    0.577215664901532861
#endif
#define M_EULER_1   0.422784335098467139 /* 1 - Euler */

void rhoinit (int, int); /* used in stage2.c */

static double *rhotable = NULL;
static int invh = 0;
static double h = 0.;
static int tablemax = 0;
#if defined(TESTDRIVE)
#define PRIME_PI_MAX 10000
#define PRIME_PI_MAP(x) (((x)+1)/2)
/* The number of primes up to i. Use prime_pi[PRIME_PI_MAP(i)].
   Only correct for i >= 2. */
static unsigned int prime_pi[PRIME_PI_MAP(PRIME_PI_MAX)+1];
#endif

#ifdef TESTDRIVE
unsigned long
gcd (unsigned long a, unsigned long b)
{
  unsigned long t;

  while (b != 0)
    {
      t = a % b;
      a = b;
      b = t;
    }

  return a;
}

unsigned long
eulerphi (unsigned long n)
{
  unsigned long phi = 1, p;

  for (p = 2; p * p <= n; p += 2)
    {
      if (n % p == 0)
        {
          phi *= p - 1;
          n /= p;
          while (n % p == 0)
            {
              phi *= p;
              n /= p;
            }
        }

      if (p == 2)
        p--;
    }

  /* now n is prime */

  return (n == 1) ? phi : phi * (n - 1);
}


/* The number of positive integers up to x that have no prime factor up to y,
   for x >= y >= 2. Uses Buchstab's identity */
unsigned long
Buchstab_Phi(unsigned long x, unsigned long y) 
{
  unsigned long p, s;
  primegen pg[1];

  if (x < 1)
    return 0;
  if (x <= y)
    return 1;
#if 0
  if (x < y^2)
    return(1 + primepi(x) - primepi (y)));
#endif

  s = 1;
  primegen_init (pg);
  primegen_skipto (pg, y + 1);
  for (p = primegen_next(pg); p <= x; p = primegen_next(pg))
    s += Buchstab_Phi(x / p, p - 1);
  return (s);
}


/* The number of positive integers up to x that have no prime factor
   greter than y, for x >= y >= 2. Uses Buchstab's identity */
unsigned long 
Buchstab_Psi(const unsigned long x, const unsigned long y) 
{
  unsigned long r, p;
  primegen pg[1];

  if (x <= y)
    return (x);

  if (y == 1UL)
    return (1);

  /* If y^2 > x, then
     Psi(x,y) = x - \sum_{y < p < x, p prime} floor(x/p)

     We separate the sum into ranges where floor(x/p) = k,
     which is x/(k+1) < p <= x/k.
     We also need to satisfy y < p, so we need k < x/y - 1,
     or k_max = ceil (x/y) - 2.
     The primes y < p <= x/(k_max + 1) are summed separately. */
  if (x <= PRIME_PI_MAX && x < y * y)
    {
      unsigned long kmax = x / y - 1;
      unsigned long s1, s2, k;
      
        s1 = (kmax + 1) * (prime_pi [PRIME_PI_MAP(x / (kmax + 1))] - 
                           prime_pi [PRIME_PI_MAP(y)]);
        s2 = 0;
        for (k = 1; k <= kmax; k++)
          s2 += prime_pi[PRIME_PI_MAP(x / k)];
        s2 -= kmax * prime_pi [PRIME_PI_MAP(x / (kmax+1))];
        return (x - s1 - s2);
    }

  r = 1;
  primegen_init (pg);
  for (p = primegen_next(pg); p <= y; p = primegen_next(pg))
    r += Buchstab_Psi (x / p, p);
  return (r);
}

#endif /* TESTDRIVE */


#ifdef HAVE_LIBGSL
static double
Li (const double x)
{
  return (- gsl_sf_expint_E1 (- log(x)));
}
#endif

/*
  Evaluate dilogarithm via the sum 
  \Li_{2}(z)=\sum_{k=1}^{\infty} \frac{z^k}{k^2}, 
  see http://mathworld.wolfram.com/Dilogarithm.html
  Assumes |z| <= 0.5, for which the sum converges quickly.
 */

static double
dilog_series (const double z)
{
  double r = 0.0, zk; /* zk = z^k */
  int k, k2; /* k2 = k^2 */
  /* Doubles have 53 bits in significand, with |z| <= 0.5 the k+1-st term
     is <= 1/(2^k k^2) of the result, so 44 terms should do */
  for (k = 1, k2 = 1, zk = z; k <= 44; k2 += 2 * k + 1, k++, zk *= z)
    r += zk / (double) k2;

  return r;
}

static double
dilog (double x)
{
  ASSERT(x <= -1.0); /* dilog(1-x) is called from rhoexact for 2 < x <= 3 */

  if (x <= -2.0)
    return -dilog_series (1./x) - M_PI_SQR_6 - 0.5 * log(-1./x) * log(-1./x);
  else /* x <= -1.0 */
    {
      /* L2(z) = -L2(1 - z) + 1/6 * Pi^2 - ln(1 - z)*ln(z) 
         L2(z) = -L2(1/z) - 1/6 * Pi^2 - 0.5*ln^2(-1/z)
         ->
         L2(z) = -(-L2(1/(1-z)) - 1/6 * Pi^2 - 0.5*ln^2(-1/(1-z))) + 1/6 * Pi^2 - ln(1 - z)*ln(z)
               = L2(1/(1-z)) - 1/6 * Pi^2 + 0.5*ln(1 - z)^2 - ln(1 - z)*ln(-z)
         z in [-1, -2) -> 1/(1-z) in [1/2, 1/3)
      */
      double log1x = log (1. - x);
      return dilog_series (1. / (1. - x)) 
             - M_PI_SQR_6 + log1x * (0.5 * log1x - log (-x));
    }
}

#if 0
static double 
L2 (double x)
{
  return log (x) * (1 - log (x-1)) + M_PI_SQR_6 - dilog (1 - x);
}
#endif

static double
rhoexact (double x)
{
  ASSERT(x <= 3.);
  if (x <= 0.)
    return 0.;
  if (x <= 1.)
    return 1.;
  if (x <= 2.)
    return 1. - log (x);
  if (x <= 3.) /* 2 < x <= 3 thus -2 <= 1-x < -1 */
    return 1. - log (x) * (1. - log (x - 1.)) + dilog (1. - x) + 0.5 * M_PI_SQR_6;
  
  return 0.; /* x > 3. and asserting not enabled: bail out with 0. */
}


/* The Buchstab omega(x) function, exact for x <= 4 where it can be 
   evaluated without numerical integration, and approximated by 
   exp(gamma) for larger x. */

static double
omega (const double x)
{
  /* magic = dilog(-1) + 1  = Pi^2/12 + 1 */
  const double magic = 1.82246703342411321824; 

  if (x < 1.) return (0.);
  if (x <= 2.) return (1. / x);
  if (x <= 3.) return ((log (x - 1.) + 1.) / x);
  if (x <= 4.)
    return ((dilog(2. - x) + (1. + log(x - 2.)) * log(x - 1.) + magic) / x);

  /* If argument is out of range, return the limiting value for 
     $x->\infty$: e^-gamma. 
     For x only a little larger than 4., this has relative error 2.2e-6,
     for larger x the error rapidly drops further */

  return 0.56145948356688516982;
}


void 
rhoinit (int parm_invh, int parm_tablemax)
{
  int i;

  if (parm_invh == invh && parm_tablemax == tablemax)
    return;

  if (rhotable != NULL)
    {
      free (rhotable);
      rhotable = NULL;
      invh = 0;
      h = 0.;
      tablemax = 0;
    }
  
  /* The integration below expects 3 * invh > 4 */
  if (parm_tablemax == 0 || parm_invh < 2)
    return;
    
  invh = parm_invh;
  h = 1. / (double) invh;
  tablemax = parm_tablemax;
  
  rhotable = (double *) malloc (parm_invh * parm_tablemax * sizeof (double));
  if (rhotable == NULL)
    {
      fprintf (stderr, "Cannot allocate memory in rhoinit\n");
      exit (1);
    }
  
  for (i = 0; i < (3 < parm_tablemax ? 3 : parm_tablemax) * invh; i++)
    rhotable[i] = rhoexact (i * h);
  
  for (i = 3 * invh; i < parm_tablemax * invh; i++)
    {
      /* rho(i*h) = 1 - \int_{1}^{i*h} rho(x-1)/x dx
                  = rho((i-4)*h) - \int_{(i-4)*h}^{i*h} rho(x-1)/x dx */
      
      rhotable[i] = rhotable[i - 4] - 2. / 45. * (
          7. * rhotable[i - invh - 4] / (double)(i - 4)
        + 32. * rhotable[i - invh - 3] / (double)(i - 3)
        + 12. * rhotable[i - invh - 2] / (double)(i - 2)
        + 32. * rhotable[i - invh - 1] / (double)(i - 1)
        + 7. * rhotable[i - invh]  / (double)i );
      if (rhotable[i] < 0.)
        {
#ifndef DEBUG_NUMINTEGRATE
          rhotable[i] = 0.;
#else
          printf (stderr, "rhoinit: rhotable[%d] = %.16f\n", i, 
                   rhotable[i]);
          exit (EXIT_FAILURE);
#endif
        }
    }
}

static double
dickmanrho (double alpha)
{
  if (alpha <= 3.)
     return rhoexact (alpha);
  if (alpha < tablemax)
    {
      int a = floor (alpha * invh);
      double rho1 = rhotable[a];
      double rho2 = (a + 1) < tablemax * invh ? rhotable[a + 1] : 0;
      return rho1 + (rho2 - rho1) * (alpha * invh - (double)a);
    }
  
  return 0.;
}

static double 
dickmanrhosigma (double alpha, double x)
{
  if (alpha <= 0.)
    return 0.;
  if (alpha <= 1.)
    return 1.;
  if (alpha < tablemax)
    return dickmanrho (alpha) + M_EULER_1 * dickmanrho (alpha - 1.) / log (x);
  
  return 0.;
}

#if 0
static double
dickmanrhosigma_i (int ai, double x)
{
  if (ai <= 0)
    return 0.;
  if (ai <= invh)
    return 1.;
  if (ai < tablemax * invh)
    return rhotable[ai] - M_EULER * rhotable[ai - invh] / log(x);
  
  return 0.;
}
#endif

static double
dickmanlocal (double alpha, double x)
{
  if (alpha <= 0.)
    return 0.;
  if (alpha <= 1.)
    return 1.;
  if (alpha < tablemax)
    return dickmanrhosigma (alpha, x) 
           - dickmanrhosigma (alpha - 1., x) / log (x);
  return 0.;
}

static double
dickmanlocal_i (int ai, double x)
{
  if (ai <= 0)
    return 0.;
  if (ai <= invh)
    return 1.;
  if (ai <= 2 * invh && ai < tablemax * invh)
    return rhotable[ai] - M_EULER / log (x);
  if (ai < tablemax * invh)
    {
      double logx = log (x);
      return rhotable[ai] - (M_EULER * rhotable[ai - invh]
             + M_EULER_1 * rhotable[ai - 2 * invh] / logx) / logx;
    }

  return 0.;
}

static double
dickmanmu (double alpha, double beta, double x)
{
  double a, b, sum;
  int ai, bi, i;
  ai = ceil ((alpha - beta) * invh);
  if (ai > tablemax * invh)
    ai = tablemax * invh;
  a = (double) ai * h;
  bi = floor ((alpha - 1.) * invh);
  if (bi > tablemax * invh)
    bi = tablemax * invh;
  b = (double) bi * h;
  sum = 0.;
  for (i = ai + 1; i < bi; i++)
    sum += dickmanlocal_i (i, x) / (alpha - i * h);
  sum += 0.5 * dickmanlocal_i (ai, x) / (alpha - a);
  sum += 0.5 * dickmanlocal_i (bi, x) / (alpha - b);
  sum *= h;
  sum += (a - alpha + beta) * 0.5 * (dickmanlocal_i (ai, x) / (alpha - a) + dickmanlocal (alpha - beta, x) / beta);
  sum += (alpha - 1. - b) * 0.5 * (dickmanlocal (alpha - 1., x) + dickmanlocal_i (bi, x) / (alpha - b));

  return sum;
}

static double
brentsuyama (double B1, double B2, double N, double nr)
{
  double a, alpha, beta, sum;
  int ai, i;
  alpha = log (N) / log (B1);
  beta = log (B2) / log (B1);
  ai = floor ((alpha - beta) * invh);
  if (ai > tablemax * invh)
    ai = tablemax * invh;
  a = (double) ai * h;
  sum = 0.;
  for (i = 1; i < ai; i++)
    sum += dickmanlocal_i (i, N) / (alpha - i * h) * (1 - exp (-nr * pow (B1, (-alpha + i * h))));
  sum += 0.5 * (1 - exp(-nr / pow (B1, alpha)));
  sum += 0.5 * dickmanlocal_i (ai, N) / (alpha - a) * (1 - exp(-nr * pow (B1, (-alpha + a))));
  sum *= h;
  sum += 0.5 * (alpha - beta - a) * (dickmanlocal_i (ai, N) / (alpha - a) + dickmanlocal (alpha - beta, N) / beta);

  return sum;
}

static double 
brsudickson (double B1, double B2, double N, double nr, int S)
{
  int i, f;
  double sum;
  sum = 0;
  f = eulerphi (S) / 2;
  for (i = 1; i <= S / 2; i++)
      if (gcd (i, S) == 1)
        sum += brentsuyama (B1, B2, N, nr * (gcd (i - 1, S) + gcd (i + 1, S) - 4) / 2);
  
  return sum / (double)f;
}

static double
brsupower (double B1, double B2, double N, double nr, int S)
{
  int i, f;
  double sum;
  sum = 0;
  f = eulerphi (S);
  for (i = 1; i < S; i++)
      if (gcd (i, S) == 1)
        sum += brentsuyama (B1, B2, N, nr * (gcd (i - 1, S) - 2));
  
  return sum / (double)f;
}

/* Assume N is as likely smooth as a number around N/delta */

static double
prob (double B1, double B2, double N, double nr, int S, double delta)
{
  double alpha, beta, stage1, stage2, brsu;
  const double effN = N / delta;

  ASSERT(rhotable != NULL);
  
  /* What to do if rhotable is not initialised and asserting is not enabled?
     For now, bail out with 0. result. Not really pretty, either */
  if (rhotable == NULL)
    return 0.;

  if (B1 < 2. || N <= 1.)
    return 0.;
  
  if (effN <= B1)
    return 1.;

#ifdef TESTDRIVE
  printf ("B1 = %f, B2 = %f, N = %.0f, nr = %f, S = %d\n", B1, B2, N, nr, S);
#endif
  
  alpha = log (effN) / log (B1);
  stage1 = dickmanlocal (alpha, effN);
  stage2 = 0.;
  if (B2 > B1)
    {
      beta = log (B2) / log (B1);
      stage2 = dickmanmu (alpha, beta, effN);
    }
  brsu = 0.;
  if (S < -1)
    brsu = brsudickson (B1, B2, effN, nr, -S * 2);
  if (S > 1)
    brsu = brsupower (B1, B2, effN, nr, S * 2);

#ifdef TESTDRIVE
  printf ("stage 1 : %f, stage 2 : %f, Brent-Suyama : %f\n", stage1, stage2, brsu);
#endif

  return (stage1 + stage2 + brsu) > 0. ? (stage1 + stage2 + brsu) : 0.;
}

double
ecmprob (double B1, double B2, double N, double nr, int S)
{
  return prob (B1, B2, N, nr, S, ECM_EXTRA_SMOOTHNESS);
}

double
pm1prob (double B1, double B2, double N, double nr, int S, const mpz_t go)
{
  mpz_t cof;
  /* A prime power q^k divides p-1, p prime, with probability 1/(q^k-q^(k-1))
     not with probability 1/q^k as for random numbers. This is taken into 
     account by the "smoothness" value here; a prime p-1 is about as likely
     smooth as a random number around (p-1)/exp(smoothness).
     smoothness = \sum_{q in Primes} log(q)/(q-1)^2 */
  double smoothness = 1.2269688;
  unsigned long i;
  
  if (go != NULL && mpz_cmp_ui (go, 1UL) > 0)
    {
      mpz_init (cof);
      mpz_set (cof, go);
      for (i = 2; i < 100; i++)
        if (mpz_divisible_ui_p (cof, i))
          {
            /* If we know that q divides p-1 with probability 1, we need to
               adjust the smoothness parameter */
            smoothness -= log ((double) i) / (double) ((i-1)*(i-1));
            /* printf ("pm1prob: Dividing out %lu\n", i); */
            while (mpz_divisible_ui_p (cof, i))
              mpz_tdiv_q_ui (cof, cof, i);
          }
      mpz_clear (cof);
      /* printf ("pm1prob: smoothness after dividing out go primes < 100: %f\n", 
               smoothness); */
      return prob (B1, B2, N, nr, S, exp(smoothness) * mpz_get_d (go));
    }

  return prob (B1, B2, N, nr, S, exp(smoothness));
}


/* Compute probability for primes p == r (mod m) */

double
pm1prob_rm (double B1, double B2, double N, double nr, int S, unsigned long r,
            unsigned long m)
{
  unsigned long cof;
  double smoothness = 1.2269688;
  unsigned long p;
  
  cof = m;
  
  for (p = 2UL; p < 100UL; p++)
    if (cof % p == 0UL) /* For each prime in m */
      {
        unsigned long cof_r, k, i;
        /* Divisibility by i is determined by r and m. We need to
           adjust the smoothness parameter. In P-1, we had estimated the 
           expected value for the exponent of p as p/(p-1)^2. Undo that. */
        smoothness -= (double)p / ((p-1)*(p-1)) * log ((double) p);
        /* The expected value for the exponent of this prime is k s.t.
           p^k || r, plus 1/(p-1) if p^k || m as well */
        cof_r = gcd (r - 1UL, m);
        for (k = 0UL; cof_r % p == 0UL; k++)
          cof_r /= p;
        smoothness += k * log ((double) p);

        cof_r = m;
        for (i = 0UL; cof_r % p == 0UL; i++)
          cof_r /= p;

        if (i == k)
          smoothness += (1./(p - 1.) * log ((double) p));
        
        while (cof % p == 0UL)
          cof /= p;
        printf ("pm1prob_rm: p = %lu, k = %lu, i = %lu, new smoothness = %f\n", 
                p, i, k, smoothness); 
      }

  return prob (B1, B2, N, nr, S, exp(smoothness));
}


/* The \Phi(x,y) function gives the number of natural numbers <= x 
   that have no prime factor <= y, see Tenenbaum, 
   "Introduction the analytical and probabilistic number theory", III.6.
   This function estimates the \Phi(x,y) function via eq. (48) of the 1st
   edition resp. equation (6.49) of the 3rd edition of Tenenbaum's book. */

#ifdef HAVE_LIBGSL

static double 
integrand1 (double x, void *params)
{
  return pow (*(double *) params, x) / x * log(x-1.);
}


double 
no_small_prime (double x, double y)
{
  double u;
  ASSERT (x >= 2.);
  ASSERT (y >= 2.);
  if (x <= y)
    return 0.;
  
  u = log(x)/log(y);

   /* If no prime factors <= sqrt(x), number must be a prime > y */
  if (u <= 2)
    return (Li(x) - Li(y));
  
  if (u <= 3)
    {
      double r, abserr;
      size_t neval;
      gsl_function f;

      f.function = &integrand1;
      f.params = &y;

      /* intnum(v=1,u,buchstab(v)/y^(u-v)) */

      /* First part: intnum(v=2, u, y^v/v*log(v-1.)) */
      gsl_integration_qng (&f, 2., u, 0., 0.01, &r, &abserr, &neval);

      /* Second part: intnum(v=1, u, y^v/v) */
      r += -gsl_sf_expint_E1(-log(y)*u) + gsl_sf_expint_E1(-log(y));
      
      r *= x / pow(y,u);
      return r;
    }

  return x * omega(u) / log(y);
}
#endif


#ifdef TESTDRIVE
int
main (int argc, char **argv)
{
  double B1, B2, N, nr, r, m;
  int S;
  unsigned long p, i, pi;
  primegen pg[1];
  
  primegen_init (pg);
  i = pi = 0;
  for (p = primegen_next (pg); p <= PRIME_PI_MAX; p = primegen_next (pg))
    {
      for ( ; i < p; i++)
        prime_pi[PRIME_PI_MAP(i)] = pi;
      pi++;
    }
  for ( ; i < p; i++)
    prime_pi[PRIME_PI_MAP(i)] = pi;
  

  if (argc < 2)
    {
      printf ("Usage: rho <B1> <B2> <N> <nr> <S> [<r> <m>]\n");
      return 1;
    }
  
  if (strcmp (argv[1], "-Buchstab_Phi") == 0)
    {
      unsigned long x, y, r;
      if (argc < 4)
        {
          printf ("-Buchstab_Phi needs x and y paramters\n");
          exit (EXIT_FAILURE);
        }
      x = strtoul (argv[2], NULL, 10);
      y = strtoul (argv[3], NULL, 10);
      r = Buchstab_Phi (x, y);
      printf ("Buchstab_Phi (%lu, %lu) = %lu\n", x, y, r);
      exit (EXIT_SUCCESS);
    }
  else if (strcmp (argv[1], "-Buchstab_Psi") == 0)
    {
      unsigned long x, y, r;
      if (argc < 4)
        {
          printf ("-Buchstab_Psi needs x and y paramters\n");
          exit (EXIT_FAILURE);
        }
      x = strtoul (argv[2], NULL, 10);
      y = strtoul (argv[3], NULL, 10);
      r = Buchstab_Psi (x, y);
      printf ("Buchstab_Psi (%lu, %lu) = %lu\n", x, y, r);
      exit (EXIT_SUCCESS);
    }
  else if (strcmp (argv[1], "-nsp") == 0)
    {
      double x, y, r;
      
      if (argc < 4)
        {
          printf ("-nsp needs x and y paramters\n");
          exit (EXIT_FAILURE);
        }
      x = atof (argv[2]);
      y = atof (argv[3]);
      r = no_small_prime (x, y);
      printf ("no_small_prime(%f, %f) = %f\n", x, y, r);
      exit (EXIT_SUCCESS);
    }


  if (argc < 6)
    {
      printf ("Need 5 or 7 arguments: B1 B2 N nr S [r m]\n");
      exit (EXIT_FAILURE);
    }
  
  B1 = atof (argv[1]);
  B2 = atof (argv[2]);
  N = atof (argv[3]);
  nr = atof (argv[4]);
  S = atoi (argv[5]);
  r = 0; m = 1;
  if (argc > 7)
    {
      r = atoi (argv[6]);
      m = atoi (argv[7]);
    }

  rhoinit (256, 10);
  if (N < 50.)
    {
      double sum;
      sum = ecmprob(B1, B2, exp2 (N), nr, S);
      sum += 4. * ecmprob(B1, B2, 3./2. * exp2 (N), nr, S);
      sum += ecmprob(B1, B2, 2. * exp2 (N), nr, S);
      sum *= 1./6.;
      printf ("ECM: %.16f\n", sum);

      sum = pm1prob_rm (B1, B2, exp2 (N), nr, S, r, m);
      sum += 4. * pm1prob_rm (B1, B2, 3./2. * exp2 (N), nr, S, r, m);
      sum += pm1prob_rm (B1, B2, 2. * exp2 (N), nr, S, r, m);
      sum *= 1./6.;
      printf ("P-1: %.16f\n", sum);
    }
  else
    {
      printf ("ECM: %.16f\n", ecmprob(B1, B2, N, nr, S));
      printf ("P-1: %.16f\n", pm1prob_rm (B1, B2, N, nr, S, r, m));
    }
  rhoinit (0, 0);
  return 0;
}
#endif
