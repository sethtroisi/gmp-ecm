/*

  Estimating the probability of success of the Elliptic Curve Method
  
  Copyright 2004, 2005 Alexander Kruppa.
  
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

  Version 0.1.1

  History: 0.1   2004.09.29
                 Initial release
           0.1.1 2005.06.14 
                 Changed extra smoothness of ECM from 12 to 23.4
                 Started on pm1prob(), but incomplete yet

  How to use this file:

  0) first start gp (http://pari.math.u-bordeaux.fr/)

  1) load the file into gp:

     ? read ("rho.gp")

  2) then call 
  
     ? ecmprob(B1,B2,N,nr,S)
     
     where N is the approximate factor, nr is the number of points evaluated
     in stage 2, i.e. nr = k*dF^2, S is the degree for Brent-Suyama. 
     Passing S>0 means S-th powers, -S means Dickson polynomials, S=0 means 
     no Brent-Suyama extension will be considered.

     The parameters can be obtained by running gmp-ecm with the -v 
     command line option. I.e. for B1=1000000
     
     $ echo 65537 | ecm5 -v 1e6-1e6
     GMP-ECM 5.0.3 [powered by GMP 4.1.2] [ECM]
     Input number is 65537 (5 digits)
     Using special division for factor of 2^16+1
     Using B1=1000000, B2=839549780, polynomial Dickson(6), sigma=2157207190
     a=51553
     starting point: x=19319
     Step 1 took 0ms
     x=19319
     B2'=948726240 k=5 b2=189560910 d=43890 dF=4320
     Initializing table of differences for F took 0ms
     Found factor while computing F[49]
     Step 2 took 20ms for 0 muls
     ********** Factor found in step 2: 65537
     Found input number N
     
     So gmp-ecm internally uses B2=948726240, k=5, dF=4320 and 
     S=-6 (Dickson(6)).
     
     Passing these paramters to ecmprob, we get
     
     ? ecmprob(1000000,948726240,sqrt(10)*10^34,5*4320^2,-6)
     %5 = 0.0009601784546838897811362317141
     ? 1./%
     %6 = 1041.473066930272122308814268
     ?
     
     Thus, the expected number of curves to find a p35 factor with 
     gmp-ecm and B1=1000000 is approximately 1041.

*/


/* Returns \int_{1}^{x} 1 - rho(t-1)/t for 1 <= x <= 2.
   This function is not actually used any more, rhoexact() needs a 
   differece L2(x)-L2(2) which can be simplified. */
L2 (x) =
{
/* \int (1 - Log(x-1))/x dx = 
   Log(x) + Log(x)*(Log(1-x) - Log(x-1)) + Dilog(x)
   Dilog(x) = Pi^2/6 - Log(x)*Log(1-x) - Dilog(1-x)
   thus:  L2(x) = Log(x) - Log(x)*Log(x-1) + Pi^2/6 - Dilog(1-x)
   L2(2) = Pi^2/4 + Log(2) */

  return (log (x) * (1 - log (x-1)) + Pi ^ 2 / 6 - real (dilog (1 - x)))
}

/*
L3(x) =
{
 \int L2(x-1)/x dx =
 \int (1 - Log[x-1] + Log[x-1]*Log[x-2] + PolyLog[2,2-x] + Pi^2/12) / x
 = Log[x] 
   + Log(x) * (Log(1-x) - Log(x-1)) + Dilog(x) 
   + (huge term, fearsome)
   - PolyLog[3,2-x]
   + Pi^2 / 12 * Log[x]
}
*/

rhoexact(x) = 
{
  if (x <= 0., return (0.));
  if (x <= 1., return (1.));

/* 1 - \int_{1}^{x} rho(t-1)/t dt = 
   1 - \int_{1}^{x} 1/t dt = 
   1 - (log(x) - log(1)) = 1 - log(x) */

  if (x <= 2., 
    return (1. - log(x)));

/* For 2 <= x <= 3, 
   1 - \int_{1}^{x} rho(t-1)/t dt =
   1 - \int_{1}^{2} 1/t dt - \int_{2}^{x} (1-log(t-1))/t dt =
   1 - log(2) - (L2(x) - L2(2))
   simplified, see L2() function. The real() is there because Pari returns a
   complex result for dilog(x), x<~-0.7, with very small imaginary part, even
   though the result should be purely real */

  if (x <= 3., 
    return (1. - log (x) * (1. - log (x - 1.)) + real (dilog (1. - x)) + Pi ^ 2 / 12.));

  error ("rhoexact: argument > 3");
}

/* With invh = 200, rho(8) = 0.000000032319, Knuth/Trapp-Pardo say ...21,
   rho(9) = 0.000000001015, Knuth/Trapp-Pardo say ...16
   With invh = 400, all digits match Knuth/Trapp-Pardo (after rounding) */

tablemax = 10;
invh = 256;
h = 1. / invh;
rhotable = listcreate (tablemax * invh);
for (i = 1, 3 * invh - 1, listput (rhotable, rhoexact (i * h), i))
/* FIXME: add listput (rhotable, rhoexact (3.), 3 * invh) here? */
/* Boole's rule. The h conveniently cancel */
for (i = 3 * invh, tablemax * invh, \
  listput (rhotable, rhotable[i - 4] - 2. / 45. * \
                    ( \
                        7. * rhotable[i - invh - 4] / (i - 4.) \
                      + 32. * rhotable[i - invh - 3] / (i - 3.) \
                      + 12. * rhotable[i - invh - 2] / (i - 2.) \
                      + 32. * rhotable[i - invh - 1] / (i - 1.) \
                      + 7. * rhotable[i - invh]  / i \
                    ), i \
  ) \
)

/* The rho function as defined by Karl Dickman, or by Knuth/Trapp-Pardo (4.1)-(4.4), 
   for alpha < tablemax. For alpha >= tablemax, returns 0. */
dickmanrho (alpha) =
{
  local (a, rho1, rho2);
  if (alpha <= 0., return (0.));
  if (alpha <= 3., return (rhoexact (alpha)));
  if (alpha < tablemax, 
    a = floor (alpha * invh);
    rho1 = rhotable[a];
    rho2 = rhotable[a + 1];
    /* Linear interpolation. Should use a better model */
    return (rho1 + (rho2 - rho1) * (alpha * invh - a));
  );
  return (0.);
}

/* The density of x^(1/alpha)-smooth positive integers below x with first 
   correction term, (4.8), (4.15) */
dickmanrhosigma (alpha, x) =
{
  if (alpha <= 0., return (0.));
  if (alpha <= 1., return (1.));
  if (alpha < tablemax, 
    return (dickmanrho (alpha) + (1. - Euler) * dickmanrho (alpha - 1.) / log (x));
  );
  return (0.);
}

/* Same, but ai is an index to rhotable, i.e. ai*h = alpha */
dickmanrhosigma_i (ai, x) =
{
  if (ai <= 0, return (0.));
  if (ai <= invh, return (1.));
  if (ai <= tablemax * invh, 
    return (rhotable[ai] + (1. - Euler) * rhotable[ai-invh] / log (x))
  );
  return (0.);
}

/* The density of x^(1/alpha)-smooth integers around x */
dickmanlocal (alpha, x) =
{
  if (alpha <= 0., return (0.));
  if (alpha <= 1., return (1.));
  /* Avoid case where alpha >= tablemax, but alpha - 1 < tablemax which 
     would give negative result */
  if (alpha < tablemax,
    return (dickmanrhosigma (alpha, x) - dickmanrhosigma (alpha - 1, x) / log (x))
  );
  return (0);
}

dickmanlocal_i (ai, x) =
{
  if (ai <= 0, return (0.));
  if (ai <= invh, return (1.));
  if (ai <= 2 * invh, 
/* dickmanrhosigma_i(ai, x) - dickmanrhosigma_i(ai-invh, x)/log(x), simplified */
    return (rhotable[ai] - Euler / log (x))
  );
  if (ai <= tablemax * invh, 
    return (
      rhotable[ai] - (Euler * rhotable[ai - invh] 
      + (1. - Euler) * rhotable[ai - 2 * invh] / log (x)) / log (x)
    )
  );
  return (0);
}

/* Probability that a number around x has all prime factors <=x^(1/alpha), 
   and exactly one >x^(1/alpha), <=x^(beta/alpha) */
dickmanmu (alpha, beta, x) =
{
  local (a, ai, b, bi);
  ai = ceil ((alpha - beta) * invh);
  a = ai * h;
  bi = floor ((alpha - 1.) * invh);
  b = bi * h;
  return (
    h * (
/* Trapezoidal rule. Could be improved */
      sum (i = ai, bi, dickmanlocal_i (i, x) / (alpha - i * h))
      - (dickmanlocal_i (ai, x) / (alpha - a) + dickmanlocal_i (bi, x) / (alpha - b)) / 2.
    )
    + (a - alpha + beta) * (dickmanlocal_i (ai, x) / (alpha - a) + dickmanlocal (alpha - beta, x) / beta) / 2.
    + (alpha - 1. - b) * (dickmanlocal (alpha - 1., x) + dickmanlocal_i (bi, x) / (alpha - b)) / 2.
  );
}

brentsuyama (B1, B2, N, nr) = 
{
  local (a, ai, i, alpha, beta);
  alpha = log (N) / log (B1);
  beta = log (B2) / log (B1);
  ai = floor ((alpha - beta) * invh);
  a = ai * h;
  return (
    h * ( 
      sum (i = 1, ai - 1, 
        dickmanlocal_i (i, N) / (alpha - i * h) * (1 - exp(-nr * B1 ^ (-alpha + i * h)))
      )
/* Between 0 and h, rho() is 1 everywhere except at 0, so we take it as 1 */
      + (1 - exp(-nr / B1 ^ alpha)) / 2.
      + dickmanlocal_i (ai, N) / (alpha - a) * (1 - exp(-nr * B1 ^ (-alpha + a))) / 2.
    )
    + (alpha - beta - a) * (
        dickmanlocal_i (ai, N) / (alpha - a) + dickmanlocal (alpha - beta, N) / beta
      ) / 2.
  );
}

/* Probability that the difference of two degree S Dickson polynomials, each
   evaluated at nr random but distinct points, includes a prime p > B2 so 
   that N / p is B1-smooth. The linear factors of Dickson_S(a)-Dickson_S(b) 
   are assumed < B2 ! */
brsudickson (B1, B2, N, nr, S) =
{
  local (i, n, phi);
  n = 0;
  phi = eulerphi (S) / 2;
  for (i = 1, S / 2, 
    if (gcd (i, S) == 1,
/* redundancy could be avoided by counting how often each gcd(i,S) value occurs */
      n = n + brentsuyama (B1, B2, N, nr * (gcd (i - 1, S) + gcd (i + 1, S) - 4) / 2)
    )
  );
  return (n / phi);
}

/* Same, but for S-th power */
brsupower (B1, B2, N, nr, S) =
{
  local (i, n, phi);
  n = 0;
  phi = eulerphi (S);
  for (i = 1, S, 
    if (gcd (i, S) == 1,
      n = n + brentsuyama (B1, B2, N, nr * (gcd (i - 1, S) - 2))
    )
  );
  return (n / phi);
}

/* The probability of ECM finding a factor near N with stage 1 parameter B1,
   stage 2 parameter B2, and evaluating nr random but distinct points in 
   stage 2 with a degree -S Dickson polynomial (if S < 0) or the 
   S-th power (S > 0) as the Brent-Suyama function */
ecmprob (B1, B2, N, nr, S) =
{
  local (alpha, beta, stage1, stage2, brsu, Nadj);
  Nadj = N / 23.4;
  alpha = log (Nadj) / log (B1);
  beta = log (B2) / log (B1);
  stage1 = dickmanlocal (alpha, Nadj);
  stage2 = 0;
  if (B2 > B1, stage2 = dickmanmu (alpha, beta, Nadj));
  brsu = 0;
  if (S < -1,
    brsu = brsudickson (B1, B2, Nadj, nr, -S * 2)
  );
  if (S > 1,
    brsu = brsupower (B1, B2, Nadj, nr, S * 2)
  );
/*  print ("ecmprob: stage 1: ", stage1, ", stage 2: ", stage2, 
         ", Brent-Suyama: ", brsu); */
  return (stage1 + stage2 + brsu)
}

/* pm1prob is incomplete! */

pm1prob (B1, B2, N, nr, S) =
{
  local (alpha, beta, stage1, stage2, brsu, divalpha);
  /* The "root properties" of a large prime minus 1 are 
     alpha = sum_{p in Primes} log(p)/(p^2-2p+1) ~= 1.2269688... */
  divalpha = exp(1.2269688056534700059656625687457626422689456478473);
  alpha = log (N / divalpha) / log (B1);
  beta = log (B2) / log (B1);
  stage1 = dickmanlocal (alpha, N / divalpha);
  stage2 = 0;
  if (B2 > B1, stage2 = dickmanmu (alpha, beta, N / divalpha));
  brsu = 0;
  if (S < -1,
    brsu = brsudickson (B1, B2, N / divalpha, nr, -S)
  );
  if (S > 1,
    brsu = brsupower (B1, B2, N / divalpha, nr, S)
  );
 print ("pm1prob: stage 1: ", stage1, ", stage 2: ", stage2, 
         ", Brent-Suyama: ", brsu);
  return (stage1 + stage2 + brsu)
}
