/* Arithmetic modulo Fermat numbers.

  Copyright 2004, 2005 Alexander Kruppa.

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
#include <stdlib.h> /* for abs if assertions enabled */
#include "ecm-impl.h"
#include "ecm-gmp.h"

#if HAVE_LIMITS_H
# include <limits.h>
#else
# ifndef UINT_MAX
#  define UINT_MAX (~(unsigned int) 0)
# endif
#endif

/*
#define DEBUG 1
#define CHECKSUM 1
*/

static mpz_t gt;
static int gt_inited = 0;
static int radix2 = 0;
unsigned int Fermat;

#define CACHESIZE 512U

/* a' <- a+b, b' <- a-b. */

#define ADDSUB_MOD(a, b) \
  mpz_sub (gt, a, b); \
  mpz_add (a, a, b);  \
  F_mod_gt (b, n);    \
  F_mod_1 (a, n);

__GMP_DECLSPEC mp_limb_t __gmpn_mod_34lsub1 (mp_limb_t*, mp_size_t);

/* compute remainder modulo 2^(mp_bits_per_limb*3/4)-1 */
#ifndef HAVE___GMPN_MOD_34LSUB1
mp_limb_t
__gmpn_mod_34lsub1 (mp_limb_t *src, mp_size_t size)
{
  mp_ptr tp;
  mp_limb_t r, d;

  ASSERT(BITS_PER_MP_LIMB % 4 == 0);
  tp = malloc (size * sizeof (mp_limb_t));
  MPN_COPY (tp, src, size);
  d = (mp_limb_t) 1 << (3 * (mp_bits_per_limb / 4)) - (mp_limb_t) 1;
  mpn_divmod_1 (&r, tp, size, d);
  free (tp);
  return r;
}
#endif

/* RS -> RS (mod 2^n+1). If input |RS| < 2^(2*n), result |RS| < 2^(n+1) */

static inline void 
F_mod_1 (mpz_t RS, unsigned int n)
{
  mp_size_t size;
  mp_limb_t v;
  
  size = mpz_size (RS);

  if ((unsigned int) size == n / GMP_NUMB_BITS + 1)
    {
      int sgn;
      sgn = mpz_sgn (RS);          /* Remember original sign */
      v = mpz_getlimbn (RS, n / GMP_NUMB_BITS);
      mpz_tdiv_r_2exp (RS, RS, n); /* Just a truncate. RS < 2^n. Can make
                                      RS zero and so change sgn(RS)! */
      if (sgn == -1)
          mpz_add_ui (RS, RS, v);
      else
          mpz_sub_ui (RS, RS, v);
    }
  else if ((unsigned int) size > n / GMP_NUMB_BITS + 1)
    {                              /* Assuming |RS| < 2^(2*n) */
      mpz_tdiv_q_2exp (gt, RS, n); /* |gt| < 2^n */
      mpz_tdiv_r_2exp (RS, RS, n); /* |RS| < 2^n */
      mpz_sub (RS, RS, gt);        /* |RS| < 2^(n+1) */
    }
}


/* R = gt (mod 2^n+1) */

static inline void 
F_mod_gt (mpz_t R, unsigned int n)
{
  mp_size_t size;
  mp_limb_t v;
  
  size = mpz_size (gt);

  ASSERT(R != gt);

  if ((unsigned int) size == n / GMP_NUMB_BITS + 1)
    {
      int sgn;
      sgn = mpz_sgn (gt);
      v = mpz_getlimbn (gt, n / GMP_NUMB_BITS);
      mpz_tdiv_r_2exp (gt, gt, n); /* Just a truncate */
      if (sgn == -1)
          mpz_add_ui (R, gt, v);
      else
          mpz_sub_ui (R, gt, v);
    }
  else if ((unsigned int) size > n / GMP_NUMB_BITS + 1)
    {
      mpz_tdiv_q_2exp (R, gt, n);
      mpz_tdiv_r_2exp (gt, gt, n); /* Just a truncate */
      mpz_sub (R, gt, R);
    }
  else 
    mpz_set (R, gt);
}


/* R = S1 * S2 (mod 2^n+1) where n is a power of 2 */
/* S1 == S2, S1 == R, S2 == R ok, but none may == gt */

static void 
F_mulmod (mpz_t R, mpz_t S1, mpz_t S2, unsigned int n)
{
  int n2 = n / __GMP_BITS_PER_MP_LIMB; /* type of _mp_size is int */
  F_mod_1 (S1, n);
  F_mod_1 (S2, n);
  while (mpz_size (S1) > (unsigned) n2)
    {
      outputf (OUTPUT_ERROR, 
               "Warning: S1 >= 2^%d after reduction, has %lu bits. "
               "Trying again\n", n, (unsigned long) mpz_sizeinbase (S1, 2));
      F_mod_1 (S1, n);
    }
  while (mpz_size (S2) > (unsigned) n2)
    {
      outputf (OUTPUT_ERROR, 
               "Warning: S2 >= 2^%d after reduction, has %lu bits. "
               "Trying again\n", n, (unsigned long) mpz_sizeinbase (S2, 2));
      F_mod_1 (S2, n);
    }
#if defined(HAVE_GWNUM) && !defined (TUNE)
  if (n >= GWTHRESHOLD)
    {
#ifdef DEBUG
      mpz_t t, t2;
      mpz_init (t);
      mpz_mul (gt, S1, S2);
      F_mod_gt (t, n);
      F_mod_1 (t, n);
#endif
      Fgwmul (R, S1, S2);
      F_mod_1 (R, n);
#ifdef DEBUG
      if (mpz_cmp (t, R) != 0)
        {
          /* Perhaps they are congruent, but of opposite sign */
          mpz_init (t2);
          mpz_sub (t2, t, R);
          F_mod_1 (t2, n);
          if (mpz_sgn (t2) != 0)
            {
              outputf (OUTPUT_ERROR, "F_mulmod: results of mpz_mul and "
                       "Fgwmul differ:\n%Zd\n%Zd\n", t, R);
              return;
            }
          mpz_clear (t2);
        }
      mpz_clear (t);
#endif
      return;
    }
#elif defined(HAVE___GMPN_MUL_FFT)
  if (n >= 32768)
    {
      unsigned long k;
      
      /* WARNING: _mpz_realloc does not keep the value!!! */
      _mpz_realloc (gt, n2 + 1);
      k = mpn_fft_best_k (n2, S1 == S2);
      mpn_mul_fft (PTR(gt), n2, PTR(S1), ABSIZ(S1), PTR(S2), ABSIZ(S2), k);
      MPN_NORMALIZE(PTR(gt), n2);
      SIZ(gt) = ((SIZ(S1) ^ SIZ(S2)) >= 0) ? n2 : -n2;
      F_mod_gt (R, n);
      return;
    }
#endif
  mpz_mul (gt, S1, S2);
  F_mod_gt (R, n);
  return;
}

/* R = S + sgn(S)*(2^e) */

static void
mpz_absadd_2exp (mpz_t RS, unsigned int e)
{
  mp_size_t siz, limb_idx, bit_idx;
  mp_limb_t cy;
  int sgn;
  
  limb_idx = e / GMP_NUMB_BITS;
  bit_idx = e % GMP_NUMB_BITS;
  siz = mpz_size (RS);
  sgn = (mpz_sgn (RS) >= 0) ? 1 : -1;
  
  if (limb_idx >= RS->_mp_alloc)
    /* WARNING: mpz_realloc2 does not keep the value!!! */
    mpz_realloc2 (RS, (limb_idx + 1) * GMP_NUMB_BITS);
  
  /* Now RS->_mp_alloc > limb_idx) */
  
  while (siz <= limb_idx)
    {
      RS->_mp_d[siz++] = 0;
      RS->_mp_size += sgn;
    }
  
  /* Now RS->_mp_alloc >= siz > limb_idx */
  
  cy = mpn_add_1 (RS->_mp_d + limb_idx, RS->_mp_d + limb_idx,
                  siz - limb_idx, ((mp_limb_t)1) << bit_idx);
  if (cy)
    {
      if (RS->_mp_alloc <= siz)
        /* WARNING: mpz_realloc2 does not keep the value!!! */
        mpz_realloc2 (RS, (siz + 1) * GMP_NUMB_BITS);

      RS->_mp_d[siz] = 1;
      RS->_mp_size += sgn;
    }
}

/* R = S / 2 (mod 2^n + 1). S == gt is ok */

static void 
F_divby2 (mpz_t R, mpz_t S, unsigned int n)
{
  int odd, sgn;
  
  odd = mpz_odd_p (S);
  sgn = mpz_sgn (S);
  mpz_tdiv_q_2exp (R, S, 1);
  
  if (odd)
    {
      /* We shifted out a set bit at the bottom. With negative wrap-around,
         that becomes -2^(n-1), so we add -2^(n-1) + 2^n+1 = 2^(n-1)+1.
         If |S| < 2^(n+1), |R| < 2^n + 2^(n-1) + 1 < 2^(n+1) for n > 1. */
      
      mpz_absadd_2exp (R, n - 1);
      if (sgn < 0)
        mpz_sub_ui (R, R, 1);
      else
        mpz_add_ui (R, R, 1);
    }
}


/* RS = RS / 3 (mod 2^n + 1). RS == gt is ok */

static void 
F_divby3_1 (mpz_t RS, unsigned int n)
{
  /* 2^2^m == 1 (mod 3) for m>0, thus F_m == 2 (mod 3) */
  int mod, sgn;
  
  sgn = mpz_sgn (RS);
  mod = __gmpn_mod_34lsub1 (RS->_mp_d, mpz_size (RS)) % 3;

  if (mod == 1)
    {
      /* Add F_m. If |RS| < 2^(n+1), |RS|+F_m < 3*2^n+1 */
      mpz_absadd_2exp (RS, n);
      if (sgn >= 0)
        mpz_add_ui (RS, RS, 1);
      else
        mpz_sub_ui (RS, RS, 1);
    }
  else if (mod == 2)
    {
      /* Add 2 * F_m.  If |RS| < 2^(n+1), |RS|+2*F_m < 4*2^n+2 */
      mpz_absadd_2exp (RS, n + 1);
      if (sgn >= 0)
        mpz_add_ui (RS, RS, 2);
      else
        mpz_sub_ui (RS, RS, 2);
    }

  mpz_divby3_1op (RS); /* |RS| < (4*2^n+2)/3 < 2^(n+1) */
}

static void 
F_divby5_1 (mpz_t RS, unsigned int n)
{
  /* 2^2^m == 1 (mod 5) for m>1, thus F_m == 2 (mod 5) */
  int mod, sgn;

  sgn = mpz_sgn (RS);
  mod = __gmpn_mod_34lsub1 (RS->_mp_d, mpz_size (RS)) % 5;

  if (mod == 1)
    {
      /* Add 2 * F_m == 4 (mod 5) */
      mpz_absadd_2exp (RS, n + 1);
      if (sgn == 1)
        mpz_add_ui (RS, RS, 2);
      else
        mpz_sub_ui (RS, RS, 2);
    }
  else if (mod == 2)
    {
      /* Add 4 * F_m == 3 (mod 5) */
      mpz_absadd_2exp (RS, n + 2);
      if (sgn == 1)
        mpz_add_ui (RS, RS, 4);
      else
        mpz_sub_ui (RS, RS, 4);
    }
  else if (mod == 3)
    {
      /* Add F_m == 3 (mod 5) */
      mpz_absadd_2exp (RS, n);
      if (sgn == 1)
        mpz_add_ui (RS, RS, 1);
      else
        mpz_sub_ui (RS, RS, 1);
    }
  else if (mod == 4)
    {
      /* Add 3 * F_m == 1 (mod 5) */
      mpz_absadd_2exp (RS, n);
      mpz_absadd_2exp (RS, n + 1);
      if (sgn == 1)
        mpz_add_ui (RS, RS, 3);
      else
        mpz_sub_ui (RS, RS, 3);
    }

  ASSERT(mpz_divisible_ui_p (RS, 5));
  mpz_divexact_ui (RS, RS, 5);
}


/* A 2^(m+2) length convolution is possible:
   (2^(3n/4) - 2^(n/4))^2 == 2 (mod 2^n+1) 
   so we have an element of order 2^(m+2) of simple enough form
   to use it as a root of unity the transform */

/* Multiply by sqrt(2)^e (mod F_m).  n = 2^m */
/* R = (S * sqrt(2)^e) % (2^n+1) */
/* R == S is ok, but neither must be == gt */
/* Assumes abs(e) < 4*n */

static void 
F_mul_sqrt2exp (mpz_t R, mpz_t S, int e, unsigned int n) 
{
  int chgsgn = 0, odd;

  ASSERT(S != gt);
  ASSERT(R != gt);
  ASSERT((unsigned) abs (e) < 4 * n);

  if (e < 0)
    e += 4 * n;
  /* 0 <= e < 4*n */
  if ((unsigned) e >= 2 * n)    /* sqrt(2)^(2*n) == -1 (mod F_m), so */
    {
      e -= 2 * n;               /* sqrt(2)^e == -sqrt(2)^(e-2*n) (mod F_m) */
      chgsgn = 1;
    }				/* Now e < 2*n */

#ifdef DEBUG_PERF
  if (e == 0)
    outputf (OUTPUT_ALWAYS, "F_mul_sqrt2exp: called for trivial case %s1\n", 
             chgsgn ? "-" : "");
#endif

  odd = e & 1;
  e >>= 1;

  if (odd)
    {
      /* Multiply by sqrt(2) == 2^(3n/4) - 2^(n/4) */
      /* S * (2^(3n/4) - 2^(n/4)) == 2^(n/4) * (S*2^(n/2) - S) */
      mpz_mul_2exp (gt, S, n / 2);
      mpz_sub (gt, gt, S);
      mpz_tdiv_q_2exp (R, gt, n / 4 * 3);
      mpz_tdiv_r_2exp (gt, gt, n / 4 * 3);
      mpz_mul_2exp (gt, gt, n / 4);
      mpz_sub (R, gt, R);
      
      if (e != 0)
        {
          mpz_tdiv_q_2exp (gt, R, n-e);
          mpz_tdiv_r_2exp (R, R, n-e);
          mpz_mul_2exp (R, R, e);
          mpz_sub (R, R, gt);
        }
    }
  else if (e != 0) 
    {
      /*  S     = a*2^(n-e) + b,   b < 2^(n-e)  */
      /*  S*2^e = a*2^n + b*2^e = b*2^e - a */
      /*  b*2^e < 2^(n-e)*2^e = 2^n */
      mpz_tdiv_q_2exp (gt, S, n - e); /* upper e bits (=a) into gt */
      mpz_tdiv_r_2exp (R, S, n - e);  /* lower n-e bits (=b) into R */
                                      /* This is simply a truncate if S == R */
      mpz_mul_2exp (R, R, e);         /* R < 2^n */
      mpz_sub (R, R, gt);
    } else 
      mpz_set (R, S);

  if (chgsgn) 
    mpz_neg (R, R);
}

/* Same, but input may be gt. Input and output must not be identical */
static void 
F_mul_sqrt2exp_2 (mpz_t R, mpz_t S, int e, unsigned int n)
{
  int chgsgn = 0, odd;

  ASSERT (S != R);
  ASSERT (R != gt);
  ASSERT ((unsigned) abs (e) < 4 * n);

  if (e < 0)
    e += 4 * n;
  if ((unsigned) e >= 2 * n)    /* sqrt(2)^(2*n) == -1 (mod F_m), so */
    {
      e -= 2 * n;               /* sqrt(2)^e == -sqrt(2)^(e-2*n) (mod F_m) */
      chgsgn = 1;
    }				/* Now e < 2*n */

#ifdef DEBUG_PERF
  if (e == 0)
    outputf (OUTPUT_ALWAYS, "F_mul_sqrt2exp_2: called for trivial case %s1\n",
	     chgsgn ? "-" : "");
#endif

  odd = e & 1;
  e >>= 1;

  if (odd != 0)
    {
      mpz_set (R, S); /* Neccessary?  n/32 mov*/
      mpz_mul_2exp (gt, S, n / 2); /* May overwrite S  n/32 mov */
      mpz_sub (gt, gt, R); /* n/32 sub*/

      mpz_tdiv_q_2exp (R, gt, n / 4 * 3); /* 3*(n/32)/4 mov */
      mpz_tdiv_r_2exp (gt, gt, n / 4 * 3); /* Just a truncate */
      mpz_mul_2exp (gt, gt, n / 4); /* 3*(n/32)/4 mov */
      mpz_sub (R, gt, R); /* (n/32)/4 sub, 3*(n/32)/4 mov */
      
      if (e != 0)
        {
          mpz_tdiv_q_2exp (gt, R, n - e);
          mpz_tdiv_r_2exp (R, R, n - e);
          mpz_mul_2exp (R, R, e);
          mpz_sub (R, R, gt);
        }
    } 
  else if (e != 0) 
    {
      mpz_tdiv_q_2exp (R, S, n - e); /* upper e bits into R */
      mpz_tdiv_r_2exp (gt, S, n - e); /* lower n-e bits into gt */
      mpz_mul_2exp (gt, gt, e);
      mpz_sub (R, gt, R);
    } else 
      mpz_set (R, S);

  if (chgsgn == -1) 
    mpz_neg (R, R);
}

#define A0s A[0]
#define A1s A[l << stride2]
#define A2s A[2 * l << stride2]
#define A3s A[3 * l << stride2]
#define A0is A[i << stride2]
#define A1is A[(i + l) << stride2]
#define A2is A[(i + 2 * l) << stride2]
#define A3is A[(i + 3 * l) << stride2]

/* Decimation-in-frequency FFT. Unscrambled input, scrambled output. */
/* Elements are (mod 2^n+1), l and n must be powers of 2, l must be <= 4*n. */
/* Performs forward transform */

static void 
F_fft_dif (mpz_t *A, int l, int stride2, int n) 
{
  int i, omega = (4 * n) / l, iomega;

  if (l <= 1)
    return;

  ASSERT((4 * n) % l == 0);

  if (l == 2)
    {
      ADDSUB_MOD(A[0], A[1<<stride2]);
      return;
    }

  if (!radix2)
    {
      l /= 4;

      mpz_sub (gt, A1s, A3s);            /* gt = a1 - a3 */
      mpz_add (A1s, A1s, A3s);           /* A1 = a1 + a3 */
      F_mul_sqrt2exp_2 (A3s, gt, n, n);  /* A3 = i * (a1 - a3) */
      
      mpz_sub (gt, A[0], A2s);           /* gt = a0 - a2 */
      mpz_add (A[0], A[0], A2s);         /* A0 = a0 + a2 */

      mpz_sub (A2s, A[0], A1s);          /* A2 = a0 - a1 + a2 - a3 */
      mpz_add (A[0], A[0], A1s);         /* A0 = a0 + a1 + a2 + a3 */
      mpz_add (A1s, gt, A3s);            /* A1 = a0 - a2 + i * (a1 - a3) */
      mpz_sub (A3s, gt, A3s);            /* A3 = a0 - a2 - i * (a1 - a3) */

      for (i = 1, iomega = omega; i < l; i++, iomega += omega)
        {
          mpz_sub (gt, A1is, A3is);
          mpz_add (A1is, A1is, A3is);
          F_mul_sqrt2exp_2 (A3is, gt, n, n);
          
          mpz_sub (gt, A0is, A2is);
          mpz_add (A0is, A0is, A2is);
          
          mpz_sub (A2is, A0is, A1is);
          mpz_add (A0is, A0is, A1is);
          mpz_add (A1is, gt, A3is);
          mpz_sub (A3is, gt, A3is);
          F_mul_sqrt2exp (A1is, A1is, iomega, n);
          F_mul_sqrt2exp (A2is, A2is, 2 * iomega, n);
          F_mul_sqrt2exp (A3is, A3is, 3 * iomega, n);
        }

      if (l > 1)
        {
          F_fft_dif (A, l, stride2, n);
          F_fft_dif (A + (l << stride2), l, stride2, n);
          F_fft_dif (A + (2 * l << stride2), l, stride2, n);
          F_fft_dif (A + (3 * l << stride2), l, stride2, n);
        }
      return;
    }

  l /= 2;

  ADDSUB_MOD(A[0], A1s);

  for (i = 1, iomega = omega; i < l; i++, iomega += omega) 
    {
      mpz_sub (gt, A0is, A1is);
      mpz_add (A0is, A0is, A1is);
      F_mul_sqrt2exp_2 (A1is, gt, iomega, n);
      F_mod_1 (A0is, n);
    }
  
  F_fft_dif (A, l, stride2, n);
  F_fft_dif (A + (l << stride2), l, stride2, n);
}


/* Decimation-in-time inverse FFT. Scrambled input, unscrambled output */
/* Does not perform divide-by-length. l, and n as in F_fft_dif() */

static void 
F_fft_dit (mpz_t *A, int l, int stride2, int n) 
{
  int i, omega = (4 * n) / l, iomega;
  
  if (l <= 1)
    return;

  ASSERT((4 * n) % l == 0);

  if (l == 2)
    {
      ADDSUB_MOD(A[0], A[1<<stride2]);
      return;
    }

  if (!radix2)
    {
      l /= 4;
      
      if (l > 1)
        {
          F_fft_dit (A, l, stride2, n);
          F_fft_dit (A + (l << stride2), l, stride2, n);
          F_fft_dit (A + (2 * l << stride2), l, stride2, n);
          F_fft_dit (A + (3 * l << stride2), l, stride2, n);
        }

      mpz_sub (gt, A3s, A1s);            /* gt = -(a1 - a3) */
      mpz_add (A1s, A1s, A3s);           /* A1 = a1 + a3 */
      F_mul_sqrt2exp_2 (A3s, gt, n, n);  /* A3 = i * -(a1 - a3) */
      
      mpz_sub (gt, A[0], A2s);           /* gt = a0 - a2 */
      mpz_add (A[0], A[0], A2s);         /* A0 = a0 + a2 */
      
      mpz_sub (A2s, A[0], A1s);          /* A2 = a0 - a1 + a2 - a3 */
      mpz_add (A[0], A[0], A1s);         /* A0 = a0 + a1 + a2 + a3 */
      mpz_add (A1s, gt, A3s);            /* A1 = a0 - a2 + i * -(a1 - a3) */
      mpz_sub (A3s, gt, A3s);            /* A3 = a0 - a2 - i * -(a1 - a3) */

      for (i = 1, iomega = omega; i < l; i++, iomega += omega)
        {
          /* Divide by omega^i. Since sqrt(2)^(4*n) == 1 (mod 2^n+1), 
             this is like multiplying by omega^(4*n-i) */
          F_mul_sqrt2exp (A1is, A1is, 4 * n - iomega, n);
          F_mul_sqrt2exp (A2is, A2is, 4 * n - 2 * iomega, n);
          F_mul_sqrt2exp (A3is, A3is, 4 * n - 3 * iomega, n);

          mpz_sub (gt, A3is, A1is);
          mpz_add (A1is, A1is, A3is);
          F_mul_sqrt2exp_2 (A3is, gt, n, n);

          mpz_sub (gt, A0is, A2is);
          mpz_add (A0is, A0is, A2is);
          
          mpz_sub (A2is, A0is, A1is);
          mpz_add (A0is, A0is, A1is);
          mpz_add (A1is, gt, A3is);
          mpz_sub (A3is, gt, A3is);
          if (1)
            {
              F_mod_1 (A0is, n);
              F_mod_1 (A1is, n);
              F_mod_1 (A2is, n);
              F_mod_1 (A3is, n);
            }
        }
      return;
    }

  l /= 2;

  F_fft_dit (A, l, stride2, n);
  F_fft_dit (A + (l << stride2), l, stride2, n);
  
  ADDSUB_MOD(A[0], A1s);
  
  for (i = 1, iomega = 4*n - omega; i < l; i++, iomega -= omega) 
    {
      F_mul_sqrt2exp (A1is, A1is, iomega, n);
      mpz_sub (gt, A0is, A1is);
      mpz_add (A0is, A0is, A1is);
      F_mod_gt (A1is, n);
      F_mod_1 (A0is, n);
    }
  
}


#define A0 A[i]
#define A1 A[l+i]
#define A2 A[2*l+i]
#define A3 A[3*l+i]
#define B0 B[i]
#define B1 B[l+i]
#define B2 B[2*l+i]
#define B3 B[3*l+i]
#define C0 C[i]
#define C1 C[l+i]
#define C2 C[2*l+i]
#define C3 C[3*l+i]
#define C4 C[4*l+i]
#define C5 C[5*l+i]
#define C6 C[6*l+i]
#define C7 C[7*l+i]
#define t0 t[i]
#define t1 t[l+i]
#define t2 t[2*l+i]
#define t3 t[3*l+i]
#define t4 t[4*l+i]
#define t5 t[5*l+i]


static unsigned int
F_toomcook4 (mpz_t *C, mpz_t *A, mpz_t *B, unsigned int len, unsigned int n, 
             mpz_t *t)
{
  unsigned int l, i, r;

  ASSERT(len % 4 == 0);

  l = len / 4;
  
  if (A == B) /* Squaring. The interpolation could probably be optimized, too */
    {
      for (i = 0; i < l; i++)
        {
          /*** Evaluate A(2), A(-2), 8*A(1/2) ***/
          mpz_mul_2exp (t0, A0, 1);
          mpz_add (t0, t0, A1);
          mpz_mul_2exp (t0, t0, 1);
          mpz_add (t0, t0, A2);
          mpz_mul_2exp (t0, t0, 1);
          mpz_add (t0, t0, A3);         /* t[0 .. l-1] = 8*A(1/2) < 15*N */
          F_mod_1 (t0, n);

          mpz_mul_2exp (t2, A3, 2);
          mpz_add (t2, t2, A1);
          mpz_mul_2exp (t2, t2, 1);     /* t[2l .. 3l-1] = 8*A_3 + 2*A_1 */
          mpz_mul_2exp (gt, A2, 2);
          mpz_add (gt, gt, A0);         /* gt = 4*A_2 + A0 */
          mpz_sub (t4, gt, t2);         /* t[4l .. 5l-1] = A(-2) */
          mpz_add (t2, t2, gt);         /* t[2l .. 3l-1] = A(2) */
          F_mod_1 (t4, n);
          F_mod_1 (t2, n);
          
          /* Evaluate A(1), A(-1) */
          mpz_add (C2, A0, A2);         /* May overwrite A2 */
          mpz_add (gt, A1, A3);
          mpz_sub (C4, C2, gt);         /* C4 = A(-1) */
          mpz_add (C2, C2, gt);         /* C2 = A(1) < 4*N */
          F_mod_1 (C2, n);
          F_mod_1 (C4, n);
        }

    /* A0  A1   A2   A3                     */
    /* A0      A(1)  A3  A(-1)              */
    /* C0  C1   C2   C3   C4    C5   C6  C7 */

      r = F_mul (t, t, t, l, DEFAULT, n, t + 6 * l);
        /* t0 = (8*A(1/2)) ^ 2 = 64*C(1/2) */
      r += F_mul (t + 2 * l, t + 2 * l, t + 2 * l, l, DEFAULT, n, t + 6 * l);
        /* t2 = A(2) ^ 2 = C(2) */
      r += F_mul (t + 4 * l, t + 4 * l, t + 4 * l, l, DEFAULT, n, t + 6 * l);
        /* t4 = A(-2) ^ 2 = C(-2) */
      r += F_mul (C, A, A, l, DEFAULT, n, t + 6 * l);
        /* C0 = A(0) ^ 2 = C(0) */
      r += F_mul (C + 6 * l, A + 3 * l, A + 3 * l, l, DEFAULT, n, t + 6 * l);
        /* C6 = A(inf) ^ 2 = C(inf) */
      r += F_mul (C + 2 * l, C + 2 * l, C + 2 * l, l, DEFAULT, n, t + 6 * l);
        /* C2 = A(1) ^ 2 = C(1). May overwrite A3 */
      r += F_mul (C + 4 * l, C + 4 * l, C + 4 * l, l, DEFAULT, n, t + 6 * l);
        /* C4 = A(-1) ^ 2 = C(-1) */
    }
  else /* Multiply */
    {
      for (i = 0; i < l; i++)
        {
          /*** Evaluate A(2), A(-2), 8*A(1/2) ***/
          mpz_mul_2exp (t0, A0, 1);
          mpz_add (t0, t0, A1);
          mpz_mul_2exp (t0, t0, 1);
          mpz_add (t0, t0, A2);
          mpz_mul_2exp (t0, t0, 1);
          mpz_add (t0, t0, A3);         /* t[0 .. l-1] = 8*A(1/2) < 15*N */
          F_mod_1 (t0, n);

          mpz_mul_2exp (t2, A3, 2);
          mpz_add (t2, t2, A1);
          mpz_mul_2exp (t2, t2, 1);     /* t[2l .. 3l-1] = 8*A_3 + 2*A_1 */

          mpz_mul_2exp (gt, A2, 2);
          mpz_add (gt, gt, A0);         /* gt = 4*A_2 + A0 */
          mpz_sub (t4, gt, t2);         /* t[4l .. 5l-1] = A(-2) */
          mpz_add (t2, t2, gt);         /* t[2l .. 3l-1] = A(2) */
          F_mod_1 (t4, n);
          F_mod_1 (t2, n);
          
          /*** Evaluate B(2), B(-2), 8*B(1/2) ***/
          mpz_mul_2exp (t1, B0, 1);
          mpz_add (t1, t1, B1);
          mpz_mul_2exp (t1, t1, 1);
          mpz_add (t1, t1, B2);
          mpz_mul_2exp (t1, t1, 1);
          mpz_add (t1, t1, B3);         /* t[l .. 2l-1] = 8*B(1/2) */
          F_mod_1 (t1, n);

          mpz_mul_2exp (t3, B3, 2);
          mpz_add (t3, t3, B1);
          mpz_mul_2exp (t3, t3, 1);     /* t[3l .. 4l-1] = 8*B_3 + 2*B_1 */

          mpz_mul_2exp (gt, B2, 2);
          mpz_add (gt, gt, B0);         /* gt = 4*B_2 + B0 */
          mpz_sub (t5, gt, t3);         /* t[5l .. 6l-1] = B(-2) */
          mpz_add (t3, t3, gt);         /* t[3l .. 4l-1] = B(2) */
          F_mod_1 (t5, n);
          F_mod_1 (t3, n);

          /* Evaluate A(1), A(-1) */
          mpz_add (C2, A0, A2);         /* May overwrite A2 */
#undef A2
          mpz_add (gt, A1, A3);
          mpz_set (C1, B0);             /* C1 = B(0) May overwrite A1 */
#undef A1
          mpz_sub (C4, C2, gt);         /* C4 = A(-1). May overwrite B0 */
#undef B0
          mpz_add (C2, C2, gt);         /* C2 = A(1) < 4*N */
          F_mod_1 (C2, n);
          F_mod_1 (C4, n);

          /* Evaluate B(1), B(-1) */
          mpz_add (gt, C1, B2);         /* B0 is in C1 */
          mpz_set (C6, A3);             /* C6 = A(inf) May overwrite B2 */
#undef B2
          mpz_add (C3, B1, B3);         /* May overwrite A3 */
#undef A3
          mpz_sub (C5, gt, C3);         /* C5 = B(-1). May overwrite B1 */
#undef B1
          mpz_add (C3, gt, C3);         /* C3 = B(1) */
          F_mod_1 (C3, n);
          F_mod_1 (C5, n);
        }

    /* A0 A1   A2   A3   B0    B1   B2 B3 */
    /* A0 B0  A(1) B(1) A(-1) B(-1) A3 B3 */
    /* C0 C1   C2   C3   C4    C5   C6 C7 */

      r = F_mul (t, t, t + l, l, DEFAULT, n, t + 6 * l);
        /* t0 = 8*A(1/2) * 8*B(1/2) = 64*C(1/2) */
      r += F_mul (t + 2 * l, t + 2 * l, t + 3 * l, l, DEFAULT, n, t + 6 * l);
        /* t2 = A(2) * B(2) = C(2) */
      r += F_mul (t + 4 * l, t + 4 * l, t + 5 * l, l, DEFAULT, n, t + 6 * l);
        /* t4 = A(-2) * B(-2) = C(-2) */
      r += F_mul (C, A, C + l, l, DEFAULT, n, t + 6 * l);
        /* C0 = A(0)*B(0) = C(0) */
      r += F_mul (C + 2 * l, C + 2 * l, C + 3 * l, l, DEFAULT, n, t + 6 * l);
        /* C2 = A(1)*B(1) = C(1) */
      r += F_mul (C + 4 * l, C + 4 * l, C + 5 * l, l, DEFAULT, n, t + 6 * l);
        /* C4 = A(-1)*B(-1) = C(-1) */
      r += F_mul (C + 6 * l, C + 6 * l, B + 3 * l, l, DEFAULT, n, t + 6 * l);
        /* C6 = A(inf)*B(inf) = C(inf) */
    }
  
/* C(0)   C(1)   C(-1)  C(inf)  64*C(1/2)  C(2)   C(-2) */
/* C0,C1  C2,C3  C4,C5  C6,C7   t0,t1      t2,t3  t4,t5 */

  for (i = 0; i < 2 * l - 1; i++)
    {
      mpz_add (t0, t0, t2);             /* t0 = 65 34 20 16 20 34 65 */

      mpz_sub (gt, C2, C4);             /* gt = 2*C_odd(1) = 0 2 0 2 0 2 0 */
      mpz_add (C2, C2, C4);             /* C2 = 2*C_even(1) = 2 0 2 0 2 0 2 */
      F_divby2 (C2, C2, n);             /* C2 = C_even(1) */

      mpz_add (C4, t2, t4);             /* C4 = 2*C_even(2) */
      F_divby2 (C4, C4, n);             /* C4 = C_even(2) */
      mpz_sub (t4, t2, t4);             /* t4 = 2*C_odd(2) */
      F_divby2 (t4, t4, n);
      F_divby2 (t4, t4, n);             /* t4 = C_odd(2)/2 = C_1 + 4*C_3 + 16*C_5 */
      F_divby2 (t2, gt, n);             /* t2 = C_odd(1) */

      mpz_sub (t0, t0, gt);             /* t0 = 65 32 20 14 20 32 65 */
      mpz_mul_2exp (gt, gt, 4);
      mpz_sub (t0, t0, gt);             /* t0 = 65 0 20 -18 20 0 65 */

      mpz_add (gt, C0, C6);             /* gt = C_0 + C_6 */
      mpz_sub (C2, C2, gt);             /* C2 = C_2 + C_4 */
      mpz_sub (t0, t0, gt);             /* t0 = 64 0 20 -18 20 0 64 */
      mpz_mul_2exp (gt, gt, 5);         /* gt = 32*C_0 + 32*C_6 */
      F_divby2 (t0, t0, n);             /* t0 = 32 0 10 -9 10 0 32 */
      mpz_sub (t0, t0, gt);             /* t0 = 0 0 10 -9 10 0 0 */
      mpz_sub (t0, t0, C2);             /* t0 = 0 0 9 -9 9 0 0 */
      F_divby3_1 (t0, n);
      F_divby3_1 (t0, n);               /* t0 = 0 0 1 -1 1 0 0 */
      mpz_sub (t0, C2, t0);             /* t0 = C_3 */
      mpz_sub (t2, t2, t0);             /* t2 = C_1 + C_5 */
      mpz_mul_2exp (gt, t0, 2);         /* gt = 4*C_3 */
      mpz_sub (t4, t4, gt);             /* t4 = C_1 + 16*C_5 */
      mpz_sub (t4, t4, t2);             /* t4 = 15*C_5 */
      F_divby3_1 (t4, n);
      F_divby5_1 (t4, n);               /* t4 = C_5 */
      mpz_sub (t2, t2, t4);             /* t2 = C_1 */

      mpz_sub (C4, C4, C0);             /* C4 = 4*C_2 + 16*C_4 + 64*C_6 */
      F_divby2 (C4, C4, n);
      F_divby2 (C4, C4, n);             /* C4 = C_2 + 4*C_4 + 16*C_6 */

      mpz_mul_2exp (gt, C6, 4);
      mpz_sub (C4, C4, gt);             /* C4 = C_2 + 4*C_4 */

      mpz_sub (C4, C4, C2);             /* C4 = 3*C_4 */
      F_divby3_1 (C4, n);               /* C4 = C_4 */
      mpz_sub (C2, C2, C4);             /* C2 = C_2 */
    }

  for (i = 0; i < l - 1; i++)
    {
      mpz_add (C1, C1, t2);
      F_mod_1 (C1, n);
    }
  mpz_set (C1, t2);
  F_mod_1 (C1, n);
  for (i = l; i < 2 * l - 1; i++)
    {
      mpz_add (C1, C1, t2);
      F_mod_1 (C1, n);
    }
  
  for (i = 0; i < l - 1; i++)
    {
      mpz_add (C3, C3, t0);
      F_mod_1 (C3, n);
    }
  mpz_set (C3, t0);
  F_mod_1 (C3, n);
  for (i = l; i < 2 * l - 1; i++)
    {
      mpz_add (C3, C3, t0);
      F_mod_1 (C3, n);
    }

  for (i = 0; i < l - 1; i++)
    {
      mpz_add (C5, C5, t4);
      F_mod_1 (C5, n);
    }
  mpz_set (C5, t4);
  F_mod_1 (C5, n);
  for (i = l; i < 2 * l - 1; i++)
    {
      mpz_add (C5, C5, t4);
      F_mod_1 (C5, n);
    }

  return r;
}


/* Karatsuba split. Calls F_mul() to multiply the three pieces. */

static unsigned int
F_karatsuba (mpz_t *R, mpz_t *A, mpz_t *B, unsigned int len, unsigned int n, 
             mpz_t *t)
{
  unsigned int i, r;

  ASSERT(len % 2 == 0);

  len /= 2;

  if (A == B) /* Squaring */
    {
      r = F_mul (t, A, A + len, len, DEFAULT, n, t + 2 * len); /* A0 * A1 */
      r += F_mul (R + 2 * len, A + len, A + len, len, DEFAULT, n, t + 2 * len); /* A1^2 */
      r += F_mul (R, A, A, len, DEFAULT, n, t + 2 * len); /* A0^2 */
      for (i = 0; i < 2 * len - 1; i++)
        {
          mpz_mul_2exp (t[i], t[i], 1);
          mpz_add (R[i + len], R[i + len], t[i]); /* i==len could be a mpz_set */
        }
      return r;
    }
  
  for (i = 0; i < len; i++) 
    {
      mpz_add (t[i],       A[i], A[i + len]); /* t0 = A0 + A1 */
      mpz_add (t[i + len], B[i], B[i + len]); /* t1 = B0 + B1 */
    }
  
  r = F_mul (t, t, t + len, len, DEFAULT, n, t + 2 * len);
  /* t[0...2*len-1] = (A0+A1) * (B0+B1) = A0*B0 + A0*B1 + A1*B0 + A1*B1 */
  
  if (R != A)
    {
      r += F_mul (R, A, B, len, DEFAULT, n, t + 2 * len);
      /* R[0...2*len-1] = A0 * B0 */
      r += F_mul (R + 2 * len, A + len, B + len, len, DEFAULT, n, t + 2 * len);
      /* R[2*len...4*len-1] = A1 * B1, may overwrite B */
    }
  else if (R + 2 * len != B)
    {
      r += F_mul (R + 2 * len, A + len, B + len, len, DEFAULT, n, t + 2 * len);
      /* R[2*len...4*len-1] = A1 * B1 */
      r += F_mul (R, A, B, len, DEFAULT, n, t + 2 * len);
      /* R[0...2*len-1] = A0 * B0, overwrites A */
    }
  else /* R == A && R + 2*len == B */
    {
      for (i = 0; i < len; i++)
        { /* mpz_swap instead? Perhaps undo later? Or interface for F_mul
             to specify separate result arrays for high/low half? */
          mpz_set (gt, A[len + i]); /* Swap A1 and B0 */
          mpz_set (A[len + i], B[i]);
          mpz_set (B[i], gt);
        }
      r += F_mul (R, R, R + len, len, DEFAULT, n, t + 2 * len);
      /* R[0...2*len-1] = A0 * B0, overwrites A */
      r += F_mul (R + 2 * len, R + 2 * len, R + 3 * len, len, DEFAULT, n, t + 2 * len);
      /* R[2*len...4*len-1] = A1 * B1, overwrites B */
    }

  /* R[0...2*len-2] == A0*B0, R[2*len-1] == 0 */
  /* R[2*len...3*len-2] == A1*B1, R[4*len-1] == 0 */
  /* t[0...2*len-2] == (A0+A1)*(B0+B1), t[2*len-1] == 0 */

  /* We're doing indices i and i+len in one loop on the assumption
     that 6 residues will probably fit into cache. After all,
     Karatsuba is only called for smallish F_m. This way, the final
     add R[i+len] += t[i] can be done inside the same loop and we need
     only one pass over main memory. */

  for (i = 0; i < len - 1; i++) 
    {
      mpz_sub (t[i], t[i], R[i]); /* t = A0*B1 + A1*B0 + A1*B1 */
      mpz_sub (t[i], t[i], R[i + 2 * len]); /* t = A0*B1 + A1*B0 */
      mpz_sub (t[i + len], t[i + len], R[i + len]);
      mpz_sub (t[i + len], t[i + len ], R[i + 3 * len]);
      
      mpz_add (R[i + len], R[i + len], t[i]);
      mpz_add (R[i + 2 * len], R[i + 2 * len], t[i + len]);
    }
  mpz_sub (t[len - 1], t[len - 1], R[len - 1]);
  mpz_sub (R[2 * len - 1], t[len - 1], R[3 * len - 1]);
  
  return r;
}

/* Multiply two polynomials with coefficients modulo 2^(2^m)+1. */
/* len is length (=degree+1) of polynomials and must be a power of 2. */
/* n=2^m */
/* Return value: number of multiplies performed, or UINT_MAX in case of error */

unsigned int 
F_mul (mpz_t *R, mpz_t *A, mpz_t *B, unsigned int len, int parameter, 
       unsigned int n, mpz_t *t)
{
  unsigned int i, r=0;
  unsigned int transformlen = (parameter == NOPAD) ? len : 2 * len;
#ifdef CHECKSUM
  mpz_t chksum1, chksum_1, chksum0, chksuminf;
#endif

  /* Handle trivial cases */
  if (len == 0)
    return 0;
  
  if (!gt_inited)
    {
      mpz_init2 (gt, 2 * n);
      gt_inited = 1;
    }
  
  if (len == 1)
    {
      if (parameter == MONIC) 
        {
          /* (x + a0)(x + b0) = x^2 + (a0 + b0)x + a0*b0 */
          mpz_add (gt, A[0], B[0]);
          F_mod_gt (t[0], n);
          F_mulmod (R[0], A[0], B[0], n); /* May overwrite A[0] */
          mpz_set (R[1], t[0]); /* May overwrite B[0] */
          /* We don't store the leading 1 monomial in the result poly */
        }
      else
        {
          F_mulmod (R[0], A[0], B[0], n); /* May overwrite A[0] */
          mpz_set_ui (R[1], 0); /* May overwrite B[0] */
        }
      
      return 1;
    }

#ifdef CHECKSUM
  mpz_init2 (chksum1, n+64);
  mpz_init2 (chksum_1, n+64);
  mpz_init2 (chksum0, n+64);
  mpz_init2 (chksuminf, n+64);

  mpz_set_ui (gt, 0);
  for (i = 0; i < len; i++) 
    {
      /* Compute A(1) and B(1) */
      mpz_add (chksum1, chksum1, A[i]);
      mpz_add (gt, gt, B[i]);

      /* Compute A(-1) and B(-1) */
      if (i % 2 == 0)
        {
          mpz_add (chksum_1, chksum_1, A[i]);
          mpz_add (chksum0, chksum0, B[i]); /* chksum0 used temporarily here */
        }
      else
        {
          mpz_sub (chksum_1, chksum_1, A[i]);
          mpz_sub (chksum0, chksum0, B[i]);
        }
    }

  if (parameter == MONIC)
    {
      mpz_add_ui (chksum1, chksum1, 1);
      mpz_add_ui (gt, gt, 1);
      mpz_add_ui (chksum_1, chksum_1, 1);
      mpz_add_ui (chksum0, chksum0, 1);
    }
  
  mpz_mul (gt, gt, chksum1);
  F_mod_gt (chksum1, n);

  mpz_mul (gt, chksum0, chksum_1);
  F_mod_gt (chksum_1, n);

  /* Compute A(0) * B(0) */
  mpz_mul (gt, A[0], B[0]);
  F_mod_gt (chksum0, n);

  /* Compute A(inf) * B(inf) */
  mpz_mul (gt, A[len - 1], B[len - 1]);
  F_mod_gt (chksuminf, n);
  if (parameter == MONIC)
    {
      mpz_add (chksuminf, chksuminf, A[len - 2]);
      mpz_add (chksuminf, chksuminf, B[len - 2]);
    }

  r += 4;
#endif /* CHECKSUM */

  /* Don't do FFT if len =< 4 (Karatsuba or Toom-Cook are faster) unless we 
     do a transform without zero padding, or if transformlen > 4*n 
     (no suitable primitive roots of 1) */
  if ((len > 4 || parameter == NOPAD) && transformlen <= 4 * n) 
    {
      unsigned int len2;
      
      /* len2 = log_2(transformlen). Assumes transformlen > 0 */
      for (i = transformlen, len2 = 0; (i&1) == 0; i >>= 1, len2++);
      
      if (i != 1) 
        {
          outputf (OUTPUT_ERROR, "F_mul: polynomial length must be power of 2, "
                           "but is %d\n", len);
          return UINT_MAX;
        }
      
      /* Are we performing a squaring or multiplication? */
      if (A != B) 
        {
          /* So it's a multiplication */
          
          /* Put transform of B into t */
          for (i = 0; i < len; i++)
            mpz_set (t[i], B[i]);
          if (parameter == MONIC)
            mpz_set_ui (t[i++], 1);
          for (; i < transformlen; i++)
            mpz_set_ui (t[i], 0);

          F_fft_dif (t, transformlen, 0, n);
        } else
          t = R; /* Do squaring */

      /* Put A into R */
      for (i = 0; i < len; i++) 
        mpz_set (R[i], A[i]);
      if (parameter == MONIC)
        mpz_set_ui (R[i++], 1); /* May overwrite B[0] */
      for (; i < transformlen; i++)
        mpz_set_ui (R[i], 0); /* May overwrite B[i - len] */

      F_fft_dif (R, transformlen, 0, n);

      for (i = 0; i < transformlen; i++) 
        {
          F_mulmod (R[i], R[i], t[i], n);
          /* Do the div-by-length. Transform length was transformlen, 
             len2 = log_2 (transformlen), so divide by 
             2^(len2) = sqrt(2)^(2*len2) */

          F_mul_sqrt2exp (R[i], R[i], - 2 * len2, n);
        }

      r += transformlen;

      F_fft_dit (R, transformlen, 0, n);

      if (parameter == MONIC)
        mpz_sub_ui (R[0], R[0], 1);
      
    } else { /* Karatsuba or Toom-Cook split */
      
      if (parameter == NOPAD)
        {
          outputf (OUTPUT_ERROR, "F_mul: cyclic/short products not supported "
                   "by Karatsuba/Toom-Cook\n");
          return UINT_MAX;
        }
      
      if (len / n == 4 || len == 2)
        r += F_karatsuba (R, A, B, len, n, t);
      else
        r += F_toomcook4 (R, A, B, len, n, t);

      if (parameter == MONIC) /* Handle the leading monomial the hard way */
        {
          /* This only works if A, B and R do not overlap */
          if (A == R || B == R + len)
            {
              outputf (OUTPUT_ERROR, "F_mul: monic polynomials with Karatsuba/"
                       "Toom-Cook and overlapping input/output not supported\n");
              return UINT_MAX;
            }
          for (i = 0; i < len; i++)
            {
              mpz_add (R[i + len], R[i + len], A[i]);
              mpz_add (R[i + len], R[i + len], B[i]);
              F_mod_1 (R[i + len], n);
            }
        }
    }
      
#ifdef DEBUG
  if (parameter != MONIC && parameter != NOPAD)
    {
      F_mod_1 (R[transformlen - 1], n);
      if (mpz_sgn (R[transformlen - 1]) != 0)
        outputf (OUTPUT_ALWAYS, "F_mul, len %d: R[%d] == %Zd != 0\n", 
		     len, transformlen - 1, R[transformlen - 1]);
    }
#endif

#ifdef CHECKSUM
  /* Compute R(1) = (A*B)(1) and subtract from chksum1 */

  for (i = 0; i < transformlen; i++) 
    mpz_sub (chksum1, chksum1, R[i]);
  
  if (parameter == MONIC)
    mpz_sub_ui (chksum1, chksum1, 1);

  while (mpz_sizeinbase (chksum1, 2) > n) 
    F_mod_1 (chksum1, n);

  if (mpz_sgn (chksum1) != 0) 
    outputf (OUTPUT_ALWAYS, "F_mul, len %d: A(1)*B(1) != R(1), difference %Zd\n", 
                 len, chksum1);

  /* Compute R(-1) = (A*B)(-1) and subtract from chksum_1 */

  for (i = 0; i < transformlen; i++) 
    if (i % 2 == 0)
      mpz_sub (chksum_1, chksum_1, R[i]);
    else
      mpz_add (chksum_1, chksum_1, R[i]);
  
  if (parameter == MONIC)
    mpz_sub_ui (chksum_1, chksum_1, 1);
  
  while (mpz_sizeinbase (chksum_1, 2) > n) 
    F_mod_1 (chksum_1, n);

  if (mpz_sgn (chksum_1) != 0) 
    outputf (OUTPUT_ALWAYS, "F_mul, len %d: A(-1)*B(-1) != R(-1), difference %Zd\n", 
		 len, chksum_1);

  if (parameter != NOPAD)
    {
      mpz_sub (chksum0, chksum0, R[0]);
      while (mpz_sizeinbase (chksum0, 2) > n) 
        F_mod_1 (chksum0, n);

      if (mpz_sgn (chksum0) != 0) 
        outputf (OUTPUT_ALWAYS, "F_mul, len %d: A(0)*B(0) != R(0), difference %Zd\n", 
                   len, chksum0);

      mpz_sub (chksuminf, chksuminf, R[transformlen - 2]);
      while (mpz_sizeinbase (chksuminf, 2) > n) 
        F_mod_1 (chksuminf, n);

      if (mpz_sgn (chksuminf) != 0) 
        outputf (OUTPUT_ALWAYS, "F_mul, len %d: A(inf)*B(inf) != R(inf), difference %Zd\n", 
                    len, chksuminf);
    }

  mpz_clear (chksum1);
  mpz_clear (chksum_1);
  mpz_clear (chksum0);
  mpz_clear (chksuminf);
#endif /* CHECKSUM */

  return r;
}

/* Transposed multiply of two polynomials with coefficients 
   modulo 2^(2^m)+1.
   len is length of polynomial B and must be a power of 2,
   the length of polynomial A is len / 2. 
   n=2^m
   t must have space for 2*len coefficients 
   Only the product coefficients [len / 2 - 1 ... len - 2] will go into 
   R[0 ... len / 2 - 1] 
   Return value: number of multiplies performed, UINT_MAX in error case. */

unsigned int 
F_mul_trans (mpz_t *R, mpz_t *A, mpz_t *B, unsigned int len, unsigned int n, 
             mpz_t *t)
{
  unsigned int i, r = 0, len2;

  /* Handle trivial cases */
  if (len < 2)
    return 0;
  
  if (!gt_inited)
    {
      mpz_init2 (gt, 2 * n);
      gt_inited = 1;
    }
  
  if (len == 2)
    {
      F_mulmod (R[0], A[0], B[0], n);
      return 1;
    }

  if (len <= 4 * n)
    {
      /* len2 = log_2(len) */
      for (i = len, len2 = 0; i > 1 && (i&1) == 0; i >>= 1, len2++);
      
      if (i != 1) 
        {
          outputf (OUTPUT_ERROR, "F_mul_trans: polynomial length must be power of 2, "
                           "but is %d\n", len);
          return UINT_MAX;
        }
      
      /* Put transform of B into t */
      for (i = 0; i < len; i++)
        mpz_set (t[i], B[i]);

      F_fft_dif (t, len, 0, n);

      /* Put transform of reversed A into t + len */
      for (i = 0; i < len / 2; i++) 
        mpz_set (t[i + len], A[len / 2 - 1 - i]);
      for (i = len / 2; i < len; i++)
        mpz_set_ui (t[i + len], 0);

      F_fft_dif (t + len, len, 0, n);

      for (i = 0; i < len; i++) 
        {
          F_mulmod (t[i], t[i], t[i + len], n);
          /* Do the div-by-length. Transform length was len, so divide by
             2^len2 = sqrt(2)^(2*len2) */
          F_mul_sqrt2exp (t[i], t[i], - 2 * len2, n);
        }

      r += len;

      F_fft_dit (t, len, 0, n);
      
      for (i = 0; i < len / 2; i++)
        mpz_set (R[i], t[i + len / 2 - 1]);

    } else { /* Only Karatsuba, no Toom-Cook here */
      unsigned int h = len / 4;
      
      /* A = a1 * x^h + a0
         B = b3 * x^3h + b2 * x^2h + b1 * x^h + b0
         mul^T(A, B) = mul^T(a0,b3) * x^4h + 
                      (mul^T(a1,b3) + mul^T(a0,b2)) * x^3h + 
                      (mul^T(a1,b2) + mul^T(a0,b1)) * x^2h + 
                      (mul^T(a1,b1) + mul^T(a0,b0)) * x + 
                       mul^T(a1,b0)
         We only want the x^h, x^2h and x^3h coefficients,
         mul^T(a1,b1) + mul^T(a0,b0)
         mul^T(a1,b2) + mul^T(a0,b1)
         mul^T(a1,b3) + mul^T(a0,b2)

         Specifically, we want R[i] = \sum_{j=0}^{len/2} A[j] * B[j+i], 0<=i<len/2
      */
      
      /* t */
      for (i = 0; i < h; i++)
        mpz_add (t[i], A[i], A[i + h]);
      F_mul_trans (t, t, B + h, 2 * h, n, t + h); /* t[0 ... h-1] = t */
                                                  /* Uses t[h ... 5h-1] as temp */

      /* t[i] = \sum_{j=0}^{h-1} (A[j]+A[j+h]) * B[j+i+h], 0 <= i < h */
      
      /* r */
      for (i = 0; i < 2 * h; i++)
        mpz_sub (t[i + h], B[i], B[h + i]);
      F_mul_trans (t + h, A, t + h, 2 * h, n, t + 3 * h); /* t[h ... 2h-1] = r */
                                                          /* Uses t[3h ... 7h-1] as temp */
      /* t[i + h] = \sum_{j=0}^{h-1} A[j] * (B[j+i] - B[j+i+h]), 0 <= i < h */
      
      for (i = 0; i < h; i++)
        mpz_add (R[i], t[i], t[i + h]); /* R[0 ... h-1] = t + r */
                                        /* r not needed anymore */
      
      /* R[i] = \sum_{j=0}^{h-1} (A[j]+A[j+h]) * B[j+i+h] +  \sum_{j=0}^{h-1} A[j] * (B[j+i] - B[j+i+h]) =
                \sum_{j=0}^{h-1} (A[j]+A[j+h]) * B[j+i+h] + A[j] * (B[j+i] - B[j+i+h]) = 
                \sum_{j=0}^{h-1} A[j]*B[j+i+h] + A[j+h]*B[j+i+h] + A[j]*B[j+i] - A[j]*B[j+i+h] =
                \sum_{j=0}^{h-1} A[j+h]*B[j+i+h] + A[j]*B[j+i] =
                \sum_{j=0}^{h-1} A[j]*B[j+i] + \sum_{j=0}^{h-1} A[j+h]*B[j+i+h] =
                \sum_{j=0}^{h-1} A[j]*B[j+i] + \sum_{j=h}^{2h-1} A[j]*B[j+i] =

                \sum_{j=0}^{2h-1} A[j]*B[j+i], 0 <= i < h   (1)
      */
      
      /* s */
      for (i = 0; i < 2 * h; i++)
        mpz_sub (t[i + h], B[i + 2 * h], B[i + h]);
      F_mul_trans (t + h, A + h, t + h, 2 * h, n, t + 3 * h); /* t[h ... 2h-1] = s */
                                                              /* Uses t[3h ... 7h - 1] as temp */

      /* t[i + h] = \sum_{j=0}^{h} A[j+h] * (B[j+i+2*h]-B[j+i+h]), 0 <= i < h */
      
      for (i = 0; i < h; i++)
        mpz_add (R[i + h], t[i], t[i + h]);
      
      /* R[i + h] = \sum_{j=0}^{h-1} (A[j]+A[j+h]) * B[j+i+h] + \sum_{j=0}^{h} A[j+h] * (B[j+i+2*h]-B[j+i+h]) =
                    \sum_{j=0}^{h-1} (A[j]+A[j+h]) * B[j+i+h] + A[j+h] * (B[j+i+2*h]-B[j+i+h]) =
                    \sum_{j=0}^{h-1} A[j]*B[j+i+h] + A[j+h]*B[j+i+h] + A[j+h]*B[j+i+2*h] - A[j+h]*B[j+i+h] =
                    \sum_{j=0}^{h-1} A[j]*B[j+i+h] + A[j+h]*B[j+i+2*h] =
                    \sum_{j=0}^{h-1} A[j]*B[j+i+h] + \sum_{j=0}^{h-1} A[j+h]*B[j+i+2*h] =
                    \sum_{j=0}^{h-1} A[j]*B[j+i+h] + \sum_{j=h}^{2h-1} A[j]*B[j+i+h] =
                    \sum_{j=0}^{2h-1} A[j]*B[j+i+h], 0 <= i < h

        R[i] = \sum_{j=0}^{2h-1} A[j]*B[j+i], h <= i < 2*h   (2)
        
        (1) and (2) : R[i] = \sum_{j=0}^{2h-1} A[j]*B[j+i], 0 <= i < 2*h
      */
          
    }
  
  return r;
}

void F_clear ()
{
  if (gt_inited)
    mpz_clear (gt);
  gt_inited = 0;
}
