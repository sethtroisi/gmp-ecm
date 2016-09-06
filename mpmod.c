/* Modular multiplication.

Copyright 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012
Paul Zimmermann, Alexander Kruppa and Cyril Bouvier.

This file is part of the ECM Library.

The ECM Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The ECM Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the ECM Library; see the file COPYING.LIB.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#include <stdio.h>
#include <stdlib.h>
#include "ecm-gmp.h"
#include "ecm-impl.h"
#include "mpmod.h"

#ifdef USE_ASM_REDC
  #include "mulredc.h"
#endif

FILE *ECM_STDOUT, *ECM_STDERR; /* define them here since needed in tune.c */

/* define WANT_ASSERT to check normalization of residues */
/* #define WANT_ASSERT 1 */
/* #define DEBUG */
/* #define WANT_ASSERT_EXPENSIVE 1 */

#define ASSERT_NORMALIZED(x) ASSERT ((modulus->repr != ECM_MOD_MODMULN && \
				      modulus->repr != ECM_MOD_REDC) || \
			     mpz_size (x) <= mpz_size (modulus->orig_modulus))
#define MPZ_NORMALIZED(x)    ASSERT (PTR(x)[ABSIZ(x)-1] != 0)



static void ecm_redc_basecase (mpz_ptr, mpz_ptr, mpmod_t) ATTRIBUTE_HOT;
static void ecm_mulredc_basecase (mpres_t, const mpres_t, const mpres_t, 
                                  mpmod_t) ATTRIBUTE_HOT;
static void base2mod (mpres_t, const mpres_t, mpres_t, mpmod_t) ATTRIBUTE_HOT;
static void REDC (mpres_t, const mpres_t, mpz_t, mpmod_t);

/* returns +/-l if n is a factor of N = 2^l +/- 1 with N <= n^threshold, 
   0 otherwise.
*/
int 
isbase2 (const mpz_t n, const double threshold)
{
  unsigned int k, lo; 
  int res = 0; 
  mpz_t u, w;

  mpz_init (u);
  mpz_init (w);
  lo = mpz_sizeinbase (n, 2) - 1; /* 2^lo <= n < 2^(lo+1) */  
  mpz_set_ui (u, 1UL);
  mpz_mul_2exp (u, u, 2UL * lo);
  mpz_mod (w, u, n); /* 2^(2lo) mod n = -/+2^(2lo-l) if m*n = 2^l+/-1 */
  if (mpz_cmp_ui (w, 1UL) == 0) /* if 2^(2lo) mod n = 1, then n divides
                                 2^(2lo)-1. If algebraic factors have been
                                 removed, n divides either 2^lo+1 or 2^lo-1.
                                 But since n has lo+1 bits, n can only divide
                                 2^lo+1. More precisely, n must be 2^lo+1. */
    {
      /* check that n equals 2^lo+1. Since n divides 2^(2lo)-1, n is odd. */
      if (mpz_scan1 (n, 1UL) != lo)
        lo = 0;
      mpz_clear (w);
      mpz_clear (u);
      return lo;
    }
  k = mpz_sizeinbase (w, 2) - 1;
  /* if w = 2^k then n divides 2^(2*lo-k)-1 */
  mpz_set_ui (u, 1UL);
  mpz_mul_2exp (u, u, k);
  if (mpz_cmp(w, u) == 0) 
    res = k - 2 * lo;
  else /* if w = -2^k then n divides 2^(2*lo-k)+1 */
    {
      mpz_neg (w, w);
      mpz_mod (w, w, n);
      k = mpz_sizeinbase (w, 2) - 1;
      mpz_set_ui (u, 1UL);
      mpz_mul_2exp (u, u, k);
      if (mpz_cmp (w, u) == 0) 
        res = 2 * lo - k;
    }
  mpz_clear (u);
  mpz_clear (w);

  if (abs (res) > (int) (threshold * (double) lo)) 
    res = 0;

  if (abs (res) < 16)
    res = 0;

  return res;
}

/* Do base-2 reduction. R must not equal S or t. */
static void
base2mod (mpres_t R, const mpres_t S, mpres_t t, mpmod_t modulus)
{
  unsigned long absbits = abs (modulus->bits);

  ASSERT (R != S && R != t);
  mpz_tdiv_q_2exp (R, S, absbits);
  mpz_tdiv_r_2exp (t, S, absbits);
  if (modulus->bits < 0)
    mpz_add (R, R, t);
  else
    mpz_sub (R, t, R);

  /* mpz_mod (R, R, modulus->orig_modulus); */
  while (mpz_sizeinbase (R, 2) > absbits)
    {
      mpz_tdiv_q_2exp (t, R, absbits);
      mpz_tdiv_r_2exp (R, R, absbits);
      if (modulus->bits < 0)
        mpz_add (R, R, t);
      else
        mpz_sub (R, R, t);
    }
}

/* Modular reduction modulo the Fermat number 2^m+1. 
   n = m / GMP_NUMB_BITS. Result is < 2^m+1.
   FIXME: this does not work with nails.
   Only copies the data to R if reduction is needed and returns 1 in that 
   case. If the value in S is reduced already, nothing is done and 0 is 
   returned. Yes, this is ugly. */
static int
base2mod_2 (mpres_t R, const mpres_t S, mp_size_t n, mpz_t modulus)
{
  mp_size_t s;

  s = ABSIZ(S);
  if (s > n)
    {
      if (s == n + 1)
        {
          mp_srcptr sp = PTR(S);
          mp_ptr rp;
          
          MPZ_REALLOC (R, s);
          rp = PTR(R);

          if ((rp[n] = mpn_sub_1 (rp, sp, n, sp[n])))
            rp[n] = mpn_add_1 (rp, rp, n, rp[n]);
          MPN_NORMALIZE(rp, s);
          ASSERT (s <= n || (s == n && rp[n] == 1));
          SIZ(R) = (SIZ(S) > 0) ? (int) s : (int) -s;
        }
      else /* should happen rarely */
        mpz_mod (R, S, modulus);

      return 1;
    }
  
  return 0;
}

/* subquadratic REDC, at mpn level.
   {orig,n} is the original modulus.
   Requires xn = 2n or 2n-1 and ABSIZ(orig_modulus)=n.
 */
static void
ecm_redc_n (mp_ptr rp, mp_srcptr x0p, mp_size_t xn,
            mp_srcptr orig, mp_srcptr invm, mp_size_t n)
{
  mp_ptr tp, up, xp;
  mp_size_t nn = n + n;
  mp_limb_t cy, cin;
  TMP_DECL;

  ASSERT((xn == nn) || (xn == nn - 1));

  TMP_MARK;
  up = TMP_ALLOC_LIMBS(nn + nn);
  if (xn < nn)
    {
      xp = TMP_ALLOC_LIMBS(nn);
      MPN_COPY (xp, x0p, xn);
      xp[nn - 1] = 0;
    }
  else
    xp = (mp_ptr) x0p;
#ifdef HAVE___GMPN_MULLO_N /* available up from GMP 5.0.0 */
  __gmpn_mullo_n (up, xp, invm, n);
#else
  ecm_mul_lo_n (up, xp, invm, n);
#endif
  tp = up + nn;
  mpn_mul_n (tp, up, orig, n);
  /* add {x, 2n} and {tp, 2n}. We know that {tp, n} + {xp, n} will give
     either 0, or a carry out. If xp[n-1] <> 0 or tp[n-1] <> 0, 
     then there is a carry. We use a binary OR, which sets the zero flag
     if and only if both operands are zero. */
  cin = (mp_limb_t) ((xp[n - 1] | tp[n - 1]) ? 1 : 0);
#ifdef HAVE___GMPN_ADD_NC
  cy = __gmpn_add_nc (rp, tp + n, xp + n, n, cin);
#else
  cy = mpn_add_n (rp, tp + n, xp + n, n);
  cy += mpn_add_1 (rp, rp, n, cin);
#endif
  /* since we add at most N-1 to the upper half of {x0p,2n},
     one adjustment is enough */
  if (cy)
    cy -= mpn_sub_n (rp, rp, orig, n);
  ASSERT (cy == 0);
  TMP_FREE;
}

/* REDC. x and t must not be identical, t has limb growth */
/* subquadratic REDC, at mpz level */
static void 
REDC (mpres_t r, const mpres_t x, mpz_t t, mpmod_t modulus)
{
  mp_size_t n = modulus->bits / GMP_NUMB_BITS;
  mp_size_t xn = ABSIZ(x);

  ASSERT (xn <= 2 * n);
  if (xn == 2 * n) /* ecm_redc_n also accepts xn=2n-1, but this seems slower
                    for now (see remark in TODO) */
    {
      mp_ptr rp;
      MPZ_REALLOC (r, n);
      rp = PTR(r);
      ecm_redc_n (rp, PTR(x), xn, PTR(modulus->orig_modulus),
                  PTR(modulus->aux_modulus), n);
      MPN_NORMALIZE(rp, n);
      SIZ(r) = (SIZ(x) > 0) ? (int) n : (int) -n;
      MPZ_NORMALIZED (r);
    }
  else
    {
      mpz_tdiv_r_2exp (t, x, modulus->bits);
      mpz_mul (t, t, modulus->aux_modulus);
      mpz_tdiv_r_2exp (t, t, modulus->bits);  /* t = (x % R) * 1/N (mod R) */
      mpz_mul (t, t, modulus->orig_modulus);  /* t <= (R-1)*N */
      mpz_add (t, t, x);                      /* t < (R-1)*N + R^2/B */
      mpz_tdiv_q_2exp (r, t, modulus->bits);  /* r = (x + m*N) / R
                                                 < N + R/B */
      if (ABSIZ (r) > n) /* this can happen in practice but is extremely
                            unlikely: in particular it requires that the
                            upper limb of N has only ones */
        mpz_sub (r, r, modulus->multiple);
    }
  ASSERT (ABSIZ(r) <= n);
}


/* Quadratic time redc for n word moduli. */
static inline void 
redc_basecase_n (mp_ptr rp, mp_ptr cp, mp_srcptr np, const mp_size_t nn, 
                 const mp_ptr invm)
{
#ifdef HAVE___GMPN_REDC_2
  REDC2(rp, cp, np, nn, invm);
#else /* HAVE___GMPN_REDC_2 is not defined */
#ifdef HAVE___GMPN_REDC_1
  REDC1(rp, cp, np, nn, invm[0]);
#else /* neither HAVE___GMPN_REDC_2 nor HAVE___GMPN_REDC_1 is defined */
  mp_limb_t cy;
  mp_size_t j;
  
  for (j = 0; j < nn; j++)
    {
      cy = mpn_addmul_1 (cp, np, nn, cp[0] * invm[0]);
      ASSERT(cp[0] == (mp_limb_t) 0);
      cp[0] = cy;
      cp++;
    }
  /* add vector of carries and shift */
  cy = mpn_add_n (rp, cp, cp - nn, nn);
  /* the result of Montgomery's REDC is less than 2^Nbits + N,
     thus at most one correction is enough */
  if (cy != 0)
    {
      mp_limb_t t;
      t = mpn_sub_n (rp, rp, np, nn); /* a borrow should always occur here */
      ASSERT (t == 1);
    }
#endif /* HAVE___GMPN_REDC_1 */
#endif /* HAVE___GMPN_REDC_2 */
}

/* r <- c/R^nn mod n, where n has nn limbs, and R=2^GMP_NUMB_BITS.
   n must be odd.
   c must have space for at least 2*nn limbs.
   r must have space for at least n limbs.
   c and r can be the same variable.
   The data in c is clobbered.
*/
static void 
ecm_redc_basecase (mpz_ptr r, mpz_ptr c, mpmod_t modulus)
{
  mp_ptr rp;
  mp_ptr cp;
  mp_srcptr np;
  mp_size_t j, nn = modulus->bits / GMP_NUMB_BITS;

  ASSERT(ABSIZ(c) <= 2 * nn);
  ASSERT(ALLOC(c) >= 2 * nn);
  ASSERT(ALLOC(r) >= nn);
  cp = PTR(c);
  rp = PTR(r);
  np = PTR(modulus->orig_modulus);
  for (j = ABSIZ(c); j < 2 * nn; j++) 
    cp[j] = 0;
  
  redc_basecase_n (rp, cp, np, nn, modulus->Nprim);

  MPN_NORMALIZE (rp, nn);
  SIZ(r) = SIZ(c) < 0 ? (int) -nn : (int) nn;
}

#ifdef USE_ASM_REDC
/* Quadratic time multiplication and REDC with nn-limb modulus.
   x and y are nn-limb residues, the nn-limb result is written to z. 
   This function merely calls the correct mulredc*() assembly function
   depending on nn, and processes any leftover carry. */

static void
mulredc (mp_ptr z, mp_srcptr x, mp_srcptr y, mp_srcptr m, 
         const mp_size_t nn, const mp_limb_t invm)
{
  mp_limb_t cy;
  switch (nn) 
    {
      case 1:
        cy = mulredc1(z, x[0], y[0], m[0], invm);
        break;
      case 2:
        cy = mulredc2(z, x, y, m, invm);
        break;
      case 3:
        cy = mulredc3(z, x, y, m, invm);
        break;
      case 4:
        cy = mulredc4(z, x, y, m, invm);
        break;
      case 5: 
        cy = mulredc5(z, x, y, m, invm);
        break;
      case 6: 
        cy = mulredc6(z, x, y, m, invm);
        break;
      case 7: 
        cy = mulredc7(z, x, y, m, invm);
        break;
      case 8:
        cy = mulredc8(z, x, y, m, invm);
        break;
      case 9:
        cy = mulredc9(z, x, y, m, invm);
        break;
      case 10:
        cy = mulredc10(z, x, y, m, invm);
        break;
      case 11:
        cy = mulredc11(z, x, y, m, invm);
        break;
      case 12:
        cy = mulredc12(z, x, y, m, invm);
        break;
      case 13:
        cy = mulredc13(z, x, y, m, invm);
        break;
      case 14:
        cy = mulredc14(z, x, y, m, invm);
        break;
      case 15:
        cy = mulredc15(z, x, y, m, invm);
        break;
      case 16:
        cy = mulredc16(z, x, y, m, invm);
        break;
      case 17:
        cy = mulredc17(z, x, y, m, invm);
        break;
      case 18:
        cy = mulredc18(z, x, y, m, invm);
        break;
      case 19:
        cy = mulredc19(z, x, y, m, invm);
        break;
      case 20:
        cy = mulredc20(z, x, y, m, invm);
        break;
      default:
        abort();
    }
  /* the result of Montgomery's REDC is less than 2^Nbits + N,
     thus at most one correction is enough */
  if (cy != 0)
    {
      ATTRIBUTE_UNUSED mp_limb_t t;
      t = mpn_sub_n (z, z, m, nn); /* a borrow should always occur here */
      ASSERT (t == 1);
    }
}

#if 0
/* {rp, n} <- {ap, n}^2/B^n mod {np, n} where B = 2^GMP_NUMB_BITS */
static void
sqrredc (mp_ptr rp, mp_srcptr ap, mp_srcptr np, const mp_size_t n,
         const mp_limb_t invm)
{
  mp_ptr cp;
  mp_size_t i;
  mp_limb_t cy, q;
  TMP_DECL;

  TMP_MARK;
  cp = TMP_ALLOC_LIMBS(2*n);
  for (i = 0; i < n; i++)
    umul_ppmm (cp[2*i+1], cp[2*i], ap[i], ap[i]);

  if (UNLIKELY(n == 1))
    {
      q = cp[0] * invm;
      rp[0] = mpn_addmul_1 (cp, np, 1, q);
      cy = mpn_add_n (rp, rp, cp + 1, 1);
      goto end_sqrredc;
    }

  if (cp[0] & (mp_limb_t) 1)
    /* cp[n] is either some ap[i]^2 mod B or floor(ap[i]^2/B),
       the latter is at most floor((B-1)^2/B) = B-2, and the former cannot be
       B-1 since -1 is not a square mod 2^n for n >1, thus there is no carry
       in cp[n] + ... below */
    cp[n] += mpn_add_n (cp, cp, np, n);
  /* now {cp, 2n} is even: divide by two */
  mpn_rshift (cp, cp, 2*n, 1);
  /* now cp[2n-1] is at most B/2-1 */

  for (i = 0; i < n - 1; i++)
    {
      q = cp[i] * invm;
      cp[i] = mpn_addmul_1 (cp + i, np, n, q);
      /* accumulate ap[i+1..n-1] * ap[i] */
      rp[i] = mpn_addmul_1 (cp + 2 * i + 1, ap + i + 1, n - 1 - i, ap[i]);
    }
  /* the last iteration did set cp[n-2] to zero, accumulated a[n-1] * a[n-2] */

  /* cp[2n-1] was untouched so far, so it is still at most B/2-1 */
  q = cp[n-1] * invm;
  rp[n-1] = mpn_addmul_1 (cp + n - 1, np, n, q);
  /* rp[n-1] <= floor((B^n-1)*(B-1)/B^n)<=B-2 */

  /* now add {rp, n}, {cp+n, n} and {cp, n-1} */
  /* cp[2n-1] still <= B/2-1 */
  rp[n-1] += mpn_add_n (rp, rp, cp, n-1); /* no overflow in rp[n-1] + ... */
  cy = mpn_add_n (rp, rp, cp + n, n);
  /* multiply by 2 */
  cy = (cy << 1) + mpn_lshift (rp, rp, n, 1);
 end_sqrredc:
  while (cy)
    cy -= mpn_sub_n (rp, rp, np, n);
  TMP_FREE;
}
#endif

#ifdef HAVE_NATIVE_MULREDC1_N
/* Multiplies y by the 1-limb value of x and does modulo reduction.
   The resulting residue may be multiplied by some constant, 
   which makes this function useful only for cases where, e.g.,
   all projective coordinates are multiplied by the same constant.
   More precisely it computes:
   {z, N} = {y, N} * x / 2^GMP_NUMB_BITS mod {m, N}
*/
static void
mulredc_1 (mp_ptr z, const mp_limb_t x, mp_srcptr y, mp_srcptr m, 
           const mp_size_t N, const mp_limb_t invm)
{
  mp_limb_t cy;

  switch (N) {
   case 1:
    cy = mulredc1(z, x, y[0], m[0], invm);
    break;
   case 2:
    cy = mulredc1_2(z, x, y, m, invm);
    break;
   case 3:
    cy = mulredc1_3(z, x, y, m, invm);
    break;
   case 4:
    cy = mulredc1_4(z, x, y, m, invm);
    break;
   case 5: 
    cy = mulredc1_5(z, x, y, m, invm);
    break;
   case 6: 
    cy = mulredc1_6(z, x, y, m, invm);
    break;
   case 7: 
    cy = mulredc1_7(z, x, y, m, invm);
    break;
   case 8:
    cy = mulredc1_8(z, x, y, m, invm);
    break;
   case 9:
    cy = mulredc1_9(z, x, y, m, invm);
    break;
   case 10:
    cy = mulredc1_10(z, x, y, m, invm);
    break;
   case 11:
    cy = mulredc1_11(z, x, y, m, invm);
    break;
   case 12:
    cy = mulredc1_12(z, x, y, m, invm);
    break;
   case 13:
    cy = mulredc1_13(z, x, y, m, invm);
    break;
   case 14:
    cy = mulredc1_14(z, x, y, m, invm);
    break;
   case 15:
    cy = mulredc1_15(z, x, y, m, invm);
    break;
   case 16:
    cy = mulredc1_16(z, x, y, m, invm);
    break;
   case 17:
    cy = mulredc1_17(z, x, y, m, invm);
    break;
   case 18:
    cy = mulredc1_18(z, x, y, m, invm);
    break;
   case 19:
    cy = mulredc1_19(z, x, y, m, invm);
    break;
   case 20:
    cy = mulredc1_20(z, x, y, m, invm);
    break;
   default:
    {
      abort ();
    }
  }
  /* the result of Montgomery's REDC is less than 2^Nbits + N,
     thus one correction (at most) is enough */
  if (cy != 0)
    {
      ATTRIBUTE_UNUSED mp_limb_t t;
      t = mpn_sub_n (z, z, m, N); /* a borrow should always occur here */
      ASSERT (t == 1);
    }
}
#endif /* ifdef HAVE_NATIVE_MULREDC1_N */
#endif

static int tune_mulredc_table[] = TUNE_MULREDC_TABLE;
static int tune_sqrredc_table[] = TUNE_SQRREDC_TABLE;

static void 
ecm_mulredc_basecase_n (mp_ptr rp, mp_srcptr s1p, mp_srcptr s2p, 
                        mp_srcptr np, mp_size_t nn, mp_ptr invm, mp_ptr tmp)
{
  mp_limb_t cy;
  mp_size_t j;

  if (nn <= MULREDC_ASSEMBLY_MAX)
    {
      switch (tune_mulredc_table[nn])
        {
        case MPMOD_MULREDC: /* use quadratic assembly mulredc */
#ifdef USE_ASM_REDC
          mulredc (rp, s1p, s2p, np, nn, invm[0]);
          break;
#endif /* otherwise go through to the next available mode */
        case MPMOD_MUL_REDC1: /* mpn_mul_n + __gmpn_redc_1 */
#ifdef HAVE___GMPN_REDC_1
          mpn_mul_n (tmp, s1p, s2p, nn);
          REDC1(rp, tmp, np, nn, invm[0]);
          break;
#endif /* otherwise go through to the next available mode */
        case MPMOD_MUL_REDC2: /* mpn_mul_n + __gmpn_redc_2 */
#ifdef HAVE___GMPN_REDC_2
          mpn_mul_n (tmp, s1p, s2p, nn);
          REDC2(rp, tmp, np, nn, invm);
          break;
#endif /* otherwise go through to the next available mode */
        case MPMOD_MUL_REDCN: /* mpn_mul_n + __gmpn_redc_n */
        /* disable redc_n for now, since it uses the opposite
           precomputed inverse wrt redc_1 and redc_2 */
#ifdef HAVE___GMPN_REDC_Nxxx
          mpn_mul_n (tmp, s1p, s2p, nn);
          __gmpn_redc_n (rp, tmp, np, nn, invm);
          break;
#endif /* otherwise go through to the next available mode */
        case MPMOD_MUL_REDC_C: /* plain C quadratic reduction */
          mpn_mul_n (tmp, s1p, s2p, nn);
          for (j = 0; j < nn; j++, tmp++)
            tmp[0] = mpn_addmul_1 (tmp, np, nn, tmp[0] * invm[0]);
          cy = mpn_add_n (rp, tmp - nn, tmp, nn);
          if (cy != 0)
            mpn_sub_n (rp, rp, np, nn); /* a borrow should always occur here */
          break;
        default:
          {
            outputf (OUTPUT_ERROR, "Invalid mulredc mode: %d\n",
                     tune_mulredc_table[nn]);
            exit (EXIT_FAILURE);
          }
        }
    }
  else /* nn > MULREDC_ASSEMBLY_MAX */
    {
      mpn_mul_n (tmp, s1p, s2p, nn);
      ecm_redc_n (rp, tmp, 2 * nn, np, invm, nn);
    }
}

static void 
ecm_sqrredc_basecase_n (mp_ptr rp, mp_srcptr s1p,
                        mp_srcptr np, mp_size_t nn, mp_ptr invm, mp_ptr tmp)
{
  mp_limb_t cy;
  mp_size_t j;

  if (nn <= MULREDC_ASSEMBLY_MAX)
    {
      switch (tune_sqrredc_table[nn])
        {
        case MPMOD_MULREDC: /* use quadratic assembly mulredc */
#ifdef USE_ASM_REDC
          mulredc (rp, s1p, s1p, np, nn, invm[0]);
          break;
#endif /* otherwise go through to the next available mode */
        case MPMOD_MUL_REDC1: /* mpn_sqr + __gmpn_redc_1 */
#ifdef HAVE___GMPN_REDC_1
          mpn_sqr (tmp, s1p, nn);
          REDC1(rp, tmp, np, nn, invm[0]);
          break;
#endif /* otherwise go through to the next available mode */
        case MPMOD_MUL_REDC2: /* mpn_sqr + __gmpn_redc_2 */
#ifdef HAVE___GMPN_REDC_2
          mpn_sqr (tmp, s1p, nn);
          REDC2(rp, tmp, np, nn, invm);
          break;
#endif /* otherwise go through to the next available mode */
        case MPMOD_MUL_REDCN: /* mpn_sqr + __gmpn_redc_n */
        /* disable redc_n for now, since it uses the opposite
           precomputed inverse wrt redc_1 and redc_2 */
#ifdef HAVE___GMPN_REDC_Nxxx
          mpn_sqr (tmp, s1p, nn);
          __gmpn_redc_n (rp, tmp, np, nn, invm);
          break;
#endif /* otherwise go through to the next available mode */
        case MPMOD_MUL_REDC_C: /* plain C quadratic reduction */
          mpn_sqr (tmp, s1p, nn);
          for (j = 0; j < nn; j++, tmp++)
            tmp[0] = mpn_addmul_1 (tmp, np, nn, tmp[0] * invm[0]);
          cy = mpn_add_n (rp, tmp - nn, tmp, nn);
          if (cy != 0)
            mpn_sub_n (rp, rp, np, nn); /* a borrow should always occur here */
          break;
        default:
          {
            outputf (OUTPUT_ERROR, "Invalid sqrredc mode: %d\n",
                     tune_sqrredc_table[nn]);
            exit (EXIT_FAILURE);
          }
        }
    }
  else /* nn > MULREDC_ASSEMBLY_MAX */
    {
      mpn_sqr (tmp, s1p, nn);
      ecm_redc_n (rp, tmp, 2 * nn, np, invm, nn);
    }
}

/* R <- S1 * S2 mod modulus
   i.e. R <- S1*S2/r^nn mod n, where n has nn limbs, and r=2^GMP_NUMB_BITS.
   Same as ecm_redc_basecase previous, but combined with mul
   Neither input argument must be in modulus->temp1
*/
static void 
ecm_mulredc_basecase (mpres_t R, const mpres_t S1, const mpres_t S2, 
                      mpmod_t modulus)
{
  mp_ptr s1p, s2p, rp = PTR(R);
  mp_size_t j, nn = modulus->bits / GMP_NUMB_BITS;

  ASSERT(ALLOC(R) >= nn);
  ASSERT(ALLOC(S1) >= nn);
  ASSERT(ALLOC(S2) >= nn);
  s1p = PTR(S1);
  s2p = PTR(S2);
  /* FIXME: S1 and S2 are input and marked const, we mustn't write to them */
  for (j = ABSIZ(S1); j < nn; j++) 
    s1p[j] = 0;
  for (j = ABSIZ(S2); j < nn; j++) 
    s2p[j] = 0;

  ecm_mulredc_basecase_n (rp, s1p, s2p, PTR(modulus->orig_modulus), nn,
                          modulus->Nprim, PTR(modulus->temp1));

  MPN_NORMALIZE (rp, nn);
  SIZ(R) = (SIZ(S1)*SIZ(S2)) < 0 ? (int) -nn : (int) nn;
}

/* R <- S1^2 mod modulus
   i.e. R <- S1^2/r^nn mod n, where n has nn limbs, and r=2^GMP_NUMB_BITS.
   Same as ecm_redc_basecase previous, but combined with sqr
   The input argument must not be in modulus->temp1 */
static void 
ecm_sqrredc_basecase (mpres_t R, const mpres_t S1, mpmod_t modulus)
{
  mp_ptr rp;
  mp_ptr s1p;
  mp_size_t j, nn = modulus->bits / GMP_NUMB_BITS;

  ASSERT(ALLOC(R) >= nn);
  ASSERT(ALLOC(S1) >= nn);
  rp = PTR(R);
  s1p = PTR(S1);
  /* FIXME: S1 is input and marked const, we mustn't write to it */
  for (j = ABSIZ(S1); j < nn; j++)
    s1p[j] = 0;

  ecm_sqrredc_basecase_n (rp, s1p, PTR(modulus->orig_modulus), nn,
                          modulus->Nprim, PTR(modulus->temp1));

  MPN_NORMALIZE (rp, nn);
  SIZ(R) = (int) nn;
}

/* If the user asked for a particular representation, always use it.
   If repr = ECM_MOD_DEFAULT, use the thresholds.
   Don't use base2 if repr = ECM_MOD_NOBASE2.
   If a value is <= -16 or >= 16, it is a base2 exponent.
   Return a non-zero value if an error occurred.
*/
int
mpmod_init (mpmod_t modulus, const mpz_t N, int repr)
{
  int base2 = 0, r = 0;
  mp_size_t n = mpz_size (N);

  switch (repr)
    {
    case ECM_MOD_DEFAULT:
      if ((base2 = isbase2 (N, BASE2_THRESHOLD)))
	{
	  repr = ECM_MOD_BASE2;
	  break;
	}
      /* else go through */
    case ECM_MOD_NOBASE2:
      if (mpz_size (N) < MPZMOD_THRESHOLD)
	repr = ECM_MOD_MODMULN;
      else if (mpz_size (N) < REDC_THRESHOLD)
        repr = ECM_MOD_MPZ;
      else
	repr = ECM_MOD_REDC;
    }

  /* now repr is {ECM_MOD_BASE2, ECM_MOD_MODMULN, ECM_MOD_MPZ, ECM_MOD_REDC},
     or |repr| >= 16. */

  switch (repr)
    {
    case ECM_MOD_MPZ:
      outputf (OUTPUT_VERBOSE, "Using mpz_mod\n");
      mpmod_init_MPZ (modulus, N);
      break;
    case ECM_MOD_MODMULN:
      outputf (OUTPUT_VERBOSE, "Using MODMULN [mulredc:%d, sqrredc:%d]\n",
               (n <= MULREDC_ASSEMBLY_MAX) ? tune_mulredc_table[n] : 4,
               (n <= MULREDC_ASSEMBLY_MAX) ? tune_sqrredc_table[n] : 4);
      mpmod_init_MODMULN (modulus, N);
      break;
    case ECM_MOD_REDC:
      outputf (OUTPUT_VERBOSE, "Using REDC\n");
      mpmod_init_REDC (modulus, N);
      break;
    default: /* base2 case: either repr=ECM_MOD_BASE2, and base2 was
		determined above, or |repr| >= 16, and we want base2 = repr */
      if (repr != ECM_MOD_BASE2)
	base2 = repr;
      r = mpmod_init_BASE2 (modulus, base2, N);
      break;
    }

  return r;
}

void
mpres_clear (mpres_t a, ATTRIBUTE_UNUSED const mpmod_t modulus) 
{
  mpz_clear (a);
  PTR(a) = NULL; /* Make sure we segfault if we access it again */
}

void 
mpmod_init_MPZ (mpmod_t modulus, const mpz_t N)
{
  size_t n;

  mpz_init_set (modulus->orig_modulus, N);
  modulus->repr = ECM_MOD_MPZ;
  
  n = mpz_size (N); /* number of limbs of N */
  modulus->bits = n * GMP_NUMB_BITS; /* Number of bits, 
					rounded up to full limb */

  mpz_init2 (modulus->temp1, 2UL * modulus->bits + GMP_NUMB_BITS);
  mpz_init2 (modulus->temp2, modulus->bits);
  mpz_init2 (modulus->aux_modulus, modulus->bits);
  mpz_set_ui (modulus->aux_modulus, 1UL);
  /* we precompute B^(n + ceil(n/2)) mod N, where B=2^GMP_NUMB_BITS */
  mpz_mul_2exp (modulus->aux_modulus, modulus->aux_modulus,
		(n + (n + 1) / 2) * GMP_NUMB_BITS);
  mpz_mod (modulus->aux_modulus, modulus->aux_modulus, N);

  return;
}

int 
mpmod_init_BASE2 (mpmod_t modulus, const int base2, const mpz_t N)
{
  int Nbits;
  
  outputf (OUTPUT_VERBOSE,
           "Using special division for factor of 2^%d%c1\n",
           abs (base2), (base2 < 0) ? '-' : '+');
  mpz_init_set (modulus->orig_modulus, N);
  modulus->repr = ECM_MOD_BASE2;
  modulus->bits = base2;

  Nbits = mpz_size (N) * GMP_NUMB_BITS; /* Number of bits, rounded
                                           up to full limb */

  mpz_init2 (modulus->temp1, 2UL * Nbits + GMP_NUMB_BITS);
  mpz_init2 (modulus->temp2, Nbits);
  
  mpz_set_ui (modulus->temp1, 1UL);
  mpz_mul_2exp (modulus->temp1, modulus->temp1, abs (base2));
  if (base2 < 0)
    mpz_sub_ui (modulus->temp1, modulus->temp1, 1UL);
  else
    mpz_add_ui (modulus->temp1, modulus->temp1, 1UL);
  if (!mpz_divisible_p (modulus->temp1, N))
    {
       outputf (OUTPUT_ERROR, "mpmod_init_BASE2: n does not divide 2^%d%c1\n",
                abs (base2), base2 < 0 ? '-' : '+');
       mpz_clear (modulus->temp2);
       mpz_clear (modulus->temp1);
       mpz_clear (modulus->orig_modulus);
       return ECM_ERROR;
    }
  
  modulus->Fermat = 0;
  if (base2 > 0)
    {
      unsigned long i;
      for (i = base2; (i & 1) == 0; i >>= 1);
      if (i == 1)
        {
          modulus->Fermat = base2;
        }
    }
  
  return 0;
}

/* initialize the following fields:
   orig_modulus - the original modulus
   bits         - # of bits of N, rounded up to a multiple of GMP_NUMB_BITS
   temp1, temp2 - auxiliary variables
   Nprim        - -1/N mod B^n where B=2^GMP_NUMB_BITS and n = #limbs(N)
   R2           - (2^bits)^2 (mod N)
   R3           - (2^bits)^3 (mod N)
   multiple     - smallest multiple of N >= 2^bits
 */
void
mpmod_init_MODMULN (mpmod_t modulus, const mpz_t N)
{
  int Nbits;

  mpz_init_set (modulus->orig_modulus, N);
  
  modulus->repr = ECM_MOD_MODMULN;
  Nbits = mpz_size (N) * GMP_NUMB_BITS; /* Number of bits, rounded
                                           up to full limb */
  modulus->bits = Nbits;

  mpz_init2 (modulus->temp1, 2UL * Nbits + GMP_NUMB_BITS);
  mpz_init2 (modulus->temp2, Nbits + 1);
  modulus->Nprim = (mp_limb_t*) malloc (mpz_size (N) * sizeof (mp_limb_t));

  mpz_init2 (modulus->R2, Nbits);
  mpz_set_ui (modulus->temp1, 1UL);
  mpz_mul_2exp (modulus->temp1, modulus->temp1, 2 * Nbits);
  mpz_mod (modulus->R2, modulus->temp1, modulus->orig_modulus);
  /* Now R2 = (2^bits)^2 (mod N) */
  
  mpz_init2 (modulus->R3, Nbits);
  mpz_mul_2exp (modulus->temp1, modulus->R2, Nbits);
  mpz_mod (modulus->R3, modulus->temp1, modulus->orig_modulus);
  /* Now R3 = (2^bits)^3 (mod N) */

  mpz_init2 (modulus->multiple, Nbits);
  mpz_set_ui (modulus->temp1, 1UL);
  mpz_mul_2exp (modulus->temp1, modulus->temp1, Nbits);
  /* compute ceil(2^bits / N) */
  mpz_cdiv_q (modulus->temp1, modulus->temp1, modulus->orig_modulus);
  mpz_mul (modulus->multiple, modulus->temp1, modulus->orig_modulus);
  /* Now multiple is the smallest multiple of N >= 2^bits */

  mpz_set_ui (modulus->temp1, 1UL);
  mpz_mul_2exp (modulus->temp1, modulus->temp1, Nbits);
  /* since we directly check even modulus in ecm/pm1/pp1,
     N is odd here, thus 1/N mod 2^Nbits always exist */
  mpz_invert (modulus->temp2, N, modulus->temp1); /* temp2 = 1/N mod B^n */
  mpz_sub (modulus->temp2, modulus->temp1, modulus->temp2);
  /* temp2 = -1/N mod B^n */
  /* ensure Nprim has all its n limbs correctly set, for ecm_redc_n */
  MPN_ZERO(modulus->Nprim, mpz_size (N));
  mpn_copyi (modulus->Nprim, PTR(modulus->temp2), ABSIZ(modulus->temp2));
}

void 
mpmod_init_REDC (mpmod_t modulus, const mpz_t N)
{
  mp_size_t n;
  int Nbits;
  
  mpz_init_set (modulus->orig_modulus, N);
  
  n = mpz_size (N);
  modulus->repr = ECM_MOD_REDC;
  Nbits = n * GMP_NUMB_BITS; /* Number of bits, rounded
                                up to full limb */
  modulus->bits = Nbits;
  
  mpz_init2 (modulus->temp1, 2 * Nbits + GMP_NUMB_BITS);
  mpz_init2 (modulus->temp2, Nbits);
  mpz_init2 (modulus->aux_modulus, Nbits);

  mpz_set_ui (modulus->temp1, 1UL);
  mpz_mul_2exp (modulus->temp1, modulus->temp1, Nbits);
  /* since we directly check even modulus in ecm/pm1/pp1,
     N is odd here, thus 1/N mod 2^Nbits always exist */
  mpz_invert (modulus->aux_modulus, N, modulus->temp1);

  mpz_sub (modulus->aux_modulus, modulus->temp1, modulus->aux_modulus);
  /* ensure aux_modulus has n allocated limbs, for ecm_redc_n */
  if (ABSIZ(modulus->aux_modulus) < n)
    {
      _mpz_realloc (modulus->aux_modulus, n);
      /* in case the reallocation fails, _mpz_realloc sets the value to 0 */
      ASSERT_ALWAYS (mpz_cmp_ui (modulus->aux_modulus, 0) != 0);
      MPN_ZERO (PTR(modulus->aux_modulus) + ABSIZ(modulus->aux_modulus),
		n - ABSIZ(modulus->aux_modulus));
    }

  mpz_init2 (modulus->R2, Nbits);
  mpz_set_ui (modulus->temp1, 1UL);
  mpz_mul_2exp (modulus->temp1, modulus->temp1, 2 * Nbits);
  mpz_mod (modulus->R2, modulus->temp1, modulus->orig_modulus);
  /* Now R2 = (2^bits)^2 (mod N) */
  
  mpz_init2 (modulus->R3, Nbits);
  mpz_mul_2exp (modulus->temp1, modulus->R2, Nbits);
  mpz_mod (modulus->R3, modulus->temp1, modulus->orig_modulus);
  /* Now R3 = (2^bits)^3 (mod N) */
  
  mpz_init (modulus->multiple);
  mpz_set_ui (modulus->temp1, 1UL);
  mpz_mul_2exp (modulus->temp1, modulus->temp1, Nbits);
  /* compute ceil(2^bits / N) */
  mpz_cdiv_q (modulus->temp1, modulus->temp1, modulus->orig_modulus);
  mpz_mul (modulus->multiple, modulus->temp1, modulus->orig_modulus);
  /* Now multiple is the largest multiple of N >= 2^bits */
}

void 
mpmod_clear (mpmod_t modulus)
{
  mpz_clear (modulus->orig_modulus);
  mpz_clear (modulus->temp1);
  mpz_clear (modulus->temp2);
  if (modulus->repr == ECM_MOD_REDC || modulus->repr == ECM_MOD_MPZ)
    mpz_clear (modulus->aux_modulus);
  if (modulus->repr == ECM_MOD_MODMULN || modulus->repr == ECM_MOD_REDC)
    {
      mpz_clear (modulus->R2);
      mpz_clear (modulus->R3);
      mpz_clear (modulus->multiple);
    }
  if (modulus->repr == ECM_MOD_MODMULN)
    free (modulus->Nprim);
  
  return;
}

/* initialize r and set all entries from those of modulus */
void
mpmod_init_set (mpmod_t r, const mpmod_t modulus)
{
  const unsigned long Nbits = abs(modulus->bits);
  const unsigned long n = mpz_size (modulus->orig_modulus);

  r->repr = modulus->repr;
  r->bits = modulus->bits;
  r->Fermat = modulus->Fermat;
  mpz_init_set (r->orig_modulus, modulus->orig_modulus);
  mpz_init2 (r->temp1, 2 * Nbits + GMP_NUMB_BITS);
  mpz_init2 (r->temp2, Nbits + GMP_NUMB_BITS);
  if (modulus->repr == ECM_MOD_MODMULN || modulus->repr == ECM_MOD_REDC)
    {
      mpz_init2 (r->multiple, Nbits);
      mpz_init2 (r->R2, Nbits);
      mpz_init2  (r->R3, Nbits);
      mpz_set (r->multiple, modulus->multiple);
      mpz_set (r->R2, modulus->R2);
      mpz_set (r->R3, modulus->R3);
    }
  if (modulus->repr == ECM_MOD_REDC || modulus->repr == ECM_MOD_MPZ)
    {
      mpz_init2 (r->aux_modulus, Nbits);
      mpz_set (r->aux_modulus, modulus->aux_modulus);
      /* for ECM_MOD_REDC, ensure r->aux_modulus has n limbs */
      _mpz_realloc (r->aux_modulus, n);
      MPN_ZERO (PTR(r->aux_modulus) + ABSIZ(r->aux_modulus),
                n - ABSIZ(r->aux_modulus));
    }
  if (modulus->repr == ECM_MOD_MODMULN)
    {
      r->Nprim = (mp_limb_t*) malloc (n * sizeof (mp_limb_t));
      mpn_copyi (r->Nprim, modulus->Nprim, n);
    }
}


void 
mpres_init (mpres_t R, const mpmod_t modulus)
{
  /* use mpz_sizeinbase since modulus->bits may not be initialized yet */
  mpz_init2 (R, mpz_sizeinbase (modulus->orig_modulus, 2) + GMP_NUMB_BITS);
}

/* realloc R so that it has at least the same number of limbs as modulus */
void
mpres_realloc (mpres_t R, const mpmod_t modulus)
{
  if (modulus->repr == ECM_MOD_MODMULN)
    MPZ_REALLOC (R, modulus->bits / GMP_NUMB_BITS);
}


/* a <- b^2 mod (modulus).
   a and b might be equal, but cannot be modulus->temp1.
   Assumes repr = ECM_MOD_BASE2 or ECM_MOD_MODMULN or ECM_MOD_REDC */
static inline void
mpres_pow_sqr (mpres_t a, const mpres_t b, mpmod_t modulus)
{
  /* Set temp2 = temp2*temp2 */
  if (modulus->repr == ECM_MOD_BASE2)
    {
      mpz_mul (modulus->temp1, b, b);
      base2mod (a, modulus->temp1, modulus->temp1, modulus);
    }
  else if (modulus->repr == ECM_MOD_MODMULN)
    ecm_sqrredc_basecase (a, b, modulus);
  else
    {
      ASSERT(modulus->repr == ECM_MOD_REDC);
      mpz_mul (modulus->temp1, b, b);
      REDC (a, modulus->temp1, modulus->temp2, modulus);
    }
}

/* a <- b*c mod (modulus)
   a, b, c must not equal modulus->temp1.
   Assumes repr = ECM_MOD_BASE2 or ECM_MOD_MODMULN or ECM_MOD_REDC */
static inline void
mpres_pow_mul (mpres_t a, const mpres_t b, const mpres_t c, mpmod_t modulus)
{
  /* Set temp2 = temp2*temp2 */
  if (modulus->repr == ECM_MOD_BASE2)
    {
      mpz_mul (modulus->temp1, b, c);
      base2mod (a, modulus->temp1, modulus->temp1, modulus);
    }
  else if (modulus->repr == ECM_MOD_MODMULN)
    ecm_mulredc_basecase (a, b, c, modulus);
  else
    {
      ASSERT(modulus->repr == ECM_MOD_REDC);
      mpz_mul (modulus->temp1, b, c);
      REDC (a, modulus->temp1, modulus->temp2, modulus);
    }
}

/* R <- BASE^EXP mod modulus.
   Warning: EXP can be negative (only occurs for P-1 in revision 2982).
 */
void 
mpres_pow (mpres_t R, const mpres_t BASE, const mpz_t EXP, mpmod_t modulus)
{
    size_t k, expnbits, K, lw, n0, i;
    mpres_t *B;
    mp_limb_t w;

    ASSERT_NORMALIZED(BASE);

  if (modulus->repr == ECM_MOD_MPZ)
    mpz_powm (R, BASE, EXP, modulus->orig_modulus);
  else if (modulus->repr == ECM_MOD_BASE2 || modulus->repr == ECM_MOD_MODMULN ||
           modulus->repr == ECM_MOD_REDC)
    {
      size_t expidx;
      mp_limb_t bitmask, expbits;
      int sgn = mpz_sgn (EXP);
      mpz_t exp;

      /* case EXP=0 */
      if (sgn == 0)
        {
          mpres_set_ui (R, 1UL, modulus); /* set result to 1 */
          ASSERT_NORMALIZED (R);
          return;
        }

      /* exp <- |EXP| */
      PTR(exp) = PTR(EXP);
      SIZ(exp) = ABSIZ(EXP);

      expidx = mpz_size (exp) - 1;         /* point at most significant limb */
      expbits = mpz_getlimbn (exp, expidx); /* get most significant limb */
      ASSERT (expbits != 0);

      /* Scan for the MSB in expbits */
      bitmask = ((mp_limb_t) 1) << (GMP_NUMB_BITS - 1);
      for (; (bitmask & expbits) == 0; bitmask >>= 1);
    
      /* here the most significant limb with any set bits is in expbits, */
      /* bitmask is set to mask in the msb of expbits */

      k = 1; /* sliding window size */
      expnbits = mpz_sizeinbase (exp, 2);
      /* the average number of multiplications is 2^(k-1) + expnbits / (k+1) */
      while ((1 << (k-1)) + expnbits / (k + 1) > (1 << k) + expnbits / (k + 2))
        k ++;
      /* precompute BASE^i for i = 2, 3, 5, ..., 2^k-1 */
      K = 1UL << (k - 1);
      B = malloc (K * sizeof (mpres_t));
      for (i = 0; i < K; i++)
        {
          mpres_init (B[i], modulus);
          mpres_realloc (B[i], modulus);
          if (i == 0) /* store BASE^2 in B[0] */
            mpres_pow_sqr (B[i], BASE, modulus);
          else if (i == 1)
            mpres_pow_mul (B[i], B[0], BASE, modulus);
          else
            mpres_pow_mul (B[i], B[i-1], B[0], modulus);
        }
      
      mpz_set (modulus->temp2, BASE);
      bitmask >>= 1;
      w = 0; /* invariant: temp2 has to be multiplied by BASE^w */
      lw = 0;   /* number of bits in w */
      expnbits --;     /* number of remaining bits to deal with */
      n0 = 0;   /* smallest bit index of current window */

      while (1) 
        {
          for ( ; bitmask != 0; bitmask >>= 1, expnbits--)
            {
              /* Set temp2 = temp2 ^ 2 */
              mpres_pow_sqr (modulus->temp2, modulus->temp2, modulus);
              w = w << 1;
              lw += (lw != 0);

              /* If bit is 1, accumulate in w */
              if (expbits & bitmask)
                {
                  w ++;
                  if (lw == 0)
                    {
                      lw = 1;
                      /* n0 points to the smallest bit of the current window */
                      n0 = (expnbits >= k) ? expnbits - k : 0;
                    }
                  /* if the current bit is the smallest one of the window,
                     we multiply by BASE^w */
                  if (mpz_scan1 (exp, n0) == expnbits - 1)
                    {
                      ASSERT(w/2 < K);
                      mpres_pow_mul (modulus->temp2, (w == 1) ? BASE : B[w/2],
                                     modulus->temp2, modulus);
                      w = lw = 0;
                    }
                }
            }
          if (expidx == 0)		/* if we just processed the least */
            break;			/* significant limb, we are done */
          expidx --;
          expbits = mpz_getlimbn (exp, expidx);
          bitmask = (mp_limb_t) 1 << (GMP_NUMB_BITS - 1);
        }
      ASSERT(w/2 < K);
      if (w != 0)
        mpres_pow_mul (modulus->temp2, (w == 1) ? BASE : B[w/2],
                       modulus->temp2, modulus);
      for (i = 0; i < K; i++)
        mpres_clear (B[i], modulus);
      free (B);
      mpz_set (R, modulus->temp2);

      /* mpz_getlimbn() ignores sign of argument, so we computed BASE^|EXP|.
         If EXP was negative, do a modular inverse */
      if (sgn < 0)
        mpres_invert (R, R, modulus);
    } /* if (modulus->repr == ECM_MOD_BASE2 || ... ) */
  ASSERT_NORMALIZED (R);
}


/* Returns 1 if S == 0 (mod modulus), 0 otherwise */

int
mpres_is_zero (const mpres_t S, mpmod_t modulus)
{
  mpz_mod (modulus->temp1, S, modulus->orig_modulus);
  /* For all currently implemented representations, a zero residue has zero
     integer representation */
  return (mpz_sgn (modulus->temp1) == 0) ? 1 : 0;
}

/* R <- BASE^EXP mod modulus */ 
void 
mpres_ui_pow (mpres_t R, const unsigned long BASE, const mpres_t EXP, 
              mpmod_t modulus)
{
  if (modulus->repr == ECM_MOD_MPZ)
    {
      mpz_set_ui (modulus->temp1, BASE);
      mpz_powm (R, modulus->temp1, EXP, modulus->orig_modulus);
    }
  else if (modulus->repr == ECM_MOD_BASE2 || modulus->repr == ECM_MOD_MODMULN ||
           modulus->repr == ECM_MOD_REDC)
    {
      size_t expidx;
      mp_limb_t bitmask, expbits;

      expidx = mpz_size (EXP) - 1;         /* point at most significant limb */
      expbits = mpz_getlimbn (EXP, expidx); /* get most significant limb */
      bitmask = (mp_limb_t) 1 << (GMP_NUMB_BITS - 1);

      /* case EXP=0 */
      if (mpz_sgn (EXP) == 0)
        {
          mpres_set_ui (R, 1UL, modulus); /* set result to 1 */
          ASSERT_NORMALIZED (R);
          return;
        }

      ASSERT (mpz_size (EXP) > 0);         /* probably redundant with _sgn() test */
      expidx = mpz_size (EXP) - 1;         /* point at most significant limb */
      expbits = mpz_getlimbn (EXP, expidx); /* get most significant limb */
      ASSERT (expbits != 0);

      /* Scan for the MSB in expbits */
      bitmask = ((mp_limb_t) 1) << (GMP_NUMB_BITS - 1);
      for (; (bitmask & expbits) == 0; bitmask >>= 1);
    
      /* here the most significant limb with any set bits is in expbits, */
      /* bitmask is set to mask in the msb of expbits */

      mpz_set_ui (modulus->temp2, BASE); /* temp2 = BASE */
      if (modulus->repr == ECM_MOD_MODMULN || modulus->repr == ECM_MOD_REDC)
        {
          mpz_mul_2exp (modulus->temp1, modulus->temp2, modulus->bits);
          mpz_mod (modulus->temp2, modulus->temp1, modulus->orig_modulus);
        }
      bitmask >>= 1;


      while (1) 
        {
          for ( ; bitmask != 0; bitmask >>= 1) 
            {
              mpres_pow_sqr (modulus->temp2, modulus->temp2, modulus);

              /* If bit is 1, set temp2 = temp2 * BASE */
              if (expbits & bitmask)
                {
                  if (BASE == 2UL)
                    {
                      mpz_mul_2exp (modulus->temp2, modulus->temp2, 1);
                      if (mpz_cmp (modulus->temp2, modulus->orig_modulus) >= 0)
                        mpz_sub (modulus->temp2, modulus->temp2, modulus->orig_modulus);
                    }
                  else
                    {
                      mpz_mul_ui (modulus->temp1, modulus->temp2, BASE);
                      mpz_mod (modulus->temp2, modulus->temp1, modulus->orig_modulus);
                    }
                }
            }
          if (expidx == 0)		/* if we just processed the least */
            break;			/* significant limb, we are done */
          expidx--;
          expbits = mpz_getlimbn (EXP, expidx);
          bitmask = (mp_limb_t) 1 << (GMP_NUMB_BITS - 1);
        }
      mpz_set (R, modulus->temp2);

      /* mpz_getlimbn() ignores sign of argument, so we computed BASE^|EXP|.
         If EXP was negative, do a modular inverse */
      if (mpz_sgn (EXP) < 0)
        {
          mpres_invert (R, R, modulus);
        }
    } /* if (modulus->repr == ECM_MOD_BASE2 || ... ) */
  ASSERT_NORMALIZED (R);
}

/* We use here the algorithm described in "Fast Modular Reduction" from
   Hasenplaugh, Gaubatz and Gobal, Arith'18, 2007: assuming N has n limbs,
   we have precomputed C = B^(n + ceil(n/2)) mod N. */
static void
mpres_mpz_mod (mpres_t R, mpz_t T, mpz_t N, mpz_t C)
{
  size_t n = mpz_size (N);
  size_t t = mpz_size (T);
  size_t m = n + (n + 1) / 2; /* n + ceil(n/2) */

  if (t > m && n > 1) /* if n=1, then m=2, thus h=0 */
    {
      size_t c = mpz_size (C);
      size_t h, l;
      mp_ptr rp;
      mp_ptr tp = PTR(T);

      /* Warning: we might have t > 2n. In that case we reduce
	 {tp+l+m, t-(m+l)} where l = t-2n. */
      l = (t > 2 * n) ? t - 2 * n : 0;

      tp += l;
      h = t - (m + l); /* since t-l <= 2n and m = n + ceil(n/2),
			  we have h <= n - ceil(n/2) = floor(n/2).
                          On the other hand, if l=0 we have h = t-m > 0;
                          if l>0, then l=t-2n, thus h=2n-m = floor(n/2) > 0
                          since n > 1. */
      mpz_realloc (R, c + h);
      rp = PTR(R);
      if (c > h)
	mpn_mul (rp, PTR(C), c, tp + m, h);
      else
	mpn_mul (rp, tp + m, h, PTR(C), c);
      /* now add {rp, c+h} to {tp, m}: we have c <= n and h <= n/2,
	 thus c + h <= m */
      if (c + h > m) abort();
      tp[m] = mpn_add (tp, tp, m, rp, c + h);
      m += l + tp[m];
      tp -= l; /* put back the low l limbs */
      MPN_NORMALIZE(tp, m);
      SIZ(T) = (SIZ(T) > 0) ? m : -m;
    }
  
  mpz_mod (R, T, N);
}

void 
mpres_mul (mpres_t R, const mpres_t S1, const mpres_t S2, mpmod_t modulus)
{
  ASSERT_NORMALIZED (S1);
  ASSERT_NORMALIZED (S2);

#ifdef WANT_ASSERT_EXPENSIVE
  mpz_t test1, test2, test_result1, test_result2;
  ASSERT_ALWAYS (S1 != modulus->temp1 && S2 != modulus->temp1 &&
                 R != modulus->temp1);
  mpz_init (test1);
  mpz_init (test2);
  mpz_init (test_result1);
  mpz_init (test_result2);
  mpres_get_z (test1, S1, modulus);
  mpres_get_z (test2, S2, modulus);
  mpz_mul (test_result1, test1, test2);
  mpz_mod (test_result1, test_result1, modulus->orig_modulus);
#endif

  if (UNLIKELY(modulus->repr == ECM_MOD_BASE2 && modulus->Fermat >= 32768))
    {
      mp_size_t n = modulus->Fermat / GMP_NUMB_BITS;
      unsigned long k;
      mp_srcptr s1p, s2p;
      mp_size_t s1s, s2s;
      
      MPZ_REALLOC (R, n + 1);
      s1p = PTR(S1);
      s1s = SIZ(S1);
      s2p = PTR(S2);
      s2s = SIZ(S2);
      
      k = mpn_fft_best_k (n, S1 == S2);
      ASSERT(mpn_fft_next_size (n, k) == n);

      if (base2mod_2 (modulus->temp1, S1, n, modulus->orig_modulus))
        {
          s1p = PTR(modulus->temp1);
          s1s = SIZ(modulus->temp1);
        }
      if (S1 == S2)
        {
          s2p = s1p;
          s2s = s1s;
        }
      else if (base2mod_2 (modulus->temp2, S2, n, modulus->orig_modulus))
        {
          s2p = PTR(modulus->temp2);
          s2s = SIZ(modulus->temp2);
        }

      /* mpn_mul_fft() computes the product modulo B^n + 1, where 
         B = 2^(machine word size in bits). So the result can be = B^n, 
         in that case R is set to zero and 1 is returned as carry-out.
         In all other cases 0 is returned. Hence the complete result is 
         R + cy * B^n, where cy is the value returned by mpn_mul_fft(). */
      PTR(R)[n] = mpn_mul_fft (PTR(R), n, s1p, ABS(s1s), s2p, ABS(s2s), k);
      n ++;
      MPN_NORMALIZE(PTR(R), n);
      SIZ(R) = ((s1s ^ s2s) >= 0) ? (int) n : (int) -n;

      return;
    }

  switch (modulus->repr)
    {
    case ECM_MOD_BASE2:
      mpz_mul (modulus->temp1, S1, S2);
      base2mod (R, modulus->temp1, modulus->temp1, modulus);
      break;
    case ECM_MOD_MODMULN:
      MPZ_REALLOC (R, modulus->bits / GMP_NUMB_BITS);
      ecm_mulredc_basecase (R, S1, S2, modulus);
      break;
    case ECM_MOD_REDC:
      mpz_mul (modulus->temp1, S1, S2);
      REDC (R, modulus->temp1, modulus->temp2, modulus);
      break;
    default: /* case ECM_MOD_MPZ */
      mpz_mul (modulus->temp1, S1, S2);
      mpres_mpz_mod (R, modulus->temp1, modulus->orig_modulus,
		     modulus->aux_modulus);
      break;
    }
  ASSERT_NORMALIZED (R);

#ifdef WANT_ASSERT_EXPENSIVE
  mpres_get_z (test_result2, R, modulus);
  if (mpz_cmp (test_result1, test_result2) != 0)
    {
      printf ("mpres_mul and mpz_mul/mpz_mod produced different results.\n");
      gmp_printf ("input 1:         %Zd\n", test1);
      gmp_printf ("input 2:         %Zd\n", test2);
      gmp_printf ("mpres_mul:       %Zd\n", test_result2);
      gmp_printf ("mpz_mul/mpz_mod: %Zd\n", test_result1);
      abort ();
    }
  mpz_clear (test1);
  mpz_clear (test2);
  mpz_clear (test_result1);
  mpz_clear (test_result2);
#endif
}

/* R <- S1^2 mod modulus */
void 
mpres_sqr (mpres_t R, const mpres_t S1, mpmod_t modulus)
{
  ASSERT_NORMALIZED (S1);

#ifdef WANT_ASSERT_EXPENSIVE
  mpz_t test1, test_result1, test_result2;
  ASSERT_ALWAYS (S1 != modulus->temp1 && R != modulus->temp1);
  mpz_init (test1);
  mpz_init (test_result1);
  mpz_init (test_result2);
  mpres_get_z (test1, S1, modulus);
  mpz_mul (test_result1, test1, test1);
  mpz_mod (test_result1, test_result1, modulus->orig_modulus);
#endif

  if (UNLIKELY(modulus->repr == ECM_MOD_BASE2 && modulus->Fermat >= 32768))
    {
      mpres_mul (R, S1, S1, modulus);
      return;
    }

  switch (modulus->repr)
    {
    case ECM_MOD_BASE2:
      mpz_mul (modulus->temp1, S1, S1);
      base2mod (R, modulus->temp1, modulus->temp1, modulus);
      break;
    case ECM_MOD_MODMULN:
      MPZ_REALLOC (R, modulus->bits / GMP_NUMB_BITS);
      ecm_sqrredc_basecase (R, S1, modulus);
      break;
    case ECM_MOD_REDC:
      mpz_mul (modulus->temp1, S1, S1);
      REDC (R, modulus->temp1, modulus->temp2, modulus);
      break;
    default: /* case ECM_MOD_MPZ */
      mpz_mul (modulus->temp1, S1, S1);
      mpres_mpz_mod (R, modulus->temp1, modulus->orig_modulus,
		     modulus->aux_modulus);
      break;
    }
  ASSERT_NORMALIZED (R);

#ifdef WANT_ASSERT_EXPENSIVE
  mpres_get_z (test_result2, R, modulus);
  if (mpz_cmp (test_result1, test_result2) != 0)
    {
      printf ("mpres_sqr and mpz_mul/mpz_mod produced different results.\n");
      gmp_printf ("input 1:         %Zd\n", test1);
      gmp_printf ("mpres_mul:       %Zd\n", test_result2);
      gmp_printf ("mpz_mul/mpz_mod: %Zd\n", test_result1);
      abort ();
    }
  mpz_clear (test1);
  mpz_clear (test_result1);
  mpz_clear (test_result2);
#endif
}

/* R <- S * n mod modulus */
void 
mpres_mul_ui (mpres_t R, const mpres_t S, const unsigned long n, 
              mpmod_t modulus)
{
  ASSERT_NORMALIZED (S);
  mpz_mul_ui (modulus->temp1, S, n);
  /* This is the same for all methods: just reduce with original modulus */
  mpz_mod (R, modulus->temp1, modulus->orig_modulus);
  ASSERT_NORMALIZED (R);
}

/* R <- S * 2^k mod modulus */
void 
mpres_mul_2exp (mpres_t R, const mpres_t S, const unsigned long k, 
              mpmod_t modulus)
{
  ASSERT_NORMALIZED (S);
  mpz_mul_2exp (modulus->temp1, S, k);
  /* This is the same for all methods: just reduce with original modulus */
  mpz_mod (R, modulus->temp1, modulus->orig_modulus);
  ASSERT_NORMALIZED (R);
}

/* This function multiplies an integer in mpres_t form with an integer in
   mpz_t form, and stores the output in mpz_t form. The advantage is that
   one REDC suffices to reduce the product and convert it to non-Montgomery
   representation. */

void 
mpres_mul_z_to_z (mpz_t R, const mpres_t S1, const mpz_t S2, mpmod_t modulus)
{
  ASSERT_NORMALIZED (S1);

  if (modulus->repr == ECM_MOD_BASE2 && modulus->Fermat >= 32768)
    {
      mp_size_t n = modulus->Fermat / GMP_NUMB_BITS;
      unsigned long k;
      mp_srcptr s1p = PTR(S1), s2p = PTR(S2);
      mp_size_t s1s = SIZ(S1), s2s = SIZ(S2);

      MPZ_REALLOC (R, n + 1);
      k = mpn_fft_best_k (n, S1 == S2);
      ASSERT(mpn_fft_next_size (n, k) == n);

      if (base2mod_2 (modulus->temp1, S1, n, modulus->orig_modulus))
        {
          s1p = PTR(modulus->temp1);
          s1s = SIZ(modulus->temp1);
        }
      if (S1 == S2)
        {
          s2p = s1p;
          s2s = s1s;
        }
      else if (base2mod_2 (modulus->temp2, S2, n, modulus->orig_modulus))
        {
          s2p = PTR(modulus->temp2);
          s2s = SIZ(modulus->temp2);
        }

      /* mpn_mul_fft() computes the product modulo B^n + 1, where 
         B = 2^(machine word size in bits). So the result can be = B^n, 
         in that case R is set to zero and 1 is returned as carry-out.
         In all other cases 0 is returned. Hence the complete result is 
         R + cy * B^n, where cy is the value returned by mpn_mul_fft(). */
      PTR(R)[n] = mpn_mul_fft (PTR(R), n, s1p, ABS(s1s), s2p, ABS(s2s), k);
      n ++;
      MPN_NORMALIZE(PTR(R), n);
      SIZ(R) = ((s1s ^ s2s) >= 0) ? (int) n : (int) -n;
      mpz_mod (R, R, modulus->orig_modulus);

      return;
    }

  switch (modulus->repr)
    {
    case ECM_MOD_BASE2:
      if (mpz_sizeinbase (S2, 2) > (unsigned) abs (modulus->bits))
	{
	  base2mod (modulus->temp2, S2, modulus->temp1, modulus);
	  mpz_mul (modulus->temp1, S1, modulus->temp2);
	}
      else
	mpz_mul (modulus->temp1, S1, S2);
      base2mod (R, modulus->temp1, modulus->temp1, modulus);
      mpz_mod (R, R, modulus->orig_modulus);
      break;
    case ECM_MOD_MODMULN:
      if (mpz_cmp (S2, modulus->orig_modulus) >= 0)
	{
	  mpz_mod (modulus->temp2, S2, modulus->orig_modulus);
          MPZ_REALLOC (R, modulus->bits / GMP_NUMB_BITS);
 	  ecm_mulredc_basecase (R, S1, modulus->temp2, modulus);
          mpz_mod (R, R, modulus->orig_modulus);
	}
      else
 	{
	  MPZ_REALLOC (R, modulus->bits / GMP_NUMB_BITS);
	  if (ABSIZ(S2) < modulus->bits / GMP_NUMB_BITS)
	    {
	      mpz_t t;
	      mpz_init2 (t, modulus->bits);
	      mpz_set (t, S2);
	      ecm_mulredc_basecase (R, S1, t, modulus);
	      mpz_clear (t);
            }
          else
            ecm_mulredc_basecase (R, S1, S2, modulus);
          mpz_mod (R, R, modulus->orig_modulus);
 	}
      break;
    case ECM_MOD_REDC:
      if (mpz_cmp (S2, modulus->orig_modulus) >= 0)
	{
	  mpz_mod (modulus->temp2, S2, modulus->orig_modulus);
	  mpz_mul (modulus->temp1, S1, modulus->temp2);
	}
      else
	mpz_mul (modulus->temp1, S1, S2);
      REDC (R, modulus->temp1, modulus->temp2, modulus);
      mpz_mod (R, R, modulus->orig_modulus);
      break;
    default:
      if (mpz_cmp (S2, modulus->orig_modulus) >= 0)
	{
	  mpz_mod (modulus->temp2, S2, modulus->orig_modulus);
	  mpz_mul (modulus->temp1, S1, modulus->temp2);
	}
      else
	mpz_mul (modulus->temp1, S1, S2);
      mpz_mod (R, modulus->temp1, modulus->orig_modulus);
      break;
    }
  ASSERT_NORMALIZED (R);
}


/* Sets R = S * c, for some constant c that is coprime to modulus.
   This is primarily useful for multiplying numbers together for a gcd with
   modulus. The advantage is that we don't need to convert the mpz_t 
   to Montgomery representation before applying REDC. */

void 
mpres_set_z_for_gcd (mpres_t R, const mpz_t S, mpmod_t modulus)
{
  mpz_mod (R, S, modulus->orig_modulus);
  ASSERT_NORMALIZED (R);  
}

/* Companion function to mpres_set_z_for_gcd(). It divides by c^n. 
   In case of MULREDC with k b-bit words, c = 1/2^(b*k), so we multiply 
   by 2^(n*b*k). The purpose is to fix products of n terms collected by 
   using mpres_set_z_for_gcd(), so that we can still get the exact residue. */

void 
mpres_set_z_for_gcd_fix (mpres_t R, const mpres_t S, const mpz_t n, mpmod_t modulus)
{
  switch (modulus->repr)
    {
      case ECM_MOD_MODMULN:
      case ECM_MOD_REDC:
        {
          mpres_t po2;
          mpz_t nb;

          mpz_init (nb);
          mpres_init (po2, modulus);

          /* Note: modulus->bits is always non-negative for MODMULN and REDC */
          mpz_mul_ui (nb, n, modulus->bits);
          mpres_set_ui (po2, 2, modulus);
          mpres_pow (po2, po2, nb, modulus);
          mpres_mul (R, S, po2, modulus);

          mpz_clear (nb);
          mpres_clear (po2, modulus);
        }
    }
}


/* R <- S / 2^n mod modulus. Does not need to be fast. */
void 
mpres_div_2exp (mpres_t R, const mpres_t S, const unsigned int n, 
                mpmod_t modulus)
{
  int i;
  ASSERT_NORMALIZED (S);

  if (n == 0)
    {
      mpres_set (R, S, modulus);
      ASSERT_NORMALIZED (R);
      return;
    }

    if (mpz_odd_p (S))
      {
        ASSERT (mpz_odd_p (modulus->orig_modulus));
        mpz_add (R, S, modulus->orig_modulus);
        mpz_tdiv_q_2exp (R, R, 1);
      }
    else
      mpz_tdiv_q_2exp (R, S, 1);

    for (i = n ; i > 1; i--)
      {
        if (mpz_odd_p (R))
          {
            ASSERT (mpz_odd_p (modulus->orig_modulus));
            mpz_add (R, R, modulus->orig_modulus);
          }
        mpz_tdiv_q_2exp (R, R, 1);
      }

    ASSERT_NORMALIZED (R);
}

void
mpres_add_ui (mpres_t R, const mpres_t S, const unsigned long n, 
              mpmod_t modulus)
{
  ASSERT_NORMALIZED (S);
  if (modulus->repr == ECM_MOD_MPZ || modulus->repr == ECM_MOD_BASE2)
    {
      mpz_add_ui (R, S, n);
      if (mpz_cmp (R, modulus->orig_modulus) > 0)
        mpz_sub (R, R, modulus->orig_modulus); /* This assumes modulus >= n */
    }
  else if (modulus->repr == ECM_MOD_MODMULN || modulus->repr == ECM_MOD_REDC)
    {
      mpz_set_ui (modulus->temp1, n);
      mpz_mul_2exp (modulus->temp1, modulus->temp1, modulus->bits);
      mpz_add (modulus->temp1, modulus->temp1, S);
      mpz_mod (R, modulus->temp1, modulus->orig_modulus);
    }
  ASSERT_NORMALIZED (R);
}

/* R <- S1 + S2 mod modulus */
void 
mpres_add (mpres_t R, const mpres_t S1, const mpres_t S2, mpmod_t modulus)
{
  ASSERT_NORMALIZED (S1);
  ASSERT_NORMALIZED (S2);
  mpz_add (R, S1, S2);
  if ((modulus->repr == ECM_MOD_MODMULN || modulus->repr == ECM_MOD_REDC) &&
      ABSIZ(R) > ABSIZ(modulus->orig_modulus))
    {
      if (SIZ(R) > 0)
	mpz_sub (R, R, modulus->multiple);
      else
	mpz_add (R, R, modulus->multiple);
      /* N <= since multiple < 2^Nbits + N, now |R| < B */
    }
  ASSERT_NORMALIZED (R);
}

/* R <- S - n mod modulus
   If repr == ECM_MOD_MODMULN or ECM_MOD_REDC, we need to convert n to
   Montgomery representation before substracting
*/
void
mpres_sub_ui (mpres_t R, const mpres_t S, const unsigned long n, 
              mpmod_t modulus)
{
  ASSERT_NORMALIZED (S);
  if (modulus->repr == ECM_MOD_MPZ || modulus->repr == ECM_MOD_BASE2)
    {
      mpz_sub_ui (R, S, n);
      if (mpz_sgn (R) < 0)
        mpz_add (R, R, modulus->orig_modulus); /* Assumes modulus >= n */
    }
  else if (modulus->repr == ECM_MOD_MODMULN || modulus->repr == ECM_MOD_REDC)
    {
      mpz_set_ui (modulus->temp1, n);
      mpz_mul_2exp (modulus->temp1, modulus->temp1, modulus->bits);
      mpz_sub (modulus->temp1, S, modulus->temp1);
      mpz_mod (R, modulus->temp1, modulus->orig_modulus);
    }
  ASSERT_NORMALIZED (R);
}

#if 0
/* R <- n - S mod modulus
   If repr == ECM_MOD_MODMULN or ECM_MOD_REDC, we need to convert n to
   Montgomery representation before substracting
*/
void
mpres_ui_sub (mpres_t R, const unsigned long n ,const mpres_t S, 
              mpmod_t modulus)
{
  ASSERT_NORMALIZED (S);
  if (modulus->repr == ECM_MOD_MPZ || modulus->repr == ECM_MOD_BASE2)
    {
      mpz_ui_sub (R, n, S);
      if (mpz_sgn (R) < 0)
        mpz_add (R, R, modulus->orig_modulus); /* Assumes modulus >= n */
    }
  else if (modulus->repr == ECM_MOD_MODMULN || modulus->repr == ECM_MOD_REDC)
    {
      mpz_set_ui (modulus->temp1, n);
      mpz_mul_2exp (modulus->temp1, modulus->temp1, modulus->bits);
      mpz_sub (modulus->temp1, modulus->temp1, S);
      mpz_mod (R, modulus->temp1, modulus->orig_modulus);
    }
  ASSERT_NORMALIZED (R);
}
#endif

/* R <- S1 - S2 mod modulus */
void 
mpres_sub (mpres_t R, const mpres_t S1, const mpres_t S2, mpmod_t modulus)
{
  ASSERT_NORMALIZED (S1);
  ASSERT_NORMALIZED (S2);
  mpz_sub (R, S1, S2);
  if ((modulus->repr == ECM_MOD_MODMULN || modulus->repr == ECM_MOD_REDC) &&
      ABSIZ(R) > ABSIZ(modulus->orig_modulus))
    {
      if (SIZ(R) > 0)
	mpz_sub (R, R, modulus->multiple);
      else
	mpz_add (R, R, modulus->multiple);
      /* N <= since multiple < 2^Nbits + N, now |R| < B */
    }
  ASSERT_NORMALIZED (R);
}

void 
mpres_set_z (mpres_t R, const mpz_t S, mpmod_t modulus)
{
  if (modulus->repr == ECM_MOD_MPZ || modulus->repr == ECM_MOD_BASE2)
    mpz_mod (R, S, modulus->orig_modulus);
  else if (modulus->repr == ECM_MOD_MODMULN)
    {
      mpz_mod (modulus->temp2, S, modulus->orig_modulus);
      ecm_mulredc_basecase (R, modulus->temp2, modulus->R2, modulus);
    }
  else if (modulus->repr == ECM_MOD_REDC)
    {
      mpz_mod (modulus->temp2, S, modulus->orig_modulus);
      mpz_mul (modulus->temp1, modulus->temp2, modulus->R2);
      REDC (R, modulus->temp1, modulus->temp2, modulus);
    }
  ASSERT_NORMALIZED (R);
}

/* R and S must not be modulus->temp1 */
void 
mpres_get_z (mpz_t R, const mpres_t S, mpmod_t modulus)
{
  ASSERT_NORMALIZED (S);
  if (modulus->repr == ECM_MOD_MPZ || modulus->repr == ECM_MOD_BASE2)
    {
      mpz_mod (R, S, modulus->orig_modulus);
    }
  else if (modulus->repr == ECM_MOD_MODMULN)
    {
      mpz_set (modulus->temp1, S);
      MPZ_REALLOC (R, modulus->bits / GMP_NUMB_BITS);
      ecm_redc_basecase (R, modulus->temp1, modulus);
      mpz_mod (R, R, modulus->orig_modulus); /* FIXME: can we avoid this? */
    }
  else if (modulus->repr == ECM_MOD_REDC)
    {
      REDC (R, S, modulus->temp1, modulus);
      mpz_mod (R, R, modulus->orig_modulus); /* FIXME: can we avoid this? */
    }
#ifdef DEBUG
  else
    {
      fprintf (ECM_STDERR, "mpres_get_z: Unexpected representation %d\n", 
               modulus->repr);
      exit (EXIT_FAILURE);
    }
#endif
}

/* R <- n mod modulus
   If repr==ECM_MOD_MPZ or ECM_MOD_BASE2, we convert n to
   Montgomery representation
 */
void 
mpres_set_ui (mpres_t R, const unsigned long n, mpmod_t modulus)
{
  if (modulus->repr == ECM_MOD_MPZ || modulus->repr == ECM_MOD_BASE2)
    {
      mpz_set_ui (R, n);
      mpz_mod (R, R, modulus->orig_modulus);
    }
  else if (modulus->repr == ECM_MOD_MODMULN || modulus->repr == ECM_MOD_REDC)
    {
      mpz_set_ui (modulus->temp1, n);
      mpz_mul_2exp (modulus->temp1, modulus->temp1, modulus->bits);
      mpz_mod (R, modulus->temp1, modulus->orig_modulus);
    }
  ASSERT_NORMALIZED (R);
}

/* same as previous but with signed long */
void 
mpres_set_si (mpres_t R, const long n, mpmod_t modulus)
{
  if (modulus->repr == ECM_MOD_MPZ || modulus->repr == ECM_MOD_BASE2)
    {
      mpz_set_si (R, n);
      mpz_mod (R, R, modulus->orig_modulus);
    }
  else if (modulus->repr == ECM_MOD_MODMULN || modulus->repr == ECM_MOD_REDC)
    {
      mpz_set_si (modulus->temp1, n);
      mpz_mul_2exp (modulus->temp1, modulus->temp1, modulus->bits);
      mpz_mod (R, modulus->temp1, modulus->orig_modulus);
    }
  ASSERT_NORMALIZED (R);
}

/* R <- -S mod modulus. Does not need to be efficient. */
void
mpres_neg (mpres_t R, const mpres_t S, ATTRIBUTE_UNUSED mpmod_t modulus)
{
  ASSERT_NORMALIZED (S);
  mpz_neg (R, S);
  ASSERT_NORMALIZED (R);
}

/* Returns non-zero if inversion succeeded, and zero if not */

int 
mpres_invert (mpres_t R, const mpres_t S, mpmod_t modulus)
{
#ifdef WANT_ASSERT_EXPENSIVE
  mpres_t test;
  mpz_t test_result;
  mpres_init (test, modulus);
  mpres_set (test, S, modulus);
#endif

  ASSERT_NORMALIZED (S);

  if (mpz_invert (modulus->temp2, S, modulus->orig_modulus) == 0)
    return 0;
  
  if (modulus->repr == ECM_MOD_MPZ || modulus->repr == ECM_MOD_BASE2)
    {
      mpz_set (R, modulus->temp2);
      ASSERT_NORMALIZED (R);
    }
  else if (modulus->repr == ECM_MOD_MODMULN)
    {
      ecm_mulredc_basecase (R, modulus->temp2, modulus->R3, modulus);
      ASSERT_NORMALIZED (R);
    }
  else if (modulus->repr == ECM_MOD_REDC)
    {
      MPZ_NORMALIZED (S);
      mpz_mul (modulus->temp1, modulus->temp2, modulus->R3);
      REDC (R, modulus->temp1, modulus->temp2, modulus);
      ASSERT_NORMALIZED (R);
    }
#ifdef DEBUG
  else
    {
      fprintf (ECM_STDERR, "mpres_invert: Unexpected representation %d\n", 
               modulus->repr);
      exit (EXIT_FAILURE);
    }
#endif

#ifdef WANT_ASSERT_EXPENSIVE
  mpres_mul (test, test, R, modulus);
  mpz_init (test_result);
  mpres_get_z (test_result, test, modulus);
  ASSERT_ALWAYS(mpz_cmp_ui (test_result, 1UL) == 0);
  mpz_clear (test_result);
  mpres_clear (test, modulus);
#endif

  return 1;
}

void 
mpres_gcd (mpz_t R, const mpres_t S, const mpmod_t modulus)
{
  /* In MODMULN and REDC form, M(x) = x*R with gcd(R, modulus) = 1 .
     Therefore gcd(M(x), modulus) = gcd(x, modulus) and we need not bother
     to convert out of Montgomery form. */
  ASSERT_NORMALIZED (S);
  mpz_gcd (R, S, modulus->orig_modulus);
}

#ifdef HAVE_ADDLAWS
/* Returns non-zero if the two residues are equal, 
   and zero if they are not */
int 
mpres_equal (const mpres_t S1, const mpres_t S2, mpmod_t modulus)
{
  mpz_mod (modulus->temp1, S1, modulus->orig_modulus);
  mpz_mod (modulus->temp2, S2, modulus->orig_modulus);
  return (mpz_cmp (modulus->temp1, modulus->temp2) == 0);
}
#endif

#if 0 /* those routines are not called in normal operation */

void 
mpres_out_str (FILE *fd, const unsigned int base, const mpres_t S, 
               mpmod_t modulus)
{
  mpres_get_z (modulus->temp2, S, modulus);
  mpz_out_str (fd, base, modulus->temp2);
}

/* Multiplies S1 by the one-limb integer S2, and does modulo reduction.
   The modulo reduction may imply multiplication of the residue class 
   by some constant, since we may not do the correct number of REDC 
   reduction passes and so fail to divide by the correct power of 2 for 
   Montgomery representation. The constant is the same for each call
   of this function with a given modulus, however. */
static void 
ecm_mulredc_1_basecase (mpres_t R, const mpres_t S1, const mp_limb_t S2, 
                        mpmod_t modulus)
{
  mp_ptr s1p;
  mp_size_t j, nn = modulus->bits / GMP_NUMB_BITS;

  ASSERT(ALLOC(R) >= nn);
  ASSERT(ALLOC(S1) >= nn);
  s1p = PTR(S1);
  for (j = ABSIZ(S1); j < nn; j++) 
    s1p[j] = 0;

#ifdef HAVE_NATIVE_MULREDC1_N
  if (nn < 20)
    {
      mp_ptr rp = PTR(R);
      mulredc_1(rp, S2, s1p, PTR(modulus->orig_modulus), nn,
                modulus->Nprim[0]);
      MPN_NORMALIZE (rp, nn);
      SIZ(R) = (SIZ(S1)) < 0 ? (int) -nn : (int) nn;
    }
  else
#endif
    {
      /* FIXME, we can do much better than this */
      mpz_mul_ui (modulus->temp1, S1, S2);
      mpz_mod(R, modulus->temp1, modulus->orig_modulus);
    }
}

/* Multiplies S by n and possibly divides by some constant. 
   Whether or not it divides depends on the modulus representation and
   the modulus size. */
void 
mpres_muldivbysomething_si (mpres_t R, const mpres_t S, const long n, 
			    mpmod_t modulus)
{
  ASSERT_NORMALIZED (S);
  if (modulus->repr == ECM_MOD_MODMULN && 
      modulus->bits / GMP_NUMB_BITS <= 20)
    /* FIXME: is the 20 here the same constant as in mulredc1_20?
       If so, it should be changed into a macro. */
    {
      MPZ_REALLOC (R, modulus->bits / GMP_NUMB_BITS);
      if (n < 0)
	{
	  ecm_mulredc_1_basecase (R, S, (mp_limb_t) -n, modulus);
	  mpres_neg (R, R, modulus);
	}
      else
	{
	  ecm_mulredc_1_basecase (R, S, (mp_limb_t) n, modulus);
	}
    }
  else
    {
      mpz_mul_si (modulus->temp1, S, n);
      /* This is the same for all methods: just reduce with original modulus */
      mpz_mod (R, modulus->temp1, modulus->orig_modulus);
    }
  ASSERT_NORMALIZED (R);
}

/* Returns 1 if successful, 0 if not */
static int
test_mpres_set_z_for_gcd_fix (const int maxk, mpmod_t modulus)
{
  mpres_t m, prod;
  mpz_t n;
  int k;

  mpres_init (prod, modulus);
  mpres_init (m, modulus);
  mpz_init (n);

  /* m = 1 * c, where c is a constant depending on mod reduction method */
  mpz_set_ui (n, 1);
  mpres_set_z_for_gcd (m, n, modulus);

  for (k = 0; k <= maxk; k++)
    {
      int i;
      mpres_set_ui (prod, 1, modulus);

      /* Compute prod = 1 * (1 * c)^k */
      for (i = 0; i < k; i++)
        mpres_mul (prod, prod, m, modulus);

      /* Divide prod by c^k */ 
      mpz_set_ui (n, k);
      mpres_set_z_for_gcd_fix (prod, prod, n, modulus);

      /* Result should be 1 */
      mpres_get_z (n, prod, modulus);
      if (mpz_cmp_ui (n, 1) != 0) {
        gmp_printf ("Error: test of mpres_set_z_for_gcd_fix() failed, k = %d, n = %Zd\n", k, n);
        return 0;
      }
    }

  mpz_clear (n);
  mpres_clear (m, modulus);
  mpres_clear (prod, modulus);
  return 1;
}

int
mpmod_selftest (const mpz_t n)
{
  mpres_t test1, test2;
  mpmod_t modulus;
  
  printf ("Performing self test of modular arithmetic\n");
  mpmod_init (modulus, n, 0);
  mpres_init (test1, modulus);
  mpres_init (test2, modulus);
  mpres_set_ui (test1, 2, modulus);
  mpres_set_ui (test2, 5, modulus);
  mpres_muldivbysomething_si (test1, test1, 5, modulus);
  mpres_muldivbysomething_si (test2, test2, 2, modulus);
  if (!mpres_equal (test1, test2, modulus))
   {
     printf ("mpres_muldivbysomething_si() wrong\n");
     fflush (stdout);
     abort();
   }
  mpres_clear (test1, modulus);
  mpres_clear (test2, modulus);

  if (!test_mpres_set_z_for_gcd_fix(10, modulus))
    abort();

  mpmod_clear (modulus);
  
  return 0;
}
#endif /* #if 0 */

/****************************************************/
/* mpresn: modular arithmetic based directly on mpn */
/****************************************************/

/* We use here a signed word-based redundant representation.

   In case N < B^n/16 (since for redc where we add to the absolute value of
   the residue), where n is the number of limbs of N in base B (2^32 or 2^64
   usually), we can prove there is no adjustment (adding or subtracting N),
   cf http://www.loria.fr/~zimmerma/papers/norm.pdf.

   However current branch predictors are quite good, thus we prefer to keep
   the tests and to allow any input N (instead of only N < B^n/16).
*/

/* ensure R has allocated space for at least n limbs,
   and if less than n limbs are used, pad with zeros,
   and set SIZ(R) to n if positive or -n if negative */
void
mpresn_pad (mpres_t R, mpmod_t N)
{
  mp_size_t n = ABSIZ(N->orig_modulus);
  mp_size_t rn;

  _mpz_realloc (R, n);
  rn = mpz_size (R);
  ASSERT_ALWAYS (rn <= n);
  if (rn < n)
    {
      MPN_ZERO (PTR(R) + rn, n - rn);
      SIZ(R) = SIZ(R) >= 0 ? n : -n;
    }
}

void
mpresn_unpad (mpres_t R)
{
  mp_size_t n = ABSIZ(R);

  while (n > 0 && PTR(R)[n-1] == 0)
    n--;
  SIZ(R) = SIZ(R) >= 0 ? n : -n;
}

/* R <- S1 * S1 mod N, used only for ECM_MOD_MODMULN */
void 
mpresn_sqr (mpres_t R, const mpres_t S1, mpmod_t modulus)
{
  mp_size_t n = ABSIZ(modulus->orig_modulus);

  ASSERT (SIZ(S1) == n || -SIZ(S1) == n);

  ecm_sqrredc_basecase_n (PTR(R), PTR(S1), PTR(modulus->orig_modulus),
                          n, modulus->Nprim, PTR(modulus->temp1));

  SIZ(R) = n;
}

/* R <- S1 * S2 mod N, used only for ECM_MOD_MODMULN */
void 
mpresn_mul (mpres_t R, const mpres_t S1, const mpres_t S2, mpmod_t modulus)
{
  mp_size_t n = ABSIZ(modulus->orig_modulus);

  ASSERT (SIZ(S1) == n || -SIZ(S1) == n);
  ASSERT (SIZ(S2) == n || -SIZ(S2) == n);

  ecm_mulredc_basecase_n (PTR(R), PTR(S1), PTR(S2), PTR(modulus->orig_modulus),
                          n, modulus->Nprim, PTR(modulus->temp1));

  SIZ(R) = SIZ(S1) == SIZ(S2) ? n : -n;
}

/* R <- S*m/B mod modulus where m fits in a mp_limb_t.
   Here S (w in dup_add_batch1) is the result of a subtraction,
   thus with the notations from http://www.loria.fr/~zimmerma/papers/norm.pdf
   we have S < 2 \alpha N.
   Then R < (2 \alpha N \beta + \beta N) = (2 \alpha + 1) N.
   This result R is used in an addition with u being the result of a squaring
   thus u < \alpha N, which gives a result < (3 \alpha + 1) N.
   Finally this result is used in a multiplication with another operand less
   than 2 \alpha N, thus we want:
   ((2 \alpha) (3 \alpha + 1) N^2 + \beta N)/\beta \leq \alpha N, i.e.,
   2 \alpha (3 \alpha + 1) \varepsilon + 1 \leq \alpha
   This implies \varepsilon \leq 7/2 - sqrt(3)/2 ~ 0.0359, in which case
   we can take \alpha = 2/3*sqrt(3)+1 ~ 2.1547.
   In that case no adjustment is needed in mpresn_mul_1.
   However we prefer to keep the adjustment here, to allow a larger set of
   inputs (\varepsilon \leq 1/16 = 0.0625 instead of 0.0359).
*/
void 
mpresn_mul_1 (mpres_t R, const mpres_t S, const mp_limb_t m, mpmod_t modulus)
{
  mp_ptr t1 = PTR(modulus->temp1);
  mp_ptr t2 = PTR(modulus->temp2);
  mp_size_t n = ABSIZ(modulus->orig_modulus);
  mp_limb_t q;

  ASSERT (SIZ(S) == n || -SIZ(S) == n);
  ASSERT (ALLOC(modulus->temp1) >= n+1);

#if defined(USE_ASM_REDC) && defined(HAVE_NATIVE_MULREDC1_N)
  if (n <= MULREDC_ASSEMBLY_MAX)
    mulredc_1 (PTR(R), m, PTR(S), PTR(modulus->orig_modulus), n,
               modulus->Nprim[0]);
  else
#endif
    {
      t1[n] = mpn_mul_1 (t1, PTR(S), n, m);
      q = t1[0] * modulus->Nprim[0];
      t2[n] = mpn_mul_1 (t2, PTR(modulus->orig_modulus), n, q);
#ifdef HAVE___GMPN_ADD_NC
      q = __gmpn_add_nc (PTR(R), t1 + 1, t2 + 1, n, t1[0] != 0);
#else
      q = mpn_add_n (PTR(R), t1 + 1, t2 + 1, n);
      q += mpn_add_1 (PTR(R), PTR(R), n, t1[0] != 0);
#endif
      while (q != 0)
        q -= mpn_sub_n (PTR(R), PTR(R), PTR(modulus->orig_modulus), n);
    }

  SIZ(R) = SIZ(S); /* sign is unchanged */
}

/* R <- S1 + S2 mod modulus */
/* we assume all numbers are allocated to n limbs, and unused most significant
   limbs are set to zero */
void
mpresn_add (mpres_t R, const mpres_t S1, const mpres_t S2, mpmod_t modulus)
{
  mp_ptr r = PTR(R);
  mp_ptr s1 = PTR(S1);
  mp_ptr s2 = PTR(S2);
  mp_size_t n = ABSIZ(modulus->orig_modulus);
  ATTRIBUTE_UNUSED mp_limb_t cy;

  ASSERT (SIZ(S1) == n || -SIZ(S1) == n);
  ASSERT (SIZ(S2) == n || -SIZ(S2) == n);

  if (SIZ(S1) == SIZ(S2)) /* S1 and S2 are of same sign */
    {
      cy = mpn_add_n (r, s1, s2, n);
      /* for N < B^n/16, the while loop will be never performed, which proves
         it will be performed a small number of times. In practice we
         observed up to 7 loops, but it happens rarely. */
#ifndef MPRESN_NO_ADJUSTMENT
      while (cy != 0)
        cy -= mpn_sub_n (r, r, PTR(modulus->orig_modulus), n);
#endif
      SIZ(R) = SIZ(S1);
    }
  else /* different signs */
    {
      if (mpn_cmp (s1, s2, n) >= 0)
        {
          mpn_sub_n (r, s1, s2, n); /* no borrow here */
          SIZ(R) = SIZ(S1);
        }
      else
        {
          mpn_sub_n (r, s2, s1, n); /* idem */
          SIZ(R) = SIZ(S2);
        }
    }
}

void
mpresn_sub (mpres_t R, const mpres_t S1, const mpres_t S2, mpmod_t modulus)
{
  mp_ptr r = PTR(R);
  mp_ptr s1 = PTR(S1);
  mp_ptr s2 = PTR(S2);
  mp_size_t n = ABSIZ(modulus->orig_modulus);
  ATTRIBUTE_UNUSED mp_limb_t cy;

  ASSERT (SIZ(S1) == n || -SIZ(S1) == n);
  ASSERT (SIZ(S2) == n || -SIZ(S2) == n);

  if (SIZ(S1) != SIZ(S2)) /* S1 and S2 are of different signs */
    {
      cy = mpn_add_n (r, s1, s2, n);
#ifndef MPRESN_NO_ADJUSTMENT
      while (cy != 0)
        cy -= mpn_sub_n (r, r, PTR(modulus->orig_modulus), n);
#endif
      SIZ(R) = SIZ(S1);
    }
  else /* same signs, it's a real subtraction */
    {
      if (mpn_cmp (s1, s2, n) >= 0)
        {
          mpn_sub_n (r, s1, s2, n); /* no borrow here */
          SIZ(R) = SIZ(S1);
        }
      else
        {
          mpn_sub_n (r, s2, s1, n); /* idem */
          SIZ(R) = -SIZ(S2);
        }
      
    }
}

/* (R, T) <- (S1 + S2, S1 - S2)
   Assume R differs from both S1 and S2.
 */
void
mpresn_addsub (mpres_t R, mpres_t T,
               const mpres_t S1, const mpres_t S2, mpmod_t modulus)
{
  mp_ptr r = PTR(R);
  mp_ptr t = PTR(T);
  mp_ptr s1 = PTR(S1);
  mp_ptr s2 = PTR(S2);
  mp_size_t n = ABSIZ(modulus->orig_modulus);
  ATTRIBUTE_UNUSED mp_limb_t cy;

  ASSERT (R != S1);
  ASSERT (R != S2);
  ASSERT (SIZ(S1) == n || -SIZ(S1) == n);
  ASSERT (SIZ(S2) == n || -SIZ(S2) == n);

  if (SIZ(S1) == SIZ(S2)) /* S1 and S2 are of same sign */
    {
      cy = mpn_add_n (r, s1, s2, n);
#ifndef MPRESN_NO_ADJUSTMENT
      while (cy != 0)
        cy -= mpn_sub_n (r, r, PTR(modulus->orig_modulus), n);
#endif
      SIZ(R) = SIZ(S1);
      if (mpn_cmp (s1, s2, n) >= 0)
        {
          mpn_sub_n (t, s1, s2, n); /* no borrow since {s1,n} >= {s2,n} */
          SIZ(T) = SIZ(S1);
        }
      else
        {
          mpn_sub_n (t, s2, s1, n); /* idem since {s2,n} >= {s1,n} */
          SIZ(T) = -SIZ(S2);
        }
    }
  else /* different signs */
    {
      if (mpn_cmp (s1, s2, n) >= 0)
        {
          mpn_sub_n (r, s1, s2, n); /* no borrow since {s1,n} >= {s2,n} */
          SIZ(R) = SIZ(S1);
        }
      else
        {
          mpn_sub_n (r, s2, s1, n); /* idem since {s2,n} >= {s1,n} */
          SIZ(R) = SIZ(S2);
        }
      cy = mpn_add_n (t, s1, s2, n);
#ifndef MPRESN_NO_ADJUSTMENT
      while (cy != 0)
        cy -= mpn_sub_n (t, t, PTR(modulus->orig_modulus), n);
#endif
      SIZ(T) = SIZ(S1);
    }
}
