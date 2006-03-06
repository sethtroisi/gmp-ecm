/* Modular multiplication.

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
#include "ecm-gmp.h"
#include "ecm-impl.h"

#ifdef NATIVE_REDC
  #include "asmredc.h"
#endif

FILE *ECM_STDOUT, *ECM_STDERR; /* define them here since needed in tune.c */

/* define WANT_ASSERT to check normalization of residues */
/* #define WANT_ASSERT 1 */
/* #define DEBUG */

#define ASSERT_NORMALIZED(x) ASSERT ((modulus->repr != ECM_MOD_MODMULN && \
				      modulus->repr != ECM_MOD_REDC) || \
			     mpz_size (x) <= mpz_size (modulus->orig_modulus))
#define MPZ_NORMALIZED(x)    ASSERT (PTR(x)[ABSIZ(x)-1] != 0)

#ifndef GMP_NUMB_BITS
#define GMP_NUMB_BITS __GMP_BITS_PER_MP_LIMB
#endif


void base2mod (mpres_t, mpres_t, mpres_t, mpmod_t);
void base2mod_1 (mpres_t, mpres_t, mpmod_t);
void REDC (mpres_t, mpres_t, mpz_t, mpmod_t);
void mod_mul2exp (mpz_t, unsigned int, mpmod_t);
void mod_div2exp (mpz_t, unsigned int, mpmod_t);

/* returns +/-l if n is a factor of N = 2^l +/- 1 with N <= n^threshold, 
   0 otherwise.
*/
int 
isbase2 (mpz_t n, double threshold)
{
  unsigned int k, lo; 
  int res = 0; 
  mpz_t u, w;

  MPZ_INIT (u);
  MPZ_INIT (w);
  lo = mpz_sizeinbase (n, 2) - 1; /* 2^lo <= n < 2^(lo+1) */  
  mpz_set_ui (u, 1);
  mpz_mul_2exp (u, u, 2 * lo);
  mpz_mod (w, u, n); /* 2^(2lo) mod n = -/+2^(2lo-l) if m*n = 2^l+/-1 */
  if (mpz_cmp_ui (w, 1) == 0) /* if 2^(2lo) mod n = 1, then n divides 2^lo+1, 
				 since n has lo+1 bits. */
    {
      mpz_clear (w);
      mpz_clear (u);
      return lo;
    }
  k = mpz_sizeinbase (w, 2) - 1;
  /* if w = 2^k then n divides 2^(2*lo-k)-1 */
  mpz_set_ui (u, 1);
  mpz_mul_2exp (u, u, k);
  if (mpz_cmp(w, u) == 0) 
    res = k - 2 * lo;
  else /* if w = -2^k then n divides 2^(2*lo-k)+1 */
    {
      mpz_neg (w, w);
      mpz_mod (w, w, n);
      k = mpz_sizeinbase (w, 2) - 1;
      mpz_set_ui (u, 1);
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
void
base2mod (mpres_t R, mpres_t S, mpres_t t, mpmod_t modulus)
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

/* Same, but source and result in same variable */
void
base2mod_1 (mpres_t RS, mpres_t t, mpmod_t modulus)
{
  unsigned long absbits = abs (modulus->bits);

  ASSERT (RS != t);
  while (mpz_sizeinbase (RS, 2) > absbits)
    {
      mpz_tdiv_q_2exp (t, RS, absbits);
      mpz_tdiv_r_2exp (RS, RS, absbits); /* Just a truncate */
      if (modulus->bits < 0)
        mpz_add (RS, RS, t);
      else
        mpz_sub (RS, RS, t);
    }
}

#ifdef HAVE___GMPN_MUL_FFT
/* Fermat-mod */
static void
base2mod_2 (mpres_t RS, mp_size_t n, mpz_t modulus)
{
  mp_size_t s;

  s = ABSIZ(RS);
  if (s > n)
    {
      if (s == n + 1)
        {
          mp_ptr rp = PTR(RS);

          if ((rp[n] = mpn_sub_1 (rp, rp, n, rp[n])))
            rp[n] = mpn_add_1 (rp, rp, n, rp[n]);
          MPN_NORMALIZE(rp, s);
          SIZ(RS) = (SIZ(RS) > 0) ? (int) s : (int) -s;
        }
      else /* should happen rarely */
        mpz_mod (RS, RS, modulus);
    }
}
#endif

/* subquadratic REDC, at mpn level.
   {orig,n} is the original modulus.
   {aux,n} is the auxiliary modulus.
   Requires ABSIZ(x) = 2n and ABSIZ(orig_modulus)=ABSIZ(aux_modulus)=n.
 */
static void
ecm_redc_n (mp_ptr rp, mp_srcptr xp, mp_srcptr orig, mp_srcptr aux, mp_size_t n)
{
  mp_ptr tp, up;
  mp_size_t nn = n + n;
  mp_limb_t cy;
  TMP_DECL(marker);

  TMP_MARK(marker);
  up = TMP_ALLOC_LIMBS(nn + nn);
  ecm_mul_lo_n (up, xp, aux, n);
  tp = up + nn;
  mpn_mul_n (tp, up, orig, n);
  /* add {x, 2n} and {tp, 2n}. We know that {tp, n} + {xp, n} will give
     either 0, or a carry out. If xp[n-1] <> 0, then there is a carry. */
#ifdef HAVE___GMPN_ADD_NC
  cy = __gmpn_add_nc (rp, tp + n, xp + n, n, (mp_limb_t) ((xp[n - 1]) ? 1 : 0));
#else
  cy = mpn_add_n (rp, tp + n, xp + n, n);
  cy += mpn_add_1 (rp, rp, n, (mp_limb_t) ((xp[n - 1]) ? 1 : 0));
#endif
  if (cy || mpn_cmp (rp, orig, n) > 0)
    cy -= mpn_sub_n (rp, rp, orig, n);
  /* ASSERT ((cy == 0) && (mpn_cmp (rp, orig, n) < 0)); */
  TMP_FREE(marker);
}

/* REDC. x and t must not be identical, t has limb growth */
/* subquadratic REDC, at mpz level */
void 
REDC (mpres_t r, mpres_t x, mpz_t t, mpmod_t modulus)
{
  mp_size_t n = modulus->bits / GMP_NUMB_BITS;

  ASSERT (ABSIZ(x) <= 2 * n);
  if (ABSIZ(x) == 2 * n)
    {
      mp_ptr rp;
      MPZ_REALLOC (r, n);
      rp = PTR(r);
      ecm_redc_n (rp, PTR(x), PTR(modulus->orig_modulus), 
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
      mpz_mul (t, t, modulus->orig_modulus);
      mpz_add (t, t, x);
      mpz_tdiv_q_2exp (r, t, modulus->bits);  /* r = (x + m*N) / R */
      if (ABSIZ (r) > n)
	mpz_sub (r, r, modulus->multiple);
    }
  ASSERT (ABSIZ(r) <= n);
}

/* multiplies c by R^k modulo n where R=2^mp_bits_per_limb 
   n is supposed odd. Does not need to be efficient. */
void 
mod_mul2exp (mpz_t c, unsigned int k, mpmod_t modulus)
{
  mpz_mul_2exp (modulus->temp1, c, k * __GMP_BITS_PER_MP_LIMB);
  mpz_mod (c, modulus->temp1, modulus->orig_modulus);
}

/* divides c by R^k modulo n where R=2^mp_bits_per_limb
   n is supposed odd. Does not need to be efficient. */
void 
mod_div2exp (mpz_t c, unsigned int k, mpmod_t modulus)
{
  mpz_set_ui (modulus->temp2, 1);
  mpz_mul_2exp (modulus->temp1, modulus->temp2, k * __GMP_BITS_PER_MP_LIMB);
  mpz_invert (modulus->temp2, modulus->temp1, modulus->orig_modulus); 
    /* temp2 = 2^(-k) (mod n) */
  mpz_mul (modulus->temp1, modulus->temp2, c);
  mpz_mod (c, modulus->temp1, modulus->orig_modulus);
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
  mp_limb_t cy;
  mp_size_t j, nn = modulus->bits / __GMP_BITS_PER_MP_LIMB;

  ASSERT(ABSIZ(c) <= 2 * nn);
  ASSERT(ALLOC(r) >= nn);
  cp = PTR(c);
  rp = PTR(r);
  np = PTR(modulus->orig_modulus);
  for (j = ABSIZ(c); j < 2 * nn; j++) 
    cp[j] = 0;
#ifndef NATIVE_REDC
  for (j = 0; j < nn; j++)
    {
      cp[0] = mpn_addmul_1 (cp, np, nn, cp[0] * modulus->Nprim);
      cp++;
    }
  /* add vector of carries and shift */
  cy = mpn_add_n (rp, cp, cp - nn, nn);
#else
  ecm_redc3 (cp, np, nn, modulus->Nprim);
  /* add vector of carries and shift */
  cy = mpn_add_n (rp, cp + nn, cp, nn);
#endif
  /* the result of Montgomery's REDC is less than 2^Nbits + N,
     thus at most one correction is enough */
  if (cy != 0)
    mpn_sub_n (rp, rp, np, nn); /* a borrow should always occur here */
  MPN_NORMALIZE (rp, nn);
  SIZ(r) = SIZ(c) < 0 ? (int) -nn : (int) nn;
}

#ifdef NATIVE_REDC
static mp_limb_t
mulredc(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y,
    const mp_limb_t *m, mp_size_t N, mp_limb_t invm, mp_limb_t *tmp)
{
  mp_limb_t cy;

  switch (N) {
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
    {
      mpn_mul_n(tmp, x, y, N);
      ecm_redc3(tmp, m, N, invm);
      cy = mpn_add_n (z, tmp + N, tmp, N);
    }
  }
  return cy;
}
#endif


/* 
 * Same as previous, but combined with mul (if in asm)
 */
static void 
ecm_mulredc_basecase (mpres_t R, mpres_t S1, mpres_t S2, mpmod_t modulus)
{
#ifndef NATIVE_REDC
  mpz_mul (modulus->temp1, S1, S2);
  ecm_redc_basecase(R, modulus->temp1, modulus);
#else
  mp_ptr rp;
  mp_ptr s1p, s2p;
  mp_srcptr np;
  mp_limb_t cy;
  mp_size_t j, nn = modulus->bits / __GMP_BITS_PER_MP_LIMB;

  ASSERT(ALLOC(R) >= nn);
  ASSERT(ALLOC(S1) >= nn);
  ASSERT(ALLOC(S2) >= nn);
  rp = PTR(R);
  s1p = PTR(S1);
  s2p = PTR(S2);
  np = PTR(modulus->orig_modulus);
  for (j = ABSIZ(S1); j < nn; j++) 
    s1p[j] = 0;
  for (j = ABSIZ(S2); j < nn; j++) 
    s2p[j] = 0;

  cy = mulredc(rp, s1p, s2p, np, nn, modulus->Nprim, PTR(modulus->temp1));

  /* the result of Montgomery's REDC is less than 2^Nbits + N,
     thus at most one correction is enough */
  if (cy != 0)
    mpn_sub_n (rp, rp, np, nn); /* a borrow should always occur here */
  MPN_NORMALIZE (rp, nn);
  SIZ(R) = (SIZ(S1)*SIZ(S2)) < 0 ? (int) -nn : (int) nn;
#endif
}


/* don't use base2 if repr == -1, i.e. -nobase2 */
void 
mpmod_init (mpmod_t modulus, mpz_t N, int repr)
{
  int base2;
  
  if ((repr != -1) && (base2 = isbase2 (N, BASE2_THRESHOLD)))
    {
      int r;
      outputf (OUTPUT_VERBOSE,
	       "Using special division for factor of 2^%d%c1\n",
	       abs (base2), (base2 < 0) ? '-' : '+');
      r = mpmod_init_BASE2 (modulus, base2, N);
      ASSERT (r == 0); /* error should not happen if isbase2 is correct */
    }
  else if (mpz_size (N) < MPZMOD_THRESHOLD)
    {
      outputf (OUTPUT_VERBOSE, "Using MODMULN\n");
      mpmod_init_MODMULN (modulus, N);
    }
  else if (mpz_size (N) < REDC_THRESHOLD)
    {
      outputf (OUTPUT_VERBOSE, "Using mpz_mod\n");
      mpmod_init_MPZ (modulus, N);
    }
  else
    {
      outputf (OUTPUT_VERBOSE, "Using REDC\n");
      mpmod_init_REDC (modulus, N);
    }
  
  return;
}

void 
mpmod_init_MPZ (mpmod_t modulus, mpz_t N)
{
  int Nbits;
  
  mpz_init_set (modulus->orig_modulus, N);
  modulus->repr = ECM_MOD_MPZ;
  
  Nbits = mpz_size (N) * __GMP_BITS_PER_MP_LIMB; /* Number of bits, rounded
                                                    up to full limb */
  MPZ_INIT2 (modulus->temp1, 2 * Nbits + __GMP_BITS_PER_MP_LIMB);
  MPZ_INIT2 (modulus->temp2, Nbits);
  
  return;
}

int 
mpmod_init_BASE2 (mpmod_t modulus, int base2, mpz_t N)
{
  int Nbits;
  
  mpz_init_set (modulus->orig_modulus, N);
  modulus->repr = ECM_MOD_BASE2;
  modulus->bits = base2;

  Nbits = mpz_size (N) * __GMP_BITS_PER_MP_LIMB; /* Number of bits, rounded
                                                    up to full limb */
  MPZ_INIT2 (modulus->temp1, 2 * Nbits + __GMP_BITS_PER_MP_LIMB);
  MPZ_INIT2 (modulus->temp2, Nbits);
  
  mpz_set_ui (modulus->temp1, 1);
  mpz_mul_2exp (modulus->temp1, modulus->temp1, abs (base2));
  if (base2 < 0)
    mpz_sub_ui (modulus->temp1, modulus->temp1, 1);
  else
    mpz_add_ui (modulus->temp1, modulus->temp1, 1);
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
#if defined(HAVE_GWNUM) && !defined (TUNE)
          if (modulus->Fermat >= GWTHRESHOLD)
            Fgwinit (modulus->Fermat);
#endif
        }
    }
  
  return 0;
}

void
mpmod_init_MODMULN (mpmod_t modulus, mpz_t N)
{
  int Nbits;

  MEMORY_TAG;
  mpz_init_set (modulus->orig_modulus, N);
  MEMORY_UNTAG;
  
  modulus->repr = ECM_MOD_MODMULN;
  Nbits = mpz_size (N) * __GMP_BITS_PER_MP_LIMB; /* Number of bits, rounded
                                                    up to full limb */
  modulus->bits = Nbits;

  MPZ_INIT2 (modulus->temp1, 2 * Nbits + __GMP_BITS_PER_MP_LIMB);
  MPZ_INIT2 (modulus->temp2, Nbits);

  mpz_set_ui (modulus->temp1, 1);
  mpz_mul_2exp (modulus->temp1, modulus->temp1, __GMP_BITS_PER_MP_LIMB);
  mpz_tdiv_r_2exp (modulus->temp2, modulus->orig_modulus, 
                   __GMP_BITS_PER_MP_LIMB);
  mpz_invert (modulus->temp2, modulus->temp2, modulus->temp1);
    /* Now temp2 = 1/n (mod 2^bits_per_limb) */
  mpz_sub (modulus->temp2, modulus->temp1, modulus->temp2);
  modulus->Nprim = mpz_getlimbn (modulus->temp2, 0);
    /* Now Nprim = -1/n (mod 2^bits_per_limb) */

  MPZ_INIT (modulus->R2);
  mpz_set_ui (modulus->temp1, 1);
  mpz_mul_2exp (modulus->temp1, modulus->temp1, 2 * Nbits);
  mpz_mod (modulus->R2, modulus->temp1, modulus->orig_modulus);
  /* Now R2 = (2^bits)^2 (mod N) */
  
  MPZ_INIT (modulus->R3);
  mpz_mul_2exp (modulus->temp1, modulus->R2, Nbits);
  mpz_mod (modulus->R3, modulus->temp1, modulus->orig_modulus);
  /* Now R3 = (2^bits)^3 (mod N) */

  MPZ_INIT (modulus->multiple);
  mpz_set_ui (modulus->temp1, 1);
  mpz_mul_2exp (modulus->temp1, modulus->temp1, Nbits);
  /* compute ceil(2^bits / N) */
  mpz_cdiv_q (modulus->temp1, modulus->temp1, modulus->orig_modulus);
  mpz_mul (modulus->multiple, modulus->temp1, modulus->orig_modulus);
  /* Now multiple is the smallest multiple of N >= 2^bits */
}

void 
mpmod_init_REDC (mpmod_t modulus, mpz_t N)
{
  mp_size_t n;
  int Nbits;
  
  mpz_init_set (modulus->orig_modulus, N);
  
  n = mpz_size (N);
  modulus->repr = ECM_MOD_REDC;
  Nbits = n * __GMP_BITS_PER_MP_LIMB; /* Number of bits, rounded
                                                    up to full limb */
  modulus->bits = Nbits;
  
  mpz_init2 (modulus->temp1, 2 * Nbits + __GMP_BITS_PER_MP_LIMB);
  mpz_init2 (modulus->temp2, Nbits);
  MPZ_INIT (modulus->aux_modulus);

  mpz_set_ui (modulus->temp1, 1);
  mpz_mul_2exp (modulus->temp1, modulus->temp1, Nbits);
  /* since we directly check even modulus in ecm/pm1/pp1,
     N is odd here, thus 1/N mod 2^Nbits always exist */
  mpz_invert (modulus->aux_modulus, N, modulus->temp1);

  mpz_sub (modulus->aux_modulus, modulus->temp1, modulus->aux_modulus);
  /* ensure aux_modulus has n allocated limbs, for ecm_redc_n */
  if (ABSIZ(modulus->aux_modulus) < n)
    {
      /* WARNING: _mpz_realloc does not keep the value!!! */
      _mpz_realloc (modulus->aux_modulus, n);
      MPN_ZERO (PTR(modulus->aux_modulus) + ABSIZ(modulus->aux_modulus),
		n - ABSIZ(modulus->aux_modulus));
    }

  MPZ_INIT (modulus->R2);
  mpz_set_ui (modulus->temp1, 1);
  mpz_mul_2exp (modulus->temp1, modulus->temp1, 2 * Nbits);
  mpz_mod (modulus->R2, modulus->temp1, modulus->orig_modulus);
  /* Now R2 = (2^bits)^2 (mod N) */
  
  MPZ_INIT (modulus->R3);
  mpz_mul_2exp (modulus->temp1, modulus->R2, Nbits);
  mpz_mod (modulus->R3, modulus->temp1, modulus->orig_modulus);
  /* Now R3 = (2^bits)^3 (mod N) */
  
  MPZ_INIT (modulus->multiple);
  mpz_set_ui (modulus->temp1, 1);
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
  if (modulus->repr == ECM_MOD_MODMULN || modulus->repr == ECM_MOD_REDC)
    {
      mpz_clear (modulus->R2);
      mpz_clear (modulus->R3);
      mpz_clear (modulus->multiple);
      if (modulus->repr == ECM_MOD_REDC)
        mpz_clear (modulus->aux_modulus);
    }
#if defined(HAVE_GWNUM) && !defined (TUNE)
  if (modulus->Fermat >= GWTHRESHOLD)
    Fgwclear ();
#endif
  
  return;
}

void 
mpmod_pausegw (ATTRIBUTE_UNUSED mpmod_t modulus)
{
#if defined(HAVE_GWNUM) && !defined (TUNE)
  if (modulus->Fermat >= GWTHRESHOLD)
    Fgwclear ();
#endif
}

void 
mpmod_contgw (ATTRIBUTE_UNUSED mpmod_t modulus)
{
#if defined(HAVE_GWNUM) && !defined (TUNE)
  if (modulus->Fermat >= GWTHRESHOLD)
    Fgwinit (modulus->Fermat);
#endif
}

void 
mpres_init (mpres_t R, mpmod_t modulus)
{
  /* use mpz_sizeinbase since modulus->bits may not be initialized yet */
  mpz_init2 (R, mpz_sizeinbase (modulus->orig_modulus, 2) + GMP_NUMB_BITS);
}

/* realloc R so that it has at least the same number of limbs as modulus */
void
mpres_realloc (mpres_t R, mpmod_t modulus)
{
  if (modulus->repr == ECM_MOD_MODMULN)
    MPZ_REALLOC (R, modulus->bits / GMP_NUMB_BITS);
}

/* R <- BASE^EXP mod modulus.
   Assume EXP >= 0.
 */
void 
mpres_pow (mpres_t R, mpres_t BASE, mpz_t EXP, mpmod_t modulus)
{
  ASSERT_NORMALIZED (BASE);
  ASSERT (mpz_sgn (EXP) >= 0);
  if (modulus->repr == ECM_MOD_MPZ)
    {
      mpz_powm (R, BASE, EXP, modulus->orig_modulus);
    }
  else if (modulus->repr == ECM_MOD_BASE2 || modulus->repr == ECM_MOD_MODMULN ||
           modulus->repr == ECM_MOD_REDC)
    {
      unsigned int expidx;
      mp_limb_t bitmask, expbits;

      /* case EXP=0 */
      if (mpz_cmp_ui (EXP, 0) == 0)
        {
          mpres_set_ui (R, 1, modulus); /* set result to 1 */
          ASSERT_NORMALIZED (R);
          return;
        }

      expidx = mpz_size (EXP) - 1;         /* point at most significant limb */
      expbits = mpz_getlimbn (EXP, expidx); /* get most significant limb */
      bitmask = ((mp_limb_t) 1) << (GMP_NUMB_BITS - 1);

      while ((bitmask & expbits) == 0)
        {
          bitmask >>= 1;
          if (bitmask == 0)                 /* no set bits in this limb */
            {
              if (expidx == 0)              /* no more limbs -> exp was 0 */
                {
                  mpres_set_ui (R, 1, modulus); /* set result to 1 */
		  ASSERT_NORMALIZED (R);
                  return;
                }
              expidx --;
              expbits = mpz_getlimbn (EXP, expidx);
              bitmask = (mp_limb_t) 1 << (GMP_NUMB_BITS - 1);
            }
        }

    /* here the most significant limb with any set bits is in expbits, */
    /* bitmask is set to mask in the msb of expbits */

      mpz_set (modulus->temp2, BASE);
      bitmask >>= 1;

      while (1) 
        {
          for ( ; bitmask != 0; bitmask >>= 1) 
            {
              mpz_mul (modulus->temp1, modulus->temp2, modulus->temp2); /* r = r^2 */

              if (modulus->repr == ECM_MOD_BASE2)
                base2mod (modulus->temp2 , modulus->temp1, modulus->temp1, modulus);
              else if (modulus->repr == ECM_MOD_MODMULN)
                {
                  ecm_redc_basecase (modulus->temp2, modulus->temp1, modulus);
                }
              else
                REDC (modulus->temp2, modulus->temp1, modulus->temp2, modulus);

              if (expbits & bitmask)
                { 
                  mpz_mul (modulus->temp1, modulus->temp2, BASE);
                  if (modulus->repr == ECM_MOD_BASE2)
                    base2mod (modulus->temp2, modulus->temp1, modulus->temp1, modulus);
                  else if (modulus->repr == ECM_MOD_MODMULN)
                    {
                      ecm_redc_basecase (modulus->temp2, modulus->temp1, modulus);
                    }
                  else
                    REDC (modulus->temp2, modulus->temp1, modulus->temp2, modulus);
                }
            }
          if (expidx == 0)		/* if we just processed the least */
            break;			/* significant limb, we are done */
          expidx --;
          expbits = mpz_getlimbn (EXP, expidx);
          bitmask = (mp_limb_t) 1 << (GMP_NUMB_BITS - 1);
        }
      mpz_set (R, modulus->temp2); /* TODO: isn't it possible to use R instead
				      of modulus->temp2 above to avoid this
				      copy? */
    } /* if (modulus->repr == ECM_MOD_BASE2 || ... ) */
  ASSERT_NORMALIZED (R);
}


/* Returns 1 if S == 0 (mod modulus), 0 otherwise */

int
mpres_is_zero (mpres_t S, mpmod_t modulus)
{
  mpz_mod (modulus->temp1, S, modulus->orig_modulus);
  /* For all currently implemented representations, a zero residue has zero
     integer representation */
  return (mpz_sgn (modulus->temp1) == 0) ? 1 : 0;
}

/* R <- BASE^EXP mod modulus */ 
void 
mpres_ui_pow (mpres_t R, unsigned int BASE, mpres_t EXP, mpmod_t modulus)
{
  if (modulus->repr == ECM_MOD_MPZ)
    {
      mpz_set_ui (modulus->temp1, BASE);
      mpz_powm (R, modulus->temp1, EXP, modulus->orig_modulus);
    }
  else if (modulus->repr == ECM_MOD_BASE2 || modulus->repr == ECM_MOD_MODMULN ||
           modulus->repr == ECM_MOD_REDC)
    {
      unsigned int expidx;
      mp_limb_t bitmask, expbits;

      expidx = mpz_size (EXP) -1;           /* point at most significant limb */
      expbits = mpz_getlimbn (EXP, expidx); /* get most significant limb */
      bitmask = (mp_limb_t) 1 << (GMP_NUMB_BITS - 1);

      while ((bitmask & expbits) == 0)
        {
          bitmask >>= 1;
          if (bitmask == 0)                 /* no set bits in this limb */
            {
              if (expidx == 0)              /* no more limbs -> exp was 0 */
                {
                  mpres_set_ui (R, 1, modulus); /* set result to 1 */
		  ASSERT_NORMALIZED (R);
                  return;
                }
              expidx --;
              expbits = mpz_getlimbn (EXP, expidx);
              bitmask = (mp_limb_t) 1 << (GMP_NUMB_BITS - 1);
            }
        }

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
              mpz_mul (modulus->temp1, modulus->temp2, modulus->temp2); /* r = r^2 */

              if (modulus->repr == ECM_MOD_BASE2)
                base2mod (modulus->temp2 , modulus->temp1, modulus->temp1, modulus);
              else if (modulus->repr == ECM_MOD_MODMULN)
                {
                  ecm_redc_basecase (modulus->temp2, modulus->temp1, modulus);
                }
              else
                REDC (modulus->temp2, modulus->temp1, modulus->temp2, modulus);

              if (expbits & bitmask)
                {
                  if (BASE == 2)
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
      mpz_set (R, modulus->temp2); /* TODO: use R instead of modulus->temp2
				      above to avoid this copy? */
    } /* if (modulus->repr == ECM_MOD_BASE2 || ... ) */
  ASSERT_NORMALIZED (R);
}

void 
mpres_mul (mpres_t R, mpres_t S1, mpres_t S2, mpmod_t modulus)
{
  ASSERT_NORMALIZED (S1);
  ASSERT_NORMALIZED (S2);

#if defined(HAVE_GWNUM) && !defined (TUNE)
  if (modulus->repr == ECM_MOD_BASE2 && modulus->Fermat >= GWTHRESHOLD)
    {
      base2mod_1 (S1, modulus->temp1, modulus);
      base2mod_1 (S2, modulus->temp1, modulus);

#ifdef DEBUG
      mpz_mul (modulus->temp1, S1, S2);
      base2mod_1 (modulus->temp1, modulus->temp2, modulus);
      mpz_mod (modulus->temp2, modulus->temp1, modulus->orig_modulus);
#endif

      ASSERT (mpz_sizeinbase (S1, 2) <= (unsigned) abs(modulus->bits));
      ASSERT (mpz_sizeinbase (S2, 2) <= (unsigned) abs(modulus->bits));
      Fgwmul (R, S1, S2);

#ifdef DEBUG
      mpz_mod (modulus->temp1, R, modulus->orig_modulus);
      if (mpz_cmp (modulus->temp1, modulus->temp2) != 0)
        {
          fprintf (stderr, "mpres_mul: results of gwmul and mpz_mul differ\n");
          gmp_fprintf (stderr, "GMP result   : %Zd\nGWNUM result : %Zd\n", 
                       modulus->temp2, modulus->temp1);
        }
#endif

      return;
    }
#elif defined(HAVE___GMPN_MUL_FFT)
  if (modulus->repr == ECM_MOD_BASE2 && modulus->Fermat >= 32768)
    {
      mp_size_t n = modulus->Fermat / __GMP_BITS_PER_MP_LIMB;
      unsigned long k;
      
      MPZ_REALLOC (R, n + 1);
      k = mpn_fft_best_k (n, S1 == S2);
      ASSERT(mpn_fft_next_size (n, k) == n);
      base2mod_2 (S1, n, modulus->orig_modulus);
      base2mod_2 (S2, n, modulus->orig_modulus);
      mpn_mul_fft (PTR(R), n, PTR(S1), ABSIZ(S1), PTR(S2), ABSIZ(S2), k);
      n ++;
      MPN_NORMALIZE(PTR(R), n);
      SIZ(R) = ((SIZ(S1) ^ SIZ(S2)) >= 0) ? (int) n : (int) -n;
      return;
    }
#endif

  if (modulus->repr != ECM_MOD_MODMULN)
    mpz_mul (modulus->temp1, S1, S2);

  switch (modulus->repr)
    {
    case ECM_MOD_BASE2:
      base2mod (R, modulus->temp1, modulus->temp1, modulus);
      break;
    case ECM_MOD_MODMULN:
      ecm_mulredc_basecase (R, S1, S2, modulus);
      break;
    case ECM_MOD_REDC:
      REDC (R, modulus->temp1, modulus->temp2, modulus);
      break;
    default:
      mpz_mod (R, modulus->temp1, modulus->orig_modulus);
      break;
    }
  ASSERT_NORMALIZED (R);
}

void 
mpres_mul_ui (mpres_t R, mpres_t S, unsigned int n, mpmod_t modulus)
{
  ASSERT_NORMALIZED (S);
  mpz_mul_ui (modulus->temp1, S, n);
  /* This is the same for all methods: just reduce with original modulus */
  mpz_mod (R, modulus->temp1, modulus->orig_modulus);
  ASSERT_NORMALIZED (R);
}

/* R <- S / 2^n mod modulus. Does not need to be fast. */
void 
mpres_div_2exp (mpres_t R, mpres_t S, unsigned int n, mpmod_t modulus)
{
  ASSERT_NORMALIZED (S);
  if (n == 0)
    {
      mpres_set (R, S, modulus);
      ASSERT_NORMALIZED (R);
      return;
    }

    if (mpz_odd_p (S))
      {
        mpz_add (R, S, modulus->orig_modulus);
        mpz_tdiv_q_2exp (R, R, 1);
      }
    else
      mpz_tdiv_q_2exp (R, S, 1);

    for ( ; n > 1; n--)
      if (mpz_odd_p (R))
        {
          mpz_add (R, R, modulus->orig_modulus);
          mpz_tdiv_q_2exp (R, R, 1);
        }
      else
        mpz_tdiv_q_2exp (R, R, 1);

    ASSERT_NORMALIZED (R);
}

void
mpres_add_ui (mpres_t R, mpres_t S, unsigned int n, mpmod_t modulus)
{
  ASSERT_NORMALIZED (S);
  if (modulus->repr == ECM_MOD_MPZ || modulus->repr == ECM_MOD_BASE2)
    {
      mpz_add_ui (R, S, n);
      if (mpz_cmp (R, modulus->orig_modulus) > 0)
        mpz_sub (R, R, modulus->orig_modulus);
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
mpres_add (mpres_t R, mpres_t S1, mpres_t S2, mpmod_t modulus)
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

void
mpres_sub_ui (mpres_t R, mpres_t S, unsigned int n, mpmod_t modulus)
{
  ASSERT_NORMALIZED (S);
  if (modulus->repr == ECM_MOD_MPZ || modulus->repr == ECM_MOD_BASE2)
    {
      mpz_sub_ui (R, S, n);
      if (mpz_sgn (R) < 0)
        mpz_add (R, R, modulus->orig_modulus);
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

/* R <- S1 - S2 mod modulus */
void 
mpres_sub (mpres_t R, mpres_t S1, mpres_t S2, mpmod_t modulus)
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
mpres_set_z (mpres_t R, mpz_t S, mpmod_t modulus)
{
  if (modulus->repr == ECM_MOD_MPZ || modulus->repr == ECM_MOD_BASE2)
    {
      mpz_mod (R, S, modulus->orig_modulus);
    }
  else if (modulus->repr == ECM_MOD_MODMULN)
    {
      mpz_mod (modulus->temp2, S, modulus->orig_modulus);
      mpz_mul (modulus->temp1, modulus->temp2, modulus->R2);
      ecm_redc_basecase (R, modulus->temp1, modulus);
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
mpres_get_z (mpz_t R, mpres_t S, mpmod_t modulus)
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

void 
mpres_set_ui (mpres_t R, unsigned int n, mpmod_t modulus)
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

/* R <- -S mod modulus. Does not need to be efficient. */
void
mpres_neg (mpres_t R, mpres_t S, ATTRIBUTE_UNUSED mpmod_t modulus)
{
  ASSERT_NORMALIZED (S);
  mpz_neg (R, S);
  ASSERT_NORMALIZED (R);
}

int 
mpres_invert (mpres_t R, mpres_t S, mpmod_t modulus)
{
  ASSERT_NORMALIZED (S);
  if (modulus->repr == ECM_MOD_MPZ || modulus->repr == ECM_MOD_BASE2)
    {
      int res = mpz_invert (R, S, modulus->orig_modulus);
      ASSERT_NORMALIZED (R);
      return res;
    }
  else if (modulus->repr == ECM_MOD_MODMULN)
    {
      if (mpz_invert (modulus->temp2, S, modulus->orig_modulus))
        {
          mpz_mul (modulus->temp1, modulus->temp2, modulus->R3);
          ecm_redc_basecase (R, modulus->temp1, modulus);
	  ASSERT_NORMALIZED (R);
          return 1;
        }
      else
        return 0;
    }
  else if (modulus->repr == ECM_MOD_REDC)
    {
      MPZ_NORMALIZED (S);
      if (mpz_invert (modulus->temp2, S, modulus->orig_modulus))
        {
          mpz_mul (modulus->temp1, modulus->temp2, modulus->R3);
          REDC (R, modulus->temp1, modulus->temp2, modulus);
	  ASSERT_NORMALIZED (R);
          return 1;
        }
      else
        return 0;
    }
#ifdef DEBUG
  else
    {
      fprintf (ECM_STDERR, "mpres_invert: Unexpected representation %d\n", 
               modulus->repr);
      exit (EXIT_FAILURE);
    }
#endif
  return 0;
}

void 
mpres_gcd (mpz_t R, mpres_t S, mpmod_t modulus)
{
  /* In MODMULN and REDC form, M(x) = x*R with gcd(R, modulus) = 1 .
     Therefore gcd(M(x), modulus) = gcd(x, modulus) and we need not bother
     to convert out of Montgomery form. */
  ASSERT_NORMALIZED (S);
  mpz_gcd (R, S, modulus->orig_modulus);
}

void 
mpres_out_str (FILE *fd, unsigned int base, mpres_t S, mpmod_t modulus)
{
  mpres_get_z (modulus->temp2, S, modulus);
  mpz_out_str (fd, base, modulus->temp2);
}
