/* Modular multiplication.

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

#include <stdlib.h>
#include <stdio.h>
#include "gmp.h"
#include "ecm.h"
#ifdef WANT_GMP_IMPL
#include "gmp-impl.h"
#else
#include "ecm-gmp.h"
#endif /* WANT_GMP_IMPL */

#ifndef MOD_PLAIN_TO_REDC_THRESHOLD
#define MOD_PLAIN_TO_REDC_THRESHOLD 20000
#endif

#define FULL_REDUCTION
/* #define DEBUG */

void base2mod (mpres_t, mpres_t, mpres_t, mpmod_t);
void REDC (mpres_t, mpres_t, mpz_t, mpmod_t);
void mpn_REDC (mp_ptr, mp_srcptr, mp_srcptr, mp_srcptr, mp_size_t);
void mod_mul2exp (mpz_t, unsigned int, mpmod_t);
void mod_div2exp (mpz_t, unsigned int, mpmod_t);
void mpz_mod_n (mpz_t, mpz_t, mpmod_t);

/* returns +/-l if n is a factor of N = 2^l +/- 1 with N <= n^threshold, 
   0 otherwise.
*/
int 
isbase2 (mpz_t n, double threshold)
{
  unsigned int k, lo; 
  int res = 0; 
  mpz_t u, w;

  mpz_init (u);
  mpz_init (w);
  lo = mpz_sizeinbase (n, 2) - 1; /* 2^lo <= n < 2^(lo+1) */  
  mpz_set_ui (u, 1);
  mpz_mul_2exp (u, u, 2 * lo);
  mpz_mod (w, u, n); /* 2^(2lo) mod n = -/+2^(2lo-l) if m*n = 2^l+/-1 */
  if (mpz_cmp_ui (w, 1) == 0) /* if 2^(2lo) mod n = 1, then n divides 2^lo+1, 
				 since n has lo+1 bits. */
    return lo;
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
  if (abs (res) > (int)(threshold * lo)) 
    res = 0;
  return res;
}

/* Do base-2 reduction. R must not equal S or t. */
void
base2mod (mpres_t R, mpres_t S, mpres_t t, mpmod_t modulus)
{
  mpz_tdiv_q_2exp (R, S, abs (modulus->bits));
  mpz_tdiv_r_2exp (t, S, abs (modulus->bits));
  if (modulus->bits < 0)
    mpz_add (R, R, t);
  else
    mpz_sub (R, t, R);
   mpz_mod (R, R, modulus->orig_modulus);
}

/* REDC. x and t must not be identical, t has limb growth */
/* subquadratic REDC, at mpz level */
void 
REDC (mpres_t r, mpres_t x, mpz_t t, mpmod_t modulus)
{
  mp_size_t n = modulus->bits / GMP_NUMB_BITS;

  if (ABSIZ(x) == 2 * n)
    {
      if (ALLOC(r) < n)
	_mpz_realloc (r, n);
      mpn_REDC (PTR(r), PTR(x), PTR(modulus->orig_modulus), 
		PTR(modulus->aux_modulus), n);
      SIZ(r) = SIZ(x) / 2;
      /* warning: r may not be normalized, i.e. high limbs may be zero */
    }
  else
    {
      mpz_tdiv_r_2exp (t, x, modulus->bits);
      mpz_mul (t, t, modulus->aux_modulus);
      mpz_tdiv_r_2exp (t, t, modulus->bits);  /* t = (x % R) * 1/N (mod R) */
      mpz_mul (t, t, modulus->orig_modulus);
      mpz_add (t, t, x);
      mpz_tdiv_q_2exp (r, t, modulus->bits);  /* r = (x + m*N) / R */
      if (mpz_cmp (r, modulus->orig_modulus) > 0)
	mpz_sub (r, r, modulus->orig_modulus);
    }
}

/* subquadratic REDC, at mpn level.
   {orig,n} is the original modulus.
   {aux,n} is the auxiliary modulus.
   Requires ABSIZ(x) = 2n and ABSIZ(orig_modulus)=ABSIZ(aux_modulus)=n.
 */
void
mpn_REDC (mp_ptr rp, mp_srcptr xp, mp_srcptr orig, mp_srcptr aux, mp_size_t n)
{
  mp_ptr tp, up;
  mp_size_t nn = n + n;
  mp_limb_t cy;
  TMP_DECL(marker);

  TMP_MARK(marker);
  up = TMP_ALLOC_LIMBS(nn + nn);
  mpn_mul_lo_n (up, xp, aux, n);
  tp = up + nn;
  mpn_mul_n (tp, up, orig, n);
  /* add {x, 2n} and {tp, 2n}. We know that {tp, n} + {xp, n} will give
     either 0, or a carry out. If xp[n-1] <> 0, then there is a carry. */
#ifdef HAVE_NATIVE_mpn_add_nc
  cy = mpn_add_nc (rp, tp + n, xp + n, n, (mp_limb_t) ((xp[n - 1]) ? 1 : 0));
#else
  cy = mpn_add_n (rp, tp + n, xp + n, n);
  cy += mpn_add_1 (rp, rp, n, (mp_limb_t) ((xp[n - 1]) ? 1 : 0));
#endif
  if (cy || mpn_cmp (rp, orig, n) > 0)
    cy -= mpn_sub_n (rp, rp, orig, n);
  /* ASSERT ((cy == 0) && (mpn_cmp (rp, orig, n) < 0)); */
  TMP_FREE(marker);
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

/* r <- c/R^nn mod n, where n are nn limbs.
   n must be odd.
   c must have space for at least 2*nn+1 limbs.
   c and r can be the same variable.
   The data in c is clobbered.
*/
void 
mpz_mod_n (mpz_ptr r, mpz_ptr c, mpmod_t modulus)
{
  mp_ptr rp;
  mp_ptr cp;
  mp_srcptr np;
  mp_limb_t cy;
  mp_ptr cys; /* vector of carries */
  mp_limb_t q;
  mp_size_t j, nn = modulus->bits / __GMP_BITS_PER_MP_LIMB;
  TMP_DECL(marker);

  if (ALLOC(r) < nn)
    _mpz_realloc (r, nn);
  cp = PTR(c);
  rp = PTR(r);
  np = PTR(modulus->orig_modulus); 
  ASSERT(ALLOC(c) >= 2 * nn + 1);
  ASSERT(ABSIZ(c) <= 2 * nn + 1);
  for (j = ABSIZ(c); j <= 2 * nn; j++) 
    cp[j] = 0;
  TMP_MARK(marker);
  cys = TMP_ALLOC_LIMBS(nn);
  for (j = 0; j < nn; j++)
    {
      q = cp[0] * modulus->Nprim;
      cys[j] = mpn_addmul_1 (cp, np, nn, q);
      cp++;
    }
  cy = cp[nn];
  /* add vector of carries and shift */
  cy += mpn_add_n (rp, cp, cys, nn);
  TMP_FREE(marker);
  if (cy)
    {
      cy -= mpn_sub_n (rp, rp, np, nn);
      while (cy > 1) /* subtract cy * {np, nn} */
        cy -= mpn_submul_1 (rp, np, nn, cy);
      while (cy) /* subtract {np, nn} */
        cy -= mpn_sub_n (rp, rp, np, nn);
    }
  MPN_NORMALIZE (rp, nn);
  SIZ(r) = SIZ(c) < 0 ? -nn : nn;
}

/* don't use base2 if repr == -1 */
void 
mpmod_init (mpmod_t modulus, mpz_t N, int repr, int verbose)
{
  int base2;
  
  if ((repr != -1) && (base2 = isbase2 (N, 2.0)))
    {
      if (verbose > 1)
	printf ("Using special division for factor of 2^%d%c1\n",
		abs (base2), (base2 < 0) ? '-' : '+');
      mpmod_init_BASE2 (modulus, base2, N);
    }
  else if (mpz_size (N) < 3 * DIV_DC_THRESHOLD / 2)
    {
      if (verbose > 1)
	printf ("Using MODMULN\n");
      mpmod_init_MODMULN (modulus, N);
    }
  else if (mpz_sizeinbase (N, 2) < MOD_PLAIN_TO_REDC_THRESHOLD)
    {
      if (verbose > 1)
	printf ("Using mpz_mod\n");
      mpmod_init_MPZ (modulus, N);
    }
  else
    {
      if (verbose > 1)
	printf("Using REDC\n");
      mpmod_init_REDC (modulus, N);
    }
  
  return;
}

void 
mpmod_init_MPZ (mpmod_t modulus, mpz_t N)
{
  int Nbits;
  
  mpz_init_set (modulus->orig_modulus, N);
  modulus->repr = MOD_PLAIN;
  
  Nbits = mpz_size (N) * __GMP_BITS_PER_MP_LIMB; /* Number of bits, rounded
                                                    up to full limb */
  mpz_init2 (modulus->temp1, 2 * Nbits + __GMP_BITS_PER_MP_LIMB);
  mpz_init2 (modulus->temp2, Nbits);
  
  return;
}

void 
mpmod_init_BASE2 (mpmod_t modulus, int base2, mpz_t N)
{
  int Nbits;
  
  mpz_init_set (modulus->orig_modulus, N);
  modulus->repr = MOD_BASE2;
  modulus->bits = base2;
  
  Nbits = mpz_size (N) * __GMP_BITS_PER_MP_LIMB; /* Number of bits, rounded
                                                    up to full limb */
  mpz_init2 (modulus->temp1, 2 * Nbits + __GMP_BITS_PER_MP_LIMB);
  mpz_init2 (modulus->temp2, Nbits);
  
  return;
}

void 
mpmod_init_MODMULN (mpmod_t modulus, mpz_t N)
{
  int Nbits;

  mpz_init_set (modulus->orig_modulus, N);
  
  modulus->repr = MOD_MODMULN;
  Nbits = mpz_size (N) * __GMP_BITS_PER_MP_LIMB; /* Number of bits, rounded
                                                    up to full limb */
  modulus->bits = Nbits;

  mpz_init2 (modulus->temp1, 2 * Nbits + __GMP_BITS_PER_MP_LIMB);
  mpz_init2 (modulus->temp2, Nbits);

  mpz_set_ui (modulus->temp1, 1);
  mpz_mul_2exp (modulus->temp1, modulus->temp1, __GMP_BITS_PER_MP_LIMB);
  mpz_tdiv_r_2exp (modulus->temp2, modulus->orig_modulus, 
                   __GMP_BITS_PER_MP_LIMB);
  mpz_invert (modulus->temp2, modulus->temp2, modulus->temp1);
    /* Now temp2 = 1/n (mod 2^bits_per_limb) */
  mpz_sub (modulus->temp2, modulus->temp1, modulus->temp2);
  modulus->Nprim = mpz_getlimbn (modulus->temp2, 0);
    /* Now Nprim = -1/n (mod 2^bits_per_limb) */

  mpz_init (modulus->R2);
  mpz_set_ui (modulus->temp1, 1);
  mpz_mul_2exp (modulus->temp1, modulus->temp1, 2 * Nbits);
  mpz_mod (modulus->R2, modulus->temp1, modulus->orig_modulus);
  /* Now R2 = (2^bits)^2 (mod N) */
  
  mpz_init (modulus->R3);
  mpz_mul_2exp (modulus->temp1, modulus->R2, Nbits);
  mpz_mod (modulus->R3, modulus->temp1, modulus->orig_modulus);
  /* Now R3 = (2^bits)^3 (mod N) */

  return;
}

void 
mpmod_init_REDC (mpmod_t modulus, mpz_t N)
{
  mp_size_t n;
  int Nbits;
  
  mpz_init_set (modulus->orig_modulus, N);
  
  n = mpz_size (N);
  modulus->repr = MOD_REDC;
  Nbits = n * __GMP_BITS_PER_MP_LIMB; /* Number of bits, rounded
                                                    up to full limb */
  modulus->bits = Nbits;
  
  mpz_init2 (modulus->temp1, 2 * Nbits + __GMP_BITS_PER_MP_LIMB);
  mpz_init2 (modulus->temp2, Nbits);
  mpz_init (modulus->aux_modulus);

  mpz_set_ui (modulus->temp1, 1);
  mpz_mul_2exp (modulus->temp1, modulus->temp1, Nbits);
  if (!mpz_invert (modulus->aux_modulus, N, modulus->temp1))
    {
      fprintf (stderr, "mpmod_init: could not invert N\n");
      exit (EXIT_FAILURE);
    }
  mpz_sub (modulus->aux_modulus, modulus->temp1, modulus->aux_modulus);
  /* ensure aux_modulus has n allocated limbs, for mpn_REDC */
  if (ABSIZ(modulus->aux_modulus) < n)
    {
      _mpz_realloc (modulus->aux_modulus, n);
      MPN_ZERO (PTR(modulus->aux_modulus) + ABSIZ(modulus->aux_modulus),
		n - ABSIZ(modulus->aux_modulus));
    }

  mpz_init (modulus->R2);
  mpz_set_ui (modulus->temp1, 1);
  mpz_mul_2exp (modulus->temp1, modulus->temp1, 2 * Nbits);
  mpz_mod (modulus->R2, modulus->temp1, modulus->orig_modulus);
  /* Now R2 = (2^bits)^2 (mod N) */
  
  mpz_init (modulus->R3);
  mpz_mul_2exp (modulus->temp1, modulus->R2, Nbits);
  mpz_mod (modulus->R3, modulus->temp1, modulus->orig_modulus);
  /* Now R3 = (2^bits)^3 (mod N) */
  
  return;
}

void 
mpmod_clear (mpmod_t modulus)
{
  mpz_clear (modulus->orig_modulus);
  mpz_clear (modulus->temp1);
  mpz_clear (modulus->temp2);
  if (modulus->repr == MOD_MODMULN || modulus->repr == MOD_REDC)
    {
      mpz_clear (modulus->R2);
      mpz_clear (modulus->R3);
    }
  if (modulus->repr == MOD_REDC)
    {
      mpz_clear (modulus->aux_modulus);
    }
  
  return;
}

void 
mpres_init (mpres_t R, mpmod_t modulus)
{
  mpz_init2 (R, mpz_sizeinbase (modulus->orig_modulus, 2));
}

void 
mpres_clear (mpres_t R, mpmod_t modulus)
{
  mpz_clear (R);
}

void 
mpres_set (mpres_t R, mpres_t S, mpmod_t modulus)
{
  mpz_set (R, S);
}

void 
mpres_swap (mpres_t R, mpres_t S, mpmod_t modulus)
{
  mpz_swap (R, S);
}

/* R <- BASE^EXP mod modulus */ 
void 
mpres_pow (mpres_t R, mpres_t BASE, mpres_t EXP, mpmod_t modulus)
{
  if (modulus->repr == MOD_PLAIN)
    {
      mpz_powm (R, BASE, EXP, modulus->orig_modulus);
    }
  else if (modulus->repr == MOD_BASE2 || modulus->repr == MOD_MODMULN ||
           modulus->repr == MOD_REDC)
    {
      unsigned int expidx;
      mp_limb_t bitmask, expbits;

      expidx = mpz_size (EXP) -1;           /* point at most significant limb */
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
                  return;
                }
              expidx--;
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

              if (modulus->repr == MOD_BASE2)
                base2mod (modulus->temp2 , modulus->temp1, modulus->temp1, modulus);
              else if (modulus->repr == MOD_MODMULN)
                {
                  mpz_mod_n (modulus->temp2, modulus->temp1, modulus);
                }
              else
                REDC (modulus->temp2, modulus->temp1, modulus->temp2, modulus);

              if (expbits & bitmask)
                { 
                  mpz_mul (modulus->temp1, modulus->temp2, BASE);
                  if (modulus->repr == MOD_BASE2)
                    base2mod (modulus->temp2, modulus->temp1, modulus->temp1, modulus);
                  else if (modulus->repr == MOD_MODMULN)
                    {
                      mpz_mod_n (modulus->temp2, modulus->temp1, modulus);
                    }
                  else
                    REDC (modulus->temp2, modulus->temp1, modulus->temp2, modulus);
                }
            }
          if (expidx == 0)		/* if we just processed the least */
            break;			/* significant limb, we are done */
          expidx--;
          expbits = mpz_getlimbn (EXP, expidx);
          bitmask = (mp_limb_t) 1 << (GMP_NUMB_BITS - 1);
        }
      mpz_set (R, modulus->temp2); /* TODO: isn't it possible to use R instead
				      of modulus->temp2 above to avoid this
				      copy? */
    } /* if (modulus->repr == MOD_BASE2 || ... ) */
}

/* R <- BASE^EXP mod modulus */ 
void 
mpres_ui_pow (mpres_t R, unsigned int BASE, mpres_t EXP, mpmod_t modulus)
{
  if (modulus->repr == MOD_PLAIN)
    {
      mpz_set_ui (modulus->temp1, BASE);
      mpz_powm (R, modulus->temp1, EXP, modulus->orig_modulus);
    }
  else if (modulus->repr == MOD_BASE2 || modulus->repr == MOD_MODMULN ||
           modulus->repr == MOD_REDC)
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
                  return;
                }
              expidx--;
              expbits = mpz_getlimbn (EXP, expidx);
              bitmask = (mp_limb_t) 1 << (GMP_NUMB_BITS - 1);
            }
        }

    /* here the most significant limb with any set bits is in expbits, */
    /* bitmask is set to mask in the msb of expbits */
    
      mpz_set_ui (modulus->temp2, BASE); /* temp2 = BASE */
      if (modulus->repr == MOD_MODMULN || modulus->repr == MOD_REDC)
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

              if (modulus->repr == MOD_BASE2)
                base2mod (modulus->temp2 , modulus->temp1, modulus->temp1, modulus);
              else if (modulus->repr == MOD_MODMULN)
                {
                  mpz_mod_n (modulus->temp2, modulus->temp1, modulus);
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
    } /* if (modulus->repr == MOD_BASE2 || ... ) */
}

void 
mpres_mul (mpres_t R, mpres_t S1, mpres_t S2, mpmod_t modulus)
{
  mpz_mul (modulus->temp1, S1, S2);

  switch (modulus->repr)
    {
    case MOD_BASE2:
      base2mod (R, modulus->temp1, modulus->temp1, modulus);
      return;
    case MOD_MODMULN:
      mpz_mod_n (R, modulus->temp1, modulus);
      return;
    case MOD_REDC:
      REDC (R, modulus->temp1, modulus->temp2, modulus);
      return;
    default:
      mpz_mod (R, modulus->temp1, modulus->orig_modulus);
      return;
    }
}

void 
mpres_mul_ui (mpres_t R, mpres_t S, unsigned int n, mpmod_t modulus)
{
  mpz_mul_ui (modulus->temp1, S, n);
  /* This is the same for all methods: just reduce with original modulus */
  mpz_mod (R, modulus->temp1, modulus->orig_modulus);
}

/* R <- S / 2^n mod modulus. Does not need to be fast. */
void 
mpres_div_2exp (mpres_t R, mpres_t S, unsigned int n, mpmod_t modulus)
{
  if (n == 0)
    {
      mpres_set (R, S, modulus);
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

  return;
}

void
mpres_add_ui (mpres_t R, mpres_t S, unsigned int n, mpmod_t modulus)
{
  if (modulus->repr == MOD_PLAIN || modulus->repr == MOD_BASE2)
    {
      mpz_add_ui (R, S, n);
      if (mpz_cmp (R, modulus->orig_modulus) > 0)
        mpz_sub (R, R, modulus->orig_modulus);
    }
  else if (modulus->repr == MOD_MODMULN || modulus->repr == MOD_REDC)
    {
      mpz_set_ui (modulus->temp1, n);
      mpz_mul_2exp (modulus->temp1, modulus->temp1, modulus->bits);
      mpz_add (modulus->temp1, modulus->temp1, S);
      mpz_mod (R, modulus->temp1, modulus->orig_modulus);
    }
}

/* R <- S1 + S2 mod modulus.
   TODO: do we really need 0 <= R < modulus? */
void 
mpres_add (mpres_t R, mpres_t S1, mpres_t S2, mpmod_t modulus)
{
  mpz_add (R, S1, S2);
#ifdef FULL_REDUCTION
  if (mpz_cmp (R, modulus->orig_modulus) > 0)
#else
  if (ABSIZ (R) > ABSIZ (modulus->orig_modulus))
#endif
    mpz_sub (R, R, modulus->orig_modulus);
}

void
mpres_sub_ui (mpres_t R, mpres_t S, unsigned int n, mpmod_t modulus)
{
  if (modulus->repr == MOD_PLAIN || modulus->repr == MOD_BASE2)
    {
      mpz_sub_ui (R, S, n);
      if (mpz_sgn (R) < 0)
        mpz_add (R, R, modulus->orig_modulus);
    }
  else if (modulus->repr == MOD_MODMULN || modulus->repr == MOD_REDC)
    {
      mpz_set_ui (modulus->temp1, n);
      mpz_mul_2exp (modulus->temp1, modulus->temp1, modulus->bits);
      mpz_sub (modulus->temp1, S, modulus->temp1);
      mpz_mod (R, modulus->temp1, modulus->orig_modulus);
    }
  else
    {
      fprintf (stderr, "mpres_sub_ui: Unexpected  representation %d\n", 
               modulus->repr);
      exit (EXIT_FAILURE);
    }
}

/* R <- S1 - S2 mod modulus.
   TODO: do we really need 0 <= R < modulus? */
void 
mpres_sub (mpres_t R, mpres_t S1, mpres_t S2, mpmod_t modulus)
{
  mpz_sub (R, S1, S2);
#ifdef FULL_REDUCTION
  if (mpz_sgn (R) < 0)
#else
  if (ABSIZ (R) > ABSIZ (modulus->orig_modulus))
#endif
    mpz_add (R, R, modulus->orig_modulus);
}

void 
mpres_set_z (mpres_t R, mpz_t S, mpmod_t modulus)
{
  if (modulus->repr == MOD_PLAIN || modulus->repr == MOD_BASE2)
    {
      mpz_mod (R, S, modulus->orig_modulus);
    }
  else if (modulus->repr == MOD_MODMULN)
    {
      mpz_mod (modulus->temp2, S, modulus->orig_modulus);
      mpz_mul (modulus->temp1, modulus->temp2, modulus->R2);
      mpz_mod_n (R, modulus->temp1, modulus);
    }
  else if (modulus->repr == MOD_REDC)
    {
      mpz_mod (modulus->temp2, S, modulus->orig_modulus);
      mpz_mul (modulus->temp1, modulus->temp2, modulus->R2);
      REDC (R, modulus->temp1, modulus->temp2, modulus);
    }
  else
    {
      fprintf (stderr, "mpres_set_z: Unexpected  representation %d\n", 
               modulus->repr);
      exit (EXIT_FAILURE);
    }
}

/* S must not be modulus->temp1 for REDC */
void 
mpres_get_z (mpz_t R, mpres_t S, mpmod_t modulus)
{
  if (modulus->repr == MOD_PLAIN || modulus->repr == MOD_BASE2)
    {
      mpz_mod (R, S, modulus->orig_modulus);
    }
  else if (modulus->repr == MOD_MODMULN)
    {
      mpz_set (modulus->temp1, S);
      mpz_mod_n (R, modulus->temp1, modulus);
    }
  else if (modulus->repr == MOD_REDC)
    {
      REDC (R, S, modulus->temp1, modulus);
    }
  else
    {
      fprintf (stderr, "mpres_get_z: Unexpected representation %d\n", 
               modulus->repr);
      exit (EXIT_FAILURE);
    }
}

void 
mpres_set_ui (mpres_t R, unsigned int n, mpmod_t modulus)
{
  if (modulus->repr == MOD_PLAIN || modulus->repr == MOD_BASE2)
    {
      mpz_set_ui (R, n);
      mpz_mod (R, R, modulus->orig_modulus);
    }
  else if (modulus->repr == MOD_MODMULN || modulus->repr == MOD_REDC)
    {
      mpz_set_ui (modulus->temp1, n);
      mpz_mul_2exp (modulus->temp1, modulus->temp1, modulus->bits);
      mpz_mod (R, modulus->temp1, modulus->orig_modulus);
    }
  else
    {
      fprintf (stderr, "mpres_set_ui: Unexpected representation %d\n", 
               modulus->repr);
      exit (EXIT_FAILURE);
    }
}

/* R <- -S mod modulus. Does not need to be efficient. */
void
mpres_neg (mpres_t R, mpres_t S, mpmod_t modulus)
{ /* TODO: R < modulus  assumes  0 < S */
  mpz_sub (R, modulus->orig_modulus, S);
}

int 
mpres_invert (mpres_t R, mpres_t S, mpmod_t modulus)
{
  if (modulus->repr == MOD_PLAIN || modulus->repr == MOD_BASE2)
    {
      return mpz_invert (R, S, modulus->orig_modulus);
    }
  else if (modulus->repr == MOD_MODMULN)
    {
      if (mpz_invert (modulus->temp2, S, modulus->orig_modulus))
        {
          mpz_mul (modulus->temp1, modulus->temp2, modulus->R3);
          mpz_mod_n (R, modulus->temp1, modulus);
          return 1;
        }
      else
        return 0;
    }
  else if (modulus->repr == MOD_REDC)
    {
      if (mpz_invert (modulus->temp2, S, modulus->orig_modulus))
        {
          mpz_mul (modulus->temp1, modulus->temp2, modulus->R3);
          REDC (R, modulus->temp1, modulus->temp2, modulus);
          return 1;
        }
      else
        return 0;
    }
  else
    {
      fprintf (stderr, "mpres_invert: Unexpected representation %d\n", 
               modulus->repr);
      exit (EXIT_FAILURE);
    }
}

void 
mpres_gcd (mpz_t R, mpres_t S, mpmod_t modulus)
{
  /* In MODMULN and REDC form, M(x) = x*R with gcd(R, modulus) = 1 .
     Therefore gcd(M(x), modulus) = gcd(x, modulus) and we need not bother
     to convert out of Montgomery form. */
  mpz_gcd (R, S, modulus->orig_modulus);
}

void 
mpres_out_str (FILE *fd, unsigned int base, mpres_t S, mpmod_t modulus)
{
  mpres_get_z (modulus->temp1, S, modulus);
  mpz_out_str (fd, base, modulus->temp1);
}
