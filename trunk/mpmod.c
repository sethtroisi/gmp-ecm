/* Modular multiplication.

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

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "ecm.h"
#include <gmp-impl.h>
#include <longlong.h>

#define DEBUG
#ifndef MOD_PLAIN_TO_MODMULN_THRESHOLD
#define MOD_PLAIN_TO_MODMULN_THRESHOLD 50
#endif
#ifndef MOD_MODMULN_TO_REDC_THRESHOLD
#define MOD_MODMULN_TO_REDC_THRESHOLD 50000
#endif

int isbase2 (mpz_t, double);
void base2mod (mpres_t, mpres_t, mpres_t, mpmod_t);
void REDC (mpres_t, mpres_t, mpz_t, mpmod_t);
void mod_mul2exp (mpz_t c, unsigned int k, mpmod_t n);
void mod_div2exp (mpz_t c, unsigned int k, mpmod_t n);
static void mpn_incr (mp_ptr p, mp_limb_t incr);
void mpz_mod_n (mpz_t, mpmod_t);

/* returns +/-k if n is a factor of N = 2^k +/- 1 with N < =n^threshold, 
   0 otherwise 
   TODO: Complete Fermat numbers (no prime factor divided out) are not
         recognized 
*/
int 
isbase2(mpz_t n, double threshold)
{
  unsigned int k, lo; 
  int res = 0; 
  mpz_t u, w;

  mpz_init (u);
  mpz_init (w);
  lo = mpz_sizeinbase (n, 2) - 1;
  mpz_set_ui (u, 1);
  mpz_mul_2exp (u, u, 2 * lo);
  mpz_mod (w, u, n); /* 2^(2lo) mod n = +/-2^(2lo-k) if m*n = 2^k+/-1 */
  k = mpz_sizeinbase (w, 2) - 1;
  /* try w = 2^k */
  mpz_set_ui (u, 1);
  mpz_mul_2exp (u, u, k);
  if (mpz_cmp(w, u) == 0) 
    res = k - 2 * lo;
  else 
    {
      /* try w = -2^k */
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
    res=0;
  return (res);
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
void 
REDC (mpres_t r, mpres_t x, mpz_t t, mpmod_t modulus)
{
  mpz_tdiv_r_2exp(t, x, modulus->bits);
  mpz_mul(t, t, modulus->aux_modulus);
  mpz_tdiv_r_2exp(t, t, modulus->bits);  /* t = (x % R) * 1/N (mod R) */
  mpz_mul(t, t, modulus->orig_modulus);
  mpz_add(t, t, x);
  mpz_tdiv_q_2exp(t, t, modulus->bits);  /* t = (x + m*N) / R */
  if (mpz_cmp(t, modulus->orig_modulus) > 0)
    mpz_sub(r, t, modulus->orig_modulus);
  else
    mpz_set(r, t);
}

/* multiplies c by R^k modulo n where R=2^mp_bits_per_limb 
   n is supposed odd. Does not need to be efficient. */
void 
mod_mul2exp (mpz_t c, unsigned int k, mpmod_t n)
{
  mpz_mul_2exp (n->temp1, c, k * __GMP_BITS_PER_MP_LIMB);
  mpz_mod (c, n->temp1, n->orig_modulus);
}

/* divides c by R^k modulo n where R=2^mp_bits_per_limb
   n is supposed odd. Does not need to be efficient. */
void 
mod_div2exp (mpz_t c, unsigned int k, mpmod_t n)
{
  mpz_set_ui (n->temp2, 1);
  mpz_mul_2exp (n->temp1, n->temp2, k * __GMP_BITS_PER_MP_LIMB);
  mpz_invert (n->temp2, n->temp1, n->orig_modulus); /* temp2 = 2^(-k) (mod n) */
  mpz_mul (n->temp1, n->temp2, c);
  mpz_mod (c, n->temp1, n->orig_modulus);
}

static void 
mpn_incr (mp_ptr p, mp_limb_t incr)
{
  mp_limb_t x;

  x = *p + incr;
  *p++ = x;
  if (x >= incr)
    return;
  while (++(*(p++)) == 0)
    ;
}

/* Computes c/R^nn mod n, where n are nn limbs
   and c has space for size(c)+1 limbs.  n must be odd.
*/
void 
mpz_mod_n (mpz_t c, mpmod_t modulus)
{
  mp_ptr cp = PTR (c), np = PTR (modulus->orig_modulus);
  mp_limb_t cy;
  mp_limb_t q;
  size_t j, nn = modulus->bits / __GMP_BITS_PER_MP_LIMB;

#ifdef DEBUG
  if (c->_mp_alloc < 2 * nn + 1)
    {
      fprintf (stderr, "mpz_mod_n: c has space for only %d limbs, needs %d\n",
               c->_mp_alloc, (int) (2 * nn + 1));
      exit (EXIT_FAILURE);
    }
#endif

  for (j = ABS (SIZ (c)); j <= 2 * nn; j++) 
    cp[j] = 0;
  for (j = 0; j < nn; j++)
    {
      q = cp[0] * modulus->Nprim;
      cy = mpn_addmul_1 (cp, np, nn, q);
      mpn_incr (cp + nn, cy);
      cp++;
    }
  cp -= nn;
  if (cp[2*nn]) 
    {
      cy = cp[2*nn] - mpn_sub_n (cp, cp + nn, np, nn);
      while (cy) cy -= mpn_sub_n (cp, cp, np, nn);
    }
  else 
    MPN_COPY (cp, cp + nn, nn);
  MPN_NORMALIZE (cp, nn);
  SIZ(c) = SIZ(c) < 0 ? -nn : nn;
}

void 
mpmod_init(mpmod_t modulus, mpz_t N)
{
  int base2, Nbits;
  mpz_t R;
  
  mpz_init_set (modulus->orig_modulus, N);
  
  Nbits = mpz_size (N) * __GMP_BITS_PER_MP_LIMB; /* Number of bits, rounded
                                                    up to full limb */
  mpz_init2 (modulus->temp1, 2 * Nbits + __GMP_BITS_PER_MP_LIMB);
  mpz_init2 (modulus->temp2, Nbits);
  
  if ((base2 = isbase2 (N, 2.0)))
    {
      printf("Using base-2: 2^%d %c 1\n", abs(base2), (base2<0)?'-':'+');
      modulus->repr = MOD_BASE2;
      modulus->bits = base2;
    }
  else if (mpz_sizeinbase (N, 2) < MOD_PLAIN_TO_MODMULN_THRESHOLD)
    {
      printf("Using plain mpz_mod\n");
      modulus->repr = MOD_PLAIN;
    }
  else if (mpz_sizeinbase (N, 2) < MOD_MODMULN_TO_REDC_THRESHOLD)
    {
      /* Init for MODMULN */
      printf("Using MODMULN\n");
      modulus->repr = MOD_MODMULN;

      modulus->bits = Nbits;

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
      mpz_mul_2exp (modulus->temp1, modulus->temp1, 2*Nbits);
      mpz_mod (modulus->R2, modulus->temp1, modulus->orig_modulus);
      /* Now R2 = (2^bits)^2 (mod N) */
      
      mpz_init (modulus->R3);
      mpz_mul_2exp (modulus->temp1, modulus->R2, Nbits);
      mpz_mod (modulus->R3, modulus->temp1, modulus->orig_modulus);
      /* Now R3 = (2^bits)^3 (mod N) */
    }
  else
    {
      /* Init for REDC */
      printf("Using REDC\n");
      
      mpz_init(modulus->aux_modulus);
      modulus->repr = MOD_REDC;

      modulus->bits = Nbits;
      
      mpz_init_set_ui(R, 1);
      mpz_mul_2exp(R, R, modulus->bits);
      if (!mpz_invert(modulus->aux_modulus, N, R))
        {
          fprintf(stderr, "mpmod_init: could not invert N\n");
          exit(EXIT_FAILURE);
        }
      mpz_sub(modulus->aux_modulus, R, modulus->aux_modulus);
      mpz_clear(R);

      mpz_set_ui (modulus->temp1, 1);
      mpz_mul_2exp (modulus->temp1, modulus->temp1, 2*Nbits);
      mpz_mod (modulus->R2, modulus->temp1, modulus->orig_modulus);
      /* Now R2 = (2^bits)^2 (mod N) */
      
      mpz_mul_2exp (modulus->temp1, modulus->R2, Nbits);
      mpz_mod (modulus->R3, modulus->temp1, modulus->orig_modulus);
      /* Now R3 = (2^bits)^3 (mod N) */
    }
  
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
  mpz_init2(R, mpz_sizeinbase (modulus->orig_modulus, 2));
}

void 
mpres_clear (mpres_t R, mpmod_t modulus)
{
  mpz_clear(R);
}

void 
mpres_set (mpres_t R, mpres_t S, mpmod_t n)
{
  mpz_set (R, S);
}

void 
mpres_swap (mpres_t R, mpres_t S, mpmod_t n)
{
  mpz_swap (R, S);
}

void 
mpres_pow(mpres_t R, mpres_t BASE, mpres_t EXP, mpmod_t modulus)
{
  if (modulus->repr == MOD_PLAIN)
    {
      mpz_powm(R, BASE, EXP, modulus->orig_modulus);
    }
  else if (modulus->repr == MOD_BASE2 || modulus->repr == MOD_MODMULN ||
           modulus->repr == MOD_REDC)
    {
      unsigned int expidx;
      mp_limb_t bitmask, expbits;

      expidx = mpz_size (EXP) -1;           /* point at most significant limb */
      expbits = mpz_getlimbn (EXP, expidx); /* get most significant limb */
      bitmask = (mp_limb_t) 1 << (__GMP_BITS_PER_MP_LIMB - 1);

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
              bitmask = (mp_limb_t) 1 << (__GMP_BITS_PER_MP_LIMB - 1);
            }
        }

    /* here the most significant limb with any set bits is in expbits, */
    /* bitmask is set to mask in the msb of expbits */

      mpz_set (modulus->temp2, BASE); /* In case R = BASE */
      bitmask >>= 1;

      while (1) 
        {
          for ( ; bitmask != 0; bitmask >>= 1) 
            {
              mpz_mul(modulus->temp1, modulus->temp2, modulus->temp2); /* r = r^2 */

              if (modulus->repr == MOD_BASE2)
                base2mod (modulus->temp2 , modulus->temp1, modulus->temp1, modulus);
              else if (modulus->repr == MOD_MODMULN)
                {
                  mpz_mod_n (modulus->temp1, modulus);
                  mpz_set (modulus->temp2, modulus->temp1); /* TODO: get rid of */
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
                      mpz_mod_n (modulus->temp1, modulus);
                      mpz_set (modulus->temp2, modulus->temp1); /* TODO: get rid of */
                    }
                  else
                    REDC (modulus->temp2, modulus->temp1, modulus->temp2, modulus);
                }
            }
          if (expidx == 0)			/* if we just processed the least */
            break;				/* significant limb, we are done */
          expidx--;
          expbits = mpz_getlimbn (EXP, expidx);
          bitmask = (mp_limb_t) 1 << (__GMP_BITS_PER_MP_LIMB - 1);
        }
      mpz_set (R, modulus->temp2);
    } /* if (modulus->repr == MOD_BASE2) */
}

void 
mpres_mul(mpres_t R, mpres_t S1, mpres_t S2, mpmod_t modulus)
{
  mpz_mul(modulus->temp1, S1, S2);
  
  if (modulus->repr == MOD_BASE2)
    base2mod (R, modulus->temp1, modulus->temp1, modulus);
  else if (modulus->repr == MOD_MODMULN)
    {
      mpz_mod_n (modulus->temp1, modulus);
      mpz_set (R, modulus->temp1);
    }
  else if (modulus->repr == MOD_REDC)
    REDC(R, modulus->temp1, modulus->temp2, modulus);
  else
    mpz_mod (R, modulus->temp1, modulus->orig_modulus);

  return;
}

void 
mpres_mul_ui (mpres_t R, mpres_t S, unsigned int n, mpmod_t modulus)
{
  mpz_mul_ui(modulus->temp1, S, n);
  /* This is the same for all methods: just reduce with original modulus */
  mpz_mod(R, modulus->temp1, modulus->orig_modulus);
}

/* TODO: make faster */
void 
mpres_div_2exp(mpres_t R, mpres_t S, unsigned int n, mpmod_t modulus)
{
  if (n == 0)
    {
      mpres_set(R, S, modulus);
      return;
    }

    if (mpz_odd_p (S))
      {
        mpz_add(R, S, modulus->orig_modulus);
        mpz_tdiv_q_2exp (R, R, 1);
      }
    else
      mpz_tdiv_q_2exp (R, S, 1);

    for ( ; n > 1; n--)
      if (mpz_odd_p (R))
        {
          mpz_add(R, R, modulus->orig_modulus);
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

void 
mpres_add (mpres_t R, mpres_t S1, mpres_t S2, mpmod_t modulus)
{
  mpz_add (R, S1, S2);
  if (mpz_cmp(R, modulus->orig_modulus) > 0)
    mpz_sub (R, R, modulus->orig_modulus);
}

void
mpres_sub_ui (mpres_t R, mpres_t S, unsigned int n, mpmod_t modulus)
{
  if (modulus->repr == MOD_PLAIN || modulus->repr == MOD_BASE2)
    {
      mpz_sub_ui(R, S, n);
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
      fprintf (stderr, "mpres_sub_ui: Unhandeled  representation %d\n", 
               modulus->repr);
      exit (EXIT_FAILURE);
    }
}

void 
mpres_sub (mpres_t R, mpres_t S1, mpres_t S2, mpmod_t modulus)
{
  mpz_sub (R, S1, S2);
  if (mpz_sgn (R) < 0)
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
      mpz_mul (modulus->temp1, S, modulus->R2);
      mpz_mod_n (modulus->temp1, modulus);
      mpz_set (R, modulus->temp1);
    }
  else if (modulus->repr == MOD_REDC)
    {
      mpz_mul (modulus->temp1, S, modulus->R2);
      REDC (R, modulus->temp1, modulus->temp2, modulus);
    }
  else
    {
      fprintf (stderr, "mpres_set_z: Unhandeled  representation %d\n", 
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
      mpz_mod_n (modulus->temp1, modulus);
      mpz_set (R, modulus->temp1);
    }
  else if (modulus->repr == MOD_REDC)
    {
      REDC(R, S, modulus->temp1, modulus);
    }
  else
    {
      fprintf (stderr, "mpres_get_z: Unhandeled representation %d\n", 
               modulus->repr);
      exit (EXIT_FAILURE);
    }
}

void 
mpres_set_ui (mpres_t R, unsigned int n, mpmod_t modulus)
{
  if (modulus->repr == MOD_PLAIN || modulus->repr == MOD_BASE2)
    {
      mpz_set_ui(R, n);
      mpz_mod(R, R, modulus->orig_modulus);
    }
  else if (modulus->repr == MOD_MODMULN || modulus->repr == MOD_REDC)
    {
      mpz_set_ui (modulus->temp1, n);
      mpz_mul_2exp (modulus->temp1, modulus->temp1, modulus->bits);
      mpz_mod(R, modulus->temp1, modulus->orig_modulus);
    }
  else
    {
      fprintf (stderr, "mpres_set_ui: Unhandeled representation %d\n", 
               modulus->repr);
      exit (EXIT_FAILURE);
    }
}

void
mpres_neg (mpres_t R, mpres_t S, mpmod_t modulus)
{ /* TODO: R < modulus  assumes  0 < S */
  mpz_sub(R, modulus->orig_modulus, S);
}

int 
mpres_invert (mpres_t R, mpres_t S, mpmod_t modulus)
{
  if (modulus->repr == MOD_PLAIN || modulus->repr == MOD_BASE2)
    {
      return mpz_invert(R, S, modulus->orig_modulus);
    }
  else if (modulus->repr == MOD_MODMULN)
    {
      if (mpz_invert (modulus->temp2, S, modulus->orig_modulus))
        {
          mpz_mul (modulus->temp1, modulus->temp2, modulus->R3);
          mpz_mod_n (modulus->temp1, modulus);
          mpz_set (R, modulus->temp1); /* TODO: get rid of */
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
      fprintf (stderr, "mpres_invert: Unhandeled representation %d\n", 
               modulus->repr);
      exit (EXIT_FAILURE);
    }
}

void 
mpres_gcd(mpz_t R, mpres_t S, mpmod_t modulus)
{
  /* In MODMULN and REDC form, M(x) = x*R with gcd(R, modulus) = 1 .
     Therefore gcd(M(x), modulus) = gcd(x, modulus) and we need not bother
     to convert out of Montgomery form. */
  mpz_gcd (R, S, modulus->orig_modulus);
}

void 
mpres_out_str (FILE *fd, unsigned int base, mpres_t S, mpmod_t modulus)
{
  mpres_get_z(modulus->temp1, S, modulus);
  mpz_out_str(fd, base, modulus->temp1);
}
