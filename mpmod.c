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

/* #define DEBUG */
#ifndef MOD_PLAIN_TO_REDC_THRESOLD
#define MOD_PLAIN_TO_REDC_THRESOLD 50000
#endif

int isbase2 (mpz_t, double);
void REDC (mpres_t, mpres_t, mpz_t, mpmod_t);
void base2mod (mpres_t, mpres_t, mpres_t, mpmod_t);

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

/* Do base-2 reduction. R must not equal S or t. S may equal t. */
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

void 
mpmod_init(mpmod_t modulus, mpz_t N)
{
  int base2;
  mpz_t R;
  
  mpz_init_set(modulus->orig_modulus, N);
  mpz_init(modulus->temp1);
  mpz_init(modulus->temp2);
  
  if ((base2 = isbase2 (N, 2.0)))
    {
      printf("Using base-2: 2^%d %c 1\n", abs(base2), (base2<0)?'-':'+');
      modulus->repr = MOD_BASE2;
      modulus->bits = base2;
    }
  else if (mpz_sizeinbase (N, 2) < MOD_PLAIN_TO_REDC_THRESOLD)
    {
      printf("Using plain mpz_mod\n");
      modulus->repr = MOD_PLAIN;
    }
  else
    {
      /* Init for REDC */
      printf("Using REDC\n");
      
      mpz_init(modulus->aux_modulus);
      mpz_init(modulus->RmodN);
      modulus->repr = MOD_REDC;

      modulus->bits = mpz_sizeinbase(N, 2) + mp_bits_per_limb;
      modulus->bits -= modulus->bits % mp_bits_per_limb;
      
      mpz_init_set_ui(R, 1);
      mpz_mul_2exp(R, R, modulus->bits);
      if (!mpz_invert(modulus->aux_modulus, N, R))
        {
          fprintf(stderr, "mpmod_init: could not invert N\n");
          exit(EXIT_FAILURE);
        }
      mpz_sub(modulus->aux_modulus, R, modulus->aux_modulus);
      mpz_mod (modulus->RmodN, R, N);
      mpz_clear(R);
    }
  
  return;
}

void 
mpmod_clear (mpmod_t modulus)
{
  mpz_clear(modulus->orig_modulus);
  mpz_clear(modulus->temp1);
  mpz_clear(modulus->temp2);
  if (modulus->repr == MOD_REDC)
    {
      mpz_clear(modulus->aux_modulus);
      mpz_clear(modulus->RmodN);
    }
  
  return;
}

void 
mpres_init (mpres_t R, mpmod_t modulus)
{
  mpz_init2(R, mpz_size(modulus->orig_modulus));
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
  else if (modulus->repr == MOD_BASE2)
    {
      unsigned int expidx;
      mp_limb_t bitmask, expbits;

      expidx = mpz_size (EXP) -1;           /* point at most significant limb */
      expbits = mpz_getlimbn (EXP, expidx); /* get most significant limb */
      bitmask = 1 << (__GMP_BITS_PER_MP_LIMB - 1);

      while ((bitmask & expbits) == 0)
        {
          bitmask >>= 1;
          if (bitmask == 0)                 /* no set bits in this limb */
            {
              if (expidx == 0)              /* no more limbs -> exp was 0 */
                {
                  mpz_set_ui(R, 1);         /* set result to 1 and return */
                  return;
                }
              expidx--;
              expbits = mpz_getlimbn (EXP, expidx);
              bitmask = 1 << (__GMP_BITS_PER_MP_LIMB - 1);
            }
        }

    /* here the most significant limb with any set bits is in expbits, */
    /* bitmask is set at the msb of expbits */

      mpz_set (modulus->temp2, BASE); /* In case R = BASE */
      mpz_set (R, BASE);
      bitmask >>= 1;

      while (1) 
        {
          for ( ; bitmask != 0; bitmask >>= 1) 
            {
              mpz_mul(modulus->temp1, R, R);    /* r = r^2 */
              base2mod (R, modulus->temp1, modulus->temp1, modulus);
              if (expbits & bitmask)
                { 
                  mpz_mul (modulus->temp1, R, modulus->temp2);
                  base2mod (R, modulus->temp1, modulus->temp1, modulus);
                }
            }
          if (expidx == 0)			/* if we just processed the least */
            break;				/* significant limb, we are done */
          expidx--;
          expbits = mpz_getlimbn (EXP, expidx);
          bitmask = 1 << (__GMP_BITS_PER_MP_LIMB - 1);
        }
    } /* if (modulus->repr == MOD_BASE2) */
  else
    { /* REDC */
      unsigned int expidx;
      mp_limb_t bitmask, expbits;

      expidx = mpz_size (EXP) -1;           /* point at most significant limb */
      expbits = mpz_getlimbn (EXP, expidx); /* get most significant limb */
      bitmask = 1 << (__GMP_BITS_PER_MP_LIMB - 1);

      while ((bitmask & expbits) == 0)
        {
          bitmask >>= 1;
          if (bitmask == 0)                 /* no set bits in this limb */
            {
              if (expidx == 0)              /* no more limbs -> EXP was 0 */
                {
                  mpz_set_ui(modulus->temp1, 1); /* set result to 1 and return */
                  mpz_mul_2exp(modulus->temp1, modulus->temp1, modulus->bits);
                  mpz_mod(R, modulus->temp1, modulus->orig_modulus);
                  return;
                }
              expidx--;
              expbits = mpz_getlimbn (EXP, expidx);
              bitmask = 1 << (__GMP_BITS_PER_MP_LIMB - 1);
            }
        }

    /* here the most significant limb with any set bits is in expbits, */
    /* bitmask is set at the msb of expbits */

      mpz_set (modulus->temp1, BASE);
      bitmask >>= 1;

      while (1) 
        {
          for ( ; bitmask != 0; bitmask >>= 1) 
            {
              mpz_mul(modulus->temp2, modulus->temp1, modulus->temp1);
              REDC(modulus->temp1, modulus->temp2, modulus->temp1, modulus);
              if (expbits & bitmask)
                { 
                  mpz_mul (modulus->temp2, modulus->temp1, BASE);
                  REDC(modulus->temp1, modulus->temp2, modulus->temp1, modulus);
                }
            }
          if (expidx == 0)			/* if we just processed the least */
            break;				/* significant limb, we are done */
          expidx--;
          expbits = mpz_getlimbn (EXP, expidx);
          bitmask = 1 << (__GMP_BITS_PER_MP_LIMB - 1);
       }
     mpz_set(R, modulus->temp1);
    }
}

void 
mpres_mul(mpres_t R, mpres_t S1, mpres_t S2, mpmod_t modulus)
{
  mpz_mul(modulus->temp1, S1, S2);
  
  if (modulus->repr == MOD_BASE2)
    base2mod (R, modulus->temp1, modulus->temp1, modulus);
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
  if (modulus->repr == MOD_REDC)
    mpz_mod(R, modulus->temp1, modulus->orig_modulus);
  else
    mpz_mod(R, modulus->temp1, modulus->orig_modulus);
}

void 
mpres_div_2exp(mpres_t R, mpres_t S, unsigned int n, mpmod_t modulus)
{
  if (n == 0)
    {
      mpres_set(R, S, modulus);
      return;
    }
  
  if (modulus->repr == MOD_REDC)
    {
      mpres_get_z(modulus->temp2, S, modulus);
      
      for ( ; n > 0; n--)
        if (mpz_odd_p (modulus->temp2))
          {
            mpz_add(modulus->temp2, modulus->temp2, modulus->orig_modulus);
            mpz_tdiv_q_2exp (modulus->temp2, modulus->temp2, 1);
          }
        else
          mpz_tdiv_q_2exp (modulus->temp2, modulus->temp2, 1);
      
      mpres_set_z(R, modulus->temp2, modulus);
    }
  else
    {
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
    }

  return;
}

void
mpres_add_ui (mpres_t R, mpres_t S, unsigned int n, mpmod_t modulus)
{
  if (modulus->repr == MOD_REDC)
    {
      mpz_mul_ui (modulus->temp1, modulus->RmodN, n);
      mpz_add (modulus->temp1, modulus->temp1, S);
      mpz_mod (R, modulus->temp1, modulus->orig_modulus);
    }
  else
    {
      mpz_add_ui(R, S, n);
      if (mpz_cmp(R, modulus->orig_modulus) > 0)
        mpz_sub(R, R, modulus->orig_modulus);
    }
}

void 
mpres_add (mpres_t R, mpres_t S1, mpres_t S2, mpmod_t modulus)
{
  mpz_add(R, S1, S2);
  if (modulus->repr == MOD_REDC)
    mpz_mod(R, R, modulus->orig_modulus);
  else
    if (mpz_cmp(R, modulus->orig_modulus) > 0)
      mpz_sub(R, R, modulus->orig_modulus);
}

void
mpres_sub_ui (mpres_t R, mpres_t S, unsigned int n, mpmod_t modulus)
{
  if (modulus->repr == MOD_REDC)
    {
      mpz_mul_ui (modulus->temp1, modulus->RmodN, n);
      mpz_sub (modulus->temp1, S, modulus->temp1);
      mpz_mod (R, modulus->temp1, modulus->orig_modulus);
    }
  else
    {
      mpz_sub_ui(R, S, n);
      if (mpz_sgn (R) < 0)
        mpz_add(R, R, modulus->orig_modulus);
    }
}

void 
mpres_sub (mpres_t R, mpres_t S1, mpres_t S2, mpmod_t modulus)
{
  mpz_sub(R, S1, S2);
  if (mpz_sgn(R) < 0)
    mpz_add(R, R, modulus->orig_modulus);
}

void 
mpres_set_z (mpres_t R, mpz_t S, mpmod_t modulus)
{
  if (modulus->repr == MOD_REDC)
    {
      mpz_mul_2exp(modulus->temp1, S, modulus->bits);
      mpz_mod(R, modulus->temp1, modulus->orig_modulus);
    }
  else
    {
      mpz_mod(R, S, modulus->orig_modulus);
    }
}

/* S must not be modulus->temp1 */
void 
mpres_get_z (mpz_t R, mpres_t S, mpmod_t modulus)
{
  if (modulus->repr == MOD_REDC)
    {
      REDC(R, S, modulus->temp1, modulus);
    }
  else
    {
      mpz_mod(R, S, modulus->orig_modulus);
    }
}

void 
mpres_set_ui (mpres_t R, unsigned int n, mpmod_t modulus)
{
  if (modulus->repr == MOD_REDC)
    {
      mpz_mul_ui(modulus->temp1, modulus->RmodN, n);
      mpz_mod(R, modulus->temp1, modulus->orig_modulus);
    }
  else
    {
      mpz_set_ui(R, n);
      mpz_mod(R, R, modulus->orig_modulus);
    }
}

void
mpres_neg (mpres_t R, mpres_t S, mpmod_t modulus)
{
  mpz_sub(R, modulus->orig_modulus, S);
}

int 
mpres_invert (mpres_t R, mpres_t S, mpmod_t modulus)
{
  if (modulus->repr == MOD_REDC)
    {
      /* Todo: make simpler */
      REDC(modulus->temp1, S, modulus->temp1, modulus);
      if (mpz_invert(modulus->temp2, modulus->temp1, modulus->orig_modulus))
        {
          mpz_mul_2exp(modulus->temp1, modulus->temp2, modulus->bits);
          mpz_mod(R, modulus->temp1, modulus->orig_modulus);
          return 1;
        }
      else
        return 0;
    }
  else
    return mpz_invert(R, S, modulus->orig_modulus);
}

void 
mpres_gcd(mpz_t R, mpres_t S, mpmod_t modulus)
{
  if (modulus->repr == MOD_REDC)
    {
      REDC (modulus->temp1, S, modulus->temp1, modulus);
      mpz_gcd (R, modulus->temp1, modulus->orig_modulus);
    }
  else
    mpz_gcd (R, S, modulus->orig_modulus);
}

void 
mpres_out_str (FILE *fd, unsigned int base, mpres_t a, mpmod_t modulus)
{
  mpres_get_z(modulus->temp1, a, modulus);
  mpz_out_str(fd, base, modulus->temp1);
}
