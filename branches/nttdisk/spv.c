/* spv.c - "small prime vector" functions for arithmetic on vectors of
   residues modulo a single small prime

  Copyright 2005, 2008 Dave Newman and Jason Papadopoulos.

  The SP Library is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or (at your
  option) any later version.

  The SP Library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
  License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with the SP Library; see the file COPYING.LIB.  If not, write to
  the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
  MA 02110-1301, USA.
*/

#include <string.h> /* for memset */
#include "ecm-impl.h"

/* Routines for vectors of integers modulo r common small prime
 * 
 * These are low-overhead routines that don't do memory allocation,
 * other than for temporary variables. Unless otherwise specified, any
 * of the input pointers can be equal. */

/* r = x */
void
spv_set (spv_t r, spv_t x, spv_size_t len)
{
#ifdef HAVE_MEMMOVE  
  /* memmove doesn't rely on the assertion below */
  memmove (r, x, len * sizeof (sp_t));
#else
  spv_size_t i;

  ASSERT (r >= x + len || x >= r);

  for (i = 0; i < len; i++)
    r[i] = x[i];
#endif
}

/* r[0 ... len - 1] = x[len - 1 ... 0]  */
void
spv_rev (spv_t r, spv_t x, spv_size_t len)
{
  spv_size_t i;
  
  ASSERT (r >= x + len || x >= r + len);

  for (i = 0; i < len; i++)
    r[i] = x[len - 1 - i];
}

/* r = [y, y, ... ] */
void
spv_set_sp (spv_t r, sp_t y, spv_size_t len)
{
  spv_size_t i;

  for (i = 0; i < len; i++)
    r[i] = y;
}

void
spv_set_zero (spv_t r, spv_size_t len)
{
  memset (r, 0, len * sizeof (sp_t));
}

int
spv_cmp (spv_t x, spv_t y, spv_size_t len)
{
  spv_size_t i;

  for (i = 0; i < len; i++)
    if (x[i] != y[i])
      return 1;

  return 0;
}

/* r = x + y */
void
spv_add (spv_t r, spv_t x, spv_t y, spv_size_t len, sp_t m)
{
  spv_size_t i;
  
  ASSERT (r >= x + len || x >= r);
  ASSERT (r >= y + len || y >= r);
  
  for (i = 0; i < len; i++)
    r[i] = sp_add (x[i], y[i], m);
}

/* r = [x[0] + y, x[1] + y, ... ] */
void
spv_add_sp (spv_t r, spv_t x, sp_t c, spv_size_t len, sp_t m)
{
  spv_size_t i;

  for (i = 0; i < len; i++)
    r[i] = sp_add (x[i], c, m);
}

/* r = x - y */
void
spv_sub (spv_t r, spv_t x, spv_t y, spv_size_t len, sp_t m)
{
  spv_size_t i;
  
  ASSERT (r >= x + len || x >= r);
  ASSERT (r >= y + len || y >= r);
   
  for (i = 0; i < len; i++)
    r[i] = sp_sub (x[i], y[i], m);
}

/* r = [x[0] - y, x[1] - y, ... ] */
void
spv_sub_sp (spv_t r, spv_t x, sp_t c, spv_size_t len, sp_t m)
{
  spv_size_t i;

  for (i = 0; i < len; i++)
    r[i] = sp_sub (x[i], c, m);
}

/* r = [-x[0], -x[1], ... ] */
void
spv_neg (spv_t r, spv_t x, spv_size_t len, sp_t m)
{
  spv_size_t i;

  for (i = 0; i < len; i++)
    r[i] = sp_sub (0, x[i], m);
}

/* Pointwise multiplication
 * r = [x[0] * y[0], x[1] * y[1], ... ] */
void
spv_pwmul (spv_t r, spv_t x, spv_t y, spv_size_t len, sp_t m, sp_t d)
{
  spv_size_t i = 0;
  
  ASSERT (r >= x + len || x >= r);
  ASSERT (r >= y + len || y >= r);

#if defined(HAVE_SSE2) && (SP_NUMB_BITS < 32)

  __m128i t0, t1, t2, t3, vm, vm2, vd;

  vm = pshufd(pcvt_i32(m), 0x00);
  vm2 = pshufd(pcvt_i32(m), 0x44);
  vd = pshufd(pcvt_i32(d), 0x00);

  for (; i < (len & (spv_size_t)(~3)); i += 4)
    {
      t0 = pload((__m128i *)(x + i));
      t1 = pload((__m128i *)(y + i));
      t2 = pshufd(t0, 0x31);
      t3 = pshufd(t1, 0x31);
      t0 = pmuludq(t0, t1);
      t2 = pmuludq(t2, t3);
      t1 = psrlq(t0, 2 * SP_NUMB_BITS - SP_TYPE_BITS);
      t3 = psrlq(t2, 2 * SP_NUMB_BITS - SP_TYPE_BITS);
      t1 = pmuludq(t1, vd);
      t3 = pmuludq(t3, vd);

      #if SP_NUMB_BITS < 31
      t1 = psrlq(t1, 33);
      t3 = psrlq(t3, 33);
      t1 = pmuludq(t1, vm);
      t3 = pmuludq(t3, vm);
      t0 = psubq(t0, t1);
      t2 = psubq(t2, t3);
      #else
      t1 = pshufd(t1, 0xf5);
      t3 = pshufd(t3, 0xf5);
      t1 = pmuludq(t1, vm);
      t3 = pmuludq(t3, vm);
      t0 = psubq(t0, t1);
      t2 = psubq(t2, t3);

      t0 = psubq(t0, vm2);
      t2 = psubq(t2, vm2);
      t1 = pshufd(t0, 0xf5);
      t3 = pshufd(t2, 0xf5);
      t1 = pand(t1, vm2);
      t3 = pand(t3, vm2);
      t0 = paddq(t0, t1);
      t2 = paddq(t2, t3);
      #endif

      t0 = pshufd(t0, 0x08);
      t1 = pshufd(t2, 0x08);
      t0 = punpcklo32(t0, t1);
      t0 = psubd(t0, vm);
      t1 = pcmpgtd(psetzero(), t0);
      t1 = pand(t1, vm);
      t0 = paddd(t0, t1);
      pstore((__m128i *)(r + i), t0);
    }

#endif

  for (; i < len; i++)
    r[i] = sp_mul (x[i], y[i], m, d);
}

/* Pointwise multiplication, second input is read in reverse
 * r = [x[0] * y[len - 1], x[1] * y[len - 2], ... x[len - 1] * y[0]] */
void
spv_pwmul_rev (spv_t r, spv_t x, spv_t y, spv_size_t len, sp_t m, sp_t d)
{
  spv_size_t i;
  
  ASSERT (r >= x + len || x >= r);
  ASSERT (r >= y + len || y >= r);

  for (i = 0; i < len; i++)
    r[i] = sp_mul (x[i], y[len - 1 - i], m, d);
}

/* dst = src * y */
void
spv_mul_sp (spv_t r, spv_t x, sp_t c, spv_size_t len, sp_t m, sp_t d)
{
  spv_size_t i = 0;
  
  ASSERT (r >= x + len || x >= r);
  
#if defined(HAVE_SSE2) && (SP_NUMB_BITS < 32)

  __m128i t0, t1, t2, t3, vm, vm2, vd, vc;

  vc = pshufd(pcvt_i32(c), 0x44);
  vm = pshufd(pcvt_i32(m), 0x00);
  vm2 = pshufd(pcvt_i32(m), 0x44);
  vd = pshufd(pcvt_i32(d), 0x00);

  for (; i < (len & (spv_size_t)(~3)); i += 4)
    {
      t0 = pload((__m128i *)(x + i));
      t2 = pshufd(t0, 0x31);
      t0 = pmuludq(t0, vc);
      t2 = pmuludq(t2, vc);
      t1 = psrlq(t0, 2 * SP_NUMB_BITS - SP_TYPE_BITS);
      t3 = psrlq(t2, 2 * SP_NUMB_BITS - SP_TYPE_BITS);
      t1 = pmuludq(t1, vd);
      t3 = pmuludq(t3, vd);

      #if SP_NUMB_BITS < 31
      t1 = psrlq(t1, 33);
      t3 = psrlq(t3, 33);
      t1 = pmuludq(t1, vm);
      t3 = pmuludq(t3, vm);
      t0 = psubq(t0, t1);
      t2 = psubq(t2, t3);
      #else
      t1 = pshufd(t1, 0xf5);
      t3 = pshufd(t3, 0xf5);
      t1 = pmuludq(t1, vm);
      t3 = pmuludq(t3, vm);
      t0 = psubq(t0, t1);
      t2 = psubq(t2, t3);

      t0 = psubq(t0, vm2);
      t2 = psubq(t2, vm2);
      t1 = pshufd(t0, 0xf5);
      t3 = pshufd(t2, 0xf5);
      t1 = pand(t1, vm2);
      t3 = pand(t3, vm2);
      t0 = paddq(t0, t1);
      t2 = paddq(t2, t3);
      #endif

      t0 = pshufd(t0, 0x08);
      t1 = pshufd(t2, 0x08);
      t0 = punpcklo32(t0, t1);
      t0 = psubd(t0, vm);
      t1 = pcmpgtd(psetzero(), t0);
      t1 = pand(t1, vm);
      t0 = paddd(t0, t1);

      pstore((__m128i *)(r + i), t0);
    }

#endif

  for (; i < len; i++)
    r[i] = sp_mul (x[i], c, m, d);
}

void
spv_random (spv_t x, spv_size_t len, sp_t m)
{
  spv_size_t i;
#if SP_TYPE_BITS == GMP_LIMB_BITS
  mpn_random ((mp_limb_t *)x, len);

#elif SP_TYPE_BITS < GMP_LIMB_BITS
  mpn_random ((mp_limb_t *)x, len / 2);
  if (len % 2)
    {
      mp_limb_t t;
      mpn_random (&t, 1);
      x[len - 1] = (sp_t)t;
    }

#else
  mpn_random ((mp_limb_t *)x, 2 * len);
#endif

  for (i = 0; i < len; i++)
#if SP_NUMB_BITS > SP_TYPE_BITS - 3
    while (x[i] >= m) 
      x[i] -= m;
#else
    x[i] %= m;
#endif
}
