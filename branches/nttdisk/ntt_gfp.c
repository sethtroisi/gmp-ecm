/* ntt_gfp.c - low-level radix-2 dif/dit ntt routines over GF(p)
   
Copyright 2005, 2006, 2007, 2008, 2009 Dave Newman, Jason Papadopoulos,
Brian Gladman, Alexander Kruppa, Paul Zimmermann.

The SP Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The SP Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the SP Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
MA 02110-1301, USA. */

#include "ecm-impl.h"

static void bfly(spv_t x0, spv_t x1, 
		spv_size_t len, sp_t p)
{
  spv_size_t i;

#if defined(HAVE_SSE2) 
  
  #if SP_NUMB_BITS < 32

  __m128i t0, t1, t2, t3, t4, t5, t6, t7;
  __m128i vp;

  vp = pshufd(pcvt_i32(p), 0x00);

  t4 = pload((__m128i *)(x0 + 0));
  t5 = pload((__m128i *)(x1 + 0));

  t0 = pload((__m128i *)(x0 + 4));
  t6 = paddd(t4, t5);
  t7 = psubd(t4, t5);
  t1 = pload((__m128i *)(x1 + 4));
  t6 = psubd(t6, vp);
  t4 = pcmpgtd(psetzero(), t6);
  t5 = pcmpgtd(psetzero(), t7);
  t4 = pand(t4, vp);
  t5 = pand(t5, vp);
  t6 = paddd(t6, t4);
  t7 = paddd(t7, t5);

  for (i = 8; i < len; i += 8)
    {
      t4 = pload((__m128i *)(x0 + i + 0));
      t2 = paddd(t0, t1);
      t3 = psubd(t0, t1);
      t5 = pload((__m128i *)(x1 + i + 0));
      t2 = psubd(t2, vp);
      t0 = pcmpgtd(psetzero(), t2);
      pstore((__m128i *)(x0 + i - 8), t6);
      t1 = pcmpgtd(psetzero(), t3);
      t0 = pand(t0, vp);
      pstore((__m128i *)(x1 + i - 8), t7);
      t1 = pand(t1, vp);
      t2 = paddd(t2, t0);
      t3 = paddd(t3, t1);

      t0 = pload((__m128i *)(x0 + i + 4));
      t6 = paddd(t4, t5);
      t7 = psubd(t4, t5);
      t1 = pload((__m128i *)(x1 + i + 4));
      t6 = psubd(t6, vp);
      t4 = pcmpgtd(psetzero(), t6);
      pstore((__m128i *)(x0 + i - 4), t2);
      t5 = pcmpgtd(psetzero(), t7);
      t4 = pand(t4, vp);
      pstore((__m128i *)(x1 + i - 4), t3);
      t5 = pand(t5, vp);
      t6 = paddd(t6, t4);
      t7 = paddd(t7, t5);
    }

  t2 = paddd(t0, t1);
  t3 = psubd(t0, t1);
  t2 = psubd(t2, vp);
  pstore((__m128i *)(x0 + i - 8), t6);
  t0 = pcmpgtd(psetzero(), t2);
  t1 = pcmpgtd(psetzero(), t3);
  pstore((__m128i *)(x1 + i - 8), t7);
  t0 = pand(t0, vp);
  t1 = pand(t1, vp);
  t2 = paddd(t2, t0);
  t3 = paddd(t3, t1);

  pstore((__m128i *)(x0 + i - 4), t2);
  pstore((__m128i *)(x1 + i - 4), t3);

  #else

  __m128i t0, t1, t2, t3, t4, t5, t6, t7;
  __m128i vp;

  vp = pshufd(pcvt_i64(p), 0x44);

  t4 = pload((__m128i *)(x0 + 0));
  t5 = pload((__m128i *)(x1 + 0));

  t0 = pload((__m128i *)(x0 + 2));
  t6 = paddq(t4, t5);
  t7 = psubq(t4, t5);
  t1 = pload((__m128i *)(x1 + 2));
  t6 = psubq(t6, vp);
  t4 = pcmpgtd(psetzero(), t6);
  t5 = pcmpgtd(psetzero(), t7);
  t4 = pshufd(t4, 0xf5);
  t5 = pshufd(t5, 0xf5);
  t4 = pand(t4, vp);
  t5 = pand(t5, vp);
  t6 = paddq(t6, t4);
  t7 = paddq(t7, t5);

  for (i = 4; i < len; i += 4)
    {
      t4 = pload((__m128i *)(x0 + i + 0));
      t2 = paddq(t0, t1);
      t3 = psubq(t0, t1);
      t5 = pload((__m128i *)(x1 + i + 0));
      t2 = psubq(t2, vp);
      t0 = pcmpgtd(psetzero(), t2);
      pstore((__m128i *)(x0 + i - 4), t6);
      t1 = pcmpgtd(psetzero(), t3);
      t0 = pshufd(t0, 0xf5);
      pstore((__m128i *)(x1 + i - 4), t7);
      t1 = pshufd(t1, 0xf5);
      t0 = pand(t0, vp);
      t1 = pand(t1, vp);
      t2 = paddq(t2, t0);
      t3 = paddq(t3, t1);

      t0 = pload((__m128i *)(x0 + i + 2));
      t6 = paddq(t4, t5);
      t7 = psubq(t4, t5);
      t1 = pload((__m128i *)(x1 + i + 2));
      t6 = psubq(t6, vp);
      t4 = pcmpgtd(psetzero(), t6);
      pstore((__m128i *)(x0 + i - 2), t2);
      t5 = pcmpgtd(psetzero(), t7);
      t4 = pshufd(t4, 0xf5);
      pstore((__m128i *)(x1 + i - 2), t3);
      t5 = pshufd(t5, 0xf5);
      t4 = pand(t4, vp);
      t5 = pand(t5, vp);
      t6 = paddq(t6, t4);
      t7 = paddq(t7, t5);
    }

  t2 = paddq(t0, t1);
  t3 = psubq(t0, t1);
  t2 = psubq(t2, vp);
  t0 = pcmpgtd(psetzero(), t2);
  pstore((__m128i *)(x0 + i - 4), t6);
  t1 = pcmpgtd(psetzero(), t3);
  t0 = pshufd(t0, 0xf5);
  pstore((__m128i *)(x1 + i - 4), t7);
  t1 = pshufd(t1, 0xf5);
  t0 = pand(t0, vp);
  t1 = pand(t1, vp);
  t2 = paddq(t2, t0);
  t3 = paddq(t3, t1);

  pstore((__m128i *)(x0 + i - 2), t2);
  pstore((__m128i *)(x1 + i - 2), t3);

  #endif

#else
  for (i = 0; i < len; i++)
    {
      sp_t t0 = x0[i];
      sp_t t1 = x1[i];
      sp_t t2, t3;
      t2 = sp_add (t0, t1, p);
      t3 = sp_sub (t0, t1, p);
      x0[i] = t2;
      x1[i] = t3;
    }
#endif
}

/*--------------------------- FORWARD NTT --------------------------------*/
static void bfly_dif(spv_t x0, spv_t x1, spv_t w,
			spv_size_t len, sp_t p, sp_t d)
{
  spv_size_t i;

#if defined(HAVE_SSE2) && (SP_NUMB_BITS < 32)

  __m128i t0, t1, t2, t3, vm, vm2, vd;

  vm = pshufd(pcvt_i32(p), 0x00);
  vm2 = pshufd(pcvt_i32(p), 0x44);
  vd = pshufd(pcvt_i32(d), 0x00);

  for (i = 0; i < len; i += 4)
    {
      t0 = pload((__m128i *)(x0 + i));
      t1 = pload((__m128i *)(x1 + i));
      t2 = paddd(t0, t1);
      t3 = psubd(t0, t1);
      t2 = psubd(t2, vm);

      t0 = pcmpgtd(psetzero(), t2);
      t1 = pcmpgtd(psetzero(), t3);
      t0 = pand(t0, vm);
      t1 = pand(t1, vm);
      t2 = paddd(t2, t0);
      t0 = pload((__m128i *)(w + i));
      t1 = paddd(t3, t1);

      pstore((__m128i *)(x0 + i), t2);

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

      pstore((__m128i *)(x1 + i), t0);
    }

#elif defined(HAVE_SSE2) && (GMP_LIMB_BITS == 32)

  __m128i t0, t1, t2, t3, t4, t5, vm, vd;

  vm = pshufd(pcvt_i64(p), 0x44);
  vd = pshufd(pcvt_i64(d), 0x44);

  for (i = 0; i < len; i += 2)
    {
      t0 = pload((__m128i *)(x0 + i));
      t1 = pload((__m128i *)(x1 + i));
      t2 = paddq(t0, t1);
      t3 = psubq(t0, t1);
      t2 = psubq(t2, vm);

      t0 = pcmpgtd(psetzero(), t2);
      t1 = pcmpgtd(psetzero(), t3);
      t0 = pshufd(t0, 0xf5);
      t1 = pshufd(t1, 0xf5);
      t0 = pand(t0, vm);
      t1 = pand(t1, vm);
      t2 = paddq(t2, t0);
      t4 = pload((__m128i *)(w + i));
      t5 = paddq(t3, t1);

      pstore((__m128i *)(x0 + i), t2);

      t3 = pshufd(t4, 0x31);
      t3 = pmuludq(t3, t5);
      t2 = pshufd(t5, 0x31);
      t2 = pmuludq(t2, t4);
      t1 = pshufd(t4, 0x31);
      t4 = pmuludq(t4, t5);
      t5 = pshufd(t5, 0x31);
      t5 = pmuludq(t5, t1);
      t3 = paddq(t2, t3);

      t0 = t4;
      t4 = psrlq(t4, 32);
      t1 = t3;
      t3 = psllq(t3, 32);
      t3 = paddq(t3, t0);
      t4 = paddq(t4, t1);
      t4 = psrlq(t4, 32);
      t4 = paddq(t4, t5);

      t0 = t3;
      t4 = psllq(t4, 2*(SP_TYPE_BITS - SP_NUMB_BITS));
      t0 = psrlq(t0, 2*SP_NUMB_BITS - SP_TYPE_BITS);
      t4 = paddq(t4, t0);

      t5 = pshufd(t4, 0x31);
      t5 = pmuludq(t5, vd);
      t2 = pshufd(vd, 0x31);
      t2 = pmuludq(t2, t4);
      t1 = pshufd(t4, 0x31);
      t4 = pmuludq(t4, vd);
      t0 = pshufd(vd, 0x31);
      t1 = pmuludq(t1, t0);
      t5 = paddq(t5, t2);

      t4 = psrlq(t4, 32);
      t5 = paddq(t5, t4);
      t5 = psrlq(t5, 32);
      t1 = paddq(t1, t5);
      t1 = psrlq(t1, 1);

      t5 = pshufd(t1, 0x31);
      t5 = pmuludq(t5, vm);
      t2 = pshufd(vm, 0x31);
      t2 = pmuludq(t2, t1);
      t1 = pmuludq(t1, vm);
      t5 = paddq(t5, t2);
      t5 = psllq(t5, 32);
      t1 = paddq(t1, t5);

      t3 = psubq(t3, t1);
      t3 = psubq(t3, vm);
      t1 = pcmpgtd(psetzero(), t3);
      t1 = pshufd(t1, 0xf5);
      t1 = pand(t1, vm);
      t3 = paddq(t3, t1);
      pstore((__m128i *)(x1 + i), t3);
    }

#else
  for (i = 0; i < len; i++)
    {
      sp_t w0 = w[i];
      sp_t t0 = x0[i];
      sp_t t1 = x1[i];
      sp_t t2, t3;
      t2 = sp_add (t0, t1, p);
      t3 = sp_sub (t0, t1, p);
      t3 = sp_mul (t3, w0, p, d);
      x0[i] = t2;
      x1[i] = t3;
    }
#endif
}

static void bfly_dif_sp(spv_t x0, spv_t x1, sp_t w,
			spv_size_t len, sp_t p, sp_t d)
{
  /* same as bfly_dif except all butterflies are scaled
     by a single value of w */

  spv_size_t i;

#if defined(HAVE_SSE2) && (SP_NUMB_BITS < 32)

  __m128i t0, t1, t2, t3, vm, vm2, vd, vw;

  vw = pshufd(pcvt_i32(w), 0x00);
  vm = pshufd(pcvt_i32(p), 0x00);
  vm2 = pshufd(pcvt_i32(p), 0x44);
  vd = pshufd(pcvt_i32(d), 0x00);

  for (i = 0; i < len; i += 4)
    {
      t0 = pload((__m128i *)(x0 + i));
      t1 = pload((__m128i *)(x1 + i));
      t2 = paddd(t0, t1);
      t3 = psubd(t0, t1);
      t2 = psubd(t2, vm);

      t0 = pcmpgtd(psetzero(), t2);
      t1 = pcmpgtd(psetzero(), t3);
      t0 = pand(t0, vm);
      t1 = pand(t1, vm);
      t2 = paddd(t2, t0);
      t0 = vw;
      t1 = paddq(t3, t1);

      pstore((__m128i *)(x0 + i), t2);

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

      pstore((__m128i *)(x1 + i), t0);
    }

#elif defined(HAVE_SSE2) && (GMP_LIMB_BITS == 32)

  __m128i t0, t1, t2, t3, t4, t5, vm, vd, vw;

  vm = pshufd(pcvt_i64(p), 0x44);
  vd = pshufd(pcvt_i64(d), 0x44);
  vw = pshufd(pcvt_i64(w), 0x44);

  for (i = 0; i < len; i += 2)
    {
      t0 = pload((__m128i *)(x0 + i));
      t1 = pload((__m128i *)(x1 + i));
      t2 = paddq(t0, t1);
      t3 = psubq(t0, t1);
      t2 = psubq(t2, vm);

      t0 = pcmpgtd(psetzero(), t2);
      t1 = pcmpgtd(psetzero(), t3);
      t0 = pshufd(t0, 0xf5);
      t1 = pshufd(t1, 0xf5);
      t0 = pand(t0, vm);
      t1 = pand(t1, vm);
      t2 = paddq(t2, t0);
      t4 = vw;
      t5 = paddq(t3, t1);

      pstore((__m128i *)(x0 + i), t2);

      t3 = pshufd(t4, 0x31);
      t3 = pmuludq(t3, t5);
      t2 = pshufd(t5, 0x31);
      t2 = pmuludq(t2, t4);
      t1 = pshufd(t4, 0x31);
      t4 = pmuludq(t4, t5);
      t5 = pshufd(t5, 0x31);
      t5 = pmuludq(t5, t1);
      t3 = paddq(t2, t3);

      t0 = t4;
      t4 = psrlq(t4, 32);
      t1 = t3;
      t3 = psllq(t3, 32);
      t3 = paddq(t3, t0);
      t4 = paddq(t4, t1);
      t4 = psrlq(t4, 32);
      t4 = paddq(t4, t5);

      t0 = t3;
      t4 = psllq(t4, 2*(SP_TYPE_BITS - SP_NUMB_BITS));
      t0 = psrlq(t0, 2*SP_NUMB_BITS - SP_TYPE_BITS);
      t4 = paddq(t4, t0);

      t5 = pshufd(t4, 0x31);
      t5 = pmuludq(t5, vd);
      t2 = pshufd(vd, 0x31);
      t2 = pmuludq(t2, t4);
      t1 = pshufd(t4, 0x31);
      t4 = pmuludq(t4, vd);
      t0 = pshufd(vd, 0x31);
      t1 = pmuludq(t1, t0);
      t5 = paddq(t5, t2);

      t4 = psrlq(t4, 32);
      t5 = paddq(t5, t4);
      t5 = psrlq(t5, 32);
      t1 = paddq(t1, t5);
      t1 = psrlq(t1, 1);

      t5 = pshufd(t1, 0x31);
      t5 = pmuludq(t5, vm);
      t2 = pshufd(vm, 0x31);
      t2 = pmuludq(t2, t1);
      t1 = pmuludq(t1, vm);
      t5 = paddq(t5, t2);
      t5 = psllq(t5, 32);
      t1 = paddq(t1, t5);

      t3 = psubq(t3, t1);
      t3 = psubq(t3, vm);
      t1 = pcmpgtd(psetzero(), t3);
      t1 = pshufd(t1, 0xf5);
      t1 = pand(t1, vm);
      t3 = paddq(t3, t1);
      pstore((__m128i *)(x1 + i), t3);
    }

#else
  for (i = 0; i < len; i++)
    {
      sp_t t0 = x0[i];
      sp_t t1 = x1[i];
      sp_t t2, t3;
      t2 = sp_add (t0, t1, p);
      t3 = sp_sub (t0, t1, p);
      t3 = sp_mul (t3, w, p, d);
      x0[i] = t2;
      x1[i] = t3;
    }
#endif
}

static void
spv_ntt_dif_core (spv_t x, spv_t w, 
		  spv_size_t log2_len, sp_t p, sp_t d)
{
  spv_size_t len;
  spv_t x0, x1;
	
  /* handle small transforms immediately */
  switch (log2_len) {
   case 0:
    return;
   case 1:
    {
      sp_t t0 = x[0];
      sp_t t1 = x[1];
      x[0] = sp_add (t0, t1, p);
      x[1] = sp_sub (t0, t1, p);
      return;
    }
   case 2:
    {
      sp_t t0 = x[0];
      sp_t t1 = x[1];
      sp_t t2 = x[2];
      sp_t t3 = x[3];
      sp_t t4, t5, t6, t7;
      t4 = sp_add (t0, t2, p);
      t6 = sp_sub (t0, t2, p);
      t5 = sp_add (t1, t3, p);
      t7 = sp_sub (t1, t3, p);
      x[0] = sp_add (t4, t5, p);
      x[1] = sp_sub (t4, t5, p);
      t7 = sp_mul (t7, w[1], p, d);
      x[2] = sp_add (t6, t7, p);
      x[3] = sp_sub (t6, t7, p);
      return;
    }
   case 3:
    {
      sp_t t0 = x[0];
      sp_t t1 = x[1];
      sp_t t2 = x[2];
      sp_t t3 = x[3];
      sp_t t4 = x[4];
      sp_t t5 = x[5];
      sp_t t6 = x[6];
      sp_t t7 = x[7];
      sp_t t8, t9, t10, t11, t12, t13, t14, t15;

      t8 = sp_add (t0, t4, p);
      t12 = sp_sub (t0, t4, p);
      t9 = sp_add (t1, t5, p);
      t13 = sp_sub (t1, t5, p);
      t13 = sp_mul (t13, w[1], p, d);
      t10 = sp_add (t2, t6, p);
      t14 = sp_sub (t2, t6, p);
      t14 = sp_mul (t14, w[2], p, d);
      t11 = sp_add (t3, t7, p);
      t15 = sp_sub (t3, t7, p);
      t15 = sp_mul (t15, w[3], p, d);

      t0 = sp_add (t8, t10, p);
      t2 = sp_sub (t8, t10, p);
      t1 = sp_add (t9, t11, p);
      t3 = sp_sub (t9, t11, p);
      t3 = sp_mul (t3, w[2], p, d);
      x[0] = sp_add (t0, t1, p);
      x[1] = sp_sub (t0, t1, p);
      x[2] = sp_add (t2, t3, p);
      x[3] = sp_sub (t2, t3, p);

      t0 = sp_add (t12, t14, p);
      t2 = sp_sub (t12, t14, p);
      t1 = sp_add (t13, t15, p);
      t3 = sp_sub (t13, t15, p);
      t3 = sp_mul (t3, w[2], p, d);
      x[4] = sp_add (t0, t1, p);
      x[5] = sp_sub (t0, t1, p);
      x[6] = sp_add (t2, t3, p);
      x[7] = sp_sub (t2, t3, p);
      return;
    }
  }

  len = 1 << (log2_len - 1);
  x0 = x;
  x1 = x + len;
  bfly_dif (x0, x1, w, len, p, d);
  spv_ntt_dif_core (x0, w + len, log2_len - 1, p, d);
  spv_ntt_dif_core (x1, w + len, log2_len - 1, p, d);
}

static void
spv_vector_ntt_dif (spv_t x, spv_t w, spv_size_t log2_len, 
			spv_size_t vsize, spv_size_t stride,
			sp_t p, sp_t d)
{
  spv_size_t i, j;
  spv_size_t num_blocks = 1;
  spv_size_t num_bfly = 1 << (log2_len - 1);
  spv_size_t bfly_stride = stride << (log2_len - 1);

  while (num_bfly > 0)
    {
      spv_t x0 = x;
      spv_t x1 = x + bfly_stride;

      for (i = 0; i < num_blocks; i++)
	{
	  bfly(x0, x1, vsize, p);

	  for (j = 1; j < num_bfly; j++)
	    {
	      x0 += stride;
	      x1 += stride;
	      bfly_dif_sp(x0, x1, w[j], vsize, p, d);
	    }

	  x0 = x1 + stride;
	  x1 = x0 + bfly_stride;
	}

      w += num_bfly;
      num_blocks *= 2;
      num_bfly /= 2;
      bfly_stride /= 2;
    }
}
	
static const int bitrev[16] =
{ 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15};

void
spv_ntt_gfp_dif (spv_t x, spv_size_t log2_len, spm_t data)
{
  sp_t p = data->sp;
  sp_t d = data->mul_c;

  if (log2_len <= NTT_GFP_TWIDDLE_DIF_BREAKOVER)
    { 
      spv_t w = data->nttdata->twiddle + 
	        data->nttdata->twiddle_size - (1 << log2_len);
      spv_ntt_dif_core (x, w, log2_len, p, d);
    }
  else
    {
      /* recursive version for data that
         doesn't fit in the L1 cache */

      spv_size_t i, j;
      spv_size_t log2_col_len = 4;
      spv_size_t col_len = 1 << log2_col_len;
      spv_size_t row_len = 1 << (log2_len - log2_col_len);
      spv_size_t block_size = MIN(row_len, MAX_NTT_BLOCK_SIZE);
      spv_size_t stride = row_len;

      sp_t root = data->nttdata->ntt_roots[log2_len];
      spv_t col_w = data->nttdata->twiddle + 
	            data->nttdata->twiddle_size - col_len;
      spv_t w0 = data->scratch1;
      spv_t w1 = data->scratch2;

      for (i = 1, w0[0] = 1; i < block_size; i++)
       	w0[i] = sp_mul (w0[i-1], root, p, d);

      root = sp_pow (root, block_size, p, d);

      for (i = 0; i < row_len; i += block_size)
	{
	  spv_t x0 = x + i;

	  spv_vector_ntt_dif(x0, col_w, log2_col_len, 
	      		block_size, stride, p, d);

	  spv_set(w1, w0, block_size);
	  for (j = 1; j < col_len; j++)
	    {
	      spv_pwmul(x0 + bitrev[j] * stride, 
		        x0 + bitrev[j] * stride, 
			w1, block_size, p, d);

	      if (j < col_len - 1)
		spv_pwmul(w1, w1, w0, block_size, p, d);
	    }

	  if (i + block_size < row_len)
	    spv_mul_sp(w0, w0, root, block_size, p, d);
	}

      for (i = 0; i < col_len; i++)
	{
          spv_ntt_gfp_dif (x, log2_len - log2_col_len, data);
	  x += stride;
	}
    }
}

/*--------------------------- INVERSE NTT --------------------------------*/
static inline void bfly_dit(spv_t x0, spv_t x1, spv_t w,
				spv_size_t len, sp_t p, sp_t d)
{
  spv_size_t i;

#if defined(HAVE_SSE2) && (SP_NUMB_BITS < 32)

  __m128i t0, t1, t2, t3, vm, vm2, vd;

  vm = pshufd(pcvt_i32(p), 0x00);
  vm2 = pshufd(pcvt_i32(p), 0x44);
  vd = pshufd(pcvt_i32(d), 0x00);

  for (i = 0; i < len; i += 4)
    {
      t0 = pload((__m128i *)(x1 + i));
      t1 = pload((__m128i *)(w + i));
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

      t1 = pload((__m128i *)(x0 + i));
      t2 = paddd(t1, t0);
      t3 = psubd(t1, t0);
      t2 = psubd(t2, vm);

      t0 = pcmpgtd(psetzero(), t2);
      t1 = pcmpgtd(psetzero(), t3);
      t0 = pand(t0, vm);
      t1 = pand(t1, vm);
      t2 = paddd(t2, t0);
      t3 = paddd(t3, t1);

      pstore((__m128i *)(x0 + i), t2);
      pstore((__m128i *)(x1 + i), t3);
    }

#elif defined(HAVE_SSE2) && (GMP_LIMB_BITS == 32)

  __m128i t0, t1, t2, t3, t4, t5, vm, vd;

  vm = pshufd(pcvt_i64(p), 0x44);
  vd = pshufd(pcvt_i64(d), 0x44);

  for (i = 0; i < len; i += 2)
    {
      t4 = pload((__m128i *)(x1 + i));
      t5 = pload((__m128i *)(w + i));

      t3 = pshufd(t4, 0x31);
      t3 = pmuludq(t3, t5);
      t2 = pshufd(t5, 0x31);
      t2 = pmuludq(t2, t4);
      t1 = pshufd(t4, 0x31);
      t4 = pmuludq(t4, t5);
      t5 = pshufd(t5, 0x31);
      t5 = pmuludq(t5, t1);
      t3 = paddq(t2, t3);

      t0 = t4;
      t4 = psrlq(t4, 32);
      t1 = t3;
      t3 = psllq(t3, 32);
      t3 = paddq(t3, t0);
      t4 = paddq(t4, t1);
      t4 = psrlq(t4, 32);
      t4 = paddq(t4, t5);

      t0 = t3;
      t4 = psllq(t4, 2*(SP_TYPE_BITS - SP_NUMB_BITS));
      t0 = psrlq(t0, 2*SP_NUMB_BITS - SP_TYPE_BITS);
      t4 = paddq(t4, t0);

      t5 = pshufd(t4, 0x31);
      t5 = pmuludq(t5, vd);
      t2 = pshufd(vd, 0x31);
      t2 = pmuludq(t2, t4);
      t1 = pshufd(t4, 0x31);
      t4 = pmuludq(t4, vd);
      t0 = pshufd(vd, 0x31);
      t1 = pmuludq(t1, t0);
      t5 = paddq(t5, t2);

      t4 = psrlq(t4, 32);
      t5 = paddq(t5, t4);
      t5 = psrlq(t5, 32);
      t1 = paddq(t1, t5);
      t1 = psrlq(t1, 1);

      t5 = pshufd(t1, 0x31);
      t5 = pmuludq(t5, vm);
      t2 = pshufd(vm, 0x31);
      t2 = pmuludq(t2, t1);
      t1 = pmuludq(t1, vm);
      t5 = paddq(t5, t2);
      t5 = psllq(t5, 32);
      t1 = paddq(t1, t5);

      t3 = psubq(t3, t1);
      t3 = psubq(t3, vm);
      t1 = pcmpgtd(psetzero(), t3);
      t1 = pshufd(t1, 0xf5);
      t1 = pand(t1, vm);
      t3 = paddq(t3, t1);

      t1 = pload((__m128i *)(x0 + i));
      t2 = paddq(t1, t3);
      t3 = psubq(t1, t3);
      t2 = psubq(t2, vm);

      t0 = pcmpgtd(psetzero(), t2);
      t1 = pcmpgtd(psetzero(), t3);
      t0 = pshufd(t0, 0xf5);
      t1 = pshufd(t1, 0xf5);
      t0 = pand(t0, vm);
      t1 = pand(t1, vm);
      t2 = paddq(t2, t0);
      t3 = paddq(t3, t1);

      pstore((__m128i *)(x0 + i), t2);
      pstore((__m128i *)(x1 + i), t3);
    }

#else
  for (i = 0; i < len; i++)
    {
      sp_t w0 = w[i];
      sp_t t0 = x0[i];
      sp_t t1 = x1[i];
      t1 = sp_mul (t1, w0, p, d);
      x0[i] = sp_add (t0, t1, p);
      x1[i] = sp_sub (t0, t1, p);
    }
#endif
}

static inline void bfly_dit_sp(spv_t x0, spv_t x1, sp_t w,
				spv_size_t len, sp_t p, sp_t d)
{
  /* same as bfly_dit except all butterflies are scaled by
     a single value of w */

  spv_size_t i;

#if defined(HAVE_SSE2) && (SP_NUMB_BITS < 32)

  __m128i t0, t1, t2, t3, vm, vm2, vd, vw;

  vw = pshufd(pcvt_i32(w), 0x00);
  vm = pshufd(pcvt_i32(p), 0x00);
  vm2 = pshufd(pcvt_i32(p), 0x44);
  vd = pshufd(pcvt_i32(d), 0x00);

  for (i = 0; i < len; i += 4)
    {
      t0 = pload((__m128i *)(x1 + i));
      t1 = vw;
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

      t1 = pload((__m128i *)(x0 + i));
      t2 = paddd(t1, t0);
      t3 = psubd(t1, t0);
      t2 = psubd(t2, vm);

      t0 = pcmpgtd(psetzero(), t2);
      t1 = pcmpgtd(psetzero(), t3);
      t0 = pand(t0, vm);
      t1 = pand(t1, vm);
      t2 = paddd(t2, t0);
      t3 = paddd(t3, t1);

      pstore((__m128i *)(x0 + i), t2);
      pstore((__m128i *)(x1 + i), t3);
    }

#elif defined(HAVE_SSE2) && (GMP_LIMB_BITS == 32)

  __m128i t0, t1, t2, t3, t4, t5, vm, vd, vw;

  vm = pshufd(pcvt_i64(p), 0x44);
  vd = pshufd(pcvt_i64(d), 0x44);
  vw = pshufd(pcvt_i64(w), 0x44);

  for (i = 0; i < len; i += 2)
    {
      t4 = pload((__m128i *)(x1 + i));
      t5 = vw;

      t3 = pshufd(t4, 0x31);
      t3 = pmuludq(t3, t5);
      t2 = pshufd(t5, 0x31);
      t2 = pmuludq(t2, t4);
      t1 = pshufd(t4, 0x31);
      t4 = pmuludq(t4, t5);
      t5 = pshufd(t5, 0x31);
      t5 = pmuludq(t5, t1);
      t3 = paddq(t2, t3);

      t0 = t4;
      t4 = psrlq(t4, 32);
      t1 = t3;
      t3 = psllq(t3, 32);
      t3 = paddq(t3, t0);
      t4 = paddq(t4, t1);
      t4 = psrlq(t4, 32);
      t4 = paddq(t4, t5);

      t0 = t3;
      t4 = psllq(t4, 2*(SP_TYPE_BITS - SP_NUMB_BITS));
      t0 = psrlq(t0, 2*SP_NUMB_BITS - SP_TYPE_BITS);
      t4 = paddq(t4, t0);

      t5 = pshufd(t4, 0x31);
      t5 = pmuludq(t5, vd);
      t2 = pshufd(vd, 0x31);
      t2 = pmuludq(t2, t4);
      t1 = pshufd(t4, 0x31);
      t4 = pmuludq(t4, vd);
      t0 = pshufd(vd, 0x31);
      t1 = pmuludq(t1, t0);
      t5 = paddq(t5, t2);

      t4 = psrlq(t4, 32);
      t5 = paddq(t5, t4);
      t5 = psrlq(t5, 32);
      t1 = paddq(t1, t5);
      t1 = psrlq(t1, 1);

      t5 = pshufd(t1, 0x31);
      t5 = pmuludq(t5, vm);
      t2 = pshufd(vm, 0x31);
      t2 = pmuludq(t2, t1);
      t1 = pmuludq(t1, vm);
      t5 = paddq(t5, t2);
      t5 = psllq(t5, 32);
      t1 = paddq(t1, t5);

      t3 = psubq(t3, t1);
      t3 = psubq(t3, vm);
      t1 = pcmpgtd(psetzero(), t3);
      t1 = pshufd(t1, 0xf5);
      t1 = pand(t1, vm);
      t3 = paddq(t3, t1);

      t1 = pload((__m128i *)(x0 + i));
      t2 = paddq(t1, t3);
      t3 = psubq(t1, t3);
      t2 = psubq(t2, vm);

      t0 = pcmpgtd(psetzero(), t2);
      t1 = pcmpgtd(psetzero(), t3);
      t0 = pshufd(t0, 0xf5);
      t1 = pshufd(t1, 0xf5);
      t0 = pand(t0, vm);
      t1 = pand(t1, vm);
      t2 = paddq(t2, t0);
      t3 = paddq(t3, t1);

      pstore((__m128i *)(x0 + i), t2);
      pstore((__m128i *)(x1 + i), t3);
    }

#else
  for (i = 0; i < len; i++)
    {
      sp_t t0 = x0[i];
      sp_t t1 = x1[i];
      t1 = sp_mul (t1, w, p, d);
      x0[i] = sp_add (t0, t1, p);
      x1[i] = sp_sub (t0, t1, p);
    }
#endif
}

static void
spv_ntt_dit_core (spv_t x, spv_t w, 
		  spv_size_t log2_len, sp_t p, sp_t d)
{
  spv_size_t len;
  spv_t x0, x1;
	
  /* handle small transforms immediately */
  switch (log2_len) {
   case 0:
    return;
   case 1:
    {
      sp_t t0 = x[0];
      sp_t t1 = x[1];
      x[0] = sp_add (t0, t1, p);
      x[1] = sp_sub (t0, t1, p);
      return;
    }
   case 2:
    {
      sp_t t0 = x[0];
      sp_t t1 = x[1];
      sp_t t2 = x[2];
      sp_t t3 = x[3];
      sp_t t4, t5, t6, t7;
      t4 = sp_add (t0, t1, p);
      t5 = sp_sub (t0, t1, p);
      t6 = sp_add (t2, t3, p);
      t7 = sp_sub (t2, t3, p);
      x[0] = sp_add (t4, t6, p);
      x[2] = sp_sub (t4, t6, p);
      t7 = sp_mul (t7, w[1], p, d);
      x[1] = sp_add (t5, t7, p);
      x[3] = sp_sub (t5, t7, p);
      return;
    }
   case 3:
    {
      sp_t t0 = x[0];
      sp_t t1 = x[1];
      sp_t t2 = x[2];
      sp_t t3 = x[3];
      sp_t t4 = x[4];
      sp_t t5 = x[5];
      sp_t t6 = x[6];
      sp_t t7 = x[7];
      sp_t t8, t9, t10, t11;

      t8 = sp_add(t0, t1, p);
      t9 = sp_sub(t0, t1, p);
      t10 = sp_add(t2, t3, p);
      t11 = sp_sub(t2, t3, p);
      t0 = sp_add(t8, t10, p);
      t2 = sp_sub(t8, t10, p);
      t11 = sp_mul (t11, w[2], p, d);
      t1 = sp_add(t9, t11, p);
      t3 = sp_sub(t9, t11, p);

      t8 = sp_add(t4, t5, p);
      t9 = sp_sub(t4, t5, p);
      t10 = sp_add(t6, t7, p);
      t11 = sp_sub(t6, t7, p);
      t4 = sp_add(t8, t10, p);
      t6 = sp_sub(t8, t10, p);
      t11 = sp_mul (t11, w[2], p, d);
      t5 = sp_add(t9, t11, p);
      t7 = sp_sub(t9, t11, p);

      x[0] = sp_add(t0, t4, p);
      x[4] = sp_sub(t0, t4, p);
      t5 = sp_mul (t5, w[1], p, d);
      x[1] = sp_add(t1, t5, p);
      x[5] = sp_sub(t1, t5, p);
      t6 = sp_mul (t6, w[2], p, d);
      x[2] = sp_add(t2, t6, p);
      x[6] = sp_sub(t2, t6, p);
      t7 = sp_mul (t7, w[3], p, d);
      x[3] = sp_add(t3, t7, p);
      x[7] = sp_sub(t3, t7, p);
      return;
    }
  }

  len = 1 << (log2_len - 1);
  x0 = x;
  x1 = x + len;
  spv_ntt_dit_core (x0, w + len, log2_len - 1, p, d);
  spv_ntt_dit_core (x1, w + len, log2_len - 1, p, d);
  bfly_dit (x0, x1, w, len, p, d);
}

static void
spv_vector_ntt_dit (spv_t x, spv_t w, spv_size_t log2_len, 
			spv_size_t vsize, spv_size_t stride,
			sp_t p, sp_t d)
{
  spv_size_t i, j;
  spv_size_t num_bfly = 1;
  spv_size_t num_blocks = 1 << (log2_len - 1);
  spv_size_t bfly_stride = stride;

  w += (1 << log2_len) - 2;
  while (num_blocks > 0)
    {
      spv_t x0 = x;
      spv_t x1 = x + bfly_stride;

      for (i = 0; i < num_blocks; i++)
	{
	  bfly(x0, x1, vsize, p);

	  for (j = 1; j < num_bfly; j++)
	    {
	      x0 += stride;
	      x1 += stride;
	      bfly_dit_sp(x0, x1, w[j], vsize, p, d);
	    }

	  x0 = x1 + stride;
	  x1 = x0 + bfly_stride;
	}

      num_blocks /= 2;
      num_bfly *= 2;
      bfly_stride *= 2;
      w -= num_bfly;
    }
}
	
void
spv_ntt_gfp_dit (spv_t x, spv_size_t log2_len, spm_t data)
{
  sp_t p = data->sp;
  sp_t d = data->mul_c;

  if (log2_len <= NTT_GFP_TWIDDLE_DIT_BREAKOVER)
    {
      spv_t w = data->inttdata->twiddle + 
	        data->inttdata->twiddle_size - (1 << log2_len);
      spv_ntt_dit_core (x, w, log2_len, p, d);
    }
  else
    {
      spv_size_t i, j;
      spv_t x0 = x;
      spv_size_t log2_col_len = 4;
      spv_size_t col_len = 1 << log2_col_len;
      spv_size_t row_len = 1 << (log2_len - log2_col_len);
      spv_size_t block_size = MIN(row_len, MAX_NTT_BLOCK_SIZE);
      spv_size_t stride = row_len;

      sp_t root = data->inttdata->ntt_roots[log2_len];
      spv_t col_w = data->inttdata->twiddle + 
	            data->inttdata->twiddle_size - col_len;
      spv_t w0 = data->scratch1;
      spv_t w1 = data->scratch2;

      for (i = 0; i < col_len; i++)
	{
	  spv_ntt_gfp_dit(x0, log2_len - log2_col_len, data);
	  x0 += stride;
	}


      for (i = 1, w0[0] = 1; i < block_size; i++)
       	w0[i] = sp_mul (w0[i-1], root, p, d);

      root = sp_pow (root, block_size, p, d);

      for (i = 0; i < row_len; i += block_size)
	{
	  x0 = x + i;

	  spv_set(w1, w0, block_size);
	  for (j = 1; j < col_len; j++)
	    {
	      spv_pwmul(x0 + bitrev[j] * stride, 
		        x0 + bitrev[j] * stride, 
			w1, block_size, p, d);

	      if (j < col_len - 1)
		spv_pwmul(w1, w1, w0, block_size, p, d);
	    }

	  if (i + block_size < row_len)
	    spv_mul_sp(w0, w0, root, block_size, p, d);

	  spv_vector_ntt_dit(x0, col_w, log2_col_len, 
	      		block_size, stride, p, d);

	}
    }
}
