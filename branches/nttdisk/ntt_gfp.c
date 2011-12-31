/* ntt_gfp.c - low-level radix-2 dif/dit ntt routines over GF(p)
   
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

#include "sp.h"
#include "ecm-impl.h"

/*--------------------------- FORWARD NTT --------------------------------*/
static void bfly_dif(spv_t x0, spv_t x1, spv_t w,
			spv_size_t len, sp_t p, sp_t d)
{
  spv_size_t i = 0;

#if defined(HAVE_SSE2) && (SP_NUMB_BITS < 32)

  __m128i t0, t1, t2, t3, vm, vm2, vd;

  vm = pshufd(pcvt_i32(p), 0x00);
  vm2 = pshufd(pcvt_i32(p), 0x44);
  vd = pshufd(pcvt_i32(d), 0x00);

  for (; i < len; i += 4)
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
      spv_size_t len = 1 << (log2_len - 1);
      spv_t x0 = x;
      spv_t x1 = x + len;
      spv_t roots = data->nttdata->ntt_roots;

        {
          spv_size_t i;
	  spv_size_t block_size = MIN(len, MAX_NTT_BLOCK_SIZE);
          sp_t root = roots[log2_len];
	  spv_t w = data->scratch;

	  w[0] = 1;
	  for (i = 1; i < block_size; i++)
	    w[i] = sp_mul (w[i-1], root, p, d);

          root = sp_pow (root, block_size, p, d);

	  for (i = 0; i < len; i += block_size)
	    {
	      if (i)
	        spv_mul_sp (w, w, root, block_size, p, d);

	      bfly_dif (x0 + i, x1 + i, w, block_size, p, d);
	    }
	}
	
      spv_ntt_gfp_dif (x0, log2_len - 1, data);
      spv_ntt_gfp_dif (x1, log2_len - 1, data);
    }
}

/*--------------------------- INVERSE NTT --------------------------------*/
static inline void bfly_dit(spv_t x0, spv_t x1, spv_t w,
				spv_size_t len, sp_t p, sp_t d)
{
  spv_size_t i = 0;

#if defined(HAVE_SSE2) && (SP_NUMB_BITS < 32)

  __m128i t0, t1, t2, t3, vm, vm2, vd;

  vm = pshufd(pcvt_i32(p), 0x00);
  vm2 = pshufd(pcvt_i32(p), 0x44);
  vd = pshufd(pcvt_i32(d), 0x00);

  for (; i < len; i += 4)
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
      spv_size_t len = 1 << (log2_len - 1);
      spv_t x0 = x;
      spv_t x1 = x + len;
      spv_t roots = data->inttdata->ntt_roots;

      spv_ntt_gfp_dit (x0, log2_len - 1, data);
      spv_ntt_gfp_dit (x1, log2_len - 1, data);

        {
          spv_size_t i;
	  spv_size_t block_size = MIN(len, MAX_NTT_BLOCK_SIZE);
          sp_t root = roots[log2_len];
	  spv_t w = data->scratch;

	  w[0] = 1;
	  for (i = 1; i < block_size; i++)
	    w[i] = sp_mul (w[i-1], root, p, d);

          root = sp_pow (root, block_size, p, d);

	  for (i = 0; i < len; i += block_size)
	    {
	      if (i)
	        spv_mul_sp (w, w, root, block_size, p, d);

	      bfly_dit (x0 + i, x1 + i, w, block_size, p, d);
	    }
	}
    }
}
