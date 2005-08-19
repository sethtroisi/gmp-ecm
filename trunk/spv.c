/* spv.c - "small prime vector" functions for arithmetic on vectors of
   residues modulo a single small prime

  Copyright 2005 Dave Newman.

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
  the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
  MA 02111-1307, USA.
*/

#include <string.h> /* for memset */
#include "sp.h"

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
  spv_size_t i;
  
  ASSERT (r >= x + len || x >= r);
  ASSERT (r >= y + len || y >= r);

  for (i = 0; i < len; i++)
    r[i] = sp_mul (x[i], y[i], m, d);
}

/* dst = src * y */
void
spv_mul_sp (spv_t r, spv_t x, sp_t c, spv_size_t len, sp_t m, sp_t d)
{
  spv_size_t i;
  
  ASSERT (r >= x + len || x >= r);
  
  for (i = 0; i < len; i++)
    r[i] = sp_mul (x[i], c, m, d);
}

#if 0
/* r += x * y */
void
spv_addmul_sp (spv_t r, spv_t x, sp_t c, spv_size_t len, sp_t m, sp_t d)
{
  spv_size_t i;
  sp_t t;
	
  ASSERT (r >= x + len || x >= r);
  
  for (i = 0; i < len; i++)
  {
    t = sp_mul (x[i], c, m, d);
    r[i] = sp_add (r[i], t, m);
  }
}

/* r -= x * y */
void
spv_submul_sp (spv_t r, spv_t x, sp_t c, spv_size_t len, sp_t m, sp_t d)
{
  spv_size_t i;
  sp_t t;
  
  ASSERT (r >= x + len || x >= r);
  
  for (i = 0; i < len; i++)
  {
    t = sp_mul (x[i], c, m, d);
    r[i] = sp_sub (x[i], t, m);
  }
}

/* r = x * y by grammar-school polynomial multiplication
 * 
 * r must be distinct from both x and y
 * x_len > 0, y_len > 0 */
void
spv_mul_basecase (spv_t r, spv_t x, spv_t y, spv_size_t x_len,
    spv_size_t y_len, sp_t m, sp_t d)
{
  spv_size_t i;
  
  ASSERT (r >= x + x_len || x >= r + x_len + y_len - 1);
  ASSERT (r >= y + y_len || y >= r + x_len + y_len - 1);
  ASSERT (x_len > 0);
  ASSERT (y_len > 0);
  
  if (x_len > y_len)
  {
    spv_mul_sp (r, x, y[0], x_len, m, d);
    spv_set_sp (r + x_len, 0, y_len - 1);
    
    for (i = 1; i < y_len; i++)
      spv_addmul_sp (r + i, x, y[i], x_len, m, d);
  }
  else
  {
    spv_mul_sp (r, y, x[0], y_len, m, d);
    spv_set_sp (r + y_len, 0, x_len - 1);
    
    for (i = 1; i < x_len; i++)
      spv_addmul_sp (r + i, y, x[i], y_len, m, d);
  }
}

/* dst = src1 * src2 by karatsuba multiplication
 * 
 * dst must be distinct from both src1 and src2
 * src1 and src2 have the same len > 0 
 * t is a temporary array which must be at least M(len) large, where
 *
 * M(1)=0, M(K) = max(3*l-1,2*l-2+M(l)) <= 2*K-1 where l = ceil(K/2).
 * 
 * Code adapted from the mpz_t karatsuba in gmp-ecm, in which multiplies
 * are substantially more expensive than additions. This is also true for
 * sp's, but to a lesser extent, so it might be the case that minimising the
 * total number of operations is more important than minimising the number
 * of muls. */

void
spv_mul_karatsuba (spv_t r, spv_t x, spv_t y, spv_t t, spv_size_t len,
    sp_t m, sp_t d)
{
  spv_size_t i, k, l;
  spv_t z;
  
  ASSERT (r >= x + len || x >= r + 2 * len - 1);
  ASSERT (r >= y + len || y >= r + 2 * len - 1);
  ASSERT (len > 0);
  /* FIXME: add assertions for t */
  
  if (len == 1)
    {
      r[0] = sp_mul (x[0], y[0], m, d);
      return;
    }
  if (len == 2)
    {
      t[0] = sp_add (x[0], x[1], m);
      r[1] = sp_add (y[0], y[1], m);
      r[1] = sp_mul (r[1], t[0], m, d);
      r[0] = sp_mul (x[0], y[0], m, d);
      r[2] = sp_mul (x[1], y[1], m, d);
      r[1] = sp_sub (r[1], r[0], m);
      r[1] = sp_sub (r[1], r[2], m);
      return;
    }
  if (len == 3)
    {
      r[0] = sp_mul (x[0], y[0], m, d);
      r[2] = sp_mul (x[1], y[1], m, d);
      r[4] = sp_mul (x[2], y[2], m, d);
      t[0] = sp_add (x[0], x[1], m);
      t[1] = sp_add (y[0], y[1], m);
      r[1] = sp_mul (t[0], t[1], m, d);
      r[1] = sp_sub (r[1], r[0], m);
      r[1] = sp_sub (r[1], r[2], m);
      t[0] = sp_add (x[1], x[2], m);
      t[1] = sp_add (y[1], y[2], m);
      r[3] = sp_mul (t[0], t[1], m, d);
      r[3] = sp_sub (r[3], r[2], m);
      r[3] = sp_sub (r[3], r[4], m);
      t[0] = sp_add (x[0], x[2], m);
      t[1] = sp_add (y[0], y[2], m);
      t[2] = sp_mul (t[0], t[1], m, d);
      t[2] = sp_sub (t[2], r[0], m);
      t[2] = sp_sub (t[2], r[4], m);
      r[2] = sp_add (r[2], t[2], m);
      return;
    }
  
  k = len / 2;
  l = len - k;

  z = t + 2 * l - 1;
  
  for (i = 0; i < k; i++)
    {
      z[i] = sp_sub (x[i], x[l + i], m);
      r[i] = sp_sub (y[i], y[l + i], m);
    }

  if (l > k)
    {
      z[k] = x[k];
      r[k] = y[k];
    }

  spv_mul_karatsuba (t, z, r, r + l, l, m, d);
       
  z = t + 2 * l - 2;
  r[2 * l - 1] = t[2 * l - 2];

  spv_mul_karatsuba (r, x, y, z, l, m, d);
  spv_mul_karatsuba (r + 2 * l, x + l, y + l, z, k, m, d);

  t[2 * l - 2] = r[2 * l - 1];
  r[2 * l - 1] = 0;

  spv_add (r + 2 * l, r + 2 * l, r + l, l - 1, m);
  if (k > 1)
    {
      spv_add (r + l, r + 2 * l, r, l, m);
      spv_add (r + 2 * l, r + 2 * l, r + 3 * l, 2 * k - 1 - l, m);
    }
  else
    {
      r[l] = sp_add (r[2 * l], r[0], m);
      if (len == 3)
        r[l + 1] = r[1];
    }

  spv_sub (r + l, r + l, t, 2 * l - 1, m);
}

/* calculate r[k], ..., r[l - 1] of the product r = x * y
 *
 * other coeffs in r are undefined
 * r has size at least the next power of two >= prod_len
 * allow l == 0 if the full product is required */
void spv_mul (spv_t r, spv_t x, spv_size_t x_len, spv_t y,
    spv_size_t y_len, spv_size_t k, spv_size_t l, int monic, spm_t spm)
{
  /* FIXME: add assertions */
  
  if (l == 0)
    l = x_len + y_len - 1 + monic;
  
  ASSERT (x_len > 0);
  ASSERT (y_len > 0);
  ASSERT (monic == 0 || monic == 1);
  ASSERT (k < l);
  ASSERT (l <= x_len + y_len - 1 + monic);
  
  x_len = MIN (l, x_len);
  y_len = MIN (l, y_len);
  
  int square = (x == y && x_len == y_len);
  int equal_op;
  spv_size_t prod_len = x_len + y_len - 1 + monic;
  spv_size_t max_len = MAX (x_len, y_len);
  
  /* ensure either x != r != y or r == x */
  if (r == y && r != x)
    {
      spv_t t = y;
      spv_size_t t_len = y_len;
      
      y = x; y_len = x_len;
      x = t; x_len = t_len;
    }
  
  equal_op = (r == x);
  
  if (prod_len < MUL_NTT_THRESHOLD)
    {
      spv_t t = (spv_t) malloc (2 * max_len * sizeof (sp_t));
      
      /* the original contents of x, y, these are respectively
       * either x, y themselves or copies */
      spv_t x0, y0;
      
      /* we cannot rely on x or y being large enough to allow us
       * to zero-pad them to max_len so if x_len != y_len then
       * we copy the smaller vector */
      
      if (equal_op || x_len < y_len)
        {
	  /* karatsuba implementation requires r != x, r != y */
	  x0 = (spv_t) malloc (max_len * sizeof (sp_t));
	  
	  spv_set (x0, x, x_len);
	  spv_set_zero (x0 + x_len, max_len - x_len);
	}
      else
        x0 = x;
      
      if (y_len < x_len)
        {
	  y0 = (spv_t) malloc (max_len * sizeof (sp_t));

	  spv_set (y0, y, y_len);
	  spv_set_zero (y0 + y_len, max_len - y_len);
	}
      else
	y0 = square ? x0 : y;
	  
      spv_mul_karatsuba (r, x0, y0, t, max_len, spm->sp, spm->mul_c);
      
      free (t);

      if (monic)
        {
	  /* FIXME crop to range [k, l) */
	  r[prod_len - 1] = 0;
	  spv_add (r + x_len, r + x_len, y0, y_len, spm->sp);
	  spv_add (r + y_len, r + y_len, x0, x_len, spm->sp);
	}
      if (equal_op || x_len < y_len)
	free (x0);
      if (y_len < x_len)
	free (y0);
      return;
    }
  
  spv_t x_ntt, y_ntt;
  spv_size_t i, ntt_size;
  
  ntt_size = 1 << ceil_log_2 (prod_len);
  
  /* threshold seems to give reasonably
   * good results but needs fine-tuning properly */
  if (prod_len < 3 * ntt_size / 4 /* && ntt_size / 2 >= x_len && ntt_size / 2 >= y_len */)
    {
      ntt_size >>= 1;
      
      if (prod_len - ntt_size > k)
	/* high/middle coeffs overflow into
	 * middle/low coeffs */
	spv_mul (r + ntt_size,
	    x, prod_len - ntt_size,
	    y, prod_len - ntt_size,
	    0, prod_len - ntt_size, 0, spm);
      else if (l > ntt_size)
	/* high/middle coeffs overflow into low coeffs */
	spv_mul (r + ntt_size,
	  x, l - ntt_size,
	  y, l - ntt_size,
	  0, l - ntt_size, 0, spm);
      /* else */
	/* high coeffs overflow into low coeffs */
    }
  
  /*  printf ("x_len = %u, y_len = %u, k = %u, l = %u, monic = %u, prod_len = %u, ntt_size = %u\n",
    x_len, y_len, k, l, monic, prod_len, ntt_size); */
  
  /* variable x_ntt only exists to help readability */
  x_ntt = r;
  
  if (!equal_op)
    spv_set (x_ntt, x, x_len);
  
  spv_set_zero (x_ntt + x_len, ntt_size - x_len);
  
  if (square)
    {
      if (monic)
	x_ntt[x_len] = 1;
      spv_sqr_ntt_gfp (r, x_ntt, ntt_size, spm);
    }
  else
    {
      y_ntt = (spv_t) malloc (ntt_size * sizeof (sp_t));
      spv_set (y_ntt, y, y_len);
      spv_set_zero (y_ntt + y_len, ntt_size - y_len);
    
      if (monic)
        x_ntt[x_len] = y_ntt[y_len] = 1;
      
      spv_mul_ntt_gfp (r, x_ntt, y_ntt, ntt_size, spm);
      free (y_ntt);
    }
  
  if (prod_len > k + ntt_size)
    {
      /* FIXME maybe deal with these cases separately */
      
      /* ####............
       *         ######## */

      /* ############....
       *     ....######## */

      /* ########........
       *     ........#### */

      /* ....########....
       *         ........ */
      
      for (i = 0; i < prod_len - ntt_size; i++)
        {
	  sp_t u = r[ntt_size + i];
	  r[ntt_size + i] = sp_sub (r[i], u, spm->sp);
	  r[i] = u;
	}
    }
  else if (l > ntt_size)
    {
      /* ########........
       *             #### */
      for (i = 0; i < l - ntt_size; i++)
        r[ntt_size + i] = sp_sub (r[i], r[ntt_size + i], spm->sp);
    }
  
  if (monic)
    /* everything is correct except for the x^{prod_len} term which may
     * have wrapped round */
    r[prod_len % ntt_size] = sp_sub (r[prod_len % ntt_size], 1, spm->sp);
}
#endif

void
spv_random (spv_t x, spv_size_t len, sp_t m)
{
  spv_size_t i;
  mpn_random (x, len);
  
  for (i = 0; i < len; i++)
    if (x[i] >= m)
      x[i] -= m;
}
