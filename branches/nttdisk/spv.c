/* spv.c - "small prime vector" functions for arithmetic on vectors of
   residues modulo a single small prime

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

#include "config.h"
#include <string.h> /* for memset */
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#ifdef USE_VALGRIND
#include <valgrind/memcheck.h>
#endif
#include "ecm-impl.h"

/* Routines for vectors of integers modulo r common small prime
 * 
 * These are low-overhead routines that don't do memory allocation,
 * other than for temporary variables. Unless otherwise specified, any
 * of the input pointers can be equal. */


/* Test that an spv_t is valid input. It tests that the memory is allocated
   and contains initialised data (if compiled with USE_VALGRIND) and that
   all the sp_t values are less than the modulus. */
int
spv_verify_in (const spv_tc x, const spv_size_t len, const sp_t m)
{
  spv_size_t i;
#ifdef USE_VALGRIND
  if (len > 0)
    {
      if (VALGRIND_CHECK_MEM_IS_DEFINED(x, len * sizeof(sp_t)) != 0)
        return 0;
    }
#endif
  for (i = 0; i < len; i++)
    {
      if (x[i] >= m)
        return 0;
    }
  
  return 1;
}

/* Test that an spv_t is valid output. It tests that the memory is allocated
   (if compiled with USE_VALGRIND). */
int
spv_verify_out (ATTRIBUTE_UNUSED const spv_tc x, 
                ATTRIBUTE_UNUSED const spv_size_t len, 
                ATTRIBUTE_UNUSED const sp_t m)
{
#ifdef USE_VALGRIND
  if (len > 0)
    {
      if (VALGRIND_CHECK_MEM_IS_ADDRESSABLE(x, len * sizeof(sp_t)) != 0)
        return 0;
    }
#endif
  
  return 1;
}


/* Wrapper function for element-wise operations. 
   Input and output operands can be on disk. */

void 
spv_elementwise(
  spv_t r, FILE *r_file, const spv_size_t r_offset, 
  spv_tc x, FILE *x_file, const spv_size_t x_offset, 
  spv_tc y, FILE *y_file, const spv_size_t y_offset, 
  const sp_t p, const sp_t d, const spv_size_t len, const int operation)
{
  spv_size_t buf_len, done = 0;
  spv_t buf_r, buf[2];
  spv_tc buf_x, buf_y;
  int xyf, rx, ry;

  ASSERT(r != NULL || r_file != NULL);
  ASSERT(r == NULL || r_file == NULL);
  ASSERT(x == NULL || x_file == NULL);
  ASSERT(y == NULL || y_file == NULL);

  /* We allow no overlap between input and output except identical data */
  ASSERT (r == NULL || x == NULL || r + r_offset == x + x_offset || 
          r + r_offset + len <= x + x_offset || 
          r + r_offset >= x + x_offset + len);
  ASSERT (r_file == NULL || x_file == NULL || r_file != x_file || 
          r_offset == x_offset || r_offset + len <= x_offset || 
          r_offset >= x_offset + len);
  ASSERT (r == NULL || y == NULL || r + r_offset == y + y_offset || 
          r + r_offset + len <= y + y_offset || 
          r + r_offset >= y + y_offset + len);
  ASSERT (r_file == NULL || y_file == NULL || r_file != y_file || 
          r_offset == y_offset || r_offset + len <= y_offset || 
          r_offset >= y_offset + len);

  if (operation != SPV_ELEMENTWISE_RANDOM)
    {
      ASSERT(x != NULL || x_file != NULL);
    }
  if (operation == SPV_ELEMENTWISE_SET || 
      operation == SPV_ELEMENTWISE_NEG || operation == SPV_ELEMENTWISE_SETSP)
    {
      ASSERT(y == NULL && y_file == NULL);
    }
  else if (operation == SPV_ELEMENTWISE_ADDSP || 
           operation == SPV_ELEMENTWISE_SUBSP)
    {
      /* The single sp_t that gets added/subtracted must be in memory */
      ASSERT(y != NULL && y_file == NULL);
    }
  else
    {
      ASSERT(y != NULL || y_file != NULL);
    }


  /*
      storage r x y:  0 = in memory, 1 = on disk
      Alias rx ry xy: whether any input pointer/file and offset of one input 
                      operand are equal to those of another  
      Alloc bufs?:    how many temp buffers we need to allocate
      Data in r x y:  Where we put the data for processing. 
                      r,x,y = input vector r,x,y, resp
                      0 = allocated buffer 0, 1 = allocated buffer 1
      xyf = (ON_DISK(x) && ON_DISK(y) && x != y). 
           This simplifies the conditions for where to store buf_y
      
      storage	Alias	 Alloc	Data in	xyf	Notes
      r x y 	rx ry xy bufs?	r x y		
      -------------------------------------------------------------------
      0 0 0	any	 no	r x y	0	work in-place

      0 0 1	0 0 0	 no	r x r	0	read y-data to r
      0 0 1	1 0 0	 1 buf	r x 0	0	Can't read y to r or we'd overwrite x-data, so read y-data to buffer 0

      0 1 0	0 0 0	 no	r r y	0	read x-data to r
      0 1 0	0 1 0	 1 buf	r 0 y	0	Can't read x-data to r or we'd overwrite y-data, so read x-data to buffer 0

      0 1 1	0 0 0	 1 buf	r r 0	1	read x-data to r, y-data to buffer 0
      0 1 1	0 0 1	 no	r r r	0	read x-data to r, use as y-data as well

      1 0 0	0 0 0	 1 buf	0 x y	0	write r-data from buffer 0
      1 0 0	0 0 1	 1 buf	0 x y	0	write r-data from buffer 0

      1 0 1	0 0 0	 1 buf	0 x 0	0	read y-data to buffer 0, write r-data from buffer 0
      1 0 1	0 1 0	 1 buf	0 x 0	0	read y-data to buffer 0, write r-data from buffer 0

      1 1 0	0 0 0	 1 buf	0 0 y	0	read x-data to buffer 0, write r-data from buffer 0
      1 1 0	1 0 0	 1 buf	0 0 y	0	read x-data to buffer 0, write r-data from buffer 0

      1 1 1	0 0 0	 2 bufs	0 0 1	1	read x-data to buffer 0, read y-data to buffer 1, write r-data from buffer 0
      1 1 1	1 0 0	 2 bufs	0 0 1	1	read x-data to buffer 0, read y-data to buffer 1, write r-data from buffer 0
      1 1 1	0 1 0	 2 bufs	0 0 1	1	read x-data to buffer 0, read y-data to buffer 1, write r-data from buffer 0
      1 1 1	0 0 1	 1 bufs	0 0 0	0	read x-data to buffer 0, use as y-data as well, write r-data from buffer 0
      1 1 1	1 1 1	 1 buf	0 0 0	0	read x-data to buffer 0, use as y-data as well, write r-data from buffer 0

      nr_buffers > 0  <==>  ON_DISK(r) || (ON_DISK(x) || ON_DISK(y)) && (r == x || r == y) || ON_DISK(y) && xyf
      nr_buffers = 2  <==>  y = buf[1]  <==>  ON_DISK(r) && ON_DISK(y) && xyf

      buf_r = (!ON_DISK(r)) ? r : buf[0];
      buf_x = (!ON_DISK(x)) ? x : (r == y) ? buf[0] : buf_r;
      buf_y = (!ON_DISK(y)) ? y : 
              Case: ON_DISK(y) == 1 && ON_DISK(r) == 0
                (x == r || xyf) ? buf[0]: buf_r
              Case: ON_DISK(y) == 1 && ON_DISK(r) == 1
               xyf ? buf[1] : buf_r;
  */

  /* xyf is non-zero iff x and y are both on disk but don't refer to the same 
     data, i.e., are in separate files or at different offsets */
  xyf = (x_file != NULL && y_file != NULL && 
         (x_file != y_file || x_offset != y_offset));
  /* Non-zero iff r and x (resp. y) are both in memory and refer 
     to the same data */
  rx = (r != NULL && x != NULL && r == x && r_offset == x_offset);
  ry = (r != NULL && y != NULL && r == y && r_offset == y_offset);

  buf_len = len;
  buf[0] = buf[1] = NULL;
  if (r_file != NULL || 
      ((x_file != NULL || y_file != NULL) && (rx || ry)) || 
      xyf)
    {
      buf_len = MIN (65536, len);
      buf[0] = (spv_t) sp_aligned_malloc (buf_len * sizeof (sp_t));
    }

  if (r_file != NULL && xyf)
    buf[1] = (spv_t) sp_aligned_malloc (buf_len * sizeof (sp_t));

  buf_r = (r_file == NULL) ? (r + r_offset) : buf[0];
  buf_x = (x_file == NULL) ? (x + x_offset) : (ry) ? buf[0] : buf_r;
  buf_y = (y_file == NULL) ? (y + y_offset) : (r_file != NULL && xyf) ? buf[1] : 
            (r_file == NULL && (rx || xyf)) ? buf[0] : buf_r;
  
  while (done < len)
    {
      spv_size_t do_now = MIN(buf_len, len - done);
      if (x_file != NULL)
        {
          /* This is stupid but buf_x is a const sp_t *, so we can't use this 
             pointer as an argument to fread(). We need to find a non-const 
             pointer which we can use instead */
          spv_t t = (buf_x == r) ? r : (buf_x == buf[0]) ? buf[0] : NULL;
          spv_seek_and_read (t, do_now, x_offset + done, x_file);
        }
      if (y_file != NULL && !(x_file == y_file && x_offset == y_offset))
        {
          spv_t t = (buf_y == r) ? r : (buf_y == buf[0]) ? buf[0] : 
                    (buf_y == buf[1]) ? buf[1] : NULL;
          spv_seek_and_read (t, do_now, y_offset + done, y_file);
        }

      /* Do operation */
      switch (operation)
        {
          case SPV_ELEMENTWISE_SET: {
                spv_set (buf_r, buf_x, do_now);
                break;
              }
          case SPV_ELEMENTWISE_ADD: {
                spv_add (buf_r, buf_x, buf_y, do_now, p);
                break;
              }
          case SPV_ELEMENTWISE_ADDSP: {
                spv_add_sp (buf_r, buf_x, y[y_offset], do_now, p);
                break;
              }
          case SPV_ELEMENTWISE_SUB: {
                spv_sub (buf_r, buf_x, buf_y, do_now, p);
                break;
              }
          case SPV_ELEMENTWISE_SUBSP: {
                spv_sub_sp (buf_r, buf_x, y[y_offset], do_now, p);
                break;
              }
          case SPV_ELEMENTWISE_NEG: {
                spv_neg (buf_r, buf_x, do_now, p);
                break;
              }
          case SPV_ELEMENTWISE_MUL: {
                spv_pwmul (buf_r, buf_x, buf_y, do_now, p, d);
                break;
              }
          case SPV_ELEMENTWISE_SETSP: {
                ASSERT_ALWAYS (x != NULL);
                spv_set_sp (buf_r, x[x_offset], do_now);
                break;
              }
          case SPV_ELEMENTWISE_RANDOM: {
                spv_random (buf_r, do_now, p);
                break;
              }
        }
      
      if (buf_r != buf[0])
        buf_r += do_now;
      else
        spv_seek_and_write (buf_r, do_now, r_offset + done, r_file);

      if (buf_x != buf[0])
        buf_x += do_now;
      if (buf_y != buf[0] && buf_y != buf[1])
        buf_y += do_now;

      done += do_now;
    }
  
  if (buf[0] != NULL)
    sp_aligned_free (buf[0]);
  if (buf[1] != NULL)
    sp_aligned_free (buf[1]);

  if (r_file != NULL)
    fflush (r_file);
}


/* r = x */
void
spv_set (spv_t r, spv_tc x, spv_size_t len)
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
  
  if (x == r)
    {
      for (i = 0; i < len - 1 - i; i++)
        {
          sp_t t = r[i];
          r[i] = r[len - 1 - i];
          r[len - 1 - i] = t;
        }
    }
  else
    {
      ASSERT (r >= x + len || x >= r + len);

      for (i = 0; i < len; i++)
        r[i] = x[len - 1 - i];
    }
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
spv_cmp (spv_tc x, spv_tc y, spv_size_t len)
{
  spv_size_t i;

  for (i = 0; i < len; i++)
    if (x[i] != y[i])
      return 1;

  return 0;
}

/* r = x + y */
void
spv_add (spv_t r, spv_tc x, spv_tc y, spv_size_t len, sp_t m)
{
  spv_size_t i;
  
  ASSERT (r >= x + len || x >= r);
  ASSERT (r >= y + len || y >= r);
  
  for (i = 0; i < len; i++)
    r[i] = sp_add (x[i], y[i], m);
}

/* r = [x[0] + y, x[1] + y, ... ] */
void
spv_add_sp (spv_t r, spv_tc x, sp_t c, spv_size_t len, sp_t m)
{
  spv_size_t i;

  for (i = 0; i < len; i++)
    r[i] = sp_add (x[i], c, m);
}

/* r = x - y */
void
spv_sub (spv_t r, spv_tc x, spv_tc y, spv_size_t len, sp_t m)
{
  spv_size_t i;
  
  ASSERT (r >= x + len || x >= r);
  ASSERT (r >= y + len || y >= r);
   
  for (i = 0; i < len; i++)
    r[i] = sp_sub (x[i], y[i], m);
}

/* r = [x[0] - y, x[1] - y, ... ] */
void
spv_sub_sp (spv_t r, spv_tc x, sp_t c, spv_size_t len, sp_t m)
{
  spv_size_t i;

  for (i = 0; i < len; i++)
    r[i] = sp_sub (x[i], c, m);
}

/* r = [-x[0], -x[1], ... ] */
void
spv_neg (spv_t r, spv_tc x, spv_size_t len, sp_t m)
{
  spv_size_t i;

  for (i = 0; i < len; i++)
    r[i] = sp_sub (0, x[i], m);
}

/* Pointwise multiplication
 * r = [x[0] * y[0], x[1] * y[1], ... ] */
void
spv_pwmul (spv_t r, spv_tc x, spv_tc y, spv_size_t len, sp_t m, sp_t d)
{
  spv_size_t i = 0;
  
  ASSERT (r >= x + len || x >= r);
  ASSERT (r >= y + len || y >= r);

#if defined(HAVE_SSE2) 
  
  #if SP_NUMB_BITS < 32

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
      pstore(t0, r + i);
    }

  #elif GMP_LIMB_BITS == 32   /* 64-bit sp_t on a 32-bit machine */

  __m128i t0, t1, t2, t3, t4, t5, vm, vd;

  vm = pshufd(pcvt_i64(m), 0x44);
  vd = pshufd(pcvt_i64(d), 0x44);

  for (; i < (len & (spv_size_t)(~1)); i += 2)
    {
      t4 = pload((__m128i *)(x + i));
      t5 = pload((__m128i *)(y + i));

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
      pstore(t3, r + i);
    }
  #endif

#endif

  for (; i < len; i++)
    r[i] = sp_mul (x[i], y[i], m, d);
}

/* Pointwise multiplication, second input is read in reverse
 * r = [x[0] * y[len - 1], x[1] * y[len - 2], ... x[len - 1] * y[0]] */
void
spv_pwmul_rev (spv_t r, spv_tc x, spv_tc y, spv_size_t len, sp_t m, sp_t d)
{
  spv_size_t i;
  
  ASSERT (r >= x + len || x >= r);
  ASSERT (r >= y + len || y >= r);

  for (i = 0; i < len; i++)
    r[i] = sp_mul (x[i], y[len - 1 - i], m, d);
}

/* dst = src * c */
void
spv_mul_sp (spv_t r, spv_tc x, sp_t c, spv_size_t len, sp_t m, sp_t d)
{
  spv_size_t i = 0;
  
  ASSERT (r >= x + len || x >= r);
  
#if defined(HAVE_SSE2) 
  
  #if SP_NUMB_BITS < 32

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

      pstore(t0, r + i);
    }

  #elif GMP_LIMB_BITS == 32   /* 64-bit sp_t on a 32-bit machine */

  __m128i t0, t1, t2, t3, t4, t5, vm, vd, vc;

  vm = pshufd(pcvt_i64(m), 0x44);
  vd = pshufd(pcvt_i64(d), 0x44);
  vc = pshufd(pcvt_i64(c), 0x44);

  for (; i < (len & (spv_size_t)(~1)); i += 2)
    {
      t4 = pload((__m128i *)(x + i));
      t5 = vc;

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
      pstore(t3, r + i);
    }
  #endif

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


spv_size_t 
spv_seek_and_read (spv_t ptr, const spv_size_t nread, const spv_size_t offset, 
                   FILE *f)
{
  spv_size_t r;
  int64_t foffset;
  
  ASSERT_ALWAYS (INT64_MAX / sizeof(sp_t) >= offset);
  foffset = (int64_t) offset * sizeof(sp_t);

#if defined(_OPENMP)
#pragma omp critical
#endif
  {
    if (aux_ftell64(f) != foffset)
      {
#ifdef TRACE_SEEK
        printf ("Seeking to file position %" PRId64 "\n", foffset);
#endif
        if (aux_fseek64 (f, foffset, SEEK_SET) != 0)
          {
            fprintf (stderr, "%s(): fseek() returned error %d\n", 
                     __func__, errno);
            abort ();
          }
      }
    
    r = fread(ptr, sizeof(sp_t), nread, f);
  }

  if (r != nread)
    {
      int saved_errno = errno;
      perror (NULL);
      fprintf (stderr, "%s(): Error reading data, r = %lu, errno = %d\n",
               __func__, (unsigned long) r, saved_errno);
      abort();
    }

  return r;
}


spv_size_t 
spv_seek_and_write (const spv_t ptr, const spv_size_t nwrite, 
                    const spv_size_t offset, FILE *f)
{
  spv_size_t r;
  int64_t foffset;

  ASSERT_ALWAYS (INT64_MAX / sizeof(sp_t) >= offset);
  foffset = (int64_t) offset * sizeof(sp_t);

#if defined(_OPENMP)
#pragma omp critical
#endif
  {
    if (aux_ftell64(f) != foffset)
      {
#ifdef TRACE_SEEK
        printf ("Seeking to file position %" PRId64 "\n", foffset);
#endif
        if (aux_fseek64 (f, foffset, SEEK_SET) != 0)
          {
            fprintf (stderr, "%s(): fseek() returned error %d\n", 
                     __func__, errno);
            abort ();
          }
      }
    
    r = fwrite(ptr, sizeof(sp_t), nwrite, f);

    if (r == nwrite)
      fflush(f); /* FIXME: Do we want this flush? */
  }
  
  if (r != nwrite)
    {
      int saved_errno = errno;
      perror (NULL);
      fprintf (stderr, "%s(): Error writing data, r = %lu, errno = %d\n",
               __func__, (unsigned long) r, saved_errno);
      abort();
    }

  return r;
}


void
spv_print_vec (const spv_t spv, const sp_t m, const spv_size_t l, 
               const char *prefix, const char *suffix)
{
  spv_size_t i;
  printf ("%s [%" PRISP, prefix, spv[0]);
  for (i = 1; i < l; i++)
    printf (", %" PRISP, spv[i]);
  printf ("] (mod %" PRISP ")%s", m, suffix);
}
