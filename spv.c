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

#if (defined(__GNUC__) || defined(__ICL)) && \
  	defined(__i386__) && defined(HAVE_SSE2)

  asm volatile (
       "movd %6, %%xmm6            \n\t"
       "pshufd $0x44, %%xmm6, %%xmm5  \n\t"
       "pshufd $0, %%xmm6, %%xmm6  \n\t"
       "movd %7, %%xmm7            \n\t"
       "pshufd $0, %%xmm7, %%xmm7  \n\t"

       "0:                         \n\t"
       "movdqa (%1,%4,4), %%xmm0   \n\t"
       "movdqa (%2,%4,4), %%xmm2   \n\t"
       "pshufd $0x31, %%xmm0, %%xmm1\n\t"
       "pshufd $0x31, %%xmm2, %%xmm3\n\t"
       "pmuludq %%xmm2, %%xmm0     \n\t"
       "pmuludq %%xmm3, %%xmm1     \n\t"

       "movdqa %%xmm0, %%xmm2      \n\t"
       "movdqa %%xmm1, %%xmm3      \n\t"
       "psrlq $" STRING((2*SP_NUMB_BITS - W_TYPE_SIZE)) ", %%xmm2  \n\t"
       "pmuludq %%xmm7, %%xmm2     \n\t"
       "psrlq $" STRING((2*SP_NUMB_BITS - W_TYPE_SIZE)) ", %%xmm3  \n\t"
       "pmuludq %%xmm7, %%xmm3     \n\t"

#if SP_NUMB_BITS < W_TYPE_SIZE - 1
       "psrlq $33, %%xmm2          \n\t"
       "pmuludq %%xmm6, %%xmm2     \n\t"
       "psrlq $33, %%xmm3          \n\t"
       "pmuludq %%xmm6, %%xmm3     \n\t"
       "psubq %%xmm2, %%xmm0       \n\t"
       "psubq %%xmm3, %%xmm1       \n\t"
#else
       "pshufd $0xf5, %%xmm2, %%xmm2 \n\t"
       "pmuludq %%xmm6, %%xmm2     \n\t"
       "pshufd $0xf5, %%xmm3, %%xmm3 \n\t"
       "pmuludq %%xmm6, %%xmm3     \n\t"
       "psubq %%xmm2, %%xmm0       \n\t"
       "psubq %%xmm3, %%xmm1       \n\t"

       "psubq %%xmm5, %%xmm0       \n\t"
       "psubq %%xmm5, %%xmm1       \n\t"
       "pshufd $0xf5, %%xmm0, %%xmm2 \n\t"
       "pshufd $0xf5, %%xmm1, %%xmm3 \n\t"
       "pand %%xmm5, %%xmm2        \n\t"
       "pand %%xmm5, %%xmm3        \n\t"
       "paddq %%xmm2, %%xmm0       \n\t"
       "paddq %%xmm3, %%xmm1       \n\t"
#endif
       "pshufd $0x8, %%xmm0, %%xmm0 \n\t"
       "pshufd $0x8, %%xmm1, %%xmm1 \n\t"
       "punpckldq %%xmm1, %%xmm0   \n\t"
       "psubd %%xmm6, %%xmm0       \n\t"

       "pxor %%xmm1, %%xmm1        \n\t"
       "pcmpgtd %%xmm0, %%xmm1     \n\t"
       "pand %%xmm6, %%xmm1        \n\t"
       "paddd %%xmm1, %%xmm0       \n\t"
       "movdqa %%xmm0, (%3,%4,4)   \n\t"

       "addl $4, %4                \n\t"  /* INC */
       "cmpl %5, %4                \n\t"
       "jne 0b                     \n\t"

       :"=r"(i)
       :"r"(x), "r"(y), "r"(r), "0"(i), 
	"g"(len & (spv_size_t)(~3)), "g"(m), "g"(d)
       :"%xmm0", "%xmm1", "%xmm2", "%xmm3",
        "%xmm5", "%xmm6", "%xmm7", "cc", "memory");
#elif defined( _MSC_VER ) && defined( SSE2)
    __asm
    {   push        esi
        push        edi
        mov         edi, x
        mov         esi, y
        mov         edx, r
        xor         ecx, ecx
        mov         eax, len
        and         eax, ~3

        movd        xmm6, m
        pshufd      xmm5, xmm6, 0x44
        pshufd      xmm6, xmm6, 0
        movd        xmm7, d
        pshufd      xmm7, xmm7, 0

    L0: movdqa      xmm0, [edi+ecx*4]
        movdqa      xmm2, [esi+ecx*4]
        pshufd      xmm1, xmm0, 0x31
        pshufd      xmm3, xmm2, 0x31
        pmuludq     xmm0, xmm2
        pmuludq     xmm1, xmm3

        movdqa      xmm2, xmm0
        movdqa      xmm3, xmm1
        psrlq       xmm2, 2*SP_NUMB_BITS - W_TYPE_SIZE
        pmuludq     xmm2, xmm7
        psrlq       xmm3, 2*SP_NUMB_BITS - W_TYPE_SIZE
        pmuludq     xmm3, xmm7

#if SP_NUMB_BITS < W_TYPE_SIZE - 1
        psrlq       xmm2, 33
        pmuludq     xmm2, xmm6
        psrlq       xmm3, 33
        pmuludq     xmm3, xmm6
        psubq       xmm0, xmm2
        psubq       xmm1, xmm3
#else
        pshufd      xmm2, xmm2, 0xf5
        pmuludq     xmm2, xmm6
        pshufd      xmm3, xmm3, 0xf5
        pmuludq     xmm3, xmm6
        psubq       xmm0, xmm2
        psubq       xmm1, xmm3

        psubq       xmm0, xmm5
        psubq       xmm1, xmm5
        pshufd      xmm2, xmm0, 0xf5
        pshufd      xmm3, xmm1, 0xf5
        pand        xmm2, xmm5
        pand        xmm3, xmm5
        paddq       xmm0, xmm2
        paddq       xmm1, xmm3
#endif
        pshufd      xmm0, xmm0, 0x8
        pshufd      xmm1, xmm1, 0x8
        punpckldq   xmm0, xmm1
        psubd       xmm0, xmm6

        pxor        xmm1, xmm1
        pcmpgtd     xmm1, xmm0
        pand        xmm1, xmm6
        paddd       xmm0, xmm1
        movdqa      [edx+ecx*4], xmm0

        add         ecx, 4
        cmp         eax, ecx
        jne         L0
        mov         i, ecx
        pop         edi
        pop         esi
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
  
#if (defined(__GNUC__) || defined(__ICL)) && \
  	defined(__i386__) && defined(HAVE_SSE2)

  asm volatile (
       "movd %2, %%xmm4            \n\t"
       "pshufd $0, %%xmm4, %%xmm4  \n\t"
       "movd %6, %%xmm6            \n\t"
       "pshufd $0x44, %%xmm6, %%xmm5  \n\t"
       "pshufd $0, %%xmm6, %%xmm6  \n\t"
       "movd %7, %%xmm7            \n\t"
       "pshufd $0, %%xmm7, %%xmm7  \n\t"

       "0:                         \n\t"
       "movdqa (%1,%4,4), %%xmm0   \n\t"
       "pshufd $0x31, %%xmm0, %%xmm1\n\t"
       "pshufd $0x31, %%xmm4, %%xmm3\n\t"
       "pmuludq %%xmm4, %%xmm0     \n\t"
       "pmuludq %%xmm3, %%xmm1     \n\t"

       "movdqa %%xmm0, %%xmm2      \n\t"
       "movdqa %%xmm1, %%xmm3      \n\t"
       "psrlq $" STRING((2*SP_NUMB_BITS - W_TYPE_SIZE)) ", %%xmm2  \n\t"
       "pmuludq %%xmm7, %%xmm2     \n\t"
       "psrlq $" STRING((2*SP_NUMB_BITS - W_TYPE_SIZE)) ", %%xmm3  \n\t"
       "pmuludq %%xmm7, %%xmm3     \n\t"

#if SP_NUMB_BITS < W_TYPE_SIZE - 1
       "psrlq $33, %%xmm2          \n\t"
       "pmuludq %%xmm6, %%xmm2     \n\t"
       "psrlq $33, %%xmm3          \n\t"
       "pmuludq %%xmm6, %%xmm3     \n\t"
       "psubq %%xmm2, %%xmm0       \n\t"
       "psubq %%xmm3, %%xmm1       \n\t"
#else
       "pshufd $0xf5, %%xmm2, %%xmm2 \n\t"
       "pmuludq %%xmm6, %%xmm2     \n\t"
       "pshufd $0xf5, %%xmm3, %%xmm3 \n\t"
       "pmuludq %%xmm6, %%xmm3     \n\t"
       "psubq %%xmm2, %%xmm0       \n\t"
       "psubq %%xmm3, %%xmm1       \n\t"

       "psubq %%xmm5, %%xmm0       \n\t"
       "psubq %%xmm5, %%xmm1       \n\t"
       "pshufd $0xf5, %%xmm0, %%xmm2 \n\t"
       "pshufd $0xf5, %%xmm1, %%xmm3 \n\t"
       "pand %%xmm5, %%xmm2        \n\t"
       "pand %%xmm5, %%xmm3        \n\t"
       "paddq %%xmm2, %%xmm0       \n\t"
       "paddq %%xmm3, %%xmm1       \n\t"
#endif
       "pshufd $0x8, %%xmm0, %%xmm0 \n\t"
       "pshufd $0x8, %%xmm1, %%xmm1 \n\t"
       "punpckldq %%xmm1, %%xmm0   \n\t"
       "psubd %%xmm6, %%xmm0       \n\t"

       "pxor %%xmm1, %%xmm1        \n\t"
       "pcmpgtd %%xmm0, %%xmm1     \n\t"
       "pand %%xmm6, %%xmm1        \n\t"
       "paddd %%xmm1, %%xmm0       \n\t"
       "movdqa %%xmm0, (%3,%4,4)   \n\t"

       "addl $4, %4                \n\t"  /* INC */
       "cmpl %5, %4                \n\t"
       "jne 0b                     \n\t"

       :"=r"(i)
       :"r"(x), "g"(c), "r"(r), "0"(i), 
	"g"(len & (spv_size_t)(~3)), "g"(m), "g"(d)
       :"%xmm0", "%xmm1", "%xmm2", "%xmm3",
        "%xmm4", "%xmm5", "%xmm6", "%xmm7", "cc", "memory");
#elif defined( _MSC_VER ) && defined( SSE2) 
    __asm
    {   push        esi
        push        edi
        xor         ecx, ecx
        mov         edi, x
        mov         esi, c
        mov         edx, r
        mov         eax, len
        and         eax, ~3

        movd        xmm4, esi
        pshufd      xmm4, xmm4, 0
        movd        xmm6, m
        pshufd      xmm5, xmm6, 0x44
        pshufd      xmm6, xmm6, 0
        movd        xmm7, d
        pshufd      xmm7, xmm7, 0

    L0: movdqa      xmm0, [edi+ecx*4]
        pshufd      xmm1, xmm0, 0x31
        pshufd      xmm3, xmm4, 0x31
        pmuludq     xmm0, xmm4
        pmuludq     xmm1, xmm3

        movdqa      xmm2, xmm0
        movdqa      xmm3, xmm1
        psrlq       xmm2, 2*SP_NUMB_BITS - W_TYPE_SIZE
        pmuludq     xmm2, xmm7
        psrlq       xmm3, 2*SP_NUMB_BITS - W_TYPE_SIZE
        pmuludq     xmm3, xmm7

#if SP_NUMB_BITS < W_TYPE_SIZE - 1
        psrlq       xmm2, 33
        pmuludq     xmm2, xmm6
        psrlq       xmm3, 33
        pmuludq     xmm3, xmm6
        psubq       xmm0, xmm2
        psubq       xmm1, xmm3
#else
        pshufd      xmm2, xmm2, 0xf5
        pmuludq     xmm2, xmm6
        pshufd      xmm3, xmm3, 0xf5
        pmuludq     xmm3, xmm6
        psubq       xmm0, xmm2
        psubq       xmm1, xmm3

        psubq       xmm0, xmm5
        psubq       xmm1, xmm5
        pshufd      xmm2, xmm0, 0xf5
        pshufd      xmm3, xmm1, 0xf5
        pand        xmm2, xmm5
        pand        xmm3, xmm5
        paddq       xmm0, xmm2
        paddq       xmm1, xmm3
#endif
        pshufd      xmm0, xmm0, 0x8
        pshufd      xmm1, xmm1, 0x8
        punpckldq   xmm0, xmm1
        psubd       xmm0, xmm6

        pxor        xmm1, xmm1
        pcmpgtd     xmm1, xmm0
        pand        xmm1, xmm6
        paddd       xmm0, xmm1
        movdqa      [edx+ecx*4], xmm0

        add         ecx, 4
        cmp         eax, ecx
        jne         L0
        mov         i, ecx
        pop         edi
        pop         esi
    }
#endif

  for (; i < len; i++)
    r[i] = sp_mul (x[i], c, m, d);
}

void
spv_random (spv_t x, spv_size_t len, sp_t m)
{
  spv_size_t i;
  mpn_random (x, len);
  
  for (i = 0; i < len; i++)
    while (x[i] >= m)
      x[i] -= m;
}
