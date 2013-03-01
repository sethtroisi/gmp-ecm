/* basicdefs.h: header with core definitions that are globally usable
 
  Copyright 2001, 2002, 2003, 2004, 2005 Paul Zimmermann and Alexander Kruppa.
 
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
  Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
  02111-1307, USA.
*/

#ifndef _BASICDEFS_H
#define _BASICDEFS_H 1

#include "config.h"

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h> /* needed for size_t */
#endif

/* name some types whose bit widths are unambiguous */

#if defined(HAVE_INTTYPES_H)
#include <inttypes.h>
#elif defined(_MSC_VER)
typedef signed __int32 int32_t;
typedef unsigned __int32 uint32_t;
typedef signed __int64 int64_t;
typedef unsigned __int64 uint64_t;

/* portable formatting of wide types */
#define PRId64 "I64d"
#define PRIu64 "I64u"
#define PRIx64 "I64x"

#endif

#include <stdio.h> /* needed for "FILE *" */
#include <limits.h>

/* Warnings about unused parameters by gcc can be suppressed by prefixing 
   parameter with ATTRIBUTE_UNUSED when parameter can't be removed, i.e. 
   for interface consistency reasons */
#ifdef __GNUC__
#if    __GNUC__ >= 3
#define ATTRIBUTE_UNUSED __attribute__ ((unused))
#else
#define ATTRIBUTE_UNUSED
#endif
#define ATTRIBUTE_CONST __attribute__ ((const))
#else
#define ATTRIBUTE_UNUSED
#define ATTRIBUTE_CONST
#endif

#ifndef LIKELY
#if defined(__GNUC__)
#define LIKELY(x) __builtin_expect ((x) != 0, 1)
#else
#define LIKELY(x) x
#endif
#endif

#ifndef UNLIKELY
#if defined(__GNUC__)
#define UNLIKELY(x) __builtin_expect ((x) != 0, 0)
#else
#define UNLIKELY(x) x
#endif
#endif

#include <assert.h>
#define ASSERT_ALWAYS(expr)   assert (expr)
#ifdef WANT_ASSERT
#define ASSERT(expr)   assert (expr)
#else
#define ASSERT(expr)   do {} while (0)
#endif

#define MAX(x,y) (((x)<(y))?(y):(x))
#define MIN(x,y) (((x)<(y))?(x):(y))

/* expanding macros and then turning them 
   into strings requires two levels of macro-izing */

#define _(x) #x
#define STRING(x) _(x)

/* include (hopefully standard) SSE2 interface, along
   with mnemonics for the small part of the ISA that 
   we need */

#ifdef HAVE_SSE2
#include <emmintrin.h>

#define pload(addr)  _mm_load_si128((__m128i const *)(addr))
#define ploadu(addr)  _mm_loadu_si128((__m128i const *)(addr))
#define pload_lo32(addr)  (__m128i)_mm_load_ss((float *)(addr))
#define pload_lo64(addr)  (__m128i)_mm_load_sd((double const *)(addr))
#define pload_hi64(x, addr)  (__m128i)_mm_loadh_pd((__m128d)x, (double const *)(addr))
#define pstore(x, addr) _mm_store_si128((__m128i *)(addr), x)
#define pstoreu(x, addr) _mm_storeu_si128((__m128i *)(addr), x)
#define pstore_lo32(x, addr)  _mm_store_ss((float *)(addr), (__m128)x)
#define pstore_lo64(x, addr)  _mm_store_sd((double *)(addr), (__m128d)x)
#define pstore_hi64(x, addr)  _mm_storeh_pd((double *)(addr), (__m128d)x)
#define pand _mm_and_si128
#define pxor _mm_xor_si128
#define psetzero() _mm_setzero_si128()
#define paddd _mm_add_epi32
#define paddq _mm_add_epi64
#define psubd _mm_sub_epi32
#define psubq _mm_sub_epi64
#define pmuludq _mm_mul_epu32
#define pslld _mm_slli_epi32
#define psllq _mm_slli_epi64
#define psrld _mm_srli_epi32
#define psrlq _mm_srli_epi64
#define pshufd _mm_shuffle_epi32
#define pcmpgtd _mm_cmpgt_epi32
#define punpcklo32 _mm_unpacklo_epi32
#define punpcklo64 _mm_unpacklo_epi64
#define pcvt_i32 _mm_cvtsi32_si128

#if defined(_WIN64) || defined(__x86_64__)
#define pcvt_i64 _mm_cvtsi64_si128
#define pstore_i64(out, in) out = _mm_cvtsi128_si64(in)
#else
#define pcvt_i64(x) _mm_loadl_epi64((__m128i const *)&(x))
#define pstore_i64(out, in) _mm_storel_epi64((__m128i *)&(out), in)
#endif

#endif

#ifndef alloca
#ifdef __GNUC__
# define alloca __builtin_alloca
#elif defined (__DECC)
# define alloca(x) __ALLOCA(x)
#elif defined (_MSC_VER)
# include <malloc.h>
# define alloca _alloca
#elif defined(HAVE_ALLOCA_H) || defined (sun)
# include <alloca.h>
#elif defined (_AIX) || defined (_IBMR2)
#pragma alloca
#else
  char *alloca ();
#endif
#endif


#endif /* _BASICDEFS_H */
