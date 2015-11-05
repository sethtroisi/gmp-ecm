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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <sys/types.h> /* needed for size_t */

#ifdef __cplusplus
extern "C" {
#endif

/* name some types whose bit widths are unambiguous */

#if defined(_MSC_VER)
typedef signed __int32 int32_t;
typedef unsigned __int32 uint32_t;
typedef signed __int64 int64_t;
typedef unsigned __int64 uint64_t;

/* portable formatting of wide types */
#define PRId64 "I64d"
#define PRIu64 "I64u"
#define PRIx64 "I64x"

#else
#include <inttypes.h>
#endif

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
#define ATTRIBUTE_ALWAYS_INLINE __attribute__ ((always_inline))
#else
#define ATTRIBUTE_UNUSED
#define ATTRIBUTE_CONST
#define ATTRIBUTE_ALWAYS_INLINE
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

#ifndef alloca
#ifdef __GNUC__
# define alloca __builtin_alloca
#elif defined (__DECC)
# define alloca(x) __ALLOCA(x)
#elif defined (_MSC_VER)
# include <malloc.h>
# define alloca _alloca
#elif HAVE_ALLOCA_H || defined (sun)
# include <alloca.h>
#elif defined (_AIX) || defined (_IBMR2)
#pragma alloca
#else
  char *alloca ();
#endif
#endif

#ifndef INLINE
#if defined(__GNUC__)
#define INLINE inline
#elif defined(_MSC_VER)
#define INLINE _inline
#else
#define INLINE /* nothing */
#endif
#endif

#ifdef __cplusplus
}
#endif


#endif /* _BASICDEFS_H */
