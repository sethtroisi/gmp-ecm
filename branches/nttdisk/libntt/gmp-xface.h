/* Part of file gmp-impl.h from GNU MP.

Copyright 1991, 1993, 1994, 1995, 1996, 1997, 1999, 2000, 2001, 2002 Free
Software Foundation, Inc.

This file contains modified code from the GNU MP Library.

The GNU MP Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The GNU MP Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MP Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
MA 02110-1301, USA. */

#ifndef _GMP_XFACE_H
#define _GMP_XFACE_H 1

#include "basicdefs.h"
#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

static INLINE void 
mpz_set_uint64 (mpz_t m, const uint64_t n)
{
  if (sizeof(uint64_t) <= sizeof(unsigned long))
    {
      mpz_set_ui(m, (unsigned long)n);
    }
  else
    {
       /* Read 1 word of 8 bytes in host endianess */
       mpz_import (m, 1, 1, 8, 0, 0, &n);
    }
}

static INLINE uint64_t
mpz_get_uint64 (mpz_t m)
{
  uint64_t n;

  ASSERT_ALWAYS (mpz_sgn(m) >= 0);
  if (sizeof(uint64_t) <= sizeof(unsigned long))
    {
      ASSERT_ALWAYS (mpz_fits_ulong_p(m));
      n = mpz_get_ui(m);
    }
  else
    {
       size_t count;
       /* Write 1 word of 8 bytes in host endianess */
       ASSERT_ALWAYS (mpz_sizeinbase (m, 2) <= 64);
       mpz_export (&n, &count, 1, 8, 0, 0, m);
       ASSERT_ALWAYS (count == 1);
    }
  return n;
}

static INLINE void 
mpz_set_int64 (mpz_t m, const int64_t n)
{
  if (n < 0)
    {
      mpz_set_uint64(m, -n);
      mpz_neg(m, m);
    }
  else
    mpz_set_uint64(m, n);
}

#ifdef __cplusplus
}
#endif

#endif /* _GMP_XFACE_H */
