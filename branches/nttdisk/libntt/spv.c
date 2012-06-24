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
