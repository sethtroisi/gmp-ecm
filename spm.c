/* spm.c - "small prime modulus" functions to precompute an inverse and a
   primitive root for a small prime

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

#include <stdlib.h>
#include "sp.h"

spm_t
spm_init (sp_t sp)
{
  sp_t a;
  spm_t spm = (spm_t) malloc (sizeof (__spm_struct));

  spm->sp = sp;
  invert_limb (spm->mul_c, sp);

  /* find generator */
  for (a = 2; sp_pow (a, (sp - 1) / 2, sp, spm->mul_c) == 1; a++);
  
  spm->prim_root = a;
  spm->inv_prim_root = sp_inv (a, sp, spm->mul_c);

  return spm;
}

void
spm_clear (spm_t spm)
{
  free (spm);
}
