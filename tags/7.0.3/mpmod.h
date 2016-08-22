/* Header for modular multiplication.

Copyright 2012, 2015 Paul Zimmermann.

This file is part of the ECM Library.

The ECM Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The ECM Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the ECM Library; see the file COPYING.LIB.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#define MPMOD_MULREDC 0    /* assembly combined mulredc */
#define MPMOD_MUL_REDC1 1  /* mpn_mul_n or mpn_sqr followed by mpn_redc_1 */
#define MPMOD_MUL_REDC2 2  /* mpn_mul_n or mpn_sqr followed by mpn_redc_2 */
#define MPMOD_MUL_REDCN 3  /* mpn_mul_n or mpn_sqr followed by mpn_redc_n */
#define MPMOD_MUL_REDC_C 4 /* mpn_mul_n or mpn_sqr followed by plain C redc */

/* Up from GMP 5.1.0, mpn_redc{1,2} do not subtract the modulus if needed,
   but return the carry of the final addition */
#ifdef HAVE___GMPN_REDC_1
#ifdef MPN_REDC12_RETURNS_CARRY
#define REDC1(rp,cp,np,nn,invm)                  \
  do {if (__gmpn_redc_1 (rp,cp,np,nn,invm))      \
    mpn_sub_n (rp, rp, np, nn);                  \
  } while(0)
#else
#define REDC1(rp,cp,np,nn,invm) __gmpn_redc_1(rp,cp,np,nn,invm)
#endif
#endif

#ifdef HAVE___GMPN_REDC_2
#ifdef MPN_REDC12_RETURNS_CARRY
#define REDC2(rp,cp,np,nn,invm)                  \
  do {if (__gmpn_redc_2 (rp,cp,np,nn,invm))      \
    mpn_sub_n (rp, rp, np, nn);                  \
  } while (0)
#else
#define REDC2(rp,cp,np,nn,invm) __gmpn_redc_2(rp,cp,np,nn,invm)
#endif
#endif

