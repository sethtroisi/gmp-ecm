/* cgbn_stage1.h: header for CGBN (GPU) based ecm stage 1.

  Copyright 2021 Seth Troisi

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

#ifndef _CGBN_STAGE1_H
#define _CGBN_STAGE1_H 1

#include <stdint.h>

#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

// define the instance structure
typedef struct {
    int n_log2;

    // Number of curves to run
    uint32_t curves;

    uint32_t num_bits;
    // Bits (malloc'ed in generate_instance)
    char    *s_bits;

    // Sigma of first curve
    uint32_t sigma;

} ecm_params_t;


int run_cgbn(mpz_t *factors, int *array_stage_found,
             const mpz_t N, const mpz_t s, float *gputime,
             ecm_params_t *ecm_params);

#ifdef __cplusplus
}
#endif


#endif /* _CGBN_STAGE1_H */
