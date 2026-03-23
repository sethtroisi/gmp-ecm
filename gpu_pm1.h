/* gpu_pm1.h: header for GPU (CGBN) P-1 stage 1.

  Copyright 2021-2026 Seth Troisi

  This program is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 3 of the License, or (at your
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

#ifndef _GPU_PM1_H
#define _GPU_PM1_H 1

#include <stdint.h>
#include <gmp.h>

#include "ecm.h"

#ifdef __cplusplus
extern "C" {
#endif

int cgbn_pm1_stage1(
             const mpz_t *numbers, const mpz_t *x0, const mpz_t *orig_x0,
             mpz_t *factors, mpz_t *residuals,
             const uint64_t B1, uint32_t curves,
             float *gputime, const ecm_params params, ecm_params mutable_params);

#ifdef __cplusplus
}
#endif


#endif /* _GPU_PM1_H */
