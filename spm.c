#include <stdlib.h>
#include "sp.h"

spm_t
spm_init (sp_t sp)
{
//  int i;
  sp_t a;
  spm_t spm = (spm_t) malloc (sizeof (__spm_struct));

  spm->sp = sp;
  invert_limb (spm->mul_c, sp);

#if 0
  /* init inverses */
  
  _spm->_sp_inv_c = (spr *) malloc (SPM_INVERSE_MAX * sizeof (spr));
  
  for (i = 2; i <= SPM_INVERSE_MAX; i++)
    _spm->_sp_inv_c[i - 2] = spr_inv (i, _sp, _spm->_sp_mul_c);
  
  /* division by 4 */
  for (i = 0; i < 4; i++)
    _spm->_sp_div_4_c[i] = spr_mul (i, _spm->_sp_inv_c[4 - 2], _sp, _spm->_sp_mul_c);
#endif
  
  /* find generator */
  for (a = 2; sp_pow (a, (sp - 1) / 2, sp, spm->mul_c) == 1; a++);
  
  spm->prim_root = a;

  return spm;
}

void
spm_clear (spm_t spm)
{
//  free (spm->_sp_inv_c);
  free (spm);
}
