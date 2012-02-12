#ifndef __ASM_REDC_H__
#define __ASM_REDC_H__

#include <gmp.h>

extern void ecm_redc3(mp_limb_t *cp, const mp_limb_t *np, mp_size_t nn, mp_limb_t Nprim);


/* WARNING: the size-1 version doesn't take pointers in input */
extern mp_limb_t mulredc1(mp_limb_t *z, mp_limb_t x, mp_limb_t y, mp_limb_t m, mp_limb_t inv_m);

extern mp_limb_t mulredc2(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y, const mp_limb_t *m, mp_limb_t inv_m);
extern mp_limb_t mulredc3(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y, const mp_limb_t *m, mp_limb_t inv_m);
extern mp_limb_t mulredc4(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y, const mp_limb_t *m, mp_limb_t inv_m);
extern mp_limb_t mulredc5(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y, const mp_limb_t *m, mp_limb_t inv_m);
extern mp_limb_t mulredc6(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y, const mp_limb_t *m, mp_limb_t inv_m);
extern mp_limb_t mulredc7(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y, const mp_limb_t *m, mp_limb_t inv_m);
extern mp_limb_t mulredc8(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y, const mp_limb_t *m, mp_limb_t inv_m);
extern mp_limb_t mulredc9(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y, const mp_limb_t *m, mp_limb_t inv_m);
extern mp_limb_t mulredc10(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y, const mp_limb_t *m, mp_limb_t inv_m);
extern mp_limb_t mulredc11(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y, const mp_limb_t *m, mp_limb_t inv_m);
extern mp_limb_t mulredc12(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y, const mp_limb_t *m, mp_limb_t inv_m);
extern mp_limb_t mulredc13(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y, const mp_limb_t *m, mp_limb_t inv_m);
extern mp_limb_t mulredc14(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y, const mp_limb_t *m, mp_limb_t inv_m);
extern mp_limb_t mulredc15(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y, const mp_limb_t *m, mp_limb_t inv_m);
extern mp_limb_t mulredc16(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y, const mp_limb_t *m, mp_limb_t inv_m);
extern mp_limb_t mulredc17(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y, const mp_limb_t *m, mp_limb_t inv_m);
extern mp_limb_t mulredc18(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y, const mp_limb_t *m, mp_limb_t inv_m);
extern mp_limb_t mulredc19(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y, const mp_limb_t *m, mp_limb_t inv_m);
extern mp_limb_t mulredc20(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y, const mp_limb_t *m, mp_limb_t inv_m);

/* 0:mulredc 1:mul+redc_1 2:mul+redc_2 3:mul+redc3 */
#define TUNE_MULREDC_TABLE {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
/* 0:mulredc 1:sqr+redc_1 2:sqr+redc_2 3:sqr+redc3 */
#define TUNE_SQRREDC_TABLE {0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,1,1,1,1,1}

#endif
