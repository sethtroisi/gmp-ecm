#ifndef __ASM_REDC_H__
#define __ASM_REDC_H__

#include <gmp.h>

/* Signals that we have assembly code for variable size redc */
#define HAVE_NATIVE_REDC3

extern void ecm_redc3(mp_limb_t *, const mp_limb_t *, mp_size_t, mp_limb_t);


/* WARNING: the size-1 version doesn't take pointers in input */
extern mp_limb_t mulredc1(mp_limb_t *, mp_limb_t, mp_limb_t, mp_limb_t, mp_limb_t);

extern mp_limb_t mulredc2(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t);
extern mp_limb_t mulredc3(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t);
extern mp_limb_t mulredc4(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t);
extern mp_limb_t mulredc5(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t);
extern mp_limb_t mulredc6(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t);
extern mp_limb_t mulredc7(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t);
extern mp_limb_t mulredc8(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t);
extern mp_limb_t mulredc9(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t);
extern mp_limb_t mulredc10(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t);
extern mp_limb_t mulredc11(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t);
extern mp_limb_t mulredc12(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t);
extern mp_limb_t mulredc13(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t);
extern mp_limb_t mulredc14(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t);
extern mp_limb_t mulredc15(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t);
extern mp_limb_t mulredc16(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t);
extern mp_limb_t mulredc17(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t);
extern mp_limb_t mulredc18(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t);
extern mp_limb_t mulredc19(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t);
extern mp_limb_t mulredc20(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t);

#endif
