#ifndef _ARI_KS_H
#define _ARI_KS_H

#include "../ecm-impl.h"
#include "auxi.h"


struct ksCstPourMul_s {
  mpres_t invZ;
  mpres_t invT;
  mpres_t x0p;
  mpres_t t0p;
  mpres_t Y0;
  mpres_t Z0;
};
typedef struct ksCstPourMul_s ksCstPourMul[1];

void ksCstPourMul_init (ksCstPourMul cMul,mpmod_t n);
void ksCstPourMul_clear (ksCstPourMul cMul,mpmod_t n);


struct ksSmallConstPourMul_s {
  long invX;
  long invY;
  long invZ;
  long invT;
  long x0p;
  long y0p;
  long z0p;
  long t0p;
  long X0;
  long Y0;
  long Z0;
  long T0;
};
typedef struct ksSmallConstPourMul_s ksSmallConstPourMul[1];


struct kspoint_s {
  mpres_t X;
  mpres_t Y;
  mpres_t Z;
  mpres_t T;
};
typedef struct kspoint_s ksPoint[1];

void ksPoint_init(ksPoint p,mpmod_t n);
void ksPoint_clear(ksPoint p,mpmod_t n);




#include "generation.h"


void mulKS(ksPoint P,ksCstPourMul cMul,mpres_t be,mpres_t ga,mpz_t k,mpmod_t n);

void doubleKS2(ksPoint P,ksCstPourMul cMul,mpres_t u,mpres_t v,mpmod_t n);

void doubleKS(ksPoint P2,const ksPoint P,ksCstPourMul cMul,mpres_t u,mpres_t v,mpmod_t n);

void loopMulKS(ksPoint Pm,ksPoint Pp,ksCstPourMul cMul,mpres_t u,mpres_t v,mpmod_t n);


void doubleKSsmallParam(ksPoint P2,const ksPoint P,ksSmallConstPourMul cMul,mpres_t u,mpres_t v,mpmod_t n);

void doubleKSsmallParam2(ksPoint P,ksSmallConstPourMul cMul,ksPoint Pi,mpres_t u,mpres_t v,mpmod_t n);

void mulKSsmallParam(ksPoint P,ksSmallConstPourMul cMul,mpres_t be,mpres_t ga,mpz_t k,mpmod_t n);

void loopMulKSsmallParam (ksPoint P,ksPoint Q,ksSmallConstPourMul cMul,ksPoint R,mpres_t u,mpres_t v,mpmod_t n);

void hadamard (ksPoint P2,const ksPoint P,mpres_t u,mpres_t v, mpmod_t n);




#endif
