#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#ifndef _MSC_VER
#include <unistd.h>
#include <sys/resource.h>
#else
#include <windows.h>
#endif

#include <gmp.h>

#ifdef _MSC_VER
#define __asm__ asm
#endif

#define VERSION_GPUECM "0.1"

#define ECM_FACTOR_FOUND 2
#define ECM_NO_FACTOR_FOUND 0


#ifndef NB_DIGITS
  #define NB_DIGITS 32 //by default
#endif

#ifndef DIGITS
  #define DIGITS 0
#endif

#if (DIGITS==0)
  #define SIZE_DIGIT 32
  typedef unsigned int digit_t;
  typedef int carry_t;
#endif

#define VOL volatile 
//#define VOL

#define MAX_BITS SIZE_DIGIT*NB_DIGITS
typedef digit_t VOL biguint_t[NB_DIGITS];
typedef digit_t VOL dbiguint_t[NB_DIGITS+1];
typedef carry_t VOL dbigint_t[NB_DIGITS+1];

typedef struct clock2_t clock2_t;
struct clock2_t
{
  long init;
  long computation;
};

