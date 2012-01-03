#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <assert.h>

//also define number of threads; problems ig bigger tjhan the size of a warp (ie 32)
//#define NB_DIGITS 32

#define CHKSUMMOD 4294967291U
#define VERSION "0.1"
#define ECM_FACTOR_FOUND 2
#define ECM_NO_FACTOR_FOUND 0

#define TWO32 4294967296 // 2^32 
#define TWO31 2147483648 // 2^31 

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
  //clock_t init;
  //clock_t computation;
  long init;
  long computation;
};


#ifdef CC20
  #define MAJOR 2
  #define MINOR 0
  #define CURVES_BY_BLOCK 32//16
  #define BLOCKS_BY_MP 1//3
  #define CURVES_BY_MP CURVES_BY_BLOCK*BLOCKS_BY_MP
#endif

#ifdef CC13
  #define MAJOR 1
  #define MINOR 3
  #define CURVES_BY_BLOCK 16
  #define BLOCKS_BY_MP 1
  #define CURVES_BY_MP CURVES_BY_BLOCK*BLOCKS_BY_MP
#endif
