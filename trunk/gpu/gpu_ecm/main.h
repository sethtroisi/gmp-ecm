#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <assert.h>

//ca définit aussi le nombre de threads; problemes si plus grand que la taille d'un warp (ie 32)
//#define SIZE_NUMBER 32

#define CHKSUMMOD 4294967291U
#define VERSION "0.1"
#define ECM_FACTOR_FOUND 2
#define ECM_NO_FACTOR_FOUND 0

#define TWO32 4294967296 // 2^32 
#define TWO31 2147483648 // 2^31 

typedef unsigned int biguint_t[SIZE_NUMBER];
typedef unsigned int dbiguint_t[SIZE_NUMBER+1];
typedef int dbigint_t[SIZE_NUMBER+1];

typedef struct clock2_t clock2_t;
struct clock2_t
{
	clock_t init;
	clock_t computation;
};


#ifdef CC20
	#define MAJOR 2
	#define MINOR 0
	#define CURVES_BY_BLOCK 16
#endif

#ifdef CC13
	#define MAJOR 1
	#define MINOR 3
	#define CURVES_BY_BLOCK 16
#endif

