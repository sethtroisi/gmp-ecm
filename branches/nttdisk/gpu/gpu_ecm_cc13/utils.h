#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <assert.h>

//ca définit aussi le nombre de threads; problemes si plus grand que la taille d'un warp (ie 32)
#define SIZE_NUMBER 32 

#define TWO32 4294967296 // 2^32 

//typedef unsigned int* biguint_t;
//typedef int* bigint_t;
typedef unsigned int biguint_t[SIZE_NUMBER];
typedef unsigned int dbiguint_t[2*SIZE_NUMBER];
typedef int dbigint_t[2*SIZE_NUMBER];

typedef struct clock2_t clock2_t;
struct clock2_t
{
	clock_t init;
	clock_t computation;
};

//#define MAX_CURVES 500

#ifdef CC20
#define MAJOR 2
#define MINOR 0
#define MAX_USE_SHARED_MEM 2560
//#define MAX_CURVES_BY_MP 20
//#define NUMBER_OF_MP 14
//#define MAX_NUMBER_OF_CURVES MAX_CURVES_BY_MP*NUMBER_OF_MP
#endif

#ifdef CC13
#define MAJOR 1
#define MINOR 3
#define MAX_USE_SHARED_MEM 2632
//#define MAX_CURVES_BY_MP 6
//#define NUMBER_OF_MP 24
//#define MAX_NUMBER_OF_CURVES MAX_CURVES_BY_MP*NUMBER_OF_MP
#endif



//functions in main.cu
unsigned int findfactor(mpz_t N, mpz_t xfin, mpz_t zfin);
void biguint_print (biguint_t a);
void bigint_print (dbigint_t a);
void mpz_to_biguint (biguint_t a, mpz_t b);
void biguint_to_mpz (mpz_t a, biguint_t b);
void calculParam (mpz_t sigma, mpz_t N, mpz_t d, mpz_t x0, mpz_t z0, mpz_t u, mpz_t v);

//functions in cudautils.cu
clock2_t cuda_Main(biguint_t h_N, biguint_t h_invmod, biguint_t *h_xarray, biguint_t *h_zarray, biguint_t *h_darray, unsigned int B1,unsigned int number_of_curves);
unsigned long getprime (unsigned long pp);

#define NV 10
#define ADD 6.0 // number of multiplications in an addition 
#define DUP 5.0 // number of multiplications in a duplicate 
double lucas_cost (unsigned long n, double v);
