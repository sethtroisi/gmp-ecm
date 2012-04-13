#ifdef WITH_GPU

#define VERSION_GPUECM "0.1"

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

/* Uncomment the next line to print the number of remaining iterations. */
//#define PRINT_REMAINING_ITER
#endif
