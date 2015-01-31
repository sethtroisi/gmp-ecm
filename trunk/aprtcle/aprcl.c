/* Copyright 2012-2015 David Cleaver
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#ifdef _MSC_VER
#  include "gettimeofday.h"
#else
#  include <sys/time.h>
#endif
#include <string.h>
#include <gmp.h>
#include "mpz_aprcl.h"

/* these are vars used to time our program */
double t0, t1, t2, t3;
int t_d;
int t_h;
int t_m;
int t_s;

double ttime(double base) {
  struct timeval tod;
  gettimeofday(&tod, NULL);
  return tod.tv_sec + tod.tv_usec/1e6 - base;
}/* method ttime */

/***********************************************************************************/

/* returns the number of decimal digits of n */
unsigned int
b10_digits (const mpz_t n)
{
  mpz_t x;
  unsigned int size;

  size = mpz_sizeinbase (n, 10);

  /* the GMP documentation says mpz_sizeinbase returns the exact value,
     or one too big, thus:
     (a) either n < 10^(size-1), and n has size-1 digits
     (b) or n >= size-1, and n has size digits
     Note: mpz_sizeinbase returns 1 for n=0, thus we always have size >= 1.
  */
				    
  mpz_init (x);
  mpz_ui_pow_ui (x, 10, size - 1);
  if (mpz_cmpabs (n, x) < 0)
    size --;
  mpz_clear (x);

  return size;
}

/***********************************************************************************/

int main(int argc, char* argv[])
{
  int ret_b = 0;
  int ret_a = 0;
  int dig = 0;
  mpz_t test;
#define MAXDIGITS 10000
  char input[MAXDIGITS];
  FILE *fp;

  if (argc < 2 || argc > 3)
  {
    printf("Usage:\n");
    printf(" %s <num>\n", argv[0]);
    printf("    Where <num> is a number to test for primality\n");
    printf(" %s -inp <file>\n", argv[0]);
    printf("    Where <file> contains one or more numbers to test for primality\n");
    printf(" Note: This program returns the APR-CL and MPZ PRP status of the input number(s).\n");
    printf(" Note: This program will also print out timing information for each check.\n");
    printf(" Note: All input numbers are assumed to be base-10 numbers.\n");
    return 0;
  }

  if (strcmp(argv[1], "-inp") == 0)
  {
    fp = fopen(argv[2], "r");
    if(!fp)
    {
      printf("Error, invalid file: %s\n", argv[2]);
      return 1; // bail out if file not found
    }

    while(fgets(input,sizeof(input),fp) != NULL)
    {
      /* make sure last character is 0 (null) */
      if(input[MAXDIGITS-1] != 0) 
         input[MAXDIGITS-1] = 0;

      if (strlen(input) == 0) continue;

      if (input[0] == 0xd || input[0] == 0xa) continue;

      if (input[0] < '0' || input[0] > '9')
      {
        printf("\n  *** Notice, not testing the following line:\n");
        printf("%s\n", input);
        continue;
      }

      if (mpz_init_set_str(test, input, 10) < 0)
      {
        printf(" *** ERROR, invalid number: %s\n", input);
        continue;
      }

      printf("===============================================================================\n");
      dig = b10_digits(test);
      gmp_printf("Testing: %Zd (%d digits)\n", test, dig); fflush(stdout);
      printf("Running the MPZ PRP test with 5 iterations...\n"); fflush(stdout);
      t0 = ttime(0);
      ret_b = mpz_probab_prime_p(test, 5);
      t1 = ttime(t0);

      printf("Running the APRCL prime test...\n"); fflush(stdout);
      t0 = ttime(0);
      ret_a = mpz_aprtcle(test, 1);
      t2 = ttime(t0);

      printf("\n MPZ PRP took %.4f seconds\n", t1);
      if (ret_b == APRTCLE_COMPOSITE)
        printf(" MPZ PRP says the number is COMPOSITE\n");
      else if (ret_b == APRTCLE_PRP)
        printf(" MPZ PRP says the number is PRP\n");
      else if (ret_b == APRTCLE_PRIME)
        printf(" MPZ PRP says the number is PRIME\n");

      printf("APRCL took %.4f seconds\n", t2);
      if (ret_a == APRTCLE_COMPOSITE)
        printf("APRCL says the number is COMPOSITE\n");
      else if (ret_a == APRTCLE_PRP)
        printf("APRCL says the number is PRP\n");
      else if (ret_a == APRTCLE_PRIME)
        printf("APRCL says the number is PRIME\n");

      if ((ret_b == APRTCLE_COMPOSITE && ret_a != APRTCLE_COMPOSITE) ||
          (ret_b != APRTCLE_COMPOSITE && ret_a == APRTCLE_COMPOSITE))
      {
        printf(" *** ATTENTION *** ATTENTION *** ATTENTION *** ATTENTION ***\n");
        printf("MPZ PRP and APRCL do not agree on the status of this number!!!\n");
        printf("Please report this to http://www.mersenneforum.org/showthread.php?t=18353\n");
        gmp_printf("N = %Zd\n", test);
      }
    }
    fclose(fp);
  }
  else
  {
    if (mpz_init_set_str(test, argv[1], 10) < 0)
    {
      mpz_clear(test);
      printf("Error, invalid number: %s\n", argv[1]);
      return 0;
    }

    dig = b10_digits(test);
    gmp_printf("Testing: %Zd (%d digits)\n", test, dig); fflush(stdout);
    printf("Running the MPZ PRP test with 5 iterations...\n"); fflush(stdout);
    t0 = ttime(0);
    ret_b = mpz_probab_prime_p(test, 5);
    t1 = ttime(t0);

    printf("Running the APRCL prime test...\n"); fflush(stdout);
    t0 = ttime(0);
    ret_a = mpz_aprtcle(test, 1);
    t2 = ttime(t0);

    printf("\n MPZ PRP took %.4f seconds\n", t1);
    if (ret_b == APRTCLE_COMPOSITE)
      printf(" MPZ PRP says the number is COMPOSITE\n");
    else if (ret_b == APRTCLE_PRP)
      printf(" MPZ PRP says the number is PRP\n");
    else if (ret_b == APRTCLE_PRIME)
      printf(" MPZ PRP says the number is PRIME\n");

    printf("APRCL took %.4f seconds\n", t2);
    if (ret_a == APRTCLE_COMPOSITE)
      printf("APRCL says the number is COMPOSITE\n");
    else if (ret_a == APRTCLE_PRP)
      printf("APRCL says the number is PRP\n");
    else if (ret_a == APRTCLE_PRIME)
      printf("APRCL says the number is PRIME\n");

    if ((ret_b == APRTCLE_COMPOSITE && ret_a != APRTCLE_COMPOSITE) ||
        (ret_b != APRTCLE_COMPOSITE && ret_a == APRTCLE_COMPOSITE))
    {
      printf(" *** ATTENTION *** ATTENTION *** ATTENTION *** ATTENTION ***\n");
      printf("MPZ PRP and APRCL do not agree on the status of this number!!!\n");
      printf("Please report this to http://www.mersenneforum.org/showthread.php?t=18353\n");
      gmp_printf("N = %Zd\n", test);
    }
  }

  mpz_clear(test);

  return 0;
}
