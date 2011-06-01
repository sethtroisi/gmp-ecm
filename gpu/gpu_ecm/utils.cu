#include "main.h"
#include "utils.h"

unsigned int findfactor(mpz_t N, mpz_t xfin, mpz_t zfin)
{
	//int probprime;
	int findfactor=0;
	mpz_t gcd;
	mpz_t factor;
	mpz_t temp;

  mpz_init (gcd);
  mpz_init (temp);
  mpz_init (factor);

	mpz_set_ui(factor,0);

	gmp_printf("  xfin=%Zd\n  zfin=%Zd\n",xfin,zfin);
	// tester si pgcd =N et =0
	mpz_gcd(gcd,zfin,N);
	
	if (mpz_cmp_ui(gcd,1)==0)
	{
		mpz_invert(zfin,zfin,N);
		mpz_mul(xfin,xfin,zfin);
		mpz_mod(xfin,xfin,N);
		gmp_printf("  xunif=%Zd\n",xfin);
		mpz_gcd(gcd,xfin,N);
			
		if (mpz_cmp_ui(gcd,1)==0)
		{
			printf("  #No factors found. You shoud try with a bigger B1.\n");
		}
		else if (mpz_cmp(gcd,N)==0)
		{
			printf("  #No factors found. You should try with a smaller B1\n");
		}
		else
		{
			mpz_set(factor,gcd);
			gmp_printf("  #Factor found : %Zd (with x/z)\n",factor);
			findfactor=1;
		}
	}
	else if (mpz_cmp(gcd,N)==0)
	{
		printf("  #No factors found. You should try with a smaller B1\n");
	}
	else //gcd !=1 gcd!=N (and gcd>0 because N>0) so we found a factor
	{
		mpz_set(factor,gcd);
		gmp_printf("  #Factor found : %Zd (with z)\n",factor);
		findfactor=1;
	}
		/*
	if (mpz_cmp_ui(factor,0)!=0)
	{
		mpz_divexact(mpztemp,N,factor);
		//gmp_sprintf (str[i],"%scofactor:=%Zd ",str[i],mpztemp);
		probprime=mpz_probab_prime_p(mpztemp,5);
		if (probprime==2)
			gmp_sprintf (str[i],"%s definitely prime\n",str[i]);
		else if (probprime==1)
			gmp_sprintf (str[i],"%s probably prime\n",str[i]);
		else if (probprime==0)
			gmp_sprintf (str[i],"%s definitely composite\n",str[i]);
		else	
			gmp_sprintf (str[i],"%s \n",str[i]);
	}
 */
 	mpz_clear(temp);
  mpz_clear(gcd);
  mpz_clear(factor);
	
	return findfactor;
}

void biguint_print (biguint_t a)
{
  unsigned int i;

  printf ("%u", a[0]);
  for (i = 1; i < SIZE_NUMBER; i++)
    if (a[i]!=0)
			printf ("+%u*2^%u", a[i], 32*i);
  //printf ("\n");
}

void bigint_print (dbigint_t a)
{
  unsigned int i;

  printf ("%d", a[0]);
  for (i = 1; i < SIZE_NUMBER; i++)
    printf ("+%d*2^%u", a[i], 32*i);
  printf ("\n");
}

void mpz_to_biguint (biguint_t a, mpz_t b)
{
	int i;

	for (i=0;i<SIZE_NUMBER;i++)
	{
		if (i%2 == 0)
			a[i]=(mpz_getlimbn(b,i/2) & 0x00000000ffffffff);
		else
			a[i]=(mpz_getlimbn(b, i/2) >> 32);	
	}
}

void biguint_to_mpz (mpz_t a, biguint_t b)
{
	int i;
	unsigned long temp;
	
	mpz_set_ui(a,0);

	for (i=SIZE_NUMBER-1;i>=0;i--)
	{
		if (i%2 == 0)
			mpz_add_ui(a,a,b[i]);
		else
		{
			temp=(unsigned long)b[i];
			mpz_add_ui(a,a,(temp<<32));
		}
		if (i!=0 && i%2==0)
			mpz_mul_2exp(a,a,64);
	}
}

unsigned long getprime (unsigned long pp)
{
 static unsigned long offset = 0;     // offset for current primes 
 static long current = -1;            // index of previous prime 
 static unsigned *primes = NULL;      // small primes up to sqrt(p) 
 static unsigned long nprimes = 0;    // length of primes[] 
 static unsigned char *sieve = NULL;  // sieving table 
 static long len = 0;                 // length of sieving table 
 static unsigned long *moduli = NULL; // offset for small primes 

 if (pp == 0) // free the tables, and reinitialize 
   {
     offset = 0.0;
     current = -1;
     free (primes);
     primes = NULL;
     nprimes = 0;
     free (sieve);
     sieve = NULL;
     len = 0;
     free (moduli);
     moduli = NULL;
     return pp;
   }

 // the following complex block is equivalent to:
 // while ((++current < len) && (sieve[current] == 0));
 // but is faster.
 
 {
   unsigned char *ptr = sieve + current;
   unsigned char *end = sieve + len;
   while ((++ptr < end) && (*ptr == 0));
   current = ptr - sieve;
 }

 if (current < len) // most calls will end here 
   return offset + 2 * current;

 // otherwise we have to sieve 
 offset += 2 * len;

 // first enlarge sieving table if too small 
 if ((unsigned long) len * len < offset)
   {
     free (sieve);
     len *= 2;
     sieve = (unsigned char *) malloc (len * sizeof (unsigned char));
     // assume this "small" malloc will not fail in normal usage 
     assert(sieve != NULL);
   }

 // now enlarge small prime table if too small 
 if ((nprimes == 0) ||
     (primes[nprimes - 1] * primes[nprimes - 1] < offset + len))
     {
       if (nprimes == 0) // initialization 
         {
           nprimes = 1;
           primes = (unsigned *) malloc (nprimes * sizeof(unsigned long));
           // assume this "small" malloc will not fail in normal usage 
           assert(primes != NULL);
           moduli = (long unsigned int *) malloc (nprimes *
                                                  sizeof(unsigned long));
           // assume this "small" malloc will not fail in normal usage 
           assert(moduli != NULL);
           len = 1;
           sieve = (unsigned char *) malloc(len * sizeof(unsigned char));//len=1 here
           // assume this "small" malloc will not fail in normal usage 
           assert(sieve != NULL);
           offset = 5.0;
           sieve[0] = 1; // corresponding to 5 
           primes[0] = 3;
           moduli[0] = 1; // next odd multiple of 3 is 7, i.e. next to 5 
           current = -1;
           return 3.0;
         }
       else
         {
           unsigned int i, p, j, ok;

           i = nprimes;
           nprimes *= 2;
           primes = (unsigned *) realloc (primes, nprimes *
                                          sizeof(unsigned long));
           moduli = (unsigned long*) realloc (moduli, nprimes *
                                                   sizeof(unsigned long));
           // assume those "small" realloc's will not fail in normal usage 
           assert(primes != NULL && moduli != NULL);
           for (p = primes[i-1]; i < nprimes; i++)
             {
               // find next (odd) prime > p 
               do
                 {
                   for (p += 2, ok = 1, j = 0; (ok != 0) && (j < i); j++)
                     ok = p % primes[j];
                 }
               while (ok == 0);
               primes[i] = p;
               // moduli[i] is the smallest m such that offset + 2*m = k*p
               j = offset % p;
               j = (j == 0) ? j : p - j; // -offset mod p 
               if ((j % 2) != 0)
                 j += p; // ensure j is even 
               moduli[i] = j / 2;
             }
         }
     }

 // now sieve for new primes 
 {
   long i;
   unsigned long j, p;

   for (i = 0; i < len; i++)
     sieve[i] = 1;
   for (j = 0; j < nprimes; j++)
     {
       p = primes[j];
       for (i = moduli[j]; i < len; i += p)
         sieve[i] = 0;
       moduli[j] = i - len; // for next sieving array 
     }
 }

 current = -1;
 while ((++current < len) && (sieve[current] == 0));
 assert(current < len);//otherwise we found a prime gap >= sqrt(x) around x
 return offset + 2 * current;
}

