#include "utils.h"

int main (int argc, char * argv[]) 
{
	unsigned int i;
#ifdef TEST
	unsigned int j;
#endif
	unsigned int B1;
	unsigned int number_of_curves;
	unsigned int nbfactor=0;

	clock2_t time_gpu;
	int device;

	char usage[]="./gpu_ecm N B1 [-s firstsigma] [-n number_of_curves] [-d device]";

	biguint_t *h_xarray;
	biguint_t *h_zarray;
	biguint_t *h_darray;
	biguint_t h_invmod;
	biguint_t h_N;
 
	//size_t size = sizeof(biguint_t);
	//size_t pitch;

  mpz_t N;
	mpz_t mpztemp;
	mpz_t mpz_invmod; //N^-1 mod (2^(32*SIZE_NUMBER))
	mpz_t mpz_Rinv; // 2^-(32*SIZE_NUMBER) mod N
  mpz_t sigma;
  mpz_t mpz_B1;
  mpz_t mpz_d;
  mpz_t xp;
  mpz_t zp;
  mpz_t xfin;
  mpz_t zfin;
  mpz_t u;
  mpz_t v;
  //mpz_t gcd;
  //mpz_t factor;

	gmp_randstate_t state;
	gmp_randinit_default(state);
	unsigned long seed=time(NULL);
	gmp_randseed_ui(state,seed);

  mpz_init (N);
  mpz_init (sigma);
  mpz_init (mpz_B1);
  mpz_init (mpz_d);
  mpz_init (xp);
  mpz_init (zp);
  mpz_init (xfin);
  mpz_init (zfin);
  mpz_init (u);
  mpz_init (v);
  mpz_init (mpztemp);
  mpz_init (mpz_invmod);
  mpz_init (mpz_Rinv);
  //mpz_init (gcd);
  //mpz_init (factor);

  if (argc < 3)
  {  
		printf ("Error in call function : not enough arguments.\n%s\n",usage);
		exit(1);
	}
	
  mpz_set_str (N, argv[1], 10); // in base 10 
  mpz_set_str (mpz_B1, argv[2], 10);

	argc-=2;
	argv+=2;

	//default values
	number_of_curves=MAX_NUMBER_OF_CURVES;
	device=-1;
	mpz_set_ui(sigma,0);
	
	while (argc > 2 && argv[1][0]=='-')
	{
		if (argv[1][1] == 's')
		{
			mpz_set_str(sigma, argv[2], 10);
			argc-=2;
			argv+=2;
		}
		else if (argv[1][1] == 'n')
		{
			sscanf(argv[2],"%u",&number_of_curves);
			argc-=2;
			argv+=2;
		}
		else if (argv[1][1] == 'd')
		{
			sscanf(argv[2],"%d",&device);
			argc-=2;
			argv+=2;
		}
	}
	
	if (argc!=1)
	{
			printf ("Error in call function : wrong number of arguments.\n%s\n",usage);
			exit(1);
	}

	if (number_of_curves>MAX_NUMBER_OF_CURVES)
		number_of_curves=MAX_NUMBER_OF_CURVES;
	if (number_of_curves==0)
		number_of_curves=1;

	if (mpz_cmp_ui(sigma,0)==0)
	{
  	mpz_urandomb(sigma,state,3);
		mpz_add_ui(sigma,sigma,6);
	}
	gmp_randclear (state);
	
  // check N is odd 
  if (mpz_divisible_ui_p (N, 2))
  {
  	fprintf (stderr, "Error, N should be odd\n");
   	exit (1);
  }

	mpz_ui_pow_ui(mpztemp,2,32*SIZE_NUMBER);	
	if (mpz_cmp(N,mpztemp) >=0)
  {
  	fprintf (stderr, "Error, N should be stricly lower than 2^%d\n",32*SIZE_NUMBER);
   	exit (1);
  }

	if (mpz_cmp_ui(mpz_B1,TWO32) >=0)
  {
  	fprintf (stderr, "Error, B1 should be stricly lower than 2^%d\n",32);
   	exit (1);
  }
	B1=(unsigned int)mpz_get_ui(mpz_B1);

	gmp_printf ("#gpu_ecm launched with :\nN=%Zd\nB1=%u\ncurves=%u\nfirstsigma=%Zd\n\n",N,B1,number_of_curves,sigma);

	h_xarray=(biguint_t *) malloc(number_of_curves*sizeof(biguint_t));
	h_zarray=(biguint_t *) malloc(number_of_curves*sizeof(biguint_t));
	h_darray=(biguint_t *) malloc(number_of_curves*sizeof(biguint_t));



	//Some precomputation
	//Compute N^-1 mod 2^(32*SIZE_NUMBER)
	mpz_invert(mpz_invmod,N,mpztemp);
	mpz_sub(mpz_invmod,mpztemp,mpz_invmod);
	//Compute  2^-(32*SIZE_NUMBER) mod N
	mpz_invert(mpz_Rinv,mpztemp,N);
	
	mpz_to_biguint(h_N,N);	
	mpz_to_biguint(h_invmod,mpz_invmod);	
	
	for(i=0;i<number_of_curves;i++)
	{
		calculParam (sigma, N, mpz_d, xp, zp, u, v);

		//Compute the Montgomery representation of x0, z0 and d
		mpz_mul_2exp(xp,xp,32*SIZE_NUMBER);
		mpz_mod(xp,xp,N);
#ifndef TEST
		mpz_to_biguint(h_xarray[i],xp);	
#endif
#ifdef TEST
		for (j=0;j<32;j++)
			h_xarray[i][j]=0;
	/*h_xarray[i][0]=987562;
	h_xarray[i][1]=655364;
	h_xarray[i][2]=25445;
	h_xarray[i][3]=6555436;
	h_xarray[i][4]=65574436;
	h_xarray[i][5]=657536;
	h_xarray[i][6]=6555536;
	h_xarray[i][7]=65536;
	h_xarray[i][8]=6555436;
	h_xarray[i][9]=675536;
	h_xarray[i][10]=655836;
	h_xarray[i][11]=6554536;
	h_xarray[i][12]=4231;
	h_xarray[i][13]=1237;
	h_xarray[i][14]=789;
	h_xarray[i][15]=6554536;
	h_xarray[i][16]=654536;
	h_xarray[i][17]=65525536;
	h_xarray[i][18]=6565536;
	h_xarray[i][19]=675536;
	h_xarray[i][20]=5536;
	h_xarray[i][21]=54536;
	*/
	h_xarray[i][0]=2;
	h_xarray[i][1]=3;
	h_xarray[i][16]=4;
	h_xarray[i][17]=5;
#endif

		mpz_mul_2exp(zp,zp,32*SIZE_NUMBER);
		mpz_mod(zp,zp,N);
		mpz_to_biguint(h_zarray[i],zp);	

		mpz_mul_2exp(mpz_d,mpz_d,32*SIZE_NUMBER);
		mpz_mod(mpz_d,mpz_d,N);
		mpz_to_biguint(h_darray[i],mpz_d);	

		mpz_add_ui(sigma,sigma,1);

	}	

	mpz_sub_ui(sigma,sigma,number_of_curves);

#ifdef TEST
	biguint_to_mpz(mpztemp,h_xarray[0]);
	biguint_print(h_xarray[0]);
	printf("\n");
#endif

	printf("#Begin GPU computation...\n");
	time_gpu=cuda_Main(h_invmod,h_N,h_xarray,h_zarray,h_darray,B1,number_of_curves,device);
	printf("#All kernels finished, analysing results...\n");

	for(i=0;i<number_of_curves;i++)
	{
		biguint_to_mpz(xfin,h_xarray[i]);	
		mpz_mul(xfin,xfin,mpz_Rinv);
		mpz_mod(xfin,xfin,N);

		biguint_to_mpz(zfin,h_zarray[i]);	
		mpz_mul(zfin,zfin,mpz_Rinv);
		mpz_mod(zfin,zfin,N);

		gmp_printf("#Looking for factors for the curves with sigma=%Zd\n",sigma);
		nbfactor+=findfactor(N,xfin,zfin);

		mpz_add_ui(sigma,sigma,1);

#ifdef TEST
		if (i==0)
		{
			printf ("a=");
			biguint_print(h_xarray[i]);
			printf ("+(");
			biguint_print(h_zarray[i]);
			printf (")*2^1024\n");
			mpz_mul(mpztemp,mpztemp,mpztemp);
			gmp_printf("b=%Zd\n",mpztemp);
			printf("print a;print b; print a-b\n");
		}
#endif

	}

	if (nbfactor==0)
		printf("#Results : No factor found\n");
	else if (nbfactor==1)
		printf("#Results : 1 factor found\n");
	else
		printf("#Results : %u factor(s) found (not necessarily different)\n",nbfactor);
	
	printf("\n#Time gpu : %.3f init&copy=%.3f computation=%.3f\n",(double)(time_gpu.init+time_gpu.computation)/CLOCKS_PER_SEC,(double)(time_gpu.init)/CLOCKS_PER_SEC,(double)(time_gpu.computation)/CLOCKS_PER_SEC);

	mpz_clear (sigma);
  mpz_clear (N);
  mpz_clear (mpztemp);
  mpz_clear (mpz_invmod);
  mpz_clear (mpz_B1);
  mpz_clear (mpz_d);
  mpz_clear (xp);
  mpz_clear (zp);
  mpz_clear (xfin);
  mpz_clear (zfin);
  mpz_clear (u);
  mpz_clear (v);
  mpz_clear (mpz_Rinv);
  
	
	free(h_xarray);
	free(h_zarray);
	free(h_darray);
	
  return 0;
}



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

	// tester si pgcd =N et =0
	mpz_gcd(gcd,zfin,N);
	
	if (mpz_cmp_ui(gcd,1)==0)
	{
		mpz_invert(zfin,zfin,N);
		mpz_mul(xfin,xfin,zfin);
		mpz_mod(xfin,xfin,N);
		mpz_gcd(gcd,xfin,N);
			
		if (mpz_cmp_ui(gcd,1)==0)
		{
			gmp_printf("  #No factors found. You shoud try with a bigger B1.\n");
		}
		else if (mpz_cmp(gcd,N))
		{
			gmp_printf("  #No factors found. You should try with a smaller B1\n");
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
		gmp_printf("  #No factors found. You should try with a smaller B1\n");
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

//calculParam computes the values of parameters of Suyama's parametrisation
//input : a random integer sigma and N
//output : parameters of Suyama's parametrisation -> a, x0, z0, u and v 
void
calculParam(mpz_t sigma, mpz_t N, mpz_t d, mpz_t x0, mpz_t z0, mpz_t u, mpz_t v)
{
	mpz_t a; 
	mpz_t tmp; 
	mpz_t bezout1;
	mpz_t bezout2; 
	mpz_t gcd; 

  mpz_init(a);
  mpz_init(tmp);
  mpz_init(bezout1);
  mpz_init(bezout2);
  mpz_init(gcd);

	// u 
	mpz_pow_ui(u,sigma,2);
	mpz_sub_ui(u,u,5);
	// v 
	mpz_mul_ui(v,sigma,4);
	// x0 
	mpz_powm_ui(x0,u,3,N);
	// z0 
	mpz_powm_ui(z0,v,3,N);
	// a 
	mpz_sub(a,v,u);
	mpz_pow_ui(a,a,3);
	mpz_mul_ui(tmp,u,3);
	mpz_add(tmp,tmp,v);
	mpz_mul(a,a,tmp);
	mpz_mod(a,a,N);
	mpz_pow_ui(tmp,u,3);
	mpz_mul(tmp,tmp,v);
	mpz_mul_ui(tmp,tmp,4);
	mpz_mod(tmp,tmp,N);
 // set gcd to the greatest common divisor of tmp and N, and in addition
 // set bezout1 and bezout2 to coefficients satisfying 
 //tmp*bezout1 + N*bezout2 = gcd -> bezout1 = (1/tmp)%N  
	mpz_gcdext(gcd,bezout1,bezout2,tmp,N);
  
	mpz_mod(tmp,bezout1,N);
	mpz_mul(a,a,tmp);
	mpz_mod(a,a,N);
	mpz_sub_ui(a,a,2);
	mpz_mod(a,a,N);

	//gmp_printf("a=%Zd\n",a);
	
	// d = (a+2)/4 mod N
	mpz_add_ui(d,a,2);
  while (mpz_divisible_ui_p (d, 4) == 0)
   	mpz_add (d, d, N);
	mpz_divexact_ui(d,d,4);
	mpz_mod(d,d,N);
	
	// calculation of the starting point x = (x0/z0)%N 
	mpz_gcdext(gcd,bezout1,bezout2,z0,N);
  // set gcd to the greatest common divisor of z0 and N, and in addition
  // set bezout1 and bezout2 to coefficients satisfying
  // z0*bezout1 + N*bezout2 = gcd -> bezout1=(1/z0)%N 
	mpz_mod(tmp,bezout1,N);
	mpz_mul(tmp,x0,tmp); // x0/z0 
	mpz_mod(tmp,tmp,N);
        
  // x0 <- x0/z0, z0 <- 1 
  mpz_set (x0, tmp);
	mpz_set_ui (z0, 1);

	mpz_clear(a);
	mpz_clear(tmp);
	mpz_clear(bezout1);
	mpz_clear(bezout2);
	mpz_clear(gcd);
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

