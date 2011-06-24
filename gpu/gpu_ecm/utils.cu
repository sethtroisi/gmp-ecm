#include "main.h"
#include "utils.h"

void usage (void)
{
	printf ("Usage: gpu_ecm [options] B1 < file\n");
	printf ("\nParameters:\n");
	printf ("  N          number to factor\n");
	printf ("  B1         stage 1 bound\n");
	printf ("\nOptions:\n");
	printf ("  -n n       compute on n curves in parallel\n");
	printf ("  -d n       compute on device n\n");
	printf ("  -s d       compute for invd from d to d+number_of_curves \n");
	printf ("  -save file save residues at end of stage 1 in file\n");
	printf ("  -h, --help prints this help and exit\n");
}

unsigned int findfactor(mpz_t N, mpz_t xfin, mpz_t zfin)
{
	//int probprime;
	mpz_t gcd;
	mpz_t factor;

  mpz_init (gcd);
  mpz_init (factor);

	mpz_set_ui(factor,0);

	mpz_gcd(gcd,zfin,N);
	
	if (mpz_cmp_ui(gcd,1)==0)
	{
		mpz_invert(zfin,zfin,N);
		mpz_mul(xfin,xfin,zfin);
		mpz_mod(xfin,xfin,N);
			
  	mpz_clear(gcd);
  	mpz_clear(factor);
		return ECM_NO_FACTOR_FOUND;
	}
	else //gcd !=1 (and gcd>0 because N>0) so we found a factor
	{
		//vérifier si factor=N
		gmp_fprintf(stdout,"********** Factor found in step 1: %Zd\n",gcd);
		if (mpz_cmp(gcd,N)==0)
			fprintf(stdout,"Found input number N\n");

		//TODO ecrire facteur cofacteur is prime facteur es prime cofacteur..
  	mpz_clear(gcd);
  	mpz_clear(factor);
		return ECM_FACTOR_FOUND;
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
 	//mpz_clear(temp);
}

void biguint_print (biguint_t a)
{
  unsigned int i;

  fprintf (stdout,"%u", a[0]);
  for (i = 1; i < SIZE_NUMBER; i++)
    if (a[i]!=0)
			fprintf (stdout,"+%u*2^%u", a[i], 32*i);
  //printf ("\n");
}

void bigint_print (dbigint_t a)
{
  unsigned int i;

  fprintf (stdout,"%d", a[0]);
  for (i = 1; i < SIZE_NUMBER; i++)
    fprintf (stdout,"+%d*2^%u", a[i], 32*i);
  fprintf (stdout,"\n");
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

#define MAX_HEIGHT 30

void compute_s (mpz_t s, unsigned int B1)
{
/*
	unsigned int PI=3,pp,ppp;
	pp=2;
	ppp=1;
	while(pp<=B1)
	{
		ppp=pp;
		pp*=2;
	}
	mpz_set_ui(s,ppp);
	//mpz_mul_2exp(s,s,power2);
	
	//for prime >=3
	PI=getprime(PI);
	while (PI<=B1)
	{
		pp=PI;
		ppp=1;
		while(pp<=B1)
		{
			ppp=pp;
			pp*=PI;

		}
		mpz_mul_ui(s,s,ppp);
		PI=getprime(PI);
	}	
*/
	mpz_t l[MAX_HEIGHT];
	mpz_t r[MAX_HEIGHT];
	unsigned int i, j;
	unsigned long pi=2, pp;
	for (i=0;i<MAX_HEIGHT;i++)
	{
		mpz_init(l[i]);
		mpz_init(r[i]);
	}

	i=0;
	while (pi<=B1)
	{
		pp=pi;
		while(pp<=B1)
			pp*=pi;
		pp/=pi;

		if (i%2==0)
			mpz_set_ui(l[0],pp);
		else
			mpz_set_ui(r[0],pp);
			
		j=0;
		while ((i&(1<<j))!=0)
		{
			if ((i&(1<<(j+1)))==0)
				mpz_mul(l[j+1],l[j],r[j]);
			else
				mpz_mul(r[j+1],l[j],r[j]);
			j++;
		}


		//printf(" %u\n",pp);
		//mpz_mul_ui(s,s,pp);

		i++;
		pi=getprime(pi);
	}

	if (i%2==0)
		mpz_set_ui(r[0],1);
	else
		mpz_set_ui(r[0],1);
		
	j=0;
	for (j=0;j<MAX_HEIGHT-1;j++)
	{
		if ((i&(1<<(j)))==0)
			mpz_set(r[j+1],r[j]);
		else
			mpz_mul(r[j+1],l[j],r[j]);
	}
	
	fprintf(stdout,"%u primes under B1\n",i);
	
	mpz_set(s,r[MAX_HEIGHT-1]);
	
	//for (i=0;i<MAX_HEIGHT;i++)
	//	gmp_fprintf("%Zd %Zd\n",l[i],r[i]);

	//gmp_fprintf("%Zd\n",s);

	//mpz_set(s,p[MAX_HEIGHT-1]);

	getprime(0);
	for (i=0;i<MAX_HEIGHT;i++)
	{
		mpz_clear(l[i]);
		mpz_clear(r[i]);
	}
}

void write_resumefile_line (FILE *file, mpz_t N, unsigned int B1, mpz_t xp, 
	unsigned int firstinvd, mpz_t mpz_d)
{
  mpz_t checksum;
  mpz_t A;
  time_t t;
  char text[256];
  char *uname, mname[32];
  
  mpz_init (checksum);
  mpz_init (A);
		
	mpz_mul_ui(A,mpz_d,firstinvd);
	mpz_mod(A,A,N);
	mpz_mul_ui(A,A,4);
	mpz_sub_ui(A,A,2);
	mpz_mod(A,A,N);


  mpz_set_ui (checksum, B1);
  
	gmp_fprintf (file, "METHOD=ECM; A=%Zd",A);
  
	mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (A, CHKSUMMOD));
  
  gmp_fprintf (file, "; B1=%u; N=%Zd", B1, N);
  
	gmp_fprintf (file, "; X=0x%Zx",xp);
  
	mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (N, CHKSUMMOD));
  mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (xp, CHKSUMMOD));
  fprintf (file, "; CHECKSUM=%lu; PROGRAM=GPU-ECM %s;",
           mpz_fdiv_ui (checksum, CHKSUMMOD), VERSION);
  mpz_clear (checksum);
  
 	fprintf (file, " X0=0x2;");
  
  /* Try to get the users and his machines name */
  /* TODO: how to make portable? */
  uname = getenv ("LOGNAME");
  if (uname == NULL)
    uname = getenv ("USERNAME");
  if (uname == NULL)
    uname = "";
  
  if (gethostname (mname, 32) != 0)
    mname[0] = 0;
  mname[31] = 0; /* gethostname() may omit trailing 0 if hostname >31 chars */
  
  if (uname[0] != 0 || mname[0] != 0)
    {
      fprintf (file, " WHO=%.233s@%.32s;", uname, mname);
    }

  t = time (NULL);
  strncpy (text, ctime (&t), 255);
  text[255] = 0;
  text[strlen (text) - 1] = 0; /* Remove newline */
  fprintf (file, " TIME=%s;", text);
  fprintf (file, "\n");
  fflush (file);
}
