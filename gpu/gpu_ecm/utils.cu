#include "main.h"
#include "utils.h"

void usage (void)
{
  printf ("Usage: gpu_ecm [options] B1 < file\n");
  printf ("\nParameters:\n");
  printf ("  B1         stage 1 bound\n");
  printf ("\nOptions:\n");
  printf ("  -n n       compute on n curves in parallel\n");
  printf ("  -d n       compute on device n\n");
  printf ("  -s d       compute for invd from d to d+number_of_curves \n");
  printf ("  -save file save residues at end of stage 1 in file"
                                                    " ( - means stdout )\n");
  printf ("  -v         verbose\n");
  printf ("  -vv        very verbose\n");
  printf ("  -h, --help prints this help and exit\n");
}

long
cputime ()
{
  struct rusage rus;
  long sec;

  getrusage (RUSAGE_SELF, &rus);
  /* This overflows a 32 bit signed int after 2147483s = 24.85 days */
  sec = (rus.ru_utime.tv_sec+rus.ru_stime.tv_sec) * 1000L;  
  return sec + (rus.ru_utime.tv_usec+rus.ru_stime.tv_usec) / 1000L;
}

unsigned int findfactor(mpz_t N, mpz_t xfin, mpz_t zfin)
{
  mpz_t gcd;
  mpz_init (gcd);

  mpz_gcd(gcd,zfin,N);
  
  if (mpz_cmp_ui(gcd,1)==0)
  {
    mpz_invert(zfin,zfin,N);
    mpz_mul(xfin,xfin,zfin);
    mpz_mod(xfin,xfin,N);
      
    mpz_clear(gcd);
    return ECM_NO_FACTOR_FOUND;
  }
  else //gcd !=1 (and gcd>0 because N>0) so we found a factor
  {
    gmp_fprintf(stdout,"********** Factor found in step 1: %Zd\n",gcd);
    if (mpz_cmp(gcd,N)==0)
      fprintf(stdout,"Found input number N\n");
    else
    {
      mpz_t cofactor;
      mpz_init (cofactor);
    
      mpz_divexact(cofactor, N, gcd);
      
      gmp_fprintf(stdout,"Found %s factor of %u digits: %Zd\n",
              mpz_probab_prime_p (gcd, 5) ? "probable prime" : "composite", 
              mpz_sizeinbase (gcd ,10) , gcd);
      
      gmp_fprintf(stdout,"%s cofactor %Zd has %u digits\n",
              mpz_probab_prime_p (cofactor, 5) ? "Probable prime" : "Composite", 
              cofactor, mpz_sizeinbase (cofactor ,10));
    
      mpz_clear(cofactor);
    }
    mpz_clear(gcd);
    return ECM_FACTOR_FOUND;
  }
}

void biguint_print (biguint_t a)
{
  unsigned int i;

  fprintf (stdout,"%u", a[0]);
  for (i = 1; i < NB_DIGITS; i++)
    if (a[i]!=0)
      fprintf (stdout,"+%u*2^%u", a[i], 32*i);
  //printf ("\n");
}

void bigint_print (dbigint_t a)
{
  unsigned int i;

  fprintf (stdout,"%d", a[0]);
  for (i = 1; i < NB_DIGITS; i++)
    fprintf (stdout,"+%d*2^%u", a[i], 32*i);
  fprintf (stdout,"\n");
}

void mpz_to_biguint (biguint_t a, mpz_t b)
{
  int i;

  for (i=0;i<NB_DIGITS;i++)
  {
    if (i%2 == 0)
      a[i]=(mpz_getlimbn(b, i/2) & 0x00000000ffffffff);
    else
      a[i]=(mpz_getlimbn(b, i/2) >> 32);  
  }
}

void biguint_to_mpz (mpz_t a, biguint_t b)
{
  int i;
  unsigned long temp;
  
  mpz_set_ui(a,0);

  for (i=NB_DIGITS-1;i>=0;i--)
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
  mpz_t acc[MAX_HEIGHT]; // To accumulate products of prime powers 
  unsigned int i, j;
  unsigned long pi = 2, pp, maxpp;

  for (j = 0; j < MAX_HEIGHT; j++)
    mpz_init (acc[j]); // sets acc[j] to 0 

  mpz_set_ui(s, 1);
  i = 0;
  while (pi <= B1)
    {
      pp = pi;
      maxpp = B1 / pi;
      while (pp <= maxpp)
          pp *= pi;

      if ((i & 1) == 0)
          mpz_set_ui (acc[0], pp);
      else
          mpz_mul_ui (acc[0], acc[0], pp);
      
      j = 0;
      while ((i & (1 << j)) != 0)
        {
          if ((i & (1 << (j + 1))) == 0) 
            mpz_swap (acc[j+1], acc[j]); 
          else
            mpz_mul (acc[j+1], acc[j+1], acc[j]); 
          mpz_set_ui (acc[j], 1);
          j++;
        }

      i++;
      pi = getprime (pi);
    }

  for (mpz_set (s, acc[0]), j = 1; mpz_cmp_ui (acc[j], 0) != 0; j++)
    mpz_mul (s, s, acc[j]);
  
  getprime (0); // free the prime tables, and reinitialize 
  
  for (i = 0; i < MAX_HEIGHT; i++)
      mpz_clear (acc[i]);
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
  //fprintf (file, " TIME=%s;", text); temporaire
  fprintf (file, "\n");
  fflush (file);
}

#define IS_NEWLINE(c) (((c) == '\n') || ((c) == '\r'))
#define MAX_LINE 1000

int read_number (mpz_t n, FILE *fd)
{
  int c;
  char line[MAX_LINE];

new_line:
  c = fgetc (fd);

  /* Skip comment lines beginning with '#' */
  if (c == '#')
    {
      do
        c = fgetc (fd);
      while (c != EOF && !IS_NEWLINE(c));
      if (IS_NEWLINE(c))
        goto new_line;
    }

  if (c == EOF)
    return 0;

  ungetc (c, fd);
  
  if (fgets(line, sizeof line, fd) == NULL || mpz_set_str(n, line, 0)==-1)
    goto new_line;

  return 1;
}

