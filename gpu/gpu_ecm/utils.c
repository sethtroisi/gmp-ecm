#include "def.h"

#define CHKSUMMOD 4294967291U

#ifdef _MSC_VER
#include "getrusage.h"
#endif

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
#if GMP_NUMB_BITS == 32
    a[i]=mpz_getlimbn(b, i);  
#else // GMP_NUMB_BITS == 64
    if (i%2 == 0)
      a[i]=(mpz_getlimbn(b, i/2) & 0x00000000ffffffff);
    else
      a[i]=(mpz_getlimbn(b, i/2) >> 32);  
#endif
  }
}

void biguint_to_mpz (mpz_t a, biguint_t b)
{
  int i;
  
  mpz_set_ui(a, 0);

  for (i=NB_DIGITS-1;i>=0;i--)
  {
    mpz_mul_2exp(a, a, 32);
	  mpz_add_ui(a , a, b[i]);
  }
}

void write_resumefile2_line (FILE *file, mpz_t N, unsigned int B1, mpz_t xp, 
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
           mpz_fdiv_ui (checksum, CHKSUMMOD), VERSION_GPUECM); 
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
#if 0
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
#endif
