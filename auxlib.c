#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "ecm-impl.h"

unsigned int
gcd (unsigned int a, unsigned int b)
{
  unsigned int t;

  while (b != 0)
    {
      t = a % b;
      a = b;
      b = t;
    }

  return a;
}

void 
mpz_sub_si (mpz_t r, mpz_t s, int i)
{
  if (i >= 0)
    mpz_sub_ui (r, s, (unsigned int) i);
  else
    mpz_add_ui (r, s, (unsigned int) (-i));
}

/* Divide RS by 3 */
void
mpz_divby3_1op (mpz_t RS)
{
  mp_size_t abssize = mpz_size (RS);
  
  if (abssize == 0)
    return;
  
  if (mpn_divexact_by3 (RS->_mp_d, RS->_mp_d, abssize) != 0)
    {
      fprintf (stderr, "mpz_divby3_1op: division by 3 left remainder");
      exit (EXIT_FAILURE);
    }

  if (RS->_mp_d[abssize - 1] == 0)
    RS->_mp_size -= mpz_sgn (RS);
}

/* returns ceil(log(n)/log(2)) */
unsigned int
ceil_log2 (unsigned int n)
{
  unsigned int k;

  /* f(1)=0, f(n)=1+f((n+1)/2) */
  for (k = 0; n > 1; n = (n + 1) / 2, k++);
  return k;
}

/* malloc with elementary error checking */
void *
xmalloc (size_t size)
{
  void *p;
  
  p = malloc (size);
  if (p == NULL)
    {
      fprintf (stderr, "Could not allocate %lu bytes\n", (unsigned long) size);
      exit (EXIT_FAILURE);
    }
  
  return p;
}

/* Return user CPU time measured in milliseconds. Thanks to Torbjorn. */
#if defined (ANSIONLY) || defined (USG) || defined (__SVR4) || defined (_UNICOS) || defined(__hpux)

#if defined (__MINGW32__) || defined (_MSC_VER)

static int
cputime_x (void)
{
  return (int) ((double) clock () * 1000 / CLOCKS_PER_SEC);
}

#include <windows.h>

int
cputime ()
{
  static int First = 1;
  static LARGE_INTEGER PF;
  LARGE_INTEGER i;
  double d;

  if (First == 1)
  {
    First = 0;
    if (!QueryPerformanceFrequency (&PF))
        First = -1;
  }
  if (First == -1)
    return cputime_x ();

  QueryPerformanceCounter (&i);
  d = (double)*(__int64*)&i;
  d /= *(__int64*)&PF;
  d *= 1000;

  /* NOTE a double converting to int is wrong!.  We need the number mod
     2^31-1 (which is correct auto type-conversion from a unsigned ) */
  /* The "other" cputime() functions probably also have this "bug" */
  return (unsigned) d;
}

#else

int
cputime ()
{
  return (int) ((double) clock () * 1000 / CLOCKS_PER_SEC);
}

#endif  /* __MINGW32___ or VC */

#else	/* ANSIONLY and others */

#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>

int
cputime ()
{
  struct rusage rus;

  getrusage (RUSAGE_SELF, &rus);
  return rus.ru_utime.tv_sec * 1000 + rus.ru_utime.tv_usec / 1000;
}
#endif



