/* example of use of libecm.a */

#include <stdio.h>
#include <stdlib.h>
#include "gmp.h" /* gmp header file */
#include "ecm.h" /* ecm header file */

int
main (int argc, char *argv[])
{
  mpz_t n, f;
  int res;
  double B1;

  if (argc != 3)
    {
      fprintf (stderr, "Usage: ecmfactor <number> <B1>\n");
      exit (1);
    }

  mpz_init (n);

  /* read number on command line */
  if (mpz_set_str (n, argv[1], 10))
    {
      fprintf (stderr, "Invalid number: %s\n", argv[1]);
      exit (1);
    }

  B1 = atof (argv[2]);

  mpz_init (f); /* for potential factor */

  printf ("Performing one curve with B1=%1.0f\n", B1);

  res = ecm_factor (f, n, B1, NULL);

  if (res > 0)
    {
      printf ("found factor in step %u: ", res);
      mpz_out_str (stdout, 10, f);
      printf ("\n");
#if 0
      printf ("lucky curve was b*y^2 = x^3 + a*x^2 + x\n");
      printf ("with a = (v-u)^3*(3*u+v)/(4*u^3*v)-2,");
      printf (" u = sigma^2-5, v = 4*sigma\n");
#endif
    }
  else if (res == ECM_NO_FACTOR_FOUND)
    printf ("found no factor\n");
  else
    printf ("error\n");

  mpz_clear (f);
  mpz_clear (n);

  return 0;
}
