/* example of use of libecm.a */

#include <stdio.h>
#include "gmp.h" /* gmp header file */
#include "ecm.h" /* ecm header file */

int
main (int argc, char *argv[])
{
  mpz_t n, f, x, sigma;
  int res;

  if (argc != 2)
    {
      fprintf (stderr, "Usage: ecmfactor <number>\n");
      exit (1);
    }

  mpz_init (n);

  /* read number on command line */
  if (mpz_set_str (n, argv[1], 10))
    {
      fprintf (stderr, "Invalid number: %s\n", argv[1]);
      exit (1);
    }

  mpz_init (f); /* for potential factor */
  mpz_init (x); /* starting point */
  mpz_init (sigma);

  /* set sigma at random */
  mpz_set_ui (sigma, getpid ());

  res = ecm (f, x, sigma, n, NULL, 1.0, 1e3, 1e3, 1e5, 1.0, ECM_DEFAULT_K,
             1, 0, 0, 0);

  if (res)
    {
      printf ("found factor ");
      mpz_out_str (stdout, 10, f);
      printf ("\n");
      printf ("lucky curve was g*y^2 = x^3 + a*x^2 + x\n");
      printf ("with a = (v-u)^3*(3*u+v),");
      printf (" u = sigma^2-5, v = 4*sigma, sigma = ");
      mpz_dump (sigma);

    }
  else
    printf ("found no factor\n");

  mpz_clear (sigma);
  mpz_clear (f);
  mpz_clear (x);
  mpz_clear (n);
}
