/* example of use of libecm.a */

#include <stdio.h>
#include "gmp.h" /* gmp header file */
#include "ecm.h" /* ecm header file */

int
main (int argc, char *argv[])
{
  mpz_t n, f, x, sigma;
  int res;
  double B1 = 1e3; /* stage 1 bound */
  double B2 = 1e5; /* stage 2 bound */
  double B2min = B1; /* lower bound for stage 2 */

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

  printf ("Performing one curve with B1=%1.0f, B2=%1.0f, sigma=%lu\n",
          B1, B2, mpz_get_ui (sigma));

  res = ecm (f, x, sigma, n, NULL, 1.0, B1, B2min, B2, 1.0, ECM_DEFAULT_K,
             1, 0, 0, 0, stdout, stderr);

  if (res)
    {
      printf ("found factor ");
      mpz_out_str (stdout, 10, f);
      printf ("\n");
      printf ("lucky curve was b*y^2 = x^3 + a*x^2 + x\n");
      printf ("with a = (v-u)^3*(3*u+v)/(4*u^3*v)-2,");
      printf (" u = sigma^2-5, v = 4*sigma\n");
    }
  else
    printf ("found no factor\n");

  mpz_clear (sigma);
  mpz_clear (f);
  mpz_clear (x);
  mpz_clear (n);
}
