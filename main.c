/* GMP-ECM -- Integer factorization with ECM and Pollard 'P-1' methods.

  Copyright (C) 2001 Paul Zimmermann,
  LORIA/INRIA Lorraine, zimmerma@loria.fr
  See http://www.loria.fr/~zimmerma/records/ecmnet.html

  This program is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 2 of the License, or (at your
  option) any later version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
  more details.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.
*/

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include "gmp.h"
#include "ecm.h"

/******************************************************************************
*                                                                             *
*                                Main program                                 *
*                                                                             *
******************************************************************************/

int
main (int argc, char *argv[])
{
  mpz_t sigma, n, p, seed;
  double B1, B2, B1cost;
  int result = 1;
  int verbose = 1; /* verbose level */
  int method = EC_METHOD;
  int k = 7; /* default number of blocks in stage 2 */
  int specific_sigma = 0; /* 1=sigma supplied by user, 0=random */
  unsigned int S = 0;
  gmp_randstate_t randstate;

  /* first look for options */
  while ((argc > 1) && (argv[1][0] == '-'))
    {
      if (strcmp (argv[1], "-pm1") == 0)
	{
	  method = PM1_METHOD;
	  argv++;
	  argc--;
	}
      else if (strcmp (argv[1], "-pp1") == 0)
	{
	  method = PP1_METHOD;
	  argv++;
	  argc--;
	}
      else if (strcmp (argv[1], "-q") == 0)
	{
	  verbose = 0;
	  argv++;
	  argc--;
	}
      else if (strcmp (argv[1], "-v") == 0)
	{
	  verbose = 2;
	  argv++;
	  argc--;
	}
      else if ((argc > 2) && (strcmp (argv[1], "-k") == 0))
	{
	  k = atoi(argv[2]);
	  argv += 2;
	  argc -= 2;
	}
      else if ((argc > 2) && (strcmp (argv[1], "-e") == 0))
	{
	  S = atoi(argv[2]);
	  argv += 2;
	  argc -= 2;
	}
      else
	{
	  fprintf (stderr, "Unknown option: %s\n", argv[1]);
	  exit (1);
	}
    }

  if (argc < 2)
    {
      fprintf (stderr, "Usage: ecm B1 [sigma] [B2] < file\n");
      fprintf (stderr, "\nParameters:\n");
      fprintf (stderr, "  B1         stage 1 bound\n");
      fprintf (stderr, "  sigma      elliptic curve seed or generator for P-1 (0 = random)\n");
      fprintf (stderr, "  B2         stage 2 bound\n");
      fprintf (stderr, "\nOptions:\n");
      fprintf (stderr, "  -k n       perform n steps in stage 2\n");
      fprintf (stderr, "  -e n       impose polynomial x^n for Brent-Suyama's extension\n");
      fprintf (stderr, "  -pm1       runs Pollard P-1 instead of ECM\n");
      fprintf (stderr, "  -pp1       runs Williams P+1 instead of ECM\n");
      fprintf (stderr, "  -q         quiet mode\n");
      fprintf (stderr, "  -v         verbose mode\n");
      exit (1);
    }

  if (verbose >= 1)
    {
      printf ("GMP-ECM 5.0 [powered by GMP %s and NTL %u.%u]\n",
	      gmp_version, NTL_major_version (), NTL_minor_version ());
    }

  /* set first stage bound B1 */
  B1 = atof (argv[1]);

  /* check B1 is not too large */
  if (B1 > MAX_B1)
    {
      fprintf (stderr, "Too large stage 1 bound, limit is %f\n", MAX_B1);
      exit (1);
    }

  NTL_init ();

  mpz_init (seed); /* starting point */
  if (method == EC_METHOD) 
    mpz_init (sigma);

  /* set initial point/sigma. If none specified, we'll set default 
     values later */
  if (argc >= 3)
    {
       /* For P-1, this is the initial seed. For ECM, it's the
          sigma or A parameter */
       if (method == PM1_METHOD || method == PP1_METHOD)
         {
           mpz_set_str (seed, argv[2], 10);
           specific_sigma = mpz_sgn (seed) != 0;
         }
       else
         {
           mpz_set_str (sigma, argv[2], 10);
           specific_sigma = mpz_sgn (sigma) != 0; /* zero sigma => random */
         }
    }

  /* We need random numbers without user-specified sigma */
  if (!specific_sigma)
    {
      gmp_randinit_default (randstate);
      gmp_randseed_ui (randstate, time (NULL) + getpid ());
      /* todo: need higher resolution */
    }

  /* set second stage bound B2: when using polynomial multiplication of
     complexity n^alpha, stage 2 has complexity about B2^(alpha/2), and
     we want stage 2 to take about half of stage 1, thus we choose
     B2 = (c*B1)^(2/alpha). Experimentally, c=1/4 seems to work well.
     For Toom-Cook 3, this gives alpha=log(5)/log(3), and B2 ~ (c*B1)^1.365.
     For Toom-Cook 4, this gives alpha=log(7)/log(4), and B2 ~ (c*B1)^1.424. */
  B1cost = (double) B1 / 6.0; /* default for P-1 */
  if (method == PP1_METHOD)
    B1cost *= 2.0;
  B2 = (argc >= 4) ? atof (argv[3]) : pow (B1cost, 1.424828748);

  /* set initial starting point for ECM */
  if (argc >= 5 && method == EC_METHOD)
    mpz_set_str (seed, argv[4], 10);
  
  /* set default Brent-Suyama's exponent */
  if (S == 0)
    S = 1;
  if (method == PP1_METHOD && S > 1)
    {
      printf ("Warning: Brent-Suyama's extension does not work with P+1, using x^1\n");
      S = 1;
    }

  if (verbose >= 1)
    {
      if (method == PM1_METHOD)
        printf ("Pollard P-1");
      else if (method == PP1_METHOD)
        printf ("Williams P+1");
      else
        printf ("Elliptic Curve");
      printf (" method with B1=%1.0f, B2=%1.0f, x^%u\n", B1, B2, S);
    }

  mpz_init (n); /* number(s) to factor */
  mpz_init (p); /* seed/factor found */
  /* loop for number in standard input or file */
  while (feof (stdin) == 0)
    {
      char c = 0;

      /* skip comment lines beginning with # */
      while ((feof (stdin) == 0) && (isdigit (c = getchar ()) == 0))
	{
	  if (c == '#') /* skip end of line */
	    while ((feof (stdin) == 0) && ((c = getchar ()) != '\n'));
	}

      ungetc (c, stdin);

      mpz_inp_str (n, stdin, 0);

      if (verbose > 0)
	{
	  if (mpz_sizeinbase (n, 10) < 1000)
	    {
	      char *str;
	      str = mpz_get_str (NULL, 10, n);
	      printf ("Input number is %s (%u digits)\n", str,
		      (unsigned) strlen (str));
	      fflush (stdout);
	      free (str);
	    }
	  else
	    printf ("Input number has around %u digits\n", (unsigned) 
		    mpz_sizeinbase (n, 10));
	}

      /* Set effective seed/sigma for factoring attempt on this number */
      mpz_set (p, seed);
      if (!specific_sigma)
        {
          if (method == EC_METHOD)
            {
              /* Make random sigma, 0 < sigma <= 2^32 */
              mpz_urandomb (sigma, randstate, 32);
              mpz_add_ui (sigma, sigma, 1);
            }
          else if (method == PP1_METHOD)
            pp1_random_seed (p, n, randstate);
          else if (method == PM1_METHOD)
            pm1_random_seed (p, n, randstate);
        }

      if (method == PM1_METHOD)
        result = pm1 (p, n, B1, B2, k, S, verbose);
      else if (method == PP1_METHOD)
        result = pp1 (p, n, B1, B2, k, S, verbose);
      else
        result = ecm (p, sigma, n, B1, B2, k, S, verbose);

      if (result != 0)
	{
          printf ("********** Factor found in step %u: ", result);
          mpz_out_str (stdout, 10, p);
          printf ("\n");
	  if (mpz_cmp (p, n))
	    {
	      /* prints factor found and cofactor on standard error. */
	      if (mpz_probab_prime_p (p, 25))
		printf ("Found probable prime factor");
	      else printf ("Found composite factor");
	      printf (" of %u digits: ", nb_digits (p));
	      mpz_out_str (stdout, 10, p);
	      printf ("\n");

	      mpz_divexact (n, n, p);
	      if (mpz_probab_prime_p (n, 25) == 0)
		printf ("Composite");
	      else
		printf ("Probable prime");
	      printf (" cofactor ");
	      if (verbose > 0)
		mpz_out_str (stdout, 10, n);
	      printf (" has %u digits", nb_digits(n));
	    }
	  else
	    printf ("Found input number N");
	  printf ("\n");
	  fflush (stdout);
	}
      
      /* if quiet, prints composite cofactors on standard output. */
      if ((verbose == 0) && (mpz_probab_prime_p (n, 25) == 0))
	{
	  mpz_out_str (stdout, 10, n);
	  putchar ('\n');
	  fflush (stdout);
	}

      while ((feof (stdin) == 0) && (isdigit (c = getchar ()) == 0));

      if (feof (stdin) != 0)
	{
	  result = (result) ? 0 : 1;
	  goto end;
	}

      ungetc (c, stdin);
    }

 end:
  NTL_clear ();
  mpz_clear (p);
  mpz_clear (n);
  if (method == EC_METHOD)
    mpz_clear (sigma);
  mpz_clear (seed);

  /* exit 0 iff a factor was found for the last input */
  return result;
}
