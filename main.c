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
#include "gmp.h"
#include "gmp-impl.h"
#include "ecm.h"

/******************************************************************************
*                                                                             *
*                                Main program                                 *
*                                                                             *
******************************************************************************/

/* maximal stage 1 bound = 2^53 + 4, the next prime being 2^53 + 5 */
#define MAX_B1 9007199254740996.0

int
main (int argc, char *argv[])
{
  mpz_t sigma, n, p;
  double B1, B2;
  int result = 1;
  int verbose = 1; /* verbose level */
  int (*factor) _PROTO((mpz_t, mpz_t, double, double, unsigned int, int)) 
    = ecm;
  int k = 3;

  /* first look for options */
  while ((argc > 1) && (argv[1][0] == '-'))
    {
      if (strcmp (argv[1], "-pm1") == 0)
	{
	  factor = pm1;
	  argv ++;
	  argc --;
	}
      else if (strcmp (argv[1], "-q") == 0)
	{
	  verbose = 0;
	  argv ++;
	  argc --;
	}
      else if (strcmp (argv[1], "-v") == 0)
	{
	  verbose = 2;
	  argv ++;
	  argc --;
	}
      else if ((argc > 2) && (strcmp (argv[1], "-k") == 0))
	{
	  k = atoi(argv[2]);
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
      fprintf (stderr, "  -pm1       runs Pollard P-1 instead of ECM\n");
      fprintf (stderr, "  -q         quiet mode\n");
      fprintf (stderr, "  -v         verbose mode\n");
      exit (1);
    }

  if (verbose >= 1)
    printf ("GMP-ECM 5.0 [powered by GMP %u.%u.%u and NTL %u.%u]\n",
            __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR,
            __GNU_MP_VERSION_PATCHLEVEL, NTL_major_version (),
            NTL_minor_version ());

  /* set first stage bound B1 */
  B1 = atof (argv[1]);

  /* check B1 is not too large */
  if (B1 > MAX_B1)
    {
      fprintf (stderr, "Too large stage 1 bound, limit is %f\n", MAX_B1);
      exit (1);
    }

  NTL_init ();

  /* set initial seed sigma */
  mpz_init (sigma);
  if (argc >= 3)
    mpz_set_str (sigma, argv[2], 10);
  if (mpz_cmp_ui (sigma, 0) == 0)
    mpz_set_ui (sigma, 17); /* default is 17 */

  /* set second stage bound B2 */
  B2 = (argc >= 4) ? atof (argv[3]) : 100.0 * (double) B1;

  if (verbose >= 1)
    {
      if (factor == pm1)
        printf ("Pollard P-1");
      else
        printf ("Elliptic Curve");
      printf (" method with B1=%1.0f, B2=%1.0f\n", B1, B2);
    }

  mpz_init (n); /* number(s) to factor */
  mpz_init (p); /* found prime */
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

      if (verbose >= 1)
        {
          printf ("Using sigma=");
          mpz_out_str (stdout, 10, sigma);
          printf ("\n");
        }

      mpz_set (p, sigma);
      if ((result = factor (p, n, B1, B2, k, verbose)))
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
  mpz_clear (sigma);

  /* exit 0 iff a factor was found for the last input */
  return result;
}
