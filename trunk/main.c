/* GMP-ECM -- Integer factorization with ECM and Pollard 'P-1' methods.

  Copyright 2001, 2002, 2003 Alexander Kruppa and Paul Zimmermann.

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

/* people keeping track of champions and corresponding url's */
unsigned int champion_digits[3] = { 51, 37, 34 };
char *champion_keeper[3] =
{ "Richard Brent <rpb@comlab.ox.ac.uk>",
  "Paul Zimmermann <zimmerma@loria.fr>",
  "Paul Zimmermann <zimmerma@loria.fr>"};
char *champion_url[3] =
{"ftp://ftp.comlab.ox.ac.uk/pub/Documents/techpapers/Richard.Brent/champs.txt",
 "http://www.loria.fr/~zimmerma/records/Pminus1.html",
 "http://www.loria.fr/~zimmerma/records/Pplus1.html"};

/******************************************************************************
*                                                                             *
*                                Main program                                 *
*                                                                             *
******************************************************************************/

int
main (int argc, char *argv[])
{
  mpz_t sigma, n, seed, f;
  mpq_t rat_seed;
  mpres_t p;
  double B1, B1done, B1cost, B2, B2min;
  int result = 1;
  int verbose = 1; /* verbose level */
  int method = EC_METHOD;
  int repr = 0;
  int k = 8; /* default number of blocks in stage 2 */
  int specific_sigma = 0; /* 1=sigma supplied by user, 0=random */
  int factor_is_prime, cofactor_is_prime;
        /* If a factor is found, indicates whether factor, cofactor are prime. 
           If no factor was found, both are zero. */
  int S = 0, use_dickson = 0;
        /* Degree for Brent-Suyama extension requested by user */
  gmp_randstate_t randstate;
  char *savefilename = NULL;
  char *endptr[1]; /* to parse B2 or B2min-B2max */
  FILE *savefile = NULL;

#ifdef MEMORY_DEBUG
  tests_memory_start ();
#endif

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
      else if (strcmp (argv[1], "-mpzmod") == 0)
        {
          repr = 1;
	  argv++;
	  argc--;
        }
      else if (strcmp (argv[1], "-modmuln") == 0)
        {
          repr = 2;
	  argv++;
	  argc--;
        }
      else if (strcmp (argv[1], "-redc") == 0)
        {
          repr = 3;
	  argv++;
	  argc--;
        }
      else if (strcmp (argv[1], "-power") == 0)
        {
          use_dickson = 0;
          S = atoi(argv[2]);
	  argv += 2;
	  argc -= 2;
        }
      else if (strcmp (argv[1], "-dickson") == 0)
        {
          use_dickson = 1;
          S = atoi(argv[2]);
	  argv += 2;
	  argc -= 2;
        }
      else if ((argc > 2) && (strcmp (argv[1], "-k") == 0))
	{
	  k = atoi(argv[2]);
	  argv += 2;
	  argc -= 2;
	}
      else if ((argc > 2) && (strcmp (argv[1], "-save") == 0))
	{
	  savefilename = argv[2];
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
      fprintf (stderr, "Usage: ecm B1 [sigma] [[B2min-]B2] < file\n");
      fprintf (stderr, "\nParameters:\n");
      fprintf (stderr, "  B1          stage 1 bound\n");
      fprintf (stderr, "  sigma       elliptic curve seed or generator for P-1 (0 = random)\n");
      fprintf (stderr, "  B2          stage 2 bound (or interval B2min-B2max)\n");
      fprintf (stderr, "\nOptions:\n");
      fprintf (stderr, "  -k n        perform n steps in stage 2\n");
      fprintf (stderr, "  -power n    use x^n for Brent-Suyama's extension\n");
      fprintf (stderr, "  -dickson n  use n-th Dickson's polynomial for Brent-Suyama's extension\n");
      fprintf (stderr, "  -pm1        runs Pollard P-1 instead of ECM\n");
      fprintf (stderr, "  -pp1        runs Williams P+1 instead of ECM\n");
      fprintf (stderr, "  -q          quiet mode\n");
      fprintf (stderr, "  -v          verbose mode\n");
      fprintf (stderr, "  -mpzmod     use GMP's mpz_mod for mod reduction\n");
      fprintf (stderr, "  -modmuln    use Montgomery's MODMULN for mod reduction\n");
      fprintf (stderr, "  -redc       use Montgomery's REDC for mod reduction\n");
      fprintf (stderr, "  -save file  save residues at end of stage 1 to file\n");
      exit (1);
    }

  if (verbose >= 1)
    {
      printf ("GMP-ECM %s [powered by GMP %s", ECM_VERSION, gmp_version);
#ifndef POLYEVAL      
      printf (" and NTL %u.%u", NTL_major_version (), NTL_minor_version ());
#endif
      printf ("]\n");
    }

  /* set first stage bound B1 */
  B1 = strtod (argv[1], &argv[1]);
  if (*argv[1] == '-')
    {
      B1done = B1;
      B1 = strtod (argv[1] + 1, NULL);
    }
  else
    B1done = 1.0;
  B2min = B1;

  if (B1 < 0 || B1done < 0)
    {
      fprintf(stderr, "Bound values must be positive\n");
      exit(EXIT_FAILURE);
    }

  /* check B1 is not too large */
  if (B1 > MAX_B1)
    {
      fprintf (stderr, "Too large stage 1 bound, limit is %f\n", MAX_B1);
      exit (1);
    }

#ifndef POLYEVAL
  NTL_init ();
#endif

  mpz_init (seed); /* starting point */
  mpq_init (rat_seed);
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
           if (mpq_set_str (rat_seed, argv[2], 10))
             {
               fprintf (stderr, "Error, invalid seed: %s\n", argv[2]);
               exit (1);
             }
           specific_sigma = mpq_sgn (rat_seed) != 0;
         }
       else
         {
           if (mpz_set_str (sigma, argv[2], 10))
             {
               fprintf (stderr, "Error, invalid sigma: %s\n", argv[2]);
               exit (1);
             }
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
  /* Since nai"ve P+1 and ECM cost respectively 2 and 11 multiplies per
     addition and duplicate, and both are optimized with PRAC, we can
     assume the ratio remains about 11/2. */
  if (method == PP1_METHOD)
    B1cost *= 2.0;
  else if (method == EC_METHOD)
    B1cost *= 11.0;

  /* parse B2 or B2min-B2max */
  if (argc >= 4)
    {
      B2 = strtod (argv[3], endptr);
      if (*endptr == argv[3])
        {
          fprintf (stderr, "Error: B2 or B2min-B2max expected: %s\n", argv[3]);
          exit (1);
        }
      if (**endptr == '-')
        {
          B2min = B2;
          B2 = atof (*endptr + 1);
        }
    }
  else
    B2 = pow (B1cost, 1.424828748);

  /* set initial starting point for ECM */
  if (argc >= 5 && method == EC_METHOD)
    mpz_set_str (seed, argv[4], 10);
  
  /* set default Brent-Suyama's exponent */
  if (S == 0)
    S = 1;
  if (use_dickson) /* Stage 2 interprets negative degree as Dickson poly */
    S = - abs (S);
  
  if (verbose >= 1)
    {
      if (method == PM1_METHOD)
        printf ("Pollard P-1");
      else if (method == PP1_METHOD)
        printf ("Williams P+1");
      else
        printf ("Elliptic Curve");
      printf (" Method with ");
      if (B1done == 1.0)
        printf("B1=%1.0f", B1);
      else
        printf("B1=%1.0f-%1.0f", B1done, B1);
      if (B2min <= B1)
        printf(", B2=%1.0f, ", B2);
      else
        printf(", B2=%1.0f-%1.0f, ", B2min, B2);
      if (S > 0)
        printf("x^%u\n", S);
      else
        printf("Dickson(%u)\n", -S);
    }

  /* Open save file for writing, if saving is requested */
  if (savefilename != NULL)
    {
      savefile = fopen (savefilename, "w");
      if (savefile == NULL)
        {
          fprintf (stderr, "Could not open file %s for writing\n", savefilename);
          exit (EXIT_FAILURE);
        }
    }

  mpz_init (n); /* number(s) to factor */
  mpz_init (f); /* factor found */
  /* loop for number in standard input or file */
  while (feof (stdin) == 0)
    {
      int c = 0;

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
	      __gmp_free_func (str, strlen (str) + 1);
	    }
	  else
	    printf ("Input number has around %u digits\n", (unsigned) 
		    mpz_sizeinbase (n, 10));
	}

      mpz_init (p); /* seed/stage 1 residue */
      factor_is_prime = cofactor_is_prime = 0;

      /* Set effective seed/sigma for factoring attempt on this number */
      if (method != EC_METHOD) /* convert rational seed to integer */
        {
          mpz_t inv;

          mpz_init (inv);
          mpz_invert (inv, mpq_denref (rat_seed), n);
          mpz_mul (inv, mpq_numref (rat_seed), inv);
          mpz_mod (seed, inv, n);
          mpz_clear (inv);
        }
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
        result = pm1 (f, p, n, B1done, B1, B2min, B2, k, S, verbose, repr);
      else if (method == PP1_METHOD)
        result = pp1 (f, p, n, B1done, B1, B2min, B2, k, S, verbose, repr);
      else
        result = ecm (f, p, sigma, n, B1done, B1, B2min, B2, k, S, verbose, repr);

      if (result != 0)
	{
          printf ("********** Factor found in step %u: ", result);
          mpz_out_str (stdout, 10, f);
          printf ("\n");
	  if (mpz_cmp (f, n))
	    {
	      /* prints factor found and cofactor on standard error. */
	      factor_is_prime = mpz_probab_prime_p (f, 25);
	      printf ("Found %s factor of %u digits: ", 
	              factor_is_prime ? "probable prime" : "composite",
	              nb_digits (f));
	      mpz_out_str (stdout, 10, f);
	      printf ("\n");

	      mpz_divexact (n, n, f);
	      cofactor_is_prime = mpz_probab_prime_p (n, 25);
	      printf ("%s cofactor ", 
	              cofactor_is_prime ? "Probable prime" : "Composite");
	      if (verbose > 0)
		mpz_out_str (stdout, 10, n);
	      printf (" has %u digits\n", nb_digits (n));
	      
              /* check for champions (top ten for each method) */
	      if (factor_is_prime && nb_digits (f) >= champion_digits[method])
                    {
                      printf ("Report your potential champion to %s\n",
                              champion_keeper[method]);
                      printf("(see %s)\n", champion_url[method]);
                    }
            }
	  else
	    printf ("Found input number N\n");
	  fflush (stdout);
	}
      
      /* if quiet, prints composite cofactors on standard output. */
      if ((verbose == 0) && !cofactor_is_prime)
	{
	  mpz_out_str (stdout, 10, n);
	  putchar ('\n');
	  fflush (stdout);
	}

#define SAVEFILE
#ifdef SAVEFILE
      /* Write composite cofactors to savefile if requested */
      /* If no factor was found, we consider cofactor composite and write it */
      if (savefile != NULL && !cofactor_is_prime)
        {
          mpz_t checksum;
          mpz_init (checksum);
          mpz_set_d (checksum, B1);
          mpz_mod (p, p, n); /* Reduce stage 1 residue wrt new cofactor, in
                                case a factor was found */
          fprintf (savefile, "METHOD=");
          if (method == PM1_METHOD)
            fprintf (savefile, "PM1");
          else if (method == PP1_METHOD)
            fprintf (savefile, "PP1");
          else 
            {
              fprintf (savefile, "ECM; SIGMA=");
              mpz_out_str (savefile, 10, sigma);
              mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (sigma, CHKSUMMOD));
            }
          
          fprintf (savefile, "; B1=%.0f; N=", B1);
          mpz_out_str (savefile, 10, n);
          fprintf (savefile, "; X=0x");
          mpz_out_str (savefile, 16, p);
          mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (n, CHKSUMMOD));
          mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (p, CHKSUMMOD));
          fprintf (savefile, "; CHECKSUM=%lu; PROGRAM=GMP-ECM ver. 5 alpha;\n", 
                   mpz_fdiv_ui (checksum, CHKSUMMOD));
          mpz_clear (checksum);
          
        }
#endif /* SAVEFILE */

      mpz_clear (p);

      while ((feof (stdin) == 0) && (isdigit (c = getchar ()) == 0));

      if (feof (stdin) != 0)
	{
	  result = (result) ? 0 : 1;
	  goto end;
	}

      ungetc (c, stdin);
    }

 end:
#ifndef POLYEVAL
  NTL_clear ();
#endif
  mpz_clear (f);
  mpz_clear (n);
  if (method == EC_METHOD)
    mpz_clear (sigma);
  mpz_clear (seed);
  mpq_clear (rat_seed);

  if (!specific_sigma)
    gmp_randclear (randstate);

#ifdef MEMORY_DEBUG
  tests_memory_end ();
#endif

  /* exit 0 iff a factor was found for the last input */
  return result;
}
