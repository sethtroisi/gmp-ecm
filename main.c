/* GMP-ECM -- Integer factorization with ECM and Pollard 'P-1' methods.

  Copyright 2001, 2002, 2003 Paul Zimmermann and Alexander Kruppa.

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
#include <unistd.h>
#include "gmp.h"
#include "ecm.h"

#if defined (__MINGW32__) || defined (_MSC_VER) /* || defined (__CYGWIN__) */
/* needed for priority setting */               /* not sure about CyGwin and windows.h and Win32 API functions */
#include <windows.h>
#endif

/* #define DEBUG */

/* people keeping track of champions and corresponding url's */
unsigned int champion_digits[3] = { 51, 38, 34 };
char *champion_keeper[3] =
{ "Richard Brent <rpb@comlab.ox.ac.uk>",
  "Paul Zimmermann <zimmerma@loria.fr>",
  "Paul Zimmermann <zimmerma@loria.fr>"};
char *champion_url[3] =
{"ftp://ftp.comlab.ox.ac.uk/pub/Documents/techpapers/Richard.Brent/champs.txt",
 "http://www.loria.fr/~zimmerma/records/Pminus1.html",
 "http://www.loria.fr/~zimmerma/records/Pplus1.html"};

/* Tries to read a number from a line from fd and stores it in r.
   Keeps reading lines until a number is found. Lines beginning with "#"
     are skipped.
   Returns 1 if a number was successfully read, 0 if no number can be read 
     (i.e. at EOF)
   Function is now simpler.  Much of the logic (other than skipping # lines
     is now contained within eval() function.
*/

int 
read_number (mpcandi_t *n, FILE *fd, int primetest)
{
  int c;
  
new_line:
  c = fgetc (fd);
  
  /* Skip comment lines beginning with '#' */
  if (c == '#')
    {
      do
        c = fgetc (fd);
      while (c != EOF && c != '\n');
      if (c == '\n')
        goto new_line;
    }

  if (c == EOF)
    return 0;

  ungetc (c, fd);
  if (!eval (n, fd, primetest))
    goto new_line;

#if 0
  /*  Code to test out eval_str function, which "appears" to work correctly. */
  {
    /* warning!! Line is pretty small, but since this is just testing code, we
       can easily control the input for this test.  This code should NEVER be
       compiled into released build, its only for testing of eval_str() */
    char Line[500], *cp;
    fgets (Line, sizeof(Line), fd);

    if (!eval_str (n, Line, primetest, &cp))
      goto new_line;
    fprintf (stderr, "\nLine is at %X cp is at %X\n", Line, cp);
  }
#endif

#if defined (DEBUG_EVALUATOR)
  if (n->cpExpr)
    fprintf (stderr, "%s\n", n->cpExpr);
  mpz_out_str (stderr, 10, n->n);
  fprintf (stderr, "\n");
#endif

  return 1;
}


/******************************************************************************
*                                                                             *
*                                Main program                                 *
*                                                                             *
******************************************************************************/

int
main (int argc, char *argv[])
{
  mpz_t x, sigma, A, f, orig_x0;
  mpcandi_t n;
  mpq_t rat_x0;
  double B1, B1done, B2, B2min;
  int result = 0;
  int verbose = 1; /* verbose level */
  int method = EC_METHOD, method1;
  int specific_x0 = 0, /* 1=starting point supplied by user, 0=random or */
                       /* compute from sigma */
      specific_sigma = 0;  /* 1=sigma from command line, 0=make random */
  int factor_is_prime;
        /* If a factor was found, indicate whether factor, cofactor are */
        /* prime. If no factor was found, both are zero. */
  int repr = 0;
#ifdef POLYEVAL
  int k = 0;
#else /* POLYGCD is more expensive -> perform more blocks */
  int k = 8; /* default number of blocks in stage 2 */
#endif
  int S = 0; /* Degree for Brent-Suyama extension requested by user */
             /* Positive value: use S-th power, */
             /* negative: use degree |S| Dickson poly */
  gmp_randstate_t randstate;
  char *savefilename = NULL, *resumefilename = NULL, *infilename = NULL;
  char *endptr[1]; /* to parse B2 or B2min-B2max */
  char rtime[256] = "", who[256] = "", comment[256] = "", program[256] = "";
  FILE *savefile = NULL, *resumefile = NULL, *infile = NULL;
  int primetest = 0;
  double autoincrementB1 = 0, startingB1;
  unsigned int autoincrementB1_calc = 0;
  unsigned int breadthfirst_maxcnt=0, breadthfirst_cnt=0;
  int breadthfirst = 0;
  unsigned int count = 1; /* number of curves for each number */
  unsigned int cnt = 0;   /* number of remaining curves for current number */
  unsigned int linenum = 0, factsfound = 0;
  mpcandi_t *pCandidates=NULL;
  unsigned int nCandidates=0, nMaxCandidates=0;
  int deep=1, trial_factor_found;
  unsigned int displayexpr=0;
  double maxtrialdiv=0;

  /* check ecm is linked with a compatible librayr */
  if (mp_bits_per_limb != GMP_NUMB_BITS)
    {
      fprintf (stderr, "Error, mp_bits_per_limb and GMP_NUMB_BITS differ\n");
      fprintf (stderr, "Please check your LD_LIBRARY_PATH variable\n");
      exit (1);
    }

#ifdef MEMORY_DEBUG
  tests_memory_start ();
#endif

  /* Init variables we might need to store options */
  mpz_init (sigma);
  mpz_init (A);
  mpq_init (rat_x0);

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
      else if (strcmp (argv[1], "-nobase2") == 0)
        {
          repr = -1;
	  argv++;
	  argc--;
        }
      else if (strcmp (argv[1], "-primetest") == 0)
        {
          primetest = 1;
	  argv++;
	  argc--;
        }
      else if (strcmp (argv[1], "-one") == 0)
        {
          deep = 0;
	  argv++;
	  argc--;
        }
      else if (strcmp (argv[1], "-b") == 0)
        {
	  breadthfirst = 1;
	  argv++;
	  argc--;
        }
      else if (strcmp (argv[1], "-d") == 0)
        {
	  /* -1 is a flag used during argv processing where a subsquent -i file will NOT change it.  Then
	     when done processing args, we change a -1 to a 0 */
	  breadthfirst = -1;
	  argv++;
	  argc--;
        }
#if defined (__MINGW32__) || (_MSC_VER)
      else if (strcmp (argv[1], "-n") == 0)
        {
/*	  fprintf (stderr, "Executing at lower priority\n"); */
	  SetPriorityClass (GetCurrentProcess (), IDLE_PRIORITY_CLASS);
	  SetThreadPriority (GetCurrentThread (), THREAD_PRIORITY_BELOW_NORMAL);
	  argv++;
	  argc--;
        }
      else if (strcmp (argv[1], "-nn") == 0)
        {
/*	  fprintf (stderr, "Executing at idle priority\n"); */
	  SetPriorityClass (GetCurrentProcess (), IDLE_PRIORITY_CLASS);
	  SetThreadPriority (GetCurrentThread (), THREAD_PRIORITY_IDLE);
	  argv++;
	  argc--;
        }
#elif __CYGWIN__
      else if (strcmp (argv[1], "-n") == 0)
        {
	  /* Not sure what is correct here */
	  nice (20);
	  argv++;
	  argc--;
        }
      else if (strcmp (argv[1], "-nn") == 0)
        {
	  /* Not sure what is correct here */
	  nice (20);
	  argv++;
	  argc--;
        }
#else
      /* The nix boys need to do this (if they want -n or -nn) */
      else if (strcmp (argv[1], "-n") == 0)
        {
	  nice (20);
	  argv++;
	  argc--;
        }
      else if (strcmp (argv[1], "-nn") == 0)
        {
	  nice (20);
	  argv++;
	  argc--;
        }
#endif
      else if ((argc > 2) && (strcmp (argv[1], "-x0")) == 0)
        {
          if (mpq_set_str (rat_x0, argv[2], 0))
            {
              fprintf (stderr, "Error, invalid starting point: %s\n", argv[2]);
              exit (EXIT_FAILURE);
            }
          specific_x0 = 1;
	  argv += 2;
	  argc -= 2;
        }
      else if ((argc > 2) && (strcmp (argv[1], "-sigma")) == 0)
        {
          if (mpz_set_str (sigma, argv[2], 0))
	    {
	      fprintf (stderr, "Error, invalid sigma value: %s\n", argv[2]);
	      exit (EXIT_FAILURE);
	    }
          specific_sigma = 1;
	  argv += 2;
	  argc -= 2;
        }
      else if ((argc > 2) && (strcmp (argv[1], "-A")) == 0)
        {
          if (mpz_set_str (A, argv[2], 0))
	    {
	      fprintf (stderr, "Error, invalid A value: %s\n", argv[2]);
              exit (EXIT_FAILURE);
	    }
	  argv += 2;
	  argc -= 2;
        }
      else if ((argc > 2) && (strcmp (argv[1], "-power")) == 0)
        {
          S = abs (atoi (argv[2]));
	  /* should this be validated? and a error/abort issued if 0 ??? */
	  argv += 2;
	  argc -= 2;
        }
      else if ((argc > 2) && (strcmp (argv[1], "-dickson") == 0))
        {
          S = - abs ( atoi (argv[2]));
	  /* should this be validated? and a error/abort issued if 0 ??? */
	  argv += 2;
	  argc -= 2;
        }
      else if ((argc > 2) && (strcmp (argv[1], "-k") == 0))
	{
	  k = atoi (argv[2]);
	  /* should this be validated? and a error/abort issued if 0 ??? */
	  argv += 2;
	  argc -= 2;
	}
      else if ((argc > 2) && (strcmp (argv[1], "-c") == 0))
	{
	  count = atoi (argv[2]);
	  /* should this be validated? and a error/abort issued if 0 ??? */
	  argv += 2;
	  argc -= 2;
	}
      else if ((argc > 2) && (strcmp (argv[1], "-save") == 0))
	{
	  savefilename = argv[2];
	  argv += 2;
	  argc -= 2;
	}
      else if ((argc > 2) && (strcmp (argv[1], "-resume") == 0))
	{
	  resumefilename = argv[2];
	  argv += 2;
	  argc -= 2;
	}
      else if ((argc > 2) && (strcmp (argv[1], "-i") == 0))
	{
	  autoincrementB1 = strtod (argv[2], NULL);
	  if (autoincrementB1 < 1)
	    {
	      fprintf (stderr, "Error, the -i command requires a whole number argument to follow it\n");
	      exit (EXIT_FAILURE);
  	    }
	  argv += 2;
	  argc -= 2;
	}
      else if ((argc > 2) && (strcmp (argv[1], "-I") == 0))
	{
	  autoincrementB1 = strtod (argv[2], NULL);
	  autoincrementB1_calc = 1;
	  if (autoincrementB1 <= 0)
	    {
	      fprintf (stderr, "Error, the -I command requires a number argument to follow it > 0\n");
	      exit (EXIT_FAILURE);
  	    }
	  argv += 2;
	  argc -= 2;
	}
      else if ((argc > 2) && (strcmp (argv[1], "-inp") == 0))
	{
	  infilename = argv[2];
	  infile = fopen (infilename, "r");
	  /* a -d depth-first switch has already been processed, so DO NOT reset to breadth-first */
	  if (breadthfirst != -1)
	    breadthfirst = 1;
	  if (!infile)
	    {
	      fprintf (stderr, "Can't find input file %s\n", infilename);
	      exit (EXIT_FAILURE);
	    }
	  argv += 2;
	  argc -= 2;
	}
      else if ((argc > 2) && (strcmp (argv[1], "-t") == 0))
	{
	  maxtrialdiv = strtod (argv[2], NULL);
	  if (maxtrialdiv == 0.)
	    {
	      fprintf (stderr, "Error, the -t command requires a number argument to follow it\n");
	      exit (EXIT_FAILURE);
  	    }
	  argv += 2;
	  argc -= 2;
	}
      else if ((argc > 2) && (strcmp (argv[1], "-ve") == 0))
        {
	  displayexpr=atoi (argv[2]);
	  if (displayexpr == 0)
	    {
	      fprintf (stderr, "Error, the -ve command requires a number argument to follow it\n");
	      exit (EXIT_FAILURE);
  	    }
	  argv += 2;
	  argc -= 2;
        }

      else
	{
	  fprintf (stderr, "Unknown option: %s\n", argv[1]);
	  exit (EXIT_FAILURE);
	}
    }

  /* check that S is even for P-1 */
  if ((method == PM1_METHOD) && (S % 2 != 0))
    {
      fprintf (stderr, "Error, S should be even for P-1\n");
      exit (EXIT_FAILURE);
    }

  /* Ok, now we can "reset" the breadthfirst switch so that we do depthfirst as requested */
  if (breadthfirst == -1)
    breadthfirst = 0;

  if (argc < 2)
    {
      fprintf (stderr, "Usage: ecm [options] B1 [[B2min-]B2] < file\n");
      fprintf (stderr, "\nParameters:\n");
      fprintf (stderr, "  B1           stage 1 bound\n");
      fprintf (stderr, "  B2           stage 2 bound (or interval B2min-B2max)\n");
      fprintf (stderr, "\nOptions:\n");
      fprintf (stderr, "  -x0 x        use x as initial point\n"); 
      fprintf (stderr, "  -sigma s     use s as curve generator [ecm]\n");
      fprintf (stderr, "  -A a         use a as curve parameter [ecm]\n");
      fprintf (stderr, "  -k n         perform >= n steps in stage 2\n");
      fprintf (stderr, "  -power n     use x^n for Brent-Suyama's extension\n");
      fprintf (stderr, "  -dickson n   use n-th Dickson's polynomial for Brent-Suyama's extension\n");
      fprintf (stderr, "  -c n         perform n runs for each input\n");
      fprintf (stderr, "  -pm1         perform P-1 instead of ECM\n");
      fprintf (stderr, "  -pp1         perform P+1 instead of ECM\n");
      fprintf (stderr, "  -q           quiet mode\n");
      fprintf (stderr, "  -v           verbose mode\n");
      fprintf (stderr, "  -mpzmod      use GMP's mpz_mod for mod reduction\n");
      fprintf (stderr, "  -modmuln     use Montgomery's MODMULN for mod reduction\n");
      fprintf (stderr, "  -redc        use Montgomery's REDC for mod reduction\n");
      fprintf (stderr, "  -nobase2     disable special base-2 code\n");
      fprintf (stderr, "  -save file   save residues at end of stage 1 to file\n");
      fprintf (stderr, "  -resume file resume residues from file, reads from stdin if file is \"-\"\n");
      fprintf (stderr, "  -primetest   perform a primality test on input\n");
      /*rintf (stderr, "  -extra functions added by JimF\n"); */
      fprintf (stderr, "\n");
      fprintf (stderr, "  Options beyond ECM 5.0  (i.e. specific to ECM 5.0c\n");
      fprintf (stderr, "  -i n         increment B1 by this constant on each run\n");
      fprintf (stderr, "  -I f         auto-calculated increment for B1 multiplied by 'f' scale factor\n");
      fprintf (stderr, "  -inp file    Use file as input (instead of redirecting stdin)\n");
      fprintf (stderr, "  -b           Use breadth-first mode of file processing (recommended)\n");
      fprintf (stderr, "  -d           Use depth-first mode of file processing\n");
      fprintf (stderr, "  -one         Stop processing a candidate when a factor is found (looping mode)\n");
      fprintf (stderr, "  -n           run ecm in \"nice\" mode (below normal priority)\n");
      fprintf (stderr, "  -nn          run ecm in \"very nice\" mode (idle priority)\n");
      fprintf (stderr, "  -t n         Trial divide candidates before P-1, P+1 or ECM up to n\n");
      fprintf (stderr, "  -ve n        Verbosely show short (< n character) expressions on each loop\n");
      exit (EXIT_FAILURE);
    }

  if (verbose >= 1)
    {
      printf ("GMP-ECM %s [powered by GMP %s", ECM_VERSION, gmp_version);
#ifdef POLYGCD
      printf (" and NTL %u.%u", NTL_major_version (), NTL_minor_version ());
#endif
      printf ("] [");
      switch (method)
	{
	case PM1_METHOD:
	  printf ("P-1");
	  break;
	case PP1_METHOD:
	  printf ("P+1");
	  break;
	default:
	  printf ("ECM");
	}
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
      fprintf (stderr, "Bound values must be positive\n");
      exit (EXIT_FAILURE);
    }

  /* check B1 is not too large */
  if (B1 > MAX_B1)
    {
      fprintf (stderr, "Too large stage 1 bound, limit is %1.0f\n", MAX_B1);
      exit (EXIT_FAILURE);
    }

#ifdef POLYGCD
  NTL_init ();
#endif

  init_expr ();

  B2 = 0.0;
  /* parse B2 or B2min-B2max */
  if (argc >= 3)
    {
      B2 = strtod (argv[2], endptr);
      if (*endptr == argv[2])
        {
          fprintf (stderr, "Error: B2 or B2min-B2max expected: %s\n", argv[2]);
          exit (EXIT_FAILURE);
        }
      if (**endptr == '-')
        {
          B2min = B2;
          B2 = atof (*endptr + 1);
        }
    }

  /* Open resume file for reading, if resuming is requested */
  if (resumefilename != NULL)
    {
      if (strcmp (resumefilename, "-") == 0)
        resumefile = stdin;
      else
        resumefile = fopen (resumefilename, "r");
      
      if (resumefile == NULL)
        {
          fprintf (stderr, "Could not open file %s for reading\n", 
                   resumefilename);
          exit (EXIT_FAILURE);
        }
    }

  /* Open save file for writing, if saving is requested */
  /* FIXME: append by default ? */
  if (savefilename != NULL)
    {
      /* Does this file already exist ? */
      if (access (savefilename, F_OK) == 0)
        {
          printf ("Save file %s already exists, will not overwrite\n", 
                  savefilename);
          exit (EXIT_FAILURE);
        }
      savefile = fopen (savefilename, "w");
      if (savefile == NULL)
        {
          fprintf (stderr, "Could not open file %s for writing\n", savefilename);
          exit (EXIT_FAILURE);
        }
    }

  if (resumefile && (specific_sigma || mpz_sgn (A) || specific_x0))
    {
      printf ("Warning: -sigma, -A and -x0 parameters are ignored when resuming from\nsave files.\n");
      mpz_set_ui (sigma, 0);
      mpz_set_ui (A, 0);
      specific_x0 = 0;
    }

  mpcandi_t_init (&n); /* number(s) to factor */
  mpz_init (f); /* factor found */
  mpz_init (x); /* stage 1 residue */
  mpz_init (orig_x0); /* starting point, for save file */

  /* We may need random numbers for sigma/starting point */
  gmp_randinit_default (randstate);
  gmp_randseed_ui (randstate, get_random_ui ());

  /* loop for number in standard input or file */

  startingB1 = B1;
  if (!infilename)
    infile = stdin;

  if (breadthfirst == 1)
    {
      breadthfirst_maxcnt=count;
      count=1;
      breadthfirst_cnt=0;
    }

BreadthFirstDoAgain:;
  if (breadthfirst == 1)
    {
      if (breadthfirst_maxcnt > breadthfirst_cnt)
        {
	  linenum = 0;
	  if (breadthfirst_cnt++)
            {
  	      B1 = calc_B1_AutoIncrement(B1, autoincrementB1, autoincrementB1_calc);
	      B2min = B1;
	    }
	  else
            {
	      /* This is ONLY entered upon the first time through.  We load the entire file here so that we can loop deep, 
		  or remove a candidate if factor found, or if in deep mode and cofactor is prp (or if original candidate
		  is prp and we are prp testing) */
	      nMaxCandidates = 100;
	      pCandidates = malloc (nMaxCandidates*sizeof(mpcandi_t));

	      while (!feof (infile))
		{
		  if (read_number (&n, infile, primetest))
		    {
		      mpcandi_t_init (&pCandidates[nCandidates]);
		      mpcandi_t_copy (&pCandidates[nCandidates++], &n);
		      if (nCandidates == nMaxCandidates)
			{
			    mpcandi_t *tmp = pCandidates;
			    pCandidates = malloc ((nMaxCandidates+100)*sizeof(mpcandi_t));
			    /*	perform a "shallow" copy, in which we do NOT need to free any of the 
				individual elements, but just the array memory */
			    if (pCandidates)
			      memcpy (pCandidates, tmp, nMaxCandidates*sizeof(mpcandi_t));
			    nMaxCandidates += 100;
			    /* Free the original "array" memory */
			    free (tmp);
			}
		    }
		}
	      /*  Now infile is at EOF, but we are in breadthfirst mode, so the main while loop will work with linenum<nCandidates */
	    }
	}
      else
	{
	  breadthfirst = 0;
	}
    }

  while ((breadthfirst && linenum < nCandidates) || feof (infile) == 0)
    {
      trial_factor_found = 0;
      if (resumefile != NULL)
        {
	  if (count != 1)
	    {
	      fprintf (stderr, "Error, option -c and -resume are incompatible\n");
	      exit (EXIT_FAILURE);
	    }

          if (!read_resumefile_line (&method, x, &n, sigma, A, orig_x0, 
                &B1done, program, who, rtime, comment, resumefile))
            break;

	  cnt = count;

          if (verbose > 0)
            {
              printf ("Resuming ");
              if (method == EC_METHOD)
                printf ("ECM");
              else if (method == PM1_METHOD)
                printf ("P-1");
              else if (method == PP1_METHOD)
                printf ("P+1");
              printf (" residue ");
              if (program[0] || who[0] || rtime[0])
                printf ("saved ");
              if (who[0])
                printf ("by %s ", who);
              if (program[0])
                printf ("with %s ", program);
              if (rtime[0])
                printf ("on %s ", rtime);
              if (comment[0])
                printf ("(%s)", comment);
              printf ("\n");
            }
        }
      else
        {
	  if (cnt) /* nothing to read: reuse old number */
	    {
	    }
	  else /* new number */
	    {
	      if (!breadthfirst && !read_number (&n, infile, primetest))
		break;
	      else if (breadthfirst)
		mpcandi_t_copy (&n,&pCandidates[linenum]);
	      linenum++;
	      cnt = count;
	      /* reset B1 value, as it could have been advanced on the prior candidate */
	      if (!breadthfirst)  /* IS THIS A BUG I INTRODUCED!!!!   Can B2min now ever be set??? */
		{
	          B1 = startingB1;
		  B2min = B1;
		}
	    }

	  /* in breadthfirst deep mode, a value of 1 is left after FULLY factoring the number, so we then skip it */
	  /* Also "blank" lines, or lines that could not be parsed correctly will leave a 1 in this value */
	  if (n.isPrp)
	  {
	    /* n is 0 or 1 (or -1 I guess) so do NOT proceed with it */
            cnt = 0;
	    continue;
	  }

          /* Set effective seed for factoring attempt on this number */

          if (specific_x0) /* convert rational value to integer */
            {
              mpz_t inv;

	      if (count != 1)
		{
		  fprintf (stderr, "Error, option -c is incompatible with -x0\n");
		  exit (EXIT_FAILURE);
		}

              mpz_init (inv);
              mpz_invert (inv, mpq_denref (rat_x0), n.n);
              mpz_mul (inv, mpq_numref (rat_x0), inv);
              mpz_mod (x, inv, n.n);
              mpz_clear (inv);
            }
          else /* Make a random starting point for P-1 and P+1. ECM will */
               /* compute a suitable value from sigma or A if x is zero */
            {
              if (method == EC_METHOD)
                mpz_set_ui (x, 0);
              if (method == PP1_METHOD)
                pp1_random_seed (x, n.n, randstate);
              if (method == PM1_METHOD)
                pm1_random_seed (x, n.n, randstate);
            }
         
          if (B1done <= 1.0)
            mpz_set (orig_x0, x);
          
          /* Make a random sigma if we have neither specific sigma nor A 
             given. Warning: sigma may still contain previous random value
             and thus be nonzero here even if no specific sigma was given */
          if (method == EC_METHOD && !specific_sigma && !mpz_sgn (A))
            {
              /* Make random sigma, 0 < sigma <= 2^32 */
              mpz_urandomb (sigma, randstate, 32);
              mpz_add_ui (sigma, sigma, 1); /* FIXME: need sigma>=5? */
            }
        }
      if (count == cnt)
	fprintf (stderr, "\r      Line=%u B1=%.0f factors=%u           \r", linenum, B1, factsfound);
      if (count != cnt)
	fprintf (stderr, "\r      Line=%u Curves=%u/%u B1=%.0f factors=%u\r", linenum, count-cnt+1, count, B1, factsfound);
      if (breadthfirst && nCandidates > 1)
	fprintf (stderr, "\r      Line=%u/%u Curves=%u/%u B1=%.0f factors=%u\r", linenum, nCandidates, breadthfirst_cnt, breadthfirst_maxcnt, B1, factsfound);
      if (verbose > 0)
	{
	  if ((!breadthfirst && cnt == count) || (breadthfirst && 1 == breadthfirst_cnt))
	    {
	      /* first time this candidate has been run (if looping more than once */
	      if (n.cpExpr && n.nexprlen < 1000)
		printf ("Input number is %s (%u digits)\n", n.cpExpr, n.ndigits);
	      else if (n.ndigits < 1000)
		{
                  char *s;
                  s = mpz_get_str (NULL, 10, n.n);
		  printf ("Input number is %s (%u digits)\n", s, n.ndigits);
                  FREE (s, n.ndigits + 1);
		}
	      else
		printf ("Input number has %u digits\n", n.ndigits);
	      if (n.isPrp)
		printf ("****** Warning: input is probably prime ******\n");
	    }
	  else
	    {
	      if (breadthfirst)
		printf ("Line=%u/%u Curves=%u/%u B1=%.0f factors=%u      \n", linenum, nCandidates, breadthfirst_cnt, breadthfirst_maxcnt, B1, factsfound);
	      else
	        printf ("Line=%u Curves=%u/%u B1=%.0f factors=%u      \n", linenum, count-cnt+1, count, B1, factsfound);
	      /* Since the expression is usally "so" short, why not just drop it out for ALL loops? */
	      if (displayexpr)
		{
		  if (n.nexprlen && n.nexprlen <= displayexpr)
		    printf ("Input number is %s (%u digits)\n", n.cpExpr, n.ndigits);
		  else if (n.ndigits <= displayexpr)
		    {
		      char *s;
		      s = mpz_get_str (NULL, 10, n.n);
		      printf ("Input number is %s (%u digits)\n", s, n.ndigits);
                      FREE (s, n.ndigits + 1);
		    }
		}
	    }
	  fflush (stdout);
	}
      /* Even in verbose==0 we should primality check if told to do so, however, we will print to stderr to keep stdout "clean"
         for verbose==0 like behavior */
      else if (((!breadthfirst && cnt == count) || (breadthfirst && breadthfirst_cnt==1)) && n.isPrp)
	{
	  char *s;
	  s = mpz_get_str (NULL, 10, n.n);
	  fprintf (stderr, "Input number is %s (%u digits)\n****** Warning: input is probably prime ******\n", s, n.ndigits);
	  FREE (s, n.ndigits + 1);
	}

      if ((!breadthfirst && cnt == count) || (breadthfirst && 1 == breadthfirst_cnt))
	{
	  int SomeFactor;
	  /*  Note, if a factors are found, then n will be adjusted "down" */
	  fprintf (stderr, "T:000 \r");
	  SomeFactor = trial_factor (&n, maxtrialdiv, deep);
	  if (SomeFactor)
	    {
	      /* should we increase factors found for trivials ??? */
	      trial_factor_found = 1;
	      factsfound += SomeFactor;
	      if (n.isPrp)
	        {
		  printf ("Probable prime cofactor ");
		  if (n.cpExpr)
		    printf ("%s", n.cpExpr);
		  else
		    mpz_out_str (stdout, 10, n.n);
		  printf (" has %u digits\n", n.ndigits);
		  /* Nothing left to do with this number, so simply continue. */
		  cnt = 0; /* no more curve to perform */
		  fflush (stdout);
		  continue;
		}
	      fflush (stdout);
	      if (!deep)
		{
		  /* Note, if we are not in deep mode, then there is no need to continue if a factor was found */
  		  factor_is_prime = 1;
  		  mpz_set_ui (f,1);
		  goto OutputFactorStuff;
		}

	    }
        }

      factor_is_prime = 0;

      cnt --; /* one more curve performed */

      if (method == PM1_METHOD)
        result = pm1 (f, x, n.n, B1done, B1, B2min, B2, k, S, verbose, repr);
      else if (method == PP1_METHOD)
        result = pp1 (f, x, n.n, B1done, B1, B2min, B2, k, S, verbose, repr);
      else /* ECM */
	{
	  if (mpz_sgn (sigma) == 0) /* If sigma is zero, then we use the A value instead */
	    result = ecm (f, x, A, n.n, B1done, B1, B2min, B2, k, S, verbose, repr, 1);
	  else
	    result = ecm (f, x, sigma, n.n, B1done, B1, B2min, B2, k, S, verbose, repr, 0);
        }

      if (result == 0)
	{
	  if (!trial_factor_found)
            fprintf (stderr, "2:100\r");
	  else
	  {
	    factor_is_prime = 1;
	    mpz_set_ui (f,1);
	    goto OutputFactorStuff;
	  }
	}
      if (result != 0)
	{
	  factsfound++;
          printf ("********** Factor found in step %u: ", ABS(result));
          mpz_out_str (stdout, 10, f);
          printf ("\n");
	  if (mpz_cmp (f, n.n))
	    {
              if (mpz_cmp_ui (f, 1) == 0)
                {
                  fprintf (stderr, "Error: factor found is 1\n");
                  /* THIS SHOULD NOT BE HERE!!! */
                  exit (1);
                }
	      /* prints factor found and cofactor on standard error. */
	      factor_is_prime = mpz_probab_prime_p (f, 25);
	      printf ("Found %s factor of %2u digits: ", 
		      factor_is_prime ? "probable prime" : "COMPOSITE     ",
		      nb_digits (f));
	      mpz_out_str (stdout, 10, f);
	      printf ("\n");

	      mpcandi_t_addfoundfactor (&n, f, 1); /* 1 for display warning if factor does not divide the current candidate */

OutputFactorStuff:;
	      if (verbose)
		{
		  printf ("%s cofactor ",
			  n.isPrp ? "Probable prime" : "Composite");
		  if (n.cpExpr)
		    printf ("%s", n.cpExpr);
		  else
		    mpz_out_str (stdout, 10, n.n);
		  printf (" has %u digits\n", n.ndigits);
		}
	      
              /* check for champions (top ten for each method) */
	      method1 = ((method == PP1_METHOD) && (result < 0))
		? PM1_METHOD : method;
	      if (factor_is_prime && nb_digits (f) >= champion_digits[method1])
                {
                  printf ("Report your potential champion to %s\n",
                          champion_keeper[method1]);
                  printf ("(see %s)\n", champion_url[method1]);
                }
	      /* Take care of fully factoring this number, in case we are in deep mode */
	      if (n.isPrp)
		  cnt = 0; /* no more curve to perform */

	      if (!deep)
	        {
		  if (breadthfirst)
		    /* I know it may not be prp, but setting this will cause all future loops to NOT 
		       check this candidate again */
		    pCandidates[linenum-1].isPrp = 1;
		  cnt = 0;
	        }
	      else if (breadthfirst)
		mpcandi_t_copy (&pCandidates[linenum-1], &n);
            }
	  else
	    {
	      if (breadthfirst)
		/* I know it may not be prp, but setting this will cause all future loops to NOT 
		   check this candidate again */
		pCandidates[linenum-1].isPrp = 1;
	      cnt = 0; /* no more curve to perform */
	      printf ("Found input number N\n");
	    }
	  fflush (stdout);
	}
      
      /* if quiet, prints composite cofactors on standard output. */
      if (!count && (verbose == 0) && !n.isPrp)
	{
	  if (n.cpExpr)
	    printf ("%s", n.cpExpr);
	  else
	    mpz_out_str (stdout, 10, n.n);
	  putchar ('\n');
	  fflush (stdout);
	}

      /* Write composite cofactors to savefile if requested */
      /* If no factor was found, we consider cofactor composite and write it */
      if (savefile != NULL && !n.isPrp)
        {
          mpz_mod (x, x, n.n); /* Reduce stage 1 residue wrt new cofactor, in
                               case a factor was found */
          write_resumefile_line (savefile, method, B1, sigma, A, x, &n, 
                                 orig_x0, comment);
        }

      /* Clean up any temp file left over.  At this time, we "assume" that if the user wants
         to resume the run, then they used -save file.  The temp save was ONLY to help in case
         of a power outage (or similar) for a long run.  It would allow finishing the current
         candidate, keeping the existing work done.   Now, we assume we "are" done. */
#if !defined (DEBUG_AUTO_SAVE)
      kill_temp_resume_file ();
#endif

      /* advance B1, if autoincrement value had been set during command line parsing */
      if (!breadthfirst && autoincrementB1)
	{
	  B1 = calc_B1_AutoIncrement(B1, autoincrementB1, autoincrementB1_calc);
	  B2min = B1;
	}
    }

  /* Allow our "breadthfirst" search to re-run the file again if enough curves have not yet been run */
  if (breadthfirst == 1)
    goto BreadthFirstDoAgain;

  /* NOTE finding a factor may have caused the loop to exit, but what is left on screen is the 
     wrong count of factors (missing the just found factor.  Update the screen to at least specify the 
     current count */
  if (breadthfirst_maxcnt)
    fprintf (stderr, "\rLine=%u Curves=%u/%u B1=%.0f factors=%u      \n", linenum, breadthfirst_cnt, breadthfirst_maxcnt, B1, factsfound);
  else if (count != 1)
    fprintf (stderr, "Line=%u Curves=%u/%u B1=%.0f factors=%u      \n", linenum, count-cnt+1, count, B1, factsfound);
  


  if (infilename)	/* note infile "might" be stdin, and don't fclose that! */
    fclose (infile);
  if (savefile)
    fclose (savefile);
  if (resumefile)
    fclose (resumefile);
  if (nCandidates)
    {
      while (nCandidates--)
	mpcandi_t_free (&pCandidates[nCandidates]);
      free (pCandidates);
    }
	  
#ifdef POLYGCD
  NTL_clear ();
#endif

  free_expr ();

  gmp_randclear (randstate);

  mpz_clear (orig_x0);
  mpz_clear (x);
  mpz_clear (f);
  mpcandi_t_free (&n);
  mpz_clear (sigma);
  mpz_clear (A);
  mpq_clear (rat_x0);

#ifdef MEMORY_DEBUG
  tests_memory_end ();
#endif

  /* exit 0 iff a factor was found for the last input */
  return result ? 0 : 1;
}
