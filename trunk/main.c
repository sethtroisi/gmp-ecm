/* GMP-ECM -- Integer factorization with ECM and Pollard 'P-1' methods.

  Copyright 2001, 2002, 2003, 2004, 2005 Jim Fougeron, Laurent Fousse, Alexander Kruppa, Paul Zimmermann.

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
  Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
  02111-1307, USA.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#ifdef _MSC_VER
#  include <winsock2.h>
#endif
#include "ecm-ecm.h"

#ifdef HAVE_UNISTD_H /* for access() */
# include <unistd.h>
#else
# define F_OK 0
# ifdef HAVE_IO_H
#  include <io.h>
# endif
#endif

#ifdef HAVE_SIGNAL_H
# include <signal.h>
#endif

#ifdef HAVE_GWNUM
/* For GWNUM_VERSION */
#include "gwnum.h"
#endif

/* Used in print_config() */
#include "ecm-params.h"

/* #define DEBUG */

#include "champions.h"

/* probab_prime_p() can get called from other modules. Instead of passing
   prpcmd to those functions, we make it static here - this variable will
   be set only in main, and read only in probab_prime_p() */
#ifdef WANT_SHELLCMD
static  char *prpcmd = NULL;
#endif

static int exit_asap_value = 0;
static int exit_asap_signalnr = 0; /* Remembers which signal we received */

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
      while (c != EOF && !IS_NEWLINE(c));
      if (IS_NEWLINE(c))
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

int
probab_prime_p (mpz_t N, int reps)
{
#ifdef WANT_SHELLCMD
  if (prpcmd != NULL)
    {
      FILE *fc;
      int r;
      fc = popen (prpcmd, "w");
      if (fc != NULL)
        {
          gmp_fprintf (fc, "%Zd\n", N);
          r = pclose (fc);
          if (r == 0) /* Exit status of 0 means success = is a PRP */
            return 1;
          else
            return 0;
        } else {
          fprintf (stderr, "Error executing the PRP command\n");
          exit (EXIT_FAILURE);
        }
    } else
#endif
      return mpz_probab_prime_p (N, reps);
}

void
signal_handler (int sig)
{
  if (sig == SIGINT || sig == SIGTERM)
    {
      exit_asap_value = 1;
      exit_asap_signalnr = sig;
      /* If one of these two signals arrives again, we'll let the default
         handler take over,  which will usually terminate the process
         immediately. */
      signal (SIGINT, SIG_DFL);
      signal (SIGTERM, SIG_DFL);
    }
  else
    {
      /* How did this happen? Let's ignore it for now, abort instead? */
    }
}

int
stop_asap_test ()
{
  return exit_asap_value;
}


static void
usage (void)
{
    printf ("Usage: ecm [options] B1 [[B2min-]B2] < file\n");
    printf ("\nParameters:\n");
    printf ("  B1           stage 1 bound\n");
    printf ("  B2           stage 2 bound (or interval B2min-B2max)\n");
    printf ("\nOptions:\n");
    printf ("  -x0 x        use x as initial point\n"); 
    printf ("  -sigma s     use s as curve generator [ecm]\n");
    printf ("  -A a         use a as curve parameter [ecm]\n");
    printf ("  -k n         perform >= n steps in stage 2\n");
    printf ("  -power n     use x^n for Brent-Suyama's extension\n");
    printf ("  -dickson n   use n-th Dickson's polynomial for Brent-Suyama's extension\n");
    printf ("  -c n         perform n runs for each input\n");
    printf ("  -pm1         perform P-1 instead of ECM\n");
    printf ("  -pp1         perform P+1 instead of ECM\n");
    printf ("  -q           quiet mode\n");
    printf ("  -v           verbose mode\n");
    printf ("  -timestamp   print a time stamp with each number\n");
    printf ("  -mpzmod      use GMP's mpz_mod for modular reduction\n");
    printf ("  -modmuln     use Montgomery's MODMULN for modular reduction\n");
    printf ("  -redc        use Montgomery's REDC for modular reduction\n");
    printf ("  -nobase2     disable special base-2 code\n");
    printf ("  -nobase2s2   disable special base-2 code in ecm stage 2 only\n");
    printf ("  -base2 n     force base 2 mode with 2^n+1 (n>0) or 2^|n|-1 (n<0)\n");
    printf ("  -ntt         enable NTT convolution routines in stage 2\n");
    printf ("  -no-ntt      disable NTT convolution routines in stage 2\n");
    printf ("  -save file   save residues at end of stage 1 to file\n");
    printf ("  -savea file  like -save, appends to existing files\n");
    printf ("  -resume file resume residues from file, reads from stdin if file is \"-\"\n");
    printf ("  -chkpnt file save periodic checkpoints during stage 1 to file\n");
    printf ("  -primetest   perform a primality test on input\n");
    printf ("  -treefile f  store product tree of F in files f.0 f.1 ... \n");
    printf ("  -maxmem n    use at most n MB of memory in stage 2\n");
    printf ("  -stage1time n add n seconds to ECM stage 1 time (for expected time est.)\n");
#ifdef WANT_SHELLCMD
    printf ("  -faccmd cmd  execute cmd when factor is found. Input number, factor\n"
            "               and cofactor are given to cmd via stdin, each on a line\n");
    printf ("  -prpcmd cmd  use shell command cmd to do prp tests (number via stdin)\n");
    printf ("  -idlecmd cmd before each curve run cmd and terminate if exit code >0\n");
#endif

    /*printf ("  -extra functions added by JimF\n"); */
    printf ("  -i n         increment B1 by this constant on each run\n");
    printf ("  -I f         auto-calculated increment for B1 multiplied by 'f' scale factor\n");
    printf ("  -inp file    Use file as input (instead of redirecting stdin)\n");
    printf ("  -b           Use breadth-first mode of file processing\n");
    printf ("  -d           Use depth-first mode of file processing (default)\n");
    printf ("  -one         Stop processing a candidate if a factor is found (looping mode)\n");
    printf ("  -n           run ecm in \"nice\" mode (below normal priority)\n");
    printf ("  -nn          run ecm in \"very nice\" mode (idle priority)\n");
    printf ("  -t n         Trial divide candidates before P-1, P+1 or ECM up to n\n");
    printf ("  -ve n        Verbosely show short (< n character) expressions on each loop\n");
    printf ("  -cofdec      Force cofactor output in decimal (even if expressions are used)\n");
    printf ("  -B2scale f   Multiplies the default B2 value by f \n");
    printf ("  -go val      Preload with group order val, which can be a simple expression,\n");
    printf ("               or can use N as a placeholder for the number being factored.\n");
    printf ("  -printconfig Print compile-time configuration and exit.\n");

    printf ("  -h, --help   Prints this help and exit.\n");
}

/* Print parameters that were used to build GMP-ECM */
static void
print_config ()
{
  printf ("Compilation options:\n");
#ifdef __GNU_MP_VERSION_PATCHLEVEL
  printf ("Included GMP header files version %d.%d.%d\n", 
          __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, 
          __GNU_MP_VERSION_PATCHLEVEL);
#else
  printf ("Included GMP header files version %d.%d\n", 
          __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR);
#endif

#ifdef GWNUM_VERSION
  printf ("Included GWNUM header files version %s\n", GWNUM_VERSION);
#else
  printf ("GWNUM_VERSION undefined\n");
#endif

#ifdef HAVE_SSE2
  printf ("HAVE_SSE2 = %d\n$", HAVE_SSE2);
#else
  printf ("HAVE_SSE2 undefined\n");
#endif

#ifdef HAVE___GMPN_ADD_NC
  printf ("HAVE___GMPN_ADD_NC = %d\n", HAVE___GMPN_ADD_NC);
#else
  printf ("HAVE___GMPN_ADD_NC undefined\n");
#endif

#ifdef HAVE___GMPN_MOD_34LSUB1
  printf ("HAVE___GMPN_MOD_34LSUB1 = %d\n", HAVE___GMPN_MOD_34LSUB1);
#else
  printf ("HAVE___GMPN_MOD_34LSUB1 undefined\n");
#endif

#ifdef HAVE___GMPN_REDC_1
  printf ("HAVE___GMPN_REDC_1 = %d\n", HAVE___GMPN_REDC_1);
#else
  printf ("HAVE___GMPN_REDC_1 undefined\n");
#endif

#ifdef MEMORY_DEBUG
  printf ("MEMORY_DEBUG = %d\n", MEMORY_DEBUG);
#else
  printf ("MEMORY_DEBUG undefined\n");
#endif

#ifdef NATIVE_REDC
  printf ("NATIVE_REDC = %d\n", NATIVE_REDC);
#ifdef TUNE_MULREDC_THRESH
  printf ("TUNE_MULREDC_THRESH = %d\n", TUNE_MULREDC_THRESH);
#else
  printf ("TUNE_MULREDC_THRESH undefined\n");
#endif
#ifdef TUNE_SQRREDC_THRESH
  printf ("TUNE_SQRREDC_THRESH = %d\n", TUNE_SQRREDC_THRESH);
#else
  printf ("TUNE_SQRREDC_THRESH undefined\n");
#endif
#ifdef WINDOWS64_ABI
  printf ("WINDOWS64_ABI = %d\n", WINDOWS64_ABI);
#else
  printf ("WINDOWS64_ABI undefined\n");
#endif
#else
  printf ("NATIVE_REDC undefined\n");
#endif

#ifdef WANT_ASSERT
  printf ("WANT_ASSERT = %d\n", WANT_ASSERT);
#else
  printf ("WANT_ASSERT undefined\n");
#endif

#ifdef WANT_SHELLCMD
  printf ("WANT_SHELLCMD = %d\n", WANT_SHELLCMD);
#else
  printf ("WANT_SHELLCMD undefined\n");
#endif

#ifdef _OPENMP
  printf ("_OPENMP = %d\n", _OPENMP);
#else
  printf ("_OPENMP undefined\n");
#endif

#ifdef MPZMOD_THRESHOLD
  printf ("MPZMOD_THRESHOLD = %d\n", MPZMOD_THRESHOLD);
#else
  printf ("MPZMOD_THRESHOLD undefined\n");
#endif

#ifdef REDC_THRESHOLD
  printf ("REDC_THRESHOLD = %d\n", REDC_THRESHOLD);
#else
  printf ("REDC_THRESHOLD undefined\n");
#endif

#ifdef MUL_NTT_THRESHOLD
  printf ("MUL_NTT_THRESHOLD = %d\n", MUL_NTT_THRESHOLD);
#else
  printf ("MUL_NTT_THRESHOLD undefined\n");
#endif

#ifdef NTT_GFP_TWIDDLE_DIF_BREAKOVER
  printf ("NTT_GFP_TWIDDLE_DIF_BREAKOVER = %d\n", 
	  NTT_GFP_TWIDDLE_DIF_BREAKOVER);
#else
  printf ("NTT_GFP_TWIDDLE_DIF_BREAKOVER undefined\n");
#endif

#ifdef NTT_GFP_TWIDDLE_DIT_BREAKOVER
  printf ("NTT_GFP_TWIDDLE_DIT_BREAKOVER = %d\n", 
	  NTT_GFP_TWIDDLE_DIT_BREAKOVER);
#else
  printf ("NTT_GFP_TWIDDLE_DIT_BREAKOVER undefined\n");
#endif

#ifdef PREREVERTDIVISION_NTT_THRESHOLD
  printf ("PREREVERTDIVISION_NTT_THRESHOLD = %d\n", 
          PREREVERTDIVISION_NTT_THRESHOLD);
#else
  printf ("PREREVERTDIVISION_NTT_THRESHOLD undefined\n");
#endif

#ifdef POLYINVERT_NTT_THRESHOLD
  printf ("POLYINVERT_NTT_THRESHOLD = %d\n", POLYINVERT_NTT_THRESHOLD);
#else
  printf ("POLYINVERT_NTT_THRESHOLD undefined\n");
#endif

#ifdef POLYEVALT_NTT_THRESHOLD
  printf ("POLYEVALT_NTT_THRESHOLD = %d\n", POLYEVALT_NTT_THRESHOLD);
#else
  printf ("POLYEVALT_NTT_THRESHOLD undefined\n");
#endif

#ifdef MPZSPV_NORMALISE_STRIDE
  printf ("MPZSPV_NORMALISE_STRIDE = %d\n", MPZSPV_NORMALISE_STRIDE);
#else
  printf ("MPZSPV_NORMALISE_STRIDE undefined\n");
#endif
}


/******************************************************************************
*                                                                             *
*                                Main program                                 *
*                                                                             *
******************************************************************************/

int
main (int argc, char *argv[])
{
  char **argv0 = argv;
  mpz_t x, sigma, A, f, orig_x0, B2, B2min, startingB2min;
  mpcandi_t n;
  mpgocandi_t go;
  mpq_t rat_x0;
  double B1, B1done;
  int result = 0, returncode = 0;
  int verbose = OUTPUT_NORMAL; /* verbose level */
  int timestamp = 0;
  int method = ECM_ECM, method1;
  int use_ntt = 1;     /* Default, use NTT if input is small enough */
  int specific_x0 = 0, /* 1=starting point supplied by user, 0=random or */
                       /* compute from sigma */
      specific_sigma = 0;  /* 1=sigma from command line, 0=make random */
  int factor_is_prime;
        /* If a factor was found, indicate whether factor, cofactor are */
        /* prime. If no factor was found, both are zero. */
  int repr = ECM_MOD_DEFAULT; /* automatic choice */
  int nobase2step2 = 0; /* flag to turn off base 2 arithmetic in ecm stage 2 */
  unsigned long k = ECM_DEFAULT_K; /* default number of blocks in stage 2 */
  int S = ECM_DEFAULT_S;
             /* Degree for Brent-Suyama extension requested by user.
                Positive value: use S-th power,
                negative: use degree |S| Dickson poly,
                default (0): automatic choice. */
  gmp_randstate_t randstate;
  char *savefilename = NULL, *resumefilename = NULL, *infilename = NULL;
  char *TreeFilename = NULL, *chkfilename = NULL;
  char rtime[256] = "", who[256] = "", comment[256] = "", program[256] = "";
  FILE *resumefile = NULL, *infile = NULL;
  mpz_t resume_lastN, resume_lastfac; /* When resuming residues from a file,
        store the last number processed and the factors found for this it */
  int resume_wasPrp = 0; /* 1 if resume_lastN/resume_lastfac is a PRP */
  int primetest = 0, saveappend = 0;
  double autoincrementB1 = 0.0, startingB1;
  unsigned int autoincrementB1_calc = 0;
  unsigned int breadthfirst_maxcnt=0, breadthfirst_cnt=0;
  int breadthfirst = 0;
  unsigned int count = 1; /* number of curves for each number */
  unsigned int cnt = 0;   /* number of remaining curves for current number */
  unsigned int linenum = 0, factsfound = 0;
  mpcandi_t *pCandidates = NULL;
  unsigned int nCandidates=0, nMaxCandidates=0;
  int deep=1, trial_factor_found;
  unsigned int displayexpr = 0;
  unsigned int decimal_cofactor = 0;
  double maxtrialdiv = 0.0;
  double B2scale = 1.0;
  double maxmem = 0.;
  double stage1time = 0.;
  ecm_params params;
#ifdef WANT_SHELLCMD
  char *faccmd = NULL;
  char *idlecmd = NULL;
#endif

  /* check ecm is linked with a compatible library */
  if (mp_bits_per_limb != GMP_NUMB_BITS)
    {
      fprintf (stderr, "Error, mp_bits_per_limb and GMP_NUMB_BITS differ\n");
      fprintf (stderr, "Please check your LD_LIBRARY_PATH variable\n");
      exit (1);
    }

#ifdef MEMORY_DEBUG
  tests_memory_start ();
#endif

  ecm_init (params);

  /* initialize the group order candidate */
  mpgocandi_t_init (&go);

  /* Init variables we might need to store options */
  MPZ_INIT (sigma);
  MPZ_INIT (A);
  MPZ_INIT (B2);
  MPZ_INIT (B2min);
  MPZ_INIT (startingB2min);
  mpq_init (rat_x0);

  /* first look for options */
  while ((argc > 1) && (argv[1][0] == '-'))
    {
      if (strcmp (argv[1], "-pm1") == 0)
	{
	  method = ECM_PM1;
	  argv++;
	  argc--;
	}
      else if (strcmp (argv[1], "-pp1") == 0)
	{
	  method = ECM_PP1;
	  argv++;
	  argc--;
	}
      else if (strcmp (argv[1], "-q") == 0)
	{
	  verbose = OUTPUT_ALWAYS;
	  argv++;
	  argc--;
	}
      else if (strcmp (argv[1], "-v") == 0)
	{
	  verbose ++;
	  argv++;
	  argc--;
	}
      else if (strcmp (argv[1], "-timestamp") == 0)
	{
	  timestamp = 1;
	  argv++;
	  argc--;
	}
      else if (strcmp (argv[1], "-mpzmod") == 0)
        {
          repr = ECM_MOD_MPZ;
	  argv++;
	  argc--;
        }
      else if (strcmp (argv[1], "-modmuln") == 0)
        {
          repr = ECM_MOD_MODMULN;
	  argv++;
	  argc--;
        }
      else if (strcmp (argv[1], "-redc") == 0)
        {
          repr = ECM_MOD_REDC;
	  argv++;
	  argc--;
        }
      else if (strcmp (argv[1], "-nobase2") == 0)
        {
          repr = ECM_MOD_NOBASE2;
	  argv++;
	  argc--;
        }
      else if (strcmp (argv[1], "-nobase2s2") == 0)
        {
          nobase2step2 = 1;
          argv++;
          argc--;
        }
      else if (strcmp (argv[1], "-ntt") == 0)
	{
	  use_ntt = 2; /* Use NTT, even for large input numbers */
	  argv++;
	  argc--;
	}
      else if (strcmp (argv[1], "-no-ntt") == 0)
	{
	  use_ntt = 0; /* Never use NTT */
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
      else if (strcmp (argv[1], "-h") == 0 || strcmp (argv[1], "--help") == 0)
        {
          usage ();
          exit (EXIT_SUCCESS);
        }
      else if (strcmp (argv[1], "-printconfig") == 0)
        {
          print_config ();
          exit (EXIT_SUCCESS);
        }
      else if (strcmp (argv[1], "-d") == 0)
        {
	  /* -1 is a flag used during argv processing where a subsquent -i file will NOT change it.  Then
	     when done processing args, we change a -1 to a 0 */
	  breadthfirst = -1;
	  argv++;
	  argc--;
        }
      else if (strcmp (argv[1], "-cofdec") == 0)
        {
	  decimal_cofactor = 1;
	  argv++;
	  argc--;
        }
      else if (strcmp (argv[1], "-n") == 0)
        {
          NICE10;
	  argv++;
	  argc--;
        }
      else if (strcmp (argv[1], "-nn") == 0)
        {
          NICE20;
	  argv++;
	  argc--;
        }
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
          if (mpz_set_str (sigma, argv[2], 0) || mpz_cmp_ui (sigma, 6) < 0)
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
	  k = atol (argv[2]);
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
	  saveappend = 0;
	  argv += 2;
	  argc -= 2;
	}
      else if ((argc > 2) && (strcmp (argv[1], "-savea") == 0))
	{
	  savefilename = argv[2];
	  saveappend = 1;
	  argv += 2;
	  argc -= 2;
	}
      else if ((argc > 2) && (strcmp (argv[1], "-resume") == 0))
	{
	  resumefilename = argv[2];
	  argv += 2;
	  argc -= 2;
	}
      else if ((argc > 2) && (strcmp (argv[1], "-chkpnt") == 0))
	{
	  chkfilename = argv[2];
	  argv += 2;
	  argc -= 2;
	}
      else if ((argc > 2) && (strcmp (argv[1], "-treefile") == 0))
	{
	  TreeFilename = argv[2];
	  argv += 2;
	  argc -= 2;
	}
      else if ((argc > 2) && (strcmp (argv[1], "-base2") == 0))
	{
	  int b = atoi (argv[2]);
	  if (abs (b) >= 16) /* |Values| < 16 are reserved for other methods */
	    repr = b;        /* keep method unchanged in that case */
	  argv += 2;
	  argc -= 2;
	}
      else if ((argc > 2) && (strcmp (argv[1], "-i") == 0))
	{
	  autoincrementB1 = strtod (argv[2], NULL);
	  if (autoincrementB1 < 1.0)
	    {
	      fprintf (stderr, "Error, the -i n option requires n >= 1\n");
	      exit (EXIT_FAILURE);
  	    }
	  argv += 2;
	  argc -= 2;
	}
      else if ((argc > 2) && (strcmp (argv[1], "-I") == 0))
	{
	  autoincrementB1 = strtod (argv[2], NULL);
	  autoincrementB1_calc = 1;
	  if (autoincrementB1 <= 0.0)
	    {
	      fprintf (stderr, "Error, the -I f option requires f > 0\n");
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
	  if (maxtrialdiv <= 0.0)
	    {
	      fprintf (stderr, "Error, the -t option requires a positive argument\n");
	      exit (EXIT_FAILURE);
  	    }
	  argv += 2;
	  argc -= 2;
	}
      else if ((argc > 2) && (strcmp (argv[1], "-ve") == 0))
        {
	  displayexpr = atoi (argv[2]);
	  if (displayexpr == 0)
	    {
	      fprintf (stderr, "Error, the -ve option requires a number argument\n");
	      exit (EXIT_FAILURE);
  	    }
	  argv += 2;
	  argc -= 2;
        }
      else if ((argc > 2) && (strcmp (argv[1], "-B2scale") == 0))
	{
	  B2scale = atof (argv[2]);
	  if (verbose >= 2)
	    printf ("Scaling B2 values by a factor of %.4f\n", B2scale);
	  argv += 2;
	  argc -= 2;
	}
      else if ((argc > 2) && (strcmp (argv[1], "-maxmem") == 0))
	{
	  maxmem = atof (argv[2]) * 1048576.;
	  argv += 2;
	  argc -= 2;
	}
      else if ((argc > 2) && (strcmp (argv[1], "-stage1time") == 0))
	{
	  stage1time = atof (argv[2]);
	  argv += 2;
	  argc -= 2;
	}
      else if ((argc > 2) && (strcmp (argv[1], "-go") == 0))
	{
	  if (go.cpOrigExpr)
            {
              fprintf (stderr, "Warning, for multiple -go options, only the last one is taken into account.\n");
              free (go.cpOrigExpr);
            }
	  go.cpOrigExpr = malloc (strlen (argv[2]) + 1);
          if (go.cpOrigExpr == NULL)
            {
              fprintf (stderr, "Cannot allocate memory in main\n");
              exit (1);
            }
	  strcpy (go.cpOrigExpr, argv[2]);
	  if (strchr (go.cpOrigExpr, 'N'))
	    {
	      go.containsN = 1;
	      go.Valid = 1;  /* we actually do not know if it is valid here,
				but we "assume" until the first time it gets
				run through */
	    }
	  else
	    { 
	      go.containsN = 0;  /* have "fully" parsed expr or number.
				    Do not recompute for each N */
	      if (eval_str (&(go.Candi), go.cpOrigExpr, 0, NULL))
		go.Valid = 1;
	    }
	  argv += 2;
	  argc -= 2;
	}
#ifdef WANT_SHELLCMD
     else if ((argc > 2) && (strcmp (argv[1], "-prpcmd") == 0))
       {
         prpcmd = argv[2];
         argv += 2;
         argc -= 2;
       }
     else if ((argc > 2) && (strcmp (argv[1], "-faccmd") == 0))
       {
         faccmd = argv[2];
         argv += 2;
         argc -= 2;
       }
     else if ((argc > 2) && (strcmp (argv[1], "-idlecmd") == 0))
       {
         idlecmd = argv[2];
         argv += 2;
         argc -= 2;
       }
#endif
      else
	{
	  fprintf (stderr, "Unknown option: %s\n", argv[1]);
	  exit (EXIT_FAILURE);
	}
    }

  /* check that S is even for old P-1 stage 2 */
  if ((method == ECM_PM1) && (S != ECM_DEFAULT_S) && (S % 2 != 0))
    {
      fprintf (stderr, "Error, S should be even for P-1\n");
      exit (EXIT_FAILURE);
    }

  /* Ok, now we can "reset" the breadthfirst switch so that we do depthfirst 
     as requested */
  if (breadthfirst == -1)
    breadthfirst = 0;

  if (argc < 2)
    {
      fprintf (stderr, "Invalid arguments. See %s --help.\n", argv0[0]);
      exit (EXIT_FAILURE);
    }

  /* start of the program */
  if (verbose >= 1)
    {
      char Gmp_version[64];
      char out0[128], *out = out0;

#ifdef __MPIR_VERSION
      sprintf (Gmp_version, "MPIR %d.%d.%d", __MPIR_VERSION,
	       __MPIR_VERSION_MINOR, __MPIR_VERSION_PATCHLEVEL);
#else /* original GMP */
      sprintf (Gmp_version, "GMP %s", gmp_version);
#endif /* __MPIR_VERSION */

      out += sprintf (out, "GMP-ECM %s [configured with %s",
                      VERSION, Gmp_version);

#ifdef HAVE_GWNUM
      out += sprintf (out, ", GWNUM %s", GWNUM_VERSION);
#endif

#ifdef NATIVE_REDC
      out += sprintf (out, ", --enable-asm-redc");
#endif

#ifdef WANT_ASSERT
      out += sprintf (out, ", --enable-assert");
#endif

      printf ("%s] [", out0);
      switch (method)
	{
	case ECM_PM1:
	  printf ("P-1");
	  break;
	case ECM_PP1:
	  printf ("P+1");
	  break;
	default:
	  printf ("ECM");
	}
      printf ("]\n");
#ifdef HAVE_GETHOSTNAME
  if (verbose >= 2)
    {
#define MNAMESIZE  64
      char mname[MNAMESIZE];
      if (gethostname (mname, MNAMESIZE) == 0)
        {
          mname[MNAMESIZE - 1] = 0; /* gethostname() may omit trailing 0 */
          printf ("Running on %s\n", mname);
        }
    }
#endif

#ifdef HAVE_GWNUM
#ifdef gwnum_is_gpl
      if (! gwnum_is_gpl())
#endif
        printf ("Due to incompatible licenses, this binary file must not "
                "be distributed.\n");
#endif
    }
  
  /* set first stage bound B1 */
  B1 = strtod (argv[1], &argv[1]);
  if (*argv[1] == '-')
    {
      B1done = B1;
      B1 = strtod (argv[1] + 1, NULL);
    }
  else
    B1done = ECM_DEFAULT_B1_DONE;
  mpz_set_si (B2min, -1); /* default, means that B2min will be set to B1 by
                             ecm(), pm1() and pp1() */

  if (B1 < 0.0 || B1done < 0.0)
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

  init_expr ();

  mpz_set_si (B2, ECM_DEFAULT_B2); /* compute it automatically from B1 */
  /* parse B2 or B2min-B2max */
  if (argc >= 3)
    {
      int c;
      double d;
      char *endptr;

      /* This is like strtok, but SunOS does not seem to have it declared in
         any header files, in spite of saying it does in the man pages... */
      for (endptr = argv[2]; *endptr != '\0' && *endptr != '-'; endptr++);
      if (*endptr == '-')
        *(endptr++) = '\0';
      else
        endptr = NULL;
      
      c = -1;
      {
	int r;

	r = gmp_sscanf (argv[2], "%Zd%n", B2, &c); /* Try parsing as integer */
	if (r <= 0)
	  {
	    /* restore original value */
	    if (endptr != NULL)
	      *(--endptr) = '-';
	    fprintf (stderr, "Invalid B2 value: %s\n", argv[2]);
	    exit (EXIT_FAILURE);
	  }
      }
#ifdef __MINGW32__
      /* MinGW scanf() returns a value 1 too high for %n */
      /* Reported to MinGW as bug number 1163607 */
      if (c > 0 && argv[2][c - 1] == 0)
        c--;
#endif

      if (c < 0 || argv[2][c] != '\0')
        {
          c = -1;
          gmp_sscanf (argv[2], "%lf%n", &d, &c); /* Try parsing scientific */
#ifdef __MINGW32__
          if (c > 0 && argv[2][c - 1] == 0)
            c--;
#endif
          mpz_set_d (B2, d);
        }
      if (c < 0 || argv[2][c] != '\0' || argv[2][0] == '\0') 
      /* If not the whole token could be parsed either way, or if there was
         no token to begin with (i.e string starting with '-') signal error */
        c = -1;
      else if (endptr != NULL) /* Did we have a '-' in there? */
        {
          mpz_set (B2min, B2);
          c = -1;
          gmp_sscanf (endptr, "%Zd%n", B2, &c);
#ifdef __MINGW32__
          if (c > 0 && endptr[c - 1] == 0)
            c--;
#endif
          if (c < 0 || endptr[c] != '\0')
            {
              gmp_sscanf (endptr, "%lf%n", &d, &c);
#ifdef __MINGW32__
              if (c > 0 && endptr[c - 1] == 0)
                c--;
#endif
              mpz_set_d (B2, d);
            }
          if (c < 0 || endptr[c] != '\0')
            c = -1;
        }
      if (c == -1)
        {
          fprintf (stderr, "Error: expected positive integer(s) B2 or "
                   "B2min-B2\n");
          exit (EXIT_FAILURE);
        }
    }

  /* set static parameters (i.e. those that don't change during the program) */
  params->verbose = verbose;
  params->method = method;
  mpz_set (params->B2, B2);
  params->k = k;
  params->S = S;
  params->repr = repr;
  params->nobase2step2 = nobase2step2;
  params->chkfilename = chkfilename;
  params->TreeFilename = TreeFilename;
  params->maxmem = maxmem;
  params->stage1time = stage1time;

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
      mpz_init (resume_lastN);
      mpz_init (resume_lastfac);
      mpz_set_ui (resume_lastfac, 1);
    }

  /* Open save file for writing, if saving is requested */
  if (savefilename != NULL)
    {
      FILE *savefile;
      /* Are we not appending and does this file already exist ? */
      if (!saveappend && access (savefilename, F_OK) == 0)
        {
          printf ("Save file %s already exists, will not overwrite\n", 
                  savefilename);
          exit (EXIT_FAILURE);
        }
      /* Test if we can open the file for writing */
      savefile = fopen (savefilename, "a");
      if (savefile == NULL)
        {
          fprintf (stderr, "Could not open file %s for writing\n", savefilename);
          exit (EXIT_FAILURE);
        }
      fclose (savefile);
    }

  if (resumefile && (specific_sigma || mpz_sgn (A) || specific_x0))
    {
      printf ("Warning: -sigma, -A and -x0 parameters are ignored when resuming from\nsave files.\n");
      mpz_set_ui (sigma, 0);
      mpz_set_ui (A, 0);
      specific_x0 = 0;
    }

  mpcandi_t_init (&n); /* number(s) to factor */
  MPZ_INIT (f); /* factor found */
  MPZ_INIT (x); /* stage 1 residue */
  MPZ_INIT (orig_x0); /* starting point, for save file */

  /* We may need random numbers for sigma/starting point */
  gmp_randinit_default (randstate);
  gmp_randseed_ui (randstate, get_random_ui ());


  /* Install signal handlers */
#ifdef HAVE_SIGNAL
  /* We catch signals only if there is a savefile. Otherwise there's nothing
     we could save by exiting cleanly, but the waiting for the code to check
     for signals may delay program end unacceptably */

  if (savefilename != NULL)
    {
      signal (SIGINT, &signal_handler);
      signal (SIGTERM, &signal_handler);
      params->stop_asap = &stop_asap_test;
    }
#endif

  /* loop for number in standard input or file */

  startingB1 = B1;
  mpz_set (startingB2min, B2min);

  if (!infilename)
    infile = stdin;

  if (breadthfirst == 1)
    {
      breadthfirst_maxcnt = count;
      count = 1;
      breadthfirst_cnt = 0;
    }

BreadthFirstDoAgain:;
  if (breadthfirst == 1)
    {
      if (breadthfirst_maxcnt > breadthfirst_cnt)
        {
	  linenum = 0;
	  if (breadthfirst_cnt++)
            {
	      double NewB1;
	      NewB1 = calc_B1_AutoIncrement (B1, autoincrementB1, autoincrementB1_calc);
	      if (mpz_cmp_d (B2min, B1) <= 0) /* floating-point equality is 
                                  unreliable, a comparison might be better */
		  mpz_set_d (B2min, NewB1);
	      B1 = NewB1;
	    }
	  else
            {
	      /* This is ONLY entered upon the first time through.  We load the entire file here so that we can loop deep, 
		  or remove a candidate if factor found, or if in deep mode and cofactor is prp (or if original candidate
		  is prp and we are prp testing) */
	      nMaxCandidates = 100;
	      pCandidates = (mpcandi_t*) malloc (nMaxCandidates *
                                                 sizeof(mpcandi_t));
              if (pCandidates == NULL)
                {
                  fprintf (stderr, "Error: not enough memory\n");
                  exit (EXIT_FAILURE);
                }

	      while (!feof (infile))
		{
		  if (read_number (&n, infile, primetest))
		    {
		      mpcandi_t_init (&pCandidates[nCandidates]);
		      mpcandi_t_copy (&pCandidates[nCandidates++], &n);
		      if (nCandidates == nMaxCandidates)
			{
			    mpcandi_t *tmp = pCandidates;
			    pCandidates = (mpcandi_t*) malloc ((nMaxCandidates
                                                 + 100) * sizeof(mpcandi_t));
                            if (pCandidates == NULL)
                              {
                                fprintf (stderr, "Error: not enough memory\n");
                                exit (EXIT_FAILURE);
                              }
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

  while (((breadthfirst && linenum < nCandidates) || cnt > 0 || 
          feof (infile) == 0) && !exit_asap_value)
    {
      trial_factor_found = 0;
      params->B1done = B1done; /* may change with resume */
      
      if (resumefile != NULL) /* resume case */
        {
	  if (count != 1)
	    {
	      fprintf (stderr, "Error, option -c and -resume are incompatible\n");
	      exit (EXIT_FAILURE);
	    }
          if (!read_resumefile_line (&method, x, &n, sigma, A, orig_x0, 
                &(params->B1done), program, who, rtime, comment, resumefile))
            break;
          
          if (mpz_cmp (n.n, resume_lastN) == 0)
            {
              /* Aha, we're trying the same number again. */
              /* We skip this attempt if: 1. the remaining cofactor after
                 the last attempt was a probable prime, or 2. if a factor
                 was found and the user gave the -one option */
              if (resume_wasPrp || 
                  (deep == 0 && mpz_cmp_ui (resume_lastfac, 1) != 0))
                  continue;
              
              /* If we found a factor in an earlier attempt, divide it out */
              if (mpz_cmp_ui (resume_lastfac, 1) > 0)
                mpcandi_t_addfoundfactor (&n, resume_lastfac, 1);
            } else {
              /* It's a different number. Set resume_lastN and 
                 resume_lastfac */
              mpz_set (resume_lastN, n.n);
              mpz_set_ui (resume_lastfac, 1);
              resume_wasPrp = n.isPrp;
            }

	  cnt = count; /* i.e. 1 */

          if (verbose >= 1)
            {
              printf ("Resuming ");
              if (method == ECM_ECM)
                printf ("ECM");
              else if (method == ECM_PM1)
                printf ("P-1");
              else if (method == ECM_PP1)
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
      else /* no-resume case */
        {
	  if (cnt) /* nothing to read: reuse old number */
	    {
              if (verbose >= OUTPUT_NORMAL)
                printf ("Run %u out of %u:\n", count - cnt + 1, count);
	    }
	  else /* new number */
	    {
	      if (!breadthfirst && !read_number (&n, infile, primetest))
		break;
	      else if (breadthfirst)
		mpcandi_t_copy (&n,&pCandidates[linenum]);
	      linenum++;
	      cnt = count;
	      /* reset B1 (and B2min) values, as they could have been advanced on the prior candidate */
	      if (!breadthfirst)
		{
	          B1 = startingB1;
		  mpz_set (B2min, startingB2min);
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

              MPZ_INIT (inv);
              mpz_invert (inv, mpq_denref (rat_x0), n.n);
              mpz_mul (inv, mpq_numref (rat_x0), inv);
              mpz_mod (x, inv, n.n);
              mpz_clear (inv);
            }
          else /* Make a random starting point for P-1 and P+1. ECM will */
               /* compute a suitable value from sigma or A if x is zero */
            {
              if (method == ECM_ECM)
                mpz_set_ui (x, 0);
              if (method == ECM_PP1)
                pp1_random_seed (x, n.n, randstate);
              if (method == ECM_PM1)
                pm1_random_seed (x, n.n, randstate);
            }
         
          if (ECM_IS_DEFAULT_B1_DONE(B1done))
            mpz_set (orig_x0, x);
          
          /* Make a random sigma if we have neither specific sigma nor A 
             given. Warning: sigma may still contain previous random value
             and thus be nonzero here even if no specific sigma was given */
          if (method == ECM_ECM && !specific_sigma && !mpz_sgn (A))
            {
              /* Make random sigma, 0 < sigma <= 2^32 */
              mpz_urandomb (sigma, randstate, 32);
              mpz_add_ui (sigma, sigma, 6); /* we need sigma >= 6 */
            }
        }
      if (verbose >= 1)
	{
	  if ((!breadthfirst && cnt == count) || (breadthfirst && 1 == breadthfirst_cnt))
	    {
	      /* first time this candidate has been run (if looping more than once */
	      if (n.cpExpr && n.nexprlen < MAX_NUMBER_PRINT_LEN)
		printf ("Input number is %s (%u digits)\n", n.cpExpr, n.ndigits);
	      else if (n.ndigits < MAX_NUMBER_PRINT_LEN)
		{
                  char *s;
                  s = mpz_get_str (NULL, 10, n.n);
		  printf ("Input number is %s (%u digits)\n", s, n.ndigits);
                  FREE (s, n.ndigits + 1);
		}
	      else
	        {
	          /* Print only first and last ten digits of the number */
	          mpz_t t, u;
	          mpz_init (t);
	          mpz_init (u);
	          mpz_ui_pow_ui (u, 5, n.ndigits - 10);
	          mpz_tdiv_q_2exp (t, n.n, n.ndigits - 10);
	          mpz_tdiv_q (t, t, u);
		  gmp_printf ("Input number is %Zd...", t);
		  mpz_ui_pow_ui (u, 10, 10);
		  mpz_tdiv_r (t, n.n, u);
		  gmp_printf ("%Zd (%u digits)\n", t, n.ndigits);
		  mpz_clear (u);
		  mpz_clear (t);
                }
              
	      if (n.isPrp)
		printf ("****** Warning: input is probably prime ******\n");
	    }
	  else /* 2nd or more try for same composite */
	    {
	      /* Since the expression is usually "so" short, why not just drop it out for ALL loops? */
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
      /* Even in verbose=0 we should primality check if told to do so, 
         however, we will print to stderr to keep stdout "clean"
         for verbose=0 like behavior */
      else if (((!breadthfirst && cnt == count) || 
                (breadthfirst && breadthfirst_cnt==1)) && n.isPrp)
	{
	  char *s;
	  s = mpz_get_str (NULL, 10, n.n);
	  fprintf (stderr, "Input number is %s (%u digits)\n"
	           "****** Warning: input is probably prime ******\n", 
	           s, n.ndigits);
	  FREE (s, n.ndigits + 1);
	}

      if ((!breadthfirst && cnt == count) || 
          (breadthfirst && 1 == breadthfirst_cnt))
	{
	  int SomeFactor;
	  /*  Note, if a factors are found, then n will be adjusted "down" */
	  SomeFactor = trial_factor (&n, maxtrialdiv, deep);
	  if (SomeFactor)
	    {
	      /* should we increase factors found for trivials ??? */
	      trial_factor_found = 1;
	      factsfound += SomeFactor;
	      if (n.isPrp)
	        {
		  printf ("Probable prime cofactor ");
		  if (n.cpExpr && !decimal_cofactor)
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
		  /* Note, if we are not in deep mode, then there is no need 
		     to continue if a factor was found */
  		  factor_is_prime = 1;
  		  mpz_set_ui (f,1);
		  goto OutputFactorStuff;
		}

	    }
        }

      factor_is_prime = 0;

      cnt --; /* one more curve performed */

      mpgocandi_fixup_with_N (&go, &n);

      /* set parameters that may change from one curve to another */
      params->method = method; /* may change with resume */
      mpz_set (params->x, x); /* may change with resume */
      /* if sigma is zero, then we use the A value instead */
      params->sigma_is_A = (mpz_sgn (sigma) == 0) ? 1 : 0;
      mpz_set (params->sigma, (params->sigma_is_A) ? A : sigma);
      mpz_set (params->go, go.Candi.n); /* may change if contains N */
      mpz_set (params->B2min, B2min); /* may change with -c */
      /* Here's an ugly hack to pass B2scale to the library somehow.
         It gets piggy-backed onto B1done */
      params->B1done = params->B1done + floor (B2scale * 128.) / 134217728.; 
      /* Default, for P-1/P+1 with old stage 2 and ECM, use NTT only 
         for small input */
      if (use_ntt == 1 && (method == ECM_ECM || S != ECM_DEFAULT_S)) 
        params->use_ntt = (mpz_size (n.n) <= NTT_SIZE_THRESHOLD);
      else 
        params->use_ntt = (use_ntt != 0);

#ifdef WANT_SHELLCMD
      /* See if the system is currently idle, if -idlecmd was given */
      if (idlecmd != NULL)
        {
          int r;
          FILE *fc;
          fc = popen (idlecmd, "r");
          if (fc == NULL)
            {
              fprintf (stderr, "Error executing idle command: %s\n",
                       idlecmd);
              exit (EXIT_FAILURE);
            }
          r = pclose (fc);
          if (r != 0) /* If exit status of idle command is non-zero */
            {
              printf ("Idle command returned %d, exiting\n", r);
              breadthfirst = 0; /* Avoid looping due to goto (ugly, FIXME!) */
              break;
            }
        }
#endif /* WANT_SHELLCMD */

      if (timestamp)
        {
          time_t t;
          
          t = time (NULL);
          printf ("[%.24s]\n", ctime (&t));
        }

#if 0
      /* Test mpres_muldivbysomething_si() which is not called in normal
         operation */
      mpmod_selftest (n.n);
#endif
      
      if (mpz_cmp_ui (n.n, 0) <= 0)
	{
	  fprintf (stderr, "Error, input number should be positive\n");
	  exit (EXIT_FAILURE);
	}

      /* now call the ecm library */
      result = ecm_factor (f, n.n, B1, params);

      if (result == ECM_ERROR)
        {
          fprintf (stderr, "Please report internal errors at <%s>.\n",
                   PACKAGE_BUGREPORT);
          exit (EXIT_FAILURE);
        }

      if (result == ECM_NO_FACTOR_FOUND)
	{
	  if (trial_factor_found)
	  {
	    factor_is_prime = 1;
	    mpz_set_ui (f, 1);
	    returncode = ECM_NO_FACTOR_FOUND;
	    goto OutputFactorStuff;
	  }
	} else {
	  factsfound++;
	  if (verbose > 0)
	    printf ("********** Factor found in step %u: ", ABS (result));
          mpz_out_str (stdout, 10, f);
	  if (verbose > 0)
            printf ("\n");

          /* Complain about non-proper factors (0, negative) */
          if (mpz_cmp_ui (f, 1) < 0)
            {
              fprintf (stderr, "Error: factor found is ");
              mpz_out_str (stderr, 10, f);
              fprintf (stderr, "\nPlease report internal errors at <%s>.\n",
                       PACKAGE_BUGREPORT);
              exit (EXIT_FAILURE);
            }
          
#ifdef WANT_SHELLCMD
	  if (faccmd != NULL)
	    {
	      FILE *fc;
	      fc = popen (faccmd, "w");
	      if (fc != NULL)
	        {
	          mpz_t cof;
	          mpz_init_set (cof, n.n);
	          mpz_divexact (cof, cof, f);
	          gmp_fprintf (fc, "%Zd\n", n.n);
	          gmp_fprintf (fc, "%Zd\n", f);
	          gmp_fprintf (fc, "%Zd\n", cof);
	          mpz_clear (cof);
	          pclose (fc);
	        }
	    }
#endif

	  if (mpz_cmp (f, n.n) != 0)
	    {
	      /* prints factor found and cofactor on standard output. */
	      factor_is_prime = probab_prime_p (f, PROBAB_PRIME_TESTS);

              if (verbose >= 1)
                {
                  printf ("Found %s factor of %2u digits: ", 
                          factor_is_prime ? "probable prime" : "composite",
                          nb_digits (f));
                  mpz_out_str (stdout, 10, f);
                  printf ("\n");
                }

	      mpcandi_t_addfoundfactor (&n, f, 1); /* 1 for display warning if factor does not divide the current candidate */

              if (resumefile != NULL)
                {
                  /* If we are resuming from a save file, add factor to the
                     discovered factors for the current number */
                  mpz_mul (resume_lastfac, resume_lastfac, f);
                  resume_wasPrp = n.isPrp;
                }

              if (factor_is_prime)
                returncode = (n.isPrp) ? ECM_PRIME_FAC_PRIME_COFAC : 
		                         ECM_PRIME_FAC_COMP_COFAC;
              else
                returncode = (n.isPrp) ? ECM_COMP_FAC_PRIME_COFAC :
		                         ECM_COMP_FAC_COMP_COFAC;

OutputFactorStuff:;
	      if (verbose >= 1)
		{
		  printf ("%s cofactor ",
			  n.isPrp ? "Probable prime" : "Composite");
		  if (n.cpExpr && !decimal_cofactor)
		    printf ("%s", n.cpExpr);
		  else
		    mpz_out_str (stdout, 10, n.n);
		  printf (" has %u digits\n", n.ndigits);
		}
              else /* quiet mode: just print a space here, remaining cofactor
                      will be printed after last curve */
                printf (" ");
	      
              /* check for champions (top ten for each method) */
	      method1 = ((method == ECM_PP1) && (result < 0))
		? ECM_PM1 : method;
	      if ((verbose > 0) && factor_is_prime && 
                  nb_digits (f) >= champion_digits[method1])
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
		/* I know it may not be prp, but setting this will cause all 
		   future loops to NOT check this candidate again */
		pCandidates[linenum-1].isPrp = 1;
	      cnt = 0; /* no more curve to perform */
              if (verbose > 0)
                printf ("Found input number N");
              printf ("\n");
              returncode = ECM_INPUT_NUMBER_FOUND;
	    }
	  fflush (stdout);
	}

      /* if quiet mode, prints remaining cofactor after last curve */
      if ((cnt == 0) && (verbose == 0))
	{
	  if (n.cpExpr && !decimal_cofactor)
	    printf ("%s", n.cpExpr);
	  else
	    mpz_out_str (stdout, 10, n.n);
	  putchar ('\n');
	  fflush (stdout);
	}

      /* Write composite cofactors to savefile if requested */
      /* If no factor was found, we consider cofactor composite and write it */
      if (savefilename != NULL && !n.isPrp)
        {
          mpz_mod (x, params->x, n.n); /* Reduce stage 1 residue wrt new co-
                                          factor, in case a factor was found */
          /* We write the B1done value to the safe file. This requires that
             a correct B1done is returned by the factoring functions */
          write_resumefile_line (savefilename, method, params->B1done, sigma, 
                                 A, x, &n, orig_x0, comment);
        }

      /* advance B1, if autoincrement value had been set during command line parsing */
      if (!breadthfirst && autoincrementB1 > 0.0)
	{
	  double NewB1;
	  NewB1 = calc_B1_AutoIncrement (B1, autoincrementB1, autoincrementB1_calc);
	  if (mpz_cmp_d (B2min, B1) <= 0) /* <= might be better than == */
	    mpz_set_d (B2min, NewB1);
	  B1 = NewB1;
	}
    }

  /* Allow our "breadthfirst" search to re-run the file again if enough curves have not yet been run */
  if (breadthfirst == 1 && !exit_asap_value)
    goto BreadthFirstDoAgain;

  /* NOTE finding a factor may have caused the loop to exit, but what is left 
     on screen is the wrong count of factors (missing the just found factor.  
     Update the screen to at least specify the current count */

  if (infilename)	/* infile might be stdin, don't fclose that! */
    fclose (infile);
  if (resumefile)
    {
      fclose (resumefile);
      mpz_clear (resume_lastN);
      mpz_clear (resume_lastfac);
    }
  if (nCandidates)
    {
      while (nCandidates--)
	mpcandi_t_free (&pCandidates[nCandidates]);
      free (pCandidates);
    }
	  
  free_expr ();

  gmp_randclear (randstate);

  mpz_clear (orig_x0);
  mpz_clear (startingB2min);
  mpz_clear (B2min);
  mpz_clear (B2);
  mpz_clear (x);
  mpz_clear (f);
  mpcandi_t_free (&n);
  mpz_clear (sigma);
  mpz_clear (A);
  mpq_clear (rat_x0);
  mpgocandi_t_free (&go);

  ecm_clear (params);

#ifdef MEMORY_DEBUG
  tests_memory_end ();
#endif

  /* exit 0 if a factor was found for the last input, except if we exit due
     to a signal */
#ifdef HAVE_SIGNAL
  if (returncode == 0 && exit_asap_value != 0)
    returncode = 143;
#endif

  return returncode;
}
