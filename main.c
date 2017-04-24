/* GMP-ECM -- Integer factorization with ECM, P-1 and P+1 methods.

Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011,
2012, 2016 Jim Fougeron, Laurent Fousse, Alexander Kruppa, Paul Zimmermann,
Cyril Bouvier, David Cleaver

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
more details.

You should have received a copy of the GNU General Public License
along with this program; see the file COPYING.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#ifdef _MSC_VER
#  include <winsock2.h>
#endif
#include "ecm-impl.h"
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

#ifdef HAVE_TORSION
#include "torsions.h" /* to benefit from more torsion groups */
#endif

/* #define DEBUG */

static int exit_asap_value = 0;
static int exit_asap_signalnr = 0; /* Remembers which signal we received */

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
    printf ("  -y0 y        use y as initial point (Weierstrass form)\n"); 
    printf ("  -param i     which parametrization should be used [ecm]\n");
    printf ("  -sigma s     use s as parameter to compute curve's coefficients"
          "\n               can use -sigma i:s to specify -param i at the same"
                                                              " time [ecm]\n");
    printf ("  -A A         use A as a curve coefficient [ecm, see README]\n");
    printf ("  -torsion T   to generate a curve with torsion group T "
	                                                "[ecm, see README]\n");
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
    printf ("  -chkpnt file save periodic checkpoints during stage 1 to file (for -param 0)\n");
    printf ("  -primetest   perform a primality test on input\n");
    printf ("  -treefile f  [ECM only] store stage 2 data in files f.0, ... \n");
    printf ("  -maxmem n    use at most n MB of memory in stage 2\n");
    printf ("  -stage1time n add n seconds to ECM stage 1 time (for expected time est.)\n");

    printf ("  -I f         increment B1 by f*sqrt(B1) on each run\n");
    printf ("  -inp file    Use file as input (instead of redirecting stdin)\n");
    printf ("  -one         Stop processing a candidate if a factor is found (looping mode)\n");
    printf ("  -go val      Preload with group order val, which can be a simple expression,\n");
    printf ("               or can use N as a placeholder for the number being factored.\n");
    printf ("  -printconfig Print compile-time configuration and exit.\n");

    printf ("  -bsaves file With -param 1-3, save stage 1 exponent in file.\n");
    printf ("  -bloads file With -param 1-3, load stage 1 exponent from file.\n");
#ifdef WITH_GPU
    printf ("  -gpu         Use GPU-ECM for stage 1.\n");
    printf ("  -gpudevice n Use device n to execute GPU code (by default, "
                                                          "CUDA chooses)\n");
    printf ("  -gpucurves n Compute on n curves in parallel on the GPU (by "
                                                  "default, CUDA chooses)\n");
#endif
    printf ("  -h, --help   Prints this help and exit.\n");
}

/* Print parameters that were used to build GMP-ECM */
static void
print_config ()
{
  printf ("Compilation options:\n");
#ifdef __MPIR_VERSION
     printf ("Included MPIR header files version %d.%d.%d\n", 
             __MPIR_VERSION, __MPIR_VERSION_MINOR, __MPIR_VERSION_PATCHLEVEL);
#else /* __MPIR_VERSION */
   #ifdef __GNU_MP_VERSION_PATCHLEVEL
     printf ("Included GMP header files version %d.%d.%d\n", 
             __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, 
             __GNU_MP_VERSION_PATCHLEVEL);
   #else
     printf ("Included GMP header files version %d.%d\n", 
             __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR);
   #endif
#endif /* __MPIR_VERSION */
  printf ("GMP_NUMB_BITS = %u\n", GMP_NUMB_BITS);

#ifdef ECM_TUNE_CASE
  printf ("Tuning parameters from %s\n", ECM_TUNE_CASE);
#else
  printf ("ECM_TUNE_CASE not defined.\n");
#endif

#ifdef GWNUM_VERSION
  printf ("Included GWNUM header files version %s\n", GWNUM_VERSION);
#else
  printf ("GWNUM_VERSION undefined\n");
#endif

#ifdef HAVE_SSE2
  printf ("HAVE_SSE2 = %d\n", HAVE_SSE2);
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

#ifdef USE_ASM_REDC
  printf ("USE_ASM_REDC = %d\n", USE_ASM_REDC);
#ifdef WINDOWS64_ABI
  printf ("WINDOWS64_ABI = %d\n", WINDOWS64_ABI);
#else
  printf ("WINDOWS64_ABI undefined\n");
#endif
#else
  printf ("USE_ASM_REDC undefined\n");
#endif

#ifdef WANT_ASSERT
  printf ("WANT_ASSERT = %d\n", WANT_ASSERT);
#else
  printf ("WANT_ASSERT undefined\n");
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

#ifdef WITH_GPU
  printf ("WITH_GPU = %d\n", WITH_GPU);
#else
  printf ("WITH_GPU undefined\n");
#endif

}

/* r <- q mod N. 
   Return value: 1 if den invertible, 0 if factor found; in this case
   gcd(den(q), N) is put in r.
 */
static int
mod_from_mpq (mpz_t r, mpq_t q, mpz_t N, int verbose)
{
    mpz_t inv, C;
    int factor_is_prime, cofactor_is_prime, ret = ECM_NO_FACTOR_FOUND;
 
    mpz_init (inv);
    if (mpz_invert (inv, mpq_denref (q), N) == 0)
      {
	mpz_gcd (r, mpq_denref (q), N);
	if (verbose > 0)
          {
            if (verbose > 1)
              printf ("Warning: factor found during initialization\n");
            printf ("********** Factor found in step 1: ");
          }
	mpz_out_str (stdout, 10, r);
	if (verbose > 0)
	    printf ("\n");
	if (mpz_cmp (r, N) == 0)
	  ret = ECM_INPUT_NUMBER_FOUND;
	else
	  {
	    factor_is_prime = mpz_probab_prime_p (r, PROBAB_PRIME_TESTS);
	    mpz_init (C);
	    mpz_tdiv_q (C, N, r);
	    cofactor_is_prime = mpz_probab_prime_p (C, PROBAB_PRIME_TESTS);
	    mpz_clear (C);
	    if (factor_is_prime)
		ret = cofactor_is_prime ? ECM_PRIME_FAC_PRIME_COFAC :
		    ECM_PRIME_FAC_COMP_COFAC;
	    else
		ret = cofactor_is_prime ? ECM_COMP_FAC_PRIME_COFAC :
		    ECM_COMP_FAC_COMP_COFAC;
	  }
      }
    else
      {
	mpz_mul (inv, mpq_numref (q), inv);
	mpz_mod (r, inv, N);
      }
    mpz_clear (inv);
    return ret;
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
  mpz_t x, y, sigma, A, f, orig_x0, orig_y0, B2, B2min, startingB2min, tmp_n;
  mpcandi_t n;
  mpgocandi_t go;
  mpq_t rat_x0, rat_y0, rat_A;
  mpz_t numer_A, denom_A; /* used in Hessian stuff */
  double B1, B1done;
  int result, returncode = 0;
  int verbose = OUTPUT_NORMAL; /* verbose level */
  int timestamp = 0;
  int method = ECM_ECM;
  int use_ntt = 1;     /* Default, use NTT if input is small enough */
  int specific_x0 = 0, /* 1=starting point supplied by user, 0=random or */
                       /* compute from sigma */
      specific_y0 = 0, /* was y0 given by the user? */
      specific_A = 0,  /* one may want its own A, including A=0 */
      specific_sigma = 0;  /*   0=make random */
                           /*   1=sigma from command line */
  int repr = ECM_MOD_DEFAULT; /* automatic choice */
  int nobase2step2 = 0; /* flag to turn off base 2 arithmetic in ecm stage 2 */
  unsigned long k = ECM_DEFAULT_K; /* default number of blocks in stage 2 */
  int S = ECM_DEFAULT_S;
             /* Degree for Brent-Suyama extension requested by user.
                Positive value: use S-th power,
                negative: use degree |S| Dickson poly,
                default (0): automatic choice. */
  char *savefilename = NULL, *resumefilename = NULL, *infilename = NULL;
  char *TreeFilename = NULL, *chkfilename = NULL;
#ifdef HAVE_TORSION
  char *torsion = NULL;
#endif
  char rtime[256] = "", who[256] = "", comment[256] = "", program[256] = "";
  FILE *resumefile = NULL, *infile = NULL;
  mpz_t resume_lastN, resume_lastfac; /* When resuming residues from a file,
        store the last number processed and the factors found for this it */
  int resume_wasPrp = 0; /* 1 if resume_lastN/resume_lastfac is a PRP */
  int primetest = 0, saveappend = 0;
  double autoincrementB1 = 0.0, startingB1;
  unsigned int count = 1; /* number of curves for each number */
  unsigned int cnt = 0;   /* number of remaining curves for current number */
  int deep=1;
  double maxmem = 0.;
  double stage1time = 0.;
  ecm_params params;
  int param = ECM_PARAM_DEFAULT; /* automatic choice */
  char *savefile_s = NULL;
  char *loadfile_s = NULL;
#ifdef HAVE_GWNUM
  double gw_k = 0.0;       /* set default values for gwnum poly k*b^n+c */
  unsigned long gw_b = 0;  /* set default values for gwnum poly k*b^n+c */
  unsigned long gw_n = 0;  /* set default values for gwnum poly k*b^n+c */
  signed long gw_c = 0;    /* set default values for gwnum poly k*b^n+c */
#endif
  int use_gpu = 0; /* Do we use the GPU for stage 1 (by default no)*/
  int gpudevice = -1; /* Which device do we use for GPU code (by default CUDA */
                      /* chooses)                                             */
  unsigned int gpucurves = 0; /* How many curves do we want for GPU code */ 
                              /* (by default CUDA chooses)               */

  /* check ecm is linked with a compatible library */
  if (mp_bits_per_limb != GMP_NUMB_BITS)
    {
      fprintf (stderr, "Error, mp_bits_per_limb and GMP_NUMB_BITS differ\n");
      fprintf (stderr, "Please check your LD_LIBRARY_PATH variable\n");
      exit (1);
    }

  ecm_init (params);

  /* initialize the group order candidate */
  mpgocandi_t_init (&go);

  /* Init variables we might need to store options */
  mpz_init (sigma);
  mpz_init (A);
  mpz_init (B2);
  mpz_init (B2min);
  mpz_init (startingB2min);
  mpq_init (rat_A);
  mpq_init (rat_x0);
  mpq_init (rat_y0);

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
      else if ((argc > 2) && (strcmp (argv[1], "-bsaves") == 0))
        {
          savefile_s = argv[2];
          argv += 2;
          argc -= 2;
        }
      else if ((argc > 2) && (strcmp (argv[1], "-bloads") == 0))
        {
          loadfile_s = argv[2];
          argv += 2;
          argc -= 2;
        }
      else if (strcmp (argv[1], "-h") == 0 || strcmp (argv[1], "--help") == 0)
        {
          usage ();
          goto free_all;
        }
      else if (strcmp (argv[1], "-printconfig") == 0)
        {
          print_config ();
          goto free_all;
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
      else if ((argc > 2) && (strcmp (argv[1], "-y0")) == 0)
        {
          if (mpq_set_str (rat_y0, argv[2], 0))
            {
              fprintf (stderr, "Error, invalid starting point: %s\n", argv[2]);
              exit (EXIT_FAILURE);
            }
          specific_y0 = 1;
	  argv += 2;
	  argc -= 2;
        }
      else if ((argc > 2) && (strcmp (argv[1], "-param")) == 0)
        {
          /* atoi does not allow use to distinct from */ 
          /* '-param 0' and '-param -otherargs' as it will always return 0 */
          if (argv[2][0] == '-') 
            {
              fprintf (stderr, "Error, invalid -param value: %s\n", argv[2]);
              exit (EXIT_FAILURE);
            }
          /* If param was already set (by -sigma i:x), we should check that */
          /* the same value is passed with -param */ 
          if (param != ECM_PARAM_DEFAULT)
            {
              if (param != atoi (argv[2]))
                {
                  fprintf (stderr, "Error, conflict between -sigma and -param "
                                   "arguments\n");
                  exit (EXIT_FAILURE);
                }
            }
          else
              param = atoi (argv[2]);
          argv += 2;
          argc -= 2;
        }
      else if ((argc > 2) && (strcmp (argv[1], "-sigma")) == 0)
        {
          /* Accept stuff like '-sigma i:s', this is equivalent to */
          /* '-param i -sigma s'. If we have only -sigma s and param is not */ 
          /* yet defined we assumed this is 0 (for compatibility reason) */
          if (argv[2][1] == ':')
            {
              argv[2][1] = '\0';
              /* If param was already set (by -param i:x), we should check */
              /* that the same value is passed with -sigma */ 
              if (param != ECM_PARAM_DEFAULT)
                {
                  if (param != atoi (argv[2]))
                    {
                      fprintf (stderr, "Error, conflict between -sigma and "
                                       "-param arguments\n");
                      exit (EXIT_FAILURE);
                    }
                }
              else
                  param = atoi (argv[2]) ;

              if (mpz_set_str (sigma, argv[2]+2, 0) || mpz_sgn (sigma) == 0) 
                {
                  fprintf (stderr, "Error, invalid sigma value: %s\n", 
                           argv[2]+2);
                  exit (EXIT_FAILURE);
                }
            }
          else
            {
              if (param == ECM_PARAM_DEFAULT)
                param = ECM_PARAM_SUYAMA; 
              if (mpz_set_str (sigma, argv[2], 0)) 
                {
                  fprintf (stderr, "Error, invalid sigma value: %s\n", argv[2]);
                  exit (EXIT_FAILURE);
                }
            }
          specific_sigma = 1;
          argv += 2;
          argc -= 2;
        }
      else if ((argc > 2) && (strcmp (argv[1], "-A")) == 0)
        {
          if (mpq_set_str (rat_A, argv[2], 0))
	    {
	      fprintf (stderr, "Error, invalid A value: %s\n", argv[2]);
              exit (EXIT_FAILURE);
	    }
	  specific_A = 1;
	  argv += 2;
	  argc -= 2;
        }
#ifdef HAVE_TORSION
      else if ((argc > 2) && (strcmp (argv[1], "-torsion")) == 0)
        {
	  torsion = argv[2];
	  argv += 2;
	  argc -= 2;
        }
#endif
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
      else if ((argc > 2) && (strcmp (argv[1], "-I") == 0))
	{
	  autoincrementB1 = strtod (argv[2], NULL);
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
	  if (!infile)
	    {
	      fprintf (stderr, "Can't find input file %s\n", infilename);
	      exit (EXIT_FAILURE);
	    }
	  argv += 2;
	  argc -= 2;
	}
      else if ((argc > 2) && (strcmp (argv[1], "-maxmem") == 0))
	{
	  maxmem = atof (argv[2]) * 1048576.;
          ASSERT_ALWAYS (maxmem >= 0.0);
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
#ifdef WITH_GPU
      else if (strcmp (argv[1], "-gpu") == 0)
        {
          use_gpu = 1;
	        argv++;
	        argc--;
        }
      else if ((argc > 2) && (strcmp (argv[1], "-gpudevice") == 0))
        {
          gpudevice = atoi (argv[2]);
          argv += 2;
          argc -= 2;
        }
      else if ((argc > 2) && (strcmp (argv[1], "-gpucurves") == 0))
        {
          gpucurves = atoi (argv[2]);
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

  if (argc < 2)
    {
      fprintf (stderr, "Invalid arguments. See %s --help.\n", argv0[0]);
      exit (EXIT_FAILURE);
    }

  if (strcmp (ecm_version (), ECM_VERSION) != 0)
    {
      fprintf (stderr, "Error, library version %s differs from header "
               "version %s with which this file was compiled\n",
               ecm_version (), ECM_VERSION);
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
                      ecm_version (), Gmp_version);

#ifdef HAVE_GWNUM
      out += sprintf (out, ", GWNUM %s", GWNUM_VERSION);
#endif

#ifdef USE_ASM_REDC
      out += sprintf (out, ", --enable-asm-redc");
#endif

#ifdef WITH_GPU
      out += sprintf (out, ", --enable-gpu");
#endif

#ifdef WANT_ASSERT
      out += sprintf (out, ", --enable-assert");
#endif

#ifdef _OPENMP
      out += sprintf (out, ", --enable-openmp");
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
#ifdef ECM_TUNE_CASE
      if (verbose >= 2)
        printf ("Tuned for %s\n", ECM_TUNE_CASE);
#endif
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
      if (! gwnum_is_gpl ())
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
          /* make sure B2min is not less than B1 */
          if (mpz_cmp_d (B2min, B1) < 0)
            mpz_set_d (B2min, B1);

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
  mpz_set (params->B2min, B2min);
  mpz_set (params->B2, B2);
  params->k = k;
  params->S = S;
  params->repr = repr;
  params->nobase2step2 = nobase2step2;
  params->chkfilename = chkfilename;
  params->TreeFilename = TreeFilename;
  params->maxmem = maxmem;
  params->stage1time = stage1time;
  params->gpu = use_gpu; /* If WITH_GPU is not defined it will always be 0 */
  params->gpu_device = gpudevice; /* If WITH_GPU is not defined or    */
                                  /* use_gpu = 0, it has no meaning   */
  params->gpu_number_of_curves = gpucurves; /* If WITH_GPU is not defined or */
                                            /* use_gpu = 0, it has no meaning*/

  /* Open resume file for reading, if resuming is requested */
  if (resumefilename != NULL)
    {
      /* -resume should not be used with -gpu */
      if (use_gpu)
        {
          fprintf (stderr, "Error, -resume not allowed with -gpu\n");
          exit (EXIT_FAILURE);
        }
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

  if (specific_sigma && (specific_x0 || specific_A))
    {
      fprintf (stderr, "Error, -sigma parameter is incompatible with "
                       "-A and -x0 parameters.\n");
      exit (EXIT_FAILURE);
    }

  if (specific_y0 && (!specific_x0 || !specific_A))
    {
      fprintf (stderr, "Error, -y0 must be used with -A and -x0 parameters.\n");
      exit (EXIT_FAILURE);
    }

  if (resumefile && (specific_sigma || param != ECM_PARAM_DEFAULT || 
                     specific_A || specific_x0))
    {
      printf ("Warning: -sigma, -param, -A and -x0 parameters are\n" 
              "ignored when resuming from save files.\n");
      mpz_set_ui (sigma, 0);
      specific_sigma = 0;
      param = ECM_PARAM_DEFAULT; 
      mpz_set_ui (A, 0);
      specific_x0 = 0;
    }

  mpcandi_t_init (&n); /* number(s) to factor */
  mpz_init (f); /* factor found */
  mpz_init (x); /* stage 1 residue */
  mpz_init (y); /* stage 1 for ECM_W or ECM_H */
  mpz_init (orig_x0); /* starting point, for save file */
  mpz_init (orig_y0); /* starting point, for save file */

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

  /* Main loop */
  while ((cnt > 0 || feof (infile) == 0) && !exit_asap_value)
    {
      result = ECM_NO_FACTOR_FOUND;
      params->B1done = B1done; /* may change with resume */
      
      if (resumefile != NULL) /* resume case */
        {
          if (count != 1)
            {
              fprintf (stderr, 
                       "Error, option -c and -resume are incompatible\n");
              exit (EXIT_FAILURE);
            }
          if (!read_resumefile_line (&method, x, y, &n, sigma, A, 
				     orig_x0, orig_y0, &(params->E->type), 
				     &(params->param), &(params->B1done), 
				     program, who, rtime, comment, resumefile))
            break;

	  if (params->E->type == ECM_EC_TYPE_WEIERSTRASS
	      || params->E->type == ECM_EC_TYPE_HESSIAN
	      || params->E->type == ECM_EC_TYPE_TWISTED_HESSIAN)
	      params->sigma_is_A = -1;
	  else
            params->sigma_is_A = mpz_sgn (sigma) == 0; /* sure? */

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
            } 
          else 
            {
              /* It's a different number. Set resume_lastN and 
                 resume_lastfac */
              mpz_set (resume_lastN, n.n);
              mpz_set_ui (resume_lastfac, 1);
              resume_wasPrp = n.isPrp;
            }
  
          cnt = count; /* i.e. 1 */

          if (verbose >= OUTPUT_NORMAL)
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
              if (!read_number (&n, infile, primetest))
                break;

              cnt = count;
              /* reset B1 (and B2min) values, as they could have been advanced 
                 on the prior candidate */
              B1 = startingB1;
              mpz_set (B2min, startingB2min);
            }

            /* "blank" lines, or lines that could not be parsed correctly 
               will leave a 1 in this value */
            if (n.isPrp)
              {
                /* n is 0 or 1 (or -1 I guess) so do NOT proceed with it */
                cnt = 0;
                continue;
              }

          /* Set effective seed for factoring attempt on this number */
	  if (specific_A)
	    {
		if (param != ECM_PARAM_TWISTED_HESSIAN)
		    {
			returncode = mod_from_mpq (A, rat_A, n.n, verbose);
		    }
		else
		    {
			mpz_init(numer_A);
			mpz_init(denom_A);
			mpz_mod(numer_A, mpq_numref(rat_A), n.n);
			mpz_mod(denom_A, mpq_denref(rat_A), n.n);
			returncode = ECM_NO_FACTOR_FOUND;
		    }
		if (returncode != ECM_NO_FACTOR_FOUND)
                  goto free_all1;
	    }
	  
          if (specific_x0) /* convert rational value to integer */
	    {
	      if (count != 1)
		{
		  fprintf (stderr, "Error, option -c is incompatible with -x0\n");
		  exit (EXIT_FAILURE);
                }

	      returncode = mod_from_mpq (x, rat_x0, n.n, verbose);
	      if (returncode != ECM_NO_FACTOR_FOUND)
                  goto free_all1;

	      if (specific_y0)
		{
		  returncode = mod_from_mpq (y, rat_y0, n.n, verbose);
		  if (returncode != ECM_NO_FACTOR_FOUND)
                  goto free_all1;
		}
            }
          else /* Set x to 0. This will tell the ecm library to choose a
                  random starting point for P-1 and P+1. ECM will
                  compute a suitable value from sigma or A if x is zero. */
            mpz_set_ui (x, 0);
         
          if (ECM_IS_DEFAULT_B1_DONE (B1done)) /* first time */
	    {
	      mpz_set (orig_x0, x);
	      if (specific_y0)
		mpz_set (orig_y0, y);  
	    }
        }
      if (verbose >= OUTPUT_NORMAL)
        {
          if (cnt == count)
            {
            /* first time this candidate has been run 
               (if looping more than once) */
              if (n.cpExpr && n.nexprlen < MAX_NUMBER_PRINT_LEN)
                  printf ("Input number is %s (%u digits)\n", n.cpExpr, 
                                                                  n.ndigits);
              else if (n.ndigits < MAX_NUMBER_PRINT_LEN)
                  gmp_printf ("Input number is %Zd (%u digits)\n", 
                              n.n, n.ndigits);
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
          fflush (stdout);
        }
      /* Even in verbose=0 we should primality check if told to do so, 
         however, we will print to stderr to keep stdout "clean"
         for verbose=0 like behavior */
      else if (cnt == count && n.isPrp)
          gmp_fprintf (stderr, "Input number is %Zd (%u digits)\n"
                           "****** Warning: input is probably prime ******\n", 
                           n.n, n.ndigits);

      mpgocandi_fixup_with_N (&go, &n);

      if (param != ECM_PARAM_DEFAULT && method != ECM_ECM)
        {
          fprintf (stderr, "Error, the -param option is only valid for ECM\n");
          exit (EXIT_FAILURE);
        }
      else if (param != ECM_PARAM_DEFAULT && !IS_BATCH_MODE (param) 
               && param != ECM_PARAM_SUYAMA && param != ECM_PARAM_WEIERSTRASS
               && param != ECM_PARAM_HESSIAN
	       && param != ECM_PARAM_TWISTED_HESSIAN)
        {
          fprintf (stderr, "Error, invalid -param value: %d\n", param);
          exit (EXIT_FAILURE);
        }
      /* params->param might be set if we resume from a file */
      if (params->param == ECM_PARAM_DEFAULT)
          params->param = param;

      /* this is a hack to produce an error in ecm() when -bsaves is used
         but we are not in batch mode */
      if (savefile_s != NULL)
        mpz_set_ui (params->batch_s, 2);

      /* load batch product s from a file */
      if (loadfile_s != NULL)
        {
          int st = cputime ();
          params->batch_last_B1_used = B1;
          if (read_s_from_file (params->batch_s, loadfile_s, B1))
            {
              fprintf (stderr, "Error while reading s from file\n");
              exit (EXIT_FAILURE);
            }
          else if (verbose >= OUTPUT_VERBOSE)
              fprintf (stdout, "Reading batch product (of %"PRIu64" bits) of "
                               "primes up to B1=%1.0f from %s took %ldms\n", 
                               mpz_sizeinbase (params->batch_s, 2), B1,
                               loadfile_s, cputime () - st);
        }

      /* set parameters that may change from one curve to another */
      params->method = method; /* may change with resume */
      mpz_set (params->x, x); /* may change with resume */
      mpz_set (params->y, y); /* may change with resume */
      /* already set when resumefile was read */
      if (resumefile == NULL)
	{
	  /* if A is not zero, we use it */
	  params->sigma_is_A = specific_A;
	  if (specific_A)
	    {
	      if (params->param == ECM_PARAM_WEIERSTRASS)
		{
		    /* compute B = y^2-x^3-A*x = y^2 - (x^2+A)*x */
		    mpz_mul (params->E->a6, y, y);
		    mpz_mul (params->E->a4, x, x);
		    mpz_add (params->E->a4, params->E->a4, A);
		    mpz_mul (params->E->a4, params->E->a4, x);
		    mpz_sub (params->E->a6, params->E->a6, params->E->a4);
		    mpz_mod (params->E->a6, params->E->a6, n.n);
		    mpz_set (params->E->a4, A);
		    params->sigma_is_A = -1;
		    params->E->type = ECM_EC_TYPE_WEIERSTRASS;
		    params->E->law = ECM_LAW_HOMOGENEOUS;
		}
	      else if (params->param == ECM_PARAM_HESSIAN)
		{
		    /* use x^3+y^3+1=3*A*x*y */
		    mpz_set (params->E->a4, A);
		    params->sigma_is_A = -1;
		    params->E->type = ECM_EC_TYPE_HESSIAN;
		    params->E->law = ECM_LAW_HOMOGENEOUS;
		}
	      else if (params->param == ECM_PARAM_TWISTED_HESSIAN)
		{
		    /* use a*x^3+y^3+1=d*x*y with A=a^3/d -- trick! */
		    mpz_set (params->E->a4, numer_A);
		    mpz_set (params->E->a6, denom_A);
		    mpz_clear(numer_A);
		    mpz_clear(denom_A);
		    params->sigma_is_A = -1;
		    params->E->type = ECM_EC_TYPE_TWISTED_HESSIAN;
		    params->E->law = ECM_LAW_HOMOGENEOUS;
		}
	    }
#ifdef HAVE_TORSION
	  else if (torsion != NULL)
	    {
	      params->param = ECM_PARAM_TORSION;
	      params->sigma_is_A = -1;
	      result = build_curves_with_torsion2 (f, n.n, params->E, 
						   params->x, params->y,
                                                   torsion, sigma);
	    }
#endif
	}
      mpz_set (params->sigma, (params->sigma_is_A) ? A : sigma);
      mpz_set (params->go, go.Candi.n); /* may change if contains N */
      mpz_set (params->B2min, B2min); /* may change with -c */
      /* Default, for P-1/P+1 with old stage 2 and ECM, use NTT only 
         for small input */
      if (use_ntt == 1 && (method == ECM_ECM || S != ECM_DEFAULT_S)) 
        params->use_ntt = (mpz_size (n.n) <= NTT_SIZE_THRESHOLD);
      else 
        params->use_ntt = use_ntt;

#ifdef HAVE_GWNUM
      /* check if the input number can be represented as k*b^n+c */
      if (kbnc_z (&gw_k, &gw_b, &gw_n, &gw_c, n.n))
        {
          params->gw_k = gw_k;
          params->gw_b = gw_b;
          params->gw_n = gw_n;
          params->gw_c = gw_c;
          if (verbose > OUTPUT_NORMAL)
            printf ("Found number: %.0f*%lu^%lu + %ld\n",
                    gw_k, gw_b, gw_n, gw_c);
        }
      else if (kbnc_str (&gw_k, &gw_b, &gw_n, &gw_c, n.cpExpr, n.n))
        {
          params->gw_k = gw_k;
          params->gw_b = gw_b;
          params->gw_n = gw_n;
          params->gw_c = gw_c;
          if (verbose > OUTPUT_NORMAL)
            printf ("Found number: %.0f*%lu^%lu + %ld\n",
                    gw_k, gw_b, gw_n, gw_c);
        }
      else
        {
          if (verbose > OUTPUT_NORMAL)
            printf ("Did not find a gwnum poly for the input number.\n");
        }
#endif

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

      /* now call the ecm library */
      if (result == ECM_NO_FACTOR_FOUND)
	  /* if torsion was used, some factor may have been found... */
	  result = ecm_factor (f, n.n, B1, params);

      if (result == ECM_ERROR)
        {
          fprintf (stderr, "Please report internal errors at <%s>.\n",
                   PACKAGE_BUGREPORT);
          exit (EXIT_FAILURE);
        }
      
      if (!params->gpu)
          cnt --; /* one more curve performed */
      else
        {
          if (cnt <= params->gpu_number_of_curves)
              cnt = 0;
          else
              cnt -= params->gpu_number_of_curves; 
        }

      /* When GPU is used we need to have the value of N before it is 
        divided by potential factor in f */
      mpz_init_set (tmp_n, n.n);

      if (result != ECM_NO_FACTOR_FOUND)
        {
          mpz_t tmp_factor;
          returncode = 0;
          mpz_init (tmp_factor);
          do 
            {
              if (params->gpu)
                  mpz_fdiv_qr (f, tmp_factor, f, tmp_n);
              else
                  mpz_set (tmp_factor, f);

              returncode = process_newfactor (tmp_factor, result, &n, method,
                                 returncode, params->gpu, &cnt, &resume_wasPrp,
                                 resume_lastfac, resumefile, verbose, deep);
            } while (params->gpu && mpz_cmp_ui (f, 0) != 0 
                                 && returncode != ECM_INPUT_NUMBER_FOUND);
          mpz_clear (tmp_factor);
        }

      /* if quiet mode, prints remaining cofactor after last curve */
      if ((cnt == 0) && (verbose == 0) && (mpz_cmp_ui (n.n, 1) > 0))
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
      if (savefilename != NULL && !n.isPrp)
        {
        /* TODO Deal with return code */
	    write_resumefile (savefilename, method, tmp_n, params, &n, 
			      orig_x0, orig_y0, comment);
        }

      mpz_clear (tmp_n);

      /* Save the batch exponent s if requested */
      if (savefile_s != NULL)
        {
          int ret = write_s_in_file (savefile_s, params->batch_s);
          if (verbose >= OUTPUT_VERBOSE && ret > 0)
              printf ("Saved batch product (of %u bytes) in %s\n", ret, 
                      savefile_s);
        }

      /* advance B1, if autoincrement value had been set during command line 
         parsing */
      if (autoincrementB1 > 0.0)
        {
          double NewB1;
          NewB1 = calc_B1_AutoIncrement (B1, autoincrementB1);
          if (mpz_cmp_d (B2min, B1) <= 0) /* <= might be better than == */
              mpz_set_d (B2min, NewB1);
          B1 = NewB1;
        }
    } /* end of main loop */

  free_all1:
  if (infilename) /* infile might be stdin, don't fclose that! */
    fclose (infile);

  if (resumefile)
    {
      fclose (resumefile);
      mpz_clear (resume_lastN);
      mpz_clear (resume_lastfac);
    }

  mpz_clear (orig_y0);
  mpz_clear (orig_x0);
  mpz_clear (y);
  mpz_clear (x);
  mpz_clear (f);
  mpcandi_t_free (&n);

 free_all:
  free_expr ();

  mpq_clear (rat_y0);
  mpq_clear (rat_x0);
  mpq_clear (rat_A);
  mpz_clear (startingB2min);
  mpz_clear (B2min);
  mpz_clear (B2);
  mpz_clear (A);
  mpz_clear (sigma);

  mpgocandi_t_free (&go);
  ecm_clear (params);

  /* exit 0 if a factor was found for the last input, except if we exit due
     to a signal */
#ifdef HAVE_SIGNAL
  if (returncode == 0 && exit_asap_value != 0)
    returncode = 143;
#endif

  /* print memory usage */
  if (verbose > 1)
    {
      long peak = PeakMemusage ();
      if (peak != -1)
        printf ("Peak memory usage: %ldMB\n", peak >> 10);
    }

  return returncode;
}
