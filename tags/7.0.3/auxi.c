/* Auxiliary functions for GMP-ECM.

Copyright 2002, 2003, 2004, 2005, 2006, 2007, 2011, 2012 Paul Zimmermann,
Alexander Kruppa, Laurent Fousse, Jim Fougeron, Cyril Bouvier.

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
#include "ecm-impl.h"
#include "ecm-ecm.h"

#ifdef HAVE_GWNUM
/* For GWNUM_VERSION */
#include "gwnum.h"
#endif

#include "champions.h"

/******************************************************************************
*                                                                             *
*                            Auxiliary functions                              *
*                                                                             *
******************************************************************************/

/* returns the number of decimal digits of n */
unsigned int
nb_digits (const mpz_t n)
{
  mpz_t x;
  unsigned int size;

  size = mpz_sizeinbase (n, 10);

  /* the GMP documentation says mpz_sizeinbase returns the exact value,
     or one too big, thus:
     (a) either n < 10^(size-1), and n has size-1 digits
     (b) or n >= size-1, and n has size digits
     Note: mpz_sizeinbase returns 1 for n=0, thus we always have size >= 1.
  */
				    
  mpz_init (x);
  mpz_ui_pow_ui (x, 10, size - 1);
  if (mpz_cmpabs (n, x) < 0)
    size --;
  mpz_clear (x);

  return size;
}

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
process_newfactor (mpz_t g, int result, mpcandi_t *n, int method, 
                   int returncode, int gpu, unsigned int *cnt, 
                   int *resume_wasPrp, mpz_t resume_lastfac, 
                   FILE *resumefile, int verbose, int deep)
{
  int factor_is_prime = 0;
        /* If a factor was found, indicate whether factor, cofactor are */
        /* prime. If no factor was found, both are zero. */
  int cofactor_is_prime = 0;
  int method1;
  int found_n;
  mpz_t f;
  
  mpz_init (f);

  mpz_gcd (f, g, n->n);

    /* When GPU is not used, the factor should divide n->n */
    if (gpu == 0)
      ASSERT_ALWAYS(mpz_cmp (g, f) == 0);
    else if (mpz_cmp (g, f) != 0 && mpz_cmp_ui (f, 1) == 0)
    /* GPU case: all factors of g were already found */
      {
        /* FIXME Maybe print something in very verbose mode */
        mpz_clear (f);
        return returncode;
      }

    /* Restoring print statement that was disabled from previous code clean-up efforts  */
    /* This print statement is important to the function of external factoring programs */
    /*    that rely on the following information being printed out */
    /* g = f (gpu or not gpu) or g != 1 with gpu */
    if (mpz_cmp (g, f) == 0 || (gpu != 0 && mpz_cmp_ui (g, 1) != 0))
      {
        if (verbose > 0)
            printf ("********** Factor found in step %u: ", ABS (result));

        mpz_out_str (stdout, 10, f);
        
        if (verbose > 0)
            printf ("\n");
      }

  /* Complain about non-proper factors (0, negative) */
  ASSERT_ALWAYS(mpz_cmp_ui (f, 1) >= 1);

  found_n = mpz_cmp (f, n->n) == 0;
  
  /* 1 for display warning if factor does not divide the current 
     candidate */
  mpcandi_t_addfoundfactor (n, f, gpu == 0);

  if (found_n == 0)
    {
      /* prints factor found and cofactor on standard output. */
      /* the second argument here tells aprtcle to not print primality proving progress info */
      factor_is_prime = mpz_aprtcle (f, 0);

      if (verbose >= 1)
        {
          if (factor_is_prime == ECM_FAC_PRIME)
            printf ("Found prime factor of %u digits: ", nb_digits (f));
          else if (factor_is_prime == ECM_FAC_PRP)
            printf ("Found probable prime factor of %u digits: ", nb_digits (f));
          else
            printf ("Found composite factor of %u digits: ", nb_digits (f));
          mpz_out_str (stdout, 10, f);
          printf ("\n");
        }
      
      if (resumefile != NULL)
        {
          /* If we are resuming from a save file, add factor to the
             discovered factors for the current number */
          mpz_mul (resume_lastfac, resume_lastfac, f);
          *resume_wasPrp = n->isPrp;
        }

      if (factor_is_prime && returncode == 0)
        returncode = (n->isPrp) ? ECM_PRIME_FAC_PRIME_COFAC : 
                      ECM_PRIME_FAC_COMP_COFAC;
      else
        returncode = (n->isPrp) ? ECM_COMP_FAC_PRIME_COFAC :
                      ECM_COMP_FAC_COMP_COFAC;

      if (verbose >= 1)
        {
          cofactor_is_prime = n->isPrp; /* from mpcandi_t_addfoundfactor */

          if (cofactor_is_prime == ECM_FAC_PRIME)
            printf ("Prime cofactor ");
          else if (cofactor_is_prime == ECM_FAC_PRP)
            printf ("Probable prime cofactor ");
          else
            printf ("Composite cofactor ");

          if (n->cpExpr)
            printf ("%s", n->cpExpr);
          else
            mpz_out_str (stdout, 10, n->n);
          printf (" has %u digits\n", n->ndigits);
        }
      else /* quiet mode: just print a space here, remaining cofactor
              will be printed after last curve */
          printf (" ");

      /* check for champions (top ten for each method) */
      method1 = ((method == ECM_PP1) && (result < 0)) ? 
                ECM_PM1 : method;
      if ((verbose > 0) && factor_is_prime && 
          nb_digits (f) >= champion_digits[method1])
        {
          printf ("Report your potential champion to %s\n",
                  champion_keeper[method1]);
          printf ("(see %s)\n", champion_url[method1]);
        }

      /* Take care of fully factoring this number, 
         in case we are in deep mode */
      if (n->isPrp)
        *cnt = 0; /* no more curve to perform */

      if (!deep)
          *cnt = 0;
    }
  else
    {
      *cnt = 0; /* no more curve to perform */
      if (verbose > 0)
          printf ("Found input number N");
      printf ("\n");
      returncode = ECM_INPUT_NUMBER_FOUND;
    }
  fflush (stdout);

  mpz_clear (f);
  return returncode;
}
