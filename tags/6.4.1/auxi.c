/* Auxiliary functions for GMP-ECM.

  Copyright 2002, 2003, 2004, 2005, 2007 Paul Zimmermann, Alexander Kruppa, Laurent Fousse, Jim Fougeron.

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

#include <gmp.h>
#include "ecm-ecm.h"

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

