/* Encapsulated candidate.  This candidate should have been a C++ class, but
   since we are using straight C for this project, I guess I can deal with it.

Copyright 2003, 2004, 2005, 2006 Jim Fougeron, Paul Zimmermann.

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
#include "ecm-ecm.h"

#define VALID_MAGIC 0x24837BF5
#define DEAD_MAGIC  0xDEADBEEF

#if defined (CANDI_DEBUG)
static void
Candi_Validate (const char *FunctionStr, const mpcandi_t *n)
{
  int abrt = 0;
  if (!FunctionStr)
    {
      fprintf (stderr, "ERROR, UNKNOWN FUNCTION, can NOT continue checks!\n");
      exit(-1);
    }
  if (!n)
    {
      abrt = fprintf (stderr, "ERROR, %s() *n was NULL, can NOT continue checks!\n", FunctionStr);
      exit(-1);
    }
  if (n->magic != VALID_MAGIC)
    abrt = fprintf (stderr, "ERROR, %s() VALID_MAGIC not valid\n", FunctionStr);
  if (n->cpExpr && n->nexprlen != strlen(n->cpExpr))
    abrt = fprintf (stderr, "ERROR, %s() Invalid cpExpr length\n", FunctionStr);
  if (n->ndigits != nb_digits(n->n))
	abrt = fprintf (stderr, "ERROR, %s() Invalid n->ndigits length\n", FunctionStr);
  if (abrt)
    exit(-1);
}
#endif

void
mpcandi_t_init (mpcandi_t *n)
{
  n->cpExpr = NULL;
  n->nexprlen = 0;
  n->ndigits = 1;
  mpz_init_set_ui (n->n, 1);
  n->isPrp = 0;
#if defined (CANDI_DEBUG)
  n->magic = VALID_MAGIC;
  Candi_Validate ("mpcandi_t_init", n);
#endif
}

void
mpcandi_t_free (mpcandi_t *n)
{
#if defined (CANDI_DEBUG)
  Candi_Validate("mpcandi_t_free", n);
#endif
  if (n->cpExpr)
    free (n->cpExpr);
  n->cpExpr = NULL;
  n->nexprlen = 0;
  n->ndigits = 0;
  mpz_clear (n->n);
  n->isPrp = 1;	  /* "default" to prp, so that if the candidate does not get
		     filled in, it will not be tested */
#if defined (CANDI_DEBUG)
  n->magic = DEAD_MAGIC;
#endif
}

int
mpcandi_t_add_candidate (mpcandi_t *n, mpz_t c, const char *cpExpr,
			 int primetest)
{
#if defined (CANDI_DEBUG)
  Candi_Validate("Pre mpcandi_t_add_candidate", n);
#endif

  if (n->cpExpr)
    free (n->cpExpr);
  n->cpExpr = NULL;
  if (cpExpr)
    {
      n->nexprlen = strlen (cpExpr);
      n->cpExpr = (char *) malloc (n->nexprlen + 1);
      ASSERT_ALWAYS (n->cpExpr != NULL);
      strcpy (n->cpExpr, cpExpr);
    }
  mpz_set (n->n, c);
  n->ndigits = nb_digits (c);
  if (primetest)
    {
      if (n->ndigits < APRCL_CUTOFF)
        n->isPrp = mpz_aprtcle (c, 0);
      else if (n->ndigits < APRCL_CUTOFF2)
        n->isPrp = mpz_aprtcle (c, 1);
      else
        n->isPrp = mpz_probab_prime_p (c, PROBAB_PRIME_TESTS);
    }
  else
    n->isPrp = 0; /* there is a candidate there now, and the user did not
		     tell us to prp it, so assume it is composite */

#if defined (CANDI_DEBUG)
  Candi_Validate("Post mpcandi_t_add_candidate", n);
#endif

  return 1;
}

int
mpcandi_t_addfoundfactor (mpcandi_t *n, mpz_t f, int displaywarning)
{
#if defined (CANDI_DEBUG)
  Candi_Validate("Pre mpcandi_t_addfoundfactor", n);
#endif
  char *cp, *cp1;

  if (!mpz_divisible_p (n->n, f))
    {
      /* ERROR was not a factor NOTE however, that this is "valid" for the 
         ui() function to call. When trial dividing, it is VERY frequent to 
         be divisible by 2^3, and we try to remove factors UNTIL */
      if (displaywarning)
        gmp_fprintf (stderr, "ECM logic ERROR. Trying to remove a "
                     "non-factor %Zd\n", f);
#if defined (CANDI_DEBUG)
      Candi_Validate("Post (no factor removed) mpcandi_t_addfoundfactor", n);
#endif
      return 0;
    }

  /* remove f from n->n */
  mpz_divexact (n->n, n->n, f);
  n->ndigits = nb_digits (n->n);
  if (n->ndigits < APRCL_CUTOFF)
    n->isPrp = mpz_aprtcle (n->n, 0);
  else if (n->ndigits < APRCL_CUTOFF2)
    n->isPrp = mpz_aprtcle (n->n, 1);
  else
    n->isPrp = mpz_probab_prime_p (n->n, PROBAB_PRIME_TESTS);
  if (n->cpExpr != NULL)
    {
      /* If there is an expression, then lets preserve it */
      cp1 = mpz_get_str (NULL, 10, f);
      cp = (char *) malloc(n->nexprlen+1 + 3 + strlen(cp1));  /* +1 for null, +3 for ()/ */
      ASSERT_ALWAYS (cp != NULL);
      sprintf (cp, "(%s)/%s", n->cpExpr, cp1);
      free(n->cpExpr);
      n->cpExpr = cp;
      n->nexprlen += (3+strlen(cp1));
      free (cp1); /* size strlen (cp1) + 1 */
    }
#if defined (CANDI_DEBUG)
  Candi_Validate("Post (removed factor) mpcandi_t_addfoundfactor", n);
#endif
  return 1;
}

/**********************************************************************
  Group order candidate functions.  These wrap the logic for the -go
  command line switch which allows the user to "insert" the proper
  group order.
**********************************************************************/
void
mpgocandi_t_init (mpgocandi_t *go)
{
  go->cpOrigExpr = NULL;
  mpcandi_t_init (&(go->Candi));
  go->containsN = 0;
  go->Valid = 0;
}

void
mpgocandi_t_free (mpgocandi_t *go)
{
  if (go->cpOrigExpr)
    free (go->cpOrigExpr);
  mpcandi_t_free (&(go->Candi));
  go->Valid = 0;
}

int
mpgocandi_fixup_with_N (mpgocandi_t *go, mpcandi_t *n)
{
  int NumNs, len;
  char *cp, *cpo, *numbuf;

  if (go->Valid == 0)
    return 0;
  if (go->containsN == 0)
    return 1;  /* a valid "normal" expression does not need updating */

  cp = strchr (go->cpOrigExpr, 'N');
  NumNs = 0;
  while (cp)
    {
      ++NumNs;
      cp = strchr (&cp[1], 'N');
    }
  /* compute size of string needed, and add some safety buffer to it */
  cp = go->cpOrigExpr;
  len = NumNs * mpz_sizeinbase (n->n, 10) + strlen (cp) + 100;
  numbuf = (char *) malloc(len);
  ASSERT_ALWAYS (numbuf != NULL);
  cpo = numbuf;
  while (*cp)
    {
      if (*cp == 'N')
	  cpo += gmp_sprintf (cpo, "%Zi", n->n);
      else
        *cpo++ = *cp;
      ++cp;
    }

  *cpo = 0; /* Null terminate the string correctly. */

  go->Valid = eval_str (&(go->Candi), numbuf, 0, NULL);
  /* for invalid inputs, eval_str() directly exit(1) */

  free (numbuf);
  return go->Valid;
}
