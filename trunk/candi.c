/* Encapsulated candidate.  This candidate should have been a C++ class, but
   since we are using straight C for this project, I guess I can deal with it.

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if !defined (_MSC_VER)
#include <unistd.h>
#endif
#include "gmp.h"
#include "ecm.h"
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

/* performs a safe "deep" copy */
int
mpcandi_t_copy (mpcandi_t *to, mpcandi_t *from)
{
#if defined (CANDI_DEBUG)
  Candi_Validate("Pre mpcandi_t_copy", to);
  Candi_Validate("Pre mpcandi_t_copy", from);
#endif
  if (to == from)
     return 1;
  if (to->cpExpr)
    free(to->cpExpr);
  to->cpExpr = NULL;
  if (from->cpExpr)
    {
      to->cpExpr = (char *) malloc(from->nexprlen+1);
      if (to->cpExpr == NULL)
        {
          fprintf (stderr, "Error: not enough memory\n");
          exit (EXIT_FAILURE);
        }
      strcpy(to->cpExpr, from->cpExpr);
    }
  to->nexprlen = from->nexprlen;
  mpz_set(to->n, from->n);
  to->isPrp = from->isPrp;
  to->ndigits = from->ndigits;

#if defined (CANDI_DEBUG)
  Candi_Validate("Post mpcandi_t_copy", to);
  Candi_Validate("Post mpcandi_t_copy", from);
#endif

  return 1;
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
      if (n->cpExpr == NULL)
        {
          fprintf (stderr, "Error: not enough memory\n");
          exit (EXIT_FAILURE);
        }
      strcpy (n->cpExpr, cpExpr);
    }
  mpz_set (n->n, c);
  if (primetest)
    n->isPrp = smart_probab_prime_p (c, PROBAB_PRIME_TESTS);
  else
    n->isPrp = 0; /* there is a candidate there now, and the user did not
		     tell us to prp it, so assume it is composite */
  n->ndigits = nb_digits (c);

#if defined (CANDI_DEBUG)
  Candi_Validate("Post mpcandi_t_add_candidate", n);
#endif

  return 1;
}

int
mpcandi_t_addfoundfactor_d (mpcandi_t *n, double f)
{
#if defined (CANDI_DEBUG)
  Candi_Validate("Pre mpcandi_t_addfoundfactor_d", n);
#endif
  int ret;
  mpz_t t;
  mpz_init_set_d(t,f);
  /* do not display a warning if this factor does not divide the remaining
     cofactor. This function is called repeatedly (until it fails) to remove
     all traces of 
     the prime factor.  It is highly likely that these smaller factors will be
     non square-free within the candidate when starting.  A return of zero is 
     exprected by the calling trial divider, as that tells it that all residue 
     of the factor has been eliminated */
  ret = mpcandi_t_addfoundfactor (n, t, 0);
  mpz_clear (t);

#if defined (CANDI_DEBUG)
  Candi_Validate("Post mpcandi_t_addfoundfactor_d", n);
#endif

  return ret;
}

int
mpcandi_t_addfoundfactor (mpcandi_t *n, mpz_t f, int displaywarning)
{
#if defined (CANDI_DEBUG)
  Candi_Validate("Pre mpcandi_t_addfoundfactor_d", n);
#endif
  mpz_t t;
  char *cp, *cp1;

  mpz_init(t);
  mpz_mod(t, n->n, f);
  if (mpz_cmp_ui(t, 0))
    {
      /* ERROR was not a factor NOTE however, that this is "valid" for the ui() function to call. When trial dividing, it 
         is VERY frequent to be divisible by 2^3, and we try to remove factors UNTIL */
      mpz_clear(t);
      if (displaywarning)
        fprintf (stderr, "\nECM logic ERROR.  Trying to remove a non-factor\n");
#if defined (CANDI_DEBUG)
      Candi_Validate("Post (no factor removed) mpcandi_t_addfoundfactor_d", n);
#endif
      return 0;
    }
  mpz_clear(t);
  /* remove f from n->n */
  mpz_divexact (n->n, n->n, f);
  n->ndigits = nb_digits (n->n);
  n->isPrp = smart_probab_prime_p (n->n, PROBAB_PRIME_TESTS);
  if (n->cpExpr != NULL)
    {
      /* If there is an expression, then lets preserve it */
      cp1 = mpz_get_str (NULL, 10, f);
      cp = (char *) malloc(n->nexprlen+1 + 3 + strlen(cp1));  /* +1 for null, +3 for ()/ */
      if (cp == NULL)
        {
          fprintf (stderr, "Error: not enough memory\n");
          exit (EXIT_FAILURE);
        }
      sprintf (cp, "(%s)/%s", n->cpExpr, cp1);
      free(n->cpExpr);
      n->cpExpr = cp;
      n->nexprlen += (3+strlen(cp1));
      FREE (cp1, strlen (cp1) + 1);
    }
#if defined (CANDI_DEBUG)
  Candi_Validate("Post (removed factor) mpcandi_t_addfoundfactor_d", n);
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
  if (numbuf == NULL)
    {
      fprintf (stderr, "Error: not enough memory\n");
      exit (EXIT_FAILURE);
    }
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

  if (eval_str (&(go->Candi), numbuf, 0, NULL))
    go->Valid = 1;
  else
    {
      static int warned = 0;
      if (!warned)
	{
	  warned = 1;
	  fprintf(stderr, "Warning, invalid expression %s for the -go option\n", go->cpOrigExpr);
	}
      go->Valid = 0;  /* it is not valid, so do not use it */
    }

  free (numbuf);
  return go->Valid;
}
