/* ecmfactor.c - example of use of libecm.a.

Copyright 2005-2016 Paul Zimmermann, Dave Newman, Pierrick Gaudry.

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
#include <pthread.h>
#include <gmp.h> /* GMP header file */
#include "ecm.h" /* ecm header file */

/* thread structure */
typedef struct
{
  mpz_t n;        /* number to factor */
  double B1;      /* stage-1 bound */
  mpz_t f;        /* potential factor */
  int ret;        /* return value */
  ecm_params q;
} __tab_struct;
typedef __tab_struct tab_t[1];

void*
one_thread (void *args)
{
  tab_t *tab = (tab_t*) args;

  tab[0]->ret = ecm_factor (tab[0]->f, tab[0]->n, tab[0]->B1, tab[0]->q);
  return NULL;
}

int
main (int argc, char *argv[])
{
  mpz_t n;
  double B1;
  unsigned long nthreads = 1, i;
  tab_t *T;
  pthread_t *tid;

  if (argc >= 3 && strcmp (argv[1], "-t") == 0)
    {
      nthreads = strtoul (argv[2], NULL, 10);
      argc -= 2;
      argv += 2;
    }

  if (argc < 3)
    {
      fprintf (stderr, "Usage: ecmfactor [-t nnn] <number> <B1>\n");
      exit (1);
    }

  /* initialize tab_t for threads */
  T = malloc (nthreads * sizeof (tab_t));
  tid = malloc (nthreads * sizeof (pthread_t));

  mpz_init (n);

  /* read number on command line */
  if (mpz_set_str (n, argv[1], 10))
    {
      fprintf (stderr, "Invalid number: %s\n", argv[1]);
      exit (1);
    }

  B1 = atof (argv[2]);

  for (i = 0; i < nthreads ; i++)
    {
      mpz_init_set (T[i]->n, n);
      T[i]->B1 = B1;
      mpz_init (T[i]->f);
      ecm_init (T[i]->q);
    }

  printf ("Performing %lu curve(s) with B1=%1.0f\n", nthreads, B1);
  for (i = 0; i < nthreads; i++)
    pthread_create (&tid[i], NULL, one_thread, (void *) (T+i));
  for (i = 0; i < nthreads; i++)
    pthread_join (tid[i], NULL);

  for (i = 0; i < nthreads; i++)
    {
      gmp_printf ("thread %lu with sigma %d:%Zd ",
                  i, T[i]->q->param, T[i]->q->sigma);
      if (T[i]->ret > 0)
        gmp_printf ("found factor in step %u: %Zd\n", T[i]->ret, T[i]->f);
      else if (T[i]->ret == ECM_NO_FACTOR_FOUND)
        printf ("found no factor\n");
      else
        printf ("gave an error\n");
      mpz_clear (T[i]->n);
      mpz_clear (T[i]->f);
      ecm_clear (T[i]->q);
    }

  mpz_clear (n);

  free (T);
  free (tid);

  return 0;
}
