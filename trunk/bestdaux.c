/* Choice of best parameters for stage 2.

  Copyright 2003 Paul Zimmermann and Alexander Kruppa.

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

#define COEFF double

typedef struct d_struct1
{
  unsigned long d;
  COEFF a, b;
  struct d_struct1 *next;
} d_struct;

COEFF toomcook3[82], *toomcook4, *PolyFromRoots, *RecursiveDivision,
  *PolyEval;

unsigned long phi (unsigned long);
void fill_toomcook3 (unsigned long);
void fill_toomcook4 (unsigned long);
COEFF list_mul (unsigned long, unsigned long);
void fill_PolyFromRoots (unsigned long);
void fill_RecursiveDivision (unsigned long);
void fill_PolyEval (unsigned long);
COEFF PolyInvert (unsigned long);
COEFF Reduce (unsigned long);
void cost (COEFF *, COEFF *, unsigned long);
d_struct* insert (d_struct *, unsigned long, COEFF, COEFF);
void print (d_struct *);
void print_d (d_struct *);

/* returns Euler's totient phi function */
unsigned long
phi (unsigned long n)
{
  unsigned int phi = 1, p;

  for (p = 2; p * p <= n; p += 2)
    {
      if (n % p == 0)
	{
	  phi *= p - 1;
	  n /= p;
	  while (n % p == 0)
	    {
	      phi *= p;
	      n /= p;
	    }
	}

      if (p == 2)
	p--;
    }

  /* now n is prime */

  return (n == 1) ? phi : phi * (n - 1);
}

void
fill_toomcook3 (unsigned long maxn)
{
  unsigned long n;
  COEFF v;

  for (n = 0; n <= maxn; n++)
    {
      if (n == 0)
        v = (COEFF) 0;
      else if (n == 1)
        v = (COEFF) 1;
      else if (n == 2)
        v = (COEFF) 3;
      else
        {
          unsigned long l, k;
          
          l = (n + 2) / 3;
          k = n - 2 * l;
          v = (COEFF) 4 * toomcook3[l] + toomcook3[k];
        }
      toomcook3[n] = v;
    }
}

void
fill_toomcook4 (unsigned long maxn)
{
  unsigned long n;
  COEFF v;

  for (n = 0; n <= maxn; n++)
    {
      if (n <= 3 || n == 5 || n == 6 || n == 9 || n == 17 || n == 18 ||
          n == 26 || n == 27 || (77 <= n && n <= 81)) /* use toomcook3 */
        {
          v = toomcook3[n];
        }
      else /* use toomcook4 */
        {
          unsigned long l, k;

          l = (n + 3) / 4;
          k = n - 3 * l;
          v = (COEFF) 6 * toomcook4[l] + toomcook4[k];
        }
      toomcook4[n] = v;
    }
}

/* assume k >= l */
COEFF list_mul
(unsigned long k, unsigned long l)
{
  if (k == l)
    return toomcook4[l];
  else if (k == l + 1)
    return toomcook4[l] + (COEFF) l;
  else
    abort ();
}

void
fill_PolyFromRoots (unsigned long n)
{
  unsigned long k, m, l;
  COEFF v;

  for (k = 0; k <= n; k++)
    {
      if (k <= 1)
        v = (COEFF) 0;
      else if (k == 2)
        v = (COEFF) 1;
      else
        {
          m = k / 2;
          l = k - m;
          v = PolyFromRoots[l] + PolyFromRoots[m] + list_mul (l, m);
        }
      PolyFromRoots[k] = v;
    }
}

void
fill_RecursiveDivision (unsigned long n)
{
  unsigned long k, m, l;
  COEFF v;

  for (k = 0; k <= n; k++)
    {
      if (k <= 1)
        v = (COEFF) k;
      else
        {
          m = k / 2;
          l = k - m;
          v = RecursiveDivision[l] + RecursiveDivision[m] + (COEFF)
            2 * toomcook4[m] + (COEFF) 2 * m * (l - m);
        }
      RecursiveDivision[k] = v;
    }
}

void
fill_PolyEval (unsigned long n)
{
  unsigned long k, m, l;
  COEFF v;

  for (k = 0; k <= n; k++)
    {
      if (k <= 1)
        v = (COEFF) 0;
      else
        {
          m = k / 2;
          l = k - m;
          v = RecursiveDivision[m];
          if (l > m)
            v += (COEFF) m;
          v += RecursiveDivision[l] + PolyEval[l] + PolyEval[m];
        }
      PolyEval[k] = v;
    }
}

COEFF
PolyInvert (unsigned long n)
{
  if (n <= 1)
    return (COEFF) 0;
  else
    {
      unsigned long k, l;
      COEFF v;
      
      k = n / 2;
      l = n - k;
      v = PolyInvert (l) + toomcook4[l] + toomcook4[k];
      if (k > 1)
        v += list_mul (l - 1, k - 1);
      return v;
    }
}

COEFF
Reduce (unsigned long n)
{
  return (COEFF) 2 * toomcook4[n-1] + (COEFF) 1;
}     

/* returns a, b such that the number of multiplies of step 2 is
   a + b * ceil(B2/(d * phi(d)/2)), i.e.
   (k+1)*PolyFromRoots(dF)+PolyInvert(dF-1)+(k-1)*toomcook4(dF)
  +(k-1)*Reduce(dF)+polyeval(dF) where dF = phi(d)/2.
*/
void
cost (COEFF *a, COEFF *b, unsigned long d)
{
  unsigned long dF;

  dF = phi (d) / 2;

  *a = PolyFromRoots[dF] + PolyInvert (dF - 1) - toomcook4[dF] - Reduce (dF)
    + PolyEval[dF];
  *b = PolyFromRoots[dF] + toomcook4[dF] + Reduce (dF);
  *b = *b / (COEFF) d / (COEFF) dF;
}

/* assumes values of a (resp. b) are increasing (resp. decreasing) in l */
d_struct*
insert (d_struct *l, unsigned long d, COEFF a, COEFF b)
{
  d_struct *oldl = NULL, *l0 = l, *newl;

  /*  printf ("inserting d=%u a=%1.0f b=%e\n", d, a, b); */

  if (l != NULL && l->a <= a && l->b <= b)
    {
      /*      printf ("d=%u a=%1.0f b=%e is better\n", l->d, l->a, l->b); */
      return l0;
    }

  /* now either l->a > a or l->b > b */

  while (l != NULL && l->a < a && l->b >= b)
    {
      /*      printf ("skipping d=%u a=%1.0f b=%e\n", l->d, l->a, l->b); */
      oldl = l; /* previous cell */
      l = l->next;
      if (l != NULL && l->a <= a && l->b <= b)
        {
          /*          printf ("d=%u a=%1.0f b=%e is better\n", l->d, l->a, l->b); */
          return l0;
        }
    }
  /* either l=NULL, or l->a > a, or l->b < b:
     (a) l->a > a, l->b <= b
     (b) l->a > a, l->b > b */

  /* if a < l->a and b < l->b, remove cell in l */
  while (l != NULL && a <= l->a && b < l->b)
    {
      /*      printf ("discard d=%u a=%1.0f b=%e\n", l->d, l->a, l->b); */
      if (oldl != NULL)
        {
          oldl->next = l->next;
          free (l);
          l = oldl->next;
        }
      else
        {
          l0 = l->next;
          free (l);
          l = l0;
        }
    }

  /* now either l=NULL, or (l->a > a and l->b <= b) */
  newl = (d_struct*) malloc (sizeof (d_struct));
  newl->d = d;
  newl->a = a;
  newl->b = b;
  newl->next = l;
  if (oldl != NULL)
    oldl->next = newl;

  return (oldl == NULL) ? newl : l0;
}

void
print (d_struct *l)
{
  while (l != NULL)
    {
      printf ("%lu %1.0f %e\n", l->d, l->a, l->b);
      l = l->next;
    }
}

void
print_d (d_struct *l)
{
  unsigned long count = 0;
  d_struct *m;

  for (m = l; m != NULL; count++, m = m->next);

  printf ("static unsigned l[%lu] = {", count);
  while (l != NULL)
    {
      printf ("{%lu,%1.0f,%1.0f}", l->d, l->a, l->b * l->d * phi(l->d) / 2.0);
      l = l->next;
      if (l != NULL)
        printf (", ");
    }
  printf ("};\n");
}

int
main (int argc, char *argv[])
{
  unsigned long dmax;
  unsigned long maxdF;
  unsigned long d;
  COEFF a, b;
  d_struct *l = NULL;

  if (argc == 1)
    {
      fprintf (stderr, "Usage: bestdaux dmax\n");
      exit (1);
    }

  dmax = atoi (argv[1]);

  if (dmax % 6 != 0)
    {
      fprintf (stderr, "Error, dmax should be multiple of 6\n");
      exit (1);
    }

  /* since d = 6*k, we have phi(d) <= phi(6)*k = 2*k <= d/3 
     so dF = phi(d)/2 <= dmax/6 */

  maxdF = dmax / 6;

  toomcook4 = (COEFF*) malloc ((maxdF + 1) * sizeof (COEFF));
  PolyFromRoots = (COEFF*) malloc ((maxdF + 1) * sizeof (COEFF));
  RecursiveDivision = (COEFF*) malloc ((maxdF + 1) * sizeof (COEFF));
  PolyEval = (COEFF*) malloc ((maxdF + 1) * sizeof (COEFF));

  fill_toomcook3 (81);
  fill_toomcook4 (maxdF);
  fill_PolyFromRoots (maxdF);
  fill_RecursiveDivision (maxdF);
  fill_PolyEval (maxdF);

  for (d = 6; d <= dmax; d += 6)
    {
      cost (&a, &b, d);
      l = insert (l, d, a, b);
    }

  print_d (l);

  free (PolyEval);
  free (RecursiveDivision);
  free (PolyFromRoots);
  free (toomcook4);

  return 0;
}
