/* Median product.

  Copyright 2003 Laurent Fousse, Paul Zimmermann and Alexander Kruppa.

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

References: 
[1] Tellegen's Principle into Practice, by A. Bostan, G. Lecerf and E. Schost,
Proc. of ISSAC'03, Philadelphia, 2003.
*/

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "ecm.h"

#ifndef MAX
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif

void list_add_wrapper (listz_t p, listz_t q, listz_t r, unsigned int n,
                       unsigned int max_r)
{
    list_add (p, q, r, MIN (n, max_r));
    if (n > max_r) 
        list_set (p + max_r, q + max_r, n - max_r);
}

void list_sub_wrapper (listz_t p, listz_t q, listz_t r, unsigned int n,
                       unsigned int max_r)
{
    list_sub (p, q, r, MIN (n, max_r));
    if (n > max_r) 
        list_set (p + max_r, q + max_r, n - max_r);
}

/* Given a[0..m] and c[0..l], puts in b[0..n] the coefficients
   of degree m to n+m of rev(a)*c, i.e.
   b[0] = a[0]*c[0] + ... + a[i]*c[i] with i = min(m, l)
   ...
   b[k] = a[0]*c[k] + ... + a[i]*c[i+k] with i = min(m, l-k)
   ...
   b[n] = a[0]*c[n] + ... + a[i]*c[i+n] with i = min(m, l-n) [=l-n].
   Using auxiliary memory in t.
   Implements algorithm TKarMul of [1].
   Assumes deg(c) = l <= m+n.
*/

unsigned int
TKarMul (listz_t b, unsigned int n,
	 listz_t a, unsigned int m, listz_t c, unsigned int l, listz_t t)
{
  unsigned int k, mu, nu, h;
  unsigned int s1;
  unsigned tot_muls = 0;
#ifdef TKARMULDEBUG
  printf ("Entrée dans TKarMul.\nm = %d\nn = %d\nl = %d\ndepth = %d\n", m, n,
                                l, depth);
  printf ("a = ");
  print_list (a, m + 1);
  printf ("\nc = ");
  print_list (c, l + 1);
  printf ("\n");
#endif

  
  if (n == 0)
    {
#ifdef TKARMULDEBUG
      printf ("Cas n = 0.\n");
#endif
      mpz_mul (b[0], a[0], c[0]);
      for (k = 1; (k <= m) && (k <= l); k++)
	mpz_addmul (b[0], a[k], c[k]);
#ifdef TKARMULDEBUG
      printf ("Sortie de TKarMul.\n");
      show_result(a, m, c, l, b, n);
#endif
      return MIN (m, l) + 1;
    }

  if (m == 0)
    {
#ifdef TKARMULDEBUG
      printf ("Cas m = 0.\n");
#endif
      for (k = 0; (k <= l) && (k <= n); k++)
	mpz_mul (b[k], a[0], c[k]);
      for (k = l + 1; k <= n; k++)
	mpz_set_ui (b[k], 0);
#ifdef TKARMULDEBUG
      printf ("Sortie de TKarMul.\n");
      show_result(a, m, c, l, b, n);
#endif
      return MIN (n, l) + 1;
    }

  mu = (m / 2) + 1;		/* 1 <= mu <= m */
  nu = (n / 2) + 1;		/* 1 <= nu <= n */
  h = MAX (mu, nu);		/* h >= 1 */

#ifdef TKARMULDEBUG
  printf ("mu = %d\nnu = %d\nh = %d\n", mu, nu, h);
#endif

  if (mu > n)
    {
#ifdef TKARMULDEBUG
      printf ("Cas mu > n.\n");
#endif

      tot_muls += TKarMul (b, n, a, mu - 1, c, l, t);
      if (l >= mu)
	{
	  /* we have to check l-mu <= n + (m-mu), i.e. l <= n+m */
	  tot_muls += TKarMul (t, n, a + mu, m - mu, c + mu, l - mu, t + n + 1);
	  list_add (b, b, t, n + 1);
	}
#ifdef TKARMULDEBUG
      printf ("Sortie de TKarMul.\n");
      show_result(a, m, c, l, b, n);
#endif
      return tot_muls;
    }

  if (nu > m)
    {
#ifdef TKARMULDEBUG
      printf ("Cas nu > m.\n");
#endif

      /* we have to check MIN(l,m+nu-1) <= nu-1+m: trivial */
      tot_muls += TKarMul (b, nu - 1, a, m, c, MIN (l, m + nu - 1), t);

      /* Description broken in reference. Should be a list
       * concatenation, not an addition.
       * Fixed now.
       */

      if (l >= nu)
	{
	  /* we have to check l-nu <= n-nu+m, i.e. l <= n+m: trivial */
	  tot_muls += TKarMul (b + nu, n - nu, a, m, c + nu, l - nu, t);
	}
      else
        list_zero (b + nu, n - nu + 1);
#ifdef TKARMULDEBUG
      printf ("Sortie de TKarMul.\n");
      show_result(a, m, c, l, b, n);
#endif
      return tot_muls;
    }

  /* We want nu = mu */

  mu = nu = h;
  
#ifdef TKARMULDEBUG
  printf ("Cas de base.\n");
#endif
  
  s1 = MIN (l + 1, n + mu);
  if (l + 1 > nu)
        list_sub_wrapper (t, c, c + nu, s1, l - nu + 1);
      else
        list_set (t, c, s1);
#ifdef TKARMULDEBUG
      printf ("DEBUG c - c[nu].\n");
      print_list (t, s1);
      printf ("On calcule (1) - (3)\n");
#endif
      tot_muls += TKarMul (b, nu - 1, a, mu - 1, t, s1 - 1, t + s1);
      /* (1) - (3) */
#ifdef TKARMULDEBUG
      print_list (b, nu);
      printf ("On calcule (2) - (4)\n");
#endif
      if (s1 >= nu + 1) { /* nu - 1 */
        tot_muls += TKarMul (b + nu, n - nu, a + mu, m - mu, 
                             t + nu, s1 - nu - 1, t + s1);
        /* (2) - (4) */
      }
      else {
          list_zero (b + nu, n - nu + 1);
      }
#ifdef TKARMULDEBUG
      print_list (b + nu, n - nu + 1);
#endif
      list_add_wrapper (t, a, a + mu, mu, m + 1 - mu);
#ifdef TKARMULDEBUG
      printf ("On calcule (2) + (3)\n");
#endif
      if (l >= nu) {
          tot_muls += TKarMul (t + mu, nu - 1, t, mu - 1, c + nu, l - nu,
                               t + mu + nu);
      }
      else
          list_zero (t + mu, nu);
      /* (2) + (3) */
#ifdef TKARMULDEBUG
      print_list (t + mu, nu);
#endif
      list_add (b, b, t + mu, nu);
      list_sub (b + nu, t + mu, b + nu, n - nu + 1);
#ifdef TKARMULDEBUG
      show_result(a, m, c, l, b, n);
#endif
      return tot_muls;
}

/* Computes the space needed for TKarMul of b[0..n],
 * a[0..m] and c[0..l]
 */

unsigned int
TKarMul_space (unsigned int n, unsigned int m, unsigned int l)
{
  unsigned int mu, nu, h;
  unsigned int s1;

  unsigned int r1, r2;
  
  if (n == 0)
      return 0;

  if (m == 0)
      return 0;

  mu = (m / 2) + 1;
  nu = (n / 2) + 1;
  h = MAX (mu, nu);

  if (mu > n)
    {
      r1 = TKarMul_space (n, mu - 1, l);
      if (l >= mu)
	{
	  r2 = TKarMul_space (n, m - mu, l - mu) + n + 1;
          r1 = MAX (r1, r2);
	}
      return r1;
    }

  if (nu > m)
    {
      r1 = TKarMul_space (nu - 1, m, MIN (l, m + nu - 1));

      if (l >= nu)
	{
	  r2 = TKarMul_space (n - nu, m,l - nu);
          r1 = MAX (r1, r2);
	}
      return r1;
    }

  mu = nu = h;
  
  s1 = MIN (l + 1, n + mu);
  r1 = TKarMul_space (nu - 1, mu - 1, s1 - 1) + s1;
  if (s1 >= nu + 1) {
    r2 = TKarMul_space (n - nu, m - mu, s1 - nu - 1) + s1;
    r1 = MAX (r1, r2);
  }
  if (l >= nu) {
    r2 = TKarMul_space (nu - 1, mu - 1, l - nu) + mu + nu;
    r1 = MAX (r1, r2);
  }
  return r1;
}

/* Returns the number of multiplication made in TKarMul
 */

unsigned int
muls_tkara (unsigned int n)
{
  unsigned int mu, nu, h;
  unsigned int m = n;
  unsigned tot_muls = 0;
  unsigned int l = n + m - 1;

  
  if (n == 0)
      return MIN(m, l)  + 1;

  if (m == 0)
      return MIN(n, l) + 1;

  mu = (m / 2) + 1;
  nu = (n / 2) + 1;
  h = MAX (mu, nu);

  mu = nu = h;
  
  if (n + 1 == 2 * nu)
      tot_muls += 3 * muls_tkara (nu - 1);
  else
  {
      tot_muls += 2 * muls_tkara (nu - 1);
      tot_muls += muls_tkara (n - nu);
  }
  return tot_muls;
}

#if 0

unsigned int
TToomCookMul (listz_t b, unsigned int n,
	      listz_t a, unsigned int m, listz_t c, unsigned int l, listz_t t)
{
    unsigned int nu, mu;

    nu = n / 3 + 1;
    mu = m / 3 + 1;

    /* First stip unnecessary trailing coefficients of c:
     */

    l = MIN(l, n + m - 1);

    /*  We should not be called with so small arguments, but
     *  treat this cases anyway.
     */

    if (n == 0)
    {
        mpz_mul (b[0], a[0], c[0]);
        for (k = 1; (k <= m) && (k <= l); k++)
            mpz_addmul (b[0], a[k], c[k]);
        return MIN(m, l);
    }

    if (m == 0)
    {
        for (k = 0; (k <= l) && (k <= n); k++)
            mpz_mul (b[k], a[0], c[k]);
        for (k = l + 1; k <= n; k++)
            mpz_set_ui (b[k], 0);
        return MIN(n, l) + 1;
    }

    /* Now the degenerate cases. We want 2 * nu < m.
     * 
     */

    if (m <= 2 * nu)
    {
       tot_muls += TToomCookMul (b, nu - 1, a, m, c, l, t);
       tot_muls += TToomCookMul (b + nu, nu - 1, a, m, c, l, t);
       tot_muls += TToomCookMul (b + 2 * nu, n - 2 * nu, a, m, c, l, t);
       return tot_muls;
    }
                  
    /* Second degenerate case. We want 2 * mu < m.
     */

    if (n <= 2 * mu)
    {
        tot_muls += TToomCookMul (b, n, a, mu - 1, c, l, t);
        tot_muls += TToomCookMul (t, n, a + mu, mu - 1, c, l, t + n + 1);
        list_add (b, b, t, n + 1);
        tot_muls += TToomCookMul (t, n, a + 2 * mu, m - 2 * mu, c, l, 
                                  t + n + 1);
        list_add (b, b, t, n + 1);
        return tot_muls;
    }

    h = MAX(nu, mu);
    nu = mu = h;


#endif
