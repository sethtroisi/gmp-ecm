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

#if !defined (_MSC_VER)
#include <sys/time.h>
#include <unistd.h>
#if !defined (__MINGW32__)
#include <sys/resource.h>
#endif
#endif

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
  printf ("Entr�e dans TKarMul.\nm = %d\nn = %d\nl = %d\ndepth = %d\n", m, n,
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

/* return number of muls when multiplying two polynomials of degree n */
unsigned int
muls_tgen (unsigned int n)
{
  return muls_toom3 (n + 1); /* put muls_toom3 in hard here, since currently
                                TMulGen directly calls TToomCookMul; put n+1
                                since n is number of terms in muls_toom3 */
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


/* list_sub with bound checking
 */

void list_sub_safe (listz_t ret, listz_t a, listz_t b,
                        unsigned int sizea, unsigned int sizeb,
                        unsigned int needed)
{
    unsigned int i;
    unsigned int safe;
    safe = MIN(sizea, sizeb);
    safe = MIN(safe, needed);

    for (i = 0; i < safe; i++)
        mpz_sub (ret[i], a[i], b[i]);

    while (i < needed)
    {
        if (i < sizea)
        {
            if (i < sizeb)
                mpz_sub (ret[i], a[i], b[i]);
            else
                mpz_set (ret[i], a[i]);
        }
        else
        {
            if (i < sizeb)
                mpz_neg (ret[i], b[i]);
            else
                mpz_set_ui (ret[i], 0);
        }
        i++;
    }
}

/* list_add with bound checking
 */

void list_add_safe (listz_t ret, listz_t a, listz_t b,
                        unsigned int sizea, unsigned int sizeb,
                        unsigned int needed)
{
    unsigned int i;
    unsigned int safe;
    safe = MIN(sizea, sizeb);
    safe = MIN(safe, needed);

    for (i = 0; i < safe; i++)
        mpz_add (ret[i], a[i], b[i]);

    while (i < needed)
    {
        if (i < sizea)
        {
            if (i < sizeb)
                mpz_add (ret[i], a[i], b[i]);
            else
                mpz_set (ret[i], a[i]);
        }
        else
        {
            if (i < sizeb)
                mpz_set (ret[i], b[i]);
            else
                mpz_set_ui (ret[i], 0);
        }
        i++;
    }
}

unsigned int
TToomCookMul (listz_t b, unsigned int n,
              listz_t a, unsigned int m, listz_t c, unsigned int l, 
              listz_t tmp)
{
    unsigned int nu, mu, h;
    unsigned int i;
    unsigned int btmp;
    mpz_t comp_tmp;
    unsigned int tot_muls = 0;

    nu = n / 3 + 1;
    mu = m / 3 + 1;

    /* ensures n + 1 > 2 * nu */
    if ((n < 2 * nu) || (m < 2 * mu))
    {
#ifdef TTCDEBUG
        printf ("Op�randes trop petites, on appelle TKara.\n");
#endif
        return TKarMul (b, n, a, m, c, l, tmp);
    }

    /* First strip unnecessary trailing coefficients of c:
     */

    l = MIN(l, n + m);

    /*  We should not be called with so small arguments, but
     *  treat this cases anyway.
     */

    if (n == 0)
    {
        mpz_mul (b[0], a[0], c[0]);
        for (i = 1; (i <= m) && (i <= l); i++)
            mpz_addmul (b[0], a[i], c[i]);
        return MIN(m, l);
    }

    if (m == 0)
    {
        for (i = 0; (i <= l) && (i <= n); i++)
            mpz_mul (b[i], a[0], c[i]);
        for (i = l + 1; i <= n; i++)
            mpz_set_ui (b[i], 0);
        return MIN(n, l) + 1;
    }

    /* Now the degenerate cases. We want 2 * nu <= m.
     * 
     */

    if (m < 2 * nu)
    {
#ifdef TTCDEBUG
        printf ("Cas d�g�n�r� 1.\n");
#endif
        tot_muls += TToomCookMul (b, nu - 1, a, m, c, l, tmp);
        if (l >= nu)
            tot_muls += TToomCookMul (b + nu, nu - 1, a, m, 
                                      c + nu, l - nu, tmp);
        else
            list_zero (b + nu, nu);
        if (l >= 2 * nu) /* n >= 2 * nu is assured. Hopefully */
            tot_muls += TToomCookMul (b + 2 * nu, n - 2 * nu, a, m, 
                                      c + 2 * nu, l - 2 * nu, tmp);
        else
            list_zero (b + 2 * nu, n - 2 * nu + 1);
        return tot_muls;
    }
                  
    /* Second degenerate case. We want 2 * mu <= n.
     */

    if (n < 2 * mu)
    {
#ifdef TTCDEBUG
        printf ("Cas d�g�n�r� 2.\n");
#endif
        tot_muls += TToomCookMul (b, n, a, mu - 1, c, l, tmp);
        if (l >= mu)
        {
            tot_muls += TToomCookMul (tmp, n, a + mu, mu - 1, 
                                      c + mu, l - mu, tmp + n + 1);
            list_add (b, b, tmp, n + 1);
        }
        if (l >= 2 * mu)
        {
            tot_muls += TToomCookMul (tmp, n, a + 2 * mu, m - 2 * mu, 
                                      c + 2 * mu, l - 2 * mu, tmp + n + 1);
            list_add (b, b, tmp, n + 1);
        }
        return tot_muls;
    }

#ifdef TTCDEBUG
    printf ("Cas de base.\n");
    printf ("a = ");
    print_list (a, m + 1);

    printf ("\nc = ");
    print_list (c, l + 1);
#endif
    h = MAX(nu, mu);
    nu = mu = h;

    mpz_init (comp_tmp);

    list_sub_safe (tmp, c + 3 * h, c + h,
                   (l + 1 > 3 * h ? l + 1 - 3 * h : 0), 
                   (l + 1 > h ? l + 1 - h : 0), 2 * h - 1);
    list_sub_safe (tmp + 2 * h - 1, c, c + 2 * h,
                   l + 1, (l + 1 > 2 * h ? l + 1 - 2 * h : 0),
                   2 * h - 1);
    for (i = 0; i < 2 * h - 1; i++)
        mpz_mul_ui (tmp[2 * h - 1 + i], tmp[2 * h - 1 + i], 2);
    
#ifdef TTCDEBUG
    print_list (tmp, 4 * h - 2);
#endif

    /* --------------------------------
     * | 0 ..  2*h-2 | 2*h-1 .. 4*h-3 |
     * --------------------------------
     * | c3 - c1     |   2(c0 - c2)   |
     * --------------------------------
     */

    list_add (tmp + 2 * h - 1, tmp + 2 * h - 1, tmp, 2 * h - 1);

    tot_muls += TToomCookMul (b, h - 1, a, h - 1, tmp + 2 * h - 1, 
                              2 * h - 2, tmp + 4 * h - 2);

    /* b[0 .. h - 1] = 2 * m0 */

#ifdef TTCDEBUG
    printf ("2 * m0 = ");
    print_list (b, h);
#endif

    for (i = 0; i < h; i++)
        mpz_add (tmp[2 * h - 1 + i], a[i], a[h + i]);
    
    for (i = 0; i < MIN(h, m + 1 - 2 * h); i++)
        mpz_add (tmp[2 * h - 1 + i], tmp[2 * h - 1 + i], a[2 * h + i]);

    /* tmp[2*h-1 .. 3*h-2] = a0 + a1 + a2 */

#ifdef TTCDEBUG
    printf ("\na0 + a1 + a2 = ");
    print_list (tmp + 2 * h - 1, h);
#endif

    list_sub_safe (tmp + 3 * h - 1, c + 2 * h, c + 3 * h, 
                   (l + 1 > 2 * h ? l + 1 - 2 * h : 0),
                   (l + 1 > 3 * h ? l + 1 - 3 * h : 0),
                   2 * h - 1);

    /* -------------------------------------------------
     * | 0 ..  2*h-2 | 2*h-1 .. 3*h-2 | 3*h-1 .. 5*h-3 |
     * -------------------------------------------------
     * | c3 - c1     |  a0 + a1 + a2  |   c2 - c3      |
     * -------------------------------------------------
     */

    btmp = (l + 1 > h ? l + 1 - h : 0);
    btmp = MIN(btmp, 2 * h - 1);
    for (i = 0; i < btmp; i++)
    {
        mpz_mul_ui (comp_tmp, c[h + i], 2);
        mpz_add (tmp[5 * h - 2 + i], comp_tmp, tmp[3 * h - 1 + i]);
    }
    while (i < 2 * h - 1)
    {
        mpz_set (tmp[5 * h - 2 + i], tmp[3 * h - 1 + i]);
        i++;
    }

    tot_muls += TToomCookMul (b + h, h - 1, tmp + 2 * h - 1, h - 1, 
                              tmp + 5 * h - 2, 2 * h - 2,
                              tmp + 7 * h - 3);

    /* b[h .. 2 * h - 1] = 2 * m1 */
#ifdef TTCDEBUG
    printf ("\n2 * m1 = ");
    print_list (b + h, h);
#endif

    /* ------------------------------------------------------------------
     * | 0 ..  2*h-2 | 2*h-1 .. 3*h-2 | 3*h-1 .. 5*h-3 | 5*h-2 .. 7*h-4 |
     * ------------------------------------------------------------------
     * | c3 - c1     |  a0 + a1 + a2  |   c2 - c3      | c2 - c3 + 2c1  |
     * ------------------------------------------------------------------
     */


    for (i = 0; i < h; i++)
    {
        mpz_add (tmp[2 * h  - 1 + i], tmp[2 * h  - 1 + i], a[i + h]);
        if (2 * h + i <= m)
        {
            mpz_mul_ui (comp_tmp, a[2 * h + i], 3);
            mpz_add (tmp[2 * h  - 1 + i], tmp[2 * h - 1 + i], comp_tmp);
        }
    }
    tot_muls += TToomCookMul (tmp + 5 * h - 2, h - 1, 
                              tmp + 2 * h - 1, h - 1,
                              tmp, 2 * h - 2, tmp + 6 * h - 2);

    /* tmp[5*h-2 .. 6*h - 3] = 6 * m2  */ 
    
#ifdef TTCDEBUG
    printf ("\n6 * m2 = ");
    print_list (tmp + 5 * h - 2, h);
#endif
    for (i = 0; i < h; i++)
    {
        mpz_sub (tmp[2 * h - 1 + i], a[i], a[h + i]);
        if (i + 2 * h <= m)
            mpz_add (tmp[2 * h - 1 + i], tmp[2 * h - 1 + i], a[2 * h + i]);
    }

    for (i = 0; i < 2 * h - 1; i++)
    {
        mpz_mul_ui (tmp[3 * h - 1 + i], tmp[3 * h - 1 + i], 3);
        mpz_mul_ui (tmp[i], tmp[i], 2);
    }

    list_add (tmp + 3 * h - 1, tmp + 3 * h - 1, tmp, 2 * h - 1);

    tot_muls += TToomCookMul (tmp + 6 * h - 2, h - 1,
                              tmp + 2 * h - 1, h - 1,
                              tmp + 3 * h - 1, 2 * h - 2, 
                              tmp + 7 * h - 2);

    /* tmp[6h-2 .. 7h - 3] = 6 * mm1 */

#ifdef TTCDEBUG
    printf ("\n6 * mm1 = ");
    print_list (tmp + 6 * h - 2, h);
#endif
    list_add_safe (tmp, tmp, c + 2 * h,
                   2 * h,
                   (l + 1 > 2 * h ? l + 1 - 2 * h : 0),
                   2 * h - 1);

    list_sub_safe (tmp, c + 4 * h, tmp,
                   (l + 1 > 4 * h ? l + 1 - 4 * h : 0),
                   2 * h - 1, 2 * h - 1);

    tot_muls += TToomCookMul (b + 2 * h, n - 2 * h, a + 2 * h, m - 2 * h,
                  tmp, 2 * h - 1, tmp + 7 * h - 2);

    /* b[2 * h .. n] = minf */

#ifdef TTCDEBUG
    printf ("\nminf = ");
    print_list (b + 2 * h, n + 1 - 2 * h);
#endif

    /* Layout of b : 
     * ---------------------------------------
     * | 0 ... h-1 | h ... 2*h-1 | 2*h ... n |
     * ---------------------------------------
     * |  2 * m0   |   2 * m1    |    minf   |
     * ---------------------------------------
     * 
     * Layout of tmp :
     * ---------------------------------------------------
     * | 0 ... 5*h-1 | 5*h-2 ... 6*h-3 | 6*h-2 ... 7*h-3 |
     * ---------------------------------------------------
     * |  ??????     |    6 * m2       |   6 * mm1       |
     * ---------------------------------------------------
     */
    
    list_add (tmp, tmp + 5 * h - 2, tmp + 6 * h - 2, h);
    for (i = 0; i < h; i++)
        mpz_divby3_1op (tmp[i]);

    /* t1 = 2 (m2 + mm1)
     * tmp[0 .. h - 1] = t1
     */
    
    list_add (b, b, b + h, h);
    list_add (b, b, tmp, h);
    for (i = 0; i < h; i++)
        mpz_tdiv_q_2exp (b[i], b[i], 1);

    /* b_{low} should be correct */

    list_add (tmp + h, b + h, tmp, h);

    /* t2 = t1 + 2 m1
     * tmp[h .. 2h - 1] = t2
     */

    list_add (b + h, tmp, tmp + h, h);
    list_sub (b + h, b + h, tmp + 6 * h - 2, h);
    for (i = 0; i < h; i++)
        mpz_tdiv_q_2exp (b[h + i], b[h + i], 1);

    /* b_{mid} should be correct */

    list_add (tmp + h, tmp + h, tmp + 5 * h - 2, n + 1 - 2 * h);
    for (i = 0; i < n + 1 - 2 * h; i++)
        mpz_tdiv_q_2exp (tmp[h + i], tmp[h + i], 1);

    list_add (b + 2 * h, b + 2 * h, tmp + h, n + 1 - 2 * h);
    /* b_{high} should be correct */

    mpz_clear (comp_tmp);
    return tot_muls;
}

/* Returns space needed by TToomCookMul */

unsigned int
TToomCookMul_space (unsigned int n, unsigned int m, unsigned int l)
              
{
    unsigned int nu, mu, h;
    unsigned int stmp1, stmp2;

    nu = n / 3 + 1;
    mu = m / 3 + 1;

    stmp1 = stmp2 = 0;

    /* ensures n + 1 > 2 * nu */
    if ((n < 2 * nu) || (m < 2 * mu))
    {
        return TKarMul_space (n, m, l);
    }

    /* First strip unnecessary trailing coefficients of c:
     */

    l = MIN(l, n + m);

    /*  We should not be called with so small arguments, but
     *  treat this cases anyway.
     */

    if (n == 0)
        return 0;

    if (m == 0)
        return 0;

    /* Now the degenerate cases. We want 2 * nu < m.
     * 
     */

    if (m <= 2 * nu)
    {
        stmp1 = TToomCookMul_space (nu - 1, m, l);
        if (l >= nu)
            stmp2 = TToomCookMul_space (nu - 1, m, l - nu);
        stmp1 = MAX(stmp1, stmp2);
        if (l >= 2 * nu)
            stmp2 = TToomCookMul_space (n - 2 * nu, m, l - 2 * nu);
        stmp1 = MAX(stmp1, stmp2);
        return stmp1;
    }
                  
    /* Second degenerate case. We want 2 * mu < n.
     */

    if (n <= 2 * mu)
    {
        stmp1 += TToomCookMul_space (n, mu - 1, l);
        if (l >= mu)
        {
            stmp2 = TToomCookMul_space (n, mu - 1, l - mu) + n + 1;
            stmp1 = MAX(stmp1, stmp2);
        }
        if (l >= 2 * mu)
        {
            stmp2 = TToomCookMul_space (n, m - 2 * mu, l - 2 * mu) + n + 1;
            stmp1 = MAX(stmp1, stmp2);
        }
        return stmp1;
    }

    h = MAX(nu, mu);

    stmp2 = TToomCookMul_space (h - 1, h - 1, 2 * h - 2) + 7 * h - 2;
    stmp1 = TToomCookMul_space (h - 1, h - 1, 2 * h - 2) + 6 * h - 2;
    stmp1 = MAX(stmp1, stmp2);
    stmp2 = TToomCookMul_space (n - 2 * h, m - 2 * h, 2 * h - 1) + 7*h-2;
    return MAX(stmp1, stmp2);
}

unsigned int
TMulGen (listz_t b, unsigned int n,
              listz_t a, unsigned int m, listz_t c, unsigned int l, 
              listz_t tmp)
{
    return TToomCookMul (b, n, a, m, c, l, tmp);
}



unsigned int
TMulGen_space (unsigned int n, unsigned int m, unsigned int l)
{
    return TToomCookMul_space (n, m, l);
}
