/* Implements algorithm POLYEVAL.

  Copyright 2003 Laurent Fousse, Alexander Kruppa, and Paul Zimmermann.

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

#include "gmp.h"
#include "ecm.h"

#ifndef MAX
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif

/* algorithm POLYEVAL from section 3.7 of Peter Montgomery's dissertation.
Input: 
   G - an array of k elements of R, G[i], 0 <= i < k
       representing the coefficients of a polynomial G(x) of degree < k
   Tree - the product tree produced by PolyFromRoots
   Tree[0][0..k-1] (degree k/2)
   Tree[1][0..k-1] (degree k/4), ...,
   Tree[lgk-1][0..k-1] (degree 1)
Output: the sequence of values of G(a[i]) are stored in G[i] for 0 <= i < k
Return value: number of multiplications used.
Remark: we need an auxiliary (k+1)-th cell G[k] in G.
The memory used is M(k) = max(3*floor(k/2)+list_mul_mem(floor(k/2)),
                              k+list_mul_mem(ceil(k/2)),
                              floor(k/2) + M(ceil(k/2))).
Since list_mul_mem(k) >= 2*k, the maximum is the 1st.
*/
unsigned int
polyeval (listz_t G, unsigned int k, listz_t *Tree, listz_t T, mpz_t n,
          int verbose, unsigned int sh)
{
  unsigned int l, m, muls;
  listz_t T0;

  if (k == 1)
    return 0;

  T0 = Tree[0] + sh;
  
  m = k / 2;
  l = k - m;

  /* divide G[0]+G[1]*x+...+G[k-1]*x^(k-1) by
            T0[l]+...+T0[k-1]*x^(m-1)+x^m,
            quotient in {T+m,l-1}, remainder in {T,m} */

  if (k == 2 * m)
    {
      /* FIXME: avoid the copy here by giving different 2nd and 3rd arguments
         to RecursiveDivision */
      list_set (T, G, k);
      /* the following needs k+m+list_mul_mem(m) in T */
      muls = RecursiveDivision (T + k, T, T0 + l, m, T + k + m, n, 1);
    }
  else /* k = 2m+1: subtract G[k-1]*x^(l-1) * T0 from G */
    {
      /* G - G[k-1] * (x^m + {T0+l,m}) * x^m */
      list_set (T, G, m);
      list_mul_z (T + m, T0 + l, G[k - 1], m, n);
      muls = m;
      list_sub (T + m, G + m, T + m, m);
      /* the following needs 3m+list_mul_mem(m) in T */
      muls += RecursiveDivision (T + 2 * m, T, T0 + l, m, T + 3 * m, n, 1);
    }
  /* in both cases we need 3*(k/2)+list_mul_mem(k/2) */

  /* right remainder is in {T,m} */

  /* k = 2l or k = 2l-1 */
  
  /* divide G[0]+G[1]*x+...+G[k-1]*x^(k-1) by
            T0[0]+...+T0[l-1]*x^(l-1)+x^l:
            quotient in {T+m,m-1}, remainder in {G,l} */

  if (k < 2 * l)
    mpz_set_ui (G[k], 0);
  /* the following needs k+list_mul_mem(l) in T */
  muls += RecursiveDivision (T + m, G, T0, l, T + k, n, 1);

  /* left remainder is in {G,l} */
  
  muls += polyeval (G, l, Tree + 1, T + m, n, verbose, sh);

  /* copy right remainder in {G+l,m} */
  list_set (G + l, T, m);
  muls += polyeval (G + l, m, Tree + 1, T, n, verbose, sh + l);

  return muls;
}

#ifdef TUPTREE_DEBUG

void print_vect (listz_t t, unsigned int l)
{
    unsigned int i;

    printf ("[");
    for (i = 0; i < l; i++)
    {
        mpz_out_str (NULL, 10, t[i]);
        if (i != l - 1)
            printf (", ");
        else
            printf ("]");
    }
}

#endif

/* Computes TUpTree as described in ref[1]. k is the degree of the
 * polynomial at the root of the tree. sh is the shift we need to
 * apply to find the actual coefficients of the polynomial at the root
 * of the tree.
 */

unsigned int TUpTree (listz_t b, listz_t *Tree, unsigned int k,
              listz_t tmp, unsigned int sh, mpz_t n)
{

    unsigned int m, l, i;
    unsigned int tot_muls = 0;
    unsigned int muls;

    m = k / 2;
    l = k - m;
    
    if (k == 1) {
        mpz_mod (b[0], b[0], n);
        return 0;
    }
   
#ifdef TUPTREE_DEBUG
    printf ("Dans TupTree, k = %d.\n", k);

    printf ("b = ");
    print_vect (b, k);
    printf ("\nLes poly de ce niveau sont : ");
    print_vect (Tree[0] + sh, k);
    printf ("\n");
#endif

    muls = TMulGen (tmp, l - 1, Tree[0] + sh + l, m - 1, b, k - 1, tmp + l);
#ifdef CHECK_MULS
    if (muls > muls_tgen (l-1)) { printf ("l=%u m=%u k=%u %u %u\n", l, m, k, muls, muls_tgen (l-1)); abort(); }
#endif
    tot_muls += muls;
    muls = TMulGen (tmp + l, m - 1, Tree[0] + sh, l - 1, b, k - 1, tmp + k);
#ifdef CHECK_MULS
    if (muls > muls_tgen (l-1)) { printf ("m=%u l=%u k=%u %u %u\n", m, l, k, muls, muls_tgen (m-1)); abort();}
#endif
    tot_muls += muls;

#ifdef TUPTREE_DEBUG
    printf ("Et le résultat à ce niveau (avant correction) est : ");
    print_vect (tmp, k);
    printf ("\n");
#endif

    /* gmp-ecm specific: leading coefficients in the product tree
     * are ones and are implied, so we need some extra work here.
     */

    for (i = 0; i < l; i++)
    {
        mpz_add (tmp[i], tmp[i], b[m + i]);
    }

    for (i = 0; i < m; i++)
    {
        mpz_add (tmp[l + i], tmp[l + i], b[l + i]);
    }

    for (i = 0; i < l; i++)
    {
        mpz_mod (b[i], tmp[i], n);
    }

    for (i = 0; i < m; i++)
    {
        mpz_mod (b[i + l], tmp[l + i], n);
    }


#ifdef TUPTREE_DEBUG
    printf ("Et le résultat à ce niveau est : ");
    print_vect (b, k);
    printf ("\n");
#endif
    
    tot_muls += TUpTree (b, Tree + 1, l, tmp, sh, n);
    tot_muls += TUpTree (b + l, Tree + 1, m, tmp, sh + l, n);
    return tot_muls;
}

unsigned int TUpTree_space (unsigned int k)
{

    unsigned int m, l;
    unsigned int r1, r2;

    m = k / 2;
    l = k - m;
    
    if (k == 1) {
        return 0;
    }
   
    r1 = TMulGen_space (l - 1, m - 1, k - 1) + l;
    r2 = TMulGen_space (m - 1, l - 1, k - 1) + k;

    r1 = MAX (r1, r2);



    r2 = TUpTree_space (l);
    r1 = MAX (r1, r2);
    r2 = TUpTree_space (m);
    r1 = MAX (r1, r2);
    return r1;
}


/* Same as polyeval. Needs invF as extra argument.
 */

unsigned int
polyeval_tellegen (listz_t b, unsigned int k, listz_t *Tree, listz_t tmp,
                   unsigned int sizeT, listz_t invF, mpz_t n, unsigned int sh)
{
    unsigned int i;
    unsigned int tupspace;
    unsigned int tkspace;
    unsigned int totmuls = 0, muls;
    int allocated = 0;
    listz_t T;

    tupspace = TUpTree_space (k) + k;
#ifndef USE_SHORT_PRODUCT
    tkspace = TMulGen_space (k - 1, k - 1, k - 1) + k;
#else
    tkspace = 2 * k - 1 + list_mul_mem (k);
#endif

    tupspace = MAX (tupspace, tkspace);

    if (sizeT >= tupspace)
        T = tmp;
    else
    {
        T = init_list (tupspace);
        allocated = 1;
    }
    
#ifdef TELLEGEN_DEBUG
    printf ("Dans polyeval_tellegen, k = %d.\n", k);
    printf ("Espace requis : %d.\n", 
            TMulGen_space (k - 1, k - 1, k - 1));
#endif
#ifndef USE_SHORT_PRODUCT
    /* revert invF for call to TMulGen below */
    for (i = 0; i < k / 2; i++)
        mpz_swap (invF[i], invF[k - 1 - i]);
    muls = TMulGen (T, k - 1, invF, k - 1, b, k - 1, T + k);
#ifdef CHECK_MULS
    if (muls != muls_tgen (k-1)) { printf ("%u %u %u\n", k, muls, muls_tgen (k-1)); abort();}
#endif
#else
    /* need space 2k-1+list_mul_mem(k) in T */
    muls = list_mul_high (T, invF, b, k, T + 2 * k - 1);
    list_set (T, T + k - 1, k);
#endif
    totmuls += muls;
#ifdef TELLEGEN_DEBUG
    printf ("\nalpha = ");
    print_list (invF, k);
    printf ("\nt = ");
    print_list (T, k);
#endif
    for (i = 0; i < k / 2; i++)
        mpz_swap (T[i], T[k - 1 - i]);
#ifdef TELLEGEN_DEBUG
    printf ("s = ");
    print_list (T, k);
#endif
    totmuls += TUpTree (T, Tree, k, T + k, sh, n);
    for (i = 0; i < k; i++)
    {
        mpz_set (b[i], T[i]);
    }

    if (allocated)
        clear_list (T, tupspace);
    return totmuls;
}

unsigned int muls_tuptree (unsigned int k)
{
    unsigned int m, l;
    unsigned int tot_muls = 0;

    m = k / 2;
    l = k - m;
    
    if (k == 1)
        return 0;

    /* bound both calls to TMulGen(l-1,m-1) and TMulGen(m-1,l-1)
       by TMulGen(l-1,l-1) */
    tot_muls += 2 * muls_tgen (l - 1);

    if (l == m)
      {
        tot_muls += 2 * muls_tuptree (l);
      }
    else
      {
        tot_muls += muls_tuptree (m);
        tot_muls += muls_tuptree (l);
      }
    return tot_muls;
}

unsigned int
muls_polyeval_tellegen (unsigned int k)
{
    unsigned int tot_muls = 0;

#ifndef USE_SHORT_PRODUCT
    tot_muls += muls_tgen (k - 1);
#else
    tot_muls += muls_gen_short (k);
#endif
    tot_muls += muls_tuptree (k);
    return tot_muls;
}
