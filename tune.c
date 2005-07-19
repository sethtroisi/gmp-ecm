/* Tune program.

  Copyright 2003, 2005 Paul Zimmermann, Alexander Kruppa, Dave Newman.

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
#include <gmp.h>
#include "ecm.h"
#include "ecm-gmp.h"
#include "ecm-impl.h"

/* we don't need any more precision */
#define GRANULARITY 1e2
#define MPZSPV_MAX_LOG2_LEN 18 /* 2 * 131072 */
#define MPZSPV_MAX_LEN (1 << MPZSPV_MAX_LOG2_LEN)

/* Throughout, each function pointer points to a function
 * 
 *   double f0 (size_t limbs, unsigned int mintime);
 *
 * that runs for at least mintime ms and then returns the number of iterations
 * performed per ms. */

#define M_str "95209938255048826235189575712705128366296557149606415206280987204268594538412191641776798249266895999715600261737863698825644292938050707507901970225804581"

mpz_t M; /* yes, global variables */
gmp_randstate_t gmp_randstate;

#undef MPZMOD_THRESHOLD
#undef REDC_THRESHOLD
extern size_t MPZMOD_THRESHOLD;
extern size_t REDC_THRESHOLD;

#ifdef HAVE_NTT
mpzspm_t mpzspm;
mpzv_t x, y, z, t;
spm_t spm;
spv_t spv;
mpzspv_t mpzspv;

#undef SPV_NTT_GFP_DIF_RECURSIVE_THRESHOLD
#undef SPV_NTT_GFP_DIT_RECURSIVE_THRESHOLD
#undef MUL_NTT_THRESHOLD
#undef PREREVERTDIVISION_NTT_THRESHOLD
#undef POLYINVERT_NTT_THRESHOLD
#undef POLYEVALT_NTT_THRESHOLD
#undef MPZSPV_NORMALISE_STRIDE
extern size_t SPV_NTT_GFP_DIF_RECURSIVE_THRESHOLD;
extern size_t SPV_NTT_GFP_DIT_RECURSIVE_THRESHOLD;
extern size_t MUL_NTT_THRESHOLD;
extern size_t PREREVERTDIVISION_NTT_THRESHOLD;
extern size_t POLYINVERT_NTT_THRESHOLD;
extern size_t POLYEVALT_NTT_THRESHOLD;
extern size_t MPZSPV_NORMALISE_STRIDE;
#endif



double
tune_mpres_mul (size_t limbs, unsigned int mintime, int repr)
{
  mpmod_t modulus;
  mpres_t x, y, z;
  mpz_t N, p, q;
  unsigned int st, k = 0;

  mpz_init (N);
  mpz_init (p);
  mpz_init (q);
  
  /* No need to generate a probable prime, just ensure N is not
     divisible by 2 or 3 */
  do
    {
      mpz_random (N, limbs);
      while (mpz_gcd_ui (NULL, N, 6) != 1)
        mpz_add_ui (N, N, 1);
    }
  while ((mp_size_t) mpz_size (N) != limbs);
  
  switch (repr)
  {
    case ECM_MOD_MPZ:
      mpmod_init_MPZ (modulus, N);
      break;
    case ECM_MOD_MODMULN:
      mpmod_init_MODMULN (modulus, N);
      break;
    case ECM_MOD_REDC:
      mpmod_init_REDC (modulus, N);
      break;
  }

  mpz_urandomm (p, gmp_randstate, N);
  mpz_urandomm (q, gmp_randstate, N);
  
  mpres_init (x, modulus);
  mpres_init (y, modulus);
  mpres_init (z, modulus);

  mpres_set_z (x, p, modulus);
  mpres_set_z (y, q, modulus);

  st = cputime ();

  do
    {
      mpres_mul (z, x, y, modulus);
      k++;
    }
  while (cputime () - st < mintime);

  st = cputime () - st;

  mpres_clear (x, modulus);
  mpres_clear (y, modulus);
  mpres_clear (z, modulus);
  mpmod_clear (modulus);
  mpz_clear (N);
  mpz_clear (p);
  mpz_clear (q);

  return (double) k / (double) st;
}

double
tune_mpres_mul_mpz (size_t n, unsigned int mintime)
{
  return tune_mpres_mul (n, mintime, ECM_MOD_MPZ);
}

double
tune_mpres_mul_modmuln (size_t n, unsigned int mintime)
{
  return tune_mpres_mul (n, mintime, ECM_MOD_MODMULN);
}

double
tune_mpres_mul_redc (size_t n, unsigned int mintime)
{
  return tune_mpres_mul (n, mintime, ECM_MOD_REDC);
}

#ifdef HAVE_NTT
double
tune_spv_ntt_gfp_dif_recursive (size_t log2_len, unsigned int mintime)
{
  unsigned int k = 0, st = cputime ();
  SPV_NTT_GFP_DIF_RECURSIVE_THRESHOLD = 1 << log2_len;

  do
    {
      spv_ntt_gfp_dif (spv, 1 << log2_len, spm->sp, spm->mul_c,
	  spm->prim_root);
      k++;
    }
  while (cputime () - st < mintime);
  
  return (double) k / (double) (cputime () - st);
}

double
tune_spv_ntt_gfp_dif_unrolled (size_t log2_len, unsigned int mintime)
{
  unsigned int k = 0, st = cputime ();
  SPV_NTT_GFP_DIF_RECURSIVE_THRESHOLD = ULONG_MAX;

  do
    {
      spv_ntt_gfp_dif (spv, 1 << log2_len, spm->sp, spm->mul_c,
	  spm->prim_root);
      k++;
    }
  while (cputime () - st < mintime);
  
  return (double) k / (double) (cputime () - st);
}


double
tune_spv_ntt_gfp_dit_recursive (size_t log2_len, unsigned int mintime)
{
  unsigned int k = 0, st = cputime ();
  SPV_NTT_GFP_DIT_RECURSIVE_THRESHOLD = 1 << log2_len;

  do
    {
      spv_ntt_gfp_dit (spv, 1 << log2_len, spm->sp, spm->mul_c,
	  spm->prim_root);
      k++;
    }
  while (cputime () - st < mintime);
  
  return (double) k / (double) (cputime () - st);
}

double
tune_spv_ntt_gfp_dit_unrolled (size_t log2_len, unsigned int mintime)
{
  unsigned int k = 0, st = cputime ();
  SPV_NTT_GFP_DIT_RECURSIVE_THRESHOLD = ULONG_MAX;

  do
    {
      spv_ntt_gfp_dit (spv, 1 << log2_len, spm->sp, spm->mul_c,
	  spm->prim_root);
      k++;
    }
  while (cputime () - st < mintime);
  
  return (double) k / (double) (cputime () - st);
}


double
tune_ntt_mul (size_t log2_len, unsigned int mintime)
{
  unsigned int k = 0, st = cputime ();
  MUL_NTT_THRESHOLD = 0;

  do
    {
      ntt_mul (z, x, y, 1 << log2_len, NULL, 1, mpzspm);
      k++;
    }
  while (cputime () - st < mintime);

  return (double) k / (double) (cputime () - st);
}

double
tune_list_mul (size_t log2_len, unsigned int mintime)
{
  unsigned int k = 0, st = cputime ();

  do
    {
      list_mul (z, x, 1 << log2_len, 1, y, 1 << log2_len, 1, t);
      k++;
    }
  while (cputime () - st < mintime);

  return (double) k / (double) (cputime () - st);
}


double
tune_ntt_PrerevertDivision (size_t log2_len, unsigned int mintime)
{
  unsigned int k = 0, st = cputime ();

  PREREVERTDIVISION_NTT_THRESHOLD = 0;

  do
    {
      ntt_PrerevertDivision (z, x, y, mpzspv, 1 << log2_len, t, mpzspm);
      k++;
    }
  while (cputime () - st < mintime);

  return (double) k / (double) (cputime () - st);
}

double
tune_PrerevertDivision (size_t log2_len, unsigned int mintime)
{
  unsigned int k = 0, st = cputime ();

  do
    {
      PrerevertDivision (z, x, y, 1 << log2_len, t, mpzspm->modulus);
      k++;
    }
  while (cputime () - st < mintime);

  return (double) k / (double) (cputime () - st);
}

double
tune_ntt_PolyInvert (size_t log2_len, unsigned int mintime)
{
  unsigned int k = 0, st = cputime ();

  POLYINVERT_NTT_THRESHOLD = 1 << log2_len;
  
  do
    {
      ntt_PolyInvert (z, x, 1 << log2_len, t, mpzspm);
      k++;
    }
  while (cputime () - st < mintime);

  return (double) k / (double) (cputime () - st);
}

double
tune_PolyInvert (size_t log2_len, unsigned int mintime)
{
  unsigned int k = 0, st = cputime ();

  do
    {
      PolyInvert (z, x, 1 << log2_len, t, mpzspm->modulus);
      k++;
    }
  while (cputime () - st < mintime);

  return (double) k / (double) (cputime () - st);
}
  

double
tune_ntt_polyevalT (size_t log2_len, unsigned int mintime)
{
  unsigned int i, k = 0, st;
  mpzv_t *Tree = (mpzv_t *) malloc ((log2_len + 1) * sizeof (mpzv_t));
  
  for (i = 0; i <= log2_len; i++)
    Tree[i] = x;

  POLYEVALT_NTT_THRESHOLD = 1 << log2_len;

  st = cputime ();
  
  do
    {
      ntt_polyevalT (z, 1 << log2_len, Tree, t, mpzspv, mpzspm, NULL);
      k++;
    }
  while (cputime () - st < mintime);

  free (Tree);
  
  return (double) k / (double) (cputime () - st);
}

double
tune_polyevalT (size_t log2_len, unsigned int mintime)
{
  unsigned int i, k = 0, st;
  mpzv_t *Tree = (mpzv_t *) malloc ((log2_len + 1) * sizeof (mpzv_t));

  for (i = 0; i <= log2_len; i++)
    Tree[i] = x;

  st = cputime ();
  
  do
    {
      polyeval_tellegen (z, 1 << log2_len, Tree, t, 3 * (1 << log2_len),
	  x, mpzspm->modulus, NULL);
      k++;
    }
  while (cputime () - st < mintime);

  free (Tree);
  
  return (double) k / (double) (cputime () - st);
}

double
tune_mpzspv_normalise (size_t log2_len, unsigned int mintime)
{
  unsigned int k = 0, st = cputime ();

  MPZSPV_NORMALISE_STRIDE = 1 << log2_len;
  
  do
    {
      mpzspv_normalise (mpzspv, 0, MPZSPV_MAX_LEN, mpzspm);
      k++;
    }
  while (cputime () - st < mintime);

  return (double) k / (double) (cputime () - st);
} 
#endif


/* Assume f0 and f1 are monotone decreasing. Return the first n in the range
 * [min_n, max_n) for which f1(n) > f0(n), or return max_n if no such n
 * exists. */
size_t
crossover (double (*f0)(size_t, unsigned int),
    double (*f1)(size_t, unsigned int), size_t min_n, size_t max_n)
{
  size_t mid_n;
  
  if (min_n == max_n)
    return min_n;

  mid_n = (max_n + min_n) / 2;
  return ((f0)(mid_n, GRANULARITY) >= (f1)(mid_n, GRANULARITY))
    ? crossover (f0, f1, mid_n + 1, max_n)
    : crossover (f0, f1, min_n, mid_n);
}

/* Return the lowest n with min_n <= n < max_n such that
 * f1(t) > f0(t) for all t in [n, n + k)
 *
 * Return max_n if no such n exists. */
size_t
crossover2 (double (*f0)(size_t, unsigned int),
    double (*f1)(size_t, unsigned int), size_t min_n, size_t max_n, size_t k)
{
  size_t n = min_n;
  size_t t;
  
  while (n < max_n)
    {
      for (t = n + k; t > n; t--)
	if ((f0)(t - 1, GRANULARITY) >= (f1)(t - 1, GRANULARITY))
	  break;

      if (t == n)
	return n;

      n = t;
    };

  return max_n;
}

/* Return the n in the range [min_n, max_n) that maximises f(n) */
size_t
maximise (double (*f)(size_t, unsigned int), size_t min_n, size_t max_n)
{
  size_t n, best_n = 0;
  double f_n, f_best_n = -1.0;

  for (n = min_n; n < max_n; n++)
    {
      f_n = f (n, GRANULARITY);
      if (f_n > f_best_n)
        {
	  f_best_n = f_n;
	  best_n = n;
	}
    }

  return best_n;
}

int main ()
{
  spv_size_t i;

  gmp_randinit_default (gmp_randstate);

#ifdef HAVE_NTT
  mpz_init_set_str (M, M_str, 10);
  
  x = init_list (MPZSPV_MAX_LEN);
  y = init_list (MPZSPV_MAX_LEN);
  z = init_list (MPZSPV_MAX_LEN);
  t = init_list (list_mul_mem (MPZSPV_MAX_LEN / 2) + 3 * MPZSPV_MAX_LEN / 2);
  
  for (i = 0; i < MPZSPV_MAX_LEN; i++)
    mpz_urandomm (x[i], gmp_randstate, M);
  for (i = 0; i < MPZSPV_MAX_LEN; i++)
    mpz_urandomm (y[i], gmp_randstate, M); 
  for (i = 0; i < MPZSPV_MAX_LEN; i++)
    mpz_urandomm (z[i], gmp_randstate, M);

  mpzspm = mpzspm_init (MPZSPV_MAX_LEN, M);

  spm = spm_init (mpzspm->spm[0].sp);
  
  spv = (spv_t) malloc (MPZSPV_MAX_LEN * sizeof (sp_t));
  spv_random (spv, MPZSPV_MAX_LEN, spm->sp);
  mpzspv = mpzspv_init (MPZSPV_MAX_LEN, mpzspm);
  mpzspv_random (mpzspv, 0, MPZSPV_MAX_LEN, mpzspm);
#endif
  
  MPZMOD_THRESHOLD = crossover2 (tune_mpres_mul_modmuln, tune_mpres_mul_mpz,
      1, 512, 10);
  
  printf ("#define MPZMOD_THRESHOLD %u\n", MPZMOD_THRESHOLD);
  
  REDC_THRESHOLD = crossover2 (tune_mpres_mul_mpz, tune_mpres_mul_redc,
      MPZMOD_THRESHOLD, 512, 10);
  
  printf ("#define REDC_THRESHOLD %u\n", REDC_THRESHOLD);

#ifdef HAVE_NTT
  SPV_NTT_GFP_DIF_RECURSIVE_THRESHOLD = 1 << crossover
    (tune_spv_ntt_gfp_dif_unrolled, tune_spv_ntt_gfp_dif_recursive, 1,
     MPZSPV_MAX_LOG2_LEN);

  printf ("#define SPV_NTT_GFP_DIF_RECURSIVE_THRESHOLD %u\n",
      SPV_NTT_GFP_DIF_RECURSIVE_THRESHOLD);
  
  SPV_NTT_GFP_DIT_RECURSIVE_THRESHOLD = 1 << crossover
    (tune_spv_ntt_gfp_dit_unrolled, tune_spv_ntt_gfp_dit_recursive, 1,
     MPZSPV_MAX_LOG2_LEN);

  printf ("#define SPV_NTT_GFP_DIT_RECURSIVE_THRESHOLD %u\n",
      SPV_NTT_GFP_DIT_RECURSIVE_THRESHOLD);
  
  MUL_NTT_THRESHOLD = 1 << crossover2 (tune_list_mul, tune_ntt_mul, 1,
      MPZSPV_MAX_LOG2_LEN - 1, 1);

  printf ("#define MUL_NTT_THRESHOLD %u\n", MUL_NTT_THRESHOLD);

  PREREVERTDIVISION_NTT_THRESHOLD = 1 << crossover2 (tune_PrerevertDivision,
      tune_ntt_PrerevertDivision, 1, MPZSPV_MAX_LOG2_LEN - 1, 1);

  printf ("#define PREREVERTDIVISION_NTT_THRESHOLD %u\n",
      PREREVERTDIVISION_NTT_THRESHOLD);

  POLYINVERT_NTT_THRESHOLD = 1 << crossover (tune_PolyInvert,
      tune_ntt_PolyInvert, 1, MPZSPV_MAX_LOG2_LEN);

  printf ("#define POLYINVERT_NTT_THRESHOLD %u\n", POLYINVERT_NTT_THRESHOLD);
  
  POLYEVALT_NTT_THRESHOLD = 1 << crossover (tune_polyevalT,
      tune_ntt_polyevalT, 1, MPZSPV_MAX_LOG2_LEN / 2);

  printf ("#define POLYEVALT_NTT_THRESHOLD %u\n", POLYEVALT_NTT_THRESHOLD);
  
  MPZSPV_NORMALISE_STRIDE = 1 << maximise (tune_mpzspv_normalise,
      7, MIN (MPZSPV_MAX_LOG2_LEN, 10));
  
  printf ("#define MPZSPV_NORMALISE_STRIDE %u\n", MPZSPV_NORMALISE_STRIDE);

  mpzspv_clear (mpzspv, mpzspm);
  free (spv);
  spm_clear (spm);
  mpzspm_clear (mpzspm);
  
  clear_list (x, MPZSPV_MAX_LEN);
  clear_list (y, MPZSPV_MAX_LEN);
  clear_list (z, MPZSPV_MAX_LEN);
  clear_list (t, list_mul_mem (MPZSPV_MAX_LEN / 2) + 3 * MPZSPV_MAX_LEN / 2);

  mpz_clear (M);
  gmp_randclear (gmp_randstate);
#endif
  
  return 0;
}

  
