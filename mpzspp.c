/* mpzspp.c - "mpz small prime polynomial" functions for arithmetic on polys
   reduced modulo a mpzspm

  Copyright 2005 Dave Newman.

  The SP Library is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or (at your
  option) any later version.

  The SP Library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
  License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with the SP Library; see the file COPYING.LIB.  If not, write to
  the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
  MA 02111-1307, USA.
*/

#include <malloc.h>
#include "sp.h"
#include <time.h> /* FIXME */
#include <stdio.h> /* FIXME */

void
mpzspp_init (mpzspp_t mpzspp, mpzspm_t mpzspm)
{
  unsigned int i;
  mpzspp->alloc_len = 0;
  mpzspp->len = 0;
  mpzspp->spv = (spv_t *) malloc (mpzspm->sp_num * sizeof (spv_t *));
  mpzspp->mpzspm = mpzspm;
  
  /* we could have used a calloc but let's allow for NULL != 0 */
  for (i = 0; i < mpzspm->sp_num; i++)
    mpzspp->spv[i] = NULL;
}

void
mpzspp_clear (mpzspp_t mpzspp)
{
  unsigned int i;
  for (i = 0; i < mpzspp->mpzspm->sp_num; i++)
    free (mpzspp->spv[i]);
  
  free (mpzspp->spv);
}

void 
mpzspp_realloc (mpzspp_t mpzspp, spv_size_t len)
{
  unsigned int i;
  for (i = 0; i < mpzspp->mpzspm->sp_num; i++)
    mpzspp->spv[i] = (spv_t) realloc (mpzspp->spv[i], len * sizeof (sp_t));
  
  mpzspp->alloc_len = len;
  
  if (mpzspp->len > len)
    mpzspp->len = len;
}

void
mpzspp_set_mpzp (mpzspp_t mpzspp, mpzp_t mpzp, spv_size_t len,
    spv_size_t offset)
{
  unsigned int i, sp_num;
  spv_size_t j;
  
  sp_num = mpzspp->mpzspm->sp_num;
  MPZSPP_REALLOC (mpzspp, len + offset);
  
  /* GMP's comments on mpn_preinv_mod_1:
   *
   * "This function used to be documented, but is now considered obsolete.  It
   * continues to exist for binary compatibility, even when not required
   * internally."
   *
   * It doesn't accept 0 as the dividend so we have to treat this case
   * separately */
  
  for (i = 0; i < len; i++)
    {
      if (SIZ(mpzp[i]) == 0)
	{
	  for (j = 0; j < sp_num; j++)
	    mpzspp->spv[j][i + offset] = 0;
	}
      else
        {
	  for (j = 0; j < sp_num; j++)
            mpzspp->spv[j][i + offset] =
              mpn_preinv_mod_1 (PTR(mpzp[i]), SIZ(mpzp[i]),
                mpzspp->mpzspm->spm[j].sp, mpzspp->mpzspm->spm[j].mul_c);
	}
    }

  mpzspp->len = MAX(len + offset, mpzspp->len);
}

void
mpzspp_set_zero (mpzspp_t mpzspp, spv_size_t len, spv_size_t offset)
{
  spv_size_t i;
}

/* B&S: ecrt mod m.
   A is clobbered and the resulting vector is stored in R. */
void
mpzspp_get_mpzp (mpzspp_t mpzspp, mpzp_t mpzp, spv_size_t len,
    spv_size_t offset)
{
  unsigned int i;
  spv_size_t j;
  unsigned int sp_num = mpzspp->mpzspm->sp_num;
  float *f = (float *) malloc (len * sizeof (float));
  float prime_recip;
  
  for (j = 0; j < len; j++)
    {
      f[j] = 0.5;
      mpz_set_ui (mpzp[j], 0);
    }
  
  for (i = 0; i < sp_num; i++)
    {
      prime_recip = 1.0 / (float) mpzspp->mpzspm->spm[i].sp;
      
      for (j = 0; j < len; j++)
        {
	  mpzspp->spv[i][j + offset] = 
	    sp_mul (mpzspp->spv[i][j + offset], mpzspp->mpzspm->crt3[i],
		mpzspp->mpzspm->spm[i].sp, mpzspp->mpzspm->spm[i].mul_c);
          
	  mpz_addmul_ui (mpzp[j], mpzspp->mpzspm->crt1[i],
	      mpzspp->spv[i][j + offset]);

	  f[j] += (float) mpzspp->spv[i][j + offset] * prime_recip;
        }
    }

  for (j = 0; j < len; j++)
      mpz_add (mpzp[j], mpzp[j], mpzspp->mpzspm->crt2[(unsigned int) f[j]]);
  
  free (f);
}  

void
mpzspp_sqr (mpzspp_t r, mpzspp_t x, int monic)
{
  unsigned int i;

  ASSERT (monic == 0 || monic == 1);

  r->mpzspm = x->mpzspm;

  if (x->len == 0)
    {
      r->len = 0;
      return;
    }

  MPZSPP_REALLOC (r, 1 << ceil_log_2 (2 * x->len - 1 + monic));
  
  for (i = 0; i < x->mpzspm->sp_num; i++)
    spv_mul (r->spv[i], x->spv[i], x->len, x->spv[i], x->len,
	monic, 0, 0, x->mpzspm->spm + i);

  r->len = 2 * x->len - 1 + monic;
}
  

void
mpzspp_mul (mpzspp_t r, mpzspp_t x, mpzspp_t y, int monic)
{
  unsigned int i;
  
  ASSERT (x->mpzspm == y->mpzspm);
  ASSERT (monic == 0 || monic == 1);

  r->mpzspm = x->mpzspm;

  if (x->len == 0 || y->len == 0)
    {
      r->len = 0;
      return;
    }
  
  MPZSPP_REALLOC (r, 1 << ceil_log_2 (x->len + y->len - 1 + monic));
  
  for (i = 0; i < x->mpzspm->sp_num; i++)
    spv_mul (r->spv[i], x->spv[i], x->len, y->spv[i], y->len,
	0, 0, monic, x->mpzspm->spm + i);
  
  r->len = x->len + y->len - 1 + monic;
}

void
mpzspp_mul_partial (mpzspp_t r, mpzspp_t x, mpzspp_t y, 
    spv_size_t k, spv_size_t l, int monic)
{
  unsigned int i;
  
  ASSERT (x->mpzspm == y->mpzspm);
  ASSERT (monic == 0 || monic == 1);

  r->mpzspm = x->mpzspm;

  if (x->len == 0 || y->len == 0)
    {
      r->len = 0;
      return;
    }
  
  MPZSPP_REALLOC (r, 1 << ceil_log_2 (x->len + y->len - 1 + monic));
  
  for (i = 0; i < x->mpzspm->sp_num; i++)
    spv_mul (r->spv[i], x->spv[i], x->len, y->spv[i], y->len, k, l, 
	monic, x->mpzspm->spm + i);
  
  r->len = l;
}

void
mpzspp_pwmul (mpzspp_t r, mpzspp_t x, mpzspp_t y)
{
  unsigned int i;

  for (i = 0; i < x->mpzspm->sp_num; i++)
    spv_pwmul (r->spv[i], x->spv[i], y->spv[i], x->len, x->mpzspm->spm[i].sp,
	x->mpzspm->spm[i].mul_c);

  r->len = x->len;
}

#define STRIDE MIN(256,len)

/* B&S: ecrt mod m mod p_j. */
void
mpzspp_normalize (mpzspp_t x, spv_size_t len, spv_size_t offset)
{
#if 0
  spv_size_t i;
  mpzp_t r = (mpzp_t) malloc (x->len * sizeof (mpz_t));

  for (i = 0; i < x->len; i++)
    mpz_init (r[i]);

  mpzspp_get_mpzp (x, r, x->len);
  mpzspp_set_mpzp (x, r, x->len, 0);

  for (i = 0; i < x->len; i++)
    mpz_clear (r[i]);
  
#else
    
  unsigned long i, j, k, l;
  unsigned int sp_num = x->mpzspm->sp_num;
  sp_t v;
  spv_t w;
  spm_t spm = x->mpzspm->spm;
  int st = clock();
  
  float *f = (float *) malloc (len * sizeof (float));

  spv_t *t = (spv_t *) malloc (sp_num * sizeof (spv_t));
  float prime_recip;
  spv_t s1 = (spv_t) calloc (3 * len, sizeof (sp_t));
  spv_t d1 = (spv_t) calloc (3 * len, sizeof (sp_t));

  for (i = 0; i < sp_num; i++)
    t[i] = (spv_t) malloc (STRIDE * sizeof (sp_t));
  
  for (i = 0; i < len; i++)
    f[i] = 0.5;
  
  /* FIXME: use B&S Theorem 2.2 */
  for (i = 0; i < sp_num; i++)
    {
      prime_recip = 1.0 / (float) spm[i].sp;
      
      for (j = 0; j < len; j++)
	{
	  x->spv[i][j + offset] = sp_mul (x->spv[i][j + offset],
	      x->mpzspm->crt3[i], spm[i].sp, spm[i].mul_c);
	  f[j] += (float) x->spv[i][j + offset] * prime_recip;
	}
    }
  
  for (l = 0; l < len; l += STRIDE)
  {
    for (i = 0; i < sp_num; i++)
      {
        for (j = 0; j < STRIDE; j++)
	  {
	    umul_ppmm (d1[3 * j + 1], d1[3 * j], x->mpzspm->crt5[i], (unsigned int) f[j + l]);
            d1[3 * j + 2] = 0;
	  }
	
        for (j = 0; j < sp_num; j++)
          {
	    w = x->spv[j];
	    v = x->mpzspm->crt4[i][j];
	    
	    for (k = 0; k < STRIDE; k++)
#if 1
	      umul_ppmm (s1[3 * k + 1], s1[3 * k], w[k + l], v);
#else
	      {
		sp_t t0, t1;
		umul_ppmm (t1, t0, w[k + l], v);
		add_sssaaaaaa (d1[3 * k + 2], d1[3 * k + 1], d1[3 * k],
                    d1[3 * k + 2], d1[3 * k + 1], d1[3 * k],
		    0, t1, t0);
	      }
#endif	      
	    /* this mpn_add_n accounts for about a third of the function's runtime */
	    mpn_add_n (d1, d1, s1, 3 * STRIDE);
          }
      
          /* FIXME: do we need to account for dividend == 0? */
          for (j = 0; j < STRIDE; j++)
	    t[i][j] = mpn_preinv_mod_1 (d1 + 3 * j, 3, spm[i].sp, spm[i].mul_c);
      }	  
  
    for (j = 0; j < sp_num; j++)
      for (k = 0; k < STRIDE; k++)
        x->spv[j][k + l + offset] = t[j][k];
  }
  
  for (i = 0; i < sp_num; i++)
    free (t[i]);
  
  free (t);
  free (s1);
  free (d1);
  free (f);

#endif
}  

/* Convert the mpzp to a mpzspp and do a fixed-length ntt on
 * the coefficients */
void
mpzspp_set_mpzp_ntt (mpzspp_t mpzspp, mpzp_t mpzp, spv_size_t mpzp_len,
    int monic, spv_size_t ntt_len)
{
  MPZSPP_REALLOC (mpzspp, MAX (mpzp_len, ntt_len));
  mpzspp_set_mpzp (mpzspp, mpzp, mpzp_len, 0);

  unsigned int i;
  spv_size_t j;
  sp_t root;
  spm_t spm;
  
  if (ntt_len > mpzp_len)
    {
      for (i = 0; i < mpzspp->mpzspm->sp_num; i++)
        {
	  spv_set_zero (mpzspp->spv[i] + mpzp_len, ntt_len - mpzp_len);
	  if (monic)
	    mpzspp->spv[i][mpzp_len] = 1;

	  spm = mpzspp->mpzspm->spm + i;
	  root = sp_pow (spm->prim_root, (spm->sp - 1) / ntt_len,
	      spm->sp, spm->mul_c);
	  spv_ntt_gfp_dif (mpzspp->spv[i], ntt_len, spm->sp, spm->mul_c, root);
	}
    }
  else
    {
      for (i = 0; i < mpzspp->mpzspm->sp_num; i++)
        {
	  for (j = ntt_len; j < mpzp_len; j += ntt_len)
	    spv_add (mpzspp->spv[i], mpzspp->spv[i], mpzspp->spv[i] + j,
	        ntt_len, mpzspp->mpzspm->spm[i].sp);
	  if (monic)
	    mpzspp->spv[i][0] = 1;

	  spm = mpzspp->mpzspm->spm + i;
	  root = sp_pow (spm->prim_root, (spm->sp - 1) / ntt_len,
	      spm->sp, spm->mul_c);
	  spv_ntt_gfp_dif (mpzspp->spv[i], ntt_len, spm->sp, spm->mul_c, root);
	}
    }
  mpzspp->len = ntt_len;
}

void mpzspp_get_mpzp_ntt (mpzspp_t mpzspp, mpzp_t mpzp, spv_size_t len)
{
  unsigned int i;
  spv_size_t j;
  sp_t root;
  spm_t spm;
  
  for (i = 0; i < mpzspp->mpzspm->sp_num; i++)
    {
      spm = mpzspp->mpzspm->spm + i;
      root = sp_pow (spm->prim_root, (spm->sp - 1) - (spm->sp - 1) / len,
	  spm->sp, spm->mul_c);
      spv_ntt_gfp_dit (mpzspp->spv[i], mpzspp->len, spm->sp, spm->mul_c, root);
      root = sp_inv (mpzspp->len, spm->sp, spm->mul_c);
      spv_mul_sp (mpzspp->spv[i], mpzspp->spv[i], root, len,
	  spm->sp, spm->mul_c);
    }
  
  mpzspp_get_mpzp_ntt (mpzspp, mpzp, len);
}
