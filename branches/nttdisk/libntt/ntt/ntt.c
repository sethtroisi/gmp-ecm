#include "ntt-impl.h"

const nttgroup_t * 
X(ntt_master_group_list)[] = 
{
#ifdef HAVE_SSE42
  & MANGLE_SSE42(X(ntt_group_simd)),
#endif
#ifdef HAVE_SSE2
  & MANGLE_SSE2(X(ntt_group_simd)),
#endif
#ifdef HAVE_AVX
  & MANGLE_AVX(X(ntt_group_simd)),
#endif
#ifdef HAVE_AVX2
  & MANGLE_AVX2(X(ntt_group_simd)),
#endif
#ifdef HAVE_FMA
  & MANGLE_FMA(X(ntt_group_simd)),
#endif
  & X(ntt_group),
};

uint32_t
X(get_master_group_list_size)()
{
  return sizeof(X(ntt_master_group_list)) /
		sizeof(X(ntt_master_group_list)[0]);
}

/*-------------------------------------------------------------------------*/
sp_t X(sp_ntt_reciprocal)(sp_t w, sp_t p)
{
    /* compute (w << SP_TYPE_BITS) / p */

#if SP_NUMB_BITS == 50  /* don't need a reciprocal */

    sp_t hi, lo;
    sp_split(hi, lo, w);
    return hi;

#elif SP_TYPE_BITS == GMP_LIMB_BITS  /* use GMP functions */
    mp_limb_t recip, dummy;

    udiv_qrnnd (recip, dummy, w, 0, p);
    return recip;

#elif SP_TYPE_BITS < GMP_LIMB_BITS  /* ordinary division */

    return ((mp_limb_t)w << SP_TYPE_BITS) / p;

#else  /* worst case: bit-at-a-time */

    sp_t r = w;
    sp_t q = 0;
    mp_limb_t i;

    for (i = 0; i < SP_TYPE_BITS; i++)
      {
	q += q;
	r += r;
	if (r >= p)
	{
	  r -= p;
	  q |= 1;
	}
      }
    return q;

#endif
}

/*-------------------------------------------------------------------------*/
#define SWAP(type, a, b) {type tmp = (a); (a) = (b); (b) = tmp; }

static void 
solve_ffmatrix(spv_t matrix, spv_size_t max_cols,
		spv_t x, spv_t b,
		sp_t p, sp_t d, spv_size_t m, spv_size_t n) 
{
  spv_size_t i, j, k;
  spv_size_t *permute = (spv_size_t *)alloca(m * sizeof(spv_size_t));
  sp_t pivot;

  for (i = 0; i < m; i++)
    permute[i] = i;

  for (i = 0; i < n; i++) 
    {

      spv_t pivot_row;

      for (j = i; j < m; j++) 
	{
	  pivot = matrix[permute[j] * max_cols + i];
	  if (pivot > 0)
	    break;
	}

      SWAP(spv_size_t, permute[i], permute[j]);
      pivot_row = matrix + permute[i] * max_cols;
      pivot = sp_inv(pivot, p, d);

      for (j = i + 1; j < m; j++) 
	{
	  spv_t curr_row = matrix + permute[j] * max_cols;
	  sp_t mult = sp_mul(curr_row[i], pivot, p, d);

	  for (k = i; k < n; k++)
	    curr_row[k] = sp_sub(curr_row[k], 
			    	sp_mul(mult, pivot_row[k], p, d), p);

	  b[permute[j]] = sp_sub(b[permute[j]], 
				sp_mul(mult, b[permute[i]], p, d), p);
	}
    }

  for (i = n - 1; (int32_t)i >= 0; i--) 
    {
      spv_t curr_row = matrix + permute[i] * max_cols;
      sp_t sum = b[permute[i]];

      for (j = i + 1; j < n; j++)
	sum = sp_sub(sum, sp_mul(x[j], curr_row[j], p, d), p);

      x[i] = sp_mul(sum, sp_inv(curr_row[i], p, d), p, d);
    }
}

/*-------------------------------------------------------------------------*/
void
X(nttdata_init_generic)(const nttconfig_t *c,
            spv_t out, sp_t p, sp_t d, 
            sp_t primroot, sp_t order, sp_t perm)
{
  /* compute the NTT constants; works for any NTT size */

  spv_size_t n = c->size;
  spv_size_t m = c->num_ntt_const;
  sp_t root = sp_pow(primroot, order / n, p, d);

  const uint8_t *must_be_unity = c->fixed_ntt_const;

  spv_t mat = (spv_t)alloca(n * n * m * sizeof(sp_t));
  spv_t rhs = (spv_t)alloca(n * n * sizeof(sp_t));
  spv_t soln = (spv_t)alloca(n * n * sizeof(sp_t));
  spv_t x = (spv_t)alloca(n * sizeof(sp_t));
  spv_t xfixed = (spv_t)alloca(n * sizeof(sp_t));
  spv_t vfixed = (spv_t)alloca(2 * m * sizeof(sp_t));
  spv_t wmult = (spv_t)alloca(n * sizeof(sp_t));

  spv_size_t i, j, k;
  spv_size_t col;

  wmult[0] = 1;
  rhs[0] = 1;
  for (i = 1; i < n; i++)
    {
      wmult[i] = sp_mul(wmult[i-1], root, p, d);
      rhs[i] = 1;
    }

  for (i = 1; i < n; i++)
    {
      for (j = 0; j < n; j++)
        rhs[n * i + j] = sp_mul(rhs[n * (i - 1) + j], wmult[j], p, d);
    }

  if (perm != 1)
    {
      for (j = n; j < n * n; j++)
        rhs[j] = sp_pow(rhs[j], perm, p, d);
    }

  for (i = 0; i < 2 * m; i++)
    vfixed[i] = 0;

  for (i = 0; i < m; i++)
    if (must_be_unity[i])
      {
	vfixed[i] = 1;
	vfixed[i+m] = X(sp_ntt_reciprocal)(1, p);
      }

  for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++)
        xfixed[j] = 0;
      xfixed[i] = 1;

      c->ntt_run(xfixed, 1, 0, 1, p, vfixed);

      #ifdef HAVE_PARTIAL_MOD
      for (j = 0; j < n; j++)
        while (xfixed[j] >= p)
          xfixed[j] -= p;
      #endif

      for (j = 0; j < n; j++)
        rhs[n * i + j] = sp_sub(rhs[n * i + j], xfixed[j], p);

      for (j = col = 0; j < m; j++)
        {
          if (must_be_unity[j])
            continue;

          for (k = 0; k < n; k++)
            x[k] = 0;
          x[i] = 1;

          vfixed[j] = 1;
          vfixed[j+m] = X(sp_ntt_reciprocal)(1, p);
          c->ntt_run(x, 1, 0, 1, p, vfixed);
          vfixed[j] = 0;
          vfixed[j+m] = 0;

          #ifdef HAVE_PARTIAL_MOD
          for (k = 0; k < n; k++)
              while(x[k] >= p)
                x[k] -= p;
          #endif

          for (k = 0; k < n; k++)
            mat[(n * i + k) * m + col] = sp_sub(x[k], xfixed[k], p);
          col++;
        }
    }

  solve_ffmatrix(mat, m, soln, rhs, p, d, n * n, col);

  for (i = j = 0; i < m; i++)
    {
      if (must_be_unity[i])
        out[i] = 1;
      else
        out[i] = soln[j++];
    }
}

/*-------------------------------------------------------------------------*/
void X(ntt_reset)(nttdata_t *d)
{
  spv_size_t i;

  for (i = 0; i < d->num_passes; i++)
    {
      sp_aligned_free(d->passes[i].codelet_const);

      if (d->passes[i].pass_type == PASS_TYPE_TWIDDLE)
	sp_aligned_free(d->passes[i].d.twiddle.w);
    }

  d->num_passes = 0;
}

/*-------------------------------------------------------------------------*/
static void
alloc_twiddle_interleaved(mpzspm_t mpzspm, nttpass_t *pass,
    		spv_size_t rows, spv_size_t cols,
		spv_size_t vsize, uint32_t order)
{
  spv_size_t i, j, k, m;
  spv_size_t num_batches = ((mpzspm->sp_num + vsize - 1) / vsize);
  spv_size_t alloc = 2 * cols * (rows - 1) * vsize;
  spv_t res = (spv_t)sp_aligned_malloc(alloc * num_batches * sizeof(sp_t));

  if (vsize * num_batches > mpzspm->sp_num)
    memset(res, 0, alloc * num_batches * sizeof(sp_t));

  for (i = 0; i < num_batches; i++)
    {
      spv_t res0 = res + alloc * i;

      for (j = 0; j < vsize; j++)
	{
	  if (i * vsize + j < mpzspm->sp_num)
	    {
	      spm_t spm = mpzspm->spm[i * vsize + j];
	      sp_t p = spm->sp;
	      sp_t d = spm->mul_c;
	      sp_t primroot = spm->primroot;
	      sp_t w = sp_pow(primroot, order / (rows * cols), p, d);
	      sp_t w_inc = 1;

	      for (k = 0; k < cols; k++)
		{
		  sp_t w0 = w_inc;

		  for (m = 0; m < rows - 1; m++)
		    {
		      res0[(2*k*(rows-1) + 2*m)*vsize + j] = w0;
		      res0[(2*k*(rows-1) + 2*m+1)*vsize + j] = X(sp_ntt_reciprocal)(w0, p);
		      w0 = sp_mul(w0, w_inc, p, d);
		    }

    		  w_inc = sp_mul(w_inc, w, p, d);
    		}
	    }
	}
    }

  pass->d.twiddle.w = res;
  pass->d.twiddle.twiddle_size = 2 * cols * (rows - 1);
}

/*-------------------------------------------------------------------------*/
static void
alloc_twiddle_packed(mpzspm_t mpzspm, nttpass_t *pass,
    		spv_size_t rows, spv_size_t cols,
		spv_size_t vsize, uint32_t order)
{
  spv_size_t i, j, k, m;
  spv_size_t num_simd = vsize * ((cols + vsize - 1) / vsize);
  spv_size_t alloc = 2 * (rows - 1) * num_simd;
  spv_t res = (spv_t)sp_aligned_malloc(alloc * 
      			mpzspm->sp_num * sizeof(sp_t));

  for (m = 0; m < mpzspm->sp_num; m++)
    {
      spm_t spm = mpzspm->spm[m];
      sp_t p = spm->sp;
      sp_t d = spm->mul_c;
      sp_t primroot = spm->primroot;
      sp_t w = sp_pow(primroot, order / (rows * cols), p, d);
      sp_t w_inc = 1;
      spv_t res0 = res + alloc * m;

      for (i = 0; i < num_simd; i += vsize)
	{
	  for (j = 0; j < vsize; j++)
	    {
	      sp_t w0 = w_inc;

	      if (i + j < cols)
		{
		  for (k = 0; k < rows - 1; k++)
		    {
		      res0[2*i*(rows-1) + (2*k)*vsize + j] = w0;
		      res0[2*i*(rows-1) + (2*k+1)*vsize + j] = X(sp_ntt_reciprocal)(w0, p);
		      w0 = sp_mul(w0, w_inc, p, d);
		    }

		  w_inc = sp_mul(w_inc, w, p, d);
		}
	      else
		{
		  for (k = 0; k < rows - 1; k++)
		    {
		      res0[2*i*(rows-1) + (2*k)*vsize + j] = 0;
		      res0[2*i*(rows-1) + (2*k+1)*vsize + j] = 0;
		    }
		}
	    }
	}
    }

  pass->d.twiddle.w = res;
  pass->d.twiddle.twiddle_size = alloc;
}

/*-------------------------------------------------------------------------*/
static void
alloc_const_interleaved(mpzspm_t mpzspm, nttpass_t * pass, 
    			uint32_t vsize, uint32_t order)
{
  spv_size_t i, j, k;
  const nttconfig_t * c = pass->codelet;
  uint32_t num_const = c->num_ntt_const;
  spv_size_t num_batches = (mpzspm->sp_num + vsize - 1) / vsize;
  spv_size_t alloc = 2 * num_const * vsize;
  spv_t tmp = (spv_t)alloca(num_const * sizeof(sp_t));

  pass->codelet_const = (spv_t)sp_aligned_malloc(
      				alloc * num_batches * sizeof(sp_t));

  if (vsize * num_batches > mpzspm->sp_num)
    memset(pass->codelet_const, 0, alloc * num_batches * sizeof(sp_t));

  for (i = 0; i < num_batches; i++)
    {
      spv_t res0 = pass->codelet_const + alloc * i;

      for (j = 0; j < vsize; j++)
	{
	  if (i * vsize + j < mpzspm->sp_num)
	    {
	      spm_t spm = mpzspm->spm[i * vsize + j];
	      sp_t p = spm->sp;
	      sp_t d = spm->mul_c;
	      sp_t primroot = spm->primroot;

	      c->nttdata_init(tmp, p, d, primroot, order, 
			    pass->pass_type == PASS_TYPE_PFA ?
      				pass->d.pfa.cofactor % pass->codelet->size :
				1);

    	      for (k = 0; k < num_const; k++)
		{
		  res0[(2*k)*vsize + j] = tmp[k];
		  res0[(2*k+1)*vsize + j] = X(sp_ntt_reciprocal)(tmp[k], p);
		}
	    }
	}
    }

  pass->const_size = 2 * num_const;
}

/*-------------------------------------------------------------------------*/
static void
alloc_const_packed(mpzspm_t mpzspm, nttpass_t * pass, uint32_t order)
{
  uint32_t i, j;
  const nttconfig_t * c = pass->codelet;
  uint32_t num_const = c->num_ntt_const;

  pass->codelet_const = (spv_t)sp_aligned_malloc(2 * mpzspm->sp_num * 
					num_const * sizeof(sp_t));

  for (i = 0; i < mpzspm->sp_num; i++)
    {
      spm_t spm = mpzspm->spm[i];
      sp_t p = spm->sp;
      sp_t d = spm->mul_c;
      sp_t primroot = spm->primroot;
      spv_t curr_const = pass->codelet_const + 2 * i * num_const;

      c->nttdata_init(curr_const, p, d, primroot, order, 
			pass->pass_type == PASS_TYPE_PFA ?
				pass->d.pfa.cofactor % pass->codelet->size :
				1);

      for (j = 0; j < num_const; j++)
  	  curr_const[num_const + j] = 
	    		X(sp_ntt_reciprocal)(curr_const[j], p);
    }

  pass->const_size = 2 * num_const;
}

/*-------------------------------------------------------------------------*/
uint32_t X(ntt_build_passes)(
    		nttplangroup_t *plans, mpzspm_t mpzspm, 
		uint32_t ntt_size, uint32_t max_ntt_size)
{
  uint32_t i, j;
  nttdata_t *data = &mpzspm->nttdata;
  nttpass_t *passes = data->passes;
  spv_size_t num_transforms = 1;

  for (i = 0; i < plans->num_plans; i++)
    {
      nttplan_t * plan = plans->plans + i;
      nttpass_t * pass = passes + i;
      const nttgroup_t * g = X(ntt_master_group_list)[plan->group_type];
      const nttconfig_t **codelets = g->get_transform_list();
      uint32_t num_codelets = g->num_transforms;
      const nttconfig_t * c = NULL;

      for (j = 0; j < num_codelets; j++)
	{
	  c = codelets[j];
	  if (c->size == plan->codelet_size)
	    break;
	}

      if (j == num_codelets)
	return 0;

      pass->pass_type = plan->pass_type;
      pass->codelet = c;
      pass->num_transforms = num_transforms;
      pass->vsize = g->vsize;

      switch (plan->pass_type)
	{
	  case PASS_TYPE_DIRECT:
	    break;

	  case PASS_TYPE_PFA:
	    pass->d.pfa.ntt_size = ntt_size;
	    pass->d.pfa.cofactor = ntt_size / plan->codelet_size;
	    break;

	  case PASS_TYPE_TWIDDLE:
	    {
	      spv_size_t cols = ntt_size / plan->codelet_size;
	      spv_size_t rows = plan->codelet_size;

	      pass->num_transforms = cols;
	      pass->d.twiddle.stride = cols;

	      if (mpzspm->interleaved)
		alloc_twiddle_interleaved(mpzspm, pass, rows, cols, 
		    			g->vsize, max_ntt_size);
	      else
		alloc_twiddle_packed(mpzspm, pass, rows, cols, 
		    			g->vsize, max_ntt_size);

	      num_transforms = rows;
	      ntt_size = cols;
	      break;
	    }
	}

      if (mpzspm->interleaved && g->vsize > 1)
	alloc_const_interleaved(mpzspm, pass, g->vsize, max_ntt_size);
      else
	alloc_const_packed(mpzspm, pass, max_ntt_size);
    }

  data->num_passes = i;
  return i;
}

/*-------------------------------------------------------------------------*/
static void 
ntt_run_interleaved(nttdata_t *d, spv_t x, 
    		spv_size_t vsize, spv_t p,
    		spv_size_t id, spv_size_t pass)
{
  spv_size_t i;
  nttpass_t *curr_pass = d->passes + pass;

  switch (curr_pass->pass_type)
    {
      case PASS_TYPE_DIRECT:
	curr_pass->codelet->ntt_run_interleaved(
	    		x, vsize, vsize * curr_pass->codelet->size,
			curr_pass->num_transforms, p,
			vsize,
			curr_pass->codelet_const + 
				id * curr_pass->const_size);
	return;

      case PASS_TYPE_PFA:
	for (i = pass; i < d->num_passes; i++)
	  {
	    curr_pass = d->passes + i;
	    curr_pass->codelet->ntt_pfa_run_interleaved(
	    		x, vsize, vsize * curr_pass->d.pfa.ntt_size,
			curr_pass->num_transforms,
			p, curr_pass->d.pfa.cofactor,
			vsize,
			curr_pass->codelet_const + 
				id * curr_pass->const_size);
	  }
	return;
    }

  curr_pass->codelet->ntt_twiddle_run_interleaved(
		x, vsize * curr_pass->d.twiddle.stride, vsize,
		curr_pass->num_transforms, p, 
		vsize,
		curr_pass->codelet_const +
			id * curr_pass->const_size,
		curr_pass->d.twiddle.w +
			id * curr_pass->d.twiddle.twiddle_size);

  /* recurse on the rows; if the next pass doesn't use twiddle
     factors, handle all the rows at once */

  if (curr_pass[1].pass_type != PASS_TYPE_TWIDDLE)
    {
      ntt_run_interleaved(d, x, vsize, p, id, pass + 1);
    }
  else
    {
      for (i = 0; i < curr_pass->codelet->size; i++)
	ntt_run_interleaved(d, x + i * vsize * curr_pass->d.twiddle.stride,
	    		vsize, p, id, pass + 1);
    }
}

/*-------------------------------------------------------------------------*/
static void 
ntt_run_packed(nttdata_t *d, spv_t x, sp_t p,
    		spv_size_t id, spv_size_t pass)
{
  spv_size_t i;
  nttpass_t *curr_pass = d->passes + pass;

  switch (curr_pass->pass_type)
    {
      case PASS_TYPE_DIRECT:
	curr_pass->codelet->ntt_run(
	    		x, 1, curr_pass->codelet->size,
			curr_pass->num_transforms,
			p, curr_pass->codelet_const + 
				id * curr_pass->const_size);
	return;

      case PASS_TYPE_PFA:
	for (i = pass; i < d->num_passes; i++)
	  {
	    curr_pass = d->passes + i;
	    curr_pass->codelet->ntt_pfa_run(
	    		x, 1, curr_pass->d.pfa.ntt_size,
			curr_pass->num_transforms,
			p, curr_pass->d.pfa.cofactor,
			curr_pass->codelet_const + 
				id * curr_pass->const_size);
	  }
	return;
    }

  curr_pass->codelet->ntt_twiddle_run(
		x, curr_pass->d.twiddle.stride, 1,
		curr_pass->num_transforms,
		p, curr_pass->codelet_const +
			id * curr_pass->const_size,
		curr_pass->d.twiddle.w +
			id * curr_pass->d.twiddle.twiddle_size);

  /* recurse on the rows; if the next pass doesn't use twiddle
     factors, handle all the rows at once */

  if (curr_pass[1].pass_type != PASS_TYPE_TWIDDLE)
    {
      ntt_run_packed(d, x, p, id, pass + 1);
    }
  else
    {
      for (i = 0; i < curr_pass->codelet->size; i++)
	ntt_run_packed(d, x + i * curr_pass->d.twiddle.stride,
	    		p, id, pass + 1);
    }
}

/*-------------------------------------------------------------------------*/
void 
X(ntt_run)(void * m, mpz_t * x, uint32_t ntt_size)
{
  uint32_t i;
  mpzspm_t mpzspm = (mpzspm_t)m;

#if SP_NUMB_BITS == 50
  fpu_precision_t old_prec = 0;
  if (!fpu_precision_is_ieee())
    old_prec = fpu_set_precision_ieee();
#endif

  if (mpzspm->interleaved)
    {
      for (i = 0; i < mpzspm->sp_num; i += mpzspm->max_vsize)
	ntt_run_interleaved(&mpzspm->nttdata, 
	    		mpzspm->work + i * ntt_size,
			mpzspm->max_vsize,
			mpzspm->sp + i, i, 0);
    }
  else
    {
      for (i = 0; i < mpzspm->sp_num; i++)
	ntt_run_packed(&mpzspm->nttdata, 
	    		mpzspm->work + i * ntt_size,
			mpzspm->spm[i]->sp, i, 0);
    }

#if SP_NUMB_BITS == 50
  if (old_prec != 0)
    fpu_clear_precision(old_prec);
#endif
}

