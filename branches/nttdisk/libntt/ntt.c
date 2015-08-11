#include <stdlib.h>
#include "ntt-impl.h"

static const nttconfig_t * ntt_config[] = 
{
  &ntt3_config,
  &ntt4_config,
  &ntt5_config,
  &ntt7_config,
  &ntt8_config,
  &ntt9_config,
  &ntt15_config,
  &ntt16_config,
  &ntt35_config,
  &ntt40_config,
};

/*-------------------------------------------------------------------------*/
const nttconfig_t ** 
ntt_master_list(void)
{
  return ntt_config;
}

uint32_t ntt_master_list_size(void)
{
  return sizeof(ntt_config) / sizeof(ntt_config[0]);
}

/*-------------------------------------------------------------------------*/
sp_t sp_ntt_reciprocal (sp_t w, sp_t p)
{
    /* compute (w << SP_TYPE_BITS) / p */

#if SP_TYPE_BITS == GMP_LIMB_BITS  /* use GMP functions */
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

	for (i = 0; i < n; i++) {

		spv_t pivot_row;

		for (j = i; j < m; j++) {
		    pivot = matrix[permute[j] * max_cols + i];
			if (pivot > 0)
                break;
		}
		SWAP(spv_size_t, permute[i], permute[j]);
		pivot_row = matrix + permute[i] * max_cols;
		pivot = sp_inv(pivot, p, d);

		for (j = i + 1; j < m; j++) {
			spv_t curr_row = matrix + permute[j] * max_cols;
			sp_t mult = sp_mul(curr_row[i], pivot, p, d);

			for (k = i; k < n; k++)
				curr_row[k] = sp_sub(curr_row[k], 
                            sp_mul(mult, pivot_row[k], p, d), p);

			b[permute[j]] = sp_sub(b[permute[j]], 
                            sp_mul(mult, b[permute[i]], p, d), p);
		}
	}

#if 0
    if (m==81)
    {
        for (j = 0; j < m; j++) {
			spv_t curr_row = matrix + permute[j] * max_cols;

            printf("%03lu: ", j);
            for (k = 0; k < n; k++)
                printf("%016llx ", curr_row[k]);
                //printf("%c ", curr_row[k] ? '*' ' ');

            printf("\t\t%016llx\n", b[permute[j]]);
        }
    }
#endif

	for (i = n - 1; (int32_t)i >= 0; i--) {

		spv_t curr_row = matrix + permute[i] * max_cols;
		sp_t sum = b[permute[i]];

		for (j = i + 1; j < n; j++) {
			sum = sp_sub(sum, 
                sp_mul(x[j], curr_row[j], p, d), p);
		}

		x[i] = sp_mul(sum, 
                sp_inv(curr_row[i], p, d), p, d);
	}
}

/*-------------------------------------------------------------------------*/
void
nttdata_init_generic(const nttconfig_t *c,
            spv_t out, sp_t p, sp_t d, 
            sp_t primroot, sp_t order, sp_t perm)
{
  /* compute the NTT constants; works for any squarefree NTT size */

  spv_size_t n = c->size;
  spv_size_t m = c->num_ntt_const;
  sp_t root = sp_pow(primroot, order / n, p, d);

  const uint8_t *must_be_unity = c->get_fixed_ntt_const();

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
      vfixed[i] = 1;

  for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++)
        xfixed[j] = 0;
      xfixed[i] = 1;

      c->ntt_run(xfixed, 1, 0, xfixed, 1, 0, 1, p, vfixed);

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
          c->ntt_run(x, 1, 0, x, 1, 0, 1, p, vfixed);
          vfixed[j] = 0;

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
void * ntt_init(sp_t size, sp_t primroot, sp_t p, sp_t recip)
{
  return (nttdata_t *)calloc(1, sizeof(nttdata_t));
}

/*-------------------------------------------------------------------------*/
void ntt_reset(void *data)
{
  spv_size_t i;
  nttdata_t *d = (nttdata_t *)data;

  if (d == NULL)
    return;

  for (i = 0; i < d->num_passes; i++)
    {
      free(d->passes[i].codelet_const);

      if (d->passes[i].pass_type == PASS_TYPE_TWIDDLE)
	sp_aligned_free(d->passes[i].d.twiddle.w);
    }

  d->num_passes = 0;
}

/*-------------------------------------------------------------------------*/
void ntt_free(void *data)
{
  spv_size_t i;
  nttdata_t *d = (nttdata_t *)data;

  if (d == NULL)
    return;

  ntt_reset(data);
  free(d);
}

/*-------------------------------------------------------------------------*/
static spv_t
make_twiddle(sp_t primroot, sp_t order, sp_t p, sp_t d,
    		spv_size_t rows, spv_size_t cols)
{
  spv_size_t size = rows * cols;
  spv_t res = (spv_t)sp_aligned_malloc(2 * (rows - 1) * cols * sizeof(sp_t));

  sp_t w = sp_pow(primroot, order / size, p, d);
  sp_t w_inc = 1;
  spv_size_t i = 0;
  spv_size_t j, k;

#ifdef HAVE_SIMD
  spv_size_t num_simd = SP_SIMD_VSIZE * (cols / SP_SIMD_VSIZE);

  for (; i < num_simd; i += SP_SIMD_VSIZE)
    {
      for (j = 0; j < SP_SIMD_VSIZE; j++)
	{
	  sp_t w0 = w_inc;

	  for (k = 0; k < rows - 1; k++)
	    {
	      res[2*i*(rows-1) + (2*k)*SP_SIMD_VSIZE + j] = w0;
	      res[2*i*(rows-1) + (2*k+1)*SP_SIMD_VSIZE + j] = sp_ntt_reciprocal(w0, p);
	      w0 = sp_mul(w0, w_inc, p, d);
	    }

	  w_inc = sp_mul(w_inc, w, p, d);
	}
    }
#endif

  for (; i < cols; i++)
    {
      sp_t w0 = w_inc;

      for (j = 0; j < rows - 1; j++)
	{
	  res[2*i*(rows-1) + 2*j] = w0;
	  res[2*i*(rows-1) + 2*j+1] = sp_ntt_reciprocal(w0, p);
	  w0 = sp_mul(w0, w_inc, p, d);
	}

      w_inc = sp_mul(w_inc, w, p, d);
    }

  return res;
}

/*-------------------------------------------------------------------------*/
uint32_t ntt_build_passes(nttdata_t *data,
    		nttplan_t *plans, uint32_t num_plans,
		sp_t size, sp_t p, sp_t primroot, sp_t order, sp_t d)
{
  uint32_t i, j;
  nttpass_t *passes = data->passes;
  uint32_t num_codelets = ntt_master_list_size();
  const nttconfig_t **codelets = ntt_master_list();

  for (i = 0; i < num_plans; i++)
    {
      nttplan_t * plan = plans + i;
      const nttconfig_t * c = NULL;

      for (j = 0; j < num_codelets; j++)
	{
	  c = codelets[j];
	  if (c->size == plan->codelet_size)
	    break;
	}

      if (j == num_codelets)
	return 0;

      passes[i].pass_type = plan->pass_type;
      passes[i].codelet = c;

      switch (plan->pass_type)
	{
	  case PASS_TYPE_DIRECT:
	    passes[i].d.direct.num_transforms = size / plan->codelet_size;
	    break;

	  case PASS_TYPE_PFA:
	    passes[i].d.pfa.cofactor = size / plan->codelet_size;
	    break;

	  case PASS_TYPE_TWIDDLE:
	    {
	      spv_size_t cols = size / plan->codelet_size;
	      spv_size_t rows = plan->codelet_size;

	      passes[i].d.twiddle.num_transforms = cols;
	      passes[i].d.twiddle.stride = cols;
	      passes[i].d.twiddle.w = make_twiddle(primroot, order, p, d, rows, cols);
	      if (plans[i+1].pass_type != PASS_TYPE_DIRECT)
		size = cols;
	      break;
	    }
	}
    }

  for (i = 0; i < num_plans; i++)
    {
      const nttconfig_t * c = passes[i].codelet;
      uint32_t num_const = c->num_ntt_const;

      passes[i].codelet_const = (sp_t *)malloc(2 * num_const * sizeof(sp_t));

      c->nttdata_init(passes[i].codelet_const, p, d, primroot, order, 
			passes[i].pass_type == PASS_TYPE_PFA ?
				passes[i].d.pfa.cofactor % 
					passes[i].codelet->size :
				1);

      for (j = 0; j < num_const; j++)
	passes[i].codelet_const[num_const + j] = sp_ntt_reciprocal(
	    				passes[i].codelet_const[j], p);
    }

  data->num_passes = i;
  return i;
}

/*-------------------------------------------------------------------------*/
static void 
ntt_run_recurse(spv_t x, sp_t p, nttdata_t *d, spv_size_t pass)
{
  spv_size_t i;
  nttpass_t *curr_pass = d->passes + pass;

  switch (curr_pass->pass_type)
    {
      case PASS_TYPE_DIRECT:
	curr_pass->codelet->ntt_run(
	    		x, 1, curr_pass->codelet->size,
	    		x, 1, curr_pass->codelet->size,
			curr_pass->d.direct.num_transforms,
			p, curr_pass->codelet_const);
	return;

      case PASS_TYPE_PFA:
	for (i = pass; i < d->num_passes; i++)
	  curr_pass[i].codelet->ntt_pfa_run(
	    		x, curr_pass[i].d.pfa.cofactor,
			p, curr_pass[i].codelet_const);
	return;
    }

  curr_pass->codelet->ntt_twiddle_run(
		x, curr_pass->d.twiddle.stride, 1,
		x, curr_pass->d.twiddle.stride, 1,
		curr_pass->d.twiddle.w,
		curr_pass->d.twiddle.num_transforms,
		p, curr_pass->codelet_const);

  /* if next pass is direct type, handle all the twiddle
     rows at once, otherwise recurse */

  if (curr_pass[1].pass_type == PASS_TYPE_DIRECT)
    {
      curr_pass[1].codelet->ntt_run(
	    		x, 1, curr_pass[1].codelet->size,
	    		x, 1, curr_pass[1].codelet->size,
			curr_pass[1].d.direct.num_transforms,
			p, curr_pass[1].codelet_const);
    }
  else
    {
      for (i = 0; i < curr_pass->codelet->size; i++)
	ntt_run_recurse(x + i * curr_pass->d.twiddle.stride,
	    		p, d, pass + 1);
    }
}

/*-------------------------------------------------------------------------*/
void 
ntt_run(spv_t x, sp_t p, void *data)
{
  nttdata_t *d = (nttdata_t *)data;

  ntt_run_recurse(x, p, d, 0);
}

