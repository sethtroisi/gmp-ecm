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

#define NUM_CODELETS (sizeof(ntt_config) / sizeof(ntt_config[0]))

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
            sp_t primroot, sp_t order)
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

      c->ntt_pfa_run(xfixed, 1, p, vfixed);

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
          c->ntt_pfa_run(x, 1, p, vfixed);
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
  uint32_t i, j, k;
  uint32_t num_const;
  spv_t curr_const;
  nttdata_t *d;

  d = (nttdata_t *)calloc(1, sizeof(nttdata_t));
  if (d == NULL)
    return d;

  for (i = num_const = 0; i < NUM_CODELETS; i++)
    {
      const nttconfig_t *c = ntt_config[i];

      if (size % c->size != 0)
	continue;

      /* allocate room for each constant and its reciprocal */
      num_const += 2 * c->num_ntt_const;
      d->num_codelets++;
    }

  if (num_const == 0)
    goto error_free;

  d->codelets = (codelet_data_t *)malloc(d->num_codelets * 
		  			sizeof(codelet_data_t));
  if (d->codelets == NULL)
    goto error_free;

  d->codelet_const = curr_const = (spv_t)malloc(num_const * sizeof(sp_t));
  if (d->codelet_const == NULL)
    goto error_free;

  for (i = j = 0; i < NUM_CODELETS; i++)
    {
      const nttconfig_t *c = ntt_config[i];

      if (size % c->size != 0)
	continue;

      d->codelets[j].config = c;
      d->codelets[j].ntt_const = curr_const;
      j++;

      /* compute the constants, then append their reciprocals */

      num_const = c->num_ntt_const;
      c->nttdata_init(curr_const, p, recip, primroot, size);

      for (k = 0; k < num_const; k++)
	{
	  curr_const[num_const + k] = sp_ntt_reciprocal(curr_const[k], p);
	}
      curr_const += 2 * num_const;
    }

  return d;

error_free:
  ntt_free(d);
  return NULL;
}

/*-------------------------------------------------------------------------*/
void ntt_free(void *data)
{
  nttdata_t *d = (nttdata_t *)data;

  if (d == NULL)
    return;

  free(d->codelets);
  free(d->codelet_const);
  free(d);
}

/*-------------------------------------------------------------------------*/
uint32_t ntt_build_passes(nttpass_t *passes,
    		codelet_data_t *codelets, uint32_t num_codelets,
    		nttplan_t *plans, uint32_t num_plans,
		sp_t size, sp_t p, sp_t primroot, sp_t d)
{
  uint32_t i, j, k;

  for (i = j = 0; i < num_plans; i++)
    {
      nttplan_t * p = plans + i;
      codelet_data_t * c = NULL;

      for (k = 0; k < num_codelets; k++)
	{
	  c = codelets + k;
	  if (c->config->size == p->codelet_size)
	    break;
	}

      if (k == num_codelets)
	return 0;

      passes[j].pass_type = p->pass_type;

      switch (p->pass_type)
	{
	  case PASS_TYPE_DIRECT:
	    passes[j++].d.direct.codelet = c;
	    break;

	  case PASS_TYPE_PFA:
	    {
	      nttpass_t * pass = passes + j;
	      pass->d.pfa.codelets[pass->d.pfa.num_codelets++] = c;

	      if (i == num_plans - 1) /* PFA transform configured */
		j++;
	      break;
	    }

	  case PASS_TYPE_TWIDDLE_PRE:
	  case PASS_TYPE_TWIDDLE:
	    return 0;
	}
    }

  return j;
}

