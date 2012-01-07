/* ECM stage 1 in batch mode, for initial point (x:z) with small coordinates,
   such that x and z fit into an unsigned long.
   For example we can start with (x=2:y=1) with the curve by^2 = x^3 + ax^2 + x
   with a = 4d-2 and b=16d+2, then we have to multiply by d=(a+2)/4 in the
   duplicates.
   With the change of variable x=b*X, y=b*Y, this curve becomes:
   Y^2 = X^3 + a/b*X^2 + 1/b^2*X.
*/

#include "ecm-impl.h"

/* (x1:z1) <- 2(x1:z1)
   (x2:z2) <- (x1:z1) + (x2:z2)
   assume (x2:z2) - (x1:z1) = (2:1)
   Uses 4 modular multiplies and 4 modular squarings.
*/
static void
dup_add (mpres_t x1, mpres_t z1, mpres_t x2, mpres_t z2,
         mpres_t q, mpres_t t, mpres_t u, mpres_t v, mpres_t w,
         unsigned long d, mpmod_t n)
{
  mpres_add (w, x1, z1, n); /* w = x1+z1 */
  mpres_sub (u, x1, z1, n); /* u = x1-z1 */
  mpres_add (t, x2, z2, n); /* t = x2+z2 */
  mpres_sub (v, x2, z2, n); /* v = x2-z2 */

  mpres_mul (t, t, u, n); /* t = (x1-z1)(x2+z2) */
  mpres_mul (v, v, w, n); /* v = (x2-z2)(x1+z1) */
  mpres_sqr (w, w, n);    /* w = (x1+z1)^2 */
  mpres_sqr (u, u, n);    /* u = (x1-z1)^2 */

  mpres_mul (x1, u, w, n); /* xdup = (x1+z1)^2 * (x1-z1)^2 */

  mpres_sub (w, w, u, n);   /* w = (x1+z1)^2 - (x1-z1)^2 */

  mpres_mul_ui (q, w, d, n); /* q = d * ((x1+z1)^2 - (x1-z1)^2) */

  mpres_add (u, u, q, n);  /* u = (x1-z1)^2 - d* ((x1+z1)^2 - (x1-z1)^2) */
  mpres_mul (z1, w, u, n); /* zdup = w * [(x1-z1)^2 - d* ((x1+z1)^2 - (x1-z1)^2)] */

  mpres_add (w, v, t, n);
  mpres_sub (v, v, t, n);

  mpres_sqr (v, v, n);
  mpres_sqr (x2, w, n);
  mpres_add (z2, v, v, n);
}


#define MAX_HEIGHT 32

#if ULONG_MAX == 4294967295
#define MAX_B1_BATCH 2977044736UL
#else
/* nth_prime(2^(MAX_HEIGHT-1)) */
#define MAX_B1_BATCH 50685770167UL
#endif

void
compute_s (mpz_t s, unsigned long B1)
{
  mpz_t acc[MAX_HEIGHT]; /* To accumulate products of prime powers */
  unsigned int i, j;
  unsigned long pi = 2, pp, maxpp;

  ASSERT_ALWAYS (B1 < MAX_B1_BATCH);

  for (j = 0; j < MAX_HEIGHT; j++)
    mpz_init (acc[j]); /* sets acc[j] to 0 */

  i = 0;
  while (pi <= B1)
    {
      pp = pi;
      maxpp = B1 / pi;
      while (pp <= maxpp)
          pp *= pi;

      if ((i & 1) == 0)
          mpz_set_ui (acc[0], pp);
      else
          mpz_mul_ui (acc[0], acc[0], pp);
			
      j = 0;
      /* We have accumulated i+1 products so far. If bits 0..j of i are all
         set, then i+1 is a multiple of 2^(j+1). */
      while ((i & (1 << j)) != 0)
        {
          /* we use acc[MAX_HEIGHT-1] as 0-sentinel below, thus we need
             j+1 < MAX_HEIGHT-1 */
          ASSERT (j + 1 < MAX_HEIGHT - 1);
          if ((i & (1 << (j + 1))) == 0) /* i+1 is not multiple of 2^(j+2),
                                            thus add[j+1] is "empty" */
            mpz_swap (acc[j+1], acc[j]); /* avoid a copy with mpz_set */
          else
            mpz_mul (acc[j+1], acc[j+1], acc[j]); /* accumulate in acc[j+1] */
          mpz_set_ui (acc[j], 1);
          j++;
        }

      i++;
      pi = getprime (pi);
    }

  for (mpz_set (s, acc[0]), j = 1; mpz_cmp_ui (acc[j], 0) != 0; j++)
    mpz_mul (s, s, acc[j]);
	
  getprime_clear (); /* free the prime tables, and reinitialize */
  
  for (i = 0; i < MAX_HEIGHT; i++)
      mpz_clear (acc[i]);
}


/* Input: x is initial point
          A is curve parameter in Montgomery's form:
          g*y^2*z = x^3 + a*x^2*z + x*z^2
          n is the number to factor
	  B1 is the stage 1 bound
   Output: If a factor is found, it is returned in x.
           Otherwise, x contains the x-coordinate of the point computed
           in stage 1 (with z coordinate normalized to 1).
	   B1done is set to B1 if stage 1 completed normally,
	   or to the largest prime processed if interrupted, but never
	   to a smaller value than B1done was upon function entry.
   Return value: ECM_FACTOR_FOUND_STEP1 if a factor, otherwise 
           ECM_NO_FACTOR_FOUND
*/
/*
For now we don't take into account go stop_asap and chkfilename
*/
int
ecm_stage1_batch (mpz_t f, mpres_t x, mpres_t A, mpmod_t n, double B1,
                  double *B1done, mpz_t s)
{
  unsigned long d;
  mpz_t x1, z1, x2, z2;
  unsigned long i;
  mpz_t q, t, u, v, w;
  int ret = ECM_NO_FACTOR_FOUND;

  MEMORY_TAG;
  mpz_init (x1);
  MEMORY_TAG;
  mpz_init (z1);
  MEMORY_TAG;
  mpz_init (x2);
  MEMORY_TAG;
  mpz_init (z2);
  MEMORY_TAG;
  mpres_init (q, n);
  MEMORY_TAG;
  mpres_init (t, n);
  MEMORY_TAG;
  mpres_init (u, n);
  MEMORY_TAG;
  mpres_init (v, n);
  MEMORY_TAG;
  mpres_init (w, n);
  MEMORY_UNTAG;

  /* initialize P */
  mpres_set (x1, x, n);
  mpres_set_ui (z1, 1, n); /* P1 <- 1P */

  /* Compute d=(A+2)/4 from A */
  mpz_add_ui (u, A, 2);
  if (mpz_fdiv_ui (u, 4) != 0)
    {
      outputf (OUTPUT_ERROR, "Error, A+2 should be divisible by 4\n");
      return ECM_ERROR;
    }
  mpz_div_2exp (u, u, 2);
  mpres_set_z_for_gcd (u, u, n);
  if (mpz_fits_ulong_p (u) == 0)
    {
      outputf (OUTPUT_ERROR,
               "Error, d=(A+2)/4 should fit in an ulong, A=%Zd\n", A);
      return ECM_ERROR;
    }
  d = mpz_get_ui (u);

  /* Compute 2P : no need to duplicate P, the coordinates are simple. */
  mpres_set_ui (x2, 9, n);
  mpres_set_ui (z2, d, n);
  mpres_mul_2exp (z2, z2, 6, n);
  mpres_add_ui (z2, z2, 8, n); /* P2 <- 2P = (9 : : 64d+8) */

  /* invariant: if j represents the upper bits of s,
     then P1 = j*P and P2=(j+1)*P */

  /* now perform the double-and-add ladder */
  for (i = mpz_sizeinbase (s, 2) - 1; i-- > 0;)
    {
      if (mpz_tstbit (s, i) == 0) /* (j,j+1) -> (2j,2j+1) */
        /* P2 <- P1+P2    P1 <- 2*P1 */
        dup_add (x1, z1, x2, z2, q, t, u, v , w, d, n);
      else /* (j,j+1) -> (2j+1,2j+2) */
          /* P1 <- P1+P2     P2 <- 2*P2 */
        dup_add (x2, z2, x1, z1, q, t, u, v, w, d, n);
    }

  *B1done=B1;

  if (!mpres_invert (u, z1, n)) /* Factor found? */
    {
      mpres_gcd (f, z1, n);
      ret = ECM_FACTOR_FOUND_STEP1;
    }
  mpres_mul (x, x1, u, n);

  mpz_clear (x1);
  mpz_clear (z1);
  mpz_clear (x2);
  mpz_clear (z2);
  mpz_clear (q);
  mpz_clear (t);
  mpz_clear (u);
  mpz_clear (v);
  mpz_clear (w);

  return ret;
}
