#include <stdio.h>
#include <limits.h>
#include "gmp.h"
#include "ecm-impl.h"

#define N 32

INLINE void
mpn_mul_lo_basecase (mp_ptr rp, mp_srcptr np, mp_srcptr mp, mp_size_t n)
{
  mpn_mul_1 (rp, np, n, mp[0]);
  for (; --n;)
    mpn_addmul_1 (++rp, np, n, (++mp)[0]);
}

mp_size_t threshold[N];

static void
mpn_mul_lo_n_tune (mp_ptr rp, mp_srcptr np, mp_srcptr mp, mp_size_t n)
{
  mp_size_t k = threshold[n];

  switch (k)
    {
    case 0:
      return mpn_mul_n (rp, np, mp, n);
    case 1:
      return mpn_mul_lo_basecase (rp, np, mp, n);
    default:
      mpn_mul_n (rp, np, mp, k);
      rp += k;
      n -= k;
      mpn_mul_lo_n_tune (rp + n, np + k, mp, n);
      mpn_add_n (rp, rp, rp + n, n);
      mpn_mul_lo_n_tune (rp + n, np, mp + k, n);
      mpn_add_n (rp, rp, rp + n, n);
    }
}

static int
test (mp_ptr cp, mp_ptr ap, mp_ptr bp, mp_size_t n, int k)
{
  int st = cputime ();

  while (k--)
    mpn_mul_lo_n_tune (cp, ap, bp, n);

  return cputime () - st;
}

#define MINTIME 1000

int
main ()
{
  mp_size_t n, t, topt;
  mp_limb_t ap[N], bp[N], cp[2*N];
  int st[3], s, k;

  printf ("n\tmul_n\tlow_bc\tlow_n\n");

  threshold[0] = 0;
  threshold[1] = 0;

  for (n = 2; n < N; n++)
    {
      printf ("%u\t", n);
      mpn_random (ap, n);
      mpn_random (bp, n);

      /* calibrate */
      threshold[n] = 0;
      for (k = 1; (st[0] = test (cp, ap, bp, n, k)) < MINTIME; k *= 2);
      k = (int) (((double) k * (double) MINTIME) / (double) st[0]);

      printf ("%u\t", st[0] = test (cp, ap, bp, n, k)); /* mpn_mul_n */

      threshold[n] = 1;
      printf ("%u\t", st[1] = test (cp, ap, bp, n, k)); /* mpn_mul_lo_basecase */
      /* find optimal threshold for mpn_mul_lo_n */
      st[2] = INT_MAX;
      for (t = (n + 1) / 2; t < n; t++)
        {
          threshold[n] = t;
          s = test (cp, ap, bp, n, k);
          if (s < st[2])
            {
              st[2] = s;
              topt = t;
            }
        }

      printf ("%u(%u)\t", st[2], topt);
      threshold[n] = topt;

      if (st[1] <= st[2]) /* mpn_mul_lo_n slower than mpn_mul_lo_basecase */
	{
	  if (st[0] <= st[1])
	    {
	      threshold[n] = 0;
	      printf ("mul_n");
	    }
	  else
	    {
	      threshold[n] = 1;
	      printf ("low_bc");
	    }
	}
      else /* mpn_mul_lo_n faster than mpn_mul_lo_basecase */
	{
	  if (st[0] <= st[2])
	    {
	      threshold[n] = 0;
	      printf ("mul_n");
	    }
	  else
	    {
	      threshold[n] = topt;
	      printf ("low_n(%u)", topt);
	    }
	}

      printf ("\n");
    }

  printf ("#define MUL_LOW_THRESHOLD_TABLE {");
  for (n = 0; n < N; n++)
    {
      if (n)
	printf (",");
      printf ("%u", threshold[n]);
    }
  printf ("}\n");

  return 0;
}
