#include <stdio.h>
#include <limits.h>
#include "gmp.h"
#include "ecm.h"

#define N 100

INLINE void
mpn_mul_lo_basecase (mp_ptr rp, mp_srcptr np, mp_srcptr mp, mp_size_t n)
{
  mpn_mul_1 (rp, np, n, mp[0]);
  for (; --n;)
    mpn_addmul_1 (++rp, np, n, (++mp)[0]);
}

mp_size_t MPN_MUL_LO_THRESHOLD = 2;

mp_size_t threshold[N];

static void
mpn_mul_lo_n_tune (mp_ptr rp, mp_srcptr np, mp_srcptr mp, mp_size_t n)
{
  if (n < MPN_MUL_LO_THRESHOLD)
    mpn_mul_lo_basecase (rp, np, mp, n);
  else
    {
      mp_size_t k = threshold[n];

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
test1 (mp_ptr cp, mp_ptr ap, mp_ptr bp, mp_size_t n, int k)
{
  int st = cputime ();
  
  while (k--)
    mpn_mul_n (cp, ap, bp, n);
  
  return cputime () - st;
}

static int
test2 (mp_ptr cp, mp_ptr ap, mp_ptr bp, mp_size_t n, int k)
{
  int st = cputime ();
  
  while (k--)
    mpn_mul_lo_basecase (cp, ap, bp, n);
  
  return cputime () - st;
}

static int
test3 (mp_ptr cp, mp_ptr ap, mp_ptr bp, mp_size_t n, int k)
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

  for (n = 2; n < N; n++)
    {
      printf ("%u\t", n);
      mpn_random (ap, n);
      mpn_random (bp, n);

      /* calibrate */
      for (k = 1; (st[0] = test1 (cp, ap, bp, n, k)) < MINTIME; k *= 2);
      k = (int) (((double) k * (double) MINTIME) / (double) st[0]);

      printf ("%u\t", st[0] = test1 (cp, ap, bp, n, k)); /* mpn_mul_n */
      printf ("%u\t", st[1] = test2 (cp, ap, bp, n, k)); /* mpn_mul_lo_basecase */
      /* find optimal threshold for mpn_mul_lo_n */
      st[2] = INT_MAX;
      for (t = (n + 1) / 2; t < n; t++)
        {
          threshold[n] = t;
          s = test3 (cp, ap, bp, n, k);
          if (s < st[2])
            {
              st[2] = s;
              topt = t;
            }
        }

      printf ("%u(%u)\n", st[2], topt);
      threshold[n] = topt;

      if (st[2] >= st[1]) /* mpn_mul_lo_n slower than mpn_mul_lo_basecase */
        MPN_MUL_LO_THRESHOLD = n + 1;
    }

  return 0;
}
