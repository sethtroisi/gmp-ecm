#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "ecm.h"
#ifdef TESTDRIVE
#include <asm/msr.h>
#include <string.h>
#endif

#define DEBUG 1
#define CHECKSUM 1

static mpz_t gt;
static int gt_inited = 0;

#define min(a,b) ((a)<(b)?(a):(b))

#define CACHESIZE 512U

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#define INLINE inline
#else
#define UNUSED
#define INLINE
#endif

/* a' <- a+b, b' <- a-b. */

#define ADDSUB_MOD(a, b) \
  mpz_sub (gt, a, b); \
  mpz_add (a, a, b);  \
  F_mod_gt (b, n);    \
  F_mod_1 (a, n);

mp_limb_t __gmpn_mod_34lsub1 (mp_limb_t *src, mp_size_t size);
void F_mod_1 (mpz_t, unsigned int);
void F_mod_gt (mpz_t, unsigned int);
void mpz_absadd_2exp (mpz_t, unsigned int);
void F_divby2 (mpz_t, mpz_t, unsigned int);
void F_divby3_1 (mpz_t, unsigned int);
void F_divby5_1 (mpz_t, unsigned int);
void F_mul_sqrt2exp (mpz_t, mpz_t, int, unsigned int);
void F_mul_sqrt2exp_2 (mpz_t, mpz_t, int, unsigned int);
void F_fft_dif (mpz_t *, int, int, int);
void F_fft_dit (mpz_t *, int, int, int);
unsigned int *make_scramble_list (int);
unsigned int *make_scramble_list_r4 (int);
void F_fft_mfa (mpz_t *, const unsigned int, const unsigned int, int, const unsigned int);
void F_fft_mfa_mul (mpz_t *, mpz_t *, const unsigned int, const unsigned int, const unsigned int);
unsigned int F_toomcook4 (mpz_t *, mpz_t *, mpz_t *, unsigned int, unsigned int, mpz_t *);
unsigned int F_karatsuba (mpz_t *, mpz_t *, mpz_t *, unsigned int, unsigned int, mpz_t *);
/*
unsigned int F_mul (mpz_t *, mpz_t *, mpz_t *, unsigned int, unsigned int, mpz_t *);
unsigned int F_mul_trans (mpz_t *, mpz_t *, mpz_t *, unsigned int, unsigned int, mpz_t *);
*/

static int radix2 = 0, do_mfa = 1;

#ifdef TESTDRIVE
int do_timing = 0, force_karatsuba = 0, verbose = 0, nofft = 0;
unsigned long long timer_start=0, timer_stop;

void * allocate_function (size_t ALLOC_SIZE)
{
  void *ptr;
  
  ptr = malloc (ALLOC_SIZE);
  if (ptr == NULL)
    {
      fprintf (stderr, "allocate_function: Failed to allocate %d bytes.\n", 
                       ALLOC_SIZE);
      exit (EXIT_FAILURE);
    }
  
  return ptr;
}

void * reallocate_function (void *PTR, UNUSED size_t OLD_SIZE, size_t NEW_SIZE)
{
  void *ptr;
  
#ifdef DEBUG_REALLOC
  printf ("reallocate_function: Reallocating %d bytes to %d bytes at %p\n",
          OLD_SIZE, NEW_SIZE, PTR);
#endif
  
  ptr = realloc (PTR, NEW_SIZE);

  if (ptr == NULL)
    {
      fprintf (stderr, "reallocate_function: Failed to reallocate %d bytes.\n", 
                       NEW_SIZE);
      exit (EXIT_FAILURE);
    }
  
  return ptr;
}

void deallocate_function (void *PTR, UNUSED size_t SIZE)
{
  free (PTR);
}

void printpoly(mpz_t *P, unsigned int l, char *t) 
{
  unsigned int i;
  
  if (verbose < 2)
    return;
  
  for (i=0; i<l; i++) {
    gmp_printf("%s[%d]: %Zd\n", t, i, P[i]);
  }
}

void print2polys(mpz_t *P, mpz_t *Q, unsigned int l, char *t, char *u) 
{
  unsigned int i;
  
  if (verbose < 2)
    return;
  
  for (i=0; i<l; i++) {
    gmp_printf("%s[%d]: %Zd   %s[%d]: %Zd\n", t, i, P[i], u, i, Q[i]);
  }
}


/* Allocates memory for a poly of degree len-1, and inits the mpz_t's.
   Each mpz_t get space for size bits allocated */
mpz_t *
poly_init (unsigned int len, unsigned int size)
{
  mpz_t *r;
  unsigned int i;
  
  r = (mpz_t *) malloc (len * sizeof (mpz_t));
  
  if (r)
    for (i = 0; i < len; i++)
      mpz_init2 (r[i], size);

  return r;
}


/* Add 2^n+1 to RS. Only used in main(), to make negative residues 
   positive for output. */

void
F_add_F (mpz_t RS, unsigned int n)
{
  mpz_set_ui (gt, 1);
  mpz_mul_2exp (gt, gt, n);
  mpz_add_ui (gt, gt, 1);
  mpz_add (gt, gt, RS); /* F_n has (n/GMP_NUMB_BITS)+1 limbs, dest gets */
  mpz_set (RS, gt);     /* allocated one extra during add, causing growth */
}

#endif /* TESTDRIVE */

/* RS -> RS (mod 2^n+1). If input |RS| < 2^(2*n), result |RS| < 2^(n+1) */

INLINE void 
F_mod_1 (mpz_t RS, unsigned int n)
{
  mp_size_t size;
  mp_limb_t v;
  
  size = mpz_size (RS);

#ifdef DEBUG_MOD
  printf ("F_mod_1: input log_2(RS) = %d, ", mpz_sizeinbase(RS, 2));

  if ((unsigned int) size <= n / GMP_NUMB_BITS)
    {
      printf ("done.\n");
      return;
    }
#endif

  if ((unsigned int) size == n / GMP_NUMB_BITS + 1)
    {
      int sgn;
      sgn = mpz_sgn (RS);          /* Remember original sign */
      v = mpz_getlimbn (RS, n / GMP_NUMB_BITS);
      mpz_tdiv_r_2exp (RS, RS, n); /* Just a truncate. RS < 2^n. Can make
                                      RS zero and so change sgn(RS)! */
      if (sgn == -1)
          mpz_add_ui (RS, RS, v);
      else
          mpz_sub_ui (RS, RS, v);
    }
  else if ((unsigned int) size > n / GMP_NUMB_BITS + 1)
    {                              /* Assuming |RS| < 2^(2*n) */
      mpz_tdiv_q_2exp (gt, RS, n); /* |gt| < 2^n */
      mpz_tdiv_r_2exp (RS, RS, n); /* |RS| < 2^n */
      mpz_sub (RS, RS, gt);        /* |RS| < 2^(n+1) */
    }
#ifdef DEBUG_MOD
  printf ("output log_2(RS) = %d\n", mpz_sizeinbase (RS, 2));
#endif
}


/* R = gt (mod 2^n+1) */

INLINE void 
F_mod_gt (mpz_t R, unsigned int n)
{
  mp_size_t size;
  mp_limb_t v;
  
  size = mpz_size (gt);

#ifdef DEBUG
  if (R == gt)
    {
      fprintf (stderr, "F_mod_gt: R == gt\n");
      exit (EXIT_FAILURE);
    }
#endif

#ifdef DEBUG_MOD
  printf ("F_mod_gt: input log_2(gt) = %d, ", mpz_sizeinbase(gt, 2));

  if ((unsigned int) size <= n / GMP_NUMB_BITS)
    {
      printf ("done.\n");
      mpz_set (R, gt);
      return;
    }
#endif

  if ((unsigned int) size == n / GMP_NUMB_BITS + 1)
    {
      int sgn;
      sgn = mpz_sgn (gt);
      v = mpz_getlimbn (gt, n / GMP_NUMB_BITS);
      mpz_tdiv_r_2exp (gt, gt, n); /* Just a truncate */
      if (sgn == -1)
          mpz_add_ui (R, gt, v);
      else
          mpz_sub_ui (R, gt, v);
    }
  else if ((unsigned int) size > n / GMP_NUMB_BITS + 1)
    {
      mpz_tdiv_q_2exp (R, gt, n);
      mpz_tdiv_r_2exp (gt, gt, n); /* Just a truncate */
      mpz_sub (R, gt, R);
    }
  else 
    mpz_set (R, gt);

#ifdef DEBUG_MOD
  printf ("output log_2(R) = %d\n", mpz_sizeinbase (R, 2));
#endif
}


/* R = S + sgn(S)*(2^e) */

void
mpz_absadd_2exp (mpz_t RS, unsigned int e)
{
  mp_size_t siz, limb_idx, bit_idx;
  mp_limb_t cy;
  int sgn;
  
  limb_idx = e / GMP_NUMB_BITS;
  bit_idx = e % GMP_NUMB_BITS;
  siz = mpz_size (RS);
  sgn = (mpz_sgn (RS) >= 0) ? 1 : -1;
  
  if (limb_idx >= RS->_mp_alloc)
    {
#ifdef DEBUG_REALLOC
      printf ("mpz_absadd_2exp: reallocating RS\n");
#endif
      mpz_realloc2 (RS, (limb_idx + 1) * GMP_NUMB_BITS);
    }
  
  /* Now RS->_mp_alloc > limb_idx) */
  
  while (siz <= limb_idx)
    {
      RS->_mp_d[siz++] = 0;
      RS->_mp_size += sgn;
    }
  
  /* Now RS->_mp_alloc >= siz > limb_idx */
  
  cy = mpn_add_1 (RS->_mp_d + limb_idx, RS->_mp_d + limb_idx,
                  siz - limb_idx, ((mp_limb_t)1) << bit_idx);
  if (cy)
    {
      if (RS->_mp_alloc <= siz)
        {
#ifdef DEBUG_REALLOC
          printf ("mpz_absadd_2exp: reallocating RS for carry\n");
#endif
          mpz_realloc2 (RS, (siz + 1) * GMP_NUMB_BITS);
        }

      RS->_mp_d[siz] = 1;
      RS->_mp_size += sgn;
    }
}

/* R = S / 2 (mod 2^n + 1). S == gt is ok */

void 
F_divby2 (mpz_t R, mpz_t S, unsigned int n)
{
  int odd, sgn;
  
  odd = mpz_odd_p (S);
  sgn = mpz_sgn (S);
  mpz_tdiv_q_2exp (R, S, 1);
  
  if (odd)
    {
      /* We shifted out a set bit at the bottom. With negative wrap-around,
         that becomes -2^(n-1), so we add -2^(n-1) + 2^n+1 = 2^(n-1)+1.
         If |S| < 2^(n+1), |R| < 2^n + 2^(n-1) + 1 < 2^(n+1) for n > 1. */
      
      mpz_absadd_2exp (R, n-1);
      if (sgn < 0)
        mpz_sub_ui (R, R, 1);
      else
        mpz_add_ui (R, R, 1);
    }
}


/* RS = RS / 3 (mod 2^n + 1). RS == gt is ok */

void 
F_divby3_1 (mpz_t RS, unsigned int n)
{
  /* 2^2^m == 1 (mod 3) for m>0, thus F_m == 2 (mod 3) */
  int mod, sgn;
  
  sgn = mpz_sgn (RS);
  mod = __gmpn_mod_34lsub1 (RS->_mp_d, mpz_size (RS)) % 3;

  if (mod == 1)
    {
      /* Add F_m. If |RS| < 2^(n+1), |RS|+F_m < 3*2^n+1 */
      mpz_absadd_2exp (RS, n);
      if (sgn >= 0)
        mpz_add_ui (RS, RS, 1);
      else
        mpz_sub_ui (RS, RS, 1);
    }
  else if (mod == 2)
    {
      /* Add 2 * F_m.  If |RS| < 2^(n+1), |RS|+2*F_m < 4*2^n+2 */
      mpz_absadd_2exp (RS, n + 1);
      if (sgn >= 0)
        mpz_add_ui (RS, RS, 2);
      else
        mpz_sub_ui (RS, RS, 2);
    }

  mpz_divby3_1op (RS); /* |RS| < (4*2^n+2)/3 < 2^(n+1) */
}

void 
F_divby5_1 (mpz_t RS, unsigned int n)
{
  /* 2^2^m == 1 (mod 5) for m>1, thus F_m == 2 (mod 5) */
  int mod, sgn;

  sgn = mpz_sgn (RS);
  mod = __gmpn_mod_34lsub1 (RS->_mp_d, mpz_size (RS)) % 5;

  if (mod == 1)
    {
      /* Add 2 * F_m == 4 (mod 5) */
      mpz_absadd_2exp (RS, n + 1);
      if (sgn == 1)
        mpz_add_ui (RS, RS, 2);
      else
        mpz_sub_ui (RS, RS, 2);
    }
  else if (mod == 2)
    {
      /* Add 4 * F_m == 3 (mod 5) */
      mpz_absadd_2exp (RS, n + 2);
      if (sgn == 1)
        mpz_add_ui (RS, RS, 4);
      else
        mpz_sub_ui (RS, RS, 4);
    }
  else if (mod == 3)
    {
      /* Add F_m == 3 (mod 5) */
      mpz_absadd_2exp (RS, n);
      if (sgn == 1)
        mpz_add_ui (RS, RS, 1);
      else
        mpz_sub_ui (RS, RS, 1);
    }
  else if (mod == 4)
    {
      /* Add 3 * F_m == 1 (mod 5) */
      mpz_absadd_2exp (RS, n);
      mpz_absadd_2exp (RS, n + 1);
      if (sgn == 1)
        mpz_add_ui (RS, RS, 3);
      else
        mpz_sub_ui (RS, RS, 3);
    }

#ifdef DEBUG
  if (!mpz_divisible_ui_p (RS, 5))
    {
      fprintf (stderr, "F_divby5_1: not divisible by 5\n");
      exit (EXIT_FAILURE);
    }
#endif

  mpz_divexact_ui (RS, RS, 5);
}


/* A 2^(m+2) length convolution is possible:
   (2^(3n/4) - 2^(n/4))^2 == 2 (mod 2^n+1) 
   so we have an element of order 2^(m+2) of simple enough form
   to use it as a root of unity the transform */

/* Multiply by sqrt(2)^e (mod F_m).  n = 2^m */
/* R = (S * sqrt(2)^e) % (2^n+1) */
/* R == S is ok, but neither must be == gt */
/* Assumes abs(e) < 4*n */

void 
F_mul_sqrt2exp (mpz_t R, mpz_t S, int e, unsigned int n) 
{
  int chgsgn = 0, odd;

#ifdef DEBUG
  if (S == gt || R == gt) 
    {
      printf ("F_mul_sqrt2exp: %c == gt\n", S==gt?'S':'R');
      exit (EXIT_FAILURE);
    }
  if (abs (e) >= 4 * n)
    {
      printf ("F_mul_sqrt2exp: e = %d, abs(e) >= 4*n\n", e);
      exit (EXIT_FAILURE);
    }
#endif

  if (e < 0)
    e += 4 * n;
  /* 0 <= e < 4*n */
  if ((unsigned) e >= 2 * n)    /* sqrt(2)^(2*n) == -1 (mod F_m), so */
    {
      e -= 2 * n;               /* sqrt(2)^e == -sqrt(2)^(e-2*n) (mod F_m) */
      chgsgn = 1;
    }				/* Now e < 2*n */

#ifdef DEBUG_PERF
  if (e == 0)
    printf ("F_mul_sqrt2exp: called for trivial case %s1\n", chgsgn ? "-" : "");
#endif

  odd = e & 1;
  e >>= 1;

  if (odd)
    {
      /* Multiply by sqrt(2) == 2^(3n/4) - 2^(n/4) */
      /* S * (2^(3n/4) - 2^(n/4)) == 2^(n/4) * (S*2^(n/2) - S) */
      mpz_mul_2exp (gt, S, n / 2);
      mpz_sub (gt, gt, S);
      mpz_tdiv_q_2exp (R, gt, n / 4 * 3);
      mpz_tdiv_r_2exp (gt, gt, n / 4 * 3);
      mpz_mul_2exp (gt, gt, n / 4);
      mpz_sub (R, gt, R);
      
      if (e != 0)
        {
          mpz_tdiv_q_2exp (gt, R, n-e);
          mpz_tdiv_r_2exp (R, R, n-e);
          mpz_mul_2exp (R, R, e);
          mpz_sub (R, R, gt);
        }
    }
  else if (e != 0) 
    {
      /*  S     = a*2^(n-e) + b,   b < 2^(n-e)  */
      /*  S*2^e = a*2^n + b*2^e = b*2^e - a */
      /*  b*2^e < 2^(n-e)*2^e = 2^n */
      mpz_tdiv_q_2exp (gt, S, n - e); /* upper e bits (=a) into gt */
      mpz_tdiv_r_2exp (R, S, n - e);  /* lower n-e bits (=b) into R */
                                      /* This is simply a truncate if S == R */
      mpz_mul_2exp (R, R, e);         /* R < 2^n */
      mpz_sub (R, R, gt);
    } else 
      mpz_set (R, S);

  if (chgsgn) 
    mpz_neg (R, R);

#ifdef DEBUG_MOD
  if (mpz_sizeinbase (R, 2) > n) 
    {
      printf ("F_mul_sqrt2exp: log_2(|R|) == %u > 2^n\n", mpz_sizeinbase(R, 2));
      F_mod_1 (R, n);
    }
#endif
}

/* Same, but input may be gt. Input and output must not be identical */
void 
F_mul_sqrt2exp_2 (mpz_t R, mpz_t S, int e, unsigned int n)
{
  int chgsgn = 0, odd;

#ifdef DEBUG
  if (S == R || R == gt) {
    printf("F_mul_sqrt2exp_2: R == %s\n", R==S ? "S" : "gt");
    exit(EXIT_FAILURE);
  }
#endif

#ifdef DEBUG_PERF
  if (abs (e) >= 4 * n)
    {
      printf ("F_mul_sqrt2exp_2: e = %d > 4*n = %d\n", e, 4*n);
    }
#endif

  if (e < 0)
    e += 4 * n;
  if ((unsigned) e >= 2 * n)    /* sqrt(2)^(2*n) == -1 (mod F_m), so */
    {
      e -= 2 * n;               /* sqrt(2)^e == -sqrt(2)^(e-2*n) (mod F_m) */
      chgsgn = 1;
    }				/* Now e < 2*n */

#ifdef DEBUG_PERF
  if (e == 0)
    printf ("F_mul_sqrt2exp_2: called for trivial case %s1\n", chgsgn ? "-" : "");
#endif

  odd = e & 1;
  e >>= 1;

  if (odd != 0)
    {
      mpz_set (R, S); /* Neccessary?  n/32 mov*/
      mpz_mul_2exp (gt, S, n / 2); /* May overwrite S  n/32 mov */
      mpz_sub (gt, gt, R); /* n/32 sub*/

      mpz_tdiv_q_2exp (R, gt, n / 4 * 3); /* 3*(n/32)/4 mov */
      mpz_tdiv_r_2exp (gt, gt, n / 4 * 3); /* Just a truncate */
      mpz_mul_2exp (gt, gt, n / 4); /* 3*(n/32)/4 mov */
      mpz_sub (R, gt, R); /* (n/32)/4 sub, 3*(n/32)/4 mov */
      
      if (e != 0)
        {
          mpz_tdiv_q_2exp (gt, R, n - e);
          mpz_tdiv_r_2exp (R, R, n - e);
          mpz_mul_2exp (R, R, e);
          mpz_sub (R, R, gt);
        }
    } 
  else if (e != 0) 
    {
      mpz_tdiv_q_2exp (R, S, n - e); /* upper e bits into R */
      mpz_tdiv_r_2exp (gt, S, n - e); /* lower n-e bits into gt */
      mpz_mul_2exp (gt, gt, e);
      mpz_sub (R, gt, R);
    } else 
      mpz_set (R, S);

  if (chgsgn == -1) 
    mpz_neg (R, R);

#ifdef DEBUG_MOD
  if (mpz_sizeinbase (R, 2) > n)
    {
      printf ("F_mul_sqrt2exp_2(%s%d): log_2(|R|) == %u > 2^n\n", 
              sgn==1?"":"-", e, mpz_sizeinbase(R, 2));
      F_mod_1 (R, n);
    }
#endif
}

#define A0s A[0]
#define A1s A[l << stride2]
#define A2s A[2 * l << stride2]
#define A3s A[3 * l << stride2]
#define A0is A[i << stride2]
#define A1is A[(i + l) << stride2]
#define A2is A[(i + 2 * l) << stride2]
#define A3is A[(i + 3 * l) << stride2]

/* Decimation-in-frequency FFT. Unscrambled input, scrambled output. */
/* Elements are (mod 2^n+1), l and n must be powers of 2, l must be <= 4*n. */
/* Performs forward transform */

void 
F_fft_dif (mpz_t *A, int l, int stride2, int n) 
{
  int i, omega = (4 * n) / l, iomega;

  if (l <= 1)
    return;
  
#ifdef DEBUG
  if ((4 * n) % l != 0)
    {
      fprintf (stderr, "F_fft_dif: l does not divide 4*n\n");
      exit (EXIT_FAILURE);
    }
#endif

  if (l == 2)
    {
      ADDSUB_MOD(A[0], A[1<<stride2]);
      return;
    }

  if (!radix2)
    {
#ifdef TESTDRIVE
      if (verbose >= 2)
        printf ("Radix-4 DIF, len=%d, omega=%d\n", l, omega);
#endif

      l /= 4;

      mpz_sub (gt, A1s, A3s);            /* gt = a1 - a3 */
      mpz_add (A1s, A1s, A3s);           /* A1 = a1 + a3 */
      F_mul_sqrt2exp_2 (A3s, gt, n, n);  /* A3 = i * (a1 - a3) */
      
      mpz_sub (gt, A[0], A2s);           /* gt = a0 - a2 */
      mpz_add (A[0], A[0], A2s);         /* A0 = a0 + a2 */

      mpz_sub (A2s, A[0], A1s);          /* A2 = a0 - a1 + a2 - a3 */
      mpz_add (A[0], A[0], A1s);         /* A0 = a0 + a1 + a2 + a3 */
      mpz_add (A1s, gt, A3s);            /* A1 = a0 - a2 + i * (a1 - a3) */
      mpz_sub (A3s, gt, A3s);            /* A3 = a0 - a2 - i * (a1 - a3) */

      for (i = 1, iomega = omega; i < l; i++, iomega += omega)
        {
          mpz_sub (gt, A1is, A3is);
          mpz_add (A1is, A1is, A3is);
          F_mul_sqrt2exp_2 (A3is, gt, n, n);
          
          mpz_sub (gt, A0is, A2is);
          mpz_add (A0is, A0is, A2is);
          
          mpz_sub (A2is, A0is, A1is);
          mpz_add (A0is, A0is, A1is);
          mpz_add (A1is, gt, A3is);
          mpz_sub (A3is, gt, A3is);
          F_mul_sqrt2exp (A1is, A1is, iomega, n);
          F_mul_sqrt2exp (A2is, A2is, 2 * iomega, n);
          F_mul_sqrt2exp (A3is, A3is, 3 * iomega, n);
        }

      if (l > 1)
        {
          F_fft_dif (A, l, stride2, n);
          F_fft_dif (A + (l << stride2), l, stride2, n);
          F_fft_dif (A + (2 * l << stride2), l, stride2, n);
          F_fft_dif (A + (3 * l << stride2), l, stride2, n);
        }
      return;
    }

  l /= 2;

  ADDSUB_MOD(A[0], A1s);

  for (i = 1, iomega = omega; i < l; i++, iomega += omega) 
    {
      mpz_sub (gt, A0is, A1is);
      mpz_add (A0is, A0is, A1is);
      F_mul_sqrt2exp_2 (A1is, gt, iomega, n);
      F_mod_1 (A0is, n);
    }
  
  F_fft_dif (A, l, stride2, n);
  F_fft_dif (A + (l << stride2), l, stride2, n);
}


/* Decimation-in-time inverse FFT. Scrambled input, unscrambled output */
/* Does not perform divide-by-length. l, and n as in F_fft_dif() */

void 
F_fft_dit (mpz_t *A, int l, int stride2, int n) 
{
  int i, omega = (4 * n) / l, iomega;
  
  if (l <= 1)
    return;

#ifdef DEBUG
  if ((4 * n) % l != 0)
    {
      fprintf (stderr, "F_fft_dit: l does not divide 4*n\n");
      exit (EXIT_FAILURE);
    }
#endif

  if (l == 2)
    {
      ADDSUB_MOD(A[0], A[1<<stride2]);
      return;
    }

  if (!radix2)
    {
#ifdef TESTDRIVE
      if (verbose >= 2)
        printf ("Radix-4 DIT, len=%d, omega=%d\n", l, omega); 
#endif

      l /= 4;
      
      if (l > 1)
        {
          F_fft_dit (A, l, stride2, n);
          F_fft_dit (A + (l << stride2), l, stride2, n);
          F_fft_dit (A + (2 * l << stride2), l, stride2, n);
          F_fft_dit (A + (3 * l << stride2), l, stride2, n);
        }

      mpz_sub (gt, A3s, A1s);            /* gt = -(a1 - a3) */
      mpz_add (A1s, A1s, A3s);           /* A1 = a1 + a3 */
      F_mul_sqrt2exp_2 (A3s, gt, n, n);  /* A3 = i * -(a1 - a3) */
      
      mpz_sub (gt, A[0], A2s);           /* gt = a0 - a2 */
      mpz_add (A[0], A[0], A2s);         /* A0 = a0 + a2 */
      
      mpz_sub (A2s, A[0], A1s);          /* A2 = a0 - a1 + a2 - a3 */
      mpz_add (A[0], A[0], A1s);         /* A0 = a0 + a1 + a2 + a3 */
      mpz_add (A1s, gt, A3s);            /* A1 = a0 - a2 + i * -(a1 - a3) */
      mpz_sub (A3s, gt, A3s);            /* A3 = a0 - a2 - i * -(a1 - a3) */

      for (i = 1, iomega = omega; i < l; i++, iomega += omega)
        {
          /* Divide by omega^i. Since sqrt(2)^(4*n) == 1 (mod 2^n+1), 
             this is like multiplying by omega^(4*n-i) */
          F_mul_sqrt2exp (A1is, A1is, 4 * n - iomega, n);
          F_mul_sqrt2exp (A2is, A2is, 4 * n - 2 * iomega, n);
          F_mul_sqrt2exp (A3is, A3is, 4 * n - 3 * iomega, n);

          mpz_sub (gt, A3is, A1is);
          mpz_add (A1is, A1is, A3is);
          F_mul_sqrt2exp_2 (A3is, gt, n, n);

          mpz_sub (gt, A0is, A2is);
          mpz_add (A0is, A0is, A2is);
          
          mpz_sub (A2is, A0is, A1is);
          mpz_add (A0is, A0is, A1is);
          mpz_add (A1is, gt, A3is);
          mpz_sub (A3is, gt, A3is);
        }
      return;
    }

  l /= 2;

  F_fft_dit (A, l, stride2, n);
  F_fft_dit (A + (l << stride2), l, stride2, n);
  
  ADDSUB_MOD(A[0], A1s);
  
  for (i = 1, iomega = 4*n - omega; i < l; i++, iomega -= omega) 
    {
      F_mul_sqrt2exp (A1is, A1is, iomega, n);
      mpz_sub (gt, A0is, A1is);
      mpz_add (A0is, A0is, A1is);
      F_mod_gt (A1is, n);
      F_mod_1 (A0is, n);
    }
  
}

unsigned int *
make_scramble_list (int n)
{
  int i, j, ldn, bit;
  unsigned int *list;
  
  if (n == 0)
    return NULL;
  
  list = (unsigned int *) malloc (n * sizeof (int));
  if (list == NULL)
    {
      printf ("make_scramble_list: could not allocate memory for list\n");
      exit (EXIT_FAILURE);
    }
    
  for (i = n, ldn = 0; (i & 1) == 0; i >>= 1, ldn++);
  
  if (i != 1)
    {
      printf ("make_scramble_list: n is not a power of 2\n");
      exit (EXIT_FAILURE);
    }
  
  j = 0;
  list[0] = 0;
  for (i = 1; i < n; i++)
    {
      bit = 1 << (ldn - 1); 
      while (j & bit)
        {
          j ^= bit;
          bit >>= 1;
        }
      j ^= bit;
      list[i] = j;
/*      printf("scramble(%d) = %d\n", i, j); */
    }

  return list;
}


unsigned int *
make_scramble_list_r4 (int n)
{
  int i, j, k, ldn, t;
  unsigned int *list;
  
  if (n == 0)
    return NULL;
  
  list = (unsigned int *) malloc (n * sizeof (int));
  if (list == NULL)
    {
      printf ("make_scramble_list_r4: could not allocate memory for list\n");
      exit (EXIT_FAILURE);
    }
    
  for (i = n, ldn = 0; (i & 1) == 0; i >>= 1, ldn++);
#ifdef DEBUG
  if (i != 1)
    {
      printf ("make_scramble_list_r4: n is not a power of 2\n");
      exit (EXIT_FAILURE);
    }
#endif
  
  list[0] = 0;
  if ((ldn & 1) == 0) /* Will we have a radix-2 step? */
    for (i = 1; i < n; i++)
      {
        j = i & 3;
        t = i >> 2;
        for (k = 2; k < ldn; k += 2)
          {
            j <<= 2;
            j |= t & 3;
            t >>= 2;
          }
        list[i] = j;
/*        printf("scramble(%d) = %d\n", i, j); */
      }
  else
    for (i = 1; i < n; i++)
      {
        j = i & 1;
        t = i >> 1;
        for (k = 1; k < ldn; k += 2)
          {
            j <<= 2;
            j |= t & 3;
            t >>= 2;
          }
        list[i] = j;
/*        printf("scramble(%d) = %d\n", i, j); */
      }

  return list;
}


/* Computes DFT (mod 2^n+1) using Matrix Fourier Algorithm.
   Matrix has r rows (each of length c) and c columns (each of length r).
   Result of forward transform is a transposed and scrambled (digit-reversal 
   of indices) matrix, inverse transform expects it like that.
   One row, that is c residues of n bits each, should fit into cache
*/

void
F_fft_mfa (mpz_t *A, const unsigned int r, const unsigned int c, int sign, const unsigned int n)
{
  unsigned int i, j, s2;
  unsigned int *scramble;

#ifdef TESTDRIVE
  if (verbose)
    printf ("MFA transform of %d rows and %d columns\n", r, c);
#endif

#ifdef DEBUG
  if ((4*n) % (c*r) != 0) 
    {
      fprintf(stderr, "F_fft_mfa: 4*n (%d) not a multiple of rows*columns (%d*%d)\n",
              4*n, r, c);
      exit(EXIT_FAILURE);
    }
#endif
  
  if (radix2)
    scramble = make_scramble_list (r);
  else
    scramble = make_scramble_list_r4 (r);
  
  /* Compute s2 = log_2(c) for use as stride */
  for (i = 1, s2 = 0; i < c; i <<= 1, s2++); 

  if (sign == 1)
    { /* Forward transform */

      /* DFT on columns */
      for (i = 0; i < c; i++)
        {
          F_fft_dif (A + i, r, s2, n);
          /* The rows are now index-scrambeled */
        }


      /* DFT on rows */

      /* Row with index 0, need no multipy by \omega */
      F_fft_dif (A, c, 0, n);

      /* The other rows */
      for (i = 1; i < r; i++)
        {
          unsigned int iomega = scramble[i] * ((4 * n) / (r * c)), ijomega;

          /* Multiply by \omega^{ij} */
          for (j = 1, ijomega = iomega; j < c; j++, ijomega += iomega)
            F_mul_sqrt2exp (A[i * c + j], A[i * c + j], ijomega, n);

          F_fft_dif (A + i * c, c, 0, n);
        }

    } else { /* Inverse transform */

      /* DFT on rows */
      F_fft_dit (A, c, 0, n);
      for (i = 1; i < r; i++)
        {
          unsigned int iomega = scramble[i] * ((4 * n) / (r * c)), ijomega;
          F_fft_dit (A + i * c, c, 0, n);

          /* Divide by \omega^{ij} */
          for (j = 1, ijomega = 4 * n - iomega; j < c; j++, ijomega -= iomega)
            F_mul_sqrt2exp (A[i * c + j], A[i * c + j], ijomega, n);
        }
      
      /* DFT on columns */
      for (i = 0; i < c; i++)
        F_fft_dit (A + i, r, s2, n);
    }

  free (scramble);
}


/* Computes convolution product A*B (mod 2^n+1) using Matrix Fourier Algorithm.
   B is assumed to be FFT'd my MFA already (rows and columns index scrambled)
*/

void
F_fft_mfa_mul (mpz_t *A, mpz_t *B, const unsigned int r, const unsigned int c, const unsigned int n)
{
  unsigned int i, j, s2, l2;
  unsigned int *scramble;

  if (radix2)
    scramble = make_scramble_list (r);
  else
    scramble = make_scramble_list_r4 (r);
  
  /* Compute s2 = log_2(c) for use as stride */
  for (i = 1, s2 = 0; i < c; i <<= 1, s2++); 

  /* Compute l2 = log_2(len) for divide-by-length */
  for (i = 1, l2 = 0; i < r; i <<= 1, l2++);
  l2 += s2;
  /* We later need to divide by (sqrt(2))^(2*l2), so we'll keep 2*l2 */
  l2 <<= 1;


  /* DFT on columns */
  for (i = 0; i < c; i++)
    {
      F_fft_dif (A + i, r, s2, n);
      /* The rows are now index-scrambeled */
    }


  /* DFT on rows, multiply, back transform of rows */

  /* Row with index 0, needs no multipy by \omega, just divide by length */
  F_fft_dif (A, c, 0, n);
  for (j = 0; j < c; j++)
    {
      F_mul_sqrt2exp (A[j], A[j], -l2, n);
      mpz_mul (gt, A[j], B[j]);
      F_mod_gt (A[j], n);
    }
  F_fft_dit (A, c, 0, n);
  
  
  /* The other rows */
  for (i = 1; i < r; i++)
    {
      int iomega = scramble[i] * ((4 * n) / (r * c)), ijomega;

      /* Multiply by \omega^{ij}, divide by transform length */
      for (j = 1, ijomega = iomega - l2; j < c; j++, ijomega += iomega)
        F_mul_sqrt2exp (A[i * c + j], A[i * c + j], ijomega, n);

      F_fft_dif (A + i * c, c, 0, n);

      for (j = 0; j < c; j++)
        {
          mpz_mul (gt, A[i * c + j], B[i * c + j]);
          F_mod_gt (A[i * c + j], n);
        }
      
      /* IDFT on rows */
      F_fft_dit (A + i * c, c, 0, n);
          
      /* Divide by \omega^{ij} */
      for (j = 1, ijomega = 4 * n - iomega; j < c; j++, ijomega -= iomega)
        F_mul_sqrt2exp (A[i * c + j], A[i * c + j], ijomega, n);
    }

      
  /* IDFT on columns */
  for (i = 0; i < c; i++)
    F_fft_dit (A + i, r, s2, n);

  free (scramble);
}

#define A0 A[i]
#define A1 A[l+i]
#define A2 A[2*l+i]
#define A3 A[3*l+i]
#define B0 B[i]
#define B1 B[l+i]
#define B2 B[2*l+i]
#define B3 B[3*l+i]
#define C0 C[i]
#define C1 C[l+i]
#define C2 C[2*l+i]
#define C3 C[3*l+i]
#define C4 C[4*l+i]
#define C5 C[5*l+i]
#define C6 C[6*l+i]
#define C7 C[7*l+i]
#define t0 t[i]
#define t1 t[l+i]
#define t2 t[2*l+i]
#define t3 t[3*l+i]
#define t4 t[4*l+i]
#define t5 t[5*l+i]


unsigned int
F_toomcook4 (mpz_t *C, mpz_t *A, mpz_t *B, unsigned int len, unsigned int n, 
             mpz_t *t)
{
  unsigned int l, i, r;

#ifdef DEBUG
  if (len % 4 != 0)
    {
      fprintf (stderr, "F_toomcook4: len not a multiple of 4.\n");
      exit (EXIT_FAILURE);
    }
#endif

#ifdef TESTDRIVE
  if (verbose)
    printf ("Toom-Cook, len = %d\n", len);
#endif
  
  l = len / 4;
  
  if (A == B) /* Squaring. The interpolation could probably be optimized, too */
    {
      for (i = 0; i < l; i++)
        {
          /*** Evaluate A(2), A(-2), 8*A(1/2) ***/
          mpz_mul_2exp (t0, A0, 1);
          mpz_add (t0, t0, A1);
          mpz_mul_2exp (t0, t0, 1);
          mpz_add (t0, t0, A2);
          mpz_mul_2exp (t0, t0, 1);
          mpz_add (t0, t0, A3);         /* t[0 .. l-1] = 8*A(1/2) < 15*N */
          F_mod_1 (t0, n);

          mpz_mul_2exp (t2, A3, 2);
          mpz_add (t2, t2, A1);
          mpz_mul_2exp (t2, t2, 1);     /* t[2l .. 3l-1] = 8*A_3 + 2*A_1 */
          mpz_mul_2exp (gt, A2, 2);
          mpz_add (gt, gt, A0);         /* gt = 4*A_2 + A0 */
          mpz_sub (t4, gt, t2);         /* t[4l .. 5l-1] = A(-2) */
          mpz_add (t2, t2, gt);         /* t[2l .. 3l-1] = A(2) */
          F_mod_1 (t4, n);
          F_mod_1 (t2, n);
          
          /* Evaluate A(1), A(-1) */
          mpz_add (C2, A0, A2);         /* May overwrite A2 */
          mpz_add (gt, A1, A3);
          mpz_sub (C4, C2, gt);         /* C4 = A(-1) */
          mpz_add (C2, C2, gt);         /* C2 = A(1) < 4*N */
          F_mod_1 (C2, n);
          F_mod_1 (C4, n);

#ifdef TESTDRIVE
          if (verbose >= 2)
            {
              gmp_printf ("A(0)[%d] = %Zd\n", i, C0);
              gmp_printf ("A(1)[%d] = %Zd\n", i, C2);
              gmp_printf ("A(-1)[%d] = %Zd\n", i, C4);
              gmp_printf ("A(inf)[%d] = %Zd\n", i, A3);
              gmp_printf ("8*A(1/2)[%d] = %Zd\n", i, t0);
              gmp_printf ("A(2)[%d] = %Zd\n", i, t2);
              gmp_printf ("A(-2)[%d] = %Zd\n", i, t4);
            }
#endif
        }

    /* A0  A1   A2   A3                     */
    /* A0      A(1)  A3  A(-1)              */
    /* C0  C1   C2   C3   C4    C5   C6  C7 */

      r = F_mul (t, t, t, l, 0, n, t + 6 * l);
        /* t0 = (8*A(1/2)) ^ 2 = 64*C(1/2) */
      r += F_mul (t + 2 * l, t + 2 * l, t + 2 * l, l, 0, n, t + 6 * l);
        /* t2 = A(2) ^ 2 = C(2) */
      r += F_mul (t + 4 * l, t + 4 * l, t + 4 * l, l, 0, n, t + 6 * l);
        /* t4 = A(-2) ^ 2 = C(-2) */
      r += F_mul (C, A, A, l, 0, n, t + 6 * l);
        /* C0 = A(0) ^ 2 = C(0) */
      r += F_mul (C + 6 * l, A + 3 * l, A + 3 * l, l, 0, n, t + 6 * l);
        /* C6 = A(inf) ^ 2 = C(inf) */
      r += F_mul (C + 2 * l, C + 2 * l, C + 2 * l, l, 0, n, t + 6 * l);
        /* C2 = A(1) ^ 2 = C(1). May overwrite A3 */
      r += F_mul (C + 4 * l, C + 4 * l, C + 4 * l, l, 0, n, t + 6 * l);
        /* C4 = A(-1) ^ 2 = C(-1) */
    }
  else /* Multiply */
    {
      for (i = 0; i < l; i++)
        {
          /*** Evaluate A(2), A(-2), 8*A(1/2) ***/
          mpz_mul_2exp (t0, A0, 1);
          mpz_add (t0, t0, A1);
          mpz_mul_2exp (t0, t0, 1);
          mpz_add (t0, t0, A2);
          mpz_mul_2exp (t0, t0, 1);
          mpz_add (t0, t0, A3);         /* t[0 .. l-1] = 8*A(1/2) < 15*N */
          F_mod_1 (t0, n);

          mpz_mul_2exp (t2, A3, 2);
          mpz_add (t2, t2, A1);
          mpz_mul_2exp (t2, t2, 1);     /* t[2l .. 3l-1] = 8*A_3 + 2*A_1 */

          mpz_mul_2exp (gt, A2, 2);
          mpz_add (gt, gt, A0);         /* gt = 4*A_2 + A0 */
          mpz_sub (t4, gt, t2);         /* t[4l .. 5l-1] = A(-2) */
          mpz_add (t2, t2, gt);         /* t[2l .. 3l-1] = A(2) */
          F_mod_1 (t4, n);
          F_mod_1 (t2, n);
          
          /*** Evaluate B(2), B(-2), 8*B(1/2) ***/
          mpz_mul_2exp (t1, B0, 1);
          mpz_add (t1, t1, B1);
          mpz_mul_2exp (t1, t1, 1);
          mpz_add (t1, t1, B2);
          mpz_mul_2exp (t1, t1, 1);
          mpz_add (t1, t1, B3);         /* t[l .. 2l-1] = 8*B(1/2) */
          F_mod_1 (t1, n);

          mpz_mul_2exp (t3, B3, 2);
          mpz_add (t3, t3, B1);
          mpz_mul_2exp (t3, t3, 1);     /* t[3l .. 4l-1] = 8*B_3 + 2*B_1 */

          mpz_mul_2exp (gt, B2, 2);
          mpz_add (gt, gt, B0);         /* gt = 4*B_2 + B0 */
          mpz_sub (t5, gt, t3);         /* t[5l .. 6l-1] = B(-2) */
          mpz_add (t3, t3, gt);         /* t[3l .. 4l-1] = B(2) */
          F_mod_1 (t5, n);
          F_mod_1 (t3, n);

          /* Evaluate A(1), A(-1) */
          mpz_add (C2, A0, A2);         /* May overwrite A2 */
    #undef A2
          mpz_add (gt, A1, A3);
          mpz_set (C1, B0);             /* C1 = B(0) May overwrite A1 */
    #undef A1
          mpz_sub (C4, C2, gt);         /* C4 = A(-1). May overwrite B0 */
    #undef B0
          mpz_add (C2, C2, gt);         /* C2 = A(1) < 4*N */
          F_mod_1 (C2, n);
          F_mod_1 (C4, n);

          /* Evaluate B(1), B(-1) */
          mpz_add (gt, C1, B2);         /* B0 is in C1 */
          mpz_set (C6, A3);             /* C6 = A(inf) May overwrite B2 */
    #undef B2
          mpz_add (C3, B1, B3);         /* May overwrite A3 */
    #undef A3
          mpz_sub (C5, gt, C3);         /* C5 = B(-1). May overwrite B1 */
    #undef B1
          mpz_add (C3, gt, C3);         /* C3 = B(1) */
          F_mod_1 (C3, n);
          F_mod_1 (C5, n);

    #ifdef TESTDRIVE
          if (verbose >= 2)
            {
              gmp_printf ("A(0)[%d] = %Zd\n", i, C0);
              gmp_printf ("B(0)[%d] = %Zd\n", i, C1);
              gmp_printf ("A(1)[%d] = %Zd\n", i, C2);
              gmp_printf ("B(1)[%d] = %Zd\n", i, C3);
              gmp_printf ("A(-1)[%d] = %Zd\n", i, C4);
              gmp_printf ("B(-1)[%d] = %Zd\n", i, C5);
              gmp_printf ("A(inf)[%d] = %Zd\n", i, C6);
              gmp_printf ("B(inf)[%d] = %Zd\n", i, B3);
              gmp_printf ("8*A(1/2)[%d] = %Zd\n", i, t0);
              gmp_printf ("8*B(1/2)[%d] = %Zd\n", i, t1);
              gmp_printf ("A(2)[%d] = %Zd\n", i, t2);
              gmp_printf ("B(2)[%d] = %Zd\n", i, t3);
              gmp_printf ("A(-2)[%d] = %Zd\n", i, t4);
              gmp_printf ("B(-2)[%d] = %Zd\n", i, t5);
            }
    #endif
        }

    /* A0 A1   A2   A3   B0    B1   B2 B3 */
    /* A0 B0  A(1) B(1) A(-1) B(-1) A3 B3 */
    /* C0 C1   C2   C3   C4    C5   C6 C7 */

      r = F_mul (t, t, t + l, l, 0, n, t + 6 * l);
        /* t0 = 8*A(1/2) * 8*B(1/2) = 64*C(1/2) */
      r += F_mul (t + 2 * l, t + 2 * l, t + 3 * l, l, 0, n, t + 6 * l);
        /* t2 = A(2) * B(2) = C(2) */
      r += F_mul (t + 4 * l, t + 4 * l, t + 5 * l, l, 0, n, t + 6 * l);
        /* t4 = A(-2) * B(-2) = C(-2) */
      r += F_mul (C, A, C + l, l, 0, n, t + 6 * l);
        /* C0 = A(0)*B(0) = C(0) */
      r += F_mul (C + 2 * l, C + 2 * l, C + 3 * l, l, 0, n, t + 6 * l);
        /* C2 = A(1)*B(1) = C(1) */
      r += F_mul (C + 4 * l, C + 4 * l, C + 5 * l, l, 0, n, t + 6 * l);
        /* C4 = A(-1)*B(-1) = C(-1) */
      r += F_mul (C + 6 * l, C + 6 * l, B + 3 * l, l, 0, n, t + 6 * l);
        /* C6 = A(inf)*B(inf) = C(inf) */
    }
  
/* C(0)   C(1)   C(-1)  C(inf)  64*C(1/2)  C(2)   C(-2) */
/* C0,C1  C2,C3  C4,C5  C6,C7   t0,t1      t2,t3  t4,t5 */

  for (i = 0; i < 2 * l - 1; i++)
    {
      mpz_add (t0, t0, t2);             /* t0 = 65 34 20 16 20 34 65 */

      mpz_sub (gt, C2, C4);             /* gt = 2*C_odd(1) = 0 2 0 2 0 2 0 */
      mpz_add (C2, C2, C4);             /* C2 = 2*C_even(1) = 2 0 2 0 2 0 2 */
      F_divby2 (C2, C2, n);             /* C2 = C_even(1) */

      mpz_add (C4, t2, t4);             /* C4 = 2*C_even(2) */
      F_divby2 (C4, C4, n);             /* C4 = C_even(2) */
      mpz_sub (t4, t2, t4);             /* t4 = 2*C_odd(2) */
      F_divby2 (t4, t4, n);
      F_divby2 (t4, t4, n);             /* t4 = C_odd(2)/2 = C_1 + 4*C_3 + 16*C_5 */
      F_divby2 (t2, gt, n);             /* t2 = C_odd(1) */

      mpz_sub (t0, t0, gt);             /* t0 = 65 32 20 14 20 32 65 */
      mpz_mul_2exp (gt, gt, 4);
      mpz_sub (t0, t0, gt);             /* t0 = 65 0 20 -18 20 0 65 */

      mpz_add (gt, C0, C6);             /* gt = C_0 + C_6 */
      mpz_sub (C2, C2, gt);             /* C2 = C_2 + C_4 */
      mpz_sub (t0, t0, gt);             /* t0 = 64 0 20 -18 20 0 64 */
      mpz_mul_2exp (gt, gt, 5);         /* gt = 32*C_0 + 32*C_6 */
      F_divby2 (t0, t0, n);             /* t0 = 32 0 10 -9 10 0 32 */
      mpz_sub (t0, t0, gt);             /* t0 = 0 0 10 -9 10 0 0 */
      mpz_sub (t0, t0, C2);             /* t0 = 0 0 9 -9 9 0 0 */
      F_divby3_1 (t0, n);
      F_divby3_1 (t0, n);               /* t0 = 0 0 1 -1 1 0 0 */
      mpz_sub (t0, C2, t0);             /* t0 = C_3 */
      mpz_sub (t2, t2, t0);             /* t2 = C_1 + C_5 */
      mpz_mul_2exp (gt, t0, 2);         /* gt = 4*C_3 */
      mpz_sub (t4, t4, gt);             /* t4 = C_1 + 16*C_5 */
      mpz_sub (t4, t4, t2);             /* t4 = 15*C_5 */
      F_divby3_1 (t4, n);
      F_divby5_1 (t4, n);               /* t4 = C_5 */
      mpz_sub (t2, t2, t4);             /* t2 = C_1 */

      mpz_sub (C4, C4, C0);             /* C4 = 4*C_2 + 16*C_4 + 64*C_6 */
      F_divby2 (C4, C4, n);
      F_divby2 (C4, C4, n);             /* C4 = C_2 + 4*C_4 + 16*C_6 */

      mpz_mul_2exp (gt, C6, 4);
      mpz_sub (C4, C4, gt);             /* C4 = C_2 + 4*C_4 */

      mpz_sub (C4, C4, C2);             /* C4 = 3*C_4 */
      F_divby3_1 (C4, n);               /* C4 = C_4 */
      mpz_sub (C2, C2, C4);             /* C2 = C_2 */

#ifdef TESTDRIVE
      if (verbose >= 2)
        {
          gmp_printf ("C0[%d] = %Zd\n", i, C0);
          gmp_printf ("C1[%d] = %Zd\n", i, t2);
          gmp_printf ("C2[%d] = %Zd\n", i, C2);
          gmp_printf ("C3[%d] = %Zd\n", i, t0);
          gmp_printf ("C4[%d] = %Zd\n", i, C4);
          gmp_printf ("C5[%d] = %Zd\n", i, t4);
          gmp_printf ("C6[%d] = %Zd\n", i, C6);
        }
#endif
    }

  for (i = 0; i < l - 1; i++)
    {
      mpz_add (C1, C1, t2);
      F_mod_1 (C1, n);
    }
  mpz_set (C1, t2);
  F_mod_1 (C1, n);
  for (i = l; i < 2 * l - 1; i++)
    {
      mpz_add (C1, C1, t2);
      F_mod_1 (C1, n);
    }
  
  for (i = 0; i < l - 1; i++)
    {
      mpz_add (C3, C3, t0);
      F_mod_1 (C3, n);
    }
  mpz_set (C3, t0);
  F_mod_1 (C3, n);
  for (i = l; i < 2 * l - 1; i++)
    {
      mpz_add (C3, C3, t0);
      F_mod_1 (C3, n);
    }

  for (i = 0; i < l - 1; i++)
    {
      mpz_add (C5, C5, t4);
      F_mod_1 (C5, n);
    }
  mpz_set (C5, t4);
  F_mod_1 (C5, n);
  for (i = l; i < 2 * l - 1; i++)
    {
      mpz_add (C5, C5, t4);
      F_mod_1 (C5, n);
    }

  return r;
}


/* Karatsuba split. Calls F_mul() to multiply the three pieces. */

unsigned int
F_karatsuba (mpz_t *R, mpz_t *A, mpz_t *B, unsigned int len, unsigned int n, 
             mpz_t *t)
{
  unsigned int i, r;
#ifdef DEBUG
  if (len % 2 != 0)
    {
      fprintf (stderr, "F_karatsuba: len is odd\n");
      exit (EXIT_FAILURE);
    }
#endif

#ifdef TESTDRIVE
  if (verbose)
    printf ("Karatsuba, len=%d\n", len);
#endif

  len /= 2;

  if (A == B) /* Squaring */
    {
      r = F_mul (t, A, A + len, len, 0, n, t + 2 * len); /* A0 * A1 */
      r += F_mul (R + 2 * len, A + len, A + len, len, 0, n, t + 2 * len); /* A1^2 */
      r += F_mul (R, A, A, len, 0, n, t + 2 * len); /* A0^2 */
      for (i = 0; i < 2 * len - 1; i++)
        {
          mpz_mul_2exp (t[i], t[i], 1);
          mpz_add (R[i + len], R[i + len], t[i]); /* i==len could be a mpz_set */
        }
      return r;
    }
  
  for (i = 0; i < len; i++) 
    {
      mpz_add (t[i],       A[i], A[i + len]); /* t0 = A0 + A1 */
      mpz_add (t[i + len], B[i], B[i + len]); /* t1 = B0 + B1 */
    }
  
  r = F_mul (t, t, t + len, len, 0, n, t + 2 * len);
  /* t[0...2*len-1] = (A0+A1) * (B0+B1) = A0*B0 + A0*B1 + A1*B0 + A1*B1 */
  
  if (R != A)
    {
      r += F_mul (R, A, B, len, 0, n, t + 2 * len);
      /* R[0...2*len-1] = A0 * B0 */
      r += F_mul (R + 2 * len, A + len, B + len, len, 0, n, t + 2 * len);
      /* R[2*len...4*len-1] = A1 * B1, may overwrite B */
    }
  else if (R + 2 * len != B)
    {
      r += F_mul (R + 2 * len, A + len, B + len, len, 0, n, t + 2 * len);
      /* R[2*len...4*len-1] = A1 * B1 */
      r += F_mul (R, A, B, len, 0, n, t + 2 * len);
      /* R[0...2*len-1] = A0 * B0, overwrites A */
    }
  else /* R == A && R + 2*len == B */
    {
      for (i = 0; i < len; i++)
        { /* mpz_swap instead? Perhaps undo later? Or interface for F_mul
             to specify separate result arrays for high/low half? (Would 
             make MFA awful!) */
          mpz_set (gt, A[len + i]); /* Swap A1 and B0 */
          mpz_set (A[len + i], B[i]);
          mpz_set (B[i], gt);
        }
      r += F_mul (R, R, R + len, len, 0, n, t + 2 * len);
      /* R[0...2*len-1] = A0 * B0, overwrites A */
      r += F_mul (R + 2 * len, R + 2 * len, R + 3 * len, len, 0, n, t + 2 * len);
      /* R[2*len...4*len-1] = A1 * B1, overwrites B */
    }

  /* R[0...2*len-2] == A0*B0, R[2*len-1] == 0 */
  /* R[2*len...3*len-2] == A1*B1, R[4*len-1] == 0 */
  /* t[0...2*len-2] == (A0+A1)*(B0+B1), t[2*len-1] == 0 */

  /* We're doing indices i and i+len in one loop on the assumption
     that 6 residues will probably fit into cache. After all,
     Karatsuba is only called for smallish F_m. This way, the final
     add R[i+len] += t[i] can be done inside the same loop and we need
     only one pass over main memory. */

  for (i = 0; i < len - 1; i++) 
    {
      mpz_sub (t[i], t[i], R[i]); /* t = A0*B1 + A1*B0 + A1*B1 */
      mpz_sub (t[i], t[i], R[i + 2 * len]); /* t = A0*B1 + A1*B0 */
      mpz_sub (t[i + len], t[i + len], R[i + len]);
      mpz_sub (t[i + len], t[i + len ], R[i + 3 * len]);
      
      mpz_add (R[i + len], R[i + len], t[i]);
      mpz_add (R[i + 2 * len], R[i + 2 * len], t[i + len]);
    }
  mpz_sub (t[len - 1], t[len - 1], R[len - 1]);
  mpz_sub (R[2 * len - 1], t[len - 1], R[3 * len - 1]);
  
  return r;
}

/* Multiply two polynomials with coefficients modulo 2^(2^m)+1. */
/* len is length (=degree+1) of polynomials and must be a power of 2. */
/* n=2^m */
/* Return value: number of multiplies performed */

unsigned int 
F_mul (mpz_t *R, mpz_t *A, mpz_t *B, unsigned int len, int monic, 
       unsigned int n, mpz_t *t) 
{
  unsigned int i, r=0, len2;
  int columns = 0;
#ifdef CHECKSUM
  mpz_t chksum1, chksum_1, chksum0, chksuminf;
#endif

  /* Handle trivial cases */
  if (len == 0)
    return 0;
  
  if (!gt_inited)
    {
      mpz_init2 (gt, 2 * n);
      gt_inited = 1;
    }
  
  if (len == 1)
    {
      mpz_mul (gt, A[0], B[0]);
      if (monic) /* (x + a0)(x + b0) = x^2 + (a0 + b0)x + a0*b0 */
        {
          mpz_add (R[1], A[0], B[0]); /* May overwrite B[0] */
          /* We don't store the leading monomial in the result poly */
        }
      else
        mpz_set_ui (R[1], 0);
      F_mod_gt (R[0], n); /* May overwrite A[0] */
      
      return 1;
    }

#ifdef CHECKSUM
  mpz_init2 (chksum1, n+64);
  mpz_init2 (chksum_1, n+64);
  mpz_init2 (chksum0, n+64);
  mpz_init2 (chksuminf, n+64);

  mpz_set_ui (gt, 0);
  for (i = 0; i < len; i++) 
    {
      /* Compute A(1) and B(1) */
      mpz_add (chksum1, chksum1, A[i]);
      mpz_add (gt, gt, B[i]);

      /* Compute A(-1) and B(-1) */
      if (i % 2 == 0)
        {
          mpz_add (chksum_1, chksum_1, A[i]);
          mpz_add (chksum0, chksum0, B[i]); /* chksum0 used temporarily here */
        }
      else
        {
          mpz_sub (chksum_1, chksum_1, A[i]);
          mpz_sub (chksum0, chksum0, B[i]);
        }
    }

  if (monic)
    {
      mpz_add_ui (chksum1, chksum1, 1);
      mpz_add_ui (gt, gt, 1);
      mpz_add_ui (chksum_1, chksum_1, 1);
      mpz_add_ui (chksum0, chksum0, 1);
    }
  
  mpz_mul (gt, gt, chksum1);
  F_mod_gt (chksum1, n);

  mpz_mul (gt, chksum0, chksum_1);
  F_mod_gt (chksum_1, n);

  /* Compute A(0) * B(0) */
  mpz_mul (gt, A[0], B[0]);
  F_mod_gt (chksum0, n);

  /* Compute A(inf) * B(inf) */
  mpz_mul (gt, A[len - 1], B[len - 1]);
  F_mod_gt (chksuminf, n);
  if (monic)
    {
      mpz_add (chksuminf, chksuminf, A[len - 2]);
      mpz_add (chksuminf, chksuminf, B[len - 2]);
    }

  r += 4;
#endif /* CHECKSUM */

  /* Don't do FFT if len=<4 (Karatsuba or Toom-Cook are faster), or if
     len > 2*n (no suitable primitive roots of 1) */
  if ((len > 4 || 1) && len <= 2 * n
#ifdef TESTDRIVE
   && !nofft
#endif
     ) 
    {
      /* len2 = log_2(len). Assumes len > 0 */
      for (i = len, len2 = 0; (i&1) == 0; i >>= 1, len2++);
      
      if (i != 1) 
        {
          fprintf (stderr, "F_mul: polynomial length must be power of 2, "
                           "but is %d\n", len);
          exit (EXIT_FAILURE);
        }
      
      /* If we use MFA, make a row use at most half the cache size, */
      /* Don't use MFA if full transform fits in half of L2 cache. */
      columns = 0;
      if (do_mfa && CACHESIZE * 8192 / 2 <= 2 * len * n)
        {
          columns = min(1U<<((len2 + 1) / 2), (CACHESIZE * 8192 / 2) / n);
        }

      /* Are we performing a squaring or multiplication? */
      if (A != B) 
        {
          /* So it's a multiplication */
          
#ifdef TESTDRIVE
          if (do_timing)
            rdtscll (timer_start);
#endif

          /* Put transform of B into t */
          for (i = 0; i < len; i++)
            mpz_set (t[i], B[i]);
          if (monic)
            mpz_set_ui (t[i++], 1);
          for (; i < 2 * len; i++)
            mpz_set_ui (t[i], 0);

          if (columns)
            F_fft_mfa (t, (2 * len) / columns, columns, 1, n);
          else
            F_fft_dif (t, 2 * len, 0, n);

#ifdef TESTDRIVE
          if (do_timing) 
            {
              rdtscll (timer_stop);
              timer_stop -= timer_start;
              printf("FFT(B), length %d: %.5fM clocks\n",  2*len, timer_stop/1000000.0);
            }
          
          if (verbose >= 2)
            printpoly (t, 2 * len, "FFT(B)");
#endif

        } else
          t = R; /* Do squaring */

#ifdef TESTDRIVE
      if (do_timing) 
        rdtscll(timer_start);
#endif

      /* Put A into R */
      for (i = 0; i < len; i++) 
        mpz_set (R[i], A[i]);
      if (monic)
        mpz_set_ui (R[i++], 1); /* May overwrite B[0] */
      for (; i < 2 * len; i++)
        mpz_set_ui (R[i], 0); /* May overwrite B[i - len] */

#ifndef MFA_MUL

      if (columns)
        F_fft_mfa (R, (2 * len) / columns, columns, 1, n);
      else
        F_fft_dif (R, 2 * len, 0, n);

#ifdef TESTDRIVE
      if (do_timing) 
        {
          rdtscll (timer_stop);
          timer_stop -= timer_start;
          printf ("FFT%s(A), length %d: %.5fM clocks\n", 
                 R==A?"":"2", 2 * len, timer_stop / 1000000.0);
        }
      
      if (verbose >= 2)
        printpoly (R, 2 * len, "FFT(A)");
      
      if (do_timing)
        rdtscll (timer_start);
#endif
      
      for (i = 0; i < 2 * len; i++) 
        {
/*          printf ("Before multiply: mpz_size (R[%d]) = %d\n", i, mpz_size (R[i])); */
          F_mod_1 (R[i], n);
/*          printf ("Before multiply, reduced: mpz_size (R[%d]) = %d\n", i, mpz_size (R[i])); */
          mpz_mul (gt, R[i], t[i]);
/*          printf ("After multiply: mpz_size (R[%d]) = %d\n", i, mpz_size (gt)); */
          F_mod_gt (R[i], n);
/*          printf ("After multiply, reduced: mpz_size (R[%d]) = %d\n", i, mpz_size (R[i])); */
          /* Do the div-by-length. Transform length was len*2, so divide by
             2^(len2+1) = sqrt(2)^(2*(len2+1)) */

          F_mul_sqrt2exp (R[i], R[i], - 2 * (len2 + 1), n);
        }

      r += 2 * len;

#ifdef TESTDRIVE
      if (do_timing) 
        {
          rdtscll (timer_stop);
          timer_stop -= timer_start;
          printf ("%d Multiplies: %.5fM clocks\n", 2 * len, timer_stop / 1000000.0);
          fflush (stdout);
          rdtscll (timer_start);
        }
#endif
      if (columns)
        F_fft_mfa (R, (2 * len) / columns, columns, -1, n);
      else
        F_fft_dit (R, 2 * len, 0, n);

#ifdef TESTDRIVE
      if (do_timing) 
        {
          rdtscll (timer_stop);
          timer_stop -= timer_start;
          printf ("IFFT, length %d: %.5fM clocks (%.0f per butterfly)\n", 
            2*len, timer_stop / 1000000.0, timer_stop / (2.*len*(len2+1)));
        }
#endif

#else /* MFA_MUL */
      
      if (columns)
        F_fft_mfa_mul (R, t, (2 * len) / columns, columns, n);
      else
        {
          F_fft_dif (R, 2 * len, 0, n);
          for (i=0; i < 2 * len; i++)
           {
             F_mod_1 (R[i], n);
             mpz_mul (gt, R[i], t[i]);
             F_mod_gt (R[i], n);
             F_mul_sqrt2exp (R[i], R[i], - 2 * (len2 + 1), n);
           }
          F_fft_dit (R, 2 * len, 0, n);
        }
      r += 2 * len;
      

#ifdef TESTDRIVE
      if (do_timing) 
        {
          rdtscll (timer_stop);
          timer_stop -= timer_start;
          printf ("Time for F_fft_mfa_mul, length %d: %.5fM clocks\n", 
            2*len, timer_stop / 1000000.0);
        }
#endif
      
#endif /* MFA_MUL */

      if (monic)
        {
          mpz_sub_ui (R[0], R[0], 1);
        }
      
    } else { /* Karatsuba or Toom-Cook split */

      if (len / n == 4 || len == 2 
#ifdef TESTDRIVE
          || force_karatsuba
#endif
         )
        r += F_karatsuba (R, A, B, len, n, t);
      else
        r += F_toomcook4 (R, A, B, len, n, t);

      if (monic) /* Handle the leading monomial the hard way */
        {
          /* This only works if A, B and R do not overlap */
          if (A == R || B == R + len)
            {
              fprintf (stderr, "F_mul: monic polynomials with Karatsuba/"
                    "Toom-Cook and overlapping input/output not supported\n");
              exit (EXIT_FAILURE);
            }
          for (i = 0; i < len; i++)
            {
              mpz_add (R[i + len], R[i + len], A[i]);
              mpz_add (R[i + len], R[i + len], B[i]);
              F_mod_1 (R[i + len], n);
            }
        }
    }
      
#ifdef DEBUG
  F_mod_1 (R[2 * len - 1], n);
  if (!monic && mpz_sgn (R[2 * len - 1]) != 0)
    {
      fprintf (stderr, "F_mul, len %d: R[%d] == ", len, 2 * len - 1);
      mpz_out_str (stderr, 10, R[2 * len - 1]);
      fprintf (stderr, " != 0\n");
    }
#endif

#ifdef CHECKSUM
  /* Compute R(1) = (A*B)(1) and subtract from chksum1 */

  for (i = 0; i < 2*len; i++) 
    mpz_sub (chksum1, chksum1, R[i]);
  
  if (monic)
    mpz_sub_ui (chksum1, chksum1, 1);

  while (mpz_sizeinbase (chksum1, 2) > n) 
    F_mod_1 (chksum1, n);

  if (mpz_sgn (chksum1) != 0) 
    {
      printf("F_mul, len %d: A(1)*B(1) != R(1), difference ", len);
      mpz_out_str(stdout, 10, chksum1);
      printf("\n");
    }

  for (i = 0; i < 2 * len; i++) 
    if (i % 2 == 0)
      mpz_sub (chksum_1, chksum_1, R[i]);
    else
      mpz_add (chksum_1, chksum_1, R[i]);
  
  if (monic)
    mpz_sub_ui (chksum_1, chksum_1, 1);
  
  while (mpz_sizeinbase (chksum_1, 2) > n) 
    F_mod_1 (chksum_1, n);

  if (mpz_sgn (chksum_1) != 0) 
    {
      printf("F_mul, len %d: A(-1)*B(-1) != R(-1), difference ", len);
      mpz_out_str(stdout, 10, chksum_1);
      printf("\n");
    }

  mpz_sub (chksum0, chksum0, R[0]);
  while (mpz_sizeinbase (chksum0, 2) > n) 
    F_mod_1 (chksum0, n);

  if (mpz_sgn (chksum0) != 0) 
    {
      printf("F_mul, len %d: A(0)*B(0) != R(0), difference ", len);
      mpz_out_str(stdout, 10, chksum0);
      printf("\n");
    }

  mpz_sub (chksuminf, chksuminf, R[2 * len - 2]);
  while (mpz_sizeinbase (chksuminf, 2) > n) 
    F_mod_1 (chksuminf, n);

  if (mpz_sgn (chksuminf) != 0) 
    {
      printf ("F_mul, len %d: A(inf)*B(inf) != R(inf), difference ", len);
      mpz_out_str (stdout, 10, chksuminf);
      printf ("\n");
    }

  mpz_clear (chksum1);
  mpz_clear (chksum_1);
  mpz_clear (chksum0);
  mpz_clear (chksuminf);
#endif /* CHECKSUM */

  return r;
}

/* Transposed multiply of two polynomials with coefficients 
   modulo 2^(2^m)+1.
   len is length of polynomial B and must be a power of 2,
   the length of polynomial A is len / 2. 
   n=2^m
   t must have space for 2*len coefficients 
   Only the product coefficients [len / 2 ... len - 1] will go into 
   R[0 ... len / 2 - 1] 
   Return value: number of multiplies performed */

unsigned int 
F_mul_trans (mpz_t *R, mpz_t *A, mpz_t *B, unsigned int len, unsigned int n, 
             mpz_t *t) 
{
  unsigned int i, r = 0, len2;
  int columns = 0;

/*  printf ("F_mul_trans: R=%p, A=%p, B=%p, len=%d, n=%d, t=%p\n", 
          R, A, B, len, n, t);
*/

  /* Handle trivial cases */
  if (len < 2)
    return 0;
  
  if (!gt_inited)
    {
      mpz_init2 (gt, 2 * n);
      gt_inited = 1;
    }
  
  if (len == 2)
    {
      mpz_mul (gt, A[0], B[0]);
      F_mod_gt (R[0], n);
      return 1;
    }

  if (len <= 4 * n)
    {
      /* len2 = log_2(len) */
      for (i = len, len2 = 0; i > 1 && (i&1) == 0; i >>= 1, len2++);
      
      if (i != 1) 
        {
          fprintf (stderr, "F_mul_trans: polynomial length must be power of 2, "
                           "but is %d\n", len);
          exit (EXIT_FAILURE);
        }
      
      /* If we use MFA, make a row use at most half the cache size */
      if (do_mfa)
        columns = min(1U<<(len2 / 2), (CACHESIZE * 4096) / n);

      /* Put transform of B into t */
      for (i = 0; i < len; i++)
        mpz_set (t[i], B[i]);

      if (do_mfa)
        F_fft_mfa (t, len / columns, columns, 1, n);
      else
        F_fft_dif (t, len, 0, n);

      /* Put transform of reversed A into t + len */
      for (i = 0; i < len / 2; i++) 
        mpz_set (t[i + len], A[len / 2 - 1 - i]);
      for (i = len / 2; i < len; i++)
        mpz_set_ui (t[i + len], 0);

      if (do_mfa)
        F_fft_mfa (t + len, len / columns, columns, 1, n);
      else
        F_fft_dif (t + len, len, 0, n);

      for (i = 0; i < len; i++) 
        {
/*          printf ("Before multiply: mpz_size (R[%d]) = %d\n", i, mpz_size (R[i])); */
          F_mod_1 (t[i], n);
          F_mod_1 (t[i + len], n);
/*          printf ("Before multiply, reduced: mpz_size (R[%d]) = %d\n", i, mpz_size (R[i])); */
          mpz_mul (gt, t[i], t[i + len]);
/*          printf ("After multiply: mpz_size (R[%d]) = %d\n", i, mpz_size (gt)); */
          F_mod_gt (t[i], n);
/*          printf ("After multiply, reduced: mpz_size (R[%d]) = %d\n", i, mpz_size (R[i])); */

          /* Do the div-by-length. Transform length was len, so divide by
             2^len2 = sqrt(2)^(2*len2) */
          F_mul_sqrt2exp (t[i], t[i], - 2 * len2, n);
        }

      r += len;

      if (do_mfa)
        F_fft_mfa (t, len / columns, columns, -1, n);
      else
        F_fft_dit (t, len, 0, n);
      
      for (i = 0; i < len / 2; i++)
        mpz_set (R[i], t[i + len / 2 - 1]);

    } else { /* Only Karatsuba, no Toom-Cook here */
      unsigned int h = len / 4;
      
/*
      fprintf (stderr, "No transposed Karatsuba yet!\n");
      exit (EXIT_FAILURE);
      printf ("Transposed Karatsuba, len = %d\n", len);
*/
      
      /* A = a1 * x^h + a0
         B = b3 * x^3h + b2 * x^2h + b1 * x^h + b0
         mul^T(A, B) = mul^T(a0,b3) * x^4h + 
                      (mul^T(a1,b3) + mul^T(a0,b2)) * x^3h + 
                      (mul^T(a1,b2) + mul^T(a0,b1)) * x^2h + 
                      (mul^T(a1,b1) + mul^T(a0,b0)) * x + 
                       mul^T(a1,b0)
         We only want the x^h, x^2h and x^3h coefficients,
         mul^T(a1,b1) + mul^T(a0,b0)
         mul^T(a1,b2) + mul^T(a0,b1)
         mul^T(a1,b3) + mul^T(a0,b2)

         Specifically, we want R[i] = \sum_{j=0}^{len/2} A[j] * B[j+i], 0<=i<len/2
      */
      
      /* t */
      for (i = 0; i < h; i++)
        mpz_add (t[i], A[i], A[i + h]);
      F_mul_trans (t, t, B + h, 2 * h, n, t + h); /* t[0 ... h-1] = t */
                                                  /* Uses t[h ... 5h-1] as temp */

      /* t[i] = \sum_{j=0}^{h-1} (A[j]+A[j+h]) * B[j+i+h], 0 <= i < h */
      
      /* r */
      for (i = 0; i < 2 * h; i++)
        mpz_sub (t[i + h], B[i], B[h + i]);
      F_mul_trans (t + h, A, t + h, 2 * h, n, t + 3 * h); /* t[h ... 2h-1] = r */
                                                          /* Uses t[3h ... 7h-1] as temp */
      /* t[i + h] = \sum_{j=0}^{h-1} A[j] * (B[j+i] - B[j+i+h]), 0 <= i < h */
      
      for (i = 0; i < h; i++)
        mpz_add (R[i], t[i], t[i + h]); /* R[0 ... h-1] = t + r */
                                        /* r not needed anymore */
      
      /* R[i] = \sum_{j=0}^{h-1} (A[j]+A[j+h]) * B[j+i+h] +  \sum_{j=0}^{h-1} A[j] * (B[j+i] - B[j+i+h]) =
                \sum_{j=0}^{h-1} (A[j]+A[j+h]) * B[j+i+h] + A[j] * (B[j+i] - B[j+i+h]) = 
                \sum_{j=0}^{h-1} A[j]*B[j+i+h] + A[j+h]*B[j+i+h] + A[j]*B[j+i] - A[j]*B[j+i+h] =
                \sum_{j=0}^{h-1} A[j+h]*B[j+i+h] + A[j]*B[j+i] =
                \sum_{j=0}^{h-1} A[j]*B[j+i] + \sum_{j=0}^{h-1} A[j+h]*B[j+i+h] =
                \sum_{j=0}^{h-1} A[j]*B[j+i] + \sum_{j=h}^{2h-1} A[j]*B[j+i] =

                \sum_{j=0}^{2h-1} A[j]*B[j+i], 0 <= i < h   (1)
      */
      
      /* s */
      for (i = 0; i < 2 * h; i++)
        mpz_sub (t[i + h], B[i + 2 * h], B[i + h]);
      F_mul_trans (t + h, A + h, t + h, 2 * h, n, t + 3 * h); /* t[h ... 2h-1] = s */
                                                              /* Uses t[3h ... 7h - 1] as temp */

      /* t[i + h] = \sum_{j=0}^{h} A[j+h] * (B[j+i+2*h]-B[j+i+h]), 0 <= i < h */
      
      for (i = 0; i < h; i++)
        mpz_add (R[i + h], t[i], t[i + h]);
      
      /* R[i + h] = \sum_{j=0}^{h-1} (A[j]+A[j+h]) * B[j+i+h] + \sum_{j=0}^{h} A[j+h] * (B[j+i+2*h]-B[j+i+h]) =
                    \sum_{j=0}^{h-1} (A[j]+A[j+h]) * B[j+i+h] + A[j+h] * (B[j+i+2*h]-B[j+i+h]) =
                    \sum_{j=0}^{h-1} A[j]*B[j+i+h] + A[j+h]*B[j+i+h] + A[j+h]*B[j+i+2*h] - A[j+h]*B[j+i+h] =
                    \sum_{j=0}^{h-1} A[j]*B[j+i+h] + A[j+h]*B[j+i+2*h] =
                    \sum_{j=0}^{h-1} A[j]*B[j+i+h] + \sum_{j=0}^{h-1} A[j+h]*B[j+i+2*h] =
                    \sum_{j=0}^{h-1} A[j]*B[j+i+h] + \sum_{j=h}^{2h-1} A[j]*B[j+i+h] =
                    \sum_{j=0}^{2h-1} A[j]*B[j+i+h], 0 <= i < h

        R[i] = \sum_{j=0}^{2h-1} A[j]*B[j+i], h <= i < 2*h   (2)
        
        (1) and (2) : R[i] = \sum_{j=0}^{2h-1} A[j]*B[j+i], 0 <= i < 2*h
      */
          
    }
  
  return r;
}

#ifdef TESTDRIVE

int 
main(int argc, char **argv) 
{
  int i;
  int length = 16;
  int fermat = 5;
  int fermat_e;
  int monic = 0;
  unsigned int len2;
  unsigned long long timer_start, timer_stop;

  mpz_t *res, *poly1, *poly2, *t;
  mpz_t F;
 
  mp_set_memory_functions (&allocate_function, &reallocate_function,
                           &deallocate_function);

  do_mfa = 0;

  while (argc > 1 && argv[1][0] == '-')
    {
      switch (argv[1][1])
        {
          case 'v' : verbose++; break;
          case 't' : do_timing = 1; break;
          case 'm' : do_mfa = 1; break;
          case 'k' : force_karatsuba = 1; break;
          case 'n' : nofft = 1; break;
          case 'o' : monic = 1; break;
          case '2' : radix2 = 1; break;
          default :
            {
              fprintf (stderr, "Unknown option: %s\n", argv[1]);
              exit (EXIT_FAILURE);
            }
        }
      argc--;
      argv++;
    }
  
  /* First command line argument is length */
  if (argc > 1)
    length = atoi(argv[1]);
  
  /* Second command line argument is m (Fermat number F_m) */
  if (argc > 2)
    fermat = atoi(argv[2]);

  fermat_e = 1<<fermat;

  for (i = length, len2 = 0; i > 1 && (i&1) == 0; i >>= 1, len2++);

  mpz_init2 (gt, 2*fermat_e + 2 * GMP_NUMB_BITS);
  mpz_init2 (F, fermat_e + 2 * GMP_NUMB_BITS);
  mpz_set_ui (F, 1);
  mpz_mul_2exp (F, F, fermat_e);
  mpz_add_ui (F, F, 1);

  /* Quick mod reduction with F_mod_* produces values in [-2^(n+1), 2^(n+1)].
     This would fit into (n/GMP_NUMB_BITS)+1 limbs.
     However, mpz_{add|sub} and others allocate max(size(S1),size(S2))+1 limbs
     to the destination, which lets the destination grow to 
     (n/GMP_NUMB_BITS)+2 limbs */
  
  poly1 = poly_init (length, fermat_e + 2 * GMP_NUMB_BITS);
  poly2 = poly_init (length, fermat_e + 2 * GMP_NUMB_BITS);
  res = poly_init (2 * length, fermat_e + 2 * GMP_NUMB_BITS);
  t = poly_init (3 * length, fermat_e + 2 * GMP_NUMB_BITS);
  
  /* Test squaring */
  printf ("Squaring, not in-place\n");
  for (i=0; i<length; i++)
    mpz_set_ui (poly1[i], 1);
  printpoly (poly1, length, "Poly1");
  
  if (do_timing)
    rdtscll (timer_start);
    
  F_mul (res, poly1, poly1, length, monic, fermat_e, t);

  if (do_timing)
    {
      rdtscll (timer_stop);
      timer_stop -= timer_start;
      printf("Total: %.5fM clocks\n",  timer_stop/1000000.0);
    }

  printpoly (res, 2*length, "Res");


  /***** Test an in-place squaring *****/
  printf ("Squaring, in-place\n");
  for (i=0; i<length; i++)
    mpz_set_ui (res[i], 1);
  printpoly (res, length, "Poly1");
  
  if (do_timing)
    rdtscll (timer_start);
    
  F_mul (res, res, res, length, monic, fermat_e, t);

  if (do_timing)
    {
      rdtscll (timer_stop);
      timer_stop -= timer_start;
      printf("Total: %.5fM clocks\n",  timer_stop/1000000.0);
    }

  printpoly (res, 2*length, "Res");

  /* Test a multiply */
  printf ("Multiply, not in-place\n");
  for (i = 0; i < length; i++) 
    {
      mpz_set_ui (poly1[i], 1);
      mpz_set_ui (poly2[i], i);
    }

  print2polys (poly1, poly2, length, "Poly1", "Poly2");
  
  if (do_timing)
    rdtscll (timer_start);
    
  F_mul (res, poly1, poly2, length, monic, fermat_e, t);

  if (do_timing)
    {
      rdtscll (timer_stop);
      timer_stop -= timer_start;
      printf("Total: %.5fM clocks\n",  timer_stop/1000000.0);
    }

  printpoly (res, 2*length, "Res");

  /* Test an in-place multiply */
  printf ("Multiply, in-place\n");
  mpz_set_ui (res[0], 3);
  mpz_set_ui (res[length], 5);
  for (i = 0; i < 20; i++) 
    {
      mpz_mul (gt, res[0], res[0]);
      F_mod_gt (res[0], fermat_e);
      mpz_mul (gt, res[length], res[length]);
      F_mod_gt (res[length], fermat_e);
    }
  
/*  printf ("Inited poly start values\n"); */
  
  for (i = 1; i < length; i++) 
    {
      mpz_mul (gt, res[i - 1], res[i - 1]);
      F_mod_gt (res[i], fermat_e);
      mpz_mul (gt, res[length + i - 1], res[length + i - 1]);
      F_mod_gt (res[length + i], fermat_e);
    }

/*  printf ("Inited polys\n"); */

  for (i = 0; i < 2 * length; i++) 
    {
      if (mpz_sizeinbase (res[i], 2) > (unsigned) fermat_e + 1)
        printf ("Init: mpz_sizeinbase(res[%d], 2) == %d\n", i, mpz_sizeinbase (res[i], 2));
    }
  
  if (4 * fermat_e >= length)
    {
      /* Try if IFFT(FFT()) is identity */

      if (verbose >= 1)
        printf ("Testing if IFFT(FFT()) is identity, plain DIT/DIF\n");

      for (i = 0; i < 2 * length; i++)
        {
          F_mod_1 (res[i], fermat_e);
          if (mpz_sgn (res[i]) < 0)
            F_add_F (res[i], fermat_e);
          mpz_set (t[i], res[i]);
        }
      
/*      printf ("Mod reduction on res done\nHere come the FFTs\n"); */

      F_fft_dif (t, length, 0, fermat_e);
      if (verbose >= 2)
        for (i = 0; i < length; i++)
            gmp_printf ("FFT_DIF(A)[%d] = %Zd\n", i, t[i]);
      F_fft_dit (t, length, 0, fermat_e);
      F_fft_dif (t + length, length, 0, fermat_e);
      F_fft_dit (t + length, length, 0, fermat_e);

/*      printf ("Done FFTs\n"); */

      for (i = 0; i < 2 * length; i++)
        {
          F_mul_sqrt2exp (t[i], t[i], - 2 * len2, fermat_e);
          F_mod_1 (t[i], fermat_e);
        }
      
/*      printf ("Done div-by-length and mod reduction\n"); */

      for (i = 0; i < 2 * length; i++)
        {
          if (mpz_sgn(t[i]) < 0)
            F_add_F (t[i], fermat_e);
          if (mpz_cmp (t[i], res[i]) != 0)
            {
              printf ("res[%d] = ", i);
              mpz_out_str (stdout, 10, res[i]);
              printf ("\nIFFT_dit(FFT_dif(res))[%d] = ", i);
              mpz_out_str (stdout, 10, t[i]);
              printf ("\n");
            }
        }

      if (verbose >= 1) 
        printf ("Testing if IFFT(FFT()) is identity, MFA\n");

      for (i = 0; i < length; i++)
        {
          mpz_set (t[i], res[i]);
          mpz_set (t[length + i], res[length + i]);
        }

      {
        unsigned int rows, columns;
        
        columns = min(1U<<((len2 + 1) / 2), (CACHESIZE * 4096) / fermat_e);
        rows = length / columns;

        F_fft_mfa (t, rows, columns, 1, fermat_e);
        F_fft_mfa (t, rows, columns, -1, fermat_e);
        F_fft_mfa (t + length, rows, columns, 1, fermat_e);
        F_fft_mfa (t + length, rows, columns, -1, fermat_e);
      }

      for (i = 0; i < 2 * length; i++)
        {
          F_mul_sqrt2exp (t[i], t[i], - 2 * len2, fermat_e);
          F_mod_1 (t[i], fermat_e);
          if (mpz_sgn (t[i]) < 0)
            F_add_F (t[i], fermat_e);
          if (mpz_cmp (t[i], res[i]) != 0)
            {
              printf ("res[%d] = ", i);
              mpz_out_str (stdout, 10, res[i]);
              printf ("\nIFFT_mfa(FFT_mfa(res))[%d] = ", i);
              mpz_out_str (stdout, 10, t[i]);
              printf ("\n");
            }
        }
      
      if (verbose >= 1)
        printf ("Finished Id test\n");

      for (i = 0; i < 2 * length; i++) 
        {
          if (mpz_sizeinbase (res[i], 2) > (unsigned) fermat_e + 1)
            printf ("After Id test: mpz_sizeinbase(res[%d], 2) == %d\n", i, mpz_sizeinbase (res[i], 2));
        }
  } else if (verbose >= 1)
    printf ("Could not do identity test, length too large for direct transform\n");

  print2polys (res, res+length, length, "Poly1", "Poly2");
  
  if (do_timing)
    rdtscll (timer_start);
    
  F_mul (res, res, res+length, length, monic, fermat_e, t);
  
  if (do_timing)
    {
      rdtscll (timer_stop);
      timer_stop -= timer_start;
      printf("Total: %.5fM clocks\n",  timer_stop/1000000.0);
    }

  printpoly(res, 2*length, "Res");

  for (i = 0; i < length; i++) 
    if (poly1[i]->_mp_alloc > fermat_e / 32 + 2)
      printf ("poly1[%d]._mp_alloc = %u\n", i, poly1[i]->_mp_alloc);

  for (i = 0; i < length; i++) 
    if (poly2[i]->_mp_alloc > fermat_e / 32 + 2)
      printf ("poly2[%d]._mp_alloc = %u\n", i, poly2[i]->_mp_alloc);

  for (i = 0; i < 2 * length; i++) 
    if (res[i]->_mp_alloc > fermat_e / 32 + 2)
      printf ("res[%d]._mp_alloc = %u\n", i, res[i]->_mp_alloc);

  for (i = 0; i < 2 * length; i++) 
    if (t[i]->_mp_alloc > fermat_e / 32 + 2)
      printf ("t[%d]._mp_alloc = %u\n", i, t[i]->_mp_alloc);

  return 0;
}
#endif
