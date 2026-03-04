/* Batched S algorithm

Copyright 2026 Seth Troisi

This file is part of the ECM Library.

The ECM Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The ECM Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the ECM Library; see the file COPYING.LIB.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#include "batched_s.h"
#include "ecm_int.h"
#include "ecm-gmp.h"
#include <math.h>
#include <gmp.h>

uint64_t SMALL_PRIME[] =
{
  2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
  73, 79, 83
};


void
batched_info_init (batched_info_t info)
{
    info->B1 = 1;
    for (int p_i = 0; p_i < SMALL_PRIMES; p_i++) {
        uint64_t p = SMALL_PRIME[p_i];
        uint64_t pp = p;
        /* by definition this shouldn't overflow in setup */
        for (int j = 1; j <= BATCHED_ITERATORS; j++)
            pp *= p;
        /* pp = p ^ (ITERATIONS+1) */
        info->small_q[p_i] = pp - 1;
    }
    for (int i = 0; i < BATCHED_ITERATORS; i++)
      {
        prime_info_init(info->iterator[i].p_i);
        info->iterator[i].current_prime = 2;
      }
}


void
batched_info_clear (batched_info_t info)
{
    for (int i = 0; i < BATCHED_ITERATORS; i++)
      {
        prime_info_clear(info->iterator[i].p_i);
      }
}

/* Compare m with t^n
   Return -1 if less than, 0 if equal, 1 if greater */
static int
cmp_power(uint64_t m, uint64_t t, uint64_t n)
{
    ASSERT_ALWAYS(n >= 2);

    uint64_t pp = t;
    for (unsigned int p = 2; p < n; p++)
      pp *= t;

    /* Avoid overflow by comparing m/t ? t^(n-1) */
    uint64_t q = m / t;
    if (q > pp) {
        return 1;
    }
    if (q < pp) {
        return -1;
    }
    /* equal if m % t == 0, else m > t^n */
    if (m - q * t == 0) {
        return 0;
    }
    return 1;
}

/* integer nth_root (m ^ (1 / n)) */
static uint64_t
floor_nth_root (uint64_t m, uint64_t n)
{
  if (n == 1)
    return m;

  uint64_t t = pow((double) m, 1.0 / n);
  // check if t^n >= m
  int sgn = cmp_power(m, t, n);
  if (sgn == 0)
      return t;
  while (1)
    {
      t += sgn;
      int new_sgn = cmp_power(m, t, n);
      if (new_sgn == 0)
          return t;
      if (new_sgn != sgn)
        {
          return t - (new_sgn == -1);
        }
    }
}

/*** Cascaded multiply ***/

#define CASCADE_THRES 3
#define CASCADE_MAX 50000000.0

typedef struct {
  unsigned int size;
  mpz_t *val;
} mul_casc;

/* return NULL if an error occurred */
static mul_casc *
mulcascade_init (void)
{
  mul_casc *t;

  t = (mul_casc *) malloc (sizeof (mul_casc));
  ASSERT_ALWAYS(t != NULL);

  t->size = 2;
  t->val = (mpz_t*) malloc (t->size * sizeof (mpz_t));
  ASSERT_ALWAYS(t->val != NULL);
  for (unsigned int i = 0; i < t->size; i++)
    mpz_init (t->val[i]);
  return t;
}

static void
mulcascade_free (mul_casc *c)
{
  unsigned int i;

  for (i = 0; i < c->size; i++)
    mpz_clear (c->val[i]);
  free (c->val);
  free (c);
}

static void
mulcascade_mul (mul_casc *c, const uint64_t n)
{
  unsigned int i;

  // Always free c->val[0] but don't continue if c->val[1] is small
  ASSERT_ALWAYS (mpz_sgn (c->val[0]) == 0);
  mpz_set_uint64 (c->val[0], n);

  if (mpz_sgn (c->val[1]) == 0)
    {
      mpz_swap(c->val[0], c->val[1]);
      return;
    }

  mpz_mul(c->val[1], c->val[1], c->val[0]);
  mpz_set_ui(c->val[0], 0);
  if (mpz_size (c->val[1]) <= CASCADE_THRES)
    return;

  for (i = 2; i < c->size; i++)
    {
      if (mpz_sgn (c->val[i]) == 0)
        {
          mpz_set (c->val[i], c->val[i-1]);
          mpz_set_ui (c->val[i-1], 0);
          return;
        }
      else
        {
          mpz_mul (c->val[i], c->val[i], c->val[i-1]);
          mpz_set_ui (c->val[i-1], 0);
        }
    }

  /* Allocate more space for cascade */

  i = c->size++;
  c->val = (mpz_t*) realloc (c->val, c->size * sizeof (mpz_t));
  ASSERT_ALWAYS(c->val != NULL);
  mpz_init (c->val[i]);
  mpz_swap (c->val[i], c->val[i-1]);
}

/* pm1.c had mulcascade_set to initialize with an existing mpz_t */
// static void mulcascade_set (mul_casc *c, mpz_t n)

static void
mulcascade_get_z_with_clear (mpz_t r, mul_casc *c)
{
  unsigned int i;

  ASSERT(c->size != 0);

  mpz_set_ui (r, 1);

  for (i = 0; i < c->size; i++)
    if (mpz_sgn (c->val[i]) != 0)
      {
        mpz_mul (r, r, c->val[i]);
        mpz_set_ui (c->val[i], 0);
      }
}

void
advance_to (batched_info_t info, uint64_t B1done)
{
    // Advance up to B1done, no need to compute s.
    for (int p_i = 0; p_i < SMALL_PRIMES; p_i++)
      {
        uint64_t pp = info->small_q[p_i];
        /* if next power of small prime is <= B1done */
        /* < is correct so that UINT_MAX can be sentinel for overflow */
        while (pp < B1done)
          {
            uint64_t p = SMALL_PRIME[p_i];
            // if pp * p will overflow set as max value.
            uint64_t t = ECM_UINT_MAX / p;
            pp += 1;
            pp = (t >= pp) ? (pp * p - 1) : ECM_UINT_MAX;
          }
        info->small_q[p_i] = pp;
      }

    for (int j = 1; j <= BATCHED_ITERATORS; j++)
      {
        struct prime_iterator_s *it = &info->iterator[j-1];
        uint64_t max_p = floor_nth_root (B1done, j);
        uint64_t p = it->current_prime;
        while (p <= max_p)
          {
              p = getprime_mt (it->p_i);
          }
        it->current_prime = p;
      }
    info->B1 = B1done;
}


void
get_batch (batched_info_t info, mpz_t s, uint64_t B1_next)
{
    mpz_set_ui(s, 1);
    /* TODO could save this between calls. */
    mul_casc *cascade = mulcascade_init ();

    for (int p_i = 0; p_i < SMALL_PRIMES; p_i++)
      {
        uint64_t pp = info->small_q[p_i];
        /* if next power of small prime is <= B1_next add a multiple of p */
        /* < is correct so that UINT_MAX can be sentinel for overflow */
        while (pp < B1_next)
          {
            uint64_t p = SMALL_PRIME[p_i];
            mulcascade_mul (cascade, p);

            // if pp * p will overflow set as max value.
            uint64_t t = ECM_UINT_MAX / p;
            pp += 1;
            pp = (t >= pp) ? (pp * p - 1) : ECM_UINT_MAX;
          }
        info->small_q[p_i] = pp;
      }

    /* Handle primes^{2,3,4,5} */
    for (int j = 2; j <= BATCHED_ITERATORS; j++)
      {
        struct prime_iterator_s *it = &info->iterator[j-1];
        uint64_t max_p = floor_nth_root (B1_next, j);
        uint64_t p = it->current_prime;
        while (p <= max_p)
          {
              mulcascade_mul (cascade, p);
              p = getprime_mt (it->p_i);
          }
        it->current_prime = p;
      }

    /* primes^1 */
    struct prime_iterator_s *it = &info->iterator[0];
    uint64_t max_p = B1_next;
    uint64_t p = it->current_prime;
    while (p <= max_p)
      {
          mulcascade_mul (cascade, p);
          p = getprime_mt (it->p_i);
      }
    it->current_prime = p;

    info->B1 = B1_next;
    mulcascade_get_z_with_clear (s, cascade);
    mulcascade_free (cascade);
}

