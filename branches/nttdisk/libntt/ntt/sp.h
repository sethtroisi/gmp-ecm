#ifndef _SP_H
#define _SP_H

#include <math.h>
#include "libntt.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Basic defs for the data types used in Number Theoretic Transforms

   The choice of basic data type has to balance computational
   efficiency with the ability to find enough small prime moduli
   to be useful. We use the GMP concept of 'nails' when choosing
   the number of bits in a prime modulus, since a modular reduction
   using B-bit arithmetic runs much faster when the modulus has
   B-2 bits or less.
   
   64-bit machines with a 64-bit OS should use 62-bit primes and 
   64-bit operations, except that Microsoft's compiler defaults 
   to 32 bits for an unsigned long even on 64-bit versions of 
   Windows, and the GMP limb size is 64 bits just to make things 
   more confusing. 
   
   For 32-bit machines the most efficient arithmetic derives from 
   using 30-bit moduli, but for large transform sizes we will run 
   out of moduli too soon. Using 31-bit moduli delays the inevitable 
   for long enough to be useful for larger problems, at a cost of 
   ~10% in slower modular reductions, but even then we will still
   run out of primes for transform lengths around 2^20, i.e. not 
   very much. So to deal with the largest problems we will need to 
   support 62-bit moduli on a 32-bit machine.

   A completely different possibility for some machines is to
   use 50-bit moduli and use the FPU for modular reductions. If 
   the instruction includes fused multiply-add (FMA) instructions 
   that perform a double-precision multiply and add with a single
   rounding operation, then the result is required to have 106 bits
   of precision before its bottom half is discarded. The upshot is 
   that we can compute a 101-bit signed product of two 51-bit signed 
   inputs in just two FPU instructions. Even without an FMA, the 
   SIMD transforms may be faster than using the ALU if the SIMD vector 
   length is long enough to amortize the cost of the extra arithmetic.

   The defines below choose 30, 31,50 or 62-63 bit prime moduli, a choice 
   that is decoupled from a 32-bit or 64-bit machine word size. */

/*********
 * TYPES *
 *********/

/* SP */

/* the type for both a small prime modulus, and a residue
 * modulo a small prime. The default type is the largest
 * supported by GMP arithmetic; redefine SP_NUMB_BITS
 * to override
 *
 * For a residue x modulo a sp p, we require 0 <= x < p */

#ifndef SP_NUMB_BITS
  #if GMP_LIMB_BITS == 64
    #define SP_NUMB_BITS 62
  #elif GMP_LIMB_BITS == 32
    #define SP_NUMB_BITS 30
  #endif
#endif

#if SP_NUMB_BITS < 30 || SP_NUMB_BITS > 63
#error "invalid choice of small prime size"
#endif

#if SP_NUMB_BITS < 32

  typedef uint32_t sp_t;
  #define SP_TYPE_BITS 32
  #define PRIdsp PRId32
  #define PRIusp PRIu32
  #define PRIxsp "08" PRIx32

  #define SP_MIN ((sp_t)1 << (SP_NUMB_BITS - 1))
  #define SP_MAX ((sp_t)(-1) >> (SP_TYPE_BITS - SP_NUMB_BITS))

#elif SP_NUMB_BITS == 50

  typedef double sp_t;
  #define SP_TYPE_BITS 64
  #define PRIdsp ".0lf"
  #define PRIusp ".0lf"
  #define PRIxsp ".0lf"

  #define SP_MIN pow(2.0, SP_NUMB_BITS - 1)
  #define SP_MAX (pow(2.0, SP_NUMB_BITS) - 1)

#else

  typedef uint64_t sp_t;
  #define SP_TYPE_BITS 64
  #define PRIdsp PRId64
  #define PRIusp PRIu64
  #define PRIxsp "016" PRIx64

  #define SP_MIN ((sp_t)1 << (SP_NUMB_BITS - 1))
  #define SP_MAX ((sp_t)(-1) >> (SP_TYPE_BITS - SP_NUMB_BITS))

#endif

#if SP_NUMB_BITS == 50

  /* The floating point routines *require* IEEE 53-bit 
     double precision, even on x86 processors that 
     support higher precision */

  #if defined(_MSC_VER)
    #include <float.h>
    typedef uint32_t fpu_precision_t;

    static INLINE fpu_precision_t fpu_set_precision_ieee(void) 
      {
	fpu_precision_t old_prec = _control87(0, 0);
	_control87(_PC_53, _MCW_PC);
	return old_prec;
      }

    static INLINE void fpu_clear_precision(fpu_precision_t old_prec) 
      {
	_control87(old_prec, 0xffffffff);
      }

    static INLINE uint32_t fpu_precision_is_ieee(void) 
      {
	fpu_precision_t prec = _control87(0, 0);
	return ((prec & _MCW_PC) == _PC_53) ? 1 : 0;
      } 

  #elif (defined(__GNUC__) || defined(__ICL)) && \
	(defined(__i386__) || defined(__x86_64__))

    #include <float.h>
    typedef uint16_t fpu_precision_t;
    
    static INLINE fpu_precision_t fpu_set_precision_ieee(void) 
      {
	fpu_precision_t old_prec, new_prec;
	asm volatile ("fnstcw %0":"=m"(old_prec)::"memory");
	new_prec = (old_prec & ~0x0300) | 0x0200;
	asm volatile ("fldcw %0": :"m"(new_prec));
	return old_prec;
      }

    static INLINE void fpu_clear_precision(fpu_precision_t old_prec) 
      {
	asm volatile ("fldcw %0": :"m"(old_prec));
      }

    static INLINE uint32_t fpu_precision_is_ieee(void) 
      {
	fpu_precision_t prec;
	asm volatile ("fnstcw %0":"=m"(prec)::"memory");
	return ((prec & ~0x0300) == 0x0200) ? 1 : 0;
      } 

  #else
    typedef uint32_t fpu_precision_t;
    
    static INLINE fpu_precision_t fpu_set_precision_ieee(void) 
      {
	return 0;
      }

    static INLINE void fpu_clear_precision(fpu_precision_t old_prec) 
      {
      }

    static INLINE uint32_t fpu_precision_is_ieee(void) 
      {
	return 1;
      } 
  #endif

#else

  /* for integer types, add GMP plumbing */

  /* types that longlong.h needs */
  typedef mp_limb_t UWtype;
  typedef unsigned int UHWtype;
  #if (defined(_PA_RISC1_1) && defined(__GNUC__))
    /* this seems to be needed, otherwise umul_ppmm() does not work properly */
    typedef mp_limb_t USItype __attribute__ ((mode (SI)));
    typedef mp_limb_t UDItype __attribute__ ((mode (DI)));
  #else
    typedef mp_limb_t USItype;
    typedef mp_limb_t UDItype;
  #endif

    #ifndef W_TYPE_SIZE
    #define W_TYPE_SIZE GMP_LIMB_BITS
    #endif

    #define LONGLONG_STANDALONE
    #include "longlong.h"
#endif

  /* name mangling for a fat binary */
#define MANGLE_NAME(name, psize, wsize) MANGLE_NAME_(name, psize, wsize)
#define MANGLE_NAME_(name, psize, wsize) name##_##sp##psize##w##wsize

#define MANGLE_NAME_STR(psize, wsize) MANGLE_NAME_STR_(psize, wsize)
#define MANGLE_NAME_STR_(psize, wsize) "sp" #psize "w" #wsize

#define MANGLE_SSE2(name) MANGLE_SSE2_(name)
#define MANGLE_SSE2_(name) name##sse2

#define MANGLE_SSE42(name) MANGLE_SSE42_(name)
#define MANGLE_SSE42_(name) name##sse42

#define MANGLE_AVX(name) MANGLE_AVX_(name)
#define MANGLE_AVX_(name) name##avx

#define MANGLE_AVX2(name) MANGLE_AVX2_(name)
#define MANGLE_AVX2_(name) name##avx2

#define MANGLE_FMA(name) MANGLE_FMA_(name)
#define MANGLE_FMA_(name) name##fma

#define X(name) MANGLE_NAME(name, SP_NUMB_BITS, GMP_LIMB_BITS)
#define SP_NAME_SUFFIX_STR MANGLE_NAME_STR(SP_NUMB_BITS, GMP_LIMB_BITS)

/* vector of residues modulo a common small prime */
typedef sp_t * spv_t;
typedef const sp_t * spv_tc;

/* the type used for array offsets */
typedef size_t spv_size_t;

/*************
 * FUNCTIONS *
 *************/

/* general */

/* Conversion functions sp_t <-> mpz_t. Using mpz_*_ui() functions is not
   portable as those take unsigned long's, but on some systems 
   (e.g. 64 bit Windows with Visual C), unsigned long has 32 bits */

static INLINE void 
mpz_set_sp (mpz_t m, const sp_t n)
{
#if SP_TYPE_BITS == 32
  mpz_set_ui(m, (unsigned long)(uint32_t)n);
#elif SP_NUMB_BITS == 50
  mpz_set_d(m, n);
#else /* 64-bit sp_t */
  mpz_set_uint64(m, n);
#endif
}

static INLINE sp_t 
mpz_get_sp (const mpz_t n)
{
#if SP_TYPE_BITS == 32

  return mpz_get_ui(n);

#elif SP_NUMB_BITS == 50

  mpz_get_d(n);

#else /* 64-bit sp_t */

  sp_t m = mpz_getlimbn(n, 0);

  #if GMP_LIMB_BITS == 32  /* 32-bit GMP limb, 64-bit sp_t */
  if (mpz_size(n) >= 2)
    m |= (sp_t)mpz_getlimbn(n, 1) << 32;
  #endif

  return m;
#endif
}

sp_t X(sp_reciprocal)(sp_t p);

/* sp */

/* Routines for arithmetic on residues modulo a small prime
 *
 * All functions return values in the range 0 <= x < p.
 *
 * The variable name of the modulus is 'p' if the input must be prime,
 *                                     'm' if we also allow composites. */


static INLINE sp_t sp_sub(sp_t a, sp_t b, sp_t m) 
{
#if (defined(__GNUC__) || defined(__ICL)) && \
    (defined(__x86_64__) || defined(__i386__)) && \
    SP_NUMB_BITS != 50

  #if SP_TYPE_BITS <= GMP_LIMB_BITS
  sp_t t = 0, tr = a;

  __asm__ (
    "sub %2, %0   # sp_sub: tr -= b\n\t"
    "cmovc %3, %1 # sp_sub: if (a < b) t = m\n\t"
    : "+&r" (tr), "+r" (t)
    : "g" (b), "g" (m)
    : "cc"
  );

  return tr + t;

  #else  /* 64-bit sp_t on 32-bit machine */

  mp_limb_t tlo = 0;
  mp_limb_t thi = 0;
  mp_limb_t alo = (uint32_t)a;
  mp_limb_t ahi = (uint32_t)(a >> 32);
  mp_limb_t blo = (uint32_t)b;
  mp_limb_t bhi = (uint32_t)(b >> 32);
  mp_limb_t mlo = (uint32_t)m;
  mp_limb_t mhi = (uint32_t)(m >> 32);

  __asm__ (
    "sub %4, %2     \n\t"
    "sbb %5, %3     \n\t"
    "cmovc %6, %0   \n\t"
    "cmovc %7, %1   \n\t"
    "add %0, %2     \n\t"
    "adc %1, %3     \n\t"
    : "+r" (tlo), "+r" (thi), "+r" (alo), "+r" (ahi)
    : "g" (blo), "g" (bhi), "g" (mlo), "g" (mhi) 
    : "cc"
  );

  return (sp_t)ahi << 32 | alo;

  #endif

#elif defined(_MSC_VER) && !defined(_WIN64) && (SP_TYPE_BITS == 32)

  __asm
    {
        mov     eax, a
        xor     edx, edx
        sub     eax, b
        cmovb   edx, m
        add     eax, edx
    }
#else
  if (a >= b)
    return a - b;
  else
    return a - b + m;
#endif
}

static INLINE sp_t sp_add(sp_t a, sp_t b, sp_t m) 
{
#if (defined(__GNUC__) || defined(__ICL)) && \
    (defined(__x86_64__) || defined(__i386__)) && \
    SP_NUMB_BITS != 50

  #if SP_TYPE_BITS <= GMP_LIMB_BITS
  sp_t t = a - m, tr = a + b;

  __asm__ (
    "add %2, %1    # sp_add: t += b\n\t"
    "cmovc %1, %0  # sp_add: if (cy) tr = t \n\t"
    : "+r" (tr), "+&r" (t)
    : "g" (b)
    : "cc"
  );

  return tr;

  #else  /* 64-bit sp_t on 32-bit machine */

  mp_limb_t tlo = 0;
  mp_limb_t thi = 0;
  mp_limb_t alo = (uint32_t)a;
  mp_limb_t ahi = (uint32_t)(a >> 32);
  mp_limb_t blo = (uint32_t)b;
  mp_limb_t bhi = (uint32_t)(b >> 32);
  mp_limb_t mlo = (uint32_t)m;
  mp_limb_t mhi = (uint32_t)(m >> 32);

  __asm__ (
    "add %4, %2     \n\t"
    "adc %5, %3     \n\t"
    "sub %6, %2     \n\t"
    "sbb %7, %3     \n\t"
    "cmovc %6, %0   \n\t"
    "cmovc %7, %1   \n\t"
    "add %0, %2     \n\t"
    "adc %1, %3     \n\t"
    : "+r" (tlo), "+r" (thi), "+r" (alo), "+r" (ahi)
    : "g" (blo), "g" (bhi), "g" (mlo), "g" (mhi) 
    : "cc"
  );

  return (sp_t)ahi << 32 | alo;

  #endif

#elif defined(_MSC_VER) && !defined(_WIN64) && (SP_TYPE_BITS == 32)
  __asm
    {
        mov     eax, a
        add     eax, b
        mov     edx, eax
        sub     edx, m
        cmovnc  eax, edx
    }
#else
  sp_t t = a + b;
  if (t >= m)
    t -= m;
  return t;
#endif
}

/* a copy of sp_add, but operating on array offsets. These
   are always size_t types, and may have a size different 
   from an sp_t */

static INLINE spv_size_t sp_array_inc(spv_size_t a, spv_size_t b, spv_size_t m) 
{
#if (defined(__GNUC__) || defined(__ICL)) && \
    (defined(__x86_64__) || defined(__i386__))

  spv_size_t t = a - m, tr = a + b;

  __asm__ (
    "add %2, %1    # sp_array_inc: t += b\n\t"
    "cmovc %1, %0  # sp_array_inc: if (cy) tr = t \n\t"
    : "+r" (tr), "+&r" (t)
    : "g" (b)
    : "cc"
  );

  return tr;

#elif defined(_MSC_VER) && !defined(_WIN64)

  __asm
    {
        mov     eax, a
        add     eax, b
        mov     edx, eax
        sub     edx, m
        cmovnc  eax, edx
    }

#else

  spv_size_t t = a + b;
  if (t >= m)
    t -= m;
  return t;

#endif
}

/* widening multiply */

#if SP_NUMB_BITS == 50 /* floating point sp_t */

  /* cribbed from DH Bailey's QD package, which references

     J. Shewchuk, "Adaptive Precision Floating-Point Arithmetic
     and Fast Robust Geometric Predicates" */

  /* find hi and lo such that a exactly equals hi + lo
     and the bottom half of each mantissa is 0 */

  #define sp_split(hi, lo, a)            \
	{				 \
	  sp_t t = 134217729.0 * a;	 \
	  hi = t - (t - a);		 \
	  lo = a - hi;			 \
	}

  /* find hi and lo such that a * b exactly equals hi+lo */

  #define sp_wide_mul(hi, lo, a, b)	 \
	{				 \
	  sp_t ahi, alo, bhi, blo;	 \
	  sp_split(ahi, alo, a);	 \
	  sp_split(bhi, blo, b);	 \
	  hi = a * b;			 \
	  lo = ((ahi * bhi - hi) + ahi * blo + alo * bhi) + alo * blo; \
	}

  #define sp_wide_sqr(hi, lo, a)	\
	{				\
	  sp_t ahi, alo;		\
	  sp_split(ahi, alo, a);	\
	  hi = a * a;			\
	  lo = ((ahi * ahi - hi) + 2.0 * alo * ahi) + alo * alo; \
	}

#elif SP_TYPE_BITS == GMP_LIMB_BITS  /* use GMP functions */

  #define sp_wide_mul(hi, lo, a, b) umul_ppmm(hi, lo, a, b)
  #define sp_wide_sqr(hi, lo, a) sp_wide_mul(hi, lo, a, a)

#elif SP_TYPE_BITS < GMP_LIMB_BITS  /* ordinary multiply */

  #define sp_wide_mul(hi, lo, a, b)		    \
  	{					    \
	  mp_limb_t prod = (mp_limb_t)(a) *	    \
	  		   (mp_limb_t)(b);	    \
	  hi = (sp_t)(prod >> (GMP_LIMB_BITS / 2)); \
	  lo = (sp_t)prod;			    \
	}

  #define sp_wide_sqr(hi, lo, a) sp_wide_mul(hi, lo, a, a)

#else  /* worst case: build product from smaller multiplies */

  #define sp_wide_mul(hi, lo, a, b)		    		\
  	{					    		\
	  mp_limb_t a1 = (mp_limb_t)((a) >> GMP_LIMB_BITS);	\
	  mp_limb_t a0 = (mp_limb_t)(a);			\
	  mp_limb_t b1 = (mp_limb_t)((b) >> GMP_LIMB_BITS);	\
	  mp_limb_t b0 = (mp_limb_t)(b);			\
	  mp_limb_t a0b0_hi, a0b0_lo;				\
	  mp_limb_t a0b1_hi, a0b1_lo;				\
	  mp_limb_t a1b0_hi, a1b0_lo;				\
	  mp_limb_t a1b1_hi, a1b1_lo;				\
								\
	  umul_ppmm(a0b0_hi, a0b0_lo, a0, b0);			\
	  umul_ppmm(a0b1_hi, a0b1_lo, a0, b1);			\
	  umul_ppmm(a1b0_hi, a1b0_lo, a1, b0);			\
	  umul_ppmm(a1b1_hi, a1b1_lo, a1, b1);			\
								\
	  /* proceed a word at a time; the high word		\
	     of each product except a0b0 have at least		\
	     two bits clear, so carries propagate one		\
	     word at most */					\
								\
	  add_ssaaaa(a0b1_hi, a0b1_lo, 				\
	             a0b1_hi, a0b1_lo,				\
		     (mp_limb_t)0, a0b0_hi);			\
	  add_ssaaaa(a0b1_hi, a0b1_lo, 				\
	             a0b1_hi, a0b1_lo,				\
		     a1b0_hi, a1b0_lo);				\
	  add_ssaaaa(a1b1_hi, a1b1_lo, 				\
	             a1b1_hi, a1b1_lo,				\
		     (mp_limb_t)0, a0b1_hi);			\
								\
	  hi = (sp_t)a1b1_hi << GMP_LIMB_BITS | a1b1_lo;	\
	  lo = (sp_t)a0b1_lo << GMP_LIMB_BITS | a0b0_lo;	\
	}

  #define sp_wide_sqr(hi, lo, a)		    		\
  	{					    		\
	  mp_limb_t a1 = (mp_limb_t)((a) >> GMP_LIMB_BITS);	\
	  mp_limb_t a0 = (mp_limb_t)(a);			\
	  mp_limb_t a0a0_hi, a0a0_lo;				\
	  mp_limb_t a0a1_hi, a0a1_lo;				\
	  mp_limb_t a1a1_hi, a1a1_lo;				\
								\
	  umul_ppmm(a0a0_hi, a0a0_lo, a0, a0);			\
	  umul_ppmm(a0a1_hi, a0a1_lo, a0, a1);			\
	  umul_ppmm(a1a1_hi, a1a1_lo, a1, a1);			\
								\
	  add_ssaaaa(a0a1_hi, a0a1_lo, 				\
	             a0a1_hi, a0a1_lo,				\
		     a0a1_hi, a0a1_lo);				\
	  add_ssaaaa(a0a1_hi, a0a1_lo, 				\
	             a0a1_hi, a0a1_lo,				\
		     (mp_limb_t)0, a0a0_hi);			\
	  add_ssaaaa(a1a1_hi, a1a1_lo, 				\
	             a1a1_hi, a1a1_lo,				\
		     (mp_limb_t)0, a0a1_hi);			\
								\
	  hi = (sp_t)a1a1_hi << GMP_LIMB_BITS | a1a1_lo;	\
	  lo = (sp_t)a0a1_lo << GMP_LIMB_BITS | a0a0_lo;	\
	}

#endif

/* modular reduction */

#if SP_NUMB_BITS == 50

static INLINE sp_t sp_udiv_rem(sp_t nhi, sp_t nlo, sp_t d, sp_t di)
{
  sp_t r, r2, q, qhi, qlo;
  ATTRIBUTE_UNUSED sp_t tmp;

  q = floor(nhi * di);
  sp_wide_mul (qhi, qlo, q, d);

  /* form (nhi:nlo) - (qhi:qlo) exactly; was assume the result
     fits into an sp_t */

  r = nhi - qhi;
  r2 = r - nhi;
  r2 = (nhi - (r - r2)) + (qhi + r2);
  r2 = (r2 + nlo) - qlo;
  r = r + r2;

  r = sp_sub(r, d, d);
  if (r < 0)
    r += d;

  return r;
}

#elif SP_NUMB_BITS <= SP_TYPE_BITS - 2

    /* having a small modulus allows the reciprocal
     * to be one bit larger, which guarantees that the
     * initial remainder fits in a word and also that at
     * most one correction is necessary */

static INLINE sp_t sp_udiv_rem(sp_t nh, sp_t nl, sp_t d, sp_t di)
{
  sp_t r, q1, q2;
  ATTRIBUTE_UNUSED sp_t tmp;

  q1 = nh << (2*(SP_TYPE_BITS - SP_NUMB_BITS)) |
	    nl >> (2*SP_NUMB_BITS - SP_TYPE_BITS);
  sp_wide_mul (q2, tmp, q1, di);
  r = nl - d * (q2 >> 1);
  return sp_sub(r, d, d);
}

#else                       /* big integer modulus; no shortcuts allowed */

static INLINE sp_t sp_udiv_rem(sp_t nh, sp_t nl, sp_t d, sp_t di)
{
  sp_t q1, q2, dqh, dql;
  ATTRIBUTE_UNUSED sp_t tmp;
  q1 = nh << (2*(SP_TYPE_BITS - SP_NUMB_BITS)) |
	    nl >> (2*SP_NUMB_BITS - SP_TYPE_BITS);
  sp_wide_mul (q2, tmp, q1, di);
  sp_wide_mul (dqh, dql, q2, d);

  #if SP_TYPE_BITS < GMP_LIMB_BITS 
  {
    mp_limb_t t = ((mp_limb_t)nh << 32 | (mp_limb_t)nl) -
                  ((mp_limb_t)dqh << 32 | (mp_limb_t)dql);
    nh = t >> SP_TYPE_BITS;
    nl = (sp_t)t;
  }
  #else
  sub_ddmmss (nh, nl, nh, nl, dqh, dql);
  #endif

  if (nh)
    nl -= d;
  nl = sp_sub(nl, d, d);
  return sp_sub(nl, d, d);
}

#endif

/* x*y mod m */
static INLINE sp_t
sp_mul (sp_t x, sp_t y, sp_t m, sp_t d)
{
  sp_t u, v;
  sp_wide_mul (u, v, x, y);
  return sp_udiv_rem (u, v, m, d);
}

/* x*y mod m */
static INLINE sp_t
sp_sqr (sp_t x, sp_t m, sp_t d)
{
  sp_t u, v;
  sp_wide_sqr (u, v, x);
  return sp_udiv_rem (u, v, m, d);
}

#define sp_neg(x,m) sp_sub(0, x, m)

/* Returns x^a % m, uses a right-to-left powering ladder */

static INLINE sp_t
sp_pow (sp_t x, sp_t e, sp_t m, sp_t d)
{
  sp_t partial = 1;
  uint64_t a = (uint64_t)e;

  while (1)
    {
      if (a & 1)
	partial = sp_mul (x, partial, m, d);

      a >>= 1;

      if (!a)
	return partial;

      x = sp_sqr (x, m, d);
    }
}

/* 1/x mod p */
#define sp_inv(x,p,d) sp_pow (x, (p) - 2, p, d)

int X(sp_prime)(sp_t);

/* spv */

void X(spv_random)(spv_t, spv_size_t, sp_t);

#ifdef __cplusplus
}
#endif

#endif /* _SP_H */
