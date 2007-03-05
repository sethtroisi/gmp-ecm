#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include "ecm-impl.h"

int 
pm1fs2(mpz_t f, mpres_t X, mpmod_t modulus, root_params_t *root_params, 
       unsigned long dF, unsigned long blocks)
{
  const unsigned long d = 2 * dF, P = root_params->d1 / 2;
  /* Our beta = 2*P, alpha = 2*P * root_params->i0 */
  mpz_t i0;
  unsigned long i, j, l, tmplen;
  listz_t F, B, C, tmp, R;
  mpres_t Xi, X2, XP, invXP, findiff[3];
  mpz_t mt;
  int youpi = 0;
  extern unsigned int Fermat;
  long timestart, timestop, muls;
  long pariline = 0;

  ASSERT  (root_params->d1 % 2 == 0);
  /*  ASSERT (P & 1 == 1); */

  if (modulus->repr == ECM_MOD_BASE2 && modulus->Fermat > 0)
    Fermat = modulus->Fermat;

  outputf (OUTPUT_VERBOSE, "P = %lu; d = %lu; /* PARI %ld */\n", 
	   P, d, pariline++);
  outputf (OUTPUT_TRACE, "compare(a, b, n) = if (a != b,print(\"In PARI \", n, \" line, \", a \" != \" b));"
	   "/* PARI %ld */\n", pariline++);

  mpz_init (i0);
  mpz_set (i0, root_params->i0);

  /* Allocate all the memory we'll need */
  /* Allocate the correct amount of space for each mpz_t or the 
     reallocations will up to double the time for stage 2! */
  mpz_init (mt);    /* All-purpose temp mpz_t */
  F = init_list2 (d, modulus->bits);
  B = init_list2 (2 * d, modulus->bits);
  C = init_list2 (d + 1, modulus->bits);
  R = init_list2 (3 * d + 1, modulus->bits);    
  tmplen = 2 * d + list_mul_mem (d / 2);
  outputf (OUTPUT_DEVVERBOSE, "tmplen = %lu\n", tmplen);
  if (TMulGen_space (d - 1, d, 2 * d - 1) > tmplen)
    tmplen = TMulGen_space (d - 1, d, 2 * d - 1) + 12;
  /* FIXME: It appears TMulGen_space() returns a too small value! */
  outputf (OUTPUT_DEVVERBOSE, "tmplen = %lu\n", tmplen);
  tmp = init_list (tmplen);
  
  mpres_get_z (mt, X, modulus); /* mpz_t copy of X for printing */
  outputf (OUTPUT_RESVERBOSE, 
	   "N = %Zd; X = Mod(%Zd, N); XP=X^P; /* PARI %ld */\n", 
	   modulus->orig_modulus, mt, pariline++);

  mpz_sub_ui (mt, i0, 1UL);
  mpz_mul_ui (mt, mt, 2UL*P);
  outputf (OUTPUT_VERBOSE, "Effective B2min = %Zd\n", mt);
  mpz_add_ui (mt, i0, blocks * d - 1UL); /* Assumes blocks*d > 0 */
  mpz_mul_ui (mt, mt, 2UL*P);
  outputf (OUTPUT_VERBOSE, "Effective B2max = %Zd\n", mt);

  outputf (OUTPUT_VERBOSE, "Computing roots of F");
  outputf (OUTPUT_TRACE, "\n"); /* So the "Computing" does not interfere */
  timestart = cputime ();
  mpres_init (X2, modulus);
  mpres_init (Xi, modulus);
  mpres_mul (X2, X, X, modulus); /* X2 = X^2 */
  mpres_set (Xi, X, modulus);    /* Xi = X^i */
  /* Prepare polynomial F(x), which is monic of degree d. The leading
     monomial is not stored. */
  /* Put in F[0 .. d-1] the values of X^i, 1<=i<2P, gcd(i, 2P) == 1 */
  ASSERT (eulerphi (2*P) <= d);
  for (i = 1, j = 0; i < 2*P; i += 2)
    {
      if (gcd (i, 2*P) == 1)
	{
	  mpres_get_z (F[j], Xi, modulus);
	  outputf (OUTPUT_TRACE, "f_%lu = X^%lu;"
		   "compare(f_%lu, Mod(%Zd, N), %ld); /* PARI %ld */\n", 
		   j, i, j, F[j], pariline, pariline); pariline++;
	  j++;
	}
      mpres_mul (Xi, Xi, X2, modulus);
    }
  
  ASSERT(j <= d);
  mpres_clear (X2, modulus);
  mpres_clear (Xi, modulus);
  timestop = cputime ();
  outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);

  outputf (OUTPUT_TRACE, "F(x) = ");
  for (j = 0; j < d - 1; j++)
    outputf (OUTPUT_TRACE, "(x - f_%lu) * ", j);
  outputf (OUTPUT_TRACE, "(x - f_%lu); /* PARI %ld */\n", d - 1, pariline++);
  
  /* Multiply all the (x - f_i) to form F(x) in monomial basis */
  outputf (OUTPUT_VERBOSE, "Building F from its roots");
  timestart = cputime ();
  PolyFromRoots (F, F, d, tmp, modulus->orig_modulus);
  timestop = cputime ();
  outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);
  
  for (j = 0; j < d; j++)
    outputf (OUTPUT_TRACE, "F_%lu = Mod(%Zd, N); /* PARI %ld */\n", 
	     j, F[j], pariline++);
  outputf (OUTPUT_TRACE, "F_%lu = 1; /* PARI %ld */\n", d, pariline++);
  
  outputf (OUTPUT_TRACE, "F(x) == x^%lu", d);
  for (j = d; j > 0; j--)
    outputf (OUTPUT_TRACE, " + F_%lu*x^%lu", j - 1, j - 1);
  outputf (OUTPUT_TRACE, " /* PARI %ld */\n", pariline++);
  
  
  mpres_init (XP, modulus);
  mpz_set_ui (mt, P);
  mpres_pow (XP, X, mt, modulus); /* XP = X^P */
  
  /* Prepare the polynomial B of degree 2*d-1, but not necessarily monic. 
     Since this is invariant over the different blocks, we need to 
     compute it only once */
  
  /* We want b_j = X^(-beta*(j-d)^2/2) = X^(-2P*(j-d)^2/2) = XP^(-(j-d)^2), 
     for j = 0 ... 2*d. This sequence is symmetric around j = d,
     so only compute one half, and copy the other. This b_j sequence
     is the same for all blocks, so we need to compute it only once.
     We can compute it with the finite differences of the exponents
     -(j+1-d)^2 - -(j-d)^2 = 2(d-j)-1,
     2(d-(j+1))-1 - 2(d-j)-1 = -2j
  */
  outputf (OUTPUT_VERBOSE, "Computing roots of B");
  timestart = cputime ();
  mpres_init (invXP, modulus);
  mpres_invert (invXP, XP, modulus);
  mpres_init (findiff[0], modulus);
  mpres_init (findiff[1], modulus);
  mpres_init (findiff[2], modulus);
  mpres_mul (findiff[0], invXP, invXP, modulus); /* findiff[0] = XP^(-2) */
  
  mpz_set_ui (mt, 2 * d - 1);
  mpres_pow (findiff[1], XP, mt, modulus); /* findiff[1] = XP^(2*d-1) */
  
  mpz_set_ui (mt, d);
  mpz_mul (mt, mt, mt); /* mt = d^2 (may exceed 2^32, so square in mpz_t) */
  mpres_pow (findiff[2], invXP, mt, modulus); /* findiff[2] = XP^(-d^2) */
  
  for (j = 0; j <= d; j++)
    {
      mpres_get_z (B[j], findiff[2], modulus);
      mpres_mul (findiff[2], findiff[2], findiff[1], modulus);
      mpres_mul (findiff[1], findiff[1], findiff[0], modulus);
    }
  
  /* B[d] = XP^(-(d-d)^2) = 1. Check that it is so */
  ASSERT(mpz_cmp_ui (B[d], 1) == 0);
  timestop = cputime ();
  outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);
  
  outputf (OUTPUT_VERBOSE, "Copying high half of B to low half");
  timestart = cputime ();
  /* Now mirror-copy the low half into the high half */
  for (j = 1; j < d; j++)
    mpz_set (B[d + j], B[d - j]);
  
  timestop = cputime ();
  outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);

  for (j = 0; j < 2*d; j++)
    {
      outputf (OUTPUT_TRACE, 
	       "B_%lu = XP^(-(%lu-d)^2); ", j, j);
      outputf (OUTPUT_TRACE, 
	       "compare(B_%lu, Mod(%Zd, N), %ld); /* PARI %ld */\n", 
	       j, B[j], pariline, pariline); pariline++;
    }
  
  outputf (OUTPUT_TRACE, "B(x) = ");
  for (j = 0; j < 2 * d; j++)
    outputf (OUTPUT_TRACE, " + B_%lu*x^%lu", j, j);
  outputf (OUTPUT_TRACE, "; /* PARI %ld */\n", pariline++);
  
  for (l = 0; l < blocks; l++)
    {
      /* Now the multipoint evaluation. We want to evaluate F(x) on
	 X^(2*P*(i0 + l*d) + 2*P*i), for successive i, or rewrite 
	 as XP^(2*l*d + 2*i0 + 2*i) with XP=X^P */
      /* Prepare polynomial C. We want 
	 c_j = f_j * X^(beta * j^2/2 + alpha * j), j = 0 ... d
	 with beta = 2P and alpha = 2P*(i0 + l*d), so we can rewrite
	 c_j = f_j * XP^(j^2 + 2*(l*d + i0)*j), j = 0 ... d
	 
	 We have 2*j + 2*(d*l+i0) + 1 and 2 for the finite differences of 
	 the exponents of XP.
      */
      outputf (OUTPUT_VERBOSE, "Computing roots of C, block %lu", l);
      timestart = cputime ();
      mpres_mul (findiff[0], XP, XP, modulus); /* fd[0] = XP^2 */
      mpz_set_ui (mt, 2*d*l + 1);
      mpz_add (mt, mt, i0);
      mpz_add (mt, mt, i0);
      mpres_pow (findiff[1], XP, mt, modulus); /* fd[1] = XP^(2*d*l + 1) */
      mpres_set_ui (findiff[2], 1, modulus); /* j=0, XP^(j^2+2 l d i0 j)=1 */
      /* Can we just init this once and let the findiff stuff continue
	 over all the blocks? */
      
      mpz_set (C[0], F[0]); /* fd[2] = 1, so c_0 = f_0 */

      for (j = 1; j < d; j++)
	{
	  mpres_mul (findiff[2], findiff[2], findiff[1], modulus);
	  mpres_mul (findiff[1], findiff[1], findiff[0], modulus);
	  mpres_get_z (mt, findiff[2], modulus);
	  mpz_mul (mt, mt, F[j]);
	  mpz_mod (C[j], mt, modulus->orig_modulus);
	}
      /* F[d] = 1 is not actually stored anywhere. Treat it separately */
      mpres_mul (findiff[2], findiff[2], findiff[1], modulus);
      /* mpres_mul (findiff[1], findiff[1], findiff[0], modulus); Needed 
	 only if we should let the findiff stuff run over several blocks*/
      mpres_get_z (C[j], findiff[2], modulus);
      timestop = cputime ();
      outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);

      for (j = 0; j <= d; j++)
	{
	  outputf (OUTPUT_TRACE, 
		   "C_%lu = F_%lu * XP^(%lu^2 + 2*(%lu*d+%Zd)*%lu); ", 
		   j, j, j, l, i0, j);
	  outputf (OUTPUT_TRACE, "compare(C_%lu, %Zd, %ld); /* PARI %ld */\n", 
		   j, C[j], pariline, pariline); pariline++;
	}
      
      outputf (OUTPUT_TRACE, "C(x) = ");
      for (j = 0; j <= d; j++)
	outputf (OUTPUT_TRACE, " + C_%lu*x^%lu", j, j);
      outputf (OUTPUT_TRACE, "; /* PARI %ld */\n", pariline++);
      
      /* Do the convolution */
#if 1
      /* Use the transposed "Middle Product" algorithm */
      /* TMulGen reverses the first input sequence, which we don't want.
	 We can fill C[] in reverse order, for now reverse it separately 
	 here. */
      outputf (OUTPUT_VERBOSE, "Swapping C\n");
      /* Remember that the C array has length d+1 ! */
      for (j = 0; j < d + 1 - 1 - j; j++)
	mpres_swap (C[j], C[d + 1 - 1 - j], modulus);

      outputf (OUTPUT_TRACE, "; /* PARI %ld */\n", pariline++);
      outputf (OUTPUT_VERBOSE, "TMulGen of B and C");
      ASSERT(tmplen >= TMulGen_space (d - 1, d, 2 * d - 1));
      muls = TMulGen (R, d - 1, C, d, B, 2 * d - 1, tmp, 
		      modulus->orig_modulus);
      list_mod (R, R, d, modulus->orig_modulus);
      timestop = cputime ();
      outputf (OUTPUT_VERBOSE, " took %lu muls, %lu ms\n", 
	       muls, timestop - timestart);
#else
      /* Use two non-transposed multiplications */
      outputf (OUTPUT_VERBOSE, "Computing B * C");
      timestart = cputime ();
      list_mul (R, C, d + 1, 0, B, d, 0, tmp);
      list_set (R, R + d, d);
      list_mul (R + d, C, d, 0, B + d, d, 0, tmp);
      list_add (R, R, R + d, d);
      list_mod (R, R, d, modulus->orig_modulus);
      timestop = cputime ();
      outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);
#endif
#if defined(WANT_ASSERT)
      if (l == 0 && mpz_cmp_ui (i0, 0) == 0)
	{
	  /* i0 == 0 so we evaluated F(X^0) = F(1). We can check that for
	     correctness easily. */
	  mpz_set_ui (mt, 0);
	  for (i = 0; i < d; i++)
	    mpz_add (mt, mt, F[i]);
	  /* Add the leading 1 term which is not stored anywhere */
	  mpz_add_ui (mt, mt, 1);
	  mpz_mod (mt, mt, modulus->orig_modulus);
	  outputf (OUTPUT_DEVVERBOSE, "Testing R[0] = F(1)\n");
	  outputf (OUTPUT_TRACE, "%Zd == %Zd /* PARI */\n", R[0], mt);
	  ASSERT (mpz_cmp (R[0], mt) == 0);
	}
#endif
      for (j = 0; j < d; j++)
	{
	  outputf (OUTPUT_TRACE, "compare(F(XP^(2*(%lu+%lu*d+%Zd)))"
		   "*XP^(-2*%lu^2/2), %Zd, %ld); /* PARI %ld */\n", 
		   j, l, i0, j, R[j], pariline, pariline); 
	  pariline++;
	}
      
      outputf (OUTPUT_VERBOSE, "Computing product of F(g_i)");
      timestart = cputime ();
#if 1
      /* Try a faster way of multiplying up the R[i] for the gcd */
      mpres_set_z_for_gcd (findiff[0], R[0], modulus);
      for (i = 1; i < d; i++)
	{
	  mpres_set_z_for_gcd (findiff[1], R[i], modulus);
	  mpres_mul (findiff[0], findiff[0], findiff[1], modulus); 
	}
      mpres_gcd (tmp[0], findiff[0], modulus);
#else
      list_mulup (R, d, modulus->orig_modulus, tmp[0]); 
      mpz_gcd (tmp[0], R[d - 1], modulus->orig_modulus);
#endif
      timestop = cputime ();
      outputf (OUTPUT_VERBOSE, " took %lu ms\n", timestop - timestart);

      if (mpz_cmp_ui (tmp[0], 1) > 0)
	{
	  outputf (OUTPUT_NORMAL, "Found factor in stage 2\n");
	  mpz_set (f, tmp[0]);
	  youpi = ECM_FACTOR_FOUND_STEP2;
	  break;
	}
      outputf (OUTPUT_RESVERBOSE, "Product of F(g_i) = %Zd\n", R[d - 1]);
    }
  
  clear_list (F, d);
  clear_list (B, 2 * d);
  clear_list (C, d + 1);
  clear_list (R, 3 * d + 1);    
  clear_list (tmp, tmplen);
  mpres_clear (findiff[0], modulus);
  mpres_clear (findiff[1], modulus);
  mpres_clear (findiff[2], modulus);
  
  return youpi;
}
