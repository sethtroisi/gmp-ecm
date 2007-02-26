#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include "ecm-impl.h"

int 
pm1fs2(mpz_t f, mpres_t X, mpmod_t modulus, unsigned long P, 
       unsigned long blocks)
{
    const unsigned long d = eulerphi (2*P);
    unsigned long i, j, l;
    listz_t F, B, C, tmp, R;
    mpres_t Xi, X2, XP, invXP, findiff[3];
    mpz_t mt;

/*    ASSERT (P & 1 == 1); */

    outputf (OUTPUT_VERBOSE, "P = %lu, d = %lu\n", P, d);

    mpz_init (mt);    /* All-purpose temp mpz_t */
    F = init_list (d);
    tmp = init_list (2 * d + list_mul_mem (d / 2));
    mpres_init (X2, modulus);
    mpres_mul (X2, X, X, modulus); /* X2 = X^2 */
    mpres_init (Xi, modulus);
    mpres_set (Xi, X, modulus);    /* Xi = X^i */
    /* Prepare polynomial F(x), which is monic of degree d. The leading
       monomial is not stored. */
    /* Put in F[0 .. d-1] the values of X^i, 1<=i<2P, gcd(i, 2P) == 1 */
    for (i = 1, j = 0; i < 2*P; i += 2)
    {
	if (gcd (i, P) == 1)
	{
	    mpres_get_z (F[j], Xi, modulus);
	    outputf (OUTPUT_RESVERBOSE, "f_%lu = %Zd\n", j, F[j]);
	    j++;
	}
	mpres_mul (Xi, Xi, X2, modulus);
    }

    ASSERT(j = d);
    mpres_clear (X2, modulus);

    /* Multiply all the (x - f_i) to form F(x) in monomial basis */
    PolyFromRoots (F, F, d, tmp, modulus->orig_modulus);

    for (j = 0; j < d; j++)
	outputf (OUTPUT_RESVERBOSE, "F[%lu] = %Zd\n", j, F[j]);
 
    mpres_init (XP, modulus);
    mpz_set_ui (mt, P);
    mpres_pow (XP, X, mt, modulus); /* XP = X^P */
    
    B = init_list (2 * d);
    C = init_list (d + 1);
    R = init_list (d);
    
    /* Prepare the polynomial B of degree 2*d-1, but not necessarily monic. 
       Since this is invariant over the different blocks, we need to 
       compute it only once */

    /* We want b_j = X^(-2P*(j-d)^2/2) = XP^(-(j-d)^2), 
       for j = 0 ... 2*d. This sequence is symmetric around j = d,
       so only compute one half, and copy the other. This b_j sequence
       is the same for all blocks, so we need to compute it only once.
       We can compute it with the finite differences of the exponents
         -(j+1-d)^2 - -(j-d)^2 = 2(d-j)-1,
         2(d-(j+1))-1 - 2(d-j)-1 = -2j
    */
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

    /* Now mirror-copy the low half into the high half */
    for (j = 1; j < d; j++)
	mpz_set (B[d + j], B[d - j]);

    for (j = 0; j < 2*d; j++)
	outputf (OUTPUT_RESVERBOSE, "B[%lu - %lu] = %Zd;\n", 2*d, j, B[j]);


    for (l = 0; l < blocks; l++)
    {
	/* Now the multipoint evaluation. We want to evaluate F(x) on
	   X^(2*P*l*d + 2*P*i), for successive i, or rewrite as
	   XP^(2*l*d + 2*i) with XP=X^P */
	/* Prepare polynomial C. We want 
	   c_j = f_j * XP^(\beta*j^2/2 + \alpha*j), j = 0 ... d
	   with \beta = 2 and \alpha = 2*l*d, so we can rewrite
	   c_j = f_j * XP^(j^2 + 2*l*d*j), j = 0 ... d
	   
	   We have 2*j + (2*d*l + 1) and 2 for the finite differences of 
	   the exponents of XP.
	*/
	
	mpres_mul (findiff[0], XP, XP, modulus); /* fd[0] = XP^2 */
	mpz_set_ui (mt, 2*d*l + 1);
	mpres_pow (findiff[1], XP, mt, modulus); /* fd[1] = XP^(2*d*l + 1) */
	mpres_set_ui (findiff[2], 1, modulus); /* j=0, XP^(j^2+2 l d j) = 1 */
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
	for (j = 0; j <= d; j++)
	    outputf (OUTPUT_RESVERBOSE, "C[%lu - %lu] = %Zd;\n", d + 1, j, C[j]);
	
	/* Do the convolution */
	/* TMulGen reverses the first input sequence, which we don't want.
	   We can fill C[] in reverse order, for now reverse it separately 
	   here. */
	printf ("Swapping C\n");
	for (j = 0; 2*j+1 < d; j++)
	    mpres_swap (C[j], C[d - 1 - j], modulus);
	printf ("TMulGen\n");
	TMulGen (R, d - 1, C, d, B, 2 * d - 1, tmp, modulus->orig_modulus);
	for (j = 0; j < d; j++)
	    mpz_mod (R[j], R[j], modulus->orig_modulus);

	for (j = 0; j < d; j++)
	    outputf (OUTPUT_RESVERBOSE, "R[%lu] = %Zd\n", j, R[j]);
    }

    return 0;
}
