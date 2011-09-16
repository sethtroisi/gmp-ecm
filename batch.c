/* ECM stage 1 in batch mode, for initial point (x:z) with small coordinates,
   such that x and z fit into an unsigned long.
   For example we can start with (x=2:y=1) with the curve by^2 = x^3 + ax^2 + x
   with a = 4d-2 and b=16d+2, then we have to multiply by d=(a+2)/4 in the
   duplicates.
*/

#include "ecm-gmp.h"
#include "ecm-impl.h"

#define USE_REDC

#ifndef __gmpn_redc_2
void __gmpn_redc_2 (mp_ptr, mp_ptr, mp_srcptr, mp_size_t, mp_srcptr);
#endif

//mpz_t u, v, w;
//mp_limb_t invm[2];

/* a <- b*c mod n */
static void
mpz_mulmod_ui (mpz_t a, mpz_t b, unsigned long c, mpz_t n)
{
  mpz_mul_ui (a, b, c);
  mpz_mod (a, a, n);
}

static unsigned long MUL = 0, SQR = 0;


/* (x3:z3) <- (x2:z2) + (x1:z1) assuming (x:z) = (2:1) = (x2:z2) - (x1:z1) */
/*
static void
add (mpres_t x3, mpres_t z3, mpres_t x2, mpres_t z2, mpres_t x1, mpres_t z1, 
		mpres_t u, mpres_t v, mpres_t w, mpmod_t n)
{
  mpres_sub (u, x2, z2, n);
  mpres_add (v, x1, z1, n);
  mpres_mul (u, u, v, n);
  mpres_add (w, x2, z2, n);
  mpres_sub (v, x1, z1, n);
  mpres_mul (v, w, v, n);
  mpres_add (w, u, v, n);
  mpres_sub (v, u, v, n);
  mpres_mul (v, v, v, n);

  mpres_mul (x3, w, w, n);
  //mpres_mul_ui (z3, v, 2, n);
  mpres_add (z3, v, v, n);
  //mpz_mul_2exp (z3, v, 1);
  //mpz_mod (z3, z3, n);
  MUL += 2;
  SQR += 2;
}
*/

/* (x2:z2) <- 2 * (x1:z1) */
/*
static void
dup (mpres_t x2, mpres_t z2, mpres_t x1, mpres_t z1, unsigned long d,
		mpres_t u, mpres_t v, mpres_t w,
     mpmod_t n)
{
  mpres_add (u, x1, z1, n);
  mpres_mul (u, u, u, n);
  mpres_sub (v, x1, z1, n);
  mpres_mul (v, v, v, n);
  mpres_mul (x2, u, v, n);
  mpres_sub (w, u, v, n);
  mpres_mul_ui (u, w, d, n);
  mpres_add (u, u, v, n);
  mpres_mul (z2, w, u, n);
  MUL += 2;
  SQR += 2;
}
*/

/* (x1:z1) <- 2(x1:z1) 
	 (x2:z2) <- (x1:z1) + (x2:z2) 
	 assume (x2:z2) - (x1:z1) = (2:1) */
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
	mpres_mul (w, w, w, n); /* w = (x1+z1)^2 */
	mpres_mul (u, u, u, n); /* u = (x1-z1)^2 */
  
	mpres_mul (x1, u, w, n); /* xdup = (x1+z1)^2 * (x1-z1)^2 */

  mpres_sub (w, w, u, n);   /* w = (x1+z1)^2 - (x1-z1)^2 */
	
	mpres_mul_ui (q, w, d, n); /* q = d* ((x1+z1)^2 - (x1-z1)^2) */
	/* use mpz_mulmod_ui?? => segfault ??*/
	//mpz_mulmod_ui(q, w, d, n);

	mpres_add (u, u, q, n);  /* u = (x1-z1)^2 - d* ((x1+z1)^2 - (x1-z1)^2) */
	mpres_mul (z1, w, u, n); /* zdup = w * [(x1-z1)^2 - d* ((x1+z1)^2 - (x1-z1)^2)] */

  mpres_add (w, v, t, n);
  mpres_sub (v, v, t, n);
  
	mpres_mul (v, v, v, n);
	mpres_mul (x2, w, w, n);
	//mpres_mul_ui (z2, v, 2, n);
	mpres_add (z2, v, v, n);
  
	MUL += 4;
  SQR += 4;
}
/*
static void
ecm_stage1_init (mpmod_t n)
{
  mpz_t u, v;

  mpz_init (u);
  mpz_init (v);

  mpz_set_ui (u, 1);
  mpz_mul_2exp (u, u, 2 * GMP_NUMB_BITS);
  mpz_neg (v, n->orig_modulus);
  mpz_invert (v, v, u);
  //invm[0] = mpz_getlimbn (v, 0);
  //invm[1] = mpz_getlimbn (v, 1);

  mpz_clear (u);
  mpz_clear (v);
}
*/
#define MAX_HEIGHT 30

void compute_s (mpz_t s, unsigned int B1)
{
	mpz_t l[MAX_HEIGHT];
	mpz_t r[MAX_HEIGHT];
	unsigned int i, j;
	unsigned long pi=2, pp;
	for (i=0;i<MAX_HEIGHT;i++)
	{
		mpz_init(l[i]);
		mpz_init(r[i]);
	}

	i=0;
	while (pi<=B1)
	{
		pp=pi;
		while(pp<=B1)
			pp*=pi;
		pp/=pi;

		if (i%2==0)
			mpz_set_ui(l[0],pp);
		else
			mpz_set_ui(r[0],pp);
			
		j=0;
		while ((i&(1<<j))!=0)
		{
			if ((i&(1<<(j+1)))==0)
				mpz_mul(l[j+1],l[j],r[j]);
			else
				mpz_mul(r[j+1],l[j],r[j]);
			j++;
		}

		i++;
		pi=getprime(pi);
	}

	if (i%2==0)
		mpz_set_ui(r[0],1);
	else
		mpz_set_ui(r[0],1);
		
	j=0;
	for (j=0;j<MAX_HEIGHT-1;j++)
	{
		if ((i&(1<<(j)))==0)
			mpz_set(r[j+1],r[j]);
		else
			mpz_mul(r[j+1],l[j],r[j]);
	}
	
	//fprintf(stdout,"%u primes under B1\n",i);
	
	mpz_set(s,r[MAX_HEIGHT-1]);
	
  getprime_clear (); /* free the prime tables, and reinitialize */
	for (i=0;i<MAX_HEIGHT;i++)
	{
		mpz_clear(l[i]);
		mpz_clear(r[i]);
	}
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
            double *B1done, mpz_t go, int (*stop_asap)(void), 
            char *chkfilename)
{
  mpz_t s; /* product of primes up to B1 */
  unsigned long d; 
  mpz_t x1, z1, x2, z2;
  long i;
  int st;
	mpz_t q, t, u, v, w;

  //double p, r;
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


  mpz_init_set_ui (s, 1);

  printf ("Batch mode : \n");

  /* construct the batch exponent */
	compute_s(s, B1);
  printf ("  Computing s took %ldms\n", cputime () - st);
	printf ("  s has %lu bits\n", mpz_sizeinbase(s,2));

  st = cputime ();
	

  /* initialize P */
  mpres_set (x1, x, n);
  mpres_set_ui (z1, 1, n); /* P1 <- 1P */

	/* Compute d from A */
  mpz_add_ui (u, A, 2);
  mpz_div_2exp (u, u, 2);
	mpz_mod (u, u, n);
  if (mpz_fits_ulong_p (u) == 0)
    {
      outputf (OUTPUT_ERROR, "Error, d=(A+2)/4 should fit in an ulong, A=%Zd\n",
		           A);
	  }
	d=mpz_get_ui(u); 

  /* Compute 2P : no need to duplicate P, the coordinates are simple. */ 
  mpres_set_ui (x2, 9, n); 
  mpres_set_ui (z2, 8+64*d, n); /* P2 <- 2P */ 
  
  /* invariant: if j represents the upper bits of s,
     then P1 = j*P and P2=(j+1)*P */

  /* now perform the double-and-add ladder */
  for (i = mpz_sizeinbase (s, 2) - 2; i>=0; i--)
  {
    if (mpz_tstbit (s, i) == 0) /* (j,j+1) -> (2j,2j+1) */
		{ /* P2 <- P1+P2    P1 <- 2*P1 */
			dup_add (x1, z1, x2, z2, q, t, u, v , w, d, n);
		}
     else /* (j,j+1) -> (2j+1,2j+2) */
		{/* P1 <- P1+P2     P2 <- 2*P2 */
			dup_add (x2, z2, x1, z1, q, t, u, v, w, d, n);
		}
	}
  
	printf ("  MUL=%lu SQR=%lu\n", MUL, SQR);
  printf ("  Stage 1 took %ldms\n", cputime () - st);

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
	mpz_clear(q);
	mpz_clear(t);
	mpz_clear(u);
	mpz_clear(v);
	mpz_clear(w);

	return ret;
}

#ifdef STANDALONE
int
main ()
{
  mpz_t n, x, z;
  double B1 = 1e6;
  unsigned long d = 1;
  
  mpz_init (n);
  mpz_init (u);
  mpz_init (v);
  mpz_init (w);
  mpz_init (x);
  mpz_init (z);

  mpz_inp_str (n, stdin, 0);
  ecm_stage1_batch (d, x, z, n, B1);

  /* the result is in (x:z) */
  mpz_invert (u, z, n);
  mpz_mul (x, x, u);
  mpz_mod (x, x, n);
  gmp_printf ("x=%Zd\n", x);

  mpz_clear (n);
  mpz_clear (u);
  mpz_clear (v);
  mpz_clear (w);
  mpz_clear (x);
  mpz_clear (z);
}
#endif
