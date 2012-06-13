#include "auxi.h"
#include <math.h>


/*
  compute k=lcm(2,..,B1)
  WARNING: slow version  (no product tree)
  NOT USED
*/
void calculk (mpz_t k,double B1) {
  double p,r;
  mpz_t q;

  mpz_set_ui (k,1);
  mpz_init (q);

  for ( p = 2.0; p <= B1; p = getprime () ) {
    mpz_set_d (q ,p);
    for (r = p; r <= B1; r *= p) {
      mpz_mul (k ,k ,q);
    }
  }

  mpz_clear (q);

}


/*
  Compute k=lcm(2,..,B1) with the use of a product tree
  This is what is used
*/
void prodTreeCalculk (mpz_t k,double B1) {
  int stop = 0;
  int n=0;
  double r;
  mpz_t t;

  mpz_init (t);

  mpz_set_ui (k,1);
  for (r = 2.0; r <= B1; r *= 2.0) {
    mpz_mul_ui (k ,k ,2);
  }

  while (stop == 0) {
    prodTreeCalculkInter (t,B1,n,&stop);
    mpz_mul (k ,k ,t);
    n++;
  }

  mpz_clear (t);
}


void prodTreeCalculkInter (mpz_t Pr,double B1, int n, int* pstop) {
  double p,r;
  mpz_t t1,q;

  mpz_init (t1);
  mpz_init (q);

  if (n==0) {
    p = getprime();
    if ( p>B1 ) {
      mpz_set_ui (Pr,1);
      *pstop = 1;
    }
    else {
      mpz_set_ui (Pr,1);
      mpz_set_ui (q,p);
      for (r = p; r <= B1; r *= p) {
	mpz_mul (Pr,Pr,q);
      }
    }
  }
  else { // n>0
    if (*pstop == 1) {
      mpz_set_ui (Pr,1);
    }
    else {
      prodTreeCalculkInter (Pr,B1,n-1,pstop);
      prodTreeCalculkInter (t1,B1,n-1,pstop);
      mpz_mul (Pr,Pr,t1);
    }
  }

  mpz_clear(q);
  mpz_clear(t1);

}






/*
  compute k="lcm(2,..,B1)" with one power more for each prime
  Use a product tree.
  NOT USED
*/

void prodTreeCalculkPlus (mpz_t k,double B1) {
  int stop = 0;
  int n=0;
  double r;
  mpz_t t;

  mpz_init (t);

  for (r = 2.0; r <= B1; r *= 2.0) {
  }
  mpz_set_d (k,r);

  while (stop == 0) {
    prodTreeCalculkInterPlus (t,B1,n,&stop);
    mpz_mul (k ,k ,t);
    n++;
  }

  mpz_clear (t);
}



//  NOT USED
void prodTreeCalculkInterPlus (mpz_t Pr,double B1, int n, int* pstop) {
  double p,r;
  mpz_t t1;

  mpz_init (t1);

  if (n==0) {
    p = getprime();
    if ( p>B1 ) {
      mpz_set_ui (Pr,1);
      *pstop = 1;
    }
    else {
      for (r = p; r <= B1; r *= p) {
      }
      mpz_set_d (Pr ,r);
    }
  }
  else { // n>0
    if (*pstop == 1) {
      mpz_set_ui (Pr,1);
    }
    else {
      prodTreeCalculkInter (Pr,B1,n-1,pstop);
      prodTreeCalculkInter (t1,B1,n-1,pstop);
      mpz_mul (Pr,Pr,t1);
    }
  }


  mpz_clear(t1);

}






// *****************************************************************************








void mpalgpol_init (mpalgpol_t pol, mpmod_t n) {
  int i;
  for (i=0; i< DEGREE_ALGEBRA; i++) {
    mpres_init (pol->coeff[i],n);
  }
  mpres_init (pol->t1,n);
  mpalgres_init (pol->aTemp1,n);
  mpalgres_init (pol->aTemp2,n);
}

void mpalgpol_clear (mpalgpol_t pol, mpmod_t n) {
  int i;
  for (i=0; i< DEGREE_ALGEBRA; i++) {
    mpres_clear (pol->coeff[i],n);
  }
  mpres_clear (pol->t1,n);
  mpalgres_clear (pol->aTemp1,n);
  mpalgres_clear (pol->aTemp2,n);
}

// give the gcd of aP[i] with n
void mpalgres_gcd (mpz_t aF[DEGREE_ALGEBRA], mpalgres_t aP, mpmod_t n) {
  int i;

  for (i=0; i< DEGREE_ALGEBRA; i++) {
    mpres_gcd (aF[i],aP[i],n);
  } 
}




void mpalgres_init (mpalgres_t aP, mpmod_t n) {
  int i;

  for (i=0; i< DEGREE_ALGEBRA; i++) {
    mpres_init (aP[i],n);
  }
}

void mpalgres_clear (mpalgres_t aP, mpmod_t n) {
  int i;

  for (i=0; i< DEGREE_ALGEBRA; i++) {
    mpres_clear (aP[i],n);
  }
}

/* aR <- aP */
void mpalgres_set (mpalgres_t aR, mpalgres_t aP,mpmod_t n) {
  int i;
  for (i=0;i<DEGREE_ALGEBRA;i++) {
    mpz_set (aR[i],aP[i]);
  }
}

void mpalgres_set_zero (mpalgres_t aP,mpmod_t n) {
  int i;
  for (i=0;i<DEGREE_ALGEBRA;i++) {
    mpres_set_ui (aP[i],0,n);
  }
}

/* aP <- u mod n
   i.e. aP[0] <- u mod n   aP[i]<- 0
 */
void mpalgres_set_ui (mpalgres_t aP,unsigned long u,mpmod_t n) {
  mpalgres_set_zero(aP,n);
  mpres_set_ui (aP[0] ,u ,n);
}   



/* aP <- p 
   i.e. aP[0] <- p   aP[i]<- 0
 */
void mpalgres_set_mpres (mpalgres_t aP,mpres_t p,mpmod_t n) {
  mpalgres_set_zero(aP,n);
  mpz_set (aP[0] ,p);
}




int mpalgres_is_zero (mpalgres_t aP ,mpalgpol_t pol ,mpmod_t n) {
  int test=1;
  int i;
  for (i=0; i<DEGREE_ALGEBRA;i++) {
    test*=mpres_is_zero(aP[i],n);
  }
  return test;
}

/*  return the degree of the element aP 
    if aP=0 return -1 */
int mpalgres_degree (mpalgres_t aP ,mpalgpol_t pol ,mpmod_t n) {
  int i=DEGREE_ALGEBRA-1;
  while (i>=0 && mpres_is_zero(aP[i],n)) {
    i--;
  }
  return i;
}


/*  aR <- -aP  */
void mpalgres_neg (mpalgres_t aR, mpalgres_t aP, mpalgpol_t pol, mpmod_t n) {
  int i;
  for (i=0; i<DEGREE_ALGEBRA; i++) {
    mpres_neg (aR[i],aP[i],n);
  }
}

/*  aR <- -aP*u mod n
    i.e.  aR[i] <- aP[i]*u mod n
  */
void mpalgres_mul_ui (mpalgres_t aR, mpalgres_t aP, unsigned long u, mpalgpol_t pol, mpmod_t n) {
  int i;
  for (i=0; i<DEGREE_ALGEBRA; i++) {
    mpres_mul_ui (aR[i],aP[i],u,n);
  }
}

/*  aR <- -aP*p mod n
    i.e.  aR[i] <- aP[i]*p mod n
  */
void mpalgres_mul_mpres (mpalgres_t aR, mpalgres_t aP, mpres_t p, mpalgpol_t pol, mpmod_t n) {
  int i;
  for (i=0; i<DEGREE_ALGEBRA; i++) {
    mpres_mul (aR[i],aP[i],p,n);
  }
}

void mpalgres_add (mpalgres_t aR, mpalgres_t aP, mpalgres_t aQ, mpalgpol_t pol, mpmod_t n) {
  int i;
  for (i=0; i<DEGREE_ALGEBRA; i++) {
    mpres_add (aR[i],aP[i],aQ[i],n);
  }
}

void mpalgres_sub (mpalgres_t aR, mpalgres_t aP, mpalgres_t aQ, mpalgpol_t pol, mpmod_t n) {
  int i;
  for (i=0; i<DEGREE_ALGEBRA; i++) {
    mpres_sub (aR[i],aP[i],aQ[i],n);
  }
}

void mpalgres_add_mpres (mpalgres_t aR, mpalgres_t aP, mpres_t q, mpalgpol_t pol, mpmod_t n) {
  mpalgres_set (aR,aP,n);
  mpres_add (aR[0] ,aR[0] ,q ,n);
}

void mpalgres_sub_mpres (mpalgres_t aR, mpalgres_t aP, mpres_t q, mpalgpol_t pol, mpmod_t n) {
  mpalgres_set (aR,aP,n);
  mpres_sub (aR[0] ,aR[0] ,q ,n);

}

void mpalgres_add_ui (mpalgres_t aR, mpalgres_t aP, unsigned long u, mpalgpol_t pol, mpmod_t n) {
  mpalgres_set (aR,aP,n);
  mpres_add_ui (aR[0] ,aR[0] ,u ,n);
}

void mpalgres_sub_ui (mpalgres_t aR, mpalgres_t aP, unsigned long u, mpalgpol_t pol, mpmod_t n) {
  mpalgres_set (aR,aP,n);
  mpres_sub_ui (aR[0] ,aR[0] ,u ,n);
}




/* multiplication by X */
static void mpalgres_shift (mpalgres_t aR, mpalgres_t aP,
			    mpalgpol_t pol, mpmod_t n) {
  int i;

  mpz_set (pol->t1,aP[DEGREE_ALGEBRA-1]);
  for (i=DEGREE_ALGEBRA-1; i > 0; i--) {
    mpres_mul (aR[i] ,pol->t1, pol->coeff[i] ,n);
    mpres_sub (aR[i] ,aP[i-1], aR[i] ,n);
  }
  mpres_mul (aR[0] ,pol->t1 ,pol->coeff[0] ,n);
  mpres_neg (aR[0] ,aR[0] ,n);

}

/* aR <- aP*aQ */
void mpalgres_mul (mpalgres_t aR, mpalgres_t aP, mpalgres_t aQ,
		   mpalgpol_t pol, mpmod_t n) {
  int i;

  if ( (DEGREE_ALGEBRA==2) && mpres_is_zero(pol->coeff[1],n) ) {
    mpres_mul (pol->t1 ,aP[1] ,aQ[1] ,n);
    mpres_mul (pol->t1 ,pol->t1 ,pol->coeff[0] ,n);
    mpres_mul (pol->aTemp1[0] ,aP[0] ,aQ[0] ,n);
    mpres_sub (pol->aTemp1[0] ,pol->aTemp1[0] ,pol->t1 ,n);

    mpres_mul (pol->aTemp1[1] ,aP[0] ,aQ[1] ,n);
    mpres_mul (pol->t1 ,aP[1] ,aQ[0] ,n);
    mpres_add (pol->aTemp1[1] ,pol->aTemp1[1] ,pol->t1 ,n);
    
    mpalgres_set (aR ,pol->aTemp1 ,n);
  }
  else {
    mpalgres_set (pol->aTemp1 ,aP ,n);
    mpalgres_mul_mpres (aR ,aP ,aQ[0] ,pol,n);
    for (i=1;i<DEGREE_ALGEBRA;i++) {
      mpalgres_shift (pol->aTemp1 ,pol->aTemp1 ,pol ,n);
      mpalgres_mul_mpres (pol->aTemp2 ,pol->aTemp1 ,aQ[i] ,pol ,n);
      mpalgres_add (aR ,aR ,pol->aTemp2 ,pol ,n);
    }
  }
}


/*
  Do the inversion of aQ in k[x]/pol(x) where k=Z/nZ

  return -1 if something failed (for instance if  gcd(aQ,aP) != 1
  return 0 if it finds a factor of n (in fact if a mpres_invert failled)
  return 1 if ok

  Use f to put a factor of n if an inversion in k failled

  TODO do the general case


*/
int mpalgres_invert (mpalgres_t aV, mpalgres_t aQ,
		     mpalgpol_t pol, mpmod_t n, mpz_t f) {

  mpres_t temp;


  if (  (DEGREE_ALGEBRA==2) && mpres_is_zero(pol->coeff[1],n) ) {
    if (mpalgres_is_zero (aQ,pol,n) ) {
      return -1;
    }


    mpz_set ( pol->aTemp1[0] ,aQ[0]);
    mpres_neg ( pol->aTemp1[1] ,aQ[1] ,n); 

    mpres_mul (pol->t1 ,aQ[0] ,aQ[0] ,n);
    mpres_mul (pol->aTemp2[0] ,aQ[1] ,aQ[1] ,n);
    mpres_mul (pol->aTemp2[0] ,pol->aTemp2[0] ,pol->coeff[0] ,n);
    mpres_add (pol->t1 ,pol->t1 ,pol->aTemp2[0] ,n);



    if ( !mpres_invert(pol->t1,pol->t1,n) ) {
      mpres_gcd(f,pol->t1,n);

      mpres_init (temp,n);
      mpres_set_z (temp ,f ,n);
      if ( mpres_is_zero (temp ,n) ) {
	mpres_clear (temp ,n);
	return -1;
      }
      else {
	mpres_clear (temp ,n);
	return 0;
      }
    }

    mpz_set ( aV[0] ,aQ[0]);
    mpres_neg ( aV[1] ,aQ[1] ,n); 
    mpalgres_mul_mpres (aV ,aV ,pol->t1, pol ,n);

    return 1;
    
  }
  else { // TODO do the general case
    return -1;
  
  }
}
