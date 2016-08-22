#include "Jacobi.h"





void coorJacobi_init (coorJacobi P, mpmod_t n ) {
  mpres_init (P->U ,n);
  mpres_init (P->V ,n);
  mpres_init (P->W ,n);
  mpres_init (P->Y ,n);
}

void coorJacobi_clear (coorJacobi P, mpmod_t n ) {
  mpres_clear (P->U ,n);
  mpres_clear (P->V ,n);
  mpres_clear (P->W ,n);
  mpres_clear (P->Y ,n);
}

void coorJacobi_set (coorJacobi P,coorJacobi Q, mpmod_t n ) {
  mpres_set (P->U ,Q->U ,n);
  mpres_set (P->V ,Q->V ,n);
  mpres_set (P->W ,Q->W ,n);
  mpres_set (P->Y ,Q->Y ,n);
}



/*
  Multiply a point by k on a Jacobi elliptic curve over ZZ (use mpz_t)
  The curve is y^2 = 1 + (-3/s^2+1/s^4)*x^2 + x^4/s^2
  The initial point is X=1, Y=1-ep  (Z=1)
  The result is put in (x,y,z)
  NOTE: for k=0 or k=1 we take k=2!
*/
void mulJacobiEntiers (mpz_t a,mpz_t b,
		       int k,
		       mpz_t x,mpz_t y,mpz_t z) {

  
  if (k<0) {
    k = -k;
  }

  if (k==0) {
    mpz_set_ui (x,0);
    mpz_set_ui (y,1);
    mpz_set_ui (z,1);
  }
  else if (k==1) {
    mpz_mul (x ,a ,a);
    mpz_mul (z ,b ,b);
    mpz_sub (y ,x ,z);

    mpz_set (x,a);
    mpz_set (z,a);

  }
  else if (k ==2) {
    mpz_mul (y ,b ,b);
    mpz_mul_ui (y ,y ,2);
    mpz_mul (z ,a ,a);
    mpz_add (y ,y ,z);    // y = a^2 + 2*b^2
    mpz_mul_ui (x ,a ,2); // x = 2*a
    mpz_set (z ,a);       // z = a


  }
  else { // general case
    mpz_t sa,sb,Yi;
    mpz_t X,Y,Z;
    mpz_t t1,t2,t3,t4,t5,t7,t8,t9;
    int i;
    
    mpz_init (sa);
    mpz_init (sb);
    mpz_init (Yi);
    mpz_init (X);
    mpz_init (Z);
    mpz_init (Y);
    mpz_init (t1);
    mpz_init (t2);
    mpz_init (t3);
    mpz_init (t4);
    mpz_init (t5);
    mpz_init (t7);
    mpz_init (t8);
    mpz_init (t9);
  
    mpz_mul (sa ,a ,a);
    mpz_mul (sb ,b ,b);
    mpz_sub (Yi ,sa ,sb);

    mpz_set (X ,a);
    mpz_set (Z ,a);
    mpz_set (Y ,Yi);

    for (i=2;i<=k;i++) {
      mpz_mul (t1 ,X ,Z);
      mpz_mul (t2 ,X ,X);
      mpz_mul (t3 ,Z ,Z);
      mpz_mul (t4 ,t3 ,sa);
      mpz_mul (t5 ,t2 ,sb);
      
      
      mpz_add (t9 ,t2 ,t3);
      mpz_mul (t9 ,t9 ,sa);
      mpz_mul (t9 ,t9 ,sb);
      mpz_mul (t9 ,t9 ,t1);
      mpz_mul_ui (t9 ,t9 ,2);
      mpz_mul (t9 ,t9 ,sa);
      
      mpz_mul (t8 ,Yi ,Y);
      mpz_mul (t8 ,t8 ,sa);
      
      mpz_mul_ui (t7 ,sa ,3);
      mpz_sub (t7 ,sb ,t7);
      mpz_mul (t7 ,t7 ,sb);
      mpz_mul (t7 ,t7 ,t1);
      
      mpz_add (t8 ,t8 ,t7);
      
      mpz_add (t7 ,t4 ,t5);
      mpz_mul (t8 ,t8 ,t7);
      mpz_add (t9 ,t9 ,t8);
      
      
      mpz_sub (t8 ,t4 ,t5);
      mpz_mul (t8 ,t8 ,a);
      
      
      mpz_mul (t7 ,t1 ,Yi);
      mpz_mul (t2 ,sa ,Y);
      mpz_add (t7 ,t7 ,t2);
      mpz_mul (t7 ,t7 ,a);
      
      mpz_set (X ,t7);
      mpz_set (Y ,t9);
      mpz_set (Z ,t8);
    }
    mpz_set (x ,X);
    mpz_set (y ,Y);
    mpz_set (z ,Z);
    

    mpz_clear (sa);
    mpz_clear (sb);
    mpz_clear (Yi);
    mpz_clear (X);
    mpz_clear (Z);
    mpz_clear (Y);
    mpz_clear (t1);
    mpz_clear (t2);
    mpz_clear (t3);
    mpz_clear (t4);
    mpz_clear (t5);
    mpz_clear (t7);
    mpz_clear (t8);
    mpz_clear (t9);
  }
   

}




/*
  Multiply a point by k on a Jacobi elliptic curve 
  The curve is y^2 = 1 + (-3/s^2+1/s^4)*x^2 + x^4/s^2
  The initial point is X=1, Y=1-ep  (Z=1)
  The result is put in x,y
  NOTE: for k=0 or k=1 we take k=2!
*/
int mulJacobi2 (mpz_t f,mpmod_t n,
		int k,
		mpres_t x,mpres_t y,
		mpres_t Y,mpres_t ep,mpres_t dep) {


  if (k <= 2) {

    doubleJacobi2DebFin(f,n,x,y,Y,ep,dep);
    return MULT_JACOBI;

  }
  else if (k == 3) {

    coorJacobi P,Q;
    mpres_t t;

    mpres_init (t,n);
    coorJacobi_init (P,n);
    coorJacobi_init (Q,n);

    mpres_set_ui (P->U,1,n);
    mpres_set_ui (P->V,1,n);
    mpres_set_ui (P->W,1,n);
    mpres_set    (P->Y,Y,n);

    doubleJacobi2Deb(n,Q,Y,ep,dep);
    addJacobi2fin(n,Q,P,x,y,t,ep,dep);

    if (!mpres_invert (t,t,n)) // t=1/Z3
    {
      mpres_gcd (f, t, n);

      mpres_clear (t,n);
      coorJacobi_clear (P,n);
      coorJacobi_clear (Q,n);

      return MULT_JACOBI_FAIL;
    }

    mpres_mul (x ,x ,t ,n);
    mpres_mul (t ,t ,t ,n);
    mpres_mul (y ,y ,t ,n);

    mpres_clear (t,n);
    coorJacobi_clear (P,n);
    coorJacobi_clear (Q,n);

    return MULT_JACOBI;

  }
  else {

    coorJacobi P,Q;
    mpres_t t;
    int mask;

    mpres_init (t,n);
    coorJacobi_init (P,n);
    coorJacobi_init (Q,n);

    mpres_set_ui (P->U,1,n);
    mpres_set_ui (P->V,1,n);
    mpres_set_ui (P->W,1,n);
    mpres_set    (P->Y,Y,n);

    
    mask = (1<<POW_MAX);// mask = 10000 with at least one more bit than k
    while ( (mask & k) == 0 ) {
      mask = (mask >> 1);
    } // mask = 100... with the same number of bits than k
    mask = (mask >> 1);// mask = 010... with the same number of bits than k


    doubleJacobi2Deb(n,Q,Y,ep,dep); // Q = 2*P
    if ( (mask & k) != 0 ) {    // case Q+2*Q
      addJacobi2(n,P,Q,Q,ep,dep);
    }
    mask = (mask >> 1);// mask = 0010... with the same number of bits than k

    while (mask > 1) {
      addJacobi2(n,Q,Q,Q,ep,dep);
      if ( (mask & k) != 0 ) {    // case Q+2*Q
	addJacobi2(n,P,Q,Q,ep,dep);
      }
      mask = (mask >> 1);
    }


    if ( (mask & k) == 0 ) { // case 2*Q
      coorJacobi_set (P ,Q ,n); 
      addJacobi2fin(n,Q,P,x,y,t,ep,dep);
    }
    else { // case 2*Q+Q
      addJacobi2(n,Q,Q,Q,ep,dep);
      addJacobi2fin(n,Q,P,x,y,t,ep,dep);
    }


    if (!mpres_invert (t,t,n)) // t=1/Z3
    {
      mpres_gcd (f, t, n);

      mpres_clear (t,n);
      coorJacobi_clear (P,n);
      coorJacobi_clear (Q,n);

      return MULT_JACOBI_FAIL;
    }

    mpres_mul (x ,x ,t ,n);
    mpres_mul (t ,t ,t ,n);
    mpres_mul (y ,y ,t ,n);

    mpres_clear (t,n);
    coorJacobi_clear (P,n);
    coorJacobi_clear (Q,n);

    return MULT_JACOBI;

  }

}




/*
  Double a point on the Jacobi curve y^2=1-dep*X^2+ep*X^4
  begin with X=1, Y=1-ep (Z=1)
  We want x,y
*/
void doubleJacobi2DebFin(mpz_t f,mpmod_t n,
			mpres_t x,mpres_t y,
			mpres_t Y,mpres_t ep,mpres_t dep) {




  mpres_set_ui (x ,2 ,n); // x=2

  mpres_mul_ui (y ,ep ,2 ,n);
  mpres_add_ui (y ,y ,1 ,n); // y = 2*ep+1 = 2/s^2+1


}




/*
  Double a point on the Jacobi curve y^2=1-dep*X^2+ep*X^4
  begin with X=1, Y=1-ep (Z=1)
  We want P2=(U3,V3,W3,Y3)
*/
void doubleJacobi2Deb(mpmod_t n,
		      coorJacobi P2,
		      mpres_t Y,mpres_t ep,mpres_t dep) {

  mpres_add_ui (P2->V ,ep, 1 ,n); // V3 = ep+1
  mpres_mul (P2->W ,Y ,Y ,n); // W3 = Y^2
  mpres_sub (P2->Y ,P2->W ,dep ,n); // Y3=Y^2-dep
  mpres_mul (P2->Y ,P2->Y ,P2->V ,n); // Y3=(1+ep)*(Y^2-dep)
  mpres_mul_ui (P2->V ,ep ,4 ,n); 
  mpres_add (P2->Y ,P2->Y ,P2->V ,n); // Y3=(1+ep)*(Y^2-dep) + 4ep

  mpres_mul_ui (P2->V ,P2->W ,2 ,n); // V3 = 2*Y^2
  mpres_mul_ui (P2->U ,P2->V ,2 ,n); // U3 = 4*Y^2

}



/*
  add two points on the jacobi curve y^2=1-dep*X^2+ep*X^4
  Initial points P1=(U1,V1,W1,Y1) and P2=(U2,V2,W2,Y2)
  P1+P2=P3=(U3,V3,W3,Y3)
*/
void addJacobi2(mpmod_t n,
		coorJacobi P1,coorJacobi P2,
		coorJacobi P3,
		mpres_t ep,mpres_t dep) {

  mpres_t t1,t3,t5,t7,t9;

  mpres_init (t1,n);
  mpres_init (t3,n);
  mpres_init (t5,n);
  mpres_init (t7,n);
  mpres_init (t9,n);



  mpres_set (t1    ,P1->U ,n);
  mpres_set (P3->U ,P2->U ,n);
  mpres_set (t3    ,P1->V ,n);
  mpres_set (P3->V ,P2->V ,n);
  mpres_set (t5    ,P1->W ,n);
  mpres_set (P3->W ,P2->W ,n);
  mpres_set (t7    ,P1->Y ,n);
  mpres_set (P3->Y ,P2->Y ,n);


  mpres_mul (t9    ,t7    ,P3->Y ,n);
  mpres_add (t7    ,t7    ,t3    ,n);
  mpres_add (P3->Y ,P3->Y ,P3->V ,n);
  mpres_mul (t3    ,t3    ,P3->V ,n);
  mpres_mul (t7    ,t7    ,P3->Y ,n);
  mpres_sub (t7    ,t7    ,t9    ,n);
  mpres_sub (t7    ,t7    ,t3    ,n); // X3
  mpres_mul (P3->V ,t1    ,P3->U ,n);
  mpres_mul (P3->Y ,t5    ,P3->W ,n);
  mpres_add (t1    ,t1    ,t5    ,n);
  mpres_add (P3->U ,P3->U ,P3->W ,n);
  mpres_mul (t5    ,t1    ,P3->U ,n);
  mpres_sub (t5    ,t5    ,P3->V ,n);
  mpres_sub (t5    ,t5    ,P3->Y ,n);
  mpres_mul (P3->V ,P3->V ,ep    ,n);
  mpres_sub (t1    ,P3->Y ,P3->V ,n); // Z3
  mpres_add (P3->U ,P3->Y ,P3->V ,n);
  mpres_mul (P3->W ,t3    ,dep   ,n);
  mpres_sub (P3->W ,t9    ,P3->W ,n);
  mpres_mul (P3->W ,P3->W ,P3->U ,n);
  mpres_mul_ui (t3 ,t3    ,2     ,n);
  mpres_mul (t3    ,t3    ,ep    ,n);
  mpres_mul (t3    ,t3    ,t5    ,n);
  mpres_add (P3->Y ,P3->W ,t3    ,n); // Y3



  mpres_mul (P3->U ,t7 ,t7 ,n); // U3
  mpres_mul (P3->V ,t1 ,t7 ,n); // V3
  mpres_mul (P3->W ,t1 ,t1 ,n); // W3


  mpres_clear (t1,n);
  mpres_clear (t3,n);
  mpres_clear (t5,n);
  mpres_clear (t7,n);
  mpres_clear (t9,n);
}



/*
  add two points on the jacobi curve y^2=1-dep*X^2+ep*X^4
  Initial points P1=(U1,V1,W1,Y1) and P2=(U2,V2,W2,Y2)
  P1+P2=P3=(U3,V3,W3,Y3)
  WARNING: value of P1 and P2 are modified
           P1 and P2 must be different
*/
void addJacobi2fin(mpmod_t n,
		   coorJacobi P1,coorJacobi P2,
		   mpres_t X3,mpres_t Y3,mpres_t Z3,
		   mpres_t ep,mpres_t dep) {

  mpres_mul (Y3    ,P1->Y ,P2->Y ,n);
  mpres_add (P1->Y ,P1->Y ,P1->V ,n);
  mpres_add (P2->Y ,P2->Y ,P2->V ,n);
  mpres_mul (P1->V ,P1->V ,P2->V ,n);
  mpres_mul (P1->Y ,P1->Y ,P2->Y ,n);
  mpres_sub (P1->Y ,P1->Y ,Y3    ,n);
  mpres_sub (X3    ,P1->Y ,P1->V ,n); // X3
  mpres_mul (P2->V ,P1->U ,P2->U ,n);
  mpres_mul (P2->Y ,P1->W ,P2->W ,n);
  mpres_add (P1->U ,P1->U ,P1->W ,n);
  mpres_add (P2->U ,P2->U ,P2->W ,n);
  mpres_mul (P1->W ,P1->U ,P2->U ,n);
  mpres_sub (P1->W ,P1->W ,P2->V ,n);
  mpres_sub (P1->W ,P1->W ,P2->Y ,n);
  mpres_mul (P2->V ,P2->V ,ep    ,n);
  mpres_sub (Z3    ,P2->Y ,P2->V ,n); // Z3
  mpres_add (P2->U ,P2->Y ,P2->V ,n);
  mpres_mul (P2->W ,P1->V ,dep   ,n);
  mpres_sub (P2->W ,Y3    ,P2->W ,n);
  mpres_mul (P2->W ,P2->W ,P2->U ,n);

  mpres_mul_ui (P1->V ,P1->V ,2 ,n);

  mpres_mul (P1->V ,P1->V ,ep    ,n);
  mpres_mul (P1->V ,P1->V ,P1->W ,n);
  mpres_add (Y3    ,P2->W ,P1->V ,n); // Y3
}


