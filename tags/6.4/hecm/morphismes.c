#include "morphismes.h"
#include "auxi.h"
#include "ariKS.h"
#include "Jacobi.h"


void DivMumfordU_init (DivMumfordU DivU, mpmod_t n ) {
  mpres_init (DivU->u0,n);
  mpres_init (DivU->u1,n);
  DivU->degree = 0;
}


void DivMumfordU_clear (DivMumfordU DivU, mpmod_t n ) {
  mpres_clear (DivU->u0,n);
  mpres_clear (DivU->u1,n);
}


void DivMumfordV_init (DivMumfordV DivV, mpmod_t n ) {
  mpres_init (DivV->V0,n);
  mpres_init (DivV->V1,n);
  mpres_init (DivV->v1v0,n);
}

void DivMumfordV_clear (DivMumfordV DivV, mpmod_t n ) {
  mpres_clear (DivV->V0,n);
  mpres_clear (DivV->V1,n);
  mpres_clear (DivV->v1v0,n);
}



/*
  Let P be a point on the Kummer surface 
  We want to obtain the Mumford polynomial of the corresponding divisor
    A divisor and its opposite have the same image on the Kummer surface
    Thus we can't have the exact v polynomial of P
    Instead we obtain V=v^2
  There is 3 cases
    if degree u=0 then we have the divisor 0
    if degree u=1 we only need the polynomial u for the following
    if degree u=2 we use the function DivMumfordDegree2
*/
int DivMumford (mpz_t f,mpmod_t n,DivMumfordU divU,DivMumfordV divV,
		ksPoint P,thetaCst th,curveHyperEll cHEll) {

  mpres_t t1,t2;
  mpres_t T13p, T14p,T16p;
  int test;

  mpres_init (t1,n);
  mpres_init (t2,n);

  mpres_init (T13p,n);
  mpres_init (T14p,n);
  mpres_init (T16p,n);

  // Construction of 3 theta functions up to constants
  //  T13p := X*Rac     - Y*t10p*Rac + Z*t7p     - T*t7p*t10p;
  //  T14p := X*t5p*Rac - Y*t6p*t7p  + Z*t5p*t7p - T*t6p*Rac;
  //  T16p := X*t6p*t7p + Y*t5p*Rac  - Z*t6p*Rac + T*t5p*t7p;
  mpres_mul (T13p ,P->Z ,th->t7p ,n); // Z*t7p
  mpres_mul (T14p ,T13p ,th->t5p ,n); // Z*t5p*t7p
  mpres_mul (t1 ,P->Y ,th->Rac ,n);
  mpres_mul (T16p ,th->t5p ,t1 ,n); // Y*t5p*Rac
  mpres_mul (t1 ,t1 ,th->t10p ,n);
  mpres_sub (T13p ,T13p ,t1 ,n); // Z*t7p - Y*t10p*Rac
  mpres_mul (t1 ,P->T ,th->t7p ,n);
  mpres_mul (t2 ,t1 ,th->t10p ,n);
  mpres_sub (T13p ,T13p ,t2 ,n); // Z*t7p - Y*t10p*Rac - T*t7p*t10p
  mpres_mul (t2 ,t1 ,th->t5p ,n);
  mpres_add (T16p ,T16p ,t2 ,n);  // T*t5p*t7p + Y*t5p*Rac
  mpres_mul (t1 ,P->X ,th->Rac ,n);
  mpres_add (T13p ,T13p ,t1 ,n); // Z*t7p - Y*t10p*Rac - T*t7p*t10p +X*Rac
  mpres_mul (t2 ,t1 ,th->t5p ,n);
  mpres_add (T14p ,T14p ,t2 ,n); // Z*t5p*t7p + X*t5p*Rac
  mpres_mul (t1 ,th->t6p ,th->Rac ,n);
  mpres_mul (t2 ,t1 ,P->T ,n);
  mpres_sub (T14p ,T14p ,t2 ,n); // Z*t5p*t7p + X*t5p*Rac - T*t6p*Rac
  mpres_mul (t2 ,t1 ,P->Z ,n);
  mpres_sub (T16p ,T16p ,t2 ,n); // T*t5p*t7p + Y*t5p*Rac - Z*t6p*Rac
  mpres_mul (t1 ,th->t6p ,th->t7p ,n);
  mpres_mul (t2 ,t1 ,P->Y ,n);
  mpres_sub (T14p ,T14p ,t2 ,n);// Z*t5p*t7p + X*t5p*Rac - T*t6p*Rac - Y*t6p*t7p
  mpres_mul (t2 ,t1 ,P->X ,n);
  mpres_sub (T16p ,T16p ,t2 ,n);// T*t5p*t7p + Y*t5p*Rac - Z*t6p*Rac + X*t6p*t7p
  // We have finished T13p,T14p et T16p
  // Note that contrary to the papper, we don't divided by t7p^2-Rac^2



  if ( !mpres_is_zero(T16p,n) ) { // T16p != 0

    divU->degree =2;

    if (!mpres_invert (T16p , T16p, n)) // T16p=1/T16p
      {
	mpres_gcd (f, T16p, n);
	

	mpres_clear (t1,n);
	mpres_clear (t2,n);

	mpres_clear (T13p,n);
	mpres_clear (T14p,n);
	mpres_clear (T16p,n);
	return MORPHISM_FAIL;
      }



    test=DivMumfordDegree2(f,n,divU,divV,P,th,cHEll,T13p,T14p, T16p);

    mpres_clear (t1,n);
    mpres_clear (t2,n);
    
    mpres_clear (T13p,n);
    mpres_clear (T14p,n);
    mpres_clear (T16p,n);
    
    return test;

  }
  else { // T16p = 0

    /* 
       denom := (la-1)*t5p*T13p - la*t8p*T14p;
       u0 := la*T14p / denom;
    */


    mpres_sub_ui (t1 ,cHEll->la ,1 ,n);
    mpres_mul (t1 ,t1 ,th->t5p ,n);
    mpres_mul (t1 ,t1 ,T13p ,n); // (la-1)*t5p*T13p
    mpres_mul (t2 ,cHEll->la ,T14p ,n); // la*t8p*t14p;
    mpres_sub (t1, t1 ,t2 ,n); // t1 = denom = (la-1)*t5p*T13p - la*t8p*T14p


    if ( !mpres_is_zero(t1,n) ) { // denom != 0

      if (!mpres_invert (t2 , t1, n)) // t2=1/( (la-1)*t5p*T13p - la*t8p*t14p )
	{
	  mpres_gcd (f, t1, n);

	  mpres_clear (t1,n);
	  mpres_clear (t2,n);
	  
	  mpres_clear (T13p,n);
	  mpres_clear (T14p,n);
	  mpres_clear (T16p,n);

	return MORPHISM_FAIL;
      }

      divU->degree =1;
      
      mpres_mul (divU->u0 ,cHEll->la ,t2 ,n); 
      mpres_mul (divU->u0 ,divU->u0 ,T14p ,n); 
      
      // we don't need V.

      mpres_clear (t1,n);
      mpres_clear (t2,n);
      
      mpres_clear (T13p,n);
      mpres_clear (T14p,n);
      mpres_clear (T16p,n);
      
      return MORPHISM;
    
    }
    else { // denom = 0

      divU->degree =0;
      
      mpres_clear (t1,n);
      mpres_clear (t2,n);
      
      mpres_clear (T13p,n);
      mpres_clear (T14p,n);
      mpres_clear (T16p,n);

    return MORPHISM;

    }

  }

}










/*
  Get the Mumford polynomial in the case of degree u=2.
  We only get V=v^2
*/
int DivMumfordDegree2 (mpz_t f,mpmod_t n,DivMumfordU divU,DivMumfordV divV,
		   ksPoint P,thetaCst th,curveHyperEll cHEll,
		   mpres_t T13p, mpres_t T14p, mpres_t T16p) {
  
  mpres_t t1,t2,t3,t4,t5;
  mpres_t T7p,T9p,T11p,T12p;

  mpres_init (t1,n);
  mpres_init (t2,n);
  mpres_init (t3,n);
  mpres_init (t4,n);
  mpres_init (t5,n);

  mpres_init (T7p,n);
  mpres_init (T9p,n);
  mpres_init (T11p,n);
  mpres_init (T12p,n);

  /* 
T7p := ( - X*t5p      + Y*t5p*t10p - Z*t6p*t10p + T*t6p      ) / (t6p^2-t5p^2);
T9p := (   X*t6p*t10p - Y*t6p      + Z*t5p      - T*t5p*t10p ) / (t8p^2-t10p^2);
T11p:= (   X*t6p      - Y*t6p*t10p + Z*t5p*t10p - T*t5p      ) / (t6p^2-t5p^2);
T12p:= ( - X*t5p*t10p + Y*t5p      - Z*t6p      + T*t6p*t10p ) / (t8p^2-t10p^2);
  */
  
  mpres_mul (T7p ,P->T ,th->t6p ,n);  // T*t6p
  mpres_mul (T12p ,T7p ,th->t10p ,n); // T*t6p*t10p
  mpres_mul (t1 ,P->Y ,th->t5p ,n); 
  mpres_add (T12p ,T12p ,t1 ,n); // T*t6p*t10p + Y*t5p
  mpres_mul (t1 ,t1 ,th->t10p ,n); 
  mpres_add (T7p ,T7p ,t1 ,n);   // T*t6p + Y*t5p*t10p 
  mpres_mul (t1 ,P->Z ,th->t6p ,n); 
  mpres_sub (T12p ,T12p ,t1 ,n); // T*t6p*t10p + Y*t5p - Z*t6p
  mpres_mul (t1 ,t1 ,th->t10p ,n); 
  mpres_sub (T7p ,T7p ,t1 ,n);   // T*t6p + Y*t5p*t10p - Z*t6p*t10p
  mpres_mul (t1 ,P->X ,th->t5p ,n); 
  mpres_sub (T7p ,T7p ,t1 ,n);   // T*t6p + Y*t5p*t10p - Z*t6p*t10p - X*t5p
  mpres_mul (t1 ,t1 ,th->t10p ,n); 
  mpres_sub (T12p ,T12p ,t1 ,n); // T*t6p*t10p + Y*t5p - Z*t6p - X*t5p*t10p
  // We finished the numerators of T7p and T12p  

  mpres_mul (T9p ,P->Z ,th->t5p ,n);  // Z*t5p
  mpres_mul (T11p ,T9p ,th->t10p ,n); // Z*t5p*t10p
  mpres_mul (t1 ,P->X ,th->t6p ,n);
  mpres_add (T11p ,T11p ,t1 ,n); // Z*t5p*t10p + X*t6p
  mpres_mul (t1 ,t1 ,th->t10p ,n);
  mpres_add (T9p ,T9p ,t1 ,n);   // Z*t5p + X*t6p*t10p
  mpres_mul (t1 ,P->T ,th->t5p ,n);
  mpres_sub (T11p ,T11p ,t1 ,n); // Z*t5p*t10p + X*t6p - T*t5p
  mpres_mul (t1 ,t1 ,th->t10p ,n);
  mpres_sub (T9p ,T9p ,t1 ,n);   // Z*t5p + X*t6p*t10p - T*t5p*t10p
  mpres_mul (t1 ,P->Y ,th->t6p ,n);
  mpres_sub (T9p ,T9p ,t1 ,n);   // Z*t5p + X*t6p*t10p - T*t5p*t10p - Y*t6p
  mpres_mul (t1 ,t1 ,th->t10p ,n);
  mpres_sub (T11p ,T11p ,t1 ,n); // Z*t5p*t10p + X*t6p - T*t5p - Y*t6p*t10p
  // We finished the numerators of T9p and T11p  

  

  mpres_mul (t1, th->t6p ,th->t6p ,n);
  mpres_mul (t2 ,th->t5p ,th->t5p ,n);
  mpres_sub (t2 ,t1 ,t2 ,n); // t6p^2-t5p^2
  if (!mpres_invert (t1 , t2, n)) // t1=1 / (t6p^2-t5p^2)
    {
      mpres_gcd (f, t2, n);
	
      mpres_clear (t1, n);
      mpres_clear (t2, n); 
      mpres_clear (t3, n);
      mpres_clear (t4, n);
      mpres_clear (t5, n);

      mpres_clear (T7p,n);
      mpres_clear (T9p,n);
      mpres_clear (T11p,n);
      mpres_clear (T12p,n);

      return MORPHISM_FAIL;
    }
  mpres_mul (T7p ,T7p ,t1 ,n);
  mpres_mul (T11p ,T11p ,t1 ,n);
  // we finished T7p et T11p
  
  mpres_mul (t1, th->t10p ,th->t10p ,n);
  mpres_ui_sub (t2 ,1 ,t1 ,n);
  if (!mpres_invert (t1 , t2, n)) // t1=1 / (t8p^2-t10p^2)
    {
      mpres_gcd (f, t2, n);
      
      mpres_clear (t1, n);
      mpres_clear (t2, n); 
      mpres_clear (t3, n);
      mpres_clear (t4, n);
      mpres_clear (t5, n);

      mpres_clear (T7p,n);
      mpres_clear (T9p,n);
      mpres_clear (T11p,n);
      mpres_clear (T12p,n);

      return MORPHISM_FAIL;
    }
  mpres_mul (T9p ,T9p ,t1 ,n);
  mpres_mul (T12p ,T12p ,t1 ,n);
  // we finished T9p et T12p

  


  
   
  // u0 := nu*be*la*T14p / T16p;
  // u1 := nu*be*(la-1)*t5p*T13p / T16p   - u0-1 ;
  mpres_mul (divU->u0 ,T16p ,cHEll->nu ,n); // we already have inversed T16p
  mpres_mul (divU->u0 ,divU->u0 ,th->be ,n); // nu*be/T16p
  mpres_mul (divU->u1 ,divU->u0 ,th->t5p ,n); // u1 = nu*be*t5p/T16p
  mpres_mul (divU->u0 ,divU->u0 ,cHEll->la ,n);
  mpres_mul (divU->u0 ,divU->u0 ,T14p ,n); // u0 = la*nu*be*T14p / T16p
  mpres_mul (divU->u1 ,divU->u1 ,T13p ,n); // nu*be*t5p*T13p/T16p
  mpres_sub_ui (t1 ,cHEll->la ,1 ,n);
  mpres_mul (divU->u1 ,divU->u1 ,t1 ,n); // nu*be*(la-1)*t5p*T13p/T16p
  mpres_sub (divU->u1 ,divU->u1 ,divU->u0 ,n);
  mpres_sub_ui (divU->u1 ,divU->u1 ,1 ,n); // u1=(la-1)*nu*be*t5p*T14p/T16p-u0-1
  




  
  // V0:= (p*ga^2*nu^3*T14p/T16p^3) * (t10p^2 + nu^2*be^2 -2)  *  ( A*D*(be*ga*(T7p*T12p*Rac-2*(XZ+YT))-t7p*t10p*T9p*T11p) + (ga+be)*((1-be*ga)*((be+ga)*(X+T)*(Y+Z)-(X^2+Y^2+Z^2+T^2))+(2-be^2-ga^2)*(XT+YZ)) );
    
  mpres_mul (t1 ,th->t10p ,th->t10p ,n);
  mpres_mul (divV->V0 ,cHEll->nu, th->be ,n);
  mpres_mul (divV->V0 ,divV->V0 ,divV->V0 ,n);
  mpres_add (divV->V0 ,t1 ,divV->V0 ,n);
  mpres_sub_ui (divV->V0 ,divV->V0 ,2 ,n); // ( t10p^2 + nu^2*be^2 -2 )
  
  mpres_mul (t1 ,cHEll->nu, T16p ,n); // On avait deja inverse T16p
  mpres_mul (t2 ,t1 ,th->ga ,n);
  mpres_mul (t2 ,t2 ,t2 ,n);
  mpres_mul (t2 ,t2 ,t1 ,n); // ga^2 * nu^3 / T16p^3
  mpres_mul (divV->V0 ,divV->V0 ,t2 ,n);
  mpres_mul (divV->V0 ,divV->V0 ,T14p ,n);
  // V0 = (t10p^2+nu^2*be^2-2)*ga^2*nu^3*T14p / T16p^3
  mpres_mul (divV->V0 ,divV->V0 ,th->p ,n);
  // V0 = p*(t10p^2+nu^2*be^2-2)*ga^2*nu^3*T14p / T16p^3


  mpres_mul (t1, P->X ,P->Z ,n);
  mpres_mul (t2, P->Y ,P->T ,n);
  mpres_add (t2 ,t2 ,t1 ,n); // XZ+YT
  mpres_mul_ui (t2 ,t2 ,2 ,n);
  mpres_mul (t1, T7p ,T12p ,n);
  mpres_mul (t1, t1 ,th->Rac ,n); // Rac*T7p*T12p
  mpres_sub (t1 ,t1 ,t2 ,n); //  Rac*T7p*T12p -2*(XZ+YT)
  mpres_mul (t1 ,t1 ,th->be ,n);
  mpres_mul (t1 ,t1 ,th->ga ,n); // be*ga* ( T7p*T12p*Rac - 2*(XZ+YT) )
  
  mpres_mul (t2 ,T9p ,T11p ,n);
  mpres_mul (t2 ,t2 ,th->t7p ,n);
  mpres_mul (t2 ,t2 ,th->t10p ,n); // al*de*t7p^2*T9p*T11p*R=t7p*t10p*T9p*T11p
  mpres_sub (t1 ,t1 ,t2 ,n);
    // t1 = be*ga* ( T7p*T12p*Rac - 2*(XZ+YT) ) - t7p*t10p*T9p*T11p
  mpres_add (t4 ,th->be ,th->ga ,n); // t4=be+ga
  mpres_add_ui (t2 ,t4 ,2 ,n); // t2 = A = 2+be+ga
  mpres_ui_sub (t3 ,2 ,t4 ,n); // t3 = D = 2-be-ga
  mpres_mul (t2 ,t3 ,t2 ,n);
  mpres_mul (t1 ,t1 ,t2 ,n);
  // t1 = A*D*( be*ga* ( T7p*T12p*Rac - 2*(XZ+YT) ) - t7p*t10p*T9p*T11p )

  mpres_add (t2 ,P->Z ,P->Y ,n);
  mpres_add (t3 ,P->X, P->T ,n);
  mpres_mul (t2 ,t2 ,t3 ,n); // (X+T)*(Y+Z)
  mpres_mul (t2 ,t2 ,t4 ,n); // t2 = (be+ga)*(X+T)*(Y+Z)
  mpres_mul (t3 ,P->X ,P->X ,n);
  mpres_mul (t4 ,P->Y ,P->Y ,n);
  mpres_add (t3 ,t3 ,t4 ,n);
  mpres_mul (t4 ,P->Z ,P->Z ,n);
  mpres_add (t3 ,t3 ,t4 ,n);
  mpres_mul (t4 ,P->T ,P->T ,n);
  mpres_add (t3 ,t3 ,t4 ,n); // t3 = X^2+Y^2+Z^2+T^2
  mpres_sub (t2 ,t2 ,t3 ,n);
  mpres_mul (t3 ,th->be ,th->ga ,n);
  mpres_ui_sub (t3 ,1 ,t3 ,n);
  mpres_mul (t2 ,t2 ,t3 ,n);
  // t2 = (1-be*ga) * ( (be+ga)*(X+T)*(Y+Z) - (X^2+Y^2+Z^2+T^2))
  
  mpres_mul (t3 ,th->be, th->be, n);
  mpres_mul (t4 ,th->ga, th->ga ,n);
  mpres_add (t3 ,t3 ,t4 ,n);
  mpres_ui_sub (t3 ,2 ,t3 ,n); // 2-be^2-ga^2
  mpres_mul (t4 ,P->X ,P->T ,n);
  mpres_mul (t5 ,P->Y ,P->Z ,n);
  mpres_add (t4 ,t4 ,t5 ,n);
  mpres_mul (t3 ,t3 ,t4 ,n); // (2-be^2-ga^2)*(XT+YZ)
  mpres_add (t2 ,t2 ,t3 ,n);
  //t2=(1-be*ga)*((be+ga)*(X+T)*(Y+Z)-(X^2+Y^2+Z^2+T^2)) + (2-be^2-ga^2)*(XT+YZ)
  mpres_add (t4 ,th->be ,th->ga ,n); // t4=be+ga
  mpres_mul (t2 ,t2 ,t4 ,n);
  //t2=(be+ga)*((1-be*ga)*((be+ga)*(X+T)*(Y+Z)-(X^2+Y^2+Z^2+T^2))+(2-be^2-ga^2)*(XT+YZ))

  mpres_add (t1 ,t1 ,t2 ,n);
  mpres_mul (divV->V0 ,divV->V0 ,t1 ,n);
  // We finished V0






  if ( mpres_is_zero (divU->u0,n) ) { // u0=0 
    // u0=0 => V0=0, v1v0=0, V1=f(-u1)/u1^2

    mpres_set_ui (divV->v1v0 ,0 ,n );

    mpres_add_ui (divV->V1 ,divU->u1 ,1 ,n);
    mpres_mul (divV->V1 ,divV->V1 ,divU->u1 ,n);
    mpres_add (t1 ,divU->u1 ,cHEll->la ,n);
    mpres_mul (divV->V1 ,divV->V1 ,t1 ,n);
    mpres_add (t1 ,divU->u1 ,cHEll->mu ,n);
    mpres_mul (divV->V1 ,divV->V1 ,t1 ,n);
    mpres_add (t1 ,divU->u1 ,cHEll->nu ,n);
    mpres_mul (divV->V1 ,divV->V1 ,t1 ,n);
    mpres_neg (divV->V1 ,divV->V1 ,n); // f(-u1)

    if (!mpres_invert (t1 , divU->u1, n)) // t1=1 / u1
      {
	mpres_gcd (f, divU->u0, n);
	
	mpres_clear (t1, n);
	mpres_clear (t2, n); 
	mpres_clear (t3, n);
	mpres_clear (t4, n);
	mpres_clear (t5, n);	  
	mpres_clear (T7p,n);
	mpres_clear (T9p,n);
	mpres_clear (T11p,n);
	mpres_clear (T12p,n);
	
	return MORPHISM_FAIL;
      }
    mpres_mul (divV->V1 ,divV->V1 ,t1 ,n);
    mpres_mul (divV->V1 ,divV->V1 ,t1 ,n); // V1=f(-u1)/u1^2


  }
  else { // u0 != 0


    if ( mpres_is_zero (divV->V0,n) ) { // V0=0
      // u0 != 0, V0 = 0  => v1v0=0
      /* 
	 V1:=-1/u0*( u0^2*Coefficient(f,4)-u0^2*u1-Coefficient(f,2)*u0+Coefficient(f,1)*u1));
      */
      
      if (!mpres_invert (t1 , divU->u0, n)) // t1=1 / u0
	{
	  mpres_gcd (f, divU->u0, n);
	  
	  mpres_clear (t1, n);
	  mpres_clear (t2, n); 
	  mpres_clear (t3, n);
	  mpres_clear (t4, n);
	  mpres_clear (t5, n);
	  
	  mpres_clear (T7p,n);
	  mpres_clear (T9p,n);
	  mpres_clear (T11p,n);
	  mpres_clear (T12p,n);
	  
	  return MORPHISM_FAIL;
	  
	}


      mpres_set_ui (divV->v1v0 ,0 ,n);


      mpres_mul (t2 ,cHEll->la ,cHEll->mu ,n);
      mpres_mul (t3 ,t2 ,cHEll->nu ,n); // t3 = la*mu*nu = Coefficient(f,1)
      mpres_add (t2 ,t2 ,t3 ,n);
      mpres_add (t5 ,cHEll->la ,cHEll->mu ,n);
      mpres_mul (t4 ,cHEll->nu ,t5 ,n);
      mpres_add (t2 ,t2 ,t4 ,n); 
      mpres_neg (t2 ,t2, n);
      // t2= - (la*mu*nu+la*mu+la*nu+mu*nu) = Coefficient(f,2)
      mpres_add(t4 ,t5 ,cHEll->nu ,n);
      mpres_add_ui (t4 ,t4 ,1 ,n);
      mpres_neg (t4 ,t4 ,n); // t4 = -(la+mu+nu+1) = Coefficient(f,4)

      mpres_mul (divV->V1 ,divU->u0 ,t4 ,n);
      mpres_mul (t5, divU->u0 ,divU->u1 ,n);
      mpres_sub (divV->V1 ,t5, divV->V1 ,n);
      mpres_add (divV->V1 ,divV->V1 ,t2 ,n);
      // -u0*Coefficient(f,4)+u0*u1+Coefficient(f,2)
      mpres_mul (t5 ,t3 ,divU->u1 ,n);
      mpres_mul (t3 ,t5 ,t1 ,n); // Coefficient(f,1)*u1/u0
      mpres_sub (divV->V1 ,divV->V1 ,t3 ,n);

    }
    else { // u0 !=0 , V0 !=0 

    /* 
      V1:=( u0^3+V0*u1 -u1^2*u0^2 -Coefficient(f,3)*u0^2 +Coefficient(f,1)*u0 + Coefficient(f,4)*u1*u0^2 )^2/( 4*V0*u0^2 );
      v1v0:=( u0^3+V0*u1 -u1^2*u0^2 -Coefficient(f,3)*u0^2 +Coefficient(f,1)*u0 + Coefficient(f,4)*u1*u0^2 )/( 2*u0 );
    */
      if (!mpres_invert (t1 , divU->u0, n)) // t1=1 / u0
	{
	  mpres_gcd (f, t1, n);
	  
	  mpres_clear (t1, n);
	  mpres_clear (t2, n); 
	  mpres_clear (t3, n);
	  mpres_clear (t4, n);
	  mpres_clear (t5, n);
	  
	  mpres_clear (T7p,n);
	  mpres_clear (T9p,n);
	  mpres_clear (T11p,n);
	  mpres_clear (T12p,n);
	  
	  return MORPHISM_FAIL;
	}


      mpres_mul (t3 ,cHEll->la ,cHEll->mu ,n);
      mpres_mul (t2 ,t3 ,cHEll->nu ,n); // t2 = la*mu*nu = Coefficient(f,1)
      mpres_add (t5 ,cHEll->la ,cHEll->mu ,n);
      mpres_mul (t4 ,cHEll->nu ,t5 ,n);
      mpres_add (t3 ,t3 ,t4 ,n); 
      mpres_add (t4 ,t5 ,cHEll->nu ,n);
      mpres_add (t3 ,t3 ,t4 ,n);
      // t3 =la*mu+la*nu+mu*nu + la+mu+nu = Coefficient(f,3)
      mpres_add_ui (t4 ,t4 ,1 ,n);
      mpres_neg (t4 ,t4 ,n); // t4 = -(la+mu+nu+1) = Coefficient(f,4)


      mpres_sub (divV->v1v0 ,t4 ,divU->u1 ,n);
      mpres_mul (divV->v1v0 ,divV->v1v0 ,divU->u1 ,n);//u1*(-u1+Coef(f,4))
      mpres_sub (divV->v1v0 ,divV->v1v0 ,t3 ,n);
      mpres_add (divV->v1v0 ,divV->v1v0 ,divU->u0 ,n);
      // u0 - u1^2 - Coefficient(f,3) + u1*Coefficient(f,4)
      mpres_mul (divV->v1v0 ,divV->v1v0 ,divU->u0 ,n);
      mpres_add (divV->v1v0 ,divV->v1v0 ,t2 ,n);
      //u0^2-u0*u1^2-u0*Coefficient(f,3)+u0*u1*Coefficient(f,4)+Coefficient(f,1)

      mpres_mul (t2 ,divV->V0 ,t1 ,n);
      mpres_mul (t1 ,divU->u1 ,t2 ,n);
      mpres_add (divV->v1v0 ,divV->v1v0 ,t1 ,n);
      
      mpres_set_ui (t1 ,2 ,n);
      mpres_invert (t2 ,t1, n); // t2 = 1/2
      mpres_mul (divV->v1v0 ,divV->v1v0 ,t2 ,n);


      if (!mpres_invert (t1 , divV->V0, n)) // t1=1 / V0
	{
	  mpres_gcd (f, divV->V0, n);
	  
	  mpres_clear (t1, n);
	  mpres_clear (t2, n); 
	  mpres_clear (t3, n);
	  mpres_clear (t4, n);
	  mpres_clear (t5, n);
	  
	  mpres_clear (T7p,n);
	  mpres_clear (T9p,n);
	  mpres_clear (T11p,n);
	  mpres_clear (T12p,n);
	  
	  return MORPHISM_FAIL;
	}
      mpres_mul (divV->V1 ,divV->v1v0 ,divV->v1v0 ,n);
      mpres_mul (divV->V1 ,divV->V1 ,t1 ,n);


    }
  }


  mpres_clear (t1, n);
  mpres_clear (t2, n); 
  mpres_clear (t3, n);
  mpres_clear (t4, n);
  mpres_clear (t5, n);

  mpres_clear (T7p, n);
  mpres_clear (T9p, n);
  mpres_clear (T11p, n);
  mpres_clear (T12p, n);


  return MORPHISM;

}











// **********************************************************************












/*
  Input:  The Mumford polynomial u of degree 1
          The square of the Mumford polynomial v // V=v^2
  Output: The points on the two underlying elliptic curves 
          These curves are in short Weierstrass form and we get A

  On the hyperelliptic curve we have the point (x::1) with x=-u0/u1
*/
int jac_to_EllW (mpz_t f, mpmod_t n,
		 curve *T1, curve *T2,
		 DivMumfordU divU, DivMumfordV divV,
		 curveHyperEll cHEll) {
  int test;

  if (divU->degree == 0) { // u=0
    return MORPHISM_FOUND_ZERO_CURVE_1_AND_2;
  }
  else if (divU->degree == 1) { // degre de u =1
    test = jac_to_EllW_Degree1 (f,n,T1,T2,divU,cHEll);
  }
  else  { // degre de u =2
    test = jac_to_EllW_Degree2 (f,n,T1,T2,divU,divV,cHEll);
  }

  return test;
}





/*
  Input:  The Mumford polynomial u of degree 1
          The square of the Mumford polynomial v // V=v^2
  Output: The points on the two underlying elliptic curves 
          These curves are in short Weierstrass form and we get A

  On the hyperelliptic curve we have the point (x::1) with x=-u0/u1
*/
int jac_to_EllW_Degree1 (mpz_t f, mpmod_t n,
			 curve *T1, curve *T2,
			 DivMumfordU divU,
			 curveHyperEll cHEll) {
  int test;
  mpres_t A2,a6;
  mpres_t x,z;

  mpres_init (A2 ,n);
  mpres_init (a6 ,n);
  mpres_init (x ,n);
  mpres_init (z ,n);

  mpres_neg (x ,divU->u0 ,n); 
  mpz_set (z ,divU->u1);   

  // for the first curve
  test = coeff_EllW (f,n,T1,A2,a6,cHEll );
  if (test ==  MORPHISM_FAIL ){
    mpres_clear (A2 ,n);
    mpres_clear (a6 ,n);
    mpres_clear (x ,n);
    mpres_clear (z ,n);
    return MORPHISM_FAIL;
  }
  test = HEll_EllW_degree1 (f,n ,x,z, T1, A2,a6 ,cHEll);
  if (test ==  MORPHISM_FAIL ){
    mpres_clear (A2 ,n);
    mpres_clear (a6 ,n);
    mpres_clear (x ,n);
    mpres_clear (z ,n);
    return test;
  }
  if (test ==  MORPHISM_FOUND_ZERO_CURVE ){
    mpres_clear (A2 ,n);
    mpres_clear (a6 ,n);
    mpres_clear (x ,n);
    mpres_clear (z ,n);
    return MORPHISM_FOUND_ZERO_CURVE_1;
  }





  // for the second curve
  mpres_neg (cHEll->q ,cHEll->q ,n);
  test = coeff_EllW (f,n,T2,A2,a6,cHEll );
  if (test ==  MORPHISM_FAIL ){
    mpres_clear (A2 ,n);
    mpres_clear (a6 ,n);
    mpres_clear (x ,n);
    mpres_clear (z ,n);
    return MORPHISM_FAIL;
  }
  test = HEll_EllW_degree1 (f,n ,x,z, T1, A2,a6 ,cHEll);
  if (test !=  MORPHISM ){
    mpres_clear (A2 ,n);
    mpres_clear (a6 ,n);
    mpres_clear (x ,n);
    mpres_clear (z ,n);
    return test;
  }
  if (test ==  MORPHISM_FOUND_ZERO_CURVE ){
    mpres_clear (A2 ,n);
    mpres_clear (a6 ,n);
    mpres_clear (x ,n);
    mpres_clear (z ,n);
    return MORPHISM_FOUND_ZERO_CURVE_1;
  }


  mpres_clear (A2 ,n);
  mpres_clear (a6 ,n);
  mpres_clear (x ,n);
  mpres_clear (z ,n);
  return MORPHISM;
}







/*
  special case of jac_to_EllW_degree2: when one of y1 ,y2 is zero
    coeff[0]=0 <=> al^2-be^2*delta <=> y1 or y2 is zero
  In this case x1=0,1,la,mu or nu
  First we find which one it is. Then it is possible to send the other point on
    the elliptic curves.
  Theoricaly we should obtain phi((x1,0))+phi((x2,y2)) on the curve with
    phi((x1,0)) zero or of 2-torsion.
  We compute 2*phi((x2,y2)) and obtain a multiple of the point
 */
int degree2_case_y1_equal_zero (mpz_t f,mpmod_t n,curve *T1, curve *T2,DivMumfordU divU, DivMumfordV divV, curveHyperEll cHEll) {
  
  int test;
  mpres_t t1,t2;
  mpres_t g;
  mpres_t x,z;
  mpres_t A2,a6;

  mpres_init (t1,n);
  mpres_init (t2,n); 
  mpres_init (g,n);
  mpres_init (x,n);

  mpres_gcd (f ,divU->u0 ,n);
  mpres_set_z (g ,f ,n);
  if ( mpres_is_zero (g,n) ) { // x1 = 0
    mpres_neg (x ,divU->u1 ,n); // x2=-u1-x1
  }
  else if ( mpz_cmp_ui (f,1) != 0 ) { // f !=0,1,n so it has a factor of n
    mpres_clear (x ,n);
    mpres_clear (g ,n);
    mpres_clear (t1,n);
    mpres_clear (t2,n);
    return MORPHISM_FAIL;
  }
  else {
    mpres_add_ui (t1 ,divU->u1 ,1 ,n);
    mpres_add (t1 ,t1 ,divU->u0 ,n); // 1^2 + u1*1 + u0
    
    mpres_gcd (f ,t1 ,n);
    mpres_set_z (g ,f ,n);
    if ( mpres_is_zero (g,n) ) { // x1 = 1
      mpres_add_ui (x ,divU->u1 ,1 ,n);
      mpres_neg (x ,x ,n); // x2=-u1-x1
    }
    else if ( mpz_cmp_ui (f,1) != 0 ) { // f !=0,1,n so it has a factor of n
      mpres_clear (g ,n);
      mpres_clear (x ,n);
      mpres_clear (t1,n);
      mpres_clear (t2,n);
      return MORPHISM_FAIL;
    }
    else {
      mpres_mul (t1 ,divU->u1 ,cHEll->la ,n);
      mpres_add (t1 ,t1 ,divU->u0 ,n);
      mpres_mul (t2 ,cHEll->la ,cHEll->la ,n);
      mpres_add (t1 ,t1 ,t2 ,n); // la^2 + u1*la + u0
      
      mpres_gcd (f ,t1 ,n);
      mpres_set_z (g ,f ,n);
      if ( mpres_is_zero (g,n) ) { // x1 = la
	mpres_add (x ,divU->u1 ,cHEll->la ,n);
	mpres_neg (x ,x ,n); // x2=-u1-x1
      }
      else if ( mpz_cmp_ui (f,1) != 0 ) {
	// f !=0,1,n so it has a factor of n
	mpres_clear (g ,n);
	mpres_clear (x ,n);
	mpres_clear (t1,n);
	mpres_clear (t2,n);
	return MORPHISM_FAIL;
      }
      else {
	mpres_mul (t1 ,divU->u1 ,cHEll->mu ,n);
	mpres_add (t1 ,t1 ,divU->u0 ,n);
	mpres_mul (t2 ,cHEll->mu ,cHEll->mu ,n);
	mpres_add (t1 ,t1 ,t2 ,n); // mu^2 + u1*mu + u0
	
	mpres_gcd (f ,t1 ,n);
	mpres_set_z (g ,f ,n);
	if ( mpres_is_zero (g,n) ) { // x1 = mu
	  mpres_add (x ,divU->u1 ,cHEll->mu ,n);
	  mpres_neg (x ,x ,n); // x2=-u1-x1
	}
	else if ( mpz_cmp_ui (f,1) != 0 ) {
	  // f !=0,1,n so it has a factor of n
	  mpres_clear (g ,n);
	  mpres_clear (x ,n);
	  mpres_clear (t1,n);
	  mpres_clear (t2,n);
	  return MORPHISM_FAIL;
	}
	else {
	  mpres_mul (t1 ,divU->u1 ,cHEll->nu ,n);
	  mpres_add (t1 ,t1 ,divU->u0 ,n);
	  mpres_mul (t2 ,cHEll->nu ,cHEll->nu ,n);
	  mpres_add (t1 ,t1 ,t2 ,n); // nu^2 + u1*nu + u0
	  
	  mpres_gcd (f ,t1 ,n);
	  mpres_set_z (g ,f ,n);
	  if ( mpres_is_zero (g,n) ) { // x1 = nu
	    mpres_add (x ,divU->u1 ,cHEll->nu ,n);
	    mpres_neg (x ,x ,n); // x2=-u1-x1
	  }
	  else if ( mpz_cmp_ui (f,1) != 0 ) {
	    // f !=0,1,n so it has a factor of n
	    mpres_clear (g ,n);
	    mpres_clear (x ,n);
	    mpres_clear (t1,n);
	    mpres_clear (t2,n);
	    return MORPHISM_FAIL;
	  }
	  else {
	    mpres_clear (g ,n);
	    mpres_clear (x ,n);
	    mpres_clear (t1,n);
	    mpres_clear (t2,n);
	    return MORPHISM_ERROR;
	  }
	}
      }
    }
  } 
  mpres_clear (g ,n);
  mpres_clear (t1,n);
  mpres_clear (t2,n);

  mpres_init( A2 ,n);
  mpres_init( a6 ,n);
  mpres_init( z ,n);

  // we have set x2 so we send it to the elliptic curve
  mpres_set_ui (z ,1 ,n);
  
  // for the first curve
  test = coeff_EllW (f,n,T1,A2,a6,cHEll );
  if (test ==  MORPHISM_FAIL ){
    mpres_clear (A2 ,n);
    mpres_clear (a6 ,n);
    mpres_clear (x ,n);
    mpres_clear (z ,n);
    return MORPHISM_FAIL;
  }
  
  test = HEll_EllW_degree1 (f,n ,x,z, T1, A2,a6 ,cHEll);
  if (test ==  MORPHISM_FAIL ){
    mpres_clear (A2 ,n);
    mpres_clear (a6 ,n);
    mpres_clear (x ,n);
    mpres_clear (z ,n);
    return MORPHISM_FAIL;
  }
  if (test ==  MORPHISM_FOUND_ZERO_CURVE ){
    mpres_clear (A2 ,n);
    mpres_clear (a6 ,n);
    mpres_clear (x ,n);
    mpres_clear (z ,n);
    return MORPHISM_FOUND_ZERO_CURVE_1;
  }
  if (test ==  MORPHISM_STEP1 ){
    mpres_clear (A2 ,n);
    mpres_clear (a6 ,n);
    mpres_clear (x ,n);
    mpres_clear (z ,n);
    return MORPHISM_STEP1;
  } 
  
  
  // for the second curve
  mpres_neg (cHEll->q ,cHEll->q ,n);
  test = coeff_EllW (f,n,T2,A2,a6,cHEll );
  if (test ==  MORPHISM_FAIL ){
    mpres_clear (A2 ,n);
    mpres_clear (a6 ,n);
    mpres_clear (x ,n);
    mpres_clear (z ,n);
    return MORPHISM_FAIL;
  }
  
  test = HEll_EllW_degree1 (f,n ,x,z, T1, A2,a6 ,cHEll);
  if (test !=  MORPHISM ){
    mpres_clear (A2 ,n);
    mpres_clear (a6 ,n);
    mpres_clear (x ,n);
    mpres_clear (z ,n);
    return test;
  }
  if (test ==  MORPHISM_FOUND_ZERO_CURVE ){
    mpres_clear (A2 ,n);
    mpres_clear (a6 ,n);
    mpres_clear (x ,n);
    mpres_clear (z ,n);
    return MORPHISM_FOUND_ZERO_CURVE_1;
  }
  if (test ==  MORPHISM_STEP1 ){
    mpres_clear (A2 ,n);
    mpres_clear (a6 ,n);
    mpres_clear (x ,n);
    mpres_clear (z ,n);
    return MORPHISM_STEP1;
  } 
  
  mpres_clear (A2 ,n);
  mpres_clear (a6 ,n);
  mpres_clear (x ,n);
  mpres_clear (z ,n);
  
      
  // double the point on the two curve
  if ( double_short_weierstrass (f,n,T1 )==0 ) {
    return MORPHISM_STEP1;
  }
  if ( double_short_weierstrass (f,n,T2 )==0 ) {
    return MORPHISM_STEP1;
  }



  
  return MORPHISM;
 }      


    
  














/*
  Input:  The Mumford polynomial u of degree 2 
          The square of the Mumford polynomial v // V=v^2
  Output: The points on the two underlying elliptic curves 
          These curves are in short Weierstrass form and we get A

	  
*/
int jac_to_EllW_Degree2 (mpz_t f, mpmod_t n,
			 curve *T1, curve *T2,
			 DivMumfordU divU, DivMumfordV divV,
			 curveHyperEll cHEll) {
  int test;
  mpalgpol_t pol;
  mpalgres_t x1,x2,y1,y2,z1,z2;
  mpalgres_t X1,X2,Y1,Y2,Z1,Z2;
  mpres_t al,be;
  mpres_t t1,t2;
  mpres_t A2,a6;
  mpalgres_t a4twist,R,B,invB;

  mpres_init (al,n);
  mpres_init (be,n);
  mpres_init (t1,n);
  mpres_init (t2,n);

  mpalgpol_init (pol,n);
  mpalgres_init (x1,n);
  mpalgres_init (x2,n);
  mpalgres_init (y1,n);
  mpalgres_init (y2,n);
  mpalgres_init (z1,n);
  mpalgres_init (z2,n);
  mpalgres_init (B,n);
  mpalgres_init (invB,n);



  mpres_set_ui (pol->coeff[1] ,0,n);
  mpres_mul_ui (divV->v1v0 ,divV->v1v0 ,4 ,n); //since we only use 4*v1v0

  // Construction of delta such that x1=-u1+sqrt(delta)
  mpres_mul (be ,divU->u1 ,divU->u1 ,n); // be=u1^2
  mpres_mul_ui (t2 ,divU->u0 ,4 ,n);
  mpres_sub (t2 ,be ,t2 ,n); // t2 = delta = u1^2-4*u0
  mpres_neg (pol->coeff[0] ,t2 ,n); // -delta



  mpres_neg (x1[0] ,divU->u1 ,n);
  mpres_neg (x2[0] ,divU->u1 ,n);
  mpres_set_ui (x1[1] ,1 ,n);  // x1 = -u1+sqrt(delta)
  mpres_set_si (x2[1] ,-1 ,n); // x2 = -u1-sqrt(delta)

  mpalgres_set_ui (z1 ,2 ,n); // z1 = 2
  mpalgres_set_ui (z2 ,2 ,n); // z2 = 2



  mpalgres_set_ui (y1 ,1 ,n); // y1' = 1



  // we have  B = y1^2 = al + be* sqrt(delta)
  //         al = 4V0 -4*v1v0*u1 + V1*(delta+u1^2)
  //         be = 4*v1v0 - 2u1*V1
  mpres_add (al ,t2 ,be ,n);
  mpres_mul (al ,al ,divV->V1 ,n); // V1*(delta+u1^2)
  mpres_mul (t1 ,divV->v1v0 ,divU->u1 ,n); // 4*v1v0*u1
  mpres_mul_ui (be ,divV->V0 ,4 ,n);  // 4V0
  mpres_sub (t1 ,be ,t1 ,n);// t1 = 4V0 -4*v1v0*u1
  mpres_add (al ,al ,t1 ,n );  // al = 4V0 -4*v1v0*u1 + V1*(delta+u1^2)

  mpres_mul (be ,divV->V1 ,divU->u1 ,n);
  mpres_mul_ui (be ,be ,2 ,n);
  mpres_sub (be ,divV->v1v0 ,be ,n); // be = 4*v1v0 - 2u1*V1

  mpz_set (B[0] ,al);
  mpz_set (B[1] ,be);


  if (mpalgres_is_zero (B ,pol,n) ) { // y1=0
    mpres_clear (al,n);
    mpres_clear (be,n);
    mpres_clear (t1,n);
    mpres_clear (t2,n);
    mpalgpol_clear (pol,n);
    mpalgres_clear (x1,n);
    mpalgres_clear (x2,n);
    mpalgres_clear (y1,n);
    mpalgres_clear (y2,n);
    mpalgres_clear (z1,n);
    mpalgres_clear (z2,n);

    mpalgres_clear (B,n);
    mpalgres_clear (invB,n);

    test = degree2_case_y1_equal_zero (f,n, T1,T2, divU,divV,cHEll);
    return test;
  }

  test = mpalgres_invert (invB ,B ,pol,n,f); // 1/B=1/y1^2
  if ( test == -1 ) {
    // This should not happen
    mpres_clear (al,n);
    mpres_clear (be,n);
    mpres_clear (t1,n);
    mpres_clear (t2,n);
    mpalgpol_clear (pol,n);
    mpalgres_clear (x1,n);
    mpalgres_clear (x2,n);
    mpalgres_clear (y1,n);
    mpalgres_clear (y2,n);
    mpalgres_clear (z1,n);
    mpalgres_clear (z2,n);

    mpalgres_clear (B,n);
    mpalgres_clear (invB,n);
    return MORPHISM_ERROR;
  }
  if ( test == 0 ) {
    mpres_clear (al,n);
    mpres_clear (be,n);
    mpres_clear (t1,n);
    mpres_clear (t2,n);
    mpalgpol_clear (pol,n);
    mpalgres_clear (x1,n);
    mpalgres_clear (x2,n);
    mpalgres_clear (y1,n);
    mpalgres_clear (y2,n);
    mpalgres_clear (z1,n);
    mpalgres_clear (z2,n);

    mpalgres_clear (B,n);
    mpalgres_clear (invB,n);
    return MORPHISM_FAIL;
  }
  // ok we have inverted B



  // y2' = y2/y1 = y2y1/y1^2 = y2y1*invB
  // y1*y2 = 4*u0*V1 - 4*v0v1*u1 + 4V0


  mpres_mul (t2 ,divU->u0 ,divV->V1 ,n);
  mpres_mul_ui(t2 ,t2 ,4 ,n);
  mpres_add (t1 ,t1 ,t2 ,n); // 4*u0*V1 - 4*v0v1*u1 + 4V0
  mpalgres_mul_mpres (y2 ,invB ,t1 ,pol ,n);



  mpres_clear (al,n);
  mpres_clear (be,n);
  mpres_clear (t1,n);
  mpres_clear (t2,n);




  // we have construct the points on the hyperelliptic curve in the good algebra
  // now go to the elliptic curves


  mpalgres_init (R,n);
  mpalgres_init (a4twist,n);
  mpres_init (A2,n);
  mpres_init (a6,n);
 
  mpalgres_init (X1,n);
  mpalgres_init (X2,n);
  mpalgres_init (Y1,n);
  mpalgres_init (Y2,n);
  mpalgres_init (Z1,n);
  mpalgres_init (Z2,n);



  // the first one:
  test= coeff_EllW_EllWtwist (f,n ,pol ,T1,a6,A2  ,R,a4twist  ,B,invB,cHEll);
  if (test == MORPHISM_FAIL) {
    mpalgres_clear (R,n);
    mpalgres_clear (a4twist,n);
    mpres_clear (a6,n);
    mpres_clear (A2,n);

    mpalgpol_clear (pol,n);
    mpalgres_clear (x1,n);
    mpalgres_clear (x2,n);
    mpalgres_clear (y1,n);
    mpalgres_clear (y2,n);
    mpalgres_clear (z1,n);
    mpalgres_clear (z2,n);

    mpalgres_clear (X1,n);
    mpalgres_clear (X2,n);
    mpalgres_clear (Y1,n);
    mpalgres_clear (Y2,n);
    mpalgres_clear (Z1,n);
    mpalgres_clear (Z2,n);

    mpalgres_clear (B,n);
    mpalgres_clear (invB,n);
    return MORPHISM_FAIL;
  }

  HEll_to_EllW(n,pol, X1,Y1,Z1, x1,y1,z1, cHEll );
  HEll_to_EllW(n,pol, X2,Y2,Z2, x2,y2,z2, cHEll ); 
  EllW_to_EllWshort(n,pol, R,A2, X1,Y1,Z1, X1,Y1,Z1); 
  EllW_to_EllWshort(n,pol, R,A2, X2,Y2,Z2, X2,Y2,Z2);
  test= addW_and_scale (f,n,pol,T1,a6,X1,Y1,Z1,X2,Y2,Z2,a4twist,R);

  if ( test != MORPHISM ) {
    mpalgres_clear (R,n);
    mpalgres_clear (a4twist,n);
    mpres_clear (a6,n);
    mpres_clear (A2,n);

    mpalgpol_clear (pol,n);
    mpalgres_clear (x1,n);
    mpalgres_clear (x2,n);
    mpalgres_clear (y1,n);
    mpalgres_clear (y2,n);
    mpalgres_clear (z1,n);
    mpalgres_clear (z2,n);

    mpalgres_clear (X1,n);
    mpalgres_clear (X2,n);
    mpalgres_clear (Y1,n);
    mpalgres_clear (Y2,n);
    mpalgres_clear (Z1,n);
    mpalgres_clear (Z2,n);

    mpalgres_clear (B,n);
    mpalgres_clear (invB,n);

    if (test == MORPHISM_FOUND_ZERO_CURVE) {
      return MORPHISM_FOUND_ZERO_CURVE_1;
    }
    else {
      return test;
    }
  }





  // the second one:
  mpres_neg (cHEll->q ,cHEll->q ,n);
  test= coeff_EllW_EllWtwist (f,n ,pol ,T2,a6,A2  ,R,a4twist  ,B,invB,cHEll);
  if (test == MORPHISM_FAIL) {
    mpalgres_clear (R,n);
    mpalgres_clear (a4twist,n);
    mpres_clear (a6,n);
    mpres_clear (A2,n);

    mpalgpol_clear (pol,n);
    mpalgres_clear (x1,n);
    mpalgres_clear (x2,n);
    mpalgres_clear (y1,n);
    mpalgres_clear (y2,n);
    mpalgres_clear (z1,n);
    mpalgres_clear (z2,n);

    mpalgres_clear (X1,n);
    mpalgres_clear (X2,n);
    mpalgres_clear (Y1,n);
    mpalgres_clear (Y2,n);
    mpalgres_clear (Z1,n);
    mpalgres_clear (Z2,n);

    mpalgres_clear (B,n);
    mpalgres_clear (invB,n);
    return MORPHISM_FAIL;
  }

  HEll_to_EllW(n,pol, X1,Y1,Z1, x1,y1,z1, cHEll );
  HEll_to_EllW(n,pol, X2,Y2,Z2, x2,y2,z2, cHEll ); 
  EllW_to_EllWshort(n,pol, R,A2, X1,Y1,Z1, X1,Y1,Z1); 
  EllW_to_EllWshort(n,pol, R,A2, X2,Y2,Z2, X2,Y2,Z2);
  test= addW_and_scale (f,n,pol,T2,a6,X1,Y1,Z1,X2,Y2,Z2,a4twist,R);

  if ( test != MORPHISM ) {
    mpalgres_clear (R,n);
    mpalgres_clear (a4twist,n);
    mpres_clear (a6,n);
    mpres_clear (A2,n);

    mpalgpol_clear (pol,n);
    mpalgres_clear (x1,n);
    mpalgres_clear (x2,n);
    mpalgres_clear (y1,n);
    mpalgres_clear (y2,n);
    mpalgres_clear (z1,n);
    mpalgres_clear (z2,n);

    mpalgres_clear (X1,n);
    mpalgres_clear (X2,n);
    mpalgres_clear (Y1,n);
    mpalgres_clear (Y2,n);
    mpalgres_clear (Z1,n);
    mpalgres_clear (Z2,n);

    mpalgres_clear (B,n);
    mpalgres_clear (invB,n);

    if (test == MORPHISM_FOUND_ZERO_CURVE) {
      return MORPHISM_FOUND_ZERO_CURVE_1;
    }
    else {
      return test;
    }
  }



  mpalgres_clear (R,n);
  mpalgres_clear (a4twist,n);
  mpres_clear (a6,n);
  mpres_clear (A2,n);
  
  mpalgpol_clear (pol,n);
  mpalgres_clear (x1,n);
  mpalgres_clear (x2,n);
  mpalgres_clear (y1,n);
  mpalgres_clear (y2,n);
  mpalgres_clear (z1,n);
  mpalgres_clear (z2,n);
  
  mpalgres_clear (X1,n);
  mpalgres_clear (X2,n);
  mpalgres_clear (Y1,n);
  mpalgres_clear (Y2,n);
  mpalgres_clear (Z1,n);
  mpalgres_clear (Z2,n);
  
  mpalgres_clear (B,n);
  mpalgres_clear (invB,n);
  return MORPHISM;
}









/*
  send a point (x,y,z) on the hyperelliptic curve to the underlying curve
  defined by cHEll->q (note that the other is defined by -q)
  Work with coordinates in k[Y]/pol(Y)
  return the result in (x,y,z)

  (x,y,z) --> ( wd*( x-(mu+q)z )^2*( x-(mu-q)z ) ,wn*y*z^2, wd*(x-(mu-q)z)^3  )
    wn = 8*q   wd = (mu-q)*( -1+(mu-q) )
*/ 
void HEll_to_EllW(mpmod_t n,mpalgpol_t pol,
		  mpalgres_t aX2,mpalgres_t aY2,mpalgres_t aZ2,
		  mpalgres_t aX ,mpalgres_t aY ,mpalgres_t aZ ,
		  curveHyperEll cHEll ) {
  mpalgres_t aT1,aT2;
  mpalgres_t x,y,z;
  mpres_t temp1,temp2;

  mpres_init (temp1 ,n);
  mpres_init (temp2 ,n);
  mpalgres_init (aT1 ,n);
  mpalgres_init (aT2 ,n);
  mpalgres_init (x,n);
  mpalgres_init (y,n);
  mpalgres_init (z,n);

  mpalgres_set (x,aX,n);
  mpalgres_set (y,aY,n);
  mpalgres_set (z,aZ,n);


  mpres_mul_ui (temp1 ,cHEll->q ,8 ,n);
  mpalgres_mul (aT1 ,z ,z ,pol ,n);
  mpalgres_mul (y ,y ,aT1 ,pol ,n);
  mpalgres_mul_mpres (y ,y ,temp1 ,pol ,n); // 8*q  *  y * z^2

  mpres_add (temp1 ,cHEll->mu ,cHEll->q ,n);
  mpalgres_mul_mpres (aT2 ,z ,temp1 ,pol ,n); // (mu+q)z
  mpres_sub (temp1 ,cHEll->mu ,cHEll->q ,n);//(mu-q)
  mpalgres_mul_mpres (aT1 ,z ,temp1 ,pol ,n); // (mu-q)z

  mpalgres_sub (z ,x ,aT1 ,pol ,n); // x-(mu-q)z
  mpalgres_sub (x ,x ,aT2 ,pol ,n); // x-(mu+q)z




  mpalgres_mul (x ,x ,x ,pol ,n);
  mpalgres_mul (x ,x ,z ,pol ,n); // ( x-(mu+q)z )^2*( x-(mu-q)z )


  mpalgres_mul (aT1 ,z ,z ,pol ,n);
  mpalgres_mul (z ,z ,aT1 ,pol ,n); // ( x-(mu-q)z )^3

  mpres_sub_ui (temp2 ,temp1 ,1 ,n);
  mpres_mul (temp1 ,temp1 ,temp2 ,n); // wd = (mu-q)*( -1+(mu-q) )

  mpalgres_mul_mpres (x ,x ,temp1 ,pol ,n);
  mpalgres_mul_mpres (z ,z ,temp1 ,pol ,n);



  mpalgres_set (aX2,x,n);
  mpalgres_set (aY2,y,n);
  mpalgres_set (aZ2,z,n);


  mpres_clear (temp1 ,n);
  mpres_clear (temp2 ,n);
  mpalgres_clear (aT1 ,n);
  mpalgres_clear (aT2 ,n);
  mpalgres_clear (x,n);
  mpalgres_clear (y,n);
  mpalgres_clear (z,n);
}








/* 
  Get many coefficient of the short weierstrass curve

  the long weierstrass curve is given by
    B*rc* y^2*z = ( x - z ) * ( x - x2^2 *z ) * ( x - x3^2 *z )
              = x^3  +  A2 * x^2  +  A4 * x  +  A6
       rc = -q*mu*(mu-1)
       x2 = (mu+q) / (mu-q)      x3 = (1-(mu+q)) / (1-(mu-q))
  first get the curve
    R * y^2*z = x^3  +  a4' * x * z^2  +  a6' * z^3
      by (x,y,z) -> (x+A2/3*z,y,z) 
      a4 = A4 - A2^2/3    a6 = A6 - A2*A4/3 + 2*A2^3/27
      R=B*rc 
  get the curve
    y^2*z = x^3  +  a4' * x * z^2  +  a6' * z^3
      by (x,y,z) -> (x,y,R*z) 
      a4' = a4/R^2   a6' = a6/R^3

  It get the coefficient a4 and a6 of the short weierstrass curve
                           a4 is put in T->A
			 A2 of the long weierstrass form
			 a4twist of the short weierstrass twisted curve
			 R
*/
int coeff_EllW_EllWtwist(mpz_t f,mpmod_t n,mpalgpol_t pol,
			 curve *T,mpres_t a6,mpres_t A2,
			 mpalgres_t R, mpalgres_t a4twist,
			 mpalgres_t B,mpalgres_t invB,curveHyperEll cHEll ) {

  mpres_t t1,t2,t3,t4;
  mpalgres_t aTemp;

  mpalgres_set_zero (a4twist ,n);
  mpalgres_set_zero (R ,n);


  mpres_init (t1,n);
  mpres_init (t2,n);
  mpres_init (t3,n);
  mpres_init (t4,n);

  mpres_set_ui (R[1] ,0 ,n);

  mpres_sub (t1, cHEll->mu, cHEll->q ,n); // t1 = mu-q 
  mpres_ui_sub (t2 ,1 ,t1 ,n); // t2 = 1-(mu-q)

  mpres_mul (t3 ,t1 ,t2 ,n);
  if ( !mpres_invert(t3,t3,n) ) { // t3 = 1 / ( (mu-q)*(1-(mu-q)) )
    mpres_gcd(f ,t3 ,n);
    mpres_clear (t1,n);
    mpres_clear (t2,n);
    mpres_clear (t3,n);
    mpres_clear (t4,n);
    return MORPHISM_FAIL;
  }


  mpres_add (t4, cHEll->mu, cHEll->q ,n); // t4 = mu+q 

  mpres_mul (t2 ,t3 ,t2 ,n); //  1 / (mu-q)
  mpres_mul (t2 ,t2 ,t4 ,n); // t2 = (mu+q) / (mu-q) = x2

  mpres_mul (t3 ,t3 ,t1 ,n); // 1 / (1-(mu-q))
  mpres_ui_sub (t4 ,1 ,t4 ,n);
  mpres_mul (t3 ,t3 ,t4 ,n); // t3 = (1-(mu+q)) / (1-(mu-q)) = x3


  mpres_mul (t2 ,t2 ,t2 ,n); // t2 = x2^2
  mpres_mul (t3 ,t3 ,t3 ,n); // t3 = x3^2


  mpres_mul (a6 ,t2 ,t3 ,n);   // x2^2*x3^2
  mpres_add (T->A ,t2 ,t3 ,n); // x2^2+x3^2




  mpres_add_ui (A2 ,T->A ,1 ,n);
  mpres_neg (A2 ,A2 ,n);         //      A2 = - (1+x2^2+x3^2)
  mpres_add (T->A ,T->A ,a6 ,n); //      A4 = x2^2*x3^2 + x2^2 + x3^2
  mpres_neg (a6 ,a6 ,n);         // a6 = A6 = - x2^2*x3^2



  mpres_mul_ui (a6 ,a6 ,27 ,n);  // 27*A6
  mpres_mul (t1 ,A2 ,T->A ,n);
  mpres_mul_ui (t1 ,t1 ,9 ,n);   // 9*A2*A4
  mpres_sub (a6 ,a6 ,t1 ,n);     // 27*A6 - 9*A2*A4
  mpres_mul (t1 ,A2 ,A2 ,n);     // A2^2
  mpres_mul (t3 ,t1 ,A2 ,n);
  mpres_mul_ui (t3 ,t3 ,2 ,n);
  mpres_add (a6 ,a6 ,t3 ,n);     // 27*A6 - 9*A2*A4 + 2*A2^3

  mpres_mul_ui (T->A ,T->A ,3 ,n);
  mpres_sub (T->A ,T->A ,t1 ,n); // 3*A4 - A2^2





  mpres_sub_ui (t1 ,cHEll->mu ,1 ,n);
  mpres_mul (t1 ,t1 ,cHEll->mu ,n);
  mpres_mul (t1 ,t1 ,cHEll->q ,n);
  mpres_neg (t1 ,t1 ,n); // t1 = rc = -q*mu*(mu-1)


  mpalgres_mul_mpres (R ,B ,t1 ,pol ,n); // R  = rc*B
  mpres_mul_ui (t2 ,t1 ,3 ,n); //  3*rc


  if ( !mpres_invert(t2,t2,n) ) { // t2 = 1 / ( 3*rc )
    mpres_gcd(f ,t2 ,n);
    mpres_clear (t1,n);
    mpres_clear (t2,n);
    mpres_clear (t3,n);
    mpres_clear (t4,n);
    return MORPHISM_FAIL; 
  }



  mpalgres_init (aTemp ,n);



  mpalgres_mul_mpres (aTemp ,invB ,t2 ,pol ,n); //  1 / ( 3*rc*B ) = 1/3R
  mpalgres_mul (aTemp ,aTemp ,aTemp ,pol ,n);   
  mpalgres_mul_ui (aTemp ,aTemp ,3 ,pol ,n);    // 1 / ( 3*R^2 )



  mpalgres_mul_mpres (a4twist ,aTemp ,T->A ,pol ,n);
  // a4' = (3*A4 - A2^2) / (3*R^2)

  mpres_mul (t1 ,t1 ,t2, n); // 1/3
  mpres_mul (T->A ,T->A ,t1 ,n); // a4 = (3*A4 - A2^2) / 3
  mpres_mul (t2 ,t1 ,t1, n); 
  mpres_mul (t1 ,t1 ,t2, n); // 1/27
  mpres_mul (a6 ,a6 ,t1 ,n); // a6 = ( 27*A6 - 9*A2*A4 + 2*A2^3 ) / 27



  mpalgres_clear (aTemp ,n);
  mpres_clear (t1,n);
  mpres_clear (t2,n);
  mpres_clear (t3,n);
  mpres_clear (t4,n);
  return MORPHISM;
}





/* 
  Get many coefficient of the short weierstrass curve

  the long weierstrass curve is given by
    rc* y^2*z = ( x - z ) * ( x - x2^2 *z ) * ( x - x3^2 *z )
              = x^3  +  A2 * x^2  +  A4 * x  +  A6
       rc = -q*mu*(mu-1)
       x2 = (mu+q) / (mu-q)      x3 = (1-(mu+q)) / (1-(mu-q))
  get the curve
    rc * y^2*z = x^3  +  a4 * x * z^2  +  a6 * z^3
      by (x,y,z) -> (x+A2/3*z,y,z) 
      a4 = A4 - A2^2/3    a6 = A6 - A2*A4/3 + 2*A2^3/27
 
  It get the coefficient a4 and a6 of the short weierstrass curve
                           (a4 is put in T->A)
                         A2 of the long weierstrass form
*/
int coeff_EllW(mpz_t f,mpmod_t n,
	       curve *T,mpres_t A2, mpres_t a6,
	       curveHyperEll cHEll ) {

  mpres_t t1,t2,t3,t4;

  mpres_init (t1,n);
  mpres_init (t2,n);
  mpres_init (t3,n);

  mpres_sub (t1, cHEll->mu, cHEll->q ,n); // t1 = mu-q 
  mpres_ui_sub (t2 ,1 ,t1 ,n); // t2 = 1-(mu-q)

  mpres_mul (t3 ,t1 ,t2 ,n);
  if ( !mpres_invert(t3,t3,n) ) { // t3 = 1 / ( (mu-q)*(1-(mu-q)) )
    mpres_gcd(f ,t3 ,n);
    mpres_clear (t1,n);
    mpres_clear (t2,n);
    mpres_clear (t3,n);
    return MORPHISM_FAIL;
  }

  mpres_init (t4,n);

  mpres_add (t4, cHEll->mu, cHEll->q ,n); // t4 = mu+q 


  mpres_mul (t2 ,t3 ,t2 ,n); //  1 / (mu-q)
  mpres_mul (t2 ,t2 ,t4 ,n); // t2 = (mu+q) / (mu-q) = x2

  mpres_mul (t3 ,t3 ,t1 ,n); // 1 / (1-(mu-q))
  mpres_ui_sub (t4 ,1 ,t4 ,n);
  mpres_mul (t3 ,t3 ,t4 ,n); // t3 = (1-(mu+q)) / (1-(mu-q)) = x3


  mpres_mul (t2 ,t2 ,t2 ,n); // t2 = x2^2
  mpres_mul (t3 ,t3 ,t3 ,n); // t3 = x3^2


  mpres_mul (a6 ,t2 ,t3 ,n); // x2^2*x3^2
  mpres_add (T->A ,t2 ,t3 ,n); // x2^2+x3^2


  mpres_clear (t4,n);


  mpres_add_ui (A2 ,T->A ,1 ,n);
  mpres_neg (A2 ,A2 ,n);         //      A2 = - (1+x2^2+x3^2)
  mpres_add (T->A ,T->A ,a6 ,n); // A  = A4 = x2^2*x3^2 + x2^2+x3^2
  mpres_neg (a6 ,a6 ,n);         // a6 = A6 = - x2^2*x3^2

  mpres_mul_ui (a6 ,a6 ,27 ,n);
  mpres_mul (t1 ,A2 ,T->A ,n);
  mpres_mul_ui (t1 ,t1 ,9 ,n);
  mpres_sub (a6 ,a6 ,t1 ,n); // 27*A6 - 9*A2*A4
  mpres_mul (t1 ,A2 ,A2 ,n); // A2^2
  mpres_mul (t3 ,t1 ,t2 ,n);
  mpres_mul_ui (t3 ,t3 ,2 ,n);
  mpres_add (a6 ,a6 ,t3 ,n); // 27*A6 - 9*A2*A4 +2*A2^3

  mpres_mul_ui (T->A ,T->A ,3 ,n);
  mpres_sub (T->A ,T->A ,t1 ,n); // 3*A4 - A2^2


  mpres_set_ui (t1 ,3 ,n);

  if ( !mpres_invert(t1,t1,n) ) { // t1 = 1 /3 
    mpres_gcd(f ,t1 ,n);
    mpres_clear (t1,n);
    mpres_clear (t2,n);
    mpres_clear (t3,n);
    return MORPHISM_FAIL; 
  }

  mpres_mul (T->A ,T->A ,t1 ,n); // A = a4 = (3*A4 - A2^2) / 3
  mpres_mul (t2 ,t1 ,t1 ,n);
  mpres_mul (t2 ,t2 ,t1 ,n);
  mpres_mul (a6 ,a6 ,t2 ,n); // a6 = (27*A6 - 9*A2*A4 +2*A2^3) / 27

  mpres_clear (t1,n);
  mpres_clear (t2,n);
  mpres_clear (t3,n);
  return MORPHISM;
}








/*
  send a point from the long Weiestrass form to the short Weierstrass twisted       form by the following morphisms:
    (x,y,z) -> (x+A2/3*z,y,z)
    (x,y,z) -> (x,y,R*z) 
  See function coeff_EllW_EllWtwist for notations
*/
void EllW_to_EllWshort(mpmod_t n,mpalgpol_t pol,
		       mpalgres_t R,mpres_t A2, 
		       mpalgres_t x,mpalgres_t y,mpalgres_t z,
		       mpalgres_t aX,mpalgres_t aY,mpalgres_t aZ) {

  mpalgres_t aT1,aT2,aT3,aT4;

  mpalgres_init (aT1 ,n);
  mpalgres_init (aT2 ,n);
  mpalgres_init (aT3 ,n);
  mpalgres_init (aT4 ,n);

  mpalgres_mul_ui (aT1 ,aX ,3 ,pol ,n);
  mpalgres_mul_mpres (aT4 ,aZ ,A2 ,pol ,n); 
  mpalgres_add (aT1 ,aT1 ,aT4 ,pol ,n);  // x =3x+A2*z
  mpalgres_mul_ui (aT2 ,aY ,3 ,pol ,n); // y =3y
  mpalgres_mul_ui (aT3 ,aZ ,3 ,pol ,n); // z =3z

  mpalgres_mul (aT3 ,aT3 ,R ,pol ,n); // z *=R


  mpalgres_set (x ,aT1 ,n);
  mpalgres_set (y ,aT2 ,n);
  mpalgres_set (z ,aT3 ,n);

  mpalgres_clear (aT1 ,n);
  mpalgres_clear (aT2 ,n);
  mpalgres_clear (aT3 ,n);
  mpalgres_clear (aT4 ,n);

}







/*
  Send the point (x::z) from the hyperelliptic curve to the elliptic curve in
  short weierstrass form defined by q and put the result in T.
  Don't touch to x,z

  (x::z) -> ( ( x-(mu+q)z )^2 :: ( x-(mu-q)z )^2 )
  (x::z) -> ( 3x + A2*z :: 3z )
  (x::z) -> (x/z :: 1)
  get B = x^3 + a4*x + a6
  Set (T->x,T->y)=(x/B,y/B) and T->A = a4/B^2
*/
int HEll_EllW_degree1 (mpz_t f,mpmod_t n,
		       mpres_t x,mpres_t z,
		       curve *T,mpres_t A2, mpres_t a6,
		       curveHyperEll cHEll ) {
  
  mpres_t t;
  mpres_t X,Z;

  mpres_init (t,n);
  mpres_init (X,n);
  mpres_init (Z,n);

  mpres_add (X ,cHEll->mu ,cHEll->q ,n);
  mpres_mul (X ,X ,z ,n);
  mpres_sub (X ,x ,X ,n);
  mpres_mul (X ,X ,X ,n); // ( x-(mu+q)z )^2

  mpres_sub (Z ,cHEll->mu ,cHEll->q ,n);
  mpres_mul (Z ,Z ,z ,n);
  mpres_sub (Z ,x ,Z ,n);
  mpres_mul (Z ,Z ,Z ,n); // ( x-(mu-q)z )^2



  mpres_mul_ui (X ,X ,3 ,n);
  mpres_mul (t ,Z ,A2 ,n);
  mpres_add (X ,X ,t ,n);   // 3x + A2*z
  mpres_mul_ui (Z ,Z ,3 ,n); // 3z

  if ( mpres_is_zero (Z ,n) ){
    mpres_gcd (f ,Z ,n);
    mpres_clear (t,n);
    mpres_clear (X,n);
    mpres_clear (Z,n);
    return MORPHISM_FOUND_ZERO_CURVE;
  }
  if ( !mpres_invert(t ,Z ,n) ){
    mpres_gcd (f ,Z ,n);
    mpres_clear (t,n);
    mpres_clear (X,n);
    mpres_clear (Z,n);
    return MORPHISM_STEP1;
  }
  
  mpres_mul (X ,X ,Z ,n);

  mpres_mul (t ,X ,X ,n);
  mpres_mul (t ,t ,X ,n);
  mpres_mul (Z ,T->A ,X ,n);
  mpres_add (Z ,Z ,t ,n);
  mpres_add (Z ,Z ,a6 ,n); // x^3 + a4*x + a6

  if ( !mpres_invert(T->y ,Z ,n) ){
    mpres_gcd (f ,Z ,n);
    mpres_clear (t,n);
    mpres_clear (X,n);
    mpres_clear (Z,n);
    return MORPHISM_FAIL;
  }

  mpres_mul (T->x ,T->x ,T->y ,n);
  mpres_mul (T->A ,T->A ,T->y ,n);
  mpres_mul (T->A ,T->A ,T->y ,n);


  mpres_clear (t,n);
  mpres_clear (X,n);
  mpres_clear (Z,n);
  return MORPHISM;
}



/*
  double the point on the weierstrass curve
  (x,y)->( la^2-2x , -la*(la^2-2x)-mu )
    la = (3x^2+a4)/(2y)    mu = y-x*la
  return 0 if failed
*/
int double_short_weierstrass (mpz_t f,mpmod_t n,curve *T) {
  mpres_t t1,t2,t3;

  mpres_init (t1 ,n);
  mpres_init (t2 ,n);
  mpres_init (t3 ,n);

  mpres_mul (t1 ,T->x ,T->x ,n);
  mpres_mul_ui (t1 ,t1 ,3 ,n);
  mpres_add (t1 ,t1 ,T->A ,n); // (3x^2+a4)

  mpres_mul_ui (t2 ,T->y ,2 ,n);
  if ( !mpres_invert(t2 ,t2 ,n) ){ // 1/2y
    mpres_gcd (f ,t2 ,n);
    mpres_clear (t1,n);
    mpres_clear (t2,n);
    mpres_clear (t3,n);
    return 0;
  }
  
  mpres_mul (t1 ,t1 ,t2 ,n); // la = (3x^2+a4)/(2y)
  mpres_mul (t2 ,T->x ,t1 ,n); 
  mpres_sub (t2 ,T->y ,t2 ,n); // mu = y-x*la

  mpres_mul_ui (T->x ,T->x ,2 ,n);
  mpres_mul (t3 ,t1 ,t1 ,n);
  mpres_sub (T->x ,t3 ,T->x ,n); // la^2-2x 

  mpres_mul (T->y ,t1 ,T->x ,n);
  mpres_add (T->y ,T->y ,t2 ,n);
  mpres_neg (T->y ,T->y ,n); // -la*(la^2-2x)-mu

  mpres_clear (t1,n);
  mpres_clear (t2,n);
  mpres_clear (t3,n);
  return 1;
}






/*
  "add" the two points on the short Weierstrass twisted curve
  We only want the (x::z) coordinate
    y^2*z = x^3  +  a4' * x * z^2  +  a6' * z^3
      a4' = a4/R^2   a6' = a6/R^3
  Then go back to the untwisted curve
    y^2*z = x^3  +  a4 * x * z^2  +  a6 * z^3
    (x,y,z) -> (x*R ,?, z);
  Scale the y coordinate and come back to non projective curve
    (x::z) -> (x/z,?)
    compute B = x^3 + a4 * x + a6   (i.e. y^2)
    (x,?)  -> (X,Y)
       X = x/B    Y=1/B   T->A = a4/B^2
*/
int addW_and_scale (mpz_t f,mpmod_t n,mpalgpol_t pol,
		    curve *T,mpres_t a6,
		    mpalgres_t X1 ,mpalgres_t Y1 ,mpalgres_t Z1,
		    mpalgres_t X2 ,mpalgres_t Y2 ,mpalgres_t Z2,
		    mpalgres_t a4twist ,mpalgres_t R) {
  int test;
  mpalgres_t aT1,aT2;
  mpres_t temp;
  mpalgres_t La,Den;
  mpalgres_t X3,Z3;
  mpz_t F[DEGREE_ALGEBRA];

  mpz_init (F[0]);
  mpz_init (F[1]);
  mpalgres_init (X3 ,n);
  mpalgres_init (Z3 ,n);

  if (mpalgres_is_zero(Z1 ,pol ,n) ) { // the first point is zero
    if (mpalgres_is_zero(Z2 ,pol ,n) ) {// the two points are zero
      mpz_clear (F[0]);
      mpz_clear (F[1]);
      mpalgres_clear (X3 ,n);
      mpalgres_clear (Z3 ,n);
      return MORPHISM_FOUND_ZERO_CURVE;
    }

    // Maybe the two points are zero modulo a factor of n
    mpalgres_gcd (F,Z2,n);
    if ( mpz_cmp_ui (F[0],1)!=0 ) {
      mpz_set (f,F[0]);
      mpz_clear (F[0]);
      mpz_clear (F[1]);
      mpalgres_clear (X3 ,n);
      mpalgres_clear (Z3 ,n);
      return MORPHISM_STEP1;
    }
    if ( mpz_cmp_ui (F[1],1)!=0 ) {
      mpz_set (f,F[1]);
      mpz_clear (F[0]);
      mpz_clear (F[1]);
      mpalgres_clear (X3 ,n);
      mpalgres_clear (Z3 ,n);
      return MORPHISM_STEP1;
    }
    
    mpalgres_set (X3,X2,n);
    mpalgres_set (Z3,Z2,n);
  
  }
  else {
    // Maybe the first points is zero modulo a factor of n
    mpalgres_gcd (F,Z1,n);
    if ( mpz_cmp_ui (F[0],1)!=0 ) {
      mpz_set (f,F[0]);
      mpz_clear (F[0]);
      mpz_clear (F[1]);
      mpalgres_clear (X3 ,n);
      mpalgres_clear (Z3 ,n);
      return MORPHISM_FAIL;
    }
    if ( mpz_cmp_ui (F[1],1)!=0 ) {
      mpz_set (f,F[1]);
      mpz_clear (F[0]);
      mpz_clear (F[1]);
      mpalgres_clear (X3 ,n);
      mpalgres_clear (Z3 ,n);
      return MORPHISM_FAIL;
    }

    // Now the first point is not the zero

    if (mpalgres_is_zero(Z2 ,pol ,n) ) { // the second point is zero
      mpalgres_set (X3,X1,n);
      mpalgres_set (Z3,Z1,n);
    }
    else {
      // Maybe the second points is zero modulo a factor of n
      mpalgres_gcd (F,Z2,n); 
      if ( mpz_cmp_ui (F[0],1)!=0 ) {
	mpz_set (f,F[0]);
	mpz_clear (F[0]);
	mpz_clear (F[1]);
	mpalgres_clear (X3 ,n);
	mpalgres_clear (Z3 ,n);
	return MORPHISM_FAIL;
      }
      if ( mpz_cmp_ui (F[1],1)!=0 ) {
	mpz_set (f,F[1]);
	mpz_clear (F[0]);
	mpz_clear (F[1]);
	mpalgres_clear (X3 ,n);
	mpalgres_clear (Z3 ,n);
	return MORPHISM_FAIL;
      }


      // general case
      
      mpalgres_init (Den,n);
      mpalgres_init (aT1,n);
      mpalgres_init (aT2,n);
  
      mpalgres_mul (Den ,Z1 ,X2 ,pol ,n);
      mpalgres_mul (aT1 ,X1 ,Z2 ,pol ,n);
      mpalgres_sub (Den ,Den ,aT1 ,pol ,n); // Z1*X2-Z2*X1

      if (mpalgres_is_zero(Den ,pol ,n) ) { // X1/Z1 = X2/Z2
	mpalgres_mul (aT1 ,Y1 ,Z2 ,pol ,n);
	mpalgres_mul (aT2 ,Z1 ,Y2 ,pol ,n);
	mpalgres_add (aT1 ,aT1 ,aT2 ,pol ,n);
	
	if (mpalgres_is_zero(aT1 ,pol ,n) ) { // Y1/Z1 = -Y2/Z2
	  // we add a point and its opposite
	  mpalgres_clear (Den,n);
	  mpalgres_clear (aT1,n);
	  mpalgres_clear (aT2,n);
	  mpz_clear (F[0]);
	  mpz_clear (F[1]);
	  mpalgres_clear (X3 ,n);
	  mpalgres_clear (Z3 ,n);
	  return MORPHISM_FOUND_ZERO_CURVE;
	}
	
	// Maybe it is the case modulo a factor of n
	mpalgres_gcd (F,aT1,n); 
	if ( mpz_cmp_ui (F[0],1)!=0 ) {
	  mpz_set (f,F[0]);
	  mpalgres_clear (Den,n);
	  mpalgres_clear (aT1,n);
	  mpalgres_clear (aT2,n);
	  mpz_clear (F[0]);
	  mpz_clear (F[1]);
	  mpalgres_clear (X3 ,n);
	  mpalgres_clear (Z3 ,n);
	  return MORPHISM_STEP1;
	}
	if ( mpz_cmp_ui (F[1],1)!=0 ) {
	  mpz_set (f,F[1]);
	  mpalgres_clear (Den,n);
	  mpalgres_clear (aT1,n);
	  mpalgres_clear (aT2,n);
	  mpz_clear (F[0]);
	  mpz_clear (F[1]);
	  mpalgres_clear (X3 ,n);
	  mpalgres_clear (Z3 ,n);
	  return MORPHISM_STEP1;
	}
	
	
	
	// Y1/Z1 = Y2/Z2 // doubling case
	
	if (mpalgres_is_zero(Y1 ,pol ,n)) { // doubling a 2-torsion point
	  mpalgres_clear (Den,n);
	  mpalgres_clear (aT1,n);
	  mpalgres_clear (aT2,n);
	  mpz_clear (F[0]);
	  mpz_clear (F[1]);
	  mpalgres_clear (X3 ,n);
	  mpalgres_clear (Z3 ,n);
	  return MORPHISM_FOUND_ZERO_CURVE;
	}

	// Maybe it is the case modulo a factor of n
	mpalgres_gcd (F,Y1,n); 
	if ( mpz_cmp_ui (F[0],1)!=0 ) {
	  mpz_set (f,F[0]);
	  mpalgres_clear (Den,n);
	  mpalgres_clear (aT1,n);
	  mpalgres_clear (aT2,n);
	  mpz_clear (F[0]);
	  mpz_clear (F[1]);
	  mpalgres_clear (X3 ,n);
	  mpalgres_clear (Z3 ,n);
	  return MORPHISM_STEP1;
	}
	if ( mpz_cmp_ui (F[1],1)!=0 ) {
	  mpz_set (f,F[1]);
	  mpalgres_clear (Den,n);
	  mpalgres_clear (aT1,n);
	  mpalgres_clear (aT2,n);
	  mpz_clear (F[0]);
	  mpz_clear (F[1]);
	  mpalgres_clear (X3 ,n);
	  mpalgres_clear (Z3 ,n);
	  return MORPHISM_STEP1;
	}

	// general doubling case
	
	mpalgres_init (La,n);
	


	mpalgres_mul (Den,Y1,Z1,pol,n);
	mpalgres_mul_ui (Den,Den,2,pol,n); // Den = 2*Y1*Z1
    
	mpalgres_mul (La ,X1 ,X1 ,pol,n);
	mpalgres_mul_ui (La ,La ,3 ,pol,n);
	mpalgres_mul (aT1 ,Z1 ,Z1 ,pol ,n);
	mpalgres_mul (aT1 ,aT1 ,a4twist ,pol ,n);
	mpalgres_add (La,La,aT1 ,pol,n); // La = 3*X1^2 + a4'*Z1^2

	mpalgres_mul (X3 ,La ,La ,pol,n);
	mpalgres_mul (X3 ,X3 ,Z1 ,pol,n);

	mpalgres_mul (Z3 ,Den ,Den ,pol,n); // Den^2
	mpalgres_mul (aT1 ,Z3 ,X1 ,pol ,n);
	mpalgres_mul_ui (aT1 ,aT1 ,2 ,pol ,n);
	mpalgres_sub (X3 ,X3 ,aT1 ,pol ,n); // X3 = La^2*Z1 - 2*X1*Den^2

	mpalgres_mul (Z3 ,Z3 ,Z1 ,pol ,n); // Z3 = Z1*Den^2
	
    

	mpalgres_clear (La,n);
	mpalgres_clear (Den,n);
	mpalgres_clear (aT1,n);
	mpalgres_clear (aT2,n);
	
      }
    
      else {

	// Maybe we have to double one point modulo a factor of n
	mpalgres_gcd (F,Den,n); 
	if ( mpz_cmp_ui (F[0],1)!=0 ) {
	  mpz_set (f,F[0]);
	  mpalgres_clear (Den,n);
	  mpalgres_clear (aT1,n);
	  mpalgres_clear (aT2,n);
	  mpz_clear (F[0]);
	  mpz_clear (F[1]);
	  mpalgres_clear (X3 ,n);
	  mpalgres_clear (Z3 ,n);
	  return MORPHISM_FAIL;
	}
	if ( mpz_cmp_ui (F[1],1)!=0 ) {
	  mpz_set (f,F[1]);
	  mpalgres_clear (Den,n);
	  mpalgres_clear (aT1,n);
	  mpalgres_clear (aT2,n);
	  mpz_clear (F[0]);
	  mpz_clear (F[1]);
	  mpalgres_clear (X3 ,n);
	  mpalgres_clear (Z3 ,n);
	  return MORPHISM_FAIL;
	}




	// adding case   

	mpalgres_init (La,n);

	

	mpalgres_mul (La ,Z1 ,Y2 ,pol ,n);
	mpalgres_mul (aT1 ,Y1 ,Z2 ,pol ,n);
	mpalgres_sub (La ,La ,aT1 ,pol ,n); // La = Z1*Y2-Y1*Z2

	
	mpalgres_mul (aT2 ,Z1 ,Z2, pol ,n);  // Z1*Z2
	mpalgres_mul (La ,La ,La ,pol ,n); 
	mpalgres_mul (La ,La ,aT2 ,pol ,n); // La^2 * Z1*Z2
	
	mpalgres_mul (X3 ,Z1 ,X2 ,pol ,n);
	mpalgres_mul (aT1 ,X1 ,Z2 ,pol ,n);
	mpalgres_add (X3 ,X3 ,aT1 ,pol ,n); // Z1*X2 + X1*Z2

	mpalgres_mul (Z3 ,Den ,Den ,pol ,n); // Den^2
	mpalgres_mul (X3 ,X3 ,Z3 ,pol,n);
	mpalgres_sub (X3 ,La ,X3 ,pol ,n);
	// X3 = La^2*Z1*Z2 - (Z1*X2+X1*Z2)*Den^2

	mpalgres_mul (Z3 ,Z3 ,aT2 ,pol ,n); // Z3 = Z1*Z2*Den^2

	
    
        mpalgres_clear (aT1,n);
	mpalgres_clear (aT2,n);         
	mpalgres_clear (La,n);
	mpalgres_clear (Den,n);

      } 
    }
  }
  
  

  // We have the point (X3::Z3) on the short Weierstrass twisted curve

  mpalgres_mul (X3, X3 ,R ,pol,n);
  



  if ( mpalgres_is_zero (Z3 ,pol ,n) ) { // Z3 =0
    mpz_clear (F[0]);
    mpz_clear (F[1]);
    mpalgres_clear (X3,n);
    mpalgres_clear (Z3,n);
    return MORPHISM_FOUND_ZERO_CURVE;
  }
  mpalgres_gcd (F,Z3,n); 
  if ( mpz_cmp_ui (F[0],1)!=0 ) { // maybe Z3 is 0 modulo a factor of n
    mpz_set (f,F[0]);
    mpz_clear (F[0]);
    mpz_clear (F[1]);
    mpalgres_clear (X3 ,n);
    mpalgres_clear (Z3 ,n);
    return MORPHISM_STEP1;
  }
  if ( mpz_cmp_ui (F[1],1)!=0 ) {
    mpz_set (f,F[1]);
    mpz_clear (F[0]);
    mpz_clear (F[1]);
    mpalgres_clear (X3 ,n);
    mpalgres_clear (Z3 ,n);
    return MORPHISM_STEP1;
  }
  mpz_clear (F[0]);
  mpz_clear (F[1]);


  test = mpalgres_invert (Z3 ,Z3 ,pol,n,f);
  if (test == -1) {
    mpalgres_clear (X3,n);
    mpalgres_clear (Z3,n);
    return MORPHISM_ERROR;
  }
  if (test == 0) {
    mpalgres_clear (X3,n);
    mpalgres_clear (Z3,n);
    return MORPHISM_FAIL;
  }

  mpalgres_mul (X3 ,X3 ,Z3 ,pol ,n); // X3/Z3


  if ( !mpres_is_zero(X3[1],n) ) {
    // This should never happen since X3 should be rationnal
    mpalgres_clear (X3,n);
    mpalgres_clear (Z3,n);
    return MORPHISM_ERROR; 
  }

  mpz_set (T->x ,X3[0]);
  


  // Compute B=x^3 + a4 * x + a6   (i.e. y^2)
  mpres_init(temp,n);

  mpres_mul (T->y ,T->x ,T->x ,n);
  mpres_mul (T->y ,T->y ,T->x ,n); // x^3
  mpres_mul (temp ,T->A ,T->x ,n); // a4 * x
  mpres_add (T->y ,T->y ,temp ,n);
  mpres_add (T->y ,T->y ,a6 ,n);   // B=x^3 + a4 * x + a6

  mpres_clear (temp,n);


  if ( !mpres_invert(T->y ,T->y ,n) ) { // 1/B
    mpres_gcd(f,T->y,n);
    mpalgres_clear (X3,n);
    mpalgres_clear (Z3,n);
    return MORPHISM_FAIL; 
  }

  mpres_mul (T->x ,T->x ,T->y ,n); // X/B
  mpres_mul (T->A ,T->A ,T->y ,n); 
  mpres_mul (T->A ,T->A ,T->y ,n); // a4/B^2

  mpalgres_clear (X3,n);
  mpalgres_clear (Z3,n);
  return MORPHISM;
}


