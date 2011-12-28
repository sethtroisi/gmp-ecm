#include "ariKS.h"

void ksPoint_init(ksPoint p,mpmod_t n) {
  mpres_init (p->X,n);
  mpres_init (p->Y,n);
  mpres_init (p->Z,n);
  mpres_init (p->T,n);
}

void ksPoint_clear(ksPoint p,mpmod_t n) {
  mpres_clear (p->X,n);
  mpres_clear (p->Y,n);
  mpres_clear (p->Z,n);
  mpres_clear (p->T,n);
}



void ksCstPourMul_init (ksCstPourMul cMul,mpmod_t n) {
  mpres_init (cMul->invZ, n);
  mpres_init (cMul->invT, n);
  mpres_init (cMul->x0p, n);
  mpres_init (cMul->t0p, n);
  mpres_init (cMul->Y0, n);
  mpres_init (cMul->Z0, n);
}

void ksCstPourMul_clear (ksCstPourMul cMul,mpmod_t n) {
  mpres_clear (cMul->invZ, n);
  mpres_clear (cMul->invT, n);
  mpres_clear (cMul->x0p, n);
  mpres_clear (cMul->t0p, n);
  mpres_clear (cMul->Y0, n);
  mpres_clear (cMul->Z0, n);
}







/*
  Multiply the point P on the Kummer surface by k and put the result in P
  The point P is such that X=-Y
  The zero of the Kummer surface is [1,be,ga,1] since al=de=1
*/
void mulKS(ksPoint P,		
	   ksCstPourMul cMul,
	   mpres_t be,mpres_t ga,
	   mpz_t k,						
	   mpmod_t n) {

  int p;
  //  printf("warning: p est un int et on y met un size_t\n");
  ksPoint Pi;
  mpres_t u,v;
  mpres_init (u,n);
  mpres_init (v,n);

  ksPoint_init (Pi,n);



  if ( !mpz_sgn(k)) { // k=0
    mpres_set_ui (P->X,1,n);
    mpres_set (P->Y,be,n); 
    mpres_set (P->Z,ga,n); 
    mpres_set_ui (P->T,1,n);
  }
 
  else if ( mpz_sgn(k) < 0 ) {  //  k < 0
    mpz_neg (k,k); // k=-k
  }

  else if ( !mpz_cmp_ui(k,2) ) { // k=2
    doubleKS2(P,cMul,u,v,n);
  }


  else if ( mpz_cmp_ui(k,2) > 0 ) { // k>2
    doubleKS(Pi,P,cMul,u,v,n);

    for (p=(mpz_sizeinbase(k,2)-2);p>=0;p--) {
      if ( mpz_tstbit (k,p) ) {
	loopMulKS(Pi,P,cMul,u,v,n);
      }
      else {
	loopMulKS(P,Pi,cMul,u,v,n);
      }
      // 
      // We could save a few operations on the last bit of k
      // We have computed k*P and (k+1)*P.
    }
  }

  ksPoint_clear (Pi,n);
  mpres_clear (u,n);
  mpres_clear (v,n);
}






/*
  Double a point P=[X,Y,Z,T] on the Kummer surface and put the result in P
*/
void doubleKS2(ksPoint P,				   
	       ksCstPourMul cMul,
	       mpres_t u,mpres_t v,
	       mpmod_t n) {

  ksPoint Q;
  ksPoint_init (Q,n);

  /* magma
  Xi := x0p*(X + Y + Z + T)^2;
  Yi :=     (X + Y - Z - T)^2;
  Zi := z0p*(X - Y + Z - T)^2;    // z0p=-1
  Ti := t0p*(X - Y - Z + T)^2;
 
  x2 :=    (Xi + Yi + Zi + Ti)^2;
  y2 := Y0*(Xi + Yi - Zi - Ti)^2;
  z2 := Z0*(Xi - Yi + Zi - Ti)^2;
  t2 := T0*(Xi - Yi - Zi + Ti)^2; // T0=1
  */

  hadamard (Q,P,u,v,n);
  // Pi->X = X + Y + Z + T 
  // Pi->Y = X + Y - Z - T 
  // Pi->Z = X - Y + Z - T
  // Pi->T = X - Y - Z + T



  mpres_mul (Q->X ,Q->X ,Q->X ,n); 
  mpres_mul (Q->Y ,Q->Y ,Q->Y ,n); // Yi
  mpres_mul (Q->Z ,Q->Z ,Q->Z ,n);
  mpres_mul (Q->T ,Q->T ,Q->T ,n);

  mpres_neg (Q->Z ,Q->Z ,n); // Zi 
  // We could change the sign of Zi later
  mpres_mul (Q->X ,Q->X ,cMul->x0p ,n);  // Xi
  mpres_mul (Q->T ,Q->T ,cMul->t0p ,n);  // Yi


  hadamard (P,Q,u,v,n);
  // P->X = Xi + Yi + Zi + Ti
  // P->Y = Xi + Yi - Zi - Ti
  // P->Z = Xi - Yi + Zi - Ti
  // P->T = Xi - Yi - Zi + Ti



  mpres_mul (P->X ,P->X ,P->X ,n);  // x2
  mpres_mul (P->Y ,P->Y ,P->Y ,n); 
  mpres_mul (P->Z ,P->Z ,P->Z ,n);
  mpres_mul (P->T ,P->T ,P->T ,n);  // t2

  mpres_mul (P->Y ,P->Y ,cMul->Y0 ,n); // y2
  mpres_mul (P->Z ,P->Z ,cMul->Z0 ,n); // z2



  ksPoint_clear (Q,n);


}






/*
  Double a point P the Kummer surface and put the result in P2
*/
void doubleKS(ksPoint P2,const ksPoint P,				
	      ksCstPourMul cMul,
	      mpres_t u,mpres_t v,
	      mpmod_t n) {

  ksPoint Pi;
  ksPoint_init (Pi,n);



  /* magma
  Xi := x0p*(X + Y + Z + T)^2;
  Yi :=     (X + Y - Z - T)^2;
  Zi := z0p*(X - Y + Z - T)^2;    // z0p=-1
  Ti := t0p*(X - Y - Z + T)^2;
 
  X2 :=    (Xi + Yi + Zi + Ti)^2;
  Y2 := Y0*(Xi + Yi - Zi - Ti)^2;
  Z2 := Z0*(Xi - Yi + Zi - Ti)^2;
  T2 := T0*(Xi - Yi - Zi + Ti)^2; // T0=1
  */



  hadamard (Pi,P,u,v,n);
  // Pi->x = (X + Y + Z + T)
  // Pi->y = (X + Y - Z - T)
  // Pi->z = (X - Y + Z - T)
  // Pi->t = (X - Y - Z + T)



  mpres_mul (Pi->X ,Pi->X ,Pi->X ,n); 
  mpres_mul (Pi->Y ,Pi->Y ,Pi->Y ,n); // Yi
  mpres_mul (Pi->Z ,Pi->Z ,Pi->Z ,n);
  mpres_mul (Pi->T ,Pi->T ,Pi->T ,n);
  mpres_neg (Pi->Z ,Pi->Z ,n); // Zi 
  // We could change the sign of Zi later
  mpres_mul (Pi->X ,Pi->X ,cMul->x0p ,n);   // Xi
  mpres_mul (Pi->T ,Pi->T ,cMul->t0p ,n);   // Ti


  hadamard (P2,Pi,u,v,n);
  // P2->X = Xi + Yi + Zi + Ti
  // P2->Y = Xi + Yi - Zi - Ti
  // P2->Z = Xi - Yi + Zi - Ti
  // P2->T = Xi - Yi - Zi + Ti



  mpres_mul (P2->X ,P2->X ,P2->X ,n); // X2
  mpres_mul (P2->Y ,P2->Y ,P2->Y ,n); 
  mpres_mul (P2->Z ,P2->Z ,P2->Z ,n);
  mpres_mul (P2->T ,P2->T ,P2->T ,n); // T2

  mpres_mul (P2->Y ,P2->Y ,cMul->Y0 ,n); // Y2
  mpres_mul (P2->Z ,P2->Z ,cMul->Z0 ,n); // Z2



  ksPoint_clear (Pi,n);
}






/*
  One step during the multiplication loop.
  We have Pm=[X,Y,Z,T] and Pp=[x,y,z,t] with Pm-Pp equal the initial point
  We want Pm <- 2*Pm=[X,Y,Z,T] and Pp <- Pm+Pp=[x,y,z,t]
*/
void loopMulKS(ksPoint Pm,ksPoint Pp,
		 ksCstPourMul cMul,
		 mpres_t u,mpres_t v,
		 mpmod_t n) {

  ksPoint Qm,Qp;
  mpres_t U,V;

  ksPoint_init (Qm,n); // Let Qm = [X2,Y2,Z2,T2]
  ksPoint_init (Qp,n); // Let Qp = [x2,y2,z2,t2]

  mpres_init (U, n); 
  mpres_init (V, n);



  hadamard (Qm,Pm,u,v,n);
  // Qm->X = X2 = X + Y + Z + T
  // Qm->Y = Y2 = X + Y - Z - T
  // Qm->Z = Z2 = X - Y + Z - T
  // Qm->T = T2 = X - Y - Z + T

  mpres_mul (U ,Qm->X ,cMul->x0p ,n);  // U = x0p * (X + Y + Z + T)
  mpres_mul (V ,Qm->T ,cMul->t0p ,n);  // V = t0p * (X - Y - Z + T)


  hadamard (Qp,Pp,u,v,n);
  // Qp->X = x2 = x + y + z + t
  // Qp->Y = y2 = x + y - z - t
  // Qp->Z = z2 = x - y + z - t
  // Qp->T = t2 = x - y - z + t



  // PseudoAdd

  mpres_mul (Qp->X ,Qp->X ,U ,n);     // Qp->X = x0p * (X+Y+Z+T) * (x+y+z+t)
  mpres_mul (Qp->Y ,Qp->Y ,Qm->Y ,n); // Qp->Y =       (X+Y-Z-T) * (x+y-z-t)
  mpres_mul (Qp->Z ,Qp->Z ,Qm->Z ,n); 
  mpres_neg (Qp->Z ,Qp->Z ,n);        // Qp->Z = z0p * (X-Y+Z-T) * (x-y+z-t)
  mpres_mul (Qp->T ,Qp->T ,V ,n);     // Qp->T = t0p * (X-Y-Z+T) * (x-y-z+t)

  hadamard (Pp,Qp,u,v,n);
  // Pp->X = x2 + y2 + z2 + t2
  // Pp->Y = x2 + y2 - z2 - t2
  // Pp->Z = x2 - y2 + z2 - t2
  // Pp->T = x2 - y2 - z2 + t2

  mpres_mul (Pp->X ,Pp->X ,Pp->X ,n);      // Pp->X =        (x2+y2+z2+t2)^2
  mpres_mul (Pp->Y ,Pp->Y ,Pp->Y ,n);
  mpres_neg (Pp->Y ,Pp->Y ,n);             // Pp->Y = invY * (x2+y2-z2-t2)^2
  mpres_mul (Pp->Z ,Pp->Z ,Pp->Z ,n);
  mpres_mul (Pp->Z ,Pp->Z ,cMul->invZ ,n); // Pp->Z = invZ * (x2-y2+z2-t2)^2
  mpres_mul (Pp->T ,Pp->T ,Pp->T ,n);
  mpres_mul (Pp->T ,Pp->T ,cMul->invT ,n); // Pp->T = invT * (x2-y2-z2+t2)^2



  // Double

  mpres_mul (Qm->X ,Qm->X ,U ,n);     // Qm->X = x0p * (X + Y + Z + T)^2
  mpres_mul (Qm->Y ,Qm->Y ,Qm->Y ,n); // Qm->Y =       (X + Y - Z - T)^2
  mpres_mul (Qm->Z ,Qm->Z ,Qm->Z ,n); 
  mpres_neg (Qm->Z ,Qm->Z ,n);        // Qm->Z = z0p * (X - Y + Z - T)^2
  mpres_mul (Qm->T ,Qm->T ,V ,n);     // Qm->T = t0p * (X - Y - Z + T)^2

  hadamard (Pm,Qm,u,v,n);
  // Pm->X = X2 + Y2 + Z2 + T2
  // Pm->Y = X2 + Y2 - Z2 - T2
  // Pm->Z = X2 - Y2 + Z2 - T2
  // Pm->T = X2 - Y2 - Z2 + T2

  mpres_mul (Pm->X ,Pm->X ,Pm->X ,n);    // Pm->X =      (X2 + Y2 + Z2 + T2)^2
  mpres_mul (Pm->Y ,Pm->Y ,Pm->Y ,n);
  mpres_mul (Pm->Y ,Pm->Y ,cMul->Y0 ,n); // Pm->Y = Y0 * (X2 + Y2 - Z2 - T2)^2
  mpres_mul (Pm->Z ,Pm->Z ,Pm->Z ,n);
  mpres_mul (Pm->Z ,Pm->Z ,cMul->Z0 ,n); // Pm->Z = Z0 * (X2 - Y2 + Z2 - T2)^2
  mpres_mul (Pm->T ,Pm->T ,Pm->T ,n);    // Pm->T =      (X2 - Y2 - Z2 + T2)^2



  ksPoint_clear (Qm,n);
  ksPoint_clear (Qp,n);

  mpres_clear (U, n);
  mpres_clear (V, n);

}














// ************* Small parameters *****************








/*
  Double a point P=[X,Y,Z,T] on the Kummer surface and put the result in P2
  We use small parameters
*/
void doubleKSsmallParam(ksPoint P2,const ksPoint P,			
			  ksSmallConstPourMul cMul,		
			  mpres_t u,mpres_t v,
			  mpmod_t n) {

  ksPoint Pi;
  ksPoint_init (Pi,n);



  /* magma
  Xi := x0p*(X + Y + Z + T)^2;
  Yi := y0p*(X + Y - Z - T)^2;
  Zi := z0p*(X - Y + Z - T)^2; 
  Ti := t0p*(X - Y - Z + T)^2;
 
  X2 := X0*(Xi + Yi + Zi + Ti)^2;
  Y2 := Y0*(Xi + Yi - Zi - Ti)^2;
  Z2 := Z0*(Xi - Yi + Zi - Ti)^2;
  T2 := T0*(Xi - Yi - Zi + Ti)^2;
  */



  hadamard (Pi,P,u,v,n);
  // Pi->x = (X + Y + Z + T)
  // Pi->y = (X + Y - Z - T)
  // Pi->z = (X - Y + Z - T)
  // Pi->t = (X - Y - Z + T)

  mpres_mul (Pi->X ,Pi->X ,Pi->X ,n); 
  mpres_mul (Pi->Y ,Pi->Y ,Pi->Y ,n);
  mpres_mul (Pi->Z ,Pi->Z ,Pi->Z ,n);
  mpres_mul (Pi->T ,Pi->T ,Pi->T ,n);

  mpres_muldivbysomething_si (Pi->X ,Pi->X ,cMul->x0p ,n);
  mpres_muldivbysomething_si (Pi->Y ,Pi->Y ,cMul->y0p ,n);
  mpres_muldivbysomething_si (Pi->Z ,Pi->Z ,cMul->z0p ,n);
  mpres_muldivbysomething_si (Pi->T ,Pi->T ,cMul->t0p ,n);



  hadamard (P2,Pi,u,v,n);
  // P2->X = Xi + Yi + Zi + Ti
  // P2->Y = Xi + Yi - Zi - Ti
  // P2->Z = Xi - Yi + Zi - Ti
  // P2->T = Xi - Yi - Zi + Ti

  mpres_mul (P2->X ,P2->X ,P2->X ,n);
  mpres_mul (P2->Y ,P2->Y ,P2->Y ,n); 
  mpres_mul (P2->Z ,P2->Z ,P2->Z ,n);
  mpres_mul (P2->T ,P2->T ,P2->T ,n);

  mpres_muldivbysomething_si (P2->X ,P2->X ,cMul->X0 ,n);
  mpres_muldivbysomething_si (P2->Y ,P2->Y ,cMul->Y0 ,n);
  mpres_muldivbysomething_si (P2->Z ,P2->Z ,cMul->Z0 ,n);
  mpres_muldivbysomething_si (P2->T ,P2->T ,cMul->T0 ,n);



  ksPoint_clear (Pi,n);
}






/*
  Double a point P=[X,Y,Z,T] on the Kummer surface and put the result in P
  We use small parameters
*/
void doubleKSsmallParam2(ksPoint P,				
			   ksSmallConstPourMul cMul,	
			   ksPoint Pi,mpres_t u,mpres_t v,
			   mpmod_t n) {

  /* magma
  Xi := x0p*(X + Y + Z + T)^2;
  Yi := y0p*(X + Y - Z - T)^2;
  Zi := z0p*(X - Y + Z - T)^2; 
  Ti := t0p*(X - Y - Z + T)^2;
 
  X := X0*(Xi + Yi + Zi + Ti)^2;
  Y := Y0*(Xi + Yi - Zi - Ti)^2;
  Z := Z0*(Xi - Yi + Zi - Ti)^2;
  T := T0*(Xi - Yi - Zi + Ti)^2;
  */


  hadamard (Pi,P,u,v,n);
  // Pi->x = (X + Y + Z + T)
  // Pi->y = (X + Y - Z - T)
  // Pi->z = (X - Y + Z - T)
  // Pi->t = (X - Y - Z + T)

  mpres_mul (Pi->X ,Pi->X ,Pi->X ,n); 
  mpres_mul (Pi->Y ,Pi->Y ,Pi->Y ,n);
  mpres_mul (Pi->Z ,Pi->Z ,Pi->Z ,n);
  mpres_mul (Pi->T ,Pi->T ,Pi->T ,n);

  mpres_muldivbysomething_si (Pi->X ,Pi->X ,cMul->x0p ,n);
  mpres_muldivbysomething_si (Pi->Y ,Pi->Y ,cMul->y0p ,n);
  mpres_muldivbysomething_si (Pi->Z ,Pi->Z ,cMul->z0p ,n);
  mpres_muldivbysomething_si (Pi->T ,Pi->T ,cMul->t0p ,n);



  hadamard (P,Pi,u,v,n);
  // P->X = Xi + Yi + Zi + Ti
  // P->Y = Xi + Yi - Zi - Ti
  // P->Z = Xi - Yi + Zi - Ti
  // P->T = Xi - Yi - Zi + Ti

  mpres_mul (P->X ,P->X ,P->X ,n);
  mpres_mul (P->Y ,P->Y ,P->Y ,n); 
  mpres_mul (P->Z ,P->Z ,P->Z ,n);
  mpres_mul (P->T ,P->T ,P->T ,n);

  mpres_muldivbysomething_si (P->X ,P->X ,cMul->X0 ,n);
  mpres_muldivbysomething_si (P->Y ,P->Y ,cMul->Y0 ,n);
  mpres_muldivbysomething_si (P->Z ,P->Z ,cMul->Z0 ,n);
  mpres_muldivbysomething_si (P->T ,P->T ,cMul->T0 ,n);


}




/*
  Multiply a point P on the Kummer surface by k and put the result in P
  We use small paramters
  The zero on the Kummer surface is [1,be,ga,1] since al=de=1
*/
void mulKSsmallParam(ksPoint P,		
		       ksSmallConstPourMul cMul,
		       mpres_t be,mpres_t ga,
		       mpz_t k,						
		       mpmod_t n) {

  int p;
  //  printf("warning: p est un int et on y met un size_t\n");
  ksPoint Q,R;
  mpres_t u,v;

  ksPoint_init (Q,n);
  ksPoint_init (R,n);
  mpres_init (u,n);
  mpres_init (v,n);


  if ( !mpz_sgn(k)) { // k=0
    mpres_set_ui (P->X,1,n);
    mpres_set (P->Y,be,n); 
    mpres_set (P->Z,ga,n); 
    mpres_set_ui (P->T,1,n);
  }
 
  else if ( mpz_sgn(k) < 0 ) {  //  k < 0
    mpz_neg (k,k); // k=-k
  }

  else if ( !mpz_cmp_ui(k,2) ) { // k=2
    doubleKSsmallParam2 (P ,cMul ,R,u,v ,n);
  }


  else if ( mpz_cmp_ui(k,2) > 0 ) { // k>2
    doubleKSsmallParam (Q ,P ,cMul ,u,v ,n);

    for (p=(mpz_sizeinbase(k,2)-2);p>=0;p--) {
      if ( mpz_tstbit (k,p) ) {
	// P  <- P+Q // Q <- 2*Q
	loopMulKSsmallParam (P,Q,cMul,R,u,v,n);

      }
      else {
	// Q <- Q+P   // P  <- 2*P
	loopMulKSsmallParam (Q,P,cMul,R,u,v,n);


      }
      // 
      // We could save a few operations on the last bit of k
      // We have computed k*P and (k+1)*P.
    }
  }

  ksPoint_clear (Q,n);
  ksPoint_clear (R,n);
  mpres_clear (u,n);
  mpres_clear (v,n);
}






/*
  The loop for multiplication on the Kummer surface with small constants
  Imput: two points P,Q
  Output P  <- P+Q  and Q <- 2*Q
*/
void loopMulKSsmallParam (ksPoint P,ksPoint Q,ksSmallConstPourMul cMul,
			  ksPoint R,mpres_t u,mpres_t v,mpmod_t n) {


  //  hadamard (R,Q,u,v,n);
  mpres_add (u ,Q->X  ,Q->Y  ,n); 
  mpres_add (v ,Q->Z  ,Q->T  ,n); 
  mpres_add (R->X ,u ,v ,n);
  mpres_sub (R->Y ,u ,v ,n); 
  mpres_sub (v ,Q->X  ,Q->Y  ,n); 
  mpres_sub (u ,Q->Z  ,Q->T  ,n); 
  mpres_add (R->Z ,v ,u ,n); 
  mpres_sub (R->T ,v ,u ,n); 
  //  hadamard (Q,P,u,v,n);
  mpres_add (u ,P->X  ,P->Y  ,n); 
  mpres_add (v ,P->Z  ,P->T  ,n); 
  mpres_add (Q->X ,u ,v ,n); 
  mpres_sub (Q->Y ,u ,v ,n); 
  mpres_sub (v ,P->X  ,P->Y  ,n); 
  mpres_sub (u ,P->Z  ,P->T  ,n); 
  mpres_add (Q->Z ,v ,u ,n);
  mpres_sub (Q->T ,v ,u ,n); 
  
  mpres_mul (P->X ,Q->X ,R->X ,n); 
  mpres_mul (P->Y ,Q->Y ,R->Y ,n); 
  mpres_mul (P->Z ,Q->Z ,R->Z ,n); 
  mpres_mul (P->T ,Q->T ,R->T ,n); 
  mpres_muldivbysomething_si (P->X ,P->X ,cMul->x0p ,n);
  mpres_muldivbysomething_si (P->Y ,P->Y ,cMul->y0p ,n);
  mpres_muldivbysomething_si (P->Z ,P->Z ,cMul->z0p ,n);
  mpres_muldivbysomething_si (P->T ,P->T ,cMul->t0p ,n);
  
  
  //  hadamard (Q,P,u,v,n);
  mpres_add (u ,P->X  ,P->Y  ,n); 
  mpres_add (v ,P->Z  ,P->T  ,n); 
  mpres_add (Q->X ,u ,v ,n); 
  mpres_sub (Q->Y ,u ,v ,n); 
  mpres_sub (v ,P->X  ,P->Y  ,n);
  mpres_sub (u ,P->Z  ,P->T  ,n); 
  mpres_add (Q->Z ,v ,u ,n); 
  mpres_sub (Q->T ,v ,u ,n); 
  
  
  
  mpres_mul (P->X ,Q->X ,Q->X ,n); 
  mpres_mul (P->Y ,Q->Y ,Q->Y ,n); 
  mpres_mul (P->Z ,Q->Z ,Q->Z ,n); 
  mpres_mul (P->T ,Q->T ,Q->T ,n); 
  mpres_muldivbysomething_si (P->X ,P->X ,cMul->invX ,n);
  mpres_muldivbysomething_si (P->Y ,P->Y ,cMul->invY ,n);
  mpres_muldivbysomething_si (P->Z ,P->Z ,cMul->invZ ,n);
  mpres_muldivbysomething_si (P->T ,P->T ,cMul->invT ,n);
	
	
  mpres_mul (Q->X ,R->X ,R->X ,n); 
  mpres_mul (Q->Y ,R->Y ,R->Y ,n); 
  mpres_mul (Q->Z ,R->Z ,R->Z ,n); 
  mpres_mul (Q->T ,R->T ,R->T ,n); 
  mpres_muldivbysomething_si (Q->X ,Q->X ,cMul->x0p ,n);
  mpres_muldivbysomething_si (Q->Y ,Q->Y ,cMul->y0p ,n);
  mpres_muldivbysomething_si (Q->Z ,Q->Z ,cMul->z0p ,n);
  mpres_muldivbysomething_si (Q->T ,Q->T ,cMul->t0p ,n);
  
  //  hadamard (R,Q,u,v,n);
  mpres_add (u ,Q->X  ,Q->Y  ,n); 
  mpres_add (v ,Q->Z  ,Q->T  ,n); 
  mpres_add (R->X ,u ,v ,n);
  mpres_sub (R->Y ,u ,v ,n); 
  mpres_sub (v ,Q->X  ,Q->Y  ,n); 
  mpres_sub (u ,Q->Z  ,Q->T  ,n); 
  mpres_add (R->Z ,v ,u ,n); 
  mpres_sub (R->T ,v ,u ,n); 
  
  mpres_mul (Q->X ,R->X ,R->X ,n); 
  mpres_mul (Q->Y ,R->Y ,R->Y ,n); 
  mpres_mul (Q->Z ,R->Z ,R->Z ,n); 
  mpres_mul (Q->T ,R->T ,R->T ,n); 
  mpres_muldivbysomething_si (Q->X ,Q->X ,cMul->X0 ,n);
  mpres_muldivbysomething_si (Q->Y ,Q->Y ,cMul->Y0 ,n);
  mpres_muldivbysomething_si (Q->Z ,Q->Z ,cMul->Z0 ,n);
  mpres_muldivbysomething_si (Q->T ,Q->T ,cMul->T0 ,n);

}






// ***************** product of Hadamard *************




/*
  Compute the product of the Hadamard matrix with the vector P
  Put the result in P
  WARNING: don't do hadamard(P,P,...)
 */
void hadamard (ksPoint Pi,const ksPoint P,
	       mpres_t u,mpres_t v,
	       mpmod_t n) {


  mpres_add (u ,P->X  ,P->Y  ,n); // u = X + Y
  mpres_add (v ,P->Z  ,P->T  ,n); // v = Z + T

  mpres_add (Pi->X ,u ,v ,n); // Pi->X = u + v = (X + Y + Z + T)
  mpres_sub (Pi->Y ,u ,v ,n); // Pi->Y = u - v = (X + Y - Z - T)


  mpres_sub (v ,P->X  ,P->Y  ,n); // v = X - Y
  mpres_sub (u ,P->Z  ,P->T  ,n); // u = Z - T

  mpres_add (Pi->Z ,v ,u ,n); // Pi->Z = v + u = (X - Y + Z - T)
  mpres_sub (Pi->T ,v ,u ,n); // Pi->T = v - u = (X - Y - Z + T)


}
