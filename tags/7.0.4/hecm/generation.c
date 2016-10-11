#include <math.h>
#include <time.h>
#include "generation.h"
#include <limits.h>


void thetaCst_init (thetaCst th, mpmod_t n ) {
  mpres_init (th->be, n);
  mpres_init (th->ga, n);
  mpres_init (th->t5p, n);
  mpres_init (th->t6p, n);
  mpres_init (th->t7p, n);
  mpres_init (th->t10p, n);
  mpres_init (th->Rac, n);
  mpres_init (th->p, n);
}


void thetaCst_clear (thetaCst th, mpmod_t n ) {
  mpres_clear (th->be, n);
  mpres_clear (th->ga, n);
  mpres_clear (th->t5p, n);
  mpres_clear (th->t6p, n);
  mpres_clear (th->t7p, n);
  mpres_clear (th->t10p, n);
  mpres_clear (th->Rac, n);
  mpres_clear (th->p, n);
}


void curveHyperEll_init (curveHyperEll cHEll, mpmod_t n ) {
  mpres_init (cHEll->la, n);
  mpres_init (cHEll->mu, n);
  mpres_init (cHEll->nu, n);
  mpres_init (cHEll->q, n);
}

void curveHyperEll_clear (curveHyperEll cHEll, mpmod_t n ) {
  mpres_clear (cHEll->la, n);
  mpres_clear (cHEll->mu, n);
  mpres_clear (cHEll->nu, n);
  mpres_clear (cHEll->q, n);
}



void paraGenCurve_init (paraGenCurve para, mpmod_t n ) {
  mpres_init (para->s,n);
  mpz_init (para->a);
  mpz_init (para->b);
  mpres_init (para->x,n);
  mpres_init (para->y,n);
}

void paraGenCurve_clear (paraGenCurve para, mpmod_t n ) {
  mpres_clear (para->s,n);
  mpz_clear (para->a);
  mpz_clear (para->b);
  mpres_clear (para->x,n);
  mpres_clear (para->y,n);
}











// ********** normal parametrization *********




/*
  Generate a correct hyperelliptic curve on Z/nZ
  Begin with the given parameter and change it if needed
    if the curve is uncorrect (i.e. one constant is zero) then do Nextparam(s)
    the program finishes since s <= n
  After having choosen s we choose nJacobi. 
    TODO this is not a good idea. we should fixed s and increase nJacobi to keep
         the advantage of small s
    N_JACOBI_MIN <= nJacobi < N_JACOBI_MAX_p1 random
*/
int generateNormalCurve (mpz_t f,mpmod_t n,paraGenCurve para,
			  thetaCst th,curveHyperEll cHEll,
			  ksPoint P,ksCstPourMul cMul,
			  optionsHECM options) { 

  int test;
  mpres_t u,g;

  mpres_init (u,n);
  mpres_init (g,n);



  do {  

    test = generateOneNormalCurve (f,n,para,th,cHEll,P,cMul);
    mpres_set_z (g,f,n);

    if (test == GENERATION_A_CURVE) {
      // We have generate a curve
      // check if it is correct

      mpres_mul (u ,th->t7p,th->t6p,n);
      mpres_mul (u ,u ,P->X ,n);
      mpres_gcd (f ,u ,n);
      mpres_set_z (g ,f ,n);


      if ( !mpz_cmp_ui (f,1) ) {
	// general case, f=1, ie X,t6p,t7p non zero modulo all factors of n
	mpres_clear (u,n);
	mpres_clear (g,n);
	return GENERATION_CORRECT_CURVE;
      }
      else if (  mpres_is_zero(g,n) != 1  ) {
	// f != 0,n   i.e. f is a true divisor of n
	mpres_clear (u,n);
	mpres_clear (g,n);
	return GENERATION_FACTOR_FOUND; 
      }
      else { // X=0 mod n or t6p =0 mod n or t7p = 0 mod n
	nextParam (f,n,para,options); // Let's try again
      }
    }
    else { // We have not generated a curve
      if (  mpres_is_zero(g,n) == 1 ) {
	// We tried to divide by 0
	nextParam (f,n,para,options); // Let's try again
      }
      else {
	// f != 0,n   i.e. f is a true divisor of n
	// ( f!=1 because otherwise the inversion would have work)
	mpres_clear (u,n);
	mpres_clear (g,n);
	return GENERATION_FACTOR_FOUND; 
      }

    }
  } while ( 1 );


}



/*
  Check if a specified curve can be used
*/
int generateNormalCurveSpecified (mpz_t f,mpmod_t n,paraGenCurve para,
				   thetaCst th,curveHyperEll cHEll,
				   ksPoint P,ksCstPourMul cMul) { 

  int test;
  mpres_t g;

  mpres_init (g,n);
  
  
  test = generateOneNormalCurve (f,n,para,th,cHEll,P,cMul);
  mpres_set_z (g,f,n);
  
  if (test == GENERATION_A_CURVE) { // We have genereted a curve
    if ( mpres_is_zero(th->t7p,n) || mpres_is_zero(th->t6p,n)) {
      // One of the theta constants t7p or t6p is zero modulo n
      // Note that with the current choice of point on the Kummer surface we 
      //   compute 1/t6p so currently t6p != 0
      mpres_clear (g,n);
      return GENERATION_FAIL;
    }
    else {
      // The genereted curve is correct modulo n
      mpres_gcd (f, P->X, n); 
      if ( !mpz_cmp_ui (f,1) ) {
	// general case, f=1, i.e. X non zero modulo all the factors of n
	mpres_clear (g,n);
	return GENERATION_CORRECT_CURVE;
      }
      else if (  mpres_is_zero(g,n) != 1  ) {
	// f is a true divisor of n
	mpres_clear (g,n);
	return GENERATION_FACTOR_FOUND; 
      }
      else { // X=0 mod n
	mpres_clear (g,n);
	return GENERATION_FAIL;
      }
    }
  }
  else { // We have not generated a curve
    if (  mpres_is_zero(g,n) == 1 ) {
      // We tried to divide by 0
      mpres_clear (g,n);
      return GENERATION_FAIL;
    }
    else {
      mpres_clear (g,n);
      return GENERATION_FACTOR_FOUND; 
    }
  }

}













/*
  generate a curve with given s and nJacobi
*/
int generateOneNormalCurve (mpz_t f,mpmod_t n,paraGenCurve para,
			     thetaCst th,curveHyperEll cHEll,
			     ksPoint P,ksCstPourMul cMul) { 

  int test;
  mpres_t  u,x,ep,dep,v,Delta;

  mpres_init (u, n);
  mpres_init (ep, n);
  mpres_init (dep, n);
  mpres_init (x, n); 
  mpres_init (v, n); 
  mpres_init (Delta, n); 



  if (!mpres_invert (th->be, para->s, n)) // be=1/s
    {
      mpres_gcd (f, para->s, n);

      mpres_clear (u, n);
      mpres_clear (ep, n);
      mpres_clear (dep, n);
      mpres_clear (v, n);
      mpres_clear (x, n);
      mpres_clear (Delta, n);

      return GENERATION_NOT_CURVE;
    }
  // be=1/s 

  mpres_mul (ep ,th->be ,th->be ,n); // ep=be^2=1/s^2
  mpres_ui_sub (dep, 3, ep, n);
  mpres_mul (dep,dep,ep,n);  // dep = (3-1/s^2)/s^2 

  mpres_ui_sub (u ,1, ep, n); // u=1-1/s^2

  test = mulJacobi2 (f,n,para->nJacobi,para->x,para->y,u,ep,dep);
  // get x and y on the Jacobi curve Y^2 = ep X^4 - dep X^2 + 1 with (1,u) as 
  //   as initial point 

  if (test==MULT_JACOBI_FAIL)
    {
      mpres_clear (u, n);
      mpres_clear (ep, n);
      mpres_clear (dep, n);
      mpres_clear (v, n);
      mpres_clear (x, n);
      mpres_clear (Delta, n);

      return GENERATION_NOT_CURVE;
   }


  mpres_mul (cHEll->q, para->s ,para->s ,n); // q=s^2
  mpres_mul (v, para->x, para->x, n); // v=x^2

  mpres_sub(cHEll->nu, v, cHEll->q, n); // x^2-s^2


  mpres_sub_ui (v ,v ,1 ,n); // v=x^2-1
  if (!mpres_invert (v, v, n)) // v=1/(x^2-1)
    {
      mpres_gcd (f, v, n);

      mpres_clear (u, n);
      mpres_clear (ep, n);
      mpres_clear (dep, n);
      mpres_clear (v, n);
      mpres_clear (x, n);
      mpres_clear (Delta, n);

      return GENERATION_NOT_CURVE;
    }

  mpres_mul (cHEll->nu ,cHEll->nu ,v, n); // nu = (x^2-s^2)/(x^2-1)  ok


  mpres_mul (cHEll->mu ,cHEll->nu ,cHEll->nu ,n);   // nu^2
  mpres_sub (Delta ,cHEll->q ,cHEll->mu ,n); // Delta=s^2-nu^2
  mpres_add (cHEll->mu ,cHEll->mu ,cHEll->q ,n);    // nu^2+s^2
  mpres_sub (cHEll->mu ,cHEll->mu ,cHEll->nu ,n);   // mu=-nu+nu^2+s^2
  mpres_mul (cHEll->mu ,cHEll->mu ,ep ,n);   // mu = (-nu+nu^2+s^2)/s^2


  mpres_sub_ui(cHEll->q ,cHEll->q ,1 ,n); // q=s^2-1
  mpres_mul (v ,v ,v ,n);   // v=1/(x^2-1)^2
  mpres_mul (cHEll->q ,cHEll->q ,v ,n);
  mpres_mul (cHEll->q ,cHEll->q ,para->y ,n); 
  mpres_mul (cHEll->q ,cHEll->q ,para->x ,n); // q = x*y*(s^2-1)/(x^2-1)^2

  // We finished be, nu, mu, q







  mpres_sub_ui (cHEll->la ,cHEll->nu ,1 ,n);
  mpres_mul (cHEll->la ,cHEll->la ,cHEll->mu, n); // mu*(nu-1)
  mpres_sub_ui (u ,cHEll->mu ,1 ,n); // u=mu-1

  if (!mpres_invert (u, u, n)) // u=1/(mu-1)
    {
      mpres_gcd (f, u, n);

      mpres_clear (u, n);
      mpres_clear (ep, n);
      mpres_clear (dep, n);
      mpres_clear (v, n);
      mpres_clear (x, n);
      mpres_clear (Delta, n);

      return GENERATION_NOT_CURVE;
    }
  mpres_mul (cHEll->la ,cHEll->la ,u ,n); // la=mu*(nu-1)/(mu-1)  ok
  // We finished la






  mpres_mul (th->ga ,cHEll->la ,th->be ,n); // ga= la*be

  mpres_sub (th->t7p ,th->ga ,th->be ,n); // t7p=ga-be   
  mpres_mul (th->Rac ,th->t7p ,cHEll->nu ,n);
  mpres_mul (th->Rac ,th->Rac ,th->be ,n);
  mpres_neg (th->Rac ,th->Rac ,n); // Rac=nu*be*(be-ga) // Rac = B/t10p

  // We finished ga, t7p et Rac



  if (!mpres_invert (u, cHEll->mu, n)) // u=1/mu
    {
      mpres_gcd (f, cHEll->mu, n);

      mpres_clear (u, n);
      mpres_clear (ep, n);
      mpres_clear (dep, n);
      mpres_clear (v, n);
      mpres_clear (x, n);
      mpres_clear (Delta, n);

      return GENERATION_NOT_CURVE;
    }


  mpres_mul (cMul->Z0 ,cHEll->nu ,u ,n);
  mpres_mul (cMul->Z0 ,cMul->Z0 ,th->be ,n); // Z0=nu/(mu*s)  // Z0=1/ga 

  mpres_mul (th->t10p ,cHEll->la ,u ,n);
  mpres_mul (th->t10p ,th->t10p ,th->be ,n); // t10p=la/(mu*s) // t10p=1/(be*nu)

  mpres_sub (th->t5p ,th->ga ,th->t10p ,n); // t5p=ga-t10p
  mpres_sub (th->t6p ,th->be ,th->t10p ,n); // t6p=be-t10p

  // We finished Z0, t10p, t5p et t6p



  mpres_add (cMul->x0p ,th->be ,th->ga ,n);
  mpres_ui_sub (cMul->t0p ,2 ,cMul->x0p ,n); // t0p= 2-be-ga = D 
  mpres_add_ui (cMul->x0p ,cMul->x0p ,2 ,n); // x0p= 2+be+ga = A 

  if (!mpres_invert (cMul->x0p,cMul->x0p, n)) // x0p=1/A
    {
      mpres_gcd (f, cMul->x0p, n);

      mpres_clear (u, n);
      mpres_clear (ep, n);
      mpres_clear (dep, n);
      mpres_clear (v, n);
      mpres_clear (x, n);
      mpres_clear (Delta, n);

      return GENERATION_NOT_CURVE;
    }

  mpres_mul (dep,th->be,th->ga,n); // dep=be*ga

  mpres_add_ui (u ,th->be ,1 ,n);
  mpres_mul (Delta ,Delta ,u ,n);
  mpres_add_ui (u ,th->ga ,1 ,n);
  mpres_mul (Delta ,Delta ,u ,n);
  mpres_mul (Delta ,Delta ,th->t10p ,n); 
  mpres_mul (Delta ,Delta ,ep, n);
 // We finished Delta

  mpres_mul (cMul->invT ,Delta ,cMul->x0p ,n); // Delta/A

  mpres_sub_ui (u ,th->be ,1 ,n);
  mpres_sub_ui (v ,th->ga ,1 ,n);
  mpres_mul (u ,u ,v ,n); // (be-1)*(ga-1)
  mpres_add_ui (v ,dep ,1 ,n); // be*ga+1
  mpres_mul (u ,u ,v ,n); // u=(be-1)*(ga-1)*(be*ga+1)

  mpres_mul_ui (v ,dep ,2 ,n);
  mpres_mul (v ,v ,cMul->t0p ,n); // v= 2*be*ga*D
  mpres_sub (u ,u ,v ,n); // u=(be-1)*(ga-1)*(be*ga+1) - 2*be*ga*D

  mpres_ui_sub (v ,1 ,dep ,n); // v= 1-be*ga
  mpres_mul (cMul->invT ,cMul->invT ,v ,n); // Delta/A*(1-be*ga)
  mpres_sub (cMul->invT ,u ,cMul->invT ,n);
  // invT = (be-1)*(ga-1)*(be*ga+1) - 2*be*ga*D - Delta/A*(1-be*ga) 

  mpres_mul_ui (v ,v ,2 ,n);
  if (!mpres_invert (v,v, n)) // x0p=1/(1-be*ga)/2
    {
      mpres_gcd (f, v, n);

      mpres_clear (u, n);
      mpres_clear (ep, n);
      mpres_clear (dep, n);
      mpres_clear (v, n);
      mpres_clear (x, n);
      mpres_clear (Delta, n);

      return GENERATION_NOT_CURVE;
    }

  mpres_mul (cMul->invT ,v ,cMul->invT, n);
  // invT = ( (be-1)*(ga-1)*(be*ga+1)+2*be*ga*D-Delta/A*(1-be*ga) )/2/(1-be*ga) 
  // We finished invT !






  mpres_mul (cMul->x0p ,cMul->x0p ,th->t7p, n);
  mpres_ui_sub (cMul->x0p ,0 ,cMul->x0p ,n); // x0p= B/A

  // We finished x0p


  if (!mpres_invert (cMul->t0p,cMul->t0p, n)) // t0p=1/D
    {
      mpres_gcd (f, cMul->t0p, n);

      mpres_clear (u, n);
      mpres_clear (ep, n);
      mpres_clear (dep, n);
      mpres_clear (v, n);
      mpres_clear (x, n);
      mpres_clear (Delta, n);

      return GENERATION_NOT_CURVE;
    }

  mpres_mul (cMul->t0p ,cMul->t0p ,th->t7p, n);
  mpres_ui_sub (cMul->t0p ,0 ,cMul->t0p ,n); // t0p= B/D

  // We finished t0p





  mpres_mul (P->T ,para->s ,cMul->Z0 ,n);
  mpres_ui_sub (P->T, 0, P->T ,n); // T=-sZ0=-1/be/ga
  mpres_mul (P->X, cMul->invT ,P->T ,n);
  mpres_neg (P->Y ,P->X ,n);


  mpres_set_ui (P->Z, 1, n);

  mpres_mul (th->p, cMul->x0p ,cMul->t0p ,n); // p=B^2/(A*D)


  mpres_set (cMul->Y0,para->s,n); // Y0 = s
  mpres_set (cMul->invZ,P->X,n); // invZ = X


  mpres_clear (x, n);
  mpres_clear (v, n);
  mpres_clear (u, n);
  mpres_clear (ep, n);
  mpres_clear (dep, n);
  mpres_clear (Delta, n);


  return GENERATION_A_CURVE;
}



















// ********** small parameters *********






int inverseCoorPointKS (mpz_t f,mpmod_t n,
			ksPoint P,ksSmallConstPourMul cMul) {

  mpres_set_si (P->X ,cMul->invX ,n);
  if (!mpres_invert (P->X ,P->X ,n)) {
    mpres_gcd (f ,P->X ,n);
    return INVERSE_NOT_COOR;
  }
  


  mpres_set_si (P->Y ,cMul->invY ,n);
  if (!mpres_invert (P->Y ,P->Y ,n)) {
    mpres_gcd (f ,P->Y ,n);
    return INVERSE_NOT_COOR;
  }


  mpres_set_si (P->Z ,cMul->invZ ,n);
  if (!mpres_invert (P->Z ,P->Z ,n)) {
    mpres_gcd (f ,P->Z ,n);
    return INVERSE_NOT_COOR;
  }
  


  mpres_set_si (P->T ,cMul->invT ,n);
  if (!mpres_invert (P->T ,P->T ,n)) {
    mpres_gcd (f ,P->T ,n);
    return INVERSE_NOT_COOR;
  }


  return INVERSE_COOR;
}





/*
  Generate an hyperelliptic curve with small parameters
  They are imposed by the choice of a,b and nJacobi
  Work with mpz_t even if the goal is that the parameters fit in  long

  We begin with creating X0,Y0,.. and x0p,y0p,...
  Then we check if they are small (i.e. if they fit in long)
  We construct the coordinate of a point on the Kummer surface
  Then we check if their inverses are small (i.e. if they fit in long)
  Now we have all small parameters and we construct the other
*/
int generateOneCurveSmallParam (mpz_t f,mpmod_t n,
				paraGenCurve para ,thetaCst th,
				curveHyperEll cHEll,ksPoint P,
				ksSmallConstPourMul cMul,
				mpz_t a, mpz_t b,
				mpz_t x, mpz_t y, mpz_t z) {

  int test;

  mpz_t thCste[3]; // list to check if the "theta constants" fit in long
  mpz_t ThCste[3]; // list to check if the "Theta constants" fit in long
  mpz_t ptCste[3]; // list to check if the "point" on KS fit in long
  mpz_t g; 
  mpz_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t0;
  mpz_t sa,sb,sx,sz;

  mpz_init (thCste[0]);
  mpz_init (thCste[1]);
  mpz_init (thCste[2]);
  mpz_init (ThCste[0]);
  mpz_init (ThCste[1]);
  mpz_init (ThCste[2]);
  mpz_init (ptCste[2]);
  mpz_init (ptCste[0]);
  mpz_init (ptCste[1]);
  mpz_init (g);

  mpz_init (t1);
  mpz_init (t2);
  mpz_init (t3);
  mpz_init (t4);
  mpz_init (t5);
  mpz_init (t6);
  mpz_init (t7);
  mpz_init (t8);
  mpz_init (t9);
  mpz_init (t0);

  mpz_init (sa);
  mpz_init (sb);
  mpz_init (sx);
  mpz_init (sz);



  mpz_mul (sa ,a ,a); // sa = a^2
  mpz_mul (sb ,b ,b); // sb = b^2
  mpz_mul (sx ,x ,x); // sx = x^2
  mpz_mul (sz ,z ,z); // sz = z^2



  mpz_mul (t8 ,a ,sz);
  mpz_mul (t7 ,t8 ,t8); // a^2*z^4
  mpz_mul (t9 ,b ,sx);
  mpz_add (t3 ,t8 ,t9);
  mpz_sub (t0 ,t8 ,t9); // t0 =   a*z^2 - b*x^2 
  mpz_mul (t3 ,t3 ,t3); // t3 = ( a*z^2 + b*x^2 )^2
  mpz_mul (t4 ,t0 ,t0); // t4 = ( a*z^2 - b*x^2 )^2

  mpz_mul (t8 ,sa ,sz);
  mpz_mul (t1 ,t8 ,t8);
  mpz_mul (t9 ,sb ,sx);
  mpz_sub (t5 ,t8 ,t9); // t5 = a^2*z^2 - b^2*x^2

  mpz_mul_ui (t6 ,sz ,2);
  mpz_sub (t6 ,sx ,t6);
  mpz_mul (t6 ,t6 ,t9); // x^2*b^2*(x^2-2*z^2)
  mpz_add (t6 ,t6 ,t7); // t6 = a^2*z^4 - 2*b^2*x^2*z^2 + b^2*x^4

  mpz_mul (t7 ,sa ,y);
  mpz_mul (t7 ,t7 ,t7);
  // t7 = b^4*x^2*z^2 - 3*a^2*b^2*x^2*z^2 + a^4*z^4 + a^2*b^2*x^4 = y^2*a^4


  mpz_mul (t1 ,a ,b); // t1 = a*b
  mpz_sub (t2 ,sx ,sz); // t2 = x^2 - z^2



  mpz_mul (thCste[0] ,b ,t7);
  // thCste[0] = b*(b^4*x^2*z^2 - 3*a^2*b^2*x^2*z^2 + a^4*z^4 + a^2*b^2*x^4)
  mpz_mul (thCste[1] ,a ,t7);
  // thCste[1] = a*(b^4*x^2*z^2 - 3*a^2*b^2*x^2*z^2 + a^4*z^4 + a^2*b^2*x^4)
  mpz_mul (thCste[2] ,a ,sb);
  mpz_mul (thCste[2] ,thCste[2] ,t5);
  mpz_mul (thCste[2] ,thCste[2] ,t2);
  mpz_neg (thCste[2] ,thCste[2]);
  // thCste[2] = - a * b^2 * (x^2-z^2) * (a^2*z^2 - b^2*x^2)

  mpz_sub (t8 ,a ,b);
  mpz_mul (ThCste[0] ,t8 ,t3);
  mpz_add (t9 ,a ,b);
  mpz_mul (ThCste[2] ,t9 ,t4);
  mpz_mul (ThCste[1] ,ThCste[0] ,ThCste[2]);
  // ThCste[1] = (a^2-b^2) * (a*z^2+b*x^2)^2 * (a*z^2-b*x^2)^2
  mpz_mul (ThCste[0] ,ThCste[0] ,t6);
  mpz_mul (ThCste[0] ,ThCste[0] ,t8);
  mpz_neg (ThCste[0] ,ThCste[0]);
  // ThCste[0] = - (a-b)^2 * (a*z^2+b*x^2)^2 * (a^4*z^4-2*b^2*x^2*z^2+b^4*z^4)
  mpz_mul (ThCste[2] ,ThCste[2] ,t6);
  mpz_mul (ThCste[2] ,ThCste[2] ,t9);
  // ThCste[2] = (a+b)^2 * (a*z^2-b*x^2)^2 * (a^4*z^4-2*b^2*x^2*z^2+b^4*z^4)




  mpz_gcd (g ,thCste[0] ,thCste[1]);
  mpz_gcd (g ,g ,thCste[2]);
  if ( mpz_sgn(g) == 0 ) { // the constants are zero

    mpz_clear (thCste[0]);
    mpz_clear (thCste[1]);
    mpz_clear (thCste[2]);
    mpz_clear (ThCste[0]);
    mpz_clear (ThCste[1]);
    mpz_clear (ThCste[2]);
    mpz_clear (ptCste[0]);
    mpz_clear (ptCste[1]);
    mpz_clear (ptCste[2]);
    mpz_clear (g);
    
    mpz_clear (t1);
    mpz_clear (t2);
    mpz_clear (t3);
    mpz_clear (t4);
    mpz_clear (t5);
    mpz_clear (t6);
    mpz_clear (t7);
    mpz_clear (t8);
    mpz_clear (t9);
    mpz_clear (t0);
    
    mpz_clear (sa);
    mpz_clear (sb);
    mpz_clear (sx);
    mpz_clear (sz);

    mpz_set_ui (f ,0);
    return GENERATION_NOT_CURVE;
  }
  mpz_divexact (thCste[0] ,thCste[0] ,g);
  mpz_divexact (thCste[1] ,thCste[1] ,g);
  mpz_divexact (thCste[2] ,thCste[2] ,g);



  mpz_gcd (g ,ThCste[0] ,ThCste[1]);
  mpz_gcd (g ,g ,ThCste[2]);
  if ( mpz_sgn(g) == 0 ) { // the constants are zero

    mpz_clear (thCste[0]);
    mpz_clear (thCste[1]);
    mpz_clear (thCste[2]);
    mpz_clear (ThCste[0]);
    mpz_clear (ThCste[1]);
    mpz_clear (ThCste[2]);
    mpz_clear (ptCste[0]);
    mpz_clear (ptCste[1]);
    mpz_clear (ptCste[2]);
    mpz_clear (g);
    
    mpz_clear (t1);
    mpz_clear (t2);
    mpz_clear (t3);
    mpz_clear (t4);
    mpz_clear (t5);
    mpz_clear (t6);
    mpz_clear (t7);
    mpz_clear (t8);
    mpz_clear (t9);
    mpz_clear (t0);
    
    mpz_clear (sa);
    mpz_clear (sb);
    mpz_clear (sx);
    mpz_clear (sz);

    mpz_set_ui (f ,0);
    return GENERATION_NOT_CURVE;
  }
  mpz_divexact (ThCste[0] ,ThCste[0] ,g);
  mpz_divexact (ThCste[1] ,ThCste[1] ,g);
  mpz_divexact (ThCste[2] ,ThCste[2] ,g);




  if ( ( mpz_cmpabs_ui (thCste[0] ,LONG_MAX ) >0 ) || ( mpz_cmpabs_ui (thCste[1] ,LONG_MAX ) >0 ) || ( mpz_cmpabs_ui (thCste[2] ,LONG_MAX ) >0 ) || ( mpz_cmpabs_ui (ThCste[0] ,LONG_MAX ) >0 ) || ( mpz_cmpabs_ui (ThCste[1] ,LONG_MAX ) >0 ) || ( mpz_cmpabs_ui (ThCste[2] ,LONG_MAX ) >0 ) ) {

    mpz_clear (thCste[0]);
    mpz_clear (thCste[1]);
    mpz_clear (thCste[2]);
    mpz_clear (ThCste[0]);
    mpz_clear (ThCste[1]);
    mpz_clear (ThCste[2]);
    mpz_clear (ptCste[0]);
    mpz_clear (ptCste[1]);
    mpz_clear (ptCste[2]);
    mpz_clear (g);
    
    mpz_clear (t1);
    mpz_clear (t2);
    mpz_clear (t3);
    mpz_clear (t4);
    mpz_clear (t5);
    mpz_clear (t6);
    mpz_clear (t7);
    mpz_clear (t8);
    mpz_clear (t9);
    mpz_clear (t0);
    
    mpz_clear (sa);
    mpz_clear (sb);
    mpz_clear (sx);
    mpz_clear (sz);

    return GENERATION_THESE_SMALL_PARAM_TOO_BIG;
  }

  // Now we know that thCste and ThCste fits in long

    cMul->X0 = mpz_get_si (thCste[0]);
    cMul->Y0 = mpz_get_si (thCste[1]);
    cMul->Z0 = mpz_get_si (thCste[2]);
    cMul->T0 = mpz_get_si (thCste[0]);


    cMul->x0p = mpz_get_si (ThCste[0]);
    cMul->y0p = mpz_get_si (ThCste[1]);
    cMul->z0p = - mpz_get_si (ThCste[1]);
    cMul->t0p = mpz_get_si (ThCste[2]);








  // Construction of a point on the Kummer surface



  mpz_mul (ptCste[0] ,t0 ,sz);
  mpz_mul (ptCste[0] ,ptCste[0] ,t5);
  mpz_mul (ptCste[0] ,ptCste[0] ,t7);
  mpz_neg (ptCste[0] ,ptCste[0]);
  // ptCste[0] = - z^2 * (a*z^2 - b*x^2) * (a^2*z^2 - b^2*x^2) * (b^4*x^2*z^2 - 3*a^2*b^2*x^2*z^2 + a^4*z^4 + a^2*b^2*x^4)

  mpz_mul (t8 ,sz ,sz);
  mpz_mul (t9 ,t8 ,a); // z^4*a
  mpz_mul (t8 ,a ,b);
  mpz_sub (t8 ,sa ,t8);
  mpz_add (t8 ,t8 ,sb); // a^2 - a*b + b^2
  mpz_mul (ptCste[1] ,t8 ,t9); // z^4 * a * (a^2 - a*b + b^2)
  mpz_mul_ui (ptCste[2] ,a ,3);
  mpz_sub (ptCste[2] ,b ,ptCste[2]); // b-3a
  mpz_mul (t8 ,ptCste[2] ,sz); // z^2 * (b-3*a)
  mpz_mul (t9 ,a ,sx);
  mpz_add (t8 ,t8 ,t9); //  z^2 * (b-3*a)  +  a*x^2
  mpz_mul (t9 ,sx ,sb);
  mpz_mul (ptCste[2] ,t8 ,t9); //  ( z^2 * (b-3*a)  +  a*x^2 ) * x^2*b^2
  mpz_add (ptCste[2] ,ptCste[2] ,ptCste[1]);
  // z^4 * a * (a^2 - a*b + b^2)  +   ( z^2 * (b-3*a)  +  a*x^2 ) * x^2*b^2
  // a*b^2*x^4 - 3*a*b^2*x^2*z^2 + b^3*x^2*z^2 + a^3*z^4 - a^2*b*z^4 + a*b^2*z^4

  mpz_mul (ptCste[1] ,ptCste[2] ,sx);
  mpz_mul (ptCste[1] ,ptCste[1] ,sb);
  mpz_mul (ptCste[1] ,ptCste[1] ,t2);
  mpz_mul (ptCste[1] ,ptCste[1] ,t5);
  mpz_neg (ptCste[1] ,ptCste[1]);
  // ptCste[1] = - x^2 * (x^2-z^2) * b^2 * (a^2*z^2 - b^2*x^2) * (a*b^2*x^4-3*a*b^2*x^2*z^2+b^3*x^2*z^2+a^3*z^4-a^2*b*z^4+a*b^2*z^4)

  mpz_mul (ptCste[2] ,ptCste[2] ,sx);
  mpz_mul (ptCste[2] ,ptCste[2] ,t7);
  mpz_neg (ptCste[2] ,ptCste[2]);
  // ptCste[2] = - x^2 * (b^4*x^2*z^2 - 3*a^2*b^2*x^2*z^2 + a^4*z^4 + a^2*b^2*x^4) * (a*b^2*x^4-3*a*b^2*x^2*z^2+b^3*x^2*z^2+a^3*z^4-a^2*b*z^4+a*b^2*z^4)



  mpz_gcd (g ,ptCste[0] ,ptCste[1]);
  mpz_gcd (g ,g ,ptCste[2]);
  if ( mpz_sgn(g) != 0 ) { // g != 0.
    mpz_divexact (ptCste[0] ,ptCste[0] ,g);
    mpz_divexact (ptCste[1] ,ptCste[1] ,g);
    mpz_divexact (ptCste[2] ,ptCste[2] ,g);
  }

  if ( ( mpz_sgn(g) != 0 ) &&  ( mpz_cmpabs_ui (ptCste[0] ,LONG_MAX ) <= 0 ) && ( mpz_cmpabs_ui (ptCste[1] ,LONG_MAX ) <= 0 ) && ( mpz_cmpabs_ui (ptCste[2] ,LONG_MAX ) <= 0 ) ) { 
    // We have genereted a correct point on the Kummer surface


    cMul->invX = mpz_get_si (ptCste[0]);
    cMul->invY = - mpz_get_si (ptCste[0]);
    cMul->invZ = mpz_get_si (ptCste[1]);
    cMul->invT = mpz_get_si (ptCste[2]);


    test = inverseCoorPointKS (f ,n ,P ,cMul);

  }
  else { // g=0 or the constants are too large

      
    mpz_mul (ptCste[0] ,t0 ,a);
    mpz_mul (ptCste[0] ,ptCste[0] ,t5);
    mpz_mul (ptCste[0] ,ptCste[0] ,t7);
    // ptCste[0] = a * (a*z^2 - b*x^2) * (a^2*z^2 - b^2*x^2) * (b^4*x^2*z^2 - 3*a^2*b^2*x^2*z^2 + a^4*z^4 + a^2*b^2*x^4)
    
    mpz_mul (t8 ,sx ,sz);
    mpz_mul (t9 ,t8 ,b); // x^2*z^2*b
    mpz_mul_ui (ptCste[1] ,sa ,2);
    mpz_mul (t8 ,a ,b);
    mpz_add (t8 ,t8 ,ptCste[1]);
    mpz_sub (t8 ,t8 ,sb); // 2*a^2 + a*b - b^2
    mpz_mul (ptCste[1] ,t9 ,t8); // x^2*z^2*b * (2*a^2 + a*b - b^2)
    mpz_mul (t8 ,sz ,sz);
    mpz_mul (t9 ,t8 ,a);
    mpz_mul (ptCste[2] ,sx ,sx);
    mpz_mul (t8 ,ptCste[2] ,b);
    mpz_add (t8 ,t8 ,t9); // a*z^4 + b*x^4
    mpz_mul (ptCste[2] ,t8 ,sa); // a^2 * (a*z^4 + b*x^4)
    mpz_sub (ptCste[1] ,ptCste[1] ,ptCste[2]);
    // x^2*z^2*b * (2*a^2 + a*b - b^2) - a^2 * (a*z^4 + b*x^4)
    // - ( a^2*b*x^4 - 2*a^2*b*x^2*z^2 - a*b^2*x^2*z^2 + b^3*x^2*z^2 + a^3*z^4 )
    mpz_mul (ptCste[2] ,ptCste[1] ,b);
    //-b*( a^2*b*x^4 - 2*a^2*b*x^2*z^2 - a*b^2*x^2*z^2 + b^3*x^2*z^2 + a^3*z^4 )
    
    mpz_mul (ptCste[1] ,ptCste[2] ,sa);
    mpz_mul (ptCste[1] ,ptCste[1] ,t2);
    mpz_mul (ptCste[1] ,ptCste[1] ,t5);
    // ptCste[1] = - a^2 * b * (x^2-z^2) * (a^2*z^2 - b^2*x^2) * (a^2*b*x^4 - 2*a^2*b*x^2*z^2 - a*b^2*x^2*z^2 + b^3*x^2*z^2 + a^3*z^4)
    
    mpz_mul (ptCste[2] ,ptCste[2] ,t7);
    // ptCste[2] = - b * (a^2*z^2 - b^2*x^2) * (a^2*b*x^4 - 2*a^2*b*x^2*z^2 - a*b^2*x^2*z^2 + b^3*x^2*z^2 + a^3*z^4) * (b^4*x^2*z^2 - 3*a^2*b^2*x^2*z^2 + a^4*z^4 + a^2*b^2*x^4)
    
    
    mpz_gcd (g ,ptCste[0] ,ptCste[1]);
    mpz_gcd (g ,g ,ptCste[2]);
    if ( mpz_sgn(g) == 0 ) { // g = 0.
      
      mpz_clear (thCste[0]);
      mpz_clear (thCste[1]);
      mpz_clear (thCste[2]);
      mpz_clear (ThCste[0]);
      mpz_clear (ThCste[1]);
      mpz_clear (ThCste[2]);
      mpz_clear (ptCste[0]);
      mpz_clear (ptCste[1]);
      mpz_clear (ptCste[2]);
      mpz_clear (g);
      
      mpz_clear (t1);
      mpz_clear (t2);
      mpz_clear (t3);
      mpz_clear (t4);
      mpz_clear (t5);
      mpz_clear (t6);
      mpz_clear (t7);
      mpz_clear (t8);
      mpz_clear (t9);
      mpz_clear (t0);
      
      mpz_clear (sa);
      mpz_clear (sb);
      mpz_clear (sx);
      mpz_clear (sz);

      mpz_set_ui (f ,0);
      return GENERATION_NOT_CURVE;
      /*
	TODO
	We could try with other points on the Kummer surface
      */
    }
    mpz_divexact (ptCste[0] ,ptCste[0] ,g);
    mpz_divexact (ptCste[1] ,ptCste[1] ,g);
    mpz_divexact (ptCste[2] ,ptCste[2] ,g);
    
    if ( ( mpz_cmpabs_ui (ptCste[0] ,LONG_MAX ) <= 0 ) && ( mpz_cmpabs_ui (ptCste[1] ,LONG_MAX ) <= 0 ) && ( mpz_cmpabs_ui (ptCste[2] ,LONG_MAX ) <= 0 ) ) {
      // We have genereted a good point on the Kummer surface


      cMul->invX = mpz_get_si (ptCste[0]);
      cMul->invY = - mpz_get_si (ptCste[0]);
      cMul->invZ = mpz_get_si (ptCste[1]);
      cMul->invT = mpz_get_si (ptCste[2]);
  
      
      test = inverseCoorPointKS (f ,n ,P ,cMul);

    }
    else { // the constants are too big
      
      mpz_clear (thCste[0]);
      mpz_clear (thCste[1]);
      mpz_clear (thCste[2]);
      mpz_clear (ThCste[0]);
      mpz_clear (ThCste[1]);
      mpz_clear (ThCste[2]);
      mpz_clear (ptCste[0]);
      mpz_clear (ptCste[1]);
      mpz_clear (ptCste[2]);
      mpz_clear (g);
      
      mpz_clear (t1);
      mpz_clear (t2);
      mpz_clear (t3);
      mpz_clear (t4);
      mpz_clear (t5);
      mpz_clear (t6);
      mpz_clear (t7);
      mpz_clear (t8);
      mpz_clear (t9);
      mpz_clear (t0);
      
      mpz_clear (sa);
      mpz_clear (sb);
      mpz_clear (sx);
      mpz_clear (sz);
      
      return GENERATION_THESE_SMALL_PARAM_TOO_BIG;
      /*
	TODO
	We could try with other points on the Kummer surface
      */
    }
  }


  if (test == 0) {
    // We had a problem during the inversion of the coordinates of the initial
    // point

    mpz_clear (thCste[0]);
    mpz_clear (thCste[1]);
    mpz_clear (thCste[2]);
    mpz_clear (ThCste[0]);
    mpz_clear (ThCste[1]);
    mpz_clear (ThCste[2]);
    mpz_clear (ptCste[0]);
    mpz_clear (ptCste[1]);
    mpz_clear (ptCste[2]);
    mpz_clear (g);
    
    mpz_clear (t1);
    mpz_clear (t2);
    mpz_clear (t3);
    mpz_clear (t4);
    mpz_clear (t5);
    mpz_clear (t6);
    mpz_clear (t7);
    mpz_clear (t8);
    mpz_clear (t9);
    mpz_clear (t0);
    
    mpz_clear (sa);
    mpz_clear (sb);
    mpz_clear (sx);
    mpz_clear (sz);
    
    return GENERATION_NOT_CURVE;
  }




  mpres_t T1,T2,T3,T4,T5,T6,T7;
  mpres_t amod,bmod;
  mpres_t u,v;

  mpres_init (T1 ,n);
  mpres_init (T2 ,n);
  mpres_init (T3 ,n);
  mpres_init (T4 ,n);
  mpres_init (T5 ,n);
  mpres_init (T6 ,n);
  mpres_init (T7 ,n);

  mpres_init (amod ,n);
  mpres_init (bmod ,n);
  mpres_init (u ,n);
  mpres_init (v ,n);

  mpres_set_z (T1 ,t1 ,n); // T1 = a*b
  mpres_set_z (T2 ,t2 ,n); // T2 = x^2 - z^2
  mpres_set_z (T3 ,t3 ,n); // T3 = ( a*z^2 + b*x^2 )^2
  mpres_set_z (T4 ,t4 ,n); // T4 = ( a*z^2 - b*x^2 )^2
  mpres_set_z (T5 ,t5 ,n); // T5 = a^2*z^2 - b^2*x^2
  mpres_set_z (T6 ,t6 ,n); // T6 = a^2*z^4 - 2*b^2*x^2*z^2 + b^2*x^4
  mpres_set_z (T7 ,t7 ,n);
  // T7 = b^4*x^2*z^2 - 3*a^2*b^2*x^2*z^2 + a^4*z^4 + a^2*b^2*x^4


  mpz_clear (thCste[0]);
  mpz_clear (thCste[1]);
  mpz_clear (thCste[2]);
  mpz_clear (ThCste[0]);
  mpz_clear (ThCste[1]);
  mpz_clear (ThCste[2]);
  mpz_clear (ptCste[0]);
  mpz_clear (ptCste[1]);
  mpz_clear (ptCste[2]);
  mpz_clear (g);

  mpz_clear (t1);
  mpz_clear (t2);
  mpz_clear (t3);
  mpz_clear (t4);
  mpz_clear (t5);
  mpz_clear (t6);
  mpz_clear (t7);
  mpz_clear (t8);
  mpz_clear (t9);
  mpz_clear (t0);
        
  mpz_clear (sa);
  mpz_clear (sb);
  mpz_clear (sx);
  mpz_clear (sz);


  mpres_set_z (para->y ,y ,n);
  mpres_set_z (para->x ,x ,n);
  mpres_set_z (u ,z ,n);
  if (!mpres_invert (u, u, n))  // u = 1/z
    {
      mpres_gcd (f, u, n);

      mpres_clear (T1 ,n);
      mpres_clear (T2 ,n);
      mpres_clear (T3 ,n);
      mpres_clear (T4 ,n);
      mpres_clear (T5 ,n);
      mpres_clear (T6 ,n);
      mpres_clear (T7 ,n);
      
      mpres_clear (amod ,n);
      mpres_clear (bmod ,n);
      mpres_clear (u ,n);
      mpres_clear (v ,n);

      return GENERATION_NOT_CURVE;
    }
  mpres_mul (para->x ,para->x ,u ,n); // x=x/z
  mpres_mul (u ,u ,u ,n); // 1/z^2
  mpres_mul (para->y ,para->y ,u ,n);

  mpres_set_z (amod ,a ,n);
  mpres_set_z (bmod ,b ,n);

  if (!mpres_invert (u, T1, n))  // u = 1/(a*b)
    {
      mpres_gcd (f, T1, n);

      mpres_clear (T1 ,n);
      mpres_clear (T2 ,n);
      mpres_clear (T3 ,n);
      mpres_clear (T4 ,n);
      mpres_clear (T5 ,n);
      mpres_clear (T6 ,n);
      mpres_clear (T7 ,n);
      
      mpres_clear (amod ,n);
      mpres_clear (bmod ,n);
      mpres_clear (u ,n);
      mpres_clear (v ,n);

      return GENERATION_NOT_CURVE;
    }
  mpres_mul (para->s ,amod ,u ,n);
  mpres_mul (para->s ,para->s ,amod ,n); // s=a/b
  mpres_mul (th->be ,bmod ,u ,n);
  mpres_mul (th->be ,th->be ,bmod ,n); // be = 1/s = b/a
  // On a fini s et be

  mpres_mul (th->t10p ,T1 ,T2 ,n); // a*b*(x^2-z^2)
  mpres_mul (th->ga ,u ,T7 ,n); // t7/(a*b)

  mpres_mul (u ,u ,u ,n); // 1/(a*b)^2

  if (!mpres_invert (cHEll->q ,T2 ,n))  // 1/(x^2-z^2)
    {
      mpres_gcd (f, T2, n);

      mpres_clear (T1 ,n);
      mpres_clear (T2 ,n);
      mpres_clear (T3 ,n);
      mpres_clear (T4 ,n);
      mpres_clear (T5 ,n);
      mpres_clear (T6 ,n);
      mpres_clear (T7 ,n);
      
      mpres_clear (amod ,n);
      mpres_clear (bmod ,n);
      mpres_clear (u ,n);
      mpres_clear (v ,n);

      return GENERATION_NOT_CURVE;
    }

  mpres_mul (th->ga ,th->ga ,cHEll->q ,n); //  t7/(a*b*(x^2-z^2))
  mpres_mul (u ,u ,cHEll->q ,n); // u = 1/( (a*b)^2 * (x^2-z^2) );
  mpres_mul (cHEll->mu ,T7 ,u ,n); // t7/( (a*b)^2 * (x^2-z^2) );

  mpres_mul (v ,amod ,amod ,n); // v=a^2
  mpres_mul (cHEll->la ,cHEll->mu ,v ,n); // t7/(b^2*(x^2-z^2));
  mpres_mul (cHEll->nu ,v ,u ,n);
  mpres_mul (cHEll->nu ,cHEll->nu ,T5 ,n);
  mpres_neg (cHEll->nu ,cHEll->nu ,n);
  // nu = -(a^2*z^2-b^2*x^2)/(b^2*(x^2-z^2))
  // We finished nu

  mpres_mul (cHEll->mu ,cHEll->mu ,cHEll->q ,n);
  // mu = t7 / ( (a*b)^2 * (x^2-z^2)^2 );
  // We finished mu


  if (!mpres_invert (v ,T5 ,n))  // 1/(a^2*z^2-b^2*x^2)
    {
      mpres_gcd (f, T5, n);

      mpres_clear (T1 ,n);
      mpres_clear (T2 ,n);
      mpres_clear (T3 ,n);
      mpres_clear (T4 ,n);
      mpres_clear (T5 ,n);
      mpres_clear (T6 ,n);
      mpres_clear (T7 ,n);
      
      mpres_clear (amod ,n);
      mpres_clear (bmod ,n);
      mpres_clear (u ,n);
      mpres_clear (v ,n);

      return GENERATION_NOT_CURVE;
    }
  mpres_neg (v ,v ,n);   // v = - 1/(a^2*z^2-b^2*x^2)

  mpres_mul (th->ga ,th->ga ,v ,n); // ga=-t7/(a*b*(x^2-z^2)*(a^2*z^2-b^2*x^2))
  mpres_mul (cHEll->la ,cHEll->la ,v ,n);
  // la = - t7/(b^2*(x^2-z^2)*(a^2*z^2-b^2*x^2));
  mpres_mul (th->t10p ,th->t10p ,v ,n);
  // t10p = -a*b*(x^2-z^2)/(a^2*z^2-b^2*x^2)
  // We finished ga, la, t10p

  mpres_set_z (u ,z ,n);
  mpres_mul (cHEll->q ,cHEll->q ,u ,n);
  mpres_mul (cHEll->q ,cHEll->q ,cHEll->q ,n);
  mpres_mul (cHEll->q ,cHEll->q ,u ,n); // z^3 / (x^2-z^2)^2
  mpres_set_z (u ,x ,n);
  mpres_mul (cHEll->q ,cHEll->q ,u ,n);
  mpres_mul (cHEll->q ,cHEll->q ,para->y ,n); // x*y*z^3 / (x^2-z^2)^2
  mpres_mul (u ,para->s ,para->s ,n);
  mpres_sub_ui (u ,u ,1 ,n); // s^2-1
  mpres_mul (cHEll->q ,cHEll->q ,u ,n);
  // q = x * y * (a^2-b^2) * z^3  /  ( b^2 * (x^2-z^2)^2 )
  // We finished q


  mpres_mul (u ,T3 ,T4 ,n);
  if (!mpres_invert (u ,u ,n)) //u = 1/( (a*z^2 + b*x^2)^2 * (a*z^2 - b*x^2)^2 )
    {
      mpres_gcd (f, u, n);

      mpres_clear (T1 ,n);
      mpres_clear (T2 ,n);
      mpres_clear (T3 ,n);
      mpres_clear (T4 ,n);
      mpres_clear (T5 ,n);
      mpres_clear (T6 ,n);
      mpres_clear (T7 ,n);
      
      mpres_clear (amod ,n);
      mpres_clear (bmod ,n);
      mpres_clear (u ,n);
      mpres_clear (v ,n);

      return GENERATION_NOT_CURVE;
    }
  mpres_mul (th->p ,T6 ,T6 ,n);
  mpres_mul (th->p ,th->p ,u ,n);
  mpres_neg (th->p ,th->p ,n);
  // p = - (a^2*z^4 - 2*b^2*x^2*z^2 + b^2*x^4) / ( (a*z^2 + b*x^2)^2 * (a*z^2 - b*x^2)^2 )
  // We finished p




  mpres_sub (th->t7p ,th->ga ,th->be ,n);
  // We finished t7p

  mpres_mul (th->Rac ,th->be ,cHEll->nu ,n);
  mpres_mul (th->Rac ,th->Rac ,th->t7p ,n);
  mpres_neg (th->Rac ,th->Rac ,n);
  // We finished Rac

  mpres_sub (th->t5p ,th->ga ,th->t10p ,n);
  mpres_sub (th->t6p ,th->be ,th->t10p ,n);
  // We finished t5p, t6p




  mpres_clear (T1 ,n);
  mpres_clear (T2 ,n);
  mpres_clear (T3 ,n);
  mpres_clear (T4 ,n);
  mpres_clear (T5 ,n);
  mpres_clear (T6 ,n);
  mpres_clear (T7 ,n);
  
  mpres_clear (amod ,n);
  mpres_clear (bmod ,n);
  mpres_clear (u ,n);
  mpres_clear (v ,n);



  return GENERATION_A_CURVE;
}












/* 
   generate curve with small parameters
   We do the special case nJacobi=2 "by hand". It produces enought curves for 
   finding factors of at least 60 digits
   In this case we only need to generate s=a/b since x,y,z are fixed.
   We generate the curve one by one such that a+b=c (for constants c) then we go
   to the next c
   We need a>0, b>0. Moreover s=2,4 doesn't work
*/
int generateCurveSmallParam (mpz_t f,mpmod_t n,
			     paraGenCurve para ,thetaCst th,
			     curveHyperEll cHEll,ksPoint P,
			     ksSmallConstPourMul cMul,
			     mpz_t a, mpz_t b,
			     int nJacobi,
			     optionsHECM options) {

  int test;
  mpz_t x,y,z;
  mpres_t u;
  mpres_t g;

  mpz_init (x);
  mpz_init (y);
  mpz_init (z);
  mpres_init (u ,n);
  mpres_init (g ,n);




  do {

    mulJacobiEntiers (a,b,nJacobi,x,y,z);
    

    test = generateOneCurveSmallParam (f,n,para,th,cHEll,P,cMul,a,b,x,y,z);

    if (test == GENERATION_NOT_CURVE) {
      mpres_set (g ,f ,n);
      if  (  mpres_is_zero(g,n) == 1  ) {// f =0 mod n
	// Let's try again
	test = nextParam (f,n,para,options);
	if ( test == NEXT_SMALL_PARAM_TOO_BIG ) { // TODO go to nJacobi=3
	  mpz_clear (x);
	  mpz_clear (y);
	  mpz_clear (z);
	  mpres_clear (u,n);
	  mpres_clear (g,n);
	  return GENERATION_PARAM_TOO_BIG; // TODO idem
	} 
      }
      else { 
	mpz_clear (x);
	mpz_clear (y);
	mpz_clear (z);
	mpres_clear (u ,n);
	mpres_clear (g ,n);
	return GENERATION_FACTOR_FOUND; 
      }
    }	
    else if (test == GENERATION_THESE_SMALL_PARAM_TOO_BIG) {
      // Let's try again
      test = nextParam (f,n,para,options);
      if ( test == NEXT_SMALL_PARAM_TOO_BIG ) { // TODO go to nJacobi=3
	mpz_clear (x);
	mpz_clear (y);
	mpz_clear (z);
	mpres_clear (u,n);
	mpres_clear (g,n);
	return GENERATION_PARAM_TOO_BIG; // TODO cf au idem
      } 
    }
    else { // We have genereted a curve

      mpres_sub_ui (u ,cHEll->mu ,1 ,n);
      mpres_mul (u ,u ,cHEll->mu ,n);
      mpres_mul (u ,u ,th->t7p ,n);
      mpres_mul (u ,u ,th->t6p ,n);

      mpres_gcd (f ,u ,n);
      mpres_set_z (g ,f ,n);

      if ( !mpz_cmp_ui (f,1) ) {
	// general case, f=1, ie non zero modulo all factors of n
	mpz_clear (x);
	mpz_clear (y);
	mpz_clear (z);
	mpres_clear (u,n);
	mpres_clear (g,n);
	return GENERATION_CORRECT_CURVE;
      }
      else if (  mpres_is_zero(g,n) != 1  ) {
	// f is a real divisor of n
	mpz_clear (x);
	mpz_clear (y);
	mpz_clear (z);
	mpres_clear (u,n);
	mpres_clear (g,n);
	return GENERATION_FACTOR_FOUND; 
      }
      else {
      // Let's try again
	test = nextParam (f,n,para,options);
	if ( test == NEXT_SMALL_PARAM_TOO_BIG ) { // TODO go to nJacobi=3
	  mpz_clear (x);
	  mpz_clear (y);
	  mpz_clear (z);
	  mpres_clear (u,n);
	  mpres_clear (g,n);
	  return GENERATION_PARAM_TOO_BIG; // TODO idem
	} 
      }
    }
  } while ( 1 );

}







int generateCurveSmallParamSpecified (mpz_t f,mpmod_t n,
				      paraGenCurve para,thetaCst th,
				      curveHyperEll cHEll,ksPoint P,
				      ksSmallConstPourMul cMul) {
  int test;
  mpz_t x,y,z;
  mpres_t u;
  mpres_t g;

  mpz_init (x);
  mpz_init (y);
  mpz_init (z);

  mpres_init (u,n);
  mpres_init (g,n);



  mulJacobiEntiers (para->a,para->b,para->nJacobi,x,y,z);




  mpres_set_z (para->s ,para->b ,n);
  if (!mpres_invert (para->s ,para->s ,n)) // s=1/b
    {
      mpres_gcd (f, para->s, n);
      mpres_set_z (g ,f ,n);
      if ( mpres_is_zero (g,n)==1 ) { // f=0 mod n
	mpres_set_ui(para->s ,0 ,n);

	mpz_clear (x);
	mpz_clear (y);
	mpz_clear (z);
	mpres_clear (u,n);
	mpres_clear (g,n);
	return GENERATION_FAIL;
      }
      else {
	mpres_set_ui(para->s ,0 ,n);

	mpz_clear (x);
	mpz_clear (y);
	mpz_clear (z);
	mpres_clear (u,n);
	mpres_clear (g,n);
	return GENERATION_FACTOR_FOUND; 
      }
    }
  mpres_set_z (u ,para->a ,n);
  mpres_mul (para->s ,para->s ,u ,n);

  test = generateOneCurveSmallParam (f ,n ,para ,th ,cHEll ,P ,cMul ,para->a ,para->b ,x ,y ,z);

  if (test == GENERATION_NOT_CURVE) {
    mpres_set (g ,f ,n);
    if  (  mpres_is_zero(g,n) == 1  ) {// f =0 mod n
      mpz_clear (x);
      mpz_clear (y);
      mpz_clear (z);
      mpres_clear (u,n);
      mpres_clear (g,n);
      return GENERATION_FAIL;
    }
    else { 
      mpz_clear (x);
      mpz_clear (y);
      mpz_clear (z);
      mpres_clear (u,n);
      mpres_clear (g,n);
      return GENERATION_FACTOR_FOUND; 
    }
  }	
  else if (test == GENERATION_THESE_SMALL_PARAM_TOO_BIG) {
    mpz_clear (x);
    mpz_clear (y);
    mpz_clear (z);
    mpres_clear (u,n);
    mpres_clear (g,n);
    return GENERATION_FAIL;
  }
  else { // We have genereted a curve

    mpres_sub_ui (u ,cHEll->mu ,1 ,n);
    mpres_mul (u ,u ,cHEll->mu ,n);
    mpres_mul (u ,u ,th->t7p ,n);
    mpres_mul (u ,u ,th->t6p ,n);

    mpres_gcd (f ,u ,n);
    mpres_set_z (g ,f ,n);

    if ( !mpz_cmp_ui (f,1) ) {
      // general case, f=1, ie non zero modulo all factors of n
      mpz_clear (x);
      mpz_clear (y);
      mpz_clear (z);
      mpres_clear (u,n);
      mpres_clear (g,n);
      return GENERATION_CORRECT_CURVE;
    }
    else if (  mpres_is_zero(g,n) != 1  ) {
      // f is a real divisor of n
      mpz_clear (x);
      mpz_clear (y);
      mpz_clear (z);
      mpres_clear (u,n);
      mpres_clear (g,n);
      return GENERATION_FACTOR_FOUND; 
    }
    else {
      mpz_clear (x);
      mpz_clear (y);
      mpz_clear (z);
      mpres_clear (u,n);
      mpres_clear (g,n);
      return GENERATION_FAIL;
    } 
  }
}






int nextSmallParam (mpz_t a, mpz_t b,
		  const long Ha,const long Hb) {

  mpz_t g;
  mpz_init (g);
  
  do {
    if ( mpz_cmp (a,b) < 0 ){
      if ( mpz_cmp_ui (a,Ha) < 0 ) {
	mpz_add_ui (a ,a ,1); // a++
      }
      else if ( mpz_cmp_ui (b,Hb) < 0 ) {
	mpz_add_ui (b ,b ,1); // b++
	mpz_set_ui (a,1);     // a=1
      }
      else {
	mpz_clear (g);
	return NEXT_SMALL_PARAM_TOO_BIG;
      }
    }
    else if ( mpz_cmp_ui (b,1) > 0) {
      mpz_sub_ui (b ,b ,1); // b--
    }
    else {
      if ( mpz_cmp_ui (a,Hb) < 0) {
	mpz_add_ui (b ,a ,1); // b=a+1
	mpz_set_ui (a ,1);    // a=1
	}
      else if ( mpz_cmp_ui (a,Ha) < 0 ) {
	mpz_add_ui (a ,a ,1); // a++
	mpz_set_ui (b ,Hb);   // b=Hb
      }
      else {
	mpz_clear (g);
	return NEXT_SMALL_PARAM_TOO_BIG;	
      }
    }
    mpz_gcd (g ,a ,b);
  } while ( mpz_cmp_ui (g,1) != 0 ); // while g != 1
  
  mpz_clear (g);
  return NEXT_PARAM_CAN_BE_USED;

}



int nextParam (mpz_t f,mpmod_t n,
	       paraGenCurve para,
	       optionsHECM options) {
  int test;

  if ( options->smallParam == TRUE ) {

    long Ha,Hb;
    mpz_t t1;
    mpz_init (t1);

    if (para->nJacobi==2) {
      mpz_set_ui (t1 ,LONG_MAX);
      mpz_mul_ui (t1 ,t1 ,4);
      mpz_root (t1 ,t1 ,5);
      Ha = mpz_get_ui (t1); // Ha = 1600
      Hb = Ha;              // Hb = 1000
    }
    else if (para->nJacobi==3) {
      mpz_set_ui (t1 ,LONG_MAX);
      mpz_mul_ui (t1 ,t1 ,3194512); // *2^15*3^8/673
      mpz_root (t1 ,t1 ,8);
      Hb = mpz_get_ui (t1);         // Hb = 200
      Ha = Hb*3; // Hb*673,2^(1/8)  // Ha = 100
    }
    else if (para->nJacobi==4) {
      mpz_set_ui (t1 ,LONG_MAX);
      mpz_mul_ui (t1 ,t1 ,65536); 
      mpz_mul_ui (t1 ,t1 ,112101); // *2^16*3^15/128
      mpz_root (t1 ,t1 ,12);
      Hb = mpz_get_ui (t1);         // Hb = 30
      Ha = Hb*2; // Hb*128^(1/8)    // Ha = 20
    }
    else {
      // TODO by default we take the same bound than for nJacobi = 4
      //      What is the real bound?
      mpz_set_ui (t1 ,LONG_MAX);
      mpz_mul_ui (t1 ,t1 ,65536); 
      mpz_mul_ui (t1 ,t1 ,112101);
      mpz_root (t1 ,t1 ,12);
      Hb = mpz_get_ui (t1);
      Ha = Hb*2;
    }


    test = nextSmallParam (para->a,para->b,Ha,Hb);

    if (test==NEXT_PARAM_CAN_BE_USED) {
      mpz_clear (t1);
      return NEXT_PARAM_CAN_BE_USED;
    }
    else { // test==NEXT_SMALL_PARAM_TOO_BIG 
      if ( para->nJacobi < N_JACOBI_SMALL_PARAM_MAX ) {
	para->nJacobi++;
	mpz_set_ui (para->a ,1);
	mpz_set_ui (para->b ,2);
	mpz_clear (t1);
	return NEXT_PARAM_CAN_BE_USED;
      }
      else {
	mpz_clear (t1);
	return NEXT_SMALL_PARAM_TOO_BIG;
      }
    }
  }
  else { // (options->smallParam == FALSE )
    para->nJacobi = (rand() % (N_JACOBI_MAX_p1 - N_JACOBI_MIN)) + N_JACOBI_MIN;
    mpres_add_ui (para->s ,para->s ,1 ,n);
    return NEXT_PARAM_CAN_BE_USED;
  }

}
