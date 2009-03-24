#include <string.h>
#include <time.h>
#include <limits.h>

#include "../ecm-impl.h"

#include "hecm.h"
#include "auxi.h"
#include "generation.h"
#include "morphismes.h"






void optionsHECM_init (optionsHECM options) {
  mpz_init (options->heightMin);
  mpz_init (options->heightMax);

  options->smallParam = TRUE;
  options->curveSpecified = FALSE;
  options->initialCurveSpecified = FALSE;
  options->verbose = OUTPUT_NORMAL;

  options->nbtests = 1;
  mpz_set_ui (options->heightMin,0);
  mpz_set_ui (options->heightMax,0);
}


void optionsHECM_clear (optionsHECM options) {
  mpz_clear (options->heightMin);
  mpz_clear (options->heightMax);
}


void printCurve (paraGenCurve para,optionsHECM options,mpmod_t n) {
  mpz_t s;
  mpz_init (s);

  mpres_get_z (s ,para->s ,n);

  if (options->smallParam == TRUE) {
    gmp_printf ("s=%Zd/%Zd",para->a,para->b);
  }
  else {
    gmp_printf("s=%Zd",s );
  }

  printf(" and nJacobi=%d\n",para->nJacobi);

  mpz_clear (s);
}



static void 
usage (void)
{
  printf ("Perform ECM with hyperelliptic curves of genus 2. The hyperelliptic curves are\n");
  printf ("(2,2)-decomposable which means that their jacobians is isogenous to the product\n");
  printf ("of two elliptic curves. Therefore one run of HECM with this kind of curves is\n");
  printf ("equivalent to two runs of ECM on the underlying curves.\n");
  printf ("Use Kummer surfaces to speed up the arithmetics in the Jacobians.\n");
  printf ("The use of small parameters make it quicker than GMP-ECM for number of\n");
  printf ("at least 300 digits.\n");


  printf ("\nParameters:\n");
  printf ("-n           number to be factored in base 10\n");
  printf ("-B1          B1 in base 10\n");
  printf ("-B2          B2 in base 10\n");
  printf ("-c           Numbers of tests to do. If unspecified run only once.\n");
  printf ("             If 0 run untill it find a factor.\n");
  printf ("-v           verbose mode\n");

  printf ("-not_low_param Work with normal parameters (not small)\n");


  printf ("-spe         We try to generate the curve we this fixed parameters.\n"); 
  printf ("             One of the following option should have been specified:\n");
  printf ("              -ab if we want small parameters\n");
  printf ("              -s  if we want normal parametrization\n");
  printf ("             If the generated curve is correct we run HECM with it else we print\n");
  printf ("             an error\n");

  printf ("-ab a b      fix the parameter s=a/b. We try to use small parameters if possible\n");
  printf ("-h_min       minimal height of s. If equal to 0 or unspecified, then we begin\n");
  printf ("             by 1\n");
  printf ("-h_max       maximal height of s. If equal to 0 or unspecified, then no limit\n");

  printf ("-nJacobi     For the parametrization, we need to have a point on a Jacobi\n");
  printf ("             curve. We begin by P=(1,*) and we use the point nJacobi*P=(x,y)\n");
  printf ("             If nJacobi<0, we use |nJacobi| (absolute value)\n");
  printf ("             If nJacobi=1, we use nJacobi=2 (first usable point)\n");
  printf ("             If nJacobi=0 or unspecified, and in the case of small parameters,\n");
  printf ("             we use nJacobi=2 and increase first a,b. If we can't increase a,b\n");
  printf ("             we increase nJacobi by 1...\n");
  printf ("             If nJacobi=0 or unspecified, and in the case of normal\n");
  printf ("             parametrization, we take random nJacobi for each run.\n");

  printf ("-s s         fix the parameter s for normal parametrization\n");
  printf ("-seed        fix the seed. Only usefull with normal parametrization\n");

 

  printf ("-h, --help   Prints this help and exit.\n");
}




double B1 = 25.0;
char * number = "4816415081";   //      58027 83003 
char * Nsigma = "2";   //      s=2                   
int seed = 0;
char * charB2 = "-1";

int main(int argc, char * argv[]) {


  clock_t
    temps_initial, // initial time in micro-seconds
    temps_final;   // final time in micro-seconds
  double 
    temps_cpu;     // total time in seconds

  temps_cpu=0;
  temps_initial = clock ();


  optionsHECM options;
  optionsHECM_init (options);

  mpmod_t n;
  mpz_t k;
  mpz_t f; // factor of N
  mpz_t s; // The parameter s in mpz_t
  paraGenCurve para;
  curve T1,T2;
  int test,test2;
  unsigned int nb;
  double B2scale = 1.0;
  mpz_t B2;

  char *a, *b, *heightMax, *heightMin;

  para->nJacobi = 2;
  a="1";b="2";
  mpz_init_set_si (B2,ECM_DEFAULT_B2); // set B2


  while ((argc > 1) && (argv[1][0] == '-')) { 
    if  ((argc > 2) && (strcmp (argv[1], "-n")) == 0) {
      number = argv[2];
      argv += 2;
      argc -= 2;
    }
    else if ((argc > 2) && (strcmp (argv[1], "-s")) == 0) {
      Nsigma = argv[2];
      options->initialCurveSpecified = TRUE;
      options->smallParam = FALSE;
      argv += 2;
      argc -= 2;
    }
    else if (strcmp (argv[1], "-h") == 0 || strcmp (argv[1], "--help") == 0) {
      usage ();
      mpz_clear(B2);
      exit (0);
    }
    else if ((argc > 3) && (strcmp (argv[1], "-ab")) == 0) {
      a = argv[2];
      b = argv[3];
      options->initialCurveSpecified = TRUE;
      options->smallParam = TRUE;
      argv += 3;
      argc -= 3;
      }
    else if ((argc > 2) && (strcmp (argv[1], "-nJacobi")) == 0) {
      para->nJacobi = atoi(argv[2]);
      options->initialCurveSpecified = TRUE;
      argv += 2;
      argc -= 2;
    }
    else if ((argc > 2) && (strcmp (argv[1], "-B1")) == 0) {
      B1 = atof(argv[2]);
      argv += 2;
      argc -= 2;
    }
    else if  ((argc > 2) && (strcmp (argv[1], "-B2")) == 0) {
      mpz_set_str (B2,argv[2],10);
      argv += 2;
      argc -= 2;
    }
    else if ((argc > 2) && (strcmp (argv[1], "-seed")) == 0) {
      seed = atoi(argv[2]);
      argv += 2;
      argc -= 2;
    }
    else if (strcmp (argv[1], "-not_low_param") == 0) {
      options->smallParam = FALSE;
      argv ++;
      argc --;
    }
    else if (strcmp (argv[1], "-spe") == 0) {
      options->initialCurveSpecified = TRUE;
      options->curveSpecified = TRUE;
      options->nbtests = 1;
      argv ++;
      argc --;
    }
    else if (strcmp (argv[1], "-v") == 0) {
      options->verbose ++;
      argv ++;
      argc --;
    }
    else if ((argc > 2) && (strcmp (argv[1], "-c")) == 0) {
      options->nbtests = atoi(argv[2]);
      argv += 2;
      argc -= 2;
    }
    else if ((argc > 2) && (strcmp (argv[1], "-h_max")) == 0) {
      heightMax = argv[2];
      mpz_set_str (options->heightMax ,heightMax ,10);
      argv += 2;
      argc -= 2;
    }
    else if ((argc > 2) && (strcmp (argv[1], "-h_min")) == 0) {
      heightMin = argv[2];
      mpz_set_str (options->heightMin ,heightMin ,10);
      argv += 2;
      argc -= 2;
    }
    else {
      fprintf(stderr, "Unexpected option: %s\n", argv[1]);
      mpz_clear(B2);
      exit(1);
    }
  }



  if (seed == 0) {
      srand(time(NULL));
  } else {
      srand(seed);
  }


  mpz_t N;
  mpz_init_set_str (N,number,10); // number in base 10
  mpmod_init(n,N,ECM_MOD_DEFAULT);
  mpz_clear (N); 

  


  if ( mpz_sgn (options->heightMax) ==0 ) { // heightMax = 0
    mpz_set (options->heightMax,n->orig_modulus); // heightMax = n
  }
  if ( mpz_cmp (options->heightMin ,options->heightMax) > 0 ) {
    gmp_printf ("We must have h_min <= h_max. Here we have %Zd > %Zd.\n",options->heightMin,options->heightMax);

    optionsHECM_clear (options);
    mpmod_clear(n);
    mpz_clear (B2);
    return 0;
  }





  paraGenCurve_init (para,n);
  mpz_init(s);


  if (options->initialCurveSpecified == TRUE) {

    if  (options->smallParam == TRUE) {
      mpz_set_str (para->a ,a ,10);
      mpz_set_str (para->b ,b ,10);

      if ( mpz_cmp_ui (para->b,0) == 0 ) { // b=0
	printf("We have s=a/b so b must be non zero.\n");

	optionsHECM_clear (options);
	mpz_clear(s);
	paraGenCurve_clear (para,n);
	mpmod_clear(n);
	mpz_clear (B2);
	return 0;
      }

      mpz_t d;
      mpz_init (d);

      mpz_gcd (d ,para->a ,para->b);
      mpz_divexact (para->a ,para->a ,d);
      mpz_divexact (para->b ,para->b ,d);

      mpz_clear (d);

      if ( ( mpz_cmp (para->a,options->heightMax) > 0 ) || ( mpz_cmp (para->b,options->heightMax) > 0 ) ) {
	printf("The height of the initial parameters is too big.\n");

	optionsHECM_clear (options);
	mpz_clear(s);
	paraGenCurve_clear (para,n);
	mpmod_clear(n);
	mpz_clear (B2);
	return 0;
      }


    }
    else {
      mpz_set_str (s,Nsigma,10);
      mpres_set_z (para->s,s,n);

      if ( ( mpz_cmp (s,options->heightMax) > 0 ) ) {
	printf("The height of the initial parameter is too big.\n");


	optionsHECM_clear (options);
	mpz_clear(s);
	paraGenCurve_clear (para,n);
	mpmod_clear(n);
	mpz_clear (B2);
	return 0;
      }

    }
  }
  else {
    if  (options->smallParam == TRUE) {
      if (mpz_cmp_ui (options->heightMax,2) < 0) {
	printf("There is no hyperelliptic curve with parameters smaller than this height.\n");

	optionsHECM_clear (options);
	mpz_clear(s);
	paraGenCurve_clear (para,n);
	mpmod_clear(n);
	mpz_clear (B2);
	return 0;
      }

      if ( mpz_cmp_ui (options->heightMin,2) < 0 ) {
	mpz_set_ui (para->a ,1);
	mpz_set_ui (para->b ,2);
      }
      else {
	mpz_set_ui (para->a ,1);
	mpz_set (para->b ,options->heightMin);
      }
      para->nJacobi = 2;

    }
    else {
      if (mpz_cmp_ui (options->heightMax,3) < 0) {
	printf("There is no hyperelliptic curve with parameters smaller than this height.\n");

	optionsHECM_clear (options);
	mpz_clear(s);
	paraGenCurve_clear (para,n);
	mpmod_clear(n);
	mpz_clear (B2);
	return 0;
      }

      para->nJacobi= (rand() % (N_JACOBI_MAX_p1 - N_JACOBI_MIN)) + N_JACOBI_MIN;
      // TODO
      if ( mpz_cmp_ui (options->heightMin,3) < 0 ) {
	mpres_set_ui (para->s ,3 ,n);
      }
      else {
	mpres_set_z (para->s ,options->heightMin ,n);
      }

    }
  }







  if (options->nbtests <= 0) {
    options->nbtests = UINT_MAX;
  }



  mpres_init (T1.x,n);
  mpres_init (T1.y,n);
  mpres_init (T1.A,n);
  mpres_init (T2.x,n);
  mpres_init (T2.y,n);
  mpres_init (T2.A,n);

  mpz_init(f);



  mpz_init(k);
  prodTreeCalculk(k,B1); // k=lcm(2,3,..,B1)
  getprime_clear ();


  if ( ECM_IS_DEFAULT_B2(B2) ) {    
    mpz_set_d (B2, B2scale * pow (ECM_COST * B1, DEFAULT_B2_EXPONENT));
    gmp_printf ("Using B1=%lu and B2=%Zd\n",(long) B1,B2);
  }






 // loop on hecm
  for(nb=1 ; nb <= options->nbtests ; nb++) {
    printf("Run %u on %u\n",nb,options->nbtests);



    // stage 1
    if ( options->smallParam == TRUE) {
      test = hecm1LowParam (f ,n ,k ,para ,&T1,&T2 ,options);
    }
    else {
      test = hecm1Normal (f ,n ,k ,para ,&T1,&T2 ,options);
    }

    if ( options->verbose >= OUTPUT_VERBOSE ) {
      // we need the parameters of the curve.
      printf ("The curve is given by ");
      printCurve (para,options,n);
    }

    if (test != HECM_NO_FACTOR_FOUND) {
      break;
    }




    // stage 2
    test = hecm2 (f,n,&T1,&T2,B2);
    if (test != HECM_NO_FACTOR_FOUND) {
      break;
    }



    // find parameters for the next curve
    test2 = nextParam (f,n,para,options);
    if (test2 == NEXT_SMALL_PARAM_TOO_BIG) {
      test = HECM_PARAM_TOO_BIG;
      break;
    }
    else if (test2 == NEXT_PARAM_ERROR) {
      test = HECM_ERROR;
      break;
    }

  }




  if (nb>options->nbtests) {
    nb = options->nbtests;
  }

  if ( ( para->nJacobi == 0 ) || ( para->nJacobi == 1 ) ) {
    para->nJacobi = 2;
  }



  if (test == HECM_ERROR) {
    printf("An error occured with the curve:\t" );
    printCurve(para,options,n);
    printf("It could be a multiplication by 0 in the projective space.\n" );
  }
  else if (test == HECM_GENERATION_FAIL) {
    printf("The generation of the curve failled. The parameters were:\t");
    printCurve(para,options,n);
  }
  else if (test == HECM_NO_FACTOR_FOUND) {
    printf("We test %u hyperelliptic curves (i.e. %u elliptic curves)\n",options->nbtests,options->nbtests*2); 
    printf("without finding factor of n.\n" );
  }
  else if (test == HECM_PARAM_TOO_BIG) {
    printf("There is no more curves with low parameters height.\n");
    printf("We test %u hyperelliptic curves (i.e. %u elliptic curves)\n",nb,nb*2); 
    printf("without finding factor of n.\n" );
  }
  else if (test == HECM_FOUND_N ) {
    printf("We divided by n\n");
    printf("We are on the curve of parameters:\t");
    printCurve(para,options,n);
    printf("We test %u hyperelliptic curves (i.e. %u elliptic curves)\n",nb,nb*2); 
  }
  else if (test == HECM_FOUND_ZERO_CURVE_1 ) {
    //TODO in that case, in general, we haven't check if we had the zero of the second curve
    printf ("k*P is the zero of the first elliptic curve modulo all the factors of n\n");
    printf("We are on the hyperelliptic curve of parameters:\t");
    printCurve(para,options,n);
    printf("We test %u hyperelliptic curves (i.e. %u elliptic curves)\n",nb,nb*2); 
  }
  else if (test == HECM_FOUND_ZERO_CURVE_2 ) {
    printf ("k*P is the zero of the second elliptic curve modulo all the factors of n\n");
    printf("We are on the hyperelliptic curve of parameters:\n");
    printf("\t");
    printCurve(para,options,n);
    printf("We test %u hyperelliptic curves (i.e. %u elliptic curves)\n",nb,nb*2); 
  }
  else if (test == HECM_FOUND_ZERO_CURVE_1_AND_2 ) {
    printf ("k*P is the zero of both elliptic curve modulo all the factors of n\n");
    printf("We are on the hyperelliptic curve of parameters:\t");
    printCurve(para,options,n);
    printf("We test %u hyperelliptic curves (i.e. %u elliptic curves)\n",nb,nb*2); 
  }
  else { // We found a real factor of n
    gmp_printf("A factor of n is %Zd\n",f);
    if ( options->verbose >= OUTPUT_VERBOSE ) {
      if (test == HECM_FACTOR_FOUND_MORPHISM) {
	printf("It was found during the computation of the morphisms.\n");
      }
      else if (test == HECM_FACTOR_FOUND_GENERATION) {
	printf("It was found during the generation of the curve.\n");
      }
      else if (test == HECM_FACTOR_FOUND_STEP1) {
	printf("It was found in stage 1.\n");
      }
      else if (test == HECM_FACTOR_FOUND_STEP2) {
	printf("It was found in stage 2.\n");
      }
    }
    printf("It was found by using the parameters:\t");
    printCurve(para,options,n);
    printf("We used %u hyperelliptic curves (i.e. %u elliptic curves).\n",nb,nb*2); 

  }


  optionsHECM_clear (options);
  mpz_clear(s);
  mpz_clear(k);
  mpz_clear(f);
  paraGenCurve_clear (para,n);
  mpres_clear (T1.x,n);
  mpres_clear (T1.y,n);
  mpres_clear (T1.A,n);
  mpres_clear (T2.x,n);
  mpres_clear (T2.y,n);
  mpres_clear (T2.A,n);
  mpmod_clear(n);
  mpz_clear (B2);

  temps_final = clock ();
  temps_cpu = (temps_final - temps_initial); 
  printf("The program took %g seconds.\n",temps_cpu/CLOCKS_PER_SEC);


  return 1;
}
