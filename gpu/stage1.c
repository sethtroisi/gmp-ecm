# include <stdio.h>
# include <math.h>
# include <gmp.h>
# include "prototype.h"

main(int argc,char*argv[])
{
	if (argc != 4)
		printf("Error in call function\n./stage1 N sigma B1\n");
	else {
		mpz_t N; mpz_init(N); mpz_set_str(N,argv[1],10); // in base 10
		mpz_t sigma; mpz_init(sigma); mpz_set_str(sigma,argv[2],10);
		mpz_t B1; mpz_init(B1); mpz_set_str(B1,argv[3],10);
		mpz_t x; mpz_init(x);
		
		stageOne(B1,sigma,N,x);
		gmp_printf("\nStep 1\nx=%Zd\n\n",x);
		
		mpz_clear(sigma);
		mpz_clear(N);
		mpz_clear(B1);
		mpz_clear(x);
	}
}
