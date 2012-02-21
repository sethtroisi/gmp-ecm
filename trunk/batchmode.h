#define BATCHMODE 1 /* use an unsigned long for 'd' */
//#define BATCHMODE 2 /* use an mpz_t for 'd' */

/* ellparam_batch.c */
int get_curve_from_ell_parametrization (mpz_t, mpres_t, mpz_t, mpmod_t);
