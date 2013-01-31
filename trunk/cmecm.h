int build_curves_with_CM(mpz_t f, int *nE, ec_curve_t *tE, ec_point_t *tP, int disc, mpmod_t n, mpz_t *sqroots);
void
set_stage2_params_CM(unsigned long *pdF, unsigned long *pk, mpz_t B2, int disc);
int ecm_rootsF_CM(mpz_t f, listz_t F, unsigned long dF, curve *C, mpmod_t modulus);
