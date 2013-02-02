#define CMECM_D_4 1 /* fast version */

int build_curves_with_CM(mpz_t f, int *nE, ec_curve_t *tE, ec_point_t *tP, int disc, mpmod_t n, mpz_t *sqroots);
void
set_stage2_params_CM(unsigned long *pdF, unsigned long *pk, mpz_t B2, int disc);
int ecm_rootsF_CM(mpz_t f, listz_t F, unsigned long dF, curve *C, mpmod_t modulus);
ecm_roots_state_t *ecm_rootsG_init_CM (mpz_t f, curve *X, root_params_t *root_params, unsigned long dF, unsigned long blocks, mpmod_t modulus);
int ecm_rootsG_CM (mpz_t f, listz_t G, unsigned long dF, ecm_roots_state_t *state, mpmod_t modulus);
