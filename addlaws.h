#define pt_is_equal(P, Q) (mpz_cmp((P)->x, (Q)->x) == 0 \
	                     && mpz_cmp((P)->y, (Q)->y) == 0 \
			     && mpz_cmp((P)->z, (Q)->z) == 0)

int
pt_is_zero(ec_point_t P, ATTRIBUTE_UNUSED mpmod_t n);
void
pt_set_to_zero(ec_point_t P, mpmod_t n);
void
pt_assign(ec_point_t Q, ec_point_t P, ATTRIBUTE_UNUSED mpmod_t n);
void
pt_neg(ec_point_t P, mpmod_t n);
void
pt_many_set_to_zero(ec_point_t *tP, int nE, mpmod_t n);
void
pt_many_neg(ec_point_t *tP, int nE, mpmod_t n);
void
pt_many_assign(ec_point_t *tQ, ec_point_t *tP, int nE, mpmod_t n);
void
print_mpz_from_mpres(mpres_t x, mpmod_t n);
int
pt_many_duplicate(ec_point_t *tQ, ec_point_t *tP, ec_curve_t *tE, int nE, 
		  mpmod_t n, 
		  mpres_t *num, mpres_t *den, mpres_t *inv, char *ok);
int
pt_many_mul(ec_point_t *tQ, ec_point_t *tP, ec_curve_t *tE, int nE,
	    mpz_t e, mpmod_t n, 
	    mpres_t *num, mpres_t *den, mpres_t *inv, char *ok);
