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
void print_mpz_from_mpres(mpres_t x, mpmod_t n);
int
pt_many_duplicate(ec_point_t *tQ, ec_point_t *tP, ec_curve_t *tE, int nE, 
		  mpmod_t n, 
		  mpres_t *num, mpres_t *den, mpres_t *inv, char *ok);
int
pt_many_mul(ec_point_t *tQ, ec_point_t *tP, ec_curve_t *tE, int nE,
	    mpz_t e, mpmod_t n, 
	    mpres_t *num, mpres_t *den, mpres_t *inv, char *ok);

int hessian_to_weierstrass(mpz_t f, mpres_t x, mpres_t y, mpres_t D, mpmod_t n);

int build_MO_chain(short *S, int Slen, mpz_t e, int w);
int build_add_sub_chain(short *S, int Slen, mpz_t e, int w);
int compute_s_4_add_sub(mpz_t s, unsigned long B1, int disc);

int mult_by_3(mpz_t f, mpres_t x, mpres_t y, mpres_t A, mpmod_t n);
void ec_point_init(ec_point_t P, ec_curve_t E, mpmod_t n);
void ec_point_clear(ec_point_t P, ATTRIBUTE_UNUSED ec_curve_t E, mpmod_t n);
void ec_point_print(ec_point_t P, ec_curve_t E, mpmod_t n);
void ec_point_set(ec_point_t Q, ec_point_t P,
		  ATTRIBUTE_UNUSED ec_curve_t E, ATTRIBUTE_UNUSED mpmod_t n);
void ec_curve_init(ec_curve_t E, int etype, mpmod_t n);
void ec_curve_init_set(ec_curve_t E, int type, mpres_t A, mpmod_t n);
void ec_curve_set_z(ec_curve_t E, ec_curve_t zE, mpmod_t n);
void ec_curve_clear(ec_curve_t E, mpmod_t n);
void ec_curve_print(ec_curve_t E, mpmod_t n);
int ec_point_is_zero(ec_point_t P, ec_curve_t E, mpmod_t n);
void ec_point_set_to_zero(ec_point_t P, ec_curve_t E, mpmod_t n);
int ec_point_add(ec_point_t R, ec_point_t P, ec_point_t Q, ec_curve_t E, mpmod_t n);
int ec_point_sub(ec_point_t R, ec_point_t P, ec_point_t Q, ec_curve_t E, mpmod_t n);
int ec_point_duplicate(ec_point_t R, ec_point_t P, ec_curve_t E, mpmod_t n);
void ec_point_negate(ec_point_t P, ec_curve_t E, mpmod_t n);
int ec_point_mul_plain (ec_point_t Q, mpz_t e, ec_point_t P, ec_curve_t E, mpmod_t n);
int get_add_sub_w(mpz_t e);
void add_sub_pack(mpz_t s, int w, short *S, int iS);
void add_sub_unpack(int *w, short **S, int *iS, mpz_t s);
int ec_point_mul_add_sub_with_S(ec_point_t Q, ec_point_t P, ec_curve_t E,
				mpmod_t n, int w, short *S, int iS);
int ec_point_mul_add_sub (ec_point_t Q, mpz_t e, ec_point_t P,
			  ec_curve_t E, mpmod_t n);
int ec_point_mul(ec_point_t Q, mpz_t e, ec_point_t P, ec_curve_t E, mpmod_t n);
