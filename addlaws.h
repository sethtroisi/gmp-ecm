#define DEBUG_ADD_LAWS 0
#define USE_ADD_SUB_CHAINS 0

#define pt_is_equal(P, Q) (mpz_cmp((P)->x, (Q)->x) == 0 \
	                     && mpz_cmp((P)->y, (Q)->y) == 0 \
			     && mpz_cmp((P)->z, (Q)->z) == 0)

void pt_print(ell_curve_t E, ell_point_t P, mpmod_t n);
int pt_is_zero(ell_point_t P, ATTRIBUTE_UNUSED mpmod_t n);
void pt_set_to_zero(ell_point_t P, mpmod_t n);
void pt_assign(ell_point_t Q, ell_point_t P, ATTRIBUTE_UNUSED mpmod_t n);
void pt_neg(ell_point_t P, mpmod_t n);

int hessian_to_weierstrass(mpz_t f, mpres_t x, mpres_t y, mpres_t D, mpmod_t n);
int
twisted_hessian_to_weierstrass(mpz_t f, mpres_t x, mpres_t y, mpres_t c, mpres_t d, mpmod_t n);

size_t build_MO_chain(short *S, size_t Slen, mpz_t e, int w);
size_t build_add_sub_chain(short *S, size_t Slen, mpz_t e, int w);
int compute_s_4_add_sub(mpz_t s, ecm_uint B1, int disc);

int mult_by_3(mpz_t f, mpres_t x, mpres_t y, mpres_t A, mpmod_t n);
void ell_point_init(ell_point_t P, ell_curve_t E, mpmod_t n);
void ell_point_clear(ell_point_t P, ATTRIBUTE_UNUSED ell_curve_t E, mpmod_t n);
void ell_point_print(ell_point_t P, ell_curve_t E, mpmod_t n);
void ell_point_set(ell_point_t Q, ell_point_t P,
		  ATTRIBUTE_UNUSED ell_curve_t E, ATTRIBUTE_UNUSED mpmod_t n);
void ell_curve_init(ell_curve_t E, int etype, int law, mpmod_t n);
void ell_curve_init_set(ell_curve_t E, int type, int law, mpres_t A, mpmod_t n);
void ell_curve_set_z(ell_curve_t E, ell_curve_t zE, mpmod_t n);
void ell_curve_clear(ell_curve_t E, mpmod_t n);
void ell_curve_print(ell_curve_t E, mpmod_t n);
int ell_point_is_on_curve(ell_point_t P, ell_curve_t E, mpmod_t n);
int ell_point_is_zero(ell_point_t P, ell_curve_t E, mpmod_t n);
void ell_point_set_to_zero(ell_point_t P, ell_curve_t E, mpmod_t n);
int ell_point_add(mpz_t f, ell_point_t R, ell_point_t P, ell_point_t Q, ell_curve_t E, mpmod_t n);
int ell_point_sub(mpz_t f, ell_point_t R, ell_point_t P, ell_point_t Q, ell_curve_t E, mpmod_t n);
int ell_point_duplicate(mpz_t f, ell_point_t R, ell_point_t P, ell_curve_t E, mpmod_t n);
void ell_point_negate(ell_point_t P, ell_curve_t E, mpmod_t n);
int ell_point_mul_plain (mpz_t f, ell_point_t Q, mpz_t e, ell_point_t P, ell_curve_t E, mpmod_t n);
int get_add_sub_w(mpz_t e);
void add_sub_pack(mpz_t s, int w, short *S, size_t iS);
void add_sub_unpack(int *w, short **S, size_t *iS, mpz_t s);
int ell_point_mul_add_sub_with_S(mpz_t f, ell_point_t Q, ell_point_t P, ell_curve_t E,mpmod_t n, int w, short *S, int iS);
int ell_point_mul_add_sub (mpz_t f, ell_point_t Q, mpz_t e, ell_point_t P,
			  ell_curve_t E, mpmod_t n);
int ell_point_mul(mpz_t f, ell_point_t Q, mpz_t e, ell_point_t P, ell_curve_t E, mpmod_t n);
