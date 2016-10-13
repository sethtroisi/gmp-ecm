void mod_div_2(mpz_t x, mpz_t n);
int mod_from_rat(mpz_t r, mpq_t q, mpz_t N);
int mod_from_rat2(mpz_t r, mpz_t num, mpz_t den, mpz_t N);
void ec_force_point(ell_curve_t E, ell_point_t P, long *x0, mpz_t n);

int
build_curves_with_torsion_Z5(mpz_t f, mpmod_t n, 
			     ell_curve_t *tE, ell_point_t *tP,
			     int smin, int smax, int nE);
int
build_curves_with_torsion_Z7(mpz_t f, mpmod_t n, 
			     ell_curve_t *tE, ell_point_t *tP,
			     int umin, int umax, int nE);
int
build_curves_with_torsion_Z9(mpz_t fac, mpmod_t n, ell_curve_t *tE, 
			     ell_point_t *tP, int umin, int umax, int nE);
int
build_curves_with_torsion_Z10(mpz_t fac, mpmod_t n, ell_curve_t *tE, 
			      ell_point_t *tP, int umin, int umax, int nE);
int
build_curves_with_torsion_Z2xZ8(mpz_t f, mpmod_t n, 
				ell_curve_t *tE, ell_point_t *tP,
				int umin, int umax, int nE);
int
build_curves_with_torsion_Z3xZ3_DuNa(mpmod_t n, ell_curve_t *tE, ell_point_t *tP,
				     int smin, int smax, int nE);
int
build_curves_with_torsion_Z3xZ3(mpz_t f, mpmod_t n, 
				ell_curve_t *tE, ell_point_t *tP,
				int umin, int umax, int nE);
int
build_curves_with_torsion_Z3xZ6(mpz_t f, mpmod_t n, 
				ell_curve_t *tE, ell_point_t *tP,
				int umin, int umax, int nE);
int build_curves_with_torsion(mpz_t f, mpmod_t n, ell_curve_t *tE, ell_point_t *tP, char *torsion, int smin, int smax, int nE);
int build_curves_with_torsion2(mpz_t f, mpz_t n, ell_curve_t E,  mpz_t x, mpz_t y, char *torsion, mpz_t sigma);
