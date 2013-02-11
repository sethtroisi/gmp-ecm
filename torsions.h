void mod_div_2(mpz_t x, mpz_t n);
int mod_from_rat(mpz_t r, mpq_t q, mpz_t N);
int mod_from_rat2(mpz_t r, mpz_t num, mpz_t den, mpz_t N);
void ec_force_point(ec_curve_t E, ec_point_t P, mpz_t B, long *x0, mpz_t n);

int
build_curves_with_torsion_Z5(mpz_t f, mpmod_t n, 
			     ec_curve_t *tE, ec_point_t *tP,
			     int smin, int smax, int nE);
int
build_curves_with_torsion_Z7(mpz_t f, mpmod_t n, 
			     ec_curve_t *tE, ec_point_t *tP,
			     int umin, int umax, int nE);
int
build_curves_with_torsion_Z9(mpz_t fac, mpmod_t n, ec_curve_t *tE, 
			     ec_point_t *tP, int umin, int umax, int nE);
int
build_curves_with_torsion_Z10(mpz_t fac, mpmod_t n, ec_curve_t *tE, 
			      ec_point_t *tP, int umin, int umax, int nE);
int
build_curves_with_torsion_Z2xZ8(mpz_t f, mpmod_t n, 
				ec_curve_t *tE, ec_point_t *tP,
				int umin, int umax, int nE);
int
build_curves_with_torsion_Z3xZ3_DuNa(mpmod_t n, ec_curve_t *tE, ec_point_t *tP,
				     int smin, int smax, int nE);
int
build_curves_with_torsion_Z3xZ3(mpz_t f, mpmod_t n, 
				ec_curve_t *tE, ec_point_t *tP,
				int umin, int umax, int nE);
int
build_curves_with_torsion_Z3xZ6(mpz_t f, mpmod_t n, 
				ec_curve_t *tE, ec_point_t *tP,
				int umin, int umax, int nE);
int
build_curves_with_torsion_Z4xZ4(mpz_t f, mpmod_t n, ec_curve_t *tE,
				ec_point_t *tP,
				int smin, int smax, int nE);
int
build_curves_with_torsion_Z5xZ5(mpmod_t n, ec_curve_t *tE,
				ec_point_t *tP,
				int smin, int smax, int nE);
int
build_curves_with_torsion_Z2xZ10(mpz_t f, mpmod_t n, ec_curve_t *tE,
				 ec_point_t *tP,
				 int smin, int smax, int nE);
int
build_curves_with_torsion_Z2xZ12(mpz_t f, mpmod_t n, ec_curve_t *tE,
				 ec_point_t *tP,
				 int smin, int smax, int nE);


