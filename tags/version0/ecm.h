double   getprime       (double);
unsigned int nb_digits  (mpz_t);
int          stage1     (mpz_t, mpz_t, double);
int          pm1        (mpz_t, mpz_t, double, double, unsigned int);
int          ecm        (mpz_t, mpz_t, double, double, unsigned int);
unsigned int phi        (unsigned int);
unsigned int bestD      (double);
double       block_size (unsigned int);
unsigned int gcd        (unsigned int, unsigned int);
int          cputime    (void);

/* stage2.c */
int          rootsF     (mpz_t *, unsigned int, mpz_t, mpz_t, mpz_t, mpz_t,
			 mpz_t);
void         rootsG     (mpz_t *, unsigned int, mpz_t, mpz_t, mpz_t, mpz_t,
			 mpz_t);
int          stage2     (mpz_t, mpz_t, double, unsigned int);

/* poly.c */
mpz_t *      init_poly  (unsigned int);
void         clear_poly (mpz_t *, unsigned int);
void         print_list (mpz_t *, unsigned int);
void         print_poly (mpz_t *, unsigned int, int);
void         polyset    (mpz_t *, mpz_t *, unsigned int);
void         polymod    (mpz_t *, mpz_t *, unsigned int, mpz_t);
void         polyadd    (mpz_t *, mpz_t *, mpz_t *, unsigned int);
void         polysub    (mpz_t *, mpz_t *, mpz_t *, unsigned int);
void         poly_zero  (mpz_t *, unsigned int);
int          poly_iszero(mpz_t *, unsigned int);
void         karatsuba  (mpz_t *, mpz_t *, mpz_t *, unsigned int, mpz_t *);
void         polymul    (mpz_t *, mpz_t *, unsigned int, mpz_t *, unsigned int,
			 mpz_t *, int);
void         polymulmod (mpz_t *, mpz_t *, mpz_t *, unsigned int, mpz_t *,
			 mpz_t);
void         buildG     (mpz_t *, unsigned int, mpz_t *, int, mpz_t, char);
void  RecursiveDivision (mpz_t *, mpz_t *, mpz_t *, unsigned int, mpz_t *,
			 mpz_t);
void         Div3by2    (mpz_t *, mpz_t *, mpz_t *, unsigned int, mpz_t *,
			 mpz_t);
int          polymod1   (mpz_t, mpz_t *, mpz_t *, unsigned int, mpz_t, mpz_t*);
int          hgcd_naive (mpz_t, mpz_t *, mpz_t *, mpz_t *, mpz_t *,
                         mpz_t *, mpz_t *, unsigned int, mpz_t, mpz_t *, int);
int          hgcd       (mpz_t, mpz_t *, mpz_t *, mpz_t *, mpz_t *,
                         mpz_t *, mpz_t *, unsigned int, mpz_t, mpz_t *, int);
int          polygcd    (mpz_t, mpz_t *, mpz_t *, unsigned int, mpz_t, mpz_t*);
