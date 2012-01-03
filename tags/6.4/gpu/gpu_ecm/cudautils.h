inline void cuda_errCheck (cudaError err, const char *file, const int line);
void cuda_Main ( biguint_t h_N, biguint_t h_invmod, biguint_t *h_xarray, 
                     biguint_t *h_zarray, biguint_t *h_x2array, 
                     biguint_t *h_z2array, mpz_t s, unsigned int firstinvd, 
                    unsigned int number_of_curves, FILE *OUTPUT_VERBOSE,
                    FILE *OUTPUT_VVERBOSE);
