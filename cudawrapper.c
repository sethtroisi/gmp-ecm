#include "ecm-ecm.h"
#include "cuda.h"

int gpu_ecm_factor()
{
  fprintf(stderr, "gpu code is not yet available via gmp-ecm.\n"
                  "See trunk/gpu/gpuecm for gpu code.\n");
  return ECM_NO_FACTOR_FOUND;
}
