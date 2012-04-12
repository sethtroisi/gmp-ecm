#include "ecm-impl.h"

#ifndef WITH_GPU

int gpu_ecm()
{
  fprintf(stderr, "This version of libecm does not contain the GPU code.\n"
                  "You should recompile it with ./configure --enable-gpu or\n"
                  "link a version of libecm which contain the GPU code.\n");
  exit(EXIT_FAILURE);
}

#else

#include "cuda.h"

int gpu_ecm()
{
  fprintf(stderr, "gpu code is not yet available via gmp-ecm.\n"
                  "See trunk/gpu/gpuecm for gpu code.\n");
  return ECM_NO_FACTOR_FOUND;
}
#endif
