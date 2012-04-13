#include "ecm-impl.h"
#include "ecm-gpu.h"

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

extern int select_and_init_GPU (int, int, FILE*);

int gpu_ecm()
{
  int number_of_curves = select_and_init_GPU (-1, 0, stdout);
  fprintf (stdout,"number_of_curves=%d\n",number_of_curves);
  fprintf(stderr, "gpu code is not yet available via gmp-ecm.\n"
                  "See trunk/gpu/gpuecm for gpu code.\n");
  return ECM_NO_FACTOR_FOUND;
}
#endif
