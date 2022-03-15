#ifndef _CUDACOMMON_H
#define _CUDACOMMON_H 1

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <stdio.h>
#ifdef _MSC_VER
#include <stdint.h>
#endif

#ifdef __cplusplus
/* cpp + CUDA only code */

#define CUDA_CHECK(action) cuda_check(action, #action, __FILE__, __LINE__)

inline void cuda_check(cudaError_t status, const char *action=NULL, const char *file=NULL, int32_t line=0) {
  if (status != cudaSuccess) {
    fprintf (stderr, "CUDA error (%d) occurred: %s\n", status, cudaGetErrorString(status));
    if (action!=NULL)
      fprintf (stderr, "While running %s   (file %s, line %d)\n", action, file, line);
    exit(EXIT_FAILURE);
  }
}


void kernel_info(const void* func, int verbose);
#endif


#ifdef __cplusplus
extern "C" {
#endif

int get_device_prop(int device, struct cudaDeviceProp *deviceProp);
int select_and_init_GPU (int, unsigned int*, int, int);

#ifdef __cplusplus
}
#endif

#endif /* _CUDACOMMON_H */
