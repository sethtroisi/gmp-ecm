/* When compiling the CUDA code, we do not want to include all ecm-impl.h*/
#define _DO_NOT_INCLUDE_ECM_IMPL_H

#include "cudacommon.h"
#include "ecm-gpu.h"

#include <stdio.h>


#ifndef __CUDACC__
#error "This file should only be compiled with nvcc"
#endif

/* First call to a global function initialize the device */
__global__ void Cuda_Init_Device ()
{
}

extern "C"
int
get_device_prop(int device, cudaDeviceProp *deviceProp)
{
  cudaError_t err;

  if (device!=-1)
    {
      err = cudaSetDevice(device);
      if (err != cudaSuccess)
        {
          fprintf (stderr, "GPU: Error: Could not use device %d\n", device);
          fprintf (stderr, "GPU: Error msg: %s\n", cudaGetErrorString(err));
          return 0;
        }
    }

  err = cudaGetDevice (&device);
  if (err != cudaSuccess)
    {
      fprintf (stderr, "GPU: Error: no active device.\n");
      fprintf (stderr, "GPU: Error msg: %s\n", cudaGetErrorString(err));
      return 0;
    }

  err = cudaGetDeviceProperties (deviceProp, device);
  if (err != cudaSuccess)
    {
      fprintf (stderr, "GPU: Error while getting device's properties.\n");
      fprintf (stderr, "GPU: Error msg: %s\n", cudaGetErrorString(err));
      return 0;
    }
  return 1;
}

extern "C"
int
select_and_init_GPU (int device, unsigned int *number_of_curves, int verbose)
{
  cudaDeviceProp deviceProp;

  if (device!=-1 && verbose)
    fprintf (stdout, "GPU: device %d is required.\n", device);

  if (!get_device_prop(device, &deviceProp))
    return -1;

  if (verbose)
    {
      printf ("GPU: will use device %d: %s, compute capability %d.%d, %d MPs.\n"
              "GPU: maxSharedPerBlock = %zu maxThreadsPerBlock = %d "
              "maxRegsPerBlock = %d\n", device, deviceProp.name,
              deviceProp.major, deviceProp.minor,
              deviceProp.multiProcessorCount, deviceProp.sharedMemPerBlock,
              deviceProp.maxThreadsPerBlock, deviceProp.regsPerBlock);
    }


  if (*number_of_curves == 0) /* if choose the number of curves */
    {
      /* Limited by the maximum number of threads per MP */
      unsigned int blocks_per_multiprocessor = 2;
      *number_of_curves = blocks_per_multiprocessor * deviceProp.multiProcessorCount
          * ECM_GPU_CURVES_BY_BLOCK;
    }
  else if (*number_of_curves % ECM_GPU_CURVES_BY_BLOCK != 0)
    {
      /* number_of_curves should be a multiple of ECM_GPU_CURVES_BY_BLOCK */
      *number_of_curves = (*number_of_curves / ECM_GPU_CURVES_BY_BLOCK + 1) *
                                                        ECM_GPU_CURVES_BY_BLOCK;
      if (verbose)
          fprintf(stderr, "GPU: the requested number of curves has been "
                          "modified to %u\n", *number_of_curves);
    }

  /* First call to a global function initialize the device */
  cuda_check (cudaSetDeviceFlags (cudaDeviceScheduleBlockingSync));
  Cuda_Init_Device<<<1, 1>>> ();
  cuda_check (cudaGetLastError());

  return 0;
}

void
kernel_info(const void* func, int verbose)
{
  if (verbose)
  {
    struct cudaFuncAttributes kernelAttr;
    cudaError_t err = cudaFuncGetAttributes (&kernelAttr, func);
    if (err == cudaSuccess)
      printf ("GPU: Using device code targeted for architecture compile_%d\n"
              "GPU: Ptx version is %d\nGPU: maxThreadsPerBlock = %d\n"
              "GPU: numRegsPerThread = %d sharedMemPerBlock = %zu bytes\n",
              kernelAttr.binaryVersion, kernelAttr.ptxVersion,
              kernelAttr.maxThreadsPerBlock, kernelAttr.numRegs,
              kernelAttr.sharedSizeBytes);
  }
}

