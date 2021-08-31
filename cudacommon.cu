/* When compiling the CUDA code, we do not want to include all ecm-impl.h*/
#define _DO_NOT_INCLUDE_ECM_IMPL_H

#include "cudacommon.h"

#include <stdio.h>


#ifndef __CUDACC__
#error "This file should only be compiled with nvcc"
#endif

/* First call to a global function initialize the device */
__global__ void Cuda_Init_Device ()
{
}

/* Given the compute compatibility (as major.minor), return the number of block
 * to be run on one multiprocessor. */
static unsigned int
getNumberOfBlockPerMultiProcessor (int major, int minor)
{
  /* For 2.0 and 2.1, limited by the maximum number of threads per MP and the
   * number of available registrer (need 23 registers per threads).
   */
  if (major == 2)
    return 1;
  /* For 3.0, 3.2, 3.5 and 3.7 limited by the maximum number of threads per MP.
   */
  else if (major == 3)
    return 2;
  /* For 5.0, 5.2, and 5.3 limited by the maximum number of threads per MP. */
  else if (major == 5)
    return 2;
  /* We assume that for newer compute capability the properties of the GPU won't
   * decrease.
   */
  else
    return 2;
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
select_and_init_GPU (int device, unsigned int *number_of_curves, int verbose, int schedule)
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
      unsigned int n, m = ECM_GPU_CURVES_BY_BLOCK;
      n = getNumberOfBlockPerMultiProcessor (deviceProp.major, deviceProp.minor);
      *number_of_curves = n * deviceProp.multiProcessorCount * m;
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
  if (schedule == 1)
    {
      cuda_check (cudaSetDeviceFlags (cudaDeviceScheduleBlockingSync));
    }
  else
    {
      cuda_check (cudaSetDeviceFlags (cudaDeviceScheduleYield));
    }
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

