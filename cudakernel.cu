/* When compiling the CUDA code, we do not want to include all ecm-impl.h*/
#define _DO_NOT_INCLUDE_ECM_IMPL_H
#include "ecm-gpu.h"
#include <gmp.h>
#include "cudakernel.h"

#ifndef __CUDACC__
#error "This file should only be compiled with nvcc"
#endif

__constant__ __device__ digit_t d_invNcst;
__device__ biguint_t d_Ncst;
__device__ biguint_t d_3Ncst;
__device__ biguint_t d_Mcst;


#define errCheck(err) cuda_errCheck (err, __FILE__, __LINE__)
#define cudaMalloc(d, size) errCheck (cudaMalloc (d, size))
#define cudaMemcpyHtoD(d, h, size) errCheck (cudaMemcpy ((void *) d, \
                                    (void *) h, size, cudaMemcpyHostToDevice))
#define cudaMemcpyDtoH(h, d, size) errCheck (cudaMemcpy ((void *) h, \
                                    (void *) d, size, cudaMemcpyDeviceToHost))
#define cudaMemcpyCst(d, h, size) errCheck (cudaMemcpyToSymbol (d, h, size))


/******************************/
/* Host code handling the GPU */
/******************************/

inline void cuda_errCheck (cudaError err, const char *file, const int line)
{
  if( err != cudaSuccess ) 
  {
    fprintf(stderr, "%s(%i) : Error cuda : %s.\n",
              file, line, cudaGetErrorString( err) );
    exit(EXIT_FAILURE);
  }
}

/* First call to a global function initialize the device */
__global__ void Cuda_Init_Device ()
{
}

/* Given the compute compatibility (as major.minor), return the number of block
 * to be run on one multiprocessor. */
extern "C"
unsigned int
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
select_and_init_GPU (int device, unsigned int *number_of_curves, int verbose)
{
  cudaDeviceProp deviceProp;
  cudaError_t err;
        
  if (device!=-1)
    {
      if (verbose)
          fprintf (stdout, "GPU: device %d is required.\n", device);

      err = cudaSetDevice(device);
      if (err != cudaSuccess)
        {
          fprintf (stderr, "GPU: Error: Could not use device %d\n", device);
          fprintf (stderr, "GPU: Error msg: %s\n", cudaGetErrorString(err));
          return -1;
        }
    }
  
  err = cudaGetDevice (&device);
  if (err != cudaSuccess)
    {
      fprintf (stderr, "GPU: Error: no active device.\n");
      fprintf (stderr, "GPU: Error msg: %s\n", cudaGetErrorString(err));
      return -1;
    }

  err = cudaGetDeviceProperties (&deviceProp, device);
  if (err != cudaSuccess)
    {
      fprintf (stderr, "GPU: Error while getting device's properties.\n");
      fprintf (stderr, "GPU: Error msg: %s\n", cudaGetErrorString(err));
      return -1;
    }

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
  errCheck (cudaSetDeviceFlags (cudaDeviceScheduleYield)); 
  Cuda_Init_Device<<<1, 1>>> ();
  errCheck (cudaGetLastError()); 

  if (verbose)
  {
    struct cudaFuncAttributes kernelAttr;
    err = cudaFuncGetAttributes (&kernelAttr, Cuda_Ell_DblAdd);
    if (err == cudaSuccess)
    {
      printf ("GPU: Using device code targeted for architecture compile_%d\n"
              "GPU: Ptx version is %d\nGPU: maxThreadsPerBlock = %d\n"
              "GPU: numRegsPerThread = %d sharedMemPerBlock = %zu bytes\n",
              kernelAttr.binaryVersion, kernelAttr.ptxVersion,
              kernelAttr.maxThreadsPerBlock, kernelAttr.numRegs,
              kernelAttr.sharedSizeBytes);
    }
  }

  return 0;
}

extern "C"
float cuda_Main (biguint_t h_N, biguint_t h_3N, biguint_t h_M, digit_t h_invN,
                 biguint_t *h_xarray, biguint_t *h_zarray,
                 biguint_t *h_x2array, biguint_t *h_z2array, mpz_t s,
                 unsigned int firstinvd, unsigned int number_of_curves,
                 int verbose)
{
  cudaEvent_t start, stop;
  cudaEventCreate (&start);
  cudaEventCreate (&stop);
  cudaEventRecord (start, 0);

  size_t j;
  int i;
  float elltime = 0.0;
  biguint_t *d_xA, *d_zA, *d_xB, *d_zB;

#define MAXEVENTS 2
#define DEPTH_EVENT 32
  cudaEvent_t event[MAXEVENTS];   // Space for some cuda Event Handles
  long nEventsRecorded = 0;   // Remember how many events are recorded
  long eventrecordix = 0;     // Remember index of next event to record
  long eventsyncix;       // Remember index of oldest recorded event

  size_t array_size = sizeof(biguint_t) * number_of_curves;

  dim3 dimBlock (ECM_GPU_NB_DIGITS, ECM_GPU_CURVES_BY_BLOCK);
  dim3 dimGrid (number_of_curves/ ECM_GPU_CURVES_BY_BLOCK);
  if (verbose)
    {
      fprintf(stdout, "GPU: Block: %ux%ux%u Grid: %ux%ux%u "
              "(%d parallel curves)\n", dimBlock.x, dimBlock.y, dimBlock.z,
              dimGrid.x, dimGrid.y, dimGrid.z, number_of_curves);
    }

  /* Create a pair of events to pace ourselves */
  for (i=0; i<MAXEVENTS; i++)
    errCheck (cudaEventCreateWithFlags (&event[i], 
                              cudaEventBlockingSync|cudaEventDisableTiming));

  cudaMalloc (&d_xA, array_size);
  cudaMalloc (&d_zA, array_size);
  cudaMalloc (&d_xB, array_size);
  cudaMalloc (&d_zB, array_size);

  /* Copy into the gpu memory */
  cudaMemcpyCst (d_invNcst, (void *) &h_invN, sizeof(digit_t));
  cudaMemcpyCst (d_Ncst, (void *) h_N, sizeof(biguint_t));
  cudaMemcpyCst (d_3Ncst, (void *) h_3N, sizeof(biguint_t));
  cudaMemcpyCst (d_Mcst, (void *) h_M, sizeof(biguint_t));

  cudaMemcpyHtoD (d_xA, h_xarray, array_size);
  cudaMemcpyHtoD (d_zA, h_zarray, array_size);
  cudaMemcpyHtoD (d_xB, h_x2array, array_size);
  cudaMemcpyHtoD (d_zB, h_z2array, array_size);

#ifdef PRINT_REMAINING_ITER
      unsigned int jmod = 100000000;
#endif

  /* Double-and-add loop: it calls the GPU for each bits of s */
  for (j = mpz_sizeinbase (s, 2) - 1; j>0; j-- )
  {
    if (mpz_tstbit (s, j-1) == 1)
      Cuda_Ell_DblAdd<<<dimGrid,dimBlock>>>(d_xB, d_zB, d_xA, d_zA, firstinvd);
    else
      Cuda_Ell_DblAdd<<<dimGrid,dimBlock>>>(d_xA, d_zA, d_xB, d_zB, firstinvd);

    /* Pace entry of events. Less overhead to enter an event every few    */
    /* iterations. But, if you exceed the depth of NVIDIA's kernel queue, */
    /* it will busy-loop!                                                 */
    /* Enter an event every DEPTH_EVENT iteration */
    if (j % DEPTH_EVENT == 0)  
    {
      cudaEventRecord(event[eventrecordix]); 
      if (nEventsRecorded == 0)     
        eventsyncix = eventrecordix; 
      nEventsRecorded += 1;          
      eventrecordix = (eventrecordix+1)%MAXEVENTS;  
    }

    if (nEventsRecorded == MAXEVENTS) 
    {
      cudaEventSynchronize(event[eventsyncix]);  
      nEventsRecorded -= 1;   
      eventsyncix = (eventsyncix+1)%MAXEVENTS; 
    }

#ifdef PRINT_REMAINING_ITER
    if (j < 100000000) jmod = 10000000;
    if (j < 10000000)  jmod =  1000000;
    if (j < 1000000)   jmod =   100000;
    if (j < 100000)    jmod =    10000;
    if (j % jmod == 0)
      printf("%lu iterations to go\n", j);
#endif
  }

  /* If an error occurs during the kernel calls in the loop */
  errCheck (cudaGetLastError()); 

  /* Await for last recorded events */
  while (nEventsRecorded != 0) 
  {
    cudaEventSynchronize(event[eventsyncix]); 
    nEventsRecorded -= 1;          
    eventsyncix = (eventsyncix+1)%MAXEVENTS; 
  }

  /* Get the results back from device memory */
  cudaMemcpyDtoH (h_xarray, d_xA, array_size);
  cudaMemcpyDtoH (h_zarray, d_zA, array_size);

  /* Clean up our events and our stream handle */
  for (i=0; i<MAXEVENTS; i++)
    errCheck (cudaEventDestroy(event[i]));


  cudaFree ((void *) d_xA);
  cudaFree ((void *) d_zA);
  cudaFree ((void *) d_xB);
  cudaFree ((void *) d_zB);

  cudaEventRecord (stop, 0);
  cudaEventSynchronize (stop);

  cudaEventElapsedTime (&elltime, start, stop);

  errCheck (cudaEventDestroy (start));
  errCheck (cudaEventDestroy (stop));

  return elltime;
}



/***************/
/* Device code */
/***************/

#if defined(_MSC_VER)
#  define ASM asm volatile
#else
#  define ASM asm __volatile__
#endif

#define __add_cc(r,a,b) ASM ("add.cc.u32 %0, %1, %2;": "=r"(r): "r"(a), "r"(b)) 
#define __addc_cc(r,a,b) ASM ("addc.cc.u32 %0, %1, %2;": "=r"(r): "r"(a), "r"(b))
#define __sub_cc(r,a,b) ASM ("sub.cc.u32 %0, %1, %2;": "=r"(r): "r"(a), "r"(b)) 

#define __addcy(carry) ASM ("addc.s32 %0, 0, 0;": "=r"(carry)) 
#define __addcy2(carry) ASM ("addc.cc.s32 %0, %0, 0;": "+r"(carry)) 

#define __subcy(carry) ASM ("subc.s32 %0, 0, 0;": "=r"(carry)) 
#define __subcy2(carry) ASM ("subc.s32 %0, %0, 0;": "+r"(carry)) 

#define __mul_lo(r,a,b) ASM("mul.lo.u32 %0, %1, %2;": "=r"(r): "r"(a),"r"(b)) 
#define __mul_hi(r,a,b) ASM("mul.hi.u32 %0, %1, %2;": "=r"(r): "r"(a),"r"(b)) 
#define __mad_lo_cc(r,a,b) ASM("mad.lo.cc.u32 %0, %1, %2, %0;":\
                                                      "+r"(r): "r"(a),"r"(b)) 
#define __madc_hi_cc(r,a,b) ASM("madc.hi.cc.u32 %0, %1, %2, %0;":\
                                                  "+r"(r):"r"(a),"r"(b)) 

#ifdef __CUDA_ARCH__
  #if __CUDA_ARCH__ >= 200
    #include "cudakernel_default.cu"
  #else
    #error "Unsupported architecture"
  #endif
#endif
