#include "def.h"
#include "cudakernel.h"

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


extern "C" 
int select_GPU (int device, int number_of_curves, FILE *OUTPUT_VERBOSE)
{
  cudaDeviceProp deviceProp;
  cudaError_t err;
        
  fprintf(OUTPUT_VERBOSE, "#Compiled for a NVIDIA GPU with " 
          "compute capability %d.%d.\n", MAJOR, MINOR);

  if (device!=-1)
  {
    fprintf(OUTPUT_VERBOSE,"#Device %d is required.\n",device);

    err= cudaSetDevice(device);
    if (err != cudaSuccess)
    {
      fprintf(stderr, "Error: Could not use device %d\n",device);
      fprintf(stderr, "Error msg: %s\n", cudaGetErrorString(err));
      exit(EXIT_FAILURE);
    }
  }
  
  err = cudaGetDevice (&device);
  if (err != cudaSuccess)
  {
    fprintf(stderr, "Error: no active device\n");
    fprintf(stderr, "Error msg: %s\n", cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }

  err = cudaGetDeviceProperties (&deviceProp, device);
  if (err != cudaSuccess)
  {
    fprintf(stderr, "Error while getting device's properties\n");
    fprintf(stderr, "Error msg: %s\n", cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }

  int minor = deviceProp.minor;
  int major = deviceProp.major;
  int MPcount = deviceProp.multiProcessorCount;

  if (10 * major + minor < 10 * MAJOR + MINOR)
  {
    fprintf(stderr, "Error: Device %d have a compute capability of %d.%d " 
                    "(required %d.%d).\n", device, major, minor, MAJOR, MINOR);
    exit(EXIT_FAILURE);
  }

  fprintf(OUTPUT_VERBOSE, "#Will use device %d : %s, compute capability %d.%d, "
          "%d MPs.\n", device, deviceProp.name, major, minor, MPcount);


  cudaSetDeviceFlags(cudaDeviceScheduleAuto); 
  //cudaSetDeviceFlags(cudaDeviceScheduleYield); 
  //cudaSetDeviceFlags(cudaDeviceScheduleSpin); //the other make performance
  //worse

  //number_of_curves should be a multiple of CURVES_BY_BLOCK
  number_of_curves=(number_of_curves/CURVES_BY_BLOCK)*CURVES_BY_BLOCK;
  if (number_of_curves==0)
    number_of_curves = MPcount * CURVES_BY_MP;

  return number_of_curves;
}

extern "C"
void cuda_Main (biguint_t h_N, biguint_t h_3N, biguint_t h_M, digit_t h_invN, 
                    biguint_t *h_xarray, biguint_t *h_zarray, 
                    biguint_t *h_x2array, biguint_t *h_z2array, mpz_t s,
                    unsigned int firstinvd, unsigned int number_of_curves, 
                    FILE *OUTPUT_VERBOSE, FILE *OUTPUT_VVERBOSE) 
{ 
  size_t j;
  biguint_t *d_xA, *d_zA, *d_xB, *d_zB;

  size_t array_size = sizeof(biguint_t) * number_of_curves;

  dim3 dimBlock (NB_DIGITS, CURVES_BY_BLOCK);
  dim3 dimGrid (number_of_curves/CURVES_BY_BLOCK);

  fprintf(OUTPUT_VVERBOSE, "Block: %ux%ux%u Grid: %ux%ux%u\n", dimBlock.x, 
                      dimBlock.y, dimBlock.z, dimGrid.x, dimGrid.y, dimGrid.z);

  cudaMalloc (&d_xA, array_size);
  cudaMalloc (&d_zA, array_size);
  cudaMalloc (&d_xB, array_size);
  cudaMalloc (&d_zB, array_size);

  /* Copy into the gpu memory */
  cudaMemcpyToSymbol (d_invNcst, (void *) &h_invN, sizeof(digit_t));
  cudaMemcpyToSymbol (d_Ncst, (void *) h_N, sizeof(biguint_t));
  cudaMemcpyToSymbol (d_3Ncst, (void *) h_3N, sizeof(biguint_t));
  cudaMemcpyToSymbol (d_Mcst, (void *) h_M, sizeof(biguint_t));

  cudaMemcpyHtoD (d_xA, h_xarray, array_size);
  cudaMemcpyHtoD (d_zA, h_zarray, array_size);
  cudaMemcpyHtoD (d_xB, h_x2array, array_size);
  cudaMemcpyHtoD (d_zB, h_z2array, array_size);

  /* Double-and-add loop: it calls the GPU for each bits of s */
  for (j = mpz_sizeinbase (s, 2) - 1; j>0; j-- )
  {
    if (mpz_tstbit (s, j-1) == 1)
      Cuda_Ell_DblAdd<<<dimGrid,dimBlock>>>(d_xB, d_zB, d_xA, d_zA, firstinvd);
    else
      Cuda_Ell_DblAdd<<<dimGrid,dimBlock>>>(d_xA, d_zA, d_xB, d_zB, firstinvd);

    //maybe only for debug mode??
    errCheck (cudaGetLastError()); 
  }

  /* Get the results back from device memory */
  cudaMemcpyDtoH (h_xarray, d_xA, array_size);
  cudaMemcpyDtoH (h_zarray, d_zA, array_size);

  cudaFree ((void *) d_xA);
  cudaFree ((void *) d_zA);
  cudaFree ((void *) d_xB);
  cudaFree ((void *) d_zB);
}



/***************/
/* Device code */
/***************/


#define __add_cc(r,a,b) __asm__("add.cc.u32 %0,%1, %2;" :"=r"(r):"r"(a),"r"(b)) 
#define __addc_cc(r,a,b) __asm__("addc.cc.u32 %0,%1, %2;":"=r"(r):"r"(a),"r"(b)) 

#define __sub_cc(r,a,b) __asm__("sub.cc.u32 %0,%1, %2;" :"=r"(r):"r"(a),"r"(b)) 

#define __addcy(carry) __asm__("addc.s32 %0, 0, 0;" :"=r"(carry)) 
#define __addcy2(carry) __asm__("addc.s32 %0, %0, 0;" :"+r"(carry)) 

#define __subcy(carry) __asm__("subc.s32 %0, 0, 0;" :"=r"(carry)) 
#define __subcy2(carry) __asm__("subc.s32 %0, %0, 0;" :"+r"(carry)) 

#define __mul(h,l,a,b) __asm__("mul.hi.u32 %0,%2,%3;\n\t" "mul.lo.u32 %1,%2,%3;"\
                                            : "=r"(h), "=r"(l) : "r"(a), "r"(b))

#define __mad_lo(r,a,b,c) __asm__("mad.lo.u32 %0,%1,%2,%3;" \
                                            : "=r"(r) : "r"(a), "r"(b), "r"(c))
#define __mad_hi(r,a,b,c) __asm__("mad.hi.u32 %0,%1,%2,%3;" \
                                            : "=r"(r) : "r"(a), "r"(b), "r"(c))




//  (A > B)?, returns 1(true), -1(false) or 0(a=b) 
//Assume A and B are normalize (no carry or borrow)
__device__ int Cuda_Cmp(const biguint_t A, const biguint_t B)
{
  int i;
  for (i = NB_DIGITS-1;i>=0;i--)
  {
    if (A[i] > B[i])
      return 1;
    else if (A[i] < B[i])
      return -1;
  }
  return 0;
}

//Assume cy[threadIdx.x] = 0,+/-1
__device__ void Cuda_Normalize(biguint_t A,dbigint_t cy)
{
  carry_t cytemp;
  cytemp = cy[threadIdx.x];
  cy[threadIdx.x]=0;
  int tmp=threadIdx.x+1 % NB_DIGITS;

  if (cytemp==1)
  {
    A[tmp]++;
    if (A[tmp]==0)
      cy[tmp]=cytemp;
  }
  else if (cytemp==-1) 
  {
    if (A[tmp]==0)
      cy[tmp]=cytemp;
    A[tmp]--;
  } 
}

__device__ void Cuda_Fully_Normalize(biguint_t A,dbigint_t cy)
{
  
  do
  {
    Cuda_Normalize(A,cy);
  }while(__any(cy[threadIdx.x])!=0);
  
}

__device__ void Cuda_Add 
(biguint_t r, dbigint_t cy ,const biguint_t a, const biguint_t b)
{
  __add_cc(r[threadIdx.x],a[threadIdx.x],b[threadIdx.x]);
  __addcy(cy[threadIdx.x]);
}

__device__ void Cuda_Subc 
(biguint_t r, dbigint_t cy, const biguint_t a, const biguint_t b)
{
  __sub_cc(r[threadIdx.x],a[threadIdx.x],b[threadIdx.x]);
  __subcy2(cy[threadIdx.x]);
}

__device__ void Cuda_Sub 
(biguint_t r, dbigint_t cy, const biguint_t a, const biguint_t b)
{
  __sub_cc(r[threadIdx.x],a[threadIdx.x],b[threadIdx.x]);
  __subcy(cy[threadIdx.x]);
}

/* Compute Rmod <- A + B */ 
/* Input: 0 <= A, B < 3*N */ 
/* Ouput: 0 <= Rmod < 6*N */ 
__device__ void Cuda_Add_mod
(biguint_t Rmod, dbigint_t cy, const biguint_t A,const biguint_t B)
{
  Cuda_Add(Rmod, cy, A, B);
  Cuda_Fully_Normalize(Rmod, cy); 
}

/* Compute Rmod <- Rmod + B */ 
/* Input: 0 <= Rmod, B < 3*N */ 
/* (except when it follows Cuda_Mulint_mod, 0 <= Rmod < 3*N, 0 < B < 7*N ) */ 
/* Ouput: 0 <= Rmod < 6*N */ 
/* (except when it follows Cuda_Mulint_mod, 0 <= Rmod < 10*N) */ 
__device__ void Cuda_Add_mod
(biguint_t Rmod, dbigint_t cy, const biguint_t A)
{
  Cuda_Add(Rmod, cy, Rmod, A);
  Cuda_Fully_Normalize(Rmod, cy);
}

/* Compute Rmod <- Rmod - B */ 
/* Input: 0 <= Rmod, B < 3*N */ 
/* Ouput: 0 <= Rmod < 6*N */ 
__device__ void Cuda_Sub_mod 
(biguint_t Rmod, dbigint_t cy, const biguint_t A)
{
  Cuda_Add (Rmod, cy, Rmod, d_3Ncst);
  Cuda_Subc (Rmod, cy, Rmod, A);
  Cuda_Fully_Normalize(Rmod, cy); 
}

/* Perform one step of REDC */ 
__device__ void Cuda_Mulmod_step
(dbiguint_t r, dbigint_t cy, unsigned int a, unsigned int b)
{
  digit_t h,l;
  int tmp;
  __mul(h,l,a,b);
  __add_cc(r[threadIdx.x],r[threadIdx.x],l);
  __addc_cc(r[threadIdx.x+1],r[threadIdx.x+1],h);
  __addcy2(cy[threadIdx.x+1]);


  __mul(h, l, d_invNcst*r[0], d_Ncst[threadIdx.x]);
  __add_cc(r[threadIdx.x],r[threadIdx.x],l);
  __addc_cc(r[threadIdx.x+1],r[threadIdx.x+1],h);
  __addcy2(cy[threadIdx.x+1]);
 
  //make one round of normalize + a right shift
  __add_cc(r[threadIdx.x],r[threadIdx.x+1],cy[threadIdx.x]);
  tmp=(threadIdx.x==NB_DIGITS-1)?cy[threadIdx.x+1]:0;
  __asm__("addc.s32 %0,%1, 0;" :"=r"(cy[threadIdx.x]): "r"(tmp)); 

  if (threadIdx.x==0)
  {
    cy[NB_DIGITS]=0;
    r[NB_DIGITS]=0;
  }
}

/* Compute r <- 2*a */ 
/* Input: 0 <= a < 3*N */ 
/* Ouput: 0 <= r < 3*N */ 
__device__ void Cuda_Dbl_mod
(biguint_t r, dbigint_t cy, biguint_t a)
{
  __add_cc(r[threadIdx.x],a[threadIdx.x],a[threadIdx.x]);
  __addcy2(r[(threadIdx.x+1)%NB_DIGITS]);

  if (Cuda_Cmp (r, d_3Ncst) >= 0) 
  {
    Cuda_Sub (r, cy, r, d_3Ncst); 
    Cuda_Fully_Normalize(r, cy);  
  }
}


/* Compute r <- A*b */ 
/* Input: 0 < b < 2^SIZE_DIGIT, 0 <= A < 6*N */ 
/* Ouput: 0 <= r < 7*N */ 
__device__ void Cuda_Mulint_mod
(dbiguint_t r, dbigint_t cy, biguint_t A, unsigned int b)
{
  digit_t h,l;
  __mul(h, r[threadIdx.x], A[threadIdx.x], b);
  __add_cc(r[threadIdx.x+1], r[threadIdx.x+1], h);
  __addcy(cy[threadIdx.x+1]);

  //h*2^32+l =A[i]*B[threadIDx.x]
  __mul(h, l, d_invNcst*r[0], d_Ncst[threadIdx.x]);
  __add_cc(r[threadIdx.x], r[threadIdx.x],l);
  __addc_cc(r[threadIdx.x+1], r[threadIdx.x+1],h);
  __addcy2(cy[threadIdx.x+1]);

  __add_cc(r[threadIdx.x], r[threadIdx.x+1], cy[threadIdx.x]);
  __addcy(cy[threadIdx.x]);
  if (threadIdx.x==0)
  {
    cy[NB_DIGITS]=0;
    r[NB_DIGITS]=0;
  }
  
  Cuda_Fully_Normalize(r,cy); 
}

/* Compute r <- A*B */ 
/* Input: 0 <= A, B < 6*N */
/* (except when it follows Cuda_Mulint_mod, 0 <= A < 6*N, 0 < B < 10*N ) */ 
/* Ouput: 0 <= r < 3*N */ 
__device__ void Cuda_Mul_mod 
(biguint_t mul, dbigint_t cy, const biguint_t A, const biguint_t B, dbiguint_t r)
{

  int i;
  digit_t temp=A[threadIdx.x];

  r[threadIdx.x]=0;
  
  for (i=0;i<NB_DIGITS;i++)
    Cuda_Mulmod_step(r,cy,temp,B[i]);

  
  Cuda_Fully_Normalize(r,cy);
  mul[threadIdx.x]=r[threadIdx.x];
}

__device__ void 
Cuda_Square_mod (biguint_t mul, dbigint_t cy, const biguint_t A, dbiguint_t r)
{
  Cuda_Mul_mod(mul,cy,A,A,r);
}

/* 
  Compute silmutaneously:
  (xarg : zarg ) <- [2](xarg : zarg) 
  (xarg2 : zarg2 ) <- (xarg : zarg) + (xarg2 : zarg2) 
*/
__global__ void 
Cuda_Ell_DblAdd (biguint_t *xarg, biguint_t *zarg, biguint_t *x2arg, 
                                       biguint_t *z2arg, unsigned int firstinvd)
{
  __shared__ VOL digit_t b_temp_r[CURVES_BY_BLOCK][NB_DIGITS+1];
  __shared__ VOL carry_t b_cy[CURVES_BY_BLOCK][NB_DIGITS+1]; 

  __shared__ VOL digit_t b_t[CURVES_BY_BLOCK][NB_DIGITS];
  __shared__ VOL digit_t b_u[CURVES_BY_BLOCK][NB_DIGITS];
  __shared__ VOL digit_t b_v[CURVES_BY_BLOCK][NB_DIGITS];
  __shared__ VOL digit_t b_w[CURVES_BY_BLOCK][NB_DIGITS];
  
  volatile unsigned int idx1=blockIdx.x*blockDim.y+threadIdx.y;
  //volatile unsigned int t1=threadIdx.x+1;
  //volatile unsigned int t2=threadIdx.x+NB_DIGITS;
  
  VOL digit_t *t=b_t[threadIdx.y];
  VOL digit_t *u=b_u[threadIdx.y];
  VOL digit_t *v=b_v[threadIdx.y];
  VOL digit_t *w=b_w[threadIdx.y];
  VOL digit_t *temp_r=b_temp_r[threadIdx.y];
  VOL carry_t *cy=b_cy[threadIdx.y];

  //init
  b_cy[threadIdx.y][threadIdx.x]=0; 
  if (threadIdx.x==0)
    b_cy[threadIdx.y][NB_DIGITS]=0; 

  v[threadIdx.x]=x2arg[idx1][threadIdx.x];
  w[threadIdx.x]=z2arg[idx1][threadIdx.x];
  temp_r[threadIdx.x]=zarg[idx1][threadIdx.x];
  u[threadIdx.x]=xarg[idx1][threadIdx.x];

  Cuda_Add_mod(t, cy, v, w);      /* C=x2+z2 */
  Cuda_Sub_mod(v, cy, w);         /* D=x2-z2 */
  Cuda_Add_mod(w, cy, u, temp_r); /* A=x+z */
  Cuda_Sub_mod(u, cy, temp_r);    /* B=x-z */

  Cuda_Mul_mod(t, cy, t, u, temp_r); /* CB=C*B=(xq+zq)(xp-zp) */
  Cuda_Mul_mod(v, cy, v, w, temp_r); /* DA=D*A=(xq-zq)(xp+zp) */

  Cuda_Square_mod(w, cy, w, temp_r); /* AA=A^2 */
  Cuda_Square_mod(u, cy, u, temp_r); /* BB=B^2 */

  Cuda_Mul_mod(temp_r, cy, u, w, temp_r); /* x2=AA*BB */
  xarg[idx1][threadIdx.x]=temp_r[threadIdx.x];

  Cuda_Sub_mod (w, cy, u); /* C= AA-BB */
  Cuda_Mulint_mod (temp_r, cy, w, idx1 + firstinvd); /* d*C */ 
  Cuda_Add_mod (u, cy, temp_r); /* BB+d*C */
 
  Cuda_Mul_mod (w, cy, w, u,temp_r); /* z2=C*(BB+d*C) */
  zarg[idx1][threadIdx.x]=w[threadIdx.x];
 
  Cuda_Add_mod(w, cy, v, t); /* DA+CB mod N */
  Cuda_Sub_mod(v, cy, t); /* DA-CB mod N */

  Cuda_Square_mod(w, cy, w, temp_r); /* (DA+CB)^2 mod N */
  Cuda_Square_mod(v, cy, v, temp_r); /* (DA-CB)^2 mod N */

  /* z0=1 so there is nothing to compute for z0*(DA+CB)^2 */
  Cuda_Dbl_mod(u, cy, v); /* x0=2 x0*(DA-CB)^2 */
  
  x2arg[idx1][threadIdx.x]=w[threadIdx.x];
  z2arg[idx1][threadIdx.x]=u[threadIdx.x];
}

