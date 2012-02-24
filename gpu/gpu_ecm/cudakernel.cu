#include "def.h"
#include "cudakernel.h"

__constant__ __device__ biguint_t d_invNcst;
__device__ biguint_t d_Ncst;


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
void cuda_Main (biguint_t h_N, biguint_t h_invN, biguint_t *h_xarray, 
                    biguint_t *h_zarray, biguint_t *h_x2array, 
                    biguint_t *h_z2array, mpz_t s, unsigned int firstinvd, 
                    unsigned int number_of_curves, FILE *OUTPUT_VERBOSE,
                    FILE *OUTPUT_VVERBOSE)
{
  cudaError_t err;
  
  biguint_t *d_xA;
  biguint_t *d_zA;
  biguint_t *d_xB;
  biguint_t *d_zB;

  size_t array_size = sizeof(biguint_t) * number_of_curves;

  dim3 dimBlock(NB_DIGITS,CURVES_BY_BLOCK);
  dim3 dimGrid(number_of_curves/CURVES_BY_BLOCK);

  fprintf(OUTPUT_VVERBOSE, "Block: %ux%ux%u Grid: %ux%ux%u\n", dimBlock.x, 
                      dimBlock.y, dimBlock.z, dimGrid.x, dimGrid.y, dimGrid.z);

  errCheck( cudaMalloc(&d_xA, array_size) );
  errCheck( cudaMalloc(&d_zA, array_size) );
  errCheck( cudaMalloc(&d_xB, array_size) );
  errCheck( cudaMalloc(&d_zB, array_size) );

  //Copy into the gpu memory
  cudaMemcpyToSymbol(d_invNcst, (void *) h_invN, sizeof(biguint_t));
  cudaMemcpyToSymbol(d_Ncst, (void *) h_N, sizeof(biguint_t));

  errCheck( cudaMemcpy((void *) d_xA, (void *) h_xarray, array_size, 
                                                    cudaMemcpyHostToDevice) );
  errCheck( cudaMemcpy((void *) d_zA, (void *) h_zarray, array_size, 
                                                    cudaMemcpyHostToDevice) );
  errCheck( cudaMemcpy((void *) d_xB, (void *) h_x2array, array_size, 
                                                    cudaMemcpyHostToDevice) );
  errCheck( cudaMemcpy((void *) d_zB, (void *) h_z2array, array_size, 
                                                    cudaMemcpyHostToDevice) );

  size_t j;

  for (j=mpz_sizeinbase(s,2)-1; j>0; j-- )
  {
    if (mpz_tstbit(s, j-1) ==1 )
      Cuda_Ell_DblAdd<<<dimGrid,dimBlock>>>(d_xB, d_zB, d_xA, d_zA, firstinvd);
    else
      Cuda_Ell_DblAdd<<<dimGrid,dimBlock>>>(d_xA, d_zA, d_xB, d_zB, firstinvd);

    //cudaThreadSynchronize();
    
    //maybe only for debug mode??
    err = cudaGetLastError(); 
    if (err != cudaSuccess )
    {
      fprintf(stderr, "%s(%i) : Error cuda : %s.\n",
              __FILE__, __LINE__, cudaGetErrorString( err) );
      exit(EXIT_FAILURE);
    }
    
  }

  errCheck( cudaMemcpy((void *) h_xarray, (void *) d_xA, array_size, 
                                                   cudaMemcpyDeviceToHost) );
  errCheck( cudaMemcpy((void *) h_zarray, (void *) d_zA, array_size, 
                                                   cudaMemcpyDeviceToHost) );


  cudaFree((void *) d_xA);
  cudaFree((void *) d_zA);
  cudaFree((void *) d_xB);
  cudaFree((void *) d_zB);

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

//Compute Rmod <- A + B [Ncst] 
__device__ void Cuda_Add_mod
(biguint_t Rmod, dbigint_t cy, const biguint_t A,const biguint_t B)
{
  Cuda_Add(Rmod, cy, A, B);
  Cuda_Fully_Normalize(Rmod, cy); 
  
  if (Cuda_Cmp (Rmod, d_Ncst) >= 0)// (Rmod >= N)? 
  //if (Cuda_Cmp2 (Rmod, cy, d_Ncst) >= 0)// (Rmod >= N)? 
  {
    Cuda_Sub (Rmod, cy, Rmod, d_Ncst); 
    Cuda_Fully_Normalize(Rmod, cy); 
  }
}

//Compute Rmod  <-Rmod + A [mod] 
__device__ void Cuda_Add_mod
(biguint_t Rmod, dbigint_t cy, const biguint_t A)
{
  Cuda_Add(Rmod, cy, Rmod, A);
  Cuda_Fully_Normalize(Rmod, cy);

  if (Cuda_Cmp (Rmod, d_Ncst) >= 0)// (Rmod >= N)? 
  {
    Cuda_Sub (Rmod, cy, Rmod, d_Ncst);  
    Cuda_Fully_Normalize(Rmod, cy); 
  }
}

//Compute Rmod <- Rmod - A [mod]
__device__ void Cuda_Sub_mod 
(biguint_t Rmod, dbigint_t cy, const biguint_t A)
{
  if (Cuda_Cmp(Rmod, A)>=0)
  {
    Cuda_Sub(Rmod, cy, Rmod, A);
    Cuda_Fully_Normalize(Rmod, cy); 
  }
  else
  {
    Cuda_Add (Rmod, cy, Rmod, d_Ncst);
    Cuda_Subc (Rmod, cy, Rmod, A);
    Cuda_Fully_Normalize(Rmod, cy); 
  }
}

__device__ void Cuda_Mulmod_step
(dbiguint_t r,dbigint_t cy, unsigned int a, unsigned int b)
{
  digit_t h,l;
  int tmp;
  __mul(h,l,a,b);
  __add_cc(r[threadIdx.x],r[threadIdx.x],l);
  __addc_cc(r[threadIdx.x+1],r[threadIdx.x+1],h);
  __addcy2(cy[threadIdx.x+1]);


/*
   h=r[threadIdx.x];
   __mad_lo(r[threadIdx.x],a,b,r[threadIdx.x]);
   cy[threadIdx.x]+=(r[threadIdx.x]<h)?1:0;
   h=r[threadIdx.x+1];
   __mad_hi(r[threadIdx.x+1],a,b,r[threadIdx.x+1]);
   cy[threadIdx.x+1]+=(r[threadIdx.x+1]<h)?1:0;
*/

  __mul(h, l, d_invNcst[0]*r[0], d_Ncst[threadIdx.x]);
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

__device__ void Cuda_Dbl_mod
(biguint_t r, dbigint_t cy, biguint_t a)
{
  //cy[threadIdx.x]=(a[threadIdx.x]>>31);
  //r[threadIdx.x]=a[threadIdx.x]<<1;
  __add_cc(r[threadIdx.x],a[threadIdx.x],a[threadIdx.x]);
  __addcy2(r[(threadIdx.x+1)%NB_DIGITS]);

  //Cuda_Normalize_test(r,cy);  

  if (Cuda_Cmp (r,d_Ncst) >= 0) 
  {
    Cuda_Sub (r, cy, r, d_Ncst); 
    Cuda_Fully_Normalize(r, cy);  
  }
}


__device__ void Cuda_Mulint_mod(dbiguint_t r,dbigint_t cy, biguint_t A, unsigned int b)
{
  digit_t h,l;
  __mul(h,r[threadIdx.x],A[threadIdx.x],b);
  __add_cc(r[threadIdx.x+1],r[threadIdx.x+1],h);
  __addcy(cy[threadIdx.x+1]);

  //h*2^32+l =A[i]*B[threadIDx.x]
  __mul(h, l, d_invNcst[0]*r[0], d_Ncst[threadIdx.x]);
  __add_cc(r[threadIdx.x],r[threadIdx.x],l);
  __addc_cc(r[threadIdx.x+1],r[threadIdx.x+1],h);
  __addcy2(cy[threadIdx.x+1]);
  /*  
  r[threadIdx.x]+=l;
  cy[threadIdx.x]+=(r[threadIdx.x]<l);
 
  r[threadIdx.x+1]+=h;
  cy[threadIdx.x+1]+=(r[1+threadIdx.x]<h);
  */
  /*
  tmp=cy[threadIdx.x];
  //cy[threadIdx.x]=0;
  r[threadIdx.x+1]+=tmp;
  //on peut mettre = au lieu de += car r<2*N<2*2^1023 = 2^1024
  cy[threadIdx.x+1]=(r[threadIdx.x+1]<tmp);
  */

  //make one round of normalize + a right shift
  /*
  asm(
    "add.cc.u32 %0, %1, %2;\n\t"
    "addc.u32 %1, 0, 0;\n\t"
    : "=r"(r[threadIdx.x]), "+r"(cy[threadIdx.x]) : "r"(r[threadIdx.x+1]));
  */
  __add_cc(r[threadIdx.x],r[threadIdx.x+1],cy[threadIdx.x]);
  __addcy(cy[threadIdx.x]);
  //r[threadIdx.x]=r[threadIdx.x+1];
  //cy[threadIdx.x]=cy[threadIdx.x+1];
  if (threadIdx.x==0)
  {
    //cy[threadIdx.x+NB_DIGITS]=0;
    //r[threadIdx.x+NB_DIGITS]=0;
    cy[NB_DIGITS]=0;
    r[NB_DIGITS]=0;
  }
  
  Cuda_Fully_Normalize(r,cy); 

  if (Cuda_Cmp (r,d_Ncst) >= 0) 
  {
    Cuda_Sub (r, cy, r, d_Ncst); 
    Cuda_Fully_Normalize(r,cy); 
  }
}

__device__ void 
Cuda_Mul_mod (biguint_t mul, dbigint_t cy, const biguint_t A, const biguint_t B, 
                                                                   dbiguint_t r)
{

  int i;
  digit_t temp=A[threadIdx.x];

  r[threadIdx.x]=0;
  //r[threadIdx.x+NB_DIGITS]=0;
  
  for (i=0;i<NB_DIGITS;i++)
    Cuda_Mulmod_step(r,cy,temp,B[i]);

  Cuda_Fully_Normalize(r,cy);

  if (Cuda_Cmp (r,d_Ncst) >= 0) // mul >= N 
  {
    Cuda_Sub (mul, cy, r, d_Ncst); 
    Cuda_Fully_Normalize(mul,cy); 
  }
  else
    mul[threadIdx.x]=r[threadIdx.x];
}

__device__ void 
Cuda_Square_mod (biguint_t mul, dbigint_t cy, const biguint_t A, dbiguint_t r)
{
  Cuda_Mul_mod(mul,cy,A,A,r);
}

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

  //C=x2+z2
  Cuda_Add_mod(t,cy,v,w);
  //D=x2-z2
  Cuda_Sub_mod(v,cy,w);
  //A=x+z
  Cuda_Add_mod(w,cy,u,temp_r);
  //B=x-z
  Cuda_Sub_mod(u,cy,temp_r);

  //CB=C*B=(xq+zq)(xp-zp)
  Cuda_Mul_mod(t,cy,t,u,temp_r);
  
  //DA=D*A=(xq-zq)(xp+zp)
  Cuda_Mul_mod(v,cy,v,w,temp_r);

  //AA=A^2
  Cuda_Square_mod(w,cy,w,temp_r);
  //BB=B^2
  Cuda_Square_mod(u,cy,u,temp_r);

  //x2=AA*BB
  Cuda_Mul_mod(temp_r,cy,u,w,temp_r);
  xarg[idx1][threadIdx.x]=temp_r[threadIdx.x];

  //C= AA-BB
  Cuda_Sub_mod(w,cy,u);
  //d*C 
  Cuda_Mulint_mod(temp_r,cy,w,idx1+firstinvd);
  //BB+d*C
  Cuda_Add_mod(u,cy,temp_r);
  
  //z2=C*(BB+d*C)
  Cuda_Mul_mod(w,cy,w,u,temp_r);

  zarg[idx1][threadIdx.x]=w[threadIdx.x];
  
  //DA+CB mod N
  Cuda_Add_mod(w,cy,v,t);
  //DA-CB mod N
  Cuda_Sub_mod(v,cy,t);

  //(DA+CB)^2 mod N
  Cuda_Square_mod(w,cy,w,temp_r);
  //(DA-CB)^2 mod N
  Cuda_Square_mod(v,cy,v,temp_r);

  //z0=1 1*(DA+CB)^2
  //x0=2 2*(DA-CB)^2
  Cuda_Dbl_mod(u,cy,v);
  
  x2arg[idx1][threadIdx.x]=w[threadIdx.x];
  z2arg[idx1][threadIdx.x]=u[threadIdx.x];
}

