/*            Default code for GPU                  */
/* A compute capability of 2.0 at least is required */


__device__ void Cuda_Fully_Normalize (biguint_t A, bigint_t cy)
{
  carry_t cytemp;
  unsigned int thm1;

  while(__any(cy[threadIdx.x])!=0)
  {
    thm1 = (threadIdx.x - 1) % ECM_GPU_NB_DIGITS;
    cytemp = cy[thm1];

    __add_cc(A[threadIdx.x], A[threadIdx.x], cytemp);
  
    if (cytemp >= 0)
      __addcy(cy[threadIdx.x]);
    else /* if (cytemp < 0) */
      __subcy(cy[threadIdx.x]);
  }
}

/* Compute Rmod <- A + B */ 
/* Input: 0 <= A, B < 3*N */ 
/* Ouput: 0 <= Rmod < 6*N */ 
__device__ void Cuda_Add_mod
(biguint_t Rmod, bigint_t cy, const biguint_t A, const biguint_t B)
{
  unsigned int thp1 = (threadIdx.x + 1) % ECM_GPU_NB_DIGITS;
  __add_cc (Rmod[threadIdx.x], A[threadIdx.x], B[threadIdx.x]);
  __addcy2(Rmod[thp1]); 
  __addcy (cy[thp1]);
  Cuda_Fully_Normalize (Rmod, cy); 
}

/* Compute Rmod <- Rmod + B */ 
/* Input: 0 <= Rmod, B < 3*N */ 
/* (except when it follows Cuda_Mulint_mod, 0 <= Rmod < 3*N, 0 < B < 7*N ) */ 
/* Ouput: 0 <= Rmod < 6*N */ 
/* (except when it follows Cuda_Mulint_mod, 0 <= Rmod < 10*N) */ 
__device__ void Cuda_Add_mod
(biguint_t Rmod, bigint_t cy, const biguint_t A)
{
  unsigned int thp1 = (threadIdx.x + 1) % ECM_GPU_NB_DIGITS;
  __add_cc (Rmod[threadIdx.x], Rmod[threadIdx.x], A[threadIdx.x]);
  //__addcy (cy[threadIdx.x]);
  __addcy2(Rmod[thp1]); 
  __addcy (cy[thp1]);
  Cuda_Fully_Normalize (Rmod, cy);
}

/* Compute Rmod <- Rmod - B */ 
/* Input: 0 <= Rmod, B < 3*N */ 
/* Ouput: 0 <= Rmod < 6*N */ 
__device__ void Cuda_Sub_mod 
(biguint_t Rmod, bigint_t cy, const biguint_t B, const digit_t N3thdx)
{
  digit_t reg_Rmod = Rmod[threadIdx.x];
  carry_t reg_cy = 0; 
  
  __add_cc (reg_Rmod, reg_Rmod, N3thdx);
  __addcy (reg_cy);
  __sub_cc (reg_Rmod, reg_Rmod, B[threadIdx.x]);
  __subcy2 (reg_cy);

  Rmod[threadIdx.x] = reg_Rmod;
  cy[threadIdx.x] = reg_cy;
  Cuda_Fully_Normalize (Rmod, cy); 
}

/* Perform one step of REDC */ 
__device__ void Cuda_Mulmod_step
(biguint_t r, bigint_t cy, digit_t a, digit_t b, const digit_t Nthdx,
 const digit_t invN)
{
  digit_t t;
  digit_t reg_hi = 0;
  unsigned int thp1= (threadIdx.x + 1) % ECM_GPU_NB_DIGITS;
  carry_t reg_cy = cy[thp1];

  __mad_lo_cc(r[threadIdx.x],a,b);
  __madc_hi_cc(reg_hi,a,b);
  __addcy2(reg_cy);

  __mul_lo(t, invN, r[0]);
  __mad_lo_cc(r[threadIdx.x],t,Nthdx);
  __madc_hi_cc(reg_hi,t,Nthdx);
  __addcy2(reg_cy);

  /* make one round of normalize + a right shift at the same time */
  __add_cc(r[threadIdx.x],r[thp1],reg_hi);
  __addc_cc(r[thp1],r[thp1],reg_cy);
  __addcy(cy[thp1]); 
}

/* Compute r <- 2*a */ 
/* Input: 0 <= a < 3*N */ 
/* Ouput: 0 <= r < 3*N */ 
__device__ void Cuda_Dbl_mod
(biguint_t r, biguint_t a)
{
  unsigned int thp1= (threadIdx.x + 1) % ECM_GPU_NB_DIGITS;
  asm ("add.cc.u32 %0, %1, %1;" : "=r"(r[threadIdx.x]) : "r"(a[threadIdx.x]));
  __addcy2(r[thp1]);
}


/* Compute r <- A*b */ 
/* Input: 0 < b < 2^SIZE_DIGIT, 0 <= A < 6*N */ 
/* Ouput: 0 <= r < 7*N */ 
__device__ void Cuda_Mulint_mod
(biguint_t r, bigint_t cy, biguint_t A, digit_t b, const digit_t Nthdx,
 const digit_t invN)
{
  digit_t t;
  digit_t reg_hi;
  unsigned int thp1= (threadIdx.x + 1) % ECM_GPU_NB_DIGITS;
  digit_t reg_A = A[threadIdx.x];
  carry_t reg_cy;

  __mul_lo(r[threadIdx.x],reg_A,b);
  __mul_hi(reg_hi,reg_A,b);

  __mul_lo(t, invN, r[0]);
  __mad_lo_cc(r[threadIdx.x],t,Nthdx);
  __madc_hi_cc(reg_hi,t,Nthdx);
  __addcy(reg_cy);

  /* make one round of normalize + a right shift at the same time */
  __add_cc(r[threadIdx.x],r[thp1],reg_hi);
  __addc_cc(r[thp1],r[thp1],reg_cy);
  __addcy(cy[thp1]); 

  Cuda_Fully_Normalize(r,cy); 
}

/* Compute r <- A*B */ 
/* Input: 0 <= A, B < 6*N */
/* (except when it follows Cuda_Mulint_mod, 0 <= A < 6*N, 0 < B < 10*N ) */ 
/* Ouput: 0 <= r < 3*N */ 
__device__ void Cuda_Mul_mod 
(biguint_t mul, bigint_t cy, const biguint_t A, const biguint_t B, biguint_t r,
 const digit_t Nthdx, const digit_t invN)
{

  int i;
  digit_t temp=A[threadIdx.x];

  r[threadIdx.x]=0;
  
  for (i=0; i<ECM_GPU_NB_DIGITS; i++)
    Cuda_Mulmod_step (r, cy, temp, B[i], Nthdx, invN);

  
  Cuda_Fully_Normalize (r, cy);
  mul[threadIdx.x]=r[threadIdx.x];
}

__device__ void Cuda_Square_mod 
(biguint_t mul, bigint_t cy, const biguint_t A, biguint_t r, 
 const digit_t Nthdx, const digit_t invN)
{
  Cuda_Mul_mod (mul, cy, A, A, r, Nthdx, invN);
}

/* 
  Compute simultaneously:
  (xarg : zarg ) <- [2](xarg : zarg) 
  (xarg2 : zarg2 ) <- (xarg : zarg) + (xarg2 : zarg2) 
*/
__global__ void 
Cuda_Ell_DblAdd (biguint_t *xAarg, biguint_t *zAarg, biguint_t *xBarg, 
                                       biguint_t *zBarg, unsigned int firstinvd)
{
  __shared__ VOL digit_t b_temp_r[ECM_GPU_CURVES_BY_BLOCK][ECM_GPU_NB_DIGITS];
  __shared__ VOL carry_t b_cy[ECM_GPU_CURVES_BY_BLOCK][ECM_GPU_NB_DIGITS]; 

  __shared__ VOL digit_t b_t[ECM_GPU_CURVES_BY_BLOCK][ECM_GPU_NB_DIGITS];
  __shared__ VOL digit_t b_u[ECM_GPU_CURVES_BY_BLOCK][ECM_GPU_NB_DIGITS];
  __shared__ VOL digit_t b_v[ECM_GPU_CURVES_BY_BLOCK][ECM_GPU_NB_DIGITS];
  __shared__ VOL digit_t b_w[ECM_GPU_CURVES_BY_BLOCK][ECM_GPU_NB_DIGITS];
  
  VOL digit_t *t=b_t[threadIdx.y];
  VOL digit_t *u=b_u[threadIdx.y];
  VOL digit_t *v=b_v[threadIdx.y];
  VOL digit_t *w=b_w[threadIdx.y];
  VOL digit_t *temp_r=b_temp_r[threadIdx.y];
  VOL carry_t *cy=b_cy[threadIdx.y];

  /* Init of shared variables */
  const unsigned int idx1=blockIdx.x*blockDim.y+threadIdx.y;
  //unsigned int t1=threadIdx.x+1;
  cy[threadIdx.x]=0; 

  w[threadIdx.x]=xBarg[idx1][threadIdx.x];
  v[threadIdx.x]=zBarg[idx1][threadIdx.x];
  temp_r[threadIdx.x]=xAarg[idx1][threadIdx.x];
  u[threadIdx.x]=zAarg[idx1][threadIdx.x];

  const digit_t Nthdx = d_Ncst[threadIdx.x]; 
  const digit_t N3thdx = d_3Ncst[threadIdx.x]; 
  const digit_t invN = d_invNcst; 

  Cuda_Add_mod(t, cy, v, w);           /* C=zB+xB */
  Cuda_Sub_mod(v, cy, w, N3thdx);      /* D=zB-xB */
  Cuda_Add_mod(w, cy, u, temp_r);      /* A=zA+xA */
  Cuda_Sub_mod(u, cy, temp_r, N3thdx); /* B=zA-xA */

  Cuda_Mul_mod(t, cy, t, u, temp_r, Nthdx, invN); /* CB=C*B=(zB+xB)(zA-xA) */
  Cuda_Mul_mod(v, cy, v, w, temp_r, Nthdx, invN); /* DA=D*A=(zB-xB)(zA+xA) */

  Cuda_Square_mod(w, cy, w, temp_r, Nthdx, invN); /* AA=A^2 */
  Cuda_Square_mod(u, cy, u, temp_r, Nthdx, invN); /* BB=B^2 */

  Cuda_Mul_mod(temp_r, cy, u, w, temp_r, Nthdx, invN); /* AA*BB */
  xAarg[idx1][threadIdx.x]=temp_r[threadIdx.x];

  Cuda_Sub_mod (w, cy, u, N3thdx); /* K= AA-BB */
  Cuda_Mulint_mod (temp_r, cy, w, idx1 + firstinvd, Nthdx, invN); /* d*K */ 
  Cuda_Add_mod (u, cy, temp_r); /* BB+d*K */
 
  Cuda_Mul_mod (w, cy, w, u, temp_r, Nthdx, invN); /* K*(BB+d*K) */
  zAarg[idx1][threadIdx.x]=w[threadIdx.x];
 
  Cuda_Add_mod(w, cy, v, t);       /* DA+CB mod N */
  Cuda_Sub_mod(v, cy, t, N3thdx);  /* DA-CB mod N */

  Cuda_Square_mod(w, cy, w, temp_r, Nthdx, invN); /* (DA+CB)^2 mod N */
  Cuda_Square_mod(v, cy, v, temp_r, Nthdx, invN); /* (DA-CB)^2 mod N */

  /* z0=1 so there is nothing to compute for z0*(DA+CB)^2 */
  Cuda_Dbl_mod(temp_r, v); /* x0=2 x0*(DA-CB)^2 */
  
  xBarg[idx1][threadIdx.x]=w[threadIdx.x];
  zBarg[idx1][threadIdx.x]=temp_r[threadIdx.x];
}

