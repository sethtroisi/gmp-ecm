#include "cudaarith.h"

//specific functions for compute capability 1.3
//#ifdef CC13
//  (A > B)?, returns 1(true), -1(false) or 0(a=b) 
//Assume A and B are normalize (no carry or borrow)
__device__ int Cuda_Cmp(const biguint_t A, const biguint_t B)
{
	unsigned int i = SIZE_NUMBER-1;
	do
	{
		if (A[i] > B[i])
			return 1;
		if (A[i] < B[i])
			return -1;
		i--;
	}while(i!=-1);
	return 0;
}
//#endif
/*
//specific functions for compute capability 20
#ifdef CC20
//  (A > B)?, returns 1(true), -1(false) or 0(a=b) 
//Assume A and B are normalize (no carry or borrow)
__device__ int Cuda_Is_Eq(const biguint_t A, const biguint_t B)
{
	return A[threadIdx.x]==B[threadIdx.x];
}
__device__ int Cuda_Is_Gt(const biguint_t A, const biguint_t B)
{
	return A[threadIdx.x]>B[threadIdx.x];
}
__device__ int Cuda_Cmp(const biguint_t A, const biguint_t B)
{
	unsigned int ballot=__ballot(Cuda_Is_Eq(A,B));
	if (ballot==0)
		return 0;
	else
		if (__ballot(Cuda_Is_Gt(A,B))>__ballot(Cuda_Is_Gt(B,A)))
			return 1;
		else
			return -1;
}
#endif
*/

//return -1 if A<cst, 0 if A=cst, 1 if A>cst
__device__ int Cuda_ge_cst(const biguint_t A, const biguint_t cst, const int cy)
{
	unsigned int i = SIZE_NUMBER-1;
	if (cy>0)
		return 1;
	do
	{
		if (A[i] > cst[i])
			return 1;
		if (A[i] < cst[i])
			return -1;
		i--;
	}while(i!=0);
	return 0;
}

//Normalise a result; 
__device__ int Cuda_Is_Normalize(dbigint_t cy)
{
	if (threadIdx.x==SIZE_NUMBER-1 || cy[threadIdx.x]==0)
		return 0;
	else
		return 1;
}

__device__ void Cuda_Normalize(biguint_t A,dbigint_t cy)
{
	//if (threadIdx.x==0)
	//	cptnormalize[blockIdx.x]++;

	int oldcy;
	if (threadIdx.x==0)
		oldcy = 0;
	else
	{
		oldcy = cy[threadIdx.x-1];
		cy[threadIdx.x-1]=0;
	}

	if (oldcy>=0)
	{
		//Cuda_Add1(A,cy,oldcy);
		A[threadIdx.x]+=oldcy;
		cy[threadIdx.x]+=(A[threadIdx.x]<oldcy);
	}
	else // oldcy < 0
	{
		//Cuda_Sub1(A,cy,-oldcy);
		cy[threadIdx.x]-=(A[threadIdx.x]< -oldcy);
		A[threadIdx.x]-= -oldcy;
	}
}

__device__ void Cuda_Fully_Normalize(biguint_t A,dbigint_t cy)
{
	do
	{
	Cuda_Normalize(A,cy);
	}while(__any(Cuda_Is_Normalize(cy))!=0);
}


__device__ void Cuda_Add (biguint_t r, dbigint_t cy ,const biguint_t b)
{
	r[threadIdx.x]+=b[threadIdx.x];
	cy[threadIdx.x]+=(r[threadIdx.x]<b[threadIdx.x]);
}

//Assume r and b are different; r and a can be the same.
__device__ void Cuda_Add 
(biguint_t r, dbigint_t cy ,const biguint_t a, const biguint_t b)
{
	r[threadIdx.x]=a[threadIdx.x];
	r[threadIdx.x]+=b[threadIdx.x];
	cy[threadIdx.x]+=(r[threadIdx.x]<b[threadIdx.x]);
}

__device__ void Cuda_Sub (biguint_t r, dbigint_t cy, const biguint_t b)
{
	cy[threadIdx.x]-=(r[threadIdx.x] < b[threadIdx.x]);
	r[threadIdx.x] -= b[threadIdx.x];
}

//Assume r and b are different; r and a can be the same.
__device__ void Cuda_Sub 
(biguint_t r, dbigint_t cy, const biguint_t a, const biguint_t b)
{
	r[threadIdx.x]=a[threadIdx.x];
	cy[threadIdx.x]-=(r[threadIdx.x] < b[threadIdx.x]);
	r[threadIdx.x] -= b[threadIdx.x];
}

//Compute Rmod <- A + B [mod] 
__device__ void Cuda_Add_mod
(biguint_t Rmod, dbigint_t cy, const biguint_t A,const biguint_t B,const biguint_t mod)
{
  Cuda_Add(Rmod,cy,A,B);
	Cuda_Fully_Normalize(Rmod,cy);	
  
  if (Cuda_ge_cst (Rmod,mod,cy[SIZE_NUMBER-1]) >= 0) // (a >= mod)? 
	{
   	Cuda_Sub (Rmod, cy, mod); // R <- R - mod 
		Cuda_Fully_Normalize(Rmod,cy);	
	}
}

//Compute Rmod  <-Rmod + A [mod] 
__device__ void Cuda_Add_mod
(biguint_t Rmod, dbigint_t cy, const biguint_t A,const biguint_t mod)
{
 	Cuda_Add(Rmod,cy,A);
	Cuda_Fully_Normalize(Rmod,cy);	

  if (Cuda_ge_cst (Rmod,mod,cy[SIZE_NUMBER-1]) >= 0) // (a >= mod)? 
	{
   	Cuda_Sub (Rmod, cy, mod);  
		Cuda_Fully_Normalize(Rmod,cy);	
	}
}


//Compute Rmod <- A - B [mod] 
__device__ void Cuda_Sub_mod 
(biguint_t Rmod, dbigint_t cy, const biguint_t A,const biguint_t B,const biguint_t mod)
{
 	Cuda_Sub(Rmod,cy,A,B);
	Cuda_Fully_Normalize(Rmod,cy);	
	
 	if (cy[SIZE_NUMBER-1] <0 ) // we should subtract 1 at smod[n] 
 	{
   	Cuda_Add (Rmod, cy, mod);
		Cuda_Fully_Normalize(Rmod,cy);	
 	}
}

//Compute Rmod <- Rmod - A [mod]
__device__ void Cuda_Sub_mod 
(biguint_t Rmod, dbigint_t cy, const biguint_t A, const biguint_t mod)
{
 	Cuda_Sub(Rmod,cy,A);
	Cuda_Fully_Normalize(Rmod,cy);	

 	if (cy[SIZE_NUMBER-1] <0 ) 
 	{
   	Cuda_Add (Rmod, cy, mod);
		Cuda_Fully_Normalize(Rmod,cy);	
 	}
}

//  Return h, l such that h*2^32 + l = A*B 
__device__ void Cuda_Mul_uint (unsigned int *h, unsigned int *l, const unsigned int A,const unsigned int B)
{
		*h=__umulhi(A,B);
		*l=A*B;
}

__device__ void Cuda_Mul
(dbiguint_t R, dbigint_t cy, const biguint_t A,const biguint_t B)
{
	int i;
	unsigned int h,l;
	
	for (i=0;i<SIZE_NUMBER;i++)
	{
		//h*2^32+l =A[i]*B[threadIDx.x]
		Cuda_Mul_uint(&h,&l,A[threadIdx.x],B[i]);
		
		R[i+threadIdx.x] +=l;
		cy[i+threadIdx.x]+=(R[i+threadIdx.x] < l);
		
		R[i+1+threadIdx.x] +=h;
		cy[i+1+threadIdx.x]+=(R[i+1+threadIdx.x]<h);
	}

	//Normalize : but R and cy are 2 * SIZE_NUMBER long
	Cuda_Fully_Normalize(R,cy);
	Cuda_Fully_Normalize(R+SIZE_NUMBER-1,cy+SIZE_NUMBER-1);

	//if (threadIdx.x==0)
	//	cptnormalizemul[blockIdx.x]++;

	if (threadIdx.x==0)
	{
		R[2*SIZE_NUMBER-1]+=cy[2*SIZE_NUMBER-2];
		cy[2*SIZE_NUMBER-2]=0;//only to let cy clean in order to re-use it
	}
}

__device__ void Cuda_Square (dbiguint_t R, dbigint_t cy, const biguint_t A)
{
	Cuda_Mul(R,cy,A,A);
}

//Assume r<N^2
__device__ void Cuda_RedMontgomery (biguint_t mul, dbigint_t cy, const biguint_t mod, dbiguint_t r, dbiguint_t temp)
{
	//temp=((r mod 2^(32*SIZE_NUMBER))*mod^-1) (mod^-1 already compute)
	Cuda_Mul(temp,cy,r,invmodcst);//pour r que la partie mod R compte : Ok
	
	//mul = temp (mod 2^(32*SIZE_NUMBER)) 
	mul[threadIdx.x]=temp[threadIdx.x];
	temp[threadIdx.x]=0;
	temp[threadIdx.x+SIZE_NUMBER]=0;
	
	//temp=mul*m
	Cuda_Mul(temp,cy,mul,mod);
	
	//r=r+temp // r and temp2 are 2*SIZE_NUMBER long
	Cuda_Add(r,cy,temp);
	Cuda_Add(r+SIZE_NUMBER,cy+SIZE_NUMBER,temp+SIZE_NUMBER);
	//Normalize : but R and cy are 2 * SIZE_NUMBER long
	Cuda_Fully_Normalize(r,cy);
	Cuda_Fully_Normalize(r+SIZE_NUMBER-1,cy+SIZE_NUMBER-1);

	if (threadIdx.x==0)
		{
		r[2*SIZE_NUMBER-1]+=cy[2*SIZE_NUMBER-2];
		cy[2*SIZE_NUMBER-2]=0;//only to let cy clean in order to re-use it
		}

	//return r/ 2^(32*SIZE_NUMBER)
	mul[threadIdx.x]=r[threadIdx.x+SIZE_NUMBER];

	if (Cuda_Cmp (mul,mod) >= 0) // (mul >= mod)? 
	{
  	Cuda_Sub (mul, cy, mod); 
		Cuda_Fully_Normalize(mul,cy);	
	}
}

//Assume A ans B are the montgomery representation
//Compute mul = A * B * 2^-(32*SIZE_NUMBER) mod[mod]
// r and temp have size 2*SIZE_NUMBER
__device__ void Cuda_Mul_mod (biguint_t mul, dbigint_t cy, const biguint_t A,const biguint_t B, const biguint_t mod, dbiguint_t r, dbiguint_t temp)
{
	temp[threadIdx.x]=0;
	temp[threadIdx.x+SIZE_NUMBER]=0;
	r[threadIdx.x]=0;
	r[threadIdx.x+SIZE_NUMBER]=0;
	
	//__syncthreads();
	//r=a*b
	Cuda_Mul(r,cy,A,B);
	Cuda_RedMontgomery (mul, cy, mod, r, temp);
}


//Assume A ans B are the montgomery representation
//Compute mul = A * A * 2^-(32*SIZE_NUMBER) mod[mod]
// r and temp have size 2*SIZE_NUMBER
__device__ void Cuda_Square_mod (biguint_t mul, dbigint_t cy, const biguint_t A, const biguint_t mod, dbiguint_t r, dbiguint_t temp)
{
	temp[threadIdx.x]=0;
	temp[threadIdx.x+SIZE_NUMBER]=0;
	r[threadIdx.x]=0;
	r[threadIdx.x+SIZE_NUMBER]=0;
	
	//__syncthreads();
	//r=a*b
	Cuda_Square(r,cy,A);
	Cuda_RedMontgomery (mul, cy, mod, r, temp);
}

__device__ void Cuda_Ell_Dbl(biguint_t x2p, biguint_t z2p, const biguint_t xp, const biguint_t zp, const biguint_t N, const biguint_t d, biguint_t temp_u, biguint_t temp_v, dbiguint_t temp_r, dbiguint_t temp_r2, dbigint_t cy)
{
	//u<-xp+zp mod N
	Cuda_Add_mod(temp_u,cy,xp,zp,N);
	//u <- u^2
	Cuda_Square_mod(temp_u,cy,temp_u,N,temp_r,temp_r2);
	//v<-xp-zp mod N
	Cuda_Sub_mod(temp_v,cy,xp,zp,N);
	//v <- v^2
	Cuda_Square_mod(temp_v,cy,temp_v,N,temp_r,temp_r2);
	//x2p=u-v mod N    x2p is used as a temporary variable here
	Cuda_Sub_mod(x2p,cy,temp_u,temp_v,N);
	//z2p<-x2p*d mod N z2p is used as a temporary variable here
	Cuda_Mul_mod(z2p,cy,x2p,d,N,temp_r,temp_r2);
	//z2p<- z2p+v mod N
	Cuda_Add_mod(z2p,cy,temp_v,N);
	//z2p <- x2p*z2p
	Cuda_Mul_mod(z2p,cy,x2p,z2p,N,temp_r,temp_r2);
	//x2p <- u*v
	Cuda_Mul_mod(x2p,cy,temp_u,temp_v,N,temp_r,temp_r2);
}

__device__ void Cuda_Ell_Add(biguint_t xplus, biguint_t zplus, const biguint_t xp, const biguint_t zp, const biguint_t xq, const biguint_t zq, const biguint_t xminus, const biguint_t zminus, const biguint_t N, biguint_t temp_u, biguint_t temp_v, biguint_t temp_w, dbiguint_t temp_r, dbiguint_t temp_r2, dbigint_t cy)
{
	unsigned int tmp;
	
	//u<-xp+zp mod N
	Cuda_Add_mod(temp_u,cy,xp,zp,N);
	//v<-xq-zq mod N
	Cuda_Sub_mod(temp_v,cy,xq,zq,N);
	//v<-u*v mod N
	Cuda_Mul_mod(temp_v,cy,temp_u,temp_v,N,temp_r,temp_r2);
	//u<-xp-zp mod N 
	Cuda_Sub_mod(temp_u,cy,xp,zp,N);
	//w<-xq+zq mod N
	Cuda_Add_mod(temp_w,cy,zq,xq,N);
	//u<-u*w mod N
	Cuda_Mul_mod(temp_u,cy,temp_u,temp_w,N,temp_r,temp_r2);

	//w<-v+u mod N
	Cuda_Add_mod(temp_w,cy,temp_v,temp_u,N);
	//v<-v-u mod N
	Cuda_Sub_mod(temp_v,cy,temp_u,N);
	//w<-w^2 mod N
	Cuda_Square_mod(temp_w,cy,temp_w,N,temp_r,temp_r2);
	//v<-v^2 mod N
	Cuda_Square_mod(temp_v,cy,temp_v,N,temp_r,temp_r2);
	
	if (xplus==xminus) //same variable : in-place variant
	{
		Cuda_Mul_mod(zplus,cy,zminus,temp_w,N,temp_r,temp_r2);
		Cuda_Mul_mod(xplus,cy,xminus,temp_v,N,temp_r,temp_r2);
		//swap
		tmp=xplus[threadIdx.x];
		xplus[threadIdx.x]=zplus[threadIdx.x];
		zplus[threadIdx.x]=tmp;
	}
	else
	{
		//xplus <- zminus*u mod N
		Cuda_Mul_mod(xplus,cy,zminus,temp_w,N,temp_r,temp_r2);
		//zplus <- xminus*v mod N
		Cuda_Mul_mod(zplus,cy,xminus,temp_v,N,temp_r,temp_r2);
	}
}


#ifdef CC13
//Compute [3^power3*2^power2]A on the elliptic curve
//prepare for prac, ie set B and C equal to A and set A to [2]A
__global__ void Cuda_Ell_Mul_2_3(biguint_t *xAarg, biguint_t *zAarg, biguint_t *xBarg, biguint_t *zBarg, biguint_t *xCarg, biguint_t *zCarg, const biguint_t *darg, const unsigned int power2, const unsigned int power3)
{
	__shared__ unsigned int temp_r[2*SIZE_NUMBER];
	__shared__ unsigned int temp_r2[2*SIZE_NUMBER];
	__shared__ unsigned int temp_u[SIZE_NUMBER];
	__shared__ unsigned int temp_v[SIZE_NUMBER];
	__shared__ unsigned int temp_w[SIZE_NUMBER];

	__shared__ int cy[2*SIZE_NUMBER]; 
	
	__shared__ unsigned int xA[SIZE_NUMBER];
	__shared__ unsigned int zA[SIZE_NUMBER];
	__shared__ unsigned int xB[SIZE_NUMBER];
	__shared__ unsigned int zB[SIZE_NUMBER];
	
	__shared__ unsigned int d[SIZE_NUMBER];
	
	unsigned int i;

	//init
	xA[threadIdx.x]=xAarg[blockIdx.x][threadIdx.x];
	zA[threadIdx.x]=zAarg[blockIdx.x][threadIdx.x];
	xB[threadIdx.x]=xAarg[blockIdx.x][threadIdx.x];
	zB[threadIdx.x]=zAarg[blockIdx.x][threadIdx.x];

	d[threadIdx.x]=darg[blockIdx.x][threadIdx.x];

	temp_r[threadIdx.x]=0;	
	temp_r[SIZE_NUMBER + threadIdx.x]=0;	
	temp_r2[threadIdx.x]=0;	
	temp_r2[SIZE_NUMBER + threadIdx.x]=0;	
	temp_u[threadIdx.x]=0;	
	temp_v[threadIdx.x]=0;	
	temp_w[threadIdx.x]=0;	
	
	cy[threadIdx.x]=0;	
	cy[SIZE_NUMBER + threadIdx.x]=0;	

 	//__syncthreads();
	//Now one can begin the computation

	for(i=0;i<power2;i++)
	{
		Cuda_Ell_Dbl(xB, zB, xB, zB, Ncst, d, temp_u, temp_v, temp_r, temp_r2, cy);
	}
	
	for(i=0;i<power3;i++)
	{
		Cuda_Ell_Dbl(xA, zA, xB, zB, Ncst, d, temp_u, temp_v, temp_r, temp_r2, cy);
	Cuda_Ell_Add(xB, zB, xA, zA, xB, zB, xB, zB, Ncst, temp_u, temp_v, temp_w, temp_r, temp_r2, cy);
	}

	Cuda_Ell_Dbl(xA, zA, xB, zB, Ncst, d, temp_u, temp_v, temp_r, temp_r2, cy);

	//End of the computation; Copy the results for the cpu
	xAarg[blockIdx.x][threadIdx.x]=xA[threadIdx.x];
	zAarg[blockIdx.x][threadIdx.x]=zA[threadIdx.x];
	xBarg[blockIdx.x][threadIdx.x]=xB[threadIdx.x];
	zBarg[blockIdx.x][threadIdx.x]=zB[threadIdx.x];
	xCarg[blockIdx.x][threadIdx.x]=xB[threadIdx.x];
	zCarg[blockIdx.x][threadIdx.x]=zB[threadIdx.x];
}

__global__ void Cuda_Ell_Add_Dbl(biguint_t *xAarg, biguint_t *zAarg, biguint_t *xBarg, biguint_t *zBarg, biguint_t *xCarg, biguint_t *zCarg, const biguint_t *darg)
{
	__shared__ unsigned int temp_r[2*SIZE_NUMBER];
	__shared__ unsigned int temp_r2[2*SIZE_NUMBER];
	__shared__ unsigned int temp_u[SIZE_NUMBER];
	__shared__ unsigned int temp_v[SIZE_NUMBER];
	__shared__ unsigned int temp_w[SIZE_NUMBER];

	__shared__ int cy[2*SIZE_NUMBER]; 
	
	__shared__ unsigned int xA[SIZE_NUMBER];
	__shared__ unsigned int zA[SIZE_NUMBER];
	__shared__ unsigned int xB[SIZE_NUMBER];
	__shared__ unsigned int zB[SIZE_NUMBER];
	__shared__ unsigned int xC[SIZE_NUMBER];
	__shared__ unsigned int zC[SIZE_NUMBER];
	
	__shared__ unsigned int d[SIZE_NUMBER];
	
	//init
	xA[threadIdx.x]=xAarg[blockIdx.x][threadIdx.x];
	zA[threadIdx.x]=zAarg[blockIdx.x][threadIdx.x];
	xB[threadIdx.x]=xBarg[blockIdx.x][threadIdx.x];
	zB[threadIdx.x]=zBarg[blockIdx.x][threadIdx.x];
	xC[threadIdx.x]=xCarg[blockIdx.x][threadIdx.x];
	zC[threadIdx.x]=zCarg[blockIdx.x][threadIdx.x];

	d[threadIdx.x]=darg[blockIdx.x][threadIdx.x];

	temp_r[threadIdx.x]=0;	
	temp_r[SIZE_NUMBER + threadIdx.x]=0;	
	temp_r2[threadIdx.x]=0;	
	temp_r2[SIZE_NUMBER + threadIdx.x]=0;	
	temp_u[threadIdx.x]=0;	
	temp_v[threadIdx.x]=0;	
	temp_w[threadIdx.x]=0;	
	
	cy[threadIdx.x]=0;	
	cy[SIZE_NUMBER + threadIdx.x]=0;	

 	//__syncthreads();
	//Now one can begin the computation

	//Add	
	Cuda_Ell_Add(xB, zB, xA, zA, xB, zB, xC, zC, Ncst, temp_u, temp_v, temp_w, temp_r, temp_r2, cy);
	
	//Dbl
	Cuda_Ell_Dbl(xA, zA, xA, zA, Ncst, d, temp_u, temp_v, temp_r, temp_r2, cy);


	//End of the computation; Copy the results for the cpu
	xAarg[blockIdx.x][threadIdx.x]=xA[threadIdx.x];
	zAarg[blockIdx.x][threadIdx.x]=zA[threadIdx.x];
	xBarg[blockIdx.x][threadIdx.x]=xB[threadIdx.x];
	zBarg[blockIdx.x][threadIdx.x]=zB[threadIdx.x];
	xCarg[blockIdx.x][threadIdx.x]=xC[threadIdx.x];
	zCarg[blockIdx.x][threadIdx.x]=zC[threadIdx.x];
}

__global__ void Cuda_Ell_3Add(biguint_t *xAarg, biguint_t *zAarg, biguint_t *xBarg, biguint_t *zBarg, biguint_t *xCarg, biguint_t *zCarg)
{
	__shared__ unsigned int temp_r[2*SIZE_NUMBER];
	__shared__ unsigned int temp_r2[2*SIZE_NUMBER];
	__shared__ unsigned int temp_u[SIZE_NUMBER];
	__shared__ unsigned int temp_v[SIZE_NUMBER];
	__shared__ unsigned int temp_w[SIZE_NUMBER];

	__shared__ int cy[2*SIZE_NUMBER]; 
	
	__shared__ unsigned int xA[SIZE_NUMBER];
	__shared__ unsigned int zA[SIZE_NUMBER];
	__shared__ unsigned int xB[SIZE_NUMBER];
	__shared__ unsigned int zB[SIZE_NUMBER];
	__shared__ unsigned int xC[SIZE_NUMBER];
	__shared__ unsigned int zC[SIZE_NUMBER];
	__shared__ unsigned int xT[SIZE_NUMBER];
	__shared__ unsigned int zT[SIZE_NUMBER];
	__shared__ unsigned int xT2[SIZE_NUMBER];
	__shared__ unsigned int zT2[SIZE_NUMBER];
	
	//init
	xA[threadIdx.x]=xAarg[blockIdx.x][threadIdx.x];
	zA[threadIdx.x]=zAarg[blockIdx.x][threadIdx.x];
	xB[threadIdx.x]=xBarg[blockIdx.x][threadIdx.x];
	zB[threadIdx.x]=zBarg[blockIdx.x][threadIdx.x];
	xC[threadIdx.x]=xCarg[blockIdx.x][threadIdx.x];
	zC[threadIdx.x]=zCarg[blockIdx.x][threadIdx.x];

	xT[threadIdx.x]=0;
	zT[threadIdx.x]=0;
	xT2[threadIdx.x]=0;
	zT2[threadIdx.x]=0;
	
	temp_r[threadIdx.x]=0;	
	temp_r[SIZE_NUMBER + threadIdx.x]=0;	
	temp_r2[threadIdx.x]=0;	
	temp_r2[SIZE_NUMBER + threadIdx.x]=0;	
	temp_u[threadIdx.x]=0;	
	temp_v[threadIdx.x]=0;	
	temp_w[threadIdx.x]=0;	
	
	cy[threadIdx.x]=0;	
	cy[SIZE_NUMBER + threadIdx.x]=0;	

 	//__syncthreads();
	//Now one can begin the computation

	//Add	
	Cuda_Ell_Add(xT, zT, xA, zA, xB, zB, xC, zC, Ncst, temp_u, temp_v, temp_w, temp_r, temp_r2, cy);
	Cuda_Ell_Add(xT2, zT2, xT, zT, xA, zA, xB, zB, Ncst, temp_u, temp_v, temp_w, temp_r, temp_r2, cy);
	Cuda_Ell_Add(xB, zB, xB, zB, xT, zT, xA, zA, Ncst, temp_u, temp_v, temp_w, temp_r, temp_r2, cy);
	
	//End of the computation; Copy the results for the cpu
	xAarg[blockIdx.x][threadIdx.x]=xT2[threadIdx.x];
	zAarg[blockIdx.x][threadIdx.x]=zT2[threadIdx.x];
	xBarg[blockIdx.x][threadIdx.x]=xB[threadIdx.x];
	zBarg[blockIdx.x][threadIdx.x]=zB[threadIdx.x];
}

__global__ void Cuda_Ell_Add_Perm(biguint_t *xAarg, biguint_t *zAarg, biguint_t *xBarg, biguint_t *zBarg, biguint_t *xCarg, biguint_t *zCarg)
{
	__shared__ unsigned int temp_r[2*SIZE_NUMBER];
	__shared__ unsigned int temp_r2[2*SIZE_NUMBER];
	__shared__ unsigned int temp_u[SIZE_NUMBER];
	__shared__ unsigned int temp_v[SIZE_NUMBER];
	__shared__ unsigned int temp_w[SIZE_NUMBER];

	__shared__ int cy[2*SIZE_NUMBER]; 
	
	__shared__ unsigned int xA[SIZE_NUMBER];
	__shared__ unsigned int zA[SIZE_NUMBER];
	__shared__ unsigned int xB[SIZE_NUMBER];
	__shared__ unsigned int zB[SIZE_NUMBER];
	__shared__ unsigned int xC[SIZE_NUMBER];
	__shared__ unsigned int zC[SIZE_NUMBER];
	
	//init
	xA[threadIdx.x]=xAarg[blockIdx.x][threadIdx.x];
	zA[threadIdx.x]=zAarg[blockIdx.x][threadIdx.x];
	xB[threadIdx.x]=xBarg[blockIdx.x][threadIdx.x];
	zB[threadIdx.x]=zBarg[blockIdx.x][threadIdx.x];
	xC[threadIdx.x]=xCarg[blockIdx.x][threadIdx.x];
	zC[threadIdx.x]=zCarg[blockIdx.x][threadIdx.x];

	temp_r[threadIdx.x]=0;	
	temp_r[SIZE_NUMBER + threadIdx.x]=0;	
	temp_r2[threadIdx.x]=0;	
	temp_r2[SIZE_NUMBER + threadIdx.x]=0;	
	temp_u[threadIdx.x]=0;	
	temp_v[threadIdx.x]=0;	
	temp_w[threadIdx.x]=0;	
	
	cy[threadIdx.x]=0;	
	cy[SIZE_NUMBER + threadIdx.x]=0;	

 	//__syncthreads();
	//Now one can begin the computation

	//Add	
	Cuda_Ell_Add(xA, zA, xB, zB, xA, zA, xC, zC, Ncst, temp_u, temp_v, temp_w, temp_r, temp_r2, cy);
	
	//End of the computation; Copy the results for the cpu
	xBarg[blockIdx.x][threadIdx.x]=xA[threadIdx.x];
	zBarg[blockIdx.x][threadIdx.x]=zA[threadIdx.x];
	xCarg[blockIdx.x][threadIdx.x]=xB[threadIdx.x];
	zCarg[blockIdx.x][threadIdx.x]=zB[threadIdx.x];
}

__global__ void Cuda_Ell_Dbl_3Add(biguint_t *xAarg, biguint_t *zAarg, biguint_t *xBarg, biguint_t *zBarg, biguint_t *xCarg, biguint_t *zCarg, const biguint_t *darg)
{
	__shared__ unsigned int temp_r[2*SIZE_NUMBER];
	__shared__ unsigned int temp_r2[2*SIZE_NUMBER];
	__shared__ unsigned int temp_u[SIZE_NUMBER];
	__shared__ unsigned int temp_v[SIZE_NUMBER];
	__shared__ unsigned int temp_w[SIZE_NUMBER];

	__shared__ int cy[2*SIZE_NUMBER]; 
	
	__shared__ unsigned int xA[SIZE_NUMBER];
	__shared__ unsigned int zA[SIZE_NUMBER];
	__shared__ unsigned int xB[SIZE_NUMBER];
	__shared__ unsigned int zB[SIZE_NUMBER];
	__shared__ unsigned int xC[SIZE_NUMBER];
	__shared__ unsigned int zC[SIZE_NUMBER];
	__shared__ unsigned int xT[SIZE_NUMBER];
	__shared__ unsigned int zT[SIZE_NUMBER];
	__shared__ unsigned int xT2[SIZE_NUMBER];
	__shared__ unsigned int zT2[SIZE_NUMBER];

	__shared__ unsigned int d[SIZE_NUMBER];
	
	//init
	xA[threadIdx.x]=xAarg[blockIdx.x][threadIdx.x];
	zA[threadIdx.x]=zAarg[blockIdx.x][threadIdx.x];
	xB[threadIdx.x]=xBarg[blockIdx.x][threadIdx.x];
	zB[threadIdx.x]=zBarg[blockIdx.x][threadIdx.x];
	xC[threadIdx.x]=xCarg[blockIdx.x][threadIdx.x];
	zC[threadIdx.x]=zCarg[blockIdx.x][threadIdx.x];

	xT[threadIdx.x]=0;
	zT[threadIdx.x]=0;
	xT2[threadIdx.x]=0;
	zT2[threadIdx.x]=0;

	d[threadIdx.x]=darg[blockIdx.x][threadIdx.x];

	temp_r[threadIdx.x]=0;	
	temp_r[SIZE_NUMBER + threadIdx.x]=0;	
	temp_r2[threadIdx.x]=0;	
	temp_r2[SIZE_NUMBER + threadIdx.x]=0;	
	temp_u[threadIdx.x]=0;	
	temp_v[threadIdx.x]=0;	
	temp_w[threadIdx.x]=0;	
	
	cy[threadIdx.x]=0;	
	cy[SIZE_NUMBER + threadIdx.x]=0;	

 	//__syncthreads();
	//Now one can begin the computation
	Cuda_Ell_Dbl(xT, zT, xA, zA, Ncst, d, temp_u, temp_v, temp_r, temp_r2, cy);
	Cuda_Ell_Add(xT2, zT2, xA, zA, xB, zB, xC, zC, Ncst, temp_u, temp_v, temp_w, temp_r, temp_r2, cy);
	Cuda_Ell_Add(xA, zA, xT, zT, xA, zA, xA, zA, Ncst, temp_u, temp_v, temp_w, temp_r, temp_r2, cy);
	Cuda_Ell_Add(xT, zT, xT, zT, xT2, zT2, xC, zC, Ncst, temp_u, temp_v, temp_w, temp_r, temp_r2, cy);
	
	//End of the computation; Copy the results for the cpu
	xAarg[blockIdx.x][threadIdx.x]=xA[threadIdx.x];
	zAarg[blockIdx.x][threadIdx.x]=zA[threadIdx.x];
	xBarg[blockIdx.x][threadIdx.x]=xT[threadIdx.x];
	zBarg[blockIdx.x][threadIdx.x]=zT[threadIdx.x];
	xCarg[blockIdx.x][threadIdx.x]=xB[threadIdx.x];
	zCarg[blockIdx.x][threadIdx.x]=zB[threadIdx.x];
}

__global__ void Cuda_Ell_2Add_Dbl_Add(biguint_t *xAarg, biguint_t *zAarg, biguint_t *xBarg, biguint_t *zBarg, biguint_t *xCarg, biguint_t *zCarg, const biguint_t *darg, const int version)
{
	__shared__ unsigned int temp_r[2*SIZE_NUMBER];
	__shared__ unsigned int temp_r2[2*SIZE_NUMBER];
	__shared__ unsigned int temp_u[SIZE_NUMBER];
	__shared__ unsigned int temp_v[SIZE_NUMBER];
	__shared__ unsigned int temp_w[SIZE_NUMBER];

	__shared__ int cy[2*SIZE_NUMBER]; 
	
	__shared__ unsigned int xA[SIZE_NUMBER];
	__shared__ unsigned int zA[SIZE_NUMBER];
	__shared__ unsigned int xB[SIZE_NUMBER];
	__shared__ unsigned int zB[SIZE_NUMBER];
	__shared__ unsigned int xC[SIZE_NUMBER];
	__shared__ unsigned int zC[SIZE_NUMBER];
	__shared__ unsigned int xT[SIZE_NUMBER];
	__shared__ unsigned int zT[SIZE_NUMBER];
	
	__shared__ unsigned int d[SIZE_NUMBER];
	
	//init
	xA[threadIdx.x]=xAarg[blockIdx.x][threadIdx.x];
	zA[threadIdx.x]=zAarg[blockIdx.x][threadIdx.x];
	xB[threadIdx.x]=xBarg[blockIdx.x][threadIdx.x];
	zB[threadIdx.x]=zBarg[blockIdx.x][threadIdx.x];
	xC[threadIdx.x]=xCarg[blockIdx.x][threadIdx.x];
	zC[threadIdx.x]=zCarg[blockIdx.x][threadIdx.x];

	xT[threadIdx.x]=0;
	zT[threadIdx.x]=0;

	d[threadIdx.x]=darg[blockIdx.x][threadIdx.x];

	temp_r[threadIdx.x]=0;	
	temp_r[SIZE_NUMBER + threadIdx.x]=0;	
	temp_r2[threadIdx.x]=0;	
	temp_r2[SIZE_NUMBER + threadIdx.x]=0;	
	temp_u[threadIdx.x]=0;	
	temp_v[threadIdx.x]=0;	
	temp_w[threadIdx.x]=0;	
	
	cy[threadIdx.x]=0;	
	cy[SIZE_NUMBER + threadIdx.x]=0;	

 	//__syncthreads();
	//Now one can begin the computation
	Cuda_Ell_Add(xT, zT, xA, zA, xB, zB, xC, zC, Ncst, temp_u, temp_v, temp_w, temp_r, temp_r2, cy);
	
	if (version==1)
		Cuda_Ell_Add(xC, zC, xT, zT, xA, zA, xB, zB, Ncst, temp_u, temp_v, temp_w, temp_r, temp_r2, cy);
	else 
	{
		Cuda_Ell_Add(xC, zC, xC, zC, xA, zA, xB, zB, Ncst, temp_u, temp_v, temp_w, temp_r, temp_r2, cy);
		xB[threadIdx.x]=xT[threadIdx.x];
		zB[threadIdx.x]=zT[threadIdx.x];
	}
	Cuda_Ell_Dbl(xT, zT, xA, zA, Ncst, d, temp_u, temp_v, temp_r, temp_r2, cy);
	Cuda_Ell_Add(xA, zA, xA, zA, xT, zT, xA, zA, Ncst, temp_u, temp_v, temp_w, temp_r, temp_r2, cy);
	
	//End of the computation; Copy the results for the cpu
	xAarg[blockIdx.x][threadIdx.x]=xA[threadIdx.x];
	zAarg[blockIdx.x][threadIdx.x]=zA[threadIdx.x];
	if (version==1)
	{
		xBarg[blockIdx.x][threadIdx.x]=xC[threadIdx.x];
		zBarg[blockIdx.x][threadIdx.x]=zC[threadIdx.x];
	}
	else
	{
		xBarg[blockIdx.x][threadIdx.x]=xB[threadIdx.x];
		zBarg[blockIdx.x][threadIdx.x]=zB[threadIdx.x];
		xCarg[blockIdx.x][threadIdx.x]=xC[threadIdx.x];
		zCarg[blockIdx.x][threadIdx.x]=zC[threadIdx.x];
		
	}

}

__global__ void Cuda_Ell_Final(biguint_t *xAarg, biguint_t *zAarg, biguint_t *xBarg, biguint_t *zBarg, biguint_t *xCarg, biguint_t *zCarg, const biguint_t *darg)
{
	__shared__ unsigned int temp_r[2*SIZE_NUMBER];
	__shared__ unsigned int temp_r2[2*SIZE_NUMBER];
	__shared__ unsigned int temp_u[SIZE_NUMBER];
	__shared__ unsigned int temp_v[SIZE_NUMBER];
	__shared__ unsigned int temp_w[SIZE_NUMBER];

	__shared__ int cy[2*SIZE_NUMBER]; 
	
	__shared__ unsigned int xA[SIZE_NUMBER];
	__shared__ unsigned int zA[SIZE_NUMBER];
	__shared__ unsigned int xB[SIZE_NUMBER];
	__shared__ unsigned int zB[SIZE_NUMBER];
	__shared__ unsigned int xC[SIZE_NUMBER];
	__shared__ unsigned int zC[SIZE_NUMBER];
	
	__shared__ unsigned int d[SIZE_NUMBER];
	
	//init
	xA[threadIdx.x]=xAarg[blockIdx.x][threadIdx.x];
	zA[threadIdx.x]=zAarg[blockIdx.x][threadIdx.x];
	xB[threadIdx.x]=xBarg[blockIdx.x][threadIdx.x];
	zB[threadIdx.x]=zBarg[blockIdx.x][threadIdx.x];
	xC[threadIdx.x]=xCarg[blockIdx.x][threadIdx.x];
	zC[threadIdx.x]=zCarg[blockIdx.x][threadIdx.x];

	d[threadIdx.x]=darg[blockIdx.x][threadIdx.x];

	temp_r[threadIdx.x]=0;	
	temp_r[SIZE_NUMBER + threadIdx.x]=0;	
	temp_r2[threadIdx.x]=0;	
	temp_r2[SIZE_NUMBER + threadIdx.x]=0;	
	temp_u[threadIdx.x]=0;	
	temp_v[threadIdx.x]=0;	
	temp_w[threadIdx.x]=0;	
	
	cy[threadIdx.x]=0;	
	cy[SIZE_NUMBER + threadIdx.x]=0;	

 	//__syncthreads();
	//Now one can begin the computation

	
	Cuda_Ell_Add(xB, zB, xA, zA, xB, zB, xC, zC, Ncst, temp_u, temp_v, temp_w, temp_r, temp_r2, cy);
	Cuda_Ell_Dbl(xA, zA, xB, zB, Ncst, d, temp_u, temp_v, temp_r, temp_r2, cy);

	//End of the computation; Copy the results for the cpu
	xAarg[blockIdx.x][threadIdx.x]=xA[threadIdx.x];
	zAarg[blockIdx.x][threadIdx.x]=zA[threadIdx.x];
	xBarg[blockIdx.x][threadIdx.x]=xB[threadIdx.x];
	zBarg[blockIdx.x][threadIdx.x]=zB[threadIdx.x];
	xCarg[blockIdx.x][threadIdx.x]=xB[threadIdx.x];
	zCarg[blockIdx.x][threadIdx.x]=zB[threadIdx.x];
	
}

#endif

#ifdef CC20
__device__ void Cuda_Swap(biguint_t A, biguint_t B)
{
	unsigned int temp=A[threadIdx.x];
	A[threadIdx.x]=B[threadIdx.x];
	B[threadIdx.x]=temp;
}

__device__ void Cuda_Circular_Perm(biguint_t A, biguint_t B, biguint_t C)
{
	unsigned int temp=A[threadIdx.x];
	A[threadIdx.x]=B[threadIdx.x];
	B[threadIdx.x]=C[threadIdx.x];
	C[threadIdx.x]=temp;
}
//Compute [3^power3*2^power2]A on the elliptic curve
//prepare for prac, ie set B and C equal to A and set A to [2]A
__global__ void Cuda_Ell_Mul_2_3(biguint_t *xAarg, biguint_t *zAarg, const biguint_t *darg, const unsigned int power2, const unsigned int power3)
{
	__shared__ unsigned int temp_r[2*SIZE_NUMBER];
	__shared__ unsigned int temp_r2[2*SIZE_NUMBER];
	__shared__ unsigned int temp_u[SIZE_NUMBER];
	__shared__ unsigned int temp_v[SIZE_NUMBER];
	__shared__ unsigned int temp_w[SIZE_NUMBER];

	__shared__ int cy[2*SIZE_NUMBER]; 
	
	__shared__ unsigned int xA[SIZE_NUMBER];
	__shared__ unsigned int zA[SIZE_NUMBER];
	__shared__ unsigned int xT[SIZE_NUMBER];
	__shared__ unsigned int zT[SIZE_NUMBER];
	
	__shared__ unsigned int dc[SIZE_NUMBER];
	
	unsigned int i;

	//init
	xA[threadIdx.x]=xAarg[blockIdx.x][threadIdx.x];
	zA[threadIdx.x]=zAarg[blockIdx.x][threadIdx.x];
	xT[threadIdx.x]=0;
	zT[threadIdx.x]=0;

	dc[threadIdx.x]=darg[blockIdx.x][threadIdx.x];

	temp_r[threadIdx.x]=0;	
	temp_r[SIZE_NUMBER + threadIdx.x]=0;	
	temp_r2[threadIdx.x]=0;	
	temp_r2[SIZE_NUMBER + threadIdx.x]=0;	
	temp_u[threadIdx.x]=0;	
	temp_v[threadIdx.x]=0;	
	temp_w[threadIdx.x]=0;	
	
	cy[threadIdx.x]=0;	
	cy[SIZE_NUMBER + threadIdx.x]=0;	

	//Now one can begin the computation

	//printf("power2=%u power3=%u\n",power2,power3);

	for(i=0;i<power2;i++)
	{
		Cuda_Ell_Dbl(xA, zA, xA, zA, Ncst, dc, temp_u, temp_v, temp_r, temp_r2, cy);
	}
	
	for(i=0;i<power3;i++)
	{
		Cuda_Ell_Dbl(xT, zT, xA, zA, dc, Ncst, temp_u, temp_v, temp_r, temp_r2, cy);
		Cuda_Ell_Add(xA, zA, xT, zT, xA, zA, xA, zA, Ncst, temp_u, temp_v, temp_w, temp_r, temp_r2, cy);
	}

	if (threadIdx.x==0 && blockIdx.x==0)
		printf("%u+%u*2^32+%u*2^64+...\n",invmodcst[0],invmodcst[1],invmodcst[2]);
	if (threadIdx.x==0 && blockIdx.x==0)
		printf("%u+%u*2^32+%u*2^64+...\n",dc[0],dc[1],dc[2]);
	if (threadIdx.x==0 && blockIdx.x==0)
		printf("%u+%u*2^32+%u*2^64+...\n",Ncst[0],Ncst[1],Ncst[2]);
	if (threadIdx.x==0 && blockIdx.x==0)
		printf("%u+%u*2^32+%u*2^64+...\n",zAarg[0],zAarg[1],zAarg[2]);
	if (threadIdx.x==0 && blockIdx.x==0)
		printf("%u+%u*2^32+%u*2^64+...\n",zA[0],zA[1],zA[2]);
	//if (threadIdx.x==0)
	//	printf("power2=%u power3=%u\n",power2,power3);
	//End of the computation; Copy the results for the cpu
	xAarg[blockIdx.x][threadIdx.x]=xA[threadIdx.x];
	zAarg[blockIdx.x][threadIdx.x]=zA[threadIdx.x];
}

__global__ void Cuda_PRAC(biguint_t *xAarg, biguint_t *zAarg, const biguint_t *darg, unsigned int PI, unsigned int B1, double val)
{
	__shared__ unsigned int temp_r[2*SIZE_NUMBER];
	__shared__ unsigned int temp_r2[2*SIZE_NUMBER];
	__shared__ unsigned int temp_u[SIZE_NUMBER];
	__shared__ unsigned int temp_v[SIZE_NUMBER];
	__shared__ unsigned int temp_w[SIZE_NUMBER];

	__shared__ int cy[2*SIZE_NUMBER]; 
	
	__shared__ unsigned int xA[SIZE_NUMBER];
	__shared__ unsigned int zA[SIZE_NUMBER];
	__shared__ unsigned int xB[SIZE_NUMBER];
	__shared__ unsigned int zB[SIZE_NUMBER];
	__shared__ unsigned int xC[SIZE_NUMBER];
	__shared__ unsigned int zC[SIZE_NUMBER];
	__shared__ unsigned int xT[SIZE_NUMBER];
	__shared__ unsigned int zT[SIZE_NUMBER];
	__shared__ unsigned int xT2[SIZE_NUMBER];
	__shared__ unsigned int zT2[SIZE_NUMBER];

	__shared__ unsigned int dc[SIZE_NUMBER];
	
	unsigned int e;
	unsigned int d;
	unsigned int r;
	unsigned int pp=PI;

	//init
	xA[threadIdx.x]=xAarg[blockIdx.x][threadIdx.x];
	zA[threadIdx.x]=zAarg[blockIdx.x][threadIdx.x];
	xB[threadIdx.x]=xA[threadIdx.x];
	zB[threadIdx.x]=zA[threadIdx.x];
	xC[threadIdx.x]=xA[threadIdx.x];
	zC[threadIdx.x]=zA[threadIdx.x];
	xT[threadIdx.x]=0;
	zT[threadIdx.x]=0;
	xT2[threadIdx.x]=0;
	zT2[threadIdx.x]=0;

	dc[threadIdx.x]=darg[blockIdx.x][threadIdx.x];

	temp_r[threadIdx.x]=0;	
	temp_r[SIZE_NUMBER + threadIdx.x]=0;	
	temp_r2[threadIdx.x]=0;	
	temp_r2[SIZE_NUMBER + threadIdx.x]=0;	
	temp_u[threadIdx.x]=0;	
	temp_v[threadIdx.x]=0;	
	temp_w[threadIdx.x]=0;	
	
	cy[threadIdx.x]=0;	
	cy[SIZE_NUMBER + threadIdx.x]=0;	

 	//__syncthreads();
	//Now one can begin the computation

	while (pp<=B1)
	{
  	d = PI;
  	r = (unsigned int) ((double) d * val + 0.5);
  	d = PI - r;
  	e = 2 * r - PI;

	
		Cuda_Ell_Dbl(xA, zA, xA, zA, Ncst, dc, temp_u, temp_v, temp_r, temp_r2, cy);
		
		while (d != e)
	  {
	 		if (d < e)
	  	{
		  	r = d;
	    	d = e;
	   		e = r;
				Cuda_Swap(xA,xB);
				Cuda_Swap(zA,zB);
	  	}
	      
			// do the first line of Table 4 whose condition qualifies 
	  	if (d - e <= e / 4 && ((d + e) % 3) == 0)
	  	{ // condition 1 
	  		d = (2 * d - e) / 3;
	    	e = (e - d) / 2;
				Cuda_Ell_Add(xT, zT, xA, zA, xB, zB, xC, zC, Ncst, temp_u, temp_v, temp_w, temp_r, temp_r2, cy);
				Cuda_Ell_Add(xT2, zT2, xT, zT, xA, zA, xB, zB, Ncst, temp_u, temp_v, temp_w, temp_r, temp_r2, cy);
				Cuda_Ell_Add(xB, zB, xB, zB, xT, zT, xA, zA, Ncst, temp_u, temp_v, temp_w, temp_r, temp_r2, cy);
				Cuda_Swap(xA,xT2);
				Cuda_Swap(zA,zT2);
			}
	
	  	else if (d - e <= e / 4 && (d - e) % 6 == 0)
	  	{ // condition 2 
	  		d = (d - e) / 2;
				Cuda_Ell_Add(xB, zB, xA, zA, xB, zB, xC, zC, Ncst, temp_u, temp_v, temp_w, temp_r, temp_r2, cy);
				Cuda_Ell_Dbl(xA, zA, xA, zA, Ncst, dc, temp_u, temp_v, temp_r, temp_r2, cy);
	  	}
	
	  	else if ((d + 3) / 4 <= e)
	  	{ // condition 3 
	  		d -= e;
				Cuda_Ell_Add(xT, zT, xB, zB, xA, zA, xC, zC, Ncst, temp_u, temp_v, temp_w, temp_r, temp_r2, cy);
				//circular permutation
				Cuda_Circular_Perm(xB,xT,xC);
				Cuda_Circular_Perm(zB,zT,zC);
			}
	
	  	else if ((d + e) % 2 == 0)
			{ // condition 4 
	  		d = (d - e) / 2;
				Cuda_Ell_Add(xB, zB, xB, zB, xA, zA, xC, zC, Ncst, temp_u, temp_v, temp_w, temp_r, temp_r2, cy);
				Cuda_Ell_Dbl(xA, zA, xA, zA, Ncst, dc, temp_u, temp_v, temp_r, temp_r2, cy);
			}
	
	  	// now d+e is odd 
	  	else if (d % 2 == 0)
			{ // condition 5 
	  		d /= 2;
				Cuda_Ell_Add(xC, zC, xC, zC, xA, zA, xB, zB, Ncst, temp_u, temp_v, temp_w, temp_r, temp_r2, cy);
				Cuda_Ell_Dbl(xA, zA, xA, zA, Ncst, dc, temp_u, temp_v, temp_r, temp_r2, cy);
			}
	
	  	// now d is odd, e is even 
	  	else if (d % 3 == 0)
			{ // condition 6 
	  		d = d / 3 - e;
				Cuda_Ell_Dbl(xT, zT, xA, zA, Ncst, dc, temp_u, temp_v, temp_r, temp_r2, cy);
				Cuda_Ell_Add(xT2, zT2, xA, zA, xB, zB, xC, zC, Ncst, temp_u, temp_v, temp_w, temp_r, temp_r2, cy);
				Cuda_Ell_Add(xA, zA, xT, zT, xA, zA, xA, zA, Ncst, temp_u, temp_v, temp_w, temp_r, temp_r2, cy);
				Cuda_Ell_Add(xT, zT, xT, zT, xT2, zT2, xC, zC, Ncst, temp_u, temp_v, temp_w, temp_r, temp_r2, cy);
				//circular permutation
				Cuda_Circular_Perm(xC,xB,xT);
				Cuda_Circular_Perm(zC,zB,zT);
			}
	
	  	else if ((d + e) % 3 == 0)
			{ // condition 7 
	  		d = (d - 2 * e) / 3;
				Cuda_Ell_Add(xT, zT, xA, zA, xB, zB, xC, zC, Ncst, temp_u, temp_v, temp_w, temp_r, temp_r2, cy);
				Cuda_Ell_Add(xB, zB, xT, zT, xA, zA, xB, zB, Ncst, temp_u, temp_v, temp_w, temp_r, temp_r2, cy);
				Cuda_Ell_Dbl(xT, zT, xA, zA, Ncst, dc, temp_u, temp_v, temp_r, temp_r2, cy);
				Cuda_Ell_Add(xA, zA, xA, zA, xT, zT, xA, zA, Ncst, temp_u, temp_v, temp_w, temp_r, temp_r2, cy);
			}
	
	  
			else if ((d - e) % 3 == 0)
			{ // condition 8 
	  		d = (d - e) / 3;
				Cuda_Ell_Add(xT, zT, xA, zA, xB, zB, xC, zC, Ncst, temp_u, temp_v, temp_w, temp_r, temp_r2, cy);
				Cuda_Ell_Add(xC, zC, xC, zC, xA, zA, xB, zB, Ncst, temp_u, temp_v, temp_w, temp_r, temp_r2, cy);
				Cuda_Swap(xB,xT);
				Cuda_Swap(zB,zT);
				Cuda_Ell_Dbl(xT, zT, xA, zA, Ncst, dc, temp_u, temp_v, temp_r, temp_r2, cy);
				Cuda_Ell_Add(xA, zA, xA, zA, xT, zT, xA, zA, Ncst, temp_u, temp_v, temp_w, temp_r, temp_r2, cy);
			}
						
			else // necessarily e is even here 
			{ // condition 9 
	  		e /= 2;
				Cuda_Ell_Add(xC, zC, xC, zC, xB, zB, xA, zA, Ncst, temp_u, temp_v, temp_w, temp_r, temp_r2, cy);
				Cuda_Ell_Dbl(xB, zB, xB, zB, Ncst, dc, temp_u, temp_v, temp_r, temp_r2, cy);
			}
		}

		Cuda_Ell_Add(xA, zA, xA, zA, xB, zB, xC, zC, Ncst, temp_u, temp_v, temp_w, temp_r, temp_r2, cy);

	pp*=PI;
	}


//End of the computation; Copy the results for the cpu
	xAarg[blockIdx.x][threadIdx.x]=xA[threadIdx.x];
	zAarg[blockIdx.x][threadIdx.x]=zA[threadIdx.x];
}
#endif

__host__ void cuda_copy_cst(biguint_t h_N, biguint_t h_invmod)
{
	cudaMemcpyToSymbol(invmodcst,h_invmod,sizeof(biguint_t));
	cudaMemcpyToSymbol(Ncst,h_N,sizeof(biguint_t));
}
