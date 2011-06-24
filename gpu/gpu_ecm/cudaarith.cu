#include "main.h"
#include "cudaarith.h"

#define __add_cc(r,a,b) __asm__("add.cc.u32 %0,%1, %2;" :"=r"(r):"r"(a),"r"(b)) 
#define __addc_cc(r,a,b) __asm__("addc.cc.u32 %0,%1, %2;":"=r"(r):"r"(a),"r"(b)) 
//#define __addc(r,a,b) __asm__("addc.u32 %0,%1, %2;" :"=r"(r): "r"(a), "r"(b)) 

#define __sub_cc(r,a,b) __asm__("sub.cc.u32 %0,%1, %2;" :"=r"(r):"r"(a),"r"(b)) 
//#define __subc_cc(r,a,b) __asm__("subc.cc.u32 %0,%1, %2;":"=r"(r):"r"(a),"r"(b)) 
//#define __subc(r,a,b) __asm__("subc.u32 %0,%1, %2;" :"=r"(r): "r"(a), "r"(b)) 

#define __addcy(carry) __asm__("addc.s32 %0, 0, 0;" :"=r"(carry)) 
#define __addcy2(carry) __asm__("addc.s32 %0, %0, 0;" :"+r"(carry)) 

#define __subcy(carry) __asm__("subc.s32 %0, 0, 0;" :"=r"(carry)) 
#define __subcy2(carry) __asm__("subc.s32 %0, %0, 0;" :"+r"(carry)) 

#define __mul(h,l,a,b) __asm__("mul.hi.u32 %0,%2,%3;\n\t" "mul.lo.u32 %1,%2,%3;" : "=r"(h), "=r"(l) : "r"(a), "r"(b))

#define __mad_lo(r,a,b,c) __asm__("mad.lo.u32 %0,%1,%2,%3;" : "=r"(r) : "r"(a), "r"(b), "r"(c))
#define __mad_hi(r,a,b,c) __asm__("mad.hi.u32 %0,%1,%2,%3;" : "=r"(r) : "r"(a), "r"(b), "r"(c))

//specific functions for compute capability 1.3
//#ifdef CC13
//cy[i]=0 or 1
__device__ int Cuda_Cmp2
(const biguint_t A, const dbigint_t cy, const biguint_t B)
{
	int i;
	int r=0;
	for (i = SIZE_NUMBER-1;i>0;i--)
	{
	  /*
		si r=0 ne rien faire
		si r=1 et A[i]=TWO32-1
			si B[i]=TWO32-1 (?=0) (?dépend pas de la retenue)
				on passe a l'itération suivante avec r=1
			sinon
				A>B return 1
		si r=1 et A[i]!=TWO32-1
			aucune retenue ne peut se propager donc égalité au niveau d'avant
			on fait cette itération normalement avec r=0
		si r=-1 et A[i]=TWO32-1
			une retenue peut se propager 
			Si B[i]==0
				suivant la retenue
			Si B[i]>0
				return -1 car si pas retenus de propagé alors c'est vrai au niveau
				précedent et si retenue propages c'est vrai pour ce niveau
		si r=-1 et A[i]!=TWO32-1
			aucune retenue ne peut se propager 
			return -1
	*/
		if (r==1 && A[i]==TWO32-1 && cy[i-1]==1)
			return 1;
		else if (r==i-1 && A[i]!=TWO32-1)
			return -1;
		else if (r==1 && A[i]==TWO32-1 && cy[i-1]==0)
		{
			//assume no propagation
			if (B[i]==TWO32-1)
				r=0;
			else
				return 1;
		}
		else if (r==-1 && A[i]==TWO32-1 && cy[i-1]==0)
		{
			//assume no propagation
			return -1;
		}
		else //r=0 or (r=1 and A[i]!=2^32-1) or (r=-1 and A[i]=2^32-1 and cy[i-1]=1)
		{
			if (A[i]+cy[i-1] > B[i])
				return 1;
			else if (A[i]+cy[i-1]==B[i])
				r=(cy[i-1]==0)?1:0;
			else if (A[i]+1 < B[i])
				return -1;
			else //it happens only for cy[i-1]=0
				r=-1;
		}
	}
	// traiter le cas i=0
	return 0;
}


//  (A > B)?, returns 1(true), -1(false) or 0(a=b) 
//Assume A and B are normalize (no carry or borrow)
__device__ int Cuda_Cmp(const biguint_t A, const biguint_t B)
{
	int i;
	for (i = SIZE_NUMBER-1;i>=0;i--)
	{
		if (A[i] > B[i])
			return 1;
		else if (A[i] < B[i])
			return -1;
	}
	return 0;
}

/*
//Assume all carries are negative
//return -1 if A+cy*2^32 < 0 else return 1
__device__ int Cuda_lt_zero(const biguint_t A, const dbigint_t cy)
{
	int i = SIZE_NUMBER-1;
	if (cy[i]<0)
		return -1;
	do
	{
		if (A[i]< -cy[i-1])
			return -1;
			i--;
	}while(i!=0);
	return 1;
}
*/
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

//Assume cy[threadIdx.x] = 0,+/-1
__device__ void Cuda_Normalize(biguint_t A,dbigint_t cy)
{
	int cytemp;
	cytemp = cy[threadIdx.x];
	cy[threadIdx.x]=0;
	int tmp=threadIdx.x+1 % SIZE_NUMBER;

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

/*
__device__ void Cuda_Normalize_64(dbiguint_t A,dbigint_t cy)
{
	__add_cc(A[threadIdx.x+1],A[threadIdx.x+1],cy[threadIdx.x]);
	__addcy(cy[(threadIdx.x+1)]); //In the end of RedMont cy[0] is always 0
}
*/
/*
__device__ void Cuda_Fully_Normalize_64(biguint_t A,dbigint_t cy)
{
	do
	{
	Cuda_Normalize_64(A,cy);
	}while(__any(cy[threadIdx.x]|cy[threadIdx.x+SIZE_NUMBER])!=0);
}
*/

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
  
	if (Cuda_Cmp (Rmod, Ncst) >= 0)// (Rmod >= N)? 
	//if (Cuda_Cmp2 (Rmod, cy, Ncst) >= 0)// (Rmod >= N)? 
	{
   	Cuda_Sub (Rmod, cy, Rmod, Ncst); 
		Cuda_Fully_Normalize(Rmod, cy);	
	}
}

//Compute Rmod  <-Rmod + A [mod] 
__device__ void Cuda_Add_mod
(biguint_t Rmod, dbigint_t cy, const biguint_t A)
{
 	Cuda_Add(Rmod, cy, Rmod, A);
	Cuda_Fully_Normalize(Rmod, cy);

	if (Cuda_Cmp (Rmod, Ncst) >= 0)// (Rmod >= N)? 
	{
   	Cuda_Sub (Rmod, cy, Rmod, Ncst);  
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
		Cuda_Add (Rmod, cy, Rmod, Ncst);
		Cuda_Subc (Rmod, cy, Rmod, A);
		Cuda_Fully_Normalize(Rmod, cy);	
	}
}
/*
__device__ void Cuda_Mul
(dbiguint_t R, dbigint_t cy, const biguint_t A,const biguint_t B)
{
	int i;
	unsigned int h,l;

	unsigned int temp=A[threadIdx.x];

	R[threadIdx.x]=0;
	R[threadIdx.x+SIZE_NUMBER]=0;

	for (i=0;i<SIZE_NUMBER;i++)
	{
		//h*2^32+l =A[i]*B[threadIDx.x]
		//Cuda_Mul_uint(&h,&l,temp,B[i]);
		__mul(h,l,temp,B[i]);
		
		//R[i+threadIdx.x] +=l;
		//cy[i+threadIdx.x]+=(R[i+threadIdx.x] < l);
		//if (R[i+threadIdx.x] < l)
		//	cy[i+threadIdx.x]++;
		
		//__add_cc(R[i+threadIdx.x],R[i+threadIdx.x],l);
		//__addcy2(cy[i+threadIdx.x]);
		
		//R[i+1+threadIdx.x] +=h;
		//cy[i+1+threadIdx.x]+=(R[i+1+threadIdx.x]<h);
		//if (R[i+1+threadIdx.x]<h)
		//	cy[i+1+threadIdx.x]++;
		
		//__add_cc(R[i+1+threadIdx.x],R[i+1+threadIdx.x],h);
		//__addcy2(cy[i+1+threadIdx.x]);

		__add_cc(R[i+threadIdx.x],R[i+threadIdx.x],l);
		__addc_cc(R[i+1+threadIdx.x],R[i+1+threadIdx.x],h);
		__addcy2(cy[i+1+threadIdx.x]);
	}
	//Cuda_Fully_Normalize_64(R,cy);
}
*/
/*
__device__ void Cuda_RedMont_Step (dbiguint_t r, dbigint_t cy)
{
	unsigned int h,l;
	unsigned int tmp; 
	//unsigned int tmp,tmp2;

	//h*2^32+l =A[i]*B[threadIDx.x]
	//Cuda_Mul_uint(&h,&l,invmodcst[0]*r[0],Ncst[threadIdx.x]);
	__mul(h,l,invmodcst[0]*r[0],Ncst[threadIdx.x]);
	
	//r[threadIdx.x]+=l;
	//cy[threadIdx.x]+=(r[threadIdx.x]<l);
 
 	//r[threadIdx.x+1]+=h;
	//cy[threadIdx.x+1]+=(r[1+threadIdx.x]<h);
		
	__add_cc(r[threadIdx.x],r[threadIdx.x],l);
	__addc_cc(r[threadIdx.x+1],r[threadIdx.x+1],h);
	__addcy2(cy[threadIdx.x+1]);

	
	//Normalize only to add the last carry
	//Cuda_Normalize(r,cy);
	
	//tmp=cy[threadIdx.x];
	//cy[threadIdx.x]=0;
	//r[threadIdx.x+1]+=tmp;
	//cy[threadIdx.x+1]+=(r[threadIdx.x+1]<tmp);//+= is mandatory, can't put only =

	//make one round of normalize + a right shift
	//__addcy(cy[threadIdx.x]);


	__add_cc(r[threadIdx.x],r[threadIdx.x+1],cy[threadIdx.x]);
	tmp=(threadIdx.x==SIZE_NUMBER-1)?cy[threadIdx.x+1]:0;
	__asm__("addc.u32 %0,%1, 0;" :"=r"(cy[threadIdx.x]): "r"(tmp)); 

	//cy[threadIdx.x]+=tmp;
	//r[threadIdx.x]=r[threadIdx.x+1];
	//cy[threadIdx.x]=cy[threadIdx.x+1];

	r[SIZE_NUMBER+threadIdx.x]=r[SIZE_NUMBER+threadIdx.x+1];
	cy[SIZE_NUMBER+threadIdx.x]=cy[SIZE_NUMBER+threadIdx.x+1];
	if (threadIdx.x==0)
	{
		cy[2*SIZE_NUMBER-1]=0;
		r[2*SIZE_NUMBER-1]=0;
	}

}
*/
__device__ void Cuda_Mulmod_step
(dbiguint_t r,dbigint_t cy, unsigned int a, unsigned int b)
{
	unsigned int h,l;
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

	__mul(h,l,invmodcst[0]*r[0],Ncst[threadIdx.x]);
	__add_cc(r[threadIdx.x],r[threadIdx.x],l);
	__addc_cc(r[threadIdx.x+1],r[threadIdx.x+1],h);
	__addcy2(cy[threadIdx.x+1]);
 
	//make one round of normalize + a right shift
	__add_cc(r[threadIdx.x],r[threadIdx.x+1],cy[threadIdx.x]);
	tmp=(threadIdx.x==SIZE_NUMBER-1)?cy[threadIdx.x+1]:0;
	__asm__("addc.s32 %0,%1, 0;" :"=r"(cy[threadIdx.x]): "r"(tmp)); 

	if (threadIdx.x==0)
	{
		cy[SIZE_NUMBER]=0;
		r[SIZE_NUMBER]=0;
	}
}

/*
//Assume r<N^2
__device__ void Cuda_RedMontgomery_V2 
(biguint_t mul, dbiguint_t r, dbigint_t cy)
{
	int i;
	for (i=0;i<SIZE_NUMBER;i++)
	{
		//__syncthreads();
		Cuda_RedMont_Step(r,cy);
		//__syncthreads();
	}

	Cuda_Fully_Normalize_64(r,cy);
	//Cuda_Fully_Normalize(r,cy);

	if (Cuda_Cmp (r,Ncst) >= 0) // mul >= N 
	{
  	Cuda_Sub (mul, cy, r, Ncst); 
		Cuda_Fully_Normalize(mul,cy);	
	}
	else
	{
		mul[threadIdx.x]=r[threadIdx.x];
	}
}
*/
__device__ void Cuda_Dbl_mod
(biguint_t r, dbigint_t cy, biguint_t a)
{
	//cy[threadIdx.x]=(a[threadIdx.x]>>31);
	//r[threadIdx.x]=a[threadIdx.x]<<1;
	__add_cc(r[threadIdx.x],a[threadIdx.x],a[threadIdx.x]);
	__addcy2(r[(threadIdx.x+1)%SIZE_NUMBER]);

	//Cuda_Normalize_test(r,cy);	

	if (Cuda_Cmp (r,Ncst) >= 0) 
	{
  	Cuda_Sub (r, cy, r, Ncst); 
		Cuda_Fully_Normalize(r, cy);	
	}
}


__device__ void Cuda_Mulint_mod(dbiguint_t r,dbigint_t cy, biguint_t A, unsigned int b)
{
	unsigned int h,l;
	__mul(h,r[threadIdx.x],A[threadIdx.x],b);
	__add_cc(r[threadIdx.x+1],r[threadIdx.x+1],h);
	__addcy(cy[threadIdx.x+1]);

	//h*2^32+l =A[i]*B[threadIDx.x]
	__mul(h,l,invmodcst[0]*r[0],Ncst[threadIdx.x]);
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
		//cy[threadIdx.x+SIZE_NUMBER]=0;
		//r[threadIdx.x+SIZE_NUMBER]=0;
		cy[SIZE_NUMBER]=0;
		r[SIZE_NUMBER]=0;
	}
	
	Cuda_Fully_Normalize(r,cy);	

	if (Cuda_Cmp (r,Ncst) >= 0) 
	{
  	Cuda_Sub (r, cy, r, Ncst); 
		Cuda_Fully_Normalize(r,cy);	
	}
}

__device__ void Cuda_Mul_mod (biguint_t mul, dbigint_t cy, const biguint_t A,const biguint_t B, dbiguint_t r)
{

	int i;
	unsigned int temp=A[threadIdx.x];

	r[threadIdx.x]=0;
	//r[threadIdx.x+SIZE_NUMBER]=0;
	
	for (i=0;i<SIZE_NUMBER;i++)
		Cuda_Mulmod_step(r,cy,temp,B[i]);

	Cuda_Fully_Normalize(r,cy);

	if (Cuda_Cmp (r,Ncst) >= 0) // mul >= N 
	{
  	Cuda_Sub (mul, cy, r, Ncst); 
		Cuda_Fully_Normalize(mul,cy);	
	}
	else
		mul[threadIdx.x]=r[threadIdx.x];
}
__device__ void Cuda_Square_mod (biguint_t mul, dbigint_t cy, const biguint_t A, dbiguint_t r)
{
	Cuda_Mul_mod(mul,cy,A,A,r);
}

#ifndef TEST
__global__ void Cuda_Ell_DblAdd(biguint_t *xarg, biguint_t *zarg, biguint_t *x2arg, biguint_t *z2arg, unsigned int firstinvd)
{
	__shared__ unsigned int b_temp_r[CURVES_BY_BLOCK][SIZE_NUMBER+1];
	__shared__ int b_cy[CURVES_BY_BLOCK][SIZE_NUMBER+1]; 

	__shared__ unsigned int b_t[CURVES_BY_BLOCK][SIZE_NUMBER];
	__shared__ unsigned int b_u[CURVES_BY_BLOCK][SIZE_NUMBER];
	__shared__ unsigned int b_v[CURVES_BY_BLOCK][SIZE_NUMBER];
	__shared__ unsigned int b_w[CURVES_BY_BLOCK][SIZE_NUMBER];
	
	volatile unsigned int idx1=blockIdx.x*blockDim.y+threadIdx.y;
	//volatile unsigned int t1=threadIdx.x+1;
	//volatile unsigned int t2=threadIdx.x+SIZE_NUMBER;
	
	unsigned int *t=b_t[threadIdx.y];
	unsigned int *u=b_u[threadIdx.y];
	unsigned int *v=b_v[threadIdx.y];
	unsigned int *w=b_w[threadIdx.y];
	unsigned int *temp_r=b_temp_r[threadIdx.y];
	int *cy=b_cy[threadIdx.y];

	//init
	b_cy[threadIdx.y][threadIdx.x]=0;	
	if (threadIdx.x==0)
		b_cy[threadIdx.y][SIZE_NUMBER]=0;	

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
#endif

#ifdef TEST
__global__ void Cuda_Test(biguint_t *Aarg, biguint_t *Barg)
{
	__shared__ int b_cy[CURVES_BY_BLOCK][2*SIZE_NUMBER]; 
	__shared__ unsigned int b_r[CURVES_BY_BLOCK][2*SIZE_NUMBER]; 
	__shared__ unsigned int b_b[CURVES_BY_BLOCK][SIZE_NUMBER]; 

	unsigned int idx1=blockIdx.x*blockDim.y+threadIdx.y;
	unsigned int *b=b_b[threadIdx.y];
	unsigned int *r=b_r[threadIdx.y];
	int *cy=b_cy[threadIdx.y];
	r[SIZE_NUMBER+threadIdx.x]=Aarg[idx1][threadIdx.x];
	b[threadIdx.x]=Barg[idx1][threadIdx.x];
		
	cy[threadIdx.x]=0;
	cy[threadIdx.x+SIZE_NUMBER]=0;

	Cuda_Mul_mod(r,cy,r+SIZE_NUMBER,b,r);

	Aarg[idx1][threadIdx.x]=r[threadIdx.x];
	Barg[idx1][threadIdx.x]=cy[threadIdx.x];
}
/*
__global__ void Cuda_Test(biguint_t *Aarg, biguint_t *Barg)
{
	__shared__ unsigned int r[CURVES_BY_BLOCK][2*SIZE_NUMBER];
	__shared__ int cy[CURVES_BY_BLOCK][2*SIZE_NUMBER]; 
	__shared__ unsigned int A[CURVES_BY_BLOCK][SIZE_NUMBER]; 
	__shared__ unsigned int B[CURVES_BY_BLOCK][SIZE_NUMBER]; 
	unsigned int idx1=blockIdx.x*blockDim.y+threadIdx.y;
	int i;

	//init
	cy[threadIdx.y][threadIdx.x]=0;	
	cy[threadIdx.y][SIZE_NUMBER + threadIdx.x]=0;	
	r[threadIdx.y][threadIdx.x]=Barg[idx1][threadIdx.x];
	r[threadIdx.y][SIZE_NUMBER + threadIdx.x]=0;	
	A[threadIdx.y][threadIdx.x]=Aarg[idx1][threadIdx.x];
	B[threadIdx.y][threadIdx.x]=Barg[idx1][threadIdx.x];

	//Cuda_Mul(r[threadIdx.y],cy[threadIdx.y],A[threadIdx.y],B[threadIdx.y]);
	//Cuda_RedMontgomery_V2(A[threadIdx.y],r[threadIdx.y],cy[threadIdx.y]);
	__syncthreads();
	Cuda_RedMont_Step(r[threadIdx.y],cy[threadIdx.y]);
	__syncthreads();
	
	
	//__syncthreads();
	//for (i=0;i<32;i++)
//	{
		//Cuda_RedMontgomery_V2(r[threadIdx.y],r[threadIdx.y],cy[threadIdx.y]);
		//Cuda_RedMont_Step(r[threadIdx.y],cy[threadIdx.y]);
		//__syncthreads();
		//Cuda_Fully_Normalize(r[threadIdx.y],cy[threadIdx.y]);
	//}
	//__syncthreads();
	//End of the computation; Copy the results for the cpu
	Aarg[idx1][threadIdx.x]=r[threadIdx.y][threadIdx.x];
	Barg[idx1][threadIdx.x]=r[threadIdx.y][threadIdx.x+32];
}*/
#endif
	
__host__ void cuda_copy_cst(biguint_t h_N, biguint_t h_invmod)
{
	cudaMemcpyToSymbol(invmodcst,h_invmod,sizeof(biguint_t));
	cudaMemcpyToSymbol(Ncst,h_N,sizeof(biguint_t));
}
