#include <stdio.h>

#define WARP_SIZE 32

#define GROUPSET 16
#define NUMFACES 3

#define fouralpha 1.82
#define fouralpha4 5.82


#define Connect(a,b,c) Connect[ a + 3 * ( b + mC * c ) ]

__device__ __forceinline__ double shfl_d(double var,int lane)
{ float lo, hi;
  asm volatile("mov.b64 {%0,%1}, %2;" : "=f"(lo), "=f"(hi) : "d"(var));
  hi = __shfl(hi, lane);
  lo = __shfl(lo, lane);
  asm volatile("mov.b64 %0, {%1,%2};" : "=d"(var) : "f"(lo), "f"(hi));
  return var;
}


extern "C"
{


  __global__ void GPU_sweep(
			    int  size_maxCorner,
			    int  size_maxcf,
			    int  nAngle,
			    int  nzones,
			    int  ncornr,
			    int  Groups,
			    int  nbelem,
			    int* AngleOrder,
			    double* soa_omega,
			    int* nextZ,
			    int* next,
			    int* soa_nCorner,
			    int* soa_nCFaces,
			    int* soa_c0,
			    double* soa_STotal,
			    double* STimeBatch,
			    double* soa_SigtInv,
			    double* soa_Volume,
			    double* soa_Sigt,
			    double* soa_A_fp,
			    double* soa_A_ez,
			    int* soa_Connect,
			    double* psic,
			    double* psib,
			    double* omega_A_fp,
			    double* omega_A_ez,
			    int* soa_Connect_ro,
			    int* passZ,
			    bool calcSTime,
			    double tau
			    )
  {

    //   double omega[3];
    int c,ig,i,icface,ifp,cez,k;
    //   double Q[Groups * size_maxCorner];
    //   double src[Groups * size_maxCorner];
    //   double SigtVol[Groups * size_maxCorner];
    //   double afpm[size_maxcf];
    //   double psifp[Groups * size_maxcf];
    //   int    ez_exit[size_maxcf];
    //   double coefpsic[size_maxcf];
    //   double tpsic[Groups * size_maxCorner];
    //   double psi_opp[Groups];
    double area_opp,area_inv,sumArea;
    double r_psifp;

    double psi_opp,tpsic,r_afpm;
    double Q[8];
    double src[8];

    //double volume[8];
    //double coefpsic_stk[3];
    //double psifp[3];
    //int ez_exit[3];
  
    //double *src;
    volatile double *volume;
    volatile double *coefpsic;
    volatile double *psifp;
    volatile int *ez_exit;
    //__shared__ volatile double sm_agg[12*128];  // 4x32 thread per tb. 8tb. 6KB
    extern __shared__ double sm_agg[];  
   
    int offset = (8+3+3*WARP_SIZE+3)*threadIdx.y;
    volume = &(sm_agg[offset]);  //8 doubles  
    offset += size_maxCorner;

    coefpsic = &(sm_agg[offset]); // 3 doubles
    offset += size_maxcf;

    psifp = &(sm_agg[offset]); // 3 x warp size doubles
    offset += size_maxcf * WARP_SIZE;

    //note ez_exit has integer type
    ez_exit = (int*) &(sm_agg[offset]); // 3 int
    //   for(int Angle=0;Angle<nAngle;Angle++)

    //   const double fouralpha = 1.82;
    //   const double fouralpha4 = 5.82;
   
#define soa_omega(a,b) soa_omega[a + 3 * b]

    //   #define tpsic(ig,c) tpsic[ (ig) + Groups * (c)]
#define EB_ListExit(a,ia) EB_ListExit[ a + 2 * (ia) ]
#define soa_A_fp(a,icface,c,zone) soa_A_fp[ a + 3 * ( icface + size_maxcf * ( c + size_maxCorner * (zone) ) )]
#define soa_A_ez(a,icface,c,zone) soa_A_ez[ a + 3 * ( icface + size_maxcf * ( c + size_maxCorner * (zone) ) )]
#define omega_A_fp(icface,c,zone) omega_A_fp[ ( icface + size_maxcf * ( c + size_maxCorner * (zone) ) )]
#define omega_A_ez(icface,c,zone) omega_A_ez[ ( icface + size_maxcf * ( c + size_maxCorner * (zone) ) )]
#define soa_Connect(a,icface,c,zone) soa_Connect[ a + 3 * ( icface + size_maxcf * ( c + size_maxCorner * (zone) ) )]
   
#define psifp(ig,jf) psifp[(ig) + WARP_SIZE * (jf)]
#define psib(ig,b,c) psib[(ig) + Groups * ((b) + nbelem * (c) )]
#define psic(ig,b,c) psic[(ig) + Groups * ((b) + ncornr *(c) )]
   
#define Q(ig,c) Q[(ig) + WARP_SIZE * (c)]
#define src(ig,c) src[c]
#define soa_Sigt(ig,zone) soa_Sigt[(ig) + Groups * (zone)]
#define soa_Volume(c,zone) soa_Volume[c + size_maxCorner * (zone)]
#define soa_SigtInv(ig,zone) soa_SigtInv[(ig) + Groups * (zone)]
#define soa_STotal(ig,c,zone) soa_STotal[ig + Groups * ( c + size_maxCorner * (zone) )]
#define STimeBatch(ig,ic,Angle) STimeBatch[ig + Groups * ( (ic) + ncornr * (Angle) ) ]
#define nextZ(a,b) nextZ[ (a) + nzones * (b) ]
#define next(a,b) next[ (a) + (ncornr+1)  * (b) ]

    //int mm = blockIdx.x;
    int Angle = AngleOrder[blockIdx.x]-1;
    ig = threadIdx.x;
 

    omega_A_fp += blockIdx.x * nzones * size_maxcf * size_maxCorner;
    omega_A_ez += blockIdx.x * nzones * size_maxcf * size_maxCorner;
    passZ      += Angle * nzones;
   
    const int group_offset=blockIdx.y * WARP_SIZE; //should be blockDim.x instead of warpsize?

    //  if (!( group_offset + threadIdx.x < Groups )) return;
    
    psib += group_offset;
    psic += group_offset;
    soa_Sigt += group_offset;
    soa_STotal += group_offset;
    soa_SigtInv += group_offset;
    
    STimeBatch += group_offset;

    int ndone = 0;
    int ndoneZ = 0;
    // hyperplane number p
    int p=0;

    while(ndoneZ < nzones)
    {
      //increment hyperplane
      p++;
      // get number of zones in this hyperplane
      int passZcnt = passZ[p] - passZ[p-1];

      // for(int ii=threadIdx.y+blockDim.y*blockIdx.z; ii<passZcnt;ii+=blockDim.y*gridDim.z)
      //for(int ii=blockIdx.z; ii<passZcnt;ii+=gridDim.z)
      for(int ii=threadIdx.y; ii<passZcnt;ii+=blockDim.y)
      {
	ndone = ( ndoneZ + ii ) * size_maxCorner;
 
	// get the zone (minus 1 so it is valid c index)
	int zone = nextZ(ndoneZ+ii,Angle) - 1;
  
  
	int nCorner   = soa_nCorner[zone];
	int nCFaces   = soa_nCFaces[zone];
	int c0        = soa_c0[zone] ;
  
	double Sigt = soa_Sigt(ig,zone);
	double r_soa_SightInv = soa_SigtInv(ig,zone);
	double r_omega_A_fp;
	double r_omega_A_ez;
	int connect0,connect1,connect2;
  

	// coallesced loads into shared memory
	if(threadIdx.x<nCorner) volume[threadIdx.x] = soa_Volume(threadIdx.x,zone);

	// different threads hold values for different icface in registers instead of shared memory
	// other threads can access the register values via a shuffle command. 
	//But now with only 16 threads (groups) there are not threads to hold nCorner*nCFaces (3*8)
	if(threadIdx.x<nCorner*nCFaces) 
	{
          int cc = size_maxcf * size_maxCorner;
          r_omega_A_fp = omega_A_fp[threadIdx.x + cc * zone];
          r_omega_A_ez = omega_A_ez[threadIdx.x + cc * zone];
          connect0     = soa_Connect_ro[threadIdx.x + cc*(0 + 3*zone)];
          connect1     = soa_Connect_ro[threadIdx.x + cc*(1 + 3*zone)];
          connect2     = soa_Connect_ro[threadIdx.x + cc*(2 + 3*zone)];
	}
   
	//if(nCorner*nCFaces>blockDim.x){printf("Error: threads are not covering nCorner*nCFaces\n");abort;}

  
	for(c=0;c<nCorner;c++)
	{
	  double source;
	  //if(!calcSTime) 
	  //{
	    source = soa_STotal(ig,c,zone) + STimeBatch(ig,c0+c,blockIdx.x);
	    //}
	    //else  // first temp and flux iteration: compute STime, use it, and zero copy back to host.
		//{
	    //double STime_temp = tau*psic(ig,c0+c,blockIdx.x);
	    //source = soa_STotal(ig,c,zone) + STime_temp;
	    //STime(ig,c0+c,Angle) = STime_temp;
	    //}

	  Q[c]       = r_soa_SightInv *source ;
	  //src(ig,c)     = soa_Volume(c,zone) *source;
	  //volume[c] = soa_Volume(c,zone);
	  src(ig,c)     = volume[c]*source; // really just src[c]
	  //SigtVol(ig,c) = soa_Sigt(ig,zone)*soa_Volume(c,zone);
	}
  
	for(i=0;i<nCorner;i++)
	{
  
	  int ic      = next(ndone+i,Angle);
	  c       = ic - c0 - 1;
  
	  sumArea = 0.0;
   
	  for(icface=0;icface<nCFaces;icface++)
	  {
	    //afpm[icface] = omega_A_fp(icface,c,zone);  
	    r_afpm = shfl_d(r_omega_A_fp,icface+size_maxcf*c);
  
	    //  if ( Angle == 1 && ig==0 && zone == 1 )
	    //    printf("a=%d,c=%d,icface=%d,afpm=%e\n",Angle,c,icface,r_afpm);
  
	    //       int icfp    = soa_Connect(0,icface,c,zone) - 1;
	    //       int ib      = soa_Connect(1,icface,c,zone) - 1;
	    int icfp= __shfl(connect0,icface+size_maxcf*c) - 1;
	    int ib= __shfl(connect1,icface+size_maxcf*c) - 1;
  
	    if ( r_afpm >= 0.0 )
	    { 
	      sumArea = sumArea + r_afpm;
	    }
	    else
	    {
	      if (icfp == -1)
	      {
		//             psifp(ig,icface) = psib(ig,ib,Angle);
		r_psifp = psib(ig,ib,blockIdx.x);
	      }
	      else
	      {
		//             psifp(ig,icface) = psic(ig,icfp,Angle);
		r_psifp = psic(ig,icfp,blockIdx.x);            
	      }
  
	      src(ig,c)  -= r_afpm*r_psifp;
	      psifp(ig,icface) = r_psifp;
	      //psifp[icface] = r_psifp;
	    }

	  }
    
	  int nxez = 0;
  
	  for(icface=0;icface<nCFaces;icface++)
	  {
  
	    //double aez = omega_A_ez(icface,c,zone);
	    double aez = shfl_d(r_omega_A_ez,icface+size_maxcf*c);
  
	    if (aez > 0.0 )
	    {
  
	      sumArea        = sumArea + aez;
	      area_opp       = .0;
	      //           cez            = soa_Connect(2,icface,c,zone) - 1;
	      cez            = __shfl(connect2,icface+size_maxcf*c) - 1;
	      ez_exit[nxez]  = cez;
	      coefpsic[nxez] = aez;
	      nxez           = nxez + 1;
  
	      if (nCFaces == 3)
	      {
  
		ifp = (icface+1)%nCFaces;
		r_afpm = shfl_d(r_omega_A_fp,ifp+size_maxcf*c);
  
		if ( r_afpm < 0.0 )
		{ 
		  area_opp   = -r_afpm;
		  psi_opp =  psifp(ig,ifp);
		  //psi_opp =  psifp[ifp];
		}
	      }
	      else
	      {
  
		ifp        = icface;
		area_opp   = 0.0;
		psi_opp = 0.0;
  
		for(k=0;k<nCFaces-2;k++)
		{
		  ifp = (ifp+1)%nCFaces;
		  r_afpm = shfl_d(r_omega_A_fp,ifp+size_maxcf*c);
		  if ( r_afpm < 0.0 )
		  {
		    area_opp   = area_opp   - r_afpm;
		    psi_opp = psi_opp - r_afpm*psifp(ig,ifp);
		    //psi_opp = psi_opp - r_afpm*psifp[ifp];
		  }
		}
  
		area_inv = 1.0/area_opp;
  
		psi_opp = psi_opp*area_inv;
  
	      }
  
	      if (area_opp > 0.0) {
  
		double aez2 = aez*aez;
  
		{
    
		  double sigv         = Sigt*volume[c];
		  double sigv2        = sigv*sigv;
		  double gnum         = aez2*( fouralpha*sigv2 +       aez*(4.0*sigv + 3.0*aez) );
		  double gtau         = gnum/( gnum + 4.0*sigv2*sigv2 + aez*sigv*(6.0*sigv2 + 2.0*aez*(2.0*sigv + aez)) ) ;
		  double sez          = gtau*sigv*( psi_opp - Q[c] ) +   0.5*aez*(1.0 - gtau)*( Q[c] - Q[cez] );
  
		  src(ig,c)    = src(ig,c)   + sez;
		  src(ig,cez)  = src(ig,cez) - sez;
  
		}
  
	      }
	      else
	      {
		double sez          = 0.5*aez*( Q[c] - Q[cez] );
		src(ig,c)    = src(ig,c)   + sez;
		src(ig,cez)  = src(ig,cez) - sez;
  
	      } 
	    }
	  }
  
  
	  tpsic = src(ig,c)/(sumArea + Sigt*volume[c]);
  
  
	  for(icface=0;icface<nxez;icface++)
	  {
	    int cez   = ez_exit[icface];
	    src(ig,cez) = src(ig,cez) + coefpsic[icface]*tpsic;
	  }
  
	  //hope that ther is no self referencing
	  psic(ig,c0+c,blockIdx.x) = tpsic;
	  //psibatch(ig,c0+c,mm)= tpsic;

	} //end of corner 
  
      } //end of zone loop 
      ndoneZ += passZcnt;
      __syncthreads();
    } //end of while

  }






  __global__ void GPU_fp_ez_hplane(
          int  size_maxCorner,
          int  size_maxcf,
          int  nzones,
          int  ncornr,
          int  Groups,
          int  nbelem,
          int* AngleOrder,
       double* soa_omega,
          int* nextZ,
          int* next,
          int* soa_nCorner,
          int* soa_nCFaces,
          int* soa_c0,
       double* soa_A_fp,
       double* soa_A_ez,
       double* omega_A_fp,
       double* omega_A_ez,
          int* soa_Connect,
          int* soa_Connect_ro,
	  int* passZ)
  {

    //int c,i,ig,icface,ii;
    int c,i,icface;


   
#define soa_omega(a,b) soa_omega[a + 3 * b]
#define omega_A_fp(icface,c,zone) omega_A_fp[  ( icface + size_maxcf * ( c + size_maxCorner * (zone) ) )]
#define omega_A_ez(icface,c,zone) omega_A_ez[  ( icface + size_maxcf * ( c + size_maxCorner * (zone) ) )]

#define soa_A_fp(a,icface,c,zone) soa_A_fp[ a + 3 * ( icface + size_maxcf * ( c + size_maxCorner * (zone) ) )]
#define soa_A_ez(a,icface,c,zone) soa_A_ez[ a + 3 * ( icface + size_maxcf * ( c + size_maxCorner * (zone) ) )]
#define soa_Connect(a,icface,c,zone) soa_Connect[ a + 3 * ( icface + size_maxcf * ( c + size_maxCorner * (zone) ) )]
#define soa_Connect_ro(a,icface,c,zone) soa_Connect_ro[ icface + size_maxcf * ( c + size_maxCorner * ( a + 3 * zone) ) ]
   

#define nextZ(a,b) nextZ[ (a) + nzones * (b) ]
#define next(a,b) next[ (a) + (ncornr+1)  * (b) ]


    //   for(int Angle=0;Angle<nAngle;Angle++)

    //int Angle = blockIdx.x;
    int Angle = AngleOrder[blockIdx.x]-1;

    double omega0, omega1, omega2;
    omega0 = soa_omega(0,Angle);
    omega1 = soa_omega(1,Angle);
    omega2 = soa_omega(2,Angle);
    
    omega_A_fp += blockIdx.x * nzones * size_maxcf * size_maxCorner;
    omega_A_ez += blockIdx.x * nzones * size_maxcf * size_maxCorner;

    int ndone = 0;
    int ndoneZ = 0;
    // hyperplane number p
    int p=0;

    while(ndoneZ < nzones)
    {
      //increment hyperplane
      p++;
    
      // get number of zones in this hyperplane
      int passZcnt = passZ[p] - passZ[p-1];

      // you can print hyperplanes for visualization
      //if( Angle == 0 && threadIdx.x==0) printf("%d \t %d\n",p,passZcnt);
      //for(int ii=threadIdx.x+blockIdx.y*blockDim.x;ii<passZcnt;ii+=blockDim.x*gridDim.y) 
      for(int ii=threadIdx.x;ii<passZcnt;ii+=blockDim.x) 
      {
	ndone = ( ndoneZ + ii ) * size_maxCorner;
    
	// get the zone (minus 1 so it is valid c index)
	int zone = nextZ(ndoneZ+ii,Angle) - 1;
 

	int nCorner   = soa_nCorner[zone];
	int nCFaces   = soa_nCFaces[zone];
	int c0        = soa_c0[zone] ;

	for(i=0;i<nCorner;i++)
	{
	  int ic      = next(ndone+i,Angle);
	  c       = ic - c0 - 1;

 
	  for(icface=0;icface<nCFaces;icface++)
	  {
	    omega_A_fp(icface,c,zone) =  omega0*soa_A_fp(0,icface,c,zone) + 
	      omega1*soa_A_fp(1,icface,c,zone) + 
	      omega2*soa_A_fp(2,icface,c,zone);
	    // could get rid of below if new order was used originally?
	    int icfp    = soa_Connect(0,icface,c,zone);
	    int ib      = soa_Connect(1,icface,c,zone);
	    int cez     = soa_Connect(2,icface,c,zone);
	    soa_Connect_ro(0,icface,c,zone) = icfp;
	    soa_Connect_ro(1,icface,c,zone) = ib  ;
	    soa_Connect_ro(2,icface,c,zone) = cez ;
	  
	  }


	  for(icface=0;icface<nCFaces;icface++)
	  {
    
	    omega_A_ez(icface,c,zone) = omega0*soa_A_ez(0,icface,c,zone) + omega1*soa_A_ez(1,icface,c,zone) + omega2*soa_A_ez(2,icface,c,zone) ;
	  }


	} // end corners 
      } // end zones in hplane

      ndoneZ += passZcnt;
      __syncthreads();

      //ndone = ndone + nCorner;

    }//end while

  }//end function


}







