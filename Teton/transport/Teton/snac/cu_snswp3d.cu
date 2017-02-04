#include<stdio.h>
#include<iostream>

#define WARP_SIZE 32

#define GROUPSET 16
#define NUMFACES 3

#define fouralpha 1.82
#define fouralpha4 5.82


#define Connect(a,b,c) Connect[ a + 3 * ( b + mC * c ) ]

extern "C" {


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
       double* STime,
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
	  int* Connect_ro,
          int* passZ,
	  bool calcSTime,
	  double tau
 );

__global__ void GPU_fp_ez_hplane(
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
       double* soa_A_fp,
       double* soa_A_ez,
       double* omega_A_fp,
       double* omega_A_ez,
          int* soa_Connect,
          int* soa_Connect_reorder,
	  int* passZ
 );


  void snswp3d_c (
		  int *anglebatch,
		  int *numzones, 
		  int *numgroups,
		  int *ncornr,
		  int *numAngles,
		  int *d_AngleOrder,
		  int *maxcorners, 
		  int *maxfaces, 
		  int *octant,  //=binRecv
		  int *NangBin, 
		  int *nbelem,
		  double *d_omega,
		  int    *d_nCorner,
		  int    *d_nCFaces,
		  int    *d_c0,
		  double *d_A_fp,
		  double *d_omega_A_fp,
		  double *d_A_ez, 
		  double* d_omega_A_ez,
		  int    *d_Connect,
		  int* d_Connect_reorder,
		  double *d_STotal,
		  double *d_STimeBatch,
		  double *d_STime,
		  double *d_Volume,
		  double *d_psic,
		  double *d_psib,
		  int *d_next, 
		  int *d_nextZ,
		  double *d_Sigt,
		  double *d_SigtInv,
		  int *d_passZ,
		  bool *calcSTime,
		  double *tau,
		  cudaStream_t streamid
		  ) 
  {
	  static int dump_cnt=0;
    //int zone,ic;
    //static double* d_omega_A_fp;
    //static double* d_omega_A_ez;
    //static    int* d_Connect_reorder;

    int nZ = *numzones;
    int nA = *numAngles;
    // will need this for large problems
    int nAbatch = *anglebatch;
    int mC = *maxcorners;
    int mF = *maxfaces;
    int nG = *numgroups;
    int nC = *ncornr;
    int nBe = *nbelem;

    {

      //printf("max faces=%d\n",mF);

      printf("*tau = %.17g\n",*tau);

      // first time being called, allocate some host and device arrays
      if ( dump_cnt == 0 )
      {

	// I guess these are synchronous now... Will need to make data members of ZData later.

	     // create device versions
	//     cudaMalloc(&d_omega_A_fp,sizeof(double)*nZ*mC*mF*nA);
	//   cudaMalloc(&d_omega_A_ez,sizeof(double)*nZ*mC*mF*nA);
	//   printf("d_omega_A_fp allocated size = %d\n",nZ*mC*mF*nA*8);

        //cudaMalloc(&d_Connect_reorder,sizeof(int)*3*nZ*mC*mF);
        
      }
      

      if( *octant == 1) {
	     // This does all the angles. Redundant when angles are done in batches.
	     // Could async copy psic or psib while doing all angles once at beginning.
	     // Actually batched works too, since this does not depend on psic or psib.
	//GPU_fp_ez<<<nA/32,32,0,streamid>>>(
	GPU_fp_ez_hplane<<<nA,32,0,streamid>>>(
				mC,                 
				mF,       //                 
				nA,                  
				nZ,                   
				nC,                          
				nG,                  
				nBe,                         
				d_AngleOrder,                
				d_omega,                     
				d_nextZ,                     
				d_next,                      
				d_nCorner,                   
				d_nCFaces,                   
				d_c0,                        
				d_A_fp,                      
				d_A_ez,                      
				d_omega_A_fp,
				d_omega_A_ez,
				d_Connect,
				d_Connect_reorder,
				d_passZ
				);
      }


      //cudaStreamSynchronize(streamid );
      
      // synchronize to be sure all streams have fully transferred psi to device. 
      //cudaDeviceSynchronize(); // remove this later.



      int nGG = ceil(nG / 32.0);
      //printf("nGG=%d\n",nGG);
      if (nG%32 != 0) {printf("current version must use groups of multiple of 32!!! sorry \n"); exit(0);}

      // shared memory needs are (8+3+3*blockDim.x+3)*blockDim.y;

      //GPU_sweep<<<dim3(*anglebatch,nGG,2),dim3(32,16,1),(8+3+3*32+3)*16*sizeof(double),streamid>>>(
      GPU_sweep<<<dim3(*anglebatch,nGG,1),dim3(32,32,1),(8+3+3*32+3)*32*sizeof(double),streamid>>>(
                       mC,                 
                       mF,       //                 
                       nA,                  
                       nZ,                   
                       nC,                          
                       nG,                  
                       nBe,                         
                       d_AngleOrder,                
                       d_omega,                     
                       d_nextZ,                     
                       d_next,                      
                       d_nCorner,                   
                       d_nCFaces,                   
                       d_c0,                        
                       d_STotal,                    
		       d_STimeBatch,
                       d_STime,                     
                       d_SigtInv,                   
                       d_Volume,                    
                       d_Sigt,                      
                       d_A_fp,                      
                       d_A_ez,                      
                       d_Connect,                   
                       d_psic,                  
                       d_psib,
                       d_omega_A_fp,
                       d_omega_A_ez,
                       d_Connect_reorder,
                       d_passZ,
		       *calcSTime,
		       *tau
                          );

      //printf("Completed a batch sweep\n");


      dump_cnt++;
      //std::cout<<"dump_cnt="<<dump_cnt<<std::endl;
    }
  } 

} // extern "C"
