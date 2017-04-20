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


  void fp_ez_c (
		  int *anglebatch,
		  int *numzones, 
		  int *numgroups,
		  int *ncornr,
		  int *numAngles,
		  int *d_AngleOrder,
		  int *maxcorners, 
		  int *maxfaces, 
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
		  int *d_next, 
		  int *d_nextZ,
		  int *d_passZ,
		  cudaStream_t streamid
		  ) 
  {

    int nZ = *numzones;
    int nA = *numAngles;
    // will need this for large problems
    int nAbatch = *anglebatch;
    int mC = *maxcorners;
    int mF = *maxfaces;
    int nG = *numgroups;
    int nC = *ncornr;
    int nBe = *nbelem;


	GPU_fp_ez_hplane<<<dim3(nAbatch,1,1),128,0,streamid>>>(
				mC,                 
				mF,       //                 
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


      int groupsize=32;

      int nGG = ceil(nG / groupsize);
      //printf("nGG=%d\n",nGG);
      if (nG%groupsize != 0) {printf("current version must use groups of multiple of %d!!! sorry \n",groupsize); exit(0);}

      // shared memory needs are (8+3+3*blockDim.x+3)*blockDim.y;

      //GPU_sweep<<<dim3(*anglebatch,nGG,2),dim3(32,16,1),(8+3+3*32+3)*16*sizeof(double),streamid>>>(
      GPU_sweep<<<dim3(*anglebatch,nGG,1),dim3(groupsize,32,1),(8+3+3*groupsize+3)*32*sizeof(double),streamid>>>(
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
