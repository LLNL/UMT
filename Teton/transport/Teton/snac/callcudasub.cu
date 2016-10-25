/***********************************************************************
!                        Version 1:  04/16 SCR                         *
!                                                                      *
!   CUDA/C optimization of:                                            *
!                                                                      *
!   SNSWP3D  - This routine calculates angular fluxes for a single     *
!              direction for an upstream corner-balance spatial        *
!              discretization in 3D.                                   *
!                                                                      *
!   Input:                                                             *
!                                                                      *
!   Output:                                                            *
!                                                                      *
!   Questions:                                                         *
!              1. Should the hyperplane pass over zones or corners?    *
!                 A: Corners                                           *
!              2.                                                      *
!                                                                      *
!**********************************************************************/

/*
  TODO:
       explicit setting of #hyperplanes

       
*/

#include<stdio.h>
//#include<omp.h>
//#include"nvToolsExt.h"

#define GROUPSET 16
#define NUMFACES 3

#define NUMSM 90

#define fouralpha 1.82
#define fouralpha4 5.82

//#include "snutil.cu"
//#include "snsweep.cu"

#define Connect(a,b,c) Connect[ a + 3 * ( b + mC * c ) ]


extern "C" {

  void callCudaHostMalloc ();
//  void CC_sweep(
__global__ void GPU_sweep(
          int  size_maxCorner,
          int  size_maxcf,
          int  nAngle,
          int  nzones,
          int  ncornr,
          int  Groups,
          int  nbelem,
          int* Angle,
       double* soa_omega,
          int* nextZ,
          int* next,
          int* soa_nCorner,
          int* soa_nCFaces,
          int* soa_c0,
       double* soa_STotal,
       double* soa_STime,
       double* soa_SigtInv,
       double* soa_Volume,
       double* soa_Sigt,
       double* soa_A_fp,
       double* soa_A_ez,
          int* soa_Connect,
       double* psic,
       double* psib );

  /*
    Simple C file to handle calling CUDA.  Easily accessible from Fortran.  
    This removes any need to compile any other part of the program with nvcc.
    Just this file needs nvcc.  Everything else just uses regular intel compilers.
    Simple.
   */

  void callcudasub (
                     int *numzones, 
		     int *numgroups,
		     int *ncornr,
		     int *numAngles,
                     int *AngleOrder,
		     int *maxcorners, 
		     int *maxfaces, 
		     int *octant,  //=binRecv
		     int *NangBin, 
                     int *nbelem,
		     double *omega,
                     int    *nCorner,
                     int    *nCFaces,
                     int    *c0,
		     double *A_fp, 
		     double *A_ez, 
		     int    *Connect,
		     double *STotal,
		     double *STime,
		     double *Volume,
                     double *psic,
                     double *psib,
		     int *next, 
                     int *nextZ,
		     double *Sigt,
		     double *SigtInv
		     ) 
  {
    //dump data to check
    static int dump_cnt=0;
    int zone,ic;
    static    int* d_AngleOrder;
    static double* d_omega;
    static    int* d_nextZ;
    static    int* d_next;
    static    int* d_nCorner;
    static    int* d_nCFaces;
    static    int* d_c0;
    static double* d_STotal;
    static double* d_STime;
    static double* d_SigtInv;
    static double* d_Volume;
    static double* d_Sigt;
    static double* d_A_fp;
    static double* d_A_ez;
    static    int* d_Connect;
    static double* d_psic;
    static double* d_psib;

    int nZ = *numzones;
    int nA = *numAngles;
    int mC = *maxcorners;
    int mF = *maxfaces;
    int nG = *numgroups;
    int nC = *ncornr;
    int nBe = *nbelem;


    if ( dump_cnt < 5 )
    {
//      for(zone=0;zone<nZ;zone++)
//      {
//        for(ic=0;ic<nCorner[zone]; ic++)
//        {
//          printf(" zone,corner,connect3 = %d,%d,%d \n",zone,ic,Connect(2,ic,zone) );
//        }
//      }

      printf("max faces=%d\n",mF);

      if ( dump_cnt == 0 )
      {
        cudaMalloc(&d_AngleOrder,sizeof(int)*8*nA);
        cudaMalloc(&d_omega,sizeof(double)*3*nA);
        cudaMalloc(&d_nextZ,sizeof(int)*nZ*nA);
        cudaMalloc(&d_next,sizeof(int)*(nC+1)*nA);
        cudaMalloc(&d_nCorner,sizeof(int)*nZ);
        cudaMalloc(&d_nCFaces,sizeof(int)*nZ);
        cudaMalloc(&d_c0,sizeof(int)*nZ);
        cudaMalloc(&d_STotal,sizeof(double)*nZ*nG*mC);
        cudaMalloc(&d_STime,sizeof(double)*nZ*nA*nG*mC);
        cudaMalloc(&d_SigtInv,sizeof(double)*nZ*nG);
        cudaMalloc(&d_Volume,sizeof(double)*nZ*mC);
        cudaMalloc(&d_Sigt,sizeof(double)*nZ*nG);
        cudaMalloc(&d_A_fp,sizeof(double)*3*nZ*mC*mF);
        cudaMalloc(&d_A_ez,sizeof(double)*3*nZ*mC*mF);
        cudaMalloc(&d_Connect,sizeof(int)*3*nZ*mC*mF);
        cudaMalloc(&d_psic,sizeof(double)*nG*nC*nA);
        cudaMalloc(&d_psib,sizeof(double)*nG*nBe*nA);
      }



      cudaMemcpy(d_AngleOrder,AngleOrder,sizeof(int)*8*nA,cudaMemcpyHostToDevice);
      cudaMemcpy(d_omega,omega,sizeof(double)*3*nA,cudaMemcpyHostToDevice);
      cudaMemcpy(d_nextZ,nextZ,sizeof(int)*nZ*nA,cudaMemcpyHostToDevice);
      cudaMemcpy(d_next,next,sizeof(int)*(nC+1)*nA,cudaMemcpyHostToDevice);
      cudaMemcpy(d_nCorner,nCorner,sizeof(int)*nZ,cudaMemcpyHostToDevice);
      cudaMemcpy(d_nCFaces,nCFaces,sizeof(int)*nZ,cudaMemcpyHostToDevice);
      cudaMemcpy(d_c0,c0,sizeof(int)*nZ,cudaMemcpyHostToDevice);
      cudaMemcpy(d_STotal,STotal,sizeof(double)*nZ*nG*mC,cudaMemcpyHostToDevice);
      cudaMemcpy(d_STime,STime,sizeof(double)*nZ*nA*nG*mC,cudaMemcpyHostToDevice);
      cudaMemcpy(d_SigtInv,SigtInv,sizeof(double)*nZ*nG,cudaMemcpyHostToDevice);
      cudaMemcpy(d_Volume,Volume,sizeof(double)*nZ*mC,cudaMemcpyHostToDevice);
      cudaMemcpy(d_Sigt,Sigt,sizeof(double)*nZ*nG,cudaMemcpyHostToDevice);
      cudaMemcpy(d_A_fp,A_fp,sizeof(double)*3*nZ*mC*mF,cudaMemcpyHostToDevice);
      cudaMemcpy(d_A_ez,A_ez,sizeof(double)*3*nZ*mC*mF,cudaMemcpyHostToDevice);
      cudaMemcpy(d_Connect,Connect,sizeof(int)*3*nZ*mC*mF,cudaMemcpyHostToDevice);
      cudaMemcpy(d_psib,psib,sizeof(double)*nG*nBe*nA,cudaMemcpyHostToDevice);
      cudaMemcpy(d_psic,psic,sizeof(double)*nG*nC*nA,cudaMemcpyHostToDevice);

//
//      
      GPU_sweep<<<nA,32>>>(
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
                       d_STime,                     
                       d_SigtInv,                   
                       d_Volume,                    
                       d_Sigt,                      
                       d_A_fp,                      
                       d_A_ez,                      
                       d_Connect,                   
                       d_psic,                  
                       d_psib        
                          );

      cudaMemcpy(psic,d_psic,sizeof(double)*nG*nC*nA,cudaMemcpyDeviceToHost);



//      CC_sweep(
//                  *maxcorners,      
//                  mF,       //*maxfaces,      
//                  *numAngles,      
//                  *numzones,      
//                  nC,      
//                  *numgroups,      
//                  nBe,      
//                  AngleOrder,      
//                  omega,      
//                  nextZ,      
//                  next,      
//                  nCorner,      
//                  nCFaces,      
//                  c0,      
//                  STotal,      
//                  STime,      
//                  SigtInv,      
//                  Volume,      
//                  Sigt,      
//                  A_fp,      
//                  A_ez,      
//                  Connect,      
//                  psic,      
//                  psib               
//             );
                  

      dump_cnt++;
    }
  } //callcuasub

}  // extern C

