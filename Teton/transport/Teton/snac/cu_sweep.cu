#include <stdio.h>

#define WARP_SIZE 32
extern "C"
{
  __global__ void cuda_test(
			    int anglebatch,
			    int numgroups,
			    int ncornr,
			    int size_maxCorner,
			    double* psicbatch)
  {

    int Groups=192;

#define psicbatch(ig,b,c) psicbatch[(ig) + Groups * ((b) + size_maxCorner *(c) )]    
    
    printf("anglebatch = %d\n", anglebatch);
    //printf("numgroups = %d\n", numgroups);
    //printf("ncornr = %d\n", ncornr);
    //printf("size_maxCorner = %d\n",size_maxCorner);
    
    for (int angle=blockIdx.x; angle<anglebatch; angle+=gridDim.x) {
    
      for(int c=0; c<size_maxCorner; c++) {
	{
      
	  for (int g=threadIdx.x; g<numgroups; g+= blockDim.x) {
	    psicbatch(g,c,angle) = 1;
	  }
	}
      }
    }
  }	     

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
       double* soa_STime,
       double* soa_SigtInv,
       double* soa_Volume,
       double* soa_Sigt,
       double* soa_A_fp,
       double* soa_A_ez,
          int* soa_Connect,
       double* psic,
       double* psib )
{

//   double omega[3];
   int c,ig,i,icface,ifp,cez,k,ii;
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

   double psi_opp,tpsic;
   double omega0, omega1, omega2;
   __shared__ double *Q, *src, *volume, *coefpsic, *afpm, *psifp;
   __shared__ int *ez_exit;
   __shared__ double sm_agg[625];
   
   int offset = 0;

   Q = &(sm_agg[0]);
   offset += size_maxCorner * WARP_SIZE;

   src = &(sm_agg[offset]);
   offset += size_maxCorner * WARP_SIZE;

   volume = &(sm_agg[offset]);
   offset += size_maxCorner;

   coefpsic = &(sm_agg[offset]);
   offset += size_maxcf;

   afpm = &(sm_agg[offset]);
   offset += size_maxcf;

   psifp = &(sm_agg[offset]);
   offset += size_maxcf * WARP_SIZE;

   //note ez_exit has integer type
   ez_exit = (int*) &(sm_agg[offset]);
   
   const double fouralpha = 1.82;
//   const double fouralpha4 = 5.82;
   
   #define soa_omega(a,b) soa_omega[a + 3 * b]

//   #define tpsic(ig,c) tpsic[ (ig) + Groups * (c)]
   #define EB_ListExit(a,ia) EB_ListExit[ a + 2 * (ia) ]
   #define soa_A_fp(a,icface,c,zone) soa_A_fp[ a + 3 * ( icface + size_maxcf * ( c + size_maxCorner * (zone) ) )]
   #define soa_A_ez(a,icface,c,zone) soa_A_ez[ a + 3 * ( icface + size_maxcf * ( c + size_maxCorner * (zone) ) )]
   #define soa_Connect(a,icface,c,zone) soa_Connect[ a + 3 * ( icface + size_maxcf * ( c + size_maxCorner * (zone) ) )]
   
   #define psifp(ig,jf) psifp[(ig) + Groups * (jf)]
   #define psib(ig,b,c) psib[(ig) + Groups * ((b) + nbelem * (c) )]
   #define psic(ig,b,c) psic[(ig) + Groups * ((b) + size_maxCorner *(c) )]
   
   #define Q(ig,c) Q[(ig) + Groups * (c)]
   #define src(ig,c) src[(ig) + Groups * (c)]
//  #define SigtVol(ig,c) SigtVol[(ig) + Groups * (c)]
   #define soa_Sigt(ig,zone) soa_Sigt[(ig) + Groups * (zone)]
   #define soa_Volume(c,zone) soa_Volume[(ig) + Groups * (zone)]
   #define soa_SigtInv(ig,zone) soa_SigtInv[(ig) + Groups * (zone)]
   #define soa_STotal(ig,c,zone) soa_STotal[ig + Groups * ( c + size_maxCorner * (zone) )]
   #define soa_STime(ig,c,Angle,zone) soa_STime[ig + Groups * ( c + size_maxCorner * ( Angle + nAngle * (zone) ) )]
   #define nextZ(a,b) nextZ[ (a) + nzones * (b) ]
   #define next(a,b) next[ (a) + (ncornr+1)  * (b) ]





//   for(int Angle=0;Angle<nAngle;Angle++)

   int Angle = blockIdx.y;
   ig = threadIdx.x;

//   if(ig==0) printf("my offset=%d\n",offset);
//   if(ig==0)
//   {
//     printf("psic=%x\n",psic);
//     printf("nextZ=%x\n",psic);
//     printf("next=%x\n",psic);
//     printf("psib=%x\n",psic);
//   }

   {

   omega0 = soa_omega(0,Angle);
   omega1 = soa_omega(1,Angle);
   omega2 = soa_omega(2,Angle);

   int ndone = 0;


   for(ii=0;ii<nzones;ii++)
   {
 
     int zone = nextZ(ii,Angle) - 1;


     int nCorner   = soa_nCorner[zone];
     int nCFaces   = soa_nCFaces[zone];
     int c0        = soa_c0[zone] ;

     double Sigt = soa_Sigt(ig,zone);

     for(c=0;c<nCorner;c++)
     {
         double source = soa_STotal(ig,c,zone) + soa_STime(ig,c,Angle,zone);
         Q(ig,c)       = soa_SigtInv(ig,zone)*source ;
         src(ig,c)     = soa_Volume(c,zone) *source;
         //SigtVol(ig,c) = soa_Sigt(ig,zone)*soa_Volume(c,zone);
         volume[c] = soa_Volume(c,zone);
     }


     for(i=0;i<nCorner;i++)
     {

       int ic      = next(ndone+i,Angle);
       c       = ic - c0 - 1;

       sumArea = 0.0;
 
       for(icface=0;icface<nCFaces;icface++)
       {
         afpm[icface] = omega0*soa_A_fp(0,icface,c,zone) + 
                        omega1*soa_A_fp(1,icface,c,zone) + 
                        omega2*soa_A_fp(2,icface,c,zone);

         int icfp    = soa_Connect(1,icface,c,zone) - 1;
         int ib      = soa_Connect(2,icface,c,zone) - 1;
                                                                                                   
         if ( afpm[icface] >= 0.0 )
         { 
           sumArea = sumArea + afpm[icface];
         }
         else
         {
           if (icfp == 0)
           {
//             psifp(ig,icface) = psib(ig,ib,Angle);
             r_psifp = psib(ig,ib,Angle);
//             r_psifp = 0.3;
           }
           else
           {
//             psifp(ig,icface) = psic(ig,icfp,Angle);
             r_psifp = psic(ig,icfp,Angle);
//             r_psifp = 0.7;
           }

           src(ig,c)  -= afpm[icface]*r_psifp;
           psifp(ig,icface) = r_psifp;
         }
       }


       int nxez = 0;

       for(icface=0;icface<nCFaces;icface++)
       {

         double aez = omega0*soa_A_ez(0,icface,c,zone) + omega1*soa_A_ez(1,icface,c,zone) + omega2*soa_A_ez(2,icface,c,zone) ;

         if (aez > 0.0 )
         {

           sumArea        = sumArea + aez;
           area_opp       = .0;
           cez            = soa_Connect(2,icface,c,zone) - 1;
           ez_exit[nxez]  = cez;
           coefpsic[nxez] = aez;
           nxez           = nxez + 1;

           if (nCFaces == 3)
           {

             ifp = (icface+1)%nCFaces;

             if ( afpm[ifp] < 0.0 )
             { 
               area_opp   = -afpm[ifp];
               psi_opp =  psifp(ig,ifp);
             }
           }
           else
           {

             ifp        = icface;
             area_opp   = 0.0;
             psi_opp = 0.0;

             for(k=0;k<nCFaces-2;k++)
             {
               ifp = ifp%nCFaces;
               if ( afpm[ifp] < 0.0 )
               {
                 area_opp   = area_opp   - afpm[ifp];
                 psi_opp = psi_opp - afpm[ifp]*psifp(ig,ifp);
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
               double sez          = gtau*sigv*( psi_opp - Q(ig,c) ) +   0.5*aez*(1.0 - gtau)*( Q(ig,c) - Q(ig,cez) );

               src(ig,c)    = src(ig,c)   + sez;
               src(ig,cez)  = src(ig,cez) - sez;

             }

           }
           else
           {
               double sez          = 0.5*aez*( Q(ig,c) - Q(ig,cez) );
               src(ig,c)    = src(ig,c)   + sez;
               src(ig,cez)  = src(ig,cez) - sez;

           } 
         }
       }

//       printf("ckim angle,zone,corner,aez_cnt %d,%d,%d,%d\n",Angle,zone,c,aez_cnt);


       tpsic = src(ig,c)/(sumArea + Sigt*volume[c]);


       for(icface=0;icface<nxez;icface++)
       {
         int cez   = ez_exit[icface];
         src(ig,cez) = src(ig,cez) + coefpsic[icface]*tpsic;
       }

       //hope that ther is no self referencing
       psic(ig,c0+c,Angle) = tpsic;
     } 

     ndone = ndone + nCorner;
   } 
   }


//   ExitBdy => getExitList(QuadSet, Angle)

//   for(i=0;i<EB_nExit;i++)
//   {
 //    int ib = EB_ListExit(1,i);
//     int ic = EB_ListExit(2,i);
//     for(ig=0;ig<Groups;ig++)
//       psib(ig,ib) = psic(ig,ic);
//   }
}



void CC_sweep(
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
       double* psib )
{

   double omega[3];
   int c,ig,i,icface,ifp,cez,k,ii;
   double Q[Groups * size_maxCorner];
   double src[Groups * size_maxCorner];
   double SigtVol[Groups * size_maxCorner];
   double afpm[size_maxcf];
   double psifp[Groups * size_maxcf];
   int    ez_exit[size_maxcf];
   double coefpsic[size_maxcf];
   double tpsic[Groups * size_maxCorner];
   double psi_opp[Groups];
   double area_opp,area_inv,sumArea;
   
   const double fouralpha = 1.82;
//   const double fouralpha4 = 5.82;
   
   #define soa_omega(a,b) soa_omega[a + 3 * b]

   #define tpsic(ig,c) tpsic[ (ig) + Groups * (c)]
   #define EB_ListExit(a,ia) EB_ListExit[ a + 2 * (ia) ]
   #define soa_A_fp(a,icface,c,zone) soa_A_fp[ a + 3 * ( icface + size_maxcf * ( c + size_maxCorner * (zone) ) )]
   #define soa_A_ez(a,icface,c,zone) soa_A_ez[ a + 3 * ( icface + size_maxcf * ( c + size_maxCorner * (zone) ) )]
   #define soa_Connect(a,icface,c,zone) soa_Connect[ a + 3 * ( icface + size_maxcf * ( c + size_maxCorner * (zone) ) )]
   
   #define psifp(ig,jf) psifp[(ig) + Groups * (jf)]
   #define psib(ig,b,c) psib[(ig) + Groups * ((b) + nbelem * (c) )]
   #define psic(ig,b,c) psic[(ig) + Groups * ((b) + size_maxCorner *(c) )]
   
   #define Q(ig,c) Q[(ig) + Groups * (c)]
   #define src(ig,c) src[(ig) + Groups * (c)]
   #define SigtVol(ig,c) SigtVol[(ig) + Groups * (c)]
   #define soa_Sigt(ig,zone) soa_Sigt[(ig) + Groups * (zone)]
   #define soa_Volume(c,zone) soa_Volume[(ig) + Groups * (zone)]
   #define soa_SigtInv(ig,zone) soa_SigtInv[(ig) + Groups * (zone)]
   #define soa_STotal(ig,c,zone) soa_STotal[ig + Groups * ( c + size_maxCorner * (zone) )]
   #define soa_STime(ig,c,Angle,zone) soa_STime[ig + Groups * ( c + size_maxCorner * ( Angle + nAngle * (zone) ) )]
   #define nextZ(a,b) nextZ[ (a) + nzones * (b) ]
   #define next(a,b) next[ (a) + (ncornr+1)  * (b) ]




   for(int Angle=0;Angle<nAngle;Angle++)
   {

   omega[0] = soa_omega(0,Angle);
   omega[1] = soa_omega(1,Angle);
   omega[2] = soa_omega(2,Angle);


   int ndone = 0;

   for(ii=0;ii<nzones;ii++)
   {
 
     int zone = nextZ(ii,Angle) - 1;

     int nCorner   = soa_nCorner[zone];
     int nCFaces   = soa_nCFaces[zone];
     int c0        = soa_c0[zone] ;

     for(c=0;c<nCorner;c++)
     {
       for(ig=0;ig<Groups;ig++)
       {
         double source = soa_STotal(ig,c,zone) + soa_STime(ig,c,Angle,zone);
         Q(ig,c)       = soa_SigtInv(ig,zone)*source ;
         src(ig,c)     = soa_Volume(c,zone) *source;
         SigtVol(ig,c) = soa_Sigt(ig,zone)*soa_Volume(c,zone);
       }
     }

     for(i=0;i<nCorner;i++)
     {

       int ic      = next(ndone+i,Angle);
       c       = ic - c0 - 1;

       sumArea = 0.0;
 
       for(icface=0;icface<nCFaces;icface++)
       {
         afpm[icface] = omega[0]*soa_A_fp(0,icface,c,zone) + 
                        omega[1]*soa_A_fp(1,icface,c,zone) + 
                        omega[2]*soa_A_fp(2,icface,c,zone);

         int icfp    = soa_Connect(1,icface,c,zone) - 1;
         int ib      = soa_Connect(2,icface,c,zone) - 1;
                                                                                                   
         if ( afpm[icface] >= 0.0 )
         { 
           sumArea = sumArea + afpm[icface];
         }
         else
         {
           if (icfp == 0)
           {
             for(ig=0;ig<Groups;ig++) psifp(ig,icface) = psib(ig,ib,Angle);
           }
           else
           {
             for(ig=0;ig<Groups;ig++) psifp(ig,icface) = psic(ig,icfp,Angle);
           }

           for(ig=0;ig<Groups;ig++) src(ig,c)  = src(ig,c) - afpm[icface]*psifp(ig,icface);
         }
       }


       int nxez = 0;

       for(icface=0;icface<nCFaces;icface++)
       {

         double aez = omega[0]*soa_A_ez(0,icface,c,zone) + omega[1]*soa_A_ez(1,icface,c,zone) + omega[2]*soa_A_ez(2,icface,c,zone) ;

         if (aez > 0.0 )
         {

           sumArea        = sumArea + aez;
           area_opp       = .0;
           cez            = soa_Connect(2,icface,c,zone) - 1;
           ez_exit[nxez]  = cez;
           coefpsic[nxez] = aez;
           nxez           = nxez + 1;

           if (nCFaces == 3)
           {

             ifp = icface%nCFaces;

             if ( afpm[ifp] < 0.0 )
             { 
               area_opp   = -afpm[ifp];
               for(ig=0;ig<Groups;ig++) psi_opp[ig] =  psifp(ig,ifp);
             }
           }
           else
           {

             ifp        = icface;
             area_opp   = 0.0;
             for(ig=0;ig<Groups;ig++) psi_opp[ig] = 0.0;

             for(k=0;k<nCFaces-2;k++)
             {
               ifp = ifp%nCFaces;
               if ( afpm[ifp] < 0.0 )
               {
                 area_opp   = area_opp   - afpm[ifp];
                 for(ig=0;ig<Groups;ig++) psi_opp[ig] = psi_opp[ig] - afpm[ifp]*psifp(ig,ifp);
               }
             }

             area_inv = 1.0/area_opp;

             for(ig=0;ig<Groups;ig++)
               psi_opp[ig] = psi_opp[ig]*area_inv;

           }

           if (area_opp > 0.0) {

             double aez2 = aez*aez;

             for(ig=0;ig<Groups;ig++)
             {
  
               double sigv         = SigtVol(ig,c);
               double sigv2        = sigv*sigv;
 
               double gnum         = aez2*( fouralpha*sigv2 +       aez*(4.0*sigv + 3.0*aez) );

               double gtau         = gnum/( gnum + 4.0*sigv2*sigv2 + aez*sigv*(6.0*sigv2 + 2.0*aez*(2.0*sigv + aez)) ) ;

               double sez          = gtau*sigv*( psi_opp[ig] - Q(ig,c) ) +   0.5*aez*(1.0 - gtau)*( Q(ig,c) - Q(ig,cez) );

               src(ig,c)    = src(ig,c)   + sez;
               src(ig,cez)  = src(ig,cez) - sez;

             }

           }
           else
           {

             for(ig=0;ig<Groups;ig++)
             {
               double sez          = 0.5*aez*( Q(ig,c) - Q(ig,cez) );
               src(ig,c)    = src(ig,c)   + sez;
               src(ig,cez)  = src(ig,cez) - sez;
             }

           } 
         }
       }

//       printf("ckim angle,zone,corner,aez_cnt %d,%d,%d,%d\n",Angle,zone,c,aez_cnt);


       for(ig=0;ig<Groups;ig++)
         tpsic(ig,c) = src(ig,c)/(sumArea + SigtVol(ig,c));


       for(icface=0;icface<nxez;icface++)
       {
         int cez   = ez_exit[icface];
         for(ig=0;ig<Groups;ig++)
           src(ig,cez) = src(ig,cez) + coefpsic[icface]*tpsic(ig,c);
       }

     } 

     ndone = ndone + nCorner;


     for(c=0;c<nCorner;c++)
     {
       for(ig=0;ig<Groups;ig++)
         psic(ig,c0+c,Angle) = tpsic(ig,c);
     }
                                                                                                
   } 
   ndone++; 
   }


//   ExitBdy => getExitList(QuadSet, Angle)

//   for(i=0;i<EB_nExit;i++)
//   {
 //    int ib = EB_ListExit(1,i);
//     int ic = EB_ListExit(2,i);
//     for(ig=0;ig<Groups;ig++)
//       psib(ig,ib) = psic(ig,ic);
//   }
}
}
