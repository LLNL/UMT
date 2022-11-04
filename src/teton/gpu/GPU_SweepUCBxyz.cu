#include <stdio.h>
#include <string.h>

#include "nvToolsExt.h"
#include <cublas_v2.h>

#include <stdlib.h>

// GPU UCB xyz Sweep
//
// Parallelization: 1 direction * 16 groups per block (SM) [directions*groups/16]
// Parallelization: 16 groups * nzones per hyperplane
// Parallelization: total = directions * groups * zones/per/hyperplane
//
//

// Macro to catch CUDA errors in CUDA runtime calls
#define CUDA_SAFE_CALL(call)                                                                                       \
   do                                                                                                              \
   {                                                                                                               \
      cudaError_t err = call;                                                                                      \
      if (cudaSuccess != err)                                                                                      \
      {                                                                                                            \
         fprintf(stderr, "Cuda error in file '%s' in line %i : %s.", __FILE__, __LINE__, cudaGetErrorString(err)); \
         exit(EXIT_FAILURE);                                                                                       \
      }                                                                                                            \
   } while (0)

#define FOURALPHA 1.82

// this should work for 32, 16 or 4 - haven't tested others.  Maybe the only other sensible number would be 8?
#define GROUPS_IN_BLOCK 4
#define THREAD_BLOCK_SIZE 128
#define MAX_CUDA_STREAMS 80

// measured to be about 570.  1024 seems like a good place to start
// RCC - TODO: find the Teton number for this (limit of cycleList, cyclePsi)
// RCC - MAX_TOTAL_CYCLES used to determine size of cycleList and cyclePsi.
//       This size overshoots actual array dimensions.
//       Need a more dynamic allocation method.
#define MAX_TOTAL_CYCLES 250000
//#define MAX_TOTAL_CYCLES 8000
#define MAX_PLANES 8192
//#define MAX_PLANES 4096

extern "C" {

__global__ void initFromCycleListKernel(const int numCycles,
                                        const int cycleOffSet,
                                        int *cycleList,
                                        double *cyclePsi,
                                        double *Psi,
                                        const int maxCorner,
                                        const int Groups)
{
   int tid = blockIdx.x * blockDim.x + threadIdx.x; // thread ID
   int myGroup = tid % Groups;
   int firstCycle = tid / Groups; // leftover
   int cycleStep = gridDim.x * blockDim.x / Groups;
   if (tid < (cycleStep * Groups))
   {
      for (int cycle = firstCycle; cycle < numCycles; cycle += cycleStep)
      {
         int mCycle = cycleOffSet + cycle;
         int c = cycleList[mCycle] - 1;
         Psi[Groups * c + myGroup] = cyclePsi[Groups * mCycle + myGroup];
      }
   }

   return;
}

__global__ void updateCycleListKernel(const int numCycles,
                                      const int cycleOffSet,
                                      int *cycleList,
                                      double *cyclePsi,
                                      double *Psi,
                                      const int maxCorner,
                                      const int Groups)
{
   int tid = blockIdx.x * blockDim.x + threadIdx.x;
   int myGroup = tid % Groups;
   int firstCycle = tid / Groups; // leftover
   int cycleStep = gridDim.x * blockDim.x / Groups;
   if (tid < (cycleStep * Groups))
   {
      for (int cycle = firstCycle; cycle < numCycles; cycle += cycleStep)
      {
         int mCycle = cycleOffSet + cycle;
         int c = cycleList[mCycle] - 1;
         cyclePsi[Groups * mCycle + myGroup] = Psi[Groups * c + myGroup];
      }
   }

   return;
}

// Leaving older implementations of helper kernels here for ideas on optimization.
//  __global__ void initFromCycleListKernelOLD ( const int numCycles, const int cycleOffSet, int *cycleList, double *cyclePsi, double *Psi, const int maxCorner, const int Groups)
//  {
//    int tid = blockIdx.x * blockDim.x + threadIdx.x;                      // thread ID
//    int icycle = tid / ( Groups * maxCorner );                            // 'initial' cycle for this thread
//    int cycleStep = (gridDim.x*blockDim.x)/( Groups*maxCorner);           // cycle step : how many cycles the entire grid of threads can handle
//    int corner = (tid%(Groups*maxCorner)) / Groups;                       // corner for this thread
//    int group = tid % Groups;                                             // group for this thread
//    int c0, nCorner, cycle;
//
//#pragma unroll
//    for ( cycle=icycle; cycle<numCycles; cycle+=cycleStep ) {             // loop over all cycles (threads get corners and groups)
//      c0 = cycleList[(cycle)*2];                                            //
//      nCorner = cycleList[(cycle)*2+1];
//      if ( corner < nCorner ) {
//        Psi[(c0+corner)* Groups + group ] = cyclePsi[(cycle) * maxCorner * Groups + corner * Groups + group ];
//      }
//    }
//
//  }
//
//  __global__ void updateCycleListKernelOLD ( const int numCycles, const int cycleOffSet, int *cycleList, double *cyclePsi, double *Psi, const int maxCorner, const int Groups)
//  {
//    int tid = blockIdx.x * blockDim.x + threadIdx.x;
//    int icycle = tid / ( Groups * maxCorner );
//    int cycleStep = (gridDim.x*blockDim.x)/( Groups*maxCorner);
//    int corner = (tid%(Groups*maxCorner)) / Groups;
//    int group = tid % Groups;
//    int c0, nCorner, cycle;
//
//#pragma unroll
//    for ( cycle=icycle; cycle<numCycles; cycle+=cycleStep ) {
//      c0 = cycleList[(cycle)*2];
//      nCorner = cycleList[(cycle)*2+1];
//      if ( corner < nCorner ) {
//        //Psi[(c0+corner)* Groups + group ] = cyclePsi[cycle * maxCorner * Groups + corner * Groups + group ];
//        cyclePsi[(cycle) * maxCorner * Groups + corner * Groups + group ] = Psi[(c0+corner)* Groups + group ];
//      }
//    }
//
//  }

// RCC - Steve believed that A_fp/A_ez did not change because mesh is fixed in this cycle, a product of fixed geometry
__global__ void compAfpKernel(
    double *omega, double *A_fp, double *d_afp, double *A_ez, double *d_aez, unsigned int nVals, unsigned int ndim)
{
   unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x; // this is the # of the aez entry being computed
   unsigned int cstride = gridDim.x * blockDim.x;

   for (unsigned int cornerCface = tid; cornerCface < nVals; cornerCface += cstride)
   {
      double afp0 = A_fp[cornerCface * ndim];
      double afp1 = A_fp[cornerCface * ndim + 1];
      double afp2 = A_fp[cornerCface * ndim + 2];

      double aez0 = A_ez[cornerCface * ndim];
      double aez1 = A_ez[cornerCface * ndim + 1];
      double aez2 = A_ez[cornerCface * ndim + 2];

      d_afp[cornerCface] = omega[0] * afp0 + omega[1] * afp1 + omega[2] * afp2;
      d_aez[cornerCface] = omega[0] * aez0 + omega[1] * aez1 + omega[2] * aez2;
   }
}

__global__ void CompSTotalTauPsi(double *d_STotal, double *d_Psi, double tau, int num)
{
   unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
   unsigned int stride = gridDim.x * blockDim.x;
   for (unsigned int i = tid; i < num; i += stride)
   {
      d_STotal[i] += tau * d_Psi[i];
   }
}

__global__ void SweepUCBxyzKernel(int Angle,
                                  int nHyperPlanes,
                                  int *nZonesInPlane,
                                  int *nextZ,
                                  int *nextC,
                                  double *STotal,
                                  double tau,
                                  const int Groups,
                                  double *Volume,
                                  double *Sigt,
                                  int *nCFacesArray,
                                  int ndim,
                                  int maxcf,
                                  int ncorner,
                                  double *omega,
                                  int *cFP,
                                  int *bdy_exit,
                                  double *Psi1,
                                  int *cEZ,
                                  double quadwt,
                                  double *Phi,
                                  double *PsiB,
                                  int maxCorner,
                                  double *d_afp,
                                  double *d_aez,
                                  int *d_GnumCorner,
                                  int *d_GcOffSet)
{
   extern __shared__ char smem[];

   //int tid = blockIdx.x * blockDim.x + threadIdx.x;                        // thread index within the grid

   // each block deals with GROUPS_IN_BLOCK groups.
   // so, I can't have more than Groups/groups16 blocks in my grid. (for now ... until we add angles)
   // so, the group handled by this thread, is threadIdx.x%GROUPS_IN_BLOCK + blockIdx.x * GROUPS_IN_BLOCK
   int g = threadIdx.x % GROUPS_IN_BLOCK + blockIdx.x * GROUPS_IN_BLOCK; // group number that the threads handles

   int zoneBatch = blockDim.x / GROUPS_IN_BLOCK; // Number of zones that can be updated simultaneously
   int ii = threadIdx.x / GROUPS_IN_BLOCK;       // index, within the chunk, of the zone handled by this thread

   int zone0;
   int zone;
   int nCorner;
   int c0;
   int nCFaces = 0;
   //int hyp_c0;
   int nxBdy;
   int cez;
   double source;
   double area_opp;
   double psi_opp;

   double *sQ = (double *) smem;
   double *Q = sQ + threadIdx.x * maxCorner;
   //double Q[MAXC];

   double *sSrc = sQ + blockDim.x * maxCorner;
   double *src = sSrc + threadIdx.x * maxCorner;
   //double src[MAXC];

   // TODO:  Replace sSigtVol with just sSigt - the replace all occurances of SigtVol with sSigt * sVol.
   double *sSigtVol = sSrc + blockDim.x * maxCorner;
   double *SigtVol = sSigtVol + threadIdx.x * maxCorner;
   //double SigtVol[MAXC];

   double *sSumArea = sSigtVol + blockDim.x * maxCorner;
   double *sumArea = sSumArea + threadIdx.x * maxCorner;
   //double sumArea[MAXC];

   double *sAfp = sSumArea + blockDim.x * maxCorner;
   double *afp = sAfp + threadIdx.x * maxcf;
   //double afp[MAX_CF];

   double *sPsifp = sAfp + blockDim.x * maxcf;
   double *psifp = sPsifp + threadIdx.x * maxcf;
   //double psifp [MAX_CF];

   double *sCoefpsi = sPsifp + blockDim.x * maxcf;
   double *coefpsi = sCoefpsi + threadIdx.x * (maxCorner * maxcf);
   //double coefpsi[MAXC * MAX_CF];

   double *sVolume = sCoefpsi + blockDim.x * (maxCorner * maxcf);
   double *lVolume = sVolume + threadIdx.x * maxCorner;

   // TODO: Volume can be divided by number of groups/block.
   // That is, all groups of the same zone share the volumes.
   int *sEz_exit = (int *) (sVolume + blockDim.x * (maxCorner));
   int *ez_exit = sEz_exit + threadIdx.x * (maxCorner * maxcf);
   //int ez_exit[MAXC * MAX_CF];

   int *sNxez = sEz_exit + blockDim.x * (maxCorner * maxcf);
   int *nxez = sNxez + threadIdx.x * maxCorner;
   //int nxez[MAXC];

   area_opp = 0;
   psi_opp = 0;

   __syncthreads();

   // Track the location of the start of the hyperplane zone in nextZ list.
   //unsigned int ndoneZ = 0;
   int ndoneZ = 0;

   int hyperPlane;

   // Loop over hyperplanes
   for (hyperPlane = 0; hyperPlane < nHyperPlanes; hyperPlane++)
   {
      // break the hyperplane up into chunks of zones of size zoneBatch
      // compute the number of chunks in this zone

      int batchNum = (nZonesInPlane[hyperPlane] + zoneBatch - 1)
                     / zoneBatch; // This is how many passes are required for this hyperplane

      for (int iBatch = 0; iBatch < batchNum; iBatch++)
      {
         // check that this zone is still in the hyperplane
         if (iBatch * zoneBatch + ii < nZonesInPlane[hyperPlane])
         {
            // Loop over zones in hyperplane - IMPLICIT

            zone0 = nextZ[(iBatch * zoneBatch + ii + ndoneZ)];
            zone = abs(zone0);
            nCorner = d_GnumCorner[zone - 1];
            c0 = d_GcOffSet[zone - 1];
            //hyp_c0  = 0;

            nxBdy = 0;

            // Contributions from volume terms
            for (int c = 0; c < nCorner; c++)
            {
               lVolume[c] = Volume[c0 + c];
               source = STotal[Groups * (c0 + c) + g]; // RCC tau*Psi already computed into STotal in CompSTotalTauPsi
               Q[c] = source;
               src[c] = lVolume[c] * source;
               SigtVol[c] = Sigt[(zone - 1) * Groups + g] * lVolume[c];
            }

            for (int qq = 0; qq < maxCorner; ++qq)
            {
               nxez[qq] = 0;
            }

            // Cornerloop
            for (int c = 0; c < nCorner; c++)
            {
               // Calculate Area CornerFace do Omega to determine the contributsions from incident fluxes across external corner faces (FP faces)

               sumArea[c] = 0;

               nCFaces = nCFacesArray[c0 + c];

               for (int cface = 0; cface < nCFaces; cface++)
               {
                  afp[cface] = d_afp[(c0 + c) * ndim + cface];

                  int cfp = cFP[(c0 + c) * maxcf + cface] - 1; // cFP is 1-indexed

                  if (afp[cface] > 0)
                  {
                     sumArea[c] += afp[cface];

                     if (cfp >= ncorner)
                     {
                        bdy_exit[2 * nxBdy + ii * 2 * maxcf * maxCorner
                                 + blockIdx.x * (zoneBatch * 2 * maxcf * maxCorner)]
                            = c;
                        bdy_exit[2 * nxBdy + 1 + ii * 2 * maxcf * maxCorner
                                 + blockIdx.x * (zoneBatch * 2 * maxcf * maxCorner)]
                            = cfp - ncorner;

                        nxBdy += 1;
                     }
                  }
                  else if (afp[cface] < 0)
                  {
                     psifp[cface] = Psi1[cfp * Groups + g];
                     src[c] -= afp[cface] * Psi1[cfp * Groups + g];

                  } // endif (afp[cface] > 0 )

               } // end for ( int cface ... )

               // Contributions from interior corner faces (EZ faces)

               for (int cface = 0; cface < nCFaces; cface++)
               {
                  double aez = d_aez[(c0 + c) * ndim + cface];

                  cez = cEZ[(c0 + c) * maxcf + cface] - 1;

                  if (cez > c)
                  {
                     if (aez > (double) 0.0)
                     {
                        nxez[c] = nxez[c] + 1;
                        ez_exit[c * maxcf + nxez[c] - 1] = cez;
                        coefpsi[c * maxcf + nxez[c] - 1] = aez;
                     }
                     else if (aez < (double) 0.0)
                     {
                        nxez[cez] = nxez[cez] + 1;
                        ez_exit[cez * maxcf + nxez[cez] - 1] = c;
                        coefpsi[cez * maxcf + nxez[cez] - 1] = -aez;
                     }
                  }

                  if (aez > 0)
                  {
                     sumArea[c] += aez;

                     area_opp = 0;

                     int ifp;

                     if (nCFaces == 3)
                     {
                        ifp = (cface + 1) % (nCFaces);
                        if (afp[ifp] < 0)
                        {
                           psi_opp = psifp[ifp];
                           area_opp = -afp[ifp];
                        }
                     }
                     else
                     {
                        ifp = cface;
                        area_opp = 0;
                        psi_opp = 0;

                        for (int k = 0; k < nCFaces - 2; k++)
                        {
                           ifp = (ifp + 1) % nCFaces;
                           if (afp[ifp] < 0)
                           {
                              area_opp -= afp[ifp];
                              psi_opp -= afp[ifp] * psifp[ifp];
                           }
                        }

                        if (area_opp > 0)
                        {
                           double area_inv = 1.0 / area_opp;
                           psi_opp *= area_inv;
                        }
                     }

                     double sez;

                     // TestOppositeFace
                     if (area_opp > 0.0)
                     {
                        double aez2 = aez * aez;

                        double sig = Sigt[(zone - 1) * Groups + g];
                        double sigv = sig * lVolume[c];
                        double sigv2 = sigv * sigv;

                        double gnum = aez2 * (FOURALPHA * sigv2 + aez * (4 * sigv + 3 * aez));
                        double gden
                            = lVolume[c] * (4.0 * sigv * sigv2 + aez * (6 * sigv2 + 2.0 * aez * (2 * sigv + aez)));
                        sez = (lVolume[c] * gnum * (sig * psi_opp - Q[c]) + 0.5 * aez * gden * (Q[c] - Q[cez]))
                              / (gnum + gden * sig);

                        src[c] += sez;
                        src[cez] -= sez;
                     }
                     else
                     {
                        double sigInv = 1.0 / Sigt[(zone - 1) * Groups + g];
                        sez = 0.5 * aez * sigInv * (Q[c] - Q[cez]);
                        src[c] += sez;
                        src[cez] -= sez;

                     } // endif TestOppositeFace

                  } // endif ( aez > 0)

               } // endloop cface

            } // endloop CornerLoop

            if (zone0 > 0)
            {
               for (int i = 0; i < nCorner; i++)
               {
                  int c = nextC[c0 + i] - 1;

                  // Corner angular flux
                  Psi1[(c0 + c) * Groups + g] = src[c] / (sumArea[c] + SigtVol[c]);

                  // Scalar Flux
                  atomicAdd(Phi + ((c0 + c) * Groups + g), quadwt * Psi1[(c0 + c) * Groups + g]);

                  // Calculate the contribution of this flux to the sources of
                  // downstream corners in this zone.  The downstream corner index is
                  // "ez_exit"

                  for (int cface = 0; cface < nxez[c]; cface++)
                  {
                     cez = ez_exit[c * maxcf + cface];
                     src[cez] += coefpsi[c * maxcf + cface] * Psi1[(c0 + c) * Groups + g];
                  }
               }
            }

            else
            {
               // Direct Solve (non-lower triangular, use old values of Psi1 for now)
               for (int c = 0; c < nCorner; c++)
               {
                  for (int cface = 0; cface < nxez[c]; cface++)
                  {
                     cez = ez_exit[c * maxcf + cface];
                     src[cez] += coefpsi[c * maxcf + cface] * Psi1[(c0 + c) * Groups + g];
                  }
               }

               for (int c = 0; c < nCorner; c++)
               {
                  Psi1[(c0 + c) * Groups + g] = src[c] / (sumArea[c] + SigtVol[c]);
                  atomicAdd(Phi + ((c0 + c) * Groups + g), quadwt * Psi1[(c0 + c) * Groups + g]);
               }
            }

            // Set exiting boundary fluxes
            for (int ib = 0; ib < nxBdy; ib++)
            {
               int c = bdy_exit[2 * ib + ii * 2 * maxcf * maxCorner + blockIdx.x * (zoneBatch * 2 * maxcf * maxCorner)];
               int b = bdy_exit[2 * ib + 1 + ii * 2 * maxcf * maxCorner
                                + blockIdx.x * (zoneBatch * 2 * maxcf * maxCorner)];
               PsiB[b * Groups + g] = Psi1[(c0 + c) * Groups + g];
            }

         } // if zone is in the hyperplane

      } // end zone batch loop

      __syncthreads();

      ndoneZ += nZonesInPlane[hyperPlane];

   } // end HyperPlaneLoop
}

static cudaStream_t sweepStream[MAX_CUDA_STREAMS];

void gpu_streamsynchronize(int *streamId)
{
   CUDA_SAFE_CALL(cudaStreamSynchronize(sweepStream[*streamId]));
}

void gpu_devicesynchronize()
{
   CUDA_SAFE_CALL(cudaDeviceSynchronize());
}

static int doonce = 0;

void gpu_sweepucbxyz(int *Angle,
                     int *nHyperPlanes,
                     int *nZonesInPlane,
                     int *nextZ,
                     int *nextC,
                     double *STotal,
                     double *tau,
                     double *Psi,
                     int *Groups,
                     double *Volume,
                     double *Sigt,
                     int *nCFacesArray,
                     int *ndim,
                     int *maxcf,
                     int *ncorner,
                     double *A_fp,
                     double *omega,
                     int *cFP,
                     double *Psi1,
                     int *nbelem,
                     double *A_ez,
                     int *cEZ,
                     int *NumAngles,
                     double *quadwt,
                     double *Phi,
                     double *PsiB,
                     int *maxCorner,
                     int *mem0solve1,
                     int *streamIdPtr,
                     int *totalStreams,
                     int *savePsi,
                     int *numCycles,
                     int *cycleOffSet,
                     double *cyclePsi,
                     int *cycleList,
                     int *b0,
                     int *nBdyElem,
                     double *PsiBMref,
                     int *Mref,
                     int *Geom_numCorner,
                     int *Geom_cOffSet)

{
   // static pointers to buffers    TODO: Size and allocate pinned and device memory at startup.
   static double *d_aez[MAX_CUDA_STREAMS];
   static double *d_afp[MAX_CUDA_STREAMS];
   //static double *h_aez[MAX_CUDA_STREAMS];
   //static double *h_afp[MAX_CUDA_STREAMS];
   //static double *h_STotal[MAX_CUDA_STREAMS];
   static double *d_Psi1[MAX_CUDA_STREAMS];
   static double *d_STotal[MAX_CUDA_STREAMS];
   static double *d_Psi[MAX_CUDA_STREAMS];
   static double *d_Phi[MAX_CUDA_STREAMS];
   static double *d_PsiB[MAX_CUDA_STREAMS];
   static double *d_Volume[MAX_CUDA_STREAMS];
   static double *h_Volume[MAX_CUDA_STREAMS];
   static double *d_Sigt[MAX_CUDA_STREAMS];
   static double *d_omega[MAX_CUDA_STREAMS];
   static double *h_omega[MAX_CUDA_STREAMS];
   static double *d_cyclePsi[MAX_CUDA_STREAMS];
   static int *d_cycleList[MAX_CUDA_STREAMS];
   static int *d_nZonesInPlane[MAX_CUDA_STREAMS];
   static int *h_nZonesInPlane[MAX_CUDA_STREAMS];
   static int *d_nextZ[MAX_CUDA_STREAMS];
   static int *h_nextZ[MAX_CUDA_STREAMS];
   static int *d_nextC[MAX_CUDA_STREAMS];
   static int *h_nextC[MAX_CUDA_STREAMS];
   static int *d_nCFacesArray[MAX_CUDA_STREAMS];
   static int *d_cFP[MAX_CUDA_STREAMS];
   static int *d_bdy_exit[MAX_CUDA_STREAMS];
   static int *d_cEZ[MAX_CUDA_STREAMS];

   static int *d_GnumCorner[MAX_CUDA_STREAMS];
   static int *d_GcOffSet[MAX_CUDA_STREAMS];

   static int streamInitialized[MAX_CUDA_STREAMS];
   //static int subStream[MAX_CUDA_STREAMS];
   static cudaEvent_t Stream_Event[MAX_CUDA_STREAMS];
   //static cudaEvent_t PsiB_Event[MAX_CUDA_STREAMS];
   //static cudaEvent_t Phi_Event[MAX_CUDA_STREAMS];

   static double *d_A_fp[MAX_CUDA_STREAMS];
   static double *d_A_ez[MAX_CUDA_STREAMS];

   static int eventId = 0;
   //static int afpAezDefined = 0;

   //cudaError_t cuerr;

   int bytes96k = 98304;

   // TRAP ERRORS
   int totalZones = 0;
   for (int i = 0; i < *nHyperPlanes; i++)
   {
      totalZones += nZonesInPlane[i];
   }

   int numberOfBlocks = (*Groups + GROUPS_IN_BLOCK - 1) / GROUPS_IN_BLOCK;
   int numZonesInThreadBlock = THREAD_BLOCK_SIZE / GROUPS_IN_BLOCK;

   // now this flag means 'first angle'
   if (*mem0solve1 == 1)
   {
   }

   // set the streamId to the thread Id
   int streamId = (*streamIdPtr) % MAX_CUDA_STREAMS;

   // COMPUTE SIZES  ====================================================================================

   if (*nHyperPlanes >= MAX_PLANES)
   {
      printf("ERROR insufficient buffer size \n\n\n\n\n");
      abort();
   }

   // initialize CUDA variables
   if (!(streamInitialized[streamId] == 789))
   {
      printf("Initialization for stream %d \n", streamId);

      // create the streams
      CUDA_SAFE_CALL(cudaStreamCreate(&sweepStream[streamId]));

      CUDA_SAFE_CALL(cudaMalloc(&(d_nZonesInPlane[streamId]), MAX_PLANES * sizeof(int)));
      printf("allocated d_nZonesInPlane  = %lu bytes\n", MAX_PLANES * sizeof(int));

      CUDA_SAFE_CALL(cudaMallocHost(&(h_nZonesInPlane[streamId]), MAX_PLANES * sizeof(int)));

      CUDA_SAFE_CALL(cudaMalloc(&(d_nextZ[streamId]), totalZones * sizeof(int)));
      printf("allocated d_nextZ          = %lu bytes\n", totalZones * 5 * sizeof(int));

      CUDA_SAFE_CALL(cudaMallocHost(&(h_nextZ[streamId]), totalZones * sizeof(int)));

      CUDA_SAFE_CALL(cudaMalloc(&(d_nextC[streamId]), *ncorner * sizeof(int)));
      printf("allocated d_nextC          = %lu bytes\n", *ncorner * sizeof(int));

      CUDA_SAFE_CALL(cudaMallocHost(&(h_nextC[streamId]), *ncorner * sizeof(int)));

      CUDA_SAFE_CALL(cudaMalloc(&(d_nCFacesArray[streamId]), *ncorner * sizeof(int)));
      printf("allocated d_nCFacesArray   = %lu bytes\n", *ncorner * sizeof(int));

      CUDA_SAFE_CALL(cudaMalloc(&(d_cFP[streamId]), *maxcf * *ncorner * sizeof(int)));
      printf("allocated d_cFP            = %lu bytes\n", *maxcf * *ncorner * sizeof(int));

      CUDA_SAFE_CALL(cudaMalloc(&(d_bdy_exit[streamId]),
                                numberOfBlocks * numZonesInThreadBlock * 2 * *maxcf * *maxCorner * sizeof(int)));
      printf("allocated d_bdy_exit       = %lu bytes\n",
             numberOfBlocks * numZonesInThreadBlock * 2 * *maxcf * *maxCorner * sizeof(int));

      CUDA_SAFE_CALL(cudaMalloc(&(d_cEZ[streamId]), *maxcf * *ncorner * sizeof(int)));
      printf("allocated d_cEZ            = %lu bytes\n", *maxcf * *ncorner * sizeof(int));

      CUDA_SAFE_CALL(cudaMalloc(&(d_Psi1[streamId]), *Groups * (*ncorner + *nbelem) * sizeof(double)));
      printf("allocated d_Psi1           = %lu bytes\n", *Groups * (*ncorner + *nbelem) * sizeof(double));

      CUDA_SAFE_CALL(cudaMalloc(&(d_STotal[streamId]), *Groups * *ncorner * sizeof(double)));
      printf("allocated d_STotal         = %lu bytes\n", *Groups * *ncorner * sizeof(double));

      CUDA_SAFE_CALL(cudaMalloc(&(d_Psi[streamId]), *Groups * *ncorner * sizeof(double) / 2));
      printf("allocated d_Psi            = %lu bytes\n", *Groups * *ncorner * sizeof(double) / 2);

      CUDA_SAFE_CALL(cudaMalloc(&(d_Phi[streamId]), *Groups * *ncorner * sizeof(double)));
      printf("allocated d_Phi            = %lu bytes\n", *Groups * *ncorner * sizeof(double));

      CUDA_SAFE_CALL(cudaMalloc(&(d_PsiB[streamId]), *Groups * *nbelem * sizeof(double)));
      printf("allocated d_PsiB           = %lu bytes\n", *Groups * *nbelem * sizeof(double));

      CUDA_SAFE_CALL(cudaMalloc(&(d_Sigt[streamId]), *Groups * totalZones * sizeof(double)));
      printf("allocated d_Sigt           = %lu bytes\n", *Groups * totalZones * sizeof(double));

      CUDA_SAFE_CALL(cudaMalloc(&(d_Volume[streamId]), *ncorner * sizeof(double)));
      printf("allocated d_Volume         = %lu bytes\n", *ncorner * sizeof(double));

      CUDA_SAFE_CALL(cudaMallocHost(&(h_Volume[streamId]), *ncorner * sizeof(double)));

      CUDA_SAFE_CALL(cudaMalloc(&(d_A_fp[streamId]), *ndim * *maxcf * *ncorner * sizeof(double)));
      printf("allocated d_A_fp           = %lu bytes\n", *ndim * *maxcf * *ncorner * sizeof(double));

      CUDA_SAFE_CALL(cudaMalloc(&(d_A_ez[streamId]), *ndim * *maxcf * *ncorner * sizeof(double)));
      printf("allocated d_A_ez           = %lu bytes\n", *ndim * *maxcf * *ncorner * sizeof(double));

      //CUDA_SAFE_CALL ( cudaMemcpy (d_A_fp[0],          A_fp,          *ndim * *maxcf * *ncorner * sizeof(double),           cudaMemcpyHostToDevice ) );
      //CUDA_SAFE_CALL ( cudaMemcpy (d_A_ez[0],          A_ez,          *ndim * *maxcf * *ncorner * sizeof(double),           cudaMemcpyHostToDevice ) );

      CUDA_SAFE_CALL(cudaMalloc(&(d_afp[streamId]), *ncorner * *maxcf * sizeof(double)));
      printf("allocated d_afp            = %lu bytes\n", *ncorner * *maxcf * sizeof(double));

      CUDA_SAFE_CALL(cudaMalloc(&(d_aez[streamId]), *ncorner * *maxcf * sizeof(double)));
      printf("allocated d_aez            = %lu bytes\n", *ncorner * *maxcf * sizeof(double));

      CUDA_SAFE_CALL(cudaMalloc(&(d_omega[streamId]), 3 * sizeof(double)));
      printf("allocated d_omega          = %lu bytes\n", 3 * sizeof(double));

      CUDA_SAFE_CALL(cudaMallocHost(&(h_omega[streamId]), 3 * sizeof(double)));

      CUDA_SAFE_CALL(cudaMalloc(&(d_cycleList[streamId]), MAX_TOTAL_CYCLES * sizeof(int)));
      printf("allocated d_cycleList      = %lu bytes\n", MAX_TOTAL_CYCLES * sizeof(int));

      CUDA_SAFE_CALL(cudaMalloc(&(d_cyclePsi[streamId]), MAX_TOTAL_CYCLES * *Groups * sizeof(double)));
      printf("allocated d_cyclePsi       = %lu bytes\n", MAX_TOTAL_CYCLES * *Groups * sizeof(double));

      CUDA_SAFE_CALL(cudaMalloc(&(d_GnumCorner[streamId]), totalZones * sizeof(int)));
      printf("allocated d_GnumCorner           = %lu bytes\n", totalZones * sizeof(int));

      CUDA_SAFE_CALL(cudaMalloc(&(d_GcOffSet[streamId]), totalZones * sizeof(int)));
      printf("allocated d_GcOffSet           = %lu bytes\n", totalZones * sizeof(int));

      if (Stream_Event[streamId] == NULL)
      {
         for (int ii = 0; ii < MAX_CUDA_STREAMS; ii++)
         {
            CUDA_SAFE_CALL(cudaEventCreate(&Stream_Event[ii]));
         }
      }

      //CUDA_SAFE_CALL ( cudaEventCreate ( &Psi1_Event[streamId]) );
      //CUDA_SAFE_CALL ( cudaEventCreate ( &PsiB_Event[streamId]) );
      //CUDA_SAFE_CALL ( cudaEventCreate ( &Phi_Event[streamId]) );

      streamInitialized[streamId] = 789;
      //subStream[streamId] = 0;

      //int sharedMemSize = 10 * THREAD_BLOCK_SIZE;
      cudaFuncSetAttribute(SweepUCBxyzKernel, cudaFuncAttributeMaxDynamicSharedMemorySize, bytes96k);
   }

   // Copy nZonesInPlane to pinned memory
   for (int qq = 0; qq < *nHyperPlanes; ++qq)
   {
      h_nZonesInPlane[streamId][qq] = nZonesInPlane[qq];
   }

   // Task - rearrange to conserve memory
   //
   //  make it like
   //
   // 1. initFromCycleListKernel
   //
   // 2. compAfp
   //
   // 3. sweepkernel
   //
   // 4. UpdateCycleListKernel
   //

   CUDA_SAFE_CALL(cudaEventSynchronize(Stream_Event[eventId]));
   eventId++;
   if (eventId == MAX_CUDA_STREAMS)
   {
      eventId = 0;
   }

   // H2D MEMCPY  ====================================================================================

   if (*numCycles >= MAX_TOTAL_CYCLES || (*numCycles + *cycleOffSet) >= MAX_TOTAL_CYCLES)
   {
      printf("ERROR - numCycles %d > MAX_TOTAL_CYCLES %d - ABORTING \n", *numCycles, MAX_TOTAL_CYCLES);
      printf("ERROR - This error was caught GPU_SweepUCBxyz.cu.\n");
      printf("ERROR - The fix it just to edit GPU_SweepUCBxyz.cu and make MAX_TOTAL_CYCLES sufficiently large. \n");
      abort();
   }

   // SNREFLECT
   CUDA_SAFE_CALL(cudaMemcpyAsync(
       d_PsiB[streamId], PsiB, *nbelem * *Groups * sizeof(double), cudaMemcpyHostToDevice, sweepStream[streamId]));

   if (*nBdyElem > 0)
   {
      CUDA_SAFE_CALL(cudaMemcpyAsync(d_PsiB[streamId] + *b0,
                                     PsiBMref + *b0,
                                     *Groups * *nBdyElem * sizeof(double),
                                     cudaMemcpyHostToDevice,
                                     sweepStream[streamId]));
   }

   // now this flag means 'first angle'
   if (*mem0solve1 == 1 && doonce == 0)
   {
      doonce = 1;
      // reset Psi1 to zero if this is the 'first angle' of an angleset
      CUDA_SAFE_CALL(
          cudaMemsetAsync(d_Psi1[streamId], 0, *Groups * (*ncorner + *nbelem) * sizeof(double), sweepStream[streamId]));

      // and Phi
      CUDA_SAFE_CALL(cudaMemsetAsync(d_Phi[streamId], 0, *Groups * *ncorner * sizeof(double), sweepStream[streamId]));
   }
   else // need to copy over Psi1 and Phi if not the 'first angle'
   {
      CUDA_SAFE_CALL(cudaMemcpyAsync(
          d_Psi1[streamId], Psi1, *ncorner * *Groups * sizeof(double), cudaMemcpyHostToDevice, sweepStream[streamId]));
      CUDA_SAFE_CALL(cudaMemcpyAsync(
          d_Phi[streamId], Phi, *ncorner * *Groups * sizeof(double), cudaMemcpyHostToDevice, sweepStream[streamId]));
   }

   CUDA_SAFE_CALL(cudaMemcpyAsync(d_Psi1[streamId] + *Groups * *ncorner,
                                  d_PsiB[streamId],
                                  *nbelem * *Groups * sizeof(double),
                                  cudaMemcpyDeviceToDevice,
                                  sweepStream[streamId]));

   CUDA_SAFE_CALL(cudaMemcpyAsync(d_cycleList[streamId],
                                  cycleList,
                                  (*numCycles + *cycleOffSet) * sizeof(int),
                                  cudaMemcpyHostToDevice,
                                  sweepStream[streamId]));

   CUDA_SAFE_CALL(cudaMemcpyAsync(d_cyclePsi[streamId],
                                  cyclePsi,
                                  (*numCycles + *cycleOffSet) * *Groups * sizeof(double),
                                  cudaMemcpyHostToDevice,
                                  sweepStream[streamId]));

   // Init cycleList
   int nICLblocks
       = (*numCycles * *maxCorner * *Groups + 1023)
         / 1024; // RCC TODO: this gives the exact number of blocks required so that cycleStep is not needed, should be numCycles+cycleOffSet, remove maxCorner, make nICLblocks up to 80 to use all available SMs

   if (*Groups > 4 * 1024) // Groups needs to be less than total number of threads
   {
      printf(
          "ERROR - Number of Groups greater than number of threads (4 blocks * 1024 threads per block), aborting before initFromCycleListKernel, please increase number of blocks. \n");
      abort();
   }

   initFromCycleListKernel<<<4, 1024, 0, sweepStream[streamId]>>>(
       *numCycles, *cycleOffSet, d_cycleList[streamId], d_cyclePsi[streamId], d_Psi1[streamId], *maxCorner, *Groups);

   h_omega[streamId][0] = omega[0];
   h_omega[streamId][1] = omega[1];
   h_omega[streamId][2] = omega[2];

   CUDA_SAFE_CALL(cudaMemcpyAsync(
       d_omega[streamId], h_omega[streamId], 3 * sizeof(double), cudaMemcpyHostToDevice, sweepStream[streamId]));
   CUDA_SAFE_CALL(cudaMemcpyAsync(d_A_fp[streamId],
                                  A_fp,
                                  *ndim * *maxcf * *ncorner * sizeof(double),
                                  cudaMemcpyHostToDevice,
                                  sweepStream[streamId]));
   CUDA_SAFE_CALL(cudaMemcpyAsync(d_A_ez[streamId],
                                  A_ez,
                                  *ndim * *maxcf * *ncorner * sizeof(double),
                                  cudaMemcpyHostToDevice,
                                  sweepStream[streamId]));
   //int nblocksafp = (*ncorner * *ndim + 511)/512;

   // compAfpKernel
   compAfpKernel<<<4, 512, 0, sweepStream[streamId]>>>(d_omega[streamId],
                                                       d_A_fp[streamId],
                                                       d_afp[streamId],
                                                       d_A_ez[streamId],
                                                       d_aez[streamId],
                                                       *ncorner * *ndim,
                                                       *ndim);

   CUDA_SAFE_CALL(cudaMemcpyAsync(
       d_STotal[streamId], STotal, *ncorner * *Groups * sizeof(double), cudaMemcpyHostToDevice, sweepStream[streamId]));

   // compSTotal
   // RCC - Steve says Psi split in half for CompSTotalTauPsi to save memory (copied into same array)
   CUDA_SAFE_CALL(cudaMemcpyAsync(
       d_Psi[streamId], Psi, *ncorner * *Groups * sizeof(double) / 2, cudaMemcpyHostToDevice, sweepStream[streamId]));

   CUDA_SAFE_CALL(cudaEventRecord(Stream_Event[eventId], sweepStream[streamId]));

   CompSTotalTauPsi<<<4, 1024, 0, sweepStream[streamId]>>>(
       d_STotal[streamId], d_Psi[streamId], *tau, *ncorner * (*Groups) / 2);

   CUDA_SAFE_CALL(cudaMemcpyAsync(d_Psi[streamId],
                                  (Psi + (*ncorner * (*Groups) / 2)),
                                  *ncorner * *Groups * sizeof(double) / 2,
                                  cudaMemcpyHostToDevice,
                                  sweepStream[streamId]));

   CompSTotalTauPsi<<<4, 1024, 0, sweepStream[streamId]>>>(
       d_STotal[streamId] + *ncorner * (*Groups) / 2, d_Psi[streamId], *tau, *ncorner * (*Groups) / 2);

   CUDA_SAFE_CALL(cudaMemcpyAsync(d_nZonesInPlane[streamId],
                                  h_nZonesInPlane[streamId],
                                  *nHyperPlanes * sizeof(int),
                                  cudaMemcpyHostToDevice,
                                  sweepStream[streamId]));
   CUDA_SAFE_CALL(cudaMemcpyAsync(
       d_nextZ[streamId], nextZ, totalZones * sizeof(int), cudaMemcpyHostToDevice, sweepStream[streamId]));
   CUDA_SAFE_CALL(cudaMemcpyAsync(
       d_nextC[streamId], nextC, *ncorner * sizeof(int), cudaMemcpyHostToDevice, sweepStream[streamId]));
   CUDA_SAFE_CALL(cudaMemcpyAsync(
       d_Volume[streamId], Volume, *ncorner * sizeof(double), cudaMemcpyHostToDevice, sweepStream[streamId]));
   CUDA_SAFE_CALL(cudaMemcpyAsync(
       d_Sigt[streamId], Sigt, totalZones * *Groups * sizeof(double), cudaMemcpyHostToDevice, sweepStream[streamId]));
   CUDA_SAFE_CALL(cudaMemcpyAsync(
       d_nCFacesArray[streamId], nCFacesArray, *ncorner * sizeof(int), cudaMemcpyHostToDevice, sweepStream[streamId]));
   CUDA_SAFE_CALL(cudaMemcpyAsync(
       d_cEZ[streamId], cEZ, *maxcf * *ncorner * sizeof(int), cudaMemcpyHostToDevice, sweepStream[streamId]));
   CUDA_SAFE_CALL(cudaMemcpyAsync(
       d_cFP[streamId], cFP, *maxcf * *ncorner * sizeof(int), cudaMemcpyHostToDevice, sweepStream[streamId]));
   CUDA_SAFE_CALL(cudaMemcpyAsync(d_GnumCorner[streamId],
                                  Geom_numCorner,
                                  totalZones * sizeof(int),
                                  cudaMemcpyHostToDevice,
                                  sweepStream[streamId]));
   CUDA_SAFE_CALL(cudaMemcpyAsync(
       d_GcOffSet[streamId], Geom_cOffSet, totalZones * sizeof(int), cudaMemcpyHostToDevice, sweepStream[streamId]));

   // Sweep kernel on the GPU
   SweepUCBxyzKernel<<<numberOfBlocks, THREAD_BLOCK_SIZE, bytes96k, sweepStream[streamId]>>>(*Angle,
                                                                                             *nHyperPlanes,
                                                                                             d_nZonesInPlane[streamId],
                                                                                             d_nextZ[streamId],
                                                                                             d_nextC[streamId],
                                                                                             d_STotal[streamId],
                                                                                             *tau,
                                                                                             *Groups,
                                                                                             d_Volume[streamId],
                                                                                             d_Sigt[streamId],
                                                                                             d_nCFacesArray[streamId],
                                                                                             *ndim,
                                                                                             *maxcf,
                                                                                             *ncorner,
                                                                                             d_omega[streamId],
                                                                                             d_cFP[streamId],
                                                                                             d_bdy_exit[streamId],
                                                                                             d_Psi1[streamId],
                                                                                             d_cEZ[streamId],
                                                                                             *quadwt,
                                                                                             d_Phi[streamId],
                                                                                             d_PsiB[streamId],
                                                                                             *maxCorner,
                                                                                             d_afp[streamId],
                                                                                             d_aez[streamId],
                                                                                             d_GnumCorner[streamId],
                                                                                             d_GcOffSet[streamId]);

   // Copy data back into ORIGINAL HOST data structures.

   // Queue async memcopies
   CUDA_SAFE_CALL(cudaMemcpyAsync(
       PsiB, d_PsiB[streamId], *nbelem * *Groups * sizeof(double), cudaMemcpyDeviceToHost, sweepStream[streamId]));
   CUDA_SAFE_CALL(cudaMemcpyAsync(
       Phi, d_Phi[streamId], *ncorner * *Groups * sizeof(double), cudaMemcpyDeviceToHost, sweepStream[streamId]));

   if (*savePsi == 1)
   {
      CUDA_SAFE_CALL(cudaMemcpyAsync(
          Psi, d_Psi1[streamId], *Groups * *ncorner * sizeof(double), cudaMemcpyDeviceToHost, sweepStream[streamId]));
   }

   // Update cycleList
   updateCycleListKernel<<<nICLblocks, 1024, 0, sweepStream[streamId]>>>(
       *numCycles, *cycleOffSet, d_cycleList[streamId], d_cyclePsi[streamId], d_Psi1[streamId], *maxCorner, *Groups);

   // Update Psi1
   CUDA_SAFE_CALL(cudaMemcpyAsync(
       Psi1, d_Psi1[streamId], *Groups * *ncorner * sizeof(double), cudaMemcpyDeviceToHost, sweepStream[streamId]));

   // TODO -Need to switch to pinned memory.
   CUDA_SAFE_CALL(cudaMemcpyAsync(cyclePsi,
                                  d_cyclePsi[streamId],
                                  (*numCycles + *cycleOffSet) * *Groups * sizeof(double),
                                  cudaMemcpyDeviceToHost,
                                  sweepStream[streamId]));
}

} // extern C
