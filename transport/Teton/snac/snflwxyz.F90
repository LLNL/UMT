!***********************************************************************
!                        Version 1:  09/96, PFN                        *
!                                                                      *
!   SNFLWXYZ - This routine, called by RSWPMD and RTACCELMD, solves    *
!              the fixed-source transport problem on an arbitrary      *
!              grid in either xyz-geometry or rz-geometry.             *
!              An upstream corner-balance spatial discretization is    *
!              used.                                                   *
!                                                                      *
!   Input:                                                             *
!                                                                      *
!   Output:                                                            *
!                                                                      *
!***********************************************************************

#define BATCHSIZE 30

   subroutine snflwxyz(ipath, PSIB, PSI, PHI, angleLoopTime)

   use, intrinsic :: iso_c_binding
   use kind_mod
   use constant_mod
   use Size_mod
   use Quadrature_mod
   use snswp3d_mod
   use cudafor

#include "assert.h"
!  Assertion checking include file for TETON

   implicit none
   include 'mpif.h'

!  Arguments

   real(adqt), intent(inout) :: psib(QuadSet%Groups,Size%nbelem,QuadSet%NumAngles)
   real(adqt), intent(inout) :: psi(QuadSet%Groups,Size%ncornr,QuadSet%NumAngles)
   real(adqt), intent(inout) :: Phi(QuadSet%Groups,Size%ncornr),angleLoopTime

   character(len=8), intent(in) :: ipath

!  Local

   integer          :: Angle, mm,mm1,mm2,anglebatch, n_cpuL, thnum
   integer          :: Groups, fluxIter, ishared
   integer          :: binSend, binRecv, NangBin, istat

   logical (kind=1) :: FluxConverged

   real(adqt)       :: maxFluxError
   real(adqt)       :: startOMPLoopTime, endOMPLoopTime, theOMPLoopTime

   type(C_DEVPTR)                    :: d_phi_p
   real(adqt), dimension(:,:), device, allocatable :: d_phi(:,:)

   real(adqt), device :: psiccache(QuadSet%Groups,Size%ncornr,BATCHSIZE)
   type(C_PTR) :: cptr
   type(C_DEVPTR) :: dptr

   integer :: OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
   integer angles, nbelem, ncornr, NumBin, myrank, info

!  Function

   istat = cudaHostGetDevicePointer(d_phi_p, C_LOC(phi(1,1)), 0)
   call c_f_pointer(d_phi_p, d_phi, [QuadSet%Groups,Size%ncornr] )

!  Set number of threads

!  n_cpuL = 1
   n_cpuL = OMP_GET_MAX_THREADS()
   theOMPLoopTime=0.0

   require(n_cpuL>0,   "Invalid Thread Count")
   require(n_cpuL<=32, "Invalid Thread Count") 

!  Mesh Constants

   Groups = QuadSet%Groups

   angles = QuadSet%NumAngles
   nbelem = Size%nbelem
   ncornr = Size%ncornr
   NangBin = maxval(QuadSet%NangBinList(:))
   NumBin = QuadSet%NumBin
   call mpi_comm_rank(mpi_comm_world, myrank, info)

   if (myrank .eq. 0) write(0,*) ' groups, ncornr, nbelem, angles, NangBin, NumBin = ', groups, ncornr, nbelem, angles, NangBin, NumBin
!  Loop over angle bins

   if (ipath == 'sweep') then
     call setIncidentFlux(psib)
   endif
                                                                                         
   FluxConverged = .FALSE.
   fluxIter      =  0

   call restoreCommOrder(QuadSet)


   FluxIteration: do

!    Post receives for all data
                                                                                                  
     call InitExchange

     fluxIter = fluxIter + 1

     if (ipath == 'sweep') then
       phi(:,:) = zero
     endif

!    Loop over angles, solving for each in turn:
     startOMPLoopTime = MPI_WTIME()

!!!$OMP PARALLEL DO  PRIVATE(binRecv,binSend,NangBin,mm1,mm2,anglebatch,thnum)
     AngleBin: do binRecv=1,QuadSet% NumBin
       binSend = QuadSet% SendOrder(binRecv)
       NangBin = QuadSet% NangBinList(binSend)
       
       ! loop over batches within the angle bin
       AngleLoop: do mm1=1,NangBin,BATCHSIZE

         mm2=min(mm1+BATCHSIZE-1,NangBin)
         anglebatch=mm2-mm1+1

         !for other bins, will begin staging in the data at the end of prev
         !iteration of the loop
         if (binRecv == 1) then
           do mm=mm1,mm2
             istat=cudaMemcpyAsync(psiccache(1,1,mm-mm1+1),                 &
                                   psi(1,1,QuadSet%AngleOrder(mm,binSend)), &
                                   QuadSet%Groups*Size%ncornr, 0)
           enddo
         endif

!        Set angular fluxes for reflected angles

         do mm=mm1,mm2
           call snreflect(QuadSet%AngleOrder(mm,binSend), PSIB)
         enddo

!        Sweep the mesh, calculating PSI for each corner; the 
!        boundary flux array PSIB is also updated here. 
!        Mesh cycles are fixed automatically.

         call snswp3d(anglebatch, QuadSet%AngleOrder(mm1,binSend), &
                      QuadSet%d_AngleOrder(mm1,binSend),           &
                      QuadSet%d_next,QuadSet%d_nextZ,              &
                      QuadSet%d_passZstart, psi, psiccache, psib)

         if (ipath == 'sweep') then
           call snmomentsD(psiccache, d_phi, QuadSet%d_Weight,     &
                           QuadSet%d_AngleOrder(mm1,binSend),      &
                           anglebatch) ! GPU version, one slice at a time
         endif

         if (binRecv < QuadSet% NumBin) then
           ! pre-stage data for next angle bin while exchange is happening
           do mm=mm1,mm2
             istat=cudaMemcpyAsync(psiccache(1,1,mm-mm1+1),        &
                                   psi(1,1,QuadSet%AngleOrder(mm,QuadSet%SendOrder(binRecv+1))), &
                                   QuadSet%Groups*Size%ncornr, 0)
           enddo
         endif

       enddo AngleLoop

!      Exchange Boundary Fluxes

       call exchange(PSIB, binSend, binRecv) 

     enddo AngleBin

     istat = cudaDeviceSynchronize()
     endOMPLoopTime = MPI_WTIME()
     theOMPLoopTime = theOMPLoopTime + (endOMPLoopTime-startOMPLoopTime)

     if (ipath == 'sweep') then
       call setIncidentFlux(psib)
       call testFluxConv(FluxConverged, fluxIter, maxFluxError)
     else
       FluxConverged = .TRUE.
     endif

     if ( FluxConverged ) then
       exit FluxIteration
     else
       call setCommOrder(QuadSet)
       cycle FluxIteration
     endif

   enddo FluxIteration

!  Update the scaler flux 

   if (ipath == 'sweep') then
     call restoreCommOrder(QuadSet)
   endif

   angleLoopTime = angleLoopTime + theOMPLoopTime


   return
   end subroutine snflwxyz


