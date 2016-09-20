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

   subroutine snflwxyz(ipath, PSIB, PSI, PHI, angleLoopTime)


   use kind_mod
   use constant_mod
   use Size_mod
   use Quadrature_mod

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

   integer          :: Angle, mm
   integer          :: Groups, fluxIter, ishared
   integer          :: binSend, binRecv, NangBin

   logical (kind=1) :: FluxConverged

   real(adqt)       :: maxFluxError
   real(adqt)       :: startOMPLoopTime, endOMPLoopTime, theOMPLoopTime

   integer angles, nbelem, ncornr, NumBin, myrank, info

   theOMPLoopTime=0.0

!  Mesh Constants

   Groups = QuadSet%Groups

   angles = QuadSet%NumAngles
   nbelem = Size%nbelem
   ncornr = Size%ncornr
   NangBin = maxval(QuadSet%NangBinList(:))
   NumBin = QuadSet%NumBin
   call mpi_comm_rank(mpi_comm_world, myrank, info)

!  Loop over angle bins

   if (ipath == 'sweep') then
     call timer_beg('_setflux')
     call setIncidentFlux(psib)
     call timer_end('_setflux')
   endif
                                                                                         
   FluxConverged = .FALSE.
   fluxIter      =  0

   call restoreCommOrder(QuadSet)


   FluxIteration: do

!    Post receives for all data
                                                                                                  
     if (myrank .eq. 0) write(0,*) 'YKT: NumBin, fluxIter = ', QuadSet% NumBin, fluxIter
     call timer_beg('_initexch')
     call InitExchange
     call timer_end('_initexch')

     fluxIter = fluxIter + 1

     AngleBin: do binRecv=1,QuadSet% NumBin
       binSend = QuadSet% SendOrder(binRecv)
       NangBin = QuadSet% NangBinList(binSend)
       
!

!    Loop over angles, solving for each in turn:
     startOMPLoopTime = MPI_WTIME()
     call timer_beg('_angleloop')
     call hpm_start("sweep")
!
!$OMP PARALLEL DO PRIVATE(Angle) schedule(static,1)
       AngleLoop: do mm=1,NangBin

         Angle = QuadSet% AngleOrder(mm,binSend)

!        Set angular fluxes for reflected angles

         call snreflect(Angle, PSIB)

!        Sweep the mesh, calculating PSI for each corner; the
!        boundary flux array PSIB is also updated here.
!        Mesh cycles are fixed automatically.

         call snswp3d(Groups, Angle,                                   &
                      QuadSet%next(1,Angle),QuadSet%nextZ(1,Angle),    &
                      PSI(1,1,Angle),PSIB(1,1,Angle))

       enddo AngleLoop
     call hpm_stop("sweep")
     call timer_end('_angleloop')
     endOMPLoopTime = MPI_WTIME()
     theOMPLoopTime = theOMPLoopTime + (endOMPLoopTime-startOMPLoopTime)

!      Exchange Boundary Fluxes

       call timer_beg('_exch')
       call exchange(PSIB, binSend, binRecv) 
       call timer_end('_exch')

     enddo AngleBin

     if (ipath == 'sweep') then
       call timer_beg('_setflux')
       call setIncidentFlux(psib)
       call timer_end('_setflux')
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
     call timer_beg('_snmoments')
     call snmoments(psi, PHI)
     call timer_end('_snmoments')
     call restoreCommOrder(QuadSet)
   endif

   angleLoopTime = angleLoopTime + theOMPLoopTime


   return
   end subroutine snflwxyz


