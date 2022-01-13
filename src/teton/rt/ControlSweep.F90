#include "macros.h"
!***********************************************************************
!                        Last Updated:  10/2016, PFN                   *
!                                                                      *
!   ControlSweep - This routine, called by the linear and GTA solvers, *
!                  solves the fixed-source transport problem on an     *
!                  arbitrary grid in either xyz-geometry or            *
!                  rz-geometry. An inner iteration is used to          *
!                  improve convergence in spatially-decomposed runs.   *
!                  An upstream corner-balance spatial discretization   *
!                  used.                                               *
!                                                                      *
!***********************************************************************

   subroutine ControlSweep(SnSweep, savePsi, PSIB, tInc)


   use kind_mod
   use constant_mod
   use Size_mod
   use QuadratureList_mod
   use GreyAcceleration_mod
   use SetData_mod
   use mpi_param_mod
   use mpif90_mod

   implicit none

!  Arguments

   real(adqt), optional, intent(inout) :: PsiB(Size%nbelem,Size%nangGTA) 
   real(adqt), optional, intent(inout) :: tInc(Size%ncornr)

   logical (kind=1),     intent(in)    :: SnSweep
   logical (kind=1),     intent(in)    :: savePsi

!  Local

   type(SetData), pointer :: Set

   integer                :: GTAsetID
   integer                :: nGTASets

   real(adqt)             :: time1
   real(adqt)             :: time2
   real(adqt)             :: dtime

   logical(kind=1)        :: useGPU
   logical(kind=1)        :: useCUDASweep

!  Constants

   nGTASets     = getNumberOfGTASets(Quad)
   useGPU       = getGPUStatus(Size)
   useCUDASweep = getCUDASweepStatus(Size)

!  Loop over all group-angle sets:

   SweepTest: if ( SnSweep ) then

     time1 = MPIWtime()

!  Hook to call optional pre-intensity iteration functions for
!  user-defined sources for test problems (i.e., MMS tests)
#if defined(TETON_ENABLE_INTENSITY_SOLVE_PREPOST_HOOKS)
     call preIntensityIterAction
#endif

     if ( useGPU ) then
       call SetSweep_GPU(savePsi)

     elseif ( useCUDASweep .and. .not. useGPU ) then

#if defined(TETON_ENABLE_CUDA)
        call SetSweep_CUDA(savePsi)
        call getPhiTotal
#else
        TETON_FATAL("Unable to use CUDA Sweep, executable was not compiled with CUDA support.")
#endif

     else ! CPU sweep

       call SetSweep(savePsi)
       call getPhiTotal

     endif

     time2 = MPIWtime()
     dtime = (time2 - time1)/sixty
     Size%SweepTimeCycle = Size%SweepTimeCycle + dtime

   else  ! GTA Solver

     if ( Size% useNewGTASolver ) then

       if ( useGPU ) then

         call GTASweep_GPU(PsiB)

       else

         GTA% PhiInc(:) = zero

         do GTAsetID=1,nGTASets
           call GTASweep(GTAsetID, PsiB)
         enddo

       endif

     else

!$omp parallel do private(GTAsetID) schedule(static)
       do GTAsetID=1,nGTASets
         call GTASweep(GTAsetID, PsiB)
       enddo
!$omp end parallel do

       tInc(:) = zero

       do GTAsetID=1,nGTASets
         Set     => getGTASetData(Quad, GTAsetID)
         tInc(:) =  tInc(:) + Set% tPhi(:)
       enddo

     endif

   endif SweepTest


   return
   end subroutine ControlSweep 
