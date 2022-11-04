#include "macros.h"
!***********************************************************************
!                       Last Update:  10/2016, PFN                     *
!                                                                      *
!   LinearSolver - Perform the inner loop intensity iteration with     *
!                  grey acceleration.                                  *
!                                                                      *
!***********************************************************************

   subroutine LinearSolver(savePsi)

   use kind_mod
   use Size_mod
   use constant_mod
   use mpi_param_mod
   use mpif90_mod
#if defined(TETON_ENABLE_CALIPER)
   use caliper_mod
#endif

   implicit none

!  Arguments

   logical(kind=1), intent(in) :: savePsi

!  Local

   integer          :: ndim
   integer          :: residualFlag

   real(adqt)       :: time1
   real(adqt)       :: time2
   real(adqt)       :: dtime

   logical (kind=1) :: useGPU

!  Constants

   ndim   = Size%ndim
   useGPU = getGPUStatus(Size)

! TODO: Add a loop around this function.   This would allow us
! to work in full LMFG mode.

!***********************************************************************
!  Initialize Collision Rate before we apply the transport operator    *
!***********************************************************************

   residualFlag = 0

   if ( useGPU ) then
     call getCollisionRate_GPU(residualFlag)
   else
     call getCollisionRate(residualFlag)
   endif

!  Sweep all angles in all groups

   START_RANGE("Teton_Sweep")

   if (ndim == 1) then
     call ControlSweep1D(savePsi)
   else
     call ControlSweep(savePsi)
   endif

   END_RANGE("Teton_Sweep")

   call PrintEnergies("LinearSolver, after call to MGSN ControlSweep (MGSN Sweep)")

!***********************************************************************
!  Update Collision Rate and GTA source
!***********************************************************************

   residualFlag = 1 

   if ( useGPU ) then
     call getCollisionRate_GPU(residualFlag)
   else
     call getCollisionRate(residualFlag)
   endif

!***********************************************************************
!     GREY ACCELERATION                                                *
!***********************************************************************
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
   time1 = MPIWtime()
   START_RANGE("Teton_Grey_Transport_Acceleration")

!  Solve for corrections and update multigroup scalar photon intensities

   if ( useGPU ) then

     if ( Size% useNewGTASolver ) then
       call GTASolver_GPU
     else
       call GTASolver
     endif

     call addGreyCorrections_GPU

   else

     if (ndim == 1) then
       call GDASolver
     else
       call GTASolver
     endif

     call addGreyCorrections

   endif


   END_RANGE("Teton_Grey_Transport_Acceleration")
   time2 = MPIWtime()

   dtime = (time2 - time1)/sixty
   Size% GTATimeCycle = Size% GTATimeCycle + dtime
#endif

   return
   end subroutine LinearSolver
