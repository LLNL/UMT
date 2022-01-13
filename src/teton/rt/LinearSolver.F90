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

   real(adqt)       :: time1
   real(adqt)       :: time2
   real(adqt)       :: dtime

   logical (kind=1) :: SnSweep

!  Constants

   parameter (SnSweep=.TRUE.)

   ndim  = Size%ndim

! TODO: Add a loop around this function.   This would allow us
! to work in full LMFG mode.

!***********************************************************************
!     BEGIN LOOP OVER BATCHES                                          *
!***********************************************************************

!  Initialize Collision Rate before we apply the transport operator
   call getCollisionRate

!  Sweep all angles in all groups in this "batch"

   START_RANGE("Teton_Sweep")

   if (ndim == 1) then
     call ControlSweep1D(savePsi)
   else
     call ControlSweep(SnSweep, savePsi)
   endif

   END_RANGE("Teton_Sweep")

   call PrintEnergies("LinearSolver, after call to MGSN ControlSweep (MGSN Sweep)")

!***********************************************************************
!     END FREQUENCY GROUP LOOP                                         *
!***********************************************************************

!  Update Collision Rate and GTA source
   call getCollisionRate

!***********************************************************************
!     GREY ACCELERATION                                                *
!***********************************************************************
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
   time1 = MPIWtime()
   START_RANGE("Teton_Grey_Transport_Acceleration")

!  Solve for corrections
   if (ndim == 1) then
     call GDASolver
   else
     call GTASolver
   endif

!  Update multigroup scalar photon intensities with grey corrections

   call addGreyCorrections

   END_RANGE("Teton_Grey_Transport_Acceleration")
   time2 = MPIWtime()

   dtime = (time2 - time1)/sixty
   Size% GTATimeCycle = Size% GTATimeCycle + dtime
#endif

   return
   end subroutine LinearSolver
