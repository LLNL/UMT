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

   subroutine ControlSweep(savePsi)


   use kind_mod
   use constant_mod
   use Size_mod
   use mpi_param_mod
   use mpif90_mod

   implicit none

!  Arguments

   logical (kind=1), intent(in) :: savePsi

!  Local

   real(adqt)       :: time1
   real(adqt)       :: time2
   real(adqt)       :: dtime

   logical(kind=1)  :: useGPU
   logical(kind=1)  :: useCUDASweep

!  Constants

   useGPU       = getGPUStatus(Size)
   useCUDASweep = getCUDASweepStatus(Size)

!  Loop over all group-angle sets:

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


   return
   end subroutine ControlSweep 
