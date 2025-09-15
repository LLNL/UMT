#include "macros.h"
!***********************************************************************
!                       Last Update:  03/2012, PFN                     *
!                                                                      *
!   RADTR  - Control program for radiation transport. It initializes   *
!            arrays, controls timestep, calls the transport package    *
!            and performs edits.                                       *
!                                                                      *
!                                                                      *
!***********************************************************************

   subroutine radtr() BIND(C,NAME="teton_radtr")

   use kind_mod
   use constant_mod
   use Size_mod
   use mpi_param_mod
   use mpif90_mod
#if defined(TETON_ENABLE_CALIPER)
   use caliper_mod
#endif

   implicit none

   real(adqt) :: time1, time2
   real(adqt) :: time3, time4

   START_RANGE("Teton_Cycle")

!***********************************************************************
!  INTERPOLATE SOURCE PROFILES                                         *
!***********************************************************************

   time1 = MPIWtime()

   call profint

!***********************************************************************
!  Initialize the radiation field and total opacity for each set       *
!*********************************************************************** 

   START_RANGE("Teton_Initialize_Sets")
   call initializeSets
   END_RANGE("Teton_Initialize_Sets")

   time2 = MPIWtime()
!***********************************************************************
!     RADIATION TRANSPORT MODULE                                       *
!***********************************************************************

   call rtmainsn

   time3 = MPIWtime()

   START_RANGE("Teton_Finalize_Sets")
   call finalizeSets
   END_RANGE("Teton_Finalize_Sets")

   time4 = MPIWtime()

!  Timings

   Size% InitTimeCycle  = (time2 - time1)/sixty
   Size% RadtrTimeCycle = (time4 - time1)/sixty
   Size% FinalTimeCycle = (time4 - time3)/sixty

   END_RANGE("Teton_Cycle")

   return
   end subroutine radtr

