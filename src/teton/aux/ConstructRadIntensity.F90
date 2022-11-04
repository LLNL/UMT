!***********************************************************************
!                        Version 0:  02/2014, PFN                      *
!                                                                      *
!   ConstructRadIntensity - Sets F90 module pointers for the angle-    *
!                           dependent (Psi) and scalar (Phi) radiation *
!                           intensity to memory allocated by the       *
!                           host code.                                 *
!                                                                      *
!***********************************************************************


   subroutine ConstructRadIntensity() &
        BIND(C,NAME="teton_constructradintensity")

!  Include

   USE ISO_C_BINDING
   use kind_mod
   use RadIntensity_mod


   implicit none

!  Construct Radiation Intensity  Module 

   allocate(Rad)

   call construct(Rad)


   return
   end subroutine ConstructRadIntensity

