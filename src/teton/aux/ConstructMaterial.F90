!***********************************************************************
!                        Last Update:  01/2012, PFN                    *
!                                                                      *
!   ConstructMaterial - Construct the material module for this         *
!                       spatial domain.                                *
!                                                                      *
!***********************************************************************


   subroutine ConstructMaterial(nonLTE, fromRestart) BIND(C,NAME="teton_constructmaterial_new")

!  Include

   use ISO_C_BINDING
   use kind_mod
   use Material_mod

   implicit none

!  Arguments
   logical(C_BOOL), intent(in) :: nonLTE
   logical(C_BOOL), intent(in) :: fromRestart

!  Construct Material Module 
   if (.not. fromRestart) then
      allocate(Mat)
   endif

   call Mat%construct(nonLTE, fromRestart)

   return
   end subroutine ConstructMaterial

