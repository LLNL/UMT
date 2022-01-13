!***********************************************************************
!                        Last Update:  01/2012, PFN                    *
!                                                                      *
!   ConstructMaterial - Construct the material module for this         *
!                       spatial domain.                                *
!                                                                      *
!***********************************************************************


   subroutine ConstructMaterial(nonLTE) BIND(C,NAME="teton_constructmaterial")

!  Include

   use ISO_C_BINDING
   use kind_mod
   use Material_mod

   implicit none

!  Arguments
   logical(C_BOOL), intent(in) :: nonLTE

!  Construct Material Module 
   allocate(Mat)
   call Mat%construct(nonLTE)

   return
   end subroutine ConstructMaterial

