!***********************************************************************
!                        Last Update:  01/2012, PFN                    *
!                                                                      *
!   ConstructBoundary - Constructor for external and shared            * 
!                       boundaries (called from C++).                  * 
!                                                                      *
!***********************************************************************


   subroutine ConstructBoundary(NumReflecting, NumVacuum,  &
                                NumSource, NumShared) &
                                BIND(C,NAME="teton_constructboundary")

!  Include

   use kind_mod
   use BoundaryList_mod

   use ISO_C_BINDING

   implicit none

!  Arguments

   integer(C_INT), intent(in)          :: NumReflecting
   integer(C_INT), intent(in)          :: NumVacuum
   integer(C_INT), intent(in)          :: NumSource
   integer(C_INT), intent(in)          :: NumShared

!  Construct the Boundary Module 

   allocate (RadBoundary)

   call construct(RadBoundary,  NumReflecting, &
                                NumVacuum,     &
                                NumSource,     &
                                NumShared) 


   return
   end subroutine ConstructBoundary

