!***********************************************************************
!                         Version 0: 11/2011 PFN                       *
!                                                                      *
!    getOpacity -  Called from host to get the total opacity           *
!                  for edit purposes only.                             *
!                                                                      *
!***********************************************************************

   subroutine getOpacity(zone, group, totalOpacity) &
        BIND(C,NAME="teton_getopacity")

   USE ISO_C_BINDING
   use kind_mod
   use constant_mod
   use Material_mod

   implicit none

!  Arguments

   integer(C_INT), intent(in)  :: zone
   integer(C_INT), intent(in)  :: group
   real(C_DOUBLE), intent(out) :: totalOpacity 

!  Local

   totalOpacity = Mat% siga(group,zone) + Mat% sigs(group,zone) + adqtEpsilon


   return
   end subroutine getOpacity 
