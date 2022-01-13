!***********************************************************************
!                         Version 0: 04/06 PFN                         *
!                                                                      *
!    InitOpacity -  Called from host to initialize opacities           *
!                   before they are set.                               *
!                                                                      *
!***********************************************************************

   subroutine initOpacity() BIND(C,NAME="teton_initopacity")

   use kind_mod
   use constant_mod
   use Size_mod 
   use Material_mod

   implicit none

!  Local

   integer  :: zone

   do zone=1,Size%nzones

     Mat% siga(:,zone) = zero 
     Mat% sigs(:,zone) = zero 

   enddo


   return
   end subroutine initOpacity 
