!***********************************************************************
!                         Last Update: 01/2012 PFN                     *
!                                                                      *
!    initMaterial -  Called from host to initialize material           *
!                    properties before they are accumulated in the     *
!                    part loop.                                        *
!                                                                      *
!***********************************************************************

   subroutine initMaterial() BIND(C,NAME="teton_initmaterial_new")

   USE ISO_C_BINDING
   use constant_mod
   use Size_mod 
   use Material_mod

   implicit none

!  Local

   integer  :: zone

   do zone=1,Size%nzones

     Mat% cve(zone) = zero 
     Mat% rho(zone) = zero 
     Mat% tez(zone) = zero 
     Mat% nez(zone) = zero 
     Mat% trz(zone) = zero 

   enddo

   return
   end subroutine initMaterial 
