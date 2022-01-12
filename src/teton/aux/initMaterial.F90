!***********************************************************************
!                         Last Update: 01/2012 PFN                     *
!                                                                      *
!    initMaterial -  Called from host to initialize material           *
!                    properties before they are accumulated in the     *
!                    part loop.                                        *
!                                                                      *
!***********************************************************************

   subroutine initMaterial(Tec) BIND(C,NAME="teton_initmaterial")

   USE ISO_C_BINDING
   use kind_mod
   use constant_mod
   use Size_mod 
   use Material_mod

   implicit none

!  Arguments

   real(C_DOUBLE), intent(in)   :: Tec(Size%ncornr)

!  Local

   integer  :: zone
   integer  :: c

   do zone=1,Size%nzones

     Mat% cve(zone) = zero 
     Mat% rho(zone) = zero 
     Mat% tez(zone) = zero 
     Mat% nez(zone) = zero 
     Mat% trz(zone) = zero 

   enddo

   do c=1,Size% ncornr
     Mat%tec(c) = Tec(c)
   enddo 


   return
   end subroutine initMaterial 
