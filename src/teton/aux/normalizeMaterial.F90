!***********************************************************************
!                         Last Update: 01/2012 PFN                     *
!                                                                      *
!    normalizeMaterial -  Called from host to "normalize" certain      *
!                         material properties.                         *
!                                                                      *
!***********************************************************************

   subroutine normalizeMaterial() BIND(C,NAME="teton_normalizematerial")

   use kind_mod
   use Size_mod 
   use Material_mod

   implicit none

!  Local

   integer  :: zone

   real(adqt), parameter :: cveFloor = 1.0e-6_adqt

   do zone=1,Size%nzones

     Mat% tez(zone) = Mat% tez(zone)/Mat% cve(zone) 
     Mat% trz(zone) = Mat% trz(zone)/Mat% rho(zone)
     Mat% cve(zone) = Mat% cve(zone)/Mat% rho(zone)
     Mat% cve(zone) = max(Mat% cve(zone), cveFloor)

   enddo


   return
   end subroutine normalizeMaterial 
