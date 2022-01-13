!***********************************************************************
!                         Version 0: 04/06 PFN                         *
!                                                                      *
!    SetMaterial -  Called from host  to set material properties       *
!                   after the material state has changed.              *
!                                                                      *
!***********************************************************************

   subroutine setMaterial(zone, cve, rho, tez, trz, nez,  &
                          stimComptonMult) &
                          BIND(C,NAME="teton_setmaterial")

   use ISO_C_BINDING
   use kind_mod
   use Material_mod

   implicit none

!  Arguments

   integer(C_INT), intent(in) :: zone
   real(C_DOUBLE), intent(in) :: cve
   real(C_DOUBLE), intent(in) :: rho
   real(C_DOUBLE), intent(in) :: tez
   real(C_DOUBLE), intent(in) :: trz
   real(C_DOUBLE), intent(in) :: nez
   real(C_DOUBLE), intent(in) :: stimComptonMult 


   Mat% cve(zone)             = Mat% cve(zone) + cve
   Mat% rho(zone)             = Mat% rho(zone) + rho
   Mat% tez(zone)             = Mat% tez(zone) + tez
   Mat% nez(zone)             = Mat% nez(zone) + nez
   Mat% trz(zone)             = Mat% trz(zone) + trz
   Mat% stimComptonMult(zone) = max(Mat% stimComptonMult(zone), stimComptonMult) 

   return
   end subroutine setMaterial 
