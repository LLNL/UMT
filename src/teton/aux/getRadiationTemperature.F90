!***********************************************************************
!                         Last Update: 01/2012 PFN                     *
!                                                                      *
!    getRadiationTemperature - Called from host to get the             *
!                              radiation temperature from the          *
!                              material module.                        *
!                                                                      *
!***********************************************************************

   subroutine getRadiationTemperature(zone, radiationTemperature) &
        BIND(C,NAME="teton_getradiationtemperature")

   USE ISO_C_BINDING
   use kind_mod
   use Material_mod

   implicit none

!  Arguments

   integer(C_INT), intent(in)  :: zone
   real(C_DOUBLE), intent(out) :: radiationTemperature


   radiationTemperature = Mat% trz(zone)


   return
   end subroutine getRadiationTemperature 
