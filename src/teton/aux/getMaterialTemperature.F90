!***********************************************************************
!                                                                      *
!    getMaterialTemperature - Called from host to get the              *
!                              radiation temperature from the          *
!                              material module.                        *
!                                                                      *
!***********************************************************************

   subroutine getMaterialTemperature(zone, materialTemperature) &
        BIND(C,NAME="teton_getmaterialtemperature")

   USE ISO_C_BINDING
   use kind_mod
   use Size_mod
   use Material_mod
   use Geometry_mod
   use Material_mod

   implicit none

!  Arguments

   integer(C_INT), intent(in)  :: zone
   real(C_DOUBLE), intent(out) :: materialTemperature

!  Local

   materialTemperature = getZoneAverage(Geom, zone, Mat%Tec)
   materialTemperature  = max( materialTemperature, Size% tfloor )


   return
   end subroutine getMaterialTemperature 
