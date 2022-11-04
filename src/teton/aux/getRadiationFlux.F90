!***********************************************************************
!                       Last Update:  03/2012, PFN                     *
!                                                                      *
!   getRadiationFlux    - Calculates the zone-average radiation        *
!                         flux vector.                                 * 
!                                                                      *
!***********************************************************************
 
   subroutine getRadiationFlux(zoneID, RadiationFlux) &
        BIND(C,NAME="teton_getradiationflux")

   USE ISO_C_BINDING
   use kind_mod
   use Size_mod
   use RadIntensity_mod

   implicit none 

!  Arguments

   integer(C_INT),    intent(in) :: zoneID

   real(C_DOUBLE), intent(inout) :: RadiationFlux(Size%ndim,Size%ngr)

!  Constants

!***********************************************************************
!  Compute the radiation flux                                          *
!***********************************************************************

   RadiationFlux(:,:) = Rad% RadiationFlux(:,:,zoneID) 


   return
   end subroutine getRadiationFlux 


