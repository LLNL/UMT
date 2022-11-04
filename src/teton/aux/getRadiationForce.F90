!***********************************************************************
!                        Last Update:  03/2012, PFN                    *
!                                                                      *
!   getRadiationForce   - Calculates the radiation force in a corner.  *
!                                                                      *
!***********************************************************************
 
   subroutine getRadiationForce(zoneID, RadiationForce) &
        BIND(C,NAME="teton_getradiationforce")

   USE ISO_C_BINDING
   use kind_mod
   use Size_mod
   use Geometry_mod
   use RadIntensity_mod

   implicit none 

!  Arguments

   integer(C_INT),    intent(in)    :: zoneID
   real(C_DOUBLE),    intent(inout) :: RadiationForce(Size%ndim,Size%maxCorner)

!  Local

   integer    :: c 
   integer    :: c0 
   integer    :: nCorner 
 
!***********************************************************************
!  Compute the radiation force on the matter                           *
!***********************************************************************

   nCorner =  Geom% numCorner(zoneID)
   c0      =  Geom% cOffSet(zoneID)


   do c=1,nCorner 
     RadiationForce(:,c) = Rad% RadiationForce(:,c0+c) 
   enddo


   return
   end subroutine getRadiationForce 


