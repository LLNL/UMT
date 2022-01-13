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
   use constant_mod
   use Size_mod
   use Geometry_mod
   use QuadratureList_mod
   use ZoneData_mod
   use RadIntensity_mod

   implicit none 

!  Arguments

   integer(C_INT),    intent(in)    :: zoneID
   real(C_DOUBLE),    intent(inout) :: RadiationForce(Size%ndim,Size%maxCorner)

!  Local

   type(RadIntensity), pointer  :: RadT

   integer    :: c 
   integer    :: c0 
   integer    :: nCorner 
   integer    :: setID
   integer    :: nSets
 
!***********************************************************************
!  Compute the radiation force on the matter                           *
!***********************************************************************

   Z       => getZoneData(Geom, zoneID)

   c0      =  Z% c0
   nCorner =  Z% nCorner
   nSets   =  getNumberOfSets(Quad)

   RadiationForce(:,:) = zero

   SetLoop: do setID=1,nSets

     RadT => getRadIntensity(Quad, setID)

     do c=1,nCorner 
       RadiationForce(:,c) = RadiationForce(:,c) + RadT% RadiationForce(:,c0+c) 
     enddo

   enddo SetLoop


   return
   end subroutine getRadiationForce 


