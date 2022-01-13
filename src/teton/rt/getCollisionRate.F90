!***********************************************************************
!                        Last Update:  03/2013, PFN                    *
!                                                                      *
!   getCollisionRate - Computes the total collision rate. This         *
!                      quantity is used to compute the residual        *
!                      source for grey-transport acceleration (GTA).   *
!                                                                      *
!***********************************************************************
   subroutine getCollisionRate 

   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use Material_mod
   use GreyAcceleration_mod

   implicit none

!  Local

   integer    :: c
   integer    :: c0
   integer    :: nCorner
   integer    :: g
   integer    :: Groups
   integer    :: zone 
   integer    :: nZones

   real(adqt) :: sumCollisionRate

!  Constants

   nZones =  Size% nZones
   Groups =  Size% ngr 

!  Calculate the total energy absorption rate density for this set 

   ZoneLoop: do zone=1,nZones
     nCorner =  Geom% numCorner(zone)
     c0      =  Geom% cOffSet(zone)

     do c=1,nCorner
       sumCollisionRate = zero
       do g=1,Groups
         sumCollisionRate = sumCollisionRate  +   &
                           (Mat%Eta(c0+c)*Mat%siga(g,zone) + Mat%sigs(g,zone))* &
                            Geom% PhiTotal(g,c0+c)
       enddo
       GTA%  GreySource(c0+c)    = sumCollisionRate - Geom% CollisionRate(c0+c)
       Geom% CollisionRate(c0+c) = sumCollisionRate
     enddo

   enddo ZoneLoop

 
   return
   end subroutine getCollisionRate 

