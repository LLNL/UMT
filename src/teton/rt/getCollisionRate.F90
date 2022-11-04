#include "macros.h"
!***********************************************************************
!                        Last Update:  03/2013, PFN                    *
!                                                                      *
!   getCollisionRate - Computes the total collision rate. This         *
!                      quantity is used to compute the residual        *
!                      source for grey-transport acceleration (GTA).   *
!                                                                      *
!***********************************************************************
   subroutine getCollisionRate(residualFlag) 

   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use RadIntensity_mod
   use Material_mod
   use GreyAcceleration_mod
   use QuadratureList_mod
   use ZoneSet_mod

   implicit none

!  Arguments

   integer,  intent(in) :: residualFlag

!  Local

   integer    :: c
   integer    :: g
   integer    :: ngr
   integer    :: zone 
   integer    :: zSetID
   integer    :: nZoneSets

   real(adqt) :: sumCollisionRate

!  Constants

   nZoneSets = getNumberOfZoneSets(Quad)
   ngr       = Size% ngr

!  Calculate the total energy absorption rate density 

   if (residualFlag == 0) then

!$omp parallel do default(none) schedule(dynamic)  &
!$omp& shared(nZoneSets, Geom, Rad, Mat, GTA, ngr)                &
!$omp& private(zone, sumCollisionRate) 

     do zSetID=1,nZoneSets
       do c=Geom% corner1(zSetID),Geom% corner2(zSetID)

         zone             = Geom% CToZone(c)
         sumCollisionRate = zero

         do g=1,ngr
           sumCollisionRate = sumCollisionRate  +   &
                             (Mat%Eta(c)*Mat%siga(g,zone) + Mat%sigs(g,zone))* &
                              Rad% PhiTotal(g,c)
         enddo

         GTA% GreySource(c) = sumCollisionRate
       enddo
     enddo
!$omp end parallel do

   else

!$omp parallel do default(none) schedule(dynamic)  &
!$omp& shared(nZoneSets, Geom, Rad, Mat, GTA, ngr)                &
!$omp& private(zone, sumCollisionRate) 

     do zSetID=1,nZoneSets
       do c=Geom% corner1(zSetID),Geom% corner2(zSetID)

         zone             = Geom% CToZone(c)
         sumCollisionRate = zero

         do g=1,ngr
           sumCollisionRate = sumCollisionRate  +   &
                             (Mat%Eta(c)*Mat%siga(g,zone) + Mat%sigs(g,zone))* &
                              Rad% PhiTotal(g,c)
         enddo

         GTA% GreySource(c) = sumCollisionRate - GTA% GreySource(c) 

       enddo
     enddo
!$omp end parallel do

   endif

 
   return
   end subroutine getCollisionRate 

