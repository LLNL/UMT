!***********************************************************************
!                        Last Update:  05/2017, PFN                    *
!                                                                      *
!   volumePWLDrz -   Calculates side and zone volumes for the          *
!                    piecewise linear discontinuous (PWLD) transport   *
!                    method in 2D cylindrical geometry.                *
!                                                                      *
!***********************************************************************
   subroutine volumePWLDrz 

   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use BoundaryList_mod
   use Boundary_mod

   implicit none

!  Local

   integer    :: zone
   integer    :: zoneOpp
   integer    :: side
   integer    :: nSides
   integer    :: nZones
   integer    :: vert1
   integer    :: vert2
   integer    :: n
   integer    :: nBoundary
   integer    :: nBdyElem
   integer    :: b
   integer    :: b0
   integer    :: side0

   real(adqt) :: Area
   real(adqt) :: r0, z0
   real(adqt) :: r1, z1
   real(adqt) :: r2, z2
   real(adqt) :: r_edge, z_edge

   real(adqt) :: zoneCenter(2)
   real(adqt) :: beta(Size% maxSides)

!  Dynamic arrays

   real(adqt),  allocatable :: A_bdy(:,:)
   real(adqt),  allocatable :: Radius(:)

   allocate( A_bdy(2,Size% nbelem) )
   allocate( Radius(Size% nbelem) )

!  Mesh Constants

   nZones    = Size%nZones
   nBoundary = getNumberOfBoundaries(RadBoundary)
   beta(:)   = one

!  Calculate the geometry factors for 2D r-z PWLD geometry;
!  all of our calculations with respect to sides 

   b = 0

   ZoneLoop: do zone=1,nZones

     zoneCenter(:) = getZoneCenter(Geom, zone)

!    In 2D nSides = nCorner
     nSides                 = Geom% numCorner(zone) 
     side0                  = Geom% cOffSet(zone)   
     Geom% VolumeZone(zone) = zero

     beta(1:nSides) = one/real(nSides, adqt)

     SideLoop: do side=1,nSides

       vert1  = side
       vert2  = mod(side,nSides) + 1

!  Load the necessary coordinates

       r0      = Geom% px(1,side0+vert1)
       z0      = Geom% px(2,side0+vert1)
       r1      = Geom% px(1,side0+vert2)
       z1      = Geom% px(2,side0+vert2)
       r2      = zoneCenter(1)
       z2      = zoneCenter(2)
       r_edge  = half*( r0 + r1 )
       z_edge  = half*( z0 + z1 )

       zoneOpp = Geom% zoneOpp(side,zone)

       if (zoneOpp < 0) then
         b          = b + 1
         A_bdy(1,b) = z0 - z1
         A_bdy(2,b) = r1 - r0
         Radius(b)  = half*(r0 + r1)
       endif

!  Compute the corner "area" (this is really an x-y area).  The
!  corner area is the sum of two half-side areas. The half-side
!  area is the one-half the dot product of the fp-outward normal
!  and the vector from point to zone-center.

       Area = abs( (r_edge - r0)*(z2 - z0) -  &
                   (z_edge - z0)*(r2 - r0) )

!  Compute the corner "volume" (this is really volume/2*pi).
!  Each half-side volume is its area multiplied by its average
!  radial coordinate.

       Geom%Volume(side0+side) = third*(r0 + r1 + r2)*Area

!  Accumulate zone volume
 
       Geom% VolumeZone(zone) = Geom% VolumeZone(zone) + Geom%Volume(side0+side)

     enddo SideLoop

   enddo ZoneLoop

!  Load components of the area vectors for the half-sides on
!  the problem boundary (only FP corner-faces live on the boundary)

   do n=1,nBoundary
     Bdy      => getBoundary(RadBoundary, n)
     nBdyElem =  getNumberOfBdyElements(Bdy)
     b0       =  getFirstBdyElement(Bdy) - 1
     do b=1,nBdyElem
       Bdy% A_bdy(:,b) = A_bdy(:,b+b0)
       Bdy% Radius(b)  = Radius(b+b0)
     enddo
   enddo

!  Release temporary arrays

   deallocate( A_bdy )
   deallocate( Radius )
 

   return
   end subroutine volumePWLDrz 

