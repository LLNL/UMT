!***********************************************************************
!                        Last Update:  05/2017, PFN                    *
!                                                                      *
!   volumeUCBrz   - Calculates corner and zone volumes for the         *
!                   upstream corner-balance (UCB) transport method     *
!                   in 2D cylindrical geometry.                        *
!                                                                      *
!***********************************************************************
   subroutine volumeUCBrz

   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use BoundaryList_mod
   use Boundary_mod

   implicit none

!  Local

   integer    :: c
   integer    :: c0
   integer    :: c1
   integer    :: c2
   integer    :: cfp1
   integer    :: cfp2
   integer    :: zone  
   integer    :: n
   integer    :: nzones
   integer    :: nCorner
   integer    :: nBoundary 
   integer    :: nBdyElem 
   integer    :: b
   integer    :: b0

   real(adqt) :: r_zone, z_zone
   real(adqt) :: r_point, z_point
   real(adqt) :: r_point1, z_point1
   real(adqt) :: r_point2, z_point2
   real(adqt) :: r_edge1, z_edge1
   real(adqt) :: r_edge2, z_edge2
   real(adqt) :: area1, area2
   real(adqt) :: rbar1, rbar2
   real(adqt) :: zoneCenter(2)

!  Dynamic arrays

   real(adqt),  allocatable :: A_bdy(:,:)
   real(adqt),  allocatable :: Radius(:)

   allocate( A_bdy(2,Size% nbelem) )
   allocate( Radius(Size% nbelem) )

!  Mesh Constants

   nzones    = Size%nzones
   nBoundary = getNumberOfBoundaries(RadBoundary)

!  Calculate the geometry factors for 2D r-z geometry;
!  all of our calculations with respect to corners 

   ZoneLoop: do zone=1,nzones

     Geom% VolumeZone(zone) = zero

!    Zone Center

     zoneCenter(:) = getZoneCenter(Geom, zone)

     nCorner = Geom% numCorner(zone)
     c0      = Geom% cOffSet(zone)

     CornerLoop: do c=1,nCorner

!  For each corner "c" find the neighboring corners (in the
!  same zone) "c1" and "c2"

       c1   = Geom% cEZ(1,c0+c)
       c2   = Geom% cEZ(2,c0+c)

       cfp1 = Geom% cFP(1,c0+c)
       cfp2 = Geom% cFP(2,c0+c)

!  Load the necessary coordinates

       r_zone   = zoneCenter(1)
       z_zone   = zoneCenter(2)
       r_point  = Geom% px(1,c0+c)
       z_point  = Geom% px(2,c0+c)
       r_point1 = Geom% px(1,c0+c1)
       z_point1 = Geom% px(2,c0+c1)
       r_point2 = Geom% px(1,c0+c2)
       z_point2 = Geom% px(2,c0+c2)

       r_edge1  = half*( r_point + r_point1 )
       z_edge1  = half*( z_point + z_point1 )
       r_edge2  = half*( r_point + r_point2 )
       z_edge2  = half*( z_point + z_point2 )

       if (cfp1 > Size%ncornr) then
         cfp1          = cfp1 - Size%ncornr
         A_bdy(1,cfp1) = half*(z_point2 - z_point) 
         A_bdy(2,cfp1) = half*(r_point  - r_point2) 
         Radius(cfp1)  = half*( r_point + r_edge2 )
       endif

       if (cfp2 > Size%ncornr) then
         cfp2          = cfp2 - Size%ncornr
         A_bdy(1,cfp2) = half*(z_point  - z_point1) 
         A_bdy(2,cfp2) = half*(r_point1 - r_point) 
         Radius(cfp2)  = half*( r_point + r_edge1 )
       endif

!  Compute the corner "area" (this is really an x-y area).  The
!  corner area is the sum of two half-side areas. The half-side
!  area is the one-half the dot product of the fp-outward normal
!  and the vector from point to zone-center.

       area1 = abs( (r_edge2 - r_point)*(z_zone - z_point) -  &
                    (z_edge2 - z_point)*(r_zone - r_point) )

       area2 = abs( (r_zone - r_point)*(z_edge1 - z_point) - &
                    (z_zone - z_point)*(r_edge1 - r_point) )

!  Compute the corner "volume" (this is really volume/2*pi).
!  Each half-side volume is its area multiplied by its average
!  radial coordinate.

       rbar1             = third*(r_point + r_edge2 + r_zone)
       rbar2             = third*(r_point + r_edge1 + r_zone)
       Geom%Volume(c0+c) = half*( rbar1*area1 + rbar2*area2 )

!  Accumulate zone volume

       Geom% VolumeZone(zone) = Geom% VolumeZone(zone) + Geom%Volume(c0+c)

     enddo CornerLoop

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
   end subroutine volumeUCBrz

