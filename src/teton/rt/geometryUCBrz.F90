!***********************************************************************
!                        Last Update:  05/2017, PFN                    *
!                                                                      *
!   geometryUCBrz - Calculates certain geometry factors for the        *
!                   upstream corner-balance (UCB) transport method     *
!                   in 2D cylindrical geometry.                        *
!                                                                      *
!***********************************************************************
   subroutine geometryUCBrz

   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod

   implicit none

!  Local

   integer    :: c
   integer    :: c0
   integer    :: c1
   integer    :: c2
   integer    :: cfp1
   integer    :: cfp2

   integer    :: zone
   integer    :: nzones
   integer    :: nCorner

   real(adqt) :: r_zone, z_zone
   real(adqt) :: r_point, z_point
   real(adqt) :: r_point1, z_point1
   real(adqt) :: r_point2, z_point2
   real(adqt) :: r_edge1, z_edge1
   real(adqt) :: r_edge2, z_edge2
   real(adqt) :: area1, area2
   real(adqt) :: rbar1, rbar2
   real(adqt) :: zoneCenter(2)

!  Mesh Constants

   nzones = Size%nzones

!  Calculate the geometry factors for 2D r-z geometry;
!  all of our calculations with respect to corners 

   ZoneLoop: do zone=1,nzones

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

!  This is one-half the average radius for the FP corner faces 

       Geom% RadiusFP(1,c0+c) = half*( r_point + r_edge2 )
       Geom% RadiusFP(2,c0+c) = half*( r_point + r_edge1 )

!  This is one-half the average radius for the EZ corner faces

       Geom% RadiusEZ(1,c0+c) = half*( r_zone  + r_edge1 )
       Geom% RadiusEZ(2,c0+c) = half*( r_zone  + r_edge2 )

!  Calculate components of the area vectors (outward normals) 

!  Outward normals on zone faces (FP)

       if (c < cfp1) then
         Geom% A_fp(1,1,c0+c) = half*(z_point2 - z_point)
         Geom% A_fp(2,1,c0+c) = half*(r_point  - r_point2)

         if (cfp1 <= Size% ncornr) then
           Geom% A_fp(1,2,cfp1) = -Geom% A_fp(1,1,c0+c)
           Geom% A_fp(2,2,cfp1) = -Geom% A_fp(2,1,c0+c)
         endif
       endif

       if (c < cfp2) then
         Geom% A_fp(1,2,c0+c) = half*(z_point  - z_point1) 
         Geom% A_fp(2,2,c0+c) = half*(r_point1 - r_point)

         if (cfp2 <= Size% ncornr) then
           Geom% A_fp(1,1,cfp2) = -Geom% A_fp(1,2,c0+c)
           Geom% A_fp(2,1,cfp2) = -Geom% A_fp(2,2,c0+c)
         endif
       endif

!  Outward normal on interior corner-faces (EZ)

       if (c < c1) then
         Geom% A_ez(1,1,c0+c)  = z_edge1 - z_zone
         Geom% A_ez(2,1,c0+c)  = r_zone  - r_edge1

         Geom% A_ez(1,2,c0+c1) = -Geom% A_ez(1,1,c0+c)
         Geom% A_ez(2,2,c0+c1) = -Geom% A_ez(2,1,c0+c)
       endif

       if (c < c2) then
         Geom% A_ez(1,2,c0+c)  = z_zone  - z_edge2
         Geom% A_ez(2,2,c0+c)  = r_edge2 - r_zone

         Geom% A_ez(1,1,c0+c2) = -Geom% A_ez(1,2,c0+c)
         Geom% A_ez(2,1,c0+c2) = -Geom% A_ez(2,2,c0+c)
       endif

!  Compute the corner "area" (this is really an x-y area).  The
!  corner area is the sum of two half-side areas. The half-side
!  area is the one-half the dot product of the fp-outward normal
!  and the vector from point to zone-center.

       area1 = abs( (r_edge2 - r_point)*(z_zone - z_point) -  &
                    (z_edge2 - z_point)*(r_zone - r_point) )

       area2 = abs( (r_zone - r_point)*(z_edge1 - z_point) - &
                    (z_zone - z_point)*(r_edge1 - r_point) )

       Geom% Area(c0+c) = half*( area1 + area2 )

!  Compute the corner "volume" (this is really volume/2*pi).
!  Each half-side volume is its area multiplied by its average
!  radial coordinate.

       rbar1             = third*(r_point + r_edge2 + r_zone)
       rbar2             = third*(r_point + r_edge1 + r_zone)
       Geom%Volume(c0+c) = half*( rbar1*area1 + rbar2*area2 )

     enddo CornerLoop

   enddo ZoneLoop


   return
   end subroutine geometryUCBrz

