!***********************************************************************
!                        Last Update:  05/2017, PFN                    *
!                                                                      *
!   geometryPWLDrz - Calculates certain geometry factors for the       *
!                    piecewise linear discontinuous (PWLD) transport   *
!                    method in 2D cylindrical geometry.                *
!                                                                      *
!***********************************************************************
   subroutine geometryPWLDrz 

   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod

   implicit none

!  Local

   integer    :: zone
   integer    :: nZones
   integer    :: ss
   integer    :: side0
   integer    :: side
   integer    :: nSides
   integer    :: vert1
   integer    :: vert2
   integer    :: j
   integer    :: k

   real(adqt) :: r0, z0
   real(adqt) :: r1, z1
   real(adqt) :: r2, z2
   real(adqt) :: r_edge, z_edge
   real(adqt) :: rbar, rbar0, rbar1, rbar2
   real(adqt) :: Area

   real(adqt) :: MM00, MM11, MM22
   real(adqt) :: MM01, MM02, MM12
   real(adqt) :: T00, T01, T10, T11
   real(adqt) :: Lx00, Lx01, Lx02
   real(adqt) :: Lx10, Lx11, Lx12
   real(adqt) :: Lx20, Lx21, Lx22
   real(adqt) :: Ly00, Ly01, Ly02
   real(adqt) :: Ly10, Ly11, Ly12
   real(adqt) :: Ly20, Ly21, Ly22

   real(adqt) :: zoneCenter(2)
   real(adqt) :: beta(Size% maxSides)

!  Mesh Constants

   nZones  = Size%nZones

   beta(:) = one

!  Calculate the geometry factors for 2D r-z PWLD geometry;
!  all of our calculations with respect to sides 

   ZoneLoop: do zone=1,nZones

     zoneCenter(:) = getZoneCenter(Geom, zone)

!    In 2D nSides = nCorner

     nSides = Geom% numCorner(zone) 
     side0  = Geom% cOffSet(zone) 

     beta(1:nSides) = one/real(nSides, adqt)

!    Intialize Matrices

     Geom% LL(:,:,:,zone) = zero
     Geom% AD(:,:,zone)   = zero

     if ( Size% useSurfaceMassLumping ) then
       Geom% MM(:,zone)    = zero
     else
       Geom% MMu(:,:,zone) = zero
     endif

     SideLoop: do side=1,nSides

       ss = side0 + side

!  Vertex 1 and 2 are the endpoints of the zone face 

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

!  Calculate components of the area vectors for the zone faces 

!       Geom% A_face(1,ss) = z0 - z1
!       Geom% A_face(2,ss) = r1 - r0

!  Compute the side "area" (this is really an x-y area)

       Area = abs( (r_edge - r0)*(z2 - z0) -  &
                   (z_edge - z0)*(r2 - r0) )

       Geom% Area(ss) = Area 
        
!  Compute the side "volume" (this is really volume/2*pi).
!  Each side volume is its area multiplied by its average
!  radial coordinate.

       rbar             = third*(r0 + r1 + r2)
       Geom% Volume(ss) = rbar*Area

!  Surface Integrals

       T00 = (three*r0 + r1)*one_12
       T01 = (r0 + r1)*one_12
       T10 = (r0 + r1)*one_12
       T11 = (r0 + three*r1)*one_12

!  Volume-gradient integral

       rbar0 = (two*r0 + r1 + r2)*one_24
       rbar1 = (r0 + two*r1 + r2)*one_24
       rbar2 = (r0 + r1 + two*r2)*one_24

       Lx00 =  rbar0*(z1 - z2)
       Lx01 = -rbar0*(z0 - z2)
       Lx02 =  rbar0*(z0 - z1)

       Lx10 =  rbar1*(z1 - z2)
       Lx11 = -rbar1*(z0 - z2)
       Lx12 =  rbar1*(z0 - z1)

       Lx20 =  rbar2*(z1 - z2)
       Lx21 = -rbar2*(z0 - z2)
       Lx22 =  rbar2*(z0 - z1)

       Ly00 = -rbar0*(r1 - r2)
       Ly01 =  rbar0*(r0 - r2)
       Ly02 = -rbar0*(r0 - r1)

       Ly10 = -rbar1*(r1 - r2)
       Ly11 =  rbar1*(r0 - r2)
       Ly12 = -rbar1*(r0 - r1)

       Ly20 = -rbar2*(r1 - r2)
       Ly21 =  rbar2*(r0 - r2)
       Ly22 =  rbar2*(r0 - r1)

       Geom% LL(1,vert1,vert1,zone) = Geom% LL(1,vert1,vert1,zone) + Lx00
       Geom% LL(1,vert1,vert2,zone) = Geom% LL(1,vert1,vert2,zone) + Lx01
       Geom% LL(1,vert2,vert1,zone) = Geom% LL(1,vert2,vert1,zone) + Lx10
       Geom% LL(1,vert2,vert2,zone) = Geom% LL(1,vert2,vert2,zone) + Lx11

       Geom% LL(2,vert1,vert1,zone) = Geom% LL(2,vert1,vert1,zone) + Ly00
       Geom% LL(2,vert1,vert2,zone) = Geom% LL(2,vert1,vert2,zone) + Ly01
       Geom% LL(2,vert2,vert1,zone) = Geom% LL(2,vert2,vert1,zone) + Ly10
       Geom% LL(2,vert2,vert2,zone) = Geom% LL(2,vert2,vert2,zone) + Ly11

       do j=1,nSides
         Geom% LL(1,vert1,j,zone) = Geom% LL(1,vert1,j,zone) + Lx02*beta(j)
         Geom% LL(1,vert2,j,zone) = Geom% LL(1,vert2,j,zone) + Lx12*beta(j)

         Geom% LL(1,j,vert1,zone) = Geom% LL(1,j,vert1,zone) + Lx20*beta(vert1)
         Geom% LL(1,j,vert2,zone) = Geom% LL(1,j,vert2,zone) + Lx21*beta(vert2)

         Geom% LL(2,vert1,j,zone) = Geom% LL(2,vert1,j,zone) + Ly02*beta(j)
         Geom% LL(2,vert2,j,zone) = Geom% LL(2,vert2,j,zone) + Ly12*beta(j)

         Geom% LL(2,j,vert1,zone) = Geom% LL(2,j,vert1,zone) + Ly20*beta(vert1)
         Geom% LL(2,j,vert2,zone) = Geom% LL(2,j,vert2,zone) + Ly21*beta(vert2)

         do k=1,nSides
           Geom% LL(1,j,k,zone) = Geom% LL(1,j,k,zone) + Lx22*beta(j)*beta(k)
           Geom% LL(2,j,k,zone) = Geom% LL(2,j,k,zone) + Ly22*beta(j)*beta(k)
         enddo
       enddo


!  Angular Derivative Matrix

       Geom% AD(vert1,vert1,zone) = Geom% AD(vert1,vert1,zone) + Area*one_6
       Geom% AD(vert1,vert2,zone) = Geom% AD(vert1,vert2,zone) + Area*one_12
       Geom% AD(vert2,vert1,zone) = Geom% AD(vert2,vert1,zone) + Area*one_12
       Geom% AD(vert2,vert2,zone) = Geom% AD(vert2,vert2,zone) + Area*one_6
 

       do j=1,nSides

         Geom% AD(vert1,j,zone) = Geom% AD(vert1,j,zone) + Area*one_12*beta(j)
         Geom% AD(vert2,j,zone) = Geom% AD(vert2,j,zone) + Area*one_12*beta(j)


         Geom% AD(j,vert1,zone) = Geom% AD(j,vert1,zone) + Area*one_12*beta(vert1)
         Geom% AD(j,vert2,zone) = Geom% AD(j,vert2,zone) + Area*one_12*beta(vert2)

         do k=1,nSides
           Geom% AD(j,k,zone) = Geom% AD(j,k,zone) + Area*one_6*beta(j)*beta(k)
         enddo
       enddo

!  Mass Matrix

       MM00 = Area*(three*r0 + r1 +r2)*one_30
       MM11 = Area*(r0 + three*r1 + r2)*one_30
       MM22 = Area*(r0 + r1 + three*r2)*one_30

       MM01 = Area*(two*r0 + two*r1 + r2)*one_60
       MM02 = Area*(two*r0 + r1 + two*r2)*one_60
       MM12 = Area*(r0 + two*r1 + two*r2)*one_60

       if ( Size% useSurfaceMassLumping ) then

         Geom% MM(vert1,zone) = Geom% MM(vert1,zone) + MM00 + MM01 + MM02
         Geom% MM(vert2,zone) = Geom% MM(vert2,zone) + MM11 + MM01 + MM12

         do j=1,nSides
           Geom% MM(j,zone) = Geom% MM(j,zone) + MM02*beta(vert1) +  &
                                               MM12*beta(vert2) +  &
                                               MM22*beta(j)
         enddo

       else

         Geom% MMu(vert1,vert1,zone) = Geom% MMu(vert1,vert1,zone) + MM00 
         Geom% MMu(vert1,vert2,zone) = Geom% MMu(vert1,vert2,zone) + MM01 
         Geom% MMu(vert2,vert1,zone) = Geom% MMu(vert2,vert1,zone) + MM01 
         Geom% MMu(vert2,vert2,zone) = Geom% MMu(vert2,vert2,zone) + MM11 

         do j=1,nSides
           Geom% MMu(vert1,j,zone) = Geom% MMu(vert1,j,zone) + MM02*beta(j)
           Geom% MMu(vert2,j,zone) = Geom% MMu(vert2,j,zone) + MM12*beta(j)

           Geom% MMu(j,vert1,zone) = Geom% MMu(j,vert1,zone) + MM02*beta(vert1)
           Geom% MMu(j,vert2,zone) = Geom% MMu(j,vert2,zone) + MM12*beta(vert2)

           do k=1,nSides
             Geom% MMu(j,k,zone)   = Geom% MMu(j,k,zone) + MM22*beta(j)*beta(k)
           enddo
         enddo

       endif

     enddo SideLoop

   enddo ZoneLoop
 


   return
   end subroutine geometryPWLDrz 

