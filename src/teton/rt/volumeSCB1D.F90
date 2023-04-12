!***********************************************************************
!                        Last Update:  05/2017, PFN                    *
!                                                                      *
!   volumeSCB1D   - Calculates corner and zone volumes for the         *
!                   simple corner-balance (SCB) transport method       *
!                   in 1D cylindrical, spherical or planar geometry.   *
!                                                                      *
!***********************************************************************
   subroutine volumeSCB1D

   use flags_mod
   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use BoundaryList_mod
   use Boundary_mod

   implicit none

!  Local

   integer    :: c0
   integer    :: zone 
   integer    :: nzones
   integer    :: n 
   integer    :: NumBoundary

   real(adqt) :: Rmin 
   real(adqt) :: Rmax 
   real(adqt) :: Rave 

!  Mesh Constants

   nzones = Size%nzones

!  Half-cell volumes need to be multiplied by a factor of 2*pi
!  to be full cylinder volumes or 4*pi for spheres.  Note that
!  for cylinders the correct dimension of volumes and areas is
!  obtained by multiplying by dz (assumed to be 1).
 
!  Cylinders
 
   select case (Size% igeom)

     case (geometry_cylinder)
 
       do zone=1,nzones
         c0                     = 2*(zone - 1)

         Rmin                   = Geom% px(1,c0+1)
         Rmax                   = Geom% px(1,c0+2)
         Rave                   = half*( Rmin + Rmax )

         Geom% Volume(c0+1)     = half*( Rave*Rave - Rmin*Rmin )
         Geom% Volume(c0+2)     = half*( Rmax*Rmax - Rave*Rave )
         Geom% VolumeZone(zone) = Geom% Volume(c0+1) + Geom% Volume(c0+2) 
       enddo

!  Spheres
 
     case (geometry_sphere)
 
       do zone=1,nzones
         c0                     = 2*(zone - 1)

         Rmin                   = Geom% px(1,c0+1) 
         Rmax                   = Geom% px(1,c0+2) 
         Rave                   = half*( Rmin + Rmax )

         Geom% Volume(c0+1)     = third*( Rave*Rave*Rave - Rmin*Rmin*Rmin )
         Geom% Volume(c0+2)     = third*( Rmax*Rmax*Rmax - Rave*Rave*Rave )
         Geom% VolumeZone(zone) = Geom% Volume(c0+1) + Geom% Volume(c0+2)
       enddo

!  Slabs
 
     case (geometry_slab)
 
       do zone=1,nzones
         c0                     = 2*(zone - 1)

         Rmin                   = Geom% px(1,c0+1) 
         Rmax                   = Geom% px(1,c0+2) 
         Rave                   = half*( Rmin + Rmax )

         Geom% Volume(c0+1)     = Rave - Rmin
         Geom% Volume(c0+2)     = Rmax - Rave
         Geom% VolumeZone(zone) = Geom% Volume(c0+1) + Geom% Volume(c0+2)
       enddo

   end select

!  Set the outward normal vector on each boundary to be consistent
!  with the 1D geometry area definition

   NumBoundary = getNumberOfBoundaries(RadBoundary)

   do n=1,NumBoundary
     Bdy  => getBoundary(RadBoundary,n)
     zone =  Bdy% BdyToZone(1)

     if (zone == 1) then
       Bdy% A_bdy(1,1) = -one
       Bdy% Radius(1)  = Geom% px(1,1)
       Bdy% BdyToC(1)  = 1
     elseif (zone == nzones) then
       c0              = 2*(nzones - 1)
       Bdy% A_bdy(1,1) = one
       Bdy% Radius(1)  = Geom% px(1,c0+2)
       Bdy% BdyToC(1)  = Size% ncornr
     else
       call f90fatal("Invalid boundary element to zone map in RTGEOM1")
     endif
   enddo


   return
   end subroutine volumeSCB1D 

