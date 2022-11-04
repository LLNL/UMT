!***********************************************************************
!                        Last Update:  05/2017, PFN                    *
!                                                                      *
!   geometrySCB1D - Calculates certain geometry factors for the        *
!                   simple corner-balance (SCB) transport method       *
!                   in 1D cylindrical, spherical or planar geometry.   *
!                                                                      *
!***********************************************************************
   subroutine geometrySCB1D

   use flags_mod
   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod

   implicit none

!  Local

   integer    :: c0, zone, nzones

   real(adqt) :: Rmin, Rmax, Rave 

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
         c0                    = 2*(zone - 1)

         Rmin                  = Geom% px(1,c0+1)
         Rmax                  = Geom% px(1,c0+2)
         Rave                  = half*( Rmin + Rmax )

         Geom% Rave(zone)      = Rave
         Geom% Rmin(zone)      = Rmin
         Geom% Rmax(zone)      = Rmax
         Geom% zoneWidth(zone) = Rmax - Rmin

         Geom% Radius(c0+1)    = Rave
         Geom% Radius(c0+2)    = Rmax
         Geom% Area(c0+1)      = Rave - Rmin
         Geom% Area(c0+2)      = Rmax - Rave
         Geom% Volume(c0+1)    = half*( Rave*Rave - Rmin*Rmin )
         Geom% Volume(c0+2)    = half*( Rmax*Rmax - Rave*Rave )
       enddo

!  Spheres
 
     case (geometry_sphere)
 
       do zone=1,nzones
         c0                    = 2*(zone - 1)

         Rmin                  = Geom% px(1,c0+1) 
         Rmax                  = Geom% px(1,c0+2) 
         Rave                  = half*( Rmin + Rmax )

         Geom% Rave(zone)      = Rave*Rave
         Geom% Rmin(zone)      = Rmin*Rmin
         Geom% Rmax(zone)      = Rmax*Rmax
         Geom% zoneWidth(zone) = Rmax - Rmin

         Geom% Radius(c0+1)    = Rave
         Geom% Radius(c0+2)    = Rmax
         Geom% Area(c0+1)      = half*( Rave*Rave - Rmin*Rmin )
         Geom% Area(c0+2)      = half*( Rmax*Rmax - Rave*Rave )
         Geom% Volume(c0+1)    = third*( Rave*Rave*Rave - Rmin*Rmin*Rmin )
         Geom% Volume(c0+2)    = third*( Rmax*Rmax*Rmax - Rave*Rave*Rave )
       enddo

!  Slabs
 
     case (geometry_slab)
 
       do zone=1,nzones
         c0                    = 2*(zone - 1)

         Rmin                  = Geom% px(1,c0+1) 
         Rmax                  = Geom% px(1,c0+2) 
         Rave                  = half*( Rmin + Rmax )

         Geom% Rave(zone)      = one 
         Geom% Rmin(zone)      = one 
         Geom% Rmax(zone)      = one 
         Geom% zoneWidth(zone) = Rmax - Rmin

         Geom% Radius(c0+1)    = Rave
         Geom% Radius(c0+2)    = Rmax
         Geom% Area(c0+1)      = zero
         Geom% Area(c0+2)      = zero
         Geom% Volume(c0+1)    = Rave - Rmin
         Geom% Volume(c0+2)    = Rmax - Rave
       enddo

   end select


   return
   end subroutine geometrySCB1D 

