!***********************************************************************
!                        Last Update:  01/2012, PFN                    *
!                                                                      *
!   ConstructGeometry - Construct the geometry module for this         *
!                       spatial domain.                                *
!                                                                      *
!***********************************************************************


   subroutine ConstructGeometry() BIND(C,NAME="teton_constructgeometry")

!  Include
   use ISO_C_BINDING
   use kind_mod
   use Geometry_mod
   use QuadratureList_mod


   implicit none

!  Construct the Geometry Module 

   allocate (Geom)

   call construct(Geom, Quad% nZoneSets)


   return
   end subroutine ConstructGeometry

