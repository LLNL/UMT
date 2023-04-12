#include "macros.h"
!***********************************************************************
!                        Last Update:  11/2022, PFN                    *
!                                                                      *
!  ConstructQuadrature  - Wrapper for module that can be called from   * 
!                         C++ to create the QuadratureList type.       *
!                                                                      *
!            QuadDef    - Quadrature Definition                        *
!                                                                      *
!                         2nd index = 1, high-order Sn quadrature      *
!                         2nd index = 2, GTA quadrature                *
!                                                                      *
!                         1 = type                                     *
!                         2 = order                                    *
!                         3 = polar angles                             *
!                         4 = azimuthal angles                         *
!                         5 = polar axis                               *
!                         6 = number of angles in the set (Output)     *
!                                                                      *
!                         "type" definitions:                          *
!                            1 = level symmetric                       *
!                            2 = product quadrature                    *
!                            3 = Lobatto (3D only)                     *
!                                                                      *
!***********************************************************************


   subroutine ConstructQuadrature(nSetsMaster, nSets, QuadDef, gnu) &
        BIND(C,NAME="teton_constructquadrature_new")

!  Include

   USE ISO_C_BINDING
   use flags_mod
   use kind_mod
   use Size_mod
   use QuadratureList_mod


   implicit none

!  Arguments

   integer(C_INT),    intent(in)    :: nSetsMaster
!   Warning: nSets may not be the value you think it should be.  See comments
!     in the QuadratureList type definition
   integer(C_INT),    intent(inout) :: nSets
   integer(C_INT),    intent(inout) :: QuadDef(6,2)
   real(C_DOUBLE),    intent(in)    :: gnu(Size%ngr+1)

!  Local

   integer :: quadID 
   integer :: quadType
   integer :: norder
   integer :: nAngles
   integer :: npolar
   integer :: nazimuthal
   integer :: nAnglesSn 
   integer :: np
   integer :: na

   integer :: nordermax

!  Construct the Quadrature Module 

   if ( .not. associated(Quad) ) then
     allocate (Quad)
   endif

!  The following loop only determines the number of angles in
!  each set based on type, order and geometry 

   do quadID=1,2

     quadType = QuadDef(1,quadID)
     norder   = QuadDef(2,quadID)
     nAngles  = 0

     ! Checks on nordermax:
     if ( quadType == 1 ) then
       if (Size%igeom == geometry_rz) then
         nordermax = 16
       else if (Size%ndim > 1) then
         nordermax = 20
       else
         nordermax = 256
       endif

       if (norder > nordermax) then

         print *, "WARNING: Quadrature order must be an even number <= ", nordermax, " for level symmetric Teton runs in this geometry.  Teton is reducing your requested quadrature order ", norder, " to ", nordermax
         norder = nordermax
         QuadDef(2,quadID) = norder

       else if (mod(norder,2) /= 0) then

         print *, "WARNING: Quadrature order must be an even number <= ", nordermax, " for level symmetric Teton runs in this geometry.  Teton is changing your requested quadrature order ", norder, " to ", max(2,norder-1)
         norder = max(2,norder-1)
         QuadDef(2,quadID) = norder

       endif
     endif

!  The number of angles depend on geometry and quadrature type

     select case (Size% igeom)

       case (geometry_cylinder)
         nAngles = norder*(norder + 4)/4

       case (geometry_sphere)
         nAngles = norder + 2

       case (geometry_slab)
         nAngles = norder

       case (geometry_rz)
         quadType   = QuadDef(1,quadID)
         npolar     = QuadDef(3,quadID)
         nazimuthal = QuadDef(4,quadID)

         if (quadType == 1) then
           nAngles = norder*(norder + 6)/2
         elseif (quadType == 2) then
           if (nazimuthal > 0) then
             nAngles = 4*npolar*(nazimuthal + 1)
           else
             nAngles = 0
             do np=1,npolar
               na = min(np - 1 + max(abs(nazimuthal),1), 32)
               ! Two symmetric polar angles, times two quadrants, plus two
               ! start+stop angles per polar level
               nAngles = nAngles + 2*(2*na + 2)
             enddo
           endif
         else
           call f90fatal("Invalid quadrature definition in Construct Quadrature")
         endif

       case (geometry_xyz)
         quadType   = QuadDef(1,quadID)
         npolar     = QuadDef(3,quadID)
         nazimuthal = QuadDef(4,quadID)

         if (quadType == 1) then
           nAngles = norder*(norder + 2)
         elseif (quadType == 2) then
           if(nazimuthal>0) then
             nAngles = 8*npolar*nazimuthal
           else
             nAngles = 0
             do np=1,npolar
               na = min(np - 1 + max(abs(nazimuthal),1), 30)
               ! Two symmetric polar angles, times 4 quadrants
               nAngles = nAngles + 8*na
             enddo
           endif
         elseif (quadType == 3) then
           nAngles = 8*npolar*nazimuthal + 2
         else
           call f90fatal("Invalid quadrature definition in Construct Quadrature")
         endif

       case default
         call f90fatal("Invalid geometry in Construct Quadrature")

     end select

     QuadDef(6,quadID) = nAngles

   enddo 

   nAnglesSn = QuadDef(6,1)

!  Construct the QuadratureList object:

   call construct(Quad, nAnglesSn, nSetsMaster, nSets) 

!  Set the high-order and GTA angle sets

   call constructAngleSets(QuadDef, gnu)


   return
   end subroutine ConstructQuadrature

