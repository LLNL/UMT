!***********************************************************************
!                        Last Update:  7/2017, PGM                     *
!                                                                      *
!  ConstructQuadrature  - Wrapper for module that can be called from   * 
!                         C++ to create the QuadratureList type.       *
!                                                                      *
!            QuadDef     - quadrature definition by group;             *
!                          entry "ngr+1" is for acceleration           *
!                                                                      *
!                          1 = type                                    *
!                          2 = order                                   *
!                          3 = polar angles                            *
!                          4 = azimuthal angles                        *
!                          5 = polar axis                              *
!                          6 = number of angles in the set (Output)    *
!                                                                      *
!                          "type" definitions:                         *
!                             1 = level symmetric                      *
!                             2 = product quadrature                   *
!                             3 = Lobatto (3D only)                    *
!                                                                      *
!***********************************************************************


   subroutine ConstructQuadrature(nSetsMaster, nSets, QuadDef, gnu) &
        BIND(C,NAME="teton_constructquadrature")

!  Include

   USE ISO_C_BINDING
   use flags_mod
   use kind_mod
   use Size_mod
   use QuadratureList_mod
   use Quadrature_mod


   implicit none

!  Arguments

   integer(C_INT),    intent(in)    :: nSetsMaster
!   Warning: nSets may not be the value you think it should be.  See comments
!     in the QuadratureList type definition
   integer(C_INT),    intent(inout) :: nSets
   integer(C_INT),    intent(inout) :: QuadDef(6,Size%ngr+1)
   real(C_DOUBLE),     intent(in)   :: gnu(Size%ngr+1)

!  Local

   integer :: g 
   integer :: NumQuadSets
   integer :: quadType
   integer :: type_set
   integer :: norder
   integer :: norder_set
   integer :: nAngles
   integer :: npolar
   integer :: nazimuthal
   integer :: set
   integer :: NumSnSets
   integer :: totalAngles
   integer :: np
   integer :: na

   integer :: nordermax

!  Construct the Quadrature Module 

   if ( .not. associated(Quad) ) then
     allocate (Quad)
   endif

!  Find the number of angle sets 

   NumQuadSets =  0 
   totalAngles =  0
   type_set    = -1 
   norder_set  = -1 

   GroupLoop: do g=1,Size%ngr+1

     quadType = QuadDef(1,g)
     norder   = QuadDef(2,g)
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
         QuadDef(2,g) = norder

       else if (mod(norder,2) /= 0) then

         print *, "WARNING: Quadrature order must be an even number <= ", nordermax, " for level symmetric Teton runs in this geometry.  Teton is changing your requested quadrature order ", norder, " to ", max(2,norder-1)
         norder = max(2,norder-1)
         QuadDef(2,g) = norder

       endif
     endif

    
     if ( quadType == type_set   .and.   &
            norder == norder_set .and.   &
            g < Size%ngr+1  ) then

!  This group belongs to the current batch
!  Last "group" is the grey set and should always be its own set

     else

       NumQuadSets = NumQuadSets + 1
       type_set    = quadType
       norder_set  = norder

     endif

!  Quadratures depend on geometry

     select case (Size% igeom)

       case (geometry_cylinder)
         nAngles = norder*(norder + 4)/4

       case (geometry_sphere)
         nAngles = norder + 2

       case (geometry_slab)
         nAngles = norder

       case (geometry_rz)
         quadtype   = QuadDef(1,g)
         npolar     = QuadDef(3,g)
         nazimuthal = QuadDef(4,g)

         if (quadtype == 1) then
           nAngles = norder*(norder + 6)/2
         elseif (quadtype == 2) then
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
         quadtype   = QuadDef(1,g)
         npolar     = QuadDef(3,g)
         nazimuthal = QuadDef(4,g)

         if (quadtype == 1) then
           nAngles = norder*(norder + 2)
         elseif (quadtype == 2) then
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
         elseif (quadtype == 3) then
           nAngles = 8*npolar*nazimuthal + 2
         else
           call f90fatal("Invalid quadrature definition in Construct Quadrature")
         endif

       case default
         call f90fatal("Invalid geometry in Construct Quadrature")

     end select

     QuadDef(6,g) = nAngles
     totalAngles  = totalAngles + nAngles 

   enddo GroupLoop

!  The last quadrature set is for GTA 

   NumSnSets = NumQuadSets - 1

   call construct(Quad, NumQuadSets, Size%ngr, totalAngles, &
                  nSetsMaster, nSets)

!  Set the unique angle sets

   call constructAngleSets(QuadDef, gnu)

   return
   end subroutine ConstructQuadrature

