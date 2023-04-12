!***********************************************************************
!                        Version 1:  05/92, PFN                        *
!                                                                      *
!   RTQUAD - Calls routines to calculate group dependent quadrature    *
!            sets for Sn radiation transport depending on geometry.    *
!                                                                      *
!   Input:   nordr  - group quadrature orders                          *
!            nangr  - total number of angles by group                  *
!            nangt  - sum of angles for all groups (with st+fin dir)   *
!            igeom  - geometry flag                                    *
!                                                                      *
!   Output:  OMEGA  - group dependent direction cosines (mu,eta,xi)    *
!            QUADWT - group dependent quadrature weights               *
!                                                                      *
!   Allowed values of "n" are:  1D: 2,4,6,8,10,12,14,16,18,20,32,64,128,256 *
!                               2D:  2, 4, 6, 8, 12, 16                *
!                               3D:  2, 4, 6, 8, 12, 16, 18, 20        *
!                                                                      *
!***********************************************************************
   subroutine rtquad(self)

   use flags_mod
   use kind_mod
   use constant_mod
   use Size_mod
   use Quadrature_mod

   implicit none

!  Arguments

   type(Quadrature)             :: self 

!  Local

   integer          :: ia, iang, NumAngles
   integer          :: ndim

   real(adqt)       :: fac, sumwgt, wtiso

   character(len=8) :: TypeName
   integer          :: igeom

!  Select the appropriate quadrature based on geometry 

   NumAngles = self% NumAngles
   TypeName  = self% TypeName
   ndim      = Size% ndim
   wtiso     = Size% wtiso
   igeom     = Size% igeom


   select case (igeom)

!  Slabs

     case (geometry_slab)

       call quadslab(self)

!  Cylinders

     case (geometry_cylinder)

       call quadcy(self)

!  Spheres

     case (geometry_sphere)

       call quadsp(self)

!  2-D RZ

     case (geometry_rz)

       call quadrz(self)

!  3-D XYZ

     case (geometry_xyz)

       if (TypeName == 'levelsym') then
         call quadxyz(self)
       elseif (TypeName == 'product') then
         call quadProduct(self)
       elseif (TypeName == 'lobatto') then
         call quadLobatto(self)
       endif

   end select
 
!  Make sure that quadrature weights are normalized correctly

   sumwgt = zero

   do ia=1,NumAngles
     sumwgt = sumwgt + self% weight(ia)
   enddo

   fac = one/(wtiso*sumwgt)

   do ia=1,NumAngles
     self% weight(ia) = fac*self% weight(ia)
   enddo

!  Identify starting and finishing directions

   self% StartingDirection(:)  = .FALSE.
   self% FinishingDirection(:) = .FALSE. 

   iang = -1
   do ia=1,NumAngles
     if (self% weight(ia) < adqtEpsilon) then
       if (iang == -1) then
         self% StartingDirection(ia) = .TRUE.
         iang = -iang
       elseif (iang == 1) then
         self% FinishingDirection(ia) = .TRUE.
         iang = -iang
       endif
     endif
   enddo

 
   return
   end subroutine rtquad


