!***********************************************************************
!                        Version 1:  08/94, PFN                        *
!                                                                      *
!   QUADSLAB - Calculates group dependent quadrature sets for Sn       *
!              radiation transport in slab geometry.                   *
!                                                                      *
!   Input:   nordr  - vector of group quadrature orders (no st. dir.)  *
!            nangt  - sum of quadrature orders for all groups (w/st d) *
!            ngr    - number of frequency groups                       *
!                                                                      *
!   Output:  OMEGA  - group dependent direction cosines (mu)           *
!            QUADWT - group dependent quadrature weights               *
!                                                                      *
!   Allowed values of "n" are: 2,4,6,8,10,12,14,16,18,20,32,64,128,256 *
!                                                                      *
!***********************************************************************

   subroutine quadslab(self)

   use kind_mod
   use Quadrature_mod
   use QuadratureData_mod

   implicit none

!  Arguments

   type(Quadrature) :: self

!  Local

   integer    :: i, NumAngles, halfAngles, offset

!  Set direction cosines and weights
!  Angular weights sum to 2 for all mu

   NumAngles  = self% NumAngles
   halfAngles = self% order/2
   offset     = iang1D(self% order) - 1
 
!  MU < 0
 
   do i=1,halfAngles
     self% omega(1,i) = -dircos1D(offset+i)
     self% weight(i)  =  weight1D(offset+i)
   enddo
 
!  MU > 0
 
   do i=halfAngles,1,-1
     self% omega(1,NumAngles-i+1) = -self% omega(1,i)
     self% weight(NumAngles-i+1)  =  self% weight(i)
   enddo
 
 
   return
   end subroutine quadslab

