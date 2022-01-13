!***********************************************************************
!                        Version 1:  05/92, PFN                        *
!                                                                      *
!   QUADSP - calculates group dependent quadrature sets for Sn         *
!            radiation transport in spheres                            *
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

   subroutine quadsp(self)

   use kind_mod
   use constant_mod
   use Quadrature_mod
   use QuadratureData_mod

   implicit none

!  Arguments

   type(Quadrature) :: self

!  Local

   integer    :: i, NumAngles, halfAngles, offset

 
!  Set direction cosines and weights
!  Angular weights sum to 2 for all mu (the starting and
!  finishing directions have zero weight)
 
   NumAngles  = self% NumAngles
   halfAngles = self% order/2
   offset     = iang1D(self% order) - 1

!  Starting Direction
 
   self% omega(1,1) = -one
   self% weight(1)  =  zero
 
!  MU < 0

   do i=1,halfAngles
     self% omega(1,i+1) = -dircos1D(offset+i)
     self% weight(i+1)  =  weight1D(offset+i)
   enddo
 
!  MU > 0

   do i=halfAngles,1,-1
     self% omega(1,NumAngles-i) = -self% omega(1,i+1)
     self% weight(NumAngles-i)  =  self% weight(i+1)
   enddo

!  Finishing Direction

   self% omega(1,NumAngles) =  one
   self% weight(NumAngles)  =  zero
 
 
   return
   end subroutine quadsp

