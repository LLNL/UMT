!***********************************************************************
!                        Version 1:  05/92, PFN                        *
!                                                                      *
!   QUADCY - calculates group dependent quadrature sets for Sn         *
!            radiation transport in cylinders                          *
!                                                                      *
!   Input:   nordr  - vector of group quadrature orders (no st. dir.)  *
!            nangt  - sum of quadrature orders for all groups (w/st d) *
!            ngr    - number of frequency groups                       *
!                                                                      *
!   Output:  OMEGA  - group dependent direction cosines (mu)           *
!            QUADWT - group dependent quadrature weights               *
!                                                                      *
!   Allowed values of "n" are:  2, 4, 6, 8, 12, 16                     *
!                                                                      *
!***********************************************************************

   subroutine quadcy(self)

   use kind_mod
   use constant_mod
   use Quadrature_mod
   use QuadratureData_mod

   implicit none

!  Arguments

   type(Quadrature) :: self

!  Local

   integer    :: i, ns, angle, offset 
   integer    :: level, nLevels, nAngLevel
   real(adqt) :: xilev


!  Set direction cosines and weights
!  Angular weights sum to 2 for all mu (starting and finishing
!  directions have zero weight)
 
   nLevels   = self% order/2
   nAngLevel = self% order/2
   ns        = iang(self% order) - 1
   offset    = icoff(self% order)
   angle     = 0
 
!  Loop over XI levels
 
   do level=1,nLevels
 
     xilev = dircos( offset + ieta(ns+1) )
 
!  Starting Direction
 
     angle                =  angle + 1
     self% omega(1,angle) = -sqrt(one - xilev*xilev)
     self% weight(angle)  =  zero
 
!  MU < 0
 
     do i=1,nAngLevel
       angle                =  angle + 1 
       self% omega(1,angle) = -dircos( offset + imu(ns+i) )
       self% weight(angle)  =  weight( offset + imu(ns+i) )
     enddo
 
!  MU > 0
 
     do i=nAngLevel,1,-1
       angle                =  angle + 1
       self% omega(1,angle) =  dircos( offset + imu(ns+i) )
       self% weight(angle)  =  weight( offset + imu(ns+i) )
     enddo

!  Finishing Direction

     angle                = angle + 1
     self% omega(1,angle) = sqrt(one - xilev*xilev)
     self% weight(angle)  = zero
 
     ns = ns + nAngLevel
     nAngLevel = nAngLevel - 1
 
   enddo
 
 
   return
   end subroutine quadcy

