!***********************************************************************
!                        Version 1:  01/93, PFN                        *
!                                                                      *
!   QUADXYZ - Calculates group dependent level-symmetric quadrature    *
!             sets for Sn radiation transport in XYZ geometry.         *
!                                                                      *
!   Input:   norder - quadrature order                                 *
!            nang   - number of angles                                 *
!                                                                      *
!   Output:  OMEGA  - group dependent direction cosines (mu,eta,xi)    *
!            QUADWT - group dependent quadrature weights               *
!                                                                      *
!   Allowed values of "n" are:  2, 4, 6, 8, 10, 12, 14, 16, 18, 20     *
!                                                                      *
!   Directions per Octant:                                             *
!                                                                      *
!                             N   N(N+2)/8                             *
!                             2       1                                *
!                             4       3                                *
!                             6       6                                *
!                             8      10                                *
!                            10      15                                *
!                            12      21                                *
!                            14      28                                *
!                            16      36                                *
!                            18      45                                *
!                            20      55                                *
!                                                                      *
!***********************************************************************

   subroutine quadxyz(self)

   use kind_mod
   use constant_mod
   use Quadrature_mod
   use QuadratureData_mod

   implicit none

!  Arguments

   type(Quadrature) :: self

!  Local

   integer    :: i,ns,nn,Nang,nangoct,jcos,norder

   real(adqt) :: halfpi

   real(adqt) :: amu,aeta,axi,awgt

!  Angular weights sum to pi/2 in each octant 

   halfpi = half*pi

!  Set the direction cosines and weights; note that the
!  angles are numbered consecutively in an octant.
!  NANGOCT is the number of angles per octant and NS is
!  an offset to the first angle and weight for the set.

   norder  = self% order
   nangoct = norder*(norder + 2)/8
   ns      = iang(norder) - 1
   jcos    = icoff(norder)
   nn      = 0

   do i=1,nangoct

     amu  = dircos( imu(ns+i)  + jcos )
     aeta = dircos( ieta(ns+i) + jcos )
     axi  = dircos( ixi(ns+i)  + jcos )
     awgt = weight( iwgt(ns+i) + jcos )

!  Octant 1  mu>0, eta>0, xsi>0

     self% omega(1,nn+1) = amu
     self% omega(2,nn+1) = aeta
     self% omega(3,nn+1) = axi
     self% weight(nn+1)  = halfpi*awgt

!  Octant 2  mu<0, eta>0, xsi>0

     self% omega(1,nn+2)   = -amu
     self% omega(2,nn+2)   =  aeta
     self% omega(3,nn+2)   =  axi
     self% weight(nn+2)    =  halfpi*awgt

!  Octant 3  mu<0, eta<0, xsi>0

     self% omega(1,nn+3) = -amu
     self% omega(2,nn+3) = -aeta
     self% omega(3,nn+3) =  axi
     self% weight(nn+3)  =  halfpi*awgt

!  Octant 4  mu>0, eta<0, xsi>0

     self% omega(1,nn+4) =  amu
     self% omega(2,nn+4) = -aeta
     self% omega(3,nn+4) =  axi
     self% weight(nn+4)  =  halfpi*awgt

!  Octant 5  mu>0, eta>0, xsi<0

     self% omega(1,nn+5) =  amu
     self% omega(2,nn+5) =  aeta
     self% omega(3,nn+5) = -axi
     self% weight(nn+5)  =  halfpi*awgt

!  Octant 6  mu<0, eta>0, xsi<0

     self% omega(1,nn+6) = -amu
     self% omega(2,nn+6) =  aeta
     self% omega(3,nn+6) = -axi
     self% weight(nn+6)  =  halfpi*awgt

!  Octant 7  mu<0, eta<0, xsi<0

     self% omega(1,nn+7) = -amu
     self% omega(2,nn+7) = -aeta
     self% omega(3,nn+7) = -axi
     self% weight(nn+7)  =  halfpi*awgt

!  Octant 8  mu>0, eta<0, xsi<0

     self% omega(1,nn+8) =  amu
     self% omega(2,nn+8) = -aeta
     self% omega(3,nn+8) = -axi
     self% weight(nn+8)  =  halfpi*awgt

     nn = nn + 8

   enddo

!  Set the Polar Angle List

   Nang = norder/2
   nn   = norder + 1

   do i=1,Nang
     self% PolarAngleList(i)    = -dircos( imu(ns+i)  + jcos )
     self% PolarAngleList(nn-i) =  dircos( imu(ns+i)  + jcos ) 
   enddo


   return
   end subroutine quadxyz


