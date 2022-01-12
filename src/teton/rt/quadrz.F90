!***********************************************************************
!                        Version 1:  07/97, PFN                        *
!                                                                      *
!   QUADRZ - Calculates group dependent quadrature sets for Sn         *
!            radiation transport in RZ geometry.                       *
!                                                                      *
!   Input:   nordr  - vector of group quadrature orders                *
!            nangt  - sum of quadrature orders for all groups          *
!            ngr    - number of frequency groups                       *
!                                                                      *
!   Output:  OMEGA  - group dependent direction cosines (mu,xi)        *
!            QUADWT - group dependent quadrature weights               *
!                                                                      *
!   Allowed values of "n" are:  2, 4, 6, 8, 12, 16                     *
!                                                                      *
!   Directions per quadrant:                                           *
!                                                                      *
!                             N   N(N+6)/8                             *
!                             2       2                                *
!                             4       5                                *
!                             6       9                                *
!                             8      14                                *
!                            12      27                                *
!                            16      44                                *
!                                                                      *
!   Angles are ordered by xi-level for parallelism:                    *
!                                                                      *
!  e.g.  S4                                                            *
!                                                                      *
!                                   xi+                                *
!                                  |                                   *
!                           17 18  |  19 20       level 4              *
!                                  |                                   *
!                                  |                                   *
!                      7  8     9  | 10    11 12  level 2              *
!                       -----------------------mu+                     *
!                      1  2     3  |  4     5  6  level 1              *
!                                  |                                   *
!                                  |                                   *
!                           13 14  | 15  16       level 3              *
!                                  |                                   *
!                                                                      *
!***********************************************************************

   subroutine quadrz(self)

   use kind_mod
   use constant_mod
   use Quadrature_mod
   use QuadratureData_mod

   implicit none

!  Arguments

   type(Quadrature) :: self

!  Local

   integer    :: i,level,ns,m,jcos,nLevels,nangLevel,norder
   integer    :: iPhi,jTheta,Phi1,Phi2,Theta1,Theta2,nazlocal

   real(adqt) :: halfpi
   real(adqt) :: xilev
   real(adqt) :: cosineTheta
   real(adqt) :: sineTheta


!  Angular weights sum to pi/2 in each quadrant

   halfpi = half*pi

!  Angles are collected by xi-level and the
!  xi-levels are ordered from longest to shortest.  The xi-levels
!  come in pairs (i.e. xi<0 and xi>0 for a given set of mu) and
!  within a pair the xi<0 level comes first (this allows us to
!  easily treat a reflecting boundary on the z-axis).  Within a
!  xi-level, angles are ordered in increasing value of mu.


   if (self% TypeName == 'levelsym') then

     norder    = self% order
     nLevels   = norder/2
     nangLevel = norder/2
     ns        = iang(norder) - 1
     jcos      = icoff(norder)
     m         = 0

!  Loop over pairs of xi-levels

     do level=1,nLevels

       xilev = dircos( ieta(ns+1) + jcos )

!    Starting direction

       m                =  m + 1
       self% omega(1,m) = -sqrt(one - xilev*xilev)
       self% omega(2,m) = -xilev
       self% weight(m)  =  zero

!    Quadrant 3  mu<0, xsi<0

       do i=1,nangLevel
         m                =  m + 1
         self% omega(1,m) = -dircos( imu(ns+i) + jcos )
         self% omega(2,m) = -xilev
         self% weight(m)  =  halfpi*weight( iwgt(ns+i) + jcos )
       enddo

!    Quadrant 4  mu>0, xsi<0

       do i=nangLevel,1,-1
         m                =  m + 1
         self% omega(1,m) =  dircos( imu(ns+i) + jcos )
         self% omega(2,m) = -xilev
         self% weight(m)  =  halfpi*weight( iwgt(ns+i) + jcos )
       enddo

!    Finishing direction

       m                =  m + 1
       self% omega(1,m) =  sqrt(one - xilev*xilev)
       self% omega(2,m) = -xilev
       self% weight(m)  =  zero

!    Starting direction

       m                =  m + 1
       self% omega(1,m) = -sqrt(one - xilev*xilev)
       self% omega(2,m) =  xilev
       self% weight(m)  =  zero

!    Quadrant 2  mu<0, xsi>0

       do i=1,nangLevel
         m                =  m + 1
         self% omega(1,m) = -dircos( imu(ns+i) + jcos )
         self% omega(2,m) =  xilev
         self% weight(m)  =  halfpi*weight( iwgt(ns+i) + jcos )
       enddo

!    Quadrant 1  mu>0, xsi>0

       do i=nangLevel,1,-1
         m                =  m + 1
         self% omega(1,m) =  dircos( imu(ns+i) + jcos )
         self% omega(2,m) =  xilev
         self% weight(m)  =  halfpi*weight( iwgt(ns+i) + jcos )
       enddo

!    Finishing direction

       m                =  m + 1
       self% omega(1,m) =  sqrt(one - xilev*xilev)
       self% omega(2,m) =  xilev
       self% weight(m)  =  zero

       ns        = ns + nangLevel
       nangLevel = nangLevel - 1

     enddo

   elseif (self% TypeName == 'product') then

!  Check for Errors

     if (self% npolar < 1 .or. self% npolar > 32) then
       call f90fatal("ERROR: npolar must be in range 1 <= npolar <= 32")
     endif

     if (self% nazimuthal < 1 .or. self% nazimuthal > 32) then
       call f90fatal("ERROR: nazimuthal must be in range 1 <= nazimuthal <= 32")
     endif
! Disable triangle sets for now, something breaks in AngleCoef2D
!     if ( abs(self% nazimuthal) > 10) then
!       call f90fatal("ERROR: nazimuthal  must be <= 10")
!     endif

     Theta1 = first(self% npolar)
     Theta2 =  last(self% npolar)
     m      = 0

     do jTheta=Theta2,Theta1,-1

       cosineTheta = cosTheta(jTheta)
       sineTheta   = sqrt( one - cosineTheta*cosineTheta )
       xilev       = cosineTheta

!  Starting direction

       m                =  m + 1
       self% omega(1,m) = -sqrt(one - xilev*xilev)
       self% omega(2,m) = -xilev
       self% weight(m)  =  zero

       if( self% nazimuthal > 0 ) then
         nazlocal = self% nazimuthal
       else
         nazlocal = min(jTheta - Theta1 + max(abs(self% nazimuthal),1), 10)
       endif

       Phi1   = first(nazlocal)
       Phi2   =  last(nazlocal)

!  Octant 3  mu<0, xsi<0

       do iPhi=Phi2,Phi1,-1
         m                =  m + 1
         self% omega(1,m) = -sineTheta*cosPhiRZ(iPhi)
         self% omega(2,m) = -xilev
         self% weight(m)  =  weightTheta(jTheta)*weightPhiRZ(iPhi)
       enddo

!  Octant 4  mu>0, xsi<0

       do iPhi=Phi1,Phi2
         m                =  m + 1
         self% omega(1,m) =  sineTheta*cosPhiRZ(iPhi)
         self% omega(2,m) = -xilev
         self% weight(m)  =  weightTheta(jTheta)*weightPhiRZ(iPhi)
       enddo

!  Finishing direction

       m                =  m + 1
       self% omega(1,m) =  sqrt(one - xilev*xilev)
       self% omega(2,m) = -xilev
       self% weight(m)  =  zero

!    Starting direction

       m                =  m + 1
       self% omega(1,m) = -sqrt(one - xilev*xilev)
       self% omega(2,m) =  xilev
       self% weight(m)  =  zero

!  Octant 2  mu<0, xsi>0

       do iPhi=Phi2,Phi1,-1
         m                = m + 1
         self% omega(1,m) = -sineTheta*cosPhiRZ(iPhi)
         self% omega(2,m) =  xilev
         self% weight(m)  =  weightTheta(jTheta)*weightPhiRZ(iPhi)
       enddo

!  Octant 1  mu>0, xsi>0

       do iPhi=Phi1,Phi2
         m                = m + 1
         self% omega(1,m) = sineTheta*cosPhiRZ(iPhi)
         self% omega(2,m) = xilev
         self% weight(m)  = weightTheta(jTheta)*weightPhiRZ(iPhi)
       enddo

!    Finishing direction

       m                =  m + 1
       self% omega(1,m) =  sqrt(one - xilev*xilev)
       self% omega(2,m) =  xilev
       self% weight(m)  =  zero

     enddo

   endif


   return
   end subroutine quadrz


