!***********************************************************************
!                        Version 1:  07/2003, PFN                      *
!                                                                      *
!   QUADPRODUCT - Computes a product quadrature set based on the       *
!                 number of polar and azimuthal directions and the     *
!                 specified polar axis.                                *
!                                                                      *
!   Input:   nang       - number of angles                             *
!            npolar     - number of polar angles                       *
!            nazimuthal - number of azimuthal angles                   *
!            polaraxis  - polar axis                                   *
!                                                                      *
!   Output:  OMEGA  - group dependent direction cosines (mu,eta,xi)    *
!            QUADWT - group dependent quadrature weights               *
!                                                                      *
!   Allowed values:  1 <= npolar <= 32, 1 <= nazimuthal <= 32          *
!                                                                      *
!                                                                      *
!***********************************************************************

   subroutine quadProduct(self)

   use kind_mod
   use constant_mod
   use Quadrature_mod
   use QuadratureData_mod

   implicit none

!  Arguments

   type(Quadrature)          :: self

!  Local

   integer    :: i,iPhi,jTheta,m,Phi1,Phi2,Theta1,Theta2,nangoct
   integer    :: nn

   integer    :: npolar,nazimuthal,polaraxis

   real(adqt) :: cosineTheta
   real(adqt) :: sineTheta
   real(adqt) :: cosinePhi

   real(adqt) :: omegaX(self% npolar*self% nazimuthal)
   real(adqt) :: omegaY(self% npolar*self% nazimuthal)
   real(adqt) :: omegaZ(self% npolar*self% nazimuthal)
   real(adqt) :: qweight(self% npolar*self% nazimuthal)

!  Constants

   npolar     = self% npolar
   nazimuthal = self% nazimuthal
   polaraxis  = self% polaraxis

!  Check for Errors

   if (npolar < 1 .or. npolar > 32) then
     call f90fatal("ERROR: npolar must be in range 1 <= npolar <= 32")
   endif

   if (nazimuthal < 1 .or. nazimuthal > 32) then
     call f90fatal("ERROR: nazimuthal must be in range 1 <= nazimuthal <= 32")
   endif

   if (polaraxis < 1 .or. polaraxis > 3) then
     call f90fatal("ERROR: polaraxis must be in range 1 <= polaraxis <= 3")
   endif

!  Set omegaX, omegaY, omegaZ based on the choice of polar axis
!  using tabulated values

   Phi1    = first(nazimuthal) 
   Phi2    =  last(nazimuthal) 
   Theta1  = first(npolar) 
   Theta2  =  last(npolar) 

   nangoct = npolar*nazimuthal

   m = 0
   do iPhi=Phi1,Phi2

     cosinePhi = cosPhiXYZ(iPhi)

!     sinePhi   = sqrt( one - cosinePhi*cosinePhi )

!  Note that we have replaced sineTheta*sinePhi below with a square-root
!  of the other two components. This gives greater accuracy in the
!  magnitude of the ordinate (which should equal one)

     do jTheta=Theta2,Theta1,-1

       cosineTheta = cosTheta(jTheta)
       sineTheta   = sqrt( one - cosineTheta*cosineTheta )
       m           = m + 1

       if (polaraxis == 1) then
         omegaX(m) = cosineTheta
         omegaY(m) = sineTheta*cosinePhi 
         omegaZ(m) = sqrt( one - omegaX(m)*omegaX(m) - omegaY(m)*omegaY(m) )
       elseif (polaraxis == 2) then
         omegaY(m) = cosineTheta
         omegaZ(m) = sineTheta*cosinePhi 
         omegaX(m) = sqrt( one - omegaY(m)*omegaY(m) - omegaZ(m)*omegaZ(m) )
       elseif (polaraxis == 3) then
         omegaZ(m) = cosineTheta
         omegaX(m) = sineTheta*cosinePhi
         omegaY(m) = sqrt( one - omegaX(m)*omegaX(m) - omegaZ(m)*omegaZ(m) )
       endif

       qweight(m) = weightTheta(jTheta)*weightPhiXYZ(iPhi)

     enddo
   enddo

!  Set the direction cosines and weights; note that the
!  angles are numbered consecutively in an octant.
!  NANGOCT is the number of angles per octant.

   nn = 0

   do i=1,nangoct

!  Octant 1  mu>0, eta>0, xsi>0

     self% omega(1,nn+1) =  omegaX(i) 
     self% omega(2,nn+1) =  omegaY(i) 
     self% omega(3,nn+1) =  omegaZ(i) 
     self% weight(nn+1)  =  qweight(i) 

!  Octant 2  mu<0, eta>0, xsi>0

     self% omega(1,nn+2) = -omegaX(i)
     self% omega(2,nn+2) =  omegaY(i) 
     self% omega(3,nn+2) =  omegaZ(i) 
     self% weight(nn+2)  =  qweight(i) 

!  Octant 3  mu<0, eta<0, xsi>0

     self% omega(1,nn+3) = -omegaX(i)
     self% omega(2,nn+3) = -omegaY(i)
     self% omega(3,nn+3) =  omegaZ(i) 
     self% weight(nn+3)  =  qweight(i) 

!  Octant 4  mu>0, eta<0, xsi>0

     self% omega(1,nn+4) =  omegaX(i) 
     self% omega(2,nn+4) = -omegaY(i)
     self% omega(3,nn+4) =  omegaZ(i) 
     self% weight(nn+4)  =  qweight(i) 

!  Octant 5  mu>0, eta>0, xsi<0

     self% omega(1,nn+5) =  omegaX(i) 
     self% omega(2,nn+5) =  omegaY(i) 
     self% omega(3,nn+5) = -omegaZ(i)
     self% weight(nn+5)  =  qweight(i) 

!  Octant 6  mu<0, eta>0, xsi<0

     self% omega(1,nn+6) = -omegaX(i)
     self% omega(2,nn+6) =  omegaY(i) 
     self% omega(3,nn+6) = -omegaZ(i)
     self% weight(nn+6)  =  qweight(i) 

!  Octant 7  mu<0, eta<0, xsi<0

     self% omega(1,nn+7) = -omegaX(i)
     self% omega(2,nn+7) = -omegaY(i)
     self% omega(3,nn+7) = -omegaZ(i)
     self% weight(nn+7)  =  qweight(i) 

!  Octant 8  mu>0, eta<0, xsi<0

     self% omega(1,nn+8) =  omegaX(i) 
     self% omega(2,nn+8) = -omegaY(i)
     self% omega(3,nn+8) = -omegaZ(i)
     self% weight(nn+8)  =  qweight(i) 

     nn = nn + 8

   enddo


   return
   end subroutine quadProduct


