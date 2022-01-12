!***********************************************************************
!                        Version 1:  05/92, PFN                        *
!                                                                      *
!   AngleCoef1D - calculates angular derivative coefficients for       *
!                 1D spheres and cylinders                             *
!                                                                      *
!   Input:   norder - quadrature order for group ig                    *
!            nangt  - sum of quadrature orders for all groups          *
!                                                                      *
!   Output:  ADWETA - alpha for angular derivative divided by weight   *
!            ETA    - angular weights for angular derivative           *
!            FALPH  - function of alpha and eta                        *
!                                                                      *
!   Local:   halfmu - half-integer mu                                  *
!            alphmu - alpha for angular derivative                     *
!                                                                      *
!   Equations:  mu(m+)    = mu(m-) + w(m)                              *
!               alpha(m+) = alpha(m-) - 2*w(m)*mu(m)                   *
!               eta(m)    = (mu(m+) - mu(m-))/(mu(m) - mu(m-))         *
!               falph(m)  = (alpha(m-)+alpha(m+)*(1-eta(m)))/w(m)      *
!                                                                      *
!***********************************************************************
   subroutine AngleCoef1D(self) 

   use flags_mod
   use kind_mod
   use constant_mod
   use Size_mod
   use AngleSet_mod

   implicit none 

!  Arguments

   type(AngleSet)         :: self

!  Local Variables

   integer    :: i,ixi,iang1,norder,nxilevels,nangxi

   real(adqt) :: phimh,phiph,wtot

   real(adqt) :: alphmu(self% NumAngles), halfmu(self% NumAngles)

!  Initialize output arrays

   self% Adweta(:) = zero
   self% Falpha(:) = zero
   self% Tau(:)    = zero
   nxilevels       = 0
 
!  If this is slab geometry, return

   if (Size% igeom == geometry_slab) then

     return

   else 
 
!  Outer loop is over xi-levels 
 
     iang1 = 1
 
     norder = self% Order 

     if (Size% igeom == geometry_sphere) then
       nxilevels = 1
     elseif (Size% igeom == geometry_cylinder) then
       nxilevels = norder/2
     endif

     XiLevelLoop: do ixi=1,nxilevels
 
       nangxi = norder - 2*(ixi-1)
 
!  Calculate angular coefficients for the current XI-level
 
       halfmu(1) = self% omega(1,iang1)
       alphmu(1) = zero
 
!  Cylinders
 
       if (Size% igeom == geometry_cylinder) then
 
         wtot = zero
 
         do i=1,nangxi
           wtot = wtot + self% weight(iang1+i)
         enddo
 
         phimh = pi
 
         do i=1,nangxi
           phiph       =  phimh + self% weight(iang1+i)*pi/wtot
           halfmu(i+1) = -self% omega(1,iang1)*cos(phiph)
           alphmu(i+1) =  alphmu(i) - self% weight(iang1+i)*self% omega(1,iang1+i)
           phimh       =  phiph
         enddo
 
!  Spheres
 
       elseif (Size% igeom == geometry_sphere) then
 
         do i=1,nangxi
           halfmu(i+1) = halfmu(i) + self% weight(iang1+i)
           alphmu(i+1) = alphmu(i) - two*self% weight(iang1+i)*self% omega(1,iang1+i)
         enddo
 
       endif
 
!  Really need alpha*eta divided by weight
 
       do i=1,nangxi
         self% Tau(iang1+i)     = (halfmu(i+1) - halfmu(i))/   &
                                  (self% omega(1,iang1+i) - halfmu(i))
         self% Falpha(iang1+i)  = (alphmu(i) + (self% Tau(iang1+i) - one)*  &
                                   alphmu(i+1))/self% weight(iang1+i)
         self% Adweta(iang1+i)  =  alphmu(i+1)*self% Tau(iang1+i)/  &
                                   self% weight(iang1+i)
       enddo
 
       iang1 = iang1 + nangxi + 1
 
     enddo XiLevelLoop
 
   endif 


   return
   end subroutine AngleCoef1D 


