!***********************************************************************
!                      Last Update:  06/2012, PFN                      *
!                                                                      *
!   AngleCoef2D - This routine computes the angle-dependent            *
!                 parameters that result from angular differencing.    *
!                 differencing.  This routine assumes that the         *
!                 angles, m=1,...,NumAngles, are grouped by xi-levels, *
!                 from xi-minimum to xi-maximum, with mu ordered       *
!                 min to max within each xi-level.                     *
!                                                                      *
!   Input:     ndir   - number of angles                               *
!                                                                      *
!   Output:    BETA   - parameters that appear in the discrete         *
!                       ordinates version of the angular derivative    *
!                       term                                           *
!              TAU    - weights that appear in the discrete-ordinates  *
!                       weighted-diamond relation in angle             *
!                                                                      *
!***********************************************************************

   subroutine AngleCoef2D(self)

   use kind_mod
   use constant_mod
   use AngleSet_mod

   implicit none

!  Arguments

   type(AngleSet)          :: self

!  Local

   integer    :: angle, angle1, angle2, level
   integer    :: NumAngles, totalAngles, nLevels, nlevelAngles

   real(adqt) :: weightLevel, Phimh, Phiph, Mumh, Muph

!  First set Alpha:  Alpha is zero for starting and finishing directions

   NumAngles = self% NumAngles

!  Now tau, which is a bit more complicated.  (The complicated
!  part is figuring out what to use for u(m+1/2) and u(m-1/2).)
 
   nLevels     = self% NumBin
   totalAngles = 0 

   XiLevelLoop: do level=1,nLevels

     nlevelAngles = self% NangBinList(level)
     angle1       = totalAngles + 1
     angle2       = totalAngles + nlevelAngles

!  Sum the angle weights for this xi-level
     weightLevel = zero
     Phimh       = pi
     Mumh        = self% omega(1,1)

     do angle=angle1,angle2
       weightLevel = weightLevel + self% weight(angle)
     enddo 

     AngleLoop: do angle=angle1,angle2

       if (self% StartingDirection(angle) ) then

         self% Alpha(angle) = zero
         self% Tau(angle)   = zero

         if ( angle == angle1 ) then
           Phimh = pi
           Mumh  = self% omega(1,angle)
         else
           call f90fatal("Mu not increasing in xi-level, AngleCoef2D")
         endif

       elseif (self% FinishingDirection(angle) ) then

         self% Alpha(angle) = zero
         self% Tau(angle)   = zero

       else

         self% Alpha(angle) = self% Alpha(angle-1) - self% weight(angle)*self% omega(1,angle)

         Phiph = Phimh - self% weight(angle)*pi/weightLevel
         Muph  = sqrt(one - self% omega(2,angle)*self% omega(2,angle))*cos(Phiph)

         if (self% omega(1,angle) < Mumh .or. self% omega(1,angle) > Muph) then
           call f90fatal("Mu not between limits, AngleCoef2D")
         endif

         self% Tau(angle) = (self% omega(1,angle) - Mumh)/(Muph - Mumh)
         Phimh            = Phiph
         Mumh             = Muph

       endif

     enddo AngleLoop

     totalAngles = totalAngles + nlevelAngles

   enddo XiLevelLoop

   if ( totalAngles /= NumAngles ) then
     call f90fatal("Incorrect number of angular coefficients in AngleCoef2D")
   endif
 

 
   return
   end subroutine AngleCoef2D 

