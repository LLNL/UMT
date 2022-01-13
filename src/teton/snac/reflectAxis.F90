!***********************************************************************
!                                                                      *
!                        Last Updated: 09/2017 by PFN                  *
!                                                                      *
!   reflectAxis - This routine, given an incident angle and outward    *
!                 normal of a boundary surface, finds the reflected    *
!                 angle in the set of angles provided.                 *
!                                                                      *
!                 This routine assumes reflections are off of          *
!                 symmetry axes or a 45-degree surface.                *
!                                                                      *
!***********************************************************************

   subroutine reflectAxis(ndim, AngleInc, NumAngles, ReflAngle,  &
                          omega, Area)

   use kind_mod
   use constant_mod

   implicit none

!  Arguments

   integer,    intent(in)    :: ndim
   integer,    intent(in)    :: AngleInc
   integer,    intent(in)    :: NumAngles
   integer,    intent(inout) :: ReflAngle 

   real(adqt), intent(in)    :: omega(ndim,NumAngles)
   real(adqt), intent(in)    :: Area(ndim)

!  Local

   integer    :: i1,i2,i3,ia,d,d2,  &
                 nzero,nmax,mref

   real(adqt) :: domega1,domega2,ratio,AreaMag 

   real(adqt) :: omega_Inc(ndim)

   logical(kind=1) :: arbitrary
 
!  Constants

   real(adqt), parameter :: fuz=1.0e-6_adqt
   real(adqt), parameter :: tol=1.0e-10_adqt

!  Most of the time, we will have reflection in one of the directions
!  of the coordinate system.  The component of OMEGA that changes sign
!  corresponds to the component of the normal to the side (or face)
!  that is non-zero ( e.g. if Area(1).ne.0 then we have reflection in 
!  mu=omega(1)  ).  First we check for the component that changes
!  sign ( omega(d) ), and then we check to see that the other
!  components are unchanged.

!  If this is a 1D problem, find the angle with the opposite sign

   mref = 0
   omega_Inc(:) = omega(:,AngleInc)

   if (ndim == 1) then
     AngleLoop1D: do ia=1,NumAngles
       if (abs(omega(1,ia) + omega_inc(1)) < fuz ) then
         mref = ia
         exit AngleLoop1D
       endif
     enddo AngleLoop1D

     ReflAngle = mref

     return
   endif

!  90 degree reflection - First check to see that all 
!  components but one are "zero"

   AreaMag = zero 

   do d=1,ndim
     AreaMag = AreaMag + Area(d)*Area(d)
   enddo

   nzero = 0
   nmax  = 0
   i1    = 2
   i2    = ndim - 1 
   i3    = 3 

   do d=1,ndim
     ratio = abs( Area(d)*Area(d)/AreaMag )
     if (ratio < tol) then
       nzero = nzero + 1
       i1    = min(i1,d)
       i2    = max(i2,d)
     else
       nmax  = d
     endif
   enddo

   if ( nzero == (ndim - 1) ) then

     mref = 0

     do ia=1,NumAngles
       if (abs(omega(nmax,ia) + omega_inc(nmax)) < fuz ) then

         domega1 = abs(omega(i1,ia) - omega_inc(i1))
         domega2 = abs(omega(i2,ia) - omega_inc(i2))

         if (domega1 < fuz .and. domega2 < fuz) then
           mref = ia 
         endif

       endif
     enddo

     if (mref == 0) then
       call f90fatal("90 degrees; No relected angle found in SNMREF")
     endif

     ReflAngle = mref

   else 
 
!  If we make it here, the next easiest case is a 45 degree angle

     arbitrary = .TRUE. 
     do d=1,ndim
       d2    = mod(ndim+d,ndim) + 1
       ratio = abs( (Area(d) + tol)/(Area(d2) + tol) )
       if (ratio < one+fuz .and. ratio > one-fuz) then
         i1        = d
         i2        = d2
         i3        = mod(ndim+i2,ndim) + 1
         arbitrary = .FALSE. 
       endif
     enddo

     if ( .not. arbitrary ) then
 
       mref = 0
       do ia=1,NumAngles
         if (abs(omega(i1,ia) - omega_inc(i2)) < fuz ) then
           if (abs(omega(i2,ia) - omega_inc(i1)) < fuz ) then
             if (abs(omega(i3,ia) - omega_inc(i3)) < fuz ) then
               mref = ia 
             endif
           endif
         endif
       enddo
 
       if (mref == 0) then
         call f90fatal("45 degrees; No reflected angle found in SNMREF")
       endif

       ReflAngle = mref
 
     endif

   endif

 

   return
   end subroutine reflectAxis 

