!***********************************************************************
!                        Last Update:  03/2012, PFN                    *
!                                                                      *
!   SNMOMENTS - This routine, called by SNFLWXYZ                       *
!               calculates the required spherical harmonic moments     *
!               [phi] of the angular flux [psi]. It uses the array     *
!               ynm(n,m), whose definition is:                         * 
!                                                                      *
!               ynm(n,m) = real part of (l,k)th spherical harmonic,    *
!                          evaluated at the mth direction, where       *
!                                                                      *
!                             n = 1 + l*(l+1)/2 + k                    *
!                                                                      *
!   Input:   psi      - angle-dependent intensity      (E/A/t/ster)    *
!            quadwt   - quadrature weights                      (0)    *
!                                                                      *
!   Output:  PHI      - scalar intensity                    (E/A/t)    *
!                                                                      *
!***********************************************************************

   subroutine snmoments(setID)

   use kind_mod
   use constant_mod
   use Geometry_mod
   use Quadrature_mod
   use QuadratureList_mod
   use Size_mod
   use SetData_mod
   use AngleSet_mod

   implicit none

!  Arguments

   integer,  intent(in)  :: setID

!  Local

   type(SetData),          pointer  :: Set
   type(AngleSet),         pointer  :: ASet

   integer    :: Angle
   integer    :: zone
   integer    :: nZones
   integer    :: nCorner
   integer    :: c
   integer    :: c0
   integer    :: NumAngles 

   real(adqt) :: quadwt0 
   real(adqt) :: quadwt1 
   real(adqt) :: quadwt2 
   real(adqt) :: quadwt3 

!  Constants

   Set       => getSetData(Quad, setID)
   ASet      => getAngleSetFromSetID(Quad, setID)
   nZones    =  Size% nZones 
   NumAngles =  Set% NumAngles

!  Add this angles contribution to the flux moments

   AngleTest: if ( mod(NumAngles,4) == 0 ) then

     ZoneLoop4: do zone=1,nZones

       nCorner = Geom% numCorner(zone) 
       c0      = Geom% cOffSet(zone) 

       do c=1,nCorner
         Set% Phi(:,c0+c) = zero
       enddo

       AngleLoop4: do Angle=1,NumAngles,4

         quadwt0 = ASet% Weight(Angle)
         quadwt1 = ASet% Weight(Angle+1)
         quadwt2 = ASet% Weight(Angle+2)
         quadwt3 = ASet% Weight(Angle+3)

         do c=1,nCorner
           Set% Phi(:,c0+c) = Set% Phi(:,c0+c)                 +  &
                              quadwt0*Set% Psi(:,c0+c,Angle)   +  &
                              quadwt1*Set% Psi(:,c0+c,Angle+1) +  &
                              quadwt2*Set% Psi(:,c0+c,Angle+2) +  &
                              quadwt3*Set% Psi(:,c0+c,Angle+3)
         enddo

       enddo AngleLoop4

     enddo ZoneLoop4

   elseif (mod(NumAngles,3) == 0 ) then

     ZoneLoop3: do zone=1,nZones

       nCorner = Geom% numCorner(zone) 
       c0      = Geom% cOffSet(zone) 

       do c=1,nCorner
         Set% Phi(:,c0+c) = zero
       enddo

       AngleLoop3: do Angle=1,NumAngles,3

         quadwt0 = ASet% Weight(Angle)
         quadwt1 = ASet% Weight(Angle+1)
         quadwt2 = ASet% Weight(Angle+2)

         do c=1,nCorner
           Set% Phi(:,c0+c) = Set% Phi(:,c0+c)                 +  &
                              quadwt0*Set% Psi(:,c0+c,Angle)   +  &
                              quadwt1*Set% Psi(:,c0+c,Angle+1) +  &
                              quadwt2*Set% Psi(:,c0+c,Angle+2)
         enddo

       enddo AngleLoop3

     enddo ZoneLoop3

   elseif (mod(NumAngles,2) == 0 ) then

     ZoneLoop2: do zone=1,nZones

       nCorner = Geom% numCorner(zone) 
       c0      = Geom% cOffSet(zone) 

       do c=1,nCorner
         Set% Phi(:,c0+c) = zero
       enddo

       AngleLoop2: do Angle=1,NumAngles,2

         quadwt0 = ASet% Weight(Angle)
         quadwt1 = ASet% Weight(Angle+1)

         do c=1,nCorner
           Set% Phi(:,c0+c) = Set% Phi(:,c0+c)               +  &
                              quadwt0*Set% Psi(:,c0+c,Angle) +  &
                              quadwt1*Set% Psi(:,c0+c,Angle+1)
         enddo

       enddo AngleLoop2

     enddo ZoneLoop2

   else

     ZoneLoop1: do zone=1,nZones

       nCorner = Geom% numCorner(zone) 
       c0      = Geom% cOffSet(zone) 

       do c=1,nCorner
         Set% Phi(:,c0+c) = zero
       enddo

       AngleLoop1: do Angle=1,NumAngles

         quadwt0 = ASet% Weight(Angle)

         do c=1,nCorner
           Set% Phi(:,c0+c) = Set% Phi(:,c0+c)               +  &
                              quadwt0*Set% Psi(:,c0+c,Angle)
         enddo

       enddo AngleLoop1

     enddo ZoneLoop1

   endif AngleTest


   return
   end subroutine snmoments


