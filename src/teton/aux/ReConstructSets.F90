!***********************************************************************
!                        Last Update:  10/2016, PFN                    *
!                                                                      *
!   ReConstructSets - Reconstructs the angle-dependent and angle-      *
!                     integrated intensities after remap.              *
!                                                                      *
!***********************************************************************
 
   subroutine ReConstructSets(setID, Groups, RadEnergyDensity)

   use kind_mod
   use constant_mod
   use radconstant_mod
   use Size_mod
   use QuadratureList_mod
   use Geometry_mod
   use SetData_mod

   implicit none 

!  Arguments

   integer,    intent(in)    :: setID
   integer,    intent(in)    :: Groups
   real(adqt), intent(in)    :: RadEnergyDensity(Size% nzones,Size% ngr)

!  Local

   type(SetData),   pointer  :: Set

   integer    :: angle 
   integer    :: g 
   integer    :: g0
   integer    :: zone 
   integer    :: NumAngles 
   integer    :: nZones 
   integer    :: nCorner 
   integer    :: c 
   integer    :: c0

   real(adqt) :: sumPhi(Groups) 
   real(adqt) :: factor(Groups) 

!  Parameters 

   real(adqt), parameter :: Phifloor = 1.0e-40_adqt
   real(adqt), parameter :: facMin   = 0.01_adqt
   real(adqt), parameter :: facMax   = 100.0_adqt


   Set  => getSetData(Quad, setID)

   nZones    = Size% nzones
   NumAngles = Set% NumAngles
   g0        = Set% g0

!***********************************************************************
!  Adjust the discrete angular intensities to conserve energy          *
!***********************************************************************

   ZoneLoop: do zone=1,nZones

     nCorner = Geom% numCorner(zone)
     c0      = Geom% cOffSet(zone)

!  Compute multiplicative correction

     sumPhi(:) = zero

     do c=1,nCorner
       do g=1,Groups
         sumPhi(g) = sumPhi(g) + Geom% Volume(c0+c)*Geom% PhiTotal(g0+g,c0+c)
       enddo
     enddo

     do g=1,Groups
       if (sumPhi(g) > Phifloor) then
         factor(g) = speed_light*RadEnergyDensity(zone,g0+g)*Geom%VolumeZone(zone)/sumPhi(g)
         factor(g) = min( factor(g), facMax )
         factor(g) = max( factor(g), facMin )
       else
         factor(g) =  one
       endif

     enddo

!  Apply correction uniformly 

     do c=1,nCorner

       do angle=1,NumAngles
         Set% Psi(:,c0+c,angle) = factor(:)*Set% Psi(:,c0+c,angle)
       enddo

     enddo

   enddo ZoneLoop



   return
   end subroutine ReConstructSets 


