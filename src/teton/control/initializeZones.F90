#include "macros.h"
!***********************************************************************
!                        Last Update:  09/2017, PFN                    *
!                                                                      *
!   InitializeZones - Save zone-average quantities from previous cycle *
!                     for delta-t calculation.  Compute internal       *
!                     scattering coefficient.                          *
!                                                                      *
!***********************************************************************
 
   subroutine initializeZones 

   use kind_mod
   use constant_mod
   use radconstant_mod
   use Size_mod
   use Geometry_mod
   use RadIntensity_mod
   use Editor_mod

   implicit none

!  Local


   integer    :: zone
   integer    :: nzones
   integer    :: c
   integer    :: c0
   integer    :: nCorner
   integer    :: g

   real(adqt) :: eradBOC 
   real(adqt) :: geometryFactor
   real(adqt) :: sumRad

   logical (kind=1) :: startCycle

!  Constants

   nzones         = Size% nzones
   geometryFactor = getGeometryFactor(Size)


!  Save the beginning of cycle radiation energy

   !$omp parallel do default(none) schedule(dynamic)  &
   !$omp& shared(Size, Geom, Rad) &
   !$omp& private(c0, nCorner, sumRad) 

   do zone=1,nzones
     nCorner              = Geom% numCorner(zone)
     c0                   = Geom% cOffSet(zone)
     Rad% radEnergy(zone) = zero

     do c=1,nCorner
       sumRad = zero

       do g=1,Size% ngr
         sumRad = sumRad + Rad% PhiTotal(g,c0+c)
       enddo

       Rad% radEnergy(zone) = Rad% radEnergy(zone) + Geom% Volume(c0+c)*sumRad
     enddo
   enddo
   !$omp end parallel do


   eradBOC = zero

   do zone=1,nzones
     eradBOC = eradBOC + Rad% radEnergy(zone)
   enddo

   RadEdit% EnergyRadBOC = geometryFactor*eradBOC/speed_light

!***********************************************************************
!  Make corner temperatures consistent with zone averages obtained     *
!  from the host code.                                                 *
!                                                                      *
!  Advance zone temperatures [set old = new]                           *
!***********************************************************************

   startCycle = .TRUE.
   call advanceMaterialProperties(startCycle)


   return
   end subroutine initializeZones 



