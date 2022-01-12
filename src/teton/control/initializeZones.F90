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
   use Editor_mod

   implicit none

!  Local


   integer         :: zone
   integer         :: nzones

   real(adqt)      :: eradBOC 
   real(adqt)      :: geometryFactor

!  Constants

   nzones         = Size% nzones
   geometryFactor = getGeometryFactor(Size)

!***********************************************************************
!  Make corner temperatures consistent with zone averages obtained     *
!  from the host code.                                                 *
!                                                                      *
!  Advance zone temperatures [set old = new]                           *
!***********************************************************************

!$omp parallel do private(zone) schedule(dynamic)

   do zone=1,nzones
     call advanceMaterialProperties(zone)
   enddo

!$omp end parallel do

!  Save the beginning of cycle radiation energy

   eradBOC = zero

   do zone=1,nzones
     eradBOC = eradBOC + Geom% radEnergy(zone)
   enddo

   RadEdit% EnergyRadBOC = geometryFactor*eradBOC/speed_light


   return
   end subroutine initializeZones 



