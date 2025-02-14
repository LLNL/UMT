#include "macros.h"
#include "omp_wrappers.h"
!***********************************************************************
!                        Last Update:  09/2017, PFN                    *
!                                                                      *
!   InitializeZones - Save zone-average quantities from previous cycle *
!                     for delta-t calculation.  Compute internal       *
!                     scattering coefficient.                          *
!                                                                      *
!***********************************************************************
 
   subroutine initializeZones 

   use cmake_defines_mod, only : omp_device_team_thread_limit
   use kind_mod
   use constant_mod
   use radconstant_mod
   use Size_mod
   use Geometry_mod
   use RadIntensity_mod
   use Editor_mod
   use QuadratureList_mod
   use ZoneSet_mod

   implicit none

!  Local


   integer    :: zone
   integer    :: nzones
   integer    :: zSetID
   integer    :: nZoneSets
   integer    :: c
   integer    :: c0
   integer    :: nCorner
   integer    :: g

   real(adqt) :: sumRad
   real(adqt) :: eradBOC 
   real(adqt) :: geometryFactor

   logical (kind=1) :: startCycle

!  Constants

   nzones         = Size% nzones
   nZoneSets      = getNumberOfZoneSets(Quad)
   geometryFactor = getGeometryFactor(Size)


!  Save the beginning of cycle radiation energy

   if ( Size% useGPU ) then

     TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none) &)
     TOMPC(shared(nZoneSets, ZSet, Geom, Rad))

     do zSetID=1,nZoneSets

       !$omp parallel do default(none) schedule(dynamic)  &
       !$omp& shared(zSetID, ZSet, Geom, Rad)

       do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
         ZSet% sumT(c) = Geom% Volume(c)*sum( Rad% PhiTotal(:,c) )
       enddo

       !$omp end parallel do

     enddo

     TOMP(end target teams distribute)

     TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none) &)
     TOMPC(shared(nZoneSets, Geom, Rad, ZSet)&)
     TOMPC(private(nCorner, c0))

     do zSetID=1,nZoneSets

       !$omp parallel do default(none) schedule(dynamic)  &
       !$omp& shared(zSetID, Geom, Rad, ZSet) &
       !$omp& private(c0, nCorner) 

       do zone=Geom% zone1(zSetID),Geom% zone2(zSetID)
         nCorner               = Geom% numCorner(zone)
         c0                    = Geom% cOffSet(zone)
         Rad% radEnergy(zone)  = zero

         do c=1,nCorner
           Rad% radEnergy(zone) = Rad% radEnergy(zone) + ZSet% sumT(c0+c)
         enddo
       enddo

       !$omp end parallel do

     enddo

     TOMP(end target teams distribute)

     TOMP(target update from(Rad% radEnergy))

   else

     !$omp parallel do default(none) schedule(dynamic)  &
     !$omp& shared(nzones, Geom, Rad, Size) &
     !$omp& private(c0, nCorner, sumRad) 

     do zone=1,nzones
       nCorner               = Geom% numCorner(zone)
       c0                    = Geom% cOffSet(zone)
       Rad% radEnergy(zone)  = zero

       do c=1,nCorner
         sumRad = zero

         do g=1,Size% ngr
           sumRad = sumRad + Rad% PhiTotal(g,c0+c)
         enddo

         Rad% radEnergy(zone) = Rad% radEnergy(zone) + Geom% Volume(c0+c)*sumRad
       enddo
     enddo

     !$omp end parallel do

   endif


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



