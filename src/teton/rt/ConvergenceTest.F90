#include "macros.h"
!***********************************************************************
!                        Version 1:  09/2011, PFN                      *
!                                                                      *
!   ConvergenceTest - Checks convergence by finding the maximum        *
!                     relative error in the radiation energy           *
!                     density and the material temperature.            *
!                                                                      *
!*********************************************************************** 
   subroutine ConvergenceTest(maxEnergyDensityError, maxTempError) 

   use kind_mod
   use mpi_param_mod
   use mpif90_mod
   use constant_mod
   use iter_control_list_mod
   use iter_control_mod
   use Size_mod
   use Geometry_mod
   use RadIntensity_mod
   use Material_mod

   implicit none

!  Arguments

   real(adqt), intent(inout) :: maxEnergyDensityError 
   real(adqt), intent(inout) :: maxTempError

!  Local

   integer    :: zoneEnergyMax 

   integer    :: myRankInGroup
   integer    :: zoneMaxError 
   integer    :: nzones
   integer    :: zone
   integer    :: RankZoneID(4)

   real(adqt) :: relerr
   real(adqt) :: cutoff
   real(adqt) :: threshold
   real(adqt) :: EnergyMax
   real(adqt) :: radEnergy
   real(adqt) :: AveTemp
   real(adqt) :: sumTempVol(2)
   real(adqt) :: maxError(2)

   type(IterControl) , pointer :: intensityControl   => NULL()
   type(IterControl) , pointer :: temperatureControl => NULL() 

!  Dynamic Arrays

   real(adqt),  allocatable :: EnergyDensity(:)

!  Constants

   parameter (cutoff=1.0d-6)

   nzones        = Size% nzones
   myRankInGroup = Size% myRankInGroup 
 
!  Compute the total energy density in a zone

   allocate( EnergyDensity(nzones) )

!  Iteration Controls

   intensityControl   => getIterationControl(IterControls,"intensity")
   temperatureControl => getIterationControl(IterControls,"temperature")

!  Compute a threshold temperature for checking convergence

   sumTempVol(:) = zero

!  Compute zone average temperature

   do zone=1,nzones
     Mat%tez(zone) = getZoneAverage(Geom, zone, Mat%Tec)

     sumTempVol(1) = sumTempVol(1) + Geom% VolumeZone(zone)*Mat%tez(zone)
     sumTempVol(2) = sumTempVol(2) + Geom% VolumeZone(zone)
   enddo

   call MPIAllReduce(sumTempVol, "sum", MY_COMM_GROUP)

   AveTemp = sumTempVol(1)/sumTempVol(2)

   ! Only zones with a temperature within 5% of the volumetric average
   ! on this spatial domain (above MPIAllReduce is over replicas)
   ! can "vote" on the zone average temperature error
   threshold = AveTemp/twenty

   !  Find maximum relative error in the zonal temperature
   !  Relative to the NEW electron temperature
   zoneMaxError = 1 
   maxTempError = zero

   do zone=1,nzones

     if (Mat%tez(zone) > threshold) then

       relerr = abs(one - Mat%tezold(zone)/Mat%tez(zone))

       if (relerr > maxTempError) then
         zoneMaxError = zone
         maxTempError = relerr
       endif

     endif

!  Save latest solution in TEZOLD for next iteration
     Mat%tezold(zone) = Mat%tez(zone)

   enddo

   maxError(1)   = maxTempError
   RankZoneID(1) = zoneMaxError

!  Set local error properties 

   call setLocalError(temperatureControl,maxTempError)
   call setZoneOfMax(temperatureControl,zoneMaxError)

!  Find the zone-average energy density 

   EnergyDensity(:) = zero

   do zone=1,nzones
     EnergyDensity(zone) = Rad% radEnergy(zone)/Geom% VolumeZone(zone)
   enddo

!  Find an energy threshold for convergence tests
   zoneEnergyMax  = maxloc( EnergyDensity, 1 )
   TETON_VERIFY(zoneEnergyMax >= 1, "Index of maximum value in EnergyDensity came back 0.  Check that the array values aren't 0 or NAN.")
   EnergyMax      = EnergyDensity( zoneEnergyMax )

   call MPIAllReduce(EnergyMax, "max", MY_COMM_GROUP)

   ! cutoff is defined above as 10^{-6} above
   threshold = cutoff*EnergyMax

!  Compute relative errors in the total energy density in
!  a zone; eliminate zones from consideration if their zonal
!  energy is less than a threshold 

   zoneMaxError          = 1 
   maxEnergyDensityError = zero
 
   do zone=1,nzones

     radEnergy = EnergyDensity(zone)

     if (RadEnergy > threshold) then

       relerr = abs( (EnergyDensity(zone) - Mat% EnergyDensityOld(zone))/ &
                      EnergyDensity(zone) )

       if (relerr > maxEnergyDensityError) then
         zoneMaxError          = zone 
         maxEnergyDensityError = relerr
       endif

     endif

     Mat% EnergyDensityOld(zone) = EnergyDensity(zone)

   enddo

   maxError(2)   = maxEnergyDensityError
   RankZoneID(2) = zoneMaxError

!  Set local error properties

   call setLocalError(intensityControl,maxEnergyDensityError)
   call setZoneOfMax(intensityControl,zoneMaxError)

!  Find the largest errors in the entire mesh and the RankZoneID

   call MPIAllReduce(maxError, "max", MY_COMM_GROUP)

   if (maxTempError == maxError(1)) then
     RankZoneID(3) = myRankInGroup    
   else
     RankZoneID(3) = -1
   endif

   if (maxEnergyDensityError == maxError(2)) then
     RankZoneID(4) = myRankInGroup 
   else
     RankZoneID(4) = -1
   endif

   ! Need to split up this reduction/sync in case there is a tie:
   call MPIAllReduce(RankZoneID(3:4), "max", MY_COMM_GROUP)
   call MPIBcast(RankZoneID(1:1), 1, RankZoneID(3), MY_COMM_GROUP)
   call MPIBcast(RankZoneID(2:2), 1, RankZoneID(4), MY_COMM_GROUP)

!  Set global error properties

   maxTempError          = maxError(1)
   maxEnergyDensityError = maxError(2)

   call setGlobalError(temperatureControl,maxError(1))
   call setZoneOfMax(temperatureControl,RankZoneID(1))
   call setProcessOfMax(temperatureControl,RankZoneID(3))

   call setGlobalError(intensityControl,maxError(2))
   call setZoneOfMax(intensityControl,RankZoneID(2))
   call setProcessOfMax(intensityControl,RankZoneID(4))


!  Release Memory

   deallocate( EnergyDensity )
 

   return
   end subroutine ConvergenceTest 

