!***********************************************************************
!                        Last Update:  02/2018, PFN                    *
!                                                                      *
!   advanceMaterialProperties: Update temperatures.                    *
!                                                                      *
!***********************************************************************
 
   subroutine advanceMaterialProperties(zone)


   use kind_mod
   use constant_mod
   use radconstant_mod
   use Size_mod
   use Geometry_mod
   use Material_mod
   use ZoneData_mod

   implicit none

!  Arguments

   integer, intent(in)      :: zone

!  Local

   type(ZoneData), pointer  :: ZT

   integer                  :: c
   integer                  :: c0
   integer                  :: nCorner
   integer                  :: g

   real(adqt)               :: tfloor
   real(adqt)               :: tr4min
   real(adqt)               :: ratio
   real(adqt)               :: Tstar
   real(adqt)               :: radEDensity 
   real(adqt)               :: sumRad

!  Constants

   tfloor  = Size% tfloor
   tr4min  = Size% tr4floor

!***********************************************************************
!  Make corner temperatures consistent with zone averages obtained     *
!  from the host code.                                                 *
!                                                                      *
!  Advance zone temperatures [set old = new]                           *
!***********************************************************************


   ZT       => getZoneData(Geom, zone)
   nCorner  =  Geom% numCorner(zone)
   c0       =  Geom% cOffSet(zone)
                                                                                                   
   Tstar          = getZoneAverage(Geom, zone, Mat%Tec)
   Mat%trz(zone)  = max( Mat%trz(zone), tfloor )
   Mat%tez(zone)  = max( Mat%tez(zone), tfloor )

   Mat%tezn(zone) = Mat%tez(zone)
   ratio          = Mat%tez(zone)/Tstar

!  Initialize the material for the next time step

   Mat% tezold(zone) = zero

   do c=1,nCorner
     Mat% Tec(c0+c)    = max( ratio*Mat%Tec(c0+c), tfloor )
     Mat% Tecn(c0+c)   = Mat%tec(c0+c)
     Mat% denec(c0+c)  = zero
     Mat% tezold(zone) = Mat% tezold(zone) + Mat% Tec(c0+c)*Geom% Volume(c0+c)
   enddo

   Mat% tezold(zone)  = Mat% tezold(zone)/Geom% VolumeZone(zone) 

!  Initialize the non-linear iteration counter

   Mat% nonLinearIterations(zone) = 0

!  Initialize the radiation energy (used for monitoring convergence)

   Geom% radEnergy(zone) = zero

   do c=1,nCorner
     sumRad = zero
     do g=1,Size% ngr
       sumRad = sumRad + Geom% PhiTotal(g,c0+c)
     enddo
     Geom% radEnergy(zone) = Geom% radEnergy(zone) + Geom% Volume(c0+c)*sumRad
   enddo

   ZT% EnergyDensityOld = Geom% radEnergy(zone)/Geom% VolumeZone(zone)
   radEDensity          = Geom% radEnergy(zone)/  &
                         (Geom% VolumeZone(zone)*rad_constant*speed_light)

   Mat%Trzn(zone) = sqrt( sqrt( max(radEDensity, tr4min) ) )



   return
   end subroutine advanceMaterialProperties 



