!***********************************************************************
!                        Last Update:  02/2018, PFN                    *
!                                                                      *
!   advanceMaterialProperties: Update temperatures.                    *
!                                                                      *
!***********************************************************************
 
   subroutine advanceMaterialProperties(startCycle)


   use kind_mod
   use constant_mod
   use radconstant_mod
   use Size_mod
   use Geometry_mod
   use RadIntensity_mod
   use Material_mod

   implicit none

!  Arguments

   logical (kind=1), intent(in)  :: startCycle

!  Local

   integer                  :: c
   integer                  :: c0
   integer                  :: nCorner
   integer                  :: zone

   real(adqt)               :: tfloor
   real(adqt)               :: tr4min
   real(adqt)               :: ratio
   real(adqt)               :: Tstar
   real(adqt)               :: radEDensity 

!  Constants

   tfloor  = Size% tfloor
   tr4min  = Size% tr4floor

!***********************************************************************
!  Make corner temperatures consistent with zone averages obtained     *
!  from the host code.                                                 *
!                                                                      *
!  Advance zone temperatures [set old = new]                           *
!***********************************************************************

!  Start of cycle

   if ( startCycle ) then

!$omp  parallel do default(none) schedule(dynamic) &
!$omp& shared(Size, Mat, Geom, Rad, tfloor, tr4min) &
!$omp& private(c0, nCorner, Tstar, ratio, radEDensity)

     do zone=1,Size% nzones

       nCorner        =  Geom% numCorner(zone)
       c0             =  Geom% cOffSet(zone)
                                                                                                   
       Tstar          =  getZoneAverage(Geom, zone, Mat%Tec)
       Mat%trz(zone)  =  max( Mat%trz(zone), tfloor )
       Mat%tez(zone)  =  max( Mat%tez(zone), tfloor )

       Mat%tezn(zone) =  Mat%tez(zone)
       ratio          =  Mat%tez(zone)/Tstar

!      Initialize the material for the next time step

       Mat% tezold(zone) = zero

       do c=1,nCorner
         Mat% Tec(c0+c)    = max( ratio*Mat%Tec(c0+c), tfloor )
         Mat% Tecn(c0+c)   = Mat%tec(c0+c)
         Mat% denec(c0+c)  = zero
         Mat% tezold(zone) = Mat% tezold(zone) + Mat% Tec(c0+c)*Geom% Volume(c0+c)
         Mat% nonLinearIterations(c0+c) = 0
       enddo

       Mat% tezold(zone)  = Mat% tezold(zone)/Geom% VolumeZone(zone) 

!      Initialize the radiation energy (used for monitoring convergence)

       Mat% EnergyDensityOld(zone) = Rad% radEnergy(zone)/Geom% VolumeZone(zone)
       radEDensity                 = Rad% radEnergy(zone)/  &
                                    (Geom% VolumeZone(zone)*rad_constant*speed_light)

       Mat%Trzn(zone) = sqrt( sqrt( max(radEDensity, tr4min) ) )

     enddo

!$omp end parallel do

   else

!  End of cycle

!$omp  parallel do default(none) schedule(dynamic) &
!$omp& shared(Size, Mat, Geom) &
!$omp& private(c0, nCorner)

     do zone=1,Size% nzones

       nCorner          = Geom% numCorner(zone)
       c0               = Geom% cOffSet(zone)
       Mat% tez(zone)   = zero
       Mat% denez(zone) = zero

       do c=1,nCorner
         Mat% tez(zone)   = Mat% tez(zone)   + Geom% Volume(c0+c)*Mat%tec(c0+c)
         Mat% denez(zone) = Mat% denez(zone) + Geom% Volume(c0+c)*Mat% denec(c0+c) 
       enddo

       Mat% tez(zone)   = Mat% tez(zone)/Geom% VolumeZone(zone)
       Mat% denez(zone) = Mat% denez(zone)/Geom% VolumeZone(zone)

     enddo

!$omp end parallel do

   endif



   return
   end subroutine advanceMaterialProperties 



