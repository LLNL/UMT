!***********************************************************************
!                        Version 1:  05/92, PFN                        *
!                                                                      *
!   ADVANCERT - Save zone-average quantities from previous cycle for   *
!               delta-t calculation.  Convert specific radiation       *
!               intensity (i.e. per unit mass) to intensity (per       *
!               unit volume) before the transport calculation.         *
!                                                                      *
!   Input:   tez,trz                                                   *
!                                                                      *
!   Output:  TEZN,TRZN                                                 *
!                                                                      *
!***********************************************************************
 
   subroutine advanceRT(dtrad, PSIR, PHI) 


   use kind_mod
   use constant_mod
   use radconstant_mod
   use Size_mod
   use Geometry_mod
   use Material_mod
   use QuadratureList_mod
   use Quadrature_mod
   use Editor_mod
   use ZoneData_mod

   implicit none

!  Arguments

   real(adqt), intent(in)    :: dtrad

   real(adqt), intent(inout) :: psir(Size%ngr,Size%ncornr,Size%nangSN), &
                                Phi(Size%ngr,Size%ncornr)

!  Local

   integer    :: ia, ic, ig
   integer    :: c, c0, nCorner, zone
   integer    :: nzones, ngr, numAngles

   real(adqt) :: factor, tfloor, tr4min, ratio, eradBOC 
   real(adqt) :: deltaVolume, aveVolume, volumeRatio
   real(adqt) :: sumRad, PhiAve, Tstar

!  Constants

   nzones     = Size%nzones
   ngr        = Size%ngr
   tfloor     = Size%tfloor

   QuadSet    => getQuadrature(Quad,1)
   numAngles  = QuadSet%NumAngles

!***********************************************************************
!  Make corner temperatures consistent with zone averages obtained     *
!  from the host code.                                                 *
!                                                                      *
!  Advance zone temperatures [set old = new]                           *
!***********************************************************************

   call timer_beg('_ZoneLoop0')

!$omp parallel do private(zone,Z,nCorner,c0,Tstar,c,ratio)
   ZoneLoop: do zone=1,nzones

     Z       => getZoneData(Geom, zone)
                                                                                                   
     nCorner = Z% nCorner
     c0      = Z% c0
                                                                                                   
     Tstar = zero
     do c=1,nCorner
       Tstar = Tstar + Z% Volume(c)*Mat%tec(c0+c)
     enddo
      Tstar = Tstar/Z% VolumeZone

     Mat%trz(zone)  = max( Mat%trz(zone), tfloor )
     Mat%tez(zone)  = max( Mat%tez(zone), tfloor )

     Mat%tezn(zone) = Mat%tez(zone)
     ratio          = Mat%tez(zone)/Tstar

     do c=1,nCorner
       Mat%tec(c0+c) = max( ratio*Mat%tec(c0+c), tfloor )
     enddo

   enddo ZoneLoop

   call timer_end('_ZoneLoop0')

!***********************************************************************
!  Set the scaler intensity                                            *
!***********************************************************************

   call timer_beg('_snmoments1')

   call snmoments(psir, PHI)

   call timer_end('_snmoments1')

!***********************************************************************
!  Compute the work done on radiation field due to volume changes.     *
!  This is an external source rate in the radiation transport equation.*
!***********************************************************************

   if (Size%radForceMultiplier > zero) then

     Mat%qext(:,:) = zero

     factor = -third*Size%radForceMultiplier/(dtrad*speed_light)

   call timer_beg('_ZoneLoop1')

     ZoneLoop1: do zone=1,nzones

       Z    => getZoneData(Geom, zone)

       nCorner = Z% nCorner
       c0      = Z% c0

       do c=1,nCorner
         deltaVolume      =       Z% Volume(c) - Z% VolumeOld(c) 
         aveVolume        = half*(Z% Volume(c) + Z% VolumeOld(c))
         volumeRatio      = factor*deltaVolume/aveVolume

         Mat%qext(:,c0+c) = volumeRatio*Phi(:,c0+c) 
       enddo

     enddo ZoneLoop1

     call timer_end('_ZoneLoop1')

   endif

!***********************************************************************
!  Scale the radiation field to account for volume changes and         *
!  tally beginning-of-cycle radiation energy                           *
!***********************************************************************
 
   eradBOC = zero
   tr4min  = tfloor*tfloor*tfloor*tfloor

   call timer_beg('_ZoneLoop2')

!$omp parallel do private(zone,PhiAve,Z,nCorner,c0,c,volumeRatio,ia,sumRad) reduction(+:eradBOC)
   ZoneLoop2: do zone=1,nzones

     PhiAve  = zero

     Z       => getZoneData(Geom, zone)
     nCorner = Z% nCorner
     c0      = Z% c0

     do c=1,nCorner
       volumeRatio = Z% VolumeOld(c)/Z% Volume(c)

       ! store volume ratio for later transfer to gpu and use with scalePsibyVol
       Z%volumeRatio(c) = volumeRatio 

       Phi(:,c0+c) = Phi(:,c0+c)*volumeRatio

       ! This scaling of psi is now done with psi on the GPU in snflwxyz.F90.
       ! or comment that out and do here:
       do ia=1,numAngles
          psir(:,c0+c,ia) = psir(:,c0+c,ia)*volumeRatio
       enddo

       sumRad = zero
       do ig=1,ngr
         sumRad = sumRad + Phi(ig,c0+c)
       enddo

       PhiAve = PhiAve + Z% Volume(c)*sumRad

     enddo


     eradBOC = eradBOC + PhiAve
     PhiAve  = PhiAve/(Z% VolumeZone*rad_constant*speed_light)

     Mat%trzn(zone) = sqrt( sqrt( max(PhiAve, tr4min) ) )

   enddo ZoneLoop2

   call timer_end('_ZoneLoop2')

   if (Size%igeom == 'rz') then
     RadEdit% EnergyRadBOC = two*pi*eradBOC/speed_light
   else
     RadEdit% EnergyRadBOC = eradBOC/speed_light
   endif


   return
   end subroutine advanceRT



