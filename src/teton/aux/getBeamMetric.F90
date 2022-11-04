!***********************************************************************
!                       Last Update:  03/2012, PFN                     *
!                                                                      *
!   getBeamMetric    - computes the zone and group averaged maximum    *
!                      diagonal of Eddington tensor.  Will be 1/3 for  *
!                      isotropic cases, 1 for beams.  Also set to 1    *
!                      if the energy density is non-positive           *
!                                                                      *
!***********************************************************************

   subroutine getBeamMetric(zoneID, beamMetric) &
        BIND(C,NAME="teton_getbeammetric")

   USE ISO_C_BINDING
   use kind_mod
   use constant_mod
   use Size_mod
   use RadIntensity_mod

   implicit none 

!  Arguments

   integer(C_INT),    intent(in) :: zoneID

   real(C_DOUBLE), intent(inout) :: beamMetric

!  Local

   real(adqt) :: EdDiag(Size%ndim)
   real(adqt) :: totRad


!  Constants

!***********************************************************************
!  Compute the beam metric                                             *
!***********************************************************************


   EdDiag(:) = Rad% EddingtonTensorDiag(:,zoneID)
   totRad    = Rad% radEnergy(zoneID)

   if( totRad > zero) then
      beamMetric = maxval(EdDiag) / totRad
   else
      beamMetric = one
   endif

   return
   end subroutine getBeamMetric


