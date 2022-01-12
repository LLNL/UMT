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
   use QuadratureList_mod
   use RadIntensity_mod
   use SetData_mod

   implicit none 

!  Arguments

   integer(C_INT),    intent(in)    :: zoneID

   real(C_DOUBLE), intent(inout) :: beamMetric

!  Local

   type(RadIntensity),     pointer  :: RadT
   type(SetData),          pointer  :: Set

   integer    :: setID
   integer    :: nSets
   integer    :: g
   integer    :: g0
   integer    :: Groups
   real(adqt) :: EdDiag(Size%ndim)
   real(adqt) :: totRad


!  Constants

!***********************************************************************
!  Compute the radiation flux                                          *
!***********************************************************************

   nSets   =  getNumberOfSets(Quad)

   EdDiag(:) = zero
   totRad = zero

   SetLoop: do setID=1,nSets

     RadT   => getRadIntensity(Quad, setID)
     Set    => getSetData(Quad, setID)

     EdDiag(:) = EdDiag(:) + RadT% EddingtonTensorDiag(:,zoneID)
     totRad = totRad + RadT% radEnergy(zoneID)

   enddo SetLoop

   if( totRad > 0.0) then
      beamMetric = maxval(EdDiag) / totRad
   else
      beamMetric = one
   endif

   return
   end subroutine getBeamMetric


