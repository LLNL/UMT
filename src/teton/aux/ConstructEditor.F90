!***********************************************************************
!                        Version 1:  01/2015, PFN                      *
!                                                                      *
!   ConstructEditor - Constructor for tallies, including the common    *
!                     radiation boundary edits.                        *
!                                                                      *
!***********************************************************************


   subroutine ConstructEditor(ngr, numSpectrumAngleBins,        &
                              spectrumAngleBinBoundaries,       &
                              RadPowerEscape, RadPowerIncident, &
                              PolarSectorPowerEscape) &
                              BIND(C,NAME="teton_constructeditor")

!  Include

   use ISO_C_BINDING
   use kind_mod
   use Editor_mod

   implicit none

!  Arguments

   integer(C_INT),    intent(in)         :: ngr
   integer(C_INT),    intent(in)         :: numSpectrumAngleBins
   real(C_DOUBLE), intent(in)         :: spectrumAngleBinBoundaries &
                                             (numSpectrumAngleBins+1) 
   real(C_DOUBLE), target, intent(in) :: RadPowerEscape(ngr)
   real(C_DOUBLE), target, intent(in) :: RadPowerIncident(ngr)
   real(C_DOUBLE), target, intent(in) :: PolarSectorPowerEscape(numSpectrumAngleBins*ngr)


!  Construct Problem Edits 

   allocate (RadEdit)

   call construct(RadEdit,                    &
                  ngr,                        &
                  numSpectrumAngleBins,       &
                  spectrumAngleBinBoundaries, &
                  RadPowerEscape,             &
                  RadPowerIncident,           &
                  PolarSectorPowerEscape)

   return
   end subroutine ConstructEditor

