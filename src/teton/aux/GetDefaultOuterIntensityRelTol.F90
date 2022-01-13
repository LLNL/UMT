!***********************************************************************
!   Created:  08/2021, TAB
!
!   GetDefaultOuterIntensityRelTol
!
!   Get default tolerance controls, for users to query before Teton setup
!   in case they need it for their parsers, etc.
!
!***********************************************************************

   subroutine GetDefaultOuterIntensityRelTol(value) &
                        BIND(C,NAME="teton_get_default_outer_intensity_reltol_internal")

   USE ISO_C_BINDING
   use default_iter_controls_mod
   use Size_mod

   implicit none

   ! Input
   real(C_DOUBLE), intent(out) :: value

   value = outer_intensity_reltol
   return
   end subroutine GetDefaultOuterIntensityRelTol
