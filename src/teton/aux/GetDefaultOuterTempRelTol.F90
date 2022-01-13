!***********************************************************************
!   Created:  08/2021, TAB
!
!   GetDefaultOuterTempRelTol
!
!   Get default tolerance controls, for users to query before Teton setup
!   in case they need it for their parsers, etc.
!
!***********************************************************************

   subroutine GetDefaultOuterTempRelTol(value) &
                        BIND(C,NAME="teton_get_default_outer_temp_reltol_internal")

   USE ISO_C_BINDING
   use default_iter_controls_mod
   use Size_mod

   implicit none

   ! Input
   real(C_DOUBLE), intent(out) :: value

   value = outer_temp_reltol
   return
   end subroutine GetDefaultOuterTempRelTol
