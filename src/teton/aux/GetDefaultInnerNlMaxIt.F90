!***********************************************************************
!   Created:  08/2021, TAB
!
!   GetDefaultInnerNlMaxIt
!
!   Get default tolerance controls, for users to query before Teton setup
!   in case they need it for their parsers, etc.
!
!***********************************************************************

   subroutine GetDefaultInnerNlMaxIt(value) &
                        BIND(C,NAME="teton_get_default_inner_nl_max_it_internal")

   USE ISO_C_BINDING
   use default_iter_controls_mod
   use Size_mod

   implicit none

   ! Input
   integer(C_INT), intent(out) :: value

   value = inner_nl_max_it
   return
   end subroutine GetDefaultInnerNlMaxIt
