!***********************************************************************
!   Created:  08/2021, TAB
!
!   GetGreyRelTol
!
!   Gets one of the iteration control values.
!
!***********************************************************************

   subroutine GetGreyRelTol(value) &
                        BIND(C,NAME="teton_get_grey_reltol_internal")

   USE ISO_C_BINDING
   use iter_control_mod
   use iter_control_list_mod
   use Size_mod

   implicit none

   ! Input
   real(C_DOUBLE), intent(out) :: value

   ! Local variables
   type(IterControl), pointer  :: greyControl => NULL()

   greyControl => getIterationControl(IterControls,"grey")

   value = getEpsilonPoint(greyControl)

   return
   end subroutine GetGreyRelTol
