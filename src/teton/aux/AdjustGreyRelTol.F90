!***********************************************************************
!   Created:  08/2021, TAB
!
!   AdjustGreyRelTol
!
!   Sets one of the iteration control values.
!
!***********************************************************************

   subroutine AdjustGreyRelTol(value) &
                        BIND(C,NAME="teton_adjust_grey_reltol")

   USE ISO_C_BINDING
   use iter_control_mod
   use iter_control_list_mod
   use Size_mod

   implicit none

   ! Input
   real(C_DOUBLE), intent(in) :: value

   ! Local variables
   type(IterControl), pointer  :: greyControl => NULL()

   greyControl => getIterationControl(IterControls,"grey")

   call setControls(greyControl, epsilonPoint=value)

   return
   end subroutine AdjustGreyRelTol
