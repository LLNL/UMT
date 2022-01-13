!***********************************************************************
!   Created:  08/2021, TAB
!
!   AdjustNonlinearRelTol
!
!   Sets one of the iteration control values.
!
!***********************************************************************

   subroutine AdjustNonlinearRelTol(value) &
                        BIND(C,NAME="teton_adjust_nonlinear_reltol")

   USE ISO_C_BINDING
   use iter_control_mod
   use iter_control_list_mod
   use Size_mod

   implicit none

   ! Input
   real(C_DOUBLE), intent(in) :: value

   ! Local variables
   type(IterControl), pointer  :: nonLinearControl => NULL()

   nonLinearControl => getIterationControl(IterControls,"nonLinear")

   call setControls(nonLinearControl, epsilonPoint=value)

   return
   end subroutine AdjustNonlinearRelTol
