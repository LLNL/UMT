!***********************************************************************
!   Created:  08/2021, TAB
!
!   AdjustNonlinearMaxIts
!
!   Sets one of the iteration control values.
!
!***********************************************************************

   subroutine AdjustNonlinearMaxIts(value) &
                        BIND(C,NAME="teton_adjust_nonlinear_maxits")

   USE ISO_C_BINDING
   use iter_control_mod
   use iter_control_list_mod
   use Size_mod

   implicit none

   ! Input
   integer(C_INT), intent(in) :: value

   ! Local variables
   type(IterControl), pointer  :: nonLinearControl => NULL()

   nonLinearControl => getIterationControl(IterControls,"nonLinear")

   call setControls(nonLinearControl, maxNumberOfIterations=value)

   return
   end subroutine AdjustNonlinearMaxIts
