!***********************************************************************
!   Created:  08/2021, TAB
!
!   GetNonlinearMaxIts
!
!   Gets one of the iteration control values.
!
!***********************************************************************

   subroutine GetNonlinearMaxIts(value) &
                        BIND(C,NAME="teton_get_nonlinear_maxits_internal")

   USE ISO_C_BINDING
   use iter_control_mod
   use iter_control_list_mod
   use Size_mod

   implicit none

   ! Input
   integer(C_INT), intent(out) :: value

   ! Local variables
   type(IterControl), pointer  :: nonLinearControl => NULL()

   nonLinearControl => getIterationControl(IterControls,"nonLinear")

   value = getMaxNumberOfIterations(nonLinearControl)

   return
   end subroutine GetNonlinearMaxIts
