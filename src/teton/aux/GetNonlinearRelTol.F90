!***********************************************************************
!   Created:  08/2021, TAB
!
!   GetNonlinearRelTol
!
!   Gets one of the iteration control values.
!
!***********************************************************************

   subroutine GetNonlinearRelTol(value) &
                        BIND(C,NAME="teton_get_nonlinear_reltol_internal")

   USE ISO_C_BINDING
   use iter_control_mod
   use iter_control_list_mod
   use Size_mod

   implicit none

   ! Input
   real(C_DOUBLE), intent(out) :: value

   ! Local variables
   type(IterControl), pointer  :: nonLinearControl => NULL()

   nonLinearControl => getIterationControl(IterControls,"nonLinear")

   value = getEpsilonPoint(nonLinearControl)

   return
   end subroutine GetNonlinearRelTol
