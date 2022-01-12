!***********************************************************************
!   Created:  08/2021, TAB
!
!   AdjustRadEnergyDensityRelTol
!
!   Sets one of the iteration control values.
!
!***********************************************************************

   subroutine AdjustRadEnergyDensityRelTol(value) &
                        BIND(C,NAME="teton_adjust_radenergydensity_reltol")

   USE ISO_C_BINDING
   use iter_control_mod
   use iter_control_list_mod
   use Size_mod

   implicit none

   ! Input
   real(C_DOUBLE), intent(in) :: value

   ! Local variables
   type(IterControl), pointer  :: intensityControl => NULL()

   intensityControl => getIterationControl(IterControls,"intensity")

   call setControls(intensityControl, epsilonPoint=value)

   return
   end subroutine AdjustRadEnergyDensityRelTol
