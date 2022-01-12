!***********************************************************************
!   Created:  08/2021, TAB
!
!   AdjustTemperatureRelTol
!
!   Sets one of the iteration control values.
!
!***********************************************************************

   subroutine AdjustTemperatureRelTol(value) &
                        BIND(C,NAME="teton_adjust_temperature_reltol")

   USE ISO_C_BINDING
   use iter_control_mod
   use iter_control_list_mod
   use Size_mod

   implicit none

   ! Input
   real(C_DOUBLE), intent(in) :: value

   ! Local variables
   type(IterControl), pointer  :: temperatureControl => NULL()

   temperatureControl => getIterationControl(IterControls,"temperature")

   call setControls(temperatureControl, epsilonPoint=value)

   return
   end subroutine AdjustTemperatureRelTol
