!***********************************************************************
!   Created:  08/2021, TAB
!
!   AdjustTemperatureMaxIts
!
!   Sets one of the iteration control values.
!
!***********************************************************************

   subroutine AdjustTemperatureMaxIts(value) &
                        BIND(C,NAME="teton_adjust_temperature_maxits")

   USE ISO_C_BINDING
   use iter_control_mod
   use iter_control_list_mod
   use Size_mod

   implicit none

   ! Input
   integer(C_INT), intent(in) :: value

   ! Local variables
   type(IterControl), pointer  :: temperatureControl => NULL()

   temperatureControl => getIterationControl(IterControls,"temperature")

   call setControls(temperatureControl, maxNumberOfIterations=value)

   return
   end subroutine AdjustTemperatureMaxIts
