!***********************************************************************
!   Created:  08/2021, TAB
!
!   GetTemperatureMaxIts
!
!   Gets one of the iteration control values.
!
!***********************************************************************

   subroutine GetTemperatureMaxIts(value) &
                        BIND(C,NAME="teton_get_temperature_maxits_internal")

   USE ISO_C_BINDING
   use iter_control_mod
   use iter_control_list_mod
   use Size_mod

   implicit none

   ! Input
   integer(C_INT), intent(out) :: value

   ! Local variables
   type(IterControl), pointer  :: temperatureControl => NULL()

   temperatureControl => getIterationControl(IterControls,"temperature")

   value = getMaxNumberOfIterations(temperatureControl)

   return
   end subroutine GetTemperatureMaxIts
