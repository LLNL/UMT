!***********************************************************************
!   Created:  08/2021, TAB
!
!   GetTemperatureRelTol
!
!   Gets one of the iteration control values.
!
!***********************************************************************

   subroutine GetTemperatureRelTol(value) &
                        BIND(C,NAME="teton_get_temperature_reltol_internal")

   USE ISO_C_BINDING
   use iter_control_mod
   use iter_control_list_mod
   use Size_mod

   implicit none

   ! Input
   real(C_DOUBLE), intent(out) :: value

   ! Local variables
   type(IterControl), pointer  :: temperatureControl => NULL()

   temperatureControl => getIterationControl(IterControls,"temperature")

   value = getEpsilonPoint(temperatureControl)

   return
   end subroutine GetTemperatureRelTol
