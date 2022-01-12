!***********************************************************************
!                       Created:  05/2018, PGM                         *
!                                                                      *
!   AdjustTemperatureGoals                                             *
!                                                                      *
!   Called to adjust the maximum number of temperature [outer]         *
!   iterations taken per time step by Teton.  Also adjusts the         *
!   pointwise error tolerance on zone temperature.                     *
!                                                                      *
!***********************************************************************

   subroutine AdjustTemperatureGoals(maxIters, goalEpsilon) &
                        BIND(C,NAME="teton_adjusttemperaturegoals")

   USE ISO_C_BINDING
   use iter_control_mod
   use iter_control_list_mod

   implicit none

   ! Input
   integer(C_INT)  ,  intent(in) :: maxIters
   real(C_DOUBLE), intent(in) :: goalEpsilon
   
   ! Local variables
   type(IterControl), pointer  :: temperatureControl => NULL()

   temperatureControl => getIterationControl(IterControls,"temperature")
   
   call setControls(temperatureControl, epsilonPoint=goalEpsilon, maxNumberOfIterations=maxIters)

   return
   end subroutine AdjustTemperatureGoals

