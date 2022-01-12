!***********************************************************************
!                       Created:  05/2018, PGM                         *
!                                                                      *
!   AdjustRadEnergyDensityGoal                                         *
!                                                                      *
!   Called to adjust pointwise tolerance of radiation energy density   *
!                                                                      *
!***********************************************************************

   subroutine AdjustRadEnergyDensityGoal(goalEpsilon) &
                        BIND(C,NAME="teton_adjustradenergydensitygoal")

   USE ISO_C_BINDING
   use iter_control_mod
   use iter_control_list_mod

   implicit none

   ! Input
   real(C_DOUBLE),  intent(in) :: goalEpsilon
   
   ! Local variables
   type(IterControl), pointer  :: intensityControl => NULL()

   intensityControl => getIterationControl(IterControls,"intensity")
   
   call setControls(intensityControl, epsilonPoint=goalEpsilon)

   return
   end subroutine AdjustRadEnergyDensityGoal

