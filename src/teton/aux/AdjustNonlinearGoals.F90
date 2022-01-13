!***********************************************************************
!                       Created:  05/2018, PGM                         *
!                                                                      *
!   AdjustNonlinearGoals                                               *
!                                                                      *
!   Called to adjust the maximum number of nonlinear solves and the    *
!   pointwise tolerance of the nonlinear solve, a.k.a the Newton       *
!   temperature update taken during each outer thermal iteration       *
!                                                                      *
!***********************************************************************

   subroutine AdjustNonlinearGoals(maxIters, goalEpsilon) &
                        BIND(C,NAME="teton_adjustnonlinearsolvegoals")

   USE ISO_C_BINDING
   use iter_control_mod
   use iter_control_list_mod

   implicit none

   ! Input
   integer(C_INT)  ,  intent(in) :: maxIters
   real(C_DOUBLE), intent(in) :: goalEpsilon
   
   ! Local variables
   type(IterControl), pointer  :: nonLinearControl => NULL()

   nonLinearControl => getIterationControl(IterControls,"nonLinear")
   
   call setControls(nonLinearControl, epsilonPoint=goalEpsilon, maxNumberOfIterations=maxIters)

   return
   end subroutine AdjustNonlinearGoals

