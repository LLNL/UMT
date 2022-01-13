!***********************************************************************
!                       Created:  05/2018, PGM                         *
!                                                                      *
!   AdjustGreyGoals                                                    *
!                                                                      *
!   Called to adjust the number of sweeps taken and the convergence    *
!   tolerance of the Grey accleration step in a LinearSolve operation  *
!   Note, this only does something in multi-d geometry                 *
!   Note2, there are two sweeps per CG iteration of the Grey solve     *
!                                                                      *
!***********************************************************************

   subroutine AdjustGreyGoals(maxIters, goalEpsilon) &
                        BIND(C,NAME="teton_adjustgreygoals")

   USE ISO_C_BINDING
   use iter_control_mod
   use iter_control_list_mod

   implicit none

   ! Input
   integer(C_INT)  ,  intent(in) :: maxIters
   real(C_DOUBLE), intent(in) :: goalEpsilon

   ! Local variables
   type(IterControl), pointer  :: greyControl => NULL()

   greyControl => getIterationControl(IterControls,"grey")

   call setControls(greyControl, epsilonPoint=goalEpsilon, maxNumberOfIterations=(2*maxIters+1))

   return
   end subroutine AdjustGreyGoals

