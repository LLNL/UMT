!***********************************************************************
!                       Created:  05/2018, PGM                         *
!                                                                      *
!   AdjustFluxExchangeGoals                                            *
!                                                                      *
!   Called to adjust the number of sweeps taken and the convergence    *
!   tolerance of the Grey accleration step in a LinearSolve operation  *
!   Note, this only does something in multi-d geometry                 *
!   Note2, there are two sweeps per CG iteration of the Grey solve     *
!                                                                      *
!***********************************************************************

   subroutine AdjustFluxExchangeGoals(maxIters, goalEpsilon) &
                        BIND(C,NAME="teton_adjustfluxexchangegoals")

   USE ISO_C_BINDING
   use iter_control_mod
   use iter_control_list_mod
   use Size_mod

   implicit none

   ! Input
   integer(C_INT)  ,  intent(in) :: maxIters
   real(C_DOUBLE), intent(in) :: goalEpsilon
   
   ! Local variables
   type(IterControl), pointer  :: incidentFluxControl => NULL()

   incidentFluxControl => getIterationControl(IterControls,"incidentFlux")
   
   call setControls(incidentFluxControl, epsilonPoint=goalEpsilon, maxNumberOfIterations=maxIters)

   return
   end subroutine AdjustFluxExchangeGoals

