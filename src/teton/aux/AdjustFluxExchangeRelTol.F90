!***********************************************************************
!   Created:  08/2021, TAB
!
!   AdjustFluxExchangeRelTol
!
!   Sets one of the iteration control values.
!
!***********************************************************************

   subroutine AdjustFluxExchangeRelTol(value) &
                        BIND(C,NAME="teton_adjust_fluxexchange_reltol")

   USE ISO_C_BINDING
   use iter_control_mod
   use iter_control_list_mod
   use Size_mod

   implicit none

   ! Input
   real(C_DOUBLE), intent(in) :: value

   ! Local variables
   type(IterControl), pointer  :: incidentFluxControl => NULL()

   incidentFluxControl => getIterationControl(IterControls,"incidentFlux")

   call setControls(incidentFluxControl, epsilonPoint=value)

   return
   end subroutine AdjustFluxExchangeRelTol
