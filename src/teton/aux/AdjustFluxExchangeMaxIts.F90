!***********************************************************************
!   Created:  08/2021, TAB
!
!   AdjustFluxExchangeMaxIts
!
!   Sets one of the iteration control values.
!
!***********************************************************************

   subroutine AdjustFluxExchangeMaxIts(value) &
                        BIND(C,NAME="teton_adjust_fluxexchange_maxits")

   USE ISO_C_BINDING
   use iter_control_mod
   use iter_control_list_mod
   use Size_mod

   implicit none

   ! Input
   integer(C_INT), intent(in) :: value

   ! Local variables
   type(IterControl), pointer  :: incidentFluxControl => NULL()

   incidentFluxControl => getIterationControl(IterControls,"incidentFlux")

   call setControls(incidentFluxControl, maxNumberOfIterations=value)

   return
   end subroutine AdjustFluxExchangeMaxIts
