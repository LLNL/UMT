!***********************************************************************
!   Created:  08/2021, TAB
!
!   GetFluxExchangeMaxIts
!
!   Gets one of the iteration control values.
!
!***********************************************************************

   subroutine GetFluxExchangeMaxIts(value) &
                        BIND(C,NAME="teton_get_fluxexchange_maxits_internal")

   USE ISO_C_BINDING
   use iter_control_mod
   use iter_control_list_mod
   use Size_mod

   implicit none

   ! Input
   integer(C_INT), intent(out) :: value

   ! Local variables
   type(IterControl), pointer  :: incidentFluxControl => NULL()

   incidentFluxControl => getIterationControl(IterControls,"incidentFlux")

   value = getMaxNumberOfIterations(incidentFluxControl)

   return
   end subroutine GetFluxExchangeMaxIts
