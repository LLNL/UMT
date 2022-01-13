!***********************************************************************
!   Created:  08/2021, TAB
!
!   GetFluxExchangeRelTol
!
!   Gets one of the iteration control values.
!
!***********************************************************************

   subroutine GetFluxExchangeRelTol(value) &
                        BIND(C,NAME="teton_get_fluxexchange_reltol_internal")

   USE ISO_C_BINDING
   use iter_control_mod
   use iter_control_list_mod
   use Size_mod

   implicit none

   ! Input
   real(C_DOUBLE), intent(out) :: value

   ! Local variables
   type(IterControl), pointer  :: incidentFluxControl => NULL()

   incidentFluxControl => getIterationControl(IterControls,"incidentFlux")

   value = getEpsilonPoint(incidentFluxControl)

   return
   end subroutine GetFluxExchangeRelTol
