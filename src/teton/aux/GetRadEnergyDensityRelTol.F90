!***********************************************************************
!   Created:  08/2021, TAB
!
!   GetRadEnergyDensityRelTol
!
!   Gets one of the iteration control values.
!
!***********************************************************************

   subroutine GetRadEnergyDensityRelTol(value) &
                        BIND(C,NAME="teton_get_radenergydensity_reltol_internal")

   USE ISO_C_BINDING
   use iter_control_mod
   use iter_control_list_mod
   use Size_mod

   implicit none

   ! Input
   real(C_DOUBLE), intent(out) :: value

   ! Local variables
   type(IterControl), pointer  :: intensityControl => NULL()

   intensityControl => getIterationControl(IterControls,"intensity")

   value = getEpsilonPoint(intensityControl)

   return
   end subroutine GetRadEnergyDensityRelTol
