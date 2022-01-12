!***********************************************************************
!   Created:  08/2021, TAB
!
!   GetGreyMaxIts
!
!   Gets one of the iteration control values.
!
!***********************************************************************

   subroutine GetGreyMaxIts(value) &
                        BIND(C,NAME="teton_get_grey_maxits_internal")

   USE ISO_C_BINDING
   use iter_control_mod
   use iter_control_list_mod
   use Size_mod

   implicit none

   ! Input
   integer(C_INT), intent(out) :: value

   ! Local variables
   type(IterControl), pointer  :: greyControl => NULL()

   greyControl => getIterationControl(IterControls,"grey")

   value = getMaxNumberOfIterations(greyControl)

   return
   end subroutine GetGreyMaxIts
