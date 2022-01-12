!***********************************************************************
!   Created:  08/2021, TAB
!
!   AdjustGreyMaxIts
!
!   Sets one of the iteration control values.
!
!***********************************************************************

   subroutine AdjustGreyMaxIts(iters) &
                        BIND(C,NAME="teton_adjust_grey_maxits")

   USE ISO_C_BINDING
   use iter_control_mod
   use iter_control_list_mod
   use Size_mod

   implicit none

   ! Input
   integer(C_INT), intent(in) :: iters

   ! Local variables
   type(IterControl), pointer  :: greyControl => NULL()

   greyControl => getIterationControl(IterControls,"grey")

   ! Really set sweeps
   call setControls(greyControl, maxNumberOfIterations=(2*iters+1))

   return
   end subroutine AdjustGreyMaxIts
