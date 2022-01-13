!***********************************************************************
!                     Version 0: 01/2005 PFN                           *
!                                                                      *
!    GETDTCONTROLINFO - API for controlling reason/process/zone.       *
!                                                                      *
!***********************************************************************

   subroutine getDtControlInfo(dtControlReason, &
                               dtControlProcess, &
                               dtControlZone) &
                           BIND(C,NAME="teton_getdtcontrolinfo")
   
   use ISO_C_BINDING
   use TimeStepControls_mod

   implicit none

!  Arguments

   integer(C_INT), intent(inout) :: dtControlReason
   integer(C_INT), intent(inout) :: dtControlProcess
   integer(C_INT), intent(inout) :: dtControlZone

!  Construct a message for the host code - append Null character for C

   dtControlReason = getDtConstraint(DtControls)
   dtControlProcess = getControlProcess(DtControls)
   dtControlZone    = getControlZone(DtControls)

   return
   end subroutine getDtControlInfo


