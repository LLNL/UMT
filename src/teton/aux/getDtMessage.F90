!***********************************************************************
!                     Version 0: 01/2005 PFN                           *
!                                                                      *
!    GETDTMESSAGE -  Called from host to get the time step control     *
!                    message.                                          *
!                                                                      *
!***********************************************************************

   subroutine getDtMessage(dtMessage) &
                           BIND(C,NAME="teton_getdtmessage")
   
   use ISO_C_BINDING
   use kind_mod
   use TimeStepControls_mod

   implicit none

!  Arguments

   type(C_PTR) :: dtMessage

!  Local

   character(len=:), allocatable, target, save :: dtString


!  Construct a message for the host code - append Null character for C

   dtString  = trim(DtControls% dtMessage)//C_NULL_CHAR
   dtMessage = C_LOC( dtString )

   return
   end subroutine getDtMessage

