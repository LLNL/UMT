!***********************************************************************
!                        Version 0:  02/02, MKN                        *
!                                                                      *
!   CINTERFACE  -   Wrapper for modules that can be called from C++    * 
!                   used to get IterControls pointer                   *
!                                                                      *
!***********************************************************************

   subroutine ConstructIterControls( ) &
         BIND(C,NAME="teton_constructitercontrols")

!  Include
   use ISO_C_BINDING
   use kind_mod
   use iter_control_list_mod
   use iter_control_mod

   implicit none

!  Construct Iteration Controls

   allocate (IterControls)
   call construct(IterControls)
   call resetNumberOfIterations(IterControls)

   return
   end subroutine ConstructIterControls

