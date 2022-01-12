!***********************************************************************
!                         Last Update: 01/2012 PFN                     *
!                                                                      *
!    setTimeStep  -  Called from host to update the radiation          *
!                    time and timestep before calling RADTR.           *
!                                                                      *
!    Input:   dtrad      - radiation timestep                          *
!             timerad    - radiation time                              *
!                                                                      *
!    Output:  DtControls - structure containing timestep information   *
!                                                                      *
!***********************************************************************
   subroutine setTimeStep(ncycle, dtrad, timerad, tfloor) &
     BIND(C,NAME="teton_settimestep")

   use kind_mod
   use Size_mod
   use TimeStepControls_mod
   use ISO_C_BINDING 
   
   implicit none


!  Arguments

   integer(C_INT), intent(in) :: ncycle

   real(C_DOUBLE), intent(in) :: dtrad 
   real(C_DOUBLE), intent(in) :: timerad
   real(C_DOUBLE), intent(in) :: tfloor

!  Update controls
                                                                                             
   call setDtControls(DtControls,         &
                      RadCycle=ncycle,    &
                      RadTimeStep=dtrad,  &
                      RadTime=timerad     )

   Size% tfloor = tfloor


   return
   end subroutine setTimeStep 



