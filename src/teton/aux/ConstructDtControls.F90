!***********************************************************************
!                        Version 0:  02/02, MKN                        *
!                                                                      *
!   CINTERFACE  -   Wrapper for modules that can be called from C++    * 
!                   used to get DtControls pointer                     *
!                                                                      *
!***********************************************************************


   subroutine ConstructDtControls(dtrad, dtrmn, dtrmx, delte, deltr) &
        BIND(C,NAME="teton_constructdtcontrols")

!  Include
   use ISO_C_BINDING
   use kind_mod
   use TimeStepControls_mod


   implicit none

!  Arguments

!  Time Step Controls
   real(C_DOUBLE), intent(in)    :: dtrad
   real(C_DOUBLE), intent(in)    :: dtrmn
   real(C_DOUBLE), intent(in)    :: dtrmx
   real(C_DOUBLE), intent(in)    :: delte
   real(C_DOUBLE), intent(in)    :: deltr

!  Construct Time Step Controls

   allocate (DtControls)

   call construct(DtControls,         & 
                  RadTimeStep=dtrad,  &
                  MaxChangeTe=delte,  &
                  MaxChangeEr=deltr,  &
                  MinTimeStep=dtrmn,  &
                  MaxTimeStep=dtrmx)


   return
   end subroutine ConstructDtControls

