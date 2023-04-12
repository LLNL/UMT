!***********************************************************************
!                        Last Update:  01/2012, PFN                    *
!                                                                      *
!   ResetSize      -  Called before each radiation cycle to set        * 
!                     parameters that may have been changed by the     *
!                     user.                                            *
!                                                                      *
!***********************************************************************

   subroutine ResetSize(functionRNLTE, tfloor,                 &
                        radForceMultiplier, betaNLTE,          &
                        gammaNLTE,                             &
                        DopplerShiftOn, useNewNonLinearSolver, &
                        useNewGTASolver,                       &
                        usePWLD, useSurfaceMassLumping)        &
                        BIND(C,NAME="teton_resetsize")

!  Include
   USE ISO_C_BINDING
   use kind_mod
   use Size_mod


   implicit none

!  Arguments

   integer(C_INT),   intent(in) :: functionRNLTE

   real(C_DOUBLE),   intent(in) :: tfloor
   real(C_DOUBLE),   intent(in) :: radForceMultiplier
   real(C_DOUBLE),   intent(in) :: betaNLTE
   real(C_DOUBLE),   intent(in) :: gammaNLTE

   logical(C_BOOL),  intent(in) :: DopplerShiftOn
   logical(C_BOOL),  intent(in) :: useNewNonLinearSolver
   logical(C_BOOL),  intent(in) :: useNewGTASolver
   logical(C_BOOL),  intent(in) :: usePWLD
   logical(C_BOOL),  intent(in) :: useSurfaceMassLumping

   Size% functionRNLTE         = functionRNLTE
   Size% tfloor                = tfloor
   Size% radForceMultiplier    = radForceMultiplier
   Size% betaNLTE              = betaNLTE
   Size% gammaNLTE             = gammaNLTE
   Size% DopplerShiftOn        = DopplerShiftOn

!  Need to phase these out of the interface.  Changing these during runtime are
!  not supported, they must be set at initial start time only via Size
!  constructor
!   Size% useNewNonLinearSolver = useNewNonLinearSolver
!   Size% useNewGTASolver       = useNewGTASolver
!   Size% usePWLD               = usePWLD
!   Size% useSurfaceMassLumping = useSurfaceMassLumping

   return
   end subroutine ResetSize
