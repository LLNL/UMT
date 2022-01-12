!***********************************************************************
!                        Last Update:  12/2016, PFN                    *
!                                                                      *
!   resetTimers    -  Called after restart to sync host code timers    * 
!                     with the internal ones.                          *
!                                                                      *
!***********************************************************************

   subroutine resetTimers(RadtrTimeTotal, MatCoupTimeTotal,  &
                          SweepTimeTotal, GPUSweepTimeTotal, &
                          GTATimeTotal,                      &
                          InitTimeTotal, FinalTimeTotal)     &
                          BIND(C,NAME="teton_resettimers")

!  Include
   USE ISO_C_BINDING
   use kind_mod
   use Size_mod


   implicit none

!  Arguments

   real(C_DOUBLE),        intent(in) :: RadtrTimeTotal 
   real(C_DOUBLE),        intent(in) :: MatCoupTimeTotal 
   real(C_DOUBLE),        intent(in) :: SweepTimeTotal
   real(C_DOUBLE),        intent(in) :: GPUSweepTimeTotal
   real(C_DOUBLE),        intent(in) :: GTATimeTotal
   real(C_DOUBLE),        intent(in) :: InitTimeTotal
   real(C_DOUBLE),        intent(in) :: FinalTimeTotal

!  Set internal timers    

   Size% RadtrTimeTotal    = RadtrTimeTotal
   Size% MatCoupTimeTotal  = MatCoupTimeTotal
   Size% SweepTimeTotal    = SweepTimeTotal
   Size% GPUSweepTimeTotal = GPUSweepTimeTotal
   Size% GTATimeTotal      = GTATimeTotal
   Size% InitTimeTotal     = InitTimeTotal
   Size% FinalTimeTotal    = FinalTimeTotal
 

   return
   end subroutine resetTimers

