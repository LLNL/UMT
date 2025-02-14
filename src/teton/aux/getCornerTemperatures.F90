!***********************************************************************
!                        Last Update:  08/2024, BCY                    *
!                                                                      *
!   Copies Mat%Tec to Tec.                                             *
!                                                                      *
!***********************************************************************
   subroutine getCornerTemperatures(Tec) BIND(C,NAME="teton_getcornertemperatures")

   USE ISO_C_BINDING

   use Material_mod
   use Size_mod

   implicit none

!  Arguments

   real(C_DOUBLE), intent(out)     :: Tec(Size%ncornr)

   Tec(:) = Mat%Tec(:)

   return
   end subroutine getCornerTemperatures

