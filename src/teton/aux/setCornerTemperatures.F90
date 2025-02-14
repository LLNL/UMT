!***********************************************************************
!                        Last Update:  08/2024, BCY                    *
!                                                                      *
!   Sets Mat%Tec to Tec, requires initialized Mat and Size modules.    *
!                                                                      *
!***********************************************************************
   subroutine setCornerTemperatures(Tec) BIND(C,NAME="teton_setcornertemperatures")

   USE ISO_C_BINDING

   use Material_mod
   use Size_mod
   use Geometry_mod

   implicit none

!  Arguments

   real(C_DOUBLE), intent(in)     :: Tec(Size%ncornr)

   Mat%Tec(:) = Tec(:)

   return
   end subroutine setCornerTemperatures

