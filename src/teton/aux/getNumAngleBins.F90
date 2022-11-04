!***********************************************************************
!                        Version 1:  04/20, BCY                        *
!                                                                      *
!   getNumAngleBins - Return the number of polar angle bins            *
!                                                                      *
!***********************************************************************
   subroutine getNumAngleBins(numAngleBins) &
                BIND(C,NAME="teton_getnumanglebins")

   USE ISO_C_BINDING
   use QuadratureList_mod

   implicit none 

!  Arguments
   integer(C_INT), intent(out)   :: numAngleBins

!  Local
   integer                  :: qset

!  There is only one higher-order SN quadrature set for all energy groups:
   qset = 1

   numAngleBins = Quad%QuadPtr(qset)%nPolarAngles

   return
   end subroutine getNumAngleBins
