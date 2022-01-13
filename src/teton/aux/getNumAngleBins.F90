!***********************************************************************
!                        Version 1:  04/20, BCY                        *
!                                                                      *
!   getNumAngleBins - Return the number of polar angle bins            *
!                                                                      *
!***********************************************************************
   subroutine getNumAngleBins(numAngleBins) &
                BIND(C,NAME="teton_getnumanglebins")

   USE ISO_C_BINDING
   use AngleSet_mod
   use QuadratureList_mod

   implicit none 

!  Arguments

   integer(C_INT), intent(out)   :: numAngleBins

!  Local
   type(AngleSet), pointer  :: ASet
   integer                  :: setID

!  There is only one higher-order SN set for all energy groups:
   setID = 1

   ASet   => getAngleSetFromSetID(Quad, setID)
   numAngleBins = ASet% nPolarAngles

   return
   end subroutine getNumAngleBins
