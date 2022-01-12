#include "macros.h"

!***********************************************************************
!                        Version 1:  04/20, BCY                        *
!                                                                      *
!   getAngleBins - Return polar angle bins                             *
!                                                                      *
!***********************************************************************
   subroutine getAngleBins(numAngleBins,       &
                           angleBinBoundaries) &
                 BIND(C,NAME="teton_getanglebins")

   USE ISO_C_BINDING
   use AngleSet_mod
   use QuadratureList_mod
   use kind_mod
   use Size_mod

   implicit none 

!  Arguments

   integer(C_INT), intent(in)   :: numAngleBins
   real(C_DOUBLE), intent(out)  :: angleBinBoundaries(numAngleBins +1)

!  Local

   type(AngleSet), pointer  :: ASet
   integer                  :: setID
   integer                  :: NumAngles
   integer                  :: angle
   integer                  :: polarAngle
   real(adqt)               :: weight
   integer                  :: bin

!  There is only one higher-order SN set for all energy groups:
   setID = 1

   ASet   => getAngleSetFromSetID(Quad, setID)

   tetonAssert(ASet% nPolarAngles == numAngleBins, &
    "numAngleBins must be # of Teton polar angles in teton_getanglebins")

   angleBinBoundaries(:) = 0.0

   NumAngles = ASet% NumAngles

   do angle=1,NumAngles
      polarAngle = ASet% polarAngle(angle)
      weight     = ASet% Weight(angle)
      angleBinBoundaries(polarAngle+1) = angleBinBoundaries(polarAngle+1) + weight
   enddo

   tetonAssert(abs(sum(angleBinBoundaries) - one) < 1.e-14_adqt, &
               "Error in getA_adqtngleBins.F90: sum(angle weights) != 1")

   angleBinBoundaries(:) = angleBinBoundaries(:)*Size%wtiso*2
   angleBinBoundaries(1) = -1.0

! Integrate the weights to get bin boundaries:
   do bin=2,numAngleBins+1
      angleBinBoundaries(bin) = angleBinBoundaries(bin)+angleBinBoundaries(bin-1)
   enddo

   tetonAssert((angleBinBoundaries(numAngleBins+1) - one) < 1.e-14_adqt, &
               "Error in getAngleBins.F90: Last bin boundary != 1")

   return
   end subroutine getAngleBins
