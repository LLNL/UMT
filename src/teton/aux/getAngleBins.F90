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
   use constant_mod
   use Size_mod

   implicit none 

!  Arguments

   integer(C_INT), intent(in)   :: numAngleBins
   real(C_DOUBLE), intent(out)  :: angleBinBoundaries(numAngleBins +1)

!  Local

   type(AngleSet), pointer  :: ASet
   integer                  :: angleSetID
   integer                  :: NumAngles
   integer                  :: numAngleSets
   integer                  :: angle
   integer                  :: polarAngle
   real(adqt)               :: weight
   integer                  :: bin

   angleBinBoundaries(:) = zero

   numAngleSets = getNumberOfAngleSets(Quad)

   do angleSetID=1,numAngleSets

      ASet   => getAngleSetData(Quad, angleSetID)

      TETON_ASSERT(ASet% nPolarAngles == numAngleBins, "numAngleBins must be # of Teton polar angles in teton_getanglebins")

      NumAngles = ASet% NumAngles

      do angle=1,NumAngles
         polarAngle = ASet% polarAngle(angle)
         weight     = ASet% Weight(angle)
         angleBinBoundaries(polarAngle+1) = angleBinBoundaries(polarAngle+1) + weight
      enddo

   enddo

   TETON_ASSERT(abs(sum(angleBinBoundaries)*Size%wtiso - one) < 1.e-12_adqt, "Error in getAngleBins.F90: sum(angle weights) != 1")

   angleBinBoundaries(:) = angleBinBoundaries(:)*Size%wtiso*2
   angleBinBoundaries(1) = -one

! Integrate the weights to get bin boundaries:
   do bin=2,numAngleBins+1
      angleBinBoundaries(bin) = angleBinBoundaries(bin)+angleBinBoundaries(bin-1)
   enddo

! Check that the final weight is close enough to one:
   TETON_ASSERT((angleBinBoundaries(numAngleBins+1) - one) < 1.e-14_adqt, "Error in getAngleBins.F90: Last bin boundary != 1")

! Then set it to one explicitly:
   angleBinBoundaries(numAngleBins+1) = one

   return
   end subroutine getAngleBins
