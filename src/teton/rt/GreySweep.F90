#include "macros.h"
#include "omp_wrappers.h"

!***********************************************************************
!                       Last Update:  10/2016, PFN                     *
!                                                                      *
!   GREYSWEEP - Performs a grey sweep with a fixed residual with       *
!               stretched TSA preconditioner. This gives the "action"  *
!               of the grey-transport operator on the input vector P.  *
!                                                                      *
!***********************************************************************

   subroutine GreySweepNEW(PsiB, P, withSource)

   use kind_mod
   use constant_mod
   use Size_mod
   use GreyAcceleration_mod

   implicit none

!  Arguments

   real(adqt), intent(inout)    :: PsiB(Size% nSurfElem,Size%nangGTA) 
   real(adqt), intent(inout)    :: P(Size%ncornr)

   logical (kind=1), intent(in) :: withSource

!  Local

   integer            :: zone 

   logical (kind=1)   :: useGPU


!  Set some constants

   useGPU = getGPUStatus(Size)

!  Perform a transport sweep to update the grey corrections, P

   if ( useGPU ) then

     TOMP_MAP(target enter data map(to: P))

     call GTASweep_GPU(P, PsiB, withSource)

     call ScalarIntensitySolve_GPU(P)

     TOMP_MAP(target exit data map(from: P))

   else

     GTA%ID = 1

     call GTASweep(P, PsiB)

!$omp parallel do default(none) schedule(static) &
!$omp& shared(Size, P, withSource)
     do zone=1,Size% nzones
       call ScalarIntensityDecompose(zone, P, withSource)
     enddo
!$omp end parallel do

   endif


   return
   end subroutine GreySweepNEW 


   subroutine GreySweep(PsiB, P)

   use kind_mod
   use constant_mod
   use Size_mod
   use GreyAcceleration_mod
   use QuadratureList_mod
   use SetData_mod

   implicit none

!  Arguments

   real(adqt),    intent(inout) :: PsiB(Size%nbelem,Size%nangGTA)
   real(adqt),    intent(inout) :: p(Size%ncornr)

!  Local

   type(SetData), pointer       :: Set

   integer                      :: GTAsetID
   integer                      :: nGTASets

   nGTASets = getNumberOfGTASets(Quad)

!  Perform a transport sweep to update the grey corrections, P
!
!  BCY: This sets P --> Linv*(S*P + GTA%GreySource)
!  Note that GTA%GreySource is zeroed out in GTASolver after the first sweep.

   GTA%OldGreyCorrection(:) = P(:)
   GTA%ID = 1

   call GTASweep(P, PsiB)

   P(:) = zero

   do GTAsetID=1,nGTASets
     Set  => getGTASetData(Quad, GTAsetID)
     P(:) =  P(:) + Set% tPhi(:)
   enddo

!  Now use "stretched" TSA as a preconditioner. 
!  The source to the TSA equations is the residual
!  from the grey transport sweep. The stretching parameter epsilon
!  has been chosen such that the scattering source vanishes.

   GTA%TsaPsib(:,:) = zero

!  Perform a "stretched" transport sweep

   GTA%ID = 2

!  BCY: With GTA%ID = 2, GTASweep applies Linv_{eps}*S_{eps} to (P-GTA%OldGreyCorrection)
   call GTASweep(P, GTA%TsaPsib)

!  Add the "stretched" TSA correction

   do GTAsetID=1,nGTASets
     Set  => getGTASetData(Quad, GTAsetID)
     P(:) =  P(:) + Set% tPhi(:)
   enddo
!  BCY: The final result is P --> Linv_{eps}S_{eps}*(Linv*(S*P + GreySource) - P) + Linv*(S*P + GreySource)
!   or equivalently,
!      P --> ( I - (I + Linv_{eps}S_{eps})(I-Linv S) ) P + (I+Linv_{eps}S_{eps})*Linv GreySource
!    When GreySource is zero, that simplifies to
!      P --> Linv_{eps}S_{eps}(Linv*S*P - P) + Linv*S*P

   return
   end subroutine GreySweep

